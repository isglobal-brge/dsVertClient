#' @title K>=3 Ring63 DCF + Beaver Gradient Loop
#' @description Pure Ring63 protocol for K>=3 binomial/Poisson GLM.
#'   Uses additive secret sharing of features to 2 DCF parties,
#'   DCF wide spline for sigmoid/exp, and Beaver matvec for gradient.
#'   Pure Ring63 additive secret sharing.
#'
#' @details Protocol:
#' \enumerate{
#'   \item One-time: all servers share features+y with 2 DCF parties (Ring63 FP)
#'   \item Per iteration: DCF parties compute eta from X_share*beta
#'   \item DCF wide spline between DCF parties -> mu shares
#'   \item Beaver matvec between DCF parties -> gradient shares
#'   \item Client aggregates Ring63 shares -> plaintext gradient
#'   \item Client L-BFGS update
#' }
#'
#' Requires: transport keys only. No CPK, Galois, or RLK.
#' Precision: ~1e-3 (DCF approximation).
#' @name k3-ring63-dcf-loop
NULL

#' @keywords internal
.k3_select_fusion_server <- function(server_list, coordinator, x_vars) {
  non_coord <- setdiff(server_list, coordinator)
  non_coord[which.max(sapply(non_coord, function(s) length(x_vars[[s]])))]
}

#' @keywords internal
.k3_ring63_gradient_loop <- function(
    datasources, server_list, server_names, x_vars,
    coordinator, coordinator_conn, non_label_servers,
    transport_pks, std_data, y_var, family,
    betas, n_obs, lambda,
    session_id, max_iter, tol, verbose,
    label_intercept, .dsAgg, .sendBlob,
    weights_active = FALSE,
    no_intercept = FALSE,
    start = NULL,
    compute_se = TRUE) {

  K <- length(server_list)
  # DCF parties: fusion (party 0) must differ from coordinator (party 1)
  # Pick the first non-coordinator server with the most features as fusion
  non_coord <- setdiff(server_list, coordinator)
  fusion_server <- .k3_select_fusion_server(server_list, coordinator, x_vars)
  fusion_conn <- which(server_names == fusion_server)
  dcf_parties <- c(fusion_server, coordinator)
  dcf_conns <- sapply(dcf_parties, function(s) which(server_names == s))

  frac_bits <- 20L
  is_gaussian <- (family == "gaussian")
  num_intervals <- if (family == "poisson") 100L else 50L
  dcf_family <- if (family == "poisson") "poisson" else "sigmoid"

  .to_b64url <- function(x) gsub("+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  # Safe beta extraction: returns NULL for empty index (avoids numeric(0) in Opal calls)
  .bsafe <- function(b, idx) { if (is.null(idx) || length(idx) == 0) NULL else b[idx] }
  .b64url_to_b64 <- function(x) {
    x <- gsub("-","+",gsub("_","/",x,fixed=TRUE),fixed=TRUE)
    pad <- nchar(x)%%4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }

  if (verbose) message(sprintf(
    "\n[Phase 3] K>=%d Ring63 DCF + Beaver L-BFGS (family=%s, n=%d, intervals=%d, lambda=%.1e)",
    K, family, n_obs, num_intervals, lambda))

  # ===========================================================================
  # Step A: Input sharing -- all servers share features with 2 DCF parties
  # ===========================================================================
  t0_share <- proc.time()[[3]]
  if (verbose) message("  [Input Sharing] Sharing features+y with DCF parties...")

  # Total features across all servers (EXCLUDING intercept -- handled separately)
  p_total <- sum(sapply(server_list, function(s) length(x_vars[[s]])))

  # Each server: split features into 2 shares, send to DCF parties
  for (server in server_list) {
    ci <- which(server_names == server)
    for (di in seq_along(dcf_parties)) {
      dcf_srv <- dcf_parties[di]
      dcf_ci <- dcf_conns[di]
      dcf_pk <- transport_pks[[dcf_srv]]

      if (server == dcf_srv) {
        # This server IS a DCF party -- just call k2ShareInputDS targeting the other DCF party
        peer_dcf <- dcf_parties[3 - di]
        peer_pk_safe <- .to_b64url(transport_pks[[peer_dcf]])
        srv_x <- x_vars[[server]]; if (length(srv_x) == 0) srv_x <- NULL
        r <- .dsAgg(datasources[ci], call(name = "k2ShareInputDS",
          data_name = std_data, x_vars = srv_x,
          y_var = if (server == coordinator) y_var else NULL,
          peer_pk = peer_pk_safe, session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]
        # Send encrypted shares to the other DCF party
        peer_ci <- dcf_conns[3 - di]
        .sendBlob(r$encrypted_x_share, "k2_peer_x_share", peer_ci)
        if (!is.null(r$encrypted_y_share))
          .sendBlob(r$encrypted_y_share, "k2_peer_y_share", peer_ci)
        break  # Only share once (with the other DCF party)
      }
    }
  }

  # Non-DCF servers: split features into 2 additive shares, one per DCF party.
  # Each server calls k2ShareInputDS ONCE (targeting fusion), then exports the
  # complement (own_share) to coordinator. This ensures share_0 + share_1 = X_k.
  non_dcf_servers <- setdiff(server_list, dcf_parties)
  for (server in non_dcf_servers) {
    ci <- which(server_names == server)
    # Split once targeting fusion as peer
    fusion_pk_safe <- .to_b64url(transport_pks[[fusion_server]])
    srv_x <- x_vars[[server]]; if (length(srv_x) == 0) srv_x <- NULL
    r <- .dsAgg(datasources[ci], call(name = "k2ShareInputDS",
      data_name = std_data, x_vars = srv_x,
      y_var = NULL,
      peer_pk = fusion_pk_safe, session_id = session_id))
    if (is.list(r) && length(r) == 1) r <- r[[1]]
    # Send peer_share to fusion
    .sendBlob(r$encrypted_x_share, paste0("k2_extra_x_share_", server), fusion_conn)

    # Export own_share (the complement) to coordinator
    coord_pk_safe <- .to_b64url(transport_pks[[coordinator]])
    r2 <- .dsAgg(datasources[ci], call(name = "glmRing63ExportOwnShareDS",
      peer_pk = coord_pk_safe, session_id = session_id))
    if (is.list(r2) && length(r2) == 1) r2 <- r2[[1]]
    .sendBlob(r2$encrypted_own_share, paste0("k2_extra_x_share_", server), coordinator_conn)
  }

  # DCF parties: receive peer's shares
  for (di in seq_along(dcf_parties)) {
    dcf_srv <- dcf_parties[di]
    dcf_ci <- dcf_conns[di]
    peer_p <- length(x_vars[[dcf_parties[3 - di]]])  # no intercept in features
    .dsAgg(datasources[dcf_ci], call(name = "k2ReceiveShareDS",
      peer_p = as.integer(peer_p), session_id = session_id))
  }

  # DCF parties: receive non-DCF feature shares
  for (server in non_dcf_servers) {
    for (di in seq_along(dcf_parties)) {
      dcf_ci <- dcf_conns[di]
      # Receive and store the extra feature shares
      # These will be assembled into the full X matrix
      .dsAgg(datasources[dcf_ci], call(name = "glmRing63ReceiveExtraShareDS",
        extra_key = paste0("k2_extra_x_share_", server),
        extra_p = as.integer(length(x_vars[[server]])),
        session_id = session_id))
    }
  }

  if (verbose) message(sprintf("  [Input Sharing] Complete (p_total=%d, %.1fs)",
    p_total, proc.time()[[3]] - t0_share))

  # ===========================================================================
  # Step B: Pre-generate DCF keys (on NON-DCF server, not client)
  #   Skipped for Gaussian (identity link, no sigmoid/exp needed)
  # ===========================================================================
  # Dealer rotation: different non-DCF server generates triples each iteration.
  # Analyst must compromise ALL dealers + 1 DCF party to extract data.
  # For K servers: security threshold = (K-1)/K.
  all_dealers <- if (length(non_dcf_servers) > 0) non_dcf_servers else c(fusion_server)
  dealer <- all_dealers[1]  # initial dealer for DCF keys
  dealer_conn <- which(server_names == dealer)

  if (!is_gaussian) {
    t0_dcf <- proc.time()[[3]]
    if (verbose) message(sprintf("  [DCF] Server %s generating keys (n=%d, %d intervals)...",
                                   dealer, n_obs, num_intervals))

    dcf_result <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenDcfKeysDS",
           dcf0_pk = transport_pks[[dcf_parties[1]]],
           dcf1_pk = transport_pks[[dcf_parties[2]]],
           family = dcf_family, n = as.integer(n_obs),
           frac_bits = frac_bits, num_intervals = num_intervals,
           session_id = session_id))
    if (is.list(dcf_result)) dcf_result <- dcf_result[[1]]

    .sendBlob(dcf_result$dcf_blob_0, "k2_dcf_keys_persistent", dcf_conns[1])
    .dsAgg(datasources[dcf_conns[1]], call(name = "k2StoreDcfKeysPersistentDS", session_id = session_id))
    .sendBlob(dcf_result$dcf_blob_1, "k2_dcf_keys_persistent", dcf_conns[2])
    .dsAgg(datasources[dcf_conns[2]], call(name = "k2StoreDcfKeysPersistentDS", session_id = session_id))

    if (verbose) message(sprintf("  [DCF] Keys distributed (%.1fs)", proc.time()[[3]] - t0_dcf))
  }

  # ===========================================================================
  # Iteration loop
  # ===========================================================================
  beta <- rep(0, p_total)
  intercept <- 0.0
  feature_order <- c(x_vars[[coordinator]])
  for (srv in server_list) {
    if (srv == coordinator) next
    feature_order <- c(feature_order, x_vars[[srv]])
  }
  if (!is.null(start)) {
    start_vec <- as.numeric(start)
    names(start_vec) <- names(start)
    if (!is.null(names(start_vec)) && "(Intercept)" %in% names(start_vec)) {
      intercept <- as.numeric(start_vec["(Intercept)"])
    }
    if (!is.null(names(start_vec))) {
      miss <- setdiff(feature_order, names(start_vec))
      if (length(miss) > 0L) {
        stop("start is missing coefficient(s): ",
             paste(miss, collapse = ", "), call. = FALSE)
      }
      beta <- as.numeric(start_vec[feature_order])
    } else if (length(start_vec) == p_total + 1L) {
      intercept <- start_vec[1L]
      beta <- start_vec[-1L]
    } else if (length(start_vec) == p_total) {
      beta <- start_vec
    } else {
      stop("start length must be p or p+1", call. = FALSE)
    }
    if (isTRUE(no_intercept)) intercept <- 0
  }
  lbfgs_s <- list(); lbfgs_y <- list()
  prev_theta <- NULL; prev_grad <- NULL
  converged <- FALSE; final_iter <- 0
  max_diff <- Inf

  # Feature counts for beta splitting (EXCLUDING intercept)
  p_coord <- length(x_vars[[coordinator]])
  p_fusion <- length(x_vars[[fusion_server]])
  p_extras <- p_total - p_coord - p_fusion
  p_nl <- p_total - p_coord

  t0_loop <- proc.time()[[3]]

  if (max_iter < 1L) {
    converged <- TRUE
    final_iter <- 0L
    max_diff <- 0
  }

  for (iter in seq_len(max_iter)) {
    t0_iter <- proc.time()[[3]]
    beta_old <- beta
    intercept_old <- intercept

    # Rotate dealer: each non-DCF server takes turns generating triples.
    # Analyst must control ALL dealers to extract data from ALL iterations.
    dealer <- all_dealers[((iter - 1L) %% length(all_dealers)) + 1L]
    dealer_conn <- which(server_names == dealer)

    # =================================================================
    # Step 1: Compute eta shares (Ring63 FP, both DCF parties)
    # =================================================================
    if (verbose && iter <= 3) message(sprintf("    [%d.1] Computing eta shares...", iter))

    # Build beta for each DCF party in the order matching X_full layout.
    # K=2 convention: X_full = [own|peer] (coordinator) or [peer|own] (non-coord)
    # Go command reorders: coord features first when is_party_zero=TRUE (coordinator),
    #                       peer first when is_party_zero=FALSE (fusion).
    # Fusion: X = [s3_peer | s2_extra | s1_own]
    # Coord:  X = [s3_own  | s1_peer  | s2_extra]
    # beta_coord = s3 betas, beta_nl = remaining in X_full order

    # Build a map: server -> beta indices
    beta_map <- list()
    idx <- 1
    if (p_coord > 0) { beta_map[[coordinator]] <- idx:(idx + p_coord - 1); idx <- idx + p_coord }
    for (srv in server_list) {
      if (srv == coordinator) next
      p_s <- length(x_vars[[srv]])
      if (p_s == 0) next
      beta_map[[srv]] <- idx:(idx + p_s - 1); idx <- idx + p_s
    }

    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]
      srv <- dcf_parties[i]
      is_coord <- (srv == coordinator)
      peer_srv <- dcf_parties[3 - i]

      if (is_coord) {
        # Coordinator: X = [own(coord) | peer(fusion) | extra...]
        # beta_coord = own betas, beta_nl = [fusion betas | extra betas]
        b_coord <- .bsafe(beta, beta_map[[coordinator]])
        b_nl <- c(.bsafe(beta, beta_map[[peer_srv]]))
        for (ns in non_dcf_servers) b_nl <- c(b_nl, .bsafe(beta, beta_map[[ns]]))
      } else {
        # Fusion: X = [peer(coord) | extra... | own(fusion)]
        # The Go cmd with is_party_zero=FALSE puts: [peer | own]
        # peer = coord_share + extras, own = fusion_share
        # So beta_coord = coord betas, beta_nl = [extra betas | own betas]
        b_coord <- .bsafe(beta, beta_map[[coordinator]])
        b_nl <- c()
        for (ns in non_dcf_servers) b_nl <- c(b_nl, .bsafe(beta, beta_map[[ns]]))
        b_nl <- c(b_nl, .bsafe(beta, beta_map[[srv]]))  # own (fusion) last
      }

      .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
        beta_coord = b_coord, beta_nl = b_nl,
        intercept = if (is_coord) intercept else 0.0,
        is_coordinator = is_coord, session_id = session_id))

      # Reorder fusion's X_full to canonical order [coord|fusion|extras]
      # so Beaver gradient works (both parties must have same column layout)
      if (!is_coord && p_extras > 0) {
        .dsAgg(datasources[ci], call(name = "glmRing63ReorderXFullDS",
          p_coord = as.integer(p_coord), p_fusion = as.integer(p_fusion),
          p_extras = as.integer(p_extras), session_id = session_id))
      }
    }

    # =================================================================
    # Step 2: Link function
    #   Gaussian: identity (mu = eta, no DCF needed)
    #   Binomial/Poisson: DCF wide spline (4-phase between DCF parties)
    # =================================================================
    if (is_gaussian) {
      # Identity link: mu_share = eta_share (stored as secure_mu_share for gradient)
      for (i in seq_along(dcf_parties)) {
        ci <- dcf_conns[i]
        .dsAgg(datasources[ci], call(name = "k2IdentityLinkDS", session_id = session_id))
      }
    } else {
      if (verbose && iter <= 3) message(sprintf("    [%d.2] DCF wide spline...", iter))

      # Generate spline Beaver triples on dealer server (not client)
      spline_t <- .dsAgg(datasources[dealer_conn],
        call(name = "glmRing63GenSplineTriplesDS",
             dcf0_pk = transport_pks[[dcf_parties[1]]],
             dcf1_pk = transport_pks[[dcf_parties[2]]],
             n = as.integer(n_obs), frac_bits = frac_bits,
             session_id = session_id))
      if (is.list(spline_t)) spline_t <- spline_t[[1]]
      .sendBlob(spline_t$spline_blob_0, "k2_spline_triples", dcf_conns[1])
      .sendBlob(spline_t$spline_blob_1, "k2_spline_triples", dcf_conns[2])

      # Phase 1-4
      ph1 <- list()
      for (i in seq_along(dcf_parties)) {
        ci <- dcf_conns[i]
        r <- .dsAgg(datasources[ci], call(name = "k2WideSplinePhase1DS",
          party_id = as.integer(i - 1), family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; ph1[[i]] <- r
      }
      .sendBlob(ph1[[1]]$dcf_masked, "k2_peer_dcf_masked", dcf_conns[2])
      .sendBlob(ph1[[2]]$dcf_masked, "k2_peer_dcf_masked", dcf_conns[1])

      ph2 <- list()
      for (i in seq_along(dcf_parties)) {
        ci <- dcf_conns[i]
        r <- .dsAgg(datasources[ci], call(name = "k2WideSplinePhase2DS",
          party_id = as.integer(i - 1), family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; ph2[[i]] <- r
      }
      for (i in seq_along(dcf_parties)) {
        peer_idx <- 3 - i; peer_ci <- dcf_conns[peer_idx]
        pk_b64 <- .b64url_to_b64(transport_pks[[dcf_parties[peer_idx]]])
        r1_json <- jsonlite::toJSON(list(
          and_xma = ph2[[i]]$and_xma, and_ymb = ph2[[i]]$and_ymb,
          had1_xma = ph2[[i]]$had1_xma, had1_ymb = ph2[[i]]$had1_ymb),
          auto_unbox = TRUE)
        sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
          data = jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk = pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", peer_ci)
      }

      ph3 <- list()
      for (i in seq_along(dcf_parties)) {
        ci <- dcf_conns[i]
        r <- .dsAgg(datasources[ci], call(name = "k2WideSplinePhase3DS",
          party_id = as.integer(i - 1), family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; ph3[[i]] <- r
      }
      for (i in seq_along(dcf_parties)) {
        peer_idx <- 3 - i; peer_ci <- dcf_conns[peer_idx]
        pk_b64 <- .b64url_to_b64(transport_pks[[dcf_parties[peer_idx]]])
        r1_json <- jsonlite::toJSON(list(
          had2_xma = ph3[[i]]$had2_xma, had2_ymb = ph3[[i]]$had2_ymb),
          auto_unbox = TRUE)
        sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
          data = jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk = pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", peer_ci)
      }

      for (i in seq_along(dcf_parties)) {
        ci <- dcf_conns[i]
        .dsAgg(datasources[ci], call(name = "k2WideSplinePhase4DS",
          party_id = as.integer(i - 1), family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          session_id = session_id))
      }
    }

    # Optional weighted GLM/IPW path. The weights vector is shared between the
    # two DCF parties; Beaver computes shares of w * (mu - y).
    if (isTRUE(weights_active)) {
      .glm_apply_shared_weight_residual(
        datasources = datasources,
        dcf_parties = dcf_parties,
        dcf_conns = dcf_conns,
        dealer_conn = dealer_conn,
        transport_pks = transport_pks,
        session_id = session_id,
        n_obs = n_obs,
        .dsAgg = .dsAgg,
        .sendBlob = .sendBlob,
        ring = 63L)
    }

    # =================================================================
    # Step 3: Beaver gradient (Ring63)
    # =================================================================
    if (verbose && iter <= 3) message(sprintf("    [%d.3] Beaver gradient...", iter))

    # Generate gradient Beaver triples on dealer server (not client)
    grad_t <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenGradTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[1]]],
           dcf1_pk = transport_pks[[dcf_parties[2]]],
           n = as.integer(n_obs), p = as.integer(p_total),
           session_id = session_id))
    if (is.list(grad_t)) grad_t <- grad_t[[1]]
    # Relay opaque blobs (client can't read)
    .sendBlob(grad_t$grad_blob_0, "k2_grad_triple_fp", dcf_conns[1])
    .sendBlob(grad_t$grad_blob_1, "k2_grad_triple_fp", dcf_conns[2])

    r1_results <- list()
    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]
      peer <- dcf_parties[3 - i]
      .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      r1_results[[i]] <- r
    }
    .sendBlob(r1_results[[1]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[2])
    .sendBlob(r1_results[[2]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[1])

    grad_results <- list()
    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR2DS",
        party_id = as.integer(i - 1), session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      grad_results[[i]] <- r
    }

    # =================================================================
    # Step 4: Aggregate gradient + L-BFGS update
    # =================================================================
    # Both DCF parties now have X_full in canonical order [coord|fusion|extras]
    # (fusion was reordered by glmRing63ReorderXFullDS)
    # Simple Ring63 aggregation -- no permutation needed.
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = grad_results[[1]]$gradient_fp,
      share_b = grad_results[[2]]$gradient_fp,
      frac_bits = frac_bits))
    gradient_canonical <- agg$values  # in canonical order [coord|fusion|extras]

    agg_res <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r1_results[[1]]$sum_residual_fp,
      share_b = r1_results[[2]]$sum_residual_fp,
      frac_bits = frac_bits))
    sum_residual <- agg_res$values[1]

    # Remap from canonical [coord | fusion | extras] to beta order [coord | non-coord-by-server-list]
    gradient <- numeric(p_total)
    gi <- 1
    if (p_coord > 0) { gradient[beta_map[[coordinator]]] <- gradient_canonical[gi:(gi+p_coord-1)]; gi <- gi + p_coord }
    if (p_fusion > 0) { gradient[beta_map[[fusion_server]]] <- gradient_canonical[gi:(gi+p_fusion-1)]; gi <- gi + p_fusion }
    for (ns in non_dcf_servers) {
      p_ns <- length(x_vars[[ns]])
      if (p_ns == 0) next
      gradient[beta_map[[ns]]] <- gradient_canonical[gi:(gi + p_ns - 1)]
      gi <- gi + p_ns
    }

    theta <- c(intercept, beta)
    full_grad <- c(sum_residual / n_obs, gradient / n_obs) + lambda * theta
    if (isTRUE(no_intercept)) {
      full_grad[1] <- 0
    }

    if (!is.null(prev_theta)) {
      sk <- theta - prev_theta; yk <- full_grad - prev_grad
      if (sum(sk * yk) > 1e-10) {
        lbfgs_s <- c(lbfgs_s, list(sk)); lbfgs_y <- c(lbfgs_y, list(yk))
        if (length(lbfgs_s) > 7) { lbfgs_s <- lbfgs_s[-1]; lbfgs_y <- lbfgs_y[-1] }
      }
    }
    prev_theta <- theta; prev_grad <- full_grad

    direction <- .lbfgs_direction(full_grad, lbfgs_s, lbfgs_y)
    # Step bound: cap |direction| at 5.0 to prevent first-iter explosion
    # (Cox K=3 Poisson trick with offset + 100 baseline dummies can produce
    # an enormous L-BFGS direction at iter 1, which without damping drove
    # eta -> inf, exp(eta) -> Inf, and a NaN interval-index slice into Go-side
    # k2WideSplinePhase1DS resulting in SIGSEGV addr=0xefffffff90a0).
    d_max <- max(abs(direction))
    if (is.finite(d_max) && d_max > 5.0) direction <- direction * (5.0 / d_max)
    if (isTRUE(no_intercept)) {
      direction[1] <- 0
    }
    # Backtracking: halve step until new theta is finite and bounded
    step_size <- if (iter <= 1) 0.3 else 1.0
    for (bt in 0:6) {
      new_theta <- theta + step_size * direction
      if (all(is.finite(new_theta)) && max(abs(new_theta)) < 50) break
      step_size <- step_size * 0.5
    }
    if (!all(is.finite(new_theta)))
      stop("L-BFGS step produced non-finite theta \u2014 divergence guard", call. = FALSE)
    intercept <- if (isTRUE(no_intercept)) 0 else new_theta[1]
    beta <- new_theta[-1]

    max_diff <- max(abs(beta - beta_old), abs(intercept - intercept_old))
    final_iter <- iter
    grad_norm <- sqrt(sum(full_grad^2))

    if (verbose) message(sprintf("  Iter %d: ||grad||=%.4f, step=%.2f, diff=%.2e, theta=[%.3f, %.3f] (%.1fs)",
      iter, grad_norm, step_size, max_diff,
      min(new_theta), max(new_theta), proc.time()[[3]] - t0_iter))

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose) message(sprintf("  Converged after %d iterations (diff = %.2e)", iter, max_diff))
      break
    }
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)", max_iter, max_diff))

  if (verbose) message(sprintf("  [Ring63-DCF-Beaver] %s after %d iters (total %.1fs)",
    if (converged) "Converged" else "Stopped",
    final_iter, proc.time()[[3]] - t0_loop))

  # Build betas in per-server format
  betas_out <- list()
  if (p_coord > 0) {
    if (label_intercept) {
      betas_out[[coordinator]] <- c(intercept, beta[1:p_coord])
    } else {
      betas_out[[coordinator]] <- beta[1:p_coord]
    }
  } else if (label_intercept) {
    betas_out[[coordinator]] <- intercept
  }
  idx <- p_coord + 1
  for (server in server_list) {
    if (server == coordinator) next
    p_s <- length(x_vars[[server]])
    if (p_s == 0) next
    betas_out[[server]] <- beta[idx:(idx + p_s - 1)]
    idx <- idx + p_s
  }

  # ===========================================================================
  # Exact standard errors via finite-difference Hessian
  # Run p+1 additional gradient evaluations with small perturbations.
  # Each gives one column of the Hessian. Exact, uses existing Beaver.
  # Disclosure: only aggregate gradients (same as iterations). SAFE.
  # ===========================================================================
  inv_hessian <- list()
  if (isTRUE(compute_se)) {
  if (verbose) message("  [SE] Computing exact Hessian (central differences)...")
  p_plus1 <- p_total + 1
  delta <- 0.01  # larger delta for better signal vs Ring63 noise
  hessian <- matrix(0, p_plus1, p_plus1)
  theta <- c(intercept, beta)

  for (j in seq_len(p_plus1)) {
    # Central difference: H_j = (g(theta+deltae_j) - g(theta-deltae_j)) / (2delta)
    # More accurate: O(delta^2) error instead of O(delta)
    theta_pert <- theta
    theta_pert[j] <- theta[j] + delta
    int_pert <- theta_pert[1]
    beta_pert <- theta_pert[-1]

    # Rotate dealer for this SE evaluation too
    dealer <- all_dealers[((final_iter + j - 1L) %% length(all_dealers)) + 1L]
    dealer_conn <- which(server_names == dealer)

    # Compute eta with perturbed beta
    for (di in seq_along(dcf_parties)) {
      ci <- dcf_conns[di]
      srv <- dcf_parties[di]
      is_coord <- (srv == coordinator)
      if (is_coord) {
        b_coord <- .bsafe(beta_pert, beta_map[[coordinator]])
        b_nl <- c(beta_pert[beta_map[[dcf_parties[1]]]])
        for (ns in non_dcf_servers) b_nl <- c(b_nl, .bsafe(beta_pert, beta_map[[ns]]))
      } else {
        b_coord <- .bsafe(beta_pert, beta_map[[coordinator]])
        b_nl <- c()
        for (ns in non_dcf_servers) b_nl <- c(b_nl, .bsafe(beta_pert, beta_map[[ns]]))
        b_nl <- c(b_nl, .bsafe(beta_pert, beta_map[[srv]]))
      }
      .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
        beta_coord = b_coord, beta_nl = b_nl,
        intercept = if (is_coord) int_pert else 0.0,
        is_coordinator = is_coord, session_id = session_id))
      if (!is_coord && p_extras > 0) {
        .dsAgg(datasources[ci], call(name = "glmRing63ReorderXFullDS",
          p_coord = as.integer(p_coord), p_fusion = as.integer(p_fusion),
          p_extras = as.integer(p_extras), session_id = session_id))
      }
    }

    # Link function (identity for Gaussian, DCF for others)
    if (is_gaussian) {
      for (di in seq_along(dcf_parties)) {
        .dsAgg(datasources[dcf_conns[di]], call(name = "k2IdentityLinkDS", session_id = session_id))
      }
    } else {
      # DCF wide spline
      spline_t <- .dsAgg(datasources[dealer_conn],
        call(name = "glmRing63GenSplineTriplesDS",
             dcf0_pk = transport_pks[[dcf_parties[1]]],
             dcf1_pk = transport_pks[[dcf_parties[2]]],
             n = as.integer(n_obs), frac_bits = frac_bits,
             session_id = session_id))
      if (is.list(spline_t)) spline_t <- spline_t[[1]]
      .sendBlob(spline_t$spline_blob_0, "k2_spline_triples", dcf_conns[1])
      .sendBlob(spline_t$spline_blob_1, "k2_spline_triples", dcf_conns[2])
      for (ph in 1:4) {
        ph_r <- list()
        for (di in seq_along(dcf_parties)) {
          fn <- paste0("k2WideSplinePhase", ph, "DS")
          r <- .dsAgg(datasources[dcf_conns[di]], call(fn,
            party_id = as.integer(di - 1), family = family,
            num_intervals = num_intervals, frac_bits = frac_bits,
            session_id = session_id))
          if (is.list(r) && length(r) == 1) r <- r[[1]]; ph_r[[di]] <- r
        }
        if (ph == 1) {
          .sendBlob(ph_r[[1]]$dcf_masked, "k2_peer_dcf_masked", dcf_conns[2])
          .sendBlob(ph_r[[2]]$dcf_masked, "k2_peer_dcf_masked", dcf_conns[1])
        } else if (ph == 2) {
          for (di in 1:2) {
            pi2 <- 3 - di; pk <- .b64url_to_b64(transport_pks[[dcf_parties[pi2]]])
            sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
              data = jsonlite::base64_enc(charToRaw(jsonlite::toJSON(
                list(and_xma=ph_r[[di]]$and_xma, and_ymb=ph_r[[di]]$and_ymb,
                     had1_xma=ph_r[[di]]$had1_xma, had1_ymb=ph_r[[di]]$had1_ymb),
                auto_unbox=TRUE))), recipient_pk=pk))
            .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", dcf_conns[pi2])
          }
        } else if (ph == 3) {
          for (di in 1:2) {
            pi2 <- 3 - di; pk <- .b64url_to_b64(transport_pks[[dcf_parties[pi2]]])
            sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
              data = jsonlite::base64_enc(charToRaw(jsonlite::toJSON(
                list(had2_xma=ph_r[[di]]$had2_xma, had2_ymb=ph_r[[di]]$had2_ymb),
                auto_unbox=TRUE))), recipient_pk=pk))
            .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", dcf_conns[pi2])
          }
        }
      }
    }

    if (isTRUE(weights_active)) {
      .glm_apply_shared_weight_residual(
        datasources = datasources,
        dcf_parties = dcf_parties,
        dcf_conns = dcf_conns,
        dealer_conn = dealer_conn,
        transport_pks = transport_pks,
        session_id = session_id,
        n_obs = n_obs,
        .dsAgg = .dsAgg,
        .sendBlob = .sendBlob,
        ring = 63L)
    }

    # Beaver gradient with perturbed beta
    grad_t <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenGradTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[1]]],
           dcf1_pk = transport_pks[[dcf_parties[2]]],
           n = as.integer(n_obs), p = as.integer(p_total),
           session_id = session_id))
    if (is.list(grad_t)) grad_t <- grad_t[[1]]
    .sendBlob(grad_t$grad_blob_0, "k2_grad_triple_fp", dcf_conns[1])
    .sendBlob(grad_t$grad_blob_1, "k2_grad_triple_fp", dcf_conns[2])
    se_r1 <- list()
    for (di in seq_along(dcf_parties)) {
      ci <- dcf_conns[di]; peer <- dcf_parties[3 - di]
      .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; se_r1[[di]] <- r
    }
    .sendBlob(se_r1[[1]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[2])
    .sendBlob(se_r1[[2]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[1])
    se_r2 <- list()
    for (di in seq_along(dcf_parties)) {
      ci <- dcf_conns[di]
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR2DS",
        party_id = as.integer(di - 1), session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; se_r2[[di]] <- r
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = se_r2[[1]]$gradient_fp, share_b = se_r2[[2]]$gradient_fp,
      frac_bits = frac_bits))
    gradient_pert <- numeric(p_total)
    gi <- 1
    if (p_coord > 0) { gradient_pert[beta_map[[coordinator]]] <- agg$values[gi:(gi+p_coord-1)]; gi <- gi + p_coord }
    if (p_fusion > 0) { gradient_pert[beta_map[[fusion_server]]] <- agg$values[gi:(gi+p_fusion-1)]; gi <- gi + p_fusion }
    for (ns in non_dcf_servers) {
      pn <- length(x_vars[[ns]])
      if (pn == 0) next
      gradient_pert[beta_map[[ns]]] <- agg$values[gi:(gi+pn-1)]; gi <- gi + pn
    }
    agg_res <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = se_r1[[1]]$sum_residual_fp, share_b = se_r1[[2]]$sum_residual_fp,
      frac_bits = frac_bits))
    grad_pert_full <- c(agg_res$values[1] / n_obs, gradient_pert / n_obs) + lambda * theta_pert

    grad_forward <- grad_pert_full

    # Backward perturbation: theta - deltae_j (same code, negative delta)
    theta_back <- theta
    theta_back[j] <- theta[j] - delta
    int_back <- theta_back[1]; beta_back <- theta_back[-1]
    dealer_b <- all_dealers[((final_iter + p_plus1 + j - 1L) %% length(all_dealers)) + 1L]
    dealer_conn_b <- which(server_names == dealer_b)
    for (di in seq_along(dcf_parties)) {
      ci <- dcf_conns[di]; srv <- dcf_parties[di]; is_coord <- (srv == coordinator)
      if (is_coord) {
        bc <- beta_back[beta_map[[coordinator]]]; bnl <- c(beta_back[beta_map[[dcf_parties[1]]]])
        for (ns in non_dcf_servers) bnl <- c(bnl, beta_back[beta_map[[ns]]])
      } else {
        bc <- beta_back[beta_map[[coordinator]]]; bnl <- c()
        for (ns in non_dcf_servers) bnl <- c(bnl, beta_back[beta_map[[ns]]])
        bnl <- c(bnl, beta_back[beta_map[[srv]]])
      }
      .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS", beta_coord=bc, beta_nl=bnl,
        intercept=if(is_coord) int_back else 0, is_coordinator=is_coord, session_id=session_id))
      if (!is_coord && p_extras > 0)
        .dsAgg(datasources[ci], call(name = "glmRing63ReorderXFullDS",
          p_coord=as.integer(p_coord), p_fusion=as.integer(p_fusion),
          p_extras=as.integer(p_extras), session_id=session_id))
    }
    if (is_gaussian) {
      for (di in seq_along(dcf_parties))
        .dsAgg(datasources[dcf_conns[di]], call(name = "k2IdentityLinkDS", session_id=session_id))
    } else {
      st2 <- .dsAgg(datasources[dealer_conn_b], call(name = "glmRing63GenSplineTriplesDS",
        dcf0_pk=transport_pks[[dcf_parties[1]]], dcf1_pk=transport_pks[[dcf_parties[2]]],
        n=as.integer(n_obs), frac_bits=frac_bits, session_id=session_id))
      if (is.list(st2)) st2 <- st2[[1]]
      .sendBlob(st2$spline_blob_0, "k2_spline_triples", dcf_conns[1])
      .sendBlob(st2$spline_blob_1, "k2_spline_triples", dcf_conns[2])
      for (ph in 1:4) {
        pr <- list()
        for (di in 1:2) {
          r <- .dsAgg(datasources[dcf_conns[di]], call(paste0("k2WideSplinePhase",ph,"DS"),
            party_id=as.integer(di-1), family=family, num_intervals=num_intervals,
            frac_bits=frac_bits, session_id=session_id))
          if (is.list(r)&&length(r)==1) r<-r[[1]]; pr[[di]]<-r
        }
        if (ph==1) { .sendBlob(pr[[1]]$dcf_masked,"k2_peer_dcf_masked",dcf_conns[2]); .sendBlob(pr[[2]]$dcf_masked,"k2_peer_dcf_masked",dcf_conns[1]) }
        else if (ph==2) { for(di in 1:2){pi2<-3-di;pk<-.b64url_to_b64(transport_pks[[dcf_parties[pi2]]]);s<-dsVert:::.callMpcTool("transport-encrypt",list(data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(and_xma=pr[[di]]$and_xma,and_ymb=pr[[di]]$and_ymb,had1_xma=pr[[di]]$had1_xma,had1_ymb=pr[[di]]$had1_ymb),auto_unbox=TRUE))),recipient_pk=pk));.sendBlob(.to_b64url(s$sealed),"k2_peer_beaver_r1",dcf_conns[pi2])} }
        else if (ph==3) { for(di in 1:2){pi2<-3-di;pk<-.b64url_to_b64(transport_pks[[dcf_parties[pi2]]]);s<-dsVert:::.callMpcTool("transport-encrypt",list(data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(had2_xma=pr[[di]]$had2_xma,had2_ymb=pr[[di]]$had2_ymb),auto_unbox=TRUE))),recipient_pk=pk));.sendBlob(.to_b64url(s$sealed),"k2_peer_had2_r1",dcf_conns[pi2])} }
      }
    }
    if (isTRUE(weights_active)) {
      .glm_apply_shared_weight_residual(
        datasources = datasources,
        dcf_parties = dcf_parties,
        dcf_conns = dcf_conns,
        dealer_conn = dealer_conn_b,
        transport_pks = transport_pks,
        session_id = session_id,
        n_obs = n_obs,
        .dsAgg = .dsAgg,
        .sendBlob = .sendBlob,
        ring = 63L)
    }
    gt2 <- .dsAgg(datasources[dealer_conn_b], call(name = "glmRing63GenGradTriplesDS",
      dcf0_pk=transport_pks[[dcf_parties[1]]], dcf1_pk=transport_pks[[dcf_parties[2]]],
      n=as.integer(n_obs), p=as.integer(p_total), session_id=session_id))
    if (is.list(gt2)) gt2 <- gt2[[1]]
    .sendBlob(gt2$grad_blob_0,"k2_grad_triple_fp",dcf_conns[1])
    .sendBlob(gt2$grad_blob_1,"k2_grad_triple_fp",dcf_conns[2])
    br1 <- list()
    for (di in 1:2) {
      ci<-dcf_conns[di]; peer<-dcf_parties[3-di]
      .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS", session_id=session_id))
      r<-.dsAgg(datasources[ci], call(name = "k2GradientR1DS", peer_pk=transport_pks[[peer]], session_id=session_id))
      if(is.list(r)&&length(r)==1) r<-r[[1]]; br1[[di]]<-r
    }
    .sendBlob(br1[[1]]$encrypted_r1,"k2_grad_peer_r1",dcf_conns[2])
    .sendBlob(br1[[2]]$encrypted_r1,"k2_grad_peer_r1",dcf_conns[1])
    br2 <- list()
    for (di in 1:2) {
      r<-.dsAgg(datasources[dcf_conns[di]], call(name = "k2GradientR2DS", party_id=as.integer(di-1), session_id=session_id))
      if(is.list(r)&&length(r)==1) r<-r[[1]]; br2[[di]]<-r
    }
    ba <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(share_a=br2[[1]]$gradient_fp, share_b=br2[[2]]$gradient_fp, frac_bits=frac_bits))
    gp2 <- numeric(p_total); gii <- 1
    if(p_coord>0){gp2[beta_map[[coordinator]]]<-ba$values[gii:(gii+p_coord-1)];gii<-gii+p_coord}
    if(p_fusion>0){gp2[beta_map[[fusion_server]]]<-ba$values[gii:(gii+p_fusion-1)];gii<-gii+p_fusion}
    for(ns in non_dcf_servers){pn<-length(x_vars[[ns]]);if(pn==0)next;gp2[beta_map[[ns]]]<-ba$values[gii:(gii+pn-1)];gii<-gii+pn}
    ar2 <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(share_a=br1[[1]]$sum_residual_fp, share_b=br1[[2]]$sum_residual_fp, frac_bits=frac_bits))
    grad_backward <- c(ar2$values[1]/n_obs, gp2/n_obs) + lambda*theta_back

    # Central difference: H_j = (g_forward - g_backward) / (2delta)
    hessian[, j] <- (grad_forward - grad_backward) / (2 * delta)
    if (verbose) message(sprintf("    [SE] Column %d/%d", j, p_plus1))
  }

  # Symmetrize (numerical noise)
  hessian <- (hessian + t(hessian)) / 2
  # Pass raw hessian for proper Fisher matrix computation in outer code
  attr(inv_hessian, "raw_hessian") <- hessian
  }

  if (isTRUE(weights_active) && family != "gaussian") {
    if (verbose) {
      message("  [Deviance] Skipped for weighted non-Gaussian K>=3 fit; canonical weighted deviance is not implemented yet.")
    }
    result <- list(betas = betas_out, converged = converged,
                   final_iter = final_iter, deviance = NA_real_,
                   inv_hessian = inv_hessian)
    if (!label_intercept) result$intercept <- intercept
    return(result)
  }

  # ===========================================================================
  # Secure canonical deviance
  # ===========================================================================
  if (verbose) message("  [Deviance] Computing canonical deviance...")

  # Recompute eta from converged beta (SE computation may have overwritten shares)
  for (di in seq_along(dcf_parties)) {
    ci <- dcf_conns[di]
    is_coord <- (dcf_parties[di] == coordinator)
    b_coord <- if (p_coord > 0) beta[1:p_coord] else NULL
    b_nl_range <- (p_coord+1):p_total
    if (p_coord >= p_total) b_nl_range <- integer(0)
    b_nl <- if (is_coord) beta[b_nl_range] else {
      if (p_fusion > 0 && length(b_nl_range) > p_fusion) {
        beta[c((p_coord+p_fusion+1):p_total, (p_coord+1):(p_coord+p_fusion))]
      } else beta[b_nl_range]
    }
    .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
      beta_coord = b_coord, beta_nl = b_nl,
      intercept = intercept, is_coordinator = is_coord,
      session_id = session_id))
  }
  if (family == "gaussian") {
    for (di in seq_along(dcf_parties)) {
      .dsAgg(datasources[dcf_conns[di]], call(name = "k2IdentityLinkDS",
        session_id = session_id))
    }
  }

  # Helper: Beaver dot-product (nx1) on both DCF parties, returns scalar
  .k3_beaver_dot <- function() {
    dt <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenGradTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[1]]],
           dcf1_pk = transport_pks[[dcf_parties[2]]],
           n = as.integer(n_obs), p = 1L, session_id = session_id))
    if (is.list(dt)) dt <- dt[[1]]
    .sendBlob(dt$grad_blob_0, "k2_grad_triple_fp", dcf_conns[1])
    .sendBlob(dt$grad_blob_1, "k2_grad_triple_fp", dcf_conns[2])
    dr1 <- list()
    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]; peer <- dcf_parties[3-i]
      .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; dr1[[i]] <- r
    }
    .sendBlob(dr1[[1]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[2])
    .sendBlob(dr1[[2]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[1])
    dr2 <- list()
    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR2DS",
        party_id = as.integer(i-1), session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; dr2[[i]] <- r
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = dr2[[1]]$gradient_fp, share_b = dr2[[2]]$gradient_fp,
      frac_bits = frac_bits))
    agg$values[1]
  }

  if (family == "gaussian") {
    # Gaussian: RSS = canonical deviance. For weighted Gaussian fits,
    # compute shares of sqrt(w) * residual before the RSS Beaver dot.
    if (isTRUE(weights_active)) {
      .glm_apply_shared_weight_residual(
        datasources = datasources,
        dcf_parties = dcf_parties,
        dcf_conns = dcf_conns,
        dealer_conn = dealer_conn,
        transport_pks = transport_pks,
        session_id = session_id,
        n_obs = n_obs,
        .dsAgg = .dsAgg,
        .sendBlob = .sendBlob,
        weight_key = "k2_sqrt_weights_share_fp",
        output_key = "k2_sqrt_weighted_residual_share_fp",
        ring = 63L)
    }
    for (i in seq_along(dcf_parties))
      .dsAgg(datasources[dcf_conns[i]], call(name = "glmRing63PrepDevianceDS",
        mode = "rss", session_id = session_id))
    secure_deviance <- .k3_beaver_dot()

  } else if (family == "binomial") {
    # Canonical binomial: D = 2*(Sumsoftplus(eta) - y^T*eta) -- future work
    sp_dcf <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenDcfKeysDS",
           dcf0_pk = transport_pks[[dcf_parties[1]]],
           dcf1_pk = transport_pks[[dcf_parties[2]]],
           family = "softplus", n = as.integer(n_obs),
           frac_bits = as.integer(frac_bits),
           num_intervals = 80L, session_id = session_id))
    if (is.list(sp_dcf)) sp_dcf <- sp_dcf[[1]]
    .sendBlob(sp_dcf$dcf_blob_0, "k2_dcf_keys_persistent", dcf_conns[1])
    .sendBlob(sp_dcf$dcf_blob_1, "k2_dcf_keys_persistent", dcf_conns[2])
    for (i in seq_along(dcf_parties))
      .dsAgg(datasources[dcf_conns[i]], call(name = "k2StoreDcfKeysPersistentDS", session_id = session_id))
    sp_t <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenSplineTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[1]]],
           dcf1_pk = transport_pks[[dcf_parties[2]]],
           n = as.integer(n_obs), frac_bits = as.integer(frac_bits),
           session_id = session_id))
    if (is.list(sp_t)) sp_t <- sp_t[[1]]
    .sendBlob(sp_t$spline_blob_0, "k2_spline_triples", dcf_conns[1])
    .sendBlob(sp_t$spline_blob_1, "k2_spline_triples", dcf_conns[2])
    for (ph in 1:4) {
      pr <- list()
      for (di in seq_along(dcf_parties)) {
        ci <- dcf_conns[di]
        r <- .dsAgg(datasources[ci], call(paste0("k2WideSplinePhase",ph,"DS"),
          party_id = as.integer(di-1), family = "softplus",
          num_intervals = 80L, frac_bits = as.integer(frac_bits),
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; pr[[di]] <- r
      }
      if (ph==1) { .sendBlob(pr[[1]]$dcf_masked,"k2_peer_dcf_masked",dcf_conns[2]);.sendBlob(pr[[2]]$dcf_masked,"k2_peer_dcf_masked",dcf_conns[1]) }
      else if (ph==2) { for(di in 1:2){pi2<-3-di;pk<-.b64url_to_b64(transport_pks[[dcf_parties[pi2]]]);s<-dsVert:::.callMpcTool("transport-encrypt",list(data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(and_xma=pr[[di]]$and_xma,and_ymb=pr[[di]]$and_ymb,had1_xma=pr[[di]]$had1_xma,had1_ymb=pr[[di]]$had1_ymb),auto_unbox=TRUE))),recipient_pk=pk));.sendBlob(.to_b64url(s$sealed),"k2_peer_beaver_r1",dcf_conns[pi2])} }
      else if (ph==3) { for(di in 1:2){pi2<-3-di;pk<-.b64url_to_b64(transport_pks[[dcf_parties[pi2]]]);s<-dsVert:::.callMpcTool("transport-encrypt",list(data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(had2_xma=pr[[di]]$had2_xma,had2_ymb=pr[[di]]$had2_ymb),auto_unbox=TRUE))),recipient_pk=pk));.sendBlob(.to_b64url(s$sealed),"k2_peer_had2_r1",dcf_conns[pi2])} }
    }
    sums <- list()
    for (di in seq_along(dcf_parties)) {
      r <- .dsAgg(datasources[dcf_conns[di]], call(name = "glmRing63DevianceSumsDS",
        family = "binomial", session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; sums[[di]] <- r
    }
    sp_agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = sums[[1]]$sum_fp, share_b = sums[[2]]$sum_fp, frac_bits = frac_bits))
    sum_softplus <- sp_agg$values[1]
    for (i in seq_along(dcf_parties))
      .dsAgg(datasources[dcf_conns[i]], call(name = "glmRing63PrepDevianceDS",
        mode = "canonical", session_id = session_id))
    y_dot_eta <- .k3_beaver_dot()
    secure_deviance <- 2 * (sum_softplus - y_dot_eta)

  } else {
    # Poisson: D = 2*(Summu - y^T*eta + C)
    sums <- list()
    for (di in seq_along(dcf_parties)) {
      r <- .dsAgg(datasources[dcf_conns[di]], call(name = "glmRing63DevianceSumsDS",
        family = "poisson", session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; sums[[di]] <- r
    }
    mu_agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = sums[[1]]$sum_fp, share_b = sums[[2]]$sum_fp, frac_bits = frac_bits))
    sum_mu <- mu_agg$values[1]
    # Only the label server returns a non-zero null_term; sum to pick it up
    null_term <- sum(sapply(sums, function(s) if(!is.null(s$null_term)) s$null_term else 0))
    for (i in seq_along(dcf_parties))
      .dsAgg(datasources[dcf_conns[i]], call(name = "glmRing63PrepDevianceDS",
        mode = "canonical", session_id = session_id))
    y_dot_eta <- .k3_beaver_dot()
    secure_deviance <- 2 * (sum_mu - y_dot_eta + null_term)
  }

  if (verbose) message(sprintf("  [Deviance] = %.4f", secure_deviance))

  result <- list(betas = betas_out, converged = converged, final_iter = final_iter,
                 deviance = secure_deviance, inv_hessian = inv_hessian)
  if (!label_intercept) result$intercept <- intercept
  result
}
