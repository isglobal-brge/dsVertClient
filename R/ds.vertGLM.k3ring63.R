#' @title K>=3 Ring63 DCF + Beaver Gradient Loop (Zero CKKS, Zero Noise)
#' @description Pure Ring63 protocol for K>=3 binomial/Poisson GLM.
#'   Uses additive secret sharing of features to 2 DCF parties,
#'   DCF wide spline for sigmoid/exp, and Beaver matvec for gradient.
#'   No CKKS, no threshold decryption, no smudging noise.
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
#' Precision: ~1e-3 (DCF approximation only, no MHE noise).
#' @name k3-ring63-dcf-loop
NULL

.lbfgs_direction_ring63 <- function(grad, s_hist, y_hist) {
  k <- length(s_hist)
  if (k == 0) return(-grad)
  q <- grad; alpha_h <- numeric(k)
  for (i in k:1) {
    rho <- 1 / sum(y_hist[[i]] * s_hist[[i]])
    if (!is.finite(rho)) rho <- 0
    alpha_h[i] <- rho * sum(s_hist[[i]] * q)
    q <- q - alpha_h[i] * y_hist[[i]]
  }
  gamma <- sum(s_hist[[k]] * y_hist[[k]]) / sum(y_hist[[k]]^2)
  if (!is.finite(gamma) || gamma <= 0) gamma <- 1.0
  r <- gamma * q
  for (i in 1:k) {
    rho <- 1 / sum(y_hist[[i]] * s_hist[[i]])
    if (!is.finite(rho)) rho <- 0
    r <- r + s_hist[[i]] * (alpha_h[i] - rho * sum(y_hist[[i]] * r))
  }
  -r
}

#' @keywords internal
.k3_ring63_gradient_loop <- function(
    datasources, server_list, server_names, x_vars,
    coordinator, coordinator_conn, non_label_servers,
    transport_pks, std_data, y_var, family,
    betas, n_obs, lambda, log_n, log_scale,
    session_id, max_iter, tol, verbose,
    label_intercept, .dsAgg, .sendBlob) {

  K <- length(server_list)
  # DCF parties: first server = fusion (party 0), coordinator = party 1
  fusion_server <- server_list[1]
  fusion_conn <- which(server_names == fusion_server)
  dcf_parties <- c(fusion_server, coordinator)
  dcf_conns <- sapply(dcf_parties, function(s) which(server_names == s))

  frac_bits <- 20L
  num_intervals <- if (family == "poisson") 100L else 50L
  dcf_family <- if (family == "poisson") "poisson" else "sigmoid"

  .to_b64url <- function(x) gsub("+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-","+",gsub("_","/",x,fixed=TRUE),fixed=TRUE)
    pad <- nchar(x)%%4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }

  if (verbose) message(sprintf(
    "\n[Phase 3] K>=%d Ring63 DCF + Beaver L-BFGS (family=%s, n=%d, intervals=%d, lambda=%.1e)",
    K, family, n_obs, num_intervals, lambda))

  # ===========================================================================
  # Step A: Input sharing — all servers share features with 2 DCF parties
  # ===========================================================================
  t0_share <- proc.time()[[3]]
  if (verbose) message("  [Input Sharing] Sharing features+y with DCF parties...")

  # Total features across all servers (EXCLUDING intercept — handled separately)
  p_total <- sum(sapply(server_list, function(s) length(x_vars[[s]])))

  # Each server: split features into 2 shares, send to DCF parties
  for (server in server_list) {
    ci <- which(server_names == server)
    for (di in seq_along(dcf_parties)) {
      dcf_srv <- dcf_parties[di]
      dcf_ci <- dcf_conns[di]
      dcf_pk <- transport_pks[[dcf_srv]]

      if (server == dcf_srv) {
        # This server IS a DCF party — just call k2ShareInputDS targeting the other DCF party
        peer_dcf <- dcf_parties[3 - di]
        peer_pk_safe <- .to_b64url(transport_pks[[peer_dcf]])
        r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
          data_name = std_data, x_vars = x_vars[[server]],
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

  # Non-DCF servers: share features with BOTH DCF parties
  # CRITICAL: Split ONCE → send each half to one DCF party.
  # k2ShareInputDS(peer=fusion) → peer_share to fusion, own_share stays on server.
  # Then export own_share → coordinator (the complement).
  non_dcf_servers <- setdiff(server_list, dcf_parties)
  for (server in non_dcf_servers) {
    ci <- which(server_names == server)
    # Split once targeting fusion as peer
    fusion_pk_safe <- .to_b64url(transport_pks[[fusion_server]])
    r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
      data_name = std_data, x_vars = x_vars[[server]],
      y_var = NULL,
      peer_pk = fusion_pk_safe, session_id = session_id))
    if (is.list(r) && length(r) == 1) r <- r[[1]]
    # Send peer_share to fusion
    .sendBlob(r$encrypted_x_share, paste0("k2_extra_x_share_", server), fusion_conn)

    # Export own_share (the complement) to coordinator
    coord_pk_safe <- .to_b64url(transport_pks[[coordinator]])
    r2 <- .dsAgg(datasources[ci], call("glmRing63ExportOwnShareDS",
      peer_pk = coord_pk_safe, session_id = session_id))
    if (is.list(r2) && length(r2) == 1) r2 <- r2[[1]]
    .sendBlob(r2$encrypted_own_share, paste0("k2_extra_x_share_", server), coordinator_conn)
  }

  # DCF parties: receive peer's shares
  for (di in seq_along(dcf_parties)) {
    dcf_srv <- dcf_parties[di]
    dcf_ci <- dcf_conns[di]
    peer_p <- length(x_vars[[dcf_parties[3 - di]]])  # no intercept in features
    .dsAgg(datasources[dcf_ci], call("k2ReceiveShareDS",
      peer_p = as.integer(peer_p), session_id = session_id))
  }

  # DCF parties: receive non-DCF feature shares
  for (server in non_dcf_servers) {
    for (di in seq_along(dcf_parties)) {
      dcf_ci <- dcf_conns[di]
      # Receive and store the extra feature shares
      # These will be assembled into the full X matrix
      .dsAgg(datasources[dcf_ci], call("glmRing63ReceiveExtraShareDS",
        extra_key = paste0("k2_extra_x_share_", server),
        extra_p = as.integer(length(x_vars[[server]])),
        session_id = session_id))
    }
  }

  if (verbose) message(sprintf("  [Input Sharing] Complete (p_total=%d, %.1fs)",
    p_total, proc.time()[[3]] - t0_share))

  # ===========================================================================
  # Step B: Pre-generate DCF keys (on NON-DCF server, not client)
  # ===========================================================================
  # Security: DCF keys generated server-side → client never sees them.
  # A malicious client cannot craft keys to leak information.
  t0_dcf <- proc.time()[[3]]
  dealer <- non_dcf_servers[1]  # first non-DCF server acts as dealer
  dealer_conn <- which(server_names == dealer)
  if (verbose) message(sprintf("  [DCF] Server %s generating keys (n=%d, %d intervals)...",
                                 dealer, n_obs, num_intervals))

  dcf_result <- .dsAgg(datasources[dealer_conn],
    call("glmRing63GenDcfKeysDS",
         dcf0_pk = transport_pks[[dcf_parties[1]]],
         dcf1_pk = transport_pks[[dcf_parties[2]]],
         family = dcf_family, n = as.integer(n_obs),
         frac_bits = frac_bits, num_intervals = num_intervals,
         session_id = session_id))
  if (is.list(dcf_result)) dcf_result <- dcf_result[[1]]

  # Relay opaque blobs to DCF parties (client can't read them)
  .sendBlob(dcf_result$dcf_blob_0, "k2_dcf_keys_persistent", dcf_conns[1])
  .dsAgg(datasources[dcf_conns[1]], call("k2StoreDcfKeysPersistentDS", session_id = session_id))
  .sendBlob(dcf_result$dcf_blob_1, "k2_dcf_keys_persistent", dcf_conns[2])
  .dsAgg(datasources[dcf_conns[2]], call("k2StoreDcfKeysPersistentDS", session_id = session_id))

  if (verbose) message(sprintf("  [DCF] Keys distributed (%.1fs)", proc.time()[[3]] - t0_dcf))

  # ===========================================================================
  # Iteration loop
  # ===========================================================================
  beta <- rep(0, p_total)
  intercept <- 0.0
  lbfgs_s <- list(); lbfgs_y <- list()
  prev_theta <- NULL; prev_grad <- NULL
  converged <- FALSE; final_iter <- 0

  # Feature counts for beta splitting (EXCLUDING intercept)
  p_coord <- length(x_vars[[coordinator]])
  p_fusion <- length(x_vars[[fusion_server]])
  p_extras <- p_total - p_coord - p_fusion
  p_nl <- p_total - p_coord

  t0_loop <- proc.time()[[3]]

  for (iter in seq_len(max_iter)) {
    t0_iter <- proc.time()[[3]]
    beta_old <- beta
    intercept_old <- intercept

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
    beta_map[[coordinator]] <- idx:(idx + p_coord - 1); idx <- idx + p_coord
    for (srv in server_list) {
      if (srv == coordinator) next
      p_s <- length(x_vars[[srv]])
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
        b_coord <- beta[beta_map[[coordinator]]]
        b_nl <- c(beta[beta_map[[peer_srv]]])
        for (ns in non_dcf_servers) b_nl <- c(b_nl, beta[beta_map[[ns]]])
      } else {
        # Fusion: X = [peer(coord) | extra... | own(fusion)]
        # The Go cmd with is_party_zero=FALSE puts: [peer | own]
        # peer = coord_share + extras, own = fusion_share
        # So beta_coord = coord betas, beta_nl = [extra betas | own betas]
        b_coord <- beta[beta_map[[coordinator]]]
        b_nl <- c()
        for (ns in non_dcf_servers) b_nl <- c(b_nl, beta[beta_map[[ns]]])
        b_nl <- c(b_nl, beta[beta_map[[srv]]])  # own (fusion) last
      }

      .dsAgg(datasources[ci], call("k2ComputeEtaShareDS",
        beta_coord = b_coord, beta_nl = b_nl,
        intercept = if (is_coord) intercept else 0.0,
        is_coordinator = is_coord, session_id = session_id))

      # Reorder fusion's X_full to canonical order [coord|fusion|extras]
      # so Beaver gradient works (both parties must have same column layout)
      if (!is_coord && p_extras > 0) {
        .dsAgg(datasources[ci], call("glmRing63ReorderXFullDS",
          p_coord = as.integer(p_coord), p_fusion = as.integer(p_fusion),
          p_extras = as.integer(p_extras), session_id = session_id))
      }
    }

    # =================================================================
    # Step 2: DCF Wide Spline (4-phase between DCF parties)
    # =================================================================
    if (verbose && iter <= 3) message(sprintf("    [%d.2] DCF wide spline...", iter))

    # Generate spline Beaver triples on dealer server (not client)
    spline_t <- .dsAgg(datasources[dealer_conn],
      call("glmRing63GenSplineTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[1]]],
           dcf1_pk = transport_pks[[dcf_parties[2]]],
           n = as.integer(n_obs), frac_bits = frac_bits,
           session_id = session_id))
    if (is.list(spline_t)) spline_t <- spline_t[[1]]
    # Relay opaque blobs (client can't read)
    .sendBlob(spline_t$spline_blob_0, "k2_spline_triples", dcf_conns[1])
    .sendBlob(spline_t$spline_blob_1, "k2_spline_triples", dcf_conns[2])

    # Phase 1-4 (same as K=2)
    ph1 <- list()
    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase1DS",
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
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase2DS",
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
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", peer_ci)
    }

    ph3 <- list()
    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase3DS",
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
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", peer_ci)
    }

    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]
      .dsAgg(datasources[ci], call("k2WideSplinePhase4DS",
        party_id = as.integer(i - 1), family = family,
        num_intervals = num_intervals, frac_bits = frac_bits,
        session_id = session_id))
    }

    # =================================================================
    # Step 3: Beaver gradient (Ring63, exact, no CKKS)
    # =================================================================
    if (verbose && iter <= 3) message(sprintf("    [%d.3] Beaver gradient...", iter))

    # Generate gradient Beaver triples on dealer server (not client)
    grad_t <- .dsAgg(datasources[dealer_conn],
      call("glmRing63GenGradTriplesDS",
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
      .dsAgg(datasources[ci], call("k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call("k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      r1_results[[i]] <- r
    }
    .sendBlob(r1_results[[1]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[2])
    .sendBlob(r1_results[[2]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[1])

    grad_results <- list()
    for (i in seq_along(dcf_parties)) {
      ci <- dcf_conns[i]
      r <- .dsAgg(datasources[ci], call("k2GradientR2DS",
        party_id = as.integer(i - 1), session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      grad_results[[i]] <- r
    }

    # =================================================================
    # Step 4: Aggregate gradient + L-BFGS update
    # =================================================================
    # Both DCF parties now have X_full in canonical order [coord|fusion|extras]
    # (fusion was reordered by glmRing63ReorderXFullDS)
    # Simple Ring63 aggregation — no permutation needed.
    agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
      share_a = grad_results[[1]]$gradient_fp,
      share_b = grad_results[[2]]$gradient_fp,
      frac_bits = frac_bits))
    gradient_canonical <- agg$values  # in canonical order [coord|fusion|extras]

    agg_res <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
      share_a = r1_results[[1]]$sum_residual_fp,
      share_b = r1_results[[2]]$sum_residual_fp,
      frac_bits = frac_bits))
    sum_residual <- agg_res$values[1]

    # Remap from canonical [coord | fusion | extras] to beta order [coord | non-coord-by-server-list]
    gradient <- numeric(p_total)
    gradient[beta_map[[coordinator]]] <- gradient_canonical[1:p_coord]
    gradient[beta_map[[fusion_server]]] <- gradient_canonical[(p_coord + 1):(p_coord + p_fusion)]
    grad_idx <- p_coord + p_fusion + 1
    for (ns in non_dcf_servers) {
      p_ns <- length(x_vars[[ns]])
      gradient[beta_map[[ns]]] <- gradient_canonical[grad_idx:(grad_idx + p_ns - 1)]
      grad_idx <- grad_idx + p_ns
    }

    theta <- c(intercept, beta)
    full_grad <- c(sum_residual / n_obs, gradient / n_obs) + lambda * theta

    if (!is.null(prev_theta)) {
      sk <- theta - prev_theta; yk <- full_grad - prev_grad
      if (sum(sk * yk) > 1e-10) {
        lbfgs_s <- c(lbfgs_s, list(sk)); lbfgs_y <- c(lbfgs_y, list(yk))
        if (length(lbfgs_s) > 7) { lbfgs_s <- lbfgs_s[-1]; lbfgs_y <- lbfgs_y[-1] }
      }
    }
    prev_theta <- theta; prev_grad <- full_grad

    direction <- .lbfgs_direction_ring63(full_grad, lbfgs_s, lbfgs_y)
    step_size <- if (iter <= 1) 0.3 else 1.0
    new_theta <- theta + step_size * direction
    intercept <- new_theta[1]
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

  # Build betas in per-server format (coordinator gets intercept prepended)
  betas_out <- list()
  betas_out[[coordinator]] <- c(intercept, beta[1:p_coord])
  idx <- p_coord + 1
  for (server in server_list) {
    if (server == coordinator) next
    p_s <- length(x_vars[[server]])
    betas_out[[server]] <- beta[idx:(idx + p_s - 1)]
    idx <- idx + p_s
  }

  list(betas = betas_out, converged = converged, final_iter = final_iter)
}
