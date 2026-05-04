#' @title K=2 Beaver L-BFGS Pipeline
#' @description Secure K=2 GLM with L-BFGS optimizer.
#'   Gaussian: identity link, constant Fisher, L-BFGS (2 rounds/iter).
#'   Binomial/Poisson: wide spline sigmoid/exp, L-BFGS (6 rounds/iter, no Fisher).
#' @name k2-beaver-lbfgs-client
NULL

#' L-BFGS two-loop recursion (Nocedal & Wright, Algorithm 7.4)
#' @keywords internal
.lbfgs_direction <- function(grad, s_hist, y_hist) {
  k <- length(s_hist)
  if (k == 0) return(-grad)  # first iteration: steepest descent

  q <- grad
  alpha_h <- numeric(k)

  # Forward loop
  for (i in k:1) {
    rho <- 1 / sum(y_hist[[i]] * s_hist[[i]])
    if (!is.finite(rho)) rho <- 0
    alpha_h[i] <- rho * sum(s_hist[[i]] * q)
    q <- q - alpha_h[i] * y_hist[[i]]
  }

  # Initial Hessian: gamma * I (scaled by most recent curvature pair)
  gamma <- sum(s_hist[[k]] * y_hist[[k]]) / sum(y_hist[[k]]^2)
  if (!is.finite(gamma) || gamma <= 0) gamma <- 1.0
  r <- gamma * q

  # Backward loop
  for (i in 1:k) {
    rho <- 1 / sum(y_hist[[i]] * s_hist[[i]])
    if (!is.finite(rho)) rho <- 0
    r <- r + s_hist[[i]] * (alpha_h[i] - rho * sum(y_hist[[i]] * r))
  }

  return(-r)
}

#' @keywords internal
.k2_strict_loop <- function(datasources, server_names, server_list,
                             coordinator, coordinator_conn,
                             non_label_servers, nl, nl_conn,
                             x_vars, y_var, std_data,
                             transport_pks, session_id,
                             family, lambda, max_iter, tol,
                             n_obs, verbose, .dsAgg, .sendBlob,
                             weights_active = FALSE,
                             no_intercept  = FALSE,
                             ring = 63L,
                             compute_se = TRUE,
                             compute_deviance = TRUE,
                             gradient_only = FALSE,
                             start = NULL) {

  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  ring_tag <- if (ring == 127L) "ring127" else "ring63"
  frac_bits <- if (ring == 127L) 50L else 20L
  is_gaussian <- (family == "gaussian")
  default_intervals <- if (family == "poisson") 100L else if (family == "binomial") 100L else 50L
  opt_name <- paste0("dsvert.glm_num_intervals_", family)
  num_intervals <- suppressWarnings(as.integer(
    getOption(opt_name, getOption("dsvert.glm_num_intervals",
                                  default_intervals))[[1L]]
  ))
  if (!is.finite(num_intervals) || num_intervals < 10L) {
    num_intervals <- default_intervals
  }
  if (is.null(lambda)) lambda <- 1e-4
  if (verbose) {
    if (is_gaussian) {
      message("  Gaussian one-shot: X^T X + X^T y via Beaver, direct solve")
    } else {
      message(sprintf("  L-BFGS: %s, %d-interval spline (6 rounds/iter, no Fisher)", family, num_intervals))
    }
  }

  p_coord <- length(x_vars[[coordinator]])
  p_nl <- length(x_vars[[nl]])
  p_total <- p_coord + p_nl

  .to_b64url <- function(x) gsub("+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-","+",gsub("_","/",x,fixed=TRUE),fixed=TRUE)
    pad <- nchar(x)%%4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }

  # === INPUT-SHARING PREAMBLE ===
  if (verbose) message("  [Input Sharing] Generating additive secret shares (2 servers)...")
  t0_share <- proc.time()[[3]]
  share_results <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    peer_pk_safe <- .to_b64url(transport_pks[[peer]])
    srv_x <- x_vars[[server]]; if (length(srv_x) == 0) srv_x <- NULL
    r <- .dsAgg(datasources[ci], call(name = "k2ShareInputDS",
      data_name = std_data, x_vars = srv_x,
      y_var = if (server == coordinator) y_var else NULL,
      peer_pk = peer_pk_safe, session_id = session_id,
      ring = ring))
    if (is.list(r) && length(r) == 1) r <- r[[1]]
    share_results[[server]] <- r
  }
  for (server in server_list) {
    peer <- setdiff(server_list, server)
    peer_ci <- which(server_names == peer)
    .sendBlob(share_results[[server]]$encrypted_x_share, "k2_peer_x_share", peer_ci)
    if (!is.null(share_results[[server]]$encrypted_y_share))
      .sendBlob(share_results[[server]]$encrypted_y_share, "k2_peer_y_share", peer_ci)
  }
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    .dsAgg(datasources[ci], call(name = "k2ReceiveShareDS",
      peer_p = as.integer(length(x_vars[[peer]])), session_id = session_id))
  }
  if (verbose) message(sprintf("  [Input Sharing] Complete: p_coord=%d, p_nl=%d (%.1fs)",
                                 p_coord, p_nl, proc.time()[[3]] - t0_share))

  # === PRE-GENERATE DCF KEYS (server-side, not client) ===
  # K=2: only 1 non-DCF server (non-label) -> fixed dealer, no rotation.
  # Security: analyst needs both servers to extract data.
  # For K>=3: dealer rotation provides stronger guarantees.
  dealer <- nl; dealer_conn <- nl_conn
  if (!is_gaussian) {
    t0_dcf <- proc.time()[[3]]
    if (verbose) message(sprintf("  [DCF] Server %s generating keys (n=%d, %d intervals)...",
                                   dealer, n_obs, num_intervals))
    dcf_result <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenDcfKeysDS",
           dcf0_pk = transport_pks[[coordinator]],
           dcf1_pk = transport_pks[[nl]],
           family = if (family == "poisson") "poisson" else "sigmoid",
           n = as.integer(n_obs), frac_bits = frac_bits,
           num_intervals = num_intervals,
           ring = ring,
           session_id = session_id))
    if (is.list(dcf_result)) dcf_result <- dcf_result[[1]]
    .sendBlob(dcf_result$dcf_blob_0, "k2_dcf_keys_persistent", coordinator_conn)
    .dsAgg(datasources[coordinator_conn], call(name = "k2StoreDcfKeysPersistentDS", session_id = session_id))
    .sendBlob(dcf_result$dcf_blob_1, "k2_dcf_keys_persistent", nl_conn)
    .dsAgg(datasources[nl_conn], call(name = "k2StoreDcfKeysPersistentDS", session_id = session_id))
    if (verbose) message(sprintf("  [DCF] Keys distributed (%.1fs)", proc.time()[[3]] - t0_dcf))
  }

  beta <- rep(0, p_total)
  intercept <- 0.0
  feature_order <- c(x_vars[[coordinator]],
                     if (p_nl > 0L) x_vars[[nl]] else character(0))
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
  beta <- unname(beta)
  intercept <- unname(intercept)
  converged <- FALSE
  final_iter <- 0
  max_diff <- Inf
  last_gradient <- NULL

  # L-BFGS state (client-side only)
  lbfgs_m <- 7L  # number of (s,y) pairs to keep
  s_hist <- list()
  y_hist <- list()
  prev_theta <- NULL
  prev_grad <- NULL

  if (verbose) message(sprintf("\n[Phase 3] K=2 L-BFGS iterations (p=%d, n=%d, lambda=%.1e)", p_total, n_obs, lambda))

  for (iter in seq_len(max_iter)) {
    t0_iter <- proc.time()[[3]]
    beta_old <- beta
    intercept_old <- intercept

    # === Step 1: Compute eta ===
    # Guard against p_nl == 0 (all predictors on coordinator, e.g.
    # sleepstudy Reaction~Days with peer holding only patient_id):
    # `beta[(p_coord+1):p_total]` would produce `beta[2:1]` = reversed
    # 2-elt vector with NAs via R's `:` semantics. Use explicit
    # seq_len-based slicing that is safe for zero-length cases.
    beta_coord_slice <- if (p_coord > 0L) unname(beta[seq_len(p_coord)]) else numeric(0)
    beta_nl_slice <- if (p_total > p_coord) unname(beta[(p_coord + 1L):p_total]) else numeric(0)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
        beta_coord = beta_coord_slice, beta_nl = beta_nl_slice,
        intercept = if (is_coord) intercept else 0.0,
        is_coordinator = is_coord, session_id = session_id))
    }

    # === Step 2: Link function ===
    if (is_gaussian) {
      # Identity link: mu = eta (copy eta share to mu share)
      for (server in server_list) {
        ci <- which(server_names == server)
        .dsAgg(datasources[ci], call(name = "k2IdentityLinkDS", session_id = session_id))
      }
    } else {
      # Spline triples generated on dealer server (not client)
      spline_t <- .dsAgg(datasources[dealer_conn],
        call(name = "glmRing63GenSplineTriplesDS",
             dcf0_pk = transport_pks[[coordinator]],
             dcf1_pk = transport_pks[[nl]],
             n = as.integer(n_obs), frac_bits = frac_bits,
             ring = ring,
             session_id = session_id))
      if (is.list(spline_t)) spline_t <- spline_t[[1]]
      .sendBlob(spline_t$spline_blob_0, "k2_spline_triples", coordinator_conn)
      .sendBlob(spline_t$spline_blob_1, "k2_spline_triples", nl_conn)

      ph1 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call(name = "k2WideSplinePhase1DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          ring = ring,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; ph1[[server]] <- r
      }
      .sendBlob(ph1[[coordinator]]$dcf_masked, "k2_peer_dcf_masked", nl_conn)
      .sendBlob(ph1[[nl]]$dcf_masked, "k2_peer_dcf_masked", coordinator_conn)

      ph2 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call(name = "k2WideSplinePhase2DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          ring = ring,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; ph2[[server]] <- r
      }
      for (server in server_list) {
        peer <- setdiff(server_list, server); peer_ci <- which(server_names == peer)
        pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
        r1_json <- jsonlite::toJSON(list(
          and_xma=ph2[[server]]$and_xma, and_ymb=ph2[[server]]$and_ymb,
          had1_xma=ph2[[server]]$had1_xma, had1_ymb=ph2[[server]]$had1_ymb),
          auto_unbox=TRUE)
        sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", peer_ci)
      }

      ph3 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call(name = "k2WideSplinePhase3DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          ring = ring,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; ph3[[server]] <- r
      }
      for (server in server_list) {
        peer <- setdiff(server_list, server); peer_ci <- which(server_names == peer)
        pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
        r1_json <- jsonlite::toJSON(list(
          had2_xma=ph3[[server]]$had2_xma, had2_ymb=ph3[[server]]$had2_ymb),
          auto_unbox=TRUE)
        sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", peer_ci)
      }

      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        .dsAgg(datasources[ci], call(name = "k2WideSplinePhase4DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          ring = ring,
          session_id = session_id))
      }
    }  # close else (non-Gaussian wide spline)

    # === Optional: apply per-patient weights before gradient ===
    # Weights are secret-shared between the two DCF parties. One Beaver
    # vecmul round computes shares of w * (mu - y), which are then fed to
    # the existing X^T r Beaver gradient path.
    if (isTRUE(weights_active)) {
      .glm_apply_shared_weight_residual(
        datasources = datasources,
        dcf_parties = c(coordinator, nl),
        dcf_conns = c(coordinator_conn, nl_conn),
        dealer_conn = dealer_conn,
        transport_pks = transport_pks,
        session_id = session_id,
        n_obs = n_obs,
        .dsAgg = .dsAgg,
        .sendBlob = .sendBlob,
        ring = ring)
    }

    # === Step 3: Gradient (Beaver matvec) ===
    # Gradient triples generated on dealer server (not client)
    grad_t <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenGradTriplesDS",
           dcf0_pk = transport_pks[[coordinator]],
           dcf1_pk = transport_pks[[nl]],
           n = as.integer(n_obs), p = as.integer(p_total),
           ring = ring,
           session_id = session_id))
    if (is.list(grad_t)) grad_t <- grad_t[[1]]
    .sendBlob(grad_t$grad_blob_0, "k2_grad_triple_fp", coordinator_conn)
    .sendBlob(grad_t$grad_blob_1, "k2_grad_triple_fp", nl_conn)
    r1_results <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      r1_results[[server]] <- r
    }
    .sendBlob(r1_results[[coordinator]]$encrypted_r1, "k2_grad_peer_r1", nl_conn)
    .sendBlob(r1_results[[nl]]$encrypted_r1, "k2_grad_peer_r1", coordinator_conn)
    grad_results <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR2DS",
        party_id = if(is_coord) 0L else 1L, session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      grad_results[[server]] <- r
    }

    # === Step 4: Aggregate gradient + L-BFGS update ===
    grad_fp_coord <- grad_results[[coordinator]]$gradient_fp
    grad_fp_nl <- grad_results[[nl]]$gradient_fp
    res_fp_coord <- r1_results[[coordinator]]$sum_residual_fp
    res_fp_nl <- r1_results[[nl]]$sum_residual_fp
    if (!is.null(grad_fp_coord) && !is.null(grad_fp_nl)) {
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a=grad_fp_coord, share_b=grad_fp_nl,
        frac_bits=frac_bits, ring=ring_tag))
      gradient <- agg$values
      agg_res <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a=res_fp_coord, share_b=res_fp_nl,
        frac_bits=frac_bits, ring=ring_tag))
      sum_residual <- agg_res$values[1]
    } else {
      gradient <- rep(0, p_total); sum_residual <- 0
      for (server in server_list) {
        gradient <- gradient + grad_results[[server]]$gradient_share
        sum_residual <- sum_residual + r1_results[[server]]$sum_residual
      }
    }

    # Full gradient with L2 regularization
    theta <- c(intercept, beta)
    full_grad <- c(sum_residual / n_obs, gradient / n_obs) + lambda * theta
    names(full_grad) <- c("(Intercept)", feature_order)
    if (isTRUE(no_intercept)) full_grad[1L] <- 0
    last_gradient <- full_grad

    if (verbose && iter <= 3) {
      message(sprintf("  [L-BFGS] sum_res=%.4f, ||grad||=%.6f, grad_range=[%.4f, %.4f]",
        sum_residual, sqrt(sum(full_grad^2)), min(full_grad), max(full_grad)))
    }

    if (isTRUE(gradient_only)) {
      converged <- TRUE
      final_iter <- iter
      max_diff <- 0
      if (verbose) {
        message(sprintf("  [Gradient] returned aggregate score at start (||grad||=%.4f)",
                        sqrt(sum(full_grad^2))))
      }
      break
    }

    # Update L-BFGS history
    if (!is.null(prev_theta)) {
      sk <- theta - prev_theta
      yk <- full_grad - prev_grad
      if (sum(sk * yk) > 1e-10) {  # curvature condition
        s_hist <- c(s_hist, list(sk))
        y_hist <- c(y_hist, list(yk))
        if (length(s_hist) > lbfgs_m) {
          s_hist <- s_hist[-1]
          y_hist <- y_hist[-1]
        }
      }
    }
    prev_theta <- theta
    prev_grad <- full_grad

    # Compute L-BFGS direction
    direction <- .lbfgs_direction(full_grad, s_hist, y_hist)

    # Cautious first step (helps Poisson convergence), full step after
    step_size <- if (iter <= 1) 0.3 else 1.0
    if (isTRUE(no_intercept)) {
      # Caller's design matrix encodes its own intercept column; we
      # freeze the auto-intercept at zero by masking its direction.
      direction[1] <- 0
    }
    new_theta <- unname(theta + step_size * direction)
    intercept <- if (isTRUE(no_intercept)) 0 else unname(new_theta[1])
    beta <- unname(new_theta[-1])

    max_diff <- max(abs(beta - beta_old), abs(intercept - intercept_old))
    final_iter <- iter

    if (verbose) message(sprintf("  Iter %d: ||grad||=%.4f, step=%.2f, diff=%.2e, theta=[%.3f, %.3f] (%.1fs)",
      iter, sqrt(sum(full_grad^2)), step_size, max_diff,
      min(new_theta), max(new_theta), proc.time()[[3]] - t0_iter))

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose) message(sprintf("  Converged after %d iterations (diff = %.2e)", iter, max_diff))
      break
    }
    if (iter %% 10 == 0) {
      for (server in server_list)
        tryCatch(.dsAgg(datasources[which(server_names==server)], call(name = "mpcGcDS")), error=function(e) NULL)
    }
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)", max_iter, max_diff))

  inv_hessian <- list()
  if (isTRUE(compute_se) && !isTRUE(gradient_only)) {
  # === Standard errors via finite-difference Hessian (K=2) ===
  if (verbose) message("  [SE] Computing Hessian (central differences)...")
  p_plus1 <- p_total + 1
  delta_se <- 0.01
  hessian_k2 <- matrix(0, p_plus1, p_plus1)
  theta_conv <- c(intercept, beta)
  grad_conv <- full_grad

  for (jj in seq_len(p_plus1)) {
    th_p <- theta_conv; th_p[jj] <- th_p[jj] + delta_se
    int_p <- th_p[1]; bet_p <- th_p[-1]
    bet_p_coord <- if (p_coord > 0L) bet_p[seq_len(p_coord)] else numeric(0)
    bet_p_nl <- if (p_total > p_coord) bet_p[(p_coord + 1L):p_total] else numeric(0)
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
        beta_coord = bet_p_coord, beta_nl = bet_p_nl,
        intercept = if (is_coord) int_p else 0.0,
        is_coordinator = is_coord, session_id = session_id))
    }
    if (is_gaussian) {
      for (server in server_list) .dsAgg(datasources[which(server_names==server)],
        call(name = "k2IdentityLinkDS", session_id=session_id))
    } else {
      st <- .dsAgg(datasources[dealer_conn], call(name = "glmRing63GenSplineTriplesDS",
        dcf0_pk=transport_pks[[coordinator]], dcf1_pk=transport_pks[[nl]],
        n=as.integer(n_obs), frac_bits=frac_bits, ring=ring, session_id=session_id))
      if(is.list(st)) st<-st[[1]]
      .sendBlob(st$spline_blob_0,"k2_spline_triples",coordinator_conn)
      .sendBlob(st$spline_blob_1,"k2_spline_triples",nl_conn)
      for(ph in 1:4){pr<-list(); for(server in server_list){ci<-which(server_names==server);is_coord<-(server==coordinator);r<-.dsAgg(datasources[ci],call(paste0("k2WideSplinePhase",ph,"DS"),party_id=if(is_coord)0L else 1L,family=family,num_intervals=num_intervals,frac_bits=frac_bits,ring=ring,session_id=session_id));if(is.list(r)&&length(r)==1)r<-r[[1]];pr[[server]]<-r};if(ph==1){.sendBlob(pr[[coordinator]]$dcf_masked,"k2_peer_dcf_masked",nl_conn);.sendBlob(pr[[nl]]$dcf_masked,"k2_peer_dcf_masked",coordinator_conn)}else if(ph==2){for(server in server_list){peer<-setdiff(server_list,server);peer_ci<-which(server_names==peer);pk_b64<-.b64url_to_b64(transport_pks[[peer]]);sealed<-dsVert:::.callMpcTool("transport-encrypt",list(data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(and_xma=pr[[server]]$and_xma,and_ymb=pr[[server]]$and_ymb,had1_xma=pr[[server]]$had1_xma,had1_ymb=pr[[server]]$had1_ymb),auto_unbox=TRUE))),recipient_pk=pk_b64));.sendBlob(.to_b64url(sealed$sealed),"k2_peer_beaver_r1",peer_ci)}}else if(ph==3){for(server in server_list){peer<-setdiff(server_list,server);peer_ci<-which(server_names==peer);pk_b64<-.b64url_to_b64(transport_pks[[peer]]);sealed<-dsVert:::.callMpcTool("transport-encrypt",list(data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(had2_xma=pr[[server]]$had2_xma,had2_ymb=pr[[server]]$had2_ymb),auto_unbox=TRUE))),recipient_pk=pk_b64));.sendBlob(.to_b64url(sealed$sealed),"k2_peer_had2_r1",peer_ci)}}}
    }
    if (isTRUE(weights_active)) {
      .glm_apply_shared_weight_residual(
        datasources = datasources,
        dcf_parties = c(coordinator, nl),
        dcf_conns = c(coordinator_conn, nl_conn),
        dealer_conn = dealer_conn,
        transport_pks = transport_pks,
        session_id = session_id,
        n_obs = n_obs,
        .dsAgg = .dsAgg,
        .sendBlob = .sendBlob,
        ring = ring)
    }
    gt <- .dsAgg(datasources[dealer_conn], call(name = "glmRing63GenGradTriplesDS",
      dcf0_pk=transport_pks[[coordinator]], dcf1_pk=transport_pks[[nl]],
      n=as.integer(n_obs), p=as.integer(p_total), ring=ring, session_id=session_id))
    if(is.list(gt)) gt<-gt[[1]]
    .sendBlob(gt$grad_blob_0,"k2_grad_triple_fp",coordinator_conn)
    .sendBlob(gt$grad_blob_1,"k2_grad_triple_fp",nl_conn)
    r1p<-list(); for(server in server_list){ci<-which(server_names==server);peer<-setdiff(server_list,server);.dsAgg(datasources[ci],call(name = "k2StoreGradTripleDS",session_id=session_id));r<-.dsAgg(datasources[ci],call(name = "k2GradientR1DS",peer_pk=transport_pks[[peer]],session_id=session_id));if(is.list(r)&&length(r)==1)r<-r[[1]];r1p[[server]]<-r}
    .sendBlob(r1p[[coordinator]]$encrypted_r1,"k2_grad_peer_r1",nl_conn)
    .sendBlob(r1p[[nl]]$encrypted_r1,"k2_grad_peer_r1",coordinator_conn)
    r2p<-list(); for(server in server_list){ci<-which(server_names==server);is_coord<-(server==coordinator);r<-.dsAgg(datasources[ci],call(name = "k2GradientR2DS",party_id=if(is_coord)0L else 1L,session_id=session_id));if(is.list(r)&&length(r)==1)r<-r[[1]];r2p[[server]]<-r}
    ag<-dsVert:::.callMpcTool("k2-ring63-aggregate",list(share_a=r2p[[coordinator]]$gradient_fp,share_b=r2p[[nl]]$gradient_fp,frac_bits=frac_bits,ring=ring_tag))
    ar<-dsVert:::.callMpcTool("k2-ring63-aggregate",list(share_a=r1p[[coordinator]]$sum_residual_fp,share_b=r1p[[nl]]$sum_residual_fp,frac_bits=frac_bits,ring=ring_tag))
    gp<-c(ar$values[1]/n_obs,ag$values/n_obs)+lambda*th_p
    grad_fwd <- gp
    # Backward perturbation for central difference
    th_m <- theta_conv; th_m[jj] <- th_m[jj] - delta_se
    int_m <- th_m[1]; bet_m <- th_m[-1]
    bet_m_coord <- if (p_coord > 0L) bet_m[seq_len(p_coord)] else numeric(0)
    bet_m_nl <- if (p_total > p_coord) bet_m[(p_coord + 1L):p_total] else numeric(0)
    for(server in server_list){ci<-which(server_names==server);is_coord<-(server==coordinator);.dsAgg(datasources[ci],call(name = "k2ComputeEtaShareDS",beta_coord=bet_m_coord,beta_nl=bet_m_nl,intercept=if(is_coord)int_m else 0,is_coordinator=is_coord,session_id=session_id))}
    if(is_gaussian){for(server in server_list).dsAgg(datasources[which(server_names==server)],call(name = "k2IdentityLinkDS",session_id=session_id))}else{stb<-.dsAgg(datasources[dealer_conn],call(name = "glmRing63GenSplineTriplesDS",dcf0_pk=transport_pks[[coordinator]],dcf1_pk=transport_pks[[nl]],n=as.integer(n_obs),frac_bits=frac_bits,ring=ring,session_id=session_id));if(is.list(stb))stb<-stb[[1]];.sendBlob(stb$spline_blob_0,"k2_spline_triples",coordinator_conn);.sendBlob(stb$spline_blob_1,"k2_spline_triples",nl_conn);for(ph in 1:4){pr<-list();for(server in server_list){ci<-which(server_names==server);is_coord<-(server==coordinator);r<-.dsAgg(datasources[ci],call(paste0("k2WideSplinePhase",ph,"DS"),party_id=if(is_coord)0L else 1L,family=family,num_intervals=num_intervals,frac_bits=frac_bits,ring=ring,session_id=session_id));if(is.list(r)&&length(r)==1)r<-r[[1]];pr[[server]]<-r};if(ph==1){.sendBlob(pr[[coordinator]]$dcf_masked,"k2_peer_dcf_masked",nl_conn);.sendBlob(pr[[nl]]$dcf_masked,"k2_peer_dcf_masked",coordinator_conn)}else if(ph==2){for(server in server_list){peer<-setdiff(server_list,server);peer_ci<-which(server_names==peer);pk_b64<-.b64url_to_b64(transport_pks[[peer]]);sealed<-dsVert:::.callMpcTool("transport-encrypt",list(data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(and_xma=pr[[server]]$and_xma,and_ymb=pr[[server]]$and_ymb,had1_xma=pr[[server]]$had1_xma,had1_ymb=pr[[server]]$had1_ymb),auto_unbox=TRUE))),recipient_pk=pk_b64));.sendBlob(.to_b64url(sealed$sealed),"k2_peer_beaver_r1",peer_ci)}}else if(ph==3){for(server in server_list){peer<-setdiff(server_list,server);peer_ci<-which(server_names==peer);pk_b64<-.b64url_to_b64(transport_pks[[peer]]);sealed<-dsVert:::.callMpcTool("transport-encrypt",list(data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(had2_xma=pr[[server]]$had2_xma,had2_ymb=pr[[server]]$had2_ymb),auto_unbox=TRUE))),recipient_pk=pk_b64));.sendBlob(.to_b64url(sealed$sealed),"k2_peer_had2_r1",peer_ci)}}}}
    if (isTRUE(weights_active)) {
      .glm_apply_shared_weight_residual(
        datasources = datasources,
        dcf_parties = c(coordinator, nl),
        dcf_conns = c(coordinator_conn, nl_conn),
        dealer_conn = dealer_conn,
        transport_pks = transport_pks,
        session_id = session_id,
        n_obs = n_obs,
        .dsAgg = .dsAgg,
        .sendBlob = .sendBlob,
        ring = ring)
    }
    gtb<-.dsAgg(datasources[dealer_conn],call(name = "glmRing63GenGradTriplesDS",dcf0_pk=transport_pks[[coordinator]],dcf1_pk=transport_pks[[nl]],n=as.integer(n_obs),p=as.integer(p_total),ring=ring,session_id=session_id));if(is.list(gtb))gtb<-gtb[[1]];.sendBlob(gtb$grad_blob_0,"k2_grad_triple_fp",coordinator_conn);.sendBlob(gtb$grad_blob_1,"k2_grad_triple_fp",nl_conn)
    r1b<-list();for(server in server_list){ci<-which(server_names==server);peer<-setdiff(server_list,server);.dsAgg(datasources[ci],call(name = "k2StoreGradTripleDS",session_id=session_id));r<-.dsAgg(datasources[ci],call(name = "k2GradientR1DS",peer_pk=transport_pks[[peer]],session_id=session_id));if(is.list(r)&&length(r)==1)r<-r[[1]];r1b[[server]]<-r}
    .sendBlob(r1b[[coordinator]]$encrypted_r1,"k2_grad_peer_r1",nl_conn);.sendBlob(r1b[[nl]]$encrypted_r1,"k2_grad_peer_r1",coordinator_conn)
    r2b<-list();for(server in server_list){ci<-which(server_names==server);is_coord<-(server==coordinator);r<-.dsAgg(datasources[ci],call(name = "k2GradientR2DS",party_id=if(is_coord)0L else 1L,session_id=session_id));if(is.list(r)&&length(r)==1)r<-r[[1]];r2b[[server]]<-r}
    agb<-dsVert:::.callMpcTool("k2-ring63-aggregate",list(share_a=r2b[[coordinator]]$gradient_fp,share_b=r2b[[nl]]$gradient_fp,frac_bits=frac_bits,ring=ring_tag))
    arb<-dsVert:::.callMpcTool("k2-ring63-aggregate",list(share_a=r1b[[coordinator]]$sum_residual_fp,share_b=r1b[[nl]]$sum_residual_fp,frac_bits=frac_bits,ring=ring_tag))
    gm<-c(arb$values[1]/n_obs,agb$values/n_obs)+lambda*th_m
    # Central difference
    hessian_k2[,jj] <- (grad_fwd - gm) / (2 * delta_se)
    if(verbose) message(sprintf("    [SE] Column %d/%d",jj,p_plus1))
  }
  hessian_k2 <- (hessian_k2 + t(hessian_k2)) / 2
  # Column order is [Intercept, x_vars[[coordinator]], x_vars[[nl]]]
  # per theta_conv layout above. Attach dimnames so downstream consumers
  # (ds.vertLASSOProximal Gram reconstruction, Wald contrasts) can
  # permute to match external x_means / x_sds ordering (which may
  # differ from server-partition order).
  hess_names <- c("(Intercept)",
                  x_vars[[coordinator]],
                  if (p_nl > 0L) x_vars[[nl]] else character(0))
  dimnames(hessian_k2) <- list(hess_names, hess_names)
  attr(inv_hessian, "raw_hessian") <- hessian_k2
  }

  if (isTRUE(weights_active) && family != "gaussian") {
    if (verbose) {
      message("  [Deviance] Skipped for weighted non-Gaussian K=2 fit; canonical weighted deviance is not implemented yet.")
    }
    betas <- list()
    betas[[coordinator]] <- if (p_coord > 0L) beta[seq_len(p_coord)] else numeric(0)
    betas[[nl]] <- if (p_total > p_coord) beta[(p_coord + 1L):p_total] else numeric(0)
    return(list(betas=betas, intercept=intercept, converged=converged,
                iterations=final_iter, max_diff=max_diff, deviance=NA_real_,
                inv_hessian=inv_hessian,
                gradient_std=last_gradient))
  }

  if (!isTRUE(compute_deviance)) {
    betas <- list()
    betas[[coordinator]] <- if (p_coord > 0L) beta[seq_len(p_coord)] else numeric(0)
    betas[[nl]] <- if (p_total > p_coord) beta[(p_coord + 1L):p_total] else numeric(0)
    return(list(betas=betas, intercept=intercept, converged=converged,
                iterations=final_iter, max_diff=max_diff, deviance=NA_real_,
                inv_hessian=inv_hessian,
                gradient_std=last_gradient))
  }

  # === Secure canonical deviance ===
  if (verbose) message("  [Deviance] Computing canonical deviance...")

  # K=2 party aliases for spline phases
  dcf_parties <- c(coordinator, nl)
  dcf_conns <- c(coordinator_conn, nl_conn)

  # Recompute eta from converged beta (SE computation may have overwritten shares)
  # Guard zero-length peer slice (sleepstudy-style p_nl == 0 case).
  for (s in server_list) {
    ci <- which(server_names == s); is_coord <- (s == coordinator)
    b_coord <- if (p_coord > 0L) beta[seq_len(p_coord)] else numeric(0)
    b_nl <- if (p_total > p_coord) beta[(p_coord + 1L):p_total] else numeric(0)
    .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
      beta_coord = b_coord, beta_nl = b_nl, intercept = intercept,
      is_coordinator = is_coord, session_id = session_id))
  }

  # Helper: run one Beaver dot-product (nx1) and return aggregated scalar
  .beaver_dot <- function() {
    dt <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenGradTriplesDS",
           dcf0_pk = transport_pks[[coordinator]],
           dcf1_pk = transport_pks[[nl]],
           n = as.integer(n_obs), p = 1L,
           ring = ring,
           session_id = session_id))
    if (is.list(dt)) dt <- dt[[1]]
    .sendBlob(dt$grad_blob_0, "k2_grad_triple_fp", coordinator_conn)
    .sendBlob(dt$grad_blob_1, "k2_grad_triple_fp", nl_conn)
    dr1 <- list()
    for (s in server_list) {
      ci <- which(server_names == s); peer <- setdiff(server_list, s)
      .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; dr1[[s]] <- r
    }
    .sendBlob(dr1[[coordinator]]$encrypted_r1, "k2_grad_peer_r1", nl_conn)
    .sendBlob(dr1[[nl]]$encrypted_r1, "k2_grad_peer_r1", coordinator_conn)
    dr2 <- list()
    for (s in server_list) {
      ci <- which(server_names == s); is_c <- (s == coordinator)
      r <- .dsAgg(datasources[ci], call(name = "k2GradientR2DS",
        party_id = if(is_c) 0L else 1L, session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; dr2[[s]] <- r
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = dr2[[coordinator]]$gradient_fp,
      share_b = dr2[[nl]]$gradient_fp, frac_bits = frac_bits, ring = ring_tag))
    agg$values[1]
  }

  if (family == "gaussian") {
    # Gaussian: RSS = canonical deviance. For weighted Gaussian fits,
    # compute shares of sqrt(w) * residual before the RSS Beaver dot.
    if (isTRUE(weights_active)) {
      .glm_apply_shared_weight_residual(
        datasources = datasources,
        dcf_parties = c(coordinator, nl),
        dcf_conns = c(coordinator_conn, nl_conn),
        dealer_conn = dealer_conn,
        transport_pks = transport_pks,
        session_id = session_id,
        n_obs = n_obs,
        .dsAgg = .dsAgg,
        .sendBlob = .sendBlob,
        weight_key = "k2_sqrt_weights_share_fp",
        output_key = "k2_sqrt_weighted_residual_share_fp",
        ring = ring)
    }
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci], call(name = "glmRing63PrepDevianceDS",
        mode = "rss", session_id = session_id))
    }
    k2_deviance <- .beaver_dot()

  } else if (family == "binomial") {
    # Binomial: D = 2*(Sumsoftplus(eta) - y^T*eta)
    # Step 1: softplus spline on eta
    sp_dcf <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenDcfKeysDS",
           dcf0_pk = transport_pks[[coordinator]],
           dcf1_pk = transport_pks[[nl]],
           family = "softplus", n = as.integer(n_obs),
           frac_bits = as.integer(frac_bits),
           num_intervals = 80L, ring = ring, session_id = session_id))
    if (is.list(sp_dcf)) sp_dcf <- sp_dcf[[1]]
    .sendBlob(sp_dcf$dcf_blob_0, "k2_dcf_keys_persistent", coordinator_conn)
    .sendBlob(sp_dcf$dcf_blob_1, "k2_dcf_keys_persistent", nl_conn)
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci], call(name = "k2StoreDcfKeysPersistentDS", session_id = session_id))
    }
    sp_triples <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenSplineTriplesDS",
           dcf0_pk = transport_pks[[coordinator]],
           dcf1_pk = transport_pks[[nl]],
           n = as.integer(n_obs), frac_bits = as.integer(frac_bits),
           ring = ring, session_id = session_id))
    if (is.list(sp_triples)) sp_triples <- sp_triples[[1]]
    .sendBlob(sp_triples$spline_blob_0, "k2_spline_triples", coordinator_conn)
    .sendBlob(sp_triples$spline_blob_1, "k2_spline_triples", nl_conn)
    # Run 4-phase softplus spline (same as sigmoid/exp but family="softplus")
    for (ph in 1:4) {
      pr <- list()
      for (di in seq_along(dcf_parties)) {
        ci <- dcf_conns[di]
        r <- .dsAgg(datasources[ci], call(paste0("k2WideSplinePhase", ph, "DS"),
          party_id = as.integer(di - 1), family = "softplus",
          num_intervals = 80L, frac_bits = as.integer(frac_bits),
          ring = ring, session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; pr[[di]] <- r
      }
      if (ph == 1) {
        # Send dcf_masked RAW (same as iteration loop -- Phase 2 expects raw, not encrypted)
        .sendBlob(pr[[1]]$dcf_masked, "k2_peer_dcf_masked", dcf_conns[2])
        .sendBlob(pr[[2]]$dcf_masked, "k2_peer_dcf_masked", dcf_conns[1])
      } else if (ph == 2) {
        for (di in 1:2) {
          pi2 <- 3 - di; pk <- .b64url_to_b64(transport_pks[[dcf_parties[pi2]]])
          s <- dsVert:::.callMpcTool("transport-encrypt", list(
            data = jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(
              and_xma=pr[[di]]$and_xma, and_ymb=pr[[di]]$and_ymb,
              had1_xma=pr[[di]]$had1_xma, had1_ymb=pr[[di]]$had1_ymb),
              auto_unbox=TRUE))), recipient_pk=pk))
          .sendBlob(.to_b64url(s$sealed), "k2_peer_beaver_r1", dcf_conns[pi2])
        }
      } else if (ph == 3) {
        for (di in 1:2) {
          pi2 <- 3 - di; pk <- .b64url_to_b64(transport_pks[[dcf_parties[pi2]]])
          s <- dsVert:::.callMpcTool("transport-encrypt", list(
            data = jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(
              had2_xma=pr[[di]]$had2_xma, had2_ymb=pr[[di]]$had2_ymb),
              auto_unbox=TRUE))), recipient_pk=pk))
          .sendBlob(.to_b64url(s$sealed), "k2_peer_had2_r1", dcf_conns[pi2])
        }
      }
    }
    # Step 2: get Sumsoftplus from both parties
    sums <- list()
    for (s in server_list) {
      ci <- which(server_names == s)
      r <- .dsAgg(datasources[ci], call(name = "glmRing63DevianceSumsDS",
        family = "binomial", session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; sums[[s]] <- r
    }
    sum_sp_agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = sums[[coordinator]]$sum_fp, share_b = sums[[nl]]$sum_fp,
      frac_bits = frac_bits, ring = ring_tag))
    sum_softplus <- sum_sp_agg$values[1]
    # Step 3: Beaver y^T*eta
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci], call(name = "glmRing63PrepDevianceDS",
        mode = "canonical", session_id = session_id))
    }
    y_dot_eta <- .beaver_dot()
    k2_deviance <- 2 * (sum_softplus - y_dot_eta)

  } else {
    # Poisson: D = 2*(Summu - y^T*eta + C) where C = Sum(y*log(y) - y)
    sums <- list()
    for (s in server_list) {
      ci <- which(server_names == s)
      r <- .dsAgg(datasources[ci], call(name = "glmRing63DevianceSumsDS",
        family = "poisson", session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; sums[[s]] <- r
    }
    sum_mu_agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = sums[[coordinator]]$sum_fp, share_b = sums[[nl]]$sum_fp,
      frac_bits = frac_bits, ring = ring_tag))
    sum_mu <- sum_mu_agg$values[1]
    # Only the label server returns non-zero null_term; sum to pick it up
    null_term <- sum(sapply(sums, function(s) if(!is.null(s$null_term)) s$null_term else 0))
    # Beaver y^T*eta
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci], call(name = "glmRing63PrepDevianceDS",
        mode = "canonical", session_id = session_id))
    }
    y_dot_eta <- .beaver_dot()
    k2_deviance <- 2 * (sum_mu - y_dot_eta + null_term)
  }

  if (verbose) message(sprintf("  [Deviance] = %.4f", k2_deviance))

  betas <- list()
  betas[[coordinator]] <- if (p_coord > 0L) beta[seq_len(p_coord)] else numeric(0)
  betas[[nl]] <- if (p_total > p_coord) beta[(p_coord + 1L):p_total] else numeric(0)
  list(betas=betas, intercept=intercept, converged=converged,
       iterations=final_iter, max_diff=max_diff, deviance=k2_deviance,
       inv_hessian=inv_hessian,
       gradient_std=last_gradient)
}
