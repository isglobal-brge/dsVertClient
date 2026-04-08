#' @title K>=3 Polynomial Sigmoid Loop (Zero Observation-Level Disclosure)
#' @description CKKS polynomial sigmoid for K>=3 binomial/Poisson GLM.
#'   Computes Enc(mu) = poly_sigmoid(Enc(eta_total)) entirely in CKKS.
#'   Nobody sees eta_total, mu, or y-mu in plaintext. Only aggregate
#'   gradient scalars (one per feature) are decrypted.
#'
#' @details Per iteration:
#' \enumerate{
#'   \item Each server encrypts X_k * beta_k under CPK
#'   \item Coordinator sums encrypted etas + evaluates degree-7 polynomial sigmoid
#'   \item Enc(mu) distributed to all servers
#'   \item Each server: encrypted gradient X_k^T * (Enc(y) - Enc(mu))
#'   \item Threshold decrypt gradients (aggregate scalars only)
#'   \item Client L-BFGS update
#' }
#'
#' Requires: log_n >= 13 (5 multiplicative levels), RLK, Galois keys.
#' @name k3-poly-sigmoid-loop
NULL

.lbfgs_direction_poly <- function(grad, s_hist, y_hist) {
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
.k3_poly_sigmoid_loop <- function(
    datasources, server_list, server_names, x_vars,
    coordinator, coordinator_conn, non_label_servers,
    transport_pks, std_data, y_var, family,
    betas, n_obs, lambda, log_n, log_scale,
    session_id, max_iter, tol, verbose,
    label_intercept, .dsAgg, .sendBlob) {

  fusion_server <- server_list[1]
  fusion_conn <- which(server_names == fusion_server)
  non_fusion_servers <- setdiff(server_list, fusion_server)

  # Polynomial coefficients for the link function
  poly_coefficients <- if (family == "poisson") {
    # Degree-4 Taylor for exp(x) on [-3, 3]
    c(1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0, 0.0, 0.0, 0.0)
  } else {
    NULL  # NULL = use built-in degree-7 sigmoid in glmHELinkStepDS
  }

  # L-BFGS state
  lbfgs_s <- list(); lbfgs_y <- list()
  prev_theta <- NULL; prev_grad <- NULL
  converged <- FALSE; final_iter <- 0

  if (verbose) message(sprintf("\n[Phase 3] K>=3 poly-sigmoid L-BFGS (family=%s, K=%d, n=%d, lambda=%.1e)",
                                 family, length(server_list), n_obs, lambda))
  t0_loop <- proc.time()[[3]]

  for (iter in seq_len(max_iter)) {
    t0_iter <- proc.time()[[3]]
    betas_old <- betas

    # =================================================================
    # Step 1: Encrypt etas + sum + polynomial sigmoid → Enc(mu)
    # =================================================================
    if (verbose) message(sprintf("    [%d.1] Encrypting etas + polynomial sigmoid...", iter))

    for (k in seq_along(server_list)) {
      server <- server_list[k]
      ci <- which(server_names == server)
      beta_for_encrypt <- betas[[server]]
      if (server == coordinator && label_intercept &&
          length(betas[[server]]) == length(x_vars[[server]]) + 1) {
        beta_for_encrypt <- betas[[server]][-1]  # strip intercept for X*beta
      }
      enc_result <- .dsAgg(datasources[ci],
        call("glmHEEncryptEtaDS", data_name = std_data,
             x_vars = x_vars[[server]], beta = beta_for_encrypt,
             clip_radius = 8.0 / length(server_list),  # adaptive: domain / K
             session_id = session_id))
      if (is.list(enc_result)) enc_result <- enc_result[[1]]
      .sendBlob(enc_result$encrypted_eta, paste0("ct_eta_", k - 1), coordinator_conn)
    }

    # Coordinator: sum etas + add intercept + polynomial sigmoid
    mws_intercept <- 0
    if (label_intercept && length(betas[[coordinator]]) > length(x_vars[[coordinator]]))
      mws_intercept <- betas[[coordinator]][1]

    link_result <- .dsAgg(datasources[coordinator_conn],
      call("glmHELinkStepDS", from_storage = TRUE,
           n_parties = as.integer(length(server_list)),
           poly_coefficients = poly_coefficients,
           skip_poly = FALSE,
           intercept = mws_intercept,
           session_id = session_id))
    if (is.list(link_result)) link_result <- link_result[[1]]

    if (verbose) message(sprintf("    [%d.2] Distributing Enc(mu)...", iter))

    # Distribute Enc(mu) to all servers
    ct_mu <- link_result$ct_mu
    for (server in server_list) {
      ci <- which(server_names == server)
      .sendBlob(ct_mu, "ct_mu", ci)
    }

    # =================================================================
    # Step 2: Compute ALL gradients, then ONE batch threshold decrypt
    # =================================================================
    if (verbose) message(sprintf("    [%d.3] Computing gradients...", iter))

    # Phase A: compute all encrypted gradients (sequential per server)
    all_enc_grads <- character(0)
    all_ct_hashes <- character(0)
    server_p_k <- integer(0)
    for (server in server_list) {
      ci <- which(server_names == server)
      use_intercept <- (server == coordinator && label_intercept &&
                        length(betas[[server]]) == length(x_vars[[server]]) + 1)
      grad_result <- .dsAgg(datasources[ci],
        call("glmHEGradientEncDS", data_name = std_data,
             x_vars = x_vars[[server]], num_obs = as.integer(n_obs),
             from_storage = TRUE, include_intercept = use_intercept,
             session_id = session_id))
      if (is.list(grad_result)) grad_result <- grad_result[[1]]
      all_enc_grads <- c(all_enc_grads, grad_result$encrypted_gradients)
      all_ct_hashes <- c(all_ct_hashes, grad_result$ct_hashes)
      server_p_k <- c(server_p_k, length(grad_result$encrypted_gradients))
    }
    total_cts <- length(all_enc_grads)

    # Phase B: ONE batch authorize + partial decrypt + fuse
    # Authorize ALL gradient CTs on non-fusion servers
    if (length(all_ct_hashes) > 0) {
      hashes_str <- paste(all_ct_hashes, collapse = ",")
      for (nf_server in non_fusion_servers) {
        nf_conn <- which(server_names == nf_server)
        .sendBlob(hashes_str, "ct_hashes", nf_conn)
        .dsAgg(datasources[nf_conn],
          call("mheAuthorizeCTDS", op_type = "he-link-gradient",
               from_storage = TRUE, session_id = session_id))
      }
      # Authorize on fusion too (for CTs from other servers)
      .sendBlob(hashes_str, "ct_hashes", fusion_conn)
      .dsAgg(datasources[fusion_conn],
        call("mheAuthorizeCTDS", op_type = "he-link-gradient",
             from_storage = TRUE, session_id = session_id))
    }

    # Send ALL CTs to non-fusion servers for partial decrypt
    for (nf_server in non_fusion_servers) {
      nf_conn <- which(server_names == nf_server)
      nf_party_id <- which(server_list == nf_server) - 1
      for (j in seq_len(total_cts))
        .sendBlob(all_enc_grads[j], paste0("ct_batch_", j), nf_conn)
      pd <- .dsAgg(datasources[nf_conn],
        call("mhePartialDecryptBatchWrappedDS",
             n_cts = as.integer(total_cts), session_id = session_id))
      if (is.list(pd)) pd <- pd[[1]]
      for (j in seq_len(total_cts))
        .sendBlob(pd$wrapped_shares[j],
                  paste0("wrapped_share_", nf_party_id, "_ct_", j), fusion_conn)
    }

    # Send ALL CTs to fusion + fuse
    for (j in seq_len(total_cts))
      .sendBlob(all_enc_grads[j], paste0("ct_batch_", j), fusion_conn)
    fuse_all <- .dsAgg(datasources[fusion_conn],
      call("mheFuseBatchDS", n_cts = as.integer(total_cts),
           n_parties = as.integer(length(server_list)),
           num_slots = 0L, session_id = session_id))
    if (is.list(fuse_all)) fuse_all <- fuse_all[[1]]

    # Split fused gradients back to per-server
    server_gradients <- list()
    idx <- 1
    for (i in seq_along(server_list)) {
      pk <- server_p_k[i]
      server_gradients[[server_list[i]]] <- fuse_all$values[idx:(idx + pk - 1)]
      idx <- idx + pk
    }

    # =================================================================
    # Step 3: L-BFGS update (client-side)
    # =================================================================
    raw_gradient <- numeric(0)
    for (server in server_list)
      raw_gradient <- c(raw_gradient, server_gradients[[server]])

    theta <- unlist(betas[server_list])
    full_grad <- -raw_gradient / n_obs + lambda * theta

    if (!is.null(prev_theta)) {
      sk <- theta - prev_theta; yk <- full_grad - prev_grad
      if (sum(sk * yk) > 1e-10) {
        lbfgs_s <- c(lbfgs_s, list(sk)); lbfgs_y <- c(lbfgs_y, list(yk))
        if (length(lbfgs_s) > 7) { lbfgs_s <- lbfgs_s[-1]; lbfgs_y <- lbfgs_y[-1] }
      }
    }
    prev_theta <- theta; prev_grad <- full_grad

    direction <- .lbfgs_direction_poly(full_grad, lbfgs_s, lbfgs_y)
    step_size <- if (iter <= 1) 0.3 else 1.0
    new_theta <- theta + step_size * direction
    max_diff <- max(abs(new_theta - unlist(betas_old[server_list])))

    idx <- 1
    for (server in server_list) {
      p_s <- length(betas[[server]])
      betas[[server]] <- new_theta[idx:(idx + p_s - 1)]
      idx <- idx + p_s
    }

    final_iter <- iter
    grad_norm <- sqrt(sum(full_grad^2))
    if (verbose) message(sprintf("  Iter %d: ||grad||=%.4f, step=%.2f, diff=%.2e, theta=[%.3f, %.3f] (%.1fs)",
                                   iter, grad_norm, step_size, max_diff,
                                   min(new_theta), max(new_theta),
                                   proc.time()[[3]] - t0_iter))

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose) message(sprintf("  Converged after %d iterations (diff = %.2e)", iter, max_diff))
      break
    }

    # Periodic garbage collection
    if (iter %% 10 == 0) {
      for (server in server_list)
        tryCatch(.dsAgg(datasources[which(server_names == server)],
          call("mheGcDS")), error = function(e) NULL)
    }
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)", max_iter, max_diff))

  if (verbose) message(sprintf("  [Poly-Sigmoid] %s after %d iters (total %.1fs)",
                                 if (converged) "Converged" else "Stopped",
                                 final_iter, proc.time()[[3]] - t0_loop))

  list(betas = betas, converged = converged, final_iter = final_iter)
}
