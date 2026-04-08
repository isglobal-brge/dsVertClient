#' @title K>=3 Secure Aggregation Loop
#' @description BCD iterations with pairwise PRG masks for K>=3 servers.
#'   Gaussian: one-shot pairwise Beaver X^T X + CKKS X^T y.
#'   Binomial/Poisson: L-BFGS with encrypted gradients + threshold decryption.
#' @name k3-secure-agg-loop
NULL

#' L-BFGS two-loop recursion (Nocedal & Wright, Algorithm 7.4)
#' Duplicated from ds.vertGLM.k2.R for standalone use.
#' @keywords internal
.lbfgs_direction_k3 <- function(grad, s_hist, y_hist) {
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
.k3_secure_agg_loop <- function(datasources, server_list, server_names, x_vars,
                                 coordinator, coordinator_conn, non_label_servers,
                                 transport_pks, std_data, y_var, family,
                                 betas, n_obs, lambda, log_n, log_scale,
                                 session_id, max_iter, tol, verbose,
                                 topology, label_intercept,
                                 .dsAgg, .sendBlob) {

  n_partitions <- length(server_list)
  non_label_count <- length(non_label_servers)
  coordinator_pk <- if (length(transport_pks) > 0) transport_pks[[coordinator]] else NULL

  converged <- FALSE
  final_iter <- 0

  non_label_pk_map <- list()
  for (s in non_label_servers) {
    non_label_pk_map[[s]] <- transport_pks[[s]]
  }

  # Sort non-label server names for canonical ordering in seed derivation
  nonlabel_sorted <- sort(non_label_servers)

  # Fusion server and non-fusion servers (for threshold decryption)
  fusion_server <- server_list[1]  # Party 0
  fusion_conn_idx <- which(server_names == fusion_server)
  non_fusion_servers <- setdiff(server_list, fusion_server)

  if (verbose) message(sprintf("\n[Phase 3] Secure aggregation (K=%d, family=%s, n=%d, lambda=%.1e)",
                                 n_partitions, family, n_obs, lambda))

  # Initialize FSM on coordinator
  .dsAgg(
    conns = datasources[coordinator_conn],
    expr = call("glmFSMInitDS",
                session_id = session_id,
                n_nonlabel = as.integer(non_label_count),
                mode = "secure_agg")
  )

  # Initialize secure aggregation on each non-label server
  for (server in non_label_servers) {
    conn_idx <- which(server_names == server)
    .dsAgg(
      conns = datasources[conn_idx],
      expr = call("glmSecureAggInitDS",
                  self_name = server,
                  session_id = session_id,
                  nonlabel_names = nonlabel_sorted,
                  scale_bits = 20L,
                  topology = topology)
    )
  }

  if (verbose) message(sprintf("  [Init] FSM + secure aggregation initialized (%d non-label servers, topology=%s)",
                                 non_label_count, topology))

  # === GAUSSIAN ONE-SHOT for K>=3 ===
  if (family == "gaussian") {
    t0_oneshot <- proc.time()[[3]]
    if (verbose) message("  [One-Shot] Gaussian: pairwise Beaver X^T X + X^T y")
    frac_bits_k3 <- 20L

    # Step 1: Each server computes local X_k^T X_k (and X_label^T y for label)
    local_results <- list()
    for (server in server_list) {
      conn_idx <- which(server_names == server)
      r <- .dsAgg(datasources[conn_idx],
        call("glmGaussianLocalXtXDS", data_name = std_data,
             x_vars = x_vars[[server]],
             y_var = if (server == coordinator) y_var else NULL,
             session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      local_results[[server]] <- r
    }
    xty_all <- list()
    xty_all[[coordinator]] <- local_results[[coordinator]]$xty_local

    # Step 2: Pairwise Beaver for X_j^T X_k cross-blocks AND X_k^T y
    # Each pair uses the k2-gaussian-oneshot Go command with zero-padded features
    pair_list <- list()
    for (i in seq_along(server_list)) {
      for (j in seq_along(server_list)) {
        if (i < j) pair_list[[length(pair_list)+1]] <- c(server_list[i], server_list[j])
      }
    }

    cross_blocks <- list()
    for (pair in pair_list) {
      s_a <- pair[1]; s_b <- pair[2]
      p_a <- length(x_vars[[s_a]]); p_b <- length(x_vars[[s_b]])
      p_combined <- p_a + p_b
      ci_a <- which(server_names == s_a); ci_b <- which(server_names == s_b)

      # Generate Beaver triples for the combined feature space
      n_elems <- as.integer((p_combined^2 + p_combined) * n_obs)
      cross_t <- dsVert:::.callMheTool("k2-gen-beaver-triples",
        list(n = n_elems, frac_bits = frac_bits_k3))

      # Distribute triples
      for (info in list(list(s=s_a, ci=ci_a, pid=0L, party="party0"),
                        list(s=s_b, ci=ci_b, pid=1L, party="party1"))) {
        pk_b64 <- gsub("-","+",gsub("_","/",transport_pks[[info$s]],fixed=TRUE),fixed=TRUE)
        pad <- nchar(pk_b64)%%4; if(pad==2) pk_b64<-paste0(pk_b64,"=="); if(pad==3) pk_b64<-paste0(pk_b64,"=")
        td <- list(a=cross_t[[paste0(info$party,"_u")]],
                   b=cross_t[[paste0(info$party,"_v")]],
                   c=cross_t[[paste0(info$party,"_w")]])
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(td,auto_unbox=TRUE))),
          recipient_pk=pk_b64))
        .sendBlob(gsub("+","-",gsub("/","_",gsub("=+$","",sealed$sealed,perl=TRUE),fixed=TRUE),fixed=TRUE), "pairwise_cross_triples", info$ci)
      }

      # Phase 1: Beaver R1 (label server also provides y for X^T y computation)
      ph1 <- list()
      for (info in list(list(s=s_a, ci=ci_a, pid=0L, peer_p=p_b),
                        list(s=s_b, ci=ci_b, pid=1L, peer_p=p_a))) {
        yv <- if (info$s == coordinator) y_var else NULL
        r <- .dsAgg(datasources[info$ci],
          call("glmPairwiseCrossProductDS", data_name = std_data,
               x_vars = x_vars[[info$s]], peer_p = as.integer(info$peer_p),
               party_id = info$pid, phase = 1L, y_var = yv,
               frac_bits = frac_bits_k3, session_id = session_id))
        if (is.list(r) && length(r)==1) r <- r[[1]]
        ph1[[info$s]] <- r
      }

      # Relay Beaver R1
      for (info in list(list(from=s_a, to=s_b, ci_to=ci_b),
                        list(from=s_b, to=s_a, ci_to=ci_a))) {
        pk_b64 <- gsub("-","+",gsub("_","/",transport_pks[[info$to]],fixed=TRUE),fixed=TRUE)
        pad <- nchar(pk_b64)%%4; if(pad==2) pk_b64<-paste0(pk_b64,"=="); if(pad==3) pk_b64<-paste0(pk_b64,"=")
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(
            list(xma=ph1[[info$from]]$xma, ymb=ph1[[info$from]]$ymb),auto_unbox=TRUE))),
          recipient_pk=pk_b64))
        .sendBlob(gsub("+","-",gsub("/","_",gsub("=+$","",sealed$sealed,perl=TRUE),fixed=TRUE),fixed=TRUE), "pairwise_cross_peer_r1", info$ci_to)
      }

      # Phase 2: Beaver close
      ph2 <- list()
      for (info in list(list(s=s_a, ci=ci_a, pid=0L, peer_p=p_b),
                        list(s=s_b, ci=ci_b, pid=1L, peer_p=p_a))) {
        yv <- if (info$s == coordinator) y_var else NULL
        r <- .dsAgg(datasources[info$ci],
          call("glmPairwiseCrossProductDS", data_name = std_data,
               x_vars = x_vars[[info$s]], peer_p = as.integer(info$peer_p),
               party_id = info$pid, phase = 2L, y_var = yv,
               frac_bits = frac_bits_k3, session_id = session_id))
        if (is.list(r) && length(r)==1) r <- r[[1]]
        ph2[[info$s]] <- r
      }

      # Reconstruct cross-product block (X^T X off-diagonal)
      cross_agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a = ph2[[s_a]]$cross_xtx_fp,
        share_b = ph2[[s_b]]$cross_xtx_fp,
        frac_bits = frac_bits_k3))
      full_mat <- matrix(cross_agg$values, p_combined, p_combined, byrow = TRUE)
      cross_block <- full_mat[1:p_a, (p_a+1):p_combined, drop = FALSE]
      cross_blocks[[paste0(s_a, "_", s_b)]] <- cross_block

      # Extract X^T y from pairs involving the label server
      if (s_a == coordinator || s_b == coordinator) {
        xty_agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
          share_a = ph2[[s_a]]$cross_xty_fp,
          share_b = ph2[[s_b]]$cross_xty_fp,
          frac_bits = frac_bits_k3))
        xty_vec <- xty_agg$values  # length p_a + p_b
        # The non-label server's X^T y is in the first p_nonlabel positions
        if (s_a == coordinator) {
          # s_b is non-label (party 0), its X^T y is positions 1:p_b
          xty_all[[s_b]] <- xty_vec[1:p_b]
        } else {
          # s_a is non-label (party 0), its X^T y is positions 1:p_a
          xty_all[[s_a]] <- xty_vec[1:p_a]
        }
      }
    }

    # Step 4: Assemble full X^T X matrix
    p_total_k3 <- sum(sapply(server_list, function(s) length(x_vars[[s]])))
    XtX <- matrix(0, p_total_k3, p_total_k3)
    XtY <- numeric(p_total_k3)
    idx <- 1
    for (server in server_list) {
      p_s <- length(x_vars[[server]])
      # Diagonal block
      XtX[idx:(idx+p_s-1), idx:(idx+p_s-1)] <- matrix(
        local_results[[server]]$xtx_local, p_s, p_s)
      # X^T y block
      XtY[idx:(idx+p_s-1)] <- xty_all[[server]]
      idx <- idx + p_s
    }
    # Off-diagonal blocks (symmetric)
    idx_i <- 1
    for (i in seq_along(server_list)) {
      p_i <- length(x_vars[[server_list[i]]])
      idx_j <- idx_i + p_i
      for (j in (i+1):min(length(server_list), length(server_list))) {
        if (j > length(server_list)) break
        p_j <- length(x_vars[[server_list[j]]])
        key <- paste0(server_list[i], "_", server_list[j])
        if (!is.null(cross_blocks[[key]])) {
          XtX[idx_i:(idx_i+p_i-1), idx_j:(idx_j+p_j-1)] <- cross_blocks[[key]]
          XtX[idx_j:(idx_j+p_j-1), idx_i:(idx_i+p_i-1)] <- t(cross_blocks[[key]])
        }
        idx_j <- idx_j + p_j
      }
      idx_i <- idx_i + p_i
    }

    # Step 5: Solve on client
    XtX_reg <- XtX / n_obs + lambda * diag(p_total_k3)
    XtY_norm <- XtY / n_obs
    beta_solved <- as.numeric(solve(XtX_reg, XtY_norm))

    if (verbose) {
      message(sprintf("  [One-Shot] X^T X: %dx%d, %d pairwise rounds, cond=%.1e (%.1fs)",
        p_total_k3, p_total_k3, length(pair_list), kappa(XtX_reg),
        proc.time()[[3]] - t0_oneshot))
      message(sprintf("  [One-Shot] ||beta||=%.4f, range=[%.4f, %.4f]",
        sqrt(sum(beta_solved^2)), min(beta_solved), max(beta_solved)))
    }

    # Distribute betas
    idx <- 1
    for (server in server_list) {
      p_s <- length(x_vars[[server]])
      betas[[server]] <- beta_solved[idx:(idx+p_s-1)]
      idx <- idx + p_s
    }
    converged <- TRUE; final_iter <- 0

  } else {
  # === L-BFGS for Binomial/Poisson K>=3 ===
  if (verbose) message(sprintf("  [L-BFGS] %s, %d non-label servers, Enc(r) threshold decrypt",
                                 family, non_label_count))
  t0_loop <- proc.time()[[3]]

  # Encrypted eta blobs (opaque to client)
  encrypted_etas <- list()

  # L-BFGS state (client-side only)
  lbfgs_s <- list(); lbfgs_y <- list(); lbfgs_prev_theta <- NULL; lbfgs_prev_grad <- NULL

  for (iter in seq_len(max_iter)) {
    t0_iter <- proc.time()[[3]]
    betas_old <- betas
    max_diff <- 0

    # --- Phase 3a: FSM check + Coordinator step (label server) ---
    # Send masked eta blobs to coordinator (from previous iteration)
    eta_blob_keys <- character(0)
    if (length(encrypted_etas) > 0) {
      for (s in names(encrypted_etas)) {
        blob <- encrypted_etas[[s]]
        if (!is.null(blob) && nzchar(blob)) {
          key <- paste0("eta_", s)
          .sendBlob(blob, key, coordinator_conn)
          eta_blob_keys <- c(eta_blob_keys, key)

          # FSM: register eta receipt
          .dsAgg(
            conns = datasources[coordinator_conn],
            expr = call("glmFSMCheckDS",
                        session_id = session_id,
                        action = "receive_eta",
                        server_name = s)
          )
        }
      }
    }

    # FSM: authorize coordinator step
    .dsAgg(
      conns = datasources[coordinator_conn],
      expr = call("glmFSMCheckDS",
                  session_id = session_id,
                  action = "coordinator_step",
                  iteration = as.integer(iter))
    )

    coord_result <- .dsAgg(
      conns = datasources[coordinator_conn],
      expr = call("glmSecureAggCoordinatorStepDS",
                  data_name = std_data,
                  y_var = y_var,
                  x_vars = x_vars[[coordinator]],
                  eta_blob_keys = if (length(eta_blob_keys) > 0) eta_blob_keys else NULL,
                  non_label_pks = non_label_pk_map,
                  family = family,
                  beta_current = betas[[coordinator]],
                  lambda = lambda,
                  intercept = label_intercept,
                  n_obs = as.integer(n_obs),
                  scale_bits = 20L,
                  session_id = session_id)
    )
    if (is.list(coord_result) && length(coord_result) == 1)
      coord_result <- coord_result[[1]]

    betas[[coordinator]] <- coord_result$beta
    mwv_blobs <- coord_result$encrypted_blobs

    diff_coord <- sum(abs(betas[[coordinator]] - betas_old[[coordinator]]))
    max_diff <- max(max_diff, diff_coord)

    # FSM: distribute_mwv
    .dsAgg(
      conns = datasources[coordinator_conn],
      expr = call("glmFSMCheckDS",
                  session_id = session_id,
                  action = "distribute_mwv")
    )

    # --- Phase 3b: Non-label servers (gradient computation + threshold decrypt) ---
    encrypted_etas <- list()
    nl_gradients <- list()

    # Send Enc(r) or (mu,w) to non-label servers for gradient computation
    enc_r <- coord_result$encrypted_residual
    use_enc_r <- !is.null(enc_r)
    if (use_enc_r) {
      # Transfer Enc(r) via chunked send — non-label servers never see mu or w
      for (server in non_label_servers) {
        ci <- which(server_names == server)
        nc <- .dsvert_adaptive_send(enc_r, function(chunk_str, chunk_idx, n_chunks) {
          .dsAgg(conns = datasources[ci],
            expr = call("mheStoreEncChunkDS", col_index = 2L,
                        chunk_index = chunk_idx, chunk = chunk_str,
                        session_id = session_id))
        })
        .dsAgg(conns = datasources[ci],
          expr = call("mheAssembleEncColumnDS", col_index = 2L,
                      n_chunks = as.integer(nc), session_id = session_id))
      }
    }

    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      vars <- x_vars[[server]]
      p_k <- length(vars)

      if (!use_enc_r) {
        .sendBlob(mwv_blobs[[server]], "mwv", conn_idx)
      }

      # Compute encrypted gradient
      grad_result <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("glmSecureGradientDS",
                    data_name = std_data,
                    x_vars = vars,
                    encrypted_mwv = NULL,
                    num_obs = as.integer(n_obs),
                    use_enc_residual = use_enc_r,
                    session_id = session_id)
      )
      if (is.list(grad_result) && length(grad_result) == 1)
        grad_result <- grad_result[[1]]

      enc_gradients <- grad_result$encrypted_gradients
      ct_hashes <- grad_result$ct_hashes

      # Protocol Firewall: authorize gradient CTs (REUSE existing pattern)
      if (!is.null(ct_hashes) && length(ct_hashes) > 0) {
        ct_hashes_blob <- paste(ct_hashes, collapse = ",")
        for (auth_server in server_list) {
          if (auth_server != server) {
            auth_conn <- which(server_names == auth_server)
            .sendBlob(ct_hashes_blob, "ct_hashes", auth_conn)
            .dsAgg(
              conns = datasources[auth_conn],
              expr = call("mheAuthorizeCTDS",
                          op_type = "glm-gradient",
                          from_storage = TRUE,
                          session_id = session_id)
            )
          }
        }
      }

      # Step 3: Batched threshold decryption of gradient components
      # Send all p_k gradient CTs as batch blobs, then one round-trip per
      # server instead of O(p_k * K) individual round-trips.
      gradient <- numeric(p_k)

      # Send all gradient CTs to each non-fusion server
      for (nf_server in non_fusion_servers) {
        nf_conn <- which(server_names == nf_server)
        nf_party_id <- which(server_list == nf_server) - 1
        for (j in seq_len(p_k)) {
          .sendBlob(enc_gradients[j], paste0("ct_batch_", j), nf_conn)
        }
        pd <- .dsAgg(
          conns = datasources[nf_conn],
          expr = call("mhePartialDecryptBatchWrappedDS",
                      n_cts = as.integer(p_k),
                      session_id = session_id)
        )
        if (is.list(pd)) pd <- pd[[1]]
        for (j in seq_len(p_k)) {
          share_key <- paste0("wrapped_share_", nf_party_id, "_ct_", j)
          .sendBlob(pd$wrapped_shares[j], share_key, fusion_conn_idx)
        }
      }

      # Send all gradient CTs to fusion server and batch fuse
      for (j in seq_len(p_k)) {
        .sendBlob(enc_gradients[j], paste0("ct_batch_", j), fusion_conn_idx)
      }
      fuse_result <- .dsAgg(
        conns = datasources[fusion_conn_idx],
        expr = call("mheFuseBatchDS",
                    n_cts = as.integer(p_k),
                    n_parties = as.integer(length(server_list)),
                    num_slots = 0L,
                    session_id = session_id)
      )
      if (is.list(fuse_result)) fuse_result <- fuse_result[[1]]
      nl_gradients[[server]] <- fuse_result$values
    }

    # --- Phase 3c: L-BFGS update (client-side) ---
    # Assemble full gradient: [coordinator block, non-label blocks]
    # All gradients are X^T(y-mu) (unnormalized). Divide by n for average.
    coord_grad <- coord_result$gradient_label
    if (is.null(coord_grad)) coord_grad <- rep(0, length(x_vars[[coordinator]]))

    raw_gradient <- coord_grad  # X_label^T(y-mu), unnormalized
    for (server in non_label_servers) {
      raw_gradient <- c(raw_gradient, nl_gradients[[server]])
    }

    # Gradient of NLL = -X^T(y-mu)/n + lambda*theta = X^T(mu-y)/n + lambda*theta
    theta <- unlist(betas[c(coordinator, non_label_servers)])
    full_grad <- -raw_gradient / n_obs + lambda * theta

    # Update L-BFGS history
    if (!is.null(lbfgs_prev_theta)) {
      sk <- theta - lbfgs_prev_theta
      yk <- full_grad - lbfgs_prev_grad
      if (sum(sk * yk) > 1e-10) {
        lbfgs_s <- c(lbfgs_s, list(sk)); lbfgs_y <- c(lbfgs_y, list(yk))
        if (length(lbfgs_s) > 7) { lbfgs_s <- lbfgs_s[-1]; lbfgs_y <- lbfgs_y[-1] }
      }
    }
    lbfgs_prev_theta <- theta; lbfgs_prev_grad <- full_grad

    # L-BFGS direction
    direction <- .lbfgs_direction_k3(full_grad, lbfgs_s, lbfgs_y)
    step_size <- if (iter <= 1) 0.3 else 1.0
    new_theta <- theta + step_size * direction

    # Distribute new betas to servers (coordinator may include intercept)
    idx <- 1
    for (server in c(coordinator, non_label_servers)) {
      p_s <- length(betas[[server]])  # use actual beta length (includes intercept if present)
      betas[[server]] <- new_theta[idx:(idx + p_s - 1)]
      idx <- idx + p_s
    }

    if (verbose && iter <= 3)
      message(sprintf("  [L-BFGS] ||grad||=%.6f, grad_range=[%.4f, %.4f]",
                       sqrt(sum(full_grad^2)), min(full_grad), max(full_grad)))

    # --- Phase 3d: Non-label servers: compute masked eta with L-BFGS beta ---
    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      vars <- x_vars[[server]]

      # Send mwv blob if available (skip_solve will clean it up)
      if (!is.null(mwv_blobs[[server]])) {
        .sendBlob(mwv_blobs[[server]], "mwv", conn_idx)
      }

      solve_result <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("glmSecureAggBlockSolveDS",
                    data_name = std_data,
                    x_vars = vars,
                    beta_current = betas[[server]],
                    gradient = rep(0, length(vars)),
                    lambda = lambda,
                    coordinator_pk = coordinator_pk,
                    iteration = as.integer(iter),
                    skip_solve = TRUE,
                    session_id = session_id)
      )
      if (is.list(solve_result) && length(solve_result) == 1)
        solve_result <- solve_result[[1]]

      encrypted_etas[[server]] <- solve_result$encrypted_masked_eta

      # FSM: block complete
      .dsAgg(
        conns = datasources[coordinator_conn],
        expr = call("glmFSMCheckDS",
                    session_id = session_id,
                    action = "block_complete")
      )
    }

    # Convergence check
    max_diff <- max(abs(new_theta - unlist(betas_old[c(coordinator, non_label_servers)])))
    final_iter <- iter

    if (iter %% 20 == 0) {
      for (server in server_list) {
        tryCatch(.dsAgg(conns = datasources[which(server_names == server)],
          expr = call("mheGcDS")), error = function(e) NULL)
      }
    }

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose)
        message(sprintf("  Converged after %d iterations (diff = %.2e)",
                        iter, max_diff))
      break
    }

    if (verbose)
      message(sprintf("  Iter %d: ||grad||=%.4f, step=%.2f, diff=%.2e, theta=[%.3f, %.3f] (%.1fs)",
                       iter, sqrt(sum(full_grad^2)), step_size, max_diff,
                       min(new_theta), max(new_theta), proc.time()[[3]] - t0_iter))
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)",
                    max_iter, max_diff))

  if (verbose) message(sprintf("  [Secure-Agg] %s after %d iters (total %.1fs)",
                                 if (converged) "Converged" else "Stopped",
                                 final_iter, proc.time()[[3]] - t0_loop))
  }  # end else (binomial/poisson L-BFGS)

  return(list(betas = betas, converged = converged, final_iter = final_iter))
}
