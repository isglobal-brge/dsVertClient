#' @title K>=3 Masked Wide Spline Loop
#' @description Non-disclosive K>=3 binomial/Poisson via masked wide spline.
#'   Combines CKKS aggregation (eta masking) with K=2 wide spline DCF
#'   (sigmoid/exp evaluation on shares). Nobody sees eta_total or mu.
#'
#' @details Per iteration:
#' \enumerate{
#'   \item CKKS encrypt eta_k, sum, mask with r → threshold decrypt
#'   \item K=2 wide spline DCF on shares (fusion + coordinator)
#'   \item Re-encrypt mu shares, sum → Enc(mu), gradient via InnerSum
#'   \item L-BFGS update on client
#' }
#'
#' @name k3-masked-wide-spline
NULL

.lbfgs_direction_mws <- function(grad, s_hist, y_hist) {
  m <- length(s_hist)
  if (m == 0L) return(-grad)
  q <- grad; alpha_v <- rho <- numeric(m)
  for (i in m:1) {
    rho[i] <- 1 / sum(y_hist[[i]] * s_hist[[i]])
    alpha_v[i] <- rho[i] * sum(s_hist[[i]] * q)
    q <- q - alpha_v[i] * y_hist[[i]]
  }
  gamma <- sum(s_hist[[m]] * y_hist[[m]]) / sum(y_hist[[m]]^2)
  r <- gamma * q
  for (i in 1:m) {
    beta_i <- rho[i] * sum(y_hist[[i]] * r)
    r <- r + (alpha_v[i] - beta_i) * s_hist[[i]]
  }
  -r
}

#' @keywords internal
.k3_masked_wide_spline_loop <- function(
    datasources, server_list, server_names, x_vars,
    coordinator, coordinator_conn, non_label_servers,
    transport_pks, std_data, y_var, family,
    betas, n_obs, lambda, log_n, log_scale,
    session_id, max_iter, tol, verbose,
    label_intercept, .dsAgg, .sendBlob) {

  # Fusion server = party 0 (s1), coordinator = label server (s3)
  # These two run the K=2 wide spline DCF
  fusion_server <- server_list[1]
  fusion_conn <- which(server_names == fusion_server)
  non_fusion_servers <- setdiff(server_list, fusion_server)

  # Wide spline parameters
  frac_bits <- 20L
  num_intervals <- if (family == "binomial") 50L else 100L

  # Helper for base64url conversion
  .to_b64url <- function(s) {
    s <- gsub("+", "-", s, fixed = TRUE)
    s <- gsub("/", "_", s, fixed = TRUE)
    gsub("=+$", "", s, perl = TRUE)
  }
  .b64url_to_b64 <- function(s) {
    s <- gsub("-", "+", s, fixed = TRUE)
    s <- gsub("_", "/", s, fixed = TRUE)
    pad <- nchar(s) %% 4
    if (pad == 2) s <- paste0(s, "==")
    if (pad == 3) s <- paste0(s, "=")
    s
  }

  # === Pre-generate DCF keys (one-time) ===
  if (verbose) message("  Pre-generating DCF keys for masked wide spline...")
  dcf_key_result <- dsVert:::.callMheTool("k2-dcf-gen-batch", list(
    n = as.integer(n_obs), num_intervals = as.integer(num_intervals)
  ))

  dcf_pair <- c(fusion_server, coordinator)
  for (srv in dcf_pair) {
    ci <- which(server_names == srv)
    is_coord <- (srv == coordinator)
    party_idx <- if (is_coord) 0L else 1L
    key_field <- paste0("party", party_idx, "_keys")
    pk_b64 <- .b64url_to_b64(transport_pks[[srv]])
    sealed <- dsVert:::.callMheTool("transport-encrypt", list(
      data = jsonlite::base64_enc(charToRaw(dcf_key_result[[key_field]])),
      recipient_pk = pk_b64
    ))
    .sendBlob(.to_b64url(sealed$sealed), "k2_dcf_keys_persistent", ci)
    .dsAgg(datasources[ci], call("k2StoreDcfKeysPersistentDS", session_id = session_id))
  }
  if (verbose) message("  DCF keys distributed to ", fusion_server, " + ", coordinator)

  # L-BFGS state
  ls <- list(); ly <- list(); pt <- NULL; pg <- NULL
  converged <- FALSE; final_iter <- 0

  for (iter in seq_len(max_iter)) {
    betas_old <- betas

    # =================================================================
    # Step 1: CKKS encrypt etas + aggregate + mask
    # =================================================================
    for (k in seq_along(server_list)) {
      server <- server_list[k]
      ci <- which(server_names == server)
      # Coordinator may have intercept in betas — separate it for encryption
      beta_for_encrypt <- betas[[server]]
      if (server == coordinator && label_intercept &&
          length(betas[[server]]) == length(x_vars[[server]]) + 1) {
        beta_for_encrypt <- betas[[server]][-1]  # drop intercept for X*beta
      }
      enc_result <- .dsAgg(datasources[ci],
        call("glmHEEncryptEtaDS", data_name = std_data,
             x_vars = x_vars[[server]], beta = beta_for_encrypt,
             session_id = session_id))
      if (is.list(enc_result)) enc_result <- enc_result[[1]]
      .sendBlob(enc_result$encrypted_eta, paste0("ct_eta_", k - 1), coordinator_conn)
    }
    # If coordinator has intercept, add it to the encrypted sum via a constant CT
    mws_intercept <- 0
    if (label_intercept && length(betas[[coordinator]]) > length(x_vars[[coordinator]])) {
      mws_intercept <- betas[[coordinator]][1]
    }

    # Coordinator: sum etas (skip polynomial — just addition)
    .dsAgg(datasources[coordinator_conn],
      call("glmHELinkStepDS", from_storage = TRUE,
           n_parties = as.integer(length(server_list)),
           skip_poly = TRUE, session_id = session_id))

    # Coordinator: mask with random r (add intercept if present)
    mask_result <- .dsAgg(datasources[coordinator_conn],
      call("glmMWSMaskEtaDS", n_obs = as.integer(n_obs),
           intercept = mws_intercept,
           frac_bits = as.integer(frac_bits), session_id = session_id))
    if (is.list(mask_result)) mask_result <- mask_result[[1]]
    ct_masked <- mask_result$ct_masked

    # Authorize ct_masked on all servers for threshold decryption
    ct_hash <- mask_result$ct_hash
    if (!is.null(ct_hash)) {
      for (server in server_list) {
        if (server != coordinator) {  # coordinator already registered it
          ci <- which(server_names == server)
          .sendBlob(ct_hash, "ct_hashes", ci)
          .dsAgg(datasources[ci],
            call("mheAuthorizeCTDS", op_type = "mws-masked-eta",
                 from_storage = TRUE, session_id = session_id))
        }
      }
    }

    # Distribute ct_masked to fusion server for threshold decryption
    .sendBlob(ct_masked, "ct_batch_1", fusion_conn)
    for (nf_server in non_fusion_servers) {
      nf_conn <- which(server_names == nf_server)
      nf_party_id <- which(server_list == nf_server) - 1
      .sendBlob(ct_masked, "ct_batch_1", nf_conn)
      pd <- .dsAgg(datasources[nf_conn],
        call("mhePartialDecryptBatchWrappedDS", n_cts = 1L, session_id = session_id))
      if (is.list(pd)) pd <- pd[[1]]
      .sendBlob(pd$wrapped_shares[1],
                paste0("wrapped_share_", nf_party_id, "_ct_1"), fusion_conn)
    }
    fuse_result <- .dsAgg(datasources[fusion_conn],
      call("mheFuseBatchDS", n_cts = 1L,
           n_parties = as.integer(length(server_list)),
           num_slots = as.integer(n_obs), session_id = session_id))
    if (is.list(fuse_result)) fuse_result <- fuse_result[[1]]
    eta_masked <- fuse_result$values  # η_total + r (plaintext on client!)

    # =================================================================
    # Step 2: Set Ring63 shares + run K=2 wide spline DCF
    # =================================================================
    # Fusion server: share_A = η_total + r (pass via blob, too large for call())
    eta_masked_json <- jsonlite::base64_enc(charToRaw(
      jsonlite::toJSON(eta_masked, digits = 17)))
    .sendBlob(eta_masked_json, "mws_eta_share", fusion_conn)
    .dsAgg(datasources[fusion_conn],
      call("glmMWSSetEtaShareDS", from_storage = TRUE,
           frac_bits = as.integer(frac_bits), session_id = session_id))

    # Coordinator already has share_B = -r (set by glmMWSMaskEtaDS)

    # Generate Beaver triples for this iteration
    triples <- lapply(1:3, function(i) dsVert:::.callMheTool("k2-gen-beaver-triples",
      list(n = as.integer(n_obs), frac_bits = as.integer(frac_bits))))

    for (srv in dcf_pair) {
      ci <- which(server_names == srv)
      is_coord <- (srv == coordinator)
      pk_b64 <- .b64url_to_b64(transport_pks[[srv]])
      party_idx <- if (is_coord) "party0" else "party1"
      td <- list()
      for (op in c("and", "had1", "had2")) {
        ti <- switch(op, and = 1, had1 = 2, had2 = 3)
        td[[paste0(op, "_a")]] <- triples[[ti]][[paste0(party_idx, "_u")]]
        td[[paste0(op, "_b")]] <- triples[[ti]][[paste0(party_idx, "_v")]]
        td[[paste0(op, "_c")]] <- triples[[ti]][[paste0(party_idx, "_w")]]
      }
      sealed_t <- dsVert:::.callMheTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(jsonlite::toJSON(td, auto_unbox = TRUE))),
        recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed_t$sealed), "k2_spline_triples", ci)
    }

    # Wide spline Phase 1: DCF masked values
    ph1 <- list()
    for (srv in dcf_pair) {
      ci <- which(server_names == srv); is_coord <- (srv == coordinator)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase1DS",
        party_id = if (is_coord) 0L else 1L, family = family,
        num_intervals = as.integer(num_intervals),
        frac_bits = as.integer(frac_bits), session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; ph1[[srv]] <- r
    }
    # Relay DCF masked values
    .sendBlob(ph1[[coordinator]]$dcf_masked, "k2_peer_dcf_masked",
              which(server_names == fusion_server))
    .sendBlob(ph1[[fusion_server]]$dcf_masked, "k2_peer_dcf_masked",
              coordinator_conn)

    # Wide spline Phase 2: DCF close + indicators + Beaver R1
    ph2 <- list()
    for (srv in dcf_pair) {
      ci <- which(server_names == srv); is_coord <- (srv == coordinator)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase2DS",
        party_id = if (is_coord) 0L else 1L, family = family,
        num_intervals = as.integer(num_intervals),
        frac_bits = as.integer(frac_bits), session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; ph2[[srv]] <- r
    }
    # Relay Beaver R1
    for (srv in dcf_pair) {
      peer <- setdiff(dcf_pair, srv)
      peer_ci <- which(server_names == peer)
      pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
      r1_json <- jsonlite::toJSON(list(
        and_xma = ph2[[srv]]$and_xma, and_ymb = ph2[[srv]]$and_ymb,
        had1_xma = ph2[[srv]]$had1_xma, had1_ymb = ph2[[srv]]$had1_ymb
      ), auto_unbox = TRUE)
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", peer_ci)
    }

    # Wide spline Phase 3: Beaver R2 for Had2
    ph3 <- list()
    for (srv in dcf_pair) {
      ci <- which(server_names == srv); is_coord <- (srv == coordinator)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase3DS",
        party_id = if (is_coord) 0L else 1L, family = family,
        num_intervals = as.integer(num_intervals),
        frac_bits = as.integer(frac_bits), session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]; ph3[[srv]] <- r
    }
    for (srv in dcf_pair) {
      peer <- setdiff(dcf_pair, srv)
      peer_ci <- which(server_names == peer)
      pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
      r1_json <- jsonlite::toJSON(list(
        had2_xma = ph3[[srv]]$had2_xma, had2_ymb = ph3[[srv]]$had2_ymb
      ), auto_unbox = TRUE)
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", peer_ci)
    }

    # Wide spline Phase 4: Final assembly → mu shares
    for (srv in dcf_pair) {
      ci <- which(server_names == srv); is_coord <- (srv == coordinator)
      .dsAgg(datasources[ci], call("k2WideSplinePhase4DS",
        party_id = if (is_coord) 0L else 1L, family = family,
        num_intervals = as.integer(num_intervals),
        frac_bits = as.integer(frac_bits), session_id = session_id))
    }

    # =================================================================
    # Step 3: Re-encrypt mu shares, sum, gradient
    # =================================================================
    # Get encrypted mu shares from both DCF parties
    mu_cts <- list()
    for (srv in dcf_pair) {
      ci <- which(server_names == srv)
      mu_result <- .dsAgg(datasources[ci],
        call("glmMWSGetMuShareDS", session_id = session_id))
      if (is.list(mu_result)) mu_result <- mu_result[[1]]
      mu_cts[[srv]] <- mu_result$encrypted_mu_share
    }

    # Sum mu shares: Enc(mu) = Enc(mu_A) + Enc(mu_B)
    ct_mu_sum <- dsVert:::.callMheTool("mhe-ct-add", list(
      ciphertext_a = .b64url_to_b64(mu_cts[[fusion_server]]),
      ciphertext_b = .b64url_to_b64(mu_cts[[coordinator]]),
      log_n = as.integer(log_n), log_scale = as.integer(log_scale)
    ))
    ct_mu <- base64_to_base64url(ct_mu_sum$ciphertext)

    # Distribute Enc(mu) to ALL servers for gradient computation
    for (server in server_list) {
      ci <- which(server_names == server)
      if (server == coordinator) {
        # Coordinator: store ct_mu locally (it generated the polynomial)
        .sendBlob(ct_mu, "ct_mu", ci)
      } else {
        # Non-coordinator: store via chunked transfer
        nc <- .dsvert_adaptive_send(ct_mu, function(chunk_str, chunk_idx, n_chunks) {
          .dsAgg(datasources[ci],
            call("mheStoreEncChunkDS", col_index = 3L,
                 chunk_index = chunk_idx, chunk = chunk_str,
                 session_id = session_id))
        })
        .dsAgg(datasources[ci],
          call("mheAssembleEncColumnDS", col_index = 3L,
               n_chunks = as.integer(nc), session_id = session_id))
      }
    }

    # Each server: gradient = X_k^T * (Enc(y) - Enc(mu))
    server_gradients <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      grad_result <- .dsAgg(datasources[ci],
        call("glmHEGradientEncDS", data_name = std_data,
             x_vars = x_vars[[server]], num_obs = as.integer(n_obs),
             from_storage = (server != coordinator),
             session_id = session_id))
      if (is.list(grad_result)) grad_result <- grad_result[[1]]

      enc_grads <- grad_result$encrypted_gradients
      p_k <- length(enc_grads)

      # Threshold decrypt each gradient component
      for (nf_server in non_fusion_servers) {
        nf_conn <- which(server_names == nf_server)
        nf_party_id <- which(server_list == nf_server) - 1
        for (j in seq_len(p_k))
          .sendBlob(enc_grads[j], paste0("ct_batch_", j), nf_conn)
        pd <- .dsAgg(datasources[nf_conn],
          call("mhePartialDecryptBatchWrappedDS",
               n_cts = as.integer(p_k), session_id = session_id))
        if (is.list(pd)) pd <- pd[[1]]
        for (j in seq_len(p_k))
          .sendBlob(pd$wrapped_shares[j],
                    paste0("wrapped_share_", nf_party_id, "_ct_", j),
                    fusion_conn)
      }
      for (j in seq_len(p_k))
        .sendBlob(enc_grads[j], paste0("ct_batch_", j), fusion_conn)
      fuse_grad <- .dsAgg(datasources[fusion_conn],
        call("mheFuseBatchDS", n_cts = as.integer(p_k),
             n_parties = as.integer(length(server_list)),
             num_slots = 0L, session_id = session_id))
      if (is.list(fuse_grad)) fuse_grad <- fuse_grad[[1]]
      server_gradients[[server]] <- fuse_grad$values
    }

    # =================================================================
    # Step 4: L-BFGS update (client-side)
    # =================================================================
    raw_gradient <- numeric(0)
    for (server in server_list)
      raw_gradient <- c(raw_gradient, server_gradients[[server]])

    theta <- unlist(betas[server_list])
    full_grad <- -raw_gradient / n_obs + lambda * theta

    if (!is.null(pt)) {
      sk <- theta - pt; yk <- full_grad - pg
      if (sum(sk * yk) > 1e-10) {
        ls <- c(ls, list(sk)); ly <- c(ly, list(yk))
        if (length(ls) > 7) { ls <- ls[-1]; ly <- ly[-1] }
      }
    }
    pt <- theta; pg <- full_grad

    direction <- .lbfgs_direction_mws(full_grad, ls, ly)
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

    if (verbose)
      message(sprintf("  Iter %d: max diff = %.2e", iter, max_diff))

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose)
        message(sprintf("  Converged after %d iterations (diff = %.2e)", iter, max_diff))
      break
    }

    if (iter %% 10 == 0) {
      for (server in server_list)
        tryCatch(.dsAgg(datasources[which(server_names == server)],
          call("mheGcDS")), error = function(e) NULL)
    }
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)", max_iter, max_diff))

  list(betas = betas, converged = converged, final_iter = final_iter)
}
