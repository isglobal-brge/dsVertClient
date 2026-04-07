#' @title K>=3 Masked Wide Spline Loop
#' @description Non-disclosive K>=3 binomial/Poisson via masked wide spline.
#'
#' Protocol per iteration:
#' \enumerate{
#'   \item CKKS: each server encrypts eta_k, coordinator sums + masks with r
#'   \item Threshold decrypt masked sum: s1 gets (eta+r), s3 keeps -r
#'   \item Ring63/DCF: K=2 wide spline between s1 and s3 on shares
#'   \item CKKS: re-encrypt mu shares, sum → Enc(mu), gradient X_k^T * Enc(y-mu)
#'   \item Threshold decrypt gradient → L-BFGS update
#' }
#'
#' Nobody sees eta_total or mu in plaintext. Precision: ~1e-2 (wide spline ~1e-4 error).
#'
#' @name k3-masked-wide-spline
NULL

# L-BFGS two-loop recursion (Nocedal & Wright Algorithm 7.4)
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

#' Run K>=3 Masked Wide Spline loop
#'
#' @param datasources DataSHIELD connections.
#' @param server_list Character vector of server names.
#' @param server_names Names of datasources.
#' @param x_vars Named list of variable names per server.
#' @param coordinator Name of label/coordinator server.
#' @param coordinator_conn Integer. Connection index for coordinator.
#' @param non_label_servers Character vector of non-label server names.
#' @param transport_pks Named list of transport public keys.
#' @param std_data Character. Name of standardized data frame.
#' @param y_var Character. Response variable name.
#' @param family Character. "binomial" or "poisson".
#' @param betas Named list of initial beta vectors per server.
#' @param n_obs Integer. Number of observations.
#' @param lambda Numeric. L2 regularization parameter.
#' @param log_n Integer. CKKS ring dimension (12 for MWS — only addition needed).
#' @param log_scale Integer. CKKS scale parameter.
#' @param session_id Character. Session ID.
#' @param max_iter Integer. Maximum iterations.
#' @param tol Numeric. Convergence tolerance.
#' @param verbose Logical. Print progress.
#' @param label_intercept Logical. Include intercept on coordinator.
#' @param .dsAgg Function. DataSHIELD aggregate wrapper.
#' @param .sendBlob Function. Blob send wrapper.
#'
#' @return List with betas, converged, final_iter.
#' @keywords internal
.k3_masked_wide_spline_loop <- function(
    datasources, server_list, server_names, x_vars,
    coordinator, coordinator_conn, non_label_servers,
    transport_pks, std_data, y_var, family,
    betas, n_obs, lambda, log_n, log_scale,
    session_id, max_iter, tol, verbose,
    label_intercept, .dsAgg, .sendBlob) {

  # The fusion server (party 0) and coordinator run the wide spline
  fusion_server <- server_list[1]
  fusion_conn <- which(server_names == fusion_server)
  non_fusion_servers <- setdiff(server_list, fusion_server)

  # Wide spline parameters (same as K=2 production)
  frac_bits <- 20L
  num_intervals <- if (family == "binomial") 50L else 100L

  # === Pre-generate DCF keys (one-time, reusable across iterations) ===
  if (verbose) message("  Pre-generating DCF keys for masked wide spline...")

  # DCF keys for fusion server and coordinator (the 2-party pair)
  dcf_pair <- list(fusion_server, coordinator)
  dcf_key_result <- dsVert:::.callMheTool("k2-dcf-gen-batch", list(
    n = as.integer(n_obs), num_intervals = as.integer(num_intervals)
  ))
  for (srv in dcf_pair) {
    ci <- which(server_names == srv)
    is_coord <- (srv == coordinator)
    party_idx <- if (is_coord) 0L else 1L
    key_field <- paste0("party", party_idx, "_keys")
    pk_b64 <- gsub("-","+", gsub("_","/", transport_pks[[srv]], fixed=TRUE), fixed=TRUE)
    pad <- nchar(pk_b64) %% 4
    if (pad == 2) pk_b64 <- paste0(pk_b64, "==")
    if (pad == 3) pk_b64 <- paste0(pk_b64, "=")
    sealed <- dsVert:::.callMheTool("transport-encrypt", list(
      data = jsonlite::base64_enc(charToRaw(dcf_key_result[[key_field]])),
      recipient_pk = pk_b64
    ))
    sealed_b64url <- gsub("+","-", gsub("/","_", gsub("=+$","", sealed$sealed, perl=TRUE), fixed=TRUE), fixed=TRUE)
    .sendBlob(sealed_b64url, "k2_dcf_keys", ci)
    .dsAgg(datasources[ci], call("k2StoreDcfKeysPersistentDS", session_id = session_id))
  }
  if (verbose) message("  DCF keys distributed")

  # L-BFGS state
  ls <- list(); ly <- list(); pt <- NULL; pg <- NULL
  converged <- FALSE; final_iter <- 0

  for (iter in seq_len(max_iter)) {
    betas_old <- betas

    # =================================================================
    # Step 1: CKKS encrypt etas + aggregate + mask
    # =================================================================
    # Each server encrypts η_k under CPK
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci], call("glmHEEncryptEtaDS",
        data_name = std_data, x_vars = x_vars[[server]],
        beta = betas[[server]], clip_radius = NULL,
        session_id = session_id))
      party_idx <- which(server_list == server) - 1
      # Send encrypted eta to coordinator
      enc_result <- .dsAgg(datasources[ci], call("glmHEEncryptEtaDS",
        data_name = std_data, x_vars = x_vars[[server]],
        beta = betas[[server]], session_id = session_id))
      if (is.list(enc_result)) enc_result <- enc_result[[1]]
      .sendBlob(enc_result$encrypted_eta, paste0("ct_eta_", party_idx), coordinator_conn)
    }

    # Coordinator: sum etas + evaluate identity (no polynomial, just addition)
    link_result <- .dsAgg(datasources[coordinator_conn],
      call("glmHELinkStepDS", from_storage = TRUE,
           n_parties = as.integer(length(server_list)),
           skip_poly = TRUE,  # Just sum, no polynomial!
           session_id = session_id))
    if (is.list(link_result)) link_result <- link_result[[1]]
    ct_eta_total <- link_result$ct_mu  # Actually ct_eta_total (skip_poly=TRUE)

    # Coordinator encrypts random mask r and adds to ct_eta_total
    # Then threshold-decrypt (η_total + r) → fusion server gets share_A
    # Coordinator keeps share_B = -r
    # TODO: implement this step (requires new server function)

    # =================================================================
    # Step 2-3: Convert shares to Ring63 + wide spline DCF
    # =================================================================
    # TODO: convert plaintext shares to Ring63 FP
    # TODO: run wide spline phases 1-4 between fusion server and coordinator
    # TODO: get μ shares

    # =================================================================
    # Step 4: Re-encrypt μ shares, gradient, threshold decrypt
    # =================================================================
    # TODO: encrypt μ shares under CPK, sum → Enc(μ)
    # TODO: Enc(y - μ) = Enc(y) - Enc(μ)
    # TODO: each server: gradient via InnerSum
    # TODO: threshold decrypt gradients

    # =================================================================
    # Step 5: L-BFGS update (client-side)
    # =================================================================
    # TODO: assemble full gradient, L-BFGS direction, update betas

    # PLACEHOLDER: use exact computation for now (will be replaced by MPC)
    break
  }

  list(betas = betas, converged = converged, final_iter = final_iter)
}
