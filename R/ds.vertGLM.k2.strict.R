#' @title K=2 Wide Spline: DCF-based Training Pipeline
#' @description Secure K=2 GLM training with wide piecewise-linear spline
#'   for sigmoid (binomial) and exp (Poisson) link functions.
#'   Fully distributed: no eta reconstruction, non-disclosive.
#'   Uses DCF comparisons + Beaver Hadamard for spline evaluation.
#' @name k2-wide-spline-client
NULL

#' @keywords internal
.k2_strict_loop <- function(datasources, server_names, server_list,
                             coordinator, coordinator_conn,
                             non_label_servers, nl, nl_conn,
                             x_vars, y_var, std_data,
                             transport_pks, session_id,
                             family, lambda, max_iter, tol,
                             n_obs, verbose, .dsAgg, .sendBlob) {

  coordinator_pk <- transport_pks[[coordinator]]
  nl_pk <- transport_pks[[nl]]
  frac_bits <- 20L

  # Wide spline configuration
  num_intervals <- if (family == "poisson") 100L else 50L
  if (verbose) message(sprintf("  Wide spline: %s, %d intervals", family, num_intervals))

  # Learning rate
  alpha <- if (family == "poisson") 0.05 else 0.3
  if (is.null(lambda)) lambda <- 1e-4

  p_coord <- length(x_vars[[coordinator]])
  p_nl <- length(x_vars[[nl]])
  p_total <- p_coord + p_nl

  # === INPUT-SHARING PREAMBLE (FixedPoint) ===
  if (verbose) message("  Input-sharing preamble (FixedPoint)...")

  share_results <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    is_label <- (server == coordinator)
    r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
      data_name = std_data, x_vars = x_vars[[server]],
      y_var = if (is_label) y_var else NULL,
      peer_pk = transport_pks[[peer]], session_id = session_id))
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
    .dsAgg(datasources[ci], call("k2ReceiveShareDS",
      peer_p = as.integer(length(x_vars[[peer]])), session_id = session_id))
  }
  if (verbose) message("  Input sharing complete")

  # Base64url helpers
  .to_b64url <- function(x) gsub("+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-", "+", gsub("_", "/", x, fixed=TRUE), fixed=TRUE)
    pad <- nchar(x) %% 4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }

  beta <- rep(0, p_total)
  intercept <- 0.0
  converged <- FALSE
  final_iter <- 0

  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    intercept_old <- intercept

    # === Step 1: Compute eta shares in FP (LOCAL) ===
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2ComputeEtaShareDS",
        beta_coord = beta[1:p_coord],
        beta_nl = beta[(p_coord+1):p_total],
        intercept = if (is_coord) intercept else 0.0,
        is_coordinator = is_coord,
        session_id = session_id))
    }

    # === Step 2: Wide spline DCF sigmoid/exp ===

    # 2a. Generate DCF keys on client
    dcf <- dsVert:::.callMheTool("k2-dcf-gen-batch", list(
      family = family, n = as.integer(n_obs),
      frac_bits = frac_bits, num_intervals = num_intervals))

    # 2b. Send keys to servers via encrypted blob
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      keys_b64 <- if (is_coord) dcf$party0_keys else dcf$party1_keys
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data = keys_b64, recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_dcf_keys", ci)
    }

    # 2c. Servers store keys + evaluate DCF phase 1
    dcf_r1 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2StoreDcfKeysDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call("k2DcfEvalDS",
        phase = 1L, party_id = if (is_coord) 0L else 1L,
        family = family, num_intervals = num_intervals,
        frac_bits = frac_bits, session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      dcf_r1[[server]] <- r
    }

    # 2d. Relay masked values
    .sendBlob(dcf_r1[[coordinator]]$masked_values, "k2_dcf_peer_masked", nl_conn)
    .sendBlob(dcf_r1[[nl]]$masked_values, "k2_dcf_peer_masked", coordinator_conn)

    # 2e. DCF phase 2 → comparison shares
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2DcfEvalDS",
        phase = 2L, party_id = if (is_coord) 0L else 1L,
        family = family, num_intervals = num_intervals,
        frac_bits = frac_bits, session_id = session_id))
    }

    # 2f. Compute spline indicators locally (no communication)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2SplineIndicatorsDS",
        party_id = if (is_coord) 0L else 1L,
        family = family, num_intervals = num_intervals,
        frac_bits = frac_bits, session_id = session_id))
    }

    # 2g. Beaver AND: I_mid = NOT(c_low) * c_high
    # Reuse k2BeaverRoundFPDS with x_key="k2_c_low_share_fp", y_key="k2_c_high_share_fp"
    triple_and <- dsVert:::.callMheTool("k2-gen-beaver-triples", list(
      n = as.integer(n_obs), frac_bits = frac_bits))
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      party_idx <- if (is_coord) "party0" else "party1"
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      triple_json <- jsonlite::toJSON(list(
        a = triple_and[[paste0(party_idx, "_u")]],
        b = triple_and[[paste0(party_idx, "_v")]],
        c = triple_and[[paste0(party_idx, "_w")]]), auto_unbox = TRUE)
      msg_b64 <- jsonlite::base64_enc(charToRaw(triple_json))
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(data = msg_b64, recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_pow_triple_1", ci)
    }
    # Phase 1
    and_r1 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2StorePowerTripleDS", triple_key = "k2_pow_triple_1", session_id = session_id))
      r <- .dsAgg(datasources[ci], call("k2BeaverRoundFPDS",
        x_key = "k2_c_low_share_fp", y_key = "k2_c_high_share_fp",
        a_fp = "NONE", b_fp = "NONE",
        party_id = if (is_coord) 0L else 1L, phase = 1L,
        use_session_triple = 1L, session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      and_r1[[server]] <- r
    }
    # Relay
    for (server in server_list) {
      peer <- setdiff(server_list, server)
      peer_ci <- which(server_names == peer)
      msg_json <- jsonlite::toJSON(list(xma = and_r1[[server]]$xma_fp, ymb = and_r1[[server]]$ymb_fp), auto_unbox = TRUE)
      msg_b64 <- jsonlite::base64_enc(charToRaw(msg_json))
      pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(data = msg_b64, recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_br_peer_1", peer_ci)
    }
    # Phase 2 → result stored as k2_i_mid_fp
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2BeaverRoundFPDS",
        x_key = "k2_c_low_share_fp", y_key = "k2_c_high_share_fp",
        a_fp = "NONE", b_fp = "NONE", c_fp = "NONE",
        peer_xma_fp = "NONE", peer_ymb_fp = "NONE",
        result_key = "k2_i_mid_fp", party_id = if (is_coord) 0L else 1L,
        phase = 2L, use_session_triple = 1L,
        peer_blob_key = "k2_br_peer_1", session_id = session_id))
    }

    # 2h. Hadamard: slope * x (reuse Beaver for Hadamard)
    triple_sx <- dsVert:::.callMheTool("k2-gen-beaver-triples", list(
      n = as.integer(n_obs), frac_bits = frac_bits))
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      party_idx <- if (is_coord) "party0" else "party1"
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      triple_json <- jsonlite::toJSON(list(
        a = triple_sx[[paste0(party_idx, "_u")]],
        b = triple_sx[[paste0(party_idx, "_v")]],
        c = triple_sx[[paste0(party_idx, "_w")]]), auto_unbox = TRUE)
      msg_b64 <- jsonlite::base64_enc(charToRaw(triple_json))
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(data = msg_b64, recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_pow_triple_2", ci)
    }
    sx_r1 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2StorePowerTripleDS", triple_key = "k2_pow_triple_2", session_id = session_id))
      r <- .dsAgg(datasources[ci], call("k2BeaverRoundFPDS",
        x_key = "k2_slope_share_fp", y_key = "secure_eta_share",
        a_fp = "NONE", b_fp = "NONE",
        party_id = if (is_coord) 0L else 1L, phase = 1L,
        use_session_triple = 1L, session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      sx_r1[[server]] <- r
    }
    for (server in server_list) {
      peer <- setdiff(server_list, server)
      peer_ci <- which(server_names == peer)
      msg_json <- jsonlite::toJSON(list(xma = sx_r1[[server]]$xma_fp, ymb = sx_r1[[server]]$ymb_fp), auto_unbox = TRUE)
      msg_b64 <- jsonlite::base64_enc(charToRaw(msg_json))
      pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(data = msg_b64, recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_br_peer_2", peer_ci)
    }
    # Phase 2 → k2_spline_value_fp = slope*x + intercept
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2BeaverRoundFPDS",
        x_key = "k2_slope_share_fp", y_key = "secure_eta_share",
        a_fp = "NONE", b_fp = "NONE", c_fp = "NONE",
        peer_xma_fp = "NONE", peer_ymb_fp = "NONE",
        result_key = "k2_slope_x_fp", party_id = if (is_coord) 0L else 1L,
        phase = 2L, use_session_triple = 1L,
        peer_blob_key = "k2_br_peer_2", session_id = session_id))
    }

    # TODO: Add intercept to spline value on server: spline = slope*x + intercept
    # For now, the k2SplineAssembleDS can handle this if we store both.
    # The intercept share was stored by k2SplineIndicatorsDS as k2_intercept_share_fp.
    # The slope*x was stored by BeaverRound as k2_slope_x_fp.
    # spline_value = k2_slope_x_fp + k2_intercept_share_fp → need a server-side add.

    # 2i. Hadamard: I_mid * spline_value (need I_mid scaled to FP first)
    # I_mid from step 2g is an integer share. Scale to FP on server, then Hadamard.
    # For now, use k2SplineAssembleDS which handles the remaining assembly.

    # 2j. Final assembly: mu = I_high + I_mid * spline
    # This needs another Beaver Hadamard for I_mid * spline.
    # But first we need spline = slope*x + intercept (server-side add).
    # TODO: Add a server-side command for this or extend k2SplineAssembleDS.

    # SIMPLIFIED for first integration: use k2SplineAssembleDS which computes
    # mu from the stored session values. (Implementation pending.)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2SplineAssembleDS",
        party_id = if (is_coord) 0L else 1L,
        family = family, frac_bits = frac_bits,
        session_id = session_id))
    }

    # === Step 3: Ring63 Beaver gradient ===
    mvt <- dsVert:::.callMheTool("k2-gen-matvec-triples", list(
      n = as.integer(n_obs), p = as.integer(p_total)))

    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      party_idx <- if (is_coord) "party0" else "party1"
      a_fp <- mvt[[paste0(party_idx, "_a")]]
      b_fp <- mvt[[paste0(party_idx, "_b")]]
      c_fp <- mvt[[paste0(party_idx, "_c")]]
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      msg_json <- jsonlite::toJSON(list(a=a_fp, b=b_fp, c=c_fp), auto_unbox=TRUE)
      msg_b64 <- jsonlite::base64_enc(charToRaw(msg_json))
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(data=msg_b64, recipient_pk=pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_grad_triple_fp", ci)
    }

    r1_results <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      .dsAgg(datasources[ci], call("k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call("k2GradientR1DS",
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
      r <- .dsAgg(datasources[ci], call("k2GradientR2DS",
        party_id = if (is_coord) 0L else 1L,
        session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      grad_results[[server]] <- r
    }

    # === Step 4: Client sums gradient shares in Ring63 ===
    grad_fp_coord <- grad_results[[coordinator]]$gradient_fp
    grad_fp_nl <- grad_results[[nl]]$gradient_fp
    res_fp_coord <- r1_results[[coordinator]]$sum_residual_fp
    res_fp_nl <- r1_results[[nl]]$sum_residual_fp

    if (!is.null(grad_fp_coord) && !is.null(grad_fp_nl)) {
      agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a = grad_fp_coord, share_b = grad_fp_nl, frac_bits = frac_bits))
      gradient <- agg$values
      agg_res <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a = res_fp_coord, share_b = res_fp_nl, frac_bits = frac_bits))
      sum_residual <- agg_res$values[1]
    } else {
      gradient <- rep(0, p_total)
      sum_residual <- 0
      for (server in server_list) {
        gradient <- gradient + grad_results[[server]]$gradient_share
        sum_residual <- sum_residual + r1_results[[server]]$sum_residual
      }
    }

    # GD update
    full_grad <- gradient / n_obs + lambda * beta
    grad_norm <- sqrt(sum(full_grad^2))
    if (grad_norm > 5.0) full_grad <- full_grad * (5.0 / grad_norm)
    beta <- beta - alpha * full_grad
    intercept <- intercept - alpha * sum_residual / n_obs

    max_diff <- max(abs(beta - beta_old), abs(intercept - intercept_old))
    final_iter <- iter

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose) message(sprintf("  Converged after %d iterations (diff = %.2e)", iter, max_diff))
      break
    }
    if (verbose && iter %% 50 == 0) message(sprintf("  Iteration %d: max diff = %.2e, ||grad||=%.2e", iter, max_diff, grad_norm))

    # Periodic GC
    if (iter %% 20 == 0) {
      for (server in server_list) {
        tryCatch(.dsAgg(datasources[which(server_names==server)], call("mheGcDS")), error=function(e) NULL)
      }
    }
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)", max_iter, max_diff))

  betas <- list()
  betas[[coordinator]] <- beta[1:p_coord]
  betas[[nl]] <- beta[(p_coord+1):p_total]

  list(betas=betas, intercept=intercept, converged=converged, iterations=final_iter, max_diff=max_diff)
}
