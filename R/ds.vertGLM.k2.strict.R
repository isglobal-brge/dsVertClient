#' @title K=2 Strict: All-Ring63 Training Pipeline
#' @description Secure K=2 training with ALL computation in FixedPoint ring.
#'   Input sharing → FP eta → Beaver poly eval → Ring63 Beaver gradient.
#'   Only final p_k gradient scalars converted to float64.
#' @name k2-strict-client
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

  poly_info <- tryCatch(
    dsVert:::.callMheTool("k2-chebyshev-coeffs", list(family = family, degree = 7L)),
    error = function(e) dsVert:::.callMheTool("mpc-get-poly-coeffs", list(family = family, degree = 7L))
  )
  poly_coeffs <- poly_info$coefficients
  poly_degree <- length(poly_coeffs) - 1

  if (verbose) message(sprintf("  Chebyshev degree-%d, max poly error: %.2e", poly_degree, poly_info$max_error))

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

  # Base64 helpers
  .b64url_to_b64 <- function(x) {
    x <- gsub("-", "+", gsub("_", "/", x, fixed=TRUE), fixed=TRUE)
    pad <- nchar(x) %% 4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }

  beta <- rep(0, p_total)
  intercept <- 0.0
  alpha <- if (family == "poisson") 0.1 else 0.3
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

    # === Step 2: Beaver polynomial eval for mu shares (FP) ===
    n_poly_triples <- poly_degree * n_obs
    triple_result <- dsVert:::.callMheTool("k2-gen-beaver-triples", list(
      n = as.integer(n_poly_triples), frac_bits = as.integer(frac_bits)))

    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      party_idx <- if (is_coord) "party0" else "party1"
      u_share <- triple_result[[paste0(party_idx, "_u")]]
      v_share <- triple_result[[paste0(party_idx, "_v")]]
      w_share <- triple_result[[paste0(party_idx, "_w")]]
      msg_json <- jsonlite::toJSON(list(u=u_share, v=v_share, w=w_share), auto_unbox=TRUE)
      msg_b64 <- jsonlite::base64_enc(charToRaw(msg_json))
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(data=msg_b64, recipient_pk=pk_b64))
      .sendBlob(gsub("+", "-", gsub("/", "_", gsub("=+$", "", sealed$sealed, perl=TRUE), fixed=TRUE), fixed=TRUE), "beaver_triples", ci)
      .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="store_triples",
             session_id=session_id))
    }

    # Power chain
    tri_offset <- 0
    for (k in seq_len(poly_degree)) {
      a_key <- if (k==1) "secure_eta_share" else paste0("secure_pow", k)
      b_key <- "secure_eta_share"
      result_key <- paste0("secure_pow", k+1)
      open_results <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        peer <- setdiff(server_list, server)
        r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="beaver_open",
          a_share_key=a_key, b_share_key=b_key,
          u_vals=tri_offset, v_vals=n_obs,
          peer_pk=transport_pks[[peer]], frac_bits=frac_bits, session_id=session_id))
        if (is.list(r)) r <- r[[1]]
        open_results[[server]] <- r
      }
      .sendBlob(open_results[[coordinator]]$peer_de_enc, "beaver_peer_de", nl_conn)
      .sendBlob(open_results[[nl]]$peer_de_enc, "beaver_peer_de", coordinator_conn)
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)
        .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="beaver_close",
          result_key=result_key, u_vals=tri_offset, v_vals=n_obs,
          party_id=if(is_coord) 0L else 1L, frac_bits=frac_bits, session_id=session_id))
      }
      tri_offset <- tri_offset + n_obs
    }

    # Polynomial eval → mu shares (FP)
    power_keys <- c("secure_eta_share", paste0("secure_pow", 2:(poly_degree+1)))
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step="poly_eval",
        power_keys=power_keys, coefficients=poly_coeffs,
        party_id=if(is_coord) 0L else 1L, frac_bits=frac_bits, session_id=session_id))
    }

    # === Step 3: Ring63 Beaver matvec gradient ===
    # Generate matvec triples (ring multiply for C)
    mvt <- dsVert:::.callMheTool("k2-gen-matvec-triples", list(
      n = as.integer(n_obs), p = as.integer(p_total)))

    # Distribute triple shares to both parties
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      party_idx <- if (is_coord) "party0" else "party1"

      # Store A, B, C shares in session
      a_fp <- mvt[[paste0(party_idx, "_a")]]
      b_fp <- mvt[[paste0(party_idx, "_b")]]
      c_fp <- mvt[[paste0(party_idx, "_c")]]

      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      msg_json <- jsonlite::toJSON(list(a=a_fp, b=b_fp, c=c_fp), auto_unbox=TRUE)
      msg_b64 <- jsonlite::base64_enc(charToRaw(msg_json))
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(data=msg_b64, recipient_pk=pk_b64))
      .sendBlob(gsub("+", "-", gsub("/", "_", gsub("=+$", "", sealed$sealed, perl=TRUE), fixed=TRUE), fixed=TRUE), "k2_grad_triple_fp", ci)
    }

    # Each party stores the triple and computes gradient R1
    r1_results <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      peer <- setdiff(server_list, server)

      # First: store the triple
      .dsAgg(datasources[ci], call("k2StoreGradTripleDS", session_id = session_id))

      # Then: compute gradient round 1
      r <- .dsAgg(datasources[ci], call("k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      r1_results[[server]] <- r
    }

    # Relay R1 messages
    .sendBlob(r1_results[[coordinator]]$encrypted_r1, "k2_grad_peer_r1", nl_conn)
    .sendBlob(r1_results[[nl]]$encrypted_r1, "k2_grad_peer_r1", coordinator_conn)

    # Gradient round 2
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

    # === Step 4: Client sums gradient shares (float64) and updates beta ===
    gradient <- rep(0, p_total)
    sum_residual <- 0
    for (server in server_list) {
      gradient <- gradient + grad_results[[server]]$gradient_share
      sum_residual <- sum_residual + r1_results[[server]]$sum_residual
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
    if (verbose && iter %% 5 == 0) message(sprintf("  Iteration %d: max diff = %.2e", iter, max_diff))

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
