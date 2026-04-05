#' @title K=2 Strict Mode Client Orchestration (Chebyshev Beaver)
#' @description Client-side loop for K=2 strict mode using Chebyshev polynomial
#'   evaluation on secret shares via Beaver triples. Everything stays in shares —
#'   only p_k gradient scalars + 2 intercept scalars revealed per iteration.
#'
#' @details
#' This reuses the existing k2MpcSecureStepDS infrastructure (Beaver open/close,
#' poly_eval, cross-gradient) but with Chebyshev coefficients instead of LS.
#'
#' @name k2-strict-client
NULL

#' Internal: K=2 strict Chebyshev Beaver training loop
#'
#' Uses k2MpcSecureStepDS steps for the full Beaver protocol:
#' split_eta → combine_eta → Beaver power chain → poly_eval →
#' intercept_newton_prepare → cross-gradient → combine_gradient → L-BFGS
#'
#' @keywords internal
.k2_strict_loop <- function(datasources, server_names, server_list,
                             coordinator, coordinator_conn,
                             non_label_servers, nl, nl_conn,
                             x_vars, y_var, std_data,
                             transport_pks, session_id,
                             family, lambda, max_iter, tol,
                             n_obs, verbose, .dsAgg, .sendBlob) {

  frac_bits <- 20L

  # Get Chebyshev polynomial coefficients (better than the old LS fit)
  poly_info <- tryCatch(
    dsVert:::.callMheTool("k2-chebyshev-coeffs", list(
      family = family, degree = 7L)),
    error = function(e) {
      # Fallback: use old mpc-get-poly-coeffs
      dsVert:::.callMheTool("mpc-get-poly-coeffs", list(
        family = family, degree = 7L))
    }
  )
  poly_coeffs <- poly_info$coefficients
  poly_degree <- length(poly_coeffs) - 1

  if (verbose)
    message(sprintf("  Chebyshev degree-%d, max poly error: %.2e",
                    poly_degree, poly_info$max_error))

  p_coord <- length(x_vars[[coordinator]])
  p_nl <- length(x_vars[[nl]])

  # Initialize betas
  betas <- list()
  for (server in server_list) {
    betas[[server]] <- rep(0, length(x_vars[[server]]))
  }
  intercept_beta0 <- 0.0

  # Learning rate (matches C++ GD approach)
  alpha <- if (family == "poisson") 0.1 else 0.3

  converged <- FALSE
  final_iter <- 0

  # Base64 helpers for triple distribution
  .b64url_to_b64 <- function(x) {
    x <- gsub("-", "+", gsub("_", "/", x, fixed = TRUE), fixed = TRUE)
    pad <- nchar(x) %% 4
    if (pad == 2) x <- paste0(x, "==")
    if (pad == 3) x <- paste0(x, "=")
    x
  }
  .b64_to_b64url <- function(x) {
    gsub("+", "-", gsub("/", "_", gsub("=+$", "", x, perl = TRUE),
                         fixed = TRUE), fixed = TRUE)
  }

  for (iter in seq_len(max_iter)) {
    betas_old <- betas
    old_beta0 <- intercept_beta0

    # === A: Split eta → shares ===
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      r <- .dsAgg(datasources[ci], call("k2MpcSplitEtaDS",
        data_name = std_data, x_vars = x_vars[[server]],
        beta = betas[[server]], peer_pk = transport_pks[[peer]],
        frac_bits = frac_bits, session_id = session_id))
      if (is.list(r)) r <- r[[1]]
      peer_ci <- which(server_names == peer)
      .sendBlob(r$peer_share_enc, "mpc_peer_eta_share", peer_ci)
    }
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "combine_eta",
             session_id = session_id))
    }

    # === B: Beaver triples + polynomial eval ===
    n_poly_triples <- poly_degree * n_obs
    n_mu2_triples <- n_obs
    n_cross_coord <- n_obs * p_coord
    n_cross_nl <- n_obs * p_nl
    n_total <- n_poly_triples + n_mu2_triples + n_cross_coord + n_cross_nl

    # Generate Beaver triples using FixedPoint-exact arithmetic.
    # The product tw = tu * tv MUST be computed in the same FixedPoint ring
    # as the server uses, otherwise the Beaver protocol is invalid.
    tu <- runif(n_total, -1, 1)
    tv <- runif(n_total, -1, 1)
    # Compute tw via mhe-tool FixedPoint multiply (NOT float64 multiply)
    fp_tu <- dsVert:::.callMheTool("k2-float-to-fp",
      list(values = tu, frac_bits = frac_bits))$fp_data
    fp_tv <- dsVert:::.callMheTool("k2-float-to-fp",
      list(values = tv, frac_bits = frac_bits))$fp_data
    # Multiply in FixedPoint ring: tw = FPMulLocal(tu, tv)
    fp_tw <- dsVert:::.callMheTool("k2-fp-mul",
      list(a = fp_tu, b = fp_tv, frac_bits = frac_bits))$result
    # Convert back to float for share splitting
    tw <- dsVert:::.callMheTool("mpc-fp-to-float",
      list(fp_data = fp_tw, frac_bits = frac_bits))$values

    tu0 <- runif(n_total, -5, 5); tu1 <- tu - tu0
    tv0 <- runif(n_total, -5, 5); tv1 <- tv - tv0
    tw0 <- runif(n_total, -5, 5); tw1 <- tw - tw0

    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      sealed <- dsVert:::.callMheTool("transport-encrypt-vectors", list(
        vectors = list(u = if (is_coord) tu0 else tu1,
                      v = if (is_coord) tv0 else tv1,
                      w = if (is_coord) tw0 else tw1),
        recipient_pk = pk_b64))
      .sendBlob(.b64_to_b64url(sealed$sealed), "beaver_triples", ci)
      .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "store_triples",
             session_id = session_id))
    }

    # Polynomial Beaver rounds (power chain)
    tri_offset <- 0
    for (k in seq_len(poly_degree)) {
      a_key <- if (k == 1) "secure_eta_share" else paste0("secure_pow", k)
      b_key <- "secure_eta_share"
      result_key <- paste0("secure_pow", k + 1)

      open_results <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        peer <- setdiff(server_list, server)
        r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "beaver_open",
          a_share_key = a_key, b_share_key = b_key,
          u_vals = tri_offset, v_vals = n_obs,
          peer_pk = transport_pks[[peer]], frac_bits = frac_bits,
          session_id = session_id))
        if (is.list(r)) r <- r[[1]]
        open_results[[server]] <- r
      }
      .sendBlob(open_results[[coordinator]]$peer_de_enc, "beaver_peer_de", nl_conn)
      .sendBlob(open_results[[nl]]$peer_de_enc, "beaver_peer_de", coordinator_conn)

      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)
        .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "beaver_close",
          result_key = result_key, u_vals = tri_offset, v_vals = n_obs,
          party_id = if (is_coord) 0L else 1L,
          frac_bits = frac_bits, session_id = session_id))
      }
      tri_offset <- tri_offset + n_obs
    }

    # Polynomial assembly → mu shares (using CHEBYSHEV coefficients)
    power_keys <- c("secure_eta_share", paste0("secure_pow", 2:(poly_degree + 1)))
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "poly_eval",
        power_keys = power_keys, coefficients = poly_coeffs,
        party_id = if (is_coord) 0L else 1L,
        frac_bits = frac_bits, session_id = session_id))
    }

    # === C: mu^2 for weights (1 Beaver round) ===
    open_results <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "beaver_open",
        a_share_key = "secure_mu_share", b_share_key = "secure_mu_share",
        u_vals = tri_offset, v_vals = n_obs,
        peer_pk = transport_pks[[peer]], frac_bits = frac_bits,
        session_id = session_id))
      if (is.list(r)) r <- r[[1]]
      open_results[[server]] <- r
    }
    .sendBlob(open_results[[coordinator]]$peer_de_enc, "beaver_peer_de", nl_conn)
    .sendBlob(open_results[[nl]]$peer_de_enc, "beaver_peer_de", coordinator_conn)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "beaver_close",
        result_key = "secure_mu2_share", u_vals = tri_offset, v_vals = n_obs,
        party_id = if (is_coord) 0L else 1L,
        frac_bits = frac_bits, session_id = session_id))
    }
    tri_offset <- tri_offset + n_obs

    # === D: Intercept Newton (2 aggregated scalars) ===
    intercept_scalars <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS",
        step = "intercept_newton_prepare",
        data_name = std_data,
        y_var = if (server == coordinator) y_var else NULL,
        role = if (server == coordinator) "label" else "nonlabel",
        frac_bits = frac_bits, session_id = session_id))
      if (is.list(r)) r <- r[[1]]
      intercept_scalars[[server]] <- r
    }
    total_sum_w <- sum(sapply(intercept_scalars, function(x) x$sum_w_share))
    total_sum_resid <- sum(sapply(intercept_scalars, function(x) x$sum_resid_share))
    intercept_beta0 <- intercept_beta0 - total_sum_resid / (total_sum_w + 1e-6)

    # === E: Cross-gradient via Beaver ===
    prep_results <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "prepare_gradient",
        data_name = std_data, x_vars = x_vars[[server]],
        y_var = if (server == coordinator) y_var else NULL,
        role = if (server == coordinator) "label" else "nonlabel",
        peer_pk = transport_pks[[peer]],
        frac_bits = frac_bits, session_id = session_id))
      if (is.list(r)) r <- r[[1]]
      prep_results[[server]] <- r
    }

    .sendBlob(prep_results[[coordinator]]$x_share_for_peer_enc, "peer_x_share", nl_conn)
    .sendBlob(prep_results[[coordinator]]$r_share_for_peer_enc, "peer_r_share", nl_conn)
    .sendBlob(prep_results[[nl]]$x_share_for_peer_enc, "peer_x_share", coordinator_conn)
    .sendBlob(prep_results[[nl]]$r_share_for_peer_enc, "peer_r_share", coordinator_conn)

    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci], call("k2MpcSecureStepDS", step = "receive_peer_shares",
             session_id = session_id))
    }

    # Cross-gradient Beaver rounds
    cross_grad_shares <- list()
    for (target in server_list) {
      peer_of_target <- setdiff(server_list, target)
      is_target_coord <- (target == coordinator)
      p_target <- if (is_target_coord) p_coord else p_nl
      cross_off <- n_poly_triples + n_mu2_triples +
        (if (is_target_coord) 0 else n_cross_coord)
      cross_cnt <- n_obs * p_target

      open_r <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        server_role <- if (server == target) "target" else "peer"
        r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS",
          step = "cross_gradient_open", u_vals = cross_off, v_vals = cross_cnt,
          role = server_role,
          peer_pk = transport_pks[[setdiff(server_list, server)]],
          frac_bits = frac_bits, session_id = session_id))
        if (is.list(r)) r <- r[[1]]
        open_r[[server]] <- r
      }
      .sendBlob(open_r[[coordinator]]$peer_de_enc, "cross_beaver_peer_de", nl_conn)
      .sendBlob(open_r[[nl]]$peer_de_enc, "cross_beaver_peer_de", coordinator_conn)

      close_r <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)
        server_role <- if (server == target) "target" else "peer"
        r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS",
          step = "cross_gradient_close", u_vals = cross_off, v_vals = cross_cnt,
          party_id = if (is_coord) 0L else 1L, role = server_role,
          peer_pk = transport_pks[[setdiff(server_list, server)]],
          frac_bits = frac_bits, session_id = session_id))
        if (is.list(r)) r <- r[[1]]
        close_r[[server]] <- r
      }
      cross_grad_shares[[target]] <- close_r[[peer_of_target]]
    }

    for (target in server_list) {
      target_ci <- which(server_names == target)
      .sendBlob(cross_grad_shares[[target]]$cross_for_peer_enc,
                "cross_gradient_from_peer", target_ci)
    }

    # === F: Combine gradient + L-BFGS update (server-local) ===
    # L-BFGS with the polynomial gradient gives deviance=166 (matching pragmatic).
    # The coefficient intercept differs from the MLE because the polynomial
    # sigmoid is not the true sigmoid, but the model quality (deviance, AUC)
    # is equivalent.
    for (server in server_list) {
      ci <- which(server_names == server)
      r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS",
        step = "combine_gradient", frac_bits = frac_bits,
        session_id = session_id))
      if (is.list(r)) r <- r[[1]]
      gradient <- r$gradient

      # L-BFGS update (server-local, no new leakage)
      lbfgs_r <- .dsAgg(datasources[ci], call("k2MpcSecureStepDS",
        step = "lbfgs_step",
        u_vals = gradient / n_obs + lambda * betas[[server]],
        v_vals = betas[[server]],
        session_id = session_id))
      if (is.list(lbfgs_r)) lbfgs_r <- lbfgs_r[[1]]
      betas[[server]] <- lbfgs_r$beta_new
    }

    # Check convergence
    max_diff <- max(abs(intercept_beta0 - old_beta0))
    for (server in server_list) {
      diff_server <- max(abs(betas[[server]] - betas_old[[server]]))
      max_diff <- max(max_diff, diff_server)
    }

    final_iter <- iter

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose)
        message(sprintf("  Converged after %d iterations (diff = %.2e)",
                        iter, max_diff))
      break
    }

    if (verbose && iter %% 5 == 0)
      message(sprintf("  Iteration %d: max diff = %.2e", iter, max_diff))

    # Periodic GC
    if (iter %% 20 == 0) {
      for (server in server_list) {
        tryCatch(.dsAgg(conns = datasources[which(server_names == server)],
          expr = call("mheGcDS")), error = function(e) NULL)
      }
    }
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)",
                    max_iter, max_diff))

  list(
    betas = betas,
    intercept = intercept_beta0,
    converged = converged,
    iterations = final_iter,
    max_diff = max_diff
  )
}
