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

  # Helper: convert base64 to base64url for passing through DataSHIELD
  .to_b64url <- function(x) gsub("+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  .from_b64url <- function(x) {
    x <- gsub("-","+",gsub("_","/",x,fixed=TRUE),fixed=TRUE)
    pad <- nchar(x) %% 4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }

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
    if (verbose && iter == 1) {
      message("  [TRACE] Step 1: compute eta shares...")
      # Diagnostic: check session state before computing eta
      for (server in server_list) {
        ci <- which(server_names == server)
        diag <- tryCatch(
          .dsAgg(datasources[ci], call("k2DiagnosticDS", session_id = session_id)),
          error = function(e) list(error = conditionMessage(e))
        )
        if (is.list(diag) && length(diag) == 1) diag <- diag[[1]]
        message("  [DIAG] ", server, ": ", paste(names(diag), "=", unlist(diag), collapse=", "))
      }
    }
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

    if (verbose && iter == 1) {
      for (server in server_list) {
        ci <- which(server_names == server)
        diag <- tryCatch(
          .dsAgg(datasources[ci], call("k2DiagnosticDS", session_id = session_id)),
          error = function(e) list(error = conditionMessage(e))
        )
        if (is.list(diag) && length(diag) == 1) diag <- diag[[1]]
        message("  [DIAG post-eta] ", server, ": has_x_full=", diag$has_x_full,
                " has_mu=", diag$has_mu, " eta_len=", diag$x_full_len,
                " n=", diag$n, " p_own=", diag$p_own, " p_peer=", diag$p_peer)
      }
    }

    # === Step 2: Beaver polynomial eval for mu shares (ALL IN FP) ===
    # Uses k2BeaverRoundFPDS instead of old store_triples + beaver_open/close.
    # Everything stays in FixedPoint base64 — no float64 conversion.

    # Generate one Beaver triple per power chain round (degree-1 triples)
    if (verbose && iter == 1) message("  [TRACE] Step 2: Beaver power chain (", poly_degree, " rounds)...")
    for (k in seq_len(poly_degree)) {
      if (verbose && iter == 1) message("  [TRACE]   Beaver round k=", k)
      a_key <- if (k == 1) "secure_eta_share" else paste0("k2_pow", k)
      b_key <- "secure_eta_share"
      result_key <- paste0("k2_pow", k + 1)

      # Generate FP triple for this round
      triple <- dsVert:::.callMheTool("k2-gen-beaver-triples", list(
        n = as.integer(n_obs), frac_bits = as.integer(frac_bits)))

      # Send triple shares as encrypted blobs (NOT as function args)
      # This avoids the DataSHIELD expression parser which can corrupt base64 data
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)
        party_idx <- if (is_coord) "party0" else "party1"

        pk_b64 <- .b64url_to_b64(transport_pks[[server]])
        triple_json <- jsonlite::toJSON(list(
          a = triple[[paste0(party_idx, "_u")]],
          b = triple[[paste0(party_idx, "_v")]],
          c = triple[[paste0(party_idx, "_w")]]), auto_unbox = TRUE)
        msg_b64 <- jsonlite::base64_enc(charToRaw(triple_json))
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(data = msg_b64, recipient_pk = pk_b64))
        .sendBlob(gsub("+","-",gsub("/","_",gsub("=+$","",sealed$sealed,perl=TRUE),fixed=TRUE),fixed=TRUE),
                  paste0("k2_pow_triple_", k), ci)
      }

      # Each party stores the triple from blob, then computes Phase 1
      if (verbose && iter == 1) message("  [TRACE]     Phase 1...")
      r1_results <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)

        # Store triple from blob
        .dsAgg(datasources[ci], call("k2StorePowerTripleDS",
          triple_key = paste0("k2_pow_triple_", k),
          session_id = session_id))

        # Phase 1: compute (X-A, Y-B) using triple from session
        r <- .dsAgg(datasources[ci], call("k2BeaverRoundFPDS",
          x_key = a_key, y_key = b_key,
          a_fp = "NONE", b_fp = "NONE",
          party_id = if (is_coord) 0L else 1L,
          phase = 1L,
          use_session_triple = 1L,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]
        r1_results[[server]] <- r
      }

      # Relay Phase 1 messages (peer's xma/ymb) as blobs
      if (verbose && iter == 1) message("  [TRACE]     Relay...")
      for (server in server_list) {
        peer <- setdiff(server_list, server)
        peer_ci <- which(server_names == peer)
        msg_json <- jsonlite::toJSON(list(
          xma = r1_results[[server]]$xma_fp,
          ymb = r1_results[[server]]$ymb_fp), auto_unbox = TRUE)
        msg_b64 <- jsonlite::base64_enc(charToRaw(msg_json))
        pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(data = msg_b64, recipient_pk = pk_b64))
        .sendBlob(gsub("+","-",gsub("/","_",gsub("=+$","",sealed$sealed,perl=TRUE),fixed=TRUE),fixed=TRUE),
                  paste0("k2_br_peer_", k), peer_ci)
      }

      # Phase 2: compute Beaver close using triple from session + peer msg from blob
      if (verbose && iter == 1) message("  [TRACE]     Phase 2...")
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)

        .dsAgg(datasources[ci], call("k2BeaverRoundFPDS",
          x_key = a_key, y_key = b_key,
          a_fp = "NONE", b_fp = "NONE", c_fp = "NONE",
          peer_xma_fp = "NONE", peer_ymb_fp = "NONE",
          result_key = result_key,
          party_id = if (is_coord) 0L else 1L,
          phase = 2L,
          use_session_triple = 1L,
          peer_blob_key = paste0("k2_br_peer_", k),
          session_id = session_id))
      }

      # Round-by-round power share validation (iter 1 only, beta=0 → eta=0 → all powers=0)
      if (verbose && iter == 1) {
        shares <- list()
        for (server in server_list) {
          ci <- which(server_names == server)
          d <- .dsAgg(datasources[ci], call("k2ReadSessionKeyDS",
            key = result_key, session_id = session_id))
          if (is.list(d) && length(d) == 1) d <- d[[1]]
          shares[[server]] <- d
        }
        # Aggregate to cleartext using k2-ring63-aggregate
        agg <- tryCatch({
          dsVert:::.callMheTool("k2-ring63-aggregate", list(
            share_a = shares[[server_list[1]]]$value,
            share_b = shares[[server_list[2]]]$value,
            frac_bits = 20L))
        }, error = function(e) list(values = paste("ERROR:", conditionMessage(e))))
        vals <- if (is.numeric(agg$values)) agg$values else NA
        message("  [VALIDATE] Round k=", k, " key=", result_key,
                " nchar=", shares[[server_list[1]]]$nchar, "/", shares[[server_list[2]]]$nchar,
                " sum=", if (is.numeric(vals)) sum(vals) else "NA",
                " first3=", paste(head(vals, 3), collapse=","),
                " expect=0 (beta=0)")
      }
    }

    if (verbose && iter == 1) message("  [TRACE] Step 2b: Polynomial eval...")
    # Polynomial eval → mu shares (FP)
    # Uses the existing mpc-secure-poly-eval (which reads FP shares correctly)
    power_keys <- c("secure_eta_share", paste0("k2_pow", 2:(poly_degree + 1)))
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2PolyEvalLocalFPDS",
        power_keys = power_keys, coefficients = poly_coeffs,
        party_id = if (is_coord) 0L else 1L,
        session_id = session_id))
    }

    # Diagnostic: reconstruct mu from shares (iter 1 only)
    if (verbose && iter == 1) {
      mu_diags <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        d <- .dsAgg(datasources[ci], call("k2DiagnosticDS", session_id = session_id))
        if (is.list(d) && length(d) == 1) d <- d[[1]]
        mu_diags[[server]] <- d
      }
      s1 <- mu_diags[[server_list[1]]]; s2 <- mu_diags[[server_list[2]]]
      message("  [DIAG] mu_len: ", s1$mu_len, " / ", s2$mu_len)
      if (!is.null(s1$eta_first5) && !is.null(s2$eta_first5)) {
        eta_sum <- s1$eta_first5 + s2$eta_first5
        message("  [DIAG] eta reconstructed (first 5): ", paste(round(eta_sum, 4), collapse=", "))
      }
      if (!is.null(s1$mu_sum_fp) && !is.null(s2$mu_sum_fp)) {
        mu_total <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
          share_a = s1$mu_sum_fp, share_b = s2$mu_sum_fp, frac_bits = frac_bits))
        message("  [DIAG] sum(mu) via Ring63 = ", mu_total$values[1], " (expected ~", n_obs*0.5, ")")
      }
      if (!is.null(s1$eta_sum_fp) && !is.null(s2$eta_sum_fp)) {
        eta_total <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
          share_a = s1$eta_sum_fp, share_b = s2$eta_sum_fp, frac_bits = frac_bits))
        message("  [DIAG] sum(eta) via Ring63 = ", eta_total$values[1], " (expected 0)")
      }
      if (!is.null(s1$y_first5) && !is.null(s2$y_first5)) {
        y_sum <- s1$y_first5 + s2$y_first5
        message("  [DIAG] y reconstructed (first 5): ", paste(round(y_sum, 4), collapse=", "))
      }
    }

    if (verbose && iter == 1) message("  [TRACE] Step 3: Gradient matvec...")
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

    # === Step 4: Client sums gradient shares in Ring63, then converts to float ===
    # Sum gradient shares in Ring63 to avoid catastrophic cancellation in float64.
    # Each party returns gradient_fp (Ring63 values as base64 FP).
    grad_fp_coord <- grad_results[[coordinator]]$gradient_fp
    grad_fp_nl <- grad_results[[nl]]$gradient_fp
    res_fp_coord <- r1_results[[coordinator]]$sum_residual_fp
    res_fp_nl <- r1_results[[nl]]$sum_residual_fp

    if (!is.null(grad_fp_coord) && !is.null(grad_fp_nl)) {
      # Use Go mhe-tool to sum Ring63 shares and convert to float
      agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a = grad_fp_coord, share_b = grad_fp_nl, frac_bits = frac_bits))
      gradient <- agg$values

      agg_res <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a = res_fp_coord, share_b = res_fp_nl, frac_bits = frac_bits))
      sum_residual <- agg_res$values[1]
    } else {
      # Fallback to float64 addition (less precise)
      gradient <- rep(0, p_total)
      sum_residual <- 0
      for (server in server_list) {
        gradient <- gradient + grad_results[[server]]$gradient_share
        sum_residual <- sum_residual + r1_results[[server]]$sum_residual
      }
    }

    # Debug: print gradient info for first iteration
    if (iter == 1 && verbose) {
      message(sprintf("  [DEBUG iter 1] gradient = %s", paste(round(gradient, 4), collapse=", ")))
      message(sprintf("  [DEBUG iter 1] sum_residual = %.4f", sum_residual))
      message(sprintf("  [DEBUG iter 1] coord grad_share = %s",
        paste(round(grad_results[[coordinator]]$gradient_share, 4), collapse=", ")))
      message(sprintf("  [DEBUG iter 1] nl grad_share = %s",
        paste(round(grad_results[[nl]]$gradient_share, 4), collapse=", ")))
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
