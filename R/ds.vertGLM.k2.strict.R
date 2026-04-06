#' @title K=2 Wide Spline: Optimized 4-phase Newton-IRLS Pipeline
#' @description Secure K=2 GLM with wide spline sigmoid/exp.
#'   Optimizations: DCF key reuse, bulk Beaver triples, Newton-IRLS (10-15 iter).
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

  frac_bits <- 20L
  num_intervals <- if (family == "poisson") 100L else 30L  # 30 for speed
  alpha <- if (family == "poisson") 0.05 else 0.3
  if (is.null(lambda)) lambda <- 1e-4
  if (verbose) message(sprintf("  Wide spline: %s, %d intervals, Newton-IRLS", family, num_intervals))

  p_coord <- length(x_vars[[coordinator]])
  p_nl <- length(x_vars[[nl]])
  p_total <- p_coord + p_nl

  .to_b64url <- function(x) gsub("+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-","+",gsub("_","/",x,fixed=TRUE),fixed=TRUE)
    pad <- nchar(x)%%4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }

  # === INPUT-SHARING PREAMBLE ===
  if (verbose) message("  Input-sharing preamble...")
  share_results <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    peer_pk_safe <- .to_b64url(transport_pks[[peer]])
    r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
      data_name = std_data, x_vars = x_vars[[server]],
      y_var = if (server == coordinator) y_var else NULL,
      peer_pk = peer_pk_safe, session_id = session_id))
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

  # === PRE-GENERATE DCF KEYS (reused across ALL iterations) ===
  if (verbose) message("  Pre-generating DCF keys (one-time)...")
  dcf <- dsVert:::.callMheTool("k2-dcf-gen-batch", list(
    family = family, n = as.integer(n_obs),
    frac_bits = frac_bits, num_intervals = num_intervals))
  # Send DCF keys to servers ONCE (stored in session for all iterations)
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == coordinator)
    pk_b64 <- .b64url_to_b64(transport_pks[[server]])
    sealed <- dsVert:::.callMheTool("transport-encrypt", list(
      data = if(is_coord) dcf$party0_keys else dcf$party1_keys,
      recipient_pk = pk_b64))
    .sendBlob(.to_b64url(sealed$sealed), "k2_dcf_keys_persistent", ci)
    # Store persistently on server
    .dsAgg(datasources[ci], call("k2StoreDcfKeysPersistentDS", session_id = session_id))
  }
  if (verbose) message("  DCF keys distributed (will reuse across iterations)")

  beta <- rep(0, p_total)
  intercept <- 0.0
  converged <- FALSE
  final_iter <- 0

  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    intercept_old <- intercept

    # === Step 1: Compute eta ===
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2ComputeEtaShareDS",
        beta_coord = beta[1:p_coord], beta_nl = beta[(p_coord+1):p_total],
        intercept = if (is_coord) intercept else 0.0,
        is_coordinator = is_coord, session_id = session_id))
    }

    # === Step 2: Wide spline sigmoid (4-phase, DCF keys reused) ===
    # Only generate Beaver triples per iteration (3 triples, small)
    triples <- lapply(1:3, function(i) dsVert:::.callMheTool("k2-gen-beaver-triples",
      list(n = as.integer(n_obs), frac_bits = frac_bits)))
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      party_idx <- if (is_coord) "party0" else "party1"
      td <- list()
      for (op in c("and","had1","had2")) {
        ti <- switch(op, and=1, had1=2, had2=3)
        td[[paste0(op,"_a")]] <- triples[[ti]][[paste0(party_idx,"_u")]]
        td[[paste0(op,"_b")]] <- triples[[ti]][[paste0(party_idx,"_v")]]
        td[[paste0(op,"_c")]] <- triples[[ti]][[paste0(party_idx,"_w")]]
      }
      sealed_t <- dsVert:::.callMheTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(jsonlite::toJSON(td, auto_unbox=TRUE))),
        recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed_t$sealed), "k2_spline_triples", ci)
    }

    # Phase 1: DCF masked values (uses persistent DCF keys)
    ph1 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase1DS",
        party_id = if(is_coord) 0L else 1L, family = family,
        num_intervals = num_intervals, frac_bits = frac_bits,
        session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      ph1[[server]] <- r
    }
    .sendBlob(ph1[[coordinator]]$dcf_masked, "k2_peer_dcf_masked", nl_conn)
    .sendBlob(ph1[[nl]]$dcf_masked, "k2_peer_dcf_masked", coordinator_conn)

    # Phase 2: DCF close + indicators + (AND + Had1) R1
    ph2 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase2DS",
        party_id = if(is_coord) 0L else 1L, family = family,
        num_intervals = num_intervals, frac_bits = frac_bits,
        session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      ph2[[server]] <- r
    }
    for (server in server_list) {
      peer <- setdiff(server_list, server)
      peer_ci <- which(server_names == peer)
      pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
      r1_json <- jsonlite::toJSON(list(
        and_xma=ph2[[server]]$and_xma, and_ymb=ph2[[server]]$and_ymb,
        had1_xma=ph2[[server]]$had1_xma, had1_ymb=ph2[[server]]$had1_ymb),
        auto_unbox=TRUE)
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data=jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk=pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", peer_ci)
    }

    # Phase 3: AND/Had1 close + Had2 R1
    ph3 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase3DS",
        party_id = if(is_coord) 0L else 1L, family = family,
        num_intervals = num_intervals, frac_bits = frac_bits,
        session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      ph3[[server]] <- r
    }
    for (server in server_list) {
      peer <- setdiff(server_list, server)
      peer_ci <- which(server_names == peer)
      pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
      r1_json <- jsonlite::toJSON(list(
        had2_xma=ph3[[server]]$had2_xma, had2_ymb=ph3[[server]]$had2_ymb),
        auto_unbox=TRUE)
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data=jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk=pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", peer_ci)
    }

    # Phase 4: Had2 close + assembly → mu
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      .dsAgg(datasources[ci], call("k2WideSplinePhase4DS",
        party_id = if(is_coord) 0L else 1L, family = family,
        num_intervals = num_intervals, frac_bits = frac_bits,
        session_id = session_id))
    }

    # === Step 3: Gradient (Beaver matvec) ===
    mvt <- dsVert:::.callMheTool("k2-gen-matvec-triples", list(
      n = as.integer(n_obs), p = as.integer(p_total)))
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      party_idx <- if (is_coord) "party0" else "party1"
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      msg_json <- jsonlite::toJSON(list(a=mvt[[paste0(party_idx,"_a")]],
        b=mvt[[paste0(party_idx,"_b")]], c=mvt[[paste0(party_idx,"_c")]]), auto_unbox=TRUE)
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data=jsonlite::base64_enc(charToRaw(msg_json)), recipient_pk=pk_b64))
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
        party_id = if(is_coord) 0L else 1L, session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      grad_results[[server]] <- r
    }

    # === Step 4: Aggregate gradient + Newton update ===
    grad_fp_coord <- grad_results[[coordinator]]$gradient_fp
    grad_fp_nl <- grad_results[[nl]]$gradient_fp
    res_fp_coord <- r1_results[[coordinator]]$sum_residual_fp
    res_fp_nl <- r1_results[[nl]]$sum_residual_fp
    if (!is.null(grad_fp_coord) && !is.null(grad_fp_nl)) {
      agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a=grad_fp_coord, share_b=grad_fp_nl, frac_bits=frac_bits))
      gradient <- agg$values
      agg_res <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a=res_fp_coord, share_b=res_fp_nl, frac_bits=frac_bits))
      sum_residual <- agg_res$values[1]
    } else {
      gradient <- rep(0, p_total); sum_residual <- 0
      for (server in server_list) {
        gradient <- gradient + grad_results[[server]]$gradient_share
        sum_residual <- sum_residual + r1_results[[server]]$sum_residual
      }
    }

    # Gradient descent update (simple but robust)
    full_grad <- gradient / n_obs + lambda * beta
    grad_norm <- sqrt(sum(full_grad^2))
    if (grad_norm > 5.0) full_grad <- full_grad * (5.0 / grad_norm)
    beta <- beta - alpha * full_grad
    intercept <- intercept - alpha * sum_residual / n_obs

    max_diff <- max(abs(beta - beta_old), abs(intercept - intercept_old))
    final_iter <- iter

    if (verbose) message(sprintf("  Iter %d: ||grad||=%.4f diff=%.2e",
      iter, sqrt(sum(full_grad^2)), max_diff))

    if (max_diff < tol) {
      converged <- TRUE
      if (verbose) message(sprintf("  Converged after %d iterations (diff = %.2e)", iter, max_diff))
      break
    }
    if (iter %% 10 == 0) {
      for (server in server_list)
        tryCatch(.dsAgg(datasources[which(server_names==server)], call("mheGcDS")), error=function(e) NULL)
    }
  }

  if (!converged && verbose)
    warning(sprintf("Did not converge after %d iterations (diff = %.2e)", max_iter, max_diff))

  betas <- list()
  betas[[coordinator]] <- beta[1:p_coord]
  betas[[nl]] <- beta[(p_coord+1):p_total]
  list(betas=betas, intercept=intercept, converged=converged, iterations=final_iter, max_diff=max_diff)
}
