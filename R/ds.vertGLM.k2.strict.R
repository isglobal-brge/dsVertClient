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
                             n_obs, verbose, .dsAgg, .sendBlob) {

  frac_bits <- 20L
  is_gaussian <- (family == "gaussian")
  num_intervals <- if (family == "poisson") 100L else 50L
  if (is.null(lambda)) lambda <- 1e-4
  if (verbose) {
    if (is_gaussian) {
      message("  L-BFGS: gaussian, identity link (2 rounds/iter)")
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

  # === PRE-GENERATE DCF KEYS (only for binomial/poisson) ===
  if (!is_gaussian) {
    if (verbose) message("  Pre-generating DCF keys (one-time)...")
    dcf <- dsVert:::.callMheTool("k2-dcf-gen-batch", list(
      family = family, n = as.integer(n_obs),
      frac_bits = frac_bits, num_intervals = num_intervals))
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == coordinator)
      pk_b64 <- .b64url_to_b64(transport_pks[[server]])
      sealed <- dsVert:::.callMheTool("transport-encrypt", list(
        data = if(is_coord) dcf$party0_keys else dcf$party1_keys,
        recipient_pk = pk_b64))
      .sendBlob(.to_b64url(sealed$sealed), "k2_dcf_keys_persistent", ci)
      .dsAgg(datasources[ci], call("k2StoreDcfKeysPersistentDS", session_id = session_id))
    }
    if (verbose) message("  DCF keys distributed (will reuse across iterations)")
  }

  beta <- rep(0, p_total)
  intercept <- 0.0
  converged <- FALSE
  final_iter <- 0

  # L-BFGS state (client-side only)
  lbfgs_m <- 7L  # number of (s,y) pairs to keep
  s_hist <- list()
  y_hist <- list()
  prev_theta <- NULL
  prev_grad <- NULL
  eta_computed <- FALSE  # track whether first eta has been computed (for x² precompute)

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

    # Gaussian: pre-compute x^2 once (first iteration, after eta creates X_full)
    if (is_gaussian && !eta_computed) {
      if (verbose) message("  Pre-computing x\u00b2 for Gaussian Fisher...")
      xsq_t <- dsVert:::.callMheTool("k2-gen-beaver-triples",
        list(n = as.integer(n_obs * p_total), frac_bits = frac_bits))
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        pk_b64 <- .b64url_to_b64(transport_pks[[server]])
        party_idx <- if (is_coord) "party0" else "party1"
        td <- list(a=xsq_t[[paste0(party_idx,"_u")]],b=xsq_t[[paste0(party_idx,"_v")]],c=xsq_t[[paste0(party_idx,"_w")]])
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(td,auto_unbox=TRUE))),recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_xsq_triples", ci)
      }
      xsq_ph1 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call("k2PrecomputeXSqPhase1DS",
          party_id=if(is_coord) 0L else 1L, frac_bits=frac_bits,
          p_total=as.integer(p_total), session_id=session_id))
        if (is.list(r) && length(r)==1) r <- r[[1]]; xsq_ph1[[server]] <- r
      }
      for (server in server_list) {
        peer <- setdiff(server_list,server); peer_ci <- which(server_names==peer)
        pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(xma=xsq_ph1[[server]]$xma,ymb=xsq_ph1[[server]]$ymb),auto_unbox=TRUE))),
          recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_xsq_peer_r1", peer_ci)
      }
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        .dsAgg(datasources[ci], call("k2PrecomputeXSqPhase2DS",
          party_id=if(is_coord) 0L else 1L, frac_bits=frac_bits,
          p_total=as.integer(p_total), session_id=session_id))
      }
      if (verbose) message("  x\u00b2 pre-computed")

      # Gaussian constant Fisher: d_j = sum(x^2_j)/n + lambda
      fisher_shares <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        r <- .dsAgg(datasources[ci], call("k2GaussianFisherDS",
          p_total = as.integer(p_total), frac_bits = frac_bits,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]
        fisher_shares[[server]] <- r
      }
      fisher_agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a = fisher_shares[[server_list[1]]]$fisher_diag_fp,
        share_b = fisher_shares[[server_list[2]]]$fisher_diag_fp,
        frac_bits = frac_bits))
      diag_fisher_gauss <- fisher_agg$values / n_obs + lambda
      if (verbose) message(sprintf("  [Fisher] constant d_j = [%s]",
        paste(round(diag_fisher_gauss, 4), collapse=", ")))
      eta_computed <- TRUE
    }
    if (!is_gaussian) eta_computed <- TRUE

    # === Step 2: Link function ===
    if (is_gaussian) {
      for (server in server_list) {
        ci <- which(server_names == server)
        .dsAgg(datasources[ci], call("k2IdentityLinkDS", session_id = session_id))
      }
    } else {
      # Binomial/Poisson: 4-phase wide spline sigmoid/exp
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

      ph1 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call("k2WideSplinePhase1DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; ph1[[server]] <- r
      }
      .sendBlob(ph1[[coordinator]]$dcf_masked, "k2_peer_dcf_masked", nl_conn)
      .sendBlob(ph1[[nl]]$dcf_masked, "k2_peer_dcf_masked", coordinator_conn)

      ph2 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call("k2WideSplinePhase2DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
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
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", peer_ci)
      }

      ph3 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call("k2WideSplinePhase3DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          session_id = session_id))
        if (is.list(r) && length(r) == 1) r <- r[[1]]; ph3[[server]] <- r
      }
      for (server in server_list) {
        peer <- setdiff(server_list, server); peer_ci <- which(server_names == peer)
        pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
        r1_json <- jsonlite::toJSON(list(
          had2_xma=ph3[[server]]$had2_xma, had2_ymb=ph3[[server]]$had2_ymb),
          auto_unbox=TRUE)
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(r1_json)), recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", peer_ci)
      }

      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        .dsAgg(datasources[ci], call("k2WideSplinePhase4DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          session_id = session_id))
      }
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

    # === Step 4: Aggregate gradient + L-BFGS update ===
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

    # Full gradient with L2 regularization
    theta <- c(intercept, beta)
    full_grad <- c(sum_residual / n_obs, gradient / n_obs) + lambda * theta

    if (verbose && iter <= 3) {
      message(sprintf("  [L-BFGS] sum_res=%.4f, ||grad||=%.6f",
        sum_residual, sqrt(sum(full_grad^2))))
      message(sprintf("  [L-BFGS] grad=[%s]", paste(round(full_grad, 6), collapse=", ")))
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
    new_theta <- theta + step_size * direction
    intercept <- new_theta[1]
    beta <- new_theta[-1]

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
