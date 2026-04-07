#' @title K=2 Beaver Newton-IRLS Pipeline
#' @description Secure K=2 GLM with Newton-IRLS and diagonal Fisher preconditioning.
#'   Gaussian: identity link, constant Fisher (2 rounds/iter, no spline/DCF).
#'   Binomial/Poisson: wide spline sigmoid/exp (9 rounds/iter).
#' @name k2-beaver-newton-client
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
  is_gaussian <- (family == "gaussian")
  num_intervals <- if (family == "poisson") 100L else 50L
  alpha <- switch(family, poisson = 0.05, gaussian = 0.5, 0.3)
  if (is.null(lambda)) lambda <- 1e-4
  if (verbose) {
    if (is_gaussian) {
      message("  Beaver Newton-IRLS: gaussian, identity link, constant Fisher")
    } else {
      message(sprintf("  Wide spline: %s, %d intervals, Newton-IRLS", family, num_intervals))
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
  xsq_precomputed <- FALSE
  intercept <- 0.0
  converged <- FALSE
  final_iter <- 0
  diag_fisher <- NULL  # For Gaussian: pre-computed once; for others: per-iteration

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

    # Pre-compute x² on first iteration (needed for Fisher, all families)
    if (!xsq_precomputed) {
      if (verbose) message("  Pre-computing x\u00b2 for diagonal Fisher...")
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
      xsq_precomputed <- TRUE
      if (verbose) message("  x\u00b2 pre-computed")

      # Gaussian: compute constant Fisher once (d_j = sum(x²_j), w=1)
      if (is_gaussian) {
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
        diag_fisher <- fisher_agg$values / n_obs + lambda
        if (verbose) message(sprintf("  [Fisher] constant d_j = [%s]",
          paste(round(diag_fisher, 4), collapse=", ")))
      }
    }

    # === Step 2: Link function ===
    if (is_gaussian) {
      # Identity link: mu = eta (no spline, no DCF, no Beaver)
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

      # Phase 1: DCF masked values
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

      # Phase 4: Had2 close + assembly -> mu
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)
        .dsAgg(datasources[ci], call("k2WideSplinePhase4DS",
          party_id = if(is_coord) 0L else 1L, family = family,
          num_intervals = num_intervals, frac_bits = frac_bits,
          session_id = session_id))
      }
    }

    # === Step 2.5: Diagonal Fisher (binomial/poisson only — Gaussian uses pre-computed) ===
    if (!is_gaussian) {
      # 3-phase protocol: w R1 -> w close + w*x^2 R1 -> w*x^2 close + d_j
      w_triple <- dsVert:::.callMheTool("k2-gen-beaver-triples",
        list(n = as.integer(n_obs), frac_bits = frac_bits))
      wx_triple <- dsVert:::.callMheTool("k2-gen-beaver-triples",
        list(n = as.integer(n_obs * p_total), frac_bits = frac_bits))

      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == coordinator)
        pk_b64 <- .b64url_to_b64(transport_pks[[server]])
        party_idx <- if (is_coord) "party0" else "party1"
        td_w <- list(a=w_triple[[paste0(party_idx,"_u")]],b=w_triple[[paste0(party_idx,"_v")]],c=w_triple[[paste0(party_idx,"_w")]])
        sealed_w <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(td_w,auto_unbox=TRUE))),recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed_w$sealed), "k2_fisher_w_triple", ci)
        td_wx <- list(a=wx_triple[[paste0(party_idx,"_u")]],b=wx_triple[[paste0(party_idx,"_v")]],c=wx_triple[[paste0(party_idx,"_w")]])
        sealed_wx <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(td_wx,auto_unbox=TRUE))),recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed_wx$sealed), "k2_fisher_wx_triple", ci)
      }

      f_ph1 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call("k2RealFisherPhase1DS", family = family,
          party_id=if(is_coord) 0L else 1L, frac_bits=frac_bits, session_id=session_id))
        if (is.list(r) && length(r)==1) r <- r[[1]]; f_ph1[[server]] <- r
      }
      for (server in server_list) {
        peer <- setdiff(server_list,server); peer_ci <- which(server_names==peer)
        pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(w_xma=f_ph1[[server]]$w_xma,w_ymb=f_ph1[[server]]$w_ymb),auto_unbox=TRUE))),
          recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_fisher_w_peer_r1", peer_ci)
      }

      f_ph2 <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call("k2RealFisherPhase2DS", family = family,
          party_id=if(is_coord) 0L else 1L, frac_bits=frac_bits,
          p_total=as.integer(p_total), session_id=session_id))
        if (is.list(r) && length(r)==1) r <- r[[1]]; f_ph2[[server]] <- r
      }
      for (server in server_list) {
        peer <- setdiff(server_list,server); peer_ci <- which(server_names==peer)
        pk_b64 <- .b64url_to_b64(transport_pks[[peer]])
        sealed <- dsVert:::.callMheTool("transport-encrypt", list(
          data=jsonlite::base64_enc(charToRaw(jsonlite::toJSON(list(wx_xma=f_ph2[[server]]$wx_xma,wx_ymb=f_ph2[[server]]$wx_ymb),auto_unbox=TRUE))),
          recipient_pk=pk_b64))
        .sendBlob(.to_b64url(sealed$sealed), "k2_fisher_wx_peer_r1", peer_ci)
      }

      fisher_shares <- list()
      for (server in server_list) {
        ci <- which(server_names == server); is_coord <- (server == coordinator)
        r <- .dsAgg(datasources[ci], call("k2RealFisherPhase3DS", family = family,
          party_id=if(is_coord) 0L else 1L, frac_bits=frac_bits,
          p_total=as.integer(p_total), session_id=session_id))
        if (is.list(r) && length(r)==1) r <- r[[1]]; fisher_shares[[server]] <- r
      }
      fisher_agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
        share_a = fisher_shares[[server_list[1]]]$fisher_diag_fp,
        share_b = fisher_shares[[server_list[2]]]$fisher_diag_fp,
        frac_bits = frac_bits))
      diag_fisher <- fisher_agg$values / n_obs + lambda
      if (verbose && iter <= 3) message(sprintf("  [Fisher] d_j = [%s]",
        paste(round(diag_fisher, 4), collapse=", ")))
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

    # Newton-IRLS: diagonal Fisher preconditioning for ALL coefficients
    full_grad <- c(sum_residual / n_obs, gradient / n_obs) + lambda * c(intercept, beta)
    intercept_d <- mean(diag_fisher)
    full_fisher <- c(intercept_d, diag_fisher)
    damping <- 0.5
    full_step <- damping * full_grad / full_fisher
    if (verbose && iter <= 3) {
      message(sprintf("  [NEWTON] sum_res=%.4f, int_grad=%.6f, int_fisher=%.6f, int_step=%.6f",
        sum_residual, full_grad[1], full_fisher[1], full_step[1]))
      message(sprintf("  [NEWTON] grad=[%s]", paste(round(full_grad, 6), collapse=", ")))
    }
    intercept <- intercept - full_step[1]
    beta <- beta - full_step[-1]

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
