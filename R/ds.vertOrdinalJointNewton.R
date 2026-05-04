.ord_joint_secure_fit <- function(formula, data = NULL, levels_ordered,
                                  cumulative_template = "%s_leq",
                                  max_outer = 8L, tol = 1e-4,
                                  verbose = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (length(datasources) < 2L) {
    stop("Strict non-disclosive ordinal joint requires at least two servers.",
         call. = FALSE)
  }
  if (!is.character(levels_ordered) || length(levels_ordered) < 3L)
    stop("levels_ordered must have >= 3 levels", call. = FALSE)

  K <- length(levels_ordered)
  K_minus_1 <- K - 1L
  thresh_levels <- head(levels_ordered, -1L)
  rhs <- attr(terms(formula), "term.labels")

  warm <- .ds_vertOrdinalWarm(formula, data = data,
                              levels_ordered = levels_ordered,
                              cumulative_template = cumulative_template,
                              verbose = FALSE,
                              datasources = datasources)
  if (is.null(warm$beta_po) || is.null(warm$thresholds)) {
    stop("Warm ordinal start did not return beta_po/thresholds",
         call. = FALSE)
  }
  beta0 <- -as.numeric(warm$beta_po)
  names(beta0) <- names(warm$beta_po)
  theta0 <- as.numeric(warm$thresholds)
  names(theta0) <- names(warm$thresholds)
  if (any(!is.finite(theta0)) || any(diff(theta0) <= 1e-4)) {
    theta0 <- stats::qlogis(seq_len(K_minus_1) / K)
    names(theta0) <- thresh_levels
  }

  y_var_char <- .ds_gee_extract_lhs(formula)
  server_names <- names(datasources)
  y_server <- .ds_gee_find_server_holding(datasources, server_names,
                                          data, y_var_char)
  if (is.null(y_server))
    stop("outcome server for ", y_var_char, " not found", call. = FALSE)

  session_id <- paste0("ordinalStrict_", as.integer(Sys.time()),
                       "_", sample.int(.Machine$integer.max, 1L))
  .dsAgg <- function(conns, expr, ...) DSI::datashield.aggregate(conns, expr, ...)
  .unwrap <- function(x) if (is.list(x) && length(x) == 1L) x[[1L]] else x
  .json_to_b64url <- function(x) {
    raw <- charToRaw(jsonlite::toJSON(x, auto_unbox = TRUE))
    b64 <- gsub("\n", "", jsonlite::base64_enc(raw), fixed = TRUE)
    chartr("+/", "-_", sub("=+$", "", b64, perl = TRUE))
  }
  .to_b64url <- function(x) chartr("+/", "-_", sub("=+$", "", x, perl = TRUE))
  .sendBlob <- function(blob, key, conn_idx) {
    if (is.null(blob) || !nzchar(blob)) return(invisible())
    conn <- datasources[conn_idx]
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        DSI::datashield.aggregate(conn,
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               session_id = session_id))
      } else {
        DSI::datashield.aggregate(conn,
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               chunk_index = chunk_idx, n_chunks = n_chunks,
               session_id = session_id))
      }
    })
  }

  x_vars_by_server <- list()
  for (srv in server_names) {
    ci <- which(server_names == srv)
    cols <- tryCatch(.unwrap(.dsAgg(datasources[ci],
      call(name = "dsvertColNamesDS", data_name = data))),
      error = function(e) NULL)
    cols_here <- if (is.list(cols) && !is.null(cols$columns)) cols$columns
                 else if (is.character(cols)) cols else character(0)
    x_vars_by_server[[srv]] <- intersect(rhs, cols_here)
  }
  if (length(server_names) == 2L) {
    nl <- setdiff(server_names, y_server)[1L]
  } else {
    if (!exists(".k3_select_fusion_server", mode = "function")) {
      stop("K>=3 strict ordinal joint requires the K>=3 DCF helper selection code",
           call. = FALSE)
    }
    nl <- .k3_select_fusion_server(server_names, y_server, x_vars_by_server)
  }
  if (is.na(nl) || is.null(nl))
    stop("non-label/fusion server not found", call. = FALSE)

  server_list <- c(y_server, nl)
  non_dcf_servers <- setdiff(server_names, server_list)
  feature_server_order <- c(y_server, nl, non_dcf_servers)
  x_vars_per_server <- x_vars_by_server[feature_server_order]
  p_coord <- length(x_vars_per_server[[y_server]])
  p_fusion <- length(x_vars_per_server[[nl]])
  p_extras <- sum(vapply(x_vars_per_server[non_dcf_servers],
                         length, integer(1L)))
  ci_os <- which(server_names == y_server)
  ci_nl <- which(server_names == nl)
  dealer_srv <- if (length(non_dcf_servers) > 0L) non_dcf_servers[[1L]] else nl
  dealer_ci <- which(server_names == dealer_srv)

  if (!setequal(unlist(x_vars_per_server, use.names = FALSE), names(beta0))) {
    stop("Strict ordinal joint requires numeric RHS variables present as ",
         "explicit columns across the participating servers.",
         call. = FALSE)
  }

  .fp_const <- function(x) {
    .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
      values = array(as.numeric(x), dim = 1L),
      frac_bits = 50L, ring = "ring127"))$fp_data)
  }
  one_fp <- .fp_const(1)
  scale_log_fp <- .fp_const(1000)
  neg_log_scale_fp <- .fp_const(-log(1000))

  transport_pks <- list()
  identity_info <- list()
  for (srv in server_names) {
    ci <- which(server_names == srv)
    r <- .unwrap(.dsAgg(datasources[ci],
      call(name = "glmRing63TransportInitDS", session_id = session_id)))
    transport_pks[[srv]] <- r$transport_pk
    if (!is.null(r$identity_pk)) {
      identity_info[[srv]] <- list(identity_pk = r$identity_pk,
                                   signature = r$signature)
    }
  }
  pk_b64 <- .json_to_b64url(transport_pks[sort(names(transport_pks))])
  id_b64 <- if (length(identity_info) > 0L) {
    .json_to_b64url(identity_info[sort(names(identity_info))])
  } else ""
  for (srv in server_names) {
    ci <- which(server_names == srv)
    .dsAgg(datasources[ci],
      call(name = "mpcStoreTransportKeysDS",
           transport_keys_b64 = pk_b64,
           identity_info_b64 = id_b64,
           session_id = session_id))
  }

  share_results <- list()
  for (srv in server_list) {
    ci <- which(server_names == srv)
    peer <- setdiff(server_list, srv)
    r <- .unwrap(.dsAgg(datasources[ci], call(name = "k2ShareInputDS",
      data_name = data, x_vars = x_vars_per_server[[srv]],
      y_var = NULL, peer_pk = .to_b64url(transport_pks[[peer]]),
      ring = 127L, session_id = session_id)))
    share_results[[srv]] <- r
  }
  for (srv in server_list) {
    peer <- setdiff(server_list, srv)
    peer_ci <- which(server_names == peer)
    .sendBlob(share_results[[srv]]$encrypted_x_share,
              "k2_peer_x_share", peer_ci)
  }
  for (srv in server_list) {
    ci <- which(server_names == srv)
    peer <- setdiff(server_list, srv)
    .dsAgg(datasources[ci], call(name = "k2ReceiveShareDS",
      peer_p = as.integer(length(x_vars_per_server[[peer]])),
      session_id = session_id))
  }
  for (srv in non_dcf_servers) {
    if (length(x_vars_per_server[[srv]]) == 0L) next
    ci <- which(server_names == srv)
    r <- .unwrap(.dsAgg(datasources[ci], call(name = "k2ShareInputDS",
      data_name = data, x_vars = x_vars_per_server[[srv]],
      y_var = NULL,
      peer_pk = .to_b64url(transport_pks[[nl]]),
      ring = 127L, session_id = session_id)))
    .sendBlob(r$encrypted_x_share, paste0("k2_extra_x_share_", srv), ci_nl)

    r2 <- .unwrap(.dsAgg(datasources[ci],
      call(name = "glmRing63ExportOwnShareDS",
           peer_pk = .to_b64url(transport_pks[[y_server]]),
           session_id = session_id)))
    .sendBlob(r2$encrypted_own_share, paste0("k2_extra_x_share_", srv), ci_os)
  }
  for (srv in non_dcf_servers) {
    extra_p <- length(x_vars_per_server[[srv]])
    if (extra_p == 0L) next
    for (dcf_srv in server_list) {
      .dsAgg(datasources[which(server_names == dcf_srv)],
        call(name = "glmRing63ReceiveExtraShareDS",
             extra_key = paste0("k2_extra_x_share_", srv),
             extra_p = as.integer(extra_p),
             session_id = session_id))
    }
  }
  n_obs <- as.integer(share_results[[y_server]]$n)

  indicator_cols_vec <- sprintf(cumulative_template, thresh_levels)
  mask_prefix <- "ord_strict_mask"
  mask_r <- .unwrap(.dsAgg(datasources[ci_os],
    call(name = "dsvertOrdinalShareClassMasksDS",
         data_name = data,
         indicator_cols = indicator_cols_vec,
         level_names = levels_ordered,
         output_prefix = mask_prefix,
         peer_pk = .to_b64url(transport_pks[[nl]]),
         session_id = session_id)))
  mask_keys <- as.character(mask_r$mask_keys)
  for (kk in seq_along(mask_keys)) {
    blob_key <- sprintf("ord_strict_mask_blob_%d", kk)
    .sendBlob(mask_r$mask_blobs[[kk]], blob_key, ci_nl)
    .dsAgg(datasources[ci_nl],
      call(name = "dsvertOrdinalReceiveClassMaskDS",
           mask_blob_key = blob_key,
           output_key = mask_keys[kk],
           session_id = session_id))
  }

  if (verbose) {
    message(sprintf("[OrdinalJointStrict] session %s n=%d p=%d class_counts=[%s]",
                    session_id, n_obs, length(beta0),
                    paste(mask_r$class_counts, collapse = ",")))
  }

  .affine <- function(a_key = NULL, b_key = NULL, sign_a = 0L, sign_b = 0L,
                      const_fp = NULL, output_key) {
    for (srv in server_list) {
      ci <- which(server_names == srv)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = a_key, b_key = b_key,
        sign_a = as.integer(sign_a), sign_b = as.integer(sign_b),
        public_const_fp = const_fp,
        is_party0 = (srv == y_server),
        output_key = output_key, n = as.numeric(n_obs),
        session_id = session_id))
    }
    invisible(output_key)
  }
  .sum_keys <- function(keys, output_key) {
    if (length(keys) < 1L) stop("no keys to sum", call. = FALSE)
    .affine(a_key = keys[1L], sign_a = 1L, output_key = output_key)
    if (length(keys) >= 2L) {
      for (kk in keys[-1L]) {
        tmp <- paste0(output_key, "_tmp_", match(kk, keys))
        .affine(a_key = output_key, b_key = kk,
                sign_a = 1L, sign_b = 1L, output_key = tmp)
        output_key_old <- output_key
        output_key <- tmp
        if (output_key_old != keys[1L]) invisible(NULL)
      }
    }
    output_key
  }
  .sum_scalar <- function(key) {
    s <- list()
    for (srv in server_list) {
      ci <- which(server_names == srv)
      s[[srv]] <- .unwrap(.dsAgg(datasources[ci],
        call(name = "k2BeaverSumShareDS",
             source_key = key, session_id = session_id,
             frac_bits = 50L, ring = "ring127")))
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = s[[y_server]]$sum_share_fp,
      share_b = s[[nl]]$sum_share_fp,
      frac_bits = 50L, ring = "ring127"))
    as.numeric(agg$values[1L])
  }
  .grad_matvec <- function(residual_key, tag) {
    for (srv in server_list) {
      ci <- which(server_names == srv)
      .dsAgg(datasources[ci], call(name = "dsvertPrepareMultinomGradDS",
        residual_key = residual_key,
        is_outcome_server = (srv == y_server),
        n = as.integer(n_obs), session_id = session_id))
    }
    p_shared <- as.integer(sum(vapply(x_vars_per_server, length, integer(1L))))
    gt <- .unwrap(.dsAgg(datasources[dealer_ci],
      call(name = "glmRing63GenGradTriplesDS",
           dcf0_pk = .to_b64url(transport_pks[[y_server]]),
           dcf1_pk = .to_b64url(transport_pks[[nl]]),
           n = as.integer(n_obs), p = p_shared,
           ring = 127L, session_id = session_id)))
    gtk <- paste0("ord_strict_grad_", tag)
    .sendBlob(gt$grad_blob_0, gtk, ci_os)
    .sendBlob(gt$grad_blob_1, gtk, ci_nl)
    r1 <- list()
    for (srv in server_list) {
      ci <- which(server_names == srv)
      peer <- setdiff(server_list, srv)
      .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS",
        session_id = session_id, grad_triple_key = gtk))
      rr <- .unwrap(.dsAgg(datasources[ci], call(name = "k2GradientR1DS",
        peer_pk = .to_b64url(transport_pks[[peer]]), session_id = session_id)))
      r1[[srv]] <- rr
    }
    .sendBlob(r1[[y_server]]$encrypted_r1, "k2_grad_peer_r1", ci_nl)
    .sendBlob(r1[[nl]]$encrypted_r1, "k2_grad_peer_r1", ci_os)
    r2 <- list()
    for (srv in server_list) {
      ci <- which(server_names == srv)
      rr <- .unwrap(.dsAgg(datasources[ci], call(name = "k2GradientR2DS",
        party_id = if (srv == y_server) 0L else 1L,
        session_id = session_id)))
      r2[[srv]] <- rr
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r2[[y_server]]$gradient_fp,
      share_b = r2[[nl]]$gradient_fp,
      frac_bits = 50L, ring = "ring127"))
    as.numeric(agg$values) / n_obs
  }

  p_beta <- length(beta0)
  q_from_par <- function(beta, theta) {
    c(as.numeric(beta), theta[1L], log(pmax(diff(theta), 1e-4)))
  }
  par_from_q <- function(q) {
    beta <- q[seq_len(p_beta)]
    names(beta) <- names(beta0)
    a1 <- q[p_beta + 1L]
    if (K_minus_1 > 1L) {
      deltas <- exp(q[p_beta + 1L + seq_len(K_minus_1 - 1L)])
      theta <- c(a1, a1 + cumsum(deltas))
    } else {
      theta <- a1
    }
    names(theta) <- thresh_levels
    list(beta = beta, theta = theta)
  }
  theta_grad_to_q <- function(g_theta, theta) {
    out <- numeric(K_minus_1)
    out[1L] <- sum(g_theta)
    if (K_minus_1 > 1L) {
      d <- diff(theta)
      for (jj in seq_len(K_minus_1 - 1L)) {
        m <- jj + 1L
        out[m] <- d[jj] * sum(g_theta[m:K_minus_1])
      }
    }
    out
  }

  eval_counter <- 0L
  cache <- new.env(parent = emptyenv())
  secure_eval <- function(q, compute_loglik = FALSE) {
    key <- paste(c(signif(q, 12), as.integer(compute_loglik)), collapse = "|")
    if (!is.null(cache$key) && identical(cache$key, key)) return(cache$value)
    eval_counter <<- eval_counter + 1L
    tag <- paste0("eval", eval_counter)
    pt <- par_from_q(q)
    beta <- pt$beta
    theta <- pt$theta
    beta_slice <- function(srv) {
      vars <- x_vars_per_server[[srv]]
      if (is.null(vars) || length(vars) == 0L) NULL
      else unname(as.numeric(beta[vars]))
    }
    eta_key <- paste0("ord_strict_eta_", tag)
    for (srv in server_list) {
      ci <- which(server_names == srv)
      is_coord <- (srv == y_server)
      if (is_coord) {
        beta_coord <- beta_slice(y_server)
        beta_nl <- c(beta_slice(nl))
        for (extra_srv in non_dcf_servers) {
          beta_nl <- c(beta_nl, beta_slice(extra_srv))
        }
      } else {
        beta_coord <- beta_slice(y_server)
        beta_nl <- c()
        for (extra_srv in non_dcf_servers) {
          beta_nl <- c(beta_nl, beta_slice(extra_srv))
        }
        beta_nl <- c(beta_nl, beta_slice(nl))
      }
      .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
        beta_coord = beta_coord,
        beta_nl = beta_nl,
        intercept = 0,
        is_coordinator = is_coord,
        session_id = session_id,
        output_key = eta_key))
      if (!is_coord && p_extras > 0L) {
        .dsAgg(datasources[ci], call(name = "glmRing63ReorderXFullDS",
          p_coord = as.integer(p_coord),
          p_fusion = as.integer(p_fusion),
          p_extras = as.integer(p_extras),
          session_id = session_id))
      }
    }

    F_keys <- f_keys <- character(K_minus_1)
    for (kk in seq_len(K_minus_1)) {
      uk <- sprintf("ord_strict_u_%s_%d", tag, kk)
      expk <- sprintf("ord_strict_exp_%s_%d", tag, kk)
      denom <- sprintf("ord_strict_den_%s_%d", tag, kk)
      Fk <- sprintf("ord_strict_F_%s_%d", tag, kk)
      one_minus_F <- sprintf("ord_strict_omF_%s_%d", tag, kk)
      fk <- sprintf("ord_strict_f_%s_%d", tag, kk)
      .affine(a_key = eta_key, sign_a = 1L,
              const_fp = .fp_const(-theta[kk]), output_key = uk)
      .ring127_exp_round_keyed_extended(uk, expk, n_obs,
        datasources, dealer_ci, server_list, server_names,
        y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
      .affine(a_key = expk, sign_a = 1L,
              const_fp = one_fp, output_key = denom)
      .ring127_recip_round_keyed(denom, Fk, n_obs,
        datasources, dealer_ci, server_list, server_names,
        y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
      .affine(a_key = Fk, sign_a = -1L,
              const_fp = one_fp, output_key = one_minus_F)
      .ring127_vecmul(Fk, one_minus_F, fk, n_obs,
        datasources, dealer_ci, server_list, server_names,
        y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
      F_keys[kk] <- Fk
      f_keys[kk] <- fk
    }

    P_keys <- N_keys <- character(K)
    P_keys[1L] <- sprintf("ord_strict_P_%s_1", tag)
    .affine(a_key = F_keys[1L], sign_a = 1L, output_key = P_keys[1L])
    N_keys[1L] <- sprintf("ord_strict_N_%s_1", tag)
    .affine(a_key = f_keys[1L], sign_a = -1L, output_key = N_keys[1L])
    if (K_minus_1 >= 2L) {
      for (jj in 2L:K_minus_1) {
        P_keys[jj] <- sprintf("ord_strict_P_%s_%d", tag, jj)
        .affine(a_key = F_keys[jj], b_key = F_keys[jj - 1L],
                sign_a = 1L, sign_b = -1L, output_key = P_keys[jj])
        N_keys[jj] <- sprintf("ord_strict_N_%s_%d", tag, jj)
        .affine(a_key = f_keys[jj - 1L], b_key = f_keys[jj],
                sign_a = 1L, sign_b = -1L, output_key = N_keys[jj])
      }
    }
    P_keys[K] <- sprintf("ord_strict_P_%s_%d", tag, K)
    .affine(a_key = F_keys[K_minus_1], sign_a = -1L,
            const_fp = one_fp, output_key = P_keys[K])
    N_keys[K] <- sprintf("ord_strict_N_%s_%d", tag, K)
    .affine(a_key = f_keys[K_minus_1], sign_a = 1L, output_key = N_keys[K])

    prodP <- prodN <- character(K)
    for (jj in seq_len(K)) {
      prodP[jj] <- sprintf("ord_strict_prodP_%s_%d", tag, jj)
      prodN[jj] <- sprintf("ord_strict_prodN_%s_%d", tag, jj)
      .ring127_vecmul(mask_keys[jj], P_keys[jj], prodP[jj], n_obs,
        datasources, dealer_ci, server_list, server_names,
        y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
      .ring127_vecmul(mask_keys[jj], N_keys[jj], prodN[jj], n_obs,
        datasources, dealer_ci, server_list, server_names,
        y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
    }
    Py_key <- .sum_keys(prodP, sprintf("ord_strict_Py_%s", tag))
    Ny_key <- .sum_keys(prodN, sprintf("ord_strict_Ny_%s", tag))
    invPy_key <- sprintf("ord_strict_invPy_%s", tag)
    T_key <- sprintf("ord_strict_T_%s", tag)
    .ring127_recip_round_keyed(Py_key, invPy_key, n_obs,
      datasources, dealer_ci, server_list, server_names,
      y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
    .ring127_vecmul(Ny_key, invPy_key, T_key, n_obs,
      datasources, dealer_ci, server_list, server_names,
      y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)

    g_beta_part <- .grad_matvec(T_key, tag)
    server_partition_names <- unlist(x_vars_per_server, use.names = FALSE)
    g_beta <- numeric(p_beta)
    names(g_beta) <- names(beta0)
    perm <- match(names(beta0), server_partition_names)
    ok <- !is.na(perm)
    g_beta[ok] <- g_beta_part[perm[ok]]

    g_theta <- numeric(K_minus_1)
    for (kk in seq_len(K_minus_1)) {
      fk_inv <- sprintf("ord_strict_f_inv_%s_%d", tag, kk)
      mdiff <- sprintf("ord_strict_mdiff_%s_%d", tag, kk)
      term <- sprintf("ord_strict_theta_term_%s_%d", tag, kk)
      .ring127_vecmul(f_keys[kk], invPy_key, fk_inv, n_obs,
        datasources, dealer_ci, server_list, server_names,
        y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
      .affine(a_key = mask_keys[kk], b_key = mask_keys[kk + 1L],
              sign_a = 1L, sign_b = -1L, output_key = mdiff)
      .ring127_vecmul(fk_inv, mdiff, term, n_obs,
        datasources, dealer_ci, server_list, server_names,
        y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
      g_theta[kk] <- .sum_scalar(term) / n_obs
    }
    names(g_theta) <- thresh_levels

    loglik <- NA_real_
    if (isTRUE(compute_loglik)) {
      scaledP <- sprintf("ord_strict_scaledP_%s", tag)
      logP <- sprintf("ord_strict_logP_%s", tag)
      logPc <- sprintf("ord_strict_logPc_%s", tag)
      for (srv in server_list) {
        ci <- which(server_names == srv)
        .dsAgg(datasources[ci], call(name = "k2Ring127LocalScaleDS",
          in_key = Py_key, scalar_fp = scale_log_fp,
          output_key = scaledP, n = as.numeric(n_obs),
          session_id = session_id,
          is_party0 = (srv == y_server)))
      }
      .ring127_log_round_keyed_nr(scaledP, logP, n_obs,
        datasources, dealer_ci, server_list, server_names,
        y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob,
        nr_iters = 4L)
      .affine(a_key = logP, sign_a = 1L,
              const_fp = neg_log_scale_fp, output_key = logPc)
      loglik <- .sum_scalar(logPc) / n_obs
    }

    grad_q <- c(as.numeric(g_beta), theta_grad_to_q(g_theta, theta))
    value <- list(loglik = loglik, grad_q = grad_q,
                  beta = beta, theta = theta,
                  grad_beta = g_beta, grad_theta = g_theta)
    cache$key <- key
    cache$value <- value
    if (verbose) {
      ll_msg <- if (is.finite(loglik)) sprintf(" avg_loglik=%.6f", loglik) else ""
      message(sprintf("[OrdinalJointStrict eval %d]%s |g|_inf=%.3e",
                      eval_counter, ll_msg, max(abs(grad_q))))
    }
    value
  }

  q0 <- q_from_par(beta0, theta0)
  q <- q0
  opt_trace <- list()
  converged <- FALSE
  Hinv <- diag(length(q0))
  if (!is.null(warm$joint_mle$covariance) &&
      nrow(warm$joint_mle$covariance) == length(q0)) {
    d0 <- abs(diag(warm$joint_mle$covariance)) * n_obs
    d0[!is.finite(d0) | d0 <= 0] <- 1
    Hinv <- diag(pmin(pmax(d0, 1e-3), 10), length(q0))
  }
  Hinv0 <- Hinv
  cur <- secure_eval(q, compute_loglik = FALSE)
  best <- cur
  best_q <- q
  best_norm <- max(abs(cur$grad_q))

  update_best <- function(val, q_val) {
    nrm <- max(abs(val$grad_q))
    if (is.finite(nrm) && nrm < best_norm) {
      best <<- val
      best_q <<- q_val
      best_norm <<- nrm
    }
    nrm
  }
  try_score_step <- function(q_cur, cur_val, step, cap, max_halves = 6L) {
    if (any(!is.finite(step))) return(list(accepted = FALSE))
    st <- max(abs(step))
    if (is.finite(st) && st > cap) step <- step * (cap / st)
    g_norm <- max(abs(cur_val$grad_q))
    for (half in 0:max_halves) {
      alpha <- 0.5^half
      q_try <- q_cur + alpha * step
      trial <- secure_eval(q_try, compute_loglik = FALSE)
      trial_norm <- update_best(trial, q_try)
      if (is.finite(trial_norm) &&
          trial_norm <= g_norm * (1 + 1e-3)) {
        return(list(accepted = TRUE, q = q_try, cur = trial,
                    norm = trial_norm))
      }
    }
    list(accepted = FALSE)
  }
  fd_newton_used <- FALSE
  fd_iters_done <- 0L
  use_fd_newton <- isTRUE(getOption("dsvert.ord_strict_fd_newton", TRUE)) &&
    length(q0) <= as.integer(getOption("dsvert.ord_strict_fd_max_dim", 8L))
  if (use_fd_newton) {
    fd_newton_used <- TRUE
    fd_iters <- min(as.integer(max_outer),
                    as.integer(getOption("dsvert.ord_strict_fd_iters", 2L)))
    fd_eps <- as.numeric(getOption("dsvert.ord_strict_fd_eps", 1e-2))
    fd_cap <- as.numeric(getOption("dsvert.ord_strict_fd_step_cap", 0.25))
    fd_ridge <- as.numeric(getOption("dsvert.ord_strict_fd_ridge", 1e-6))
    for (fd_iter in seq_len(fd_iters)) {
      fd_iters_done <- fd_iter
      g <- cur$grad_q
      g_norm <- update_best(cur, q)
      opt_trace[[length(opt_trace) + 1L]] <-
        list(grad_inf = g_norm, q = q, optimizer = "fd_newton")
      if (is.finite(g_norm) && g_norm < tol) {
        converged <- TRUE
        break
      }
      p_dim <- length(q)
      H <- matrix(0, p_dim, p_dim)
      for (jj in seq_len(p_dim)) {
        h <- fd_eps * max(1, abs(q[jj]))
        qp <- q
        qm <- q
        qp[jj] <- qp[jj] + h
        qm[jj] <- qm[jj] - h
        ep <- secure_eval(qp, compute_loglik = FALSE)
        em <- secure_eval(qm, compute_loglik = FALSE)
        H[, jj] <- (ep$grad_q - em$grad_q) / (2 * h)
      }
      H <- (H + t(H)) / 2
      step <- tryCatch({
        ee <- eigen(H, symmetric = TRUE)
        vals <- -pmax(abs(ee$values), fd_ridge)
        -as.numeric(ee$vectors %*%
          (as.numeric(crossprod(ee$vectors, g)) / vals))
      }, error = function(e) {
        tryCatch(-as.numeric(solve(H - fd_ridge * diag(p_dim), g)),
                 error = function(e2) rep(NA_real_, p_dim))
      })
      fd_trial <- try_score_step(q, cur, step, fd_cap, max_halves = 8L)
      if (!isTRUE(fd_trial$accepted)) {
        q <- best_q
        cur <- best
        break
      }
      q <- fd_trial$q
      cur <- fd_trial$cur
      if (is.finite(fd_trial$norm) && fd_trial$norm < tol) {
        converged <- TRUE
        break
      }
    }
    Hinv <- Hinv0
  }

  remaining_outer <- max(0L, as.integer(max_outer) - length(opt_trace))
  for (iter in seq_len(remaining_outer)) {
    g <- cur$grad_q
    g_norm <- max(abs(g))
    opt_trace[[length(opt_trace) + 1L]] <-
      list(grad_inf = g_norm, q = q, optimizer = "bfgs_score")
    update_best(cur, q)
    if (is.finite(g_norm) && g_norm < tol) {
      converged <- TRUE
      break
    }
    step <- as.numeric(Hinv %*% g)
    if (any(!is.finite(step))) step <- g
    cap <- getOption("dsvert.ord_strict_step_cap", 0.01)
    st <- max(abs(step))
    if (is.finite(st) && st > cap) step <- step * (cap / st)
    accepted <- FALSE
    trial <- NULL
    for (half in 0:5) {
      alpha <- 0.5^half
      q_try <- q + alpha * step
      trial <- secure_eval(q_try, compute_loglik = FALSE)
      trial_norm <- update_best(trial, q_try)
      if (is.finite(trial_norm) &&
          trial_norm <= g_norm * (1 + 1e-3)) {
        accepted <- TRUE
        q_new <- q_try
        break
      }
    }
    if (!accepted || is.null(trial)) {
      q <- best_q
      cur <- best
      break
    }
    s <- q_new - q
    y <- trial$grad_q - g
    sy <- sum(s * y)
    # For maximisation, the local Hessian of the score is negative
    # definite. Use a damped positive BFGS update on -score curvature only
    # when the secant has usable orientation; otherwise keep the diagonal
    # warm-start metric.
    if (is.finite(sy) && sy < -1e-10) {
      y_pos <- -y
      sy_pos <- sum(s * y_pos)
      Hy <- Hinv %*% y_pos
      denom <- as.numeric(t(y_pos) %*% Hy)
      if (is.finite(sy_pos) && sy_pos > 1e-10 &&
          is.finite(denom) && denom > 1e-10) {
        Hinv <- Hinv + (s %*% t(s)) / sy_pos -
          (Hy %*% t(Hy)) / denom
      }
    }
    q <- q_new
    cur <- trial
  }
  q <- best_q
  final <- best
  if (isTRUE(getOption("dsvert.ord_strict_loglik", FALSE))) {
    final <- secure_eval(q, compute_loglik = TRUE)
  }

  for (srv in server_list) {
    ci <- which(server_names == srv)
    try(.dsAgg(datasources[ci],
               call(name = "mpcCleanupDS", session_id = session_id)),
        silent = TRUE)
  }

  out <- warm
  out$beta_po_joint <- final$beta
  out$thresholds_joint <- final$theta
  out$strict_non_disclosive <- TRUE
  out$strict_optimizer <- list(trace = opt_trace, evals = eval_counter,
                               converged = converged,
                               fd_newton = fd_newton_used,
                               fd_iters = fd_iters_done)
  out$strict_avg_loglik <- final$loglik
  out$strict_grad_beta <- final$grad_beta
  out$strict_grad_theta <- final$grad_theta
  out$outer_iter <- length(opt_trace)
  out$converged <- converged
  out$family <- "ordinal_joint_po_ring127_strict"
  out$session_id <- session_id
  class(out) <- c("ds.vertOrdinalJointNewton", class(out))
  out
}


#' @title Federated joint proportional-odds ordinal regression via
#'   Ring127 MPC-orchestrated Newton iteration
#' @description Strict K=2 share-domain proportional-odds Newton route by
#'   default. Patient-level linear predictors, class probabilities,
#'   reciprocals, and score terms remain Ring127 shares; the client receives
#'   only guarded class counts, aggregate score/Hessian probes, and final
#'   model parameters. The older per-patient reconstruction route has been
#'   removed from the product package.
#'
#'   For a K-level ordered outcome
#'   the PO log-likelihood is
#'     \deqn{\ell(\beta, \theta) = \sum_i \log[F(\theta_{y_i} - \eta_i)
#'                                            - F(\theta_{y_i-1} - \eta_i)],}
#'   with \eqn{F = \mathrm{sigmoid}}, \eqn{\eta_i = X_i \beta},
#'   \eqn{\theta_0 = -\infty}, \eqn{\theta_K = +\infty}. Score:
#'     \deqn{\partial \ell / \partial \beta_j = -\sum_i x_{ij} \cdot
#'           \frac{f(\theta_{y_i}-\eta_i) - f(\theta_{y_i-1}-\eta_i)}
#'                {F(\theta_{y_i}-\eta_i) - F(\theta_{y_i-1}-\eta_i)},}
#'   with \eqn{f(u) = F(u)(1-F(u))}.
#'
#'   MPC pipeline per outer Newton iter:
#'   \enumerate{
#'     \item Compute \eqn{\eta} share via \code{k2ComputeEtaShareDS}.
#'     \item For each threshold \eqn{k}: compute \eqn{\theta_k - \eta}
#'           share via affine-combine (plaintext \eqn{\theta_k} public).
#'     \item Apply exp + recip on \eqn{\theta_k - \eta} -> \eqn{F_k} share
#'           (\eqn{F(u) = 1/(1+\exp(-u)) = e^u/(1+e^u)} -- evaluated as
#'           \code{exp(u) * (1/(1+exp(u)))} via existing primitives).
#'     \item \eqn{f_k = F_k (1 - F_k)} via Beaver vecmul.
#'     \item Per-patient residual numerator/denominator built from
#'           indicator-weighted differences (done on outcome server
#'           since it holds \eqn{y_i} plaintext).
#'     \item \code{.ring127_recip_round_keyed} on the F-difference share.
#'     \item Beaver vecmul (f-diff) * (1/F-diff) -> \eqn{T_i} share.
#'     \item Beaver matvec \eqn{X^\top T} -> aggregate score for \eqn{\beta}.
#'     \item Client Newton on stacked \eqn{(\beta, \theta)} using an
#'           aggregate finite-difference Hessian over share-domain scores.
#'   }
#'
#' @param formula Ordered outcome on LHS.
#' @param data Aligned data name.
#' @param levels_ordered Character vector of ordered levels (low -> high).
#' @param cumulative_template e.g. \code{"\%s_leq"} for Y <= k indicator.
#' @param max_outer Outer Newton iterations.
#' @param tol Convergence tolerance on \eqn{\|\Delta (\beta, \theta)\|_\infty}.
#' @param verbose Logical.
#' @param datasources DataSHIELD connections.
#' @return \code{ds.vertOrdinalJointNewton} object.
#' @export
ds.vertOrdinalJointNewton <- function(formula, data = NULL, levels_ordered,
                                      cumulative_template = "%s_leq",
                                      max_outer = 8L, tol = 1e-4,
                                      verbose = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  .ord_joint_secure_fit(
    formula = formula, data = data,
    levels_ordered = levels_ordered,
    cumulative_template = cumulative_template,
    max_outer = max_outer, tol = tol,
    verbose = verbose, datasources = datasources)
}

#' @export
print.ds.vertOrdinalJointNewton <- function(x, ...) {
  cat("dsVert joint PO ordinal (Ring127 Newton)\n")
  cat(sprintf("  N = %d  levels = %s  outer_iter = %d  converged = %s\n",
              x$n_obs, paste(x$levels, collapse = ","),
              x$outer_iter, x$converged))
  invisible(x)
}
