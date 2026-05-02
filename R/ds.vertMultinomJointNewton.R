#' @title Federated joint-softmax multinomial logistic regression via
#'   Ring127 MPC-orchestrated Newton iteration
#' @description Full per-patient softmax Newton: orchestrates K-1 parallel
#'   exp(eta_k) shares, sums to denominator D = 1+Sumexp(eta_k), computes 1/D
#'   via Ring127 Chebyshev + Newton-Raphson, multiplies per class to get
#'   p_k(x_i) share per patient, builds residual y_k - p_k on outcome
#'   server, and aggregates X^T(y_k - p_k) via existing Beaver matvec
#'   pipeline for each class. Client-side Bohning-Hessian-bounded Newton
#'   step on stacked beta.
#'
#'   All per-patient quantities stay as Ring127 additive shares; only the
#'   final p(K-1)-dim aggregate gradient per iter is revealed -- same
#'   privacy class as the single-class gradient of ds.vertGLM. **P3 delta:
#'   zero new reveal types.**
#'
#' @param formula R formula with categorical outcome on LHS.
#' @param data Aligned data frame name.
#' @param levels Character vector of outcome levels (first = reference).
#' @param indicator_template sprintf template for class-indicator columns
#'   on the outcome server, e.g. \code{"\%s_ind"}. Must exist server-side.
#' @param max_outer Outer Newton iterations (default 8).
#' @param tol Convergence tolerance on max |Deltabeta| (default 1e-4).
#' @param verbose Logical.
#' @param datasources DataSHIELD connections.
#' @return \code{ds.vertMultinomJointNewton} object.
#' @export
ds.vertMultinomJointNewton <- function(formula, data = NULL, levels,
                                        indicator_template = "%s_ind",
                                        max_outer = 8L, tol = 1e-4,
                                        verbose = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (!inherits(formula, "formula"))
    stop("formula required", call. = FALSE)
  if (!is.character(levels) || length(levels) < 3L)
    stop("levels must have >= 3 ordered levels (first = reference)",
         call. = FALSE)
  ref <- levels[1L]
  non_ref <- levels[-1L]
  K_minus_1 <- length(non_ref)

  # Warm start via ds.vertMultinom (OVR + softmax-anchor).
  warm <- ds.vertMultinom(formula, data = data, classes = levels,
                           reference = ref,
                           indicator_template = indicator_template,
                           verbose = verbose, datasources = datasources)
  beta_mat <- warm$coefficients
  if (is.null(beta_mat))
    stop("warm start produced no coefficients", call. = FALSE)
  cnames <- rownames(beta_mat)
  p <- length(cnames)

  # Setup a fresh Ring127 session compatible with the joint pipeline.
  # We reuse the internal `.glm_mpc_setup` helper from ds.vertGLM.setup
  # to get {y_server, nl, server_list, transport_pks, session_id}.
  y_var_char <- .ds_gee_extract_lhs(formula)
  server_names <- names(datasources)
  y_server <- .ds_gee_find_server_holding(datasources, server_names,
                                           data, y_var_char)
  if (is.null(y_server)) stop("outcome server for ", y_var_char,
                               " not found", call. = FALSE)
  nl <- setdiff(server_names, y_server)[1L]
  if (is.na(nl) || is.null(nl))
    stop("non-label server not found", call. = FALSE)
  server_list <- c(y_server, nl)
  y_server_ci <- which(server_names == y_server)
  dealer_ci <- which(server_names == nl)

  session_id <- paste0("multinomJoint_", as.integer(Sys.time()),
                       "_", sample.int(.Machine$integer.max, 1L))
  transport_pks <- list()
  for (srv in server_list) {
    ci <- which(server_names == srv)
    r <- DSI::datashield.aggregate(datasources[ci],
      call(name = "glmRing63TransportInitDS", session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    transport_pks[[srv]] <- r$transport_pk
  }
  .json_to_b64url <- function(x) {
    raw <- charToRaw(jsonlite::toJSON(x, auto_unbox = TRUE))
    b64 <- gsub("\n", "", jsonlite::base64_enc(raw), fixed = TRUE)
    chartr("+/", "-_", sub("=+$", "", b64, perl = TRUE))
  }
  # base64url for already-base64 strings (peer transport PKs). Required
  # to avoid Opal DSL parser eating "=" / "+" / "/" in the call expr.
  .to_b64url <- function(x) chartr("+/", "-_", sub("=+$", "", x, perl = TRUE))
  for (srv in server_list) {
    ci <- which(server_names == srv)
    peer_srv <- setdiff(server_list, srv)
    peers <- setNames(list(transport_pks[[peer_srv]]), peer_srv)
    DSI::datashield.aggregate(datasources[ci],
      call(name = "mpcStoreTransportKeysDS",
           transport_keys_b64 = .json_to_b64url(peers),
           session_id = session_id))
  }

  # PSI + share input (uses ds.vertCox-like pattern). For simplicity we
  # delegate input sharing to a quick pre-pass and REUSE warm$fits state
  # when the implementation is tested locally. For federated correctness
  # with K=2 servers, we invoke k2ShareInputDS once per server at
  # ring = 127L using the SAME X columns as warm.
  rhs <- attr(terms(formula), "term.labels")
  x_vars_y <- intersect(rhs, names(warm$fits[[1L]]$x_means))
  # Partition by which server holds each feature (query each server's
  # columns)
  x_vars_per_server <- list()
  for (srv in server_list) {
    ci <- which(server_names == srv)
    r <- tryCatch(DSI::datashield.aggregate(datasources[ci],
      call(name = "dsvertColNamesDS", data_name = data)),
      error = function(e) NULL)
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    cols_here <- if (is.list(r) && !is.null(r$columns)) r$columns
                 else if (is.character(r)) r else character(0)
    x_vars_per_server[[srv]] <- intersect(rhs, cols_here)
  }

  .dsAgg <- function(conns, expr, ...)
    DSI::datashield.aggregate(conns, expr, ...)
  # Adaptive chunking (same pattern as ds.psiAlign). Canonical server
  # signature: mpcStoreBlobDS(key, chunk, chunk_index, n_chunks, session_id).
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

  # Share input at Ring127
  share_results <- list()
  for (srv in server_list) {
    ci <- which(server_names == srv)
    peer <- setdiff(server_list, srv)
    peer_pk <- .to_b64url(transport_pks[[peer]])
    # Categorical outcome (bp_cls) is not shared via k2ShareInputDS
    # (fails Go float64 unmarshal). Per-class indicator columns are
    # read directly by the label server inside later aggregates (see
    # indicator_template usage at ~line 240). Pass y_var=NULL for all
    # servers in the multinom joint path.
    r <- .dsAgg(datasources[ci], call(name = "k2ShareInputDS",
      data_name = data, x_vars = x_vars_per_server[[srv]],
      y_var = NULL,
      peer_pk = peer_pk, ring = 127L, session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    share_results[[srv]] <- r
  }
  for (srv in server_list) {
    peer <- setdiff(server_list, srv)
    peer_ci <- which(server_names == peer)
    .sendBlob(share_results[[srv]]$encrypted_x_share, "k2_peer_x_share", peer_ci)
    if (!is.null(share_results[[srv]]$encrypted_y_share))
      .sendBlob(share_results[[srv]]$encrypted_y_share, "k2_peer_y_share", peer_ci)
  }
  for (srv in server_list) {
    ci <- which(server_names == srv)
    peer <- setdiff(server_list, srv)
    .dsAgg(datasources[ci], call(name = "k2ReceiveShareDS",
      peer_p = as.integer(length(x_vars_per_server[[peer]])),
      session_id = session_id))
  }
  n_obs <- share_results[[y_server]]$n
  if (verbose) message(sprintf("[MultinomJointNewton] session %s  n=%d  K=%d  p=%d",
                                session_id, n_obs, K_minus_1 + 1L, p))

  converged <- FALSE
  final_iter <- max_outer
  # Best-beta tracking: retain argmin_k |g(beta_k)|_L2. Required because the
  # MPC-approximated gradient + step-cap interaction induces late-iter
  # oscillation (empirically iter 7 and iter 10 spikes on NHANES
  # Mnl trace 20260423-112816). Returning best-so-far mirrors Cox Path
  # B's revert-on-grad-grow (memory project_ring127_p1_progress).
  best_beta <- NULL
  best_g_norm <- Inf
  best_iter <- 0L
  best_step_beta <- NULL
  best_step_norm <- Inf
  best_step_iter <- 0L

  # Bohning (1992) constant upper-bound Hessian: H* = (1/2)*(I_{K-1} -
  # (1/K) 1 1^T) (x) (X^T X / n). PSD + beta-independent -> monotone Newton
  # descent (Bohning 1992 Ann Inst Stat Math 44:197-200, Theorem 2;
  # Krishnapuram et al 2005 IEEE PAMI 27(6)).
  #
  # XtX_over_n must be in FORMULA order to match beta_mat rows. The
  # warm fit's $covariance is in fit-internal (server-partition) order
  # per the LASSO permutation bug (same class as 4ce55a3 -- see
  # dsVertGLM.k2.R theta_conv layout). Reconstruct via hessian_std +
  # x_means + x_sds from the warm fit (those fields ARE in formula
  # order), applying the same Gram-from-hessian formula as LASSO:
  #   G[j,k] = x_j x_k + x_sd_j x_sd_k * H_std[perm(j), perm(k)]  (slopes)
  #   G[0,j] = x_j  G[0,0] = 1
  w0 <- warm$fits[[non_ref[1L]]]
  XtX_over_n <- NULL
  if (!is.null(w0$hessian_std) && !is.null(w0$x_means) &&
      !is.null(w0$x_sds) && is.matrix(w0$hessian_std) &&
      all(dim(w0$hessian_std) == c(p, p))) {
    lam_ridge <- if (!is.null(w0$lambda) && is.finite(w0$lambda)) w0$lambda else 0
    H_std <- w0$hessian_std - lam_ridge * diag(p)
    if (!is.null(dimnames(H_std)) && !is.null(rownames(H_std))) {
      perm <- match(cnames, rownames(H_std))
      if (all(!is.na(perm))) H_std <- H_std[perm, perm, drop = FALSE]
    }
    x_m <- as.numeric(w0$x_means[cnames]); x_m[is.na(x_m)] <- 0
    x_s <- as.numeric(w0$x_sds[cnames]);    x_s[is.na(x_s)] <- 0
    XtX_over_n <- matrix(0, p, p, dimnames = list(cnames, cnames))
    int_j <- which(cnames == "(Intercept)")
    if (length(int_j) != 1L) int_j <- NA_integer_
    for (jj in seq_len(p)) for (kk in seq_len(p)) {
      if (!is.na(int_j) && jj == int_j && kk == int_j) XtX_over_n[jj,kk] <- 1
      else if (!is.na(int_j) && jj == int_j) XtX_over_n[jj,kk] <- x_m[kk]
      else if (!is.na(int_j) && kk == int_j) XtX_over_n[jj,kk] <- x_m[jj]
      else XtX_over_n[jj,kk] <- x_m[jj]*x_m[kk] + x_s[jj]*x_s[kk]*H_std[jj,kk]
    }
  } else {
    # Fallback for external warm starts without hessian_std -- use covariance.
    cov_k0 <- w0$covariance
    sigma2 <- if (!is.null(w0$deviance))
      w0$deviance / max(n_obs - p, 1L) else 1
    XtX_over_n <- tryCatch(sigma2 * solve(cov_k0) / n_obs,
                           error = function(e) diag(p))
  }
  B <- matrix(0, p * K_minus_1, p * K_minus_1)
  for (i in seq_len(K_minus_1))
    for (j in seq_len(K_minus_1)) {
      row_i <- ((i-1L)*p + 1L):(i*p)
      col_j <- ((j-1L)*p + 1L):(j*p)
      coef_ij <- if (i == j) 0.5 * (1 - 1/(K_minus_1 + 1L))
                 else -0.5 / (K_minus_1 + 1L)
      B[row_i, col_j] <- coef_ij * XtX_over_n
    }
  B_reg <- B + 1e-6 * max(abs(diag(B))) * diag(p * K_minus_1)

  coord <- y_server   # label server
  for (outer in seq_len(max_outer)) {
    t_iter <- proc.time()[[3L]]
    # Step 1-2: for each class, compute eta_k share + exp(eta_k) share
    exp_eta_keys <- character(K_minus_1)
    for (ki in seq_along(non_ref)) {
      k <- non_ref[ki]
      b_k <- beta_mat[, k]
      b_coord_vec <- b_k[intersect(names(b_k), x_vars_per_server[[coord]])]
      b_nl_vec    <- b_k[intersect(names(b_k), x_vars_per_server[[nl]])]
      intercept_k <- unname(b_k["(Intercept)"])
      if (is.na(intercept_k)) intercept_k <- 0
      eta_key <- paste0("eta_class_", ki)
      for (srv in server_list) {
        ci <- which(server_names == srv)
        .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
          beta_coord = b_coord_vec, beta_nl = b_nl_vec,
          intercept = if (srv == coord) intercept_k else 0,
          is_coordinator = (srv == coord),
          session_id = session_id, output_key = eta_key))
      }
      exp_key <- paste0("exp_eta_class_", ki)
      .ring127_exp_round_keyed_extended(eta_key, exp_key, n_obs,
                               datasources, dealer_ci, server_list,
                               server_names, y_server, nl, transport_pks,
                               session_id, .dsAgg, .sendBlob)
      exp_eta_keys[ki] <- exp_key
    }

    # Step 3: D = 1 + sum(exp(eta_k))
    for (srv in server_list) {
      ci <- which(server_names == srv)
      .dsAgg(datasources[ci], call(name = "dsvertSoftmaxDenominatorDS",
        exp_eta_keys = exp_eta_keys, output_key = "D_share",
        is_party0 = (srv == coord), n = as.integer(n_obs),
        session_id = session_id))
    }

    # Step 4: 1/D via Chebyshev + NR
    .ring127_recip_round_keyed("D_share", "inv_D", n_obs,
                               datasources, dealer_ci, server_list,
                               server_names, y_server, nl, transport_pks,
                               session_id, .dsAgg, .sendBlob)

    # Step 5: p_k = exp(eta_k) * inv_D per class
    p_keys <- character(K_minus_1)
    for (ki in seq_along(non_ref)) {
      p_keys[ki] <- paste0("p_class_", ki)
      .ring127_vecmul(exp_eta_keys[ki], "inv_D", p_keys[ki], n_obs,
                      datasources, dealer_ci, server_list, server_names,
                      y_server, nl, transport_pks, session_id,
                      .dsAgg, .sendBlob)
    }

    # Step 6-7: residual r_k = y_k - p_k + Beaver matvec X^T r_k.
    # Dimension handling: X is shared via k2ShareInputDS with p_shared =
    # sum(lengths(x_vars_per_server)) columns (slopes only, no intercept
    # column). The Beaver matvec therefore returns p_shared-length slope
    # gradient. The intercept gradient = Sum_i r_i is read off the SAME
    # Beaver round as a side-product of k2GradientR1DS (its
    # `sum_residual_fp` field; identical pattern to ds.vertGLM K=2
    # in dsVertGLM.k2.R:309). One Beaver round per (iter x class) with
    # two k2-ring63-aggregate reveals (slope vector + intercept scalar)
    # -- threat model isomorphic to GLM K=2 audit OK non-disclosive.
    # The previous separate k2BeaverSumShareDS round (which consumed an
    # extra Beaver triple per iter per class) was removed per reviewer
    # directive 2026-04-26. Cites: Bohning 1992 Ann Inst Stat Math
    # 44:197-200 Theorem 2 (constant majorant Hessian); Krishnapuram
    # -Carin-Figueiredo-Hartemink 2005 IEEE PAMI 27(6) (peer use of
    # Bohning H* for MPC-friendly multinomial logistic regression).
    p_shared <- sum(vapply(x_vars_per_server, length, integer(1L)))
    int_row <- which(cnames == "(Intercept)")
    slope_rows <- setdiff(seq_len(p), int_row)
    gradients <- matrix(0, p, K_minus_1,
                         dimnames = list(cnames, non_ref))
    for (ki in seq_along(non_ref)) {
      k <- non_ref[ki]
      r_key <- paste0("r_class_", ki)
      for (srv in server_list) {
        ci <- which(server_names == srv)
        .dsAgg(datasources[ci], call(name = "dsvertComputeResidualShareDS",
          p_key = p_keys[ki],
          indicator_col = if (srv == coord) sprintf(indicator_template, k) else NULL,
          data_name = if (srv == coord) data else NULL,
          output_key = r_key,
          is_outcome_server = (srv == coord),
          n = as.integer(n_obs), session_id = session_id))
      }
      # Prep the standard gradient machinery: pipeline computes X^T(mu-y),
      # we set mu=0 and y=-r_k so X^T(mu-y) = X^T r_k. As a side-product
      # k2GradientR1DS emits Sum(mu-y) = Sum_i r_i (per-server share) in
      # `sum_residual_fp`, which we will aggregate below for the
      # intercept gradient -- same pattern as dsVertGLM.k2.R:309.
      for (srv in server_list) {
        ci <- which(server_names == srv)
        .dsAgg(datasources[ci], call(name = "dsvertPrepareMultinomGradDS",
          residual_key = r_key,
          is_outcome_server = (srv == coord),
          n = as.integer(n_obs), session_id = session_id))
      }
      # Beaver matvec: X^T (mu - y) -> per-class SLOPE gradient (length p_shared)
      grad_t <- .dsAgg(datasources[dealer_ci],
        call(name = "glmRing63GenGradTriplesDS",
             dcf0_pk = .to_b64url(transport_pks[[coord]]),
             dcf1_pk = .to_b64url(transport_pks[[nl]]),
             n = as.integer(n_obs), p = as.integer(p_shared),
             ring = 127L, session_id = session_id))
      if (is.list(grad_t) && length(grad_t) == 1L) grad_t <- grad_t[[1L]]
      # Per-class blob-key namespacing (defensive: ABY3 Sec.IV.D pool
      # isolation; MP-SPDZ Multiplications.hpp). Eliminates cross-class
      # blob-key collision under the shared key "k2_grad_triple_fp" that
      # the K-1 classes within this outer iter (and the next iter) all
      # used to share. Hypothesised root cause of intermittent iter-2
      # NPE on s2 (3/9 approx 33% empirical rate, see paper Sec.VIII bullet #4).
      grad_triple_key <- sprintf("k2_grad_triple_fp_iter%d_class%d", outer, ki)
      .sendBlob(grad_t$grad_blob_0, grad_triple_key, y_server_ci)
      .sendBlob(grad_t$grad_blob_1, grad_triple_key, dealer_ci)
      r1 <- list()
      for (srv in server_list) {
        ci <- which(server_names == srv)
        peer <- setdiff(server_list, srv)
        .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS",
          session_id = session_id,
          grad_triple_key = grad_triple_key))
        rr <- .dsAgg(datasources[ci], call(name = "k2GradientR1DS",
          peer_pk = .to_b64url(transport_pks[[peer]]), session_id = session_id))
        if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
        r1[[srv]] <- rr
      }
      .sendBlob(r1[[coord]]$encrypted_r1, "k2_grad_peer_r1", dealer_ci)
      .sendBlob(r1[[nl]]$encrypted_r1, "k2_grad_peer_r1", y_server_ci)
      r2 <- list()
      for (srv in server_list) {
        ci <- which(server_names == srv)
        is_c <- (srv == coord)
        rr <- .dsAgg(datasources[ci], call(name = "k2GradientR2DS",
          party_id = if (is_c) 0L else 1L, session_id = session_id))
        if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
        r2[[srv]] <- rr
      }
      # Intercept gradient: Sum_i r_i is the side-product of the SAME
      # k2GradientR1DS round that produced the slope share (its
      # `sum_residual_fp` field). Aggregate the two server shares to
      # plaintext Sum r_i; divide by n for the intercept-row gradient.
      # Same disclosure pattern as dsVertGLM.k2.R:309 (audit OK).
      int_res_agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = r1[[coord]]$sum_residual_fp,
        share_b = r1[[nl]]$sum_residual_fp,
        frac_bits = 50L, ring = "ring127"))
      int_val <- as.numeric(int_res_agg$values)[1L]
      if (!is.finite(int_val)) {
        stop(sprintf("[mnl_joint] iter=%d class=%d: intercept grad NA from r1$sum_residual_fp; coord=%s nl=%s",
                     outer, ki,
                     substr(r1[[coord]]$sum_residual_fp %||% "NULL", 1, 12),
                     substr(r1[[nl]]$sum_residual_fp %||% "NULL", 1, 12)),
             call. = FALSE)
      }
      gradients[int_row, ki] <- int_val / n_obs

      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = r2[[coord]]$gradient_fp,
        share_b = r2[[nl]]$gradient_fp,
        frac_bits = 50L, ring = "ring127"))
      # Slope gradient X^T r_k (length p_shared); intercept already filled.
      slope_vals <- as.numeric(agg$values)
      if (length(slope_vals) != length(slope_rows)) {
        stop(sprintf("slope gradient length mismatch: got %d, expected %d",
                     length(slope_vals), length(slope_rows)), call. = FALSE)
      }
      if (any(!is.finite(slope_vals))) {
        stop(sprintf("[debug bug#9] iter=%d class=%d: slope grad NA at positions [%s]",
                     outer, ki, paste(which(!is.finite(slope_vals)), collapse=",")),
             call. = FALSE)
      }
      # AUDITORIA root cause of multinom Newton divergence (|g|_L2
      # INCREASING monotonically at iters 1->10 per grad-trace 9e1e8fc):
      # slope_vals comes back from the Beaver matvec in SERVER-PARTITION
      # order (concatenation of x_vars_per_server[[srv]] over
      # server_list). slope_rows indexes cnames in FORMULA order (from
      # rownames(beta_mat)). Writing slope_vals directly was pairing
      # gradient of variable A with beta of variable B. Exactly the same
      # class of bug as LASSO 4ce55a3 (federated hessian_std column
      # permutation). Permute to formula order before assigning.
      server_partition_names <- unlist(x_vars_per_server, use.names = FALSE)
      formula_slope_names <- cnames[slope_rows]
      perm_grad <- match(formula_slope_names, server_partition_names)
      if (any(is.na(perm_grad))) {
        stop(sprintf("slope gradient permutation mismatch: formula names [%s] not all present in server partition [%s]",
                     paste(formula_slope_names, collapse=","),
                     paste(server_partition_names, collapse=",")),
             call. = FALSE)
      }
      gradients[slope_rows, ki] <- slope_vals[perm_grad] / n_obs
    }

    # H_emp activation gate. Default FALSE per reviewer architectural
    # call 2026-04-29: Bohning B_reg ships as the production Hessian.
    #
    # Diagnostic history. The first L3 sweep of the H_emp pipeline at
    # SHA fe03cf1 produced max|Deltapi| = 1.19e-1 (regression vs the prior
    # Bohning baseline 4.86e-2) with iter-4 H_emp_full containing a
    # single sign-flipped diagonal entry of magnitude approx10^6 -- non-
    # physical for X^T diag(p(1-p)) X which is PSD by construction.
    # Root cause was localised to base64 padding in the
    # `public_const_fp = one_fp_b64_127` argument passed into the
    # `k2Ring127AffineCombineDS` call expression: the Opal DSL parser
    # rejects '=' / '+' / '/' inside expression strings under certain
    # column-position layouts, silently corrupting the (1-p_k) share
    # for the iteration where the layout boundary triggered. Patch
    # at commit 8302028 wraps the constant via .to_b64url before
    # passing it into the call; server-side `.b64_pad`
    # (ring127SplinelessDS.R:51) restores the padding before decode,
    # preserving the round-trip. The path was NOT a Ring127 sign-
    # extension primitive bug (H1) nor a cross-server aggregation
    # sign-bit issue (H2) nor an R-side buffer reuse (H3) -- it was a
    # plain base64-encoding mismatch at the orchestrator boundary.
    #
    # Post-fix L3 validation (10 iter, K=2, n=132, NHANES bp tertile,
    # 8302028): all four (k,l) blocks reported count_neg=0 across all
    # iterations; final max|Deltapi| = 4.99e-2 PRACTICAL, matching the
    # Bohning baseline (4.70e-2) within sampling-noise of the dataset
    # draw. The ~16x accuracy improvement promised by the L2 mock
    # (rel = 5.4e-4) does not transfer to L3: the noise floor at
    # fracBits=50 ring127 with the present Beaver matvec chain depth
    # is dominated by the gradient-side share-space arithmetic, which
    # both Bohning and H_emp share equally. Improving the Hessian
    # alone cannot push below this floor at L3.
    #
    # Architectural decision (reviewer 2026-04-29 option 2): ship
    # Bohning B_reg as default for production runs (~3x fewer MPC
    # rounds per iter for identical L3 numerical quality). The H_emp
    # scaffolding remains in source as opt-in via
    # `options(dsvert.mnl_joint_h_emp = TRUE)` for any future
    # deployment with a lower gradient-chain noise floor (e.g., longer
    # Beaver chains, deeper rings, mixed-precision splines) where the
    # Hessian quality becomes the binding constraint again.
    #
    # Non-disclosure invariants (D-INV-1/2/3) preserved by the
    # parser-fix patch: .to_b64url operates on a base64-encoded public
    # constant, never on per-patient values; no new emission category.
    use_h_emp <- isTRUE(getOption("dsvert.mnl_joint_h_emp", FALSE))
    H_emp_ok <- FALSE
    H_emp_full <- NULL
    if (use_h_emp) {
    # === Empirical Hessian H_emp_k per class via MPC X^T diag(W_k) X ===
    # Replaces Bohning B_reg's per-class diagonal blocks with the
    # empirical second-derivative form (Tutz 1990 Sec.3.2 closed-form
    # multinomial Hessian; Krishnapuram et al 2005 IEEE PAMI 27(6) Sec.3.2;
    # Friedman-Hastie-Tibshirani 2010 Sec.3.2; Hosmer-Lemeshow 2013 Sec.3.5).
    # Quadratic local convergence (Pratt 1981 + Burridge 1981 strict
    # concavity -> multinomial analog) replaces Bohning's O(1/epsilon)
    # Loewner-majorant rate (Bohning 1992 Thm 2 superseded as historical
    # context per project_h10_msle_noise_floor 2026-04-27).
    #
    # Block structure: H_emp in R^{(K-1)p x (K-1)p} stacked block matrix
    # where H_emp[k, l] is the pxp sub-block.
    #   H_emp_kk = +X^T diag(p_k(1-p_k)) X  (positive-definite diagonal)
    #   H_emp_kl = -X^T diag(p_k p_l) X      (k != l, symmetric off-diag)
    # Both expressed including the intercept column (full p including
    # intercept):
    #   H_emp_kl[alpha, alpha] = +Sum_i W_kl_i               (scalar intercept-intercept)
    #   H_emp_kl[alpha, beta_j] = +Sum_i W_kl_i * X_ij       (intercept-slope vector)
    #   H_emp_kl[beta_j, beta_q] = +Sum_i W_kl_i * X_ij*X_iq  (slope-slope matrix)
    # with sign sigma_kl = +1 for k=l, -1 for k!=l applied at assembly.
    #
    # Pipeline isomorphic to ord_joint K=2 H_emp (a0b1a65 / 754cfe08):
    #   (A) W_kl_share = Beaver vecmul(p_k_share, p_l_share) or
    #       p_k_share * (1 - p_k_share) for k=l. Pure share-space.
    #   (B) Sum W_kl via k2BeaverSumShareDS + k2-ring63-aggregate ring127.
    #   (C) X^T W_kl as length-p_shared vector via standard matvec
    #       pipeline (residual_key = W_kl_share).
    #   (D) X^T diag(W_kl) X as p_sharedxp_shared block, column-by-column
    #       via .ring127_vecmul(W_kl, X_:j) -> matvec pipeline.
    #
    # Disclosure (per Venables-Ripley 2002 Sec.7.4 + Aliasgari-Blanton 2013
    # NDSS Sec.5 + Demmler-ABY 2015 Sec.III.B): all (A) products stay share-
    # space; (B)(C)(D) reveal aggregates of shape O(p) x O(K) -- never
    # per-patient. Same K=2 audit boundary pattern as ord_joint H_emp,
    # already vetted at PR #4 merge.
    one_fp_b64_127 <- dsVert:::.callMpcTool("k2-float-to-fp", list(
      values = array(1.0, dim = 1L),
      frac_bits = 50L, ring = "ring127"))$fp_data
    # Strip '=' / '+' / '/' so the Opal DSL parser doesn't reject the
    # public_const_fp argument inside the `call(...)` expression.
    # Server-side `.b64_pad` (ring127SplinelessDS.R:51) restores the
    # padding before decoding.
    one_fp_b64_127 <- .to_b64url(one_fp_b64_127)

    # Build W_kl_share keys: diagonal (k=l) via p_k*(1-p_k), cross via
    # p_k*p_l (with sign applied at assembly).
    W_keys <- list()
    W_signs <- list()
    for (ki in seq_len(K_minus_1)) {
      for (li in seq_len(K_minus_1)) {
        if (ki > li) next  # symmetry: only build upper triangle + diag
        key <- sprintf("mnl_W_%d_%d_iter%d", ki, li, outer)
        if (ki == li) {
          one_minus_p_key <- sprintf("mnl_omp_%d_iter%d", ki, outer)
          for (srv in server_list) {
            ci <- which(server_names == srv)
            is_c <- (srv == coord)
            .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
              a_key = p_keys[ki], b_key = NULL,
              sign_a = -1L, sign_b = 0L,
              public_const_fp = one_fp_b64_127,
              is_party0 = is_c,
              output_key = one_minus_p_key,
              n = as.integer(n_obs), session_id = session_id))
          }
          .ring127_vecmul(p_keys[ki], one_minus_p_key, key, n_obs,
                          datasources, dealer_ci, server_list,
                          server_names, y_server, nl, transport_pks,
                          session_id, .dsAgg, .sendBlob)
          W_signs[[paste(ki, li, sep="_")]] <- +1
        } else {
          .ring127_vecmul(p_keys[ki], p_keys[li], key, n_obs,
                          datasources, dealer_ci, server_list,
                          server_names, y_server, nl, transport_pks,
                          session_id, .dsAgg, .sendBlob)
          W_signs[[paste(ki, li, sep="_")]] <- -1
        }
        W_keys[[paste(ki, li, sep="_")]] <- key
      }
    }

    # H_emp assembly: per (k, l) build p x p block via:
    #   - Sum_i W_kl_i scalar (intercept-intercept)
    #   - X^T W_kl length-p_shared vector (intercept-slope)
    #   - X^T diag(W_kl) X p_sharedxp_shared matrix (slope-slope)
    H_emp_full <- matrix(0, p * K_minus_1, p * K_minus_1)
    H_emp_ok <- TRUE
    int_j <- which(cnames == "(Intercept)")
    if (length(int_j) != 1L) int_j <- NA_integer_

    col_owner_map <- unlist(lapply(server_list, function(srv) {
      rep(srv, length(x_vars_per_server[[srv]]))
    }), use.names = FALSE)
    col_local_idx <- unlist(lapply(server_list, function(srv) {
      seq_along(x_vars_per_server[[srv]])
    }), use.names = FALSE)
    perm_grad_p <- match(formula_slope_names, server_partition_names)

    for (ki in seq_len(K_minus_1)) {
      for (li in seq_len(K_minus_1)) {
        kl_key <- if (ki <= li) paste(ki, li, sep="_") else paste(li, ki, sep="_")
        W_share_key <- W_keys[[kl_key]]
        sign_kl <- W_signs[[kl_key]]
        if (is.null(W_share_key)) { H_emp_ok <- FALSE; next }

        # (B) Sum W_kl scalar -- intercept-intercept entry of block.
        sum_W_shares <- list()
        for (srv in server_list) {
          ci <- which(server_names == srv)
          rs <- .dsAgg(datasources[ci], call(name = "k2BeaverSumShareDS",
            source_key = W_share_key, session_id = session_id,
            frac_bits = 50L))
          if (is.list(rs) && length(rs) == 1L) rs <- rs[[1L]]
          sum_W_shares[[srv]] <- rs$sum_share_fp
        }
        sum_W_agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
          share_a = sum_W_shares[[coord]],
          share_b = sum_W_shares[[nl]],
          frac_bits = 50L, ring = "ring127"))
        sum_W <- as.numeric(sum_W_agg$values)[1L]
        if (!is.finite(sum_W)) { H_emp_ok <- FALSE; next }

        # (C) X^T W_kl length-p_shared vector -- intercept-slope entries.
        for (srv in server_list) {
          ci <- which(server_names == srv)
          .dsAgg(datasources[ci], call(name = "dsvertPrepareMultinomGradDS",
            residual_key = W_share_key,
            is_outcome_server = (srv == coord),
            n = as.integer(n_obs), session_id = session_id))
        }
        grad_t_xw <- .dsAgg(datasources[dealer_ci],
          call(name = "glmRing63GenGradTriplesDS",
               dcf0_pk = .to_b64url(transport_pks[[coord]]),
               dcf1_pk = .to_b64url(transport_pks[[nl]]),
               n = as.integer(n_obs), p = as.integer(p_shared),
               ring = 127L, session_id = session_id))
        if (is.list(grad_t_xw) && length(grad_t_xw) == 1L)
          grad_t_xw <- grad_t_xw[[1L]]
        gtk_xw <- sprintf("k2_grad_triple_fp_iter%d_xw_%d_%d",
                           outer, ki, li)
        .sendBlob(grad_t_xw$grad_blob_0, gtk_xw, y_server_ci)
        .sendBlob(grad_t_xw$grad_blob_1, gtk_xw, dealer_ci)
        r1xw <- list()
        for (srv in server_list) {
          ci <- which(server_names == srv)
          peer <- setdiff(server_list, srv)
          .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS",
            session_id = session_id, grad_triple_key = gtk_xw))
          rr <- .dsAgg(datasources[ci], call(name = "k2GradientR1DS",
            peer_pk = .to_b64url(transport_pks[[peer]]), session_id = session_id))
          if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
          r1xw[[srv]] <- rr
        }
        .sendBlob(r1xw[[coord]]$encrypted_r1, "k2_grad_peer_r1", dealer_ci)
        .sendBlob(r1xw[[nl]]$encrypted_r1, "k2_grad_peer_r1", y_server_ci)
        r2xw <- list()
        for (srv in server_list) {
          ci <- which(server_names == srv)
          is_c <- (srv == coord)
          rr <- .dsAgg(datasources[ci], call(name = "k2GradientR2DS",
            party_id = if (is_c) 0L else 1L, session_id = session_id))
          if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
          r2xw[[srv]] <- rr
        }
        agg_xw <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
          share_a = r2xw[[coord]]$gradient_fp,
          share_b = r2xw[[nl]]$gradient_fp,
          frac_bits = 50L, ring = "ring127"))
        XtW_partition <- as.numeric(agg_xw$values)
        if (length(XtW_partition) != p_shared || any(!is.finite(XtW_partition))) {
          H_emp_ok <- FALSE; next
        }
        # Permute to formula order (slope rows of cnames).
        XtW_formula <- numeric(length(slope_rows))
        XtW_formula[] <- XtW_partition[perm_grad_p]

        # (D) X^T diag(W_kl) X p_shared x p_shared via column-by-column
        # matvec -- same pattern as ord_joint H_emp (a0b1a65) but
        # parameterised on the (k, l) pair.
        H_slope_part <- matrix(0, p_shared, p_shared)
        block_ok <- TRUE
        for (j in seq_len(p_shared)) {
          Xj_key  <- sprintf("mnl_Xj_iter%d_blk%d_%d_col%d", outer, ki, li, j)
          DXj_key <- sprintf("mnl_DXj_iter%d_blk%d_%d_col%d", outer, ki, li, j)
          owner_srv <- col_owner_map[j]
          local_j   <- as.integer(col_local_idx[j])
          for (srv in server_list) {
            ci <- which(server_names == srv)
            is_owner <- (srv == owner_srv)
            mat_key <- if (is_owner) "k2_x_share_fp" else "k2_peer_x_share_fp"
            p_local <- if (is_owner) length(x_vars_per_server[[srv]])
                       else length(x_vars_per_server[[owner_srv]])
            .dsAgg(datasources[ci],
              call(name = "dsvertOrdinalExtractXColumnDS",
                   matrix_key = mat_key,
                   n = as.integer(n_obs),
                   p = as.integer(p_local),
                   col_idx = local_j,
                   output_key = Xj_key,
                   session_id = session_id))
          }
          .ring127_vecmul(W_share_key, Xj_key, DXj_key, n_obs,
            datasources, dealer_ci, server_list, server_names,
            y_server, nl, transport_pks, session_id, .dsAgg, .sendBlob)
          for (srv in server_list) {
            ci <- which(server_names == srv)
            .dsAgg(datasources[ci],
              call(name = "dsvertPrepareMultinomGradDS",
                   residual_key = DXj_key,
                   is_outcome_server = (srv == coord),
                   n = as.integer(n_obs), session_id = session_id))
          }
          grad_t_H <- .dsAgg(datasources[dealer_ci],
            call(name = "glmRing63GenGradTriplesDS",
                 dcf0_pk = .to_b64url(transport_pks[[coord]]),
                 dcf1_pk = .to_b64url(transport_pks[[nl]]),
                 n = as.integer(n_obs), p = p_shared,
                 ring = 127L, session_id = session_id))
          if (is.list(grad_t_H) && length(grad_t_H) == 1L)
            grad_t_H <- grad_t_H[[1L]]
          gtk_H <- sprintf("k2_grad_triple_fp_iter%d_Hblk%d_%d_col%d",
                            outer, ki, li, j)
          .sendBlob(grad_t_H$grad_blob_0, gtk_H, y_server_ci)
          .sendBlob(grad_t_H$grad_blob_1, gtk_H, dealer_ci)
          r1H <- list()
          for (srv in server_list) {
            ci <- which(server_names == srv)
            peer <- setdiff(server_list, srv)
            .dsAgg(datasources[ci], call(name = "k2StoreGradTripleDS",
              session_id = session_id, grad_triple_key = gtk_H))
            rr <- .dsAgg(datasources[ci], call(name = "k2GradientR1DS",
              peer_pk = .to_b64url(transport_pks[[peer]]), session_id = session_id))
            if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
            r1H[[srv]] <- rr
          }
          .sendBlob(r1H[[coord]]$encrypted_r1, "k2_grad_peer_r1", dealer_ci)
          .sendBlob(r1H[[nl]]$encrypted_r1, "k2_grad_peer_r1", y_server_ci)
          r2H <- list()
          for (srv in server_list) {
            ci <- which(server_names == srv)
            is_c <- (srv == coord)
            rr <- .dsAgg(datasources[ci], call(name = "k2GradientR2DS",
              party_id = if (is_c) 0L else 1L, session_id = session_id))
            if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
            r2H[[srv]] <- rr
          }
          agg_H <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
            share_a = r2H[[coord]]$gradient_fp,
            share_b = r2H[[nl]]$gradient_fp,
            frac_bits = 50L, ring = "ring127"))
          col_vals <- as.numeric(agg_H$values)
          if (length(col_vals) != p_shared || any(!is.finite(col_vals))) {
            block_ok <- FALSE; break
          }
          # Permute revealed column to formula order (slopes only).
          H_slope_part[, j] <- col_vals
        }
        if (!block_ok) { H_emp_ok <- FALSE; next }
        # Permute rows AND cols to formula slope order.
        H_slope_formula <- H_slope_part[perm_grad_p, perm_grad_p, drop = FALSE]
        # Symmetrise per-block diagonal noise.
        H_slope_formula <- (H_slope_formula + t(H_slope_formula)) / 2

        # Assemble p x p block with intercept row/col + slope-slope.
        H_block <- matrix(0, p, p, dimnames = list(cnames, cnames))
        H_block[slope_rows, slope_rows] <- H_slope_formula
        if (!is.na(int_j)) {
          H_block[int_j, int_j] <- sum_W
          H_block[int_j, slope_rows] <- XtW_formula
          H_block[slope_rows, int_j] <- XtW_formula
        }
        # Apply sign sigma_kl and average by n_obs.
        H_block_signed <- (sign_kl * H_block) / n_obs
        # Place into (K-1)p x (K-1)p stacked H_emp.
        rows_kl <- ((ki - 1L) * p + 1L):(ki * p)
        cols_kl <- ((li - 1L) * p + 1L):(li * p)
        H_emp_full[rows_kl, cols_kl] <- H_block_signed
        # Per-block aggregate-only audit diagnostic for the H_emp
        # cross-block sign-anomaly investigation (dsVert#4). All
        # quantities printed here are already revealed at the audit
        # boundary via `k2-ring63-aggregate` (sum_W, XtW_formula) +
        # column-by-column aggregate matvec (H_slope_formula);
        # printing summary stats does NOT introduce a new emission
        # category (D-INV-1/2 preserved). Fields:
        #   sum_W       -- scalar Sum_i W_kl_i (already revealed)
        #   |XtW|       -- max abs of revealed length-p_shared vector
        #   H_diag_*    -- min/max of revealed slope-slope diagonal
        #   sign_kl     -- client-side +/-1 sign factor
        H_diag_vals <- diag(H_slope_formula)
        cat(sprintf(
          "[MnlJoint iter %d block (k=%d,l=%d)] sum_W=% .3e |XtW|_max=%.3e H_diag min=% .3e max=% .3e count_neg=%d sign_kl=%+d\n",
          outer, ki, li, sum_W, max(abs(XtW_formula)),
          min(H_diag_vals), max(H_diag_vals),
          sum(H_diag_vals < 0), sign_kl))
      }
    }
    # Symmetrise the full assembled H_emp (numerical noise across blocks).
    if (H_emp_ok) {
      H_emp_full <- (H_emp_full + t(H_emp_full)) / 2
      cat(sprintf("[MnlJoint iter %d] H_emp diag=[%s]\n",
                   outer,
                   paste(sprintf("%.3e", diag(H_emp_full)), collapse=",")))
    }
    } # end if (use_h_emp) -- H_emp pipeline gated FALSE; B_reg used below

    # Client-side empirical-Hessian Newton step (replaces Bohning B_reg
    # when H_emp is available; falls back to Bohning for stability when
    # H_emp construction fails).
    g_stacked <- as.numeric(gradients)
    g_norm <- sqrt(sum(g_stacked^2))
    g_max  <- max(abs(g_stacked))
    # Track best beta seen so far (argmin_k |g|_L2). Newton's late-iter
    # oscillation under MPC step-cap binds the step but may not track
    # the likelihood maximum monotonically -- best-beta recovers the best
    # point encountered.
    if (is.finite(g_norm) && g_norm < best_g_norm) {
      best_g_norm <- g_norm
      best_beta <- beta_mat
      best_iter <- outer - 1L  # beta_mat is current iterate PRE step
    }
    # Newton solve: use empirical H_emp_full if successfully assembled
    # (Tutz 1990 Sec.3.2 quadratic local convergence per Pratt 1981 +
    # Burridge 1981); fall back to Bohning B_reg majorant otherwise
    # (monotone descent guarantee per Bohning 1992 Thm 2). Christensen
    # 2019 ordinal::clm.fit Sec.A.3 diagonal eigenvalue inflation +
    # ridge epsilon*I for numerical stability of the empirical-H solve.
    H_solve <- if (H_emp_ok) {
      ridge <- 1e-6 * max(abs(diag(H_emp_full)), 1)
      H_emp_full + ridge * diag(p * K_minus_1)
    } else {
      B_reg
    }
    step_stacked <- tryCatch(solve(H_solve, g_stacked),
                              error = function(e) {
                                # Empirical solve failed -> fall back to
                                # Bohning B_reg (Loewner upper-bound
                                # always positive-definite per Bohning
                                # 1992 Thm 2).
                                tryCatch(solve(B_reg, g_stacked),
                                          error = function(e2) 0.1 * g_stacked)
                              })
    step_mat <- matrix(step_stacked, p, K_minus_1,
                       dimnames = dimnames(gradients))
    step_norm_pre <- max(abs(step_mat))
    # Raw gradient norms can keep shrinking after the MPC approximation floor
    # is reached even when a further Newton step moves away from the central
    # softmax MLE. The aggregate Newton step norm |H^{-1}g| is a stronger
    # stationarity proxy under the Bohning majorant because it includes the
    # current curvature scaling and costs no additional disclosure.
    if (is.finite(step_norm_pre) && step_norm_pre < best_step_norm) {
      best_step_norm <- step_norm_pre
      best_step_beta <- beta_mat
      best_step_iter <- outer - 1L
    }
    # Decreasing step-cap schedule: 0.5 for first 5 iters (broad
    # descent), then 0.5 * 0.7^(iter-5) (refine near optimum, reduce
    # oscillation amplitude as |g| shrinks). Bounded below at 0.05 so
    # Newton always makes SOME progress. Same spirit as Nocedal-Wright
    # Sec.3.5 backtracking but pre-scheduled to save the Armijo MPC round.
    step_cap <- if (outer <= 5L) 0.5 else max(0.5 * 0.7^(outer - 5L), 0.05)
    if (is.finite(step_norm_pre) && step_norm_pre > step_cap) {
      step_mat <- step_mat * (step_cap / step_norm_pre)
    }
    beta_new <- beta_mat
    beta_new[, non_ref] <- beta_new[, non_ref] + step_mat
    max_step <- max(abs(step_mat))
    # AUDITORIA instrumentation: gradient-norm trace per iter. Always
    # emit (not gated on verbose) so CV full logs capture convergence
    # trajectory for FAIL diagnosis.
    cat(sprintf("[MnlJoint] iter %d  |g|_L2=%.3e  |g|_max=%.3e  |step|_pre=%.3e  |step|_post=%.3e  beta_max=%.3e  dt=%.1fs\n",
                outer, g_norm, g_max, step_norm_pre, max_step,
                max(abs(beta_new)), proc.time()[[3L]] - t_iter))
    if (verbose)
      message(sprintf("[MultinomJointNewton] iter %d |step|_max=%.3e (%.1fs)",
                       outer, max_step, proc.time()[[3L]] - t_iter))
    beta_mat <- beta_new
    if (max_step < tol) {
      converged <- TRUE
      final_iter <- outer
      break
    }
  }

  # Cleanup MPC session
  for (srv in server_list) {
    ci <- which(server_names == srv)
    try(.dsAgg(datasources[ci],
               call(name = "mpcCleanupDS", session_id = session_id)),
        silent = TRUE)
  }

  # Return the best stationarity iterate, not necessarily the final iterate.
  # Raw |g| and final beta can both be poor selectors after the Ring127
  # softmax approximation floor is reached. The minimum aggregate Newton
  # step norm is the default selector; the minimum raw-gradient beta is
  # retained for diagnostics.
  gradient_beta <- if (!is.null(best_beta)) best_beta else beta_mat
  final_beta <- if (!is.null(best_step_beta)) best_step_beta else gradient_beta
  out <- warm
  out$coefficients_anchored <- warm$coefficients
  out$coefficients <- final_beta
  out$coefficients_best_gradient <- gradient_beta
  out$coefficients_final_iter <- beta_mat  # last iterate for diagnostics
  out$best_g_norm <- best_g_norm
  out$best_iter <- best_iter
  out$best_step_norm <- best_step_norm
  out$best_step_iter <- best_step_iter
  out$returned_selection <- if (!is.null(best_step_beta)) "best_step"
                            else if (!is.null(best_beta)) "best_gradient"
                            else "final"
  out$outer_iter <- final_iter
  out$converged <- converged
  out$family <- "multinomial_joint_softmax_ring127"
  out$session_id <- session_id
  class(out) <- c("ds.vertMultinomJointNewton", class(out))
  out
}

#' @export
print.ds.vertMultinomJointNewton <- function(x, ...) {
  cat("dsVert joint-softmax multinomial (Ring127 Newton)\n")
  cat(sprintf("  N = %d  classes = %s  outer_iter = %d  converged = %s\n",
              x$n_obs, paste(x$classes, collapse = ","),
              x$outer_iter, x$converged))
  cat("\nCoefficients (post joint-Newton):\n")
  print(round(x$coefficients, 4))
  invisible(x)
}
