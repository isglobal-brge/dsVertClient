#' @title Federated joint-softmax multinomial logistic regression via
#'   Ring127 MPC-orchestrated Newton iteration
#' @description Full per-patient softmax Newton: orchestrates K−1 parallel
#'   exp(η_k) shares, sums to denominator D = 1+Σexp(η_k), computes 1/D
#'   via Ring127 Chebyshev + Newton-Raphson, multiplies per class to get
#'   p_k(x_i) share per patient, builds residual y_k − p_k on outcome
#'   server, and aggregates X^T(y_k − p_k) via existing Beaver matvec
#'   pipeline for each class. Client-side Bohning-Hessian-bounded Newton
#'   step on stacked β.
#'
#'   All per-patient quantities stay as Ring127 additive shares; only the
#'   final p(K−1)-dim aggregate gradient per iter is revealed — same
#'   privacy class as the single-class gradient of ds.vertGLM. **P3 delta:
#'   zero new reveal types.**
#'
#' @param formula R formula with categorical outcome on LHS.
#' @param data Aligned data frame name.
#' @param levels Character vector of outcome levels (first = reference).
#' @param indicator_template sprintf template for class-indicator columns
#'   on the outcome server, e.g. \code{"\%s_ind"}. Must exist server-side.
#' @param max_outer Outer Newton iterations (default 8).
#' @param tol Convergence tolerance on max |Δβ| (default 1e-4).
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
      call("glmRing63TransportInitDS", session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    transport_pks[[srv]] <- r$transport_pk
  }
  .json_to_b64url <- function(x) {
    raw <- charToRaw(jsonlite::toJSON(x, auto_unbox = TRUE))
    b64 <- gsub("\n", "", jsonlite::base64_enc(raw), fixed = TRUE)
    chartr("+/", "-_", sub("=+$", "", b64, perl = TRUE))
  }
  for (srv in server_list) {
    ci <- which(server_names == srv)
    peer_srv <- setdiff(server_list, srv)
    peers <- setNames(list(transport_pks[[peer_srv]]), peer_srv)
    DSI::datashield.aggregate(datasources[ci],
      call("mpcStoreTransportKeysDS",
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
      call("dsvertColNamesDS", data_name = data)),
      error = function(e) NULL)
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    cols_here <- if (is.list(r) && !is.null(r$columns)) r$columns
                 else if (is.character(r)) r else character(0)
    x_vars_per_server[[srv]] <- intersect(rhs, cols_here)
  }

  .dsAgg <- function(conns, expr, ...)
    DSI::datashield.aggregate(conns, expr, ...)
  .sendBlob <- function(blob, key, conn_idx) {
    if (is.null(blob) || !nzchar(blob)) return(invisible())
    DSI::datashield.aggregate(datasources[conn_idx],
      call("mpcStoreBlobDS", blob_b64 = blob, blob_key = key,
           session_id = session_id))
  }

  # Share input at Ring127
  share_results <- list()
  for (srv in server_list) {
    ci <- which(server_names == srv)
    peer <- setdiff(server_list, srv)
    peer_pk <- transport_pks[[peer]]
    r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
      data_name = data, x_vars = x_vars_per_server[[srv]],
      y_var = if (srv == y_server) y_var_char else NULL,
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
    .dsAgg(datasources[ci], call("k2ReceiveShareDS",
      peer_p = as.integer(length(x_vars_per_server[[peer]])),
      session_id = session_id))
  }
  n_obs <- share_results[[y_server]]$n
  if (verbose) message(sprintf("[MultinomJointNewton] session %s  n=%d  K=%d  p=%d",
                                session_id, n_obs, K_minus_1 + 1L, p))

  converged <- FALSE
  final_iter <- max_outer

  # Bohning Hessian upper bound (client-side, constant): use any warm$fits
  # covariance to approximate X^T X / n, then form block kronecker.
  cov_k0 <- warm$fits[[non_ref[1L]]]$covariance
  sigma2 <- if (!is.null(warm$fits[[non_ref[1L]]]$deviance))
    warm$fits[[non_ref[1L]]]$deviance /
      max(n_obs - p, 1L) else 1
  XtX_over_n <- tryCatch(sigma2 * solve(cov_k0) / n_obs,
                         error = function(e) diag(p))
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
        .dsAgg(datasources[ci], call("k2ComputeEtaShareDS",
          beta_coord = b_coord_vec, beta_nl = b_nl_vec,
          intercept = if (srv == coord) intercept_k else 0,
          is_coordinator = (srv == coord),
          session_id = session_id, output_key = eta_key))
      }
      exp_key <- paste0("exp_eta_class_", ki)
      .ring127_exp_round_keyed(eta_key, exp_key, n_obs,
                               datasources, dealer_ci, server_list,
                               server_names, y_server, nl, transport_pks,
                               session_id, .dsAgg, .sendBlob)
      exp_eta_keys[ki] <- exp_key
    }

    # Step 3: D = 1 + sum(exp(eta_k))
    for (srv in server_list) {
      ci <- which(server_names == srv)
      .dsAgg(datasources[ci], call("dsvertSoftmaxDenominatorDS",
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

    # Step 6-7: residual r_k = y_k - p_k + Beaver matvec X^T r_k
    gradients <- matrix(0, p, K_minus_1,
                         dimnames = list(cnames, non_ref))
    for (ki in seq_along(non_ref)) {
      k <- non_ref[ki]
      r_key <- paste0("r_class_", ki)
      for (srv in server_list) {
        ci <- which(server_names == srv)
        .dsAgg(datasources[ci], call("dsvertComputeResidualShareDS",
          p_key = p_keys[ki],
          indicator_col = if (srv == coord) sprintf(indicator_template, k) else NULL,
          data_name = if (srv == coord) data else NULL,
          output_key = r_key,
          is_outcome_server = (srv == coord),
          n = as.integer(n_obs), session_id = session_id))
      }
      # Prep the standard gradient machinery: pipeline computes X^T(mu-y),
      # we set mu=0 and y=-r_k so X^T(mu-y) = X^T r_k.
      for (srv in server_list) {
        ci <- which(server_names == srv)
        .dsAgg(datasources[ci], call("dsvertPrepareMultinomGradDS",
          residual_key = r_key,
          is_outcome_server = (srv == coord),
          n = as.integer(n_obs), session_id = session_id))
      }
      # Beaver matvec: X^T (mu - y) → per-class gradient
      grad_t <- .dsAgg(datasources[dealer_ci],
        call("glmRing63GenGradTriplesDS",
             dcf0_pk = transport_pks[[coord]],
             dcf1_pk = transport_pks[[nl]],
             n = as.integer(n_obs), p = as.integer(p),
             ring = 127L, session_id = session_id))
      if (is.list(grad_t) && length(grad_t) == 1L) grad_t <- grad_t[[1L]]
      .sendBlob(grad_t$grad_blob_0, "k2_grad_triple_fp", y_server_ci)
      .sendBlob(grad_t$grad_blob_1, "k2_grad_triple_fp", dealer_ci)
      r1 <- list()
      for (srv in server_list) {
        ci <- which(server_names == srv)
        peer <- setdiff(server_list, srv)
        .dsAgg(datasources[ci], call("k2StoreGradTripleDS",
          session_id = session_id))
        rr <- .dsAgg(datasources[ci], call("k2GradientR1DS",
          peer_pk = transport_pks[[peer]], session_id = session_id))
        if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
        r1[[srv]] <- rr
      }
      .sendBlob(r1[[coord]]$encrypted_r1, "k2_grad_peer_r1", dealer_ci)
      .sendBlob(r1[[nl]]$encrypted_r1, "k2_grad_peer_r1", y_server_ci)
      r2 <- list()
      for (srv in server_list) {
        ci <- which(server_names == srv)
        is_c <- (srv == coord)
        rr <- .dsAgg(datasources[ci], call("k2GradientR2DS",
          party_id = if (is_c) 0L else 1L, session_id = session_id))
        if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
        r2[[srv]] <- rr
      }
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = r2[[coord]]$gradient_fp,
        share_b = r2[[nl]]$gradient_fp,
        frac_bits = 50L, ring = "ring127"))
      # Aggregate gradient is X^T(mu-y) = X^T(0 - (-r_k)) = X^T r_k → good
      gradients[, ki] <- as.numeric(agg$values) / n_obs
    }

    # Client-side Bohning Newton step
    g_stacked <- as.numeric(gradients)
    step_stacked <- tryCatch(solve(B_reg, g_stacked),
                              error = function(e) 0.1 * g_stacked)
    step_mat <- matrix(step_stacked, p, K_minus_1,
                       dimnames = dimnames(gradients))
    beta_new <- beta_mat
    beta_new[, non_ref] <- beta_new[, non_ref] + step_mat
    max_step <- max(abs(step_mat))
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
               call("mpcCleanupDS", session_id = session_id)),
        silent = TRUE)
  }

  out <- warm
  out$coefficients_anchored <- warm$coefficients
  out$coefficients <- beta_mat
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
