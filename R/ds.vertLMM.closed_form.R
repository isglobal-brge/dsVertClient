# LMM is retained as a quarantined compatibility route while its statistical
# release contract is redesigned.  Even so, every remote phase must use the
# same fail-closed DSI boundary as promoted routes: named results, transport
# validation and typed identity/session failures. This helper adds no LMM
# endpoint to the replay allowlist; generic typed producers replay only when
# already classified as idempotent by the shared transport contract.
.dsvert_lmm_aggregate_strict <- function(conns, expr) {
  endpoints <- unique(.dsvert_dsi_call_names(expr))
  endpoints <- endpoints[!is.na(endpoints) & nzchar(endpoints)]
  operation <- if (length(endpoints)) {
    paste0("LMM ", paste(endpoints, collapse = "+"))
  } else {
    "LMM remote phase"
  }
  .dsvert_aggregate_strict(
    conns = conns, expr = expr, operation = operation)
}

# A few historical LMM optional branches deliberately convert ordinary
# capability/numeric failures into a local fallback. Identity mismatch and an
# unresolved DSI session are never optional and must cross those boundaries
# with their typed class intact.
.dsvert_lmm_transport_fallback <- function(condition, fallback = NULL) {
  if (inherits(condition, c("dsvert_peer_not_recognized",
                            "dsvert_dsi_poisoned_session"))) {
    stop(condition)
  }
  fallback
}

# Complete the public fixed-effect fit from the protected GLS sufficient
# statistics. The closed-form path replaces the inner Gaussian fit, so its
# covariance must come from the same transformed information matrix rather
# than from the preliminary OLS fit.
.dsvert_lmm_closed_form_fit <- function(cf_fit, coefficients, sigma2) {
  coefficient_names <- names(coefficients)
  information <- as.matrix(cf_fit$XtX)
  if (!is.numeric(coefficients) || !length(coefficients) ||
      is.null(coefficient_names) || any(!nzchar(coefficient_names)) ||
      anyDuplicated(coefficient_names) || any(!is.finite(coefficients)) ||
      !is.numeric(information) || length(dim(information)) != 2L ||
      nrow(information) != ncol(information) ||
      is.null(rownames(information)) || is.null(colnames(information)) ||
      !all(coefficient_names %in% rownames(information)) ||
      !all(coefficient_names %in% colnames(information)) ||
      length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0) {
    .dsvert_stop_non_identifiable(
      "The protected LMM fixed-effect information is incomplete.",
      reason = "invalid_lmm_fixed_effect_information")
  }
  information <- information[
    coefficient_names, coefficient_names, drop = FALSE]
  covariance <- sigma2 * .dsvert_solve_identifiable(
    information,
    context = "The protected LMM fixed-effect information",
    reason = "singular_lmm_fixed_effect_information",
    symmetric = TRUE)
  covariance <- (covariance + t(covariance)) / 2
  dimnames(covariance) <- list(coefficient_names, coefficient_names)
  standard_errors <- sqrt(diag(covariance))
  if (any(!is.finite(standard_errors)) || any(standard_errors <= 0)) {
    .dsvert_stop_non_identifiable(
      "The protected LMM fixed-effect covariance is invalid.",
      reason = "invalid_lmm_fixed_effect_covariance")
  }
  list(
    coefficients = coefficients,
    covariance = covariance,
    std_errors = stats::setNames(standard_errors, coefficient_names))
}

#' @title LMM closed-form GLS driver
#' @description Computes the exact Laird-Ware generalised-least-squares
#'   estimate for a random-intercept LMM without relying on
#'   \code{ds.vertGLM}'s L-BFGS iteration. Directly assembles the
#'   Gram matrix X'X and right-hand side X'y over the cluster-mean-
#'   centred design via:
#'     - local per-server block computations
#'     - Beaver vecmul + FP-sum for cross-server entries
#'   and solves \code{beta = solve(XtX, Xty)} client-side. Matches
#'   \code{lme4::lmer} to FP precision (~1e-5 on small n).
#'
#'   This internal compatibility driver reconstructs the exact Gram matrix
#'   and right-hand side at the analyst client. Its paired per-cluster share
#'   consumers likewise reconstruct exact cluster-level residual moments.
#'   Transport encryption protects peer-to-peer operands; it does not make
#'   those client-side exact openings differentially private or
#'   non-reconstructible. The driver is quarantined and must not be treated as
#'   a promoted biomedical release route.
#'
#'   Assumptions:
#'     - Cluster IDs have been broadcast via
#'       \code{dsvertLMMBroadcastClusterIDsDS} /
#'       \code{dsvertLMMReceiveClusterIDsDS}.
#'     - Transport keys are initialised in the session.
#'     - Exactly K=2 servers (outcome + one peer).
#'
#' @param conns Named list of DS connections.
#' @param server_names Character vector of server names (names of \code{conns}).
#' @param y_srv,peer_srv Character. Outcome server name and peer server name.
#' @param data Aligned data-frame name on both servers.
#' @param y_var Response column (on \code{y_srv}).
#' @param x_ysrv,x_peer Character vectors of predictor names on each
#'   server (after formula parsing).
#' @param lambda_i Numeric vector of length n_clusters.
#' @param transport_pks Named list of public keys per server.
#' @param session_id Active MPC session (cluster IDs stored within).
#' @return list(coefficients, XtX, Xty, yty, n)
#' @keywords internal
.ds_vertLMM_closed_form <- function(conns, server_names, y_srv, peer_srv,
                                     data, y_var, x_ysrv, x_peer,
                                     lambda_i, transport_pks, session_id,
                                     verbose = FALSE,
                                     share_scale = 1.0,
                                     column_scales = NULL,
                                     standardize = FALSE,
                                     ring = "ring63") {
  # Task #121 LMM Ring127 migration: when ring=="ring127", Beaver
  # element-wise products use fracBits=50 (Uint128 shares) -- a ~2^30
  # reduction in per-op truncation noise vs Ring63 fracBits=20. All
  # Beaver/gen-triples/aggregate DS ops already dispatch on `ring`.
  ring_tag <- if (is.character(ring)) ring else "ring63"
  if (!(ring_tag %in% c("ring63", "ring127"))) {
    stop(".ds_vertLMM_closed_form: ring must be 'ring63' or 'ring127'",
         call. = FALSE)
  }
  fb <- if (ring_tag == "ring127") 50L else 20L
  # The LMMGram DS wrappers accept the string form; the OT-Beaver
  # preprocessing bridge uses the integer ring width.
  ring_int <- if (ring_tag == "ring127") 127L else 63L
  if (is.null(peer_srv))
    stop("Closed-form LMM solver currently requires K=2 (peer_srv != NULL)",
         call. = FALSE)

  ysrv_ci <- which(server_names == y_srv)
  peer_ci <- which(server_names == peer_srv)

  int_col <- "dsvertlmmint"

  # SNR-boost: pre-multiply shared columns by `share_scale` before Beaver
  # cross products, then de-scale the final XtX/Xty by c^2 client-side.
  # Because beta = solve(c^2XtX, c^2Xty) = solve(XtX, Xty) the GLS solution is
  # invariant but the Beaver product entries have c^2x larger magnitude
  # vs the fixed absolute Ring63 FP noise -- relative precision improves
  # by a factor c^2. For c=8, this turns ~1.8e-4 relative on X4 into
  # ~3e-6, well below the plan's rel<1e-4 bar for |beta|>1.
  # See docs/acceptance/path_b_targets.md Sec.LMM iterative-refinement
  # band-aid for the full rationale and escape hatch if this isn't
  # enough (escalate to task #111 B_full).
  sc <- as.numeric(share_scale)
  if (!is.finite(sc) || sc <= 0) sc <- 1.0

  # --- Phase 1: local Gram blocks + share transformed columns. -------
  if (verbose) message(sprintf(
    "[closed_form] Phase 1: local gram + shares (share_scale=%.2f)", sc))

  # Post-centering L2 standardization is enabled by default (the
  # Codex-approved structural fix 2026-04-19 for kappa=5.57e5 Gram
  # ill-conditioning). Server computes L2 of each centered column and
  # divides by it, then returns `l2_scales` so the client can unscale beta.
  # Default ON; caller can pass `standardize = FALSE` for debugging.
  std_flag <- isTRUE(standardize)

  # Outcome and peer transforms are independent and launch in one named DSI
  # fan-out phase. The outcome creates the intercept/y shares; the peer only
  # transforms its predictor block.
  gram_conns <- conns[c(ysrv_ci, peer_ci)]
  gram_expressions <- stats::setNames(list(
    call(name = "dsvertLMMLocalGramDS",
         data_name = data,
         columns = as.character(x_ysrv),
         y_var = y_var,
         lambda_per_cluster = as.numeric(lambda_i),
         create_intercept = TRUE,
         intercept_col = int_col,
         peer_pk = transport_pks[[peer_srv]],
         session_id = session_id,
         share_scale = sc,
         standardize = std_flag,
         ring = ring_tag),
    call(name = "dsvertLMMLocalGramDS",
         data_name = data,
         columns = as.character(x_peer),
         y_var = NULL,
         lambda_per_cluster = as.numeric(lambda_i),
         create_intercept = FALSE,
         peer_pk = transport_pks[[y_srv]],
         session_id = session_id,
         share_scale = sc,
         standardize = std_flag,
         ring = ring_tag)), c(y_srv, peer_srv))
  gram_results <- .dsvert_fanout_by_site(
    gram_conns, gram_expressions, operation = "LMM local Gram preparation")
  local_y <- gram_results[[y_srv]]
  local_p <- gram_results[[peer_srv]]

  # Relay peer blobs to the opposite party.
  .dsvert_store_blob(
    local_y$peer_blob, "k2_lmm_gram_peer_shares",
    conns[peer_ci], session_id)
  .dsvert_store_blob(
    local_p$peer_blob, "k2_lmm_gram_peer_shares",
    conns[ysrv_ci], session_id)
  # Each party ingests the other's column shares in one fan-out phase.
  .dsvert_fanout_by_site(
    gram_conns,
    stats::setNames(rep(list(call(
      name = "dsvertLMMReceiveGramSharesDS", session_id = session_id)), 2L),
      c(y_srv, peer_srv)),
    operation = "LMM Gram-share receipt")

  y_key <- local_y$y_key
  ysrv_cols <- as.character(local_y$column_names)   # int_col + x_ysrv
  peer_cols <- as.character(local_p$column_names)   # x_peer
  p_y <- length(ysrv_cols); p_p <- length(peer_cols)
  n <- as.integer(local_y$n)
  if (n != as.integer(local_p$n))
    stop("server row counts disagree: ", n, " vs ", local_p$n, call. = FALSE)

  # --- Phase 2: cross-server Beaver dot products. --------------------
  # pair list: all (a in peer_cols) x (b in ysrv_cols) + (a in peer_cols, y).
  pairs <- list()
  for (a in peer_cols) {
    for (b in ysrv_cols) pairs[[length(pairs) + 1L]] <- list(a = a, b = b)
  }
  for (a in peer_cols) pairs[[length(pairs) + 1L]] <- list(a = a, b = y_key)
  if (verbose) message(sprintf("[closed_form] Phase 2: %d Beaver dot products",
                                length(pairs)))

  cross_results <- vector("numeric", length(pairs))

  send_blob <- function(blob, key, conn_idx) {
    .dsvert_store_blob(blob, key, conns[conn_idx], session_id)
  }
  ds_agg <- function(ds, expr) .dsvert_lmm_aggregate_strict(ds, expr)

  # One Beaver vecmul per pair (sequential; each pair uses a fresh
  # OT-generated triple and then the existing online R1/R2 path).
  for (k in seq_along(pairs)) {
    pr <- pairs[[k]]
    a_col <- pr$a; b_col <- pr$b
    .ot_beaver_prepare_vecmul(
      datasources = conns,
      party_conns = c(ysrv_ci, peer_ci),
      party_names = c(y_srv, peer_srv),
      transport_pks = transport_pks,
      session_id = session_id,
      n = n,
      ring = ring_int,
      .dsAgg = ds_agg,
      .sendBlob = send_blob)
    # R1 on both parties.
    r1 <- .dsvert_fanout_by_site(
      gram_conns,
      stats::setNames(list(
        call(name = "dsvertLMMGramR1DS",
             peer_pk = transport_pks[[peer_srv]],
             x_col = a_col, y_col = b_col,
             session_id = session_id, frac_bits = fb,
             ring = ring_tag),
        call(name = "dsvertLMMGramR1DS",
             peer_pk = transport_pks[[y_srv]],
             x_col = a_col, y_col = b_col,
             session_id = session_id, frac_bits = fb,
             ring = ring_tag)), c(y_srv, peer_srv)),
      operation = "LMM Gram R1")
    r1_y <- r1[[y_srv]]
    r1_p <- r1[[peer_srv]]
    # Relay masks between parties.
    .dsvert_store_blob(
      r1_y$peer_blob, "k2_beaver_vecmul_peer_masked",
      conns[peer_ci], session_id)
    .dsvert_store_blob(
      r1_p$peer_blob, "k2_beaver_vecmul_peer_masked",
      conns[ysrv_ci], session_id)
    # R2 on both parties reduces to scalar shares.
    r2 <- .dsvert_fanout_by_site(
      gram_conns,
      stats::setNames(list(
        call(name = "dsvertLMMGramR2DS",
             is_party0 = TRUE,
             x_col = a_col, y_col = b_col,
             session_id = session_id, frac_bits = fb,
             ring = ring_tag),
        call(name = "dsvertLMMGramR2DS",
             is_party0 = FALSE,
             x_col = a_col, y_col = b_col,
             session_id = session_id, frac_bits = fb,
             ring = ring_tag)), c(y_srv, peer_srv)),
      operation = "LMM Gram R2")
    r2_y <- r2[[y_srv]]
    r2_p <- r2[[peer_srv]]
    # Aggregate two scalar shares into the true dot product.
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r2_y$scalar_share,
      share_b = r2_p$scalar_share, frac_bits = fb,
      ring = ring_tag))
    cross_results[k] <- as.numeric(agg$values[1L])
  }

  # --- Phase 3: assemble full Gram matrix. ---------------------------
  all_cols <- c(ysrv_cols, peer_cols)
  p_total <- p_y + p_p
  XtX <- matrix(0, p_total, p_total)
  rownames(XtX) <- colnames(XtX) <- all_cols
  # y_srv block.
  XtX[seq_len(p_y), seq_len(p_y)] <- as.matrix(local_y$XtX_local)
  # peer block.
  XtX[p_y + seq_len(p_p), p_y + seq_len(p_p)] <- as.matrix(local_p$XtX_local)
  # Cross block (peer x ysrv).
  idx <- 1L
  for (a in peer_cols) {
    for (b in ysrv_cols) {
      XtX[p_y + which(peer_cols == a), which(ysrv_cols == b)] <-
        cross_results[idx]
      XtX[which(ysrv_cols == b), p_y + which(peer_cols == a)] <-
        cross_results[idx]
      idx <- idx + 1L
    }
  }
  # Build Xty: y_srv part is local_y$Xty_local, peer part from the last
  # p_p entries of cross_results (pairs (a, y_key) for a in peer_cols).
  Xty <- numeric(p_total)
  Xty[seq_len(p_y)] <- as.numeric(local_y$Xty_local)
  for (j in seq_len(p_p)) {
    Xty[p_y + j] <- cross_results[length(pairs) - p_p + j]
  }
  names(Xty) <- all_cols

  # --- Phase 4: direct solve. ----------------------------------------
  beta_hat <- drop(.dsvert_solve_identifiable(
    XtX, Xty,
    context = "The protected Gaussian/GLS normal equations",
    reason = "rank_deficient_design_matrix",
    symmetric = TRUE))
  names(beta_hat) <- all_cols

  # Rename int_col -> (Intercept) for downstream callers.
  nm_out <- names(beta_hat)
  nm_out[nm_out == int_col] <- "(Intercept)"
  names(beta_hat) <- nm_out

  # Un-standardization: server divided each centered column X~_j by its
  # L2 norm s_j before Gram. Model: y = beta_0 + Sum beta_j x_j produces a
  # standardized system where beta_std_j = beta_raw_j x s_j. Unscale:
  # beta_raw_j = beta_std_j / s_j. Same for XtX and Xty to return raw-basis
  # quantities for downstream sigma^2, SE consumers:
  #   XtX_raw[i,j] = XtX_std[i,j] / (s_i x s_j)
  #   Xty_raw[j]   = Xty_std[j]   / s_j
  #   yty unchanged (y was not standardized).
  # The server returns l2_scales per column; the intercept column lives
  # on the outcome server only, so its scale is in local_y$l2_scales.
  l2_y <- as.list(local_y$l2_scales %||% list())
  l2_p <- as.list(local_p$l2_scales %||% list())
  # Build a name -> scale map, with "(Intercept)" mapped from int_col.
  scales_raw <- c(l2_y, l2_p)
  nm_beta <- names(beta_hat)
  scale_per_coef <- setNames(rep(1.0, length(nm_beta)), nm_beta)
  for (nm in nm_beta) {
    key <- if (nm == "(Intercept)") int_col else nm
    sj <- scales_raw[[key]]
    if (!is.null(sj) && length(sj) == 1L && is.finite(as.numeric(sj)) &&
        as.numeric(sj) > 0) {
      scale_per_coef[nm] <- as.numeric(sj)
    }
  }
  # Unscale beta.
  beta_unscaled <- beta_hat / scale_per_coef[nm_beta]
  # Unscale XtX and Xty to raw basis.
  S_inv <- 1.0 / scale_per_coef
  raw_nm <- rownames(XtX)
  raw_nm[raw_nm == int_col] <- "(Intercept)"
  rownames(XtX) <- colnames(XtX) <- raw_nm
  names(Xty) <- raw_nm
  S_inv_vec <- S_inv[raw_nm]
  XtX_raw <- XtX * outer(S_inv_vec, S_inv_vec)
  Xty_raw <- Xty * S_inv_vec

  # De-scale by 1/c^2 (SNR-boost legacy path; sc=1.0 default is a no-op).
  inv_sc2 <- 1.0 / (sc * sc)
  list(coefficients = beta_unscaled,
       XtX = XtX_raw * inv_sc2,
       Xty = Xty_raw * inv_sc2,
       yty = as.numeric(local_y$yty) * inv_sc2,
       n = n,
       share_scale = sc,
       standardize_applied = isTRUE(standardize),
       l2_scales_applied = scale_per_coef)
}

# Null-coalescing helper (base R has no %||%)
`%||%` <- function(a, b) if (is.null(a)) b else a
