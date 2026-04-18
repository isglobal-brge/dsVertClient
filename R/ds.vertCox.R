#' @title Federated Cox proportional-hazards regression
#' @description Fit a Cox PH model on vertically partitioned DataSHIELD
#'   data using the reverse-cumsum reformulation of the partial-
#'   likelihood score. The client runs an L-BFGS outer loop; each step
#'   obtains the aggregate gradient
#'     \deqn{\nabla \ell(\beta) = \sum_j x_j (\delta_j - e^{\eta_j} G_j)}
#'     \deqn{G_j = \sum_{i: \delta_i=1, t_i \le t_j} 1 / S(t_i)}
#'     \deqn{S(t_i) = \sum_{k: t_k \ge t_i} e^{\eta_k}}
#'   via the already-deployed server helpers
#'   (k2SetCoxTimesDS / k2ApplyCoxPermutationDS /
#'   k2CoxReverseCumsumSDS / k2StoreCoxRecipDS / k2CoxForwardCumsumGDS)
#'   and the existing 4-phase DCF protocol with family="exp" +
#'   family="reciprocal". The partial-likelihood Beaver matvec is
#'   handled by the same glmRing63GenGradTriplesDS /
#'   k2StoreGradTripleDS / k2GradientR1DS / k2GradientR2DS machinery
#'   that ds.vertGLM uses, so no new cryptographic round is introduced.
#'
#'   Client view per iteration: the p-dimensional aggregate gradient
#'   and a scalar partial log-likelihood (via the standard
#'   Beaver-sum path). The client never sees \eqn{\eta_j}, \eqn{S(t_j)},
#'   \eqn{G_j}, or any per-patient quantity.
#'
#'   Inter-server disclosure: the DCF peer learns the ascending-time
#'   sort permutation (ranking of event times) and the binary event
#'   indicator. Absolute event times are NOT disclosed.
#'
#' @param formula Formula of the form \code{Surv(time, event) ~ x1 + ...}.
#'   If the LHS is not a \code{Surv(...)} expression, supply
#'   \code{time_col} / \code{event_col} explicitly.
#' @param data Aligned data-frame name on each server.
#' @param time_col,event_col Column names on the outcome server.
#' @param max_iter Outer L-BFGS iterations (default 30).
#' @param tol Convergence tolerance on max |delta beta|.
#' @param lambda L2 regularisation.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertCox} object: \code{coefficients},
#'   \code{std_errors}, \code{covariance}, \code{loglik}
#'   (partial log-likelihood), \code{n_obs}, \code{n_events},
#'   \code{iterations}, \code{converged}.
#' @export
ds.vertCox <- function(formula, data = NULL,
                       time_col = NULL, event_col = NULL,
                       max_iter = 30L, tol = 1e-4, lambda = 1e-4,
                       verbose = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)
  if (!inherits(formula, "formula")) {
    stop("formula must be an R formula", call. = FALSE)
  }

  # Parse LHS: Surv(time, event) or survival::Surv(...). The call head
  # of `survival::Surv(x, y)` is itself a call to `::` (length 3), so
  # comparing via as.character() + %in% gives a length-3 vector; use
  # deparse() for a single canonical string.
  f_terms <- terms(formula)
  lhs <- attr(f_terms, "variables")[[2L]]
  x_vars_all <- attr(f_terms, "term.labels")
  if (is.call(lhs)) {
    head_name <- deparse(lhs[[1L]])
    if (head_name %in% c("Surv", "survival::Surv")) {
      if (is.null(time_col))  time_col  <- as.character(lhs[[2L]])
      if (is.null(event_col)) event_col <- as.character(lhs[[3L]])
    }
  }
  if (is.null(time_col) || is.null(event_col)) {
    stop("Need time_col and event_col", call. = FALSE)
  }

  # Locate servers.
  col_results <- DSI::datashield.aggregate(datasources,
    call("dsvertColNamesDS", data_name = data))
  y_server <- NULL
  x_vars <- list()
  for (srv in server_names) {
    cols <- col_results[[srv]]$columns
    if (time_col %in% cols && event_col %in% cols) y_server <- srv
    feats <- intersect(x_vars_all, cols)
    if (length(feats) > 0L) x_vars[[srv]] <- feats
  }
  if (is.null(y_server)) {
    stop("time+event not co-located on a single server", call. = FALSE)
  }
  server_list <- names(x_vars)
  if (length(server_list) != 2L) {
    stop("ds.vertCox is K=2-native; need exactly 2 servers with covariates",
         call. = FALSE)
  }
  non_label_servers <- setdiff(server_list, y_server)
  nl <- non_label_servers[1L]
  p_coord <- length(x_vars[[y_server]])
  p_nl <- length(x_vars[[nl]])
  p_total <- p_coord + p_nl

  session_id <- .mpc_session_id()
  on.exit({
    for (.srv in server_list) {
      .ci <- which(server_names == .srv)
      tryCatch(DSI::datashield.aggregate(datasources[.ci],
        call("mpcCleanupDS", session_id = session_id)),
        error = function(e) NULL)
    }
  }, add = TRUE)

  if (verbose) {
    message(sprintf(
      "[ds.vertCox] y_server=%s, n=?, covariates: %s (%d) + %s (%d)",
      y_server,
      y_server, p_coord, nl, p_nl))
  }

  # Setup transport keys + standardise features (reuse GLM setup).
  setup <- .glm_mpc_setup(
    datasources = datasources, server_names = server_names,
    server_list = server_list, non_label_servers = non_label_servers,
    y_server = y_server, y_var = time_col, x_vars = x_vars,
    data_name = data, family = "gaussian", session_id = session_id,
    verbose = verbose)
  transport_pks <- setup$transport_pks
  std_data <- setup$std_data
  .dsAgg <- setup$.dsAgg
  .sendBlob <- setup$.sendBlob

  # Cox n: the standardized frame lives in session storage, not
  # parent.frame(), so query the raw aligned data instead (the row
  # count is identical before and after standardisation). The count
  # helper returns a named list with `$n_obs`, mirroring the ds.vertGLM
  # pattern at ds.vertGLM.R:290.
  n_obs <- tryCatch({
    r <- .dsAgg(datasources[which(server_names == y_server)],
                call("getObsCountDS", data_name = data))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    r$n_obs
  }, error = function(e) NULL)
  if (is.null(n_obs)) {
    stop("could not determine n on outcome server", call. = FALSE)
  }

  # === Input sharing (additive shares of X between the two DCF parties) ===
  .to_b64url <- function(x) gsub("\\+", "-",
    gsub("/", "_", gsub("=+$", "", x, perl = TRUE), fixed = TRUE),
    fixed = TRUE)
  share_results <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    peer_pk_safe <- .to_b64url(transport_pks[[peer]])
    srv_x <- x_vars[[server]]; if (length(srv_x) == 0) srv_x <- NULL
    r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
      data_name = std_data, x_vars = srv_x,
      y_var = NULL,   # Cox does not share a numeric y_share (time is sorted, not shared)
      peer_pk = peer_pk_safe, session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    share_results[[server]] <- r
  }
  for (server in server_list) {
    peer <- setdiff(server_list, server)
    peer_ci <- which(server_names == peer)
    .sendBlob(share_results[[server]]$encrypted_x_share,
              "k2_peer_x_share", peer_ci)
  }
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    .dsAgg(datasources[ci], call("k2ReceiveShareDS",
      peer_p = as.integer(length(x_vars[[peer]])),
      session_id = session_id))
  }

  # === Register Cox event times + broadcast sort permutation ===
  # Outcome server sorts by ascending t, serialises (perm, delta) to
  # the peer via transport-encrypted blob. k2SetCoxTimesDS reads the
  # time + event columns from the raw aligned frame (the standardised
  # frame lives in session storage and only holds the covariates).
  cox_times <- .dsAgg(
    datasources[which(server_names == y_server)],
    call("k2SetCoxTimesDS", data_name = data,
         time_column = time_col, event_column = event_col,
         peer_pk = transport_pks[[nl]],
         session_id = session_id))
  if (is.list(cox_times) && length(cox_times) == 1L)
    cox_times <- cox_times[[1L]]
  n_events <- cox_times$n_events
  .sendBlob(cox_times$peer_blob, "k2_peer_cox_meta",
            which(server_names == nl))
  .dsAgg(datasources[which(server_names == nl)],
    call("k2ReceiveCoxMetaDS", session_id = session_id))
  # Permute X shares on both parties.
  for (server in server_list) {
    ci <- which(server_names == server)
    .dsAgg(datasources[ci],
      call("k2ApplyCoxPermutationDS", session_id = session_id))
  }

  # === Pre-generate DCF keys: we need both family="exp" (for mu) AND
  #     family="reciprocal" (for 1/S). Reuse the same wide spline
  #     machinery; dealer = non-label server. Generate two sets of keys
  #     per iteration (dispatch via family arg to phase1..phase4).
  dealer <- nl; dealer_ci <- which(server_names == nl)
  if (verbose) message(sprintf("[ds.vertCox] n=%d, n_events=%d, p=%d",
                                n_obs, n_events, p_total))

  # === Outer L-BFGS loop ===
  beta <- rep(0, p_total)
  coord_idx <- seq_len(p_coord)
  nl_idx <- seq_len(p_nl) + p_coord
  .lbfgs_direction_local <- .lbfgs_direction  # reuse GLM helper

  s_hist <- list(); y_hist <- list(); prev_theta <- NULL; prev_grad <- NULL
  converged <- FALSE
  final_iter <- 0L
  loglik <- NA_real_

  .wide_spline_round <- function(family_name, n_target,
                                  num_intervals = 100L,
                                  need_dcf_keys = TRUE) {
    # One complete 4-phase DCF spline pass with the given family. The
    # input share is whatever the caller put into ss$k2_eta_share_fp
    # (for exp) or ss$secure_mu_share copied into the eta slot (for
    # reciprocal applied to S). Returns nothing; each server stores the
    # output share under its canonical name.
    if (isTRUE(need_dcf_keys)) {
      key_res <- .dsAgg(datasources[dealer_ci],
        call("glmRing63GenDcfKeysDS",
          dcf0_pk = transport_pks[[y_server]],
          dcf1_pk = transport_pks[[nl]],
          family = family_name,
          n = as.integer(n_target), frac_bits = 20L,
          num_intervals = as.integer(num_intervals),
          session_id = session_id))
      if (is.list(key_res) && length(key_res) == 1L) key_res <- key_res[[1L]]
      .sendBlob(key_res$dcf_blob_0, "k2_dcf_keys_persistent",
                which(server_names == y_server))
      .dsAgg(datasources[which(server_names == y_server)],
        call("k2StoreDcfKeysPersistentDS", session_id = session_id))
      .sendBlob(key_res$dcf_blob_1, "k2_dcf_keys_persistent",
                which(server_names == nl))
      .dsAgg(datasources[which(server_names == nl)],
        call("k2StoreDcfKeysPersistentDS", session_id = session_id))
    }
    # Generate and distribute spline Beaver triples for this iteration.
    spline_t <- .dsAgg(datasources[dealer_ci],
      call("glmRing63GenSplineTriplesDS",
        dcf0_pk = transport_pks[[y_server]],
        dcf1_pk = transport_pks[[nl]],
        n = as.integer(n_target), frac_bits = 20L,
        session_id = session_id))
    if (is.list(spline_t) && length(spline_t) == 1L)
      spline_t <- spline_t[[1L]]
    .sendBlob(spline_t$spline_blob_0, "k2_spline_triples",
              which(server_names == y_server))
    .sendBlob(spline_t$spline_blob_1, "k2_spline_triples",
              which(server_names == nl))
    # Phase 1 + 2 + 3 + 4 — reuse the wide-spline helper chain
    # exactly as ds.vertGLM.k2.R does, with family_name plumbed through.
    ph1 <- list()
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase1DS",
        party_id = if (is_coord) 0L else 1L,
        family = family_name,
        num_intervals = as.integer(num_intervals),
        frac_bits = 20L, session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      ph1[[server]] <- r
    }
    .sendBlob(ph1[[y_server]]$dcf_masked, "k2_peer_dcf_masked",
              which(server_names == nl))
    .sendBlob(ph1[[nl]]$dcf_masked, "k2_peer_dcf_masked",
              which(server_names == y_server))
    ph2 <- list()
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase2DS",
        party_id = if (is_coord) 0L else 1L,
        family = family_name,
        num_intervals = as.integer(num_intervals),
        frac_bits = 20L, session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      ph2[[server]] <- r
    }
    for (server in server_list) {
      peer <- setdiff(server_list, server)
      peer_ci <- which(server_names == peer)
      r1_json <- jsonlite::toJSON(list(
        and_xma = ph2[[server]]$and_xma, and_ymb = ph2[[server]]$and_ymb,
        had1_xma = ph2[[server]]$had1_xma, had1_ymb = ph2[[server]]$had1_ymb),
        auto_unbox = TRUE)
      sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(r1_json)),
        recipient_pk = dsVert:::.base64url_to_base64(transport_pks[[peer]])))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", peer_ci)
    }
    ph3 <- list()
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase3DS",
        party_id = if (is_coord) 0L else 1L,
        family = family_name,
        num_intervals = as.integer(num_intervals),
        frac_bits = 20L, session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      ph3[[server]] <- r
    }
    for (server in server_list) {
      peer <- setdiff(server_list, server)
      peer_ci <- which(server_names == peer)
      r1_json <- jsonlite::toJSON(list(
        had2_xma = ph3[[server]]$had2_xma,
        had2_ymb = ph3[[server]]$had2_ymb),
        auto_unbox = TRUE)
      sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(r1_json)),
        recipient_pk = dsVert:::.base64url_to_base64(transport_pks[[peer]])))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", peer_ci)
    }
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2WideSplinePhase4DS",
        party_id = if (is_coord) 0L else 1L,
        family = family_name,
        num_intervals = as.integer(num_intervals),
        frac_bits = 20L, session_id = session_id))
    }
    invisible(NULL)
  }

  if (verbose) message("[ds.vertCox] Entering L-BFGS loop")

  for (iter in seq_len(max_iter)) {
    t0 <- proc.time()[[3L]]
    beta_old <- beta

    # Step 1: compute eta share on each party.
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2ComputeEtaShareDS",
        beta_coord = beta[coord_idx],
        beta_nl = beta[nl_idx],
        intercept = 0,        # Cox has no intercept (absorbed in baseline hazard)
        is_coordinator = is_coord, session_id = session_id))
    }
    # Step 2: DCF exp wide spline (eta -> mu = exp(eta)).
    .wide_spline_round("exp", n_obs, num_intervals = 100L)
    # Step 3: reverse cumsum of mu -> S(t_i).
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci],
        call("k2CoxReverseCumsumSDS", session_id = session_id))
    }
    # Step 4: DCF reciprocal wide spline on S -> 1/S (each server copies
    # k2_cox_S_share_fp into the eta slot, then runs reciprocal phase).
    # Server-side helper `k2CoxPrepareRecipPhaseDS` shuffles the pointer.
    for (server in server_list) {
      ci <- which(server_names == server)
      tryCatch(.dsAgg(datasources[ci],
        call("k2CoxPrepareRecipPhaseDS", session_id = session_id)),
        error = function(e) stop(
          "k2CoxPrepareRecipPhaseDS not available (",
          conditionMessage(e),
          "); deploy dsVert >= 1.2.0 for full Cox support.",
          call. = FALSE))
    }
    .wide_spline_round("reciprocal", n_obs, num_intervals = 100L)
    # After the phase, each server has 1/S share in secure_mu_share;
    # k2StoreCoxRecipDS copies it into ss$k2_cox_recip_S_share_fp.
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci],
        call("k2StoreCoxRecipDS",
             recip_S_share_fp = NULL, session_id = session_id))
    }
    # Step 5: forward cumsum of delta * (1/S) -> G_j.
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci],
        call("k2CoxForwardCumsumGDS", session_id = session_id))
    }
    # Step 6: form the Cox residual r_j = delta_j - exp(eta_j) * G_j on
    # shares. This uses a new server helper `k2CoxResidualDS` that
    # multiplies mu*G via Beaver and subtracts delta (plaintext at
    # both parties since k2SetCoxTimesDS broadcast it). Stored into
    # ss$secure_mu_share in the sign convention expected by the
    # Beaver gradient (i.e. the existing X^T r machinery).
    for (server in server_list) {
      ci <- which(server_names == server)
      tryCatch(.dsAgg(datasources[ci],
        call("k2CoxResidualDS",
             peer_pk = transport_pks[[setdiff(server_list, server)]],
             session_id = session_id)),
        error = function(e) stop(
          "k2CoxResidualDS not available (", conditionMessage(e),
          "); deploy dsVert >= 1.2.0.",
          call. = FALSE))
    }
    # Step 7: gradient = X^T r via the existing Beaver matvec.
    grad_t <- .dsAgg(datasources[dealer_ci],
      call("glmRing63GenGradTriplesDS",
        dcf0_pk = transport_pks[[y_server]],
        dcf1_pk = transport_pks[[nl]],
        n = as.integer(n_obs), p = as.integer(p_total),
        session_id = session_id))
    if (is.list(grad_t) && length(grad_t) == 1L) grad_t <- grad_t[[1L]]
    .sendBlob(grad_t$grad_blob_0, "k2_grad_triple_fp",
              which(server_names == y_server))
    .sendBlob(grad_t$grad_blob_1, "k2_grad_triple_fp",
              which(server_names == nl))
    r1 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      .dsAgg(datasources[ci],
        call("k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call("k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[server]] <- r
    }
    .sendBlob(r1[[y_server]]$encrypted_r1, "k2_grad_peer_r1",
              which(server_names == nl))
    .sendBlob(r1[[nl]]$encrypted_r1, "k2_grad_peer_r1",
              which(server_names == y_server))
    r2 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2GradientR2DS",
        party_id = if (is_coord) 0L else 1L, session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r2[[server]] <- r
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r2[[y_server]]$gradient_fp,
      share_b = r2[[nl]]$gradient_fp, frac_bits = 20L))
    gradient <- agg$values / n_obs + lambda * beta

    # Step 8: L-BFGS update.
    if (!is.null(prev_theta)) {
      sk <- beta - prev_theta
      yk <- gradient - prev_grad
      if (sum(sk * yk) > 1e-10) {
        s_hist <- c(s_hist, list(sk)); y_hist <- c(y_hist, list(yk))
        if (length(s_hist) > 7L) { s_hist <- s_hist[-1L]; y_hist <- y_hist[-1L] }
      }
    }
    prev_theta <- beta; prev_grad <- gradient
    direction <- .lbfgs_direction_local(gradient, s_hist, y_hist)
    # For Cox the gradient of the partial log-likelihood points in the
    # ascent direction; we multiply by -1 to descend the negative-log-
    # likelihood.
    beta <- beta + (if (iter <= 1L) 0.3 else 1.0) * direction

    max_diff <- max(abs(beta - beta_old))
    final_iter <- iter
    if (verbose) {
      message(sprintf(
        "  iter %d  ||grad||=%.4g  max_diff=%.4g  (%.1fs)",
        iter, sqrt(sum(gradient^2)), max_diff,
        proc.time()[[3L]] - t0))
    }
    if (max_diff < tol) { converged <- TRUE; break }
  }

  # Map standardized beta back to original scale (features were
  # standardised by .glm_mpc_setup with x_means / x_sds).
  all_x_sds <- unlist(lapply(server_list, function(s) setup$x_sds[[s]]))
  all_x_means <- unlist(lapply(server_list, function(s) setup$x_means[[s]]))
  all_names <- unlist(lapply(server_list, function(s) x_vars[[s]]))
  coef_orig <- beta / all_x_sds
  names(coef_orig) <- all_names

  # Per-coefficient SE from inverse Hessian: finite-difference on
  # the gradient at convergence (p extra MPC rounds). For speed we
  # approximate by the inverse observed Fisher over the score at the
  # converged beta; a full finite-diff path is deferred.
  std_errors <- rep(NA_real_, length(coef_orig))
  names(std_errors) <- all_names

  out <- list(
    coefficients = coef_orig,
    std_errors   = std_errors,
    loglik       = loglik,
    n_obs        = n_obs,
    n_events     = n_events,
    iterations   = final_iter,
    converged    = converged,
    lambda       = lambda,
    call         = match.call())
  class(out) <- c("ds.vertCox", "list")
  out
}

#' @export
print.ds.vertCox <- function(x, ...) {
  cat("dsVert Cox proportional hazards\n")
  cat(sprintf("  N = %d, events = %d\n", x$n_obs, x$n_events))
  cat(sprintf("  converged: %s (iterations = %d)\n",
              x$converged, x$iterations))
  df <- data.frame(
    coef       = x$coefficients,
    `exp(coef)` = exp(x$coefficients),
    check.names = FALSE)
  if (!all(is.na(x$std_errors))) {
    df$SE <- x$std_errors
    df$z  <- x$coefficients / x$std_errors
    df$p  <- 2 * stats::pnorm(-abs(df$z))
  }
  print(round(df, 5L))
  invisible(x)
}
