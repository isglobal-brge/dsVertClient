#' @title Quarantined linear-mixed-model compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and computes or returns no mixed model, variance component,
#'   cluster statistic, or diagnostic. Retained LMM implementation code after
#'   the gate is unreachable through this public frontdoor and carries no
#'   disclosure, DP, accuracy, or availability claim.
#' @details Promotion requires contribution-bounded cluster statistics,
#'   private cluster handling, certified ML/REML and random-effects semantics,
#'   covariance and identifiability certificates, and independent multi-host
#'   validation.
#' @param formula,data,cluster_col,random_slopes,reml,max_iter,inner_iter,tol,exact_cross_server,sigma_b2_override,ring,verbose,datasources
#'   Retained compatibility arguments. They are not evaluated because the
#'   public frontdoor fails locally.
#' @return No fitted object. The function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertLMM <- function(formula, data = NULL, cluster_col,
                       random_slopes = NULL,
                       reml = TRUE, max_iter = 30L, inner_iter = 50L,
                       tol = 1e-4,
                       exact_cross_server = TRUE,
                       sigma_b2_override = NULL,
                       ring = c("ring63", "ring127"),
                       verbose = TRUE, datasources = NULL) {
  .dsvert_block_retired_remote_route("lmm")
  # Task #121: Ring127 LMM migration. ring="ring127" switches the
  # Beaver vecmul pipeline (LocalGram shares, R1, R2, Aggregate) to
  # Ring127 fracBits=50 to drive the per-Gram-entry noise from
  # rel~1e-4 down to rel~1e-8, closing the X4 rel STRICT gap.
  ring <- match.arg(ring)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)
  if (missing(cluster_col) || !is.character(cluster_col) ||
      length(cluster_col) != 1L) {
    stop("cluster_col: single character column name required",
         call. = FALSE)
  }

  # Locate outcome + cluster servers.
  y_var <- .ds_gee_extract_lhs(formula)
  column_results <- .dsvert_lmm_aggregate_strict(
    datasources, call(name = "dsvertColNamesDS", data_name = data))
  server_columns <- lapply(column_results, function(value) value$columns)
  locate_server <- function(variable) {
    present <- server_names[vapply(
      server_columns, function(columns) variable %in% columns, logical(1L))]
    if (length(present) > 1L) {
      stop("variable '", variable,
           "' is present on more than one server; choose an unambiguous vertical partition",
           call. = FALSE)
    }
    if (length(present) == 1L) present[[1L]] else NULL
  }
  y_srv <- locate_server(y_var)
  clust_srv <- locate_server(cluster_col)
  if (is.null(y_srv)) stop("y '", y_var, "' not located", call. = FALSE)
  if (is.null(clust_srv)) stop("cluster_col '", cluster_col,
                                "' not located", call. = FALSE)
  if (clust_srv != y_srv) {
    stop("cluster_col must live on the outcome server (y is on '",
         y_srv, "', cluster_col is on '", clust_srv, "'). ",
         "This LMM frontdoor requires the cluster column on the outcome ",
         "server.",
         call. = FALSE)
  }

  # Ask the outcome server for the cluster-size vector (aggregate).
  clust_info <- tryCatch(
    .dsvert_lmm_aggregate_strict(
      datasources[which(server_names == y_srv)],
      call(name = "dsvertClusterSizesDS", data_name = data,
           cluster_col = cluster_col)),
    error = function(e) {
      .dsvert_lmm_transport_fallback(e)
      stop("dsvertClusterSizesDS not available (", conditionMessage(e),
           "); deploy dsVert >= 1.2.0 for LMM support.",
           call. = FALSE)
    })
  if (is.list(clust_info) && length(clust_info) == 1L)
    clust_info <- clust_info[[1]]
  n_per_cluster <- as.integer(clust_info$sizes)
  n_clusters <- length(n_per_cluster)
  n_total <- sum(n_per_cluster)

  if (verbose) {
    message(sprintf(
      "[ds.vertLMM] %d clusters, n_total=%d (sizes: median=%.1f, max=%d)",
      n_clusters, n_total,
      as.numeric(stats::median(n_per_cluster)),
      as.integer(max(n_per_cluster))))
  }

  # Initial unweighted fit to prime beta and sigma^2.
  fit0 <- ds.vertGLM(formula = formula, data = data, family = "gaussian",
                     max_iter = inner_iter, tol = tol,
                     verbose = FALSE, datasources = datasources)
  if (is.null(fit0$covariance)) {
    stop("inner ds.vertGLM must expose Cov(beta); refit with ",
         "dsVert >= 8bb7902.", call. = FALSE)
  }
  sigma2 <- fit0$deviance / max(n_total - length(fit0$coefficients), 1L)
  if (!is.finite(sigma2) || sigma2 <= 0) {
    .dsvert_stop_non_identifiable(
      "The LMM residual variance is non-positive; variance components are not identifiable.",
      reason = "degenerate_lmm_residual_variance")
  }

  # Discover which predictors live on the outcome server so we can
  # pass the right slice of betahat and fitted-value predictors to the
  # server-side residual helper. Cross-server predictors get absorbed
  # into the intercept correction term below.
  y_srv_cols <- server_columns[[y_srv]]
  x_all <- attr(terms(formula), "term.labels")
  x_local_ysrv <- intersect(x_all, y_srv_cols)
  x_remote <- setdiff(x_all, x_local_ysrv)
  if (length(x_remote) > 0L && verbose) {
    message("[ds.vertLMM] non-outcome-server predictors (",
            paste(x_remote, collapse = ","),
            ") are absorbed into the intercept for residual SS; ",
            "ICC estimate is on the outcome-server projection only")
  }

  # Exact cross-server path: use the shipped Beaver-based residual
  # pipeline (dsvertLMMPeerFittedShareDS + dsvertLMMCoordResidualShareDS
  # + dsvertLMMPeerResidualFinaliseDS + k2BeaverVecmul chain +
  # dsvertLMMPerClusterSumDS) to compute exact per-cluster SS without
  # the intercept-absorption approximation. If this requested backend is
  # unavailable it fails closed; the approximation is only an explicit
  # exact_cross_server = FALSE estimand.
  peer_servers <- setdiff(server_names, y_srv)
  # K-dispatch. This closed-form solver is the K=2 backend (exact GLS via the
  # Beaver Gram chain). K>=3 fits the same one-way random-intercept estimand
  # with the variance-ratio profile algorithm; delegate to it so a single
  # ds.vertLMM entry point serves every K (mirroring ds.vertGLM / ds.vertCox).
  # K=1 has no peer and cannot fit a cross-server model.
  if (length(peer_servers) == 0L)
    stop("ds.vertLMM requires at least one non-outcome (peer) server (K>=2).",
         call. = FALSE)
  if (length(peer_servers) >= 2L) {
    if (!is.null(random_slopes))
      stop("random_slopes (q>1 random effects) are supported only at K=2; the ",
           "K>=3 path fits a random-intercept model. Drop random_slopes for ",
           "K>=3.", call. = FALSE)
    return(.ds_vertLMM_k3_impl(formula = formula, data = data,
                               cluster_col = cluster_col,
                               tol = tol, max_outer = max_iter,
                               ring = "ring127", verbose = verbose,
                               datasources = datasources))
  }
  peer_srv <- peer_servers[1L]

  get_cluster_resids_exact <- function(beta_hat, session_id_active,
                                         transport_pks_active) {
    # 1. Peer computes fitted_peer + splits shares.
    x_remote_vars <- setdiff(names(beta_hat), c(x_local_ysrv, "(Intercept)"))
    b_remote <- as.numeric(beta_hat[x_remote_vars])
    r <- .dsvert_lmm_aggregate_strict(
      datasources[which(server_names == peer_srv)],
      call(name = "dsvertLMMPeerFittedShareDS",
           data_name = data, x_names = x_remote_vars,
           betahat = b_remote,
           peer_pk = transport_pks_active[[y_srv]],
           session_id = session_id_active))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    # 2. Relay peer's blob to outcome server.
    .dsvert_store_blob(
      r$peer_blob, "k2_lmm_exact_peer_blob",
      datasources[which(server_names == y_srv)], session_id_active)
    # 3. Outcome: compute r_share_0.
    .dsvert_lmm_aggregate_strict(
      datasources[which(server_names == y_srv)],
      call(name = "dsvertLMMCoordResidualShareDS",
           data_name = data, y_var = y_var,
           x_names = x_local_ysrv,
           betahat_local = as.numeric(beta_hat[x_local_ysrv]),
           intercept = as.numeric(beta_hat["(Intercept)"]),
           session_id = session_id_active))
    # 4. Peer finalises its share.
    .dsvert_lmm_aggregate_strict(
      datasources[which(server_names == peer_srv)],
      call(name = "dsvertLMMPeerResidualFinaliseDS",
           session_id = session_id_active))
    # (The subsequent Beaver vecmul + per-cluster sum would follow;
    # for the outer REML update we only need per-cluster sums which
    # we get via the approximate helper for now. Full Beaver vecmul
    # integration in the outer loop is shipped in ds.vertLMM v3.)
    list(exact = TRUE)
  }

  get_cluster_resids <- function(beta_hat) {
    b_local <- as.numeric(beta_hat[x_local_ysrv])
    tryCatch(
      .dsvert_lmm_aggregate_strict(
        datasources[which(server_names == y_srv)],
        call(name = "dsvertClusterResidualsDS",
             data_name = data,
             y_var = y_var,
             x_names = x_local_ysrv,
             intercept = as.numeric(beta_hat["(Intercept)"]),
             betahat = b_local,
             cluster_col = cluster_col)),
      error = function(e) {
        .dsvert_lmm_transport_fallback(e)
        stop("dsvertClusterResidualsDS failure: ",
             conditionMessage(e), call. = FALSE)
      })
  }

  # Random-slopes path: if random_slopes is non-empty, fetch per-cluster
  # Z^T Z matrices from the outcome server for a q x q Woodbury inverse
  # per cluster. The random-effects design per cluster is
  #   Z_i = [1, slope_var_1, slope_var_2, ...]_{j in cluster i}
  # and the covariance matrix Omega is q x q (q = 1 + length(random_slopes)).
  q <- 1L + length(random_slopes)
  Z_info <- NULL
  if (!is.null(random_slopes) && length(random_slopes) > 0L) {
    missing_slopes <- setdiff(random_slopes, y_srv_cols)
    if (length(missing_slopes) > 0L) {
      stop("random_slopes not on outcome server: ",
           paste(missing_slopes, collapse = ","),
           ". Cross-server random-slope columns are not supported by ",
           "this LMM frontdoor.",
           call. = FALSE)
    }
    Z_info <- tryCatch(
      .dsvert_lmm_aggregate_strict(
        datasources[which(server_names == y_srv)],
        call(name = "dsvertClusterZtZDS",
             data_name = data,
             cluster_col = cluster_col,
             slope_columns = random_slopes)),
      error = function(e) {
        .dsvert_lmm_transport_fallback(e)
        stop("dsvertClusterZtZDS unavailable: ",
             conditionMessage(e), call. = FALSE)
      })
    if (is.list(Z_info) && length(Z_info) == 1L) Z_info <- Z_info[[1L]]
    if (verbose) {
      message(sprintf("[ds.vertLMM] random effects: intercept + %s (q=%d)",
                       paste(random_slopes, collapse = "+"), q))
    }
  }

  # ==== Exact cross-server residual pipeline orchestration ====
  # Runs the full Beaver-based per-cluster r / r^2 aggregation; returns
  # (rsum_per_cluster, rss_total) that we plug into the REML update
  # INSTEAD of the intercept-absorption approximation. Requires all the
  # new dsvertLMM* helpers (commits ad0df6c + this one).
  get_cluster_resids_full_exact <- function(beta_hat) {
    if (is.null(peer_srv)) {
      return(NULL)  # single-server case falls back to approximate path
    }
    # CRITICAL: use actual row count on outcome server, not
    # sum(n_per_cluster) which may be smaller due to privacy
    # suppression of small clusters. All Beaver ops need n_actual.
    n_actual <- tryCatch({
      r <- .dsvert_lmm_aggregate_strict(datasources[which(server_names == y_srv)],
        call(name = "getObsCountDS", data_name = data))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      as.integer(r$n_obs)
    }, error = function(e) .dsvert_lmm_transport_fallback(e))
    if (is.null(n_actual) || n_actual <= 0) return(NULL)
    # Reuse the already-live session (ds.vertGLM ran with fit-time session
    # but cleaned up; we open a fresh one for the exact pipeline).
    sess <- .mpc_session_id()
    y_srv_ci <- which(server_names == y_srv)
    peer_ci <- which(server_names == peer_srv)
    residual_conns <- datasources[c(y_srv_ci, peer_ci)]
    on.exit({
      for (.ci in c(y_srv_ci, peer_ci)) {
        .dsvert_cleanup_best_effort(
          datasources[.ci], call(name = "mpcCleanupDS", session_id = sess))
      }
    }, add = TRUE)
    # Transport keys.
    pks <- .dsvert_setup_peer_transport(datasources, server_names,
                                        c(y_srv, peer_srv), sess)
    x_remote_vars <- setdiff(names(beta_hat),
                              c(x_local_ysrv, "(Intercept)"))
    b_remote <- as.numeric(beta_hat[x_remote_vars])
    # 1. Peer fitted share.
    r1 <- tryCatch(.dsvert_lmm_aggregate_strict(datasources[peer_ci],
      call(name = "dsvertLMMPeerFittedShareDS",
           data_name = data, x_names = x_remote_vars,
           betahat = b_remote, peer_pk = pks[[y_srv]],
           session_id = sess)),
      error = function(e) .dsvert_lmm_transport_fallback(e))
    if (is.null(r1)) return(NULL)
    if (is.list(r1) && length(r1) == 1L) r1 <- r1[[1L]]
    .dsvert_store_blob(
      r1$peer_blob, "k2_lmm_exact_peer_blob", datasources[y_srv_ci], sess)
    # 2. Coord residual share.
    tryCatch(.dsvert_lmm_aggregate_strict(datasources[y_srv_ci],
      call(name = "dsvertLMMCoordResidualShareDS",
           data_name = data, y_var = y_var,
           x_names = x_local_ysrv,
           betahat_local = as.numeric(beta_hat[x_local_ysrv]),
           intercept = as.numeric(beta_hat["(Intercept)"]),
           session_id = sess)),
      error = function(e) {
        .dsvert_lmm_transport_fallback(e)
        message("[LMM exact] coord: ", conditionMessage(e))
        NULL
      })
    # 3. Peer residual finalise. Needs k2_x_n; we attempt to set it via a
    # no-op that also caches length. Fallback: skip if unavailable.
    tryCatch({
      # Set peer's k2_x_n from the row count
      n_check <- .dsvert_lmm_aggregate_strict(datasources[peer_ci],
        call(name = "getObsCountDS", data_name = data))[[1L]]$n_obs
      # Manually stamp via a helper; if not available, finalise may
      # fail but we can still proceed with partial info.
    }, error = function(e) .dsvert_lmm_transport_fallback(e))
    fin_ok <- tryCatch({
      .dsvert_lmm_aggregate_strict(datasources[peer_ci],
        call(name = "dsvertLMMPeerResidualFinaliseDS",
             n = as.integer(n_actual),
             session_id = sess))
      TRUE
    }, error = function(e) {
      .dsvert_lmm_transport_fallback(e)
      message("[LMM exact] peer finalise: ", conditionMessage(e))
      FALSE
    })
    if (!fin_ok) return(NULL)
    # 4. Cluster-ID broadcast.
    cb <- tryCatch(.dsvert_lmm_aggregate_strict(datasources[y_srv_ci],
      call(name = "dsvertLMMBroadcastClusterIDsDS",
           data_name = data, cluster_col = cluster_col,
           peer_pk = pks[[peer_srv]], session_id = sess)),
      error = function(e) .dsvert_lmm_transport_fallback(e))
    if (is.null(cb)) return(NULL)
    if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
    .dsvert_store_blob(
      cb$peer_blob, "k2_lmm_cluster_ids_blob", datasources[peer_ci], sess)
    .dsvert_lmm_aggregate_strict(datasources[peer_ci],
      call(name = "dsvertLMMReceiveClusterIDsDS", session_id = sess))
    # 5. Per-cluster rsum (both parties) + aggregate client-side.
    rs <- .dsvert_lmm_aggregate_strict(residual_conns,
      call(name = "dsvertLMMPerClusterSumDS",
           share_key = "k2_lmm_exact_r_share", session_id = sess))
    rs_y <- rs[[y_srv]]
    rs_p <- rs[[peer_srv]]
    K <- length(rs_y$per_cluster_fp)
    rsum_cluster <- numeric(K)
    for (ck in seq_len(K)) {
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = rs_y$per_cluster_fp[[ck]],
        share_b = rs_p$per_cluster_fp[[ck]], frac_bits = 20L))
      rsum_cluster[ck] <- as.numeric(agg$values[1L])
    }
    # 6. Beaver vecmul r x r -> r^2 share.
    send_blob <- function(blob, key, conn_idx) {
      .dsvert_store_blob(blob, key, datasources[conn_idx], sess)
    }
    ds_agg <- function(ds, expr) .dsvert_lmm_aggregate_strict(ds, expr)
    .ot_beaver_prepare_vecmul(
      datasources = datasources,
      party_conns = c(y_srv_ci, peer_ci),
      party_names = c(y_srv, peer_srv),
      transport_pks = pks,
      session_id = sess,
      n = n_actual,
      ring = 63L,
      .dsAgg = ds_agg,
      .sendBlob = send_blob)
    r1 <- .dsvert_fanout_by_site(
      residual_conns,
      stats::setNames(list(
        call(name = "k2BeaverVecmulR1DS",
             peer_pk = pks[[peer_srv]],
             x_key = "k2_lmm_exact_r_share",
             y_key = "k2_lmm_exact_r_share",
             n = as.integer(n_actual), session_id = sess, frac_bits = 20L),
        call(name = "k2BeaverVecmulR1DS",
             peer_pk = pks[[y_srv]],
             x_key = "k2_lmm_exact_r_share",
             y_key = "k2_lmm_exact_r_share",
             n = as.integer(n_actual), session_id = sess, frac_bits = 20L)),
        c(y_srv, peer_srv)),
      operation = "LMM residual vecmul R1")
    r1a <- r1[[y_srv]]
    r1b <- r1[[peer_srv]]
    .dsvert_store_typed_blob(
      r1a$peer_blob, r1a$peer_transfer, datasources[peer_ci], sess,
      producer_conn = datasources[y_srv_ci])
    .dsvert_store_typed_blob(
      r1b$peer_blob, r1b$peer_transfer, datasources[y_srv_ci], sess,
      producer_conn = datasources[peer_ci])
    .dsvert_fanout_by_site(
      residual_conns,
      stats::setNames(lapply(c(y_srv, peer_srv), function(srv) {
        call(name = "k2BeaverVecmulR2DS",
             is_party0 = (srv == y_srv),
             peer_name = setdiff(c(y_srv, peer_srv), srv),
             x_key = "k2_lmm_exact_r_share",
             y_key = "k2_lmm_exact_r_share",
             output_key = "k2_lmm_exact_r2_share",
             n = as.integer(n_actual), session_id = sess,
             frac_bits = 20L)
      }), c(y_srv, peer_srv)),
      operation = "LMM residual vecmul R2")
    # 7. Global sum r^2 (both parties) -> aggregate.
    gs <- .dsvert_lmm_aggregate_strict(residual_conns,
      call(name = "dsvertLMMGlobalSumDS",
           share_key = "k2_lmm_exact_r2_share", session_id = sess))
    gs_y <- gs[[y_srv]]
    gs_p <- gs[[peer_srv]]
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = gs_y$sum_fp, share_b = gs_p$sum_fp, frac_bits = 20L))
    rss_total <- as.numeric(agg$values[1L])
    # 7b. EXACT per-cluster rss: reuse dsvertLMMPerClusterSumDS on the
    # r^2 share (same helper we just used on r). The two DCF parties
    # already hold the additive shares of r^2 under k2_lmm_exact_r2_share
    # after step 6's Beaver vecmul; per-cluster sum is a linear op that
    # preserves additive sharing, so we aggregate client-side to get the
    # exact per-cluster r^2 vector that the REML profile likelihood needs
    # to detect a weak sigma_b^2 signal. Closes the LMM intercept-
    # absorption approximation shipped in v2.
    r2s <- .dsvert_lmm_aggregate_strict(residual_conns,
      call(name = "dsvertLMMPerClusterSumDS",
           share_key = "k2_lmm_exact_r2_share", session_id = sess))
    r2s_y <- r2s[[y_srv]]
    r2s_p <- r2s[[peer_srv]]
    rss_cluster <- numeric(K)
    for (ck in seq_len(K)) {
      agg2 <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = r2s_y$per_cluster_fp[[ck]],
        share_b = r2s_p$per_cluster_fp[[ck]], frac_bits = 20L))
      rss_cluster[ck] <- as.numeric(agg2$values[1L])
    }
    list(rsum_per_cluster = rsum_cluster,
         rss_per_cluster  = rss_cluster,
         rss_total = rss_total, n_per_cluster = rs_y$cluster_sizes)
  }

  # Outer REML loop: update variance components (Omega q x q + sigma^2).
  # Initialise Omega = 0.1 sigma^2 I, covariance matrix of random effects.
  Omega <- diag(q) * sigma2 * 0.1
  sigma_b2 <- Omega[1L, 1L]  # keep legacy slot for intercept variance
  converged <- FALSE
  rho_prev <- Inf
  fit <- fit0
  # HYBRID (y) Aitken acceleration attempted 2026-04-19 late but
  # empirically DESTABILIZED the outer REML loop: determinism |Delta|
  # degraded from 2.2e-5 (V2 alone) to 11.35 units across runs, a
  # catastrophic regression. Aitken extrapolation off the contraction
  # path induces oscillation; with max_iter=30 the iterate wanders.
  # Reverted to plain Picard (V2 state). Task #115 remains open for a
  # stability-preserving acceleration scheme. The sigma^2-outer-loop
  # coupling floor of 6e-4 rel (-> beta_X4 rel ~1.8e-4) is accepted as
  # the current LMM precision limit.
  for (iter in seq_len(max_iter)) {
    # The exact route and the explicitly requested approximation are distinct
    # estimands. Never substitute the latter after an exact-route failure.
    cl_exact <- NULL
    if (isTRUE(exact_cross_server) && !is.null(peer_srv) &&
        length(setdiff(attr(terms(formula), "term.labels"),
                         x_local_ysrv)) > 0L) {
      cl_exact <- tryCatch(
        get_cluster_resids_full_exact(fit$coefficients),
        error = function(e) {
          if (inherits(e, "non_identifiable")) stop(e)
          .dsvert_lmm_transport_fallback(e)
          .dsvert_stop_numeric(
            "numeric_backend_unavailable",
            paste0("The exact LMM residual backend failed: ",
                   conditionMessage(e),
                   ". No approximate residual fallback was applied."),
            reason = "lmm_exact_residual_backend_failure")
        })
      if (is.null(cl_exact)) {
        .dsvert_stop_numeric(
          "numeric_backend_unavailable",
          paste0("The exact LMM residual backend returned no result. ",
                 "No approximate residual fallback was applied."),
          reason = "lmm_exact_residual_backend_failure")
      }
    }
    if (!is.null(cl_exact)) {
      rsum <- as.numeric(cl_exact$rsum_per_cluster)
      # EXACT per-cluster rss (from the second reuse of
      # dsvertLMMPerClusterSumDS on the r^2 share). Required for the
      # REML profile likelihood to detect weak sigma_b^2 signals --
      # the rss_total/n_clusters approximation previously shipped here
      # washed out within-cluster variance heterogeneity and caused the
      # LMM FAIL (|Delta| ~ 0.07 vs lme4) on datasets with small ICC.
      if (!is.null(cl_exact$rss_per_cluster) &&
          length(cl_exact$rss_per_cluster) == n_clusters) {
        rss <- as.numeric(cl_exact$rss_per_cluster)
      } else {
        .dsvert_stop_numeric(
          "numeric_backend_unavailable",
          "The exact LMM backend did not return per-cluster residual sums of squares.",
          reason = "incomplete_lmm_exact_statistics")
      }
    } else {
      cl <- get_cluster_resids(fit$coefficients)
      if (is.list(cl) && length(cl) == 1L) cl <- cl[[1]]
      rss <- as.numeric(cl$rss_per_cluster)
      rsum <- as.numeric(cl$rsum_per_cluster)
    }
    stopifnot(length(rss) == n_clusters)

    # Variance components: sigma^2 by MoM (matches lme4 exactly for
    # balanced data); sigma_b^2 by a 1-D profile ML optimiser over
    # the per-cluster log-likelihood. The profile approach avoids the
    # MoM floor at 0 when the signal is weak (sigma_b^2 < 0.2) --
    # lme4's REML picks it up and so should we.
    #
    # Per-cluster marginal log-lik in r_i = y_i - X_i beta:
    #   log p(r_i | sigma^2, sigma_b^2) = -0.5 * [
    #     log|V_i| + r_i^T V_i^{-1} r_i + n_i * log(2 pi)
    #   ]
    #   where V_i = sigma^2 I + sigma_b^2 1 1^T and closed-form:
    #     log|V_i| = (n_i - 1) log sigma^2 + log(sigma^2 + n_i sigma_b^2)
    #     r_i^T V_i^{-1} r_i = rss_i / sigma^2
    #         - (sigma_b^2 / (sigma^2 (sigma^2 + n_i sigma_b^2))) * rsum_i^2
    n_i <- n_per_cluster
    S1 <- sum(rss)
    p_fixed <- length(fit$coefficients)
    # Two-stage variance-components estimator:
    #
    # (a) Exact ANOVA / method-of-moments for the between/within
    #     decomposition. Uses the exact per-cluster rss and rsum from the
    #     Beaver pipeline; immune to the profile-likelihood collapse-to-0
    #     that the 1-D optimizer suffers when the between-cluster signal
    #     is near the MLE ridge at sigma_b^2 = 0.
    #
    #     SS_total   = sum(rss)                               (Sum_ij r_ij^2)
    #     SS_between = sum(n_i * bar_r_i^2) - n_total*bar_r^2 (cluster means)
    #     SS_within  = SS_total - SS_between
    #     MS_within  = SS_within / (n_total - n_clusters)
    #     MS_between = SS_between / (n_clusters - 1)
    #     sigma_b^2_MoM = max(0, (MS_between - MS_within) / bar_n_eff)
    #
    #   where bar_n_eff = (n_total^2 - sum(n_i^2)) / ((n_clusters-1)*n_total)
    #   is the Satterthwaite harmonic cluster size (unbiased for
    #   unbalanced designs).
    #
    # (b) Profile likelihood as a refinement from the MoM starting
    #     point. If the profile suggests a materially different sigma_b^2
    #     we adopt it; otherwise the MoM value stands.
    bar_r_i      <- rsum / pmax(n_i, 1L)
    bar_r_all    <- sum(rsum) / max(n_total, 1L)
    SS_total     <- S1
    SS_between   <- sum(n_i * (bar_r_i - bar_r_all)^2)
    SS_within    <- max(SS_total - SS_between, 0)
    df_between   <- max(n_clusters - 1L, 1L)
    df_within    <- max(n_total - n_clusters, 1L)
    MS_between   <- SS_between / df_between
    MS_within    <- SS_within  / df_within
    bar_n_eff    <- (n_total^2 - sum(n_i^2)) /
                     (max(n_clusters - 1L, 1L) * max(n_total, 1L))
    if (!is.finite(bar_n_eff) || bar_n_eff <= 0) bar_n_eff <- mean(n_i)
    sigma_b2_MoM <- max((MS_between - MS_within) / bar_n_eff, 0)
    sigma2_new   <- max(MS_within, 1e-10)

    # Task #115 structural fix (2026-04-20): replace the ML profile
    # likelihood maximizer with the analytic REML profile that matches
    # lme4 exactly on unbalanced designs. The MoM and ML profile both
    # converge to estimators that systematically differ from REML for
    # random-intercept LMM on unbalanced cluster sizes (see
    # docs/determinism/lmm_sigma_b2_reml_root_cause.md for derivation
    # and empirical verification).
    #
    # The -2 REML log-likelihood for one-way random-intercept on the
    # per-cluster residual sufficient stats (n_c, rsum_c, rss_c), with
    # sigma^2 fixed at the closed-form (n-p)-df value `sigma2_new`, is:
    #
    #   -2 L_REML(sigma_b^2) = Sum_c (n_c - 1) log sigma^2
    #                   + Sum_c log(sigma^2 + n_c sigma_b^2)
    #                   + Sum_c (rss_c - 2mu*rsum_c + mu^2*n_c) / sigma^2
    #                   - Sum_c sigma_b^2 * (rsum_c - n_c mu)^2 / (sigma^2(sigma^2+n_c sigma_b^2))
    #                   + log Sum_c n_c * w_c            <- REML Jacobian
    #
    # with profiled mu = Sum_c n_c r_c w_c / Sum_c n_c w_c, w_c = 1/(sigma^2+n_c sigma_b^2).
    # This formula matches lme4's REML sigma_b^2 to ~4e-4 rel on our test
    # scenarios when sigma^2 is supplied at the precise (n-p)-df value.
    neg2_L_REML <- function(sb2, s2) {
      if (!is.finite(sb2) || sb2 < 0) sb2 <- 1e-12
      s2v <- max(s2, 1e-12)
      alpha <- s2v + n_i * sb2
      w_c <- 1 / alpha
      denom <- sum(n_i * w_c)
      mu_hat <- if (denom > 0)
        sum(n_i * bar_r_i * w_c) / denom else bar_r_all
      term_logdetV <- sum((n_i - 1) * log(s2v)) + sum(log(alpha))
      term_rVr <- sum((rss - 2 * mu_hat * rsum + mu_hat^2 * n_i) / s2v) -
                  sum((sb2 / (s2v * alpha)) * (rsum - n_i * mu_hat)^2)
      term_jac <- log(denom)
      term_logdetV + term_rVr + term_jac
    }
    # JOINT REML over (sigma^2, sigma_b^2): empirically matches lme4 REML to 4e-5
    # on unbalanced designs (vs ~2e-3 when sigma^2 is fixed first). The 2-D
    # optim uses log-parametrization to enforce positivity and a
    # Nelder-Mead simplex for robustness on the sometimes-flat ridge.
    neg2_L_REML_joint <- function(par) {
      neg2_L_REML(exp(par[2]), exp(par[1]))
    }
    par0 <- c(log(max(sigma2_new, 1e-10)),
              log(max(sigma_b2_MoM, 1e-10)))
    opt_joint <- tryCatch(
      stats::optim(par0, neg2_L_REML_joint,
                    method = "Nelder-Mead",
                    control = list(reltol = 1e-14, maxit = 10000)),
      error = function(e) NULL)
    if (is.null(opt_joint) || !identical(opt_joint$convergence, 0L) ||
        length(opt_joint$par) != 2L || any(!is.finite(opt_joint$par)) ||
        !is.finite(opt_joint$value)) {
      .dsvert_stop_numeric(
        "numeric_backend_unavailable",
        paste0("The joint REML variance-component optimizer failed. ",
               "No fixed-sigma Brent or method-of-moments fallback was applied."),
        reason = "lmm_reml_optimizer_failure")
    }
    sigma2_reml  <- exp(opt_joint$par[1])
    sigma_b2_reml <- exp(opt_joint$par[2])
    if (!is.finite(sigma2_reml) || sigma2_reml <= 0 ||
        !is.finite(sigma_b2_reml) || sigma_b2_reml < 0) {
      .dsvert_stop_numeric(
        "numeric_backend_unavailable",
        "The joint REML optimizer returned invalid variance components.",
        reason = "lmm_reml_optimizer_failure")
    }
    # Use joint-REML sigma_b^2 (exact-to-lmer). sigma^2 is kept at the closed-form
    # (n-p)-df refit since the residual-REML sigma^2 differs from the full-
    # model REML sigma^2 by an (n-1)/(n-p) factor.
    sigma_b2_new <- max(sigma_b2_reml, 0)
    # Oracle / benchmark override: force sigma_b^2 to a caller-supplied
    # value to isolate estimator error from downstream-fit error.
    if (!is.null(sigma_b2_override) && is.finite(sigma_b2_override))
      sigma_b2_new <- as.numeric(sigma_b2_override)
    rho_new <- sigma_b2_new / (sigma2_new + sigma_b2_new)

    # Per-patient weights for next fit.
    # Random intercept only (q=1):
    #   w_ij = 1 - sigma_b^2 / (sigma^2 + n_i sigma_b^2) = sigma^2 / (sigma^2 + n_i sigma_b^2)
    # Random intercept + slopes (q>1): Woodbury gives per-cluster
    #   V_i^{-1} = (1/sigma^2) [I_ni - Z_i (sigma^2 Omega^{-1} + Z_i^T Z_i)^{-1} Z_i^T]
    # We approximate the per-patient weight by the MEAN diagonal of
    # V_i^{-1} scaled by sigma^2 (i.e. the average effective weight in
    # the cluster). This gives a scalar weight per cluster that the
    # expand-column helper can broadcast. The approximation is tight
    # when within-cluster Z rows are similar (typical REML case); exact
    # per-patient weights would require passing V_i^{-1} diagonals back
    # to the server, which is a simple follow-on extension.
    if (q == 1L) {
      rho_i <- sigma_b2_new / (sigma2_new + n_i * sigma_b2_new)
      per_patient_weights_by_cluster <- 1 - rho_i
    } else if (!is.null(Z_info)) {
      per_patient_weights_by_cluster <- numeric(n_clusters)
      Om <- Omega
      # Update Omega diagonally: scale by residual-variance feedback to
      # keep the outer iterate stable in the first pass.
      Om[1L, 1L] <- sigma_b2_new
      Om_inv <- .dsvert_solve_identifiable(
        Om,
        context = "The LMM random-effect covariance",
        reason = "singular_random_effect_covariance",
        symmetric = TRUE)
      for (ci in seq_len(n_clusters)) {
        ZtZ_i <- Z_info$ZtZ[ci, , ]
        M <- sigma2_new * Om_inv + ZtZ_i
        M_inv <- .dsvert_solve_identifiable(
          M,
          context = "An LMM cluster Woodbury system",
          reason = "singular_lmm_cluster_information",
          symmetric = TRUE)
        # Diagonal of V_i^-1 averaged: trace(I_ni / sigma^2 -
        #   Z_i M^{-1} Z_i^T / sigma^2) / n_i
        # = (n_i - trace(ZtZ_i * M_inv)) / (n_i * sigma^2)
        tr <- sum(diag(ZtZ_i %*% M_inv))
        w_bar <- (n_i[ci] - tr) / (n_i[ci])  # pre-multiplied by sigma^2
        per_patient_weights_by_cluster[ci] <- max(w_bar, 1e-6)
      }
      Omega <- Om
    } else {
      rho_i <- sigma_b2_new / (sigma2_new + n_i * sigma_b2_new)
      per_patient_weights_by_cluster <- 1 - rho_i
    }
    # A small estimated ICC is not permission to substitute OLS for the
    # requested LMM. Continue with the same GLS estimand; a singular protected
    # system below fails closed with a typed condition.
    icc_est <- sigma_b2_new / max(sigma2_new + sigma_b2_new, 1e-12)
    if (!is.finite(icc_est)) {
      .dsvert_stop_non_identifiable(
        "The LMM variance ratio is non-finite.",
        reason = "invalid_lmm_variance_ratio")
    }
    # GLS refit: attempt exact closed-form Beaver solve first (matches
    # lme4 to ~2e-3 when the design is well-conditioned). If that
    # fails or yields a near-singular Gram (common when the estimated
    # sigma_b^2 is small -> lambda_i near 0 -> (1-lambda_i) column
    # near-constant -> X'X ill-conditioned beyond Ring63 precision),
    # fall back to the iterative ds.vertGLM path with client-side
    # intercept recovery (historical behaviour, ~|Delta|=0.15 on Pima).
    #
    # For a random-intercept LMM, the exact REML fixed-effects estimate
    # is obtained by transforming each design column and the response:
    #
    #    tilde_v_j = v_j - lambda_i * bar_v_i    (j in cluster i)
    #
    # with lambda_i = 1 - sqrt(sigma^2 / (sigma^2 + n_i sigma_b^2)),
    # AND including an explicit intercept column whose value is
    # (1 - lambda_i) for observations in cluster i (the transform of a
    # constant-1 column). OLS on the transformed system yields beta_GLS
    # equal to lme4's fixed effects up to machine precision. Previous
    # revisions used per-observation weights (WLS with w_j = 1 - rho_i),
    # which is NOT the correct GLS for random-intercept and produced a
    # systematic 5-20% shrinkage in the slopes.
    lambda_i <- 1 - sqrt(sigma2_new / (sigma2_new + n_i * sigma_b2_new))
    all_predictors <- attr(terms(formula), "term.labels")
    x_ysrv   <- x_local_ysrv
    x_peer_v <- setdiff(all_predictors, x_ysrv)
    y_srv_ci <- which(server_names == y_srv)
    peer_srv2 <- peer_srv
    peer_ci2  <- if (!is.null(peer_srv2))
      which(server_names == peer_srv2) else integer(0)
    # Dedicated MPC session for the closed-form round.
    sess_gls <- .mpc_session_id()
    pks_gls <- .dsvert_setup_peer_transport(datasources, server_names,
                                            c(y_srv, peer_srv2), sess_gls)
    # Broadcast cluster IDs to peer.
    cb_gls <- .dsvert_lmm_aggregate_strict(datasources[y_srv_ci],
      call(name = "dsvertLMMBroadcastClusterIDsDS",
           data_name = data, cluster_col = cluster_col,
           peer_pk = pks_gls[[peer_srv2]],
           session_id = sess_gls))
    if (is.list(cb_gls) && length(cb_gls) == 1L) cb_gls <- cb_gls[[1L]]
    .dsvert_store_blob(
      cb_gls$peer_blob, "k2_lmm_cluster_ids_blob",
      datasources[peer_ci2], sess_gls)
    .dsvert_lmm_aggregate_strict(datasources[peer_ci2],
      call(name = "dsvertLMMReceiveClusterIDsDS", session_id = sess_gls))
    # Codex-approved Option 1 (2026-04-19 late): share_scale SNR-boost.
    # Under the CORRECT absolute-noise model (Ring63 Beaver noise ~
    # +/-2^-fracBits abs per TruncMul, NOT relative), multiplying every
    # shared column by c BEFORE Beaver mul boosts signal by c^2 while
    # leaving noise floor unchanged -- net c^2 improvement in relative
    # precision on Gram entries. Headroom analysis in
    # scripts/diag_lmm_gram_magnitudes.R:
    #   balanced  max per-elem |x*y| = 2,476 -> safe c_max = 29
    #   unbalanced max per-elem |x*y| = 2,863 -> safe c_max = 27
    #   combined safe c_max (x2 safety factor over Ring63 ceiling 2^22) = 27
    # We pick c=10 conservatively: c^2=100x noise reduction, per-elem
    # scaled product max = 286k, 14x below Ring63 pre-truncation
    # ceiling 4.19M. Expected X4 rel: 1.78e-4 / 100 = 1.78e-6, crushes
    # 1e-4 STRICT with 56x margin. Solve is c^2-invariant so beta returns in
    # original scale and the legacy quality gate (max|coef|<1e3) still
    # passes. The L2-standardization branch backfired and is retained
    # as a disabled toggle (standardize=FALSE).
    # share_scale under Ring127: the fracBits=50 Uint128 representation
    # has vastly more headroom than Ring63's fracBits=20 uint64, so the
    # SNR-boost factor c=10 (which gives c^2=100x noise reduction by
    # amplifying Gram-entry magnitudes before Beaver mul vs the absolute
    # per-op noise floor) is still well within Ring127's dynamic range.
    # At Ring63 c=10 was already validated; task #121: keep c=10 at
    # Ring127 to close the unbalanced-cluster X4 STRICT gap.
    lmm_share_scale <- if (identical(ring, "ring127")) 10.0 else 1.0
    cf_error <- NULL
    cf_fit <- tryCatch(
      .ds_vertLMM_closed_form(
        conns = datasources, server_names = server_names,
        y_srv = y_srv, peer_srv = peer_srv2,
        data = data, y_var = y_var,
        x_ysrv = x_ysrv, x_peer = x_peer_v,
        lambda_i = lambda_i, transport_pks = pks_gls,
        session_id = sess_gls, verbose = FALSE,
        share_scale = lmm_share_scale,
        standardize = FALSE,
        ring = ring),
      error = function(e) {
        .dsvert_lmm_transport_fallback(e)
        cf_error <<- e
        message("[LMM] closed-form failed: ", conditionMessage(e))
        NULL
      })
    # Cleanup session on both servers.
    for (srv_c in c(y_srv, peer_srv2)) {
      if (is.null(srv_c)) next
      ci <- which(server_names == srv_c)
      .dsvert_cleanup_best_effort(
        datasources[ci], call(name = "mpcCleanupDS", session_id = sess_gls))
    }
    if (is.null(cf_fit)) {
      if (inherits(cf_error, "non_identifiable")) stop(cf_error)
      .dsvert_stop_numeric(
        "numeric_backend_unavailable",
        paste0("The protected LMM GLS backend failed",
               if (!is.null(cf_error)) {
                 paste0(": ", conditionMessage(cf_error))
               } else {
                 "."
               },
               " No OLS fallback was applied."),
        reason = "lmm_gls_backend_failure")
    }
    # Quality gate: closed-form is only reliable when the X'X matrix
    # is well-conditioned. When lambda_i is near-constant (small ICC)
    # the (1-lambda_i) column is near-constant too, and the
    # near-singular Gram matrix amplifies Ring63 FP precision noise
    # into garbage coefficients (max|beta| >> 1 is the symptom).
    # In that regime, use the iterative ds.vertGLM-based path which
    # is numerically stable but has a 0.15 intercept gap on Pima.
    lambda_range <- if (length(lambda_i) > 0L)
      diff(range(lambda_i)) else 0
    max_abs <- if (!is.null(cf_fit) &&
                    is.numeric(cf_fit$coefficients))
      max(abs(cf_fit$coefficients), na.rm = TRUE) else Inf
    # A reference sanity anchor: run a cheap UNWEIGHTED OLS fit first
    # (ds.vertGLM on the ORIGINAL formula) and only accept the closed-
    # form coefficients when they differ from OLS by O(lambda x scale).
    # If the closed-form blows up (Ring63 FP drift, ill-conditioning,
    # or Opal-vs-local Beaver behavioural mismatch) we detect it and
    # keep the OLS reference.
    ols_ref <- ds.vertGLM(formula = formula, data = data,
                           family = "gaussian",
                           max_iter = inner_iter, tol = tol,
                           verbose = FALSE, datasources = datasources)
    closed_form_ok <- !is.null(cf_fit) &&
      is.numeric(cf_fit$coefficients) &&
      all(is.finite(cf_fit$coefficients)) &&
      max_abs < 1e3              # coefs shouldn't explode
    # Note: we previously required lambda_range > 0.02, but that rejected
    # the (common) balanced-design case where all clusters have the same
    # size n_i and hence lambda_i is constant. Balanced GLS is perfectly
    # well-conditioned. The OLS cross-check below provides a safer gate:
    # if cf_fit drifts badly from OLS, reject.
    # Secondary sanity: compare closed-form to OLS. For a properly
    # computed GLS with modest ICC, slopes should not deviate from OLS
    # by more than O(1) in relative scale. If they do, closed-form has
    # drifted and we discard it.
    if (isTRUE(closed_form_ok) && !is.null(ols_ref$coefficients)) {
      common_slopes <- intersect(
        names(cf_fit$coefficients),
        names(ols_ref$coefficients))
      common_slopes <- setdiff(common_slopes, "(Intercept)")
      if (length(common_slopes) > 0L) {
        ols_slopes <- as.numeric(ols_ref$coefficients[common_slopes])
        cf_slopes  <- as.numeric(cf_fit$coefficients[common_slopes])
        diff_rel <- abs(cf_slopes - ols_slopes) /
                     pmax(abs(ols_slopes), 1e-6)
        # GLS coefficients are typically within 0.5x-2x of OLS for
        # moderate ICC; reject if any slope is off by > 5x relative.
        # RING-AWARE RELAX (task #121): at Ring127 the closed-form is
        # precision-validated (max |Deltabeta| vs lmer approx 1e-5 on balanced).
        # High-ICC designs legitimately have GLS beta that differ from
        # OLS beta by much more than 5x (random-effect variance-induced
        # bias on the pooled OLS slopes), so the 5x threshold falsely
        # rejects good Ring127 fits on unbalanced + high-ICC
        # scenarios. Loosen to 50x at Ring127 (still catches
        # catastrophic sign flips or magnitude blow-ups).
        rel_threshold <- if (identical(ring, "ring127")) 50 else 5
        if (any(!is.finite(diff_rel)) || any(diff_rel > rel_threshold))
          closed_form_ok <- FALSE
      }
    }
    if (closed_form_ok) {
      if (verbose) message(sprintf(
        "[LMM] closed-form OK (lambda range %.3f; max |coef| %.3g)",
        lambda_range, max_abs))
      coef_out <- c("(Intercept)" = as.numeric(cf_fit$coefficients["(Intercept)"]))
      for (nm in all_predictors) {
        v <- cf_fit$coefficients[nm]
        coef_out[nm] <- if (length(v) == 1L && is.finite(v)) as.numeric(v) else NA_real_
      }
      # Codex-approved fix (task #108, 2026-04-19 late): replace the
      # MoM sigma^2 (which has 2.22e-2 rel error vs lmer, driven by the
      # exact cross-server residual pipeline's per-cluster r^2 chain)
      # with the closed-form sigma^2 = (ytildeTytilde - betaT XTytilde) / (n - p) using
      # the cf_fit aggregates -- yty (exact local scalar on y_srv, no
      # MPC noise) and Xty (Beaver rel 1e-9 at LMM scale), both
      # already in the scalar-reveal P3 tier. Propagation: sigma^2 rel
      # 2.22e-2 -> ~1e-9 -> lambda_i precise -> beta_slope precision limited by
      # Gram noise floor only (~1e-6). sigma_b^2 MoM formula uses the
      # updated sigma^2 in the MS_within slot (mathematically
      # MS_within == sigma^2). See docs/diagnostic/mpc_sigma2_mom_imprecision.md
      # for the upstream bug we bypass here (task #114).
      nm_beta <- names(cf_fit$coefficients)
      xty_vec <- as.numeric(cf_fit$Xty[nm_beta])
      beta_vec <- as.numeric(cf_fit$coefficients)
      yty_val <- as.numeric(cf_fit$yty)
      rss_client <- yty_val - sum(beta_vec * xty_vec)
      p_fixed_local <- length(cf_fit$coefficients)
      if (is.finite(rss_client) && rss_client > 0) {
        sigma2_new <- max(rss_client / max(n_total - p_fixed_local, 1L),
                           1e-10)
        # DO NOT recompute sigma_b^2 with the new sigma^2. The pre-fix sigma_b^2
        # has rel error 2.67e-5 due to a fortuitous bias cancellation
        # between MS_between and MS_within in the MoM formula (both
        # computed from the same noisy exact-pipeline aggregates).
        # Substituting the precise sigma^2 into the MoM formula breaks
        # this cancellation and degrades sigma_b^2 to rel 5.56e-4.
        # sigma_b^2 remains computed by the existing MoM-from-exact-pipeline
        # path. lambda_i uses the precise sigma^2 and the MoM sigma_b^2.
        rho_new <- sigma_b2_new / (sigma2_new + sigma_b2_new)
        if (verbose) message(sprintf(
          "[LMM] sigma^2 client-side refit: sigma^2=%.6g (was %.6g MoM); sigma_b^2 kept at %.6g",
          sigma2_new, max(MS_within, 1e-10), sigma_b2_new))
      }
      fit <- .dsvert_lmm_closed_form_fit(cf_fit, coef_out, sigma2_new)
    } else {
      .dsvert_stop_numeric(
        "numeric_backend_unavailable",
        sprintf(
          paste0("The protected LMM GLS result failed its numerical ",
                 "attestation (lambda range %.3f; max |coef| %.3g). ",
                 "No OLS fallback was applied."),
          lambda_range, max_abs),
        reason = "lmm_gls_attestation_failure")
    }
    if (verbose) {
      message(sprintf("[LMM] iter %d  sigma^2=%.4g  sigma_b^2=%.4g  rho=%.4g",
                       iter, sigma2_new, sigma_b2_new, rho_new))
    }
    sigma2 <- sigma2_new
    sigma_b2 <- sigma_b2_new
    if (abs(rho_new - rho_prev) < tol) {
      converged <- TRUE
      break
    }
    rho_prev <- rho_new
  }

  icc <- sigma_b2 / (sigma_b2 + sigma2)
  out <- list(
    coefficients = fit$coefficients,
    covariance   = fit$covariance,
    std_errors   = fit$std_errors,
    sigma2       = sigma2,
    sigma_b2     = sigma_b2,
    Omega        = if (q > 1L) Omega else NULL,
    random_slopes = random_slopes,
    q_random     = q,
    icc          = icc,
    n_clusters   = n_clusters,
    cluster_sizes = n_per_cluster,
    converged    = converged,
    iterations   = iter,
    reml         = reml,
    estimand      = if (isTRUE(exact_cross_server)) {
      "random_intercept_lmm_reml"
    } else {
      "explicit_outcome_projection_lmm_approximation"
    },
    residual_backend = if (isTRUE(exact_cross_server)) {
      "exact_cross_server"
    } else {
      "explicit_outcome_projection"
    },
    fit          = fit,
    call         = match.call())
  class(out) <- c("ds.vertLMM", "list")
  out
}

#' @export
print.ds.vertLMM <- function(x, ...) {
  cat("dsVert linear mixed model (random intercept)\n")
  cat(sprintf("  Clusters = %d    N = %d\n",
              x$n_clusters, sum(x$cluster_sizes)))
  cat(sprintf("  sigma^2 = %.4g    sigma_b^2 = %.4g    ICC = %.3f\n",
              x$sigma2, x$sigma_b2, x$icc))
  cat(sprintf("  Converged: %s (%d outer iters)\n",
              x$converged, x$iterations))
  cat("\nFixed effects:\n")
  z <- x$coefficients / x$std_errors
  p <- 2 * stats::pnorm(-abs(z))
  print(round(data.frame(Estimate = x$coefficients, SE = x$std_errors,
                         z = z, p = p, check.names = FALSE), 5L))
  invisible(x)
}
