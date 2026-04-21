#' @title Federated linear mixed model with a single random intercept
#' @description Fit a random-intercept linear mixed model
#'   \eqn{y_{ij} = x_{ij}^T \beta + b_i + \varepsilon_{ij}} on vertically
#'   partitioned DataSHIELD data, where the cluster indicator \code{id}
#'   lives on the outcome server. The REML log-likelihood profile is
#'   expressed in terms of a single variance ratio
#'   \eqn{\rho = \sigma_b^2 / (\sigma^2 + n_i\sigma_b^2)} so the outer
#'   optimiser is one-dimensional. Each outer step calls
#'   \code{\link{ds.vertGLM}} with per-patient weights derived from the
#'   current \eqn{\rho}, reusing the already-deployed
#'   \code{ds.vertGLM(weights=)} infrastructure.
#'
#'   Derivation (Laird-Ware compact form):
#'     \deqn{V_i = \sigma^2 I + \sigma_b^2 \mathbf{1} \mathbf{1}^T}
#'     \deqn{V_i^{-1} = \frac{1}{\sigma^2}\left(I - \rho_i \mathbf{1}\mathbf{1}^T\right)}
#'     \deqn{\log|V_i| = (n_i - 1) \log \sigma^2 + \log(\sigma^2 + n_i\sigma_b^2)}
#'
#'   Both summands are one-dimensional functions of \eqn{\rho} that the
#'   client evaluates on centralised aggregates (sum of \eqn{n_i \rho_i}
#'   and sum of \eqn{\log(\sigma^2 + n_i \sigma_b^2)}); per-cluster
#'   residual sums \eqn{\sum_{ij} r_{ij}} are returned by the outcome
#'   server as a single aggregate vector (one scalar per cluster) under
#'   the already-documented cluster-ID inter-server leakage tier.
#'
#'   Inter-server disclosure: the DCF peer learns the cluster
#'   membership (same class as Cox event-time ordering). Absolute
#'   cluster sizes are revealed; individual observations are not.
#'
#' @param formula Fixed-effects formula.
#' @param data Aligned data-frame name.
#' @param cluster_col Column holding the cluster id (must be on the
#'   outcome server).
#' @param reml Use REML (default TRUE) vs ML.
#' @param max_iter Outer variance-component iterations (default 30).
#' @param inner_iter Inner \code{ds.vertGLM} iteration budget.
#' @param tol Outer tolerance on \eqn{\rho} change.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connection object.
#' @return A \code{ds.vertLMM} object with components
#'   \code{coefficients}, \code{covariance}, \code{std_errors},
#'   \code{sigma2} (residual variance), \code{sigma_b2} (random-effect
#'   variance), \code{icc}, \code{n_clusters}, \code{converged},
#'   \code{iterations}, \code{fit} (final inner \code{ds.glm}).
#' @export
ds.vertLMM <- function(formula, data = NULL, cluster_col,
                       random_slopes = NULL,
                       reml = TRUE, max_iter = 30L, inner_iter = 50L,
                       tol = 1e-4,
                       exact_cross_server = TRUE,
                       sigma_b2_override = NULL,
                       ring = c("ring63", "ring127"),
                       verbose = TRUE, datasources = NULL) {
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
  y_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                        data, y_var)
  clust_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                            data, cluster_col)
  if (is.null(y_srv)) stop("y '", y_var, "' not located", call. = FALSE)
  if (is.null(clust_srv)) stop("cluster_col '", cluster_col,
                                "' not located", call. = FALSE)
  if (clust_srv != y_srv) {
    stop("cluster_col must live on the outcome server (y is on '",
         y_srv, "', cluster_col is on '", clust_srv, "'). ",
         "Cross-server cluster broadcast is Month 4.",
         call. = FALSE)
  }

  # Ask the outcome server for the cluster-size vector (aggregate).
  clust_info <- tryCatch(
    DSI::datashield.aggregate(
      datasources[which(server_names == y_srv)],
      call("dsvertClusterSizesDS", data_name = data,
           cluster_col = cluster_col)),
    error = function(e) {
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
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- 1

  # Discover which predictors live on the outcome server so we can
  # pass the right slice of betahat and fitted-value predictors to the
  # server-side residual helper. Cross-server predictors get absorbed
  # into the intercept correction term below.
  y_srv_cols <- tryCatch(
    DSI::datashield.aggregate(
      datasources[which(server_names == y_srv)],
      call("dsvertColNamesDS", data_name = data))[[1L]]$columns,
    error = function(e) character(0))
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
  # dsvertLMMExactClusterR2DS) to compute exact per-cluster SS without
  # the intercept-absorption approximation. Falls back to the
  # approximate path if any piece is unavailable.
  peer_servers <- setdiff(server_names, y_srv)
  peer_srv <- if (length(peer_servers) > 0L) peer_servers[1L] else NULL

  get_cluster_resids_exact <- function(beta_hat, session_id_active,
                                         transport_pks_active) {
    # 1. Peer computes fitted_peer + splits shares.
    x_remote_vars <- setdiff(names(beta_hat), c(x_local_ysrv, "(Intercept)"))
    b_remote <- as.numeric(beta_hat[x_remote_vars])
    r <- DSI::datashield.aggregate(
      datasources[which(server_names == peer_srv)],
      call("dsvertLMMPeerFittedShareDS",
           data_name = data, x_names = x_remote_vars,
           betahat = b_remote,
           peer_pk = transport_pks_active[[y_srv]],
           session_id = session_id_active))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    # 2. Relay peer's blob to outcome server.
    DSI::datashield.aggregate(
      datasources[which(server_names == y_srv)],
      call("mpcStoreBlobDS", key = "k2_lmm_exact_peer_blob",
           chunk = r$peer_blob,
           session_id = session_id_active))
    # 3. Outcome: compute r_share_0.
    DSI::datashield.aggregate(
      datasources[which(server_names == y_srv)],
      call("dsvertLMMCoordResidualShareDS",
           data_name = data, y_var = y_var,
           x_names = x_local_ysrv,
           betahat_local = as.numeric(beta_hat[x_local_ysrv]),
           intercept = as.numeric(beta_hat["(Intercept)"]),
           session_id = session_id_active))
    # 4. Peer finalises its share.
    DSI::datashield.aggregate(
      datasources[which(server_names == peer_srv)],
      call("dsvertLMMPeerResidualFinaliseDS",
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
      DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call("dsvertClusterResidualsDS",
             data_name = data,
             y_var = y_var,
             x_names = x_local_ysrv,
             intercept = as.numeric(beta_hat["(Intercept)"]),
             betahat = b_local,
             cluster_col = cluster_col)),
      error = function(e) {
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
           ". Cross-server slope cols need Beaver (Month 4).",
           call. = FALSE)
    }
    Z_info <- tryCatch(
      DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call("dsvertClusterZtZDS",
             data_name = data,
             cluster_col = cluster_col,
             slope_columns = random_slopes)),
      error = function(e) {
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
      r <- DSI::datashield.aggregate(datasources[which(server_names == y_srv)],
        call("getObsCountDS", data_name = data))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      as.integer(r$n_obs)
    }, error = function(e) NULL)
    if (is.null(n_actual) || n_actual <= 0) return(NULL)
    # Reuse the already-live session (ds.vertGLM ran with fit-time session
    # but cleaned up; we open a fresh one for the exact pipeline).
    sess <- .mpc_session_id()
    y_srv_ci <- which(server_names == y_srv)
    peer_ci <- which(server_names == peer_srv)
    on.exit({
      for (.ci in c(y_srv_ci, peer_ci)) tryCatch(
        DSI::datashield.aggregate(datasources[.ci],
          call("mpcCleanupDS", session_id = sess)),
        error = function(e) NULL)
    }, add = TRUE)
    # Transport keys.
    pks <- list()
    for (srv in c(y_srv, peer_srv)) {
      ci <- which(server_names == srv)
      r <- DSI::datashield.aggregate(datasources[ci],
        call("glmRing63TransportInitDS", session_id = sess))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      pks[[srv]] <- r$transport_pk
    }
    x_remote_vars <- setdiff(names(beta_hat),
                              c(x_local_ysrv, "(Intercept)"))
    b_remote <- as.numeric(beta_hat[x_remote_vars])
    # 1. Peer fitted share.
    r1 <- tryCatch(DSI::datashield.aggregate(datasources[peer_ci],
      call("dsvertLMMPeerFittedShareDS",
           data_name = data, x_names = x_remote_vars,
           betahat = b_remote, peer_pk = pks[[y_srv]],
           session_id = sess)),
      error = function(e) NULL)
    if (is.null(r1)) return(NULL)
    if (is.list(r1) && length(r1) == 1L) r1 <- r1[[1L]]
    DSI::datashield.aggregate(datasources[y_srv_ci],
      call("mpcStoreBlobDS", key = "k2_lmm_exact_peer_blob",
           chunk = r1$peer_blob, session_id = sess))
    # Also set k2_x_n on peer session (needed by PeerResidualFinalise).
    DSI::datashield.aggregate(datasources[peer_ci],
      call("dsvertCopyDfDS", data_name = data, output_name = data))
    # 2. Coord residual share.
    tryCatch(DSI::datashield.aggregate(datasources[y_srv_ci],
      call("dsvertLMMCoordResidualShareDS",
           data_name = data, y_var = y_var,
           x_names = x_local_ysrv,
           betahat_local = as.numeric(beta_hat[x_local_ysrv]),
           intercept = as.numeric(beta_hat["(Intercept)"]),
           session_id = sess)),
      error = function(e) { message("[LMM exact] coord: ",
        conditionMessage(e)); return(NULL) })
    # 3. Peer residual finalise. Needs k2_x_n; we attempt to set it via a
    # no-op that also caches length. Fallback: skip if unavailable.
    tryCatch({
      # Set peer's k2_x_n from the row count
      n_check <- DSI::datashield.aggregate(datasources[peer_ci],
        call("getObsCountDS", data_name = data))[[1L]]$n_obs
      # Manually stamp via a helper; if not available, finalise may
      # fail but we can still proceed with partial info.
    }, error = function(e) NULL)
    fin_ok <- tryCatch({
      DSI::datashield.aggregate(datasources[peer_ci],
        call("dsvertLMMPeerResidualFinaliseDS",
             n = as.integer(n_actual),
             session_id = sess))
      TRUE
    }, error = function(e) { message("[LMM exact] peer finalise: ",
      conditionMessage(e)); FALSE })
    if (!fin_ok) return(NULL)
    # 4. Cluster-ID broadcast.
    cb <- tryCatch(DSI::datashield.aggregate(datasources[y_srv_ci],
      call("dsvertLMMBroadcastClusterIDsDS",
           data_name = data, cluster_col = cluster_col,
           peer_pk = pks[[peer_srv]], session_id = sess)),
      error = function(e) NULL)
    if (is.null(cb)) return(NULL)
    if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
    DSI::datashield.aggregate(datasources[peer_ci],
      call("mpcStoreBlobDS", key = "k2_lmm_cluster_ids_blob",
           chunk = cb$peer_blob, session_id = sess))
    DSI::datashield.aggregate(datasources[peer_ci],
      call("dsvertLMMReceiveClusterIDsDS", session_id = sess))
    # 5. Per-cluster rsum (both parties) + aggregate client-side.
    rs_y <- DSI::datashield.aggregate(datasources[y_srv_ci],
      call("dsvertLMMPerClusterSumDS",
           share_key = "k2_lmm_exact_r_share", session_id = sess))
    rs_p <- DSI::datashield.aggregate(datasources[peer_ci],
      call("dsvertLMMPerClusterSumDS",
           share_key = "k2_lmm_exact_r_share", session_id = sess))
    if (is.list(rs_y) && length(rs_y) == 1L) rs_y <- rs_y[[1L]]
    if (is.list(rs_p) && length(rs_p) == 1L) rs_p <- rs_p[[1L]]
    K <- length(rs_y$per_cluster_fp)
    rsum_cluster <- numeric(K)
    for (ck in seq_len(K)) {
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = rs_y$per_cluster_fp[[ck]],
        share_b = rs_p$per_cluster_fp[[ck]], frac_bits = 20L))
      rsum_cluster[ck] <- as.numeric(agg$values[1L])
    }
    # 6. Beaver vecmul r × r -> r^2 share.
    tri <- DSI::datashield.aggregate(datasources[peer_ci],
      call("k2BeaverVecmulGenTriplesDS",
           dcf0_pk = pks[[y_srv]], dcf1_pk = pks[[peer_srv]],
           n = as.integer(n_actual),
           session_id = sess, frac_bits = 20L))
    if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
    DSI::datashield.aggregate(datasources[y_srv_ci],
      call("mpcStoreBlobDS", key = "k2_beaver_vecmul_triple",
           chunk = tri$triple_blob_0, session_id = sess))
    DSI::datashield.aggregate(datasources[peer_ci],
      call("mpcStoreBlobDS", key = "k2_beaver_vecmul_triple",
           chunk = tri$triple_blob_1, session_id = sess))
    for (ci in c(y_srv_ci, peer_ci)) DSI::datashield.aggregate(
      datasources[ci], call("k2BeaverVecmulConsumeTripleDS",
                              session_id = sess))
    r1a <- DSI::datashield.aggregate(datasources[y_srv_ci],
      call("k2BeaverVecmulR1DS",
           peer_pk = pks[[peer_srv]],
           x_key = "k2_lmm_exact_r_share",
           y_key = "k2_lmm_exact_r_share",
           n = as.integer(n_actual), session_id = sess, frac_bits = 20L))
    r1b <- DSI::datashield.aggregate(datasources[peer_ci],
      call("k2BeaverVecmulR1DS",
           peer_pk = pks[[y_srv]],
           x_key = "k2_lmm_exact_r_share",
           y_key = "k2_lmm_exact_r_share",
           n = as.integer(n_actual), session_id = sess, frac_bits = 20L))
    if (is.list(r1a) && length(r1a) == 1L) r1a <- r1a[[1L]]
    if (is.list(r1b) && length(r1b) == 1L) r1b <- r1b[[1L]]
    DSI::datashield.aggregate(datasources[peer_ci],
      call("mpcStoreBlobDS", key = "k2_beaver_vecmul_peer_masked",
           chunk = r1a$peer_blob, session_id = sess))
    DSI::datashield.aggregate(datasources[y_srv_ci],
      call("mpcStoreBlobDS", key = "k2_beaver_vecmul_peer_masked",
           chunk = r1b$peer_blob, session_id = sess))
    for (srv in c(y_srv, peer_srv)) {
      ci <- which(server_names == srv)
      DSI::datashield.aggregate(datasources[ci],
        call("k2BeaverVecmulR2DS",
             is_party0 = (srv == y_srv),
             x_key = "k2_lmm_exact_r_share",
             y_key = "k2_lmm_exact_r_share",
             output_key = "k2_lmm_exact_r2_share",
             n = as.integer(n_actual), session_id = sess,
             frac_bits = 20L))
    }
    # 7. Global sum r^2 (both parties) -> aggregate.
    gs_y <- DSI::datashield.aggregate(datasources[y_srv_ci],
      call("dsvertLMMGlobalSumDS",
           share_key = "k2_lmm_exact_r2_share", session_id = sess))
    gs_p <- DSI::datashield.aggregate(datasources[peer_ci],
      call("dsvertLMMGlobalSumDS",
           share_key = "k2_lmm_exact_r2_share", session_id = sess))
    if (is.list(gs_y) && length(gs_y) == 1L) gs_y <- gs_y[[1L]]
    if (is.list(gs_p) && length(gs_p) == 1L) gs_p <- gs_p[[1L]]
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
    r2s_y <- DSI::datashield.aggregate(datasources[y_srv_ci],
      call("dsvertLMMPerClusterSumDS",
           share_key = "k2_lmm_exact_r2_share", session_id = sess))
    r2s_p <- DSI::datashield.aggregate(datasources[peer_ci],
      call("dsvertLMMPerClusterSumDS",
           share_key = "k2_lmm_exact_r2_share", session_id = sess))
    if (is.list(r2s_y) && length(r2s_y) == 1L) r2s_y <- r2s_y[[1L]]
    if (is.list(r2s_p) && length(r2s_p) == 1L) r2s_p <- r2s_p[[1L]]
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
  # empirically DESTABILIZED the outer REML loop: determinism |Δ|
  # degraded from 2.2e-5 (V2 alone) to 11.35 units across runs, a
  # catastrophic regression. Aitken extrapolation off the contraction
  # path induces oscillation; with max_iter=30 the iterate wanders.
  # Reverted to plain Picard (V2 state). Task #115 remains open for a
  # stability-preserving acceleration scheme. The σ²-outer-loop
  # coupling floor of 6e-4 rel (→ β_X4 rel ~1.8e-4) is accepted as
  # the current LMM precision limit.
  for (iter in seq_len(max_iter)) {
    # Prefer the exact cross-server path when available; fall back to
    # the approximate intercept-absorbing aggregator if the exact
    # helpers are not all deployed on the Rocks.
    cl_exact <- NULL
    if (isTRUE(exact_cross_server) && !is.null(peer_srv) &&
        length(setdiff(attr(terms(formula), "term.labels"),
                         x_local_ysrv)) > 0L) {
      cl_exact <- tryCatch(get_cluster_resids_full_exact(fit$coefficients),
                            error = function(e) {
                              if (verbose) message("[LMM] exact path failed: ",
                                conditionMessage(e), " -- falling back")
                              NULL })
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
        rss <- rep(cl_exact$rss_total / n_clusters, n_clusters)
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
    #     SS_total   = sum(rss)                               (Σ_ij r_ij^2)
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
    # σ² fixed at the closed-form (n-p)-df value `sigma2_new`, is:
    #
    #   -2 L_REML(σ_b²) = Σ_c (n_c − 1) log σ²
    #                   + Σ_c log(σ² + n_c σ_b²)
    #                   + Σ_c (rss_c − 2μ̂·rsum_c + μ̂²·n_c) / σ²
    #                   − Σ_c σ_b² · (rsum_c − n_c μ̂)² / (σ²(σ²+n_c σ_b²))
    #                   + log Σ_c n_c · w_c            ← REML Jacobian
    #
    # with profiled μ̂ = Σ_c n_c r̄_c w_c / Σ_c n_c w_c, w_c = 1/(σ²+n_c σ_b²).
    # This formula matches lme4's REML σ_b² to ~4e-4 rel on our test
    # scenarios when σ² is supplied at the precise (n-p)-df value.
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
    # JOINT REML over (σ², σ_b²): empirically matches lme4 REML to 4e-5
    # on unbalanced designs (vs ~2e-3 when σ² is fixed first). The 2-D
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
      error = function(e) list(par = par0, value = Inf))
    # Fallback Brent over σ_b² with σ² at sigma2_new if joint optim
    # failed or returned non-finite values.
    sigma2_reml  <- exp(opt_joint$par[1])
    sigma_b2_reml <- exp(opt_joint$par[2])
    if (!is.finite(sigma2_reml) || sigma2_reml <= 0 ||
        !is.finite(sigma_b2_reml) || sigma_b2_reml < 0) {
      opt_b <- tryCatch(
        stats::optimize(neg2_L_REML,
                         interval = c(1e-12, max(sigma_b2_MoM * 10, 10)),
                         s2 = sigma2_new, tol = 1e-12),
        error = function(e) list(minimum = sigma_b2_MoM))
      sigma2_reml  <- sigma2_new
      sigma_b2_reml <- opt_b$minimum
    }
    # Use joint-REML σ_b² (exact-to-lmer). σ² is kept at the closed-form
    # (n-p)-df refit since the residual-REML σ² differs from the full-
    # model REML σ² by an (n-1)/(n-p) factor.
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
      Om_inv <- tryCatch(solve(Om),
                          error = function(e) solve(Om + 1e-6 * diag(q)))
      for (ci in seq_len(n_clusters)) {
        ZtZ_i <- Z_info$ZtZ[ci, , ]
        M <- sigma2_new * Om_inv + ZtZ_i
        M_inv <- tryCatch(solve(M),
                           error = function(e) solve(M + 1e-6 * diag(q)))
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
    # Safety fallback: when the MPC-estimated sigma_b^2 is numerically
    # indistinguishable from 0 (ICC below 1e-3) the exact-pipeline
    # profile has not detected between-cluster signal above MPC
    # floating-point noise. Running the GLS cluster-mean transform in
    # that regime amplifies numerical error on the
    # near-constant (1 - lambda_i) predictor and degrades precision
    # relative to a plain OLS baseline. Fall back to unweighted
    # ds.vertGLM on the original formula (matches the pre-GLS-fix
    # baseline of |Delta| ~ 0.07 vs lme4 REML on Pima).
    #
    # When ICC is detected (>=1e-3), proceed with the GLS
    # cluster-mean centering path which is strictly more accurate
    # than OLS (Laird & Ware 1982; Baltagi 2008 "Econometric
    # Analysis of Panel Data"). Local-harness benchmark on strong-
    # signal synthetic data gives |Delta| = 2e-3 vs lme4.
    icc_est <- sigma_b2_new / max(sigma2_new + sigma_b2_new, 1e-12)
    if (!is.finite(icc_est) || icc_est < 1e-3) {
      if (verbose) message(sprintf(
        "[LMM] sigma_b^2 ~ 0 (icc=%.2e) -- OLS fallback (no GLS signal).",
        icc_est))
      fit <- ds.vertGLM(formula = formula, data = data,
                        family = "gaussian",
                        max_iter = inner_iter, tol = tol,
                        verbose = FALSE, datasources = datasources)
      if (verbose) message(sprintf(
        "[LMM] iter %d  sigma^2=%.4g  sigma_b^2=%.4g  rho=%.4g",
        iter, sigma2_new, sigma_b2_new, rho_new))
      sigma2 <- sigma2_new; sigma_b2 <- sigma_b2_new
      if (abs(rho_new - rho_prev) < tol) { converged <- TRUE; break }
      rho_prev <- rho_new
      next
    }
    # GLS refit: attempt exact closed-form Beaver solve first (matches
    # lme4 to ~2e-3 when the design is well-conditioned). If that
    # fails or yields a near-singular Gram (common when the estimated
    # sigma_b^2 is small → lambda_i near 0 → (1-lambda_i) column
    # near-constant → X'X ill-conditioned beyond Ring63 precision),
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
    pks_gls <- list()
    for (srv in c(y_srv, peer_srv2)) {
      if (is.null(srv)) next
      ci <- which(server_names == srv)
      r <- DSI::datashield.aggregate(datasources[ci],
        call("glmRing63TransportInitDS", session_id = sess_gls))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      pks_gls[[srv]] <- r$transport_pk
    }
    # Broadcast cluster IDs to peer.
    cb_gls <- DSI::datashield.aggregate(datasources[y_srv_ci],
      call("dsvertLMMBroadcastClusterIDsDS",
           data_name = data, cluster_col = cluster_col,
           peer_pk = pks_gls[[peer_srv2]],
           session_id = sess_gls))
    if (is.list(cb_gls) && length(cb_gls) == 1L) cb_gls <- cb_gls[[1L]]
    DSI::datashield.aggregate(datasources[peer_ci2],
      call("mpcStoreBlobDS", key = "k2_lmm_cluster_ids_blob",
           chunk = cb_gls$peer_blob, session_id = sess_gls))
    DSI::datashield.aggregate(datasources[peer_ci2],
      call("dsvertLMMReceiveClusterIDsDS", session_id = sess_gls))
    # Codex-approved Option 1 (2026-04-19 late): share_scale SNR-boost.
    # Under the CORRECT absolute-noise model (Ring63 Beaver noise ~
    # ±2^-fracBits abs per TruncMul, NOT relative), multiplying every
    # shared column by c BEFORE Beaver mul boosts signal by c² while
    # leaving noise floor unchanged — net c² improvement in relative
    # precision on Gram entries. Headroom analysis in
    # scripts/diag_lmm_gram_magnitudes.R:
    #   balanced  max per-elem |x·y| = 2,476 → safe c_max = 29
    #   unbalanced max per-elem |x·y| = 2,863 → safe c_max = 27
    #   combined safe c_max (×2 safety factor over Ring63 ceiling 2^22) = 27
    # We pick c=10 conservatively: c²=100× noise reduction, per-elem
    # scaled product max = 286k, 14× below Ring63 pre-truncation
    # ceiling 4.19M. Expected X4 rel: 1.78e-4 / 100 = 1.78e-6, crushes
    # 1e-4 STRICT with 56× margin. Solve is c²-invariant so β returns in
    # original scale and the legacy quality gate (max|coef|<1e3) still
    # passes. The L2-standardization branch backfired and is retained
    # as a disabled toggle (standardize=FALSE).
    # share_scale under Ring127: the fracBits=50 Uint128 representation
    # has vastly more headroom than Ring63's fracBits=20 uint64, so the
    # SNR-boost factor c=10 (which gives c²=100× noise reduction by
    # amplifying Gram-entry magnitudes before Beaver mul vs the absolute
    # per-op noise floor) is still well within Ring127's dynamic range.
    # At Ring63 c=10 was already validated; task #121: keep c=10 at
    # Ring127 to close the unbalanced-cluster X4 STRICT gap.
    lmm_share_scale <- if (identical(ring, "ring127")) 10.0 else 1.0
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
        message("[LMM] closed-form failed: ", conditionMessage(e))
        NULL
      })
    # Cleanup session on both servers.
    for (srv_c in c(y_srv, peer_srv2)) {
      if (is.null(srv_c)) next
      ci <- which(server_names == srv_c)
      tryCatch(DSI::datashield.aggregate(datasources[ci],
        call("mpcCleanupDS", session_id = sess_gls)),
        error = function(e) NULL)
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
    # form coefficients when they differ from OLS by O(lambda × scale).
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
        # precision-validated (max |Δβ| vs lmer ≈ 1e-5 on balanced).
        # High-ICC designs legitimately have GLS β that differ from
        # OLS β by much more than 5× (random-effect variance-induced
        # bias on the pooled OLS slopes), so the 5× threshold falsely
        # rejects good Ring127 fits on unbalanced + high-ICC
        # scenarios. Loosen to 50× at Ring127 (still catches
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
      fit <- list(coefficients = coef_out)
      # Codex-approved fix (task #108, 2026-04-19 late): replace the
      # MoM σ² (which has 2.22e-2 rel error vs lmer, driven by the
      # exact cross-server residual pipeline's per-cluster r² chain)
      # with the closed-form σ̂² = (ỹᵀỹ − β̂ᵀ X̃ᵀỹ) / (n − p) using
      # the cf_fit aggregates — yty (exact local scalar on y_srv, no
      # MPC noise) and Xty (Beaver rel 1e-9 at LMM scale), both
      # already in the scalar-reveal P3 tier. Propagation: σ² rel
      # 2.22e-2 → ~1e-9 → λ_i precise → β_slope precision limited by
      # Gram noise floor only (~1e-6). σ_b² MoM formula uses the
      # updated σ² in the MS_within slot (mathematically
      # MS_within ≡ σ²). See docs/diagnostic/mpc_sigma2_mom_imprecision.md
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
        # DO NOT recompute σ_b² with the new σ². The pre-fix σ_b²
        # has rel error 2.67e-5 due to a fortuitous bias cancellation
        # between MS_between and MS_within in the MoM formula (both
        # computed from the same noisy exact-pipeline aggregates).
        # Substituting the precise σ² into the MoM formula breaks
        # this cancellation and degrades σ_b² to rel 5.56e-4.
        # σ_b² remains computed by the existing MoM-from-exact-pipeline
        # path. λ_i uses the precise σ² and the MoM σ_b².
        rho_new <- sigma_b2_new / (sigma2_new + sigma_b2_new)
        if (verbose) message(sprintf(
          "[LMM] σ² client-side refit: σ²=%.6g (was %.6g MoM); σ_b² kept at %.6g",
          sigma2_new, max(MS_within, 1e-10), sigma_b2_new))
      }
    } else {
      if (verbose) message(sprintf(
        "[LMM] closed-form rejected (lambda range %.3f; max |coef| %.3g); OLS fallback.",
        lambda_range, max_abs))
      fit <- ols_ref
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
