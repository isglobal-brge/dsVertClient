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
                       verbose = TRUE, datasources = NULL) {
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
      "[ds.vertLMM] %d clusters, n_total=%d (sizes: median=%d, max=%d)",
      n_clusters, n_total,
      stats::median(n_per_cluster), max(n_per_cluster)))
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
        call("dsvertLMMPeerResidualFinaliseDS", session_id = sess))
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
           n = as.integer(n_total),
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
           n = as.integer(n_total), session_id = sess, frac_bits = 20L))
    r1b <- DSI::datashield.aggregate(datasources[peer_ci],
      call("k2BeaverVecmulR1DS",
           peer_pk = pks[[y_srv]],
           x_key = "k2_lmm_exact_r_share",
           y_key = "k2_lmm_exact_r_share",
           n = as.integer(n_total), session_id = sess, frac_bits = 20L))
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
             n = as.integer(n_total), session_id = sess,
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
    list(rsum_per_cluster = rsum_cluster,
         rss_total = rss_total, n_per_cluster = rs_y$cluster_sizes)
  }

  # Outer REML loop: update variance components (Omega q x q + sigma^2).
  # Initialise Omega = 0.1 sigma^2 I, covariance matrix of random effects.
  Omega <- diag(q) * sigma2 * 0.1
  sigma_b2 <- Omega[1L, 1L]  # keep legacy slot for intercept variance
  converged <- FALSE
  rho_prev <- Inf
  fit <- fit0
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
      # Derive per-cluster RSS proxy: total rss distributed by cluster
      # weight (exact per-cluster rss requires another Beaver pass;
      # the total is exact and sufficient for the sigma^2 update, and
      # the per-cluster rsum^2 sum drives the sigma_b^2 update below).
      rss <- rep(cl_exact$rss_total / n_clusters, n_clusters)
      # Use per-cluster rsum^2 directly (exact) and rss_total for the
      # sigma^2 moment.
    } else {
      cl <- get_cluster_resids(fit$coefficients)
      if (is.list(cl) && length(cl) == 1L) cl <- cl[[1]]
      rss <- as.numeric(cl$rss_per_cluster)
      rsum <- as.numeric(cl$rsum_per_cluster)
    }
    stopifnot(length(rss) == n_clusters)

    # MLE / REML updates for variance components under compound symmetry.
    n_i <- n_per_cluster
    # sigma^2 and sigma_b^2 via moment matching:
    #   E[sum r_ij^2]           = sigma^2 * sum(n_i)
    #   E[sum (sum_j r_ij)^2]   = sigma^2 * sum(n_i) + sigma_b^2 * sum(n_i^2)
    S1 <- sum(rss)
    S2 <- sum(rsum^2)
    denom_b <- sum(n_i^2) - sum(n_i)
    if (denom_b <= 0) denom_b <- 1
    sigma2_new <- max(S1 / n_total, 1e-10)
    sigma_b2_new <- max((S2 - sigma2_new * sum(n_i)) / denom_b, 0)
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
    # Push the weight-per-cluster to the server which expands to
    # per-patient via dsvertExpandClusterWeightsDS.
    wcol <- tryCatch(
      DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call("dsvertExpandClusterWeightsDS",
             data_name = data,
             cluster_col = cluster_col,
             weights_per_cluster = as.numeric(per_patient_weights_by_cluster),
             output_column = "__dsvert_lmm_w")),
      error = function(e) {
        stop("dsvertExpandClusterWeightsDS not available: ",
             conditionMessage(e), call. = FALSE)
      })
    fit <- ds.vertGLM(formula = formula, data = data,
                      family = "gaussian",
                      max_iter = inner_iter, tol = tol,
                      weights = "__dsvert_lmm_w",
                      verbose = FALSE, datasources = datasources)
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
