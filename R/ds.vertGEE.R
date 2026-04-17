#' @title Federated generalised estimating equations (skeleton)
#' @description Placeholder for a GEE fit with exchangeable or AR(1)
#'   working correlation structure on vertically partitioned data.
#'
#'   Key property: uses the full Cov(beta_hat) matrix already exposed
#'   by ds.vertGLM (commit 8bb7902) as the "bread" in the
#'   Liang-Zeger sandwich estimator. The "meat" -- sum over clusters
#'   of D_i^T V_i^{-1} r_i r_i^T V_i^{-1} D_i -- is the new MPC
#'   aggregate that this function would compute via per-cluster Beaver
#'   reductions on the share-encoded residuals.
#'
#'   Requires:
#'     - Cluster-ID broadcast (same as LMM; Month 3 plan).
#'     - Per-cluster outer products on shares (Beaver matrix-matrix,
#'       straightforward extension of the existing gradient Beaver).
#'
#'   This file reserves the signature; full implementation is Month 3.
#' @export
ds.vertGEE <- function(formula, data = NULL, id_col,
                       family = "gaussian",
                       corstr = c("independence", "exchangeable", "ar1"),
                       max_iter = 50L, tol = 1e-4,
                       verbose = TRUE, datasources = NULL) {
  corstr <- match.arg(corstr)
  stop("ds.vertGEE: GEE is scheduled for Month 3. The sandwich SE ",
       "reuses fit$covariance (already exposed by ds.vertGLM in ",
       "commit 8bb7902); the per-cluster outer-product 'meat' requires ",
       "the cluster-ID broadcast and a new server helper for per-cluster ",
       "Beaver reductions. See V2_PROGRESS.md.",
       call. = FALSE)
}
