#' @title Bounded random-intercept LMM compatibility frontdoor
#' @description Admits only the signed \code{y ~ 1} random-intercept
#'   method-of-moments estimand. It uses one existing canonical DP Synopsis;
#'   it neither calls the retired mixed-model endpoints nor exposes cluster
#'   statistics.
#' @details \code{analysis_id} selects the custodian-signed artifact. Random
#'   slopes, ML/REML, covariance estimates, standard errors and classical
#'   inference are unavailable. The returned variances are DP-projected
#'   method-of-moments components, not an ML/REML fit.
#' @param formula A formula exactly of the form \code{outcome ~ 1}.
#' @param data Signed protected dataset name or federation.
#' @param cluster_col Cluster column required to match the signed artifact.
#' @param analysis_id Custodian-configured signed random-intercept artifact id.
#' @param random_slopes Must be \code{NULL}.
#' @param reml Must be \code{FALSE}; ML/REML is not implemented.
#' @param max_iter,inner_iter,tol,ring,verbose Retained compatibility controls;
#'   they do not alter the signed Synopsis estimand.
#' @param exact_cross_server Must be \code{TRUE}; outcome projection is not
#'   available.
#' @param sigma_b2_override Must be \code{NULL}.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertLMM} / \code{ds.vertDPLMM} object containing only the
#'   certified public DP moment projection.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertLMM <- function(formula, data = NULL, cluster_col,
                       analysis_id = NULL,
                       random_slopes = NULL,
                       reml = FALSE, max_iter = 30L, inner_iter = 50L,
                       tol = 1e-4,
                       exact_cross_server = TRUE,
                       sigma_b2_override = NULL,
                       ring = c("ring63", "ring127"),
                       verbose = TRUE, datasources = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]]) ||
      length(attr(stats::terms(formula), "term.labels")) != 0L ||
      !identical(attr(stats::terms(formula), "intercept"), 1L)) {
    stop("ds.vertLMM currently supports only an outcome ~ 1 formula",
         call. = FALSE)
  }
  if (!is.character(cluster_col) || length(cluster_col) != 1L ||
      !nzchar(cluster_col) || is.null(analysis_id) ||
      !is.character(analysis_id) || length(analysis_id) != 1L ||
      !nzchar(analysis_id)) {
    stop("ds.vertLMM requires cluster_col and signed analysis_id strings",
         call. = FALSE)
  }
  if (!is.null(random_slopes) || !identical(reml, FALSE) ||
      !identical(exact_cross_server, TRUE) || !is.null(sigma_b2_override)) {
    stop(paste(
      "ds.vertLMM supports only the signed random-intercept",
      "method-of-moments route: random_slopes=NULL, reml=FALSE,",
      "exact_cross_server=TRUE and sigma_b2_override=NULL"), call. = FALSE)
  }
  signed <- ds.vertDPLMM(
    data_name = data, analysis_id = analysis_id, datasources = datasources)
  outcome <- as.character(formula[[2L]])
  artifact <- signed$signed_artifact
  if (!identical(artifact$outcome$column, outcome) ||
      !identical(artifact$cluster$column, cluster_col)) {
    stop("formula or cluster_col does not match the signed LMM artifact",
         call. = FALSE)
  }
  signed$call <- match.call()
  signed$reml <- FALSE
  signed$random_slopes <- NULL
  signed$converged <- identical(signed$status, "ok")
  signed$iterations <- 0L
  signed$n_clusters <- signed$cluster_count
  signed$cluster_sizes <- NULL
  signed$estimation_scope <- artifact$estimation_scope
  signed$legacy_fallback_called <- FALSE
  class(signed) <- unique(c("ds.vertLMM", class(signed)))
  return(signed)
}

#' @export
print.ds.vertLMM <- function(x, ...) {
  cat("dsVert bounded random-intercept DP moment model\n")
  cat(sprintf("  Clusters = %d    N = %d\n",
              x$cluster_count, x$n_obs))
  cat(sprintf("  sigma^2 = %.4g    sigma_b^2 = %.4g    ICC = %.3f\n",
              x$sigma2, x$sigma_b2, x$icc))
  cat("  Estimation: signed DP method of moments; no ML/REML inference\n")
  cat("\nFixed effects:\n")
  print(round(data.frame(Estimate = x$coefficients,
                         check.names = FALSE), 5L))
  invisible(x)
}
