#' @title Bounded random-intercept LMM compatibility frontdoor
#' @description Admits a signed random-intercept artifact: either the legacy
#'   \code{y ~ 1} moment estimand or a fixed-effect finite-grid GLS estimand.
#'   It uses one existing canonical DP Synopsis;
#'   it neither calls the retired mixed-model endpoints nor exposes cluster
#'   statistics.
#' @details \code{analysis_id} selects the custodian-signed artifact. Random
#'   slopes, unrestricted profile optimisation, covariance estimates, standard
#'   errors and classical inference are unavailable. The intercept-only
#'   artifact returns a method-of-moments projection; the fixed-effect artifact
#'   evaluates its signed finite variance-ratio grid by the signed ML or REML
#'   profile.
#' @param formula A formula matching the signed intercept-only or fixed-effect
#'   random-intercept artifact.
#' @param data Signed protected dataset name or federation.
#' @param cluster_col Cluster column required to match the signed artifact.
#' @param analysis_id Custodian-configured signed random-intercept artifact id.
#' @param random_slopes Must be \code{NULL}.
#' @param reml Must match the signed fixed-effect artifact profile. The
#'   intercept-only moment artifact requires \code{FALSE}.
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
  terms <- if (inherits(formula, "formula")) stats::terms(formula) else NULL
  labels <- if (is.null(terms)) character() else attr(terms, "term.labels")
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]]) ||
      !identical(attr(terms, "intercept"), 1L) ||
      any(grepl("[:*^|()]", labels))) {
    stop("ds.vertLMM requires an intercept-only or additive fixed-effect formula",
         call. = FALSE)
  }
  if (!is.character(cluster_col) || length(cluster_col) != 1L ||
      !nzchar(cluster_col) || is.null(analysis_id) ||
      !is.character(analysis_id) || length(analysis_id) != 1L ||
      !nzchar(analysis_id)) {
    stop("ds.vertLMM requires cluster_col and signed analysis_id strings",
         call. = FALSE)
  }
  if (!is.null(random_slopes) ||
      !is.logical(reml) || length(reml) != 1L || is.na(reml) ||
      !identical(exact_cross_server, TRUE) || !is.null(sigma_b2_override)) {
    stop(paste(
      "ds.vertLMM supports only the signed random-intercept",
      "random-intercept route: random_slopes=NULL,",
      "exact_cross_server=TRUE and sigma_b2_override=NULL"), call. = FALSE)
  }
  signed <- ds.vertDPLMM(
    data_name = data, analysis_id = analysis_id, datasources = datasources)
  outcome <- as.character(formula[[2L]])
  artifact <- signed$signed_artifact
  fixed <- identical(
    artifact$version,
    .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION) ||
    identical(artifact$version,
              .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION)
  artifact_reml <- identical(artifact$estimation_profile, "reml")
  if (!identical(isTRUE(reml), artifact_reml)) {
    stop(if (artifact_reml) {
      "reml=TRUE is required by the signed LMM artifact"
    } else {
      "reml=FALSE is required by the signed LMM artifact"
    }, call. = FALSE)
  }
  expected_labels <- if (fixed) artifact$predictor_order else character()
  if (!identical(sort(labels, method = "radix"), expected_labels)) {
    stop(if (fixed) {
      "formula does not match the signed fixed-effect LMM artifact"
    } else {
      "ds.vertLMM currently supports only an outcome ~ 1 formula"
    }, call. = FALSE)
  }
  if (!identical(artifact$outcome$column, outcome) ||
      !identical(artifact$cluster$column, cluster_col)) {
    stop("formula or cluster_col does not match the signed LMM artifact",
         call. = FALSE)
  }
  signed$call <- match.call()
  signed$reml <- artifact_reml
  signed$random_slopes <- NULL
  signed$converged <- identical(signed$status, "ok")
  signed$iterations <- if (fixed) 1L else 0L
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
  profile <- if (isTRUE(x$reml)) "REML" else "ML/moment"
  cat("  Estimation: signed DP moment projection or ", profile,
      " signed finite grid\n", sep = "")
  cat("\nFixed effects:\n")
  print(round(data.frame(Estimate = x$coefficients,
                         check.names = FALSE), 5L))
  invisible(x)
}
