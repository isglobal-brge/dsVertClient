#' @title Signed LMM K>=3 compatibility frontdoor
#' @description This historical K>=3 name delegates to the signed
#'   LMM Synopsis route, including a custodian-signed finite-grid
#'   ML/REML profile or Gaussian random-slope candidate grid.
#' @param formula An intercept-only or additive formula matching the signed
#'   LMM artifact.
#' @param data Signed protected dataset name or federation.
#' @param cluster_col Cluster column required to match the signed artifact.
#' @param analysis_id Custodian-configured signed LMM artifact id.
#' @param random_slopes Optional exact signed random-slope set.
#' @param rho_lo,rho_hi,tol,max_outer,ring,verbose Retained compatibility
#'   controls; they do not alter the signed Synopsis estimand.
#' @param reml Must match the signed fixed-effect artifact profile.
#' @param datasources At least three DataSHIELD connections.
#' @return A \code{ds.vertLMM} / \code{ds.vertDPLMM} object.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertLMM.k3 <- function(formula, data, cluster_col,
                           analysis_id = NULL,
                           random_slopes = NULL,
                           rho_lo = 0.001, rho_hi = 0.999,
                           tol = 1e-4, max_outer = 30L, reml = FALSE,
                           ring = c("ring127", "ring63"),
                           verbose = TRUE, datasources = NULL) {
  datasources <- .dsvert_datasources(datasources)
  if (length(datasources) < 3L) {
    stop("ds.vertLMM.k3 requires at least three DataSHIELD connections",
         call. = FALSE)
  }
  fit <- ds.vertLMM(
    formula = formula, data = data, cluster_col = cluster_col,
    analysis_id = analysis_id, random_slopes = random_slopes,
    reml = reml, max_iter = max_outer,
    tol = tol, ring = ring, verbose = verbose, datasources = datasources)
  .dsvert_set_frontdoor(fit, "ds.vertLMM.k3", "ds.vertLMM", length(datasources))
}

#' @export
print.ds.vertLMM.k3 <- function(x, ...) {
  if (!is.null(x$random_effect_covariance)) {
    class(x) <- unique(c("ds.vertLMM", setdiff(class(x), "ds.vertLMM.k3")))
    return(print.ds.vertLMM(x, ...))
  }
  cat("dsVert signed random-intercept LMM compatibility result (K>=3)\n")
  if (!identical(x$status, "ok")) {
    cat("  Status: ", x$status %||% "unavailable", "\n", sep = "")
    return(invisible(x))
  }
  cat(sprintf("  Clusters = %d    N = %d\n",
              x$cluster_count, x$n_obs))
  cat(sprintf("  sigma^2 = %.4g    sigma_b^2 = %.4g    ICC = %.3f\n",
              x$sigma2, x$sigma_b2, x$icc))
  cat("  Estimation: signed DP moment projection or fixed signed grid\n")
  cat("\nFixed effects:\n")
  print(round(x$coefficients, 5L))
  invisible(x)
}
