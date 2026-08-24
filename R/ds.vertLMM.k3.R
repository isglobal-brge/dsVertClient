#' @title Signed random-intercept LMM K>=3 compatibility frontdoor
#' @description This historical K>=3 name delegates only to the signed
#'   random-intercept method-of-moments Synopsis route. It does not activate
#'   the retired K>=3 profile/REML implementation.
#' @param formula A formula exactly of the form \code{outcome ~ 1}.
#' @param data Signed protected dataset name or federation.
#' @param cluster_col Cluster column required to match the signed artifact.
#' @param analysis_id Custodian-configured signed random-intercept artifact id.
#' @param rho_lo,rho_hi,tol,max_outer,ring,verbose Retained compatibility
#'   controls; they do not alter the signed Synopsis estimand.
#' @param datasources At least three DataSHIELD connections.
#' @return A \code{ds.vertLMM} / \code{ds.vertDPLMM} object.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertLMM.k3 <- function(formula, data, cluster_col,
                           analysis_id = NULL,
                           rho_lo = 0.001, rho_hi = 0.999,
                           tol = 1e-4, max_outer = 30L,
                           ring = c("ring127", "ring63"),
                           verbose = TRUE, datasources = NULL) {
  datasources <- .dsvert_datasources(datasources)
  if (length(datasources) < 3L) {
    stop("ds.vertLMM.k3 requires at least three DataSHIELD connections",
         call. = FALSE)
  }
  fit <- ds.vertLMM(
    formula = formula, data = data, cluster_col = cluster_col,
    analysis_id = analysis_id, reml = FALSE, max_iter = max_outer,
    tol = tol, ring = ring, verbose = verbose, datasources = datasources)
  .dsvert_set_frontdoor(fit, "ds.vertLMM.k3", "ds.vertLMM", length(datasources))
}

#' @export
print.ds.vertLMM.k3 <- function(x, ...) {
  cat("dsVert signed random-intercept LMM compatibility result (K>=3)\n")
  if (!identical(x$status, "ok")) {
    cat("  Status: ", x$status %||% "unavailable", "\n", sep = "")
    return(invisible(x))
  }
  cat(sprintf("  Clusters = %d    N = %d\n",
              x$cluster_count, x$n_obs))
  cat(sprintf("  sigma^2 = %.4g    sigma_b^2 = %.4g    ICC = %.3f\n",
              x$sigma2, x$sigma_b2, x$icc))
  cat("  Estimation: signed DP method of moments; no K>=3 REML/profile fit\n")
  cat("\nFixed effects:\n")
  print(round(x$coefficients, 5L))
  invisible(x)
}
