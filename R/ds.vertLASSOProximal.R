#' @title Gaussian LASSO from one signed DP Synopsis
#' @description Fits one deterministic, KKT-checked Gaussian LASSO solution
#'   from the projected moments in a validated sticky
#'   \code{ds.vertDPGaussian} release. It performs no new DataSHIELD, MPC or
#'   privacy operation. Unauthenticated historical \code{ds.glm} normal-
#'   equation inputs are rejected.
#'
#' @param fit A validated same-owner \code{ds.vertDPGaussian} object.
#' @param lambda One finite non-negative L1 penalty on the signed-bound
#'   normalized Gaussian objective.
#' @param max_iter Positive bound on coordinate-descent iterations.
#' @param tol Positive KKT and convergence tolerance.
#' @param keep_intercept If TRUE (default), the intercept is not penalised.
#' @param warm_start Optional original-scale coefficient vector with the
#'   signed intercept and predictor names.
#' @param accelerate Retained compatibility argument; the deterministic
#'   coordinate-descent solver does not use acceleration.
#' @param trusted_pinset Optional trusted peer pinset for an offline release.
#' @return A \code{ds.vertLASSOProximal} / \code{ds.vertDPLASSO} object with
#'   original- and normalized-scale coefficients, a KKT certificate and no
#'   sampling inference or coefficient confidence regions.
#'
#' @section Disclosure budget:
#' The result is deterministic post-processing of one authenticated sticky
#' Synopsis: additional privacy cost is \code{(epsilon, delta) = (0, 0)} and
#' no additional server calls occur.
#'
#' @seealso \code{\link{ds.vertLASSO}}, \code{\link{ds.vertLASSOCV}}
#' @export
ds.vertLASSOProximal <- function(fit, lambda,
                                 max_iter = 2000L, tol = 1e-9,
                                 keep_intercept = TRUE,
                                 warm_start = NULL,
                                 accelerate = TRUE,
                                 trusted_pinset = NULL) {
  if (!inherits(fit, "ds.vertDPGaussian")) {
    stop(
      "fit must be a validated ds.vertDPGaussian release; the legacy ds.glm normal-equation route is unavailable",
      call. = FALSE)
  }
  .dsvert_lasso_dp_proximal(
    fit = fit, lambda = lambda, max_iter = max_iter, tol = tol,
    keep_intercept = keep_intercept, warm_start = warm_start,
    accelerate = accelerate, trusted_pinset = trusted_pinset)
}

#' @export
print.ds.vertLASSOProximal <- function(x, ...) {
  cat("dsVert Gaussian LASSO from a signed DP Synopsis\n")
  cat(sprintf("  lambda = %.4g  converged = %s (iter %d)\n",
              x$lambda, x$converged, x$iterations))
  cat(sprintf("  objective = %.6g  |support| = %d / %d\n",
              x$objective, length(x$support),
              length(x$coefficients_normalized)))
  cat("  extra DP cost = (0, 0); no sampling inference\n")
  cat("\nOriginal-scale coefficients:\n")
  print(round(x$coefficients, 5L))
  cat("\nNormalized coefficients and KKT certificate:\n")
  print(round(x$coefficients_normalized, 5L))
  cat(sprintf("  max KKT violation = %.6g (tolerance %.6g)\n",
              x$kkt$max_violation, x$kkt$tolerance))
  invisible(x)
}
