#' @title DP-projected information-criterion selection for Gaussian LASSO
#' @description Selects an L1 penalty from the signed projected moments in one
#'   validated same-owner \code{ds.vertDPGaussian} Synopsis. Candidates are
#'   ranked by an explicitly labelled DP-projected pseudo-AIC/BIC/EBIC; the
#'   historical \code{CV} suffix does not mean cross-validation or resampling.
#'   Unauthenticated historical \code{ds.glm} covariance inputs are rejected.
#'
#' @param fit A validated same-owner \code{ds.vertDPGaussian} object.
#' @param lambda_grid Numeric vector of candidate non-negative lambda values,
#'   or \code{NULL} for the signed-source default grid.
#' @param criterion One of \code{"BIC"} (default), \code{"AIC"}, or
#'   \code{"EBIC"}.
#' @param ebic_gamma Extended-BIC gamma parameter, used only for EBIC.
#' @param keep_intercept If TRUE (default), the intercept is not penalised.
#' @param se_threshold Relative pseudo-IC tolerance for the parsimonious
#'   selection. The historical name does not denote a sampling standard error.
#' @param trusted_pinset Optional trusted peer pinset for an offline release.
#' @return A \code{ds.vertLASSOCV} / \code{ds.vertDPLASSOSelect} object with
#'   the selected path, KKT certificates and metadata explicitly identifying
#'   information-criterion selection without cross-validation.
#'
#' @section Disclosure budget:
#' Candidate paths are deterministic post-processing of one authenticated
#' sticky Synopsis: additional privacy cost is \code{(epsilon, delta) =
#' (0, 0)} and no additional server calls occur.
#' @export
ds.vertLASSOCV <- function(fit, lambda_grid = NULL,
                           criterion = c("BIC", "AIC", "EBIC"),
                           ebic_gamma = 0.5,
                           keep_intercept = TRUE,
                           se_threshold = 0.02,
                           trusted_pinset = NULL) {
  if (!inherits(fit, "ds.vertDPGaussian")) {
    stop(
      "fit must be a validated ds.vertDPGaussian release; the legacy ds.glm quadratic-surrogate selector is unavailable",
      call. = FALSE)
  }
  .dsvert_lasso_dp_select(
    fit = fit, lambda_grid = lambda_grid, criterion = match.arg(criterion),
    ebic_gamma = ebic_gamma, keep_intercept = keep_intercept,
    se_threshold = se_threshold, trusted_pinset = trusted_pinset)
}

#' @export
print.ds.vertLASSOCV <- function(x, ...) {
  cat("dsVert DP-projected pseudo-IC LASSO selection (no cross-validation)\n")
  cat(sprintf("  Criterion : %s\n", x$criterion))
  if (identical(x$selection_available, FALSE)) {
    cat("  Selection unavailable: ", x$selection_unavailable_reason,
        "\n", sep = "")
    cat("  No classical IC, CV, or sampling-inference claim is made.\n")
    return(invisible(x))
  }
  cat(sprintf("  lambda.min: %.4g  (IC = %.4g, df = %d)\n",
              x$lambda.min, min(x$ic), x$df[which.min(x$ic)]))
  tolerance <- x$relative_ic_tolerance %||% 0.02
  cat(sprintf("  lambda.parsimonious: %.4g  (sparsest within %g%% of IC_min)\n",
              x$lambda.parsimonious %||% x$lambda.1se,
              100 * tolerance))
  cat("  lambda.1se is a historical alias, not a sampling standard-error rule.\n")
  cat("  Input: validated sticky Gaussian DP Synopsis; extra DP cost = (0, 0)\n")
  cat("\nbeta at lambda.min:\n")
  print(round(x$beta.min, 5L))
  if (!isTRUE(all.equal(x$lambda.min, x$lambda.1se))) {
    cat("\nbeta at lambda.1se:\n")
    print(round(x$beta.1se, 5L))
  }
  invisible(x)
}
