#' @title Post-hoc soft-thresholded GLM coefficients (naive LASSO)
#' @description Apply the LASSO soft-threshold operator
#'   \eqn{\mathrm{sign}(\beta_j) \cdot \max(|\beta_j| - \alpha \lambda_1, 0)}
#'   to the coefficient vector of a fitted \code{ds.glm} object, without
#'   re-running the MPC iteration loop. This is the simplest possible
#'   approximation to LASSO: it zeros coefficients whose magnitude falls
#'   below the threshold but does not iteratively re-optimise under the
#'   sparsity constraint.
#'
#'   Useful for quick variable-selection sketches and for checking
#'   whether a subsequent proper proximal-gradient implementation would
#'   pay off. Proper LASSO / elastic net via client-side proximal
#'   gradient with per-iteration MPC gradient calls is a planned
#'   Month 3 deliverable (see V2_PROGRESS.md).
#'
#' @param fit          A \code{ds.glm} object from \code{ds.vertGLM}.
#' @param lambda_1     L1 penalty magnitude.
#' @param alpha_grid   Step-size multipliers to sweep (default 1, 0.5,
#'   0.25, 0.125, 0.0625).
#' @param keep_intercept  If TRUE (default) never threshold the intercept.
#' @return List of class \code{ds.vertLASSO} containing the thresholded
#'   coefficient vectors indexed by effective lambda = alpha * lambda_1.
#' @export
ds.vertLASSO <- function(fit, lambda_1,
                         alpha_grid = c(1, 0.5, 0.25, 0.125, 0.0625),
                         keep_intercept = TRUE) {
  if (!inherits(fit, "ds.glm")) {
    stop("`fit` must be a ds.glm object", call. = FALSE)
  }
  if (!is.numeric(lambda_1) || length(lambda_1) != 1L || lambda_1 < 0) {
    stop("lambda_1 must be a non-negative scalar", call. = FALSE)
  }

  soft_threshold <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

  coef <- fit$coefficients
  int_idx <- which(names(coef) == "(Intercept)")

  paths <- list()
  for (a in alpha_grid) {
    thresh <- a * lambda_1
    shrunk <- soft_threshold(coef, thresh)
    if (keep_intercept && length(int_idx) == 1L) {
      shrunk[int_idx] <- coef[int_idx]
    }
    paths[[sprintf("%.4f", thresh)]] <- shrunk
  }

  out <- list(
    lambda_grid = alpha_grid * lambda_1,
    paths = paths,
    original = coef,
    fit = fit,
    notes = paste0(
      "Naive post-hoc soft-thresholding: zero coefficients below ",
      "|beta| < alpha*lambda_1 without iterative re-optimisation. ",
      "Use proper proximal-gradient LASSO (Month 3) for true convex ",
      "sparsity."))
  class(out) <- c("ds.vertLASSO", "list")
  out
}

#' @export
print.ds.vertLASSO <- function(x, ...) {
  cat("dsVert naive post-hoc LASSO sketch\n")
  cat("  (note: ", x$notes, ")\n\n", sep = "")
  cat(sprintf("Threshold grid: %s\n",
              paste(sprintf("%.4f", x$lambda_grid), collapse = " ")))
  cat("\nCoefficient paths:\n")
  m <- do.call(cbind, x$paths)
  print(round(m, 5L))
  invisible(x)
}
