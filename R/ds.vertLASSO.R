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


#' @title One-step LASSO via quadratic-surrogate proximal gradient
#' @description Fit a proper L1-regularised coefficient vector by
#'   expanding the log-likelihood at the converged GLM point as a local
#'   quadratic and then running proximal coordinate descent with soft-
#'   thresholding to minimise
#'
#'     0.5 (beta - betahat)^T H (beta - betahat) + lambda ||beta||_1
#'
#'   where H = solve(Cov(betahat)) is the Fisher information matrix
#'   already exposed by \code{ds.vertGLM} (\code{fit$covariance}). This
#'   yields a consistent, efficient L1-regularised estimator in large
#'   samples without a second MPC round: it uses only the already-
#'   returned coefficient vector and full covariance matrix.
#'
#'   Useful for:
#'     - Regularised post-selection inference.
#'     - Model compression while preserving the statistical geometry.
#'     - A baseline for the planned Month 3 proper proximal-gradient
#'       LASSO that fits from scratch with per-iteration MPC gradients.
#'
#' @param fit      A \code{ds.glm} object with \code{fit$covariance}
#'                 populated (commits >= 8bb7902).
#' @param lambda   Numeric vector of L1 penalty values (a regularisation
#'                 path).
#' @param keep_intercept Never penalise the intercept.
#' @param max_iter Coordinate-descent iterations per lambda.
#' @param tol      Convergence tolerance on max |Delta beta|.
#' @return A ds.vertLASSO1Step object: per-lambda coefficient vectors,
#'         the penalty path, the quadratic-surrogate objective at each
#'         lambda, and the fit.
#' @export
ds.vertLASSO1Step <- function(fit, lambda,
                              keep_intercept = TRUE,
                              max_iter = 500L,
                              tol = 1e-8) {
  if (!inherits(fit, "ds.glm")) {
    stop("`fit` must be a ds.glm object", call. = FALSE)
  }
  if (is.null(fit$covariance)) {
    stop("fit does not expose the full covariance matrix; refit with a
          dsVert version >= 8bb7902.", call. = FALSE)
  }
  if (!is.numeric(lambda) || any(lambda < 0)) {
    stop("lambda must be a non-negative numeric vector", call. = FALSE)
  }

  cov <- as.matrix(fit$covariance)
  # Fisher = inverse of covariance. Symmetrise for numerical stability.
  H <- tryCatch(solve((cov + t(cov)) / 2), error = function(e) NULL)
  if (is.null(H)) {
    stop("Covariance matrix is singular; cannot form H = Cov^{-1}",
         call. = FALSE)
  }
  H <- (H + t(H)) / 2
  diag_H <- diag(H)
  if (any(diag_H <= 0)) {
    stop("Non-positive Hessian diagonal; check fit quality", call. = FALSE)
  }

  betahat <- as.numeric(fit$coefficients)
  names(betahat) <- names(fit$coefficients)
  p <- length(betahat)
  int_idx <- which(names(betahat) == "(Intercept)")

  soft_threshold <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

  # For each lambda in the path, run coordinate descent with warm-start
  # initialised from betahat. Coordinate update for coordinate j:
  #   r_j = betahat_j - (1/H_jj) * sum_{k != j} H_jk (beta_k - betahat_k)
  #   beta_j <- soft_threshold(r_j, lambda / H_jj)
  # Intercept (if requested) is never thresholded: beta_j <- r_j.

  paths <- list()
  objective_at_lambda <- numeric(length(lambda))
  for (li in seq_along(lambda)) {
    lam <- lambda[li]
    beta <- betahat
    for (iter in seq_len(max_iter)) {
      max_delta <- 0
      for (j in seq_len(p)) {
        delta_k <- beta - betahat
        delta_k[j] <- 0  # exclude j from the correction sum
        r_j <- betahat[j] - sum(H[j, ] * delta_k) / diag_H[j]
        if (keep_intercept && length(int_idx) == 1L && j == int_idx) {
          new_j <- r_j
        } else {
          new_j <- soft_threshold(r_j, lam / diag_H[j])
        }
        d <- abs(new_j - beta[j])
        if (d > max_delta) max_delta <- d
        beta[j] <- new_j
      }
      if (max_delta < tol) break
    }
    # Final objective value at this lambda
    delta <- beta - betahat
    quad <- 0.5 * drop(t(delta) %*% H %*% delta)
    l1 <- lam * sum(abs(beta[setdiff(seq_len(p), if (keep_intercept) int_idx else integer(0))]))
    objective_at_lambda[li] <- quad + l1
    paths[[sprintf("%.6g", lam)]] <- beta
  }

  out <- list(
    lambda = lambda,
    paths = paths,
    original = betahat,
    covariance = cov,
    H = H,
    objective = objective_at_lambda,
    fit = fit)
  class(out) <- c("ds.vertLASSO1Step", "list")
  out
}

#' @export
print.ds.vertLASSO1Step <- function(x, ...) {
  cat("dsVert one-step LASSO (quadratic-surrogate proximal)\n")
  cat(sprintf("  Fit family : %s\n", x$fit$family))
  cat(sprintf("  N          : %d\n", x$fit$n_obs))
  cat(sprintf("  p          : %d\n", length(x$original)))
  cat(sprintf("  Lambda grid: %s\n",
              paste(sprintf("%.4g", x$lambda), collapse = " ")))
  m <- do.call(cbind, x$paths)
  cat("\nCoefficient path (rows = coefficients, columns = lambda):\n")
  print(round(m, 5L))
  cat("\nObjective (quadratic surrogate + L1) at each lambda:\n")
  print(round(x$objective, 5L))
  invisible(x)
}
