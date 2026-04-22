#' @title Proper LASSO via client-side proximal gradient on normal equations
#' @description Cierra el gap OLS-soft-threshold → proper-LASSO usando solo
#'   cantidades ya expuestas por ds.vertGLM (β̂, covariance, n), sin añadir
#'   ninguna ronda MPC adicional.
#'
#'   Normal equations: para Gaussian \eqn{y = X\beta + \varepsilon} el
#'   minimizador LASSO es
#'     \deqn{\beta^* = \arg\min_{\beta} \tfrac{1}{2n} \|y - X\beta\|^2 + \lambda \|\beta\|_1.}
#'   El gradiente de la parte cuadrática es
#'     \deqn{\nabla f(\beta) = \tfrac{1}{n} X^\top X \cdot (\beta - \hat\beta_{OLS}).}
#'   ds.vertGLM ya expone \eqn{\hat\beta_{OLS}} y \eqn{\mathrm{Cov}(\hat\beta) = \sigma^2 (X^\top X)^{-1}}.
#'   De ahí el cliente reconstruye
#'     \eqn{X^\top X / n = \mathrm{Cov}^{-1} \cdot \hat\sigma^2 / n}
#'   e itera proximal-gradient puramente en memoria:
#'     \eqn{\beta_{t+1} = S_{\lambda/L}(\beta_t - \eta \nabla f(\beta_t))}
#'   con \eqn{L \ge \lambda_{\max}(X^\top X / n)} (upper bound local) y
#'   \eqn{S_t(x) = \mathrm{sign}(x) \max(|x|-t, 0)} el operador soft-threshold.
#'
#' @param fit A \code{ds.glm} object from \code{ds.vertGLM(family="gaussian")}.
#'   Must have \code{$covariance} and \code{$n_obs} populated.
#' @param lambda Numeric. L1 penalty magnitude (on the 1/n-normalised objective).
#' @param max_iter Integer. Outer proximal-gradient iterations (default 200).
#' @param tol Numeric. Convergence tolerance on \eqn{\|\beta_{t+1} - \beta_t\|}
#'   (default 1e-7).
#' @param keep_intercept Logical. If TRUE, do NOT penalise the intercept.
#' @param warm_start Numeric vector. Optional β_0 (default = β̂_OLS).
#'
#' @return An object of class \code{ds.vertLASSOProximal} with the proximal-MLE
#'   coefficients, number of iterations, convergence flag, support, final
#'   objective value, and the reconstructed Gram matrix used. The slot
#'   \code{$comparison$coefficients_soft} reports the naive post-hoc
#'   soft-thresholded OLS for comparison.
#'
#' @section P3 disclosure budget:
#'   Zero new MPC rounds beyond the initial \code{ds.vertGLM} call. All
#'   iteration is client-side on quantities already in the \code{fit}
#'   object. The optimisation volume argument is therefore moot — no
#'   additional reveals happen regardless of iteration count.
#'
#' @seealso \code{\link{ds.vertLASSO}}, \code{\link{ds.vertLASSOCV}}
#' @export
ds.vertLASSOProximal <- function(fit, lambda,
                                 max_iter = 200L, tol = 1e-7,
                                 keep_intercept = TRUE,
                                 warm_start = NULL) {
  if (!inherits(fit, "ds.glm")) {
    stop("`fit` must be a ds.glm object", call. = FALSE)
  }
  if (is.null(fit$covariance) || is.null(fit$n_obs)) {
    stop("fit must have $covariance and $n_obs", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0) {
    stop("lambda must be a non-negative scalar", call. = FALSE)
  }

  beta_ols <- as.numeric(fit$coefficients)
  names(beta_ols) <- names(fit$coefficients)
  p <- length(beta_ols)
  int_idx <- if (keep_intercept) which(names(beta_ols) == "(Intercept)") else integer(0L)

  # Reconstruct X^T X / n = Cov^{-1} * σ² / n
  # For Gaussian OLS: Cov(β̂) = σ² (X^T X)^{-1}, so (X^T X) = σ² Cov^{-1}
  # and (X^T X)/n = σ² Cov^{-1} / n. σ² is estimated as deviance/(n-p).
  sigma2_hat <- if (!is.null(fit$deviance) && is.finite(fit$deviance))
    fit$deviance / max(fit$n_obs - p, 1L) else 1
  cov_inv <- tryCatch(solve(fit$covariance),
                      error = function(e)
                        solve(fit$covariance +
                              1e-8 * diag(p) * mean(abs(diag(fit$covariance)))))
  XtX_over_n <- sigma2_hat * cov_inv / fit$n_obs

  # Step size: 1/L where L ≥ λ_max(XtX/n). Use power iteration upper bound.
  L <- max(abs(eigen(XtX_over_n, only.values = TRUE)$values))
  if (!is.finite(L) || L <= 0) L <- max(abs(diag(XtX_over_n))) * p
  eta <- 1 / L

  soft <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

  beta <- if (is.null(warm_start)) beta_ols else as.numeric(warm_start)
  names(beta) <- names(beta_ols)
  converged <- FALSE
  final_iter <- max_iter
  for (t in seq_len(max_iter)) {
    # ∇f(β) = (X^T X / n) · (β - β̂_OLS)
    grad <- as.numeric(XtX_over_n %*% (beta - beta_ols))
    beta_new <- beta - eta * grad
    # Soft-threshold (skip intercept if keep_intercept)
    beta_new <- soft(beta_new, eta * lambda)
    if (length(int_idx) == 1L) beta_new[int_idx] <- beta[int_idx] - eta * grad[int_idx]
    if (max(abs(beta_new - beta)) < tol) {
      converged <- TRUE
      final_iter <- t
      beta <- beta_new
      break
    }
    beta <- beta_new
  }

  # Final objective: 0.5 * (β - β_OLS)^T (XtX/n) (β - β_OLS) + λ||β||_1
  resid <- beta - beta_ols
  quad_part <- 0.5 * as.numeric(t(resid) %*% XtX_over_n %*% resid)
  l1_part <- lambda * sum(abs(if (length(int_idx)) beta[-int_idx] else beta))
  obj <- quad_part + l1_part

  # Comparison: naive post-hoc soft-threshold
  beta_soft <- soft(beta_ols, lambda)
  if (length(int_idx)) beta_soft[int_idx] <- beta_ols[int_idx]
  names(beta_soft) <- names(beta_ols)

  out <- list(
    coefficients = beta,
    coefficients_ols = beta_ols,
    support = which(abs(beta) > tol),
    lambda = lambda,
    converged = converged,
    iterations = final_iter,
    L = L, step_size = eta,
    XtX_over_n = XtX_over_n,
    sigma2_hat = sigma2_hat,
    objective = obj,
    n_obs = fit$n_obs,
    comparison = list(coefficients_soft = beta_soft))
  class(out) <- c("ds.vertLASSOProximal", "list")
  out
}

#' @export
print.ds.vertLASSOProximal <- function(x, ...) {
  cat(sprintf("dsVert proper LASSO (proximal gradient on normal equations)\n"))
  cat(sprintf("  N = %d  λ = %.4g  L = %.4g  converged = %s (iter %d)\n",
              x$n_obs, x$lambda, x$L, x$converged, x$iterations))
  cat(sprintf("  objective = %.6g  |support| = %d / %d\n",
              x$objective, length(x$support), length(x$coefficients)))
  df <- data.frame(
    OLS     = round(x$coefficients_ols, 5),
    softthr = round(x$comparison$coefficients_soft, 5),
    proximal= round(x$coefficients, 5),
    check.names = FALSE)
  print(df)
  invisible(x)
}
