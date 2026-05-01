#' @title Proper LASSO via client-side proximal gradient on normal equations
#' @description Cierra el gap OLS-soft-threshold -> proper-LASSO usando solo
#'   cantidades ya expuestas por ds.vertGLM (beta, covariance, n), sin anadir
#'   ninguna ronda MPC adicional.
#'
#'   Normal equations: para Gaussian \eqn{y = X\beta + \varepsilon} el
#'   minimizador LASSO es
#'     \deqn{\beta^* = \arg\min_{\beta} \tfrac{1}{2n} \|y - X\beta\|^2 + \lambda \|\beta\|_1.}
#'   El gradiente de la parte cuadratica es
#'     \deqn{\nabla f(\beta) = \tfrac{1}{n} X^\top X \cdot (\beta - \hat\beta_{OLS}).}
#'   ds.vertGLM ya expone \eqn{\hat\beta_{OLS}} y \eqn{\mathrm{Cov}(\hat\beta) = \sigma^2 (X^\top X)^{-1}}.
#'   De ahi el cliente reconstruye
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
#' @param warm_start Numeric vector. Optional beta_0 (default = beta_OLS).
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
#'   object. The optimisation volume argument is therefore moot -- no
#'   additional reveals happen regardless of iteration count.
#'
#' @seealso \code{\link{ds.vertLASSO}}, \code{\link{ds.vertLASSOCV}}
#' @export
ds.vertLASSOProximal <- function(fit, lambda,
                                 max_iter = 2000L, tol = 1e-9,
                                 keep_intercept = TRUE,
                                 warm_start = NULL,
                                 accelerate = TRUE) {
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

  # Reconstruct X^T X / n. Preferred route: use the standardized-space
  # Hessian exposed by ds.vertGLM as $hessian_std (X_std^T X_std / n +
  # lambdaI), plus $x_means, $x_sds. This avoids the double matrix
  # inversion (Cov -> inv -> rescale) that loses precision on ill-
  # conditioned designs.
  #
  # For standardized features X_std[,j] = (X_raw[,j] - x_j)/x_sd_j:
  #   G[j,k] = x_j x_k + x_sd_j * x_sd_k * H_std[j,k]   (slope j,k)
  #   G[0,j] = x_j                                       (intercept-slope)
  #   G[0,0] = 1
  #
  # Fallback for ds.glm-like inputs without $hessian_std (e.g. lm()
  # mocks): reconstruct via sigma^2_epsilon * Cov^{-1} / n (lm() convention) or
  # y_sd^2 * Cov^{-1} / n (ds.vertGLM convention).
  sigma2_hat <- if (!is.null(fit$deviance) && is.finite(fit$deviance))
    fit$deviance / max(fit$n_obs - p, 1L) else 1

  if (!is.null(fit$hessian_std) && !is.null(fit$x_means) &&
      !is.null(fit$x_sds) && is.matrix(fit$hessian_std) &&
      all(dim(fit$hessian_std) == c(p, p))) {
    lambda_ridge <- if (!is.null(fit$lambda) && is.finite(fit$lambda))
      fit$lambda else 0
    H_std <- fit$hessian_std - lambda_ridge * diag(p)
    # hessian_std may be in fit-internal (server-partition) column order,
    # which differs from the formula-order used by beta_ols / x_means /
    # x_sds. Permute to match beta_ols ordering when dimnames are set.
    if (!is.null(dimnames(H_std)) && !is.null(rownames(H_std))) {
      perm <- match(names(beta_ols), rownames(H_std))
      if (all(!is.na(perm))) {
        H_std <- H_std[perm, perm, drop = FALSE]
      }
    }
    x_means_vec <- as.numeric(fit$x_means[names(beta_ols)])
    x_sds_vec <- as.numeric(fit$x_sds[names(beta_ols)])
    # Intercept row: x_j = sum(X_ij)/n for slope columns; 1 for
    # intercept. Slope-slope block: x_j x_k + x_sd_j x_sd_k H_std[j,k].
    XtX_over_n <- matrix(0, p, p,
                          dimnames = list(names(beta_ols), names(beta_ols)))
    int_j <- if (length(int_idx) == 1L) int_idx else NA_integer_
    for (jj in seq_len(p)) {
      for (kk in seq_len(p)) {
        if (!is.na(int_j) && jj == int_j && kk == int_j) {
          XtX_over_n[jj, kk] <- 1
        } else if (!is.na(int_j) && jj == int_j) {
          XtX_over_n[jj, kk] <- x_means_vec[kk]
        } else if (!is.na(int_j) && kk == int_j) {
          XtX_over_n[jj, kk] <- x_means_vec[jj]
        } else {
          XtX_over_n[jj, kk] <- x_means_vec[jj] * x_means_vec[kk] +
            x_sds_vec[jj] * x_sds_vec[kk] * H_std[jj, kk]
        }
      }
    }
  } else {
    cov_inv <- tryCatch(solve(fit$covariance),
                        error = function(e)
                          solve(fit$covariance +
                                1e-8 * diag(p) * mean(abs(diag(fit$covariance)))))
    scale_factor <- if (!is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
      fit$y_sd^2
    } else {
      sigma2_hat
    }
    XtX_over_n <- scale_factor * cov_inv / fit$n_obs
  }

  # Step size: 1/L where L >= lambda_max(XtX/n). Use power iteration upper bound.
  L <- max(abs(eigen(XtX_over_n, only.values = TRUE)$values))
  if (!is.finite(L) || L <= 0) L <- max(abs(diag(XtX_over_n))) * p
  eta <- 1 / L

  soft <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

  beta <- if (is.null(warm_start)) beta_ols else as.numeric(warm_start)
  names(beta) <- names(beta_ols)
  converged <- FALSE
  final_iter <- max_iter
  # Coordinate descent (Friedman-Hastie-Tibshirani 2010 JSS Sec.2.4) --
  # the canonical LASSO solver (same as glmnet). Performs exact
  # minimization along one coordinate at a time, which is crucial for
  # unstandardized designs where proximal gradient has L/mu-limited
  # convergence (on NHANES-like p=4 n=132 glu in [50,200] even FISTA
  # needs >50K iters for 1e-3 accuracy; CD converges in <100 passes).
  #
  # Per-coordinate update:
  #   beta_j <- S_{lambda/G_jj}(beta_OLS_j - Sum_{k!=j} G_jk (beta_k - beta_OLS_k) / G_jj)
  # where G = XtX/n. Intercept (if keep_intercept) updated without
  # soft-threshold.
  G <- XtX_over_n
  G_diag <- diag(G)
  for (t in seq_len(max_iter)) {
    beta_old <- beta
    for (j in seq_along(beta)) {
      Gjj <- G_diag[j]
      if (!is.finite(Gjj) || Gjj < 1e-12) next
      off_j <- sum(G[j, -j] * (beta[-j] - beta_ols[-j]))
      base <- beta_ols[j] - off_j / Gjj
      if (length(int_idx) == 1L && j == int_idx) {
        beta[j] <- base
      } else {
        beta[j] <- soft(base, lambda / Gjj)
      }
    }
    if (max(abs(beta - beta_old)) < tol) {
      converged <- TRUE
      final_iter <- t
      break
    }
  }

  # Final objective: 0.5 * (beta - beta_OLS)^T (XtX/n) (beta - beta_OLS) + lambda||beta||_1
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
  cat(sprintf("  N = %d  lambda = %.4g  L = %.4g  converged = %s (iter %d)\n",
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
