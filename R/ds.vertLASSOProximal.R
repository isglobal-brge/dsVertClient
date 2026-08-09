#' @title Gaussian LASSO via client-side coordinate descent
#' @description Solve the Gaussian LASSO objective entirely client-side. The
#'   preferred route accepts a validated \code{ds.vertDPGaussian} release and
#'   optimises its signed, projected DP sufficient statistics without another
#'   DSI call or privacy cost. The historical \code{ds.glm} route is retained
#'   for compatibility with already-created fits.
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
#' @param fit Preferably a \code{ds.vertDPGaussian} object. A historical
#'   unpenalised Gaussian \code{ds.glm} object with \code{$covariance} and
#'   \code{$n_obs} remains accepted for compatibility.
#' @param lambda Numeric. L1 penalty magnitude (on the 1/n-normalised objective).
#' @param max_iter Positive integer. Coordinate-descent passes (default 2000).
#' @param tol Numeric. Convergence tolerance on \eqn{\|\beta_{t+1} - \beta_t\|}
#'   (default 1e-9).
#' @param keep_intercept Logical. If TRUE, do NOT penalise the intercept.
#' @param warm_start Numeric vector. For a \code{ds.vertDPGaussian} input this
#'   is an optional original-scale coefficient vector containing the intercept
#'   and every signed predictor. For a legacy \code{ds.glm}, it is beta_0.
#' @param accelerate Compatibility argument retained for historical callers.
#'   The current solver is coordinate descent, so this value has no effect.
#'
#' @return An object of class \code{ds.vertLASSOProximal} with the proximal-MLE
#'   coefficients, number of iterations, convergence flag, support, final
#'   objective value, and the reconstructed Gram matrix used. The slot
#'   \code{$comparison$coefficients_soft} reports the naive post-hoc
#'   soft-thresholded OLS for comparison.
#'
#' @section Disclosure budget:
#'   Zero new DSI or MPC rounds. For a \code{ds.vertDPGaussian} input, every
#'   lambda and iteration is deterministic post-processing of the same sticky
#'   release and has additional privacy cost \code{(epsilon, delta) = (0, 0)}.
#'   No classical sampling inference or coefficient confidence region is
#'   implied by the projected noisy moments.
#'
#' @seealso \code{\link{ds.vertLASSO}}, \code{\link{ds.vertLASSOCV}}
#' @export
ds.vertLASSOProximal <- function(fit, lambda,
                                 max_iter = 2000L, tol = 1e-9,
                                 keep_intercept = TRUE,
                                 warm_start = NULL,
                                 accelerate = TRUE) {
  if (inherits(fit, "ds.vertDPGaussian")) {
    return(.dsvert_lasso_dp_proximal(
      fit = fit, lambda = lambda, max_iter = max_iter, tol = tol,
      keep_intercept = keep_intercept, warm_start = warm_start,
      accelerate = accelerate))
  }
  if (!inherits(fit, "ds.glm")) {
    stop("`fit` must be a ds.vertDPGaussian or ds.glm object",
         call. = FALSE)
  }
  if (!is.character(fit$family) || length(fit$family) != 1L ||
      is.na(fit$family) || !identical(tolower(fit$family), "gaussian")) {
    stop("fit must use family='gaussian'; this solver implements only the Gaussian LASSO objective",
         call. = FALSE)
  }
  if (is.null(fit$lambda) || !is.numeric(fit$lambda) ||
      length(fit$lambda) != 1L || !is.finite(fit$lambda) ||
      fit$lambda != 0) {
    stop("fit must be unpenalised (fit$lambda = 0); refit the Gaussian model with lambda = 0",
         call. = FALSE)
  }
  if (is.null(fit$covariance) || is.null(fit$n_obs)) {
    stop("fit must have $covariance and $n_obs", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one finite non-negative number", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1 || max_iter != floor(max_iter)) {
    stop("max_iter must be one positive integer", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive number", call. = FALSE)
  }
  if (!is.logical(keep_intercept) || length(keep_intercept) != 1L ||
      is.na(keep_intercept)) {
    stop("keep_intercept must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(accelerate) || length(accelerate) != 1L ||
      is.na(accelerate)) {
    stop("accelerate must be TRUE or FALSE", call. = FALSE)
  }

  if (!is.numeric(fit$coefficients)) {
    stop("fit$coefficients must be a finite, uniquely named numeric vector",
         call. = FALSE)
  }
  beta_ols <- as.numeric(fit$coefficients)
  names(beta_ols) <- names(fit$coefficients)
  if (!length(beta_ols) || any(!is.finite(beta_ols)) ||
      is.null(names(beta_ols)) || any(!nzchar(names(beta_ols))) ||
      anyDuplicated(names(beta_ols))) {
    stop("fit$coefficients must be a finite, uniquely named numeric vector",
         call. = FALSE)
  }
  p <- length(beta_ols)
  intercept_idx <- which(names(beta_ols) == "(Intercept)")
  unpenalised_idx <- if (keep_intercept) intercept_idx else integer(0L)
  if (!is.numeric(fit$n_obs) || length(fit$n_obs) != 1L ||
      !is.finite(fit$n_obs) || fit$n_obs < 1 ||
      fit$n_obs != floor(fit$n_obs)) {
    stop("fit$n_obs must be one positive integer", call. = FALSE)
  }
  n_obs <- as.integer(fit$n_obs)
  covariance <- fit$covariance
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      !identical(dim(covariance), c(p, p)) || any(!is.finite(covariance))) {
    stop("fit$covariance must be a finite numeric p-by-p matrix",
         call. = FALSE)
  }
  if (!is.null(rownames(covariance)) || !is.null(colnames(covariance))) {
    if (is.null(rownames(covariance)) || is.null(colnames(covariance)) ||
        !setequal(rownames(covariance), names(beta_ols)) ||
        !setequal(colnames(covariance), names(beta_ols))) {
      stop("fit$covariance dimnames must match fit$coefficients",
           call. = FALSE)
    }
    covariance <- covariance[names(beta_ols), names(beta_ols), drop = FALSE]
  } else {
    dimnames(covariance) <- list(names(beta_ols), names(beta_ols))
  }
  covariance <- (covariance + t(covariance)) / 2

  if (!is.null(warm_start)) {
    if (!is.numeric(warm_start) || length(warm_start) != p ||
        any(!is.finite(warm_start))) {
      stop("warm_start must be NULL or a finite numeric vector with one value per coefficient",
           call. = FALSE)
    }
    if (!is.null(names(warm_start))) {
      if (any(!nzchar(names(warm_start))) || anyDuplicated(names(warm_start)) ||
          !setequal(names(warm_start), names(beta_ols))) {
        stop("named warm_start values must match fit$coefficients",
             call. = FALSE)
      }
      warm_start <- warm_start[names(beta_ols)]
    }
  }

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
  # Fallback for ds.glm-like inputs without $hessian_std: invert the
  # information covariance directly when available, or reconstruct via
  # sigma^2_epsilon * Cov(beta)^{-1} / n.
  sigma2_hat <- if (!is.null(fit$deviance) &&
                    is.numeric(fit$deviance) &&
                    length(fit$deviance) == 1L &&
                    is.finite(fit$deviance) && fit$deviance >= 0) {
    fit$deviance / max(n_obs - p, 1L)
  } else {
    NA_real_
  }

  if (!is.null(fit$hessian_std) && !is.null(fit$x_means) &&
      !is.null(fit$x_sds) && is.matrix(fit$hessian_std) &&
      all(dim(fit$hessian_std) == c(p, p))) {
    lambda_ridge <- if (!is.null(fit$lambda) && is.finite(fit$lambda))
      fit$lambda else 0
    H_std <- fit$hessian_std - lambda_ridge * diag(p)
    if (!is.numeric(H_std) || any(!is.finite(H_std))) {
      stop("fit$hessian_std must contain only finite numeric values",
           call. = FALSE)
    }
    # hessian_std may be in fit-internal (server-partition) column order,
    # which differs from the formula-order used by beta_ols / x_means /
    # x_sds. Permute to match beta_ols ordering when dimnames are set.
    if (!is.null(rownames(H_std)) || !is.null(colnames(H_std))) {
      if (is.null(rownames(H_std)) || is.null(colnames(H_std)) ||
          !setequal(rownames(H_std), names(beta_ols)) ||
          !setequal(colnames(H_std), names(beta_ols))) {
        stop("fit$hessian_std dimnames must match fit$coefficients",
             call. = FALSE)
      }
      H_std <- H_std[names(beta_ols), names(beta_ols), drop = FALSE]
    }
    x_means_vec <- as.numeric(fit$x_means[names(beta_ols)])
    x_sds_vec <- as.numeric(fit$x_sds[names(beta_ols)])
    slope_idx <- setdiff(seq_len(p), intercept_idx)
    if (any(!is.finite(x_means_vec[slope_idx])) ||
        any(!is.finite(x_sds_vec[slope_idx])) ||
        any(x_sds_vec[slope_idx] <= 0)) {
      stop("fit$x_means and fit$x_sds must cover every slope with finite values and positive scales",
           call. = FALSE)
    }
    # Intercept row: x_j = sum(X_ij)/n for slope columns; 1 for
    # intercept. Slope-slope block: x_j x_k + x_sd_j x_sd_k H_std[j,k].
    XtX_over_n <- matrix(0, p, p,
                          dimnames = list(names(beta_ols), names(beta_ols)))
    int_j <- if (length(intercept_idx) == 1L) intercept_idx else NA_integer_
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
  } else if (!is.null(fit$covariance_information)) {
    cov_info <- fit$covariance_information
    if (!is.matrix(cov_info) || !is.numeric(cov_info) ||
        !identical(dim(cov_info), c(p, p)) || any(!is.finite(cov_info))) {
      stop("fit$covariance_information must be a finite numeric p-by-p matrix",
           call. = FALSE)
    }
    if (!is.null(rownames(cov_info)) || !is.null(colnames(cov_info))) {
      if (is.null(rownames(cov_info)) || is.null(colnames(cov_info)) ||
          !setequal(rownames(cov_info), names(beta_ols)) ||
          !setequal(colnames(cov_info), names(beta_ols))) {
        stop("fit$covariance_information dimnames must match fit$coefficients",
             call. = FALSE)
      }
      cov_info <- cov_info[names(beta_ols), names(beta_ols), drop = FALSE]
    }
    cov_info <- (cov_info + t(cov_info)) / 2
    XtX_over_n <- .dsvert_solve_identifiable(
      cov_info,
      context = "The Gaussian LASSO source information",
      reason = "rank_deficient_lasso_design",
      symmetric = TRUE) / n_obs
  } else {
    if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
      stop("a finite positive Gaussian deviance is required when the fit exposes neither hessian_std nor covariance_information",
           call. = FALSE)
    }
    cov_inv <- .dsvert_solve_identifiable(
      covariance,
      context = "The Gaussian LASSO source covariance",
      reason = "rank_deficient_lasso_design",
      symmetric = TRUE)
    XtX_over_n <- sigma2_hat * cov_inv / n_obs
  }

  XtX_over_n <- (XtX_over_n + t(XtX_over_n)) / 2
  gram_values <- eigen(XtX_over_n, symmetric = TRUE,
                       only.values = TRUE)$values
  gram_tolerance <- 1e-8 * max(1, max(abs(gram_values)))
  if (min(gram_values) < -gram_tolerance) {
    stop("the reconstructed Gaussian Gram matrix is not positive semidefinite",
         call. = FALSE)
  }
  if (min(gram_values) <= gram_tolerance) {
    .dsvert_stop_non_identifiable(
      paste0(
        "The Gaussian LASSO coefficient vector is not uniquely identifiable ",
        "under the reconstructed rank-deficient design. No elastic-net or ",
        "ridge fallback was applied."),
      reason = "rank_deficient_lasso_design")
  }
  # Retained in the result for compatibility with previous output objects.
  L <- max(gram_values)
  if (!is.finite(L) || L <= 0) {
    stop("the reconstructed Gaussian Gram matrix has no positive curvature",
         call. = FALSE)
  }
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
      if (j %in% unpenalised_idx) {
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
  penalised_idx <- setdiff(seq_along(beta), unpenalised_idx)
  l1_part <- lambda * sum(abs(beta[penalised_idx]))
  obj <- quad_part + l1_part

  # Comparison: naive post-hoc soft-threshold
  beta_soft <- soft(beta_ols, lambda)
  if (length(unpenalised_idx)) {
    beta_soft[unpenalised_idx] <- beta_ols[unpenalised_idx]
  }
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
    n_obs = n_obs,
    family = "gaussian",
    estimand = "gaussian_lasso",
    solver = "coordinate_descent_normal_equations",
    source_fit_penalty = fit$lambda,
    compatibility = list(accelerate_argument = accelerate,
                         accelerate_effective = FALSE),
    comparison = list(coefficients_soft = beta_soft))
  class(out) <- c("ds.vertLASSOProximal", "list")
  out
}

#' @export
print.ds.vertLASSOProximal <- function(x, ...) {
  cat("dsVert Gaussian LASSO (coordinate descent on normal equations)\n")
  cat(sprintf("  N = %.6g  lambda = %.4g  L = %.4g  converged = %s (iter %d)\n",
              x$n_obs, x$lambda, x$L, x$converged, x$iterations))
  cat(sprintf("  objective = %.6g  |support| = %d / %d\n",
              x$objective, length(x$support),
              if (!is.null(x$coefficients_normalized)) {
                length(x$coefficients_normalized)
              } else {
                length(x$coefficients)
              }))
  if (identical(x$input_provenance, "signed_dp_gaussian_capsule")) {
    cat("  source = validated sticky Gaussian DP capsule; extra DP cost = (0, 0)\n")
    cat("\nOriginal-scale coefficients:\n")
    print(round(x$coefficients, 5L))
    cat("\nNormalized coefficients and KKT certificate:\n")
    print(round(x$coefficients_normalized, 5L))
    cat(sprintf("  max KKT violation = %.6g (tolerance %.6g)\n",
                x$kkt$max_violation, x$kkt$tolerance))
    return(invisible(x))
  }
  df <- data.frame(
    OLS     = round(x$coefficients_ols, 5),
    softthr = round(x$comparison$coefficients_soft, 5),
    proximal= round(x$coefficients, 5),
    check.names = FALSE)
  print(df)
  invisible(x)
}
