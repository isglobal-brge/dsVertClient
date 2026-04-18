#' @title Information-criterion lambda selection for one-step LASSO
#' @description Select the L1 penalty on the one-step quadratic-surrogate
#'   LASSO path (\code{ds.vertLASSO1Step}) by minimising AIC / BIC /
#'   extended-BIC of the quadratic-surrogate criterion. Because the
#'   surrogate uses only the already-computed coefficient vector and the
#'   full covariance matrix, selection is entirely client-side: no new
#'   MPC rounds are spent beyond the single \code{ds.vertGLM} fit that
#'   produced \code{fit}.
#'
#'   For each lambda we solve the one-step LASSO, then score
#'       IC(lambda) = surrogate_misfit(beta_lambda) + penalty * df
#'   where \code{df} is the number of nonzero coefficients and
#'   \code{penalty} is 2 (AIC), log(n) (BIC, default), or
#'   log(n) + 2 gamma log(p) (extended BIC). The selected
#'   \code{lambda.min} is the global minimiser; \code{lambda.1se} is
#'   reported as a more parsimonious alternative that preserves at
#'   least \code{(1 - se_threshold) * IC_min} of the fit.
#'
#'   This is the standard selector for one-step / SCAD-style penalised
#'   maximum likelihood and is defensible in large samples without the
#'   need for full K-fold refitting (which would require rerunning the
#'   MPC GLM loop K times).
#'
#' @param fit A \code{ds.glm} object (with \code{fit$covariance}).
#' @param lambda_grid Numeric vector of candidate lambda values
#'   (default: a 50-point log-spaced grid from \code{lambda_max} to
#'   \code{lambda_max / 1000}).
#' @param criterion One of \code{"BIC"} (default), \code{"AIC"}, or
#'   \code{"EBIC"}.
#' @param ebic_gamma Extended-BIC gamma parameter (default 0.5;
#'   effective only when \code{criterion = "EBIC"}).
#' @param keep_intercept Never penalise the intercept.
#' @param se_threshold For \code{lambda.1se}, retain the sparsest lambda
#'   whose IC is within this fraction of \code{IC_min} (default 0.02,
#'   i.e. 2\%).
#' @return A \code{ds.vertLASSOCV} object: \code{lambda}, \code{ic},
#'   \code{df}, \code{lambda.min}, \code{lambda.1se}, \code{beta.min},
#'   \code{beta.1se}, the original fit.
#' @export
ds.vertLASSOCV <- function(fit, lambda_grid = NULL,
                           criterion = c("BIC", "AIC", "EBIC"),
                           ebic_gamma = 0.5,
                           keep_intercept = TRUE,
                           se_threshold = 0.02) {
  criterion <- match.arg(criterion)
  if (!inherits(fit, "ds.glm")) {
    stop("fit must be a ds.glm object", call. = FALSE)
  }
  if (is.null(fit$covariance)) {
    stop("fit does not expose covariance; refit with dsVert >= 8bb7902",
         call. = FALSE)
  }

  cov <- as.matrix(fit$covariance)
  H <- tryCatch(solve((cov + t(cov)) / 2), error = function(e) NULL)
  if (is.null(H)) stop("Cov(beta) is singular", call. = FALSE)
  H <- (H + t(H)) / 2
  diag_H <- diag(H)
  if (any(diag_H <= 0)) stop("Non-positive Hessian diagonal", call. = FALSE)

  betahat <- as.numeric(fit$coefficients)
  names(betahat) <- names(fit$coefficients)
  p <- length(betahat)
  int_idx <- which(names(betahat) == "(Intercept)")
  penalisable <- setdiff(seq_len(p), if (keep_intercept) int_idx else integer(0))

  soft_threshold <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

  # Default grid: scan from roughly where all penalisable coeffs are
  # zeroed down to lambda_max/1000.
  if (is.null(lambda_grid)) {
    lambda_max <- max(abs(diag_H[penalisable] * betahat[penalisable]))
    lambda_grid <- exp(seq(log(max(lambda_max, 1e-10)),
                           log(max(lambda_max, 1e-10) / 1000),
                           length.out = 50L))
  }

  solve_lambda <- function(lam) {
    beta <- betahat
    for (iter in seq_len(500L)) {
      max_delta <- 0
      for (j in seq_len(p)) {
        delta_k <- beta - betahat
        delta_k[j] <- 0
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
      if (max_delta < 1e-10) break
    }
    beta
  }

  n <- if (!is.null(fit$n_obs)) fit$n_obs else 1
  pen_const <- switch(criterion,
    AIC  = 2,
    BIC  = log(max(n, 2)),
    EBIC = log(max(n, 2)) + 2 * ebic_gamma * log(max(p, 2)))

  ic_vals <- numeric(length(lambda_grid))
  df_vals <- integer(length(lambda_grid))
  paths <- list()
  for (li in seq_along(lambda_grid)) {
    lam <- lambda_grid[li]
    beta <- solve_lambda(lam)
    delta <- beta - betahat
    surrogate <- drop(t(delta) %*% H %*% delta)  # 2 * misfit
    # Effective df = number of nonzero coords (excluding intercept).
    nz <- sum(abs(beta[penalisable]) > 1e-6) +
          (if (keep_intercept && length(int_idx) == 1L) 1L else 0L)
    df_vals[li] <- nz
    ic_vals[li] <- surrogate + pen_const * nz
    paths[[sprintf("%.6g", lam)]] <- beta
  }

  li_min <- which.min(ic_vals)
  lambda.min <- lambda_grid[li_min]
  beta.min <- paths[[li_min]]

  lambda.1se <- lambda.min
  beta.1se <- beta.min
  ic_min <- ic_vals[li_min]
  thresh <- ic_min + abs(ic_min) * se_threshold
  viable <- which(ic_vals <= thresh)
  if (length(viable) > 0L) {
    # Choose sparsest viable lambda (highest lambda value among viable).
    li_1se <- viable[which.max(lambda_grid[viable])]
    lambda.1se <- lambda_grid[li_1se]
    beta.1se <- paths[[li_1se]]
  }
  names(beta.min) <- names(beta.1se) <- names(betahat)

  out <- list(
    lambda = lambda_grid,
    ic = ic_vals,
    df = df_vals,
    criterion = criterion,
    lambda.min = lambda.min,
    lambda.1se = lambda.1se,
    beta.min = beta.min,
    beta.1se = beta.1se,
    paths = paths,
    fit = fit)
  class(out) <- c("ds.vertLASSOCV", "list")
  out
}

#' @export
print.ds.vertLASSOCV <- function(x, ...) {
  cat("dsVert information-criterion LASSO path selection\n")
  cat(sprintf("  Criterion : %s\n", x$criterion))
  cat(sprintf("  lambda.min: %.4g  (IC = %.4g, df = %d)\n",
              x$lambda.min, min(x$ic),
              x$df[which.min(x$ic)]))
  cat(sprintf("  lambda.1se: %.4g  (sparsest within %g%% of IC_min)\n",
              x$lambda.1se, 2))
  cat("\nbeta at lambda.min:\n")
  print(round(x$beta.min, 5L))
  if (!isTRUE(all.equal(x$lambda.min, x$lambda.1se))) {
    cat("\nbeta at lambda.1se:\n")
    print(round(x$beta.1se, 5L))
  }
  invisible(x)
}
