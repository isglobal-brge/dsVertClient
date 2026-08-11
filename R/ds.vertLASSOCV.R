#' @title Client-side information-criterion selection for Gaussian LASSO
#' @description Select an L1 penalty entirely client-side. With a
#'   \code{ds.vertDPGaussian} input, each candidate is solved from the same
#'   validated projected DP moments and ranked by an explicitly labelled DP
#'   projected pseudo-AIC/BIC/EBIC. With a historical \code{ds.glm} input,
#'   the existing quadratic-surrogate selector is retained. Despite the
#'   historical function/class suffix \code{CV}, neither route performs
#'   cross-validation or resampling.
#'
#'   For a DP source, the score uses the projected noisy residual moment and
#'   noisy effective count. It is therefore a deterministic model-selection
#'   heuristic, not a classical likelihood information criterion; selection
#'   is returned as unavailable when those released quantities cannot define
#'   it. For a legacy source, the score is the quadratic-surrogate misfit plus
#'   the requested degrees-of-freedom penalty. The preferred metadata name
#'   for the more parsimonious alternative is \code{lambda.parsimonious}; the
#'   historical \code{lambda.1se} slot is retained as an exact alias, but it
#'   is not a one-standard-error rule. It uses a relative IC tolerance set by
#'   \code{se_threshold} and involves no estimated sampling standard error.
#'
#'   Candidate paths are reusable post-processing and never trigger K-fold
#'   refitting or repeated private releases.
#'
#' @param fit Preferably a \code{ds.vertDPGaussian} object. A historical
#'   \code{ds.glm} object with \code{fit$covariance} remains accepted.
#' @param lambda_grid Numeric vector of candidate lambda values
#'   (default: a 50-point log-spaced grid from \code{lambda_max} to
#'   \code{lambda_max / 1000}).
#' @param criterion One of \code{"BIC"} (default), \code{"AIC"}, or
#'   \code{"EBIC"}.
#' @param ebic_gamma Extended-BIC gamma parameter (default 0.5;
#'   effective only when \code{criterion = "EBIC"}).
#' @param keep_intercept Never penalise the intercept.
#' @param se_threshold Relative IC tolerance for the parsimonious selection:
#'   retain the sparsest lambda whose IC is no more than
#'   \code{abs(IC_min) * se_threshold} above \code{IC_min}. The historical
#'   argument name is retained for compatibility; it is not a standard error.
#' @param trusted_pinset Optional trusted named peer-to-Ed25519-public-key map
#'   used to authenticate a saved or offline \code{ds.vertDPGaussian} release.
#'   Online releases created in the current session use their transport-
#'   anchored pinset cache. It is invalid for the legacy \code{ds.glm} route.
#' @return A \code{ds.vertLASSOCV} object: \code{lambda}, \code{ic},
#'   \code{df}, \code{lambda.min}, \code{lambda.parsimonious},
#'   compatibility aliases \code{lambda.1se}/\code{beta.1se}, and metadata
#'   explicitly identifying information-criterion selection without
#'   cross-validation.
#' @export
ds.vertLASSOCV <- function(fit, lambda_grid = NULL,
                           criterion = c("BIC", "AIC", "EBIC"),
                           ebic_gamma = 0.5,
                           keep_intercept = TRUE,
                           se_threshold = 0.02,
                           trusted_pinset = NULL) {
  criterion <- match.arg(criterion)
  if (inherits(fit, "ds.vertDPGaussian")) {
    return(.dsvert_lasso_dp_select(
      fit = fit, lambda_grid = lambda_grid, criterion = criterion,
      ebic_gamma = ebic_gamma, keep_intercept = keep_intercept,
      se_threshold = se_threshold, trusted_pinset = trusted_pinset))
  }
  if (!is.null(trusted_pinset)) {
    stop("trusted_pinset applies only to a ds.vertDPGaussian fit",
         call. = FALSE)
  }
  if (!inherits(fit, "ds.glm")) {
    stop("fit must be a ds.vertDPGaussian or ds.glm object",
         call. = FALSE)
  }
  if (is.null(fit$covariance)) {
    stop("fit does not expose covariance; refit with dsVert >= 8bb7902",
         call. = FALSE)
  }

  cov <- as.matrix(fit$covariance)
  H <- .dsvert_solve_identifiable(
    cov,
    context = "The LASSO-CV source covariance",
    reason = "singular_lasso_source_covariance",
    symmetric = TRUE)
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
    lambda.parsimonious = lambda.1se,
    lambda.1se = lambda.1se,
    beta.min = beta.min,
    beta.parsimonious = beta.1se,
    beta.1se = beta.1se,
    selection_method = "information_criterion_quadratic_surrogate",
    cross_validation = FALSE,
    one_standard_error_rule = FALSE,
    relative_ic_tolerance = se_threshold,
    paths = paths,
    fit = fit)
  class(out) <- c("ds.vertLASSOCV", "list")
  out
}

#' @export
print.ds.vertLASSOCV <- function(x, ...) {
  cat("dsVert information-criterion LASSO path selection (no cross-validation)\n")
  cat(sprintf("  Criterion : %s\n", x$criterion))
  if (identical(x$selection_available, FALSE)) {
    cat("  Selection unavailable: ", x$selection_unavailable_reason,
        "\n", sep = "")
    cat("  No classical IC, CV, or sampling-inference claim is made.\n")
    return(invisible(x))
  }
  cat(sprintf("  lambda.min: %.4g  (IC = %.4g, df = %d)\n",
              x$lambda.min, min(x$ic),
              x$df[which.min(x$ic)]))
  tolerance <- x$relative_ic_tolerance %||% 0.02
  cat(sprintf("  lambda.parsimonious: %.4g  (sparsest within %g%% of IC_min)\n",
              x$lambda.parsimonious %||% x$lambda.1se,
              100 * tolerance))
  cat("  Legacy alias lambda.1se is retained; this is not a sampling standard-error rule.\n")
  if (identical(x$input_provenance, "signed_dp_gaussian_capsule")) {
    cat("  Scores are DP-projected pseudo-IC values; extra DP cost = (0, 0).\n")
  }
  cat("\nbeta at lambda.min:\n")
  print(round(x$beta.min, 5L))
  if (!isTRUE(all.equal(x$lambda.min, x$lambda.1se))) {
    cat("\nbeta at lambda.1se:\n")
    print(round(x$beta.1se, 5L))
  }
  invisible(x)
}
