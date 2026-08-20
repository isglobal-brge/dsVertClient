#' @title Gaussian LASSO path from one signed DP Synopsis
#' @description Fits a deterministic, KKT-checked Gaussian LASSO path from
#'   the projected moments in one validated sticky \code{ds.vertDPGaussian}
#'   release. It performs no new DataSHIELD, MPC or privacy operation. The
#'   retired \code{ds.glm} soft-thresholding sketch is rejected rather than
#'   being represented as an L1-optimised estimator.
#'
#' @param fit A validated \code{ds.vertDPGaussian} object.
#' @param lambda_1 Base non-negative L1 penalty on the signed-bound normalized
#'   Gaussian objective.
#' @param alpha_grid Non-negative multipliers defining
#'   \code{lambda_1 * alpha_grid} (default 1, 0.5, 0.25, 0.125, 0.0625).
#' @param keep_intercept If TRUE (default), the intercept is not penalised.
#' @param max_iter Positive bound on coordinate-descent iterations per path
#'   member.
#' @param tol Positive KKT and convergence tolerance.
#' @param trusted_pinset Optional trusted named peer-to-Ed25519-public-key map
#'   required to authenticate an offline saved Gaussian release; online
#'   session-anchored releases use their cached pinset.
#' @return A \code{ds.vertLASSO} / \code{ds.vertDPLASSOPath} object with one
#'   KKT certificate per penalty, original- and normalized-scale paths and no
#'   sampling inference or coefficient confidence regions.
#' @section Disclosure budget:
#' Every path member is deterministic post-processing of one authenticated
#' sticky Synopsis: additional privacy cost is \code{(epsilon, delta) =
#' (0, 0)} and no additional server calls occur.
#' @export
ds.vertLASSO <- function(fit, lambda_1,
                         alpha_grid = c(1, 0.5, 0.25, 0.125, 0.0625),
                         keep_intercept = TRUE,
                         max_iter = 2000L, tol = 1e-9,
                         trusted_pinset = NULL) {
  if (!inherits(fit, "ds.vertDPGaussian")) {
    stop(
      "fit must be a validated ds.vertDPGaussian release; the retired ds.glm thresholding sketch is unavailable",
      call. = FALSE)
  }
  if (!is.numeric(lambda_1) || length(lambda_1) != 1L ||
      !is.finite(lambda_1) || lambda_1 < 0) {
    stop("lambda_1 must be one finite non-negative number", call. = FALSE)
  }
  if (!is.numeric(alpha_grid) || !length(alpha_grid) ||
      any(!is.finite(alpha_grid)) || any(alpha_grid < 0)) {
    stop("alpha_grid must contain finite non-negative values", call. = FALSE)
  }
  if (!is.logical(keep_intercept) || length(keep_intercept) != 1L ||
      is.na(keep_intercept)) {
    stop("keep_intercept must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1L || max_iter != floor(max_iter)) {
    stop("max_iter must be one positive integer", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive number", call. = FALSE)
  }

  lambda_grid <- as.numeric(lambda_1) * as.numeric(alpha_grid)
  source <- .dsvert_lasso_dp_source(fit, trusted_pinset = trusted_pinset)
  if (!identical(source$artifact$implementation_state,
                 "same_owner_materialized") ||
      !identical(source$artifact$cross_owner_state,
                 "reserved_not_materialized")) {
    stop(
      "ds.vertLASSO currently requires a validated same-owner Gaussian Synopsis",
      call. = FALSE)
  }
  solutions <- .dsvert_lasso_dp_lambda_path(
    source$moments$gram, source$moments$cross, lambda_grid,
    max_iter = as.integer(max_iter), tol = as.numeric(tol),
    keep_intercept = keep_intercept)
  paths_normalized <- lapply(solutions, `[[`, "coefficients")
  paths <- lapply(paths_normalized, .dsvert_lasso_dp_original_scale,
                  artifact = source$artifact)
  names(paths) <- names(paths_normalized) <- names(solutions)
  zero_index <- which(lambda_grid == 0)
  original_normalized <- if (length(zero_index)) {
    paths_normalized[[zero_index[[1L]]]]
  } else {
    .dsvert_lasso_dp_solver(
      source$moments$gram, source$moments$cross, lambda = 0,
      max_iter = as.integer(max_iter), tol = as.numeric(tol),
      keep_intercept = keep_intercept)$coefficients
  }

  out <- list(
    lambda_grid = lambda_grid,
    paths = paths,
    paths_normalized = paths_normalized,
    path_certificates = lapply(solutions, function(solution) {
      list(kkt = solution$kkt, curvature = solution$curvature)
    }),
    original = .dsvert_lasso_dp_original_scale(
      original_normalized, source$artifact),
    fit = fit,
    family = "gaussian",
    estimand = "bounded_normalized_DP_gaussian_lasso_path",
    input_provenance = "signed_dp_gaussian_capsule",
    source_certificate_validation = source$verification,
    source_query_contract_sha256 = fit$query_contract_sha256,
    source_release_contract_hash = fit$release_contract_hash,
    inference = list(
      classical_standard_errors = NULL, p_values = NULL,
      confidence_intervals = NULL, sampling_inference_available = FALSE),
    coefficient_regions_available = FALSE,
    additional_server_calls_after_capsule = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    disclosure_guard = list(
      satisfied = TRUE,
      basis = "deterministic_postprocessing_of_one_validated_DP_capsule"))
  if (identical(source$verification$certificate$version,
                .DSVERT_DP_GAUSSIAN_SYNOPSIS_CERTIFICATE_VERSION)) {
    out$input_provenance <- "signed_dp_gaussian_synopsis"
    out$source_query_contract_sha256 <- NULL
    out$source_release_contract_hash <- NULL
    out$source_contract_sha256 <- fit$contract_sha256
    out$additional_server_calls_after_capsule <- NULL
    out$additional_server_calls_after_synopsis <- 0L
    out$disclosure_guard$basis <-
      "deterministic_postprocessing_of_one_validated_DP_synopsis"
  }
  class(out) <- c("ds.vertLASSO", "ds.vertDPLASSOPath", "list")
  out
}

#' @export
print.ds.vertLASSO <- function(x, ...) {
  cat("dsVert Gaussian LASSO path from a signed DP Synopsis\n")
  cat(sprintf("  Penalty grid: %s\n",
              paste(sprintf("%.4f", x$lambda_grid), collapse = " ")))
  cat("  Extra DP cost: (0, 0); no sampling inference\n\n")
  cat("Coefficient paths:\n")
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
#'     - A fast surrogate baseline for comparison with the product
#'       proximal / iterative LASSO routes.
#'
#' @param fit      A \code{ds.glm} object with \code{fit$covariance}
#'                 populated.
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
  H <- .dsvert_solve_identifiable(
    cov,
    context = "The LASSO source covariance",
    reason = "singular_lasso_source_covariance",
    symmetric = TRUE)
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
