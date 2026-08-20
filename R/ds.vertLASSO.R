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
  objective <- vapply(solutions, `[[`, numeric(1L), "objective_without_constant")
  names(objective) <- names(solutions)
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
    objective = objective,
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


#' @title Gaussian LASSO path under the historical one-step name
#' @description Compatibility frontdoor for the ordinary one-step LASSO name.
#'   It now delegates to the KKT-certified path over one validated sticky
#'   \code{ds.vertDPGaussian} Synopsis. The retired \code{ds.glm}
#'   quadratic-surrogate approximation is unavailable.
#'
#' @param fit A validated \code{ds.vertDPGaussian} object.
#' @param lambda Numeric vector of finite non-negative L1 penalties.
#' @param keep_intercept If TRUE (default), the intercept is not penalised.
#' @param max_iter Positive bound on coordinate-descent iterations per path
#'   member.
#' @param tol Positive KKT and convergence tolerance.
#' @param trusted_pinset Optional trusted peer pinset for an offline release.
#' @return A \code{ds.vertLASSO1Step} / \code{ds.vertDPLASSOPath} object with
#'   original- and normalized-scale paths and KKT certificates.
#' @export
ds.vertLASSO1Step <- function(fit, lambda,
                              keep_intercept = TRUE,
                              max_iter = 500L,
                              tol = 1e-8,
                              trusted_pinset = NULL) {
  out <- ds.vertLASSO(
    fit = fit, lambda_1 = 1, alpha_grid = lambda,
    keep_intercept = keep_intercept, max_iter = max_iter, tol = tol,
    trusted_pinset = trusted_pinset)
  out$lambda <- out$lambda_grid
  class(out) <- c("ds.vertLASSO1Step", "ds.vertDPLASSOPath", "list")
  out
}

#' @export
print.ds.vertLASSO1Step <- function(x, ...) {
  cat("dsVert Gaussian LASSO path (historical one-step frontdoor)\n")
  cat(sprintf("  Penalty grid: %s\n",
              paste(sprintf("%.4g", x$lambda), collapse = " ")))
  cat("  Extra DP cost: (0, 0); no sampling inference\n\n")
  print(round(do.call(cbind, x$paths), 5L))
  invisible(x)
}
