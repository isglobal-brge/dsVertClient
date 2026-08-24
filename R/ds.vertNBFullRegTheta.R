#' @title Sticky-DP intercept-only NB2 frontdoor
#' @description With a released, validated \code{ds.vertDPFrequency} object,
#'   or an explicit \code{server} from which it can be read,
#'   whose signed domain is bounded non-negative integer counts, this frontdoor
#'   fits \code{y ~ 1} by deterministic NB2 method-of-moments post-processing.
#'   It never starts a new analysis or DataSHIELD request.
#' @details This is deliberately narrower than NB2 regression. Predictors,
#'   likelihood optimization, covariance, standard errors and inference remain
#'   unavailable until a protected beta/theta score-and-information artifact is
#'   implemented. A frequency release with variance no greater than its mean is
#'   reported as the Poisson-limit boundary rather than as a finite NB2 fit.
#' @param formula,data,theta,joint,theta_max_iter,theta_tol,variant,beta_max_iter,beta_tol,compute_covariance,verbose,datasources
#'   Retained compatibility arguments. With \code{frequency}, only
#'   \code{formula}, \code{data}, \code{verbose}, and \code{frequency} are
#'   accepted; all legacy NB2 controls fail locally.
#' @param server Required source-owner when \code{frequency} is absent. The
#'   frontdoor reads the canonical signed Frequency artifact for the outcome.
#' @param frequency A released, validated \code{ds.vertDPFrequency} object for
#'   a bounded count outcome. It enables only an intercept-only, no-inference
#'   NB2 method-of-moments result.
#' @param ... Retained compatibility arguments; not evaluated.
#' @return With \code{frequency} or \code{server}, a coefficient-only
#'   \code{ds.vertNBFullRegTheta} object. Otherwise the function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertNBFullRegTheta <- function(formula, data = NULL, theta = NULL,
                                  joint = TRUE, theta_max_iter = 5L,
                                  theta_tol = 1e-3, variant = "full_reg_nd",
                                  beta_max_iter = 2L, beta_tol = 1e-4,
                                  compute_covariance = TRUE,
                                  verbose = TRUE, datasources = NULL, ...,
                                  frequency = NULL, server = NULL) {
  explicit_arguments <- names(match.call())[-1L]
  if (is.null(frequency) && (!is.null(server) ||
                             (!missing(formula) &&
                              .dsvert_dp_frequency_intercept_formula(formula)))) {
    frequency <- .dsvert_dp_frequency_intercept_artifact(
      formula = if (missing(formula)) NULL else formula, data = data,
      server = server, datasources = datasources, method = "NB2")
    return(.dsvert_formal_nb2_frequency_adapter(
      explicit_arguments = setdiff(explicit_arguments, c("server", "datasources")),
      formula = if (missing(formula)) NULL else formula,
      data = data, frequency = frequency))
  }
  if (!is.null(frequency)) {
    return(.dsvert_formal_nb2_frequency_adapter(
      explicit_arguments = explicit_arguments,
      formula = if (missing(formula)) NULL else formula,
      data = data, frequency = frequency))
  }
  if (!identical(variant, "full_reg_nd")) {
    stop("variant must be 'full_reg_nd'", call. = FALSE)
  }
  .dsvert_block_retired_remote_route("negative_binomial")
}
#' @export
print.ds.vertNBFullRegTheta <- function(x, ...) {
  if (inherits(x, "dsvert_dp_frequency_nb2")) {
    cat("dsVert sticky-DP NB2 intercept-only method-of-moments fit\n")
    cat(sprintf("  DP effective count = %.6g   mean = %.6g   variance = %.6g\n",
                x$effective_count_dp, x$mean, x$variance))
    if (is.infinite(x$theta)) {
      cat("  Dispersion: Poisson-limit (finite NB2 theta not identified)\n")
    } else {
      cat(sprintf("  theta = %.6g\n", x$theta))
    }
    print(round(x$coefficients, 5L))
    return(invisible(x))
  }
  cat("dsVert NB regression (full-reg theta: variance-corrected profile MLE)\n")
  cat(sprintf("  N = %d   theta = %.4g   theta_iid = %.4g   variant = %s\n",
              x$n_obs, x$theta, x$theta_iid, x$variant))
  cat(sprintf("  Var(mu) estimate = %.4g   var-inflation = %.3f\n",
              x$variance_correction, x$var_inflation))
  if (!is.null(x$quality$status)) {
    cat(sprintf("  Quality: %s\n", x$quality$status))
    if (length(x$quality$warnings)) {
      for (w in x$quality$warnings) cat("  - ", w, "\n", sep = "")
    }
  }
  df <- data.frame(
    Estimate = x$coefficients,
    SE       = x$std_errors,
    z        = x$z_values,
    p        = x$p_values,
    check.names = FALSE)
  print(round(df, 5L))
  invisible(x)
}

.dsvert_formal_nb2_frequency_adapter <- function(
    explicit_arguments, formula, data, frequency) {
  allowed <- c("formula", "data", "verbose", "frequency")
  unexpected <- setdiff(explicit_arguments, allowed)
  if (length(unexpected)) {
    stop(paste(
      "Frequency-backed NB2 does not accept legacy controls:",
      paste(sort(unexpected, method = "radix"), collapse = ", ")),
      call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]])) {
    stop("Frequency-backed NB2 requires a simple outcome formula",
         call. = FALSE)
  }
  terms <- stats::terms(formula)
  if (!identical(as.integer(attr(terms, "intercept")), 1L) ||
      length(attr(terms, "term.labels")) != 0L) {
    stop("Frequency-backed NB2 supports only an intercept-only y ~ 1 formula",
         call. = FALSE)
  }
  frequency <- .dsvert_dp_frequency_contract(frequency)
  levels <- frequency$levels
  counts <- frequency$counts
  if (!is.character(levels) || length(levels) < 2L || anyNA(levels) ||
      any(!nzchar(levels)) || anyDuplicated(levels) ||
      !is.numeric(counts) || length(counts) != length(levels) ||
      is.null(names(counts)) || !identical(names(counts), levels) ||
      any(!is.finite(counts)) || any(counts < 0)) {
    stop("Frequency-backed NB2 requires a non-empty signed count release",
         call. = FALSE)
  }
  outcome <- as.character(formula[[2L]])
  if (!identical(outcome, frequency$variable)) {
    stop("Frequency-backed NB2 outcome does not match the signed frequency variable",
         call. = FALSE)
  }
  descriptor <- frequency$coordinate_descriptor
  dataset <- if (is.list(descriptor)) descriptor$dataset else NULL
  if (!is.null(data) && (!is.character(data) || length(data) != 1L ||
                         is.na(data) || !identical(data, dataset))) {
    stop("Frequency-backed NB2 data does not match the signed frequency release",
         call. = FALSE)
  }
  if (any(!grepl("^(0|[1-9][0-9]*)$", levels)) ||
      any(nchar(levels) > 16L) ||
      any(nchar(levels) == 16L & levels > "9007199254740991")) {
    stop(paste(
      "Frequency-backed NB2 requires non-negative integer count levels",
      "no larger than 2^53 - 1"), call. = FALSE)
  }
  support <- suppressWarnings(as.numeric(levels))
  if (any(!is.finite(support))) {
    stop("Frequency-backed NB2 requires non-negative integer count levels",
         call. = FALSE)
  }
  total <- sum(counts)
  if (!is.finite(total) || total <= 0) {
    stop("Frequency-backed NB2 requires a non-empty signed count release",
         call. = FALSE)
  }
  probabilities <- counts / total
  mean_dp <- sum(probabilities * support)
  variance_dp <- sum(probabilities * (support - mean_dp)^2)
  if (!is.finite(mean_dp) || !is.finite(variance_dp) || mean_dp <= 0) {
    stop("Frequency-backed NB2 requires a positive finite DP mean", call. = FALSE)
  }
  if (variance_dp > mean_dp) {
    theta <- mean_dp^2 / (variance_dp - mean_dp)
    dispersion_status <- "overdispersed_nb2_moment"
  } else {
    theta <- Inf
    dispersion_status <- "poisson_limit"
  }
  result <- list(
    status = "public_certified_intercept_only_nb2",
    family = "negative_binomial_intercept_only_frequency_postprocessing",
    coefficients = stats::setNames(log(mean_dp), "(Intercept)"),
    theta = theta,
    dispersion_status = dispersion_status,
    mean = mean_dp,
    variance = variance_dp,
    outcome_support = stats::setNames(support, levels),
    dp_counts = counts,
    effective_count_dp = total,
    frequency_release_sha256 = frequency$release_sha256,
    sticky_noise = TRUE,
    sticky_replay = TRUE,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    covariance = NULL,
    std_errors = NULL,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    production_ready = FALSE,
    inference = "unavailable_without_a_protected_nb2_score_information_artifact",
    called_via = "ds.vertNBFullRegTheta_frequency")
  class(result) <- c("dsvert_dp_frequency_nb2", "ds.vertNBFullRegTheta", "list")
  result
}

# Null-coalescing helper -- internal, may already exist elsewhere but
# redefined here for standalone loading safety.
`%||%` <- function(a, b) if (is.null(a)) b else a
