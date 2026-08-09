.ds_epi_validate_level_scale <- function(level, scale) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) ||
      scale <= 0) {
    stop("scale must be one finite, strictly positive number", call. = FALSE)
  }
}

.ds_epi_validate_vector <- function(x, name, non_negative = TRUE,
                                    whole_number = FALSE,
                                    strictly_positive = FALSE) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop(name, " must be a non-empty numeric vector of finite values",
         call. = FALSE)
  }
  if (non_negative && any(x < 0)) {
    stop(name, " must contain non-negative values", call. = FALSE)
  }
  if (strictly_positive && any(x <= 0)) {
    stop(name, " must contain strictly positive values", call. = FALSE)
  }
  if (whole_number &&
      any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop(name, " must contain whole-number aggregate counts", call. = FALSE)
  }
  as.numeric(x)
}

.ds_epi_align_strata <- function(vectors, labels) {
  lengths <- vapply(vectors, length, integer(1L))
  if (length(unique(lengths)) != 1L) {
    stop(paste(labels, collapse = ", "), " must have the same length",
         call. = FALSE)
  }

  vector_names <- lapply(vectors, names)
  has_names <- vapply(vector_names, function(x) !is.null(x), logical(1L))
  if (any(has_names) && !all(has_names)) {
    stop(paste(labels, collapse = ", "),
         " must all be named or all be unnamed", call. = FALSE)
  }
  if (!all(has_names)) {
    return(vectors)
  }

  complete_unique <- vapply(vector_names, function(x) {
    length(x) > 0L && all(nzchar(x)) && !anyDuplicated(x)
  }, logical(1L))
  if (!all(complete_unique)) {
    stop("stratum names must be non-empty and unique", call. = FALSE)
  }
  reference <- vector_names[[1L]]
  same_sets <- vapply(vector_names[-1L], function(x) setequal(reference, x),
                      logical(1L))
  if (!all(same_sets)) {
    stop(paste(labels, collapse = ", "),
         " must contain the same stratum names", call. = FALSE)
  }
  lapply(vectors, function(x) x[reference])
}

#' Direct standardization of aggregate stratum-specific rates
#'
#' Computes a directly standardized rate using caller-supplied aggregate event
#' counts, person-time, and an external standard population. The function is
#' entirely client-side: it makes no DSI or DataSHIELD calls and exposes no new
#' server-side values.
#'
#' The default confidence interval is the Poisson-gamma interval of Fay and
#' Feuer, including its one-event upper-tail adjustment. A normal interval
#' using the independent-Poisson variance is also available; its lower bound is
#' truncated at zero and the truncation is reported in \code{correction}.
#'
#' @param cases Numeric vector of non-negative, whole-number aggregate event
#'   counts by stratum.
#' @param person_time Numeric vector of strictly positive aggregate person-time
#'   by stratum.
#' @param standard_population Numeric vector of non-negative external standard
#'   population weights. Values are normalized internally to sum to one.
#' @param scale Positive multiplier for reported rates, commonly \code{1e5}.
#' @param level Confidence level.
#' @param ci_method \code{"poisson_gamma"} (default) for the Fay--Feuer
#'   interval or \code{"normal"} for a bounded Wald interval.
#'
#' @return An object of class \code{ds.vertDirectStandardization}. Its estimate,
#'   standard error, confidence interval, and stratum rates are on the scale
#'   selected by \code{scale}. The object records aggregate provenance,
#'   assumptions, and any interval correction.
#'
#' @references Fay, M. P. and Feuer, E. J. (1997). Confidence intervals for
#'   directly standardized rates: a method based on the gamma distribution.
#'   \emph{Statistics in Medicine}, 16, 791--801.
#' @export
ds.vertDirectStandardization <- function(
    cases, person_time, standard_population, scale = 1e5, level = 0.95,
    ci_method = c("poisson_gamma", "normal")) {
  ci_method <- match.arg(ci_method)
  .ds_epi_validate_level_scale(level, scale)

  original_names <- list(names(cases), names(person_time),
                         names(standard_population))
  cases <- .ds_epi_validate_vector(
    cases, "cases", whole_number = TRUE)
  person_time <- .ds_epi_validate_vector(
    person_time, "person_time", strictly_positive = TRUE)
  standard_population <- .ds_epi_validate_vector(
    standard_population, "standard_population")
  names(cases) <- original_names[[1L]]
  names(person_time) <- original_names[[2L]]
  names(standard_population) <- original_names[[3L]]
  aligned <- .ds_epi_align_strata(
    list(cases, person_time, standard_population),
    c("cases", "person_time", "standard_population"))
  cases <- aligned[[1L]]
  person_time <- aligned[[2L]]
  standard_population <- aligned[[3L]]
  standard_total <- sum(standard_population)
  if (!is.finite(standard_total) || standard_total <= 0) {
    stop("standard_population must have a finite positive total",
         call. = FALSE)
  }

  weights <- standard_population / standard_total
  stratum_rates <- cases / person_time
  estimate <- sum(weights * stratum_rates)
  variance <- sum(weights^2 * cases / person_time^2)
  if (!is.finite(estimate) || !is.finite(variance)) {
    stop("The supplied aggregates produce a non-finite rate or variance",
         call. = FALSE)
  }
  std_error <- sqrt(variance)
  alpha <- 1 - level

  if (identical(ci_method, "normal")) {
    raw_interval <- estimate + c(lower = -1, upper = 1) *
      stats::qnorm(1 - alpha / 2) * std_error
    lower_truncated <- raw_interval[["lower"]] < 0
    interval <- c(lower = max(0, raw_interval[["lower"]]),
                  upper = raw_interval[["upper"]])
    correction <- list(
      name = "non-negative lower bound",
      applied = lower_truncated,
      value = if (lower_truncated) -raw_interval[["lower"]] * scale else 0)
  } else {
    one_event <- max(weights / person_time)
    if (estimate > 0 && variance > 0) {
      lower <- stats::qgamma(
        alpha / 2, shape = estimate^2 / variance,
        scale = variance / estimate)
    } else {
      lower <- 0
    }
    upper_estimate <- estimate + one_event
    upper_variance <- variance + one_event^2
    upper <- stats::qgamma(
      1 - alpha / 2,
      shape = upper_estimate^2 / upper_variance,
      scale = upper_variance / upper_estimate)
    interval <- c(lower = lower, upper = upper)
    correction <- list(
      name = "Fay-Feuer one-event upper-tail adjustment",
      applied = TRUE,
      value = one_event * scale)
  }

  scaled_estimate <- estimate * scale
  scaled_std_error <- std_error * scale
  scaled_interval <- interval * scale
  scaled_stratum_rates <- stratum_rates * scale
  if (any(!is.finite(c(scaled_estimate, scaled_std_error,
                       scaled_interval, scaled_stratum_rates,
                       correction$value)))) {
    stop("The requested reported scale produces non-finite results",
         call. = FALSE)
  }

  out <- list(
    call = match.call(),
    estimate = scaled_estimate,
    std_error = scaled_std_error,
    conf_int = scaled_interval,
    level = level,
    ci_method = ci_method,
    scale = scale,
    cases = cases,
    person_time = person_time,
    standard_population = standard_population,
    standard_weights = weights,
    stratum_rates = scaled_stratum_rates,
    correction = correction,
    assumptions = c(
      "Inputs are disclosure-authorized aggregates supplied by the caller; no DSI request is made.",
      "Stratum event counts are independent Poisson observations.",
      "Rates are constant within strata and the external standard population is appropriate for the target comparison."),
    input_provenance = "caller_supplied_aggregate",
    method_status = "promoted",
    privacy = list(dsi_calls = 0L, new_server_outputs = FALSE),
    quality = .dsvert_quality(status = "ok"))
  class(out) <- c("ds.vertDirectStandardization", "list")
  out
}

#' @export
print.ds.vertDirectStandardization <- function(x, ...) {
  cat("Directly standardized rate\n")
  cat(sprintf("  Estimate: %.6g per %.6g\n", x$estimate, x$scale))
  cat(sprintf("  %.1f%% %s CI: %.6g to %.6g\n",
              100 * x$level, x$ci_method,
              x$conf_int[["lower"]], x$conf_int[["upper"]]))
  invisible(x)
}

#' Indirect standardization from aggregate observed and expected events
#'
#' Computes a standardized mortality ratio (SMR) or standardized incidence
#' ratio (SIR) from caller-supplied aggregate observed and expected events.
#' The two labels share the same estimator; \code{measure} records the intended
#' scientific interpretation. The confidence interval is the exact central
#' Garwood interval for a Poisson count with fixed expected events.
#'
#' This function is entirely client-side. It makes no DSI or DataSHIELD calls
#' and exposes no new server-side values.
#'
#' @param observed Numeric vector of non-negative, whole-number aggregate
#'   observed event counts by stratum.
#' @param expected Numeric vector of non-negative expected events by stratum.
#' @param measure Either \code{"SMR"} or \code{"SIR"}.
#' @param scale Positive multiplier for the reported ratio. The conventional
#'   value \code{100} reports percent of expected events; use \code{1} for an
#'   unscaled ratio.
#' @param level Confidence level.
#'
#' @return An object of class \code{ds.vertIndirectStandardization} containing
#'   the observed-to-expected ratio, its scaled estimate, an exact Poisson
#'   confidence interval, provenance, and assumptions.
#'
#' @references Garwood, F. (1936). Fiducial limits for the Poisson
#'   distribution. \emph{Biometrika}, 28, 437--442.
#' @export
ds.vertIndirectStandardization <- function(
    observed, expected, measure = c("SMR", "SIR"), scale = 100,
    level = 0.95) {
  measure <- match.arg(measure)
  .ds_epi_validate_level_scale(level, scale)

  original_names <- list(names(observed), names(expected))
  observed <- .ds_epi_validate_vector(
    observed, "observed", whole_number = TRUE)
  expected <- .ds_epi_validate_vector(expected, "expected")
  names(observed) <- original_names[[1L]]
  names(expected) <- original_names[[2L]]
  aligned <- .ds_epi_align_strata(
    list(observed, expected), c("observed", "expected"))
  observed <- aligned[[1L]]
  expected <- aligned[[2L]]

  total_observed <- sum(observed)
  total_expected <- sum(expected)
  if (!is.finite(total_observed)) {
    stop("observed must have a finite total", call. = FALSE)
  }
  if (!is.finite(total_expected) || total_expected <= 0) {
    stop("expected must have a finite positive total", call. = FALSE)
  }
  alpha <- 1 - level
  lower <- if (total_observed == 0) {
    0
  } else {
    0.5 * stats::qchisq(alpha / 2, 2 * total_observed) / total_expected
  }
  upper <- 0.5 * stats::qchisq(
    1 - alpha / 2, 2 * (total_observed + 1)) / total_expected
  ratio <- total_observed / total_expected
  scaled_estimate <- ratio * scale
  scaled_interval <- c(lower = lower, upper = upper) * scale
  if (any(!is.finite(c(scaled_estimate, scaled_interval)))) {
    stop("The requested reported scale produces non-finite results",
         call. = FALSE)
  }

  out <- list(
    call = match.call(),
    measure = measure,
    ratio = ratio,
    estimate = scaled_estimate,
    conf_int = scaled_interval,
    level = level,
    ci_method = "exact_poisson_garwood",
    scale = scale,
    observed = observed,
    expected = expected,
    total_observed = total_observed,
    total_expected = total_expected,
    correction = list(name = "none", applied = FALSE, value = 0),
    assumptions = c(
      "Inputs are disclosure-authorized aggregates supplied by the caller; no DSI request is made.",
      "The observed total is Poisson and the expected total is fixed and externally valid."),
    input_provenance = "caller_supplied_aggregate",
    method_status = "promoted",
    privacy = list(dsi_calls = 0L, new_server_outputs = FALSE),
    quality = .dsvert_quality(status = "ok"))
  class(out) <- c("ds.vertIndirectStandardization", "list")
  out
}

#' @export
print.ds.vertIndirectStandardization <- function(x, ...) {
  cat(sprintf("Indirect standardization (%s)\n", x$measure))
  cat(sprintf("  Observed / expected: %.6g / %.6g\n",
              x$total_observed, x$total_expected))
  cat(sprintf("  %s: %.6g (%.1f%% exact Poisson CI %.6g to %.6g)\n",
              x$measure, x$estimate, 100 * x$level,
              x$conf_int[["lower"]], x$conf_int[["upper"]]))
  invisible(x)
}
