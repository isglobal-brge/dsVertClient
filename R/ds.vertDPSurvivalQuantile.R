.dsvert_dp_survival_quantile_endpoint <- function(
    time_grid, survival, probability) {
  valid_curve <- is.numeric(time_grid) && is.numeric(survival) &&
    length(time_grid) > 0L && length(survival) == length(time_grid) &&
    !anyNA(time_grid) && !anyNA(survival) &&
    all(is.finite(time_grid)) && all(is.finite(survival)) &&
    all(diff(time_grid) > 0) && all(survival >= 0 & survival <= 1) &&
    all(diff(survival) <= 64 * .Machine$double.eps)
  valid_probability <- is.numeric(probability) &&
    length(probability) == 1L && !is.na(probability) &&
    is.finite(probability) && probability > 0 && probability < 1
  if (!isTRUE(valid_curve) || !isTRUE(valid_probability)) {
    stop("Invalid fixed-grid survival quantile inputs", call. = FALSE)
  }
  crossing <- which(survival <= 1 - probability)
  if (length(crossing)) time_grid[[crossing[[1L]]]] else Inf
}

.dsvert_dp_survival_quantile_region <- function(
    time_grid, survival_lower, survival_upper, probability) {
  valid <- is.numeric(survival_lower) && is.numeric(survival_upper) &&
    length(survival_lower) == length(time_grid) &&
    length(survival_upper) == length(time_grid) &&
    !anyNA(survival_lower) && !anyNA(survival_upper) &&
    all(is.finite(survival_lower)) && all(is.finite(survival_upper)) &&
    all(survival_lower <= survival_upper)
  if (!isTRUE(valid)) {
    stop("Invalid fixed-grid survival quantile band", call. = FALSE)
  }
  c(
    lower = .dsvert_dp_survival_quantile_endpoint(
      time_grid, survival_lower, probability),
    upper = .dsvert_dp_survival_quantile_endpoint(
      time_grid, survival_upper, probability))
}

#' Fixed-grid survival quantiles from one DP survival release
#'
#' Inverts the Kaplan--Meier curve and its simultaneous DP-mechanism band from
#' one validated `ds.vertDPSurvival` object. The finite-snapshot quantile is the
#' first public interval endpoint at which the product-limit survival is no
#' greater than `1 - probability`. If the target is not reached by the public
#' horizon, the point is `NA` and the mechanism-region endpoint is represented
#' by `Inf` with an explicit status flag. This is zero-call post-processing and
#' consumes no additional privacy. The limits exclude sampling uncertainty and
#' public-grid discretisation error.
#'
#' @param x A validated released `ds.vertDPSurvival` object.
#' @param probabilities One or more finite event-distribution probabilities
#'   strictly between zero and one.
#' @return A `ds.vertDPSurvivalQuantile` data frame containing fixed-grid
#'   quantiles, inverted simultaneous DP-mechanism limits, estimability flags,
#'   and complete source-release provenance.
#' @export
ds.vertDPSurvivalQuantile <- function(x, probabilities = 0.5) {
  x <- .dsvert_dp_survival_object(x)
  if (!is.numeric(probabilities) || !length(probabilities) ||
      anyNA(probabilities) || any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1)) {
    stop("probabilities must be finite numbers strictly between zero and one",
         call. = FALSE)
  }
  probabilities <- as.numeric(probabilities)
  point <- vapply(probabilities, function(probability) {
    .dsvert_dp_survival_quantile_endpoint(
      x$time_grid, x$curve$kaplan_meier, probability)
  }, numeric(1L))
  regions <- t(vapply(probabilities, function(probability) {
    .dsvert_dp_survival_quantile_region(
      x$time_grid,
      x$curve$kaplan_meier_mechanism_lower_95,
      x$curve$kaplan_meier_mechanism_upper_95,
      probability)
  }, numeric(2L)))
  finite_lower <- is.finite(regions[, "lower"])
  finite_upper <- is.finite(regions[, "upper"])
  result <- data.frame(
    probability = probabilities,
    survival_probability = 1 - probabilities,
    quantile = ifelse(is.finite(point), point, NA_real_),
    point_status = ifelse(
      is.finite(point), "reached_on_public_grid", "not_reached_by_grid_end"),
    quantile_mechanism_lower_95 = regions[, "lower"],
    quantile_mechanism_upper_95 = regions[, "upper"],
    mechanism_region_has_finite_values = finite_lower,
    mechanism_region_includes_beyond_grid = !finite_upper,
    status = ifelse(
      !finite_lower,
      "not_reached_by_any_curve_in_mechanism_band",
      ifelse(!finite_upper, "mechanism_region_extends_beyond_grid", "ok")),
    stringsAsFactors = FALSE)
  attr(result, "uncertainty_scope") <- x$uncertainty_scope
  attr(result, "mechanism_band_scope") <- x$mechanism_band_scope
  attr(result, "mechanism_band_tightness") <- x$mechanism_band_tightness
  attr(result, "mechanism_band_confidence") <-
    x$mechanism_band_confidence
  attr(result, "mechanism_band_method") <- x$mechanism_band_method
  attr(result, "grid_error_scope") <- x$grid_error_scope
  attr(result, "statistical_inference") <- x$statistical_inference
  attr(result, "discretisation_error") <- x$discretisation_error
  attr(result, "grid_quantile_definition") <- paste(
    "first public endpoint at which Kaplan-Meier survival is no greater",
    "than one minus the requested event-distribution probability")
  attr(result, "beyond_grid_representation") <- paste(
    "an infinite mechanism endpoint means the target quantile can lie",
    "beyond the signed public horizon; the point estimate is NA when the",
    "released product-limit curve itself does not cross")
  attr(result, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  attr(result, "additional_server_calls") <- 0L
  provenance_fields <- c(
    "analysis_id", "analysis_version", "server", "capsule_id",
    "final_vector_root", "coordinate_order_sha256", "privacy_epoch",
    "noise_key_id", "mechanism", "implementation", "sampler", "epsilon",
    "delta", "implementation_delta", "adjacency",
    "composition_partitions", "time_grid", "time_lower_bound",
    "time_upper_bound", "security_claim")
  attr(result, "source_release_provenance") <- c(
    list(source_class = "ds.vertDPSurvival"), x[provenance_fields])
  class(result) <- c("ds.vertDPSurvivalQuantile", class(result))
  result
}

#' Median survival from one DP survival release
#'
#' Exact convenience view of `ds.vertDPSurvivalQuantile(x, 0.5)`. It makes no
#' server call, draws no new noise, and retains the same fixed-grid and
#' uncertainty qualifications.
#'
#' @param x A validated released `ds.vertDPSurvival` object.
#' @return A one-row `ds.vertDPMedianSurvival` data frame retaining the full
#'   survival-quantile provenance and mechanism-band contract.
#' @export
ds.vertDPMedianSurvival <- function(x) {
  result <- ds.vertDPSurvivalQuantile(x, 0.5)
  class(result) <- c("ds.vertDPMedianSurvival", class(result))
  result
}
