.dsvert_dp_logrank_counts <- function(x, cause) {
  time_count <- length(x$time_grid)
  outcome_count <- length(x$causes) + 1L
  exit_start <- 1L + if (isTRUE(x$delayed_entry)) time_count else 0L
  exit_end <- exit_start + time_count * outcome_count - 1L
  exit_lower <- matrix(
    x$histogram_mechanism_lower_95[exit_start:exit_end],
    nrow = time_count, ncol = outcome_count)
  exit_upper <- matrix(
    x$histogram_mechanism_upper_95[exit_start:exit_end],
    nrow = time_count, ncol = outcome_count)
  event_columns <- if (is.null(cause)) {
    seq_along(x$causes) + 1L
  } else {
    match(cause, x$causes) + 1L
  }
  list(
    event = rowSums(x$exit_histogram[, event_columns, drop = FALSE]),
    event_lower = rowSums(exit_lower[, event_columns, drop = FALSE]),
    event_upper = rowSums(exit_upper[, event_columns, drop = FALSE]),
    at_risk = x$curve$at_risk_dp,
    at_risk_lower = x$mechanism_bands$at_risk$lower,
    at_risk_upper = x$mechanism_bands$at_risk$upper)
}

.dsvert_dp_logrank_score_box <- function(comparison, reference) {
  fields <- c(
    "event", "event_lower", "event_upper", "at_risk",
    "at_risk_lower", "at_risk_upper")
  valid <- function(value) {
    is.numeric(value) && length(value) && !anyNA(value) &&
      all(is.finite(value)) && all(value >= 0)
  }
  if (!all(vapply(comparison[fields], valid, logical(1L))) ||
      !all(vapply(reference[fields], valid, logical(1L))) ||
      length(unique(c(
        vapply(comparison[fields], length, integer(1L)),
        vapply(reference[fields], length, integer(1L))))) != 1L ||
      any(comparison$event_lower > comparison$event_upper) ||
      any(reference$event_lower > reference$event_upper) ||
      any(comparison$at_risk_lower > comparison$at_risk_upper) ||
      any(reference$at_risk_lower > reference$at_risk_upper)) {
    stop("The DP log-rank count box is invalid", call. = FALSE)
  }

  d1 <- comparison$event
  d0 <- reference$event
  n1 <- comparison$at_risk
  n0 <- reference$at_risk
  n <- n1 + n0
  point_terms <- ifelse(n > 0, (d1 * n0 - d0 * n1) / n, 0)

  lower_denominator <- comparison$at_risk_upper +
    reference$at_risk_lower
  upper_denominator <- comparison$at_risk_lower +
    reference$at_risk_upper
  lower_terms <- ifelse(
    lower_denominator > 0,
    (comparison$event_lower * reference$at_risk_lower -
       reference$event_upper * comparison$at_risk_upper) /
      lower_denominator,
    -reference$event_upper)
  upper_terms <- ifelse(
    upper_denominator > 0,
    (comparison$event_upper * reference$at_risk_upper -
       reference$event_lower * comparison$at_risk_lower) /
      upper_denominator,
    comparison$event_upper)

  # This is only a non-negative plug-in diagnostic.  It is deliberately not
  # used to form a test, a standard error, or a sampling interval.
  variance_terms <- ifelse(
    n > 1,
    n1 * n0 * (d1 + d0) * pmax(0, n - d1 - d0) /
      (n^2 * (n - 1)),
    0)
  variance <- sum(variance_terms)
  score <- sum(point_terms)
  list(
    score = score, lower = sum(lower_terms), upper = sum(upper_terms),
    null_variance_plugin = variance,
    standardised_score_plugin = if (variance > 0) score / sqrt(variance) else NA_real_,
    contributing_time_points = sum(n > 0),
    zero_risk_time_points = sum(n <= 0),
    status = if (any(n > 0)) "ok" else "no_positive_combined_risk")
}

.dsvert_dp_logrank_attributes <- function(result, inputs, cause) {
  attr(result, "comparison_label") <- unname(inputs$labels[["comparison"]])
  attr(result, "reference_label") <- unname(inputs$labels[["reference"]])
  attr(result, "contrast_estimand") <- if (is.null(cause)) {
    "fixed-grid all-cause log-rank score from the declared release populations"
  } else {
    paste(
      "fixed-grid cause-specific log-rank score for", cause,
      "from the declared release populations")
  }
  attr(result, "uncertainty_scope") <- paste(
    "simultaneous DP mechanism noise only; sampling uncertainty and",
    "public-grid discretisation error excluded")
  attr(result, "joint_mechanism_confidence") <- inputs$joint_confidence
  attr(result, "joint_event") <- if (inputs$same_artifact) {
    "same_signed_survival_artifact"
  } else {
    "bonferroni_across_distinct_releases"
  }
  attr(result, "joint_confidence_method") <- if (inputs$same_artifact) {
    paste(
      "the two inputs identify the same signed survival artifact, so its",
      "original simultaneous mechanism event is retained")
  } else {
    paste(
      "Bonferroni lower bound max(0, confidence_comparison +",
      "confidence_reference - 1); no independence assumption")
  }
  attr(result, "statistical_inference") <- paste(
    "mechanism-only score region; null variance and standardised score are",
    "descriptive plug-in diagnostics, not a standard error, null-hypothesis",
    "test, sampling confidence interval, or p-value")
  attr(result, "population_comparison_contract") <- paste(
    "Input authentication establishes signed release provenance but cannot",
    "establish that declared comparison and reference populations are",
    "disjoint; that scientific cohort contract remains external")
  attr(result, "grid_error_scope") <- paste(
    "continuous-time and between-grid discretisation error is not included")
  attr(result, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  attr(result, "additional_server_calls") <- 0L
  attr(result, "source_release_provenance") <- list(
    comparison = .dsvert_dp_survival_source_provenance(inputs$comparison),
    reference = .dsvert_dp_survival_source_provenance(inputs$reference))
  result
}

#' Fixed-grid differentially private log-rank score
#'
#' Computes an all-cause or cause-specific fixed-grid log-rank score from two
#' compatible, already released `ds.vertDPSurvival` objects. It makes no
#' server call, draws no new noise and consumes no additional privacy. Its
#' conservative region uses the two releases' simultaneous DP histogram boxes;
#' distinct releases use a Bonferroni joint-confidence lower bound without an
#' independence assumption. The null-variance and standardised-score columns
#' are descriptive plug-in diagnostics only. No null-hypothesis test, p-value,
#' sampling confidence interval or proportional-hazards estimate is provided.
#'
#' The authenticated releases do not prove that the two declared populations
#' are disjoint. That cohort definition, and whether a cause-specific score is
#' clinically meaningful when competing events are present, remain scientific
#' preconditions outside this client-side post-processing.
#'
#' @param comparison A validated released `ds.vertDPSurvival` object for the
#'   comparison group.
#' @param reference A validated released `ds.vertDPSurvival` object for the
#'   reference group on the identical public time grid.
#' @param cause `NULL` for all released event causes, or one released cause
#'   label for a cause-specific score.
#' @param comparison_label Non-empty label for the comparison group.
#' @param reference_label Distinct non-empty label for the reference group.
#' @return A one-row `ds.vertDPLogRank` data frame containing the fixed-grid
#'   score, conservative joint DP-mechanism limits, descriptive plug-in
#'   diagnostics and source-release provenance.
#' @export
ds.vertDPLogRank <- function(
    comparison, reference, cause = NULL,
    comparison_label = "comparison", reference_label = "reference") {
  inputs <- .dsvert_dp_survival_contrast_inputs(
    comparison, reference, comparison_label, reference_label)
  if (!is.null(cause) &&
      (!is.character(cause) || length(cause) != 1L || is.na(cause) ||
       !cause %in% inputs$comparison$causes)) {
    stop("cause must be NULL or name one released event cause", call. = FALSE)
  }
  score <- .dsvert_dp_logrank_score_box(
    .dsvert_dp_logrank_counts(inputs$comparison, cause),
    .dsvert_dp_logrank_counts(inputs$reference, cause))
  result <- data.frame(
    cause = if (is.null(cause)) "all_causes" else cause,
    logrank_score = score$score,
    logrank_score_mechanism_lower = score$lower,
    logrank_score_mechanism_upper = score$upper,
    null_variance_plugin = score$null_variance_plugin,
    standardised_score_plugin = score$standardised_score_plugin,
    contributing_time_points = score$contributing_time_points,
    zero_risk_time_points = score$zero_risk_time_points,
    score_status = score$status,
    stringsAsFactors = FALSE)
  result <- .dsvert_dp_logrank_attributes(result, inputs, cause)
  class(result) <- c("ds.vertDPLogRank", class(result))
  result
}
