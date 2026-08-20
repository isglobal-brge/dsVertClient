.dsvert_dp_survival_contrast_labels <- function(comparison_label,
                                                 reference_label) {
  valid <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value)
  }
  if (!valid(comparison_label) || !valid(reference_label) ||
      identical(comparison_label, reference_label)) {
    stop("comparison_label and reference_label must be distinct non-empty strings",
         call. = FALSE)
  }
  c(comparison = comparison_label, reference = reference_label)
}

.dsvert_dp_survival_same_signed_artifact <- function(comparison, reference) {
  hex256 <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      grepl("^[0-9a-f]{64}$", value)
  }
  synopsis <- .dsvert_dp_survival_is_synopsis(comparison) &&
    .dsvert_dp_survival_is_synopsis(reference)
  identity_fields <- if (isTRUE(synopsis)) {
    c(
      "analysis_id", "analysis_version", "server", "artifact_key",
      "execution_id", "contract_sha256", "attempt_sha256",
      "source_contract_sha256", "result_set_sha256", "final_vector_root",
      "coordinate_order_sha256")
  } else {
    c(
      "analysis_id", "analysis_version", "server", "capsule_id",
      "final_vector_root", "coordinate_order_sha256", "privacy_epoch")
  }
  required_hashes <- if (isTRUE(synopsis)) {
    c(
      "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
      "source_contract_sha256", "result_set_sha256", "final_vector_root",
      "coordinate_order_sha256")
  } else {
    c("capsule_id", "final_vector_root", "coordinate_order_sha256")
  }
  # Provenance fields are mutable after materialisation in R. Retain the
  # single-event coverage claim only for identical validated release objects.
  all(vapply(required_hashes, function(field) {
    hex256(comparison[[field]]) && hex256(reference[[field]])
  }, logical(1L))) &&
    all(vapply(identity_fields, function(field) {
      identical(comparison[[field]], reference[[field]])
    }, logical(1L))) &&
    identical(comparison, reference)
}

.dsvert_dp_survival_contrast_inputs <- function(
    comparison, reference, comparison_label, reference_label) {
  comparison <- .dsvert_dp_survival_object(comparison)
  reference <- .dsvert_dp_survival_object(reference)
  labels <- .dsvert_dp_survival_contrast_labels(
    comparison_label, reference_label)
  compatible <- identical(comparison$time_grid, reference$time_grid) &&
    identical(comparison$time_lower_bound, reference$time_lower_bound) &&
    identical(comparison$time_upper_bound, reference$time_upper_bound) &&
    identical(comparison$interval_semantics, reference$interval_semantics) &&
    identical(comparison$unit_collapse, reference$unit_collapse) &&
    identical(comparison$censor_level, reference$censor_level) &&
    identical(comparison$causes, reference$causes) &&
    identical(comparison$delayed_entry, reference$delayed_entry)
  if (!isTRUE(compatible)) {
    stop(
      paste(
        "comparison and reference must use the same signed public time grid,",
        "event/censor contract, delayed-entry semantics, and survival estimand"),
      call. = FALSE)
  }

  same_artifact <- .dsvert_dp_survival_same_signed_artifact(
    comparison, reference)
  comparison_confidence <- comparison$mechanism_band_confidence
  reference_confidence <- reference$mechanism_band_confidence
  joint_confidence <- if (same_artifact) {
    min(comparison_confidence, reference_confidence)
  } else {
    max(0, comparison_confidence + reference_confidence - 1)
  }
  list(
    comparison = comparison, reference = reference, labels = labels,
    same_artifact = same_artifact, joint_confidence = joint_confidence)
}

.dsvert_dp_nonnegative_interval_ratio <- function(
    numerator, denominator, numerator_lower, numerator_upper,
    denominator_lower, denominator_upper) {
  point <- ifelse(denominator > 0, numerator / denominator, NA_real_)
  lower <- ifelse(
    denominator_upper > 0, numerator_lower / denominator_upper, NA_real_)
  upper <- ifelse(
    denominator_lower > 0,
    numerator_upper / denominator_lower,
    ifelse(numerator_upper > 0, Inf, NA_real_))
  status <- ifelse(
    denominator <= 0,
    "reference_point_zero_not_estimable",
    ifelse(
      denominator_lower <= 0,
      "point_estimable_mechanism_upper_unbounded_at_zero",
      "ok"))
  list(
    point = point, lower = lower, upper = upper, status = status,
    denominator_interval_includes_zero = denominator_lower <= 0)
}

.dsvert_dp_survival_release_provenance <- function(x) {
  .dsvert_dp_survival_source_provenance(x)
}

.dsvert_dp_survival_contrast_attributes <- function(
    result, inputs, estimand) {
  attr(result, "comparison_label") <- unname(inputs$labels[["comparison"]])
  attr(result, "reference_label") <- unname(inputs$labels[["reference"]])
  attr(result, "contrast_estimand") <- estimand
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
    "mechanism-only contrast region; no sampling confidence interval,",
    "standard error, null-hypothesis test or p-value")
  attr(result, "grid_error_scope") <- paste(
    "continuous-time and between-grid discretisation error is not included")
  attr(result, "ratio_denominator_contract") <- paste(
    "the point ratio is unavailable when the reference point is zero;",
    "the mechanism upper limit is infinite when its reference interval",
    "includes zero and a positive numerator remains possible")
  attr(result, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  attr(result, "additional_server_calls") <- 0L
  attr(result, "source_release_provenance") <- list(
    comparison = .dsvert_dp_survival_release_provenance(
      inputs$comparison),
    reference = .dsvert_dp_survival_release_provenance(inputs$reference))
  result
}

#' Contrast two differentially private survival releases
#'
#' Computes fixed-grid survival differences and ratios from two already
#' released, validated `ds.vertDPSurvival` objects. It performs no server call
#' and consumes no additional privacy. Difference and ratio limits use
#' conservative interval arithmetic on the inputs' simultaneous
#' DP-mechanism bands. If the inputs are not the exact same signed survival
#' artifact, joint coverage is bounded by Bonferroni without assuming
#' independent DP noise. These are mechanism-only regions, not sampling
#' confidence intervals or hypothesis tests.
#'
#' @param comparison A validated released `ds.vertDPSurvival` object for the
#'   numerator/comparison group.
#' @param reference A validated released `ds.vertDPSurvival` object for the
#'   denominator/reference group on the identical public time grid.
#' @param comparison_label Non-empty label for the comparison group.
#' @param reference_label Distinct non-empty label for the reference group.
#' @return A `ds.vertDPSurvivalContrast` data frame with point contrasts,
#'   conservative joint DP-mechanism limits, denominator status and complete
#'   source-release provenance.
#' @export
ds.vertDPSurvivalContrast <- function(
    comparison, reference, comparison_label = "comparison",
    reference_label = "reference") {
  inputs <- .dsvert_dp_survival_contrast_inputs(
    comparison, reference, comparison_label, reference_label)
  comparison <- inputs$comparison
  reference <- inputs$reference
  comparison_curve <- comparison$curve
  reference_curve <- reference$curve
  ratio <- .dsvert_dp_nonnegative_interval_ratio(
    comparison_curve$kaplan_meier,
    reference_curve$kaplan_meier,
    comparison_curve$kaplan_meier_mechanism_lower_95,
    comparison_curve$kaplan_meier_mechanism_upper_95,
    reference_curve$kaplan_meier_mechanism_lower_95,
    reference_curve$kaplan_meier_mechanism_upper_95)

  result <- data.frame(
    time = comparison$time_grid,
    comparison_survival = comparison_curve$kaplan_meier,
    reference_survival = reference_curve$kaplan_meier,
    comparison_survival_mechanism_lower =
      comparison_curve$kaplan_meier_mechanism_lower_95,
    comparison_survival_mechanism_upper =
      comparison_curve$kaplan_meier_mechanism_upper_95,
    reference_survival_mechanism_lower =
      reference_curve$kaplan_meier_mechanism_lower_95,
    reference_survival_mechanism_upper =
      reference_curve$kaplan_meier_mechanism_upper_95,
    survival_difference =
      comparison_curve$kaplan_meier - reference_curve$kaplan_meier,
    survival_difference_mechanism_lower =
      comparison_curve$kaplan_meier_mechanism_lower_95 -
        reference_curve$kaplan_meier_mechanism_upper_95,
    survival_difference_mechanism_upper =
      comparison_curve$kaplan_meier_mechanism_upper_95 -
        reference_curve$kaplan_meier_mechanism_lower_95,
    survival_ratio = ratio$point,
    survival_ratio_mechanism_lower = ratio$lower,
    survival_ratio_mechanism_upper = ratio$upper,
    survival_ratio_status = ratio$status,
    reference_mechanism_interval_includes_zero =
      ratio$denominator_interval_includes_zero,
    stringsAsFactors = FALSE)
  result <- .dsvert_dp_survival_contrast_attributes(
    result, inputs,
    paste(
      "comparison minus/reference-divided fixed-grid Kaplan-Meier",
      "survival from the declared release populations"))
  class(result) <- c("ds.vertDPSurvivalContrast", class(result))
  result
}

#' Contrast restricted mean survival time from two DP releases
#'
#' Computes comparison-minus-reference RMST and comparison/reference RMST from
#' two compatible, already released DP survival objects. It reuses
#' `ds.vertDPRMST()` and conservative interval arithmetic, makes no server
#' call, draws no noise and spends no additional privacy. The joint
#' mechanism confidence is retained only for the exact same signed artifact;
#' otherwise it uses a Bonferroni lower bound. No sampling inference is
#' implied.
#'
#' @inheritParams ds.vertDPSurvivalContrast
#' @param tau One or more common finite restriction times within the signed
#'   public time interval. `NULL` uses the common upper bound.
#' @return A `ds.vertDPRMSTContrast` data frame with RMST differences, ratios,
#'   conservative joint mechanism limits and source-release provenance.
#' @export
ds.vertDPRMSTContrast <- function(
    comparison, reference, tau = NULL,
    comparison_label = "comparison", reference_label = "reference") {
  inputs <- .dsvert_dp_survival_contrast_inputs(
    comparison, reference, comparison_label, reference_label)
  comparison_rmst <- ds.vertDPRMST(inputs$comparison, tau)
  reference_rmst <- ds.vertDPRMST(inputs$reference, tau)
  ratio <- .dsvert_dp_nonnegative_interval_ratio(
    comparison_rmst$rmst, reference_rmst$rmst,
    comparison_rmst$rmst_mechanism_lower_95,
    comparison_rmst$rmst_mechanism_upper_95,
    reference_rmst$rmst_mechanism_lower_95,
    reference_rmst$rmst_mechanism_upper_95)

  result <- data.frame(
    time_lower_bound = comparison_rmst$time_lower_bound,
    tau = comparison_rmst$tau,
    comparison_rmst = comparison_rmst$rmst,
    reference_rmst = reference_rmst$rmst,
    comparison_rmst_mechanism_lower =
      comparison_rmst$rmst_mechanism_lower_95,
    comparison_rmst_mechanism_upper =
      comparison_rmst$rmst_mechanism_upper_95,
    reference_rmst_mechanism_lower =
      reference_rmst$rmst_mechanism_lower_95,
    reference_rmst_mechanism_upper =
      reference_rmst$rmst_mechanism_upper_95,
    rmst_difference = comparison_rmst$rmst - reference_rmst$rmst,
    rmst_difference_mechanism_lower =
      comparison_rmst$rmst_mechanism_lower_95 -
        reference_rmst$rmst_mechanism_upper_95,
    rmst_difference_mechanism_upper =
      comparison_rmst$rmst_mechanism_upper_95 -
        reference_rmst$rmst_mechanism_lower_95,
    rmst_ratio = ratio$point,
    rmst_ratio_mechanism_lower = ratio$lower,
    rmst_ratio_mechanism_upper = ratio$upper,
    rmst_ratio_status = ratio$status,
    reference_mechanism_interval_includes_zero =
      ratio$denominator_interval_includes_zero,
    stringsAsFactors = FALSE)
  result <- .dsvert_dp_survival_contrast_attributes(
    result, inputs,
    paste(
      "comparison minus/reference-divided fixed-grid restricted mean",
      "survival time through common tau"))
  attr(result, "integration_rule") <- attr(
    comparison_rmst, "integration_rule")
  class(result) <- c("ds.vertDPRMSTContrast", class(result))
  result
}
