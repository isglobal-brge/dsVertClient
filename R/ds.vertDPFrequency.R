# One-way categorical epidemiology from the reusable joint-DP capsule.

.dsvert_dp_frequency_radius <- function(x, confidence) {
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      is.na(confidence) || !is.finite(confidence) ||
      confidence <= 0 || confidence >= 1) {
    stop("confidence must be one finite number in (0, 1)", call. = FALSE)
  }
  total_coordinates <- x$capsule_coordinate_count
  if (!.dsvert_dp_is_integer(total_coordinates, 1,
                             .DSVERT_DP_MAX_COORDINATES)) {
    stop("The frequency release has an invalid capsule dimension",
         call. = FALSE)
  }
  pseudo_release <- list(
    epsilon = x$epsilon, implementation_delta = x$implementation_delta,
    mechanism = x$mechanism, backend = x$implementation,
    mechanism_plan = x$mechanism_plan, plan_sha256 = x$plan_sha256,
    backend_selection = x$backend_selection,
    backend_assessment = x$backend_assessment,
    manifest_sha256 = x$manifest_sha256)
  pseudo_manifest <- list(workload = list(
    coordinate_count = total_coordinates,
    capsule_mechanism = x$capsule_mechanism,
    mechanism_selection = x$mechanism_selection,
    release_lattice = list(
      output_lattice_bits = x$output_lattice_bits,
      output_lattice_scale = x$output_lattice_scale,
      natural_l1_sensitivity = x$l1_sensitivity,
      integer_l1_sensitivity_steps =
        x$l1_sensitivity * x$output_lattice_scale,
      natural_l2_sensitivity = x$l2_sensitivity,
      integer_l2_sensitivity_steps =
        x$l2_sensitivity * x$output_lattice_scale)))
  .dsvert_dp_vector_accuracy_radius(
    pseudo_release, pseudo_manifest, coordinate_count = length(x$counts),
    confidence = confidence, maximum_error = x$coordinate_maximum)$radius
}

.dsvert_dp_frequency_regions <- function(counts, radius, capacity) {
  if (!is.numeric(counts) || !length(counts) || anyNA(counts) ||
      any(!is.finite(counts)) || any(counts < 0) ||
      !is.numeric(radius) || length(radius) != 1L || is.na(radius) ||
      !is.finite(radius) || radius < 0 ||
      !is.numeric(capacity) || length(capacity) != 1L || is.na(capacity) ||
      !is.finite(capacity) || capacity < 1) {
    stop("Invalid frequency-region inputs", call. = FALSE)
  }
  lower <- pmax(0, counts - radius)
  upper <- pmin(capacity, counts + radius)
  intervals <- t(vapply(seq_along(counts), function(index) {
    other <- setdiff(seq_along(counts), index)
    other_lower <- sum(lower[other])
    other_upper <- min(sum(upper[other]), capacity - lower[[index]])
    denominator_lower <- lower[[index]] + max(0, other_upper)
    denominator_upper <- upper[[index]] + other_lower
    c(
      lower = if (denominator_lower > 0) {
        lower[[index]] / denominator_lower
      } else 0,
      upper = if (denominator_upper > 0) {
        upper[[index]] / denominator_upper
      } else 1)
  }, numeric(2L)))
  rownames(intervals) <- names(counts)
  list(
    count = cbind(lower = lower, upper = upper),
    proportion = intervals,
    includes_zero_effective_count = sum(lower) == 0,
    has_positive_effective_count = sum(upper) > 0)
}

.dsvert_dp_frequency_contract <- function(x) {
  valid_number <- function(value, lower = 0, upper = Inf,
                           lower_open = FALSE) {
    .dsvert_dp_is_number(value, lower, upper, lower_open)
  }
  descriptor <- if (is.list(x)) x$coordinate_descriptor else NULL
  levels <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    descriptor$levels, "categorical levels"), error = function(error) NULL)
  expected_l1 <- if (is.list(x) && identical(
      x$adjacency, "replace_one_fixed_cohort")) 2 else 1
  expected_l2 <- sqrt(expected_l1)
  valid <- inherits(x, "ds.vertDPFrequency") && is.list(x) &&
    isTRUE(x$released) &&
    identical(x$backend, "exact_signed_Ring128_global_vector") &&
    identical(x$coordinate_family, "categorical_marginals") &&
    is.character(levels) && length(levels) > 0L && !anyNA(levels) &&
    !anyDuplicated(levels) &&
    is.numeric(x$counts) && length(x$counts) == length(levels) &&
    identical(names(x$counts), levels) && !anyNA(x$counts) &&
    all(is.finite(x$counts)) && all(x$counts >= 0) &&
    valid_number(x$coordinate_maximum, 1) &&
    all(x$counts <= x$coordinate_maximum) &&
    identical(x$levels, levels) &&
    identical(descriptor$repeated_record_policy,
              "consistent_level_else_exclude_v1") &&
    identical(descriptor$missingness_policy,
              "missing_or_out_of_domain_rows_are_ignored") &&
    x$adjacency %in% c(
      "add_remove_patient", "replace_one_fixed_cohort") &&
    .dsvert_dp_num_equal(x$artifact_l1_sensitivity, expected_l1) &&
    .dsvert_dp_num_equal(x$artifact_l2_sensitivity, expected_l2, 2048) &&
    valid_number(x$epsilon, 0, .DSVERT_DP_MAXIMUM_EPSILON, TRUE) &&
    valid_number(x$delta, 0, 1) && x$delta < 1 &&
    valid_number(x$l1_sensitivity, 0, Inf, TRUE) &&
    valid_number(x$l2_sensitivity, 0, Inf, TRUE) &&
    identical(x$sticky_noise, "immutable_capsule_durable_replay_v3") &&
    identical(x$sticky_replay, TRUE) &&
    identical(x$history_gate, TRUE) && identical(x$request_limit, FALSE) &&
    identical(x$operation_limit, TRUE) &&
    identical(x$source_values_exposed, FALSE) &&
    identical(x$intermediate_values_exposed, FALSE) &&
    identical(x$additional_privacy_cost, c(epsilon = 0, delta = 0))
  if (isTRUE(valid)) {
    radius <- tryCatch(.dsvert_dp_frequency_radius(x, 0.95),
                       error = function(error) NA_real_)
    valid <- is.finite(radius) &&
      .dsvert_dp_num_equal(radius, x$accuracy_simultaneous_95_abs, 2048)
  }
  if (!isTRUE(valid)) {
    stop("x must be a released, validated ds.vertDPFrequency object",
         call. = FALSE)
  }
  x
}

#' Differentially private one-way frequency distribution
#'
#' Extracts one signed categorical-marginal block from the single reusable
#' joint-DP capsule. Each admitted patient contributes to at most one public
#' level after the custodian-defined repeated-record collapse. The returned
#' count and proportion regions cover mechanism noise only; no exact marginal,
#' missingness count, or patient-level value is requested.
#'
#' @param data_name Signed protected dataset name.
#' @param variable Signed fixed-domain categorical variable.
#' @param server Optional expected owner server.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPFrequency` object.
#' @export
ds.vertDPFrequency <- function(data_name, variable, server = NULL,
                               datasources = NULL) {
  .dsvert_dp_frequency_impl(
    data_name, variable, server, datasources, DSI::datashield.aggregate)
}

.dsvert_dp_frequency_impl <- function(
    data_name, variable, server = NULL, datasources = NULL, .aggregate) {
  values <- list(data_name = data_name, variable = variable)
  if (any(!vapply(values, function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
          nzchar(value)
      }, logical(1L)))) {
    stop("data_name and variable must be non-empty strings", call. = FALSE)
  }
  datasources <- .dsvert_dp_datasources(datasources)
  owner <- .dsvert_dp_vector_server_filter(server, datasources)
  run <- .dsvert_dp_capsule_vector_run(
    datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run)
  block <- .dsvert_dp_capsule_single_block(
    context$layout, "categorical_marginals", dataset = data_name,
    owner_peer = owner,
    predicate = function(candidate) {
      identical(candidate$descriptor$column, variable)
    },
    description = paste0("signed categorical-marginal block for '",
                         variable, "'"))
  descriptor <- block$descriptor
  levels <- .dsvert_dp_capsule_manifest_strings(
    descriptor$levels, "categorical levels")
  if (!length(levels) || anyDuplicated(levels) ||
      !identical(descriptor$repeated_record_policy,
                 "consistent_level_else_exclude_v1") ||
      !identical(descriptor$missingness_policy,
                 "missing_or_out_of_domain_rows_are_ignored")) {
    stop("The signed categorical-marginal descriptor is invalid",
         call. = FALSE)
  }
  maximum <- .dsvert_dp_capsule_manifest_numbers(
    descriptor$statistic_maximum, "categorical-marginal maximum")
  if (length(maximum) != 1L || maximum[[1L]] < 1 ||
      maximum[[1L]] != floor(maximum[[1L]]) ||
      maximum[[1L]] > 2^53 - 1) {
    stop("The signed categorical-marginal bound is invalid", call. = FALSE)
  }
  counts <- .dsvert_dp_capsule_vector_values(context$release, block)
  if (length(counts) != length(levels) || any(counts < 0) ||
      any(counts > maximum[[1L]])) {
    stop("The released categorical marginal violates its signed bounds",
         call. = FALSE)
  }
  names(counts) <- levels
  simultaneous <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest, coordinate_count = length(counts),
    confidence = 0.95, maximum_error = maximum[[1L]])
  regions <- .dsvert_dp_frequency_regions(
    counts, simultaneous$radius, maximum[[1L]])
  total <- sum(counts)
  proportions <- if (total > 0) counts / total else
    stats::setNames(rep(NA_real_, length(counts)), levels)
  artifact_l1 <- if (identical(
      context$adjacency, "replace_one_fixed_cohort")) 2 else 1
  result <- c(.dsvert_dp_vector_public_metadata(context), list(
    status = if (total > 0) "ok" else "dp_effective_count_zero",
    server = block$owner_peer, data_name = data_name, variable = variable,
    coordinate_family = "categorical_marginals",
    coordinate_descriptor = descriptor, levels = levels,
    counts = unname(counts), effective_count_dp = total,
    proportions = proportions,
    coordinate_maximum = maximum[[1L]],
    unit_aggregation_policy = "consistent_level_else_exclude_v1",
    artifact_l1_sensitivity = artifact_l1,
    artifact_l2_sensitivity = sqrt(artifact_l1),
    accuracy_simultaneous_95_abs = simultaneous$radius,
    accuracy_simultaneous_confidence = simultaneous$confidence,
    accuracy_simultaneous_method = simultaneous$method,
    accuracy_implementation_tv_upper_bound =
      simultaneous$implementation_tv_upper_bound,
    mechanism_regions = regions,
    uncertainty_scope = paste(
      "Simultaneous regions cover DP mechanism noise only; clipping,",
      "category collapse and population sampling uncertainty are excluded"),
    inferential_scope = paste(
      "Finite-dataset DP frequency distribution; use",
      "ds.vertDPFrequencyInference() for conservative iid multinomial",
      "sampling regions"),
    additional_server_calls_after_capsule = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0)))
  names(result$counts) <- levels
  class(result) <- c("ds.vertDPFrequency", "list")
  .dsvert_dp_frequency_contract(result)
}

#' Conservative sampling regions for a DP frequency distribution
#'
#' Combines one simultaneous DP count-box event with exact Clopper--Pearson
#' intervals for every public level. A Bonferroni union bound supplies the
#' advertised joint coverage. This is pure post-processing and makes no DSI
#' call.
#'
#' @param x A validated `ds.vertDPFrequency` object.
#' @param level Requested joint coverage.
#' @param dp_fraction Fraction of the total error probability assigned to the
#'   simultaneous DP mechanism event; the remainder is split over categories.
#' @return A `ds.vertDPFrequencyInference` object.
#' @export
ds.vertDPFrequencyInference <- function(x, level = 0.95,
                                        dp_fraction = 0.5) {
  x <- .dsvert_dp_frequency_contract(x)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1 ||
      !is.numeric(dp_fraction) || length(dp_fraction) != 1L ||
      is.na(dp_fraction) || !is.finite(dp_fraction) ||
      dp_fraction <= 0 || dp_fraction >= 1) {
    stop("level and dp_fraction must be finite numbers in (0, 1)",
         call. = FALSE)
  }
  alpha <- 1 - level
  dp_alpha <- alpha * dp_fraction
  sampling_alpha <- alpha - dp_alpha
  per_level_alpha <- sampling_alpha / length(x$counts)
  radius <- .dsvert_dp_frequency_radius(x, 1 - dp_alpha)
  regions <- .dsvert_dp_frequency_regions(
    x$counts, radius, x$coordinate_maximum)
  integer_lower <- vapply(seq_along(x$counts), function(index) {
    .dsvert_dp_integer_count_range(
      regions$count[index, "lower"],
      regions$count[index, "upper"])[["lower"]]
  }, numeric(1L))
  integer_upper <- vapply(seq_along(x$counts), function(index) {
    .dsvert_dp_integer_count_range(
      regions$count[index, "lower"],
      regions$count[index, "upper"])[["upper"]]
  }, numeric(1L))
  intervals <- t(vapply(seq_along(x$counts), function(index) {
    other <- setdiff(seq_along(x$counts), index)
    failure_lower <- sum(integer_lower[other])
    failure_upper <- min(
      sum(integer_upper[other]),
      x$coordinate_maximum - integer_lower[[index]])
    .dsvert_dp_cp_union_over_box(
      integer_lower[[index]], integer_upper[[index]],
      failure_lower, max(failure_lower, failure_upper),
      per_level_alpha)
  }, numeric(2L)))
  rownames(intervals) <- x$levels
  result <- list(
    status = "ok", levels = x$levels, counts_dp = x$counts,
    proportions_dp = x$proportions,
    intervals = intervals, level = level,
    coverage_lower_bound = level,
    coverage_method = paste(
      "simultaneous DP mechanism count box plus Bonferroni union of",
      "Clopper-Pearson binomial intervals"),
    dp_event_confidence = 1 - dp_alpha,
    sampling_familywise_confidence = 1 - sampling_alpha,
    base_sampling_interval_confidence = 1 - per_level_alpha,
    mechanism_radius = radius,
    sampling_model = paste(
      "iid multinomial privacy units conditional on a valid consistent",
      "public-domain category"),
    selection_warning = paste(
      "Missing, out-of-domain and conflicting repeated records are excluded;",
      "a population interpretation requires that this selection be",
      "scientifically ignorable"),
    p_values = NULL, hypothesis_tests = NULL,
    additional_server_calls = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    capsule_id = x$capsule_id, manifest_sha256 = x$manifest_sha256)
  class(result) <- c("ds.vertDPFrequencyInference", "list")
  result
}

#' @export
print.ds.vertDPFrequency <- function(x, ...) {
  print(data.frame(
    level = x$levels, count_dp = unname(x$counts),
    proportion_dp = unname(x$proportions), row.names = NULL), ...)
  cat("simultaneous 95% mechanism radius:",
      x$accuracy_simultaneous_95_abs, "\n")
  invisible(x)
}

#' @export
print.ds.vertDPFrequencyInference <- function(x, ...) {
  print(x$intervals, ...)
  cat("Conservative joint coverage >= ", format(100 * x$level), "%\n",
      sep = "")
  invisible(x)
}
