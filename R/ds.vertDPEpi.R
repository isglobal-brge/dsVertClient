.dsvert_dp_table_vector_profile <- function(x) {
  capsule_mechanism <- x$capsule_mechanism
  if (is.null(capsule_mechanism)) {
    capsule_mechanism <- if (identical(
        x$mechanism, .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM)) {
      "discrete-laplace"
    } else if (identical(
        x$mechanism, .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM
    } else {
      return(NULL)
    }
  }
  tryCatch(
    .dsvert_vector_profile(
      capsule_mechanism, x$mechanism_selection,
      backend = if (.dsvert_dp_is_string(x$implementation)) {
        x$implementation
      } else NULL),
    error = function(error) NULL)
}

.dsvert_dp_table_vector_accuracy_method_is_valid <- function(x, profile) {
  if (!is.character(x$accuracy_simultaneous_method) ||
      length(x$accuracy_simultaneous_method) != 1L ||
      is.na(x$accuracy_simultaneous_method)) return(FALSE)
  expected <- if (isTRUE(profile$gaussian)) {
    paste(
      "signed fixed-work dyadic discrete-Gaussian plan v2 simultaneous",
      "95% bound; tail and CDF TV transfers already charged;",
      "fixed-clamp range applied")
  } else if (isTRUE(profile$exact_gc)) {
    paste(
      "exact ideal one-draw two-sided-geometric tail with union bound;",
      "signed vector sampler TV deducted once; clamp inside exact GC applied")
  } else {
    paste(
      "exact ideal two-sided-geometric convolution tail with union bound;",
      "two-peer finite-sampler TV deducted; fixed-clamp range applied")
  }
  identical(x$accuracy_simultaneous_method, expected)
}

.dsvert_dp_table_synopsis_contract <- function(x) {
  finite_scalar <- function(value, lower = -Inf, upper = Inf,
                            lower_open = FALSE) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) &&
      (if (lower_open) value > lower else value >= lower) &&
      value <= upper
  }
  table <- x$table
  profile <- .dsvert_dp_table_vector_profile(x)
  valid_table <- is.matrix(table) && is.numeric(table) &&
    length(table) > 0L && !anyNA(table) && all(is.finite(table)) &&
    all(table >= 0) && finite_scalar(x$coordinate_maximum, 1) &&
    all(table <= x$coordinate_maximum) &&
    identical(as.integer(nrow(table)), x$nrow) &&
    identical(as.integer(ncol(table)), x$ncol) &&
    identical(unname(rownames(table)), unname(x$row_levels)) &&
    identical(unname(colnames(table)), unname(x$col_levels)) &&
    identical(unname(as.numeric(table)), unname(x$counts))
  roots <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  valid_roots <- all(vapply(roots, function(field) {
    .dsvert_vector_hex(x[[field]]) &&
      identical(x$release_provenance[[field]], x[[field]])
  }, logical(1L)))
  privacy <- x$privacy
  valid_privacy <- is.list(privacy) && identical(names(privacy), c(
    "version", "per_artifact_epsilon", "per_artifact_delta",
    "sticky_noise", "unlimited_replay", "replay_is_postprocessing",
    "public_openings", "distinct_artifacts_compose",
    "finite_global_composition_claim")) &&
    identical(privacy$version, "dsvert-per-synopsis-dp-v1") &&
    .dsvert_dp_num_equal(privacy$per_artifact_epsilon, x$epsilon, 128) &&
    .dsvert_dp_num_equal(privacy$per_artifact_delta, x$delta, 128) &&
    identical(privacy$sticky_noise, TRUE) &&
    identical(privacy$unlimited_replay, TRUE) &&
    identical(privacy$replay_is_postprocessing, TRUE) &&
    identical(privacy$public_openings, 1L) &&
    identical(privacy$distinct_artifacts_compose, TRUE) &&
    identical(privacy$finite_global_composition_claim, FALSE)
  recursive_names <- function(value) {
    if (!is.list(value)) return(character())
    c(names(value), unlist(lapply(value, recursive_names), use.names = FALSE))
  }
  forbidden <- c(
    "capsule_id", "capsule_coordinate_count", "privacy_epoch",
    "privacy_epochs", "noise_key_id", "noise_key_ids", "history_gate",
    "request_limit", "operation_limit", "lifetime_budget",
    "lifetime_composition", "privacy_accountant", "release_instance",
    "release_instance_id", "allocation_certificate", "ledger",
    "reservation", "rate_limit", "catalog_limit", "quota")
  no_legacy <- !any(recursive_names(x) %in% forbidden)
  radius <- tryCatch(
    .dsvert_dp_vector_table_radius(x, 0.95),
    error = function(error) NA_real_)
  selected_sensitivity <- if (is.list(profile) && isTRUE(profile$gaussian)) {
    x$l2_sensitivity
  } else x$l1_sensitivity
  descriptor <- x$coordinate_descriptor
  cross_owner <- isTRUE(x$cross_owner)
  cross_owners <- if (cross_owner && is.list(descriptor)) sort(unique(c(
    descriptor$left$owner_peer, descriptor$right$owner_peer
  )), method = "radix") else character()
  cross_datasets <- if (cross_owner && is.list(descriptor)) sort(unique(c(
    descriptor$left$dataset, descriptor$right$dataset
  )), method = "radix") else character()
  valid_unit_policy <- if (cross_owner) {
    identical(
      descriptor$version,
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION) &&
      identical(
        x$unit_aggregation_policy,
        .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_UNIT_POLICY) &&
      length(cross_owners) == 2L &&
      identical(sort(unique(x$servers), method = "radix"), cross_owners) &&
      identical(sort(unique(x$datasets), method = "radix"), cross_datasets)
  } else {
    identical(x$cross_owner, FALSE) &&
      identical(x$unit_aggregation_policy,
                "consistent_joint_cell_else_exclude_v1")
  }
  valid <- isTRUE(x$released) && valid_table && valid_roots &&
    valid_privacy && no_legacy && is.list(profile) &&
    identical(x$implementation, profile$backend) &&
    identical(x$sampler, profile$sampler) &&
    identical(x$mechanism, profile$release_mechanism) &&
    identical(
      x$randomness,
      paste("independent pinned-peer HKDF-SHA256/ChaCha20 streams;",
            "no analyst-controlled seed")) &&
    identical(x$sticky_noise,
              "one_immutable_canonical_synopsis_artifact") &&
    identical(x$sticky_replay, TRUE) &&
    identical(x$unlimited_replay, TRUE) &&
    identical(x$sensitivity_scope,
              "complete_signed_canonical_synopsis_vector") &&
    finite_scalar(x$epsilon, 0, .DSVERT_DP_MAXIMUM_EPSILON, TRUE) &&
    finite_scalar(x$delta, 0, 1) && x$delta < 1 &&
    finite_scalar(x$sensitivity, 0, Inf, TRUE) &&
    identical(x$sensitivity_norm,
              if (isTRUE(profile$gaussian)) "l2" else "l1") &&
    finite_scalar(x$l1_sensitivity, 0, Inf, TRUE) &&
    finite_scalar(x$l2_sensitivity, 0, Inf, TRUE) &&
    .dsvert_dp_num_equal(selected_sensitivity, x$sensitivity, 2048) &&
    x$adjacency %in% c(
      "add_remove_patient", "replace_one_fixed_cohort") &&
    .dsvert_dp_num_equal(
      x$artifact_l1_sensitivity,
      if (identical(x$adjacency, "add_remove_patient")) 1 else 2) &&
    .dsvert_dp_num_equal(
      x$artifact_l2_sensitivity,
      if (identical(x$adjacency, "add_remove_patient")) 1 else sqrt(2),
      2048) &&
    isTRUE(valid_unit_policy) &&
    identical(x$source_values_exposed, FALSE) &&
    identical(x$intermediate_values_exposed, FALSE) &&
    identical(x$clipped_coordinates, NA_integer_) &&
    identical(x$clamp_activation_disclosed, FALSE) && is.finite(radius) &&
    .dsvert_dp_num_equal(x$accuracy_simultaneous_95_abs, radius) &&
    identical(x$accuracy_simultaneous_confidence, 0.95) &&
    .dsvert_dp_table_vector_accuracy_method_is_valid(x, profile) &&
    identical(x$accuracy_additional_privacy_cost, c(epsilon = 0, delta = 0)) &&
    is.list(x$release_provenance) &&
    identical(x$release_provenance$version,
              "dsvert-stateless-synopsis-public-provenance-v1") &&
    identical(x$release_provenance$durable_replay, TRUE) &&
    identical(x$composition_rule,
              "one_sticky_release_per_canonical_signed_artifact")
  if (!isTRUE(valid)) {
    stop("x must be a released, validated ds.vertDPContingency object",
         call. = FALSE)
  }
  x
}

.dsvert_dp_table_contract <- function(x) {
  table <- if (is.list(x)) x$table else NULL
  if (inherits(x, "ds.vertDPContingency") && is.list(x) &&
      identical(
        x$backend, "exact_signed_Ring128_canonical_synopsis_vector")) {
    return(.dsvert_dp_table_synopsis_contract(x))
  }
  if (inherits(x, "ds.vertDPContingency") && is.list(x) &&
      identical(x$backend, "exact_signed_Ring128_global_vector")) {
    profile <- .dsvert_dp_table_vector_profile(x)
    finite_scalar <- function(value, lower = -Inf, upper = Inf,
                              lower_open = FALSE) {
      is.numeric(value) && length(value) == 1L && !is.na(value) &&
        is.finite(value) &&
        (if (lower_open) value > lower else value >= lower) &&
        value <= upper
    }
    valid_table <- is.matrix(table) && is.numeric(table) &&
      length(table) > 0L && !anyNA(table) && all(is.finite(table)) &&
      all(table >= 0) && finite_scalar(x$coordinate_maximum, 1) &&
      all(table <= x$coordinate_maximum) &&
      identical(as.integer(nrow(table)), x$nrow) &&
      identical(as.integer(ncol(table)), x$ncol) &&
      identical(unname(rownames(table)), unname(x$row_levels)) &&
      identical(unname(colnames(table)), unname(x$col_levels)) &&
      identical(unname(as.numeric(table)), unname(x$counts))
    valid_epochs <- is.numeric(x$privacy_epochs) &&
      length(x$privacy_epochs) >= 2L && !anyNA(x$privacy_epochs) &&
      all(is.finite(x$privacy_epochs)) && all(x$privacy_epochs >= 1) &&
      all(x$privacy_epochs == floor(x$privacy_epochs)) &&
      is.character(x$noise_key_ids) &&
      length(x$noise_key_ids) == length(x$privacy_epochs) &&
      !anyNA(x$noise_key_ids) && all(nzchar(x$noise_key_ids))
    radius <- tryCatch(
      .dsvert_dp_vector_table_radius(x, 0.95),
      error = function(error) NA_real_)
    selected_sensitivity <- if (isTRUE(profile$gaussian)) {
      x$l2_sensitivity
    } else {
      x$l1_sensitivity
    }
    valid <- isTRUE(x$released) && valid_table && valid_epochs &&
      is.list(profile) &&
      identical(x$implementation, profile$backend) &&
      identical(x$sampler, profile$sampler) &&
      identical(x$mechanism, profile$release_mechanism) &&
      identical(
        x$randomness,
        paste("independent pinned-peer HKDF-SHA256/ChaCha20 streams;",
              "no analyst-controlled seed")) &&
      identical(x$sticky_noise, "immutable_capsule_durable_replay_v3") &&
      identical(x$sticky_replay, TRUE) &&
      identical(x$sensitivity_scope,
                "complete_signed_biomedical_capsule_vector") &&
      finite_scalar(x$epsilon, 0, .DSVERT_DP_MAXIMUM_EPSILON, TRUE) &&
      finite_scalar(x$delta, 0, 1) && x$delta < 1 &&
      finite_scalar(x$sensitivity, 0, Inf, TRUE) &&
      identical(x$sensitivity_norm,
                if (isTRUE(profile$gaussian)) "l2" else "l1") &&
      finite_scalar(x$l1_sensitivity, 0, Inf, TRUE) &&
      finite_scalar(x$l2_sensitivity, 0, Inf, TRUE) &&
      .dsvert_dp_num_equal(selected_sensitivity, x$sensitivity, 2048) &&
      x$adjacency %in% c(
        "add_remove_patient", "replace_one_fixed_cohort") &&
      .dsvert_dp_num_equal(
        x$artifact_l1_sensitivity,
        if (identical(x$adjacency, "add_remove_patient")) 1 else 2) &&
      .dsvert_dp_num_equal(
        x$artifact_l2_sensitivity,
        if (identical(x$adjacency, "add_remove_patient")) 1 else sqrt(2),
        2048) &&
      identical(x$unit_aggregation_policy,
                "consistent_joint_cell_else_exclude_v1") &&
      identical(x$history_gate, TRUE) &&
      identical(x$request_limit, FALSE) &&
      identical(x$operation_limit, TRUE) &&
      identical(x$source_values_exposed, FALSE) &&
      identical(x$intermediate_values_exposed, FALSE) &&
      identical(x$clipped_coordinates, NA_integer_) &&
      identical(x$clamp_activation_disclosed, FALSE) &&
      is.finite(radius) &&
      .dsvert_dp_num_equal(x$accuracy_simultaneous_95_abs, radius) &&
      identical(x$accuracy_simultaneous_confidence, 0.95) &&
      .dsvert_dp_table_vector_accuracy_method_is_valid(x, profile) &&
      identical(x$accuracy_additional_privacy_cost,
                c(epsilon = 0, delta = 0))
    if (!isTRUE(valid)) {
      stop("x must be a released, validated ds.vertDPContingency object",
           call. = FALSE)
    }
    return(x)
  }
  stop("x must be a released, validated ds.vertDPContingency object",
       call. = FALSE)
}

.dsvert_dp_vector_table_radius <- function(x, level) {
  profile <- .dsvert_dp_table_vector_profile(x)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1 ||
      !is.numeric(x$output_lattice_scale) ||
      length(x$output_lattice_scale) != 1L ||
      is.na(x$output_lattice_scale) || !is.finite(x$output_lattice_scale) ||
      x$output_lattice_scale <= 0 ||
      !is.numeric(x$l1_sensitivity) || length(x$l1_sensitivity) != 1L ||
      is.na(x$l1_sensitivity) || !is.finite(x$l1_sensitivity) ||
      x$l1_sensitivity <= 0 ||
      !is.numeric(x$l2_sensitivity) || length(x$l2_sensitivity) != 1L ||
      is.na(x$l2_sensitivity) || !is.finite(x$l2_sensitivity) ||
      x$l2_sensitivity <= 0 || !is.list(profile)) {
    stop("Invalid vector-table uncertainty contract", call. = FALSE)
  }
  total_coordinates <- if (.dsvert_dp_is_integer(
      x$synopsis_coordinate_count, 1, .DSVERT_DP_MAX_COORDINATES)) {
    x$synopsis_coordinate_count
  } else x$capsule_coordinate_count
  if (!.dsvert_dp_is_integer(total_coordinates, 1,
                             .DSVERT_DP_MAX_COORDINATES)) {
    request <- tryCatch(
      x$mechanism_selection$gaussian_calibration_request,
      error = function(error) NULL)
    total_coordinates <- if (is.list(request)) {
      request$total_coordinate_count
    } else {
      length(x$table)
    }
  }
  pseudo_release <- list(
    epsilon = x$epsilon, implementation_delta = x$implementation_delta,
    mechanism = x$mechanism, backend = x$implementation,
    mechanism_plan = x$mechanism_plan, plan_sha256 = x$plan_sha256,
    backend_selection = x$backend_selection,
    backend_assessment = x$backend_assessment,
    manifest_sha256 = x$manifest_sha256)
  capsule_mechanism <- if (!is.null(x$synopsis_mechanism)) {
    x$synopsis_mechanism
  } else x$capsule_mechanism
  if (is.null(capsule_mechanism)) {
    capsule_mechanism <- if (isTRUE(profile$gaussian)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM
    } else {
      "discrete-laplace"
    }
  }
  pseudo_manifest <- list(workload = list(
    coordinate_count = total_coordinates,
    capsule_mechanism = capsule_mechanism,
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
    pseudo_release, pseudo_manifest, coordinate_count = length(x$table),
    confidence = level, maximum_error = x$coordinate_maximum)$radius
}

.dsvert_dp_table_simultaneous_radius <- function(x, level) {
  .dsvert_dp_vector_table_radius(x, level)
}

.dsvert_dp_table_coverage_method <- function(x) {
  profile <- .dsvert_dp_table_vector_profile(x)
  if (isTRUE(profile$gaussian)) {
    return(paste(
      "Signed fixed-work dyadic discrete-Gaussian plan v2 confidence-specific",
      "subgaussian bound; signed tail and CDF TV transfers are deducted from",
      "the requested accuracy probability; the published 95% radius is used",
      "unchanged at level 0.95"))
  }
  if (isTRUE(profile$exact_gc)) {
    return(paste(
      "Exact ideal one-draw two-sided-geometric tail with union bound;",
      "the signed vector sampler TV bound is deducted once"))
  }
  paste(
    "Exact ideal two-sided-geometric convolution tail with union bound;",
    "both pinned peers' finite-sampler TV bounds are deducted")
}

.dsvert_dp_table_published_accuracy_matches <- function(x, radius) {
  profile <- .dsvert_dp_table_vector_profile(x)
  valid_method <- is.list(profile) &&
    .dsvert_dp_table_vector_accuracy_method_is_valid(x, profile)
  is.numeric(x$accuracy_simultaneous_95_abs) &&
    length(x$accuracy_simultaneous_95_abs) == 1L &&
    !is.na(x$accuracy_simultaneous_95_abs) &&
    is.finite(x$accuracy_simultaneous_95_abs) &&
    .dsvert_dp_num_equal(x$accuracy_simultaneous_95_abs, radius) &&
    isTRUE(valid_method)
}

.dsvert_dp_dimension_index <- function(value, labels, size, what) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    if (is.null(labels) || !value %in% labels) {
      stop("Unknown ", what, " level: '", value, "'", call. = FALSE)
    }
    return(match(value, labels))
  }
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) ||
      value < 1 || value > size) {
    stop(what, " must identify one table dimension by name or index",
         call. = FALSE)
  }
  as.integer(value)
}

.dsvert_dp_count_box <- function(table, radius) {
  lower <- table - radius
  lower[lower < 0] <- 0
  list(lower = lower, upper = table + radius)
}

.dsvert_dp_risk_bounds <- function(event_lower, event_upper,
                                   nonevent_lower, nonevent_upper) {
  lower_total <- event_lower + nonevent_upper
  upper_total <- event_upper + nonevent_lower
  lower <- if (lower_total > 0) event_lower / lower_total else 0
  upper <- if (upper_total > 0) event_upper / upper_total else 1
  c(lower = lower, upper = upper)
}

.dsvert_dp_attributable_fraction_exposed_region <- function(
    exposed_risk, unexposed_risk) {
  valid <- function(value) {
    is.numeric(value) && length(value) == 2L &&
      identical(names(value), c("lower", "upper")) &&
      !anyNA(value) && all(is.finite(value)) &&
      value[["lower"]] >= 0 && value[["upper"]] <= 1 &&
      value[["lower"]] <= value[["upper"]]
  }
  if (!valid(exposed_risk) || !valid(unexposed_risk)) {
    stop("The DP risk regions are invalid", call. = FALSE)
  }
  exposed_has_positive <- exposed_risk[["upper"]] > 0
  includes_non_estimable <- exposed_risk[["lower"]] == 0
  interval <- if (!exposed_has_positive) {
    c(lower = -Inf, upper = 1)
  } else {
    c(
      lower = if (exposed_risk[["lower"]] > 0) {
        1 - unexposed_risk[["upper"]] / exposed_risk[["lower"]]
      } else if (unexposed_risk[["upper"]] > 0) {
        -Inf
      } else {
        1
      },
      upper = 1 - unexposed_risk[["lower"]] /
        exposed_risk[["upper"]])
  }
  list(
    interval = interval,
    includes_non_estimable = includes_non_estimable,
    has_estimable_values = exposed_has_positive)
}

.dsvert_dp_number_needed_region <- function(risk_difference) {
  if (!is.numeric(risk_difference) || length(risk_difference) != 2L ||
      !identical(names(risk_difference), c("lower", "upper")) ||
      anyNA(risk_difference) || any(!is.finite(risk_difference)) ||
      risk_difference[["lower"]] > risk_difference[["upper"]]) {
    stop("The DP risk-difference region is invalid", call. = FALSE)
  }
  lower <- risk_difference[["lower"]]
  upper <- risk_difference[["upper"]]
  regions <- list()
  if (lower < 0) {
    minimum_effect <- if (upper < 0) abs(upper) else 0
    regions$benefit <- c(
      lower = 1 / abs(lower),
      upper = if (minimum_effect > 0) 1 / minimum_effect else Inf)
  }
  if (upper > 0) {
    minimum_effect <- if (lower > 0) lower else 0
    regions$harm <- c(
      lower = 1 / upper,
      upper = if (minimum_effect > 0) 1 / minimum_effect else Inf)
  }
  list(
    regions = regions,
    includes_infinite = lower <= 0 && upper >= 0,
    possible_directions = names(regions))
}

.dsvert_dp_diagnostic_values <- function(cells) {
  expected <- c(
    "true_positive", "false_negative", "false_positive", "true_negative")
  if (!is.numeric(cells) || length(cells) != 4L ||
      is.null(names(cells)) || anyNA(names(cells)) ||
      anyDuplicated(names(cells)) || !setequal(names(cells), expected) ||
      anyNA(cells) || any(!is.finite(cells)) || any(cells < 0)) {
    stop("The diagnostic 2-by-2 cells are invalid", call. = FALSE)
  }
  cells <- unname(cells[expected])
  names(cells) <- expected
  tp <- cells[["true_positive"]]
  fn <- cells[["false_negative"]]
  fp <- cells[["false_positive"]]
  tn <- cells[["true_negative"]]
  disease_total <- tp + fn
  nondisease_total <- fp + tn
  test_positive_total <- tp + fp
  test_negative_total <- fn + tn
  total <- disease_total + nondisease_total
  proportion <- function(numerator, denominator) {
    if (denominator > 0) numerator / denominator else NA_real_
  }
  extended_ratio <- function(numerator, denominator) {
    if (is.na(numerator) || is.na(denominator)) return(NA_real_)
    if (denominator > 0) return(numerator / denominator)
    if (numerator > 0) Inf else NA_real_
  }
  sensitivity <- proportion(tp, disease_total)
  specificity <- proportion(tn, nondisease_total)
  false_positive_rate <- proportion(fp, nondisease_total)
  false_negative_rate <- proportion(fn, disease_total)
  f1_denominator <- 2 * tp + fp + fn
  odds_numerator <- tp * tn
  odds_denominator <- fp * fn
  c(
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = proportion(tp, test_positive_total),
    npv = proportion(tn, test_negative_total),
    prevalence = proportion(disease_total, total),
    accuracy = proportion(tp + tn, total),
    balanced_accuracy = if (disease_total > 0 && nondisease_total > 0) {
      (sensitivity + specificity) / 2
    } else {
      NA_real_
    },
    f1_score = proportion(2 * tp, f1_denominator),
    lr_positive = extended_ratio(sensitivity, false_positive_rate),
    lr_negative = extended_ratio(false_negative_rate, specificity),
    diagnostic_odds_ratio = extended_ratio(
      odds_numerator, odds_denominator))
}

.dsvert_dp_diagnostic_regions <- function(lower, upper) {
  expected <- c(
    "true_positive", "false_negative", "false_positive", "true_negative")
  valid <- function(value) {
    is.numeric(value) && length(value) == 4L &&
      !is.null(names(value)) && !anyNA(names(value)) &&
      !anyDuplicated(names(value)) && setequal(names(value), expected) &&
      !anyNA(value) && all(is.finite(value)) && all(value >= 0)
  }
  if (!valid(lower) || !valid(upper)) {
    stop("The diagnostic DP count box is invalid", call. = FALSE)
  }
  lower <- unname(lower[expected])
  upper <- unname(upper[expected])
  names(lower) <- names(upper) <- expected
  if (any(lower > upper)) {
    stop("The diagnostic DP count box is invalid", call. = FALSE)
  }
  cell_range <- function(name) {
    c(lower = lower[[name]], upper = upper[[name]])
  }
  add_ranges <- function(...) {
    values <- list(...)
    c(lower = sum(vapply(values, `[[`, numeric(1L), "lower")),
      upper = sum(vapply(values, `[[`, numeric(1L), "upper")))
  }
  probability <- function(numerator, complement) {
    list(
      interval = .dsvert_dp_risk_bounds(
        numerator[["lower"]], numerator[["upper"]],
        complement[["lower"]], complement[["upper"]]),
      zero = numerator[["lower"]] == 0 &&
        complement[["upper"]] > 0,
      positive = numerator[["upper"]] > 0,
      non_estimable = numerator[["lower"]] == 0 &&
        complement[["lower"]] == 0)
  }
  ratio_interval <- function(numerator, denominator) {
    c(
      lower = if (denominator[["upper"]] > 0) {
        numerator[["lower"]] / denominator[["upper"]]
      } else if (numerator[["lower"]] > 0) Inf else 0,
      upper = if (denominator[["lower"]] > 0) {
        numerator[["upper"]] / denominator[["lower"]]
      } else if (numerator[["upper"]] > 0) Inf else 0)
  }
  ratio <- function(numerator, denominator) {
    list(
      interval = ratio_interval(
        numerator$interval, denominator$interval),
      zero = numerator$zero && denominator$positive,
      infinite = numerator$positive && denominator$zero,
      non_estimable = numerator$non_estimable ||
        denominator$non_estimable ||
        (numerator$zero && denominator$zero))
  }

  tp <- cell_range("true_positive")
  fn <- cell_range("false_negative")
  fp <- cell_range("false_positive")
  tn <- cell_range("true_negative")
  diseased <- add_ranges(tp, fn)
  nondiseased <- add_ranges(fp, tn)
  correct <- add_ranges(tp, tn)
  incorrect <- add_ranges(fp, fn)

  sensitivity <- probability(tp, fn)
  specificity <- probability(tn, fp)
  ppv <- probability(tp, fp)
  npv <- probability(tn, fn)
  prevalence <- probability(diseased, nondiseased)
  accuracy <- probability(correct, incorrect)
  balanced_accuracy <- list(
    interval = c(
      lower = (sensitivity$interval[["lower"]] +
                 specificity$interval[["lower"]]) / 2,
      upper = (sensitivity$interval[["upper"]] +
                 specificity$interval[["upper"]]) / 2),
    zero = tp[["lower"]] == 0 && tn[["lower"]] == 0 &&
      fn[["upper"]] > 0 && fp[["upper"]] > 0,
    positive = tp[["upper"]] > 0 || tn[["upper"]] > 0,
    non_estimable = sensitivity$non_estimable ||
      specificity$non_estimable)
  f1_score <- probability(
    c(lower = 2 * tp[["lower"]], upper = 2 * tp[["upper"]]),
    add_ranges(fp, fn))
  lr_positive <- ratio(sensitivity, probability(fp, tn))
  lr_negative <- ratio(probability(fn, tp), specificity)

  odds_numerator <- c(
    lower = tp[["lower"]] * tn[["lower"]],
    upper = tp[["upper"]] * tn[["upper"]])
  odds_denominator <- c(
    lower = fp[["lower"]] * fn[["lower"]],
    upper = fp[["upper"]] * fn[["upper"]])
  diagnostic_odds_ratio <- list(
    interval = ratio_interval(odds_numerator, odds_denominator),
    zero = (tp[["lower"]] == 0 || tn[["lower"]] == 0) &&
      fp[["upper"]] > 0 && fn[["upper"]] > 0,
    infinite = tp[["upper"]] > 0 && tn[["upper"]] > 0 &&
      (fp[["lower"]] == 0 || fn[["lower"]] == 0),
    non_estimable =
      (tp[["lower"]] == 0 || tn[["lower"]] == 0) &&
      (fp[["lower"]] == 0 || fn[["lower"]] == 0))

  probabilities <- list(
    sensitivity = sensitivity, specificity = specificity,
    ppv = ppv, npv = npv, prevalence = prevalence, accuracy = accuracy,
    balanced_accuracy = balanced_accuracy, f1_score = f1_score)
  ratios <- list(
    lr_positive = lr_positive, lr_negative = lr_negative,
    diagnostic_odds_ratio = diagnostic_odds_ratio)
  all_metrics <- c(probabilities, ratios)
  regions <- lapply(all_metrics, `[[`, "interval")
  flags <- data.frame(
    includes_zero = vapply(all_metrics, `[[`, logical(1L), "zero"),
    includes_infinite = c(
      rep(FALSE, length(probabilities)),
      vapply(ratios, `[[`, logical(1L), "infinite")),
    includes_non_estimable = vapply(
      all_metrics, `[[`, logical(1L), "non_estimable"),
    row.names = names(all_metrics), check.names = FALSE)
  region_types <- vapply(seq_len(nrow(flags)), function(index) {
    if (flags$includes_non_estimable[[index]] &&
        flags$includes_infinite[[index]]) {
      "includes_non_estimable_and_infinite_boundary"
    } else if (flags$includes_non_estimable[[index]]) {
      "includes_non_estimable"
    } else if (flags$includes_infinite[[index]] &&
               is.infinite(regions[[index]][["lower"]])) {
      "infinite_boundary"
    } else if (flags$includes_infinite[[index]]) {
      "unbounded_above"
    } else {
      "finite"
    }
  }, character(1L))
  names(region_types) <- names(all_metrics)
  list(regions = regions, flags = flags, region_types = region_types)
}

.dsvert_dp_diagnostic_point_status <- function(cells, estimates) {
  tp <- cells[["true_positive"]]
  fn <- cells[["false_negative"]]
  fp <- cells[["false_positive"]]
  tn <- cells[["true_negative"]]
  missing_reason <- c(
    sensitivity = "non_estimable_zero_disease_total",
    specificity = "non_estimable_zero_nondisease_total",
    ppv = "non_estimable_zero_test_positive_total",
    npv = "non_estimable_zero_test_negative_total",
    prevalence = "non_estimable_zero_total",
    accuracy = "non_estimable_zero_total",
    balanced_accuracy = "non_estimable_zero_reference_total",
    f1_score = "non_estimable_zero_f1_denominator",
    lr_positive = if (tp + fn == 0 || fp + tn == 0) {
      "non_estimable_zero_reference_total"
    } else "non_estimable_undefined_zero_over_zero",
    lr_negative = if (tp + fn == 0 || fp + tn == 0) {
      "non_estimable_zero_reference_total"
    } else "non_estimable_undefined_zero_over_zero",
    diagnostic_odds_ratio = "non_estimable_undefined_zero_over_zero")
  status <- vapply(names(estimates), function(name) {
    value <- estimates[[name]]
    if (is.na(value)) return(missing_reason[[name]])
    if (is.infinite(value)) return("boundary_infinite")
    if (value == 0) return("boundary_zero")
    if (name %in% c(
          "sensitivity", "specificity", "ppv", "npv",
          "prevalence", "accuracy", "balanced_accuracy", "f1_score") &&
        value == 1) {
      return("boundary_one")
    }
    "ok"
  }, character(1L))
  status
}

#' Diagnostic-accuracy measures with simultaneous DP-noise regions
#'
#' Purely post-processes one validated DP 2-by-2 histogram. Rows must encode
#' disease status and columns must encode the diagnostic-test result. The
#' positive disease row and positive test column are mandatory explicit
#' choices. No continuity correction, server call, p-value, or classical
#' sampling confidence interval is used. Measures include sensitivity,
#' specificity, predictive values, prevalence, accuracy, balanced accuracy,
#' F1 score, likelihood ratios, and the diagnostic odds ratio.
#'
#' @param x A released `ds.vertDPContingency` with a 2-by-2 table.
#' @param disease_positive Positive disease-status row by name or index.
#' @param test_positive Positive diagnostic-test column by name or index.
#' @param level Simultaneous coverage for DP mechanism noise.
#' @return A `ds.vertDPDiagnostic2x2` object containing typed estimates and
#'   simultaneous DP-mechanism regions. No server call is made.
#' @export
ds.vertDPDiagnostic2x2 <- function(
    x, disease_positive, test_positive, level = 0.95) {
  x <- .dsvert_dp_table_contract(x)
  if (!identical(dim(x$table), c(2L, 2L))) {
    stop("x must contain exactly a 2-by-2 DP table", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }
  disease_positive <- .dsvert_dp_dimension_index(
    disease_positive, rownames(x$table), 2L, "disease_positive")
  test_positive <- .dsvert_dp_dimension_index(
    test_positive, colnames(x$table), 2L, "test_positive")
  disease_negative <- setdiff(1:2, disease_positive)
  test_negative <- setdiff(1:2, test_positive)
  cells <- c(
    true_positive = x$table[disease_positive, test_positive],
    false_negative = x$table[disease_positive, test_negative],
    false_positive = x$table[disease_negative, test_positive],
    true_negative = x$table[disease_negative, test_negative])
  oriented_table <- matrix(
    cells, nrow = 2L, byrow = TRUE,
    dimnames = list(
      disease = c("positive", "negative"),
      test = c("positive", "negative")))

  radius <- .dsvert_dp_table_simultaneous_radius(x, level)
  if (identical(level, 0.95) &&
      !.dsvert_dp_table_published_accuracy_matches(x, radius)) {
    stop("x does not carry a valid simultaneous DP accuracy certificate",
         call. = FALSE)
  }
  box <- .dsvert_dp_count_box(oriented_table, radius)
  lower <- c(
    true_positive = box$lower[1L, 1L],
    false_negative = box$lower[1L, 2L],
    false_positive = box$lower[2L, 1L],
    true_negative = box$lower[2L, 2L])
  upper <- c(
    true_positive = box$upper[1L, 1L],
    false_negative = box$upper[1L, 2L],
    false_positive = box$upper[2L, 1L],
    true_negative = box$upper[2L, 2L])
  estimates <- .dsvert_dp_diagnostic_values(cells)
  point_status <- .dsvert_dp_diagnostic_point_status(cells, estimates)
  region <- .dsvert_dp_diagnostic_regions(lower, upper)
  axis_level <- function(labels, index) {
    if (is.null(labels)) NA_character_ else labels[[index]]
  }
  non_estimable <- grepl("^non_estimable", point_status)
  boundary <- grepl("^boundary", point_status)
  status <- if (all(non_estimable)) {
    "non_estimable"
  } else if (any(non_estimable)) {
    "partially_non_estimable"
  } else if (any(boundary)) {
    "boundary"
  } else {
    "ok"
  }
  result <- list(
    status = status,
    orientation = list(
      row_role = "disease_status",
      column_role = "test_result",
      disease_positive = list(
        index = disease_positive,
        level = axis_level(rownames(x$table), disease_positive)),
      disease_negative = list(
        index = disease_negative,
        level = axis_level(rownames(x$table), disease_negative)),
      test_positive = list(
        index = test_positive,
        level = axis_level(colnames(x$table), test_positive)),
      test_negative = list(
        index = test_negative,
        level = axis_level(colnames(x$table), test_negative))),
    oriented_table = oriented_table,
    estimates = estimates,
    point_status = point_status,
    count_lower = box$lower,
    count_upper = box$upper,
    mechanism_regions = region$regions,
    mechanism_region_flags = region$flags,
    mechanism_region_types = region$region_types,
    level = level,
    simultaneous_radius = radius,
    coverage_method = .dsvert_dp_table_coverage_method(x),
    continuity_correction = "none",
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    inferential_scope = paste(
      "Finite-dataset diagnostic accuracy from one DP-noised table;",
      "no sampling confidence interval is provided; no hypothesis test or",
      "p-value is provided"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPDiagnostic2x2", "list")
  result
}

#' Epidemiological measures with simultaneous DP-noise regions
#'
#' Post-processes one already released 2-by-2 DP histogram. The regions cover
#' the certified granular-Laplace or approximate-Gaussian mechanism noise with
#' at least `level` probability by a union bound over all cells. They are not
#' sampling confidence intervals and do not make an ordinary chi-square or
#' Fisher p-value valid for noisy counts. Returned finite-snapshot measures
#' include group and population risks, risk difference, risk and odds ratios,
#' attributable fractions among the exposed and in the population, and a
#' direction-typed number needed to benefit or harm.
#'
#' @param x A released `ds.vertDPContingency` with a 2-by-2 table.
#' @param exposed,event Exposed row and event column by index or level name.
#' @param level Simultaneous coverage for DP mechanism noise.
#' @return A `ds.vertDPEpi2x2` object. No server call is made.
#' @export
ds.vertDPEpi2x2 <- function(x, exposed = 2L, event = 2L, level = 0.95) {
  x <- .dsvert_dp_table_contract(x)
  if (!identical(dim(x$table), c(2L, 2L))) {
    stop("x must contain exactly a 2-by-2 DP table", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }
  exposed <- .dsvert_dp_dimension_index(
    exposed, rownames(x$table), 2L, "exposed")
  event <- .dsvert_dp_dimension_index(
    event, colnames(x$table), 2L, "event")
  unexposed <- setdiff(1:2, exposed)
  nonevent <- setdiff(1:2, event)

  radius <- .dsvert_dp_table_simultaneous_radius(x, level)
  if (identical(level, 0.95) &&
      !.dsvert_dp_table_published_accuracy_matches(x, radius)) {
    stop("x does not carry a valid simultaneous DP accuracy certificate",
         call. = FALSE)
  }
  box <- .dsvert_dp_count_box(x$table, radius)
  exposed_risk <- .dsvert_dp_risk_bounds(
    box$lower[exposed, event], box$upper[exposed, event],
    box$lower[exposed, nonevent], box$upper[exposed, nonevent])
  unexposed_risk <- .dsvert_dp_risk_bounds(
    box$lower[unexposed, event], box$upper[unexposed, event],
    box$lower[unexposed, nonevent], box$upper[unexposed, nonevent])
  population_risk <- .dsvert_dp_risk_bounds(
    sum(box$lower[, event]), sum(box$upper[, event]),
    sum(box$lower[, nonevent]), sum(box$upper[, nonevent]))

  risk_ratio <- c(
    lower = if (unexposed_risk[["upper"]] > 0) {
      exposed_risk[["lower"]] / unexposed_risk[["upper"]]
    } else 0,
    upper = if (unexposed_risk[["lower"]] > 0) {
      exposed_risk[["upper"]] / unexposed_risk[["lower"]]
    } else Inf)
  a <- c(box$lower[exposed, event], box$upper[exposed, event])
  b <- c(box$lower[exposed, nonevent], box$upper[exposed, nonevent])
  c0 <- c(box$lower[unexposed, event], box$upper[unexposed, event])
  d <- c(box$lower[unexposed, nonevent], box$upper[unexposed, nonevent])
  odds_ratio <- c(
    lower = if (b[[2L]] > 0 && c0[[2L]] > 0) {
      a[[1L]] * d[[1L]] / (b[[2L]] * c0[[2L]])
    } else 0,
    upper = if (b[[1L]] > 0 && c0[[1L]] > 0) {
      a[[2L]] * d[[2L]] / (b[[1L]] * c0[[1L]])
    } else Inf)

  a_point <- x$table[exposed, event]
  b_point <- x$table[exposed, nonevent]
  c_point <- x$table[unexposed, event]
  d_point <- x$table[unexposed, nonevent]
  exposed_total <- a_point + b_point
  unexposed_total <- c_point + d_point
  population_total <- sum(x$table)
  population_risk_point <- if (population_total > 0) {
    (a_point + c_point) / population_total
  } else NULL
  risk_exposed_point <- if (exposed_total > 0) {
    a_point / exposed_total
  } else NULL
  risk_unexposed_point <- if (unexposed_total > 0) {
    c_point / unexposed_total
  } else NULL
  groups_estimable <- !is.null(risk_exposed_point) &&
    !is.null(risk_unexposed_point)
  risk_difference_point <- if (groups_estimable) {
    risk_exposed_point - risk_unexposed_point
  } else NULL
  risk_ratio_point <- if (!groups_estimable) {
    NULL
  } else if (risk_unexposed_point > 0) {
    risk_exposed_point / risk_unexposed_point
  } else if (risk_exposed_point > 0) {
    Inf
  } else NULL
  odds_numerator <- a_point * d_point
  odds_denominator <- b_point * c_point
  odds_ratio_point <- if (!groups_estimable) {
    NULL
  } else if (odds_denominator > 0) {
    odds_numerator / odds_denominator
  } else if (odds_numerator > 0) {
    Inf
  } else NULL
  attributable_fraction_point <- if (!groups_estimable ||
      risk_exposed_point <= 0) {
    NULL
  } else {
    risk_difference_point / risk_exposed_point
  }
  population_attributable_fraction_point <- if (
      is.null(population_risk_point) || population_risk_point <= 0 ||
      is.null(risk_unexposed_point)) {
    NULL
  } else {
    (population_risk_point - risk_unexposed_point) /
      population_risk_point
  }
  point_values <- list(
    risk_exposed = risk_exposed_point,
    risk_unexposed = risk_unexposed_point,
    population_risk = population_risk_point,
    risk_difference = risk_difference_point,
    risk_ratio = risk_ratio_point,
    odds_ratio = odds_ratio_point,
    attributable_fraction_exposed = attributable_fraction_point,
    population_attributable_fraction =
      population_attributable_fraction_point)
  point_status <- list(
    risk_exposed = if (is.null(risk_exposed_point)) {
      "non_estimable_zero_group_total"
    } else "ok",
    risk_unexposed = if (is.null(risk_unexposed_point)) {
      "non_estimable_zero_group_total"
    } else "ok",
    population_risk = if (is.null(population_risk_point)) {
      "non_estimable_zero_population_total"
    } else "ok",
    risk_difference = if (is.null(risk_difference_point)) {
      "non_estimable_zero_group_total"
    } else "ok",
    risk_ratio = if (!groups_estimable) {
      "non_estimable_zero_group_total"
    } else if (is.null(risk_ratio_point)) {
      "non_estimable_undefined_ratio"
    } else if (is.infinite(risk_ratio_point)) {
      "boundary_infinite"
    } else "ok",
    odds_ratio = if (!groups_estimable) {
      "non_estimable_zero_group_total"
    } else if (is.null(odds_ratio_point)) {
      "non_estimable_undefined_ratio"
    } else if (is.infinite(odds_ratio_point)) {
      "boundary_infinite"
    } else "ok",
    attributable_fraction_exposed = if (!groups_estimable) {
      "non_estimable_zero_group_total"
    } else if (is.null(attributable_fraction_point)) {
      "non_estimable_zero_exposed_risk"
    } else {
      "ok"
    },
    population_attributable_fraction = if (is.null(risk_unexposed_point)) {
      "non_estimable_zero_unexposed_total"
    } else if (is.null(population_risk_point)) {
      "non_estimable_zero_population_total"
    } else if (is.null(population_attributable_fraction_point)) {
      "non_estimable_zero_population_risk"
    } else {
      "ok"
    })
  point_estimates <- Map(function(value, status) {
    if (identical(status, "ok")) value else NULL
  }, point_values, point_status)
  risk_difference_region <- c(
    lower = exposed_risk[["lower"]] - unexposed_risk[["upper"]],
    upper = exposed_risk[["upper"]] - unexposed_risk[["lower"]])
  attributable_fraction_region <-
    .dsvert_dp_attributable_fraction_exposed_region(
      exposed_risk, unexposed_risk)
  population_attributable_fraction_region <-
    .dsvert_dp_attributable_fraction_exposed_region(
      population_risk, unexposed_risk)
  unexposed_total_upper <-
    box$upper[unexposed, event] + box$upper[unexposed, nonevent]
  if (unexposed_total_upper <= 0) {
    population_attributable_fraction_region <- list(
      interval = c(lower = -Inf, upper = 1),
      includes_non_estimable = TRUE,
      has_estimable_values = FALSE)
  }
  regions <- list(
    risk_exposed = exposed_risk,
    risk_unexposed = unexposed_risk,
    population_risk = population_risk,
    risk_difference = risk_difference_region,
    risk_ratio = risk_ratio,
    odds_ratio = odds_ratio,
    attributable_fraction_exposed =
      attributable_fraction_region$interval,
    population_attributable_fraction =
      population_attributable_fraction_region$interval)
  region_types <- vapply(regions, function(interval) {
    if (!is.numeric(interval) || length(interval) != 2L ||
        anyNA(interval) || interval[["lower"]] > interval[["upper"]]) {
      stop("DP mechanism interval construction failed", call. = FALSE)
    }
    if (all(is.finite(interval))) return("finite")
    if (is.finite(interval[["lower"]]) &&
        is.infinite(interval[["upper"]]) && interval[["upper"]] > 0) {
      return(if (interval[["lower"]] == 0) {
        "unbounded_nonnegative"
      } else {
        "unbounded_above"
      })
    }
    "unbounded"
  }, character(1L))
  all_point_ok <- all(vapply(
    point_status, identical, logical(1L), "ok"))
  any_boundary <- any(vapply(
    point_status, identical, logical(1L), "boundary_infinite"))
  number_needed_region <-
    .dsvert_dp_number_needed_region(risk_difference_region)
  number_needed_point <- if (is.null(risk_difference_point)) {
    NA_real_
  } else if (risk_difference_point == 0) {
    Inf
  } else {
    1 / abs(risk_difference_point)
  }
  number_needed_direction <- if (is.null(risk_difference_point)) {
    "non_estimable"
  } else if (risk_difference_point < 0) {
    "benefit"
  } else if (risk_difference_point > 0) {
    "harm"
  } else {
    "none"
  }
  result <- list(
    status = if (all_point_ok) {
      "ok"
    } else if (any_boundary) {
      "dp_point_on_ratio_boundary"
    } else {
      "dp_point_non_estimable"
    },
    noisy_table = x$table,
    count_lower = box$lower,
    count_upper = box$upper,
    point_estimates = point_estimates,
    point_status = point_status,
    mechanism_regions = regions,
    mechanism_region_types = region_types,
    attributable_fraction_exposed_region_includes_non_estimable =
      attributable_fraction_region$includes_non_estimable,
    attributable_fraction_exposed_region_has_estimable_values =
      attributable_fraction_region$has_estimable_values,
    population_attributable_fraction_region_includes_non_estimable =
      population_attributable_fraction_region$includes_non_estimable,
    population_attributable_fraction_region_has_estimable_values =
      population_attributable_fraction_region$has_estimable_values,
    number_needed = list(
      point_estimate = number_needed_point,
      point_direction = number_needed_direction,
      mechanism_regions = number_needed_region$regions,
      mechanism_region_includes_infinite =
        number_needed_region$includes_infinite,
      mechanism_region_possible_directions =
        number_needed_region$possible_directions,
      definition = paste(
        "absolute reciprocal of the finite-snapshot risk difference;",
        "benefit means exposure lowers event risk and harm means it raises",
        "event risk")),
    level = level,
    simultaneous_radius = radius,
    coverage_method = .dsvert_dp_table_coverage_method(x),
    uncertainty_scope = "DP mechanism noise only; sampling uncertainty excluded",
    inferential_scope = paste(
      "Finite-dataset effects from one DP-noised table; no sampling",
      "confidence interval or ordinary hypothesis test is provided"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPEpi2x2", "list")
  result
}

#' Direct standardisation with simultaneous DP-noise bounds
#'
#' Treats rows of an already released DP contingency table as strata and one
#' column as the event. Reference weights are public client-side inputs.
#'
#' @param x A released `ds.vertDPContingency` (strata by outcome).
#' @param standard_weights Non-negative public weights, optionally named by
#'   table row levels.
#' @param event Event column by name or index.
#' @param level Simultaneous coverage for all DP table cells.
#' @return A `ds.vertDPStandardization` object. No server call is made.
#' @export
ds.vertDPDirectStandardization <- function(
    x, standard_weights, event = 2L, level = 0.95) {
  x <- .dsvert_dp_table_contract(x)
  if (ncol(x$table) != 2L || nrow(x$table) < 1L) {
    stop("x must be a strata-by-binary-outcome DP table", call. = FALSE)
  }
  if (!is.numeric(standard_weights) ||
      length(standard_weights) != nrow(x$table) ||
      anyNA(standard_weights) || any(!is.finite(standard_weights)) ||
      any(standard_weights < 0) || !any(standard_weights > 0)) {
    stop("standard_weights must be finite, non-negative, and match the strata",
         call. = FALSE)
  }
  if (!is.null(names(standard_weights))) {
    if (is.null(rownames(x$table)) || anyDuplicated(names(standard_weights)) ||
        !setequal(names(standard_weights), rownames(x$table))) {
      stop("Named standard_weights must match the DP stratum levels",
           call. = FALSE)
    }
    standard_weights <- standard_weights[rownames(x$table)]
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }
  event <- .dsvert_dp_dimension_index(
    event, colnames(x$table), 2L, "event")
  nonevent <- setdiff(1:2, event)
  # Scale first so finite but very large public weights cannot overflow their
  # sum and silently normalise to all zeros.
  weight_scale <- max(standard_weights)
  scaled_weights <- as.numeric(standard_weights) / weight_scale
  scaled_total <- sum(scaled_weights)
  if (!is.finite(scaled_total) || scaled_total <= 0) {
    stop("standard_weights cannot be normalised safely", call. = FALSE)
  }
  weights <- scaled_weights / scaled_total
  radius <- .dsvert_dp_table_simultaneous_radius(x, level)
  if (identical(level, 0.95) &&
      !.dsvert_dp_table_published_accuracy_matches(x, radius)) {
    stop("x does not carry a valid simultaneous DP accuracy certificate",
         call. = FALSE)
  }
  box <- .dsvert_dp_count_box(x$table, radius)
  bounds <- t(vapply(seq_len(nrow(x$table)), function(i) {
    .dsvert_dp_risk_bounds(
      box$lower[i, event], box$upper[i, event],
      box$lower[i, nonevent], box$upper[i, nonevent])
  }, numeric(2L)))
  noisy_denominator <- rowSums(x$table)
  point_rates <- ifelse(
    noisy_denominator > 0, x$table[, event] / noisy_denominator, NA_real_)
  positive_weight <- weights > 0
  point <- if (any(!is.finite(point_rates[positive_weight]))) {
    NA_real_
  } else {
    sum(weights[positive_weight] * point_rates[positive_weight])
  }
  result <- list(
    status = if (is.finite(point)) "ok" else "dp_point_non_estimable",
    estimate = point,
    mechanism_region = c(
      lower = sum(weights * bounds[, "lower"]),
      upper = sum(weights * bounds[, "upper"])),
    stratum_estimates = point_rates,
    stratum_regions = bounds,
    count_lower = box$lower,
    count_upper = box$upper,
    weights = weights,
    level = level,
    simultaneous_radius = radius,
    coverage_method = .dsvert_dp_table_coverage_method(x),
    uncertainty_scope = "DP mechanism noise only; sampling uncertainty excluded",
    inferential_scope = paste(
      "Finite-dataset directly standardized risk from one DP-noised table;",
      "no sampling confidence interval or hypothesis test is provided"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPStandardization", "list")
  result
}

.dsvert_dp_indirect_box_region <- function(
    event_lower, event_upper, nonevent_lower, nonevent_upper,
    expected_rates) {
  values <- list(
    event_lower, event_upper, nonevent_lower, nonevent_upper,
    expected_rates)
  lengths <- vapply(values, length, integer(1L))
  if (!length(event_lower) || length(unique(lengths)) != 1L ||
      any(!vapply(values, is.numeric, logical(1L))) ||
      any(vapply(values, anyNA, logical(1L))) ||
      any(!is.finite(unlist(values, use.names = FALSE))) ||
      any(unlist(values[1:4], use.names = FALSE) < 0) ||
      any(event_lower > event_upper) ||
      any(nonevent_lower > nonevent_upper) ||
      any(expected_rates < 0 | expected_rates > 1) ||
      !any(expected_rates > 0)) {
    stop("The indirect-standardisation DP count box is invalid",
         call. = FALSE)
  }

  denominator_lower <- sum(
    expected_rates * (event_lower + nonevent_lower))
  denominator_upper <- sum(
    expected_rates * (event_upper + nonevent_upper))
  numerator_upper <- sum(event_upper)
  positive_rate <- expected_rates > 0
  includes_non_estimable <-
    all(event_lower == 0) && all(nonevent_lower[positive_rate] == 0)
  includes_infinite <-
    all(event_lower[positive_rate] == 0) &&
    all(nonevent_lower[positive_rate] == 0) &&
    any(event_upper[!positive_rate] > 0)

  if (denominator_upper <= 0) {
    return(list(
      interval = c(
        lower = 0,
        upper = if (numerator_upper > 0) Inf else 0),
      includes_non_estimable = includes_non_estimable,
      includes_infinite = numerator_upper > 0,
      estimable_values = FALSE,
      tightness = "no positive expected denominator exists in the count box"))
  }

  # For q >= 0, extrema of N - qD over the rectangular count box are
  # separable. Bisection of those monotone envelopes gives the exact extrema
  # over the continuous relaxation, which conservatively contains every
  # attainable integer table. Nonevents use their upper bound for the minimum
  # ratio and their lower bound for the maximum ratio.
  minimum_envelope <- function(q) {
    coefficients <- 1 - q * expected_rates
    selected_event <- ifelse(
      coefficients >= 0, event_lower, event_upper)
    sum(coefficients * selected_event) -
      q * sum(expected_rates * nonevent_upper)
  }
  maximum_envelope <- function(q) {
    coefficients <- 1 - q * expected_rates
    selected_event <- ifelse(
      coefficients >= 0, event_upper, event_lower)
    sum(coefficients * selected_event) -
      q * sum(expected_rates * nonevent_lower)
  }
  bisect <- function(lower, upper, envelope, lower_when_nonnegative) {
    for (iteration in seq_len(256L)) {
      middle <- lower / 2 + upper / 2
      if (!is.finite(middle) || middle <= lower || middle >= upper) break
      nonnegative <- envelope(middle) >= 0
      if (identical(nonnegative, lower_when_nonnegative)) {
        lower <- middle
      } else {
        upper <- middle
      }
    }
    c(lower = lower, upper = upper)
  }

  feasible_ratio <- numerator_upper / denominator_upper
  minimum_bracket <- bisect(
    0, feasible_ratio, minimum_envelope,
    lower_when_nonnegative = TRUE)
  minimum <- max(
    0, minimum_bracket[["lower"]] -
      64 * .Machine$double.eps *
        max(1, abs(minimum_bracket[["lower"]])))

  maximum <- if (denominator_lower > 0) {
    upper <- numerator_upper / denominator_lower
    bracket <- bisect(
      0, upper, maximum_envelope,
      lower_when_nonnegative = TRUE)
    bracket[["upper"]] + 64 * .Machine$double.eps *
      max(1, abs(bracket[["upper"]]))
  } else if (any(event_upper[!positive_rate] > 0)) {
    Inf
  } else if (numerator_upper == 0) {
    0
  } else {
    upper <- 1 / min(expected_rates[positive_rate])
    bracket <- bisect(
      0, upper, maximum_envelope,
      lower_when_nonnegative = TRUE)
    bracket[["upper"]] + 64 * .Machine$double.eps *
      max(1, abs(bracket[["upper"]]))
  }
  list(
    interval = c(lower = minimum, upper = maximum),
    includes_non_estimable = includes_non_estimable,
    includes_infinite = includes_infinite || is.infinite(maximum),
    estimable_values = TRUE,
    tightness = paste(
      "outward-bracketed linear-fractional extrema over the continuous",
      "rectangular count relaxation; conservative for integer tables"))
}

#' Differentially private indirect standardisation
#'
#' Computes an observed-to-expected ratio (for example, an SMR or SIR) from
#' one validated DP strata-by-binary-outcome table and public stratum-specific
#' expected event probabilities. It makes no server call and consumes no
#' additional privacy budget. The returned region propagates simultaneous DP
#' mechanism noise; it is not a Garwood, sampling, or population confidence
#' interval.
#'
#' @param x A released `ds.vertDPContingency` whose rows are strata and whose
#'   two columns are event/non-event outcomes.
#' @param expected_rates Public expected event probabilities, one per stratum.
#'   Named values are reordered to the released row levels.
#' @param event Event column by name or index.
#' @param level Simultaneous coverage for DP mechanism noise.
#' @return A `ds.vertDPIndirectStandardization` object. No server call is made.
#' @export
ds.vertDPIndirectStandardization <- function(
    x, expected_rates, event = 2L, level = 0.95) {
  x <- .dsvert_dp_table_contract(x)
  if (ncol(x$table) != 2L || nrow(x$table) < 1L) {
    stop("x must be a strata-by-binary-outcome DP table", call. = FALSE)
  }
  if (!is.numeric(expected_rates) ||
      length(expected_rates) != nrow(x$table) ||
      anyNA(expected_rates) || any(!is.finite(expected_rates)) ||
      any(expected_rates < 0 | expected_rates > 1) ||
      !any(expected_rates > 0)) {
    stop(
      "expected_rates must be public probabilities in [0, 1] matching the strata",
      call. = FALSE)
  }
  if (!is.null(names(expected_rates))) {
    if (is.null(rownames(x$table)) || anyNA(names(expected_rates)) ||
        anyDuplicated(names(expected_rates)) ||
        !setequal(names(expected_rates), rownames(x$table))) {
      stop("Named expected_rates must match the DP stratum levels",
           call. = FALSE)
    }
    expected_rates <- expected_rates[rownames(x$table)]
  }
  expected_rates <- as.numeric(expected_rates)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }
  event <- .dsvert_dp_dimension_index(
    event, colnames(x$table), 2L, "event")
  nonevent <- setdiff(1:2, event)
  radius <- .dsvert_dp_table_simultaneous_radius(x, level)
  if (identical(level, 0.95) &&
      !.dsvert_dp_table_published_accuracy_matches(x, radius)) {
    stop("x does not carry a valid simultaneous DP accuracy certificate",
         call. = FALSE)
  }
  box <- .dsvert_dp_count_box(x$table, radius)
  region <- .dsvert_dp_indirect_box_region(
    box$lower[, event], box$upper[, event],
    box$lower[, nonevent], box$upper[, nonevent], expected_rates)
  observed <- sum(x$table[, event])
  expected <- sum(expected_rates * rowSums(x$table))
  estimate <- if (expected > 0) {
    observed / expected
  } else if (observed > 0) {
    Inf
  } else {
    NA_real_
  }
  status <- if (is.na(estimate)) {
    "dp_point_non_estimable"
  } else if (is.infinite(estimate)) {
    "boundary_infinite"
  } else if (estimate == 0) {
    "boundary_zero"
  } else {
    "ok"
  }
  result <- list(
    status = status,
    estimate = estimate,
    observed_events_dp = observed,
    expected_events_dp = expected,
    expected_rates = expected_rates,
    mechanism_region = region$interval,
    mechanism_region_includes_non_estimable =
      region$includes_non_estimable,
    mechanism_region_includes_infinite = region$includes_infinite,
    mechanism_region_has_estimable_values = region$estimable_values,
    mechanism_region_tightness = region$tightness,
    count_lower = box$lower,
    count_upper = box$upper,
    level = level,
    simultaneous_radius = radius,
    coverage_method = .dsvert_dp_table_coverage_method(x),
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    inferential_scope = paste(
      "Finite-dataset observed-to-expected ratio from one DP-noised table;",
      "no Garwood or sampling confidence interval, hypothesis test, causal",
      "interpretation, or population transportability claim is provided"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPIndirectStandardization", "list")
  result
}

#' @export
print.ds.vertDPEpi2x2 <- function(x, ...) {
  cat("dsVert DP 2x2 post-processing:", x$status, "\n")
  if (!is.null(x$point_estimates)) print(x$point_estimates, ...)
  cat("Simultaneous DP-noise regions (", format(100 * x$level), "%):\n",
      sep = "")
  print(do.call(rbind, x$mechanism_regions), ...)
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}

#' @export
print.ds.vertDPDiagnostic2x2 <- function(x, ...) {
  cat("dsVert DP diagnostic 2x2 post-processing:", x$status, "\n")
  cat("rows=disease status; columns=test result; positive levels: ",
      format(x$orientation$disease_positive$level), " / ",
      format(x$orientation$test_positive$level), "\n", sep = "")
  summary <- data.frame(
    estimate = unname(x$estimates),
    point_status = unname(x$point_status),
    lower = vapply(x$mechanism_regions, `[[`, numeric(1L), "lower"),
    upper = vapply(x$mechanism_regions, `[[`, numeric(1L), "upper"),
    row.names = names(x$estimates), check.names = FALSE)
  print(summary, ...)
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}

#' @export
print.ds.vertDPStandardization <- function(x, ...) {
  cat("dsVert DP direct standardisation:", x$status, "\n")
  cat("estimate:", format(x$estimate), " | DP-noise region: [",
      format(x$mechanism_region[["lower"]]), ", ",
      format(x$mechanism_region[["upper"]]), "]\n", sep = "")
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}

#' @export
print.ds.vertDPIndirectStandardization <- function(x, ...) {
  cat("dsVert DP indirect standardisation:", x$status, "\n")
  cat("observed/expected:", format(x$estimate), " | DP-noise region: [",
      format(x$mechanism_region[["lower"]]), ", ",
      format(x$mechanism_region[["upper"]]), "]\n", sep = "")
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}
