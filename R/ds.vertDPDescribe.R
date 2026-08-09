.dsvert_dp_validate_describe_release <- function(result, analysis_id,
                                                  server, status) {
  invalid <- function() {
    stop("Server '", server, "' returned an invalid DP describe release",
         call. = FALSE)
  }
  if (is.list(result) && length(result) == 1L &&
      is.list(result[[1L]]) && !"released" %in% names(result)) {
    if (is.null(names(result)) || !identical(names(result), server)) invalid()
    result <- result[[1L]]
  }
  if (!is.list(result) || !is.list(status) || !is.list(status$policy)) {
    invalid()
  }
  if (identical(result$released, FALSE)) {
    expected <- c(
      "released", "value", "reason", "epsilon", "delta", "adjacency",
      "composition_partitions")
    if (!.dsvert_dp_has_exact_names(result, expected) ||
        !is.null(result$value) ||
        !identical(result$reason,
                   "privacy_allocation_not_representable")) invalid()
    tryCatch(
      .dsvert_dp_validate_release_common(result, status, server),
      error = function(e) invalid())
    result$server <- server
    return(result)
  }

  expected <- c(
    "released", "analysis_id", "analysis_version", "variables",
    "variable_count", "lower_bounds", "upper_bounds", "grid_lengths",
    "grid_values", "histogram_semantics", "unit_collapse",
    "count_definition", "invalid_unit_rule", "statistics",
    "coordinate_count", "coordinate_layout", "numeric_grid_bits",
    "numeric_grid_scale", "quantization",
    "max_abs_normalized_quantization_per_unit", "allocation_names",
    "allocation_weights", "epsilon_per_variable_family",
    "epsilon_allocation_sum", "submechanism_count", "composition_rule",
    "calibration_rule", "family_sensitivity",
    "sum_family_l1_sensitivity_bound", "l2_sensitivity_bound",
    "noise_selection", "mechanism",
    "implementation", "sampler", "randomness", "postprocessing",
    "clipped_coordinates", "accuracy_95_abs_by_variable_family",
    "accuracy_simultaneous_95_abs_by_variable_family",
    "accuracy_simultaneous_confidence",
    "accuracy_simultaneous_method", "uncertainty_scope", "privacy_epoch",
    "noise_key_id", "sticky_noise", "epsilon", "delta", "adjacency",
    "composition_partitions")
  if (!identical(result$released, TRUE) ||
      !.dsvert_dp_has_exact_names(result, expected)) invalid()
  tryCatch(
    .dsvert_dp_validate_release_common(result, status, server),
    error = function(e) invalid())

  variables <- result$variables
  if (!identical(result$analysis_id, analysis_id) ||
      !.dsvert_dp_is_string(result$analysis_version) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$",
             result$analysis_version) ||
      !is.character(variables) || !length(variables) || anyNA(variables) ||
      any(!nzchar(variables)) || anyDuplicated(variables) ||
      !.dsvert_dp_is_integer(result$variable_count, 1,
                             .DSVERT_DP_MAX_COORDINATES) ||
      result$variable_count != length(variables)) invalid()
  variable_count <- length(variables)

  numeric_vector <- function(value, size) {
    is.numeric(value) && length(value) == size && !anyNA(value) &&
      all(is.finite(value))
  }
  if (!numeric_vector(result$lower_bounds, variable_count) ||
      !numeric_vector(result$upper_bounds, variable_count) ||
      any(result$lower_bounds >= result$upper_bounds) ||
      !is.numeric(result$grid_lengths) ||
      length(result$grid_lengths) != variable_count ||
      anyNA(result$grid_lengths) || any(!is.finite(result$grid_lengths)) ||
      any(result$grid_lengths < 1) ||
      any(result$grid_lengths != floor(result$grid_lengths)) ||
      any(result$grid_lengths > .DSVERT_DP_MAX_COORDINATES)) invalid()
  grid_lengths <- as.integer(result$grid_lengths)
  if (sum(grid_lengths) > .DSVERT_DP_MAX_COORDINATES ||
      !numeric_vector(result$grid_values, sum(grid_lengths))) invalid()

  policy_bounds <- status$policy$numeric_bounds
  grid_cursor <- 1L
  for (index in seq_along(variables)) {
    bounds <- policy_bounds[[variables[[index]]]]
    if (!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds) ||
        any(!is.finite(bounds)) ||
        !.dsvert_dp_num_equal(result$lower_bounds[[index]], bounds[[1L]]) ||
        !.dsvert_dp_num_equal(result$upper_bounds[[index]], bounds[[2L]])) {
      invalid()
    }
    length_i <- grid_lengths[[index]]
    grid <- result$grid_values[grid_cursor:(grid_cursor + length_i - 1L)]
    if (any(diff(grid) <= 0) || result$lower_bounds[[index]] >= grid[[1L]] ||
        !.dsvert_dp_num_equal(
          result$upper_bounds[[index]], grid[[length_i]])) invalid()
    grid_cursor <- grid_cursor + length_i
  }
  coordinate_count <- sum(grid_lengths + 4L)
  if (coordinate_count > .DSVERT_DP_MAX_COORDINATES ||
      !.dsvert_dp_is_integer(result$coordinate_count, 1,
                             .DSVERT_DP_MAX_COORDINATES) ||
      result$coordinate_count != coordinate_count ||
      !numeric_vector(result$statistics, coordinate_count) ||
      any(result$statistics < 0) ||
      any(result$statistics != floor(result$statistics)) ||
      any(result$statistics > 2^53 - 1)) invalid()

  grid_bits <- status$policy$numeric_grid_bits
  grid_scale <- 2^grid_bits
  allocation_names <- c("count", "sum", "sumsq", "histogram")
  if (!.dsvert_dp_is_integer(result$numeric_grid_bits, 8, 18) ||
      result$numeric_grid_bits != grid_bits ||
      !.dsvert_dp_is_integer(result$numeric_grid_scale, 2^8, 2^18) ||
      result$numeric_grid_scale != grid_scale ||
      !identical(result$allocation_names, allocation_names) ||
      !numeric_vector(result$allocation_weights, 4L) ||
      any(result$allocation_weights <= 0) ||
      abs(sum(result$allocation_weights) - 1) >
        4096 * .Machine$double.eps ||
      !numeric_vector(result$epsilon_per_variable_family, 4L)) invalid()
  expected_epsilon <-
    result$epsilon * result$allocation_weights / variable_count
  if (any(!vapply(seq_len(4L), function(index) {
        .dsvert_dp_num_equal(
          result$epsilon_per_variable_family[[index]],
          expected_epsilon[[index]], multiplier = 4096)
      }, logical(1L))) ||
      !.dsvert_dp_num_equal(
        result$epsilon_allocation_sum,
        variable_count * sum(result$epsilon_per_variable_family),
        multiplier = 4096) ||
      !.dsvert_dp_num_equal(
        result$epsilon_allocation_sum, result$epsilon,
        multiplier = 4096) ||
      !.dsvert_dp_is_integer(result$submechanism_count, 1L,
                             4L * .DSVERT_DP_MAX_COORDINATES)) invalid()

  family_sensitivity <- c(
    1, grid_scale, grid_scale,
    if (identical(
      status$policy$adjacency, "add_remove_patient")) 1 else 2)
  sum_family_l1_sensitivity_bound <-
    variable_count * sum(family_sensitivity)
  l2_sensitivity_bound <- sqrt(variable_count * (
    2 * grid_scale^2 + if (identical(
      status$policy$adjacency, "add_remove_patient")) 2 else 3))
  ratios <- family_sensitivity / result$epsilon_per_variable_family
  worst <- which.max(ratios)
  local_total_delta <- status$policy$local_total_delta
  if (is.null(local_total_delta)) local_total_delta <- 0
  scheduled_delta <- result$epsilon * local_total_delta /
    status$policy$local_total_epsilon
  selection <- .dsvert_dp_validate_noise_selection_certificate(
    result$noise_selection, coordinate_count, "simultaneous_95_abs",
    result$epsilon_per_variable_family[[worst]],
    family_sensitivity[[worst]], result$epsilon, scheduled_delta,
    l2_sensitivity_bound)
  gaussian_selected <- !is.null(selection) &&
    identical(selection$winner, "gaussian")
  laplace_marginal <- vapply(seq_len(4L), function(index) {
    .dsvert_dp_google_ci_radius(
      result$epsilon_per_variable_family[[index]],
      family_sensitivity[[index]])
  }, numeric(1L))
  laplace_simultaneous <- vapply(seq_len(4L), function(index) {
    .dsvert_dp_google_ci_radius(
      result$epsilon_per_variable_family[[index]],
      family_sensitivity[[index]], alpha = 0.05 / coordinate_count)
  }, numeric(1L))
  expected_marginal <- rep(
    if (gaussian_selected) selection$selected$marginal_95_abs else
      laplace_marginal,
    times = variable_count)
  expected_simultaneous <- rep(
    if (gaussian_selected) selection$selected$simultaneous_95_abs else
      laplace_simultaneous,
    times = variable_count)
  accuracy_length <- 4L * variable_count
  accuracy_ok <- function(value, expected_value) {
    numeric_vector(value, accuracy_length) &&
      all(value >= 0) && all(value == floor(value)) &&
      all(vapply(seq_len(accuracy_length), function(index) {
        .dsvert_dp_num_equal(value[[index]], expected_value[[index]],
                             multiplier = 4096)
      }, logical(1L)))
  }
  if (!identical(
        result$histogram_semantics,
        "(previous_endpoint,current_endpoint] plus fixed invalid bin") ||
      !identical(
        result$unit_collapse,
        "mean_of_finite_rows_after_public_bound_clipping") ||
      !identical(
        result$count_definition,
        "DP-noisy effective units with at least one finite bounded value") ||
      !identical(result$invalid_unit_rule,
                 "invalid_patient_ids_rejected_by_admission") ||
      !identical(
        result$coordinate_layout,
        paste("per_variable[count,qsum,qsumsq,",
              "histogram[grid_bins+invalid]]")) ||
      !identical(
        result$quantization,
        "round(z*scale) and round(z^2*scale) independently") ||
      !.dsvert_dp_is_number(
        result$max_abs_normalized_quantization_per_unit, 0,
        lower_open = TRUE) ||
      !.dsvert_dp_num_equal(
        result$max_abs_normalized_quantization_per_unit,
        0.5 / grid_scale) ||
      !identical(result$submechanism_count,
                 if (gaussian_selected) 1L else 4L * variable_count) ||
      !identical(
        result$composition_rule,
        if (gaussian_selected) {
          "single_joint_l2_gaussian_vector"
        } else {
          paste("sequential composition across variables and count/sum/sumsq/",
                "histogram families; each histogram is one L1 vector mechanism")
        }) ||
      !identical(
        result$calibration_rule,
        if (gaussian_selected) {
          "one analytic Gaussian calibrated to the joint L2 bound"
        } else {
          paste("each family is calibrated with its own epsilon and sensitivity;",
                "the summed sensitivity bound is descriptive, not a shared scale")
        }) ||
      !numeric_vector(result$family_sensitivity, 4L) ||
      !identical(as.numeric(result$family_sensitivity), family_sensitivity) ||
      !.dsvert_dp_is_integer(
        result$sum_family_l1_sensitivity_bound, 1) ||
      result$sum_family_l1_sensitivity_bound !=
        sum_family_l1_sensitivity_bound ||
      !.dsvert_dp_num_equal(
        result$l2_sensitivity_bound, l2_sensitivity_bound,
        multiplier = 2048) ||
      is.null(selection) ||
      !.dsvert_dp_selected_sampler_metadata_is_valid(
        result, selection, coordinate_count) ||
      !.dsvert_dp_num_equal(result$delta, selection$selected$delta) ||
      !identical(result$postprocessing,
                 "coordinatewise_nonnegative_integer") ||
      !accuracy_ok(
        result$accuracy_95_abs_by_variable_family,
        expected_marginal) ||
      !accuracy_ok(
        result$accuracy_simultaneous_95_abs_by_variable_family,
        expected_simultaneous) ||
      any(result$accuracy_simultaneous_95_abs_by_variable_family <
          result$accuracy_95_abs_by_variable_family) ||
      !.dsvert_dp_num_equal(
        result$accuracy_simultaneous_confidence, 0.95) ||
      !identical(result$accuracy_simultaneous_method, "union_bound") ||
      !identical(
        result$uncertainty_scope,
        "DP mechanism noise only; sampling uncertainty excluded")) invalid()
  result$server <- server
  result
}

.dsvert_dp_describe_quantiles <- function(counts, grid, lower_bound, probs,
                                           simultaneous_radius) {
  bin_count <- length(grid)
  if (!is.numeric(counts) || length(counts) != bin_count ||
      anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      !is.numeric(simultaneous_radius) ||
      length(simultaneous_radius) != 1L || is.na(simultaneous_radius) ||
      !is.finite(simultaneous_radius) || simultaneous_radius < 0) {
    stop("The DP histogram quantile inputs are invalid", call. = FALSE)
  }
  total <- sum(counts)
  lower_count <- pmax(0, counts - simultaneous_radius)
  upper_count <- counts + simultaneous_radius
  lower_prefix <- cumsum(lower_count)
  upper_prefix <- cumsum(upper_count)
  lower_cdf <- upper_cdf <- numeric(bin_count)
  for (index in seq_len(bin_count)) {
    upper_after <- sum(upper_count) - upper_prefix[[index]]
    lower_after <- sum(lower_count) - lower_prefix[[index]]
    lower_denominator <- lower_prefix[[index]] + upper_after
    upper_denominator <- upper_prefix[[index]] + lower_after
    lower_cdf[[index]] <- if (lower_denominator > 0) {
      lower_prefix[[index]] / lower_denominator
    } else 0
    upper_cdf[[index]] <- if (upper_denominator > 0) {
      upper_prefix[[index]] / upper_denominator
    } else 0
  }
  bin_lower <- c(lower_bound, head(grid, -1L))
  result <- lapply(probs, function(probability) {
    point_index <- if (total > 0) {
      which(cumsum(counts) >= probability * total)[[1L]]
    } else {
      NA_integer_
    }
    lower_candidates <- which(upper_cdf >= probability)
    upper_candidates <- which(lower_cdf >= probability)
    lower_index <- if (length(lower_candidates)) {
      lower_candidates[[1L]]
    } else 1L
    upper_index <- if (length(upper_candidates)) {
      upper_candidates[[1L]]
    } else bin_count
    if (upper_index < lower_index) upper_index <- lower_index
    data.frame(
      probability = probability,
      estimate = if (is.na(point_index)) NA_real_ else grid[[point_index]],
      dp_grid_lower = bin_lower[[lower_index]],
      dp_grid_upper = grid[[upper_index]],
      status = if (is.na(point_index)) {
        "dp_histogram_empty_after_postprocessing"
      } else {
        "ok"
      },
      stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}

.dsvert_dp_describe_moment_region <- function(
    n_dp, qsum_dp, qsumsq_dp, count_radius, sum_radius, sumsq_radius,
    grid_scale, lower_bound, upper_bound) {
  values <- c(
    n_dp = n_dp, qsum_dp = qsum_dp, qsumsq_dp = qsumsq_dp,
    count_radius = count_radius, sum_radius = sum_radius,
    sumsq_radius = sumsq_radius, grid_scale = grid_scale,
    lower_bound = lower_bound, upper_bound = upper_bound)
  if (anyNA(values) || any(!is.finite(values)) ||
      any(values[c("n_dp", "qsum_dp", "qsumsq_dp", "count_radius",
                   "sum_radius", "sumsq_radius")] < 0) ||
      grid_scale <= 0 || lower_bound >= upper_bound) {
    stop("The DP moment-region inputs are invalid", call. = FALSE)
  }

  n_lower <- max(0, ceiling(n_dp - count_radius))
  n_upper <- max(0, floor(n_dp + count_radius))
  width <- upper_bound - lower_bound
  full <- function(status) list(
    status = status,
    effective_count = c(lower = n_lower, upper = n_upper),
    mean = c(lower = lower_bound, upper = upper_bound),
    variance = c(lower = 0, upper = width^2 / 4))
  if (n_upper < 1) {
    return(full("dp_region_has_no_positive_effective_count"))
  }

  positive_n_lower <- max(1, n_lower)
  qsum_lower <- max(0, qsum_dp - sum_radius)
  qsum_upper <- qsum_dp + sum_radius
  qsumsq_lower <- max(0, qsumsq_dp - sumsq_radius)
  qsumsq_upper <- qsumsq_dp + sumsq_radius
  quantization <- 0.5 / grid_scale
  guard <- 64 * .Machine$double.eps

  mean_lower <- max(
    0, qsum_lower / grid_scale / n_upper - quantization - guard)
  mean_upper <- min(
    1, qsum_upper / grid_scale / positive_n_lower + quantization + guard)
  second_lower <- max(
    0, qsumsq_lower / grid_scale / n_upper - quantization - guard)
  second_upper <- min(
    1, qsumsq_upper / grid_scale / positive_n_lower + quantization + guard)
  if (mean_lower > mean_upper || second_lower > second_upper) {
    return(full("dp_simultaneous_noise_box_moment_geometry_infeasible"))
  }

  variance_lower <- max(0, second_lower - mean_upper^2 - guard)
  bounded_variance_upper <- if (mean_lower <= 0.5 && mean_upper >= 0.5) {
    0.25
  } else {
    max(mean_lower * (1 - mean_lower),
        mean_upper * (1 - mean_upper))
  }
  variance_upper <- min(
    bounded_variance_upper,
    max(0, second_upper - mean_lower^2 + guard))
  if (variance_lower > variance_upper) {
    return(full("dp_simultaneous_noise_box_moment_geometry_infeasible"))
  }

  value_guard <- 64 * .Machine$double.eps *
    max(1, abs(lower_bound), abs(upper_bound), width^2)
  list(
    status = if (n_lower > 0) {
      "ok"
    } else {
      "dp_region_includes_zero_effective_count"
    },
    effective_count = c(lower = n_lower, upper = n_upper),
    mean = c(
      lower = max(lower_bound,
                  lower_bound + width * mean_lower - value_guard),
      upper = min(upper_bound,
                  lower_bound + width * mean_upper + value_guard)),
    variance = c(
      lower = max(0, width^2 * variance_lower - value_guard),
      upper = min(width^2 / 4,
                  width^2 * variance_upper + value_guard)))
}

.dsvert_dp_describe_postprocess <- function(result, probs) {
  probs <- sort(unique(as.numeric(probs)), method = "radix")
  if (!isTRUE(result$released)) {
    result$descriptives <- NULL
    result$histograms <- NULL
    result$quantiles <- NULL
    return(result)
  }
  accuracy <- matrix(
    result$accuracy_95_abs_by_variable_family,
    ncol = 4L, byrow = TRUE,
    dimnames = list(result$variables, result$allocation_names))
  simultaneous_accuracy <- matrix(
    result$accuracy_simultaneous_95_abs_by_variable_family,
    ncol = 4L, byrow = TRUE,
    dimnames = list(result$variables, result$allocation_names))
  cursor <- grid_cursor <- 1L
  descriptions <- vector("list", result$variable_count)
  histograms <- setNames(vector("list", result$variable_count),
                         result$variables)
  quantiles <- vector("list", result$variable_count)
  for (index in seq_len(result$variable_count)) {
    variable <- result$variables[[index]]
    grid_length <- as.integer(result$grid_lengths[[index]])
    grid <- result$grid_values[
      grid_cursor:(grid_cursor + grid_length - 1L)]
    n_dp <- result$statistics[[cursor]]
    qsum_dp <- result$statistics[[cursor + 1L]]
    qsumsq_dp <- result$statistics[[cursor + 2L]]
    histogram <- result$statistics[
      (cursor + 3L):(cursor + 3L + grid_length)]
    cursor <- cursor + grid_length + 4L
    grid_cursor <- grid_cursor + grid_length
    valid_histogram <- histogram[seq_len(grid_length)]
    invalid_dp <- histogram[[grid_length + 1L]]
    count_lower <- max(0, n_dp - simultaneous_accuracy[index, "count"])
    point_estimable <- n_dp > 0
    count_certified <- count_lower > 0
    mean_value <- variance_value <- sd_value <- NA_real_
    if (point_estimable) {
      normalized_sum <- qsum_dp / result$numeric_grid_scale
      normalized_sumsq <- qsumsq_dp / result$numeric_grid_scale
      sum_projected <- min(n_dp, max(0, normalized_sum))
      sumsq_projected <- min(
        sum_projected,
        max(sum_projected^2 / n_dp, normalized_sumsq))
      mean_normalized <- sum_projected / n_dp
      variance_normalized <- max(
        0, sumsq_projected / n_dp - mean_normalized^2)
      width <- result$upper_bounds[[index]] - result$lower_bounds[[index]]
      mean_value <- result$lower_bounds[[index]] + width * mean_normalized
      variance_value <- width^2 * variance_normalized
      sd_value <- sqrt(variance_value)
    }
    moment_region <- .dsvert_dp_describe_moment_region(
      n_dp = n_dp, qsum_dp = qsum_dp, qsumsq_dp = qsumsq_dp,
      count_radius = simultaneous_accuracy[index, "count"],
      sum_radius = simultaneous_accuracy[index, "sum"],
      sumsq_radius = simultaneous_accuracy[index, "sumsq"],
      grid_scale = result$numeric_grid_scale,
      lower_bound = result$lower_bounds[[index]],
      upper_bound = result$upper_bounds[[index]])
    descriptions[[index]] <- data.frame(
      variable = variable,
      status = if (count_certified) {
        "ok"
      } else if (point_estimable) {
        "dp_point_available_count_not_certified_positive"
      } else {
        "dp_effective_count_zero_after_postprocessing"
      },
      n_dp = n_dp,
      n_simultaneous_95_lower = count_lower,
      mean = mean_value,
      variance = variance_value,
      sd = sd_value,
      mean_mechanism_grid_lower_95 =
        moment_region$mean[["lower"]],
      mean_mechanism_grid_upper_95 =
        moment_region$mean[["upper"]],
      variance_mechanism_grid_lower_95 =
        moment_region$variance[["lower"]],
      variance_mechanism_grid_upper_95 =
        moment_region$variance[["upper"]],
      moment_region_status = moment_region$status,
      lower_bound = result$lower_bounds[[index]],
      upper_bound = result$upper_bounds[[index]],
      count_noise_radius_95 = accuracy[index, "count"],
      count_noise_radius_simultaneous_95 =
        simultaneous_accuracy[index, "count"],
      normalized_sum_noise_radius_95 =
        accuracy[index, "sum"] / result$numeric_grid_scale,
      normalized_sumsq_noise_radius_95 =
        accuracy[index, "sumsq"] / result$numeric_grid_scale,
      max_abs_normalized_quantization_per_unit =
        result$max_abs_normalized_quantization_per_unit,
      invalid_dp = invalid_dp,
      stringsAsFactors = FALSE)
    histograms[[index]] <- list(
      grid = grid,
      counts = valid_histogram,
      invalid_dp = invalid_dp,
      cell_noise_radius_95 = accuracy[index, "histogram"],
      cell_noise_radius_simultaneous_95 =
        simultaneous_accuracy[index, "histogram"])
    quantile <- .dsvert_dp_describe_quantiles(
      valid_histogram, grid, result$lower_bounds[[index]], probs,
      simultaneous_accuracy[index, "histogram"])
    quantile$variable <- variable
    quantile <- quantile[c(
      "variable", "probability", "estimate", "dp_grid_lower",
      "dp_grid_upper", "status")]
    quantiles[[index]] <- quantile
  }
  result$status <- if (isTRUE(result$clipped_coordinates > 0)) {
    "dp_sampler_coordinates_clipped"
  } else if (is.na(result$clipped_coordinates)) {
    "fixed_public_clamp_applied_preclamp_state_not_released"
  } else {
    "ok"
  }
  result$descriptives <- do.call(rbind, descriptions)
  rownames(result$descriptives) <- NULL
  result$histograms <- histograms
  result$quantiles <- do.call(rbind, quantiles)
  rownames(result$quantiles) <- NULL
  result$quantile_band_confidence <- 0.95
  result$quantile_band_scope <- paste(
    "simultaneous DP histogram-noise box plus public grid interval;",
    "sampling uncertainty excluded")
  result$moment_region_confidence <- 0.95
  result$moment_region_method <- paste(
    "simultaneous coordinate-box propagation with bounded-moment geometry",
    "and deterministic per-unit quantization error")
  result$moment_region_scope <- paste(
    "simultaneous DP mechanism-noise box plus deterministic quantization",
    "bounds, conditional on a positive effective count;",
    "sampling uncertainty excluded")
  result$statistical_inference <- paste(
    "DP-noised point estimates; no sampling confidence intervals,",
    "standard errors, extrema, or p-values are provided")
  result
}

.dsvert_dp_describe_vector_result <- function(
    capsule, data_name, analysis_id, server = NULL) {
  capsule <- .dsvert_dp_vector_context(capsule)
  release <- capsule$release
  manifest <- release$manifest
  families <- manifest$workload$families
  artifact <- families$describe_artifacts[[analysis_id]]
  variables <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$variables, "describe variable list", sorted = TRUE),
    error = function(error) character())
  if (!is.list(artifact) || !identical(artifact$dataset, data_name) ||
      !length(variables) ||
      !is.character(artifact$owner_peer) ||
      length(artifact$owner_peer) != 1L ||
      !artifact$owner_peer %in% names(capsule$status)) {
    stop("The signed capsule does not contain the requested describe artifact",
         call. = FALSE)
  }
  owner <- artifact$owner_peer
  if (!is.null(server)) {
    if (!is.character(server) || length(server) != 1L || is.na(server) ||
        !identical(server, owner)) {
      stop("server does not own the signed describe artifact", call. = FALSE)
    }
  }
  scale <- as.numeric(manifest$workload$release_lattice$output_lattice_scale)
  bits <- as.integer(manifest$workload$release_lattice$output_lattice_bits)
  capacity <- as.numeric(capsule$status[[owner]]$policy$unit_capacity)
  marginal <- .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = 1L, maximum_error = capacity)

  moment_blocks <- histogram_blocks <- vector("list", length(variables))
  total_coordinates <- 0L
  for (index in seq_along(variables)) {
    variable <- variables[[index]]
    moment_blocks[[index]] <- .dsvert_dp_capsule_single_block(
      capsule$layout, "numeric_moments", dataset = data_name,
      owner_peer = owner,
      predicate = function(block) {
        identical(block$descriptor$column, variable)
      }, description = paste0("numeric moment block for '", variable, "'"))
    reference <- artifact$histogram_references[[index]]
    primitive_id <- if (is.list(reference)) reference$primitive_id else NULL
    histogram_blocks[[index]] <- .dsvert_dp_capsule_single_block(
      capsule$layout, "fixed_numeric_histograms", dataset = data_name,
      owner_peer = owner,
      predicate = function(block) identical(block$key, primitive_id),
      description = paste0("fixed histogram block for '", variable, "'"))
    descriptor <- histogram_blocks[[index]]$descriptor
    descriptor_grid <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
      descriptor$grid, "histogram grid"), error = function(error) numeric())
    artifact_grid <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
      artifact$histogram_grids[[index]], "describe histogram grid"),
      error = function(error) numeric())
    if (!identical(descriptor$column, variable) ||
        !identical(descriptor_grid, artifact_grid) ||
        histogram_blocks[[index]]$length != length(descriptor_grid) + 1L) {
      stop("A signed describe histogram reference is inconsistent",
           call. = FALSE)
    }
    total_coordinates <- total_coordinates + 3L +
      histogram_blocks[[index]]$length
  }
  simultaneous <- .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = total_coordinates,
    maximum_error = capacity)

  statistics <- numeric()
  grids <- numeric()
  lower <- upper <- numeric(length(variables))
  grid_lengths <- integer(length(variables))
  marginal_accuracy <- simultaneous_accuracy <- numeric()
  for (index in seq_along(variables)) {
    moments <- .dsvert_dp_capsule_vector_values(
      release, moment_blocks[[index]])
    histogram <- .dsvert_dp_capsule_vector_values(
      release, histogram_blocks[[index]])
    # The released common lattice is natural: count and histogram coordinates
    # are counts, while numeric sums are normalized values. The established
    # postprocessor consumes exact integer qsum/qsumsq coordinates.
    qmoments <- c(moments[[1L]], moments[[2L]] * scale,
                  moments[[3L]] * scale)
    if (any(abs(qmoments - round(qmoments)) >
            64 * .Machine$double.eps * pmax(1, abs(qmoments)))) {
      stop("A released numeric moment is not on its signed lattice",
           call. = FALSE)
    }
    statistics <- c(statistics, round(qmoments), histogram)
    descriptor <- moment_blocks[[index]]$descriptor
    lower[[index]] <- descriptor$lower
    upper[[index]] <- descriptor$upper
    grid <- .dsvert_dp_capsule_manifest_numbers(
      histogram_blocks[[index]]$descriptor$grid, "histogram grid")
    grids <- c(grids, grid)
    grid_lengths[[index]] <- length(grid)
    marginal_accuracy <- c(
      marginal_accuracy, marginal$radius,
      marginal$radius * scale, marginal$radius * scale, marginal$radius)
    simultaneous_accuracy <- c(
      simultaneous_accuracy, simultaneous$radius,
      simultaneous$radius * scale, simultaneous$radius * scale,
      simultaneous$radius)
  }
  noise_root <- capsule$status[[owner]]$noise_root
  profile <- .dsvert_vector_profile(
    manifest$workload$capsule_mechanism,
    manifest$workload$mechanism_selection)
  list(
    released = TRUE, analysis_id = analysis_id,
    analysis_version = artifact$version, variables = variables,
    variable_count = length(variables), lower_bounds = lower,
    upper_bounds = upper, grid_lengths = grid_lengths, grid_values = grids,
    histogram_semantics =
      "(previous_endpoint,current_endpoint] plus fixed invalid bin",
    unit_collapse =
      "mean_of_finite_rows_after_public_bound_clipping",
    count_definition =
      "DP-noisy effective units with at least one finite bounded value",
    invalid_unit_rule = "invalid_patient_ids_rejected_by_admission",
    statistics = statistics, coordinate_count = total_coordinates,
    coordinate_layout = paste(
      "referenced capsule blocks per variable[count,qsum,qsumsq,",
      "histogram[grid_bins+invalid]]"),
    numeric_grid_bits = bits, numeric_grid_scale = scale,
    quantization = "round(z*scale) and round(z^2*scale) independently",
    max_abs_normalized_quantization_per_unit = 0.5 / scale,
    allocation_names = c("count", "sum", "sumsq", "histogram"),
    allocation_weights = .dsvert_dp_capsule_manifest_numbers(
      artifact$allocation_weights, "describe allocation weights"),
    mechanism = release$mechanism,
    implementation = "two pinned peers; Ring128 exact signed finalizer",
    sampler = profile$sampler,
    randomness = "independent HMAC/HKDF/ChaCha20 streams at both peers",
    postprocessing = "fixed public per-coordinate clamp",
    clipped_coordinates = NA_integer_, clipping_observable = FALSE,
    accuracy_95_abs_by_variable_family = marginal_accuracy,
    accuracy_simultaneous_95_abs_by_variable_family = simultaneous_accuracy,
    accuracy_simultaneous_confidence = simultaneous$confidence,
    accuracy_simultaneous_method = simultaneous$method,
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    privacy_epoch = noise_root$privacy_epoch,
    noise_key_id = noise_root$key_id,
    sticky_noise = "one immutable capsule vector; unlimited replay",
    epsilon = release$epsilon, delta = release$delta,
    implementation_delta = release$implementation_delta,
    adjacency = capsule$status[[owner]]$policy$adjacency,
    composition_partitions = 1L, capsule_id = release$capsule_id,
    final_vector_root = release$final_vector_root,
    coordinate_order_sha256 = release$coordinate_order_sha256,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    security_claim = .dsvert_dp_capsule_security_claim(),
    server = owner)
}

#' Differentially private fixed-grid descriptive statistics
#'
#' Makes one sticky release for a custodian-owned variable set and returns
#' noisy effective counts, bounded mean/variance/SD, fixed histograms, and
#' histogram-CDF quantiles. It propagates the simultaneous coordinate-noise
#' box and deterministic quantisation bound into conservative mean/variance
#' regions. Observed extrema and adaptive bins are never used.
#'
#' @param data_name Name of the registered protected data frame.
#' @param analysis_id Custodian-owned describe specification id.
#' @param probs Public quantile probabilities strictly between zero and one.
#' @param server Optional datasource name. If omitted, the lexicographically
#'   first datasource is selected deterministically.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPDescribe` object. Mechanism/grid regions exclude
#'   sampling error.
#' @export
ds.vertDPDescribe <- function(data_name, analysis_id,
                              probs = c(0.25, 0.5, 0.75),
                              server = NULL, datasources = NULL) {
  .dsvert_dp_describe_impl(
    data_name, analysis_id, probs, server, datasources,
    DSI::datashield.aggregate)
}

.dsvert_dp_describe_impl <- function(data_name, analysis_id, probs,
                                     server = NULL, datasources = NULL,
                                     .aggregate) {
  for (value in list(data_name, analysis_id)) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop("data_name and analysis_id must be non-empty strings",
           call. = FALSE)
    }
  }
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("probs must contain finite probabilities strictly inside (0,1)",
         call. = FALSE)
  }
  probs <- sort(unique(as.numeric(probs)), method = "radix")
  datasources <- .dsvert_dp_datasources(datasources)
  capsule <- .dsvert_dp_capsule_vector_run(
    datasources, .aggregate = .aggregate)
  result <- .dsvert_dp_describe_vector_result(
    capsule, data_name, analysis_id, server)
  result <- .dsvert_dp_describe_postprocess(result, probs)
  class(result) <- c("ds.vertDPDescribe", "list")
  result
}

#' @export
print.ds.vertDPDescribe <- function(x, ...) {
  if (!isTRUE(x$released)) {
    cat("dsVert DP describe: suppressed (", x$reason, ")\n", sep = "")
  } else {
    cat("dsVert DP describe:", x$analysis_id, "[", x$server, "]\n")
    print(x$descriptives, row.names = FALSE)
    cat("epsilon:", format(x$epsilon), "| variables:",
        x$variable_count, "| status:", x$status, "\n")
    cat("Quantile bands cover DP histogram noise and grid intervals,",
        "not sampling uncertainty.\n")
  }
  invisible(x)
}
