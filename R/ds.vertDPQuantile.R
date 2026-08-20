.dsvert_dp_describe_provenance_contract <- function(x) {
  invalid <- function() {
    stop(
      "x must be an intact released ds.vertDPDescribe capsule object",
      call. = FALSE)
  }
  legacy_fields <- c(
    "capsule_id", "privacy_epoch", "noise_key_id", "history_gate",
    "request_limit", "operation_limit")
  synopsis_fields <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "release_provenance",
    "privacy", "composition_rule")
  if (!inherits(x, "ds.vertDPDescribe") || !is.list(x) ||
      !identical(x$released, TRUE) || is.null(names(x)) ||
      anyNA(names(x)) || anyDuplicated(names(x)) ||
      !.dsvert_vector_hex(x$final_vector_root) ||
      !.dsvert_vector_hex(x$coordinate_order_sha256) ||
      !.dsvert_dp_is_number(
        x$epsilon, 0, .DSVERT_DP_MAXIMUM_EPSILON, lower_open = TRUE) ||
      !.dsvert_dp_is_number(x$delta, 0, 1) || x$delta >= 1) {
    invalid()
  }

  legacy <- all(legacy_fields %in% names(x)) &&
    !any(synopsis_fields %in% names(x))
  synopsis <- all(synopsis_fields %in% names(x)) &&
    !any(legacy_fields %in% names(x))
  if (identical(legacy, synopsis)) invalid()

  if (isTRUE(legacy)) {
    if (!identical(x$history_gate, TRUE) ||
        !identical(x$request_limit, FALSE) ||
        !identical(x$operation_limit, TRUE) ||
        !.dsvert_vector_whole(x$privacy_epoch, 1, 2^53 - 1) ||
        !.dsvert_vector_string(x$noise_key_id) ||
        !.dsvert_vector_hex(x$capsule_id)) {
      invalid()
    }
    return("legacy")
  }

  anchors <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  provenance <- x$release_provenance
  if (!is.list(provenance) || is.null(names(provenance)) ||
      anyNA(names(provenance)) || anyDuplicated(names(provenance)) ||
      !all(c("version", anchors) %in% names(provenance)) ||
      any(legacy_fields %in% names(provenance)) ||
      !identical(
        provenance$version,
        "dsvert-stateless-synopsis-public-provenance-v1") ||
      !all(vapply(anchors, function(field) {
        .dsvert_vector_hex(x[[field]]) &&
          identical(x[[field]], provenance[[field]])
      }, logical(1L)))) {
    invalid()
  }
  privacy_fields <- c(
    "version", "per_artifact_epsilon", "per_artifact_delta",
    "sticky_noise", "public_openings", "distinct_artifacts_compose",
    "finite_global_composition_claim")
  privacy <- x$privacy
  if (!.dsvert_dp_has_exact_names(privacy, privacy_fields) ||
      !identical(privacy$version, "dsvert-per-synopsis-dp-v1") ||
      !.dsvert_dp_is_number(
        privacy$per_artifact_epsilon, 0, .DSVERT_DP_MAXIMUM_EPSILON,
        lower_open = TRUE) ||
      !.dsvert_dp_is_number(privacy$per_artifact_delta, 0, 1) ||
      privacy$per_artifact_delta >= 1 ||
      !identical(as.numeric(privacy$per_artifact_epsilon),
                 as.numeric(x$epsilon)) ||
      !identical(as.numeric(privacy$per_artifact_delta),
                 as.numeric(x$delta)) ||
      !identical(privacy$sticky_noise, TRUE) ||
      !.dsvert_vector_whole(privacy$public_openings, 1, 1) ||
      !identical(privacy$distinct_artifacts_compose, TRUE) ||
      !identical(privacy$finite_global_composition_claim, FALSE) ||
      !identical(
        x$composition_rule,
        "one_sticky_release_per_canonical_signed_artifact")) {
    invalid()
  }
  "synopsis"
}

.dsvert_dp_describe_postprocess_contract <- function(x) {
  invalid <- function() {
    stop(
      "x must be an intact released ds.vertDPDescribe capsule object",
      call. = FALSE)
  }
  numeric_vector <- function(value, size, non_negative = FALSE) {
    is.numeric(value) && length(value) == size && !anyNA(value) &&
      all(is.finite(value)) &&
      (!isTRUE(non_negative) || all(value >= 0))
  }
  scalar_string <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value) && nchar(value, type = "bytes") <= 1024L
  }

  .dsvert_dp_describe_provenance_contract(x)
  if (!inherits(x, "ds.vertDPDescribe") || !is.list(x) ||
      !identical(x$released, TRUE) ||
      !identical(
        x$status,
        "fixed_public_clamp_applied_preclamp_state_not_released") ||
      !scalar_string(x$analysis_id) || !scalar_string(x$analysis_version) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", x$analysis_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$",
             x$analysis_version) ||
      !scalar_string(x$server) ||
      !is.character(x$variables) || !length(x$variables) ||
      anyNA(x$variables) || any(!nzchar(x$variables)) ||
      any(nchar(x$variables, type = "bytes") > 1024L) ||
      anyDuplicated(x$variables) ||
      !.dsvert_dp_is_integer(
        x$variable_count, 1L, .DSVERT_DP_MAX_COORDINATES) ||
      x$variable_count != length(x$variables)) {
    invalid()
  }
  variable_count <- length(x$variables)

  if (!numeric_vector(x$lower_bounds, variable_count) ||
      !numeric_vector(x$upper_bounds, variable_count) ||
      any(x$lower_bounds >= x$upper_bounds) ||
      !numeric_vector(x$grid_lengths, variable_count, TRUE) ||
      any(x$grid_lengths < 1) ||
      any(x$grid_lengths != floor(x$grid_lengths)) ||
      any(x$grid_lengths > .DSVERT_DP_MAX_COORDINATES)) {
    invalid()
  }
  grid_lengths <- as.integer(x$grid_lengths)
  grid_total <- sum(grid_lengths)
  coordinate_count <- sum(grid_lengths + 4L)
  if (grid_total > .DSVERT_DP_MAX_COORDINATES ||
      coordinate_count > .DSVERT_DP_MAX_COORDINATES ||
      !numeric_vector(x$grid_values, grid_total) ||
      !.dsvert_dp_is_integer(
        x$coordinate_count, 1L, .DSVERT_DP_MAX_COORDINATES) ||
      x$coordinate_count != coordinate_count ||
      !numeric_vector(x$statistics, coordinate_count, TRUE) ||
      any(x$statistics > 2^53 - 1)) {
    invalid()
  }

  if (!.dsvert_dp_is_integer(x$numeric_grid_bits, 8L, 18L) ||
      !.dsvert_dp_is_number(x$numeric_grid_scale, 0, Inf,
                            lower_open = TRUE) ||
      !identical(as.numeric(x$numeric_grid_scale),
                 2^as.numeric(x$numeric_grid_bits)) ||
      !identical(x$allocation_names,
                 c("count", "sum", "sumsq", "histogram")) ||
      !numeric_vector(x$allocation_weights, 4L, TRUE) ||
      any(x$allocation_weights <= 0) ||
      abs(sum(x$allocation_weights) - 1) >
        4096 * .Machine$double.eps ||
      !.dsvert_dp_num_equal(
        x$max_abs_normalized_quantization_per_unit,
        0.5 / x$numeric_grid_scale, multiplier = 4096)) {
    invalid()
  }

  accuracy_length <- 4L * variable_count
  if (!numeric_vector(
        x$accuracy_95_abs_by_variable_family,
        accuracy_length, TRUE) ||
      !numeric_vector(
        x$accuracy_simultaneous_95_abs_by_variable_family,
        accuracy_length, TRUE) ||
      any(x$accuracy_simultaneous_95_abs_by_variable_family <
          x$accuracy_95_abs_by_variable_family) ||
      !.dsvert_dp_num_equal(x$accuracy_simultaneous_confidence, 0.95) ||
      !scalar_string(x$accuracy_simultaneous_method) ||
      !identical(
        x$uncertainty_scope,
        "DP mechanism noise only; sampling uncertainty excluded") ||
      !.dsvert_dp_num_equal(x$quantile_band_confidence, 0.95) ||
      !scalar_string(x$quantile_band_scope) ||
      !grepl("sampling uncertainty excluded", x$quantile_band_scope,
             fixed = TRUE)) {
    invalid()
  }

  profile <- if (identical(
      x$mechanism, .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM)) {
    list(sampler = .DSVERT_CLIENT_VECTOR_SAMPLER)
  } else if (identical(
      x$mechanism, .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)) {
    list(sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER)
  } else {
    NULL
  }
  if (is.null(profile) || !identical(x$sampler, profile$sampler) ||
      !identical(
        x$implementation,
        "two pinned peers; Ring128 exact signed finalizer") ||
      !identical(
        x$randomness,
        "independent HMAC/HKDF/ChaCha20 streams at both peers") ||
      !identical(x$postprocessing, "fixed public per-coordinate clamp") ||
      !identical(x$clipped_coordinates, NA_integer_) ||
      !identical(x$clipping_observable, FALSE) ||
      !identical(
        x$sticky_noise,
        "one immutable capsule vector; unlimited replay") ||
      !.dsvert_dp_is_number(
        x$epsilon, 0, .DSVERT_DP_MAXIMUM_EPSILON, lower_open = TRUE) ||
      !.dsvert_dp_is_number(x$delta, 0, 1) || x$delta >= 1 ||
      !scalar_string(x$implementation_delta) ||
      !scalar_string(x$adjacency) ||
      !x$adjacency %in% c(
        "add_remove_patient", "replace_one_fixed_cohort") ||
      !.dsvert_vector_hex(x$final_vector_root) ||
      !.dsvert_vector_hex(x$coordinate_order_sha256) ||
      !identical(
        x$histogram_semantics,
        "(previous_endpoint,current_endpoint] plus fixed invalid bin") ||
      !identical(
        x$unit_collapse,
        "mean_of_finite_rows_after_public_bound_clipping") ||
      !identical(
        x$count_definition,
        "DP-noisy effective units with at least one finite bounded value") ||
      !identical(
        x$invalid_unit_rule,
        "invalid_patient_ids_rejected_by_admission") ||
      !identical(
        x$coordinate_layout,
        paste(
          "referenced capsule blocks per variable[count,qsum,qsumsq,",
          "histogram[grid_bins+invalid]]")) ||
      !identical(
        x$quantization,
        "round(z*scale) and round(z^2*scale) independently")) {
    invalid()
  }

  required_descriptives <- c(
    "variable", "status", "n_dp", "lower_bound", "upper_bound",
    "invalid_dp")
  if (!is.data.frame(x$descriptives) ||
      !all(required_descriptives %in% names(x$descriptives)) ||
      nrow(x$descriptives) != variable_count ||
      !identical(as.character(x$descriptives$variable), x$variables) ||
      !is.list(x$histograms) ||
      !identical(names(x$histograms), x$variables)) {
    invalid()
  }

  cursor <- grid_cursor <- 1L
  for (index in seq_len(variable_count)) {
    grid_length <- grid_lengths[[index]]
    grid <- x$grid_values[
      grid_cursor:(grid_cursor + grid_length - 1L)]
    histogram <- x$histograms[[index]]
    block <- x$statistics[cursor:(cursor + grid_length + 3L)]
    expected_counts <- block[4L:(grid_length + 3L)]
    expected_invalid <- block[[grid_length + 4L]]
    lattice_values <- c(expected_counts, expected_invalid) *
      x$numeric_grid_scale
    lattice_ok <- all(
      abs(lattice_values - round(lattice_values)) <=
        64 * .Machine$double.eps * pmax(1, abs(lattice_values)))
    expected_description_status <- if (
        block[[1L]] >
          x$accuracy_simultaneous_95_abs_by_variable_family[
            4L * (index - 1L) + 1L]) {
      "ok"
    } else if (block[[1L]] > 0) {
      "dp_point_available_count_not_certified_positive"
    } else {
      "dp_effective_count_zero_after_postprocessing"
    }
    accuracy_index <- 4L * (index - 1L) + 4L
    if (any(diff(grid) <= 0) ||
        x$lower_bounds[[index]] >= grid[[1L]] ||
        !.dsvert_dp_num_equal(
          x$upper_bounds[[index]], grid[[grid_length]]) ||
        !is.list(histogram) ||
        !identical(
          names(histogram),
          c("grid", "counts", "invalid_dp", "cell_noise_radius_95",
            "cell_noise_radius_simultaneous_95")) ||
        !numeric_vector(histogram$grid, grid_length) ||
        !identical(as.numeric(histogram$grid), as.numeric(grid)) ||
        !numeric_vector(histogram$counts, grid_length, TRUE) ||
        !identical(as.numeric(histogram$counts),
                   as.numeric(expected_counts)) ||
        !lattice_ok ||
        !.dsvert_dp_num_equal(
          block[[1L]] * x$numeric_grid_scale,
          round(block[[1L]] * x$numeric_grid_scale), multiplier = 64) ||
        any(block[2:3] != floor(block[2:3])) ||
        !.dsvert_dp_is_number(histogram$invalid_dp, 0, 2^53 - 1) ||
        !identical(as.numeric(histogram$invalid_dp),
                   as.numeric(expected_invalid)) ||
        !.dsvert_dp_num_equal(
          histogram$cell_noise_radius_95,
          x$accuracy_95_abs_by_variable_family[[accuracy_index]]) ||
        !.dsvert_dp_num_equal(
          histogram$cell_noise_radius_simultaneous_95,
          x$accuracy_simultaneous_95_abs_by_variable_family[
            accuracy_index]) ||
        !.dsvert_dp_num_equal(
          x$descriptives$n_dp[[index]], block[[1L]]) ||
        !identical(
          as.character(x$descriptives$status[[index]]),
          expected_description_status) ||
        !.dsvert_dp_num_equal(
          x$descriptives$lower_bound[[index]], x$lower_bounds[[index]]) ||
        !.dsvert_dp_num_equal(
          x$descriptives$upper_bound[[index]], x$upper_bounds[[index]]) ||
        !.dsvert_dp_num_equal(
          x$descriptives$invalid_dp[[index]], expected_invalid)) {
      invalid()
    }
    cursor <- cursor + grid_length + 4L
    grid_cursor <- grid_cursor + grid_length
  }
  x
}

.dsvert_dp_quantile_band_indices <- function(counts, probability, radius) {
  bin_count <- length(counts)
  lower <- pmax(0, counts - radius)
  upper <- counts + radius
  possible <- which(upper > 0)
  forced <- which(lower > 0)
  if (!length(possible)) return(c(lower = 1L, upper = bin_count))

  if (probability == 0) {
    upper_index <- if (length(forced)) min(forced) else max(possible)
    return(c(lower = min(possible), upper = upper_index))
  }
  if (probability == 1) {
    lower_index <- if (length(forced)) max(forced) else min(possible)
    return(c(lower = lower_index, upper = max(possible)))
  }

  lower_prefix <- cumsum(lower)
  upper_prefix <- cumsum(upper)
  lower_cdf <- upper_cdf <- numeric(bin_count)
  for (index in seq_len(bin_count)) {
    lower_after <- sum(lower) - lower_prefix[[index]]
    upper_after <- sum(upper) - upper_prefix[[index]]
    lower_denominator <- lower_prefix[[index]] + upper_after
    upper_denominator <- upper_prefix[[index]] + lower_after
    lower_cdf[[index]] <- if (lower_denominator > 0) {
      lower_prefix[[index]] / lower_denominator
    } else 0
    upper_cdf[[index]] <- if (upper_denominator > 0) {
      upper_prefix[[index]] / upper_denominator
    } else 0
  }
  tolerance <- 64 * .Machine$double.eps
  lower_candidates <- which(upper_cdf >= probability - tolerance)
  upper_candidates <- which(lower_cdf >= probability + tolerance)
  lower_index <- if (length(lower_candidates)) {
    lower_candidates[[1L]]
  } else 1L
  upper_index <- if (length(upper_candidates)) {
    upper_candidates[[1L]]
  } else bin_count
  c(lower = lower_index, upper = max(lower_index, upper_index))
}

.dsvert_dp_quantile_rows <- function(x, probs) {
  rows <- vector("list", length(x$variables))
  for (index in seq_along(x$variables)) {
    variable <- x$variables[[index]]
    histogram <- x$histograms[[index]]
    counts <- as.numeric(histogram$counts)
    grid <- as.numeric(histogram$grid)
    bin_lower <- c(x$lower_bounds[[index]], head(grid, -1L))
    total <- sum(counts)
    cumulative <- cumsum(counts)
    radius <- histogram$cell_noise_radius_simultaneous_95
    values <- lapply(probs, function(probability) {
      point_index <- if (total <= 0) {
        NA_integer_
      } else if (probability == 0) {
        which(counts > 0)[[1L]]
      } else if (probability == 1) {
        tail(which(counts > 0), 1L)
      } else {
        which(cumulative / total >= probability)[[1L]]
      }
      band <- .dsvert_dp_quantile_band_indices(
        counts, probability, radius)
      data.frame(
        variable = variable,
        probability = probability,
        status = if (is.na(point_index)) {
          "dp_projected_histogram_empty"
        } else {
          "ok_binned_postprocessed_estimate"
        },
        bin_id = if (is.na(point_index)) {
          NA_character_
        } else {
          paste0(variable, "::bin_", point_index)
        },
        bin_index = point_index,
        bin_lower = if (is.na(point_index)) {
          NA_real_
        } else bin_lower[[point_index]],
        bin_upper = if (is.na(point_index)) {
          NA_real_
        } else grid[[point_index]],
        bin_left_closed = if (is.na(point_index)) {
          NA
        } else point_index == 1L,
        bin_right_closed = if (is.na(point_index)) NA else TRUE,
        projected_histogram_mass_dp = total,
        cell_noise_radius_simultaneous_95 = radius,
        mechanism_grid_lower_95 = bin_lower[[band[["lower"]]]],
        mechanism_grid_upper_95 = grid[[band[["upper"]]]],
        stringsAsFactors = FALSE)
    })
    rows[[index]] <- do.call(rbind, values)
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Binned DP quantiles from one validated describe release
#'
#' Computes fixed-grid quantiles by deterministic client-side
#' post-processing of one already released `ds.vertDPDescribe()` object. No
#' server call is made and the additional privacy cost is `(0, 0)`. The result
#' identifies a public histogram bin and its interval; it is not an exact
#' sample quantile and does not interpolate within bins.
#'
#' @param x An intact released `ds.vertDPDescribe` object with either legacy
#'   capsule provenance or stateless synopsis provenance.
#' @param probs Finite public probabilities in `[0, 1]`. Duplicates are
#'   removed and the result is sorted by probability within each variable.
#' @return A `ds.vertDPQuantile` data frame with the selected public bin,
#'   interval, projected DP histogram mass, and a simultaneous 95-percent
#'   mechanism/grid region. Sampling uncertainty is excluded.
#' @export
ds.vertDPQuantile <- function(x, probs = c(0.25, 0.5, 0.75)) {
  x <- .dsvert_dp_describe_postprocess_contract(x)
  provenance_kind <- if ("release_provenance" %in% names(x)) {
    "synopsis"
  } else {
    "legacy"
  }
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs < 0 | probs > 1)) {
    stop("probs must contain finite probabilities in [0, 1]",
         call. = FALSE)
  }
  probs <- sort(unique(as.numeric(probs)), method = "radix")
  result <- .dsvert_dp_quantile_rows(x, probs)
  attr(result, "source_provenance") <- if (identical(
      provenance_kind, "legacy")) {
    list(
      source_class = "ds.vertDPDescribe",
      analysis_id = x$analysis_id,
      analysis_version = x$analysis_version,
      server = x$server,
      capsule_id = x$capsule_id,
      final_vector_root = x$final_vector_root,
      coordinate_order_sha256 = x$coordinate_order_sha256,
      mechanism = x$mechanism,
      epsilon = x$epsilon,
      delta = x$delta,
      implementation_delta = x$implementation_delta,
      adjacency = x$adjacency,
      privacy_epoch = x$privacy_epoch,
      noise_key_id = x$noise_key_id)
  } else {
    list(
      source_class = "ds.vertDPDescribe",
      analysis_id = x$analysis_id,
      analysis_version = x$analysis_version,
      server = x$server,
      artifact_key = x$artifact_key,
      execution_id = x$execution_id,
      contract_sha256 = x$contract_sha256,
      attempt_sha256 = x$attempt_sha256,
      source_contract_sha256 = x$source_contract_sha256,
      result_set_sha256 = x$result_set_sha256,
      final_vector_root = x$final_vector_root,
      coordinate_order_sha256 = x$coordinate_order_sha256,
      release_provenance = x$release_provenance,
      privacy = x$privacy,
      composition_rule = x$composition_rule,
      security_claim = x$security_claim,
      mechanism = x$mechanism,
      epsilon = x$epsilon,
      delta = x$delta,
      implementation_delta = x$implementation_delta,
      adjacency = x$adjacency)
  }
  attr(result, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  attr(result, "additional_server_calls") <- 0L
  attr(result, "postprocessing_only") <- TRUE
  attr(result, "quantile_definition") <- paste(
    "first public histogram bin whose cumulative coordinatewise-",
    "nonnegative DP-projected mass reaches probability times total mass;",
    "probability zero selects the first positive-mass bin")
  attr(result, "bin_interval_semantics") <- paste(
    "the first bin is [public_lower, endpoint]; subsequent bins are",
    "(previous_endpoint, endpoint]")
  attr(result, "estimate_scope") <- paste(
    "bin-identified fixed-grid postprocessed estimate; no exact sample",
    "quantile and no within-bin interpolation")
  attr(result, "histogram_mass_scope") <- paste(
    "coordinatewise-nonnegative DP-projected valid-value histogram mass;",
    "the fixed invalid bin is excluded and the mass is not silently forced",
    "to equal the separately noised effective count")
  attr(result, "mechanism_band_confidence") <- 0.95
  attr(result, "mechanism_band_scope") <- x$quantile_band_scope
  attr(result, "uncertainty_scope") <- x$uncertainty_scope
  class(result) <- c("ds.vertDPQuantile", class(result))
  result
}

#' Binned DP median from one validated describe release
#'
#' A release-only convenience wrapper around `ds.vertDPQuantile()` at
#' probability `0.5`. It performs no DSI operation and has additional privacy
#' cost `(0, 0)`.
#'
#' @param x An intact released `ds.vertDPDescribe` object with either legacy
#'   capsule provenance or stateless synopsis provenance.
#' @return A `ds.vertDPMedian` data frame with one binned median per released
#'   variable and the same mechanism/grid metadata as `ds.vertDPQuantile()`.
#' @export
ds.vertDPMedian <- function(x) {
  result <- ds.vertDPQuantile(x, probs = 0.5)
  class(result) <- c("ds.vertDPMedian", class(result))
  result
}
