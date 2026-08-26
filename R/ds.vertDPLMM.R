# Internal validator and deterministic post-processing for the signed
# random-intercept LMM capsule artifact.  The caller supplies only public DP
# coordinates after synopsis verification; this code has no DSI path.

.DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION <-
  "bounded-normalized-random-intercept-moments-v1"
.DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION <-
  "bounded-normalized-random-intercept-fixed-sufficient-statistics-v2"
.DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION <-
  "bounded-normalized-random-intercept-fixed-sufficient-statistics-v3"

.dsvert_dp_lmm_fixed_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  version <- if (is.list(artifact)) artifact$version else NULL
  fixed_ml <- identical(
    version, .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION)
  fixed_reml <- identical(
    version, .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "cluster", "predictors", "predictor_order", "intercept",
    "design_terms", "observation_capacity", "max_patients_per_cluster",
    "variance_ratio_grid", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "quantization_contract", "statistic_maximum",
    "source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
    "natural_l1_sensitivity", "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis", "estimation_scope",
    "implementation_state", "cross_owner_state")
  if (fixed_reml) required <- c(required, "estimation_profile")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    (fixed_ml || fixed_reml) &&
    identical(artifact$spec_version, if (fixed_reml) {
      "random_intercept_fixed_v3"
    } else {
      "random_intercept_fixed_v2"
    }) &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid fixed-effect random-intercept LMM ",
      "artifact '", analysis_id, "' for dataset '", data_name, "'"))
  }
  estimation_profile <- if (fixed_reml) artifact$estimation_profile else "ml"
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "LMM outcome")
  cluster <- artifact$cluster
  levels <- cluster$levels
  if (is.character(levels) && length(levels) &&
      all(vapply(levels, .dsvert_dp_is_string, logical(1L)))) {
    levels <- enc2utf8(unname(levels))
  } else {
    levels <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
      levels, "LMM cluster levels"), error = function(error) character())
  }
  cluster_valid <- .dsvert_dp_has_exact_names(cluster, c("column", "levels")) &&
    .dsvert_dp_is_string(cluster$column) && length(levels) >= 2L &&
    !anyDuplicated(levels) && all(nzchar(trimws(levels))) &&
    !identical(cluster$column, outcome$column)
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
    artifact$predictor_order, "LMM predictor order"),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
    artifact$design_terms, "LMM design terms"),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order) &&
    !outcome$column %in% predictor_order && !cluster$column %in% predictor_order
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) {
        .dsvert_dp_gaussian_bound(predictors[[variable]], "LMM predictor")
      }), predictor_order), error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(
      predictors, function(value) !identical(value$column, outcome$column) &&
        !identical(value$column, cluster$column), logical(1L)))
  }
  grid <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
    artifact$variance_ratio_grid, "LMM fixed variance-ratio grid"),
    error = function(error) numeric())
  grid_valid <- length(grid) &&
    !anyNA(grid) && all(is.finite(grid)) && all(grid >= 0) &&
    identical(grid[[1L]], 0) && all(diff(grid) > 0)
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(
    as.numeric(artifact$max_patients_per_cluster))
  design_count <- 1L + length(predictor_order)
  gram_count <- design_count * (design_count + 1L) / 2L
  summary_count <- gram_count + design_count + 1L
  coordinate_count <- (cluster_capacity + 1L) * (summary_count + 1L)
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "LMM statistic maxima"),
    error = function(error) numeric())
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) {
    2
  } else if (identical(adjacency, "add_remove_patient")) {
    1
  } else {
    NA_real_
  }
  raw_l1 <- if (is.finite(cluster_capacity) && is.finite(bits)) {
    multiplier * (3 + summary_count * (1 + 2 * cluster_capacity^2)) * scale
  } else NA_real_
  raw_l2_lower <- if (is.finite(cluster_capacity) && is.finite(bits)) {
    multiplier * sqrt(3 + summary_count * (1 + 4 * cluster_capacity^4)) *
      scale
  } else NA_real_
  expected_maximum <- if (is.finite(cluster_capacity)) c(
    capacity, rep(capacity * scale, summary_count),
    rep(c(capacity, rep(capacity * cluster_capacity * scale,
                       summary_count)), cluster_capacity)) else numeric()
  quantization <- artifact$quantization_contract
  valid <- cluster_valid && isTRUE(predictors_valid) &&
    is.character(estimation_profile) && length(estimation_profile) == 1L &&
    !is.na(estimation_profile) && estimation_profile %in% c("ml", "reml") &&
    .dsvert_dp_is_integer(bits, 8L, 18L) && identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    .dsvert_dp_is_integer(cluster_capacity, 2L, capacity) &&
    isTRUE(grid_valid) && isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count,
                          coordinate_count, coordinate_count) &&
    identical(artifact$coordinate_order, paste(
      "n_then_global_xtx_upper_xty_yty_then_each_cluster_size_from_1",
      "through_C_as_count_xtx_upper_xty_yty_v2", sep = "_")) &&
    identical(artifact$source_coordinate_scaling,
              "counts_left_shifted_to_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "clip_finite_complete_outcome_predictor_rows_then_mean_once_per",
      "admitted_patient_and_require_one_consistent_public_cluster_level_v2",
      sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "missing_or_nonfinite_outcome_predictor_or_missing_or_inconsistent",
      "cluster_excludes_the_patient_from_every_LMM_coordinate_v2",
      sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient_vector_and_one_consistent_cluster_with",
      "public_cluster_size_cap_v2", sep = "_")) &&
    .dsvert_dp_has_exact_names(quantization, c(
      "version", "input_rounding", "common_lattice")) &&
    identical(quantization$version,
              "random-intercept-fixed-common-lattice-quantization-v1") &&
    identical(quantization$input_rounding,
              "nearest_integer_ties_to_even_r_v1") &&
    identical(quantization$common_lattice, "numeric_grid_v1") &&
    identical(maximum, expected_maximum) &&
    identical(as.numeric(artifact$source_raw_l1_sensitivity), raw_l1) &&
    is.finite(as.numeric(artifact$source_raw_l2_sensitivity)) &&
    as.numeric(artifact$source_raw_l2_sensitivity) >= raw_l2_lower &&
    as.numeric(artifact$source_raw_l2_sensitivity) <= raw_l2_lower *
      (1 + 1e-8) &&
    identical(as.numeric(artifact$natural_l1_sensitivity), raw_l1 / scale) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity),
                     as.numeric(artifact$source_raw_l2_sensitivity) / scale,
                     tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis, paste(
      "one_patient_changes_n_and_at_most_two_cluster_size_bins_and",
      "bounded_squared_grid_cluster_summaries_with_replace_one_as",
      "two_add_remove_changes_v2", sep = "_")) &&
    identical(artifact$estimation_scope, paste(
      "bounded_random_intercept_GLS_fixed_effects_finite_signed",
      "variance_ratio_grid", toupper(estimation_profile),
      "profile_v1", sep = "_"))
  if (!isTRUE(valid)) {
    stop("The signed fixed-effect random-intercept LMM descriptor is invalid",
         call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$cluster <- list(column = enc2utf8(cluster$column), levels = levels)
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$variance_ratio_grid <- grid
  artifact$estimation_profile <- estimation_profile
  artifact$coordinate_count <- as.integer(coordinate_count)
  artifact
}

.dsvert_dp_lmm_fixed_coordinate_upper <- function(artifact) {
  dimension <- length(artifact$design_terms)
  summary_count <- dimension * (dimension + 1L) / 2L + dimension + 1L
  capacity <- as.numeric(artifact$observation_capacity)
  cluster_capacity <- as.numeric(artifact$max_patients_per_cluster)
  c(capacity, rep(capacity, summary_count),
    rep(c(capacity, rep(capacity * cluster_capacity, summary_count)),
        cluster_capacity))
}

.dsvert_dp_lmm_fixed_counts <- function(n, observed, cluster_capacity) {
  n <- as.integer(round(n))
  if (!is.finite(n) || n < 2L) return(NULL)
  counts <- pmax(0L, as.integer(floor(observed + 0.5)))
  if (length(counts) != cluster_capacity || anyNA(counts)) return(NULL)
  total <- sum(counts * seq_len(cluster_capacity))
  if (total > n) {
    excess <- total - n
    for (size in rev(seq_len(cluster_capacity))) {
      removable <- min(counts[[size]], excess %/% size)
      counts[[size]] <- counts[[size]] - removable
      excess <- excess - removable * size
    }
    if (excess != 0L) counts[] <- 0L
  }
  total <- sum(counts * seq_len(cluster_capacity))
  if (total > n) return(NULL)
  counts[[1L]] <- counts[[1L]] + n - total
  counts
}

.dsvert_dp_lmm_fixed_psd_global <- function(global) {
  dimension <- nrow(global$xtx)
  augmented <- rbind(
    cbind(global$xtx, global$xty),
    c(global$xty, global$yty))
  augmented <- (augmented + t(augmented)) / 2
  decomposition <- eigen(augmented, symmetric = TRUE)
  scale <- max(1, max(abs(augmented)))
  floor <- 128 * .Machine$double.eps * scale
  values <- pmax(decomposition$values, floor)
  projected <- tcrossprod(
    decomposition$vectors * rep(values, each = nrow(decomposition$vectors)),
    decomposition$vectors)
  projected <- (projected + t(projected)) / 2
  terms <- colnames(global$xtx) %||% paste0("x", seq_len(dimension))
  dimnames(projected) <- list(c(terms, ".outcome"),
                              c(terms, ".outcome"))
  list(
    n = global$n,
    xtx = projected[seq_len(dimension), seq_len(dimension), drop = FALSE],
    xty = stats::setNames(
      projected[seq_len(dimension), dimension + 1L], terms),
    yty = projected[[dimension + 1L, dimension + 1L]],
    applied = !isTRUE(all.equal(augmented, projected, tolerance = 0)))
}

.dsvert_dp_lmm_fixed_moments <- function(coordinates, artifact) {
  upper <- .dsvert_dp_lmm_fixed_coordinate_upper(artifact)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates))) {
    stop("The released fixed-effect random-intercept LMM coordinates are invalid",
         call. = FALSE)
  }
  projected <- pmin(upper, pmax(0, coordinates))
  dimension <- length(artifact$design_terms)
  summary_count <- dimension * (dimension + 1L) / 2L + dimension + 1L
  n <- as.integer(round(projected[[1L]]))
  n <- min(as.integer(artifact$observation_capacity), max(0L, n))
  count_positions <- 1L + summary_count +
    (seq_len(artifact$max_patients_per_cluster) - 1L) *
      (summary_count + 1L) + 1L
  counts <- .dsvert_dp_lmm_fixed_counts(
    n, projected[count_positions], artifact$max_patients_per_cluster)
  if (is.null(counts)) {
    return(.dsvert_lmm_gls_non_identifiable("inconsistent_cluster_counts"))
  }
  unpack <- function(values) {
    cursor <- 1L
    xtx <- matrix(0, dimension, dimension,
                  dimnames = list(artifact$design_terms,
                                  artifact$design_terms))
    for (right in seq_len(dimension)) {
      for (left in seq_len(right)) {
        xtx[left, right] <- xtx[right, left] <- values[[cursor]]
        cursor <- cursor + 1L
      }
    }
    list(xtx = xtx,
         xty = stats::setNames(values[seq.int(cursor, cursor + dimension - 1L)],
                               artifact$design_terms),
         yty = values[[cursor + dimension]])
  }
  global <- unpack(projected[seq.int(2L, summary_count + 1L)])
  global$n <- n
  by_size <- vector("list", artifact$max_patients_per_cluster)
  for (size in seq_along(by_size)) {
    start <- count_positions[[size]] + 1L
    item <- unpack(projected[seq.int(start, start + summary_count - 1L)])
    by_size[[size]] <- c(list(count = counts[[size]]), item)
  }
  fit <- .dsvert_lmm_random_intercept_gls(
    global[c("n", "xtx", "xty", "yty")], by_size,
    artifact$variance_ratio_grid, artifact$estimation_profile)
  algebraic_projection <- FALSE
  if (!identical(fit$status, "ok")) {
    stabilized <- .dsvert_dp_lmm_fixed_psd_global(
      global[c("n", "xtx", "xty", "yty")])
    fallback <- .dsvert_lmm_random_intercept_gls(
      stabilized[c("n", "xtx", "xty", "yty")], by_size, 0,
      artifact$estimation_profile)
    if (identical(fallback$status, "ok")) {
      fallback$reason <- "dp_psd_projected_zero_variance_profile"
      fit <- fallback
      algebraic_projection <- isTRUE(stabilized$applied)
    }
  }
  if (!identical(fit$status, "ok")) return(fit)
  y_span <- artifact$outcome$upper - artifact$outcome$lower
  coefficient <- fit$coefficients
  original <- coefficient
  original[["(Intercept)"]] <- artifact$outcome$lower + y_span *
    coefficient[["(Intercept)"]]
  for (variable in artifact$predictor_order) {
    bound <- artifact$predictors[[variable]]
    x_span <- bound$upper - bound$lower
    slope <- y_span * coefficient[[variable]] / x_span
    original[["(Intercept)"]] <- original[["(Intercept)"]] -
      slope * bound$lower
    original[[variable]] <- slope
  }
  fit$normalized_coefficients <- coefficient
  fit$coefficients <- original
  fit$sigma2 <- fit$sigma2 * y_span^2
  fit$sigma_b2 <- fit$sigma_b2 * y_span^2
  fit$n_obs <- n
  fit$cluster_count <- sum(counts)
  fit$projected_summary <- list(
    n = n, cluster_counts = counts,
    coordinate_projection_applied = !isTRUE(all.equal(
      projected, coordinates, tolerance = 0)),
    algebraic_psd_projection_applied = algebraic_projection)
  fit$projection_applied <- isTRUE(
    fit$projected_summary$coordinate_projection_applied) ||
    algebraic_projection
  fit
}

.dsvert_dp_lmm_any_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  if (is.list(artifact) && artifact$version %in% c(
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION)) {
    return(.dsvert_dp_lmm_fixed_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
  if (is.list(artifact) && identical(
        artifact$version,
        .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    return(.dsvert_dp_lmm_random_slope_grid_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
  .dsvert_dp_lmm_artifact(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity)
}

.dsvert_dp_lmm_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "cluster", "observation_capacity",
    "max_patients_per_cluster", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "quantization_contract", "statistic_maximum",
    "source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
    "natural_l1_sensitivity", "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis", "estimation_scope",
    "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version,
              .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "random_intercept_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid same-owner random-intercept LMM ",
      "artifact '", analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "LMM outcome")
  cluster <- artifact$cluster
  cluster_levels <- cluster$levels
  if (is.character(cluster_levels) && length(cluster_levels) > 0L &&
      all(vapply(cluster_levels, .dsvert_dp_is_string, logical(1L)))) {
    cluster_levels <- enc2utf8(unname(cluster_levels))
  } else {
    cluster_levels <- tryCatch(
      .dsvert_dp_capsule_manifest_string_array(
        cluster_levels, "LMM cluster levels"),
      error = function(error) character())
  }
  cluster_valid <- .dsvert_dp_has_exact_names(cluster, c("column", "levels")) &&
    .dsvert_dp_is_string(cluster$column) && length(cluster_levels) >= 2L &&
    !anyDuplicated(cluster_levels) && all(nzchar(trimws(cluster_levels))) &&
    !identical(cluster$column, outcome$column)
  maximum <- tryCatch(
    .dsvert_dp_capsule_manifest_numbers(
      artifact$statistic_maximum, "LMM statistic maxima"),
    error = function(error) numeric())
  grid_bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(
    as.numeric(artifact$max_patients_per_cluster))
  expected_multiplier <- if (identical(adjacency,
                                        "replace_one_fixed_cohort")) {
    2
  } else if (identical(adjacency, "add_remove_patient")) {
    1
  } else {
    NA_real_
  }
  raw_l1 <- if (is.finite(cluster_capacity) && is.finite(grid_bits)) {
    expected_multiplier * (2 * cluster_capacity + 2 + 3 * scale)
  } else {
    NA_real_
  }
  raw_l2_lower <- if (is.finite(cluster_capacity) && is.finite(grid_bits)) {
    expected_multiplier * sqrt(
      2 + (2 * cluster_capacity - 1)^2 + 2 * scale^2 +
        (scale + 1)^2)
  } else {
    NA_real_
  }
  quantization <- artifact$quantization_contract
  quantization_valid <- .dsvert_dp_has_exact_names(quantization, c(
    "version", "input_rounding", "sum_y_max_abs_error_normalized",
    "sum_y_sq_max_abs_error_normalized",
    "cluster_mean_sq_max_abs_error_normalized")) &&
    identical(quantization$version,
              "random-intercept-unit-moment-quantization-v1") &&
    identical(quantization$input_rounding,
              "nearest_integer_ties_to_even_r_v1") &&
    isTRUE(all.equal(
      as.numeric(quantization$sum_y_max_abs_error_normalized),
      capacity / (2 * scale), tolerance = 1e-12)) &&
    isTRUE(all.equal(
      as.numeric(quantization$sum_y_sq_max_abs_error_normalized),
      capacity / (2 * scale), tolerance = 1e-12)) &&
    isTRUE(all.equal(
      as.numeric(quantization$cluster_mean_sq_max_abs_error_normalized),
      3 * capacity / (2 * scale) + capacity / (4 * scale^2),
      tolerance = 1e-12))
  valid <- cluster_valid && .dsvert_dp_is_integer(grid_bits, 8, 18) &&
    identical(2^grid_bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    .dsvert_dp_is_integer(cluster_capacity, 2, capacity) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, 6, 6) &&
    identical(artifact$coordinate_order, paste(
      "n_then_cluster_count_then_cluster_size_sq_then_quantized_sum_y",
      "then_quantized_sum_y_sq_then_quantized_cluster_mean_sq_v1",
      sep = "_")) &&
    identical(artifact$source_coordinate_scaling,
              "three_counts_then_three_common_lattice_moments_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "clip_finite_outcome_then_mean_once_per_admitted_patient_and",
      "require_one_consistent_public_cluster_level_v1", sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "missing_or_nonfinite_outcome_or_missing_or_inconsistent_cluster",
      "excludes_the_patient_from_all_six_coordinates_v1", sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient_outcome_and_one_consistent_cluster_level",
      "with_public_cluster_size_cap_v1", sep = "_")) &&
    isTRUE(quantization_valid) && length(maximum) == 6L &&
    identical(maximum, c(
      capacity, capacity, capacity * cluster_capacity,
      rep(capacity * scale, 3L))) &&
    identical(as.numeric(artifact$source_raw_l1_sensitivity), raw_l1) &&
    is.finite(as.numeric(artifact$source_raw_l2_sensitivity)) &&
    as.numeric(artifact$source_raw_l2_sensitivity) >= raw_l2_lower &&
    as.numeric(artifact$source_raw_l2_sensitivity) <= raw_l2_lower *
      (1 + 1e-8) &&
    identical(as.numeric(artifact$natural_l1_sensitivity), raw_l1 / scale) &&
    isTRUE(all.equal(
      as.numeric(artifact$natural_l2_sensitivity),
      as.numeric(artifact$source_raw_l2_sensitivity) / scale,
      tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis, paste(
      "one_patient_changes_three_counts_by_1_1_and_at_most",
      "2C_minus_1_and_three_quantized_moments_by_S_S_and_S_plus_1",
      "with_replace_one_as_two_add_remove_changes_v1", sep = "_")) &&
    identical(artifact$estimation_scope,
              "bounded_random_intercept_method_of_moments_no_fixed_covariates_v1")
  if (!isTRUE(valid)) {
    stop("The signed random-intercept LMM descriptor is invalid",
         call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$cluster <- list(
    column = enc2utf8(cluster$column),
    levels = enc2utf8(unname(cluster_levels)))
  artifact$coordinate_count <- 6L
  artifact
}

.dsvert_dp_lmm_moments <- function(coordinates, artifact) {
  if (!is.numeric(coordinates) || length(coordinates) != 6L ||
      anyNA(coordinates) || any(!is.finite(coordinates))) {
    stop("The released random-intercept LMM coordinates are invalid",
         call. = FALSE)
  }
  .dsvert_lmm_random_intercept_moments(
    list(
      n = coordinates[[1L]], clusters = coordinates[[2L]],
      cluster_size_sq = coordinates[[3L]],
      sum_y_normalized = coordinates[[4L]],
      sum_y_sq_normalized = coordinates[[5L]],
      sum_cluster_mean_sq_normalized = coordinates[[6L]]),
    outcome_lower = artifact$outcome$lower,
    outcome_upper = artifact$outcome$upper,
    observation_capacity = artifact$observation_capacity,
    cluster_capacity = artifact$max_patients_per_cluster)
}

.dsvert_dp_lmm_synopsis_release <- function(
    data_name, analysis_id, server = NULL, datasources = NULL, .aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  if (!is.null(server)) server <- .dsvert_dp_server(server, datasources)
  run <- .dsvert_dp_synopsis_vector_run(
    datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  metadata <- .dsvert_dp_vector_public_metadata(context)
  scale <- as.numeric(context$lattice$output_lattice_scale)
  count_block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count",
    description = "signed admitted-count capacity block")
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  artifact <- .dsvert_dp_lmm_any_artifact(
    context$manifest, data_name, analysis_id, server,
    context$adjacency, scale, capacity)
  blocks <- .dsvert_dp_capsule_vector_blocks(
    context$layout, "gaussian_models", dataset = data_name,
    owner_peer = artifact$owner_peer)
  blocks <- blocks[vapply(
    blocks, function(block) identical(block$key, analysis_id), logical(1L))]
  signed_descriptor <- tryCatch(
    context$manifest$workload$families$gaussian_models$
      artifacts[[analysis_id]], error = function(error) NULL)
  if (length(blocks) != 1L || !identical(
        .dsvert_joint_dp_client_json(blocks[[1L]]$descriptor),
        .dsvert_joint_dp_client_json(signed_descriptor))) {
    stop("The signed random-intercept LMM artifact does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(
    context$release, blocks[[1L]])
  coordinate_upper <- if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION)) {
    .dsvert_dp_lmm_fixed_coordinate_upper(artifact)
  } else if (identical(artifact$version,
                       .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else c(
    capacity, capacity, capacity * artifact$max_patients_per_cluster,
    rep(capacity, 3L))
  if (length(coordinates) != artifact$coordinate_count ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > coordinate_upper)) {
    stop("The released random-intercept LMM block violates its signed bounds",
         call. = FALSE)
  }
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity,
                 "session_transport_anchored") ||
      !verification$artifact$version %in% c(
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    stop("The random-intercept LMM Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(
    context = context, metadata = metadata,
    artifact = verification$artifact, block = blocks[[1L]],
    coordinates = verification$coordinates,
    moment = verification$validated_moment,
    certificate = certificate, verification = verification,
    scale = scale, capacity = capacity)
}

.dsvert_dp_lmm_impl <- function(
    data_name, analysis_id, server = NULL, datasources = NULL, .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(
    analysis_id, "analysis_id")
  released <- .dsvert_dp_lmm_synopsis_release(
    data_name, analysis_id, server, datasources, .aggregate)
  moment <- released$moment
  artifact <- released$artifact
  coefficients <- if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    moment$coefficients
  } else moment$coefficient %||% moment$coefficients
  result <- c(released$metadata, list(
    status = moment$status, analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    family = if (identical(artifact$version,
                           .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
      "gaussian_random_slope"
    } else "gaussian_random_intercept",
    estimand = artifact$estimation_scope,
    coefficient = coefficients,
    coefficients = coefficients,
    sigma2 = moment$sigma2 %||% NULL,
    sigma_b2 = moment$sigma_b2 %||% NULL,
    random_effect_covariance = moment$random_effect_covariance %||% NULL,
    random_effect_order = moment$random_effect_order %||% NULL,
    selected_candidate = moment$selected_candidate %||% NULL,
    selected_dp_negative_log_likelihood =
      moment$selected_dp_negative_log_likelihood %||% NULL,
    icc = moment$icc %||% NULL,
    effective_cluster_size = moment$effective_cluster_size %||% NULL,
    n_obs = moment$n_obs %||% moment$projected_summary[["n"]],
    cluster_count = moment$cluster_count %||%
      moment$projected_summary[["clusters"]],
    projected_moments = moment$projected_summary,
    moment_projection_applied = moment$projection_applied,
    identifiability_reason = moment$reason %||% NULL,
    accuracy = list(
      confidence = released$verification$accuracy_simultaneous_95$confidence,
      simultaneous_abs_mechanism_radius =
        released$verification$accuracy_simultaneous_95$radius,
      coordinate_count = artifact$coordinate_count,
      max_abs_quantization_normalized = artifact$quantization_contract,
      additional_privacy_cost = c(epsilon = 0, delta = 0)),
    inference = list(
      classical_standard_errors = NULL, p_values = NULL,
      confidence_intervals = NULL, sampling_inference_available = FALSE),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE,
    provenance_certificate = released$certificate,
    disclosure_guard = list(
      satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("ds.vertDPLMM", "list")
  result
}

#' Bounded random-intercept model from a canonical signed DP Synopsis
#'
#' Fits the documented intercept-only Gaussian random-intercept
#' method-of-moments estimand from six bounded, jointly released coordinates.
#' It never calls a legacy mixed-model endpoint or receives cluster-level
#' statistics.
#'
#' @param data_name Signed protected dataset name.
#' @param analysis_id Custodian-configured signed random-intercept artifact id.
#' @param server Optional expected signed outcome-owner server name.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPLMM` object with DP variance components when the
#'   released projected moments are identifiable.
#' @export
ds.vertDPLMM <- function(
    data_name, analysis_id, server = NULL, datasources = NULL) {
  resolved <- .dsvert_federation_argument(data_name, datasources)
  .dsvert_dp_lmm_impl(
    resolved$value, analysis_id, server, resolved$datasources,
    DSI::datashield.aggregate)
}

#' Verify a random-intercept LMM Synopsis certificate without DSI
#'
#' @param x A `ds.vertDPLMM` result or its provenance certificate.
#' @param trusted_pinset Optional trusted named peer-to-Ed25519-public-key map.
#' @return The verified public Synopsis provenance and projected LMM moments.
#' @export
ds.validateDPLMMCertificate <- function(x, trusted_pinset = NULL) {
  certificate <- if (is.list(x) && identical(
        x$version, .DSVERT_DP_GAUSSIAN_SYNOPSIS_CERTIFICATE_VERSION)) {
    x
  } else if (is.list(x)) {
    x$provenance_certificate
  } else {
    NULL
  }
  verified <- ds.validateDPGaussianCertificate(
    certificate, trusted_pinset = trusted_pinset)
  if (!verified$artifact$version %in% c(
      .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION,
      .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
      .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION,
      .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    stop("The certificate is not a signed LMM artifact",
         call. = FALSE)
  }
  verified
}
