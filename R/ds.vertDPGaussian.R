# Bounded Gaussian regression from one signed sticky capsule.
# The released sufficient statistics are DP values. No exact statistic,
# protected row, discovery endpoint, or legacy Gaussian route is reachable.

.dsvert_dp_gaussian_identifier <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)) {
    stop(what, " must be one non-empty capsule identifier", call. = FALSE)
  }
  enc2utf8(value)
}

.dsvert_dp_gaussian_reserved <- function(message) {
  stop("reserved_not_materialized: ", message, call. = FALSE)
}

.dsvert_dp_gaussian_reference <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    return(NULL)
  }
  value <- enc2utf8(value)
  separator <- regexpr("$", value, fixed = TRUE)[[1L]]
  if (separator > 0L && grepl(
        "$", substring(value, separator + 1L), fixed = TRUE)) {
    return(NULL)
  }
  parts <- if (separator < 0L) value else c(
    substring(value, 1L, separator - 1L),
    substring(value, separator + 1L))
  if (!length(parts) %in% 1:2 ||
      any(!grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", parts))) {
    return(NULL)
  }
  list(reference = value,
       server = if (length(parts) == 2L) parts[[1L]] else NULL,
       column = parts[[length(parts)]])
}

.dsvert_dp_gaussian_reference_matches <- function(value, descriptor,
                                                   default_owner = NULL) {
  reference <- .dsvert_dp_gaussian_reference(value)
  owner <- descriptor$owner_peer %||% default_owner
  !is.null(reference) && is.list(descriptor) &&
    identical(reference$column, descriptor$column) &&
    (is.null(reference$server) || identical(reference$server, owner))
}

.dsvert_dp_gaussian_bound <- function(value, what) {
  required <- c("column", "lower", "upper")
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      anyDuplicated(names(value)) || !setequal(names(value), required) ||
      !.dsvert_dp_is_string(value$column) ||
      !.dsvert_dp_is_number(value$lower) ||
      !.dsvert_dp_is_number(value$upper) || value$lower >= value$upper) {
    stop("The signed Gaussian ", what, " bound is invalid", call. = FALSE)
  }
  list(
    column = enc2utf8(value$column), lower = as.numeric(value$lower),
    upper = as.numeric(value$upper))
}

.dsvert_dp_gaussian_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  family <- tryCatch(
    manifest$workload$families$gaussian_models,
    error = function(error) NULL)
  artifact <- if (is.list(family$artifacts)) {
    family$artifacts[[analysis_id]]
  } else NULL
  if (is.list(artifact) && identical(
        artifact$version,
        .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION)) {
    return(.dsvert_dp_gaussian_cross_artifact(
      artifact, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
  if (is.list(artifact) && identical(
        artifact$version,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION)) {
    return(.dsvert_dp_lmm_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
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
  if (is.list(artifact) && identical(
        artifact$version,
        .DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION)) {
    return(.dsvert_dp_gee_ar1_grid_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
  if (is.list(artifact) && artifact$version %in% c(
        .DSVERT_CLIENT_DP_GLMM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_POISSON_GLMM_GRID_ARTIFACT_VERSION)) {
    return(.dsvert_dp_glmm_grid_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
  if (is.list(artifact) && artifact$version %in%
      unname(.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS)) {
    family_name <- names(.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS)[[
      match(artifact$version,
            unname(.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS))]]
    return(.dsvert_dp_glm_grid_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity, family_name))
  }
  if (is.list(artifact) && identical(
        artifact$version, .DSVERT_CLIENT_DP_NB_GRID_ARTIFACT_VERSION)) {
    return(.dsvert_dp_nb_grid_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
  if (is.list(artifact) && identical(
        artifact$version, .DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION)) {
    return(.dsvert_dp_multinom_grid_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
  if (is.list(artifact) && identical(
        artifact$version, .DSVERT_CLIENT_DP_ORDINAL_GRID_ARTIFACT_VERSION)) {
    return(.dsvert_dp_ordinal_grid_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale,
      capacity))
  }
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "predictors", "predictor_order", "intercept",
    "design_terms", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "repeated_record_policy", "missingness_policy",
    "contribution_domain", "count_gram_intercept_policy",
    "statistic_maximum",
    "source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
    "natural_l1_sensitivity", "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis", "regularization_policy",
    "implementation_state", "cross_owner_state")
  basic <- is.list(artifact) && !is.null(names(artifact)) &&
    !anyNA(names(artifact)) && !anyDuplicated(names(artifact)) &&
    setequal(names(artifact), required) &&
    identical(
      artifact$version,
      "bounded-normalized-gaussian-sufficient-statistics-v1") &&
    .dsvert_dp_is_string(artifact$spec_version) &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    is.logical(artifact$intercept) && length(artifact$intercept) == 1L &&
    !is.na(artifact$intercept) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid same-owner Gaussian artifact '",
      analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "outcome")
  predictor_order <- .dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "Gaussian predictor order", sorted = TRUE)
  if (!length(predictor_order) || outcome$column %in% predictor_order ||
      !is.list(artifact$predictors) ||
      !identical(names(artifact$predictors), predictor_order)) {
    stop("The signed Gaussian predictor contract is invalid", call. = FALSE)
  }
  predictors <- lapply(predictor_order, function(variable) {
    bound <- .dsvert_dp_gaussian_bound(
      artifact$predictors[[variable]], "predictor")
    if (!identical(bound$column, variable)) {
      stop("The signed Gaussian predictor order changed", call. = FALSE)
    }
    bound
  })
  names(predictors) <- predictor_order
  design_terms <- c(
    if (isTRUE(artifact$intercept)) "(Intercept)" else character(),
    predictor_order)
  observed_terms <- .dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "Gaussian design terms", sorted = FALSE)
  q <- length(design_terms)
  coordinate_count <- q * (q + 1) / 2 + q + 2
  maxima <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "Gaussian maxima"),
  error = function(error) numeric())
  expected_raw_l1 <- 1 + (coordinate_count - 1) * scale
  expected_raw_l2 <- sqrt(1 + (coordinate_count - 1) * scale^2)
  valid <- identical(observed_terms, design_terms) &&
    .dsvert_dp_is_integer(artifact$numeric_grid_bits, 8, 18) &&
    identical(2^as.numeric(artifact$numeric_grid_bits), scale) &&
    .dsvert_dp_is_integer(
      artifact$coordinate_count, coordinate_count, coordinate_count) &&
    identical(
      artifact$coordinate_order,
      paste0(
        "n_then_xtx_upper_column_major_then_xty_design_order_then_yty_",
        "v1")) &&
    identical(
      artifact$repeated_record_policy,
      paste0(
        "clip_finite_rows_then_mean_each_variable_once_per_admitted_unit_",
        "v1")) &&
    identical(
      artifact$missingness_policy,
      "complete_case_across_outcome_and_all_predictors_v1") &&
    identical(
      artifact$contribution_domain,
      "one_vector_of_normalized_monomials_in_closed_unit_interval_v1") &&
    identical(
      artifact$count_gram_intercept_policy,
      paste0(
        "n_is_complete_case_count_and_moment_upper_bound_gram11_governs_",
        "the_solve_no_averaging_v1")) &&
    length(maxima) == coordinate_count &&
    identical(maxima, c(capacity, rep(
      capacity * scale, coordinate_count - 1L))) &&
    isTRUE(all.equal(
      as.numeric(artifact$source_raw_l1_sensitivity), expected_raw_l1,
      tolerance = 1e-12)) &&
    isTRUE(all.equal(
      as.numeric(artifact$source_raw_l2_sensitivity), expected_raw_l2,
      tolerance = 1e-12)) &&
    identical(
      as.numeric(artifact$natural_l1_sensitivity), coordinate_count) &&
    isTRUE(all.equal(
      as.numeric(artifact$natural_l2_sensitivity), sqrt(coordinate_count),
      tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(
      artifact$adjacency_sensitivity_basis,
      paste0(
        "zero_missing_complete_case_vs_all_one_complete_unit_is_worst_",
        "case_for_add_remove_and_replace_one")) &&
    identical(
      artifact$regularization_policy,
      "none_in_release_explicit_client_postprocessing_only_v1")
  if (!isTRUE(valid)) {
    stop("The signed Gaussian sufficient-statistic descriptor is invalid",
         call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- design_terms
  artifact$coordinate_count <- as.integer(coordinate_count)
  artifact
}

.dsvert_dp_gaussian_cross_artifact <- function(
    artifact, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "predictors", "predictor_order", "input_variable_order",
    "participating_peers", "computation_peers", "intercept",
    "design_terms", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling",
    "private_input_layout", "repeated_record_policy", "missingness_policy",
    "contribution_domain", "count_gram_intercept_policy",
    "statistic_maximum", "source_raw_l1_sensitivity",
    "source_raw_l2_sensitivity", "natural_l1_sensitivity",
    "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis", "quantization_contract",
    "numeric_certificate", "transcript", "alignment_contract",
    "regularization_policy", "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version,
              .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "v2") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    is.logical(artifact$intercept) && length(artifact$intercept) == 1L &&
    !is.na(artifact$intercept) &&
    identical(artifact$implementation_state,
              "cross_owner_exact_gc_materialized") &&
    identical(artifact$cross_owner_state,
              "exact_gc_to_joint_dp_vector_v1")
  if (!isTRUE(basic)) {
    stop("The signed capsule has no valid cross-owner Gaussian artifact '",
         analysis_id, "' for outcome dataset '", data_name, "'",
         call. = FALSE)
  }
  bound <- function(value, what) {
    fields <- c("column", "dataset", "owner_peer", "lower", "upper")
    if (!.dsvert_dp_has_exact_names(value, fields) ||
        !.dsvert_dp_is_string(value$column) ||
        !.dsvert_dp_is_string(value$dataset) ||
        !.dsvert_dp_is_string(value$owner_peer) ||
        !.dsvert_dp_is_number(value$lower) ||
        !.dsvert_dp_is_number(value$upper) || value$lower >= value$upper) {
      stop("The signed cross-owner Gaussian ", what, " bound is invalid",
           call. = FALSE)
    }
    list(
      column = enc2utf8(value$column), dataset = enc2utf8(value$dataset),
      owner_peer = enc2utf8(value$owner_peer),
      lower = as.numeric(value$lower), upper = as.numeric(value$upper))
  }
  outcome <- bound(artifact$outcome, "outcome")
  predictor_order <- .dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "Gaussian predictor order", sorted = TRUE)
  if (!length(predictor_order) || outcome$column %in% predictor_order ||
      !is.list(artifact$predictors) ||
      !identical(names(artifact$predictors), predictor_order)) {
    stop("The signed cross-owner Gaussian predictor contract is invalid",
         call. = FALSE)
  }
  predictors <- lapply(predictor_order, function(variable) {
    value <- bound(artifact$predictors[[variable]], "predictor")
    if (!identical(value$column, variable)) {
      stop("The signed cross-owner Gaussian predictor order changed",
           call. = FALSE)
    }
    value
  })
  names(predictors) <- predictor_order
  input_order <- .dsvert_dp_capsule_manifest_strings(
    artifact$input_variable_order, "Gaussian private input order",
    sorted = FALSE)
  participants <- .dsvert_dp_gaussian_cross_names_client(
    artifact$participating_peers, "participant list")
  computation <- .dsvert_dp_gaussian_cross_names_client(
    artifact$computation_peers, "computation-peer list")
  descriptor_owners <- c(
    vapply(predictors, `[[`, character(1L), "owner_peer"),
    outcome$owner_peer)
  design_terms <- c(
    if (isTRUE(artifact$intercept)) "(Intercept)" else character(),
    predictor_order)
  observed_terms <- .dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "Gaussian design terms", sorted = FALSE)
  q <- length(design_terms)
  coordinate_count <- q * (q + 1) / 2 + q + 2
  maxima <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "Gaussian maxima"),
  error = function(error) numeric())
  grid_bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  required_signed_bits <- if (is.finite(grid_bits)) {
    as.integer(ceiling(log2(capacity + 1))) + 2L * as.integer(grid_bits) + 3L
  } else {
    NA_integer_
  }
  product_error <- 2 + 1 / (4 * scale)
  quantization <- artifact$quantization_contract
  certificate <- artifact$numeric_certificate
  transcript <- artifact$transcript
  alignment <- artifact$alignment_contract
  valid <- identical(outcome$dataset, data_name) &&
    identical(outcome$owner_peer, artifact$owner_peer) &&
    identical(input_order, c(predictor_order, outcome$column)) &&
    !anyDuplicated(participants) &&
    identical(participants, sort(participants, method = "radix")) &&
    setequal(participants, unique(descriptor_owners)) &&
    length(computation) == 2L && !anyDuplicated(computation) &&
    identical(computation, sort(computation, method = "radix")) &&
    identical(observed_terms, design_terms) &&
    .dsvert_dp_is_integer(artifact$numeric_grid_bits, 8, 18) &&
    identical(2^grid_bits, scale) &&
    .dsvert_dp_is_integer(
      artifact$coordinate_count, coordinate_count, coordinate_count) &&
    identical(artifact$coordinate_order,
              paste0("n_then_xtx_upper_column_major_then_xty_design_order_",
                     "then_yty_v2")) &&
    identical(artifact$source_coordinate_scaling,
              "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$private_input_layout,
              paste0("capacity_padded_value_then_validity_per_signed_",
                     "variable_manifest_order_v1")) &&
    identical(artifact$repeated_record_policy,
              paste0("clip_finite_rows_then_mean_each_variable_once_per_",
                     "admitted_unit_v1")) &&
    identical(artifact$missingness_policy,
              paste0("complete_case_mask_remains_secret_shared_through_",
                     "joint_noise_v1")) &&
    identical(artifact$contribution_domain,
              paste0("round_normalized_inputs_then_exact_floor_ring128_",
                     "products_on_closed_unit_interval_v1")) &&
    identical(artifact$count_gram_intercept_policy,
              paste0("n_and_all_moments_share_one_secret_complete_case_",
                     "mask_and_are_released_only_after_joint_dp_v1")) &&
    length(maxima) == coordinate_count &&
    identical(maxima, rep(capacity * scale, coordinate_count)) &&
    isTRUE(all.equal(
      as.numeric(artifact$source_raw_l1_sensitivity),
      coordinate_count * scale, tolerance = 1e-12)) &&
    isTRUE(all.equal(
      as.numeric(artifact$source_raw_l2_sensitivity),
      sqrt(coordinate_count) * scale, tolerance = 1e-12)) &&
    identical(as.numeric(artifact$natural_l1_sensitivity),
              as.numeric(coordinate_count)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity),
                     sqrt(coordinate_count), tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis,
              paste0("zero_missing_complete_case_vs_all_one_complete_unit_",
                     "is_worst_case_for_add_remove_and_replace_one")) &&
    .dsvert_dp_has_exact_names(quantization, c(
      "input_rounding", "product_rounding",
      "per_product_max_abs_error_lattice_steps",
      "per_product_max_abs_error_normalized",
      "per_sum_max_abs_error_lattice_steps",
      "same_owner_v1_numerically_identical")) &&
    identical(quantization$input_rounding,
              "nearest_integer_ties_to_even_r_v1") &&
    identical(quantization$product_rounding,
              "exact_signed_floor_after_division_by_scale_v1") &&
    isTRUE(all.equal(
      as.numeric(quantization$per_product_max_abs_error_lattice_steps),
      product_error, tolerance = 1e-12)) &&
    isTRUE(all.equal(
      as.numeric(quantization$per_product_max_abs_error_normalized),
      product_error / scale, tolerance = 1e-12)) &&
    isTRUE(all.equal(
      as.numeric(quantization$per_sum_max_abs_error_lattice_steps),
      capacity * product_error, tolerance = 1e-12)) &&
    identical(quantization$same_owner_v1_numerically_identical, FALSE) &&
    .dsvert_dp_has_exact_names(certificate, c(
      "version", "ring_bits", "frac_bits", "required_signed_bits",
      "operand_maximum", "raw_product_maximum",
      "accumulated_coordinate_maximum", "truncation", "comparison",
      "modular_wrap_proved_absent", "overflow_behavior")) &&
    identical(certificate$version,
              "dsvert-cross-gaussian-numeric-certificate-v1") &&
    identical(as.numeric(certificate$ring_bits), 128) &&
    identical(as.numeric(certificate$frac_bits), grid_bits) &&
    identical(as.numeric(certificate$required_signed_bits),
              as.numeric(required_signed_bits)) &&
    identical(as.numeric(certificate$operand_maximum), scale) &&
    identical(as.numeric(certificate$raw_product_maximum), scale^2) &&
    identical(as.numeric(certificate$accumulated_coordinate_maximum),
              capacity * scale) &&
    identical(certificate$truncation,
              "exact_signed_floor_gc_ot_or_direct_wide_v1") &&
    identical(certificate$comparison,
              "not_used_after_custodian_bound_clipping") &&
    identical(certificate$modular_wrap_proved_absent, TRUE) &&
    identical(certificate$overflow_behavior, "typed_abort_before_commit") &&
    .dsvert_dp_has_exact_names(transcript, c(
      "version", "padded_units", "variable_count",
      "validity_product_rounds", "masked_value_rounds",
      "moment_product_rounds", "data_dependent_branches",
      "exact_intermediate_release_count")) &&
    identical(transcript$version,
              "dsvert-cross-gaussian-fixed-transcript-v1") &&
    identical(as.numeric(transcript$padded_units), as.numeric(capacity)) &&
    identical(as.numeric(transcript$variable_count),
              as.numeric(length(input_order))) &&
    identical(as.numeric(transcript$validity_product_rounds),
              as.numeric(length(input_order) - 1L)) &&
    identical(as.numeric(transcript$masked_value_rounds), 1) &&
    identical(as.numeric(transcript$moment_product_rounds), 1) &&
    identical(as.numeric(transcript$data_dependent_branches), 0) &&
    identical(as.numeric(transcript$exact_intermediate_release_count), 0) &&
    .dsvert_dp_has_exact_names(alignment, c(
      "version", "public_patient_dependent_hash", "mismatch_behavior")) &&
    identical(alignment$version,
              "private-psi-ordered-manifest-consensus-v1") &&
    identical(alignment$public_patient_dependent_hash, FALSE) &&
    identical(alignment$mismatch_behavior,
              "typed_non_prealigned_cohort_failure") &&
    identical(artifact$regularization_policy,
              "none_in_release_explicit_client_postprocessing_only_v1")
  if (!isTRUE(valid)) {
    stop("The signed cross-owner Gaussian descriptor is invalid",
         call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$input_variable_order <- input_order
  artifact$participating_peers <- participants
  artifact$computation_peers <- computation
  artifact$design_terms <- design_terms
  artifact$coordinate_count <- as.integer(coordinate_count)
  artifact
}

.dsvert_dp_gaussian_unpack <- function(coordinates, artifact, capacity) {
  q <- length(artifact$design_terms)
  n <- min(capacity, max(0, coordinates[[1L]]))
  if (!is.finite(n) || n <= 0) {
    .dsvert_stop_non_identifiable(
      "The DP complete-case count is non-positive.",
      reason = "non_positive_dp_complete_case_count")
  }
  cursor <- 2L
  gram <- matrix(
    0, q, q,
    dimnames = list(artifact$design_terms, artifact$design_terms))
  gram_released <- gram
  for (right in seq_len(q)) {
    for (left in seq_len(right)) {
      released_value <- coordinates[[cursor]]
      value <- min(n, max(0, released_value))
      gram_released[left, right] <- gram_released[right, left] <-
        released_value
      gram[left, right] <- gram[right, left] <- value
      cursor <- cursor + 1L
    }
  }
  cross <- pmin(n, pmax(0, coordinates[cursor:(cursor + q - 1L)]))
  names(cross) <- artifact$design_terms
  cursor <- cursor + q
  outcome_square <- min(n, max(0, coordinates[[cursor]]))
  augmented <- rbind(
    cbind(gram, cross),
    c(cross, outcome_square))
  augmented <- (augmented + t(augmented)) / 2
  decomposition <- eigen(augmented, symmetric = TRUE)
  clipped <- pmax(0, decomposition$values)
  projected <- sweep(
    decomposition$vectors, 2L, clipped, `*`) %*%
    t(decomposition$vectors)
  projected <- (projected + t(projected)) / 2
  projected_gram <- projected[seq_len(q), seq_len(q), drop = FALSE]
  dimnames(projected_gram) <- dimnames(gram)
  projected_cross <- projected[seq_len(q), q + 1L]
  names(projected_cross) <- artifact$design_terms
  tolerance <- 256 * .Machine$double.eps * max(
    1, max(abs(decomposition$values))) * (q + 1L)
  list(
    n = n, gram_released = gram_released, gram_raw_clamped = gram,
    cross_raw_clamped = cross,
    outcome_square_raw_clamped = outcome_square,
    augmented_projected = projected, gram = projected_gram,
    cross = projected_cross, outcome_square = projected[q + 1L, q + 1L],
    projection = list(
      method = "augmented_second_moment_eigenvalue_clipping_v1",
      exact_nearest_with_all_moment_constraints = FALSE,
      input_min_eigenvalue = min(decomposition$values),
      output_min_eigenvalue = min(clipped),
      clipped_eigenvalues = as.integer(sum(decomposition$values < 0)),
      frobenius_distance = sqrt(sum((projected - augmented)^2)),
      numerical_tolerance = tolerance))
}

.dsvert_dp_gaussian_solve <- function(moment, artifact, ridge) {
  q <- nrow(moment$gram)
  penalty <- rep(1, q)
  if (isTRUE(artifact$intercept)) penalty[[1L]] <- 0
  gram_values <- eigen(
    moment$gram, symmetric = TRUE, only.values = TRUE)$values
  gram_tolerance <- 256 * .Machine$double.eps *
    max(1, max(abs(gram_values))) * q
  rank <- sum(gram_values > gram_tolerance)
  condition <- if (rank == q) {
    max(gram_values) / min(gram_values)
  } else {
    Inf
  }
  if (ridge == 0 && rank < q) {
    .dsvert_stop_non_identifiable(
      paste(
        "The released DP Gaussian design is singular; an explicit positive",
        "ridge value requests a different estimand."),
      reason = "singular_dp_gaussian_design")
  }
  system <- moment$gram + diag(ridge * penalty, q)
  decomposition <- eigen(system, symmetric = TRUE)
  tolerance <- 256 * .Machine$double.eps *
    max(1, max(abs(decomposition$values))) * q
  if (min(decomposition$values) <= tolerance) {
    .dsvert_stop_non_identifiable(
      "The requested Gaussian system is not invertible.",
      reason = "noninvertible_dp_gaussian_system")
  }
  coefficients <- as.numeric(decomposition$vectors %*%
    ((t(decomposition$vectors) %*% moment$cross) /
       decomposition$values))
  names(coefficients) <- artifact$design_terms
  rss <- as.numeric(
    moment$outcome_square - 2 * crossprod(coefficients, moment$cross) +
      crossprod(coefficients, moment$gram %*% coefficients))
  rss <- max(0, rss)
  list(
    coefficients = coefficients, residual_second_moment = rss,
    identifiability = list(
      scope = "released_DP_projected_normalized_design",
      design_dimension = q, numerical_rank = as.integer(rank),
      full_rank = rank == q, minimum_eigenvalue = min(gram_values),
      maximum_eigenvalue = max(gram_values),
      condition_number = condition, tolerance = gram_tolerance),
    system = list(
      ridge = ridge, intercept_penalized = FALSE,
      minimum_eigenvalue = min(decomposition$values),
      condition_number = max(decomposition$values) /
        min(decomposition$values)))
}

.dsvert_dp_gaussian_original_coefficients <- function(
    coefficients, artifact) {
  outcome_range <- artifact$outcome$upper - artifact$outcome$lower
  predictor_ranges <- vapply(artifact$predictors, function(bound) {
    bound$upper - bound$lower
  }, numeric(1L))
  predictor_lowers <- vapply(
    artifact$predictors, `[[`, numeric(1L), "lower")
  normalized_slopes <- coefficients[artifact$predictor_order]
  slopes <- outcome_range * normalized_slopes / predictor_ranges
  normalized_intercept <- if (isTRUE(artifact$intercept)) {
    coefficients[["(Intercept)"]]
  } else {
    0
  }
  offset <- artifact$outcome$lower +
    outcome_range * normalized_intercept -
    sum(slopes * predictor_lowers)
  c(`(Intercept)` = offset, slopes)
}

.dsvert_dp_gaussian_synopsis_release <- function(
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
  artifact <- .dsvert_dp_gaussian_artifact(
    context$manifest, data_name, analysis_id, server,
    context$adjacency, scale, capacity)
  blocks <- .dsvert_dp_capsule_vector_blocks(
    context$layout, "gaussian_models", dataset = data_name,
    owner_peer = artifact$owner_peer)
  blocks <- blocks[vapply(
    blocks, function(block) identical(block$key, analysis_id), logical(1L))]
  signed_descriptor <- tryCatch(
    context$manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  if (length(blocks) != 1L ||
      !identical(
        .dsvert_joint_dp_client_json(blocks[[1L]]$descriptor),
        .dsvert_joint_dp_client_json(signed_descriptor))) {
    stop("The signed Gaussian artifact does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(
    context$release, blocks[[1L]])
  if (length(coordinates) != artifact$coordinate_count ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > capacity)) {
    stop("The released Gaussian Synopsis block violates its signed bounds",
         call. = FALSE)
  }
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity,
                 "session_transport_anchored")) {
    stop("The Gaussian Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(
    context = context, metadata = metadata,
    artifact = verification$artifact,
    block = blocks[[1L]], coordinates = verification$coordinates,
    moment = verification$validated_moment,
    certificate = certificate, verification = verification,
    scale = scale, capacity = capacity)
}

#' Bounded Gaussian regression from a canonical signed DP Synopsis
#'
#' Fits a descriptive Gaussian linear model from one signed sufficient-
#' statistic artifact. Outcome and predictors are clipped to custodian-owned
#' public bounds, collapsed once per admitted patient, normalized to `[0,1]`,
#' and restricted to a fixed complete-case estimand. The function performs one
#' canonical Synopsis retrieval and then only deterministic client
#' post-processing.
#'
#' @param data_name Signed protected dataset name.
#' @param analysis_id Custodian-configured signed Gaussian artifact id.
#' @param ridge Explicit non-negative ridge penalty on normalized predictors.
#'   The default zero preserves the unpenalized bounded-statistic estimand.
#' @param server Optional expected signed outcome-owner server name.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPGaussian` object. It contains no classical standard
#'   errors, p-values, individual fitted values, residuals, or scores.
#' @export
ds.vertDPGaussian <- function(
    data_name, analysis_id, ridge = 0, server = NULL,
    datasources = NULL) {
  resolved <- .dsvert_federation_argument(data_name, datasources)
  .dsvert_dp_gaussian_impl(
    resolved$value, analysis_id, ridge, server, resolved$datasources,
    DSI::datashield.aggregate)
}

.dsvert_dp_gaussian_impl <- function(
    data_name, analysis_id, ridge = 0, server = NULL,
    datasources = NULL, .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(
    analysis_id, "analysis_id")
  if (!is.numeric(ridge) || length(ridge) != 1L || is.na(ridge) ||
      !is.finite(ridge) || ridge < 0) {
    stop("ridge must be one finite non-negative number", call. = FALSE)
  }
  ridge <- as.numeric(ridge)
  released <- .dsvert_dp_gaussian_synopsis_release(
    data_name, analysis_id, server, datasources, .aggregate)
  context <- released$context
  artifact <- released$artifact
  coordinates <- released$coordinates
  scale <- released$scale
  capacity <- released$capacity
  provenance_certificate <- released$certificate
  provenance_verification <- released$verification
  if (identical(artifact$version,
                .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION) ||
      artifact$version %in% c(
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    stop("The signed artifact is an LMM; use ds.vertDPLMM",
         call. = FALSE)
  }
  simultaneous <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest,
    coordinate_count = artifact$coordinate_count,
    confidence = 0.95, maximum_error = capacity)
  n_interval <- c(
    lower = max(0, coordinates[[1L]] - simultaneous$radius),
    upper = min(capacity, coordinates[[1L]] + simultaneous$radius))
  quantization <- if (identical(
      artifact$version,
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION)) {
    as.numeric(
      artifact$quantization_contract$per_sum_max_abs_error_lattice_steps) /
      scale
  } else {
    0.5 * n_interval[["upper"]] / scale
  }
  moment_intervals <- cbind(
    lower = pmax(
      0, coordinates[-1L] - simultaneous$radius - quantization),
    upper = pmin(
      n_interval[["upper"]],
      coordinates[-1L] + simultaneous$radius + quantization))
  moment <- released$moment
  fit <- .dsvert_dp_gaussian_solve(moment, artifact, ridge)
  original <- .dsvert_dp_gaussian_original_coefficients(
    fit$coefficients, artifact)
  result <- c(released$metadata, list(
    status = "ok", analysis_id = analysis_id,
    cohort_id = provenance_verification$cohort_id,
    logical_snapshot = provenance_verification$logical_snapshot,
    certificate_sha256 = provenance_certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    participating_servers = if (is.null(artifact$participating_peers)) {
      artifact$owner_peer
    } else {
      artifact$participating_peers
    },
    computation_servers = if (is.null(artifact$computation_peers)) {
      NULL
    } else {
      artifact$computation_peers
    },
    family = "gaussian",
    estimand = paste(
      "patient-collapsed, public-bound-clipped, normalized complete-case",
      "Gaussian sufficient-statistic regression"),
    coefficients_normalized = fit$coefficients,
    coefficients_original_scale = original,
    coefficients = original,
    normalized_intercept_included = artifact$intercept,
    original_scale_intercept_semantics = if (isTRUE(artifact$intercept)) {
      "fitted_intercept_after_inverse_bound_normalization"
    } else {
      "deterministic_bound_transform_offset_not_a_fitted_intercept"
    },
    ridge = ridge,
    regularized = ridge > 0,
    regularization_estimand = if (ridge > 0) {
      "explicit_predictor_only_ridge_on_normalized_design"
    } else {
      "unpenalized_normalized_least_squares"
    },
    n_obs = moment$n,
    n_obs_definition = paste0(
      "DP_noisy_complete_case_count_coordinate_not_averaged_with_",
      "Gram11"),
    residual_second_moment_dp = fit$residual_second_moment,
    sufficient_statistics_dp = list(
      gram_projected = moment$gram, cross_projected = moment$cross,
      outcome_square_projected = moment$outcome_square),
    count_gram_intercept_reconciliation = if (isTRUE(
      artifact$intercept)) {
      list(
        policy = artifact$count_gram_intercept_policy,
        n_coordinate_dp = moment$n,
        gram_intercept_intercept_dp =
          moment$gram_released[[1L, 1L]],
        discrepancy = moment$gram_released[[1L, 1L]] - moment$n,
        silently_averaged = FALSE)
    } else {
      list(
        policy = artifact$count_gram_intercept_policy,
        applicable = FALSE, silently_averaged = FALSE)
    },
    augmented_moment_projection = moment$projection,
    identifiability = fit$identifiability,
    solve_diagnostics = fit$system,
    accuracy = list(
      confidence = simultaneous$confidence,
      simultaneous_abs_mechanism_radius = simultaneous$radius,
      coordinate_count = artifact$coordinate_count,
      complete_case_count_interval = n_interval,
      non_count_coordinate_intervals = moment_intervals,
      max_abs_quantization_per_non_count_sum = quantization,
      method = simultaneous$method,
      implementation_tv_upper_bound =
        simultaneous$implementation_tv_upper_bound,
      additional_privacy_cost = simultaneous$additional_privacy_cost,
      coefficient_regions_available = FALSE),
    uncertainty_scope = paste(
      "The certificate covers simultaneous Synopsis-mechanism noise and",
      "deterministic quantization of bounded finite-snapshot sufficient",
      "statistics; nonlinear coefficient regions and population sampling",
      "inference are not claimed"),
    inference = list(
      classical_standard_errors = NULL, p_values = NULL,
      confidence_intervals = NULL, sampling_inference_available = FALSE),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    individual_fitted_values_available = FALSE,
    individual_residuals_available = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE,
    provenance_certificate = provenance_certificate,
    disclosure_guard = list(
      satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("ds.vertDPGaussian", "list")
  verified <- ds.validateDPGaussianCertificate(result)
  result$provenance_integrity <- verified$integrity_valid
  result$provenance_authenticity <- verified$authenticity
  result
}

.dsvert_dp_gaussian_glm_adapter <- function(
    explicit_arguments, formula, data, x_vars, y_server, family, lambda,
    no_intercept, data_name, y_var, missing, verbose, datasources,
    analysis_id) {
  analysis_id <- .dsvert_dp_gaussian_identifier(
    analysis_id, "dp_analysis_id")
  if (!identical(family, "gaussian")) {
    stop("dp_analysis_id is available only for family='gaussian'",
         call. = FALSE)
  }
  unsupported <- intersect(explicit_arguments, c(
    "max_iter", "tol", "log_n", "offset", "weights", "ring",
    "binomial_sigmoid_intervals", "eta_privacy", "keep_session",
    "std_mode", "start", "compute_se", "compute_deviance",
    "gradient_only", "numeric_backend"))
  if (length(unsupported)) {
    stop(
      "The explicit capsule Gaussian route does not accept: ",
      paste(unsupported, collapse = ", "),
      "; these arguments belong to the legacy iterative estimand",
      call. = FALSE)
  }
  if ("missing" %in% explicit_arguments &&
      !identical(missing, "complete_case_capsule")) {
    stop(paste(
      "The signed capsule Gaussian estimand has fixed complete-case",
      "missingness; use missing='complete_case_capsule' or omit missing"),
    call. = FALSE)
  }
  if (is.null(formula) ||
      !(inherits(formula, "formula") ||
        (is.character(formula) && length(formula) == 1L &&
         !is.na(formula) && grepl("~", formula, fixed = TRUE)))) {
    stop(paste(
      "The capsule Gaussian ds.vertGLM adapter requires an explicit",
      "additive formula and data name"), call. = FALSE)
  }
  formula_spec <- .dsvert_plain_formula(formula)
  if (!is.null(y_var) && !identical(y_var, formula_spec$response)) {
    stop("y_var disagrees with the explicit capsule formula", call. = FALSE)
  }
  data_name <- if (!is.null(data_name)) data_name else data
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data")
  predictors <- if (is.null(x_vars)) formula_spec$predictors else x_vars
  requested_predictor_owners <- NULL
  if (is.list(predictors)) {
    if (is.null(names(predictors)) || anyNA(names(predictors)) ||
        any(!nzchar(names(predictors))) || anyDuplicated(names(predictors)) ||
        any(!vapply(predictors, function(value) {
          is.character(value) && !anyNA(value) && all(nzchar(value)) &&
            !anyDuplicated(value)
        }, logical(1L)))) {
      stop("x_vars must be a valid named owner-to-predictor list",
           call. = FALSE)
    }
    nonempty <- names(predictors)[lengths(predictors) > 0L]
    requested_predictor_owners <- unlist(lapply(nonempty, function(owner) {
      stats::setNames(rep(owner, length(predictors[[owner]])),
                      predictors[[owner]])
    }), use.names = TRUE)
    predictors <- unname(unlist(predictors[nonempty], use.names = FALSE))
    if (anyDuplicated(predictors)) {
      stop("A Gaussian predictor cannot be assigned to multiple owners",
           call. = FALSE)
    }
  }
  if (!is.character(predictors) || !length(predictors) ||
      anyNA(predictors) || any(!nzchar(predictors)) ||
      anyDuplicated(predictors)) {
    stop("The capsule Gaussian formula needs unique numeric predictors",
         call. = FALSE)
  }
  if (!setequal(predictors, formula_spec$predictors)) {
    stop("x_vars disagrees with the explicit capsule formula", call. = FALSE)
  }
  if (!is.null(y_server)) {
    y_server <- .dsvert_dp_gaussian_identifier(y_server, "y_server")
  }
  intercept <- isTRUE(formula_spec$intercept) && !isTRUE(no_intercept)
  if (isTRUE(verbose)) {
    message(paste(
      "Using the signed bounded complete-case Gaussian DP artifact;",
      "this is not the legacy iterative GLM estimand."))
  }
  result <- ds.vertDPGaussian(
    data_name = data_name, analysis_id = analysis_id, ridge = lambda,
    server = y_server, datasources = datasources)
  artifact <- result$signed_artifact
  owner_match <- is.null(requested_predictor_owners) || all(vapply(
    names(requested_predictor_owners), function(variable) {
      identical(artifact$predictors[[variable]]$owner_peer,
                unname(requested_predictor_owners[[variable]])) ||
        (is.null(artifact$predictors[[variable]]$owner_peer) &&
         identical(artifact$owner_peer,
                   unname(requested_predictor_owners[[variable]])))
    }, logical(1L)))
  if (!identical(artifact$outcome$column, formula_spec$response) ||
      !identical(artifact$predictor_order,
                 sort(enc2utf8(predictors), method = "radix")) ||
      !identical(artifact$intercept, intercept) || !isTRUE(owner_match)) {
    stop(paste(
      "The explicit formula does not match the signed Gaussian artifact;",
      "no legacy fallback was attempted"), call. = FALSE)
  }
  result$called_via <- "ds.vertGLM_explicit_dp_analysis_id"
  result$requested_formula <- paste(deparse(formula), collapse = " ")
  result$legacy_glm_estimand <- FALSE
  result
}

#' @export
print.ds.vertDPGaussian <- function(x, digits = 4, ...) {
  cat("Bounded Gaussian Regression (Sticky Joint DP Capsule)\n")
  cat("=====================================================\n\n")
  cat("Signed artifact:", x$analysis_id, "\n")
  cat("Outcome owner:", x$server, "\n")
  if (length(x$participating_servers) > 1L) {
    cat("Participating owners:",
        paste(x$participating_servers, collapse = ", "), "\n")
  }
  cat("DP complete-case count:", round(x$n_obs, digits), "\n")
  cat("Regularization:", x$regularization_estimand, "\n")
  cat("Released-design rank:", x$identifiability$numerical_rank, "/",
      x$identifiability$design_dimension, "\n\n")
  print(round(x$coefficients_original_scale, digits), ...)
  cat("\n", x$uncertainty_scope, "\n", sep = "")
  invisible(x)
}
