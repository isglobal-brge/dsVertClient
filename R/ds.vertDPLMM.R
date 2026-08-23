# Internal validator and deterministic post-processing for the signed
# random-intercept LMM capsule artifact.  The caller supplies only public DP
# coordinates after synopsis verification; this code has no DSI path.

.DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION <-
  "bounded-normalized-random-intercept-moments-v1"

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
  artifact <- .dsvert_dp_lmm_artifact(
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
  coordinate_upper <- c(
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
      !identical(verification$artifact$version,
                 .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION)) {
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
  result <- c(released$metadata, list(
    status = moment$status, analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    family = "gaussian_random_intercept",
    estimand = "bounded_random_intercept_method_of_moments",
    coefficient = moment$coefficient,
    coefficients = moment$coefficient,
    sigma2 = moment$sigma2 %||% NULL,
    sigma_b2 = moment$sigma_b2 %||% NULL,
    icc = moment$icc %||% NULL,
    effective_cluster_size = moment$effective_cluster_size %||% NULL,
    n_obs = moment$projected_summary[["n"]],
    cluster_count = moment$projected_summary[["clusters"]],
    projected_moments = moment$projected_summary,
    moment_projection_applied = moment$projection_applied,
    identifiability_reason = moment$reason %||% NULL,
    accuracy = list(
      confidence = released$verification$accuracy_simultaneous_95$confidence,
      simultaneous_abs_mechanism_radius =
        released$verification$accuracy_simultaneous_95$radius,
      coordinate_count = artifact$coordinate_count,
      max_abs_quantization_normalized = artifact$quantization_contract[c(
        "sum_y_max_abs_error_normalized",
        "sum_y_sq_max_abs_error_normalized",
        "cluster_mean_sq_max_abs_error_normalized")],
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
  if (!identical(verified$artifact$version,
                .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION)) {
    stop("The certificate is not a random-intercept LMM artifact",
         call. = FALSE)
  }
  verified
}
