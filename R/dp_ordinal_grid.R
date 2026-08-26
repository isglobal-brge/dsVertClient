# Internal validator and post-processing for a custodian-signed finite ordinal
# cumulative-logit likelihood grid. It consumes only one public sticky Synopsis.

.DSVERT_CLIENT_DP_ORDINAL_GRID_ARTIFACT_VERSION <-
  "bounded-ordinal-likelihood-grid-v1"

.dsvert_dp_ordinal_grid_loss_bounds <- function(candidate_grid) {
  vapply(candidate_grid, function(candidate) {
    score_bound <- sum(abs(candidate$beta))
    gap <- min(diff(candidate$thresholds))
    tail <- 8 + score_bound
    probability_lower <- min(
      stats::plogis(-tail),
      stats::plogis(tail + gap) - stats::plogis(tail))
    if (!is.finite(probability_lower) || probability_lower <= 0) {
      stop("The signed ordinal loss bound is invalid", call. = FALSE)
    }
    -log(probability_lower)
  }, numeric(1L))
}

.dsvert_dp_ordinal_grid_outcome <- function(value) {
  required <- c("column", "levels", "ordered_levels")
  normalize <- function(levels) {
    if (is.list(levels) && is.null(names(levels))) {
      levels <- unlist(levels, use.names = FALSE)
    }
    levels
  }
  levels <- if (is.list(value)) normalize(value$levels) else NULL
  ordered_levels <- if (is.list(value)) normalize(value$ordered_levels) else NULL
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      anyDuplicated(names(value)) || !setequal(names(value), required) ||
      !.dsvert_dp_is_string(value$column) || !is.character(levels) ||
      !is.character(ordered_levels) || length(levels) < 3L ||
      anyNA(levels) || anyNA(ordered_levels) || any(!nzchar(levels)) ||
      any(!nzchar(ordered_levels)) || anyDuplicated(levels) ||
      anyDuplicated(ordered_levels) || !identical(levels, sort(levels,
        method = "radix")) || !setequal(levels, ordered_levels)) {
    stop("The signed ordinal outcome domain is invalid", call. = FALSE)
  }
  list(column = enc2utf8(value$column), levels = enc2utf8(levels),
       ordered_levels = enc2utf8(ordered_levels))
}

.dsvert_dp_ordinal_grid_candidates <- function(value, levels, dimension) {
  if (!is.list(value) || !length(value) || !is.null(names(value))) {
    stop("The signed ordinal candidate grid is invalid", call. = FALSE)
  }
  candidates <- lapply(value, function(candidate) {
    if (!is.list(candidate) || is.null(names(candidate)) ||
        !setequal(names(candidate), c("thresholds", "beta"))) {
      return(NULL)
    }
    thresholds <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$thresholds, "ordinal threshold row"),
      error = function(error) numeric())
    beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$beta, "ordinal beta row"),
      error = function(error) numeric())
    if (length(thresholds) != length(levels) - 1L ||
        length(beta) != dimension || any(!is.finite(thresholds)) ||
        any(!is.finite(beta)) || any(abs(thresholds) > 8) ||
        any(abs(beta) > 8) || any(diff(thresholds) < 1 / 256)) {
      return(NULL)
    }
    list(thresholds = as.numeric(thresholds), beta = as.numeric(beta))
  })
  keys <- if (length(candidates) && all(vapply(candidates, is.list,
                                               logical(1L)))) {
    vapply(candidates, .dsvert_joint_dp_client_json, character(1L))
  } else character()
  if (length(candidates) > 256L || !length(keys) || anyDuplicated(keys) ||
      !identical(keys, sort(keys, method = "radix"))) {
    stop("The signed ordinal candidate grid is invalid", call. = FALSE)
  }
  candidates
}

.dsvert_dp_ordinal_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "predictors", "predictor_order", "intercept",
    "design_terms", "observation_capacity", "candidate_grid",
    "candidate_order", "candidate_loss_bounds", "numeric_grid_bits",
    "coordinate_count", "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "statistic_maximum", "source_raw_l1_sensitivity",
    "source_raw_l2_sensitivity", "natural_l1_sensitivity",
    "natural_l2_sensitivity", "adjacency", "adjacency_sensitivity_basis",
    "estimation_scope", "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version, .DSVERT_CLIENT_DP_ORDINAL_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "ordinal_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) && .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid ordinal likelihood grid '",
      analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_ordinal_grid_outcome(artifact$outcome)
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "ordinal predictor order", sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "ordinal design terms", sorted = FALSE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order) &&
    !outcome$column %in% predictor_order
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) .dsvert_dp_gaussian_bound(
        predictors[[variable]], "ordinal predictor")), predictor_order),
      error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(
      predictors, function(value) !identical(value$column, outcome$column),
      logical(1L)))
  }
  candidates <- tryCatch(.dsvert_dp_ordinal_grid_candidates(
    artifact$candidate_grid, outcome$levels, 1L + length(predictor_order)),
    error = function(error) list())
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  loss_bounds <- if (length(candidates)) {
    tryCatch(.dsvert_dp_ordinal_grid_loss_bounds(candidates),
             error = function(error) numeric())
  } else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "ordinal candidate loss bounds"),
    error = function(error) numeric())
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "ordinal statistic maxima"),
    error = function(error) numeric())
  raw_per_candidate <- ceiling(loss_bounds * scale)
  expected_maximum <- capacity * raw_per_candidate
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) 2 else if (
    identical(adjacency, "add_remove_patient")) 1 else NA_real_
  raw_l1 <- multiplier * sum(raw_per_candidate)
  raw_l2 <- multiplier * sqrt(sum(raw_per_candidate^2))
  valid <- isTRUE(predictors_valid) && .dsvert_dp_is_integer(bits, 8L, 18L) &&
    identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    length(candidates) && isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$candidate_order,
              "canonical_ordinal_cumulative_logit_grid_v1") &&
    length(observed_loss_bounds) == length(candidates) &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, length(candidates),
                          length(candidates)) &&
    identical(artifact$coordinate_order, paste(
      "canonical_ordinal_candidate_grid_cumulative_logit_negative_log_likelihood_v1",
      sep = "_")) &&
    identical(artifact$source_coordinate_scaling,
              "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "require_one_categorical_outcome_and_mean_once_per_admitted",
      "patient_v1", sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "missing_outcome_or_missing_or_nonfinite_predictor_excludes_patient",
      "and_unknown_or_conflicting_nonmissing_outcome_rejects_v1", sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient_ordinal_cumulative_logit_negative_log_likelihood",
      "contribution_for_every_signed_candidate_v1", sep = "_")) &&
    identical(maximum, expected_maximum) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l1_sensitivity), raw_l1,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l2_sensitivity), raw_l2,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l1_sensitivity), raw_l1 / scale,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity), raw_l2 / scale,
                     tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis, paste(
      "one_patient_changes_one_candidate_loss_by_at_most_its_signed",
      "ordinal_cumulative_logit_loss_bound_v1", sep = "_")) &&
    identical(artifact$estimation_scope, paste(
      "bounded_ordinal_cumulative_logit_fixed_covariates_finite_signed",
      "candidate_grid_v1", sep = "_"))
  if (!isTRUE(valid)) stop("The signed ordinal grid descriptor is invalid",
                           call. = FALSE)
  artifact$outcome <- outcome
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$candidate_grid <- candidates
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$coordinate_count <- as.integer(length(candidates))
  artifact$statistic_maximum <- expected_maximum
  artifact
}

.dsvert_dp_ordinal_grid_moment <- function(coordinates, artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released ordinal grid violates its signed bounds", call. = FALSE)
  }
  candidate <- which.min(coordinates)[[1L]]
  selected <- artifact$candidate_grid[[candidate]]
  normalized_beta <- selected$beta
  ranges <- vapply(artifact$predictors, function(bound) {
    bound$upper - bound$lower
  }, numeric(1L))
  lowers <- vapply(artifact$predictors, `[[`, numeric(1L), "lower")
  beta <- normalized_beta[-1L] / ranges
  names(beta) <- artifact$predictor_order
  thresholds <- selected$thresholds - normalized_beta[[1L]] +
    sum(beta * lowers)
  names(thresholds) <- head(artifact$outcome$ordered_levels, -1L)
  list(status = "ok", coefficients = beta, normalized_coefficients = normalized_beta,
       thresholds = thresholds, normalized_thresholds = selected$thresholds,
       ordered_levels = artifact$outcome$ordered_levels,
       selected_candidate = as.integer(candidate),
       selected_dp_negative_log_likelihood = coordinates[[candidate]] /
         (2^as.numeric(artifact$numeric_grid_bits)),
       candidate_selection =
         "minimum_signed_finite_grid_dp_postprocessing_v1")
}

.dsvert_dp_ordinal_grid_synopsis_release <- function(
    data_name, analysis_id, server = NULL, datasources = NULL, .aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  if (!is.null(server)) server <- .dsvert_dp_server(server, datasources)
  run <- .dsvert_dp_synopsis_vector_run(datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  metadata <- .dsvert_dp_vector_public_metadata(context)
  scale <- as.numeric(context$lattice$output_lattice_scale)
  count_block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count", description = "signed admitted-count capacity block")
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  artifact <- .dsvert_dp_ordinal_grid_artifact(
    context$manifest, data_name, analysis_id, server, context$adjacency,
    scale, capacity)
  blocks <- .dsvert_dp_capsule_vector_blocks(
    context$layout, "gaussian_models", dataset = data_name,
    owner_peer = artifact$owner_peer)
  blocks <- blocks[vapply(blocks, function(block) {
    identical(block$key, analysis_id)
  }, logical(1L))]
  signed_descriptor <- tryCatch(context$manifest$workload$families$
    gaussian_models$artifacts[[analysis_id]], error = function(error) NULL)
  if (length(blocks) != 1L || !identical(
        .dsvert_joint_dp_client_json(blocks[[1L]]$descriptor),
        .dsvert_joint_dp_client_json(signed_descriptor))) {
    stop("The signed ordinal grid artifact does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  moment <- .dsvert_dp_ordinal_grid_moment(coordinates, artifact)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !identical(verification$artifact$version,
                 .DSVERT_CLIENT_DP_ORDINAL_GRID_ARTIFACT_VERSION)) {
    stop("The ordinal Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(context = context, metadata = metadata, artifact = verification$artifact,
       block = blocks[[1L]], coordinates = verification$coordinates,
       moment = verification$validated_moment %||% moment,
       certificate = certificate, verification = verification,
       scale = scale, capacity = capacity)
}

.dsvert_dp_ordinal_grid_impl <- function(
    formula, data_name, analysis_id, server = NULL, datasources = NULL,
    .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(analysis_id, "analysis_id")
  outcome <- as.character(formula[[2L]])
  predictors <- attr(stats::terms(formula), "term.labels")
  released <- .dsvert_dp_ordinal_grid_synopsis_release(
    data_name, analysis_id, server, datasources, .aggregate)
  artifact <- released$artifact
  if (!identical(artifact$outcome$column, outcome) ||
      !setequal(artifact$predictor_order, predictors)) {
    stop("formula must match the signed ordinal grid artifact", call. = FALSE)
  }
  moment <- released$moment
  result <- c(released$metadata, list(
    status = moment$status, analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    family = "ordinal_finite_grid", estimand = artifact$estimation_scope,
    coefficients = moment$coefficients,
    normalized_coefficients = moment$normalized_coefficients,
    thresholds = moment$thresholds,
    normalized_thresholds = moment$normalized_thresholds,
    ordered_levels = moment$ordered_levels,
    selected_candidate = moment$selected_candidate,
    selected_dp_negative_log_likelihood = moment$selected_dp_negative_log_likelihood,
    candidate_selection = moment$candidate_selection,
    covariance = NULL, std_errors = NULL, standard_errors = NULL,
    p_values = NULL, inference = "unavailable_for_finite_grid_ordinal",
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE, provenance_certificate = released$certificate,
    disclosure_guard = list(satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("dsvert_dp_ordinal_grid", "ds.vertOrdinal", "list")
  result
}
