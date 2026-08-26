# Internal validator and post-processing for a custodian-signed finite
# multinomial likelihood grid. It consumes only one public sticky Synopsis.

.DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION <-
  "bounded-multinomial-likelihood-grid-v1"

.dsvert_dp_multinom_grid_loss_bounds <- function(beta_grid, levels, reference) {
  non_reference <- setdiff(levels, reference)
  dimension <- length(beta_grid[[1L]]) / length(non_reference)
  vapply(beta_grid, function(beta) {
    score_bound <- max(colSums(abs(matrix(
      beta, nrow = dimension, ncol = length(non_reference)))))
    2 * score_bound + log(length(levels))
  }, numeric(1L))
}

.dsvert_dp_multinom_grid_outcome <- function(value) {
  required <- c("column", "levels", "reference")
  levels <- if (is.list(value) && is.list(value$levels) &&
      is.null(names(value$levels))) {
    unlist(value$levels, use.names = FALSE)
  } else if (is.list(value)) {
    value$levels
  } else {
    NULL
  }
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      anyDuplicated(names(value)) || !setequal(names(value), required) ||
      !.dsvert_dp_is_string(value$column) || !is.character(levels) ||
      length(levels) < 3L || anyNA(levels) ||
      any(!nzchar(levels)) || anyDuplicated(levels) ||
      !identical(levels, sort(levels, method = "radix")) ||
      !.dsvert_dp_is_string(value$reference) ||
      !value$reference %in% levels) {
    stop("The signed multinomial outcome domain is invalid", call. = FALSE)
  }
  list(column = enc2utf8(value$column), levels = enc2utf8(levels),
       reference = enc2utf8(value$reference))
}

.dsvert_dp_multinom_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "predictors", "predictor_order", "intercept",
    "design_terms", "observation_capacity", "beta_grid", "candidate_order",
    "candidate_loss_bounds", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "statistic_maximum", "source_raw_l1_sensitivity",
    "source_raw_l2_sensitivity", "natural_l1_sensitivity",
    "natural_l2_sensitivity", "adjacency", "adjacency_sensitivity_basis",
    "estimation_scope", "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version, .DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "multinomial_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) && .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid multinomial likelihood grid '",
      analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_multinom_grid_outcome(artifact$outcome)
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "multinomial predictor order", sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "multinomial design terms", sorted = FALSE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order) &&
    !outcome$column %in% predictor_order
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) .dsvert_dp_gaussian_bound(
        predictors[[variable]], "multinomial predictor")), predictor_order),
      error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(
      predictors, function(value) !identical(value$column, outcome$column),
      logical(1L)))
  }
  beta_grid <- artifact$beta_grid
  if (!is.list(beta_grid) || !length(beta_grid) || !is.null(names(beta_grid))) {
    beta_grid <- list()
  } else {
    beta_grid <- lapply(beta_grid, function(beta) tryCatch(
      .dsvert_dp_capsule_manifest_number_array(
        beta, "multinomial beta grid row"),
      error = function(error) numeric()))
  }
  beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
    .dsvert_joint_dp_client_json(as.list(beta))
  }, character(1L)) else character()
  dimension <- (length(outcome$levels) - 1L) * (1L + length(predictor_order))
  beta_valid <- length(beta_grid) && all(vapply(beta_grid, function(beta) {
    length(beta) == dimension && !anyNA(beta) && all(is.finite(beta)) &&
      all(abs(beta) <= 8)
  }, logical(1L))) && !anyDuplicated(beta_keys) &&
    identical(beta_keys, sort(beta_keys, method = "radix"))
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  loss_bounds <- if (isTRUE(beta_valid)) {
    .dsvert_dp_multinom_grid_loss_bounds(beta_grid, outcome$levels,
                                         outcome$reference)
  } else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "multinomial candidate loss bounds"),
    error = function(error) numeric())
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "multinomial statistic maxima"),
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
    isTRUE(beta_valid) && length(beta_grid) <= 256L && isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$candidate_order, "canonical_beta_grid_softmax_v1") &&
    length(observed_loss_bounds) == length(beta_grid) &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, length(beta_grid),
                          length(beta_grid)) &&
    identical(artifact$coordinate_order, paste(
      "canonical_beta_grid_multinomial_softmax_negative_log_likelihood_v1",
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
      "one_bounded_patient_multinomial_softmax_negative_log_likelihood",
      "contribution_for_every_signed_candidate_v1", sep = "_")) &&
    identical(maximum, expected_maximum) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l1_sensitivity), raw_l1,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l2_sensitivity), raw_l2,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l1_sensitivity),
                     raw_l1 / scale, tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity),
                     raw_l2 / scale, tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis, paste(
      "one_patient_changes_one_candidate_loss_by_at_most_its_signed",
      "multinomial_softmax_loss_bound_v1", sep = "_")) &&
    identical(artifact$estimation_scope, paste(
      "bounded_multinomial_fixed_covariates_finite_signed_beta_grid_v1",
      sep = "_"))
  if (!isTRUE(valid)) stop("The signed multinomial grid descriptor is invalid",
                           call. = FALSE)
  artifact$outcome <- outcome
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$beta_grid <- beta_grid
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$coordinate_count <- as.integer(length(beta_grid))
  artifact$statistic_maximum <- expected_maximum
  artifact
}

.dsvert_dp_multinom_grid_moment <- function(coordinates, artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released multinomial grid violates its signed bounds",
         call. = FALSE)
  }
  candidate <- which.min(coordinates)[[1L]]
  beta <- artifact$beta_grid[[candidate]]
  non_reference <- setdiff(artifact$outcome$levels, artifact$outcome$reference)
  terms <- artifact$design_terms
  normalized <- matrix(beta, nrow = length(terms), ncol = length(non_reference),
                       dimnames = list(terms, non_reference))
  ranges <- vapply(artifact$predictors, function(bound) {
    bound$upper - bound$lower
  }, numeric(1L))
  lowers <- vapply(artifact$predictors, `[[`, numeric(1L), "lower")
  coefficients <- normalized
  coefficients[-1L, ] <- sweep(normalized[-1L, , drop = FALSE], 1L, ranges, "/")
  coefficients[1L, ] <- normalized[1L, ] - colSums(
    sweep(coefficients[-1L, , drop = FALSE], 1L, lowers, "*"))
  list(status = "ok", coefficients = coefficients,
       normalized_coefficients = normalized, reference = artifact$outcome$reference,
       outcome_levels = artifact$outcome$levels,
       selected_candidate = as.integer(candidate),
       selected_dp_negative_log_likelihood = coordinates[[candidate]] /
         (2^as.numeric(artifact$numeric_grid_bits)),
       candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1")
}

.dsvert_dp_multinom_grid_synopsis_release <- function(
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
  artifact <- .dsvert_dp_multinom_grid_artifact(
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
    stop("The signed multinomial grid artifact does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  moment <- .dsvert_dp_multinom_grid_moment(coordinates, artifact)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !identical(verification$artifact$version,
                 .DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION)) {
    stop("The multinomial Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(context = context, metadata = metadata, artifact = verification$artifact,
       block = blocks[[1L]], coordinates = verification$coordinates,
       moment = verification$validated_moment %||% moment,
       certificate = certificate, verification = verification,
       scale = scale, capacity = capacity)
}

.dsvert_dp_multinom_grid_impl <- function(
    formula, data_name, analysis_id, server = NULL, datasources = NULL,
    .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(analysis_id, "analysis_id")
  outcome <- as.character(formula[[2L]])
  predictors <- attr(stats::terms(formula), "term.labels")
  released <- .dsvert_dp_multinom_grid_synopsis_release(
    data_name, analysis_id, server, datasources, .aggregate)
  artifact <- released$artifact
  if (!identical(artifact$outcome$column, outcome) ||
      !setequal(artifact$predictor_order, predictors)) {
    stop("formula must match the signed multinomial grid artifact",
         call. = FALSE)
  }
  moment <- released$moment
  result <- c(released$metadata, list(
    status = moment$status, analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    family = "multinomial_finite_grid", estimand = artifact$estimation_scope,
    coefficients = moment$coefficients,
    normalized_coefficients = moment$normalized_coefficients,
    reference = moment$reference, outcome_levels = moment$outcome_levels,
    selected_candidate = moment$selected_candidate,
    selected_dp_negative_log_likelihood = moment$selected_dp_negative_log_likelihood,
    candidate_selection = moment$candidate_selection,
    covariance = NULL, std_errors = NULL, standard_errors = NULL,
    p_values = NULL, inference = "unavailable_for_finite_grid_multinomial",
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE, provenance_certificate = released$certificate,
    disclosure_guard = list(satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("dsvert_dp_multinom_grid", "ds.vertMultinom", "list")
  result
}
