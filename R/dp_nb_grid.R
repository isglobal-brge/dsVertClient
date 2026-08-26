# Internal validator and post-processing for a custodian-signed finite NB2
# likelihood grid.  It consumes only public DP coordinates from one sticky
# Synopsis; no profile-score RPC or client-side optimisation is performed.

.DSVERT_CLIENT_DP_NB_GRID_ARTIFACT_VERSION <-
  "bounded-negative-binomial-likelihood-grid-v1"

.dsvert_dp_nb_grid_loss_bounds <- function(beta_grid, theta_grid, max_outcome) {
  log1pexp <- function(value) pmax(value, 0) + log1p(exp(-abs(value)))
  unlist(lapply(theta_grid, function(theta) vapply(beta_grid, function(beta) {
    eta_bound <- sum(abs(beta))
    y <- 0:max_outcome
    max(0, max(outer(y, c(-eta_bound, eta_bound), function(count, eta) {
      lgamma(theta) + lgamma(count + 1) - lgamma(count + theta) +
        theta * log1pexp(eta - log(theta)) - count * eta
    })))
  }, numeric(1L))), use.names = FALSE)
}

.dsvert_dp_nb_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "predictors", "predictor_order", "intercept",
    "design_terms", "observation_capacity", "max_outcome", "beta_grid",
    "theta_grid", "candidate_order", "candidate_loss_bounds",
    "numeric_grid_bits", "coordinate_count", "coordinate_order",
    "source_coordinate_scaling", "repeated_record_policy",
    "missingness_policy", "contribution_domain", "statistic_maximum",
    "source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
    "natural_l1_sensitivity", "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis", "estimation_scope",
    "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version, .DSVERT_CLIENT_DP_NB_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "negative_binomial_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) && .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid negative-binomial likelihood grid '",
      analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "NB2 outcome")
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "NB2 predictor order", sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "NB2 design terms", sorted = FALSE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order) &&
    !outcome$column %in% predictor_order
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) .dsvert_dp_gaussian_bound(
        predictors[[variable]], "NB2 predictor")), predictor_order),
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
      .dsvert_dp_capsule_manifest_number_array(beta, "NB2 beta grid row"),
      error = function(error) numeric()))
  }
  beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
    .dsvert_joint_dp_client_json(as.list(beta))
  }, character(1L)) else character()
  theta_grid <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
    artifact$theta_grid, "NB2 theta grid"), error = function(error) numeric())
  max_outcome <- suppressWarnings(as.numeric(artifact$max_outcome))
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  dimension <- 1L + length(predictor_order)
  beta_valid <- length(beta_grid) && all(vapply(beta_grid, function(beta) {
    length(beta) == dimension && !anyNA(beta) && all(is.finite(beta)) &&
      all(abs(beta) <= 8)
  }, logical(1L))) && !anyDuplicated(beta_keys) &&
    identical(beta_keys, sort(beta_keys, method = "radix"))
  theta_valid <- length(theta_grid) && !anyNA(theta_grid) &&
    all(is.finite(theta_grid)) && all(theta_grid > 0) &&
    all(theta_grid <= 64) && all(diff(theta_grid) > 0)
  candidate_count <- length(beta_grid) * length(theta_grid)
  loss_bounds <- if (isTRUE(beta_valid) && isTRUE(theta_valid) &&
      .dsvert_dp_is_integer(max_outcome, 1L, 1024L)) {
    .dsvert_dp_nb_grid_loss_bounds(beta_grid, theta_grid, as.integer(max_outcome))
  } else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "NB2 candidate loss bounds"),
    error = function(error) numeric())
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "NB2 statistic maxima"),
    error = function(error) numeric())
  raw_per_candidate <- ceiling(loss_bounds * scale)
  expected_maximum <- capacity * raw_per_candidate
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) 2 else if (
    identical(adjacency, "add_remove_patient")) 1 else NA_real_
  raw_l1 <- multiplier * sum(raw_per_candidate)
  raw_l2 <- multiplier * sqrt(sum(raw_per_candidate^2))
  valid <- isTRUE(predictors_valid) && isTRUE(all.equal(outcome$lower, 0)) &&
    .dsvert_dp_is_integer(max_outcome, 1L, 1024L) &&
    isTRUE(all.equal(outcome$upper, max_outcome)) &&
    .dsvert_dp_is_integer(bits, 8L, 18L) && identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    isTRUE(beta_valid) && isTRUE(theta_valid) && candidate_count <= 256L &&
    isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$candidate_order, "theta_grid_then_beta_grid_v1") &&
    length(observed_loss_bounds) == candidate_count &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, candidate_count,
                          candidate_count) &&
    identical(artifact$coordinate_order, paste(
      "theta_grid_then_beta_grid_negative_binomial_log_likelihood_v1",
      sep = "_")) &&
    identical(artifact$source_coordinate_scaling,
              "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "require_one_bounded_count_outcome_and_mean_once_per_admitted",
      "patient_v1", sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "noninteger_or_out_of_range_or_missing_outcome_or_missing_or",
      "nonfinite_predictor_excludes_patient_v1", sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient_negative_binomial_log_likelihood",
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
      "negative_binomial_loss_bound_v1", sep = "_")) &&
    identical(artifact$estimation_scope, paste(
      "bounded_negative_binomial_fixed_covariates_finite_signed",
      "beta_theta_grid_v1", sep = "_"))
  if (!isTRUE(valid)) stop("The signed negative-binomial grid descriptor is invalid",
                           call. = FALSE)
  artifact$outcome <- outcome
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$beta_grid <- beta_grid
  artifact$theta_grid <- theta_grid
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$coordinate_count <- as.integer(candidate_count)
  artifact$statistic_maximum <- expected_maximum
  artifact
}

.dsvert_dp_nb_grid_moment <- function(coordinates, artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released negative-binomial grid violates its signed bounds",
         call. = FALSE)
  }
  candidate <- which.min(coordinates)[[1L]]
  beta_count <- length(artifact$beta_grid)
  theta_index <- as.integer((candidate - 1L) %/% beta_count) + 1L
  beta_index <- as.integer((candidate - 1L) %% beta_count) + 1L
  beta <- artifact$beta_grid[[beta_index]]
  ranges <- vapply(artifact$predictors, function(bound) {
    bound$upper - bound$lower
  }, numeric(1L))
  lowers <- vapply(artifact$predictors, `[[`, numeric(1L), "lower")
  slopes <- beta[-1L] / ranges
  names(slopes) <- artifact$predictor_order
  intercept <- beta[[1L]] - sum(slopes * lowers)
  list(status = "ok", coefficients = c(`(Intercept)` = intercept, slopes),
       normalized_coefficients = stats::setNames(beta, artifact$design_terms),
       theta = artifact$theta_grid[[theta_index]],
       selected_candidate = as.integer(candidate),
       selected_beta_index = as.integer(beta_index),
       selected_theta_index = as.integer(theta_index),
       selected_dp_negative_log_likelihood = coordinates[[candidate]] /
         (2^as.numeric(artifact$numeric_grid_bits)),
       candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1")
}

.dsvert_dp_nb_grid_synopsis_release <- function(
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
  artifact <- .dsvert_dp_nb_grid_artifact(
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
    stop("The signed negative-binomial grid artifact does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  moment <- .dsvert_dp_nb_grid_moment(coordinates, artifact)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !identical(verification$artifact$version,
                 .DSVERT_CLIENT_DP_NB_GRID_ARTIFACT_VERSION)) {
    stop("The negative-binomial Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(context = context, metadata = metadata, artifact = verification$artifact,
       block = blocks[[1L]], coordinates = verification$coordinates,
       moment = verification$validated_moment %||% moment,
       certificate = certificate, verification = verification,
       scale = scale, capacity = capacity)
}

.dsvert_dp_nb_grid_impl <- function(
    formula, data_name, analysis_id, server = NULL, datasources = NULL,
    .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(analysis_id, "analysis_id")
  outcome <- as.character(formula[[2L]])
  predictors <- attr(stats::terms(formula), "term.labels")
  released <- .dsvert_dp_nb_grid_synopsis_release(
    data_name, analysis_id, server, datasources, .aggregate)
  artifact <- released$artifact
  if (!identical(artifact$outcome$column, outcome) ||
      !setequal(artifact$predictor_order, predictors)) {
    stop("formula must match the signed negative-binomial grid artifact",
         call. = FALSE)
  }
  moment <- released$moment
  result <- c(released$metadata, list(
    status = moment$status, analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    family = "negative_binomial_finite_grid", estimand = artifact$estimation_scope,
    coefficients = moment$coefficients,
    normalized_coefficients = moment$normalized_coefficients,
    theta = moment$theta, selected_candidate = moment$selected_candidate,
    selected_beta_index = moment$selected_beta_index,
    selected_theta_index = moment$selected_theta_index,
    selected_dp_negative_log_likelihood = moment$selected_dp_negative_log_likelihood,
    candidate_selection = moment$candidate_selection,
    covariance = NULL, std_errors = NULL, standard_errors = NULL,
    p_values = NULL, inference = "unavailable_for_finite_grid_nb2",
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE, provenance_certificate = released$certificate,
    disclosure_guard = list(satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("dsvert_dp_nb2_grid", "ds.vertNBFullRegTheta", "list")
  result
}
