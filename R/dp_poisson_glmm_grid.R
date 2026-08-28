# Internal validator for a signed finite Poisson random-intercept GLMM grid.
# It consumes only public DP coordinates and does not expose an optimizer.

.DSVERT_CLIENT_DP_POISSON_GLMM_GRID_ARTIFACT_VERSION <-
  "bounded-poisson-random-intercept-likelihood-grid-v1"

.dsvert_dp_poisson_glmm_grid_loss_bounds <- function(
    beta_grid, variance_grid, max_outcome) {
  if (!is.list(beta_grid) || !length(beta_grid) ||
      !is.numeric(variance_grid) || !length(variance_grid) ||
      !.dsvert_dp_is_integer(max_outcome, 1L, 1024L)) {
    stop("The signed Poisson GLMM loss bound is invalid", call. = FALSE)
  }
  unlist(lapply(variance_grid, function(variance) {
    vapply(beta_grid, function(beta) {
      eta_bound <- sum(abs(beta)) + sqrt(2 * variance) *
        .DSVERT_CLIENT_DP_GLMM_GRID_MAX_GH_NODE
      exp(eta_bound) + max_outcome * eta_bound + lgamma(max_outcome + 1)
    }, numeric(1L))
  }), use.names = FALSE)
}

.dsvert_dp_poisson_glmm_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "max_outcome", "cluster", "predictors", "predictor_order",
    "intercept", "design_terms", "observation_capacity",
    "max_patients_per_cluster", "beta_grid", "variance_grid",
    "quadrature_rule", "candidate_order", "candidate_loss_bounds",
    "numeric_grid_bits", "coordinate_count", "coordinate_order",
    "source_coordinate_scaling", "repeated_record_policy",
    "missingness_policy", "contribution_domain", "statistic_maximum",
    "source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
    "natural_l1_sensitivity", "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis", "estimation_scope",
    "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version,
              .DSVERT_CLIENT_DP_POISSON_GLMM_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "poisson_random_intercept_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid Poisson random-intercept likelihood",
      " grid artifact '", analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "Poisson GLMM outcome")
  max_outcome <- artifact$max_outcome
  cluster <- artifact$cluster
  levels <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
    cluster$levels, "Poisson GLMM cluster levels"),
    error = function(error) character())
  cluster_valid <- .dsvert_dp_has_exact_names(cluster, c("column", "levels")) &&
    .dsvert_dp_is_string(cluster$column) && length(levels) >= 2L &&
    !anyDuplicated(levels) && all(nzchar(trimws(levels))) &&
    !identical(cluster$column, outcome$column)
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "Poisson GLMM predictor order", sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "Poisson GLMM design terms", sorted = FALSE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order) &&
    !outcome$column %in% predictor_order &&
    !cluster$column %in% predictor_order
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) .dsvert_dp_gaussian_bound(
        predictors[[variable]], "Poisson GLMM predictor")), predictor_order),
      error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(
      predictors, function(value) !identical(value$column, outcome$column) &&
        !identical(value$column, cluster$column), logical(1L)))
  }
  beta_grid <- artifact$beta_grid
  if (!is.list(beta_grid) || !length(beta_grid) || !is.null(names(beta_grid))) {
    beta_grid <- list()
  } else {
    beta_grid <- lapply(beta_grid, function(beta) tryCatch(
      .dsvert_dp_capsule_manifest_number_array(
        beta, "Poisson GLMM beta grid row"),
      error = function(error) numeric()))
  }
  beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
    .dsvert_joint_dp_client_json(as.list(beta))
  }, character(1L)) else character()
  variance_grid <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
    artifact$variance_grid, "Poisson GLMM variance grid"),
    error = function(error) numeric())
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(
    as.numeric(artifact$max_patients_per_cluster))
  dimension <- 1L + length(predictor_order)
  beta_valid <- length(beta_grid) && all(vapply(beta_grid, function(beta) {
    length(beta) == dimension && !anyNA(beta) && all(is.finite(beta)) &&
      all(abs(beta) <= 8) && sum(abs(beta)) <= 16
  }, logical(1L))) && !anyDuplicated(beta_keys) &&
    identical(beta_keys, sort(beta_keys, method = "radix"))
  variance_valid <- length(variance_grid) && !anyNA(variance_grid) &&
    all(is.finite(variance_grid)) && all(variance_grid >= 0) &&
    all(variance_grid <= 16) && all(diff(variance_grid) > 0)
  max_outcome_valid <- .dsvert_dp_is_integer(max_outcome, 1L, 1024L) &&
    isTRUE(all.equal(outcome$lower, 0, tolerance = 0)) &&
    isTRUE(all.equal(outcome$upper, as.numeric(max_outcome), tolerance = 0))
  candidate_count <- length(beta_grid) * length(variance_grid)
  loss_bounds <- if (isTRUE(beta_valid) && isTRUE(variance_valid) &&
      isTRUE(max_outcome_valid)) {
    tryCatch(.dsvert_dp_poisson_glmm_grid_loss_bounds(
      beta_grid, variance_grid, as.integer(max_outcome)),
      error = function(error) numeric())
  } else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "Poisson GLMM candidate loss bounds"),
    error = function(error) numeric())
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "Poisson GLMM statistic maxima"),
    error = function(error) numeric())
  raw_per_candidate <- ceiling(loss_bounds * scale)
  expected_maximum <- capacity * raw_per_candidate
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) {
    2
  } else if (identical(adjacency, "add_remove_patient")) {
    1
  } else {
    NA_real_
  }
  raw_l1 <- multiplier * sum(raw_per_candidate)
  raw_l2 <- multiplier * sqrt(sum(raw_per_candidate^2))
  valid <- cluster_valid && isTRUE(predictors_valid) &&
    isTRUE(max_outcome_valid) && .dsvert_dp_is_integer(bits, 8L, 18L) &&
    identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    .dsvert_dp_is_integer(cluster_capacity, 2L, capacity) &&
    isTRUE(beta_valid) && isTRUE(variance_valid) && candidate_count <= 256L &&
    isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$quadrature_rule, "gauss_hermite_9_standard_normal_v1") &&
    identical(artifact$candidate_order, "variance_grid_then_beta_grid_v1") &&
    length(observed_loss_bounds) == candidate_count &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count,
                          candidate_count, candidate_count) &&
    identical(artifact$coordinate_order, paste(
      "variance_grid_then_beta_grid_cluster_marginal_negative_log",
      "likelihood_v1", sep = "_")) &&
    identical(artifact$source_coordinate_scaling,
              "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "require_one_bounded_poisson_outcome_and_mean_once_per_admitted",
      "patient_with_one_consistent_public_cluster_level_v1", sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "noninteger_or_out_of_range_or_missing_outcome_or_missing_or",
      "nonfinite_predictor_or_missing_or_inconsistent_cluster_excludes",
      "patient_v1", sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient_poisson_log_likelihood_contribution_in_one",
      "consistent_cluster_for_every_signed_candidate_v1", sep = "_")) &&
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
      "one_patient_changes_one_cluster_marginal_log_likelihood_by_at",
      "most_its_signed_poisson_loss_bound_v1", sep = "_")) &&
    identical(artifact$estimation_scope, paste(
      "bounded_poisson_random_intercept_marginal_likelihood_fixed",
      "covariates_finite_signed_parameter_grid_v1", sep = "_"))
  if (!isTRUE(valid)) {
    stop("The signed Poisson random-intercept GLMM descriptor is invalid",
         call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$max_outcome <- as.integer(max_outcome)
  artifact$cluster <- list(column = enc2utf8(cluster$column), levels = levels)
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$beta_grid <- beta_grid
  artifact$variance_grid <- variance_grid
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$coordinate_count <- as.integer(candidate_count)
  artifact$statistic_maximum <- expected_maximum
  artifact
}
