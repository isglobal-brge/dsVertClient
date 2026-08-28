# Internal validator for a signed finite Poisson random-intercept GLMM grid.
# It consumes only public DP coordinates and does not expose an optimizer.

.DSVERT_CLIENT_DP_POISSON_GLMM_GRID_ARTIFACT_VERSION <-
  "bounded-poisson-random-intercept-likelihood-grid-v1"
.DSVERT_CLIENT_DP_POISSON_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION <-
  "bounded-poisson-random-slope-likelihood-grid-v1"

.dsvert_dp_poisson_glmm_random_slope_candidates <- function(
    candidate_grid, dimension, random_effect_order, max_outcome) {
  if (!is.list(candidate_grid) || !length(candidate_grid) ||
      !is.null(names(candidate_grid)) ||
      !is.character(random_effect_order) || length(random_effect_order) < 2L ||
      length(random_effect_order) > 4L ||
      !identical(random_effect_order[[1L]], "(Intercept)") ||
      !.dsvert_dp_is_integer(max_outcome, 1L, 1024L)) return(list())
  candidates <- lapply(candidate_grid, function(candidate) {
    if (!.dsvert_dp_has_exact_names(candidate, c("beta", "covariance"))) {
      return(NULL)
    }
    beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$beta, "Poisson GLMM random-slope beta"),
      error = function(error) numeric())
    covariance <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$covariance, "Poisson GLMM random-slope covariance"),
      error = function(error) numeric())
    effect_count <- length(random_effect_order)
    if (length(beta) != dimension || anyNA(beta) || any(!is.finite(beta)) ||
        any(abs(beta) > 8) || sum(abs(beta)) > 16 ||
        length(covariance) != effect_count^2 || anyNA(covariance) ||
        any(!is.finite(covariance)) || any(abs(covariance) > 16)) {
      return(NULL)
    }
    covariance <- matrix(covariance, effect_count, effect_count, byrow = TRUE)
    if (!isTRUE(all.equal(covariance, t(covariance), tolerance = 0)) ||
        any(diag(covariance) < 0) ||
        any(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values < -1e-10)) {
      return(NULL)
    }
    decomposition <- eigen(covariance, symmetric = TRUE)
    root <- decomposition$vectors %*% diag(sqrt(pmax(0, decomposition$values)))
    eta_bound <- sum(abs(beta)) + sqrt(2) *
      .DSVERT_CLIENT_DP_GLMM_GRID_MAX_GH_NODE * sum(abs(root))
    loss_bound <- exp(eta_bound) + max_outcome * eta_bound +
      lgamma(max_outcome + 1)
    if (!is.finite(loss_bound) || loss_bound <= 0) return(NULL)
    list(beta = beta, covariance = covariance, loss_bound = loss_bound)
  })
  if (any(vapply(candidates, is.null, logical(1L)))) return(list())
  keys <- vapply(candidates, function(candidate) .dsvert_joint_dp_client_json(list(
    beta = as.list(candidate$beta),
    covariance = as.list(as.vector(t(candidate$covariance))))), character(1L))
  if (anyDuplicated(keys) || !identical(keys, sort(keys, method = "radix"))) {
    return(list())
  }
  candidates
}

.dsvert_dp_poisson_glmm_random_slope_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "max_outcome", "cluster", "predictors", "predictor_order",
    "intercept", "design_terms", "random_effect_order",
    "observation_capacity", "max_patients_per_cluster", "candidate_grid",
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
              .DSVERT_CLIENT_DP_POISSON_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "poisson_random_slope_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) && .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid Poisson random-slope GLMM grid artifact '",
      analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "Poisson GLMM outcome")
  max_outcome <- artifact$max_outcome
  cluster <- artifact$cluster
  levels <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
    cluster$levels, "Poisson GLMM random-slope cluster levels"),
    error = function(error) character())
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "Poisson GLMM random-slope predictor order",
    sorted = TRUE), error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "Poisson GLMM random-slope design terms",
    sorted = FALSE), error = function(error) character())
  effects <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$random_effect_order, "Poisson GLMM random-slope effects",
    sorted = FALSE), error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order)
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(predictor_order, function(x) {
      .dsvert_dp_gaussian_bound(predictors[[x]],
                                "Poisson GLMM random-slope predictor")
    }), predictor_order), error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(predictors, function(x) {
      !identical(x$column, outcome$column) && !identical(x$column, cluster$column)
    }, logical(1L)))
  }
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(as.numeric(artifact$max_patients_per_cluster))
  max_outcome_valid <- .dsvert_dp_is_integer(max_outcome, 1L, 1024L) &&
    isTRUE(all.equal(outcome$lower, 0, tolerance = 0)) &&
    isTRUE(all.equal(outcome$upper, as.numeric(max_outcome), tolerance = 0))
  effects_valid <- length(effects) >= 2L && length(effects) <= 4L &&
    identical(effects[[1L]], "(Intercept)") &&
    all(effects[-1L] %in% predictor_order) &&
    !anyDuplicated(effects)
  candidates <- if (isTRUE(effects_valid) && isTRUE(max_outcome_valid)) {
    .dsvert_dp_poisson_glmm_random_slope_candidates(
      artifact$candidate_grid, 1L + length(predictor_order), effects,
      as.integer(max_outcome))
  } else list()
  loss_bounds <- if (length(candidates)) {
    vapply(candidates, `[[`, numeric(1L), "loss_bound")
  } else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "Poisson GLMM random-slope candidate loss bounds"),
    error = function(error) numeric())
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "Poisson GLMM random-slope maxima"),
    error = function(error) numeric())
  raw <- ceiling(loss_bounds * scale)
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) 2 else if (
    identical(adjacency, "add_remove_patient")) 1 else NA_real_
  estimation_scope <- if (length(effects) == 2L) {
    "bounded_poisson_random_intercept_and_one_random_slope_marginal_likelihood_finite_signed_parameter_grid_v1"
  } else {
    "bounded_poisson_random_intercept_and_one_to_three_random_slopes_marginal_likelihood_finite_signed_parameter_grid_v1"
  }
  valid <- .dsvert_dp_has_exact_names(cluster, c("column", "levels")) &&
    .dsvert_dp_is_string(cluster$column) && length(levels) >= 2L &&
    !anyDuplicated(levels) && all(nzchar(trimws(levels))) &&
    isTRUE(predictors_valid) && isTRUE(max_outcome_valid) &&
    .dsvert_dp_is_integer(bits, 8L, 18L) && identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    .dsvert_dp_is_integer(cluster_capacity, 2L, capacity) &&
    isTRUE(effects_valid) && length(candidates) && length(candidates) <= 64L &&
    isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$quadrature_rule,
              .dsvert_dp_glmm_random_slope_quadrature_rule_v1(length(effects))) &&
    identical(artifact$candidate_order, "canonical_signed_candidate_grid_v1") &&
    length(observed_loss_bounds) == length(candidates) &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count,
                          length(candidates), length(candidates)) &&
    identical(artifact$coordinate_order,
      "signed_candidate_grid_cluster_marginal_poisson_negative_log_likelihood_v1") &&
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
    identical(maximum, capacity * raw) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l1_sensitivity),
                     multiplier * sum(raw), tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l2_sensitivity),
                     multiplier * sqrt(sum(raw^2)), tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l1_sensitivity),
                     multiplier * sum(raw) / scale, tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity),
                     multiplier * sqrt(sum(raw^2)) / scale, tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis, paste(
      "one_patient_changes_one_cluster_marginal_log_likelihood_by_at",
      "most_its_signed_poisson_loss_bound_v1", sep = "_")) &&
    identical(artifact$estimation_scope, estimation_scope)
  if (!isTRUE(valid)) {
    stop("The signed Poisson random-slope GLMM descriptor is invalid",
         call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$max_outcome <- as.integer(max_outcome)
  artifact$cluster <- list(column = enc2utf8(cluster$column), levels = levels)
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$random_effect_order <- effects
  artifact$candidate_grid <- lapply(candidates, function(candidate) list(
    beta = candidate$beta, covariance = as.vector(t(candidate$covariance))))
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$statistic_maximum <- capacity * raw
  artifact$coordinate_count <- as.integer(length(candidates))
  artifact
}

.dsvert_dp_poisson_glmm_random_slope_grid_moment <- function(coordinates,
                                                               artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released Poisson random-slope GLMM grid violates its signed bounds",
         call. = FALSE)
  }
  selected <- which.min(coordinates)[[1L]]
  candidate <- .dsvert_dp_poisson_glmm_random_slope_candidates(
    artifact$candidate_grid, 1L + length(artifact$predictor_order),
    artifact$random_effect_order, artifact$max_outcome)[[selected]]
  slopes <- candidate$beta[-1L] / vapply(
    artifact$predictors, function(x) x$upper - x$lower, numeric(1L))
  names(slopes) <- artifact$predictor_order
  intercept <- candidate$beta[[1L]] - sum(slopes * vapply(
    artifact$predictors, `[[`, numeric(1L), "lower"))
  random_slopes <- artifact$random_effect_order[-1L]
  spans <- vapply(random_slopes, function(slope) {
    artifact$predictors[[slope]]$upper - artifact$predictors[[slope]]$lower
  }, numeric(1L))
  lowers <- vapply(random_slopes, function(slope) {
    artifact$predictors[[slope]]$lower
  }, numeric(1L))
  effect_count <- length(artifact$random_effect_order)
  transform <- diag(effect_count)
  transform[1L, -1L] <- -lowers / spans
  transform[cbind(seq.int(2L, effect_count), seq.int(2L, effect_count))] <-
    1 / spans
  covariance <- transform %*% candidate$covariance %*% t(transform)
  dimnames(covariance) <- list(artifact$random_effect_order,
                               artifact$random_effect_order)
  list(
    status = "ok", coefficients = c(`(Intercept)` = intercept, slopes),
    normalized_coefficients = stats::setNames(candidate$beta, artifact$design_terms),
    random_effect_covariance = covariance,
    random_effect_order = artifact$random_effect_order,
    selected_candidate = as.integer(selected),
    selected_dp_negative_log_likelihood = coordinates[[selected]] /
      (2^as.numeric(artifact$numeric_grid_bits)),
    candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1")
}

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
