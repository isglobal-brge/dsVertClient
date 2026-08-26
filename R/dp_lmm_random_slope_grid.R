# Internal validation and deterministic post-processing for the signed
# finite-grid Gaussian random-slope LMM artifact.  The artifact contains only
# public DP-clipped cluster marginal losses; it never contains a per-cluster
# design, residual, covariance estimate, or fitted random effect.

.DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION <-
  "bounded-gaussian-random-slope-likelihood-grid-v1"

.dsvert_dp_lmm_random_slope_candidates <- function(candidate_grid,
                                                    dimension,
                                                    random_effect_order,
                                                    cluster_capacity) {
  effect_count <- length(random_effect_order)
  if (!is.list(candidate_grid) || !length(candidate_grid) ||
      !is.null(names(candidate_grid))) return(list())
  candidates <- lapply(candidate_grid, function(candidate) {
    required <- c("beta", "sigma2", "covariance")
    if (!.dsvert_dp_has_exact_names(candidate, required)) return(NULL)
    beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$beta, "LMM random-slope beta"), error = function(error) numeric())
    covariance <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$covariance, "LMM random-slope covariance"),
      error = function(error) numeric())
    sigma2 <- suppressWarnings(as.numeric(candidate$sigma2))
    if (length(beta) != dimension || anyNA(beta) || any(!is.finite(beta)) ||
        any(abs(beta) > 8) || length(sigma2) != 1L || !is.finite(sigma2) ||
        sigma2 < 0.25 || sigma2 > 16 ||
        length(covariance) != effect_count^2 || anyNA(covariance) ||
        any(!is.finite(covariance)) || any(abs(covariance) > 16)) return(NULL)
    covariance <- matrix(covariance, effect_count, effect_count, byrow = TRUE)
    if (!isTRUE(all.equal(covariance, t(covariance), tolerance = 0)) ||
        any(diag(covariance) < 0) ||
        any(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values <
            -1e-10)) return(NULL)
    residual_bound <- 1 + sum(abs(beta))
    variance_upper <- sigma2 + 16 * cluster_capacity * effect_count^2
    loss_bound <- 0.5 * cluster_capacity * (
      log(2 * pi) + log(variance_upper) + residual_bound^2 / sigma2)
    if (!is.finite(loss_bound) || loss_bound <= 0) return(NULL)
    list(beta = beta, sigma2 = sigma2, covariance = covariance,
         loss_bound = loss_bound)
  })
  if (any(vapply(candidates, is.null, logical(1L)))) return(list())
  keys <- vapply(candidates, function(candidate) .dsvert_joint_dp_client_json(
    list(beta = as.list(candidate$beta), sigma2 = candidate$sigma2,
         covariance = as.list(as.vector(t(candidate$covariance))))), character(1L))
  if (anyDuplicated(keys) || !identical(keys, sort(keys, method = "radix"))) {
    return(list())
  }
  candidates
}

.dsvert_dp_lmm_random_slope_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "cluster", "predictors", "predictor_order", "intercept",
    "design_terms", "random_effect_order", "observation_capacity",
    "max_patients_per_cluster", "candidate_grid", "candidate_order",
    "candidate_loss_bounds", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "statistic_maximum", "source_raw_l1_sensitivity",
    "source_raw_l2_sensitivity", "natural_l1_sensitivity",
    "natural_l2_sensitivity", "adjacency", "adjacency_sensitivity_basis",
    "estimation_scope", "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version,
              .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "gaussian_random_slope_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) && .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) .dsvert_dp_gaussian_reserved(paste0(
    "the signed capsule has no valid Gaussian random-slope LMM grid artifact '",
    analysis_id, "' for dataset '", data_name, "'"))
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "LMM outcome")
  cluster <- artifact$cluster
  levels <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
    cluster$levels, "LMM random-slope cluster levels"),
    error = function(error) character())
  cluster_valid <- .dsvert_dp_has_exact_names(cluster, c("column", "levels")) &&
    .dsvert_dp_is_string(cluster$column) && length(levels) >= 2L &&
    !anyDuplicated(levels) && all(nzchar(trimws(levels))) &&
    !identical(cluster$column, outcome$column)
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "LMM random-slope predictor order", sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "LMM random-slope design terms", sorted = FALSE),
    error = function(error) character())
  random_effect_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$random_effect_order, "LMM random-slope effects", sorted = FALSE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order)
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(predictor_order, function(x) {
      .dsvert_dp_gaussian_bound(predictors[[x]], "LMM random-slope predictor")
    }), predictor_order), error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(
      predictors, function(x) !identical(x$column, outcome$column) &&
        !identical(x$column, cluster$column), logical(1L)))
  }
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(as.numeric(artifact$max_patients_per_cluster))
  dimension <- 1L + length(predictor_order)
  effects_valid <- length(random_effect_order) >= 2L &&
    identical(random_effect_order[[1L]], "(Intercept)") &&
    all(random_effect_order[-1L] %in% predictor_order) &&
    !anyDuplicated(random_effect_order)
  candidates <- if (isTRUE(effects_valid) && is.finite(cluster_capacity)) {
    .dsvert_dp_lmm_random_slope_candidates(
      artifact$candidate_grid, dimension, random_effect_order, cluster_capacity)
  } else list()
  loss_bounds <- if (length(candidates)) vapply(candidates, `[[`, numeric(1L), "loss_bound") else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "LMM random-slope candidate loss bounds"),
    error = function(error) numeric())
  raw_per_candidate <- ceiling(loss_bounds * scale)
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "LMM random-slope statistic maxima"),
    error = function(error) numeric())
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) 2 else if (
    identical(adjacency, "add_remove_patient")) 1 else NA_real_
  raw_l1 <- multiplier * sum(raw_per_candidate)
  raw_l2 <- multiplier * sqrt(sum(raw_per_candidate^2))
  valid <- cluster_valid && isTRUE(predictors_valid) &&
    .dsvert_dp_is_integer(bits, 8L, 18L) && identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    .dsvert_dp_is_integer(cluster_capacity, 2L, capacity) &&
    isTRUE(effects_valid) && length(candidates) && length(candidates) <= 128L &&
    isTRUE(artifact$intercept) && identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$candidate_order, "canonical_signed_candidate_grid_v1") &&
    length(observed_loss_bounds) == length(candidates) &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, length(candidates), length(candidates)) &&
    identical(artifact$coordinate_order,
              "signed_candidate_grid_clipped_cluster_gaussian_negative_log_likelihood_v1") &&
    identical(artifact$source_coordinate_scaling,
              "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy,
              "require_one_complete_bounded_row_per_admitted_patient_with_one_consistent_public_cluster_level_v1") &&
    identical(artifact$missingness_policy,
              "missing_or_nonfinite_outcome_predictor_or_missing_or_inconsistent_cluster_excludes_patient_v1") &&
    identical(artifact$contribution_domain,
              "one_bounded_patient_changes_one_clipped_cluster_gaussian_loss_per_signed_candidate_v1") &&
    identical(maximum, capacity * raw_per_candidate) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l1_sensitivity), raw_l1,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l2_sensitivity), raw_l2,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l1_sensitivity), raw_l1 / scale,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity), raw_l2 / scale,
                     tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis,
              "one_patient_can_change_one_entire_clipped_cluster_loss_by_at_most_its_signed_bound_v1") &&
    identical(artifact$estimation_scope,
              "bounded_gaussian_random_slope_marginal_likelihood_finite_signed_parameter_grid_v1")
  if (!isTRUE(valid)) stop("The signed Gaussian random-slope LMM descriptor is invalid",
                           call. = FALSE)
  artifact$outcome <- outcome
  artifact$cluster <- list(column = enc2utf8(cluster$column), levels = levels)
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$random_effect_order <- random_effect_order
  artifact$candidate_grid <- lapply(candidates, function(candidate) list(
    beta = candidate$beta, sigma2 = candidate$sigma2,
    covariance = as.vector(t(candidate$covariance))))
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$statistic_maximum <- capacity * raw_per_candidate
  artifact$coordinate_count <- as.integer(length(candidates))
  artifact
}

.dsvert_dp_lmm_random_slope_grid_moment <- function(coordinates, artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released Gaussian random-slope LMM grid violates its signed bounds",
         call. = FALSE)
  }
  selected <- which.min(coordinates)[[1L]]
  candidates <- .dsvert_dp_lmm_random_slope_candidates(
    artifact$candidate_grid, 1L + length(artifact$predictor_order),
    artifact$random_effect_order, artifact$max_patients_per_cluster)
  if (!length(candidates)) {
    stop("The Gaussian random-slope LMM candidate grid is invalid", call. = FALSE)
  }
  candidate <- candidates[[selected]]
  y_span <- artifact$outcome$upper - artifact$outcome$lower
  beta <- candidate$beta
  slopes <- beta[-1L] * y_span / vapply(
    artifact$predictors, function(x) x$upper - x$lower, numeric(1L))
  names(slopes) <- artifact$predictor_order
  intercept <- artifact$outcome$lower + y_span * beta[[1L]] -
    sum(slopes * vapply(artifact$predictors, `[[`, numeric(1L), "lower"))
  transform <- diag(length(artifact$random_effect_order))
  for (index in seq_along(artifact$random_effect_order)[-1L]) {
    variable <- artifact$random_effect_order[[index]]
    span <- artifact$predictors[[variable]]$upper -
      artifact$predictors[[variable]]$lower
    transform[1L, index] <- -artifact$predictors[[variable]]$lower / span
    transform[index, index] <- 1 / span
  }
  covariance <- y_span^2 * transform %*% candidate$covariance %*% t(transform)
  dimnames(covariance) <- list(artifact$random_effect_order,
                               artifact$random_effect_order)
  list(status = "ok", coefficients = c(`(Intercept)` = intercept, slopes),
       normalized_coefficients = stats::setNames(beta, artifact$design_terms),
       sigma2 = candidate$sigma2 * y_span^2,
       random_effect_covariance = covariance,
       random_effect_order = artifact$random_effect_order,
       selected_candidate = as.integer(selected),
       selected_dp_negative_log_likelihood = coordinates[[selected]] /
         (2^as.numeric(artifact$numeric_grid_bits)),
       candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1")
}
