# Internal validator and post-processing for the signed finite-grid binary
# random-intercept GLMM artifact.  It consumes only public DP coordinates.

.DSVERT_CLIENT_DP_GLMM_GRID_ARTIFACT_VERSION <-
  "bounded-binary-random-intercept-likelihood-grid-v1"
.DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION <-
  "bounded-binary-random-slope-likelihood-grid-v1"
.DSVERT_CLIENT_DP_GLMM_GRID_MAX_GH_NODE <- 3.19099320178153

.dsvert_dp_glmm_random_slope_quadrature_rule_v1 <- function(effect_count) {
  if (!is.numeric(effect_count) || length(effect_count) != 1L ||
      is.na(effect_count) || effect_count != floor(effect_count) ||
      effect_count < 2L || effect_count > 4L) return(NA_character_)
  paste0("gauss_hermite_", paste(rep.int("9", effect_count), collapse = "x"),
         "_standard_normal_v1")
}

.dsvert_dp_glmm_random_slope_candidates <- function(candidate_grid, dimension,
                                                     random_effect_order) {
  if (!is.list(candidate_grid) || !length(candidate_grid) ||
      !is.null(names(candidate_grid)) || length(random_effect_order) < 2L ||
      length(random_effect_order) > 4L ||
      !identical(random_effect_order[[1L]], "(Intercept)")) return(list())
  effect_count <- length(random_effect_order)
  candidates <- lapply(candidate_grid, function(candidate) {
    if (!.dsvert_dp_has_exact_names(candidate, c("beta", "covariance"))) return(NULL)
    beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$beta, "GLMM random-slope beta"), error = function(error) numeric())
    covariance <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$covariance, "GLMM random-slope covariance"),
      error = function(error) numeric())
    if (length(beta) != dimension || anyNA(beta) || any(!is.finite(beta)) ||
        any(abs(beta) > 8) || length(covariance) != effect_count^2 ||
        anyNA(covariance) ||
        any(!is.finite(covariance)) || any(abs(covariance) > 16)) return(NULL)
    covariance <- matrix(covariance, effect_count, effect_count, byrow = TRUE)
    if (!isTRUE(all.equal(covariance, t(covariance), tolerance = 0)) ||
        any(diag(covariance) < 0) ||
        any(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values < -1e-10)) {
      return(NULL)
    }
    decomposition <- eigen(covariance, symmetric = TRUE)
    root <- decomposition$vectors %*% diag(sqrt(pmax(0, decomposition$values)))
    loss_bound <- log1p(exp(sum(abs(beta)) + sqrt(2) *
      .DSVERT_CLIENT_DP_GLMM_GRID_MAX_GH_NODE * sum(abs(root))))
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

.dsvert_dp_glmm_random_slope_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "cluster", "predictors", "predictor_order", "intercept",
    "design_terms", "random_effect_order", "observation_capacity",
    "max_patients_per_cluster", "candidate_grid", "quadrature_rule",
    "candidate_order", "candidate_loss_bounds", "numeric_grid_bits",
    "coordinate_count", "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "statistic_maximum", "source_raw_l1_sensitivity",
    "source_raw_l2_sensitivity", "natural_l1_sensitivity",
    "natural_l2_sensitivity", "adjacency", "adjacency_sensitivity_basis",
    "estimation_scope", "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version,
              .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "binary_random_slope_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) && identical(artifact$dataset, data_name) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) .dsvert_dp_gaussian_reserved(paste0(
    "the signed capsule has no valid binary random-slope GLMM grid artifact '",
    analysis_id, "' for dataset '", data_name, "'"))
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "GLMM outcome")
  cluster <- artifact$cluster
  levels <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
    cluster$levels, "GLMM random-slope cluster levels"),
    error = function(error) character())
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "GLMM random-slope predictor order", sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "GLMM random-slope design terms", sorted = FALSE),
    error = function(error) character())
  effects <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$random_effect_order, "GLMM random-slope effects", sorted = FALSE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order)
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(predictor_order, function(x) {
      .dsvert_dp_gaussian_bound(predictors[[x]], "GLMM random-slope predictor")
    }), predictor_order), error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(predictors, function(x) {
      !identical(x$column, outcome$column) && !identical(x$column, cluster$column)
    }, logical(1L)))
  }
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(as.numeric(artifact$max_patients_per_cluster))
  effects_valid <- length(effects) >= 2L && length(effects) <= 4L &&
    identical(effects[[1L]], "(Intercept)") &&
    all(effects[-1L] %in% predictor_order) && !anyDuplicated(effects)
  candidates <- if (isTRUE(effects_valid)) .dsvert_dp_glmm_random_slope_candidates(
    artifact$candidate_grid, 1L + length(predictor_order), effects) else list()
  loss_bounds <- if (length(candidates)) vapply(candidates, `[[`, numeric(1L), "loss_bound") else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "GLMM random-slope candidate loss bounds"),
    error = function(error) numeric())
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "GLMM random-slope maxima"), error = function(error) numeric())
  raw <- ceiling(loss_bounds * scale)
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) 2 else if (
    identical(adjacency, "add_remove_patient")) 1 else NA_real_
  valid <- .dsvert_dp_has_exact_names(cluster, c("column", "levels")) &&
    .dsvert_dp_is_string(cluster$column) && length(levels) >= 2L &&
    !anyDuplicated(levels) && all(nzchar(trimws(levels))) &&
    isTRUE(predictors_valid) && .dsvert_dp_is_integer(bits, 8L, 18L) &&
    identical(2^bits, scale) && .dsvert_dp_is_integer(artifact$observation_capacity,
      capacity, capacity) && .dsvert_dp_is_integer(cluster_capacity, 2L, capacity) &&
    outcome$lower == 0 && outcome$upper == 1 && isTRUE(artifact$intercept) &&
    isTRUE(effects_valid) && length(candidates) && length(candidates) <= 128L &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$quadrature_rule,
              .dsvert_dp_glmm_random_slope_quadrature_rule_v1(length(effects))) &&
    identical(artifact$candidate_order, "canonical_signed_candidate_grid_v1") &&
    length(observed_loss_bounds) == length(candidates) &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, length(candidates), length(candidates)) &&
    identical(artifact$coordinate_order,
      "signed_candidate_grid_cluster_marginal_binary_negative_log_likelihood_v1") &&
    identical(artifact$source_coordinate_scaling,
      "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "require_one_binary_outcome_and_mean_once_per_admitted_patient",
      "with_one_consistent_public_cluster_level_v1", sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "nonbinary_or_missing_outcome_or_missing_or_nonfinite_predictor",
      "or_missing_or_inconsistent_cluster_excludes_patient_v1", sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient_binary_log_likelihood_contribution_in_one",
      "consistent_cluster_for_every_signed_candidate_v1", sep = "_")) &&
    identical(maximum, capacity * raw) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l1_sensitivity), multiplier * sum(raw), tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l2_sensitivity), multiplier * sqrt(sum(raw^2)), tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l1_sensitivity), multiplier * sum(raw) / scale, tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity), multiplier * sqrt(sum(raw^2)) / scale, tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis,
      "one_patient_changes_one_cluster_marginal_log_likelihood_by_at_most_its_signed_candidate_loss_bound_v1") &&
    identical(artifact$estimation_scope,
      "bounded_binary_random_slope_marginal_likelihood_finite_signed_parameter_grid_v1")
  if (!isTRUE(valid)) stop("The signed binary random-slope GLMM descriptor is invalid",
                           call. = FALSE)
  artifact$outcome <- outcome
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

.dsvert_dp_glmm_random_slope_grid_moment <- function(coordinates, artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released binary random-slope GLMM grid violates its signed bounds",
         call. = FALSE)
  }
  candidate <- .dsvert_dp_glmm_random_slope_candidates(
    artifact$candidate_grid, 1L + length(artifact$predictor_order),
    artifact$random_effect_order)[[which.min(coordinates)[[1L]]]]
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
  selected <- which.min(coordinates)[[1L]]
  list(status = "ok", coefficients = c(`(Intercept)` = intercept, slopes),
       normalized_coefficients = stats::setNames(candidate$beta, artifact$design_terms),
       random_effect_covariance = covariance,
       random_effect_order = artifact$random_effect_order,
       selected_candidate = as.integer(selected),
       selected_dp_negative_log_likelihood = coordinates[[selected]] /
         (2^as.numeric(artifact$numeric_grid_bits)),
       candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1")
}

.dsvert_dp_glmm_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale,
    capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  if (is.list(artifact) && identical(
        artifact$version, .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    return(.dsvert_dp_glmm_random_slope_grid_artifact(
      manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity))
  }
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "cluster", "predictors", "predictor_order", "intercept",
    "design_terms", "observation_capacity", "max_patients_per_cluster",
    "beta_grid", "variance_grid", "quadrature_rule", "candidate_order",
    "candidate_loss_bounds", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "statistic_maximum", "source_raw_l1_sensitivity",
    "source_raw_l2_sensitivity", "natural_l1_sensitivity",
    "natural_l2_sensitivity", "adjacency", "adjacency_sensitivity_basis",
    "estimation_scope", "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version, .DSVERT_CLIENT_DP_GLMM_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "binary_random_intercept_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid binary random-intercept likelihood",
      " grid artifact '", analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_gaussian_bound(artifact$outcome, "GLMM outcome")
  cluster <- artifact$cluster
  levels <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
    cluster$levels, "GLMM cluster levels"), error = function(error) character())
  cluster_valid <- .dsvert_dp_has_exact_names(cluster, c("column", "levels")) &&
    .dsvert_dp_is_string(cluster$column) && length(levels) >= 2L &&
    !anyDuplicated(levels) && all(nzchar(trimws(levels))) &&
    !identical(cluster$column, outcome$column)
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "GLMM predictor order", sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "GLMM design terms", sorted = FALSE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order) &&
    !outcome$column %in% predictor_order &&
    !cluster$column %in% predictor_order
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) {
        .dsvert_dp_gaussian_bound(predictors[[variable]], "GLMM predictor")
      }), predictor_order), error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(
      predictors, function(value) !identical(value$column, outcome$column) &&
        !identical(value$column, cluster$column), logical(1L)))
  }
  beta_grid <- artifact$beta_grid
  if (!is.list(beta_grid) || !length(beta_grid) || !is.null(names(beta_grid))) {
    beta_grid <- list()
  } else {
    beta_grid <- lapply(beta_grid, function(beta) tryCatch(
      .dsvert_dp_capsule_manifest_number_array(beta, "GLMM beta grid row"),
      error = function(error) numeric()))
  }
  beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
    .dsvert_joint_dp_client_json(as.list(beta))
  }, character(1L)) else character()
  variance_grid <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
    artifact$variance_grid, "GLMM variance grid"),
    error = function(error) numeric())
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(
    as.numeric(artifact$max_patients_per_cluster))
  dimension <- 1L + length(predictor_order)
  beta_valid <- length(beta_grid) && all(vapply(beta_grid, function(beta) {
    length(beta) == dimension && !anyNA(beta) && all(is.finite(beta)) &&
      all(abs(beta) <= 8)
  }, logical(1L))) && !anyDuplicated(beta_keys) &&
    identical(beta_keys, sort(beta_keys, method = "radix"))
  variance_valid <- length(variance_grid) && !anyNA(variance_grid) &&
    all(is.finite(variance_grid)) && all(variance_grid >= 0) &&
    all(variance_grid <= 16) && all(diff(variance_grid) > 0)
  candidate_count <- length(beta_grid) * length(variance_grid)
  candidates <- if (isTRUE(beta_valid) && isTRUE(variance_valid)) unlist(
    lapply(variance_grid, function(variance) lapply(beta_grid, function(beta) {
      eta_bound <- sum(abs(beta)) + sqrt(2 * variance) *
        .DSVERT_CLIENT_DP_GLMM_GRID_MAX_GH_NODE
      list(beta = beta, variance = variance,
           loss_bound = log1p(exp(eta_bound)))
    })), recursive = FALSE) else list()
  loss_bounds <- if (length(candidates)) vapply(
    candidates, `[[`, numeric(1L), "loss_bound") else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "GLMM candidate loss bounds"),
    error = function(error) numeric())
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "GLMM statistic maxima"),
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
    isTRUE(all.equal(outcome$lower, 0, tolerance = 0)) &&
    isTRUE(all.equal(outcome$upper, 1, tolerance = 0)) &&
    .dsvert_dp_is_integer(bits, 8L, 18L) && identical(2^bits, scale) &&
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
      "require_one_binary_outcome_and_mean_once_per_admitted_patient",
      "with_one_consistent_public_cluster_level_v1", sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "nonbinary_or_missing_outcome_or_missing_or_nonfinite_predictor",
      "or_missing_or_inconsistent_cluster_excludes_patient_v1", sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient_binary_log_likelihood_contribution_in_one",
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
      "most_its_signed_candidate_loss_bound_v1", sep = "_")) &&
    identical(artifact$estimation_scope, paste(
      "bounded_binary_random_intercept_marginal_likelihood_fixed",
      "covariates_finite_signed_parameter_grid_v1", sep = "_"))
  if (!isTRUE(valid)) {
    stop("The signed binary random-intercept GLMM descriptor is invalid",
         call. = FALSE)
  }
  artifact$outcome <- outcome
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

.dsvert_dp_glmm_grid_moment <- function(coordinates, artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released binary random-intercept GLMM grid violates its signed bounds",
         call. = FALSE)
  }
  candidate <- which.min(coordinates)[[1L]]
  beta_count <- length(artifact$beta_grid)
  variance_index <- as.integer((candidate - 1L) %/% beta_count) + 1L
  beta_index <- as.integer((candidate - 1L) %% beta_count) + 1L
  beta <- artifact$beta_grid[[beta_index]]
  predictor_ranges <- vapply(artifact$predictors, function(bound) {
    bound$upper - bound$lower
  }, numeric(1L))
  predictor_lowers <- vapply(artifact$predictors, `[[`, numeric(1L), "lower")
  slopes <- beta[-1L] / predictor_ranges
  names(slopes) <- artifact$predictor_order
  intercept <- beta[[1L]] - sum(slopes * predictor_lowers)
  list(
    status = "ok",
    coefficients = c(`(Intercept)` = intercept, slopes),
    normalized_coefficients = stats::setNames(beta, artifact$design_terms),
    random_intercept_variance = artifact$variance_grid[[variance_index]],
    selected_candidate = as.integer(candidate),
    selected_beta_index = as.integer(beta_index),
    selected_variance_index = as.integer(variance_index),
    selected_dp_negative_log_likelihood = coordinates[[candidate]] /
      (2^as.numeric(artifact$numeric_grid_bits)),
    candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1")
}

.dsvert_dp_glmm_grid_synopsis_release <- function(
    data_name, analysis_id, server = NULL, datasources = NULL, .aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  if (!is.null(server)) server <- .dsvert_dp_server(server, datasources)
  run <- .dsvert_dp_synopsis_vector_run(datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  metadata <- .dsvert_dp_vector_public_metadata(context)
  scale <- as.numeric(context$lattice$output_lattice_scale)
  count_block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count",
    description = "signed admitted-count capacity block")
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  artifact <- .dsvert_dp_glmm_grid_artifact(
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
    stop("The signed GLMM grid artifact does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  moment <- if (identical(
        artifact$version, .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    .dsvert_dp_glmm_random_slope_grid_moment(coordinates, artifact)
  } else .dsvert_dp_glmm_grid_moment(coordinates, artifact)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !verification$artifact$version %in% c(
        .DSVERT_CLIENT_DP_GLMM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    stop("The binary GLMM Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(context = context, metadata = metadata,
       artifact = verification$artifact, block = blocks[[1L]],
       coordinates = verification$coordinates,
       moment = verification$validated_moment %||% moment,
       certificate = certificate, verification = verification,
       scale = scale, capacity = capacity)
}

.dsvert_dp_glmm_grid_impl <- function(
    formula, data_name, cluster_col, analysis_id, server = NULL,
    datasources = NULL, .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(analysis_id, "analysis_id")
  outcome <- as.character(formula[[2L]])
  predictors <- attr(stats::terms(formula), "term.labels")
  released <- .dsvert_dp_glmm_grid_synopsis_release(
    data_name, analysis_id, server, datasources, .aggregate)
  artifact <- released$artifact
  if (!identical(artifact$outcome$column, outcome) ||
      !identical(artifact$cluster$column, cluster_col) ||
      !setequal(artifact$predictor_order, predictors)) {
    stop("formula and cluster_col must match the signed GLMM grid artifact",
         call. = FALSE)
  }
  moment <- released$moment
  random_slope <- identical(
    artifact$version, .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)
  result <- c(released$metadata, list(
    status = moment$status, analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    family = if (isTRUE(random_slope)) "binomial_random_slope" else
      "binomial_random_intercept",
    estimand = artifact$estimation_scope,
    coefficients = moment$coefficients,
    normalized_coefficients = moment$normalized_coefficients,
    sigma_b2 = if (isTRUE(random_slope)) moment$random_effect_covariance[
      "(Intercept)", "(Intercept)"] else moment$random_intercept_variance,
    random_effect_covariance = if (isTRUE(random_slope))
      moment$random_effect_covariance else NULL,
    random_effect_order = if (isTRUE(random_slope)) moment$random_effect_order else NULL,
    selected_candidate = moment$selected_candidate,
    selected_beta_index = moment$selected_beta_index,
    selected_variance_index = moment$selected_variance_index,
    selected_dp_negative_log_likelihood =
      moment$selected_dp_negative_log_likelihood,
    candidate_selection = moment$candidate_selection,
    standard_errors = NULL, p_values = NULL,
    inference = list(classical_standard_errors = NULL, p_values = NULL,
      confidence_intervals = NULL, sampling_inference_available = FALSE),
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE,
    provenance_certificate = released$certificate,
    disclosure_guard = list(satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("ds.vertGLMM", "list")
  result
}
