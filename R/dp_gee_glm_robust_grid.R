# Internal validator and post-processing for a same-owner signed binomial or
# Poisson independence GEE grid with a componentwise-clipped cluster-score
# sandwich. It consumes one authenticated Synopsis vector and exposes neither
# rows, cluster scores nor an optimizer.

.DSVERT_CLIENT_DP_GEE_GLM_ROBUST_GRID_ARTIFACT_VERSIONS <- c(
  binomial = "bounded-binomial-robust-independence-gee-grid-v1",
  poisson = "bounded-poisson-robust-independence-gee-grid-v1")

.dsvert_dp_gee_glm_robust_grid_candidates <- function(
    family, beta_grid, max_outcome = NULL) {
  if (!family %in% names(.DSVERT_CLIENT_DP_GEE_GLM_ROBUST_GRID_ARTIFACT_VERSIONS) ||
      !is.list(beta_grid) || !length(beta_grid) || !is.null(names(beta_grid))) {
    return(list())
  }
  poisson <- identical(family, "poisson")
  if (isTRUE(poisson) && !.dsvert_dp_is_integer(max_outcome, 1L, 1024L)) {
    return(list())
  }
  softplus <- function(value) pmax(value, 0) + log1p(exp(-abs(value)))
  candidates <- lapply(beta_grid, function(value) {
    beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      value, paste(family, "robust independence GEE beta grid row")),
      error = function(error) numeric())
    if (!length(beta) || anyNA(beta) || any(!is.finite(beta)) ||
        any(abs(beta) > 8) || sum(abs(beta)) > 8) return(NULL)
    eta_bound <- sum(abs(beta))
    loss <- if (isTRUE(poisson)) {
      outer(0:as.integer(max_outcome), c(-eta_bound, eta_bound),
            function(count, eta) exp(eta) - count * eta + lgamma(count + 1))
    } else {
      outer(c(0, 1), c(-eta_bound, eta_bound),
            function(observed, eta) softplus(eta) - observed * eta)
    }
    loss_bound <- max(0, max(loss))
    bread_bound <- if (isTRUE(poisson)) exp(eta_bound) else 0.25
    if (!is.finite(loss_bound) || !is.finite(bread_bound) ||
        loss_bound <= 0 || bread_bound <= 0) return(NULL)
    list(beta = beta, loss_bound = loss_bound, bread_bound = bread_bound)
  })
  if (any(vapply(candidates, is.null, logical(1L)))) return(list())
  keys <- vapply(candidates, function(candidate) .dsvert_joint_dp_client_json(
    as.list(candidate$beta)), character(1L))
  if (anyDuplicated(keys) || !identical(keys, sort(keys, method = "radix"))) {
    return(list())
  }
  candidates
}

.dsvert_dp_gee_glm_robust_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity,
    family) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  if (!family %in% names(.DSVERT_CLIENT_DP_GEE_GLM_ROBUST_GRID_ARTIFACT_VERSIONS)) {
    stop("The robust independence GEE family is invalid", call. = FALSE)
  }
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "family", "outcome", "max_outcome", "cluster", "predictors",
    "predictor_order", "intercept", "design_terms", "observation_capacity",
    "max_patients_per_cluster", "beta_grid", "candidate_order",
    "candidate_loss_bounds", "score_clip", "candidate_bread_bounds",
    "candidate_meat_bounds", "candidate_loss_sensitivity_bounds",
    "candidate_bread_sensitivity_bounds", "candidate_meat_sensitivity_bounds",
    "public_cluster_levels", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling", "repeated_record_policy",
    "missingness_policy", "contribution_domain", "statistic_maximum",
    "source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
    "natural_l1_sensitivity", "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis", "estimation_scope", "implementation_state",
    "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version,
              .DSVERT_CLIENT_DP_GEE_GLM_ROBUST_GRID_ARTIFACT_VERSIONS[[family]]) &&
    identical(artifact$spec_version,
              paste0(family, "_robust_independence_gee_grid_v1")) &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) && identical(artifact$family, family) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed ", family, " robust independence GEE grid '", analysis_id,
      "' is unavailable for dataset '", data_name, "'"))
  }
  outcome <- tryCatch(.dsvert_dp_gaussian_bound(
    artifact$outcome, paste(family, "robust independence GEE outcome")),
    error = function(error) NULL)
  cluster <- artifact$cluster
  cluster_levels <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    cluster$levels, "robust independence GEE cluster levels", sorted = TRUE),
    error = function(error) character())
  cluster_valid <- is.list(cluster) && .dsvert_dp_has_exact_names(
    cluster, c("column", "levels")) && .dsvert_dp_is_string(cluster$column) &&
    length(cluster_levels) >= 2L && length(cluster_levels) <= 64L &&
    !anyDuplicated(cluster_levels) && all(nzchar(trimws(cluster_levels)))
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, paste(family, "robust independence GEE predictor order"),
    sorted = TRUE), error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, paste(family, "robust independence GEE design terms"),
    sorted = FALSE), error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order) &&
    length(predictor_order) <= 3L
  if (isTRUE(predictors_valid) && !is.null(outcome) && isTRUE(cluster_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) .dsvert_dp_gaussian_bound(
        predictors[[variable]], paste(family, "robust independence GEE predictor"))),
      predictor_order), error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(predictors, function(value) {
      !identical(value$column, outcome$column) &&
        !identical(value$column, cluster$column)
    }, logical(1L)))
  }
  maximum <- artifact$max_outcome
  poisson <- identical(family, "poisson")
  maximum_valid <- if (isTRUE(poisson)) {
    .dsvert_dp_is_integer(maximum, 1L, 1024L) && !is.null(outcome) &&
      isTRUE(all.equal(outcome$lower, 0, tolerance = 0)) &&
      isTRUE(all.equal(outcome$upper, as.numeric(maximum), tolerance = 0))
  } else {
    is.null(maximum) && !is.null(outcome) &&
      isTRUE(all.equal(outcome$lower, 0, tolerance = 0)) &&
      isTRUE(all.equal(outcome$upper, 1, tolerance = 0))
  }
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(as.numeric(artifact$max_patients_per_cluster))
  score_clip <- suppressWarnings(as.numeric(artifact$score_clip))
  candidates <- if (isTRUE(maximum_valid)) {
    .dsvert_dp_gee_glm_robust_grid_candidates(
      family, artifact$beta_grid, if (isTRUE(poisson)) as.integer(maximum) else NULL)
  } else list()
  dimension <- 1L + length(predictor_order)
  candidate_valid <- length(candidates) && length(candidates) <= 32L &&
    all(vapply(candidates, function(candidate) length(candidate$beta) == dimension,
               logical(1L)))
  loss_bounds <- if (isTRUE(candidate_valid)) {
    vapply(candidates, `[[`, numeric(1L), "loss_bound")
  } else numeric()
  bread_bounds <- if (isTRUE(candidate_valid)) {
    vapply(candidates, `[[`, numeric(1L), "bread_bound")
  } else numeric()
  meat_bounds <- if (length(score_clip) == 1L && is.finite(score_clip)) {
    rep(score_clip^2, length(candidates))
  } else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, paste(family, "robust GEE loss bounds")),
    error = function(error) numeric())
  observed_bread_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_bread_bounds, paste(family, "robust GEE bread bounds")),
    error = function(error) numeric())
  observed_meat_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_meat_bounds, paste(family, "robust GEE meat bounds")),
    error = function(error) numeric())
  observed_loss_sensitivity_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_sensitivity_bounds,
    paste(family, "robust GEE loss sensitivity bounds")),
    error = function(error) numeric())
  observed_bread_sensitivity_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_bread_sensitivity_bounds,
    paste(family, "robust GEE bread sensitivity bounds")),
    error = function(error) numeric())
  observed_meat_sensitivity_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_meat_sensitivity_bounds,
    paste(family, "robust GEE meat sensitivity bounds")),
    error = function(error) numeric())
  maximums <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, paste(family, "robust GEE statistic maxima")),
    error = function(error) numeric())
  triangle_count <- dimension * (dimension + 1L) / 2L
  loss_raw <- ceiling(loss_bounds * scale)
  bread_raw <- ceiling(bread_bounds * scale)
  meat_raw <- ceiling(2 * meat_bounds * scale)
  loss_sensitivity_raw <- loss_raw
  bread_sensitivity_raw <- bread_raw
  meat_sensitivity_raw <- meat_raw + 2
  public_cluster_levels <- suppressWarnings(as.numeric(artifact$public_cluster_levels))
  expected_maximum <- if (length(loss_raw) && length(public_cluster_levels) == 1L &&
      is.finite(public_cluster_levels)) {
    unlist(lapply(seq_along(candidates), function(index) {
      c(capacity * loss_raw[[index]],
        rep(capacity * bread_raw[[index]], triangle_count),
        rep(public_cluster_levels * meat_raw[[index]], triangle_count))
    }), use.names = FALSE)
  } else numeric()
  raw_coordinate_bounds <- unlist(lapply(seq_along(candidates), function(index) {
    c(loss_sensitivity_raw[[index]],
      rep(bread_sensitivity_raw[[index]], triangle_count),
      rep(meat_sensitivity_raw[[index]], triangle_count))
  }), use.names = FALSE)
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) 2 else if (
    identical(adjacency, "add_remove_patient")) 1 else NA_real_
  raw_l1 <- multiplier * sum(raw_coordinate_bounds)
  raw_l2 <- multiplier * sqrt(sum(raw_coordinate_bounds^2))
  valid <- !is.null(outcome) && isTRUE(cluster_valid) &&
    isTRUE(predictors_valid) && isTRUE(maximum_valid) &&
    .dsvert_dp_is_integer(bits, 8L, 18L) && identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    .dsvert_dp_is_integer(cluster_capacity, 2L, min(capacity, 32L)) &&
    length(score_clip) == 1L && is.finite(score_clip) &&
    score_clip >= 0.25 && score_clip <= 32 &&
    .dsvert_dp_is_integer(public_cluster_levels, length(cluster_levels),
                          length(cluster_levels)) &&
    isTRUE(candidate_valid) && isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$candidate_order, "canonical_beta_grid_glm_v1") &&
    length(observed_loss_bounds) == length(loss_bounds) &&
    length(observed_bread_bounds) == length(bread_bounds) &&
    length(observed_meat_bounds) == length(meat_bounds) &&
    length(observed_loss_sensitivity_bounds) == length(loss_bounds) &&
    length(observed_bread_sensitivity_bounds) == length(bread_bounds) &&
    length(observed_meat_sensitivity_bounds) == length(meat_bounds) &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    isTRUE(all.equal(observed_bread_bounds, bread_bounds, tolerance = 1e-12)) &&
    isTRUE(all.equal(observed_meat_bounds, meat_bounds, tolerance = 1e-12)) &&
    isTRUE(all.equal(observed_loss_sensitivity_bounds, loss_bounds,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(observed_bread_sensitivity_bounds, bread_bounds,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(observed_meat_sensitivity_bounds, 2 * meat_bounds,
                     tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, length(raw_coordinate_bounds),
                          length(raw_coordinate_bounds)) &&
    identical(artifact$coordinate_order, paste(
      "signed_candidate_grid_cluster", family,
      "independence_loss_bread_meat_upper_v1", sep = "_")) &&
    identical(artifact$source_coordinate_scaling,
              "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "require_one_complete_bounded", family,
      "outcome_and_mean_once_per_admitted_patient_with_one",
      "consistent_public_cluster_level_v1", sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "noninteger_or_out_of_range_or_missing_outcome_or_missing_or",
      "nonfinite_predictor_or_missing_or_inconsistent_cluster_excludes",
      "patient_v1", sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient_can_change_one", family,
      "loss_and_bread_term_and_one_componentwise_clipped_cluster_score",
      "outer_product_per_signed_candidate_v1", sep = "_")) &&
    identical(maximums, expected_maximum) &&
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
      "one_patient_removal_insertion_or_replacement_changes_one",
      "bounded_loss_and_bread_term_and_one_clipped_cluster_score_outer",
      "product_with_quantization_slack_v1", sep = "_")) &&
    identical(artifact$estimation_scope, paste(
      "bounded", family,
      "independence_gee_finite_signed_beta_grid_with_componentwise",
      "clipped_cluster_score_sandwich_v1", sep = "_"))
  if (!isTRUE(valid)) {
    stop("The signed robust independence GEE descriptor is invalid", call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$cluster <- list(column = cluster$column, levels = cluster_levels)
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$max_outcome <- if (isTRUE(poisson)) as.integer(maximum) else NULL
  artifact$beta_grid <- lapply(candidates, `[[`, "beta")
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$candidate_bread_bounds <- bread_bounds
  artifact$candidate_meat_bounds <- meat_bounds
  artifact$candidate_loss_sensitivity_bounds <- loss_bounds
  artifact$candidate_bread_sensitivity_bounds <- bread_bounds
  artifact$candidate_meat_sensitivity_bounds <- 2 * meat_bounds
  artifact$public_cluster_levels <- as.integer(public_cluster_levels)
  artifact$coordinate_count <- as.integer(length(raw_coordinate_bounds))
  artifact$statistic_maximum <- expected_maximum
  artifact
}

.dsvert_dp_gee_glm_robust_grid_moment <- function(coordinates, artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released robust independence GEE grid violates its signed bounds",
         call. = FALSE)
  }
  scale <- 2^as.numeric(artifact$numeric_grid_bits)
  dimension <- 1L + length(artifact$predictor_order)
  triangle_count <- dimension * (dimension + 1L) / 2L
  width <- 1L + 2L * triangle_count
  loss_coordinates <- coordinates[seq.int(1L, length(coordinates), by = width)]
  selected <- which.min(loss_coordinates)[[1L]]
  beta <- artifact$beta_grid[[selected]]
  ranges <- vapply(artifact$predictors, function(bound) {
    bound$upper - bound$lower
  }, numeric(1L))
  lowers <- vapply(artifact$predictors, function(bound) bound$lower, numeric(1L))
  slopes <- beta[-1L] / ranges
  names(slopes) <- artifact$predictor_order
  intercept <- beta[[1L]] - sum(slopes * lowers)
  coefficients <- stats::setNames(c(intercept, unname(slopes)),
                                  c("(Intercept)", names(slopes)))
  start <- (selected - 1L) * width
  upper_index <- which(upper.tri(matrix(0, dimension, dimension), diag = TRUE))
  inflate <- function(values) {
    value <- matrix(0, dimension, dimension)
    value[upper_index] <- values
    value[lower.tri(value)] <- t(value)[lower.tri(value)]
    value
  }
  bread <- inflate(coordinates[start + seq_len(triangle_count) + 1L] / scale)
  meat <- inflate(coordinates[start + triangle_count + seq_len(triangle_count) +
    1L] / scale - artifact$public_cluster_levels *
    artifact$candidate_meat_bounds[[selected]])
  bread <- (bread + t(bread)) / 2
  meat <- (meat + t(meat)) / 2
  bread_values <- tryCatch(eigen(bread, symmetric = TRUE,
                                 only.values = TRUE)$values,
                            error = function(error) numeric())
  bread_tolerance <- if (length(bread_values)) {
    max(1e-10, max(abs(bread_values)) * 1e-8)
  } else Inf
  transform <- diag(c(1, 1 / ranges))
  transform[1L, -1L] <- -lowers / ranges
  dimnames(transform) <- list(names(coefficients), artifact$design_terms)
  covariance <- NULL
  meat_projected <- FALSE
  status <- "non_identifiable_dp_bread"
  if (length(bread_values) == dimension && all(is.finite(bread_values)) &&
      min(bread_values) > bread_tolerance) {
    inverse_bread <- tryCatch(solve(bread), error = function(error) NULL)
    meat_decomposition <- tryCatch(eigen(meat, symmetric = TRUE),
                                   error = function(error) NULL)
    if (!is.null(inverse_bread) && !is.null(meat_decomposition) &&
        all(is.finite(meat_decomposition$values))) {
      meat_projected <- any(meat_decomposition$values < -bread_tolerance)
      meat_psd <- meat_decomposition$vectors %*%
        (pmax(meat_decomposition$values, 0) * t(meat_decomposition$vectors))
      covariance <- inverse_bread %*% meat_psd %*% inverse_bread
      covariance <- transform %*% covariance %*% t(transform)
      covariance <- (covariance + t(covariance)) / 2
      dimnames(covariance) <- list(names(coefficients), names(coefficients))
      status <- "ok"
    }
  }
  list(
    status = status, coefficients = coefficients,
    normalized_coefficients = stats::setNames(beta, artifact$design_terms),
    selected_candidate = as.integer(selected),
    selected_dp_negative_log_likelihood = loss_coordinates[[selected]] / scale,
    candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1",
    robust_covariance = covariance, robust_covariance_status = status,
    meat_psd_projection_applied = meat_projected,
    inference = "no_standard_errors_or_p_values_for_clipped_score_sandwich_v1")
}

.dsvert_dp_gee_glm_robust_grid_synopsis_release <- function(
    data_name, analysis_id, family, server = NULL, datasources = NULL,
    .aggregate) {
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
  artifact <- .dsvert_dp_gee_glm_robust_grid_artifact(
    context$manifest, data_name, analysis_id, server, context$adjacency,
    scale, capacity, family)
  blocks <- .dsvert_dp_capsule_vector_blocks(
    context$layout, "gaussian_models", dataset = data_name,
    owner_peer = artifact$owner_peer)
  blocks <- blocks[vapply(blocks, function(block) identical(block$key, analysis_id),
                          logical(1L))]
  signed_descriptor <- tryCatch(context$manifest$workload$families$
    gaussian_models$artifacts[[analysis_id]], error = function(error) NULL)
  if (length(blocks) != 1L || !identical(
        .dsvert_joint_dp_client_json(blocks[[1L]]$descriptor),
        .dsvert_joint_dp_client_json(signed_descriptor))) {
    stop("The signed robust independence GEE grid does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  moment <- .dsvert_dp_gee_glm_robust_grid_moment(coordinates, artifact)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !identical(verification$artifact$version,
                 .DSVERT_CLIENT_DP_GEE_GLM_ROBUST_GRID_ARTIFACT_VERSIONS[[family]])) {
    stop("The robust independence GEE Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(context = context, metadata = metadata, artifact = verification$artifact,
       block = blocks[[1L]], coordinates = verification$coordinates,
       moment = verification$validated_moment %||% moment,
       certificate = certificate, verification = verification)
}

.dsvert_dp_gee_glm_robust_grid_impl <- function(
    formula, data_name, analysis_id, family, id_col, server = NULL,
    datasources = NULL, .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(analysis_id, "analysis_id")
  terms <- stats::terms(formula)
  outcome <- as.character(formula[[2L]])
  predictors <- attr(terms, "term.labels")
  released <- .dsvert_dp_gee_glm_robust_grid_synopsis_release(
    data_name, analysis_id, family, server, datasources, .aggregate)
  artifact <- released$artifact
  if (!identical(artifact$outcome$column, outcome) ||
      !setequal(artifact$predictor_order, predictors) ||
      !identical(artifact$cluster$column, id_col)) {
    stop("formula and id_col must match the signed robust independence GEE artifact",
         call. = FALSE)
  }
  moment <- released$moment
  result <- c(released$metadata, list(
    status = paste0("public_certified_", family,
                    "_independence_clipped_score_sandwich_finite_grid"),
    family = family, corstr = "independence", analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    estimand = artifact$estimation_scope, coefficients = moment$coefficients,
    normalized_coefficients = moment$normalized_coefficients,
    selected_candidate = moment$selected_candidate,
    selected_dp_negative_log_likelihood =
      moment$selected_dp_negative_log_likelihood,
    candidate_selection = moment$candidate_selection,
    robust_covariance = moment$robust_covariance,
    covariance = moment$robust_covariance, vcov = moment$robust_covariance,
    std_errors = NULL, standard_errors = NULL, p_values = NULL,
    cluster_correlation_estimated = FALSE, cluster_columns = id_col,
    robust_covariance_status = moment$robust_covariance_status,
    meat_psd_projection_applied = moment$meat_psd_projection_applied,
    inference = moment$inference, source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE, production_ready = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE, provenance_certificate = released$certificate,
    disclosure_guard = list(satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("dsvert_dp_glm_robust_gee", "ds.vertGEE", "list")
  result
}

.dsvert_dp_gee_glm_robust_grid_adapter <- function(
    explicit_arguments, formula, data, family, id_col, order_col, corstr,
    verbose, datasources, analysis_id) {
  if (!family %in% c("binomial", "poisson") ||
      !identical(corstr, "independence")) {
    stop("robust finite-grid GEE supports binomial or Poisson independence only",
         call. = FALSE)
  }
  if (!is.character(id_col) || length(id_col) != 1L || is.na(id_col) ||
      !nzchar(id_col) || !is.null(order_col)) {
    stop("robust finite-grid GEE requires one id_col and no order_col",
         call. = FALSE)
  }
  allowed <- c("formula", "data", "family", "id_col", "corstr", "verbose",
               "datasources", "analysis_id")
  unexpected <- setdiff(explicit_arguments, allowed)
  if (length(unexpected)) {
    stop("robust finite-grid GEE does not accept legacy controls: ",
         paste(sort(unexpected, method = "radix"), collapse = ", "), call. = FALSE)
  }
  terms <- if (inherits(formula, "formula")) stats::terms(formula) else NULL
  predictors <- if (is.null(terms)) character() else attr(terms, "term.labels")
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]]) || !identical(attr(terms, "intercept"), 1L) ||
      !length(predictors) || length(predictors) > 3L ||
      any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", predictors))) {
    stop(paste(
      "robust finite-grid GEE requires an intercept and one to three additive",
      "bare predictor names"), call. = FALSE)
  }
  resolved <- .dsvert_federation_argument(data, datasources)
  result <- .dsvert_dp_gee_glm_robust_grid_impl(
    formula, resolved$value, analysis_id, family, id_col,
    datasources = resolved$datasources, .aggregate = DSI::datashield.aggregate)
  result$called_via <- "ds.vertGEE_robust_independence_analysis_id"
  result
}
