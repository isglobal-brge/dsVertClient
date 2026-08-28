# Internal validator and post-processing for a custodian-signed Gaussian AR(1)
# working-GLS candidate grid.  It consumes one authenticated Synopsis vector;
# it never opens the legacy clustered session helpers.

.DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION <-
  "bounded-gaussian-ar1-working-gls-grid-v1"

.DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION <-
  "bounded-gaussian-ar1-robust-working-gls-grid-v1"

.dsvert_dp_gee_ar1_grid_candidates <- function(candidate_grid, cluster_capacity) {
  if (!is.list(candidate_grid) || !length(candidate_grid) ||
      !is.finite(cluster_capacity) || cluster_capacity < 2) return(list())
  candidates <- lapply(candidate_grid, function(candidate) {
    if (!is.list(candidate) || !setequal(names(candidate), c("beta", "rho"))) {
      return(NULL)
    }
    beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$beta, "GEE AR1 beta grid row"), error = function(error) numeric())
    rho <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
      candidate$rho, "GEE AR1 rho grid value"), error = function(error) numeric())
    if (!length(beta) || any(!is.finite(beta)) || any(abs(beta) > 8) ||
        sum(abs(beta)) > 16 || length(rho) != 1L || !is.finite(rho) ||
        abs(rho) > 0.8) return(NULL)
    residual_bound <- 1 + sum(abs(beta))
    loss_bound <- 0.5 * cluster_capacity * residual_bound^2 *
      (1 + abs(rho)) / (1 - abs(rho))
    if (!is.finite(loss_bound) || loss_bound <= 0) return(NULL)
    list(beta = as.numeric(beta), rho = as.numeric(rho),
         loss_bound = as.numeric(loss_bound))
  })
  if (any(vapply(candidates, is.null, logical(1L)))) return(list())
  keys <- vapply(candidates, function(candidate) .dsvert_joint_dp_client_json(list(
    beta = as.list(candidate$beta), rho = candidate$rho)), character(1L))
  if (anyDuplicated(keys) || !identical(keys, sort(keys, method = "radix"))) {
    return(list())
  }
  candidates
}

.dsvert_dp_gee_ar1_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity) {
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  robust_sandwich <- is.list(artifact) && identical(
    artifact$version, .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "cluster", "order", "predictors", "predictor_order",
    "intercept", "design_terms", "observation_capacity",
    "max_patients_per_cluster", "candidate_grid", "candidate_order",
    "candidate_loss_bounds", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "statistic_maximum", "source_raw_l1_sensitivity",
    "source_raw_l2_sensitivity", "natural_l1_sensitivity",
    "natural_l2_sensitivity", "adjacency", "adjacency_sensitivity_basis",
    "estimation_scope", "implementation_state", "cross_owner_state")
  if (isTRUE(robust_sandwich)) {
    required <- c(required, "score_clip", "candidate_bread_bounds",
                  "candidate_meat_bounds", "candidate_loss_sensitivity_bounds",
                  "candidate_bread_sensitivity_bounds",
                  "candidate_meat_sensitivity_bounds", "public_cluster_levels")
  }
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    artifact$version %in% c(.DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION,
                            .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, if (isTRUE(robust_sandwich)) {
      "gaussian_ar1_robust_working_gls_grid_v1"
    } else {
      "gaussian_ar1_working_gls_grid_v1"
    }) &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) && .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed capsule has no valid Gaussian AR1 working-GLS artifact '",
      analysis_id, "' for dataset '", data_name, "'"))
  }
  outcome <- tryCatch(.dsvert_dp_gaussian_bound(artifact$outcome, "GEE outcome"),
                      error = function(error) NULL)
  order <- tryCatch(.dsvert_dp_gaussian_bound(artifact$order, "GEE order"),
                    error = function(error) NULL)
  cluster <- artifact$cluster
  cluster_levels <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    cluster$levels, "GEE cluster levels"), error = function(error) character())
  cluster_valid <- is.list(cluster) && .dsvert_dp_has_exact_names(
    cluster, c("column", "levels")) && .dsvert_dp_is_string(cluster$column) &&
    length(cluster_levels) >= 2L && !anyDuplicated(cluster_levels) &&
    all(nzchar(trimws(cluster_levels)))
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "GEE predictor order", sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, "GEE design terms"), error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order)
  if (isTRUE(predictors_valid) && !is.null(outcome) && !is.null(order) &&
      isTRUE(cluster_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) .dsvert_dp_gaussian_bound(
        predictors[[variable]], "GEE predictor")), predictor_order),
      error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(
      predictors, function(value) !identical(value$column, outcome$column) &&
        !identical(value$column, order$column) &&
        !identical(value$column, cluster$column), logical(1L)))
  }
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  cluster_capacity <- suppressWarnings(as.numeric(artifact$max_patients_per_cluster))
  candidates <- .dsvert_dp_gee_ar1_grid_candidates(
    artifact$candidate_grid, cluster_capacity)
  candidate_valid <- length(candidates) && length(candidates) <=
    (if (isTRUE(robust_sandwich)) 64L else 128L) &&
    all(vapply(candidates, function(candidate) {
      length(candidate$beta) == 1L + length(predictor_order)
    }, logical(1L)))
  loss_bounds <- if (isTRUE(candidate_valid)) {
    vapply(candidates, `[[`, numeric(1L), "loss_bound")
  } else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "GEE AR1 candidate loss bounds"),
    error = function(error) numeric())
  maximums <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "GEE AR1 statistic maxima"),
    error = function(error) numeric())
  score_clip <- suppressWarnings(as.numeric(artifact$score_clip))
  bread_bounds <- if (isTRUE(robust_sandwich) && isTRUE(candidate_valid)) {
    vapply(candidates, function(candidate) {
      cluster_capacity * (1 + abs(candidate$rho)) / (1 - abs(candidate$rho))
    }, numeric(1L))
  } else numeric()
  meat_bounds <- if (isTRUE(robust_sandwich) && length(score_clip) == 1L &&
      is.finite(score_clip)) rep(score_clip^2, length(candidates)) else numeric()
  local_ar1_change <- if (isTRUE(robust_sandwich) && isTRUE(candidate_valid)) {
    vapply(candidates, function(candidate) {
      2 * (1 + candidate$rho^2 + 6 * abs(candidate$rho)) /
        (1 - candidate$rho^2)
    }, numeric(1L))
  } else numeric()
  loss_sensitivity_bounds <- if (length(local_ar1_change)) {
    local_ar1_change * vapply(candidates, function(candidate) {
      (1 + sum(abs(candidate$beta)))^2
    }, numeric(1L))
  } else numeric()
  bread_sensitivity_bounds <- local_ar1_change
  meat_sensitivity_bounds <- if (length(meat_bounds)) 4 * meat_bounds else numeric()
  observed_bread_bounds <- if (isTRUE(robust_sandwich)) tryCatch(
    .dsvert_dp_capsule_manifest_numbers(
      artifact$candidate_bread_bounds, "GEE AR1 candidate bread bounds"),
    error = function(error) numeric()) else numeric()
  observed_meat_bounds <- if (isTRUE(robust_sandwich)) tryCatch(
    .dsvert_dp_capsule_manifest_numbers(
      artifact$candidate_meat_bounds, "GEE AR1 candidate meat bounds"),
    error = function(error) numeric()) else numeric()
  observed_loss_sensitivity_bounds <- if (isTRUE(robust_sandwich)) tryCatch(
    .dsvert_dp_capsule_manifest_numbers(
      artifact$candidate_loss_sensitivity_bounds,
      "GEE AR1 candidate loss sensitivity bounds"),
    error = function(error) numeric()) else numeric()
  observed_bread_sensitivity_bounds <- if (isTRUE(robust_sandwich)) tryCatch(
    .dsvert_dp_capsule_manifest_numbers(
      artifact$candidate_bread_sensitivity_bounds,
      "GEE AR1 candidate bread sensitivity bounds"),
    error = function(error) numeric()) else numeric()
  observed_meat_sensitivity_bounds <- if (isTRUE(robust_sandwich)) tryCatch(
    .dsvert_dp_capsule_manifest_numbers(
      artifact$candidate_meat_sensitivity_bounds,
      "GEE AR1 candidate meat sensitivity bounds"),
    error = function(error) numeric()) else numeric()
  if (isTRUE(robust_sandwich)) {
    triangle_count <- (1L + length(predictor_order)) *
      (2L + length(predictor_order)) / 2L
    loss_raw <- ceiling(loss_bounds * scale)
    bread_raw <- ceiling(2 * bread_bounds * scale)
    meat_raw <- ceiling(2 * meat_bounds * scale)
    public_cluster_levels <- suppressWarnings(as.numeric(
      artifact$public_cluster_levels))
    expected_maximum <- unlist(lapply(seq_along(candidates), function(index) {
      c(capacity * loss_raw[[index]],
        rep(public_cluster_levels * bread_raw[[index]], triangle_count),
        rep(public_cluster_levels * meat_raw[[index]], triangle_count))
    }), use.names = FALSE)
    loss_sensitivity_raw <- ceiling(loss_sensitivity_bounds * scale) + 2
    bread_sensitivity_raw <- ceiling(bread_sensitivity_bounds * scale) + 2
    meat_sensitivity_raw <- ceiling(meat_sensitivity_bounds * scale) + 2
    raw_coordinate_bounds <- unlist(lapply(seq_along(candidates), function(index) {
      c(loss_sensitivity_raw[[index]],
        rep(bread_sensitivity_raw[[index]], triangle_count),
        rep(meat_sensitivity_raw[[index]], triangle_count))
    }), use.names = FALSE)
    expected_coordinate_count <- length(raw_coordinate_bounds)
    raw_l1_base <- sum(raw_coordinate_bounds)
    raw_l2_base <- sqrt(sum(raw_coordinate_bounds^2))
  } else {
    raw_per_candidate <- ceiling(loss_bounds * scale)
    expected_maximum <- capacity * raw_per_candidate
    expected_coordinate_count <- length(candidates)
    raw_l1_base <- sum(raw_per_candidate)
    raw_l2_base <- sqrt(sum(raw_per_candidate^2))
    public_cluster_levels <- numeric()
  }
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) {
    2
  } else if (identical(adjacency, "add_remove_patient")) {
    1
  } else NA_real_
  raw_l1 <- multiplier * raw_l1_base
  raw_l2 <- multiplier * raw_l2_base
  robust_contract <- !isTRUE(robust_sandwich) ||
    (length(predictor_order) <= 3L &&
     .dsvert_dp_is_integer(cluster_capacity, 2L, min(capacity, 32L)) &&
     length(score_clip) == 1L && is.finite(score_clip) &&
     score_clip >= 0.25 && score_clip <= 32 &&
     .dsvert_dp_is_integer(public_cluster_levels, length(cluster_levels),
                           length(cluster_levels)) &&
     public_cluster_levels <= cluster_capacity &&
     length(observed_bread_bounds) == length(bread_bounds) &&
     length(observed_meat_bounds) == length(meat_bounds) &&
     length(observed_loss_sensitivity_bounds) == length(loss_sensitivity_bounds) &&
     length(observed_bread_sensitivity_bounds) == length(bread_sensitivity_bounds) &&
     length(observed_meat_sensitivity_bounds) == length(meat_sensitivity_bounds) &&
     isTRUE(all.equal(observed_bread_bounds, bread_bounds, tolerance = 1e-12)) &&
     isTRUE(all.equal(observed_meat_bounds, meat_bounds, tolerance = 1e-12)) &&
     isTRUE(all.equal(observed_loss_sensitivity_bounds,
                      loss_sensitivity_bounds, tolerance = 1e-12)) &&
     isTRUE(all.equal(observed_bread_sensitivity_bounds,
                      bread_sensitivity_bounds, tolerance = 1e-12)) &&
     isTRUE(all.equal(observed_meat_sensitivity_bounds,
                      meat_sensitivity_bounds, tolerance = 1e-12)))
  valid <- !is.null(outcome) && !is.null(order) && cluster_valid &&
    !identical(outcome$column, order$column) &&
    !outcome$column %in% c(cluster$column, predictor_order) &&
    !order$column %in% c(cluster$column, predictor_order) &&
    isTRUE(predictors_valid) && .dsvert_dp_is_integer(bits, 8L, 18L) &&
    identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    .dsvert_dp_is_integer(cluster_capacity, 2L, capacity) &&
    isTRUE(candidate_valid) && isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$candidate_order, "canonical_beta_rho_grid_v1") &&
    length(observed_loss_bounds) == length(loss_bounds) &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, expected_coordinate_count,
                          expected_coordinate_count) &&
    identical(artifact$coordinate_order, if (isTRUE(robust_sandwich)) {
      "signed_candidate_grid_cluster_gaussian_ar1_loss_bread_meat_upper_v1"
    } else {
      "signed_candidate_grid_cluster_gaussian_ar1_working_gls_loss_v1"
    }) &&
    identical(artifact$source_coordinate_scaling,
              "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "require_one_complete_bounded_row_per_admitted_patient_with_one",
      "consistent_public_cluster_level_and_strict_within_cluster_order_v1")) &&
    identical(artifact$missingness_policy, paste(
      "missing_or_nonfinite_outcome_predictor_or_order_or_missing_or",
      "inconsistent_cluster_excludes_patient_and_order_ties_reject_v1")) &&
    identical(artifact$contribution_domain, if (isTRUE(robust_sandwich)) {
      paste("one_bounded_patient_can_change_local_ar1_loss_and_bread",
            "terms_and_at_most_two_componentwise_clipped_cluster_score",
            "meat_products_per_signed_candidate_v1", sep = "_")
    } else {
      paste("one_bounded_patient_can_change_one_clipped_cluster_gaussian_ar1",
            "working_gls_loss_per_signed_candidate_v1")
    }) &&
    identical(maximums, expected_maximum) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l1_sensitivity), raw_l1,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l2_sensitivity), raw_l2,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l1_sensitivity), raw_l1 / scale,
                     tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity), raw_l2 / scale,
                     tolerance = 1e-12)) && identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis, if (isTRUE(robust_sandwich)) {
      paste("one_patient_removal_insertion_or_order_change_affects_at",
            "most_two_local_ar1_neighborhoods_and_two_clipped_score",
            "outer_products_with_quantization_slack_v1", sep = "_")
    } else {
      paste("one_patient_can_change_one_entire_clipped_cluster_ar1_working",
            "gls_loss_by_at_most_its_signed_bound_v1")
    }) &&
    identical(artifact$estimation_scope, if (isTRUE(robust_sandwich)) {
      paste("bounded_gaussian_ar1_working_gls_finite_signed_beta_rho_grid",
            "with_componentwise_clipped_cluster_score_sandwich_v1", sep = "_")
    } else {
      "bounded_gaussian_ar1_working_gls_finite_signed_beta_rho_grid_v1"
    }) && isTRUE(robust_contract)
  if (!isTRUE(valid)) {
    stop("The signed Gaussian AR1 working-GLS descriptor is invalid", call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$order <- order
  artifact$cluster <- list(column = cluster$column, levels = cluster_levels)
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$candidate_grid <- lapply(candidates, function(candidate) {
    list(beta = candidate$beta, rho = candidate$rho)
  })
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$coordinate_count <- as.integer(expected_coordinate_count)
  artifact$statistic_maximum <- expected_maximum
  if (isTRUE(robust_sandwich)) {
    artifact$score_clip <- score_clip
    artifact$candidate_bread_bounds <- bread_bounds
    artifact$candidate_meat_bounds <- meat_bounds
    artifact$candidate_loss_sensitivity_bounds <- loss_sensitivity_bounds
    artifact$candidate_bread_sensitivity_bounds <- bread_sensitivity_bounds
    artifact$candidate_meat_sensitivity_bounds <- meat_sensitivity_bounds
    artifact$public_cluster_levels <- as.integer(public_cluster_levels)
  }
  artifact
}

.dsvert_dp_gee_ar1_grid_moment <- function(coordinates, artifact) {
  scale <- 2^as.numeric(artifact$numeric_grid_bits)
  upper <- as.numeric(artifact$statistic_maximum) / scale
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper)) {
    stop("The released Gaussian AR1 working-GLS grid violates its signed bounds",
         call. = FALSE)
  }
  robust_sandwich <- identical(
    artifact$version, .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION)
  dimension <- 1L + length(artifact$predictor_order)
  triangle_count <- dimension * (dimension + 1L) / 2L
  candidate_width <- if (isTRUE(robust_sandwich)) 1L + 2L * triangle_count else 1L
  loss_coordinates <- coordinates[seq.int(1L, length(coordinates),
                                          by = candidate_width)]
  candidate_index <- which.min(loss_coordinates)[[1L]]
  candidate <- artifact$candidate_grid[[candidate_index]]
  ranges <- vapply(artifact$predictors, function(bound) bound$upper - bound$lower,
                   numeric(1L))
  lowers <- vapply(artifact$predictors, function(bound) bound$lower,
                   numeric(1L))
  outcome_range <- artifact$outcome$upper - artifact$outcome$lower
  slopes <- outcome_range * candidate$beta[-1L] / ranges
  names(slopes) <- artifact$predictor_order
  intercept <- artifact$outcome$lower + outcome_range * candidate$beta[[1L]] -
    sum(slopes * lowers)
  coefficients <- stats::setNames(c(intercept, unname(slopes)),
                                  c("(Intercept)", names(slopes)))
  if (isTRUE(robust_sandwich)) {
    start <- (candidate_index - 1L) * candidate_width
    upper_index <- which(upper.tri(matrix(0, dimension, dimension), diag = TRUE))
    inflate <- function(values) {
      matrix_value <- matrix(0, dimension, dimension)
      matrix_value[upper_index] <- values
      matrix_value[lower.tri(matrix_value)] <- t(matrix_value)[lower.tri(matrix_value)]
      matrix_value
    }
    bread_shifted <- coordinates[start + seq_len(triangle_count) + 1L]
    meat_shifted <- coordinates[start + triangle_count + seq_len(triangle_count) +
      1L]
    bread <- inflate(bread_shifted - artifact$public_cluster_levels *
      artifact$candidate_bread_bounds[[candidate_index]])
    meat <- inflate(meat_shifted - artifact$public_cluster_levels *
      artifact$candidate_meat_bounds[[candidate_index]])
    bread <- (bread + t(bread)) / 2
    meat <- (meat + t(meat)) / 2
    bread_values <- tryCatch(eigen(bread, symmetric = TRUE,
                                   only.values = TRUE)$values,
                              error = function(error) numeric())
    bread_tolerance <- if (length(bread_values)) {
      max(1e-10, max(abs(bread_values)) * 1e-8)
    } else {
      Inf
    }
    source_transform <- diag(c(outcome_range, outcome_range / ranges))
    source_transform[1L, -1L] <- -outcome_range * lowers / ranges
    dimnames(source_transform) <- list(names(coefficients), artifact$design_terms)
    robust_covariance <- NULL
    meat_projected <- FALSE
    bread_status <- "non_identifiable_dp_bread"
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
        normalized_covariance <- inverse_bread %*% meat_psd %*% inverse_bread
        robust_covariance <- source_transform %*% normalized_covariance %*%
          t(source_transform)
        robust_covariance <- (robust_covariance + t(robust_covariance)) / 2
        dimnames(robust_covariance) <- list(names(coefficients), names(coefficients))
        bread_status <- "ok"
      }
    }
    return(list(
      status = bread_status,
      coefficients = coefficients,
      normalized_coefficients = stats::setNames(candidate$beta,
                                                 artifact$design_terms),
      working_correlation = candidate$rho,
      selected_candidate = as.integer(candidate_index),
      selected_dp_working_gls_loss = loss_coordinates[[candidate_index]],
      candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1",
      robust_covariance = robust_covariance,
      robust_covariance_status = bread_status,
      meat_psd_projection_applied = meat_projected,
      inference = "no_standard_errors_or_p_values_for_clipped_score_sandwich_v1"))
  }
  list(
    status = "ok",
    coefficients = coefficients,
    normalized_coefficients = stats::setNames(candidate$beta,
                                               artifact$design_terms),
    working_correlation = candidate$rho,
    selected_candidate = as.integer(candidate_index),
    selected_dp_working_gls_loss = loss_coordinates[[candidate_index]],
    candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1")
}

.dsvert_dp_gee_ar1_grid_synopsis_release <- function(
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
  artifact <- .dsvert_dp_gee_ar1_grid_artifact(
    context$manifest, data_name, analysis_id, server, context$adjacency,
    scale, capacity)
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
    stop("The signed Gaussian AR1 working-GLS grid does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  moment <- .dsvert_dp_gee_ar1_grid_moment(coordinates, artifact)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !verification$artifact$version %in% c(
        .DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION)) {
    stop("The Gaussian AR1 working-GLS Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(context = context, metadata = metadata, artifact = verification$artifact,
       block = blocks[[1L]], coordinates = verification$coordinates,
       moment = verification$validated_moment %||% moment,
       certificate = certificate, verification = verification)
}

.dsvert_dp_gee_ar1_grid_impl <- function(
    formula, data_name, analysis_id, id_col, order_col, server = NULL,
    datasources = NULL, .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(analysis_id, "analysis_id")
  terms <- stats::terms(formula)
  outcome <- as.character(formula[[2L]])
  predictors <- attr(terms, "term.labels")
  released <- .dsvert_dp_gee_ar1_grid_synopsis_release(
    data_name, analysis_id, server, datasources, .aggregate)
  artifact <- released$artifact
  if (!identical(artifact$outcome$column, outcome) ||
      !setequal(artifact$predictor_order, predictors) ||
      !identical(artifact$cluster$column, id_col) ||
      !identical(artifact$order$column, order_col)) {
    stop("formula, id_col and order_col must match the signed Gaussian AR1 working-GLS artifact",
         call. = FALSE)
  }
  moment <- released$moment
  result <- c(released$metadata, list(
    status = if (identical(
        artifact$version, .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION)) {
      "public_certified_gaussian_ar1_clipped_score_sandwich_finite_grid"
    } else {
      "public_certified_gaussian_ar1_working_gls_finite_grid"
    },
    family = "gaussian", corstr = "ar1", analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    estimand = artifact$estimation_scope, coefficients = moment$coefficients,
    normalized_coefficients = moment$normalized_coefficients,
    working_correlation = moment$working_correlation,
    selected_candidate = moment$selected_candidate,
    selected_dp_working_gls_loss = moment$selected_dp_working_gls_loss,
    candidate_selection = moment$candidate_selection,
    robust_covariance = moment$robust_covariance %||% NULL,
    covariance = moment$robust_covariance %||% NULL,
    vcov = moment$robust_covariance %||% NULL,
    std_errors = NULL, standard_errors = NULL, p_values = NULL,
    cluster_correlation_estimated = TRUE, cluster_columns = id_col,
    order_column = order_col,
    robust_covariance_status = moment$robust_covariance_status %||% NULL,
    meat_psd_projection_applied = moment$meat_psd_projection_applied %||% NULL,
    inference = moment$inference %||%
      "unavailable_without_protected_cluster_score_and_meat",
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    production_ready = FALSE, additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE, provenance_certificate = released$certificate,
    disclosure_guard = list(satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("dsvert_dp_gaussian_ar1_gee", "ds.vertGEE", "list")
  result
}

.dsvert_dp_gee_ar1_grid_adapter <- function(
    explicit_arguments, formula, data, id_col, order_col, corstr, verbose,
    datasources, analysis_id) {
  if (!identical(corstr, "ar1")) {
    stop("Gaussian AR1 working-GLS requires corstr='ar1'", call. = FALSE)
  }
  if (!is.character(id_col) || length(id_col) != 1L || is.na(id_col) ||
      !nzchar(id_col) || !is.character(order_col) || length(order_col) != 1L ||
      is.na(order_col) || !nzchar(order_col) || identical(id_col, order_col)) {
    stop("Gaussian AR1 working-GLS requires distinct id_col and order_col",
         call. = FALSE)
  }
  allowed <- c("formula", "data", "family", "id_col", "order_col", "corstr",
               "verbose", "datasources", "analysis_id")
  unexpected <- setdiff(explicit_arguments, allowed)
  if (length(unexpected)) {
    stop("Gaussian AR1 working-GLS does not accept legacy controls: ",
         paste(sort(unexpected, method = "radix"), collapse = ", "), call. = FALSE)
  }
  terms <- if (inherits(formula, "formula")) stats::terms(formula) else NULL
  predictors <- if (is.null(terms)) character() else attr(terms, "term.labels")
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]]) || !identical(attr(terms, "intercept"), 1L) ||
      !length(predictors) || any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", predictors))) {
    stop("Gaussian AR1 working-GLS requires an intercept and additive bare column names",
         call. = FALSE)
  }
  resolved <- .dsvert_federation_argument(data, datasources)
  result <- .dsvert_dp_gee_ar1_grid_impl(
    formula = formula, data_name = resolved$value, analysis_id = analysis_id,
    id_col = id_col, order_col = order_col,
    datasources = resolved$datasources, .aggregate = DSI::datashield.aggregate)
  result$called_via <- "ds.vertGEE_gaussian_ar1_analysis_id"
  result
}
