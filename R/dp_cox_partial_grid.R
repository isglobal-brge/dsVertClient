# Internal post-processing for a same-owner signed Breslow Cox partial-
# likelihood grid.  This opens no legacy Cox session and selects only one
# custodian-signed beta candidate from an authenticated Synopsis vector.

.DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION <-
  "bounded-cox-partial-likelihood-grid-v1"

.dsvert_dp_cox_partial_grid_candidates <- function(value, dimension, capacity) {
  if (!is.list(value) || !length(value) || !is.numeric(capacity) ||
      length(capacity) != 1L || !is.finite(capacity) || capacity < 2L) {
    return(list())
  }
  candidates <- lapply(value, function(beta) {
    beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      beta, "Cox beta grid row"), error = function(error) numeric())
    if (length(beta) != dimension || any(!is.finite(beta)) ||
        any(abs(beta) > 4) || sum(abs(beta)) > 8) return(NULL)
    loss_bound <- capacity * (log(capacity) + 2 * sum(abs(beta)))
    if (!is.finite(loss_bound) || loss_bound <= 0) return(NULL)
    list(beta = as.numeric(beta), loss_bound = as.numeric(loss_bound))
  })
  if (any(vapply(candidates, is.null, logical(1L)))) return(list())
  keys <- vapply(candidates, function(candidate) .dsvert_joint_dp_client_json(
    as.list(candidate$beta)), character(1L))
  if (anyDuplicated(keys) || !identical(keys, sort(keys, method = "radix"))) {
    return(list())
  }
  candidates
}

.dsvert_dp_cox_partial_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity) {
  artifact <- tryCatch(
    manifest$workload$families$survival_artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "time", "event", "time_grid", "predictors", "predictor_order",
    "intercept", "design_terms", "observation_capacity", "candidate_grid",
    "candidate_order", "candidate_loss_bounds", "numeric_grid_bits",
    "coordinate_count", "coordinate_order", "source_coordinate_scaling",
    "repeated_record_policy", "missingness_policy", "contribution_domain",
    "statistic_maximum", "source_raw_l1_sensitivity",
    "source_raw_l2_sensitivity", "natural_l1_sensitivity",
    "natural_l2_sensitivity", "adjacency", "adjacency_sensitivity_basis",
    "estimation_scope", "implementation_state", "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version, .DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "cox_partial_likelihood_grid_v1") &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) && .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) stop("The signed Cox partial-likelihood descriptor is invalid",
                            call. = FALSE)
  time <- tryCatch(.dsvert_dp_gaussian_bound(artifact$time, "Cox time"),
                   error = function(error) NULL)
  event <- artifact$event
  event_levels <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    event$levels, "Cox event levels"), error = function(error) character())
  event_valid <- is.list(event) && .dsvert_dp_has_exact_names(
    event, c("column", "levels", "censor", "event_level")) &&
    .dsvert_dp_is_string(event$column) && .dsvert_dp_is_string(event$censor) &&
    .dsvert_dp_is_string(event$event_level) && !identical(event$censor,
      event$event_level) && length(event_levels) == 2L && !anyDuplicated(event_levels) &&
    setequal(event_levels, c(event$censor, event$event_level))
  grid <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$time_grid, "Cox time grid"), error = function(error) numeric())
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, "Cox predictor order", sorted = TRUE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order)
  if (isTRUE(predictors_valid) && !is.null(time) && isTRUE(event_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(predictor_order, function(name) {
      .dsvert_dp_gaussian_bound(predictors[[name]], "Cox predictor")
    }), predictor_order), error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(predictors,
      function(value) !value$column %in% c(time$column, event$column), logical(1L)))
  }
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  candidates <- .dsvert_dp_cox_partial_grid_candidates(
    artifact$candidate_grid, length(predictor_order), capacity)
  loss_bounds <- if (length(candidates)) vapply(candidates, `[[`, numeric(1L),
                                                 "loss_bound") else numeric()
  observed_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, "Cox candidate loss bounds"),
    error = function(error) numeric())
  maximum <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, "Cox statistic maxima"),
    error = function(error) numeric())
  raw <- ceiling(loss_bounds * scale)
  multiplier <- if (identical(adjacency, "add_remove_patient")) 1 else if (
      identical(adjacency, "replace_one_fixed_cohort")) 2 else NA_real_
  valid <- !is.null(time) && isTRUE(event_valid) && length(grid) &&
    all(diff(grid) > 0) && time$lower < grid[[1L]] &&
    identical(time$upper, grid[[length(grid)]]) && isTRUE(predictors_valid) &&
    .dsvert_dp_is_integer(bits, 8L, 18L) && identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    length(candidates) && length(candidates) <= 128L &&
    identical(artifact$intercept, FALSE) &&
    identical(artifact$design_terms, predictor_order) &&
    identical(artifact$candidate_order, "canonical_beta_grid_v1") &&
    length(observed_bounds) == length(loss_bounds) &&
    isTRUE(all.equal(observed_bounds, loss_bounds, tolerance = 1e-12)) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, length(candidates),
                          length(candidates)) &&
    identical(artifact$coordinate_order,
      "signed_candidate_grid_breslow_cox_partial_likelihood_loss_v1") &&
    identical(artifact$source_coordinate_scaling,
      "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy,
      "require_one_complete_bounded_row_per_admitted_patient_v1") &&
    identical(artifact$missingness_policy,
      "missing_or_nonfinite_time_or_predictor_or_event_excludes_patient_v1") &&
    identical(artifact$contribution_domain,
      "one_bounded_patient_can_change_the_clipped_cox_cohort_loss_per_signed_candidate_v1") &&
    identical(maximum, raw) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l1_sensitivity),
      multiplier * sum(raw), tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$source_raw_l2_sensitivity),
      multiplier * sqrt(sum(raw^2)), tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l1_sensitivity),
      multiplier * sum(raw) / scale, tolerance = 1e-12)) &&
    isTRUE(all.equal(as.numeric(artifact$natural_l2_sensitivity),
      multiplier * sqrt(sum(raw^2)) / scale, tolerance = 1e-12)) &&
    identical(artifact$adjacency, adjacency) &&
    identical(artifact$adjacency_sensitivity_basis,
      "one_patient_can_change_one_entire_clipped_breslow_cox_loss_by_at_most_its_signed_bound_v1") &&
    identical(artifact$estimation_scope,
      "bounded_same_owner_breslow_cox_partial_likelihood_finite_signed_beta_grid_v1")
  if (!isTRUE(valid)) stop("The signed Cox partial-likelihood descriptor is invalid",
                            call. = FALSE)
  artifact$time <- time
  artifact$event <- list(column = event$column, levels = event_levels,
                         censor = event$censor, event_level = event$event_level)
  artifact$time_grid <- grid
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$candidate_grid <- lapply(candidates, `[[`, "beta")
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$statistic_maximum <- raw
  artifact
}

.dsvert_dp_cox_partial_grid_moment <- function(coordinates, artifact) {
  upper <- as.numeric(artifact$statistic_maximum)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) || any(coordinates < 0) ||
      any(coordinates > upper)) stop("The released Cox grid violates its signed bounds",
                                       call. = FALSE)
  selected <- which.min(coordinates)[[1L]]
  beta <- artifact$candidate_grid[[selected]]
  ranges <- vapply(artifact$predictors, function(value) value$upper - value$lower,
                   numeric(1L))
  coefficients <- beta / ranges
  names(coefficients) <- artifact$predictor_order
  list(coefficients = coefficients, hazard_ratio = exp(coefficients),
       selected_candidate = as.integer(selected),
       selected_dp_partial_loss = coordinates[[selected]] /
         (2^as.numeric(artifact$numeric_grid_bits)))
}

.dsvert_dp_cox_partial_grid_synopsis_release <- function(
    data_name, analysis_id, datasources = NULL, .aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  run <- .dsvert_dp_synopsis_vector_run(datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  metadata <- .dsvert_dp_vector_public_metadata(context)
  count_block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count",
    description = "signed admitted-count capacity block")
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  artifact <- .dsvert_dp_cox_partial_grid_artifact(
    context$manifest, data_name, analysis_id, NULL, context$adjacency,
    as.numeric(context$lattice$output_lattice_scale), capacity)
  blocks <- .dsvert_dp_capsule_vector_blocks(
    context$layout, "survival_artifacts", dataset = data_name,
    owner_peer = artifact$owner_peer)
  blocks <- blocks[vapply(blocks, function(block) identical(
    block$key, analysis_id), logical(1L))]
  signed_descriptor <- tryCatch(context$manifest$workload$families$
    survival_artifacts[[analysis_id]], error = function(error) NULL)
  if (length(blocks) != 1L || !identical(.dsvert_joint_dp_client_json(
      blocks[[1L]]$descriptor), .dsvert_joint_dp_client_json(
        signed_descriptor))) {
    stop("The signed Cox grid does not match its Synopsis layout", call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  moment <- .dsvert_dp_cox_partial_grid_moment(coordinates, artifact)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !identical(verification$artifact$version,
                 .DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION)) {
    stop("The Cox partial-likelihood Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(context = context, metadata = metadata, artifact = verification$artifact,
       block = blocks[[1L]], coordinates = verification$coordinates,
       moment = verification$validated_moment %||% moment,
       certificate = certificate, verification = verification)
}

.dsvert_dp_cox_partial_grid_impl <- function(
    formula, data_name, analysis_id, datasources = NULL, .aggregate) {
  released <- .dsvert_dp_cox_partial_grid_synopsis_release(
    data_name, analysis_id, datasources, .aggregate)
  artifact <- released$artifact
  lhs <- formula[[2L]]
  terms <- stats::terms(formula)
  if (!is.call(lhs) || !identical(as.character(lhs[[1L]]), "Surv") ||
      length(lhs) != 3L || !is.symbol(lhs[[2L]]) || !is.symbol(lhs[[3L]]) ||
      !identical(as.character(lhs[[2L]]), artifact$time$column) ||
      !identical(as.character(lhs[[3L]]), artifact$event$column) ||
      !setequal(attr(terms, "term.labels"), artifact$predictor_order)) {
    stop("formula must match the signed Cox partial-likelihood artifact", call. = FALSE)
  }
  moment <- released$moment
  result <- c(released$metadata, list(
    status = "public_certified_breslow_cox_partial_likelihood_finite_grid",
    analysis_id = analysis_id, family = "cox", coefficients = moment$coefficients,
    hazard_ratio = moment$hazard_ratio, selected_candidate = moment$selected_candidate,
    selected_dp_partial_loss = moment$selected_dp_partial_loss,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    estimand = artifact$estimation_scope,
    covariance = NULL, std_errors = NULL, p_values = NULL, baseline_hazard = NULL,
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    production_ready = FALSE,
    inference = "unavailable_without_protected_cox_score_information_artifact",
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    cross_owner_state = artifact$cross_owner_state,
    legacy_fallback_called = FALSE,
    provenance_certificate = released$certificate,
    disclosure_guard = list(satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("dsvert_dp_cox_partial_grid", "ds.vertCox", "list")
  result
}

.dsvert_dp_cox_partial_grid_adapter <- function(
    explicit_arguments, formula, data, verbose, datasources, analysis_id) {
  allowed <- c("formula", "data", "verbose", "datasources", "analysis_id")
  unsupported <- setdiff(explicit_arguments, allowed)
  if (length(unsupported) || !is.logical(verbose) || length(verbose) != 1L ||
      is.na(verbose) || !inherits(formula, "formula") || is.null(data)) {
    stop(paste(
      "analysis_id requires an explicit additive Surv(time, event) formula,",
      "data, one logical verbose value and no legacy Cox controls"),
      call. = FALSE)
  }
  lhs <- formula[[2L]]
  terms <- stats::terms(formula)
  labels <- attr(terms, "term.labels")
  if (!is.call(lhs) || !identical(as.character(lhs[[1L]]), "Surv") ||
      length(lhs) != 3L || !is.symbol(lhs[[2L]]) || !is.symbol(lhs[[3L]]) ||
      !length(labels) || any(grepl("[:*^|()]", labels))) {
    stop("analysis_id requires Surv(time, event) with additive bare predictors",
         call. = FALSE)
  }
  resolved <- .dsvert_federation_argument(data, datasources)
  result <- .dsvert_dp_cox_partial_grid_impl(
    formula = formula, data_name = resolved$value, analysis_id = analysis_id,
    datasources = resolved$datasources, .aggregate = DSI::datashield.aggregate)
  result$called_via <- "ds.vertCox_analysis_id"
  result
}
