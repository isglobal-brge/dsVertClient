# Internal protected finite L1 paths for binomial and Poisson models.  Every
# beta/lambda pair is custodian-signed; no iterative score worker is reopened.

.DSVERT_CLIENT_DP_LASSO_GRID_ARTIFACT_VERSIONS <- c(
  binomial = "bounded-binomial-lasso-grid-v1",
  poisson = "bounded-poisson-lasso-grid-v1")

.dsvert_dp_lasso_grid_candidates <- function(value, dimension, family) {
  if (!is.list(value) || !length(value) || !is.null(names(value))) return(list())
  candidates <- lapply(value, function(candidate) {
    if (!is.list(candidate) || is.null(names(candidate)) ||
        !setequal(names(candidate), c("lambda", "beta"))) return(NULL)
    lambda <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$lambda, paste(family, "L1 lambda")),
      error = function(error) numeric())
    beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
      candidate$beta, paste(family, "L1 beta row")),
      error = function(error) numeric())
    if (length(lambda) != 1L || anyNA(lambda) || any(!is.finite(lambda)) ||
        lambda < 0 || lambda > 8 || length(beta) != dimension ||
        anyNA(beta) || any(!is.finite(beta)) || any(abs(beta) > 8) ||
        sum(abs(beta)) > 16) return(NULL)
    list(lambda = unname(lambda), beta = unname(beta))
  })
  if (!all(vapply(candidates, is.list, logical(1L)))) return(list())
  candidates
}

.dsvert_dp_lasso_grid_artifact <- function(
    manifest, data_name, analysis_id, owner_peer, adjacency, scale, capacity,
    family) {
  if (!family %in% names(.DSVERT_CLIENT_DP_LASSO_GRID_ARTIFACT_VERSIONS)) {
    stop("The finite L1 family is invalid", call. = FALSE)
  }
  artifact <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts[[analysis_id]],
    error = function(error) NULL)
  required <- c(
    "version", "spec_version", "analysis_id", "dataset", "owner_peer",
    "outcome", "predictors", "predictor_order", "intercept", "design_terms",
    "observation_capacity", "max_outcome", "candidate_grid",
    "penalty_normalizer", "candidate_l1_penalty_raw", "candidate_order",
    "candidate_loss_bounds", "numeric_grid_bits", "coordinate_count",
    "coordinate_order", "source_coordinate_scaling", "repeated_record_policy",
    "missingness_policy", "contribution_domain", "statistic_maximum",
    "source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
    "natural_l1_sensitivity", "natural_l2_sensitivity", "adjacency",
    "adjacency_sensitivity_basis", "estimation_scope", "implementation_state",
    "cross_owner_state")
  basic <- .dsvert_dp_has_exact_names(artifact, required) &&
    identical(artifact$version,
              .DSVERT_CLIENT_DP_LASSO_GRID_ARTIFACT_VERSIONS[[family]]) &&
    identical(artifact$spec_version, paste0(family, "_lasso_grid_v1")) &&
    identical(artifact$analysis_id, analysis_id) &&
    identical(artifact$dataset, data_name) &&
    .dsvert_dp_is_string(artifact$owner_peer) &&
    (is.null(owner_peer) || identical(artifact$owner_peer, owner_peer)) &&
    identical(artifact$implementation_state, "same_owner_materialized") &&
    identical(artifact$cross_owner_state, "reserved_not_materialized")
  if (!isTRUE(basic)) {
    .dsvert_dp_gaussian_reserved(paste0(
      "the signed ", family, " L1 path '", analysis_id,
      "' is unavailable for dataset '", data_name, "'"))
  }
  outcome <- .dsvert_dp_gaussian_bound(
    artifact$outcome, paste(family, "L1 outcome"))
  predictor_order <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$predictor_order, paste(family, "L1 predictor order"), sorted = TRUE),
    error = function(error) character())
  design_terms <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$design_terms, paste(family, "L1 design terms"), sorted = FALSE),
    error = function(error) character())
  predictors <- artifact$predictors
  predictors_valid <- is.list(predictors) && !is.null(names(predictors)) &&
    !anyNA(names(predictors)) && !anyDuplicated(names(predictors)) &&
    identical(names(predictors), predictor_order) && length(predictor_order) &&
    !outcome$column %in% predictor_order
  if (isTRUE(predictors_valid)) {
    predictors <- tryCatch(stats::setNames(lapply(
      predictor_order, function(variable) .dsvert_dp_gaussian_bound(
        predictors[[variable]], paste(family, "L1 predictor"))),
      predictor_order), error = function(error) NULL)
    predictors_valid <- !is.null(predictors) && all(vapply(
      predictors, function(value) !identical(value$column, outcome$column),
      logical(1L)))
  }
  dimension <- 1L + length(predictor_order)
  candidates <- .dsvert_dp_lasso_grid_candidates(
    artifact$candidate_grid, dimension, family)
  lambda <- if (length(candidates)) vapply(
    candidates, function(candidate) candidate$lambda, numeric(1L)) else numeric()
  beta_grid <- if (length(candidates)) lapply(
    candidates, function(candidate) candidate$beta) else list()
  beta_keys <- if (length(beta_grid)) vapply(
    beta_grid, function(beta) .dsvert_joint_dp_client_json(as.list(beta)),
    character(1L)) else character()
  candidate_keys <- if (length(candidates)) vapply(seq_along(candidates),
    function(index) .dsvert_joint_dp_client_json(list(
      lambda = lambda[[index]], beta = beta_grid[[index]])), character(1L)) else
    character()
  canonical <- length(candidates) && !anyDuplicated(candidate_keys) &&
    identical(seq_along(candidates), order(-lambda, beta_keys, method = "radix"))
  poisson <- identical(family, "poisson")
  maximum <- artifact$max_outcome
  outcome_valid <- if (isTRUE(poisson)) {
    .dsvert_dp_is_integer(maximum, 1L, 1024L) &&
      isTRUE(all.equal(outcome$lower, 0, tolerance = 0)) &&
      isTRUE(all.equal(outcome$upper, as.numeric(maximum), tolerance = 0))
  } else is.null(maximum) &&
    isTRUE(all.equal(outcome$lower, 0, tolerance = 0)) &&
    isTRUE(all.equal(outcome$upper, 1, tolerance = 0))
  bits <- suppressWarnings(as.numeric(artifact$numeric_grid_bits))
  loss_bounds <- if (length(candidates) && isTRUE(outcome_valid)) {
    tryCatch(.dsvert_dp_glm_grid_loss_bounds(
      family, beta_grid, if (isTRUE(poisson)) as.integer(maximum) else NULL),
      error = function(error) numeric())
  } else numeric()
  observed_loss_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_loss_bounds, paste(family, "L1 loss bounds")),
    error = function(error) numeric())
  observed_penalties <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$candidate_l1_penalty_raw, paste(family, "L1 penalties")),
    error = function(error) numeric())
  maximums <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$statistic_maximum, paste(family, "L1 maxima")),
    error = function(error) numeric())
  raw_per_candidate <- ceiling(loss_bounds * scale)
  expected_maximum <- capacity * raw_per_candidate
  expected_penalty <- ceiling(capacity * scale * lambda * vapply(
    beta_grid, function(beta) sum(abs(beta[-1L])), numeric(1L)))
  multiplier <- if (identical(adjacency, "replace_one_fixed_cohort")) 2 else if (
    identical(adjacency, "add_remove_patient")) 1 else NA_real_
  raw_l1 <- multiplier * sum(raw_per_candidate)
  raw_l2 <- multiplier * sqrt(sum(raw_per_candidate^2))
  valid <- isTRUE(predictors_valid) && isTRUE(outcome_valid) &&
    .dsvert_dp_is_integer(bits, 8L, 18L) && identical(2^bits, scale) &&
    .dsvert_dp_is_integer(artifact$observation_capacity, capacity, capacity) &&
    length(candidates) >= 2L && length(candidates) <= 256L &&
    length(unique(lambda)) >= 2L && isTRUE(canonical) &&
    isTRUE(artifact$intercept) &&
    identical(design_terms, c("(Intercept)", predictor_order)) &&
    identical(artifact$penalty_normalizer, "observation_capacity_v1") &&
    identical(artifact$candidate_order,
              "canonical_lambda_descending_beta_grid_lasso_v1") &&
    length(observed_loss_bounds) == length(candidates) &&
    isTRUE(all.equal(observed_loss_bounds, loss_bounds, tolerance = 1e-12)) &&
    length(observed_penalties) == length(candidates) &&
    isTRUE(all.equal(observed_penalties, expected_penalty, tolerance = 0)) &&
    all(expected_maximum + expected_penalty <= 2^53) &&
    .dsvert_dp_is_integer(artifact$coordinate_count, length(candidates),
                          length(candidates)) &&
    identical(artifact$coordinate_order, paste(
      "canonical_lambda_descending_beta_grid", family,
      "negative_log_likelihood_lasso_path_v1", sep = "_")) &&
    identical(artifact$source_coordinate_scaling,
              "all_coordinates_already_on_common_numeric_lattice_v1") &&
    identical(artifact$repeated_record_policy, paste(
      "require_one_bounded", family, "outcome_and_mean_once_per_admitted",
      "patient_v1", sep = "_")) &&
    identical(artifact$missingness_policy, paste(
      "noninteger_or_out_of_range_or_missing_outcome_or_missing_or",
      "nonfinite_predictor_excludes_patient_v1", sep = "_")) &&
    identical(artifact$contribution_domain, paste(
      "one_bounded_patient", family, "negative_log_likelihood",
      "contribution_for_every_signed_candidate_v1", sep = "_")) &&
    identical(maximums, expected_maximum) &&
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
      family, "loss_bound_v1", sep = "_")) &&
    identical(artifact$estimation_scope, paste(
      "bounded", family, "fixed_covariates_finite_signed_l1",
      "candidate_path_capacity_normalized_v1", sep = "_"))
  if (!isTRUE(valid)) {
    stop("The signed finite L1 path descriptor is invalid", call. = FALSE)
  }
  artifact$outcome <- outcome
  artifact$predictors <- predictors
  artifact$predictor_order <- predictor_order
  artifact$design_terms <- c("(Intercept)", predictor_order)
  artifact$candidate_grid <- candidates
  artifact$candidate_loss_bounds <- loss_bounds
  artifact$candidate_l1_penalty_raw <- expected_penalty
  artifact$coordinate_count <- as.integer(length(candidates))
  artifact$statistic_maximum <- expected_maximum
  artifact
}

.dsvert_dp_lasso_grid_moment <- function(coordinates, artifact, family) {
  upper <- as.numeric(artifact$statistic_maximum)
  penalties <- as.numeric(artifact$candidate_l1_penalty_raw)
  if (!is.numeric(coordinates) || length(coordinates) != length(upper) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > upper) ||
      anyNA(penalties) || any(!is.finite(penalties)) || any(penalties < 0)) {
    stop("The released finite L1 path violates its signed bounds", call. = FALSE)
  }
  candidates <- artifact$candidate_grid
  lambda <- vapply(candidates, function(candidate) candidate$lambda, numeric(1L))
  selected_lambda <- sort(unique(lambda), decreasing = TRUE)
  source_scale <- 2^as.numeric(artifact$numeric_grid_bits)
  transform_beta <- function(beta) {
    ranges <- vapply(artifact$predictors, function(bound) {
      bound$upper - bound$lower
    }, numeric(1L))
    lowers <- vapply(artifact$predictors, function(bound) bound$lower,
                     numeric(1L))
    slopes <- beta[-1L] / ranges
    names(slopes) <- artifact$predictor_order
    intercept <- beta[[1L]] - sum(slopes * lowers)
    stats::setNames(c(intercept, unname(slopes)),
                    c("(Intercept)", names(slopes)))
  }
  selected <- lapply(selected_lambda, function(value) {
    indexes <- which(lambda == value)
    index <- indexes[[which.min(coordinates[indexes] + penalties[indexes])]]
    list(lambda = value, candidate = as.integer(index),
         coefficients = transform_beta(candidates[[index]]$beta),
         normalized_coefficients = stats::setNames(
           candidates[[index]]$beta, artifact$design_terms),
         dp_negative_log_likelihood = coordinates[[index]] / source_scale,
         l1_penalty = penalties[[index]] / source_scale,
         dp_objective = (coordinates[[index]] + penalties[[index]]) / source_scale)
  })
  names(selected) <- format(selected_lambda, scientific = FALSE, trim = TRUE)
  list(status = "ok", lambda = selected_lambda, selected = selected,
       candidate_selection =
         "minimum_signed_finite_capacity_normalized_l1_dp_objective_v1",
       family = family)
}

.dsvert_dp_lasso_grid_synopsis_release <- function(
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
  artifact <- .dsvert_dp_lasso_grid_artifact(
    context$manifest, data_name, analysis_id, server, context$adjacency,
    scale, capacity, family)
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
    stop("The signed finite L1 path does not match its Synopsis layout",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  moment <- .dsvert_dp_lasso_grid_moment(coordinates, artifact, family)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !identical(verification$artifact$version,
                 .DSVERT_CLIENT_DP_LASSO_GRID_ARTIFACT_VERSIONS[[family]])) {
    stop("The finite L1 Synopsis certificate is not transport-anchored",
         call. = FALSE)
  }
  list(context = context, metadata = metadata, artifact = verification$artifact,
       block = blocks[[1L]], coordinates = verification$coordinates,
       moment = verification$validated_moment %||% moment,
       certificate = certificate, verification = verification)
}

.dsvert_dp_lasso_grid_impl <- function(
    formula, data_name, analysis_id, family, lambda = NULL, server = NULL,
    datasources = NULL, .aggregate) {
  data_name <- .dsvert_dp_gaussian_identifier(data_name, "data_name")
  analysis_id <- .dsvert_dp_gaussian_identifier(analysis_id, "analysis_id")
  outcome <- as.character(formula[[2L]])
  predictors <- attr(stats::terms(formula), "term.labels")
  released <- .dsvert_dp_lasso_grid_synopsis_release(
    data_name, analysis_id, family, server, datasources, .aggregate)
  artifact <- released$artifact
  if (!identical(artifact$outcome$column, outcome) ||
      !setequal(artifact$predictor_order, predictors)) {
    stop("formula must match the signed finite L1 path artifact", call. = FALSE)
  }
  moment <- released$moment
  requested <- if (is.null(lambda)) moment$lambda else {
    if (!is.numeric(lambda) || !length(lambda) || anyNA(lambda) ||
        any(!is.finite(lambda)) || anyDuplicated(lambda) ||
        !all(lambda %in% moment$lambda)) {
      stop("lambda must select unique signed L1 path penalty values exactly",
           call. = FALSE)
    }
    moment$lambda[moment$lambda %in% lambda]
  }
  labels <- format(requested, scientific = FALSE, trim = TRUE)
  selected <- moment$selected[labels]
  result <- c(released$metadata, list(
    status = moment$status, analysis_id = analysis_id,
    cohort_id = released$verification$cohort_id,
    logical_snapshot = released$verification$logical_snapshot,
    certificate_sha256 = released$certificate$certificate_sha256,
    signed_artifact = artifact, server = artifact$owner_peer,
    family = paste0(family, "_finite_lasso_grid"),
    estimand = artifact$estimation_scope, lambda = requested,
    paths = lapply(selected, function(value) value$coefficients),
    normalized_paths = lapply(selected, function(value) value$normalized_coefficients),
    selected_candidates = vapply(selected, function(value) value$candidate,
                                 integer(1L)),
    selected_dp_negative_log_likelihood = vapply(
      selected, function(value) value$dp_negative_log_likelihood, numeric(1L)),
    l1_penalty = vapply(selected, function(value) value$l1_penalty, numeric(1L)),
    selected_dp_objective = vapply(
      selected, function(value) value$dp_objective, numeric(1L)),
    candidate_selection = moment$candidate_selection,
    covariance = NULL, vcov = NULL, std_errors = NULL,
    standard_errors = NULL, p_values = NULL,
    inference = "unavailable_for_finite_l1_path",
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    legacy_iterative_estimand = FALSE,
    cross_owner_state = artifact$cross_owner_state,
    provenance_certificate = released$certificate,
    disclosure_guard = list(satisfied = TRUE,
      basis = "formal_canonical_sticky_DP_synopsis_postprocessing")))
  class(result) <- c("ds.vertLASSOIter", "ds.vertDPLASSOPath", "list")
  result
}

.dsvert_dp_lasso_grid_adapter <- function(
    explicit_arguments, formula, data, family, lambda, verbose, datasources,
    analysis_id) {
  allowed <- c("formula", "data", "family", "lambda", "verbose",
               "datasources", "analysis_id")
  unsupported <- setdiff(explicit_arguments, allowed)
  if (length(unsupported)) {
    stop("The signed finite L1 path does not accept legacy controls: ",
         paste(unsupported, collapse = ", "), call. = FALSE)
  }
  if (!is.character(analysis_id) || length(analysis_id) != 1L ||
      is.na(analysis_id) || !nzchar(analysis_id) ||
      !family %in% c("binomial", "poisson")) {
    stop("The signed finite L1 path requires binomial or Poisson analysis_id",
         call. = FALSE)
  }
  terms <- if (inherits(formula, "formula")) stats::terms(formula) else NULL
  predictors <- if (is.null(terms)) character() else attr(terms, "term.labels")
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]]) || !identical(attr(terms, "intercept"), 1L) ||
      !length(predictors) ||
      any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", predictors))) {
    stop(paste(
      "The signed finite L1 path requires an intercept and additive bare",
      "column names"), call. = FALSE)
  }
  resolved <- .dsvert_federation_argument(data, datasources)
  result <- .dsvert_dp_lasso_grid_impl(
    formula, resolved$value, analysis_id, family, lambda,
    datasources = resolved$datasources, .aggregate = DSI::datashield.aggregate)
  result$called_via <- "ds.vertLASSOIter_analysis_id"
  result
}

