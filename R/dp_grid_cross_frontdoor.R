# Public front doors discover only contracts already embedded in an authenticated
# Synopsis. Discovery never creates an unsigned candidate set or opens a source.
.dsvert_dp_grid_cross_formula_requested <- function(formula) {
  inherits(formula, "formula") && any(all.names(formula) == "$")
}

.dsvert_dp_grid_cross_formula <- function(formula, spec = NULL) {
  reference <- function(x) {
    if (!is.call(x) || length(x) != 3L || !identical(x[[1L]], as.name("$")) ||
        !is.symbol(x[[2L]]) || !is.symbol(x[[3L]])) {
      stop("Cross-owner grids require additive owner-qualified references", call. = FALSE)
    }
    paste0(as.character(x[[2L]]), "$", as.character(x[[3L]]))
  }
  terms <- function(x) {
    if (is.call(x) && identical(x[[1L]], as.name("+")) && length(x) == 3L) {
      return(c(terms(x[[2L]]), terms(x[[3L]])))
    }
    reference(x)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("Cross-owner grids require an additive formula", call. = FALSE)
  }
  outcome <- reference(formula[[2L]])
  predictors <- terms(formula[[3L]])
  if (anyDuplicated(predictors) || outcome %in% predictors) {
    stop("Cross-owner grids require distinct outcome and predictor references", call. = FALSE)
  }
  predictors <- .dsvert_dp_glm_grid_cross_references(predictors)
  if (!is.null(spec) && (!identical(outcome, spec$outcome$reference) ||
      !setequal(predictors, unlist(spec$predictor_order, use.names = FALSE)))) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  list(outcome = outcome, predictors = predictors)
}

.dsvert_dp_grid_cross_bootstrap <- function(datasources,
    .aggregate = DSI::datashield.aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  bootstrap <- .dsvert_dp_synopsis_bootstrap_build_v1(datasources, .aggregate = .aggregate)
  trusted <- .dsvert_dp_synopsis_client_bundle(bootstrap$manifest_bundle, bootstrap$status)
  manifest <- trusted$manifest
  policy <- list(peer_pinset = trusted$context$pinset,
    peer_pinset_sha256 = .dsvert_dp_capsule_source_hash(as.list(trusted$context$pinset)),
    designated_noise_peers = trusted$context$designated,
    unit_capacity = manifest$admission$unit_capacity,
    numeric_grid_bits = manifest$bounds$numeric_grid_bits,
    adjacency = manifest$admission$adjacency)
  schema <- .dsvert_joint_dp_client_decode(bootstrap$manifest_bundle$schema_json,
    "signed cross-grid frontdoor schema", .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES)
  list(manifest = manifest, policy = policy, schema = schema, datasources = datasources)
}

.dsvert_dp_grid_cross_frontdoor <- function(family, formula, data, analysis_id,
    datasources, classes = NULL, reference = NULL, levels_ordered = NULL) {
  .dsvert_dp_grid_cross_formula(formula)
  resolved <- .dsvert_federation_argument(data, datasources)
  trusted <- .dsvert_dp_grid_cross_bootstrap(resolved$datasources)
  artifact <- trusted$manifest$workload$families$gaussian_models$artifacts[[analysis_id]]
  if (!is.list(artifact)) .dsvert_dp_glm_grid_cross_fail()
  contract <- .dsvert_dp_glm_grid_profile_admit(
    .dsvert_dp_glm_grid_cross_embedded_contract(artifact), trusted$policy, trusted$schema)
  spec <- contract$spec
  if (!identical(spec$family, family) || !identical(spec$analysis_id, analysis_id) ||
      !identical(spec$dataset, resolved$value)) .dsvert_dp_glm_grid_cross_fail()
  .dsvert_dp_glm_grid_cross_client_artifact(artifact, resolved$value, analysis_id,
    NULL, trusted$policy$adjacency, 2^trusted$policy$numeric_grid_bits,
    trusted$policy$unit_capacity, family)
  .dsvert_dp_grid_cross_formula(formula, spec)
  if ((!is.null(classes) && !identical(classes, unlist(spec$class_order, use.names = FALSE))) ||
      (!is.null(reference) && !identical(reference, spec$reference)) ||
      (!is.null(levels_ordered) &&
       !identical(levels_ordered, unlist(spec$class_order, use.names = FALSE)))) {
    stop("Categorical controls must exactly match the signed cross-owner domain", call. = FALSE)
  }
  switch(family,
    nb = dp_nb_grid(formula, resolved$value, analysis_id, contract,
      trusted$policy, trusted$schema, trusted$datasources),
    multinomial = dp_multinomial_grid(formula, resolved$value, analysis_id, contract,
      trusted$policy, trusted$schema, trusted$datasources),
    ordinal = dp_ordinal_grid(formula, resolved$value, analysis_id, contract,
      trusted$policy, trusted$schema, trusted$datasources),
    .dsvert_dp_glm_grid_cross_fail())
}

.dsvert_dp_lasso_cross_frontdoor <- function(explicit_arguments, formula, data,
    family, lambda, datasources, analysis_id, cross_contract, base_contract) {
  allowed <- c("formula", "data", "family", "lambda", "verbose", "datasources",
    "analysis_id", "cross_contract", "base_contract")
  if (length(setdiff(explicit_arguments, allowed))) {
    stop("The signed cross-owner L1 path does not accept legacy controls", call. = FALSE)
  }
  if (is.null(cross_contract) || is.null(base_contract) ||
      !family %in% c("binomial", "poisson") ||
      !.dsvert_dp_is_string(analysis_id)) {
    stop("Cross-owner LASSO requires binomial/Poisson analysis_id and both signed contracts",
      call. = FALSE)
  }
  refs <- .dsvert_dp_grid_cross_formula(formula)
  resolved <- .dsvert_federation_argument(data, datasources)
  trusted <- .dsvert_dp_grid_cross_bootstrap(resolved$datasources)
  base <- .dsvert_dp_glm_grid_profile_admit(base_contract, trusted$policy, trusted$schema)
  validated <- .dsvert_dp_lasso_cross_contract_validate(cross_contract,
    trusted$policy, trusted$schema, base,
    .base_validator = function(value, policy, schema) {
      .dsvert_dp_lasso_cross_base(.dsvert_dp_glm_grid_profile_admit(value, policy, schema), policy)
    })
  if (!identical(validated$spec$analysis_id, analysis_id) ||
      !identical(base$spec$family, family) ||
      !identical(base$spec$dataset, resolved$value)) .dsvert_dp_glm_grid_cross_fail()
  .dsvert_dp_grid_cross_formula(formula, base$spec)
  signed_lambda <- unique(vapply(validated$spec$candidate_grid, `[[`, numeric(1L), "lambda"))
  if (!is.null(lambda) && (!is.numeric(lambda) || !length(lambda) || anyNA(lambda) ||
      any(!is.finite(lambda)) || anyDuplicated(lambda) || any(!lambda %in% signed_lambda))) {
    stop("lambda must select unique exact signed penalty levels", call. = FALSE)
  }
  result <- dp_lasso_grid(refs$predictors, validated, base,
    trusted$policy, trusted$schema, trusted$datasources)
  if (!is.null(lambda)) {
    keep <- which(result$lambda %in% lambda)
    for (field in c("lambda", "selected_candidates", "selected_objectives",
                    "normalized_paths", "paths")) result[[field]] <- result[[field]][keep]
  }
  labels <- format(result$lambda, scientific = FALSE, trim = TRUE)
  names(result$paths) <- names(result$normalized_paths) <- labels
  result$status <- "ok"
  result$selected_dp_objective <- result$selected_objectives
  result$inference <- "unavailable_for_finite_l1_path"
  result$source_values_exposed <- FALSE
  result$intermediate_values_exposed <- FALSE
  result$lambda_grid <- result$lambda
  result$coefficients <- result$paths
  result$analysis_id <- analysis_id
  result$data <- resolved$value
  result$method <- "signed_dp_cross_owner_finite_l1_path"
  result$called_via <- "ds.vertLASSOIter_cross_contract"
  result$legacy_iterative_estimand <- FALSE
  result$additional_server_calls_after_synopsis <- 0L
  class(result) <- c("ds.vertLASSOIter", "ds.vertDPLASSOPath", "list")
  result
}
