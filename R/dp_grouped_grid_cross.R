# Public grouped-grid entry points are separate from the sealed same-owner
# adapters. Only the integration-owned authenticated release reader may supply
# coordinates; no option, callback argument, or plaintext fallback enables it.
.dsvert_dp_grouped_grid_cross_release <- function(contract, policy, schema, datasources,
    .aggregate = DSI::datashield.aggregate) {
  if (!contract$spec$family %in% c("lmm", "binomial_glmm", "poisson_glmm")) .dsvert_dp_glm_grid_cross_fail()
  contract <- .dsvert_dp_glm_grid_profile_admit(contract, policy, schema)
  datasources <- .dsvert_dp_datasources(datasources)
  bootstrap <- .dsvert_dp_synopsis_bootstrap_build_v1(datasources, .aggregate = .aggregate)
  trusted <- .dsvert_dp_synopsis_client_bundle(bootstrap$manifest_bundle, bootstrap$status)
  .dsvert_dp_glm_grid_cross_equal(as.list(policy$peer_pinset), as.list(trusted$context$pinset))
  if (!setequal(unname(policy$designated_noise_peers), trusted$context$designated)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  .dsvert_dp_glm_grid_cross_equal(schema,
    .dsvert_joint_dp_client_decode(bootstrap$manifest_bundle$schema_json,
      "signed grouped-grid schema", .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES))
  spec <- contract$spec
  check <- function(manifest) {
    artifact <- manifest$workload$families$gaussian_models$artifacts[[spec$analysis_id]]
    if (!is.list(artifact)) .dsvert_dp_glm_grid_cross_fail()
    .dsvert_dp_glm_grid_cross_equal(.dsvert_dp_glm_grid_cross_embedded_contract(artifact), contract)
    .dsvert_dp_grouped_cross_client_artifact(artifact, spec$dataset,
      spec$analysis_id, NULL, manifest$admission$adjacency,
      2^manifest$bounds$numeric_grid_bits, manifest$admission$unit_capacity)
  }
  check(trusted$manifest)
  run <- .dsvert_dp_synopsis_vector_run(datasources, .aggregate = .aggregate,
    .request_check = check)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  artifact <- check(context$manifest)
  blocks <- .dsvert_dp_capsule_vector_blocks(context$layout, "gaussian_models",
    dataset = spec$dataset, owner_peer = artifact$owner_peer)
  blocks <- blocks[vapply(blocks, function(x) identical(x$key, spec$analysis_id), logical(1L))]
  if (length(blocks) != 1L) .dsvert_dp_glm_grid_cross_fail()
  .dsvert_dp_glm_grid_cross_equal(blocks[[1L]]$descriptor,
    context$manifest$workload$families$gaussian_models$artifacts[[spec$analysis_id]])
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !identical(as.numeric(verification$output_lattice_scale), 2^spec$numeric_grid_bits)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  result <- .dsvert_dp_grouped_grid_cross_moment(
    unname(verification$coordinates * verification$output_lattice_scale), spec, contract$artifact)
  c(.dsvert_dp_vector_public_metadata(context), result, list(
    provenance_certificate = certificate, certificate_sha256 = certificate$certificate_sha256,
    signed_contract = contract))
}

.dsvert_dp_grouped_grid_cross_moment <- function(coordinates, spec, artifact) {
  .dsvert_dp_grouped_grid_cross_artifact_validate(artifact, spec)
  maxima <- unlist(spec$sensitivity$maximum_coordinates, use.names = FALSE)
  if (!is.numeric(coordinates) || length(coordinates) != length(maxima) ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates != floor(coordinates)) || any(coordinates < 0) ||
      any(coordinates > maxima)) .dsvert_dp_glm_grid_cross_fail()
  dimension <- length(spec$predictor_order) + 1L
  triangle <- dimension * (dimension + 1L) / 2L
  gee <- grepl("_gee$", spec$family)
  width <- if (gee) 1L + 2L * triangle else 1L
  ml <- identical(spec$family, "lmm") && identical(spec$parameters$objective, "ml")
  glmm_grid <- spec$family %in% c("binomial_glmm", "poisson_glmm") &&
    !is.null(spec$parameters$variance_grid)
  count <- if (ml || glmm_grid) length(spec$candidate_grid) else length(spec$beta_grid)
  if (length(coordinates) != count * width) .dsvert_dp_glm_grid_cross_fail()
  losses <- coordinates[seq.int(1L, length(coordinates), by = width)]
  selected <- which.min(losses)
  candidate <- if (ml || glmm_grid) spec$candidate_grid[[selected]] else list(beta_index = selected)
  beta <- unlist(spec$beta_grid[[candidate$beta_index]], use.names = FALSE)
  predictors <- unlist(spec$predictor_order, use.names = FALSE)
  spans <- vapply(spec$predictors, function(x) x$upper - x$lower, numeric(1L))
  lower <- vapply(spec$predictors, `[[`, numeric(1L), "lower")
  slopes <- beta[-1L] / spans
  coefficients <- stats::setNames(c(beta[1L] - sum(slopes * lower), slopes),
                                  c("(Intercept)", predictors))
  if (identical(spec$family, "lmm")) {
    coefficients <- coefficients * (spec$outcome$upper - spec$outcome$lower)
    coefficients[1L] <- coefficients[1L] + spec$outcome$lower
  }
  result <- list(status = "ok", family = spec$family,
       analysis_id = spec$analysis_id, coefficients = coefficients,
       normalized_coefficients = stats::setNames(beta,
         c("(Intercept)", predictors)), selected_candidate = as.integer(selected),
       selected_dp_loss = losses[[selected]] / 2^spec$numeric_grid_bits,
       candidate_selection = "minimum_signed_finite_grid_dp_postprocessing_v1",
       parameters = spec$parameters, dp_unit = "patient",
       inference = "no_standard_errors_or_p_values_v1",
       standard_errors = NULL, p_values = NULL,
       implementation_state = artifact$implementation_state,
       cross_owner_state = artifact$cross_owner_state)
  if (ml) {
    result$parameters <- c(list(objective = "ml"),
      spec$parameters$variance_grid[[candidate$variance_index]])
    result$loss_objective <- spec$numeric_contract$objective
  }
  if (glmm_grid) {
    result$parameters <- list(
      random_intercept_variance = spec$parameters$variance_grid[[candidate$variance_index]],
      quadrature = spec$parameters$quadrature)
    result$loss_objective <- spec$numeric_contract$objective
  }
  if (identical(spec$family, "poisson_glmm"))
    result$loss_objective <- "factorial_free_gh5_selection_plus_2_per_live_row_v1"
  result
}

.dsvert_dp_grouped_grid_cross_impl <- function(
    families, outcome, predictors, contract, policy, schema_manifest,
    datasources = NULL) {
  tryCatch({
    contract <- .dsvert_dp_grouped_grid_cross_contract_validate(
      contract, policy, schema_manifest)
    spec <- contract$spec
    if (!spec$family %in% families || !is.character(outcome) ||
        length(outcome) != 1L || !identical(outcome, spec$outcome$reference) ||
        !identical(.dsvert_dp_glm_grid_cross_references(predictors),
                   unlist(spec$predictor_order, use.names = FALSE))) {
      .dsvert_dp_glm_grid_cross_fail()
    }
    .dsvert_dp_grouped_grid_cross_release(contract, policy, schema_manifest, datasources)
  }, error = .dsvert_dp_glm_grid_cross_transcript_stop)
}

#' Select a signed cross-owner grouped finite-grid candidate
#'
#' These typed entry points accept owner-qualified covariates and a contract
#' signed by both compute-and-noise authorities. They never fit an optimiser
#' to protected data. LMM ML and binomial GH5 grids read authenticated sticky
#' joint-DP releases; the remaining grouped readers fail closed.
#'
#' @param outcome One owner-qualified outcome, such as `site_a$y`.
#' @param predictors Owner-qualified covariates in the signed canonical order.
#' @param contract The complete doubly signed grouped-grid contract.
#' @param policy Public custodian policy with the two pinned authorities.
#' @param schema_manifest Doubly signed aligned dataset schema.
#' @param datasources Optional DataSHIELD connections for the integrated reader.
#' @return The DP-selected fixed-effect candidate,
#'   selected signed parameters and provenance. No standard errors or p-values.
#' @details LMM accepts a signed fixed covariance or an explicit finite ML
#'   variance grid, ordered by variance then coefficient candidate. The ML
#'   objective includes the private log determinant and a candidate-independent
#'   count shift; it does not provide REML. The authenticated LMM reader requires
#'   an explicit ML grid. GLMM uses signed fixed variance or, for binomial and Poisson,
#'   an explicit variance grid over zero and one quarter, ordered by variance
#'   then coefficient candidate, with the nonadaptive GH5 surrogate.
#'   The authenticated binomial and Poisson GLMM readers require an explicit variance grid.
#'   GEE selects using independent likelihood; its complete DP workload also
#'   includes bread and clipped cluster-score meat. Its correlation structure
#'   and parameter are signed. All families protect one admitted patient with
#'   exactly one analysis row; the entire affected cluster enters the bound.
#' @name dp_grouped_grid_cross
NULL

#' @rdname dp_grouped_grid_cross
#' @export
dp_lmm_grid <- function(outcome, predictors, contract, policy, schema_manifest,
                        datasources = NULL) {
  .dsvert_dp_grouped_grid_cross_impl("lmm", outcome, predictors, contract,
                                    policy, schema_manifest, datasources)
}

#' @rdname dp_grouped_grid_cross
#' @export
dp_glmm_grid <- function(outcome, predictors, contract, policy, schema_manifest,
                         datasources = NULL) {
  .dsvert_dp_grouped_grid_cross_impl(c("binomial_glmm", "poisson_glmm"),
    outcome, predictors, contract, policy, schema_manifest, datasources)
}

#' @rdname dp_grouped_grid_cross
#' @export
dp_gee_grid <- function(outcome, predictors, contract, policy, schema_manifest,
                        datasources = NULL) {
  .dsvert_dp_grouped_grid_cross_impl(c("binomial_gee", "poisson_gee"),
    outcome, predictors, contract, policy, schema_manifest, datasources)
}

# One integration surface; obtaining these callbacks does not promote a reader.
.dsvert_dp_grouped_grid_cross_register <- function() {
  list(versions = .DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_SPEC_VERSIONS,
       validate = .dsvert_dp_grouped_grid_cross_contract_validate,
       artifact_validate = .dsvert_dp_grouped_grid_cross_artifact_validate,
       postprocess = .dsvert_dp_grouped_grid_cross_moment,
       release = .dsvert_dp_grouped_grid_cross_release)
}
