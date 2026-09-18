# Public grouped-grid entry points are separate from the sealed same-owner
# adapters. Only the integration-owned authenticated release reader may supply
# coordinates; no option, callback argument, or plaintext fallback enables it.
.dsvert_dp_grouped_grid_cross_release <- function(contract, datasources) {
  .dsvert_dp_glm_grid_cross_fail()
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
  count <- length(spec$beta_grid)
  if (length(coordinates) != count * width) .dsvert_dp_glm_grid_cross_fail()
  losses <- coordinates[seq.int(1L, length(coordinates), by = width)]
  selected <- which.min(losses)
  beta <- unlist(spec$beta_grid[[selected]], use.names = FALSE)
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
  list(status = "ok", family = spec$family,
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
    # The release reader must bind sticky semantic identity, signed mechanism,
    # authenticated MPC result evidence and both noise-authority contributions.
    # It intentionally fails closed until that integration is promoted.
    released <- .dsvert_dp_grouped_grid_cross_release(contract, datasources)
    .dsvert_dp_grouped_grid_cross_moment(
      released$coordinates, spec, contract$artifact)
  }, error = .dsvert_dp_glm_grid_cross_transcript_stop)
}

#' Select a signed cross-owner grouped finite-grid candidate
#'
#' These typed entry points accept owner-qualified covariates and a contract
#' signed by both compute-and-noise authorities. They never fit an optimiser
#' to protected data. The production reader fails closed until the grouped
#' producer is connected to authenticated sticky joint-DP releases.
#'
#' @param outcome One owner-qualified outcome, such as `site_a$y`.
#' @param predictors Owner-qualified covariates in the signed canonical order.
#' @param contract The complete doubly signed grouped-grid contract.
#' @param policy Public custodian policy with the two pinned authorities.
#' @param schema_manifest Doubly signed aligned dataset schema.
#' @param datasources Optional DataSHIELD connections for the integrated reader.
#' @return After producer integration, the DP-selected fixed-effect candidate,
#'   public fixed parameters and provenance. No standard errors or p-values.
#' @details Random-intercept variance is fixed and signed. GLMM uses the signed
#'   nonadaptive GH5 surrogate, not adaptive quadrature or a fitted variance.
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
