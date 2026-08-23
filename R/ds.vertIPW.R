#' @title Sticky intercept-only inverse-probability weighting
#' @description Computes the binary, intercept-only ATE IPW identity from one
#' signed sticky DP contingency release. With an intercept-only propensity
#' model, IPW is exactly the treated-minus-control risk difference, so no
#' weights, source rows, propensity fit, or second private computation is
#' exposed. Covariate-adjusted IPW, ATT/ATC, continuous outcomes and weighted
#' outcome regression remain unavailable.
#'
#' @param outcome_formula A formula exactly of the form code{outcome ~ treatment}.
#' @param propensity_formula A formula exactly of the form
#'   code{treatment ~ 1}.
#' @param data Protected data frame name or a version-2 code{ds.vertFederation}.
#' @param weights_column Retained for compatibility and must be code{"ipw"}.
#' @param outcome_family Must be code{"binomial"}.
#' @param treated Explicit treated level from the signed binary treatment domain.
#' @param event Event outcome column name or index in the signed binary outcome
#'   domain.
#' @param level DP-mechanism and combined-region coverage level.
#' @param server Optional same-owner assertion for the signed table artifact.
#' @param verbose Print a short route message.
#' @param datasources DataSHIELD connections.
#' @param ... No additional controls are accepted.
#' @return A code{ds.vertIPW} object with the intercept-only ATE, conservative
#'   DP-mechanism and combined DP-plus-sampling regions, and no released
#'   individual weights or propensity fit.
#' @details This is the exact intercept-only IPW special case, not a
#' covariate-adjusted IPW implementation. The release is one canonical
#' treatment-by-outcome DP table; unlimited retries return the same sticky
#' release and do not create another draw.
#' @seealso code{link{ds.vertDPCausalStandardization}}
#' @export
ds.vertIPW <- function(outcome_formula, propensity_formula, data = NULL,
                      weights_column = "ipw", outcome_family = "binomial",
                      treated = NULL, event = 2L, level = 0.95,
                      server = NULL, verbose = TRUE, datasources = NULL, ...) {
  dots <- list(...)
  if (length(dots)) {
    stop("intercept-only IPW does not accept additional controls", call. = FALSE)
  }
  if (!identical(weights_column, "ipw")) {
    stop("weights_column must be 'ipw' for intercept-only IPW", call. = FALSE)
  }
  if (!identical(outcome_family, "binomial")) {
    stop("intercept-only IPW requires outcome_family = 'binomial'", call. = FALSE)
  }
  if (!is.character(treated) || length(treated) != 1L || is.na(treated) ||
      !nzchar(treated)) {
    stop("treated must name one level in the signed treatment domain", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be in (0, 1)", call. = FALSE)
  }
  if (is.null(data)) {
    stop("data must select one signed treatment-outcome artifact", call. = FALSE)
  }
  spec <- .dsvert_ipw_intercept_only_spec(
    outcome_formula, propensity_formula)
  table <- ds.vertDPContingency(
    data, spec$treatment, spec$outcome, server = server,
    datasources = datasources)
  treatment <- table$row_levels
  if (!is.character(treatment) || length(treatment) != 2L || anyNA(treatment) ||
      any(!nzchar(treatment)) || anyDuplicated(treatment) ||
      !treated %in% treatment) {
    stop("The signed IPW treatment table must have the requested two-level domain",
         call. = FALSE)
  }
  causal_arguments <- list(
    x = table, strata = rep("overall", 2L), treatment = treatment,
    treated = treated, standard_weights = c(overall = 1), event = event,
    level = level)
  causal <- do.call(ds.vertDPCausalStandardization, causal_arguments)
  combined <- do.call(ds.vertDPCausalStandardizationInference, causal_arguments)
  if (isTRUE(verbose)) {
    message("[ds.vertIPW] intercept-only ATE from one signed sticky DP table")
  }
  out <- list(
    status = causal$status,
    estimand = "ATE",
    estimate = causal$point_estimates$risk_difference,
    risk_treated = causal$point_estimates$risk_treated,
    risk_control = causal$point_estimates$risk_control,
    mechanism_region = causal$mechanism_regions$risk_difference,
    confidence_region = combined$combined_regions$risk_difference,
    coverage_level = level,
    uncertainty_scope = combined$uncertainty_scope,
    propensity_model = "intercept_only",
    outcome_model = "saturated_binary_treatment_means",
    weights_released = FALSE,
    result_contract = "intercept_only_ipw_equivalence_from_sticky_dp_table_v1",
    source_artifact_key = table$artifact_key,
    final_vector_root = table$final_vector_root,
    sticky_replay = TRUE,
    server_calls_for_artifact = 1L,
    additional_server_calls_after_artifact = 0L,
    additional_privacy_cost_after_artifact = c(epsilon = 0, delta = 0),
    limitations = paste(
      "No covariates, ATT, ATC, propensity weights, outcome regression,",
      "standard errors or p-values are available."))
  class(out) <- c("ds.vertIPW", "list")
  out
}

.dsvert_ipw_intercept_only_spec <- function(outcome_formula,
                                             propensity_formula) {
  symbol <- function(value, label) {
    if (!is.symbol(value)) {
      stop(label, " must be one untransformed column name", call. = FALSE)
    }
    as.character(value)
  }
  valid_formula <- function(value, label) {
    if (!inherits(value, "formula") || length(value) != 3L) {
      stop(label, " must be a two-sided formula", call. = FALSE)
    }
  }
  valid_formula(outcome_formula, "outcome_formula")
  valid_formula(propensity_formula, "propensity_formula")
  outcome <- symbol(outcome_formula[[2L]], "outcome_formula response")
  treatment <- symbol(outcome_formula[[3L]], "outcome_formula")
  propensity_treatment <- symbol(
    propensity_formula[[2L]], "propensity_formula response")
  intercept <- propensity_formula[[3L]]
  if (!identical(propensity_treatment, treatment) ||
      !is.numeric(intercept) || length(intercept) != 1L ||
      is.na(intercept) || !identical(as.numeric(intercept), 1)) {
    stop(paste(
      "intercept-only IPW requires outcome ~ treatment and",
      "the matching propensity_formula treatment ~ 1"), call. = FALSE)
  }
  if (identical(outcome, treatment)) {
    stop("outcome and treatment must be distinct columns", call. = FALSE)
  }
  list(outcome = outcome, treatment = treatment)
}

#' @export
print.ds.vertIPW <- function(x, ...) {
  cat("dsVert intercept-only IPW ATE from a sticky DP table\n")
  cat("  risk difference = ", format(x$estimate, digits = 6L), "\n", sep = "")
  if (is.numeric(x$confidence_region) && length(x$confidence_region) == 2L) {
    cat("  combined region = [", format(x$confidence_region[[1L]], digits = 6L),
        ", ", format(x$confidence_region[[2L]], digits = 6L), "]\n", sep = "")
  }
  cat("  No individual weights, covariates, standard errors or p-values are released.\n")
  invisible(x)
}
