#' @title Sticky saturated inverse-probability weighting
#' @description Computes a binary ATE from one signed sticky DP contingency
#' release. With an intercept-only propensity model, IPW is exactly the
#' treated-minus-control risk difference. With one categorical propensity
#' stratum, a signed stratum-treatment-by-outcome table and fixed public target
#' weights, the saturated IPW identity is exactly the corresponding
#' stratum-standardised g-formula. Neither route releases weights, source rows,
#' a propensity fit, or performs a second private computation. ATT/ATC,
#' continuous outcomes and weighted outcome regression remain unavailable.
#'
#' @param outcome_formula A formula exactly of the form
#'   \code{outcome ~ treatment}.
#' @param propensity_formula A formula exactly of the form \code{treatment ~ 1}
#'   or \code{treatment ~ stratum} for one categorical saturated stratum.
#' @param data Protected data frame name or a version-2
#'   \code{ds.vertFederation}.
#' @param weights_column Retained for compatibility and must be \code{"ipw"}.
#' @param outcome_family Must be \code{"binomial"}.
#' @param treated Explicit treated level from the signed binary treatment domain.
#' @param event Event outcome column name or index in the signed binary outcome
#'   domain.
#' @param arm_column For \code{treatment ~ stratum}, the signed categorical
#'   column whose levels are the complete public stratum-treatment rows.
#' @param arm_strata,arm_treatment Named public vectors, with names exactly the
#'   signed \code{arm_column} levels, assigning every row to its stratum and
#'   treatment level respectively.
#' @param standard_weights Named, non-negative public target-population weights
#'   for every declared stratum. Required for \code{treatment ~ stratum}.
#' @param level DP-mechanism and combined-region coverage level.
#' @param server Optional same-owner assertion for the signed table artifact.
#' @param verbose Print a short route message.
#' @param datasources DataSHIELD connections.
#' @param ... No additional controls are accepted.
#' @return A \code{ds.vertIPW} object with a conservative DP-mechanism and
#'   combined DP-plus-sampling region, and no released individual weights or
#'   propensity fit.
#' @details The categorical stratum route is an exact saturated IPW/g-formula
#' identity, not a fitted propensity-score model. Its release is one canonical
#' stratum-treatment-by-outcome DP table; unlimited retries return the same
#' sticky release and do not create another draw. Arbitrary covariates,
#' continuous strata, ATT/ATC, propensity fitting, trimming and outcome
#' regression remain unavailable.
#' @seealso \code{\link{ds.vertDPCausalStandardization}}
#' @export
ds.vertIPW <- function(outcome_formula, propensity_formula, data = NULL,
                      weights_column = "ipw", outcome_family = "binomial",
                      treated = NULL, event = 2L, level = 0.95,
                      arm_column = NULL, arm_strata = NULL,
                      arm_treatment = NULL, standard_weights = NULL,
                      server = NULL, verbose = TRUE, datasources = NULL, ...) {
  dots <- list(...)
  if (length(dots)) {
    stop("saturated IPW does not accept additional controls", call. = FALSE)
  }
  if (!identical(weights_column, "ipw")) {
    stop("weights_column must be 'ipw' for saturated IPW", call. = FALSE)
  }
  if (!identical(outcome_family, "binomial")) {
    stop("saturated IPW requires outcome_family = 'binomial'", call. = FALSE)
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
    stop("data must select one signed IPW contingency artifact", call. = FALSE)
  }
  spec <- .dsvert_ipw_formula_spec(
    outcome_formula, propensity_formula)
  if (identical(spec$route, "intercept_only")) {
    if (!is.null(arm_column) || !is.null(arm_strata) ||
        !is.null(arm_treatment) || !is.null(standard_weights)) {
      stop(paste(
        "intercept-only IPW does not accept arm_column, arm_strata,",
        "arm_treatment or standard_weights"), call. = FALSE)
    }
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
    propensity_model <- "intercept_only"
    outcome_model <- "saturated_binary_treatment_means"
    result_contract <- "intercept_only_ipw_equivalence_from_sticky_dp_table_v1"
    strata_variable <- NULL
  } else {
    .dsvert_ipw_validate_stratified_controls(
      arm_column, arm_strata, arm_treatment, standard_weights)
    if (arm_column %in% c(spec$outcome, spec$treatment, spec$stratum)) {
      stop(paste(
        "saturated categorical IPW arm_column must differ from outcome,",
        "treatment and stratum"), call. = FALSE)
    }
    table <- ds.vertDPContingency(
      data, arm_column, spec$outcome, server = server,
      datasources = datasources)
    mapping <- .dsvert_ipw_stratified_mapping(
      table$row_levels, arm_strata, arm_treatment, standard_weights, treated)
    causal_arguments <- list(
      x = table, strata = mapping$strata, treatment = mapping$treatment,
      treated = treated, standard_weights = mapping$standard_weights,
      event = event, level = level)
    propensity_model <- "saturated_categorical_treatment_given_stratum"
    outcome_model <- "saturated_stratum_treatment_binary_means"
    result_contract <-
      "saturated_categorical_ipw_g_formula_equivalence_from_sticky_dp_table_v1"
    strata_variable <- spec$stratum
  }
  causal <- do.call(ds.vertDPCausalStandardization, causal_arguments)
  combined <- do.call(ds.vertDPCausalStandardizationInference, causal_arguments)
  if (isTRUE(verbose)) {
    message(if (identical(spec$route, "intercept_only")) {
      "[ds.vertIPW] intercept-only ATE from one signed sticky DP table"
    } else {
      "[ds.vertIPW] saturated categorical stratum ATE from one signed sticky DP table"
    })
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
    propensity_model = propensity_model,
    outcome_model = outcome_model,
    strata_variable = strata_variable,
    standard_weights = causal_arguments$standard_weights,
    weights_released = FALSE,
    result_contract = result_contract,
    source_artifact_key = table$artifact_key,
    final_vector_root = table$final_vector_root,
    sticky_replay = TRUE,
    server_calls_for_artifact = 1L,
    additional_server_calls_after_artifact = 0L,
    additional_privacy_cost_after_artifact = c(epsilon = 0, delta = 0),
    limitations = paste(
      "Only the intercept-only or one categorical saturated-stratum ATE",
      "identity is available. No ATT, ATC, individual propensity weights,",
      "continuous or multiple covariates, outcome regression, standard errors",
      "or p-values are available."))
  class(out) <- c("ds.vertIPW", "list")
  out
}

.dsvert_ipw_formula_spec <- function(outcome_formula, propensity_formula) {
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
  propensity_rhs <- propensity_formula[[3L]]
  route <- if (is.numeric(propensity_rhs) && length(propensity_rhs) == 1L &&
      !is.na(propensity_rhs) && identical(as.numeric(propensity_rhs), 1)) {
    "intercept_only"
  } else if (is.symbol(propensity_rhs)) {
    "saturated_stratum"
  } else {
    NA_character_
  }
  if (!identical(propensity_treatment, treatment) || is.na(route)) {
    stop(paste(
      "IPW requires outcome ~ treatment and the matching propensity_formula",
      "treatment ~ 1 or treatment ~ stratum"), call. = FALSE)
  }
  if (identical(outcome, treatment)) {
    stop("outcome and treatment must be distinct columns", call. = FALSE)
  }
  stratum <- if (identical(route, "saturated_stratum")) {
    as.character(propensity_rhs)
  } else NULL
  if (!is.null(stratum) && (identical(stratum, outcome) ||
      identical(stratum, treatment))) {
    stop("the IPW stratum must differ from outcome and treatment", call. = FALSE)
  }
  list(route = route, outcome = outcome, treatment = treatment, stratum = stratum)
}

#' @keywords internal
.dsvert_ipw_validate_stratified_controls <- function(
    arm_column, arm_strata, arm_treatment, standard_weights) {
  if (!is.character(arm_column) || length(arm_column) != 1L ||
      is.na(arm_column) || !nzchar(arm_column)) {
    stop("saturated categorical IPW requires one arm_column", call. = FALSE)
  }
  for (value in list(arm_strata, arm_treatment)) {
    if (!is.character(value) || length(value) < 4L || anyNA(value) ||
        any(!nzchar(value)) || is.null(names(value)) ||
        anyNA(names(value)) || any(!nzchar(names(value))) ||
        anyDuplicated(names(value))) {
      stop(paste(
        "saturated categorical IPW requires named arm_strata and",
        "arm_treatment vectors"), call. = FALSE)
    }
  }
  if (!identical(names(arm_strata), names(arm_treatment)) ||
      length(unique(arm_strata)) < 2L || length(unique(arm_treatment)) != 2L ||
      any(table(arm_strata, arm_treatment) != 1L)) {
    stop(paste(
      "saturated categorical IPW requires one complete binary treatment",
      "pair for every named stratum"), call. = FALSE)
  }
  if (!is.numeric(standard_weights) || length(standard_weights) < 2L ||
      anyNA(standard_weights) || any(!is.finite(standard_weights)) ||
      any(standard_weights < 0) || sum(standard_weights) <= 0 ||
      is.null(names(standard_weights)) || anyNA(names(standard_weights)) ||
      any(!nzchar(names(standard_weights))) || anyDuplicated(names(standard_weights))) {
    stop(paste(
      "saturated categorical IPW requires named non-negative",
      "standard_weights for at least two strata"), call. = FALSE)
  }
  if (!setequal(names(standard_weights), unique(arm_strata))) {
    stop(paste(
      "saturated categorical IPW standard_weights must name exactly the",
      "declared strata"), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.dsvert_ipw_stratified_mapping <- function(
    row_levels, arm_strata, arm_treatment, standard_weights, treated) {
  valid_levels <- is.character(row_levels) && length(row_levels) >= 4L &&
    !anyNA(row_levels) && all(nzchar(row_levels)) && !anyDuplicated(row_levels)
  if (!isTRUE(valid_levels) || !identical(names(arm_strata), row_levels) ||
      !identical(names(arm_treatment), row_levels)) {
    stop(paste(
      "signed arm levels must exactly match the names of arm_strata and",
      "arm_treatment"), call. = FALSE)
  }
  strata <- unname(arm_strata[row_levels])
  treatment <- unname(arm_treatment[row_levels])
  if (length(unique(treatment)) != 2L || !treated %in% treatment ||
      !setequal(names(standard_weights), unique(strata))) {
    stop(paste(
      "signed arm mapping must contain the requested binary treatment and",
      "standard_weights for exactly its strata"), call. = FALSE)
  }
  list(strata = strata, treatment = treatment,
       standard_weights = standard_weights)
}

#' @export
print.ds.vertIPW <- function(x, ...) {
  title <- if (identical(x$propensity_model,
                         "saturated_categorical_treatment_given_stratum")) {
    "dsVert saturated categorical-stratum IPW-equivalent ATE from a sticky DP table"
  } else {
    "dsVert intercept-only IPW ATE from a sticky DP table"
  }
  cat(title, "\n")
  cat("  risk difference = ", format(x$estimate, digits = 6L), "\n", sep = "")
  if (is.numeric(x$confidence_region) && length(x$confidence_region) == 2L) {
    cat("  combined region = [", format(x$confidence_region[[1L]], digits = 6L),
        ", ", format(x$confidence_region[[2L]], digits = 6L), "]\n", sep = "")
  }
  cat("  No individual weights, propensity fit, standard errors or p-values are released.\n")
  invisible(x)
}
