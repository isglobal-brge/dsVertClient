#' @title Sticky saturated inverse-probability weighting
#' @description Computes a binary ATE from one signed sticky DP contingency
#' release. With an intercept-only propensity model, IPW is exactly the
#' treated-minus-control risk difference. With one categorical propensity
#' stratum, a signed stratum-treatment-by-outcome table and fixed public ATE
#' target weights, the saturated IPW identity is exactly the corresponding
#' stratum-standardised g-formula. ATT and ATC use the same signed table and
#' derive their target-arm stratum weights deterministically from that sticky
#' release. Neither route releases weights, source rows, a propensity fit, or
#' performs a second private computation.
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
#' @param standard_weights Named, non-negative public ATE target-population
#'   weights for every declared stratum. Required for the stratified
#'   \code{estimand = "ATE"} route and forbidden for ATT/ATC.
#' @param level DP-mechanism and combined-region coverage level.
#' @param server Optional same-owner assertion for the signed table artifact.
#' @param verbose Print a short route message.
#' @param datasources DataSHIELD connections.
#' @param estimand One of \code{"ATE"}, \code{"ATT"} or \code{"ATC"}. ATT and
#'   ATC require \code{treatment ~ stratum}; their target weights are derived
#'   from the signed released target arm, so they provide a DP-mechanism region
#'   only, not sampling inference.
#' @param ... No additional controls are accepted.
#' @return A \code{ds.vertIPW} object with a conservative DP-mechanism and
#'   combined DP-plus-sampling region for ATE, or a DP-mechanism region for
#'   ATT/ATC, and no released individual weights or propensity fit.
#' @details The categorical stratum route is an exact saturated IPW/g-formula
#' identity, not a fitted propensity-score model. Its release is one canonical
#' stratum-treatment-by-outcome DP table; unlimited retries return the same
#' sticky release and do not create another draw. Arbitrary covariates,
#' continuous strata, propensity fitting, trimming and outcome regression
#' remain unavailable. ATT/ATC target-arm weights are release-adaptive, so the
#' reported outer region accounts for the simultaneous DP count box but does
#' not claim superpopulation sampling coverage.
#' @seealso \code{\link{ds.vertDPCausalStandardization}}
#' @export
ds.vertIPW <- function(outcome_formula, propensity_formula, data = NULL,
                      weights_column = "ipw", outcome_family = "binomial",
                      treated = NULL, event = 2L, level = 0.95,
                      arm_column = NULL, arm_strata = NULL,
                      arm_treatment = NULL, standard_weights = NULL,
                      server = NULL, verbose = TRUE, datasources = NULL,
                      estimand = c("ATE", "ATT", "ATC"), ...) {
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
  estimand <- match.arg(estimand)
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
    if (!identical(estimand, "ATE")) {
      stop("ATT and ATC require treatment ~ one categorical stratum",
           call. = FALSE)
    }
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
    if (!identical(estimand, "ATE") && !is.null(standard_weights)) {
      stop("ATT and ATC derive target weights from the signed table and do not accept standard_weights",
           call. = FALSE)
    }
    .dsvert_ipw_validate_stratified_controls(
      arm_column, arm_strata, arm_treatment, standard_weights,
      require_standard_weights = identical(estimand, "ATE"))
    if (arm_column %in% c(spec$outcome, spec$treatment, spec$stratum)) {
      stop(paste(
        "saturated categorical IPW arm_column must differ from outcome,",
        "treatment and stratum"), call. = FALSE)
    }
    table <- ds.vertDPContingency(
      data, arm_column, spec$outcome, server = server,
      datasources = datasources)
    mapping_weights <- if (identical(estimand, "ATE")) {
      standard_weights
    } else {
      stats::setNames(rep(1, length(unique(arm_strata))), unique(arm_strata))
    }
    mapping <- .dsvert_ipw_stratified_mapping(
      table$row_levels, arm_strata, arm_treatment, mapping_weights, treated)
    causal_arguments <- list(
      x = table, strata = mapping$strata, treatment = mapping$treatment,
      treated = treated, standard_weights = mapping$standard_weights,
      event = event, level = level)
    propensity_model <- "saturated_categorical_treatment_given_stratum"
    outcome_model <- "saturated_stratum_treatment_binary_means"
    result_contract <- if (identical(estimand, "ATE")) {
      "saturated_categorical_ipw_g_formula_equivalence_from_sticky_dp_table_v1"
    } else {
      paste0("saturated_categorical_ipw_", tolower(estimand),
             "_from_sticky_dp_table_v1")
    }
    strata_variable <- spec$stratum
  }
  if (identical(estimand, "ATE")) {
    causal <- do.call(ds.vertDPCausalStandardization, causal_arguments)
    combined <- do.call(ds.vertDPCausalStandardizationInference, causal_arguments)
  } else {
    causal <- .dsvert_ipw_target_arm_postprocess(
      x = table, strata = causal_arguments$strata,
      treatment = causal_arguments$treatment, treated = treated,
      estimand = estimand, event = event, level = level)
    causal_arguments$standard_weights <- causal$standard_weights
    combined <- NULL
  }
  if (isTRUE(verbose)) {
    message(if (identical(spec$route, "intercept_only")) {
      "[ds.vertIPW] intercept-only ATE from one signed sticky DP table"
    } else if (identical(estimand, "ATE")) {
      "[ds.vertIPW] saturated categorical stratum ATE from one signed sticky DP table"
    } else {
      paste0("[ds.vertIPW] saturated categorical stratum ", estimand,
             " from one signed sticky DP table")
    })
  }
  out <- list(
    status = causal$status,
    estimand = estimand,
    estimate = causal$point_estimates$risk_difference,
    risk_treated = causal$point_estimates$risk_treated,
    risk_control = causal$point_estimates$risk_control,
    mechanism_region = causal$mechanism_regions$risk_difference,
    confidence_region = if (is.null(combined)) NULL else
      combined$combined_regions$risk_difference,
    coverage_level = level,
    uncertainty_scope = if (is.null(combined)) causal$uncertainty_scope else
      combined$uncertainty_scope,
    propensity_model = propensity_model,
    outcome_model = outcome_model,
    strata_variable = strata_variable,
    standard_weights = causal_arguments$standard_weights,
    target_weight_source = if (identical(estimand, "ATE")) {
      "fixed_public_standard_weights"
    } else {
      "target_arm_weights_derived_from_signed_sticky_dp_table"
    },
    sampling_inference_available = !is.null(combined),
    weights_released = FALSE,
    result_contract = result_contract,
    source_artifact_key = table$artifact_key,
    final_vector_root = table$final_vector_root,
    sticky_replay = TRUE,
    server_calls_for_artifact = 1L,
    additional_server_calls_after_artifact = 0L,
    additional_privacy_cost_after_artifact = c(epsilon = 0, delta = 0),
    limitations = if (identical(estimand, "ATE")) paste(
      "Only the intercept-only or one categorical saturated-stratum ATE",
      "identity is available. No individual propensity weights, continuous",
      "or multiple covariates, outcome regression, standard errors or p-values",
      "are available.") else paste(
      "ATT/ATC is limited to one categorical saturated stratum. Target weights",
      "are derived from the signed sticky DP target arm and its region is",
      "DP-mechanism-only: sampling inference, individual propensity weights,",
      "continuous or multiple covariates, outcome regression, standard errors",
      "and p-values are unavailable."))
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
    arm_column, arm_strata, arm_treatment, standard_weights,
    require_standard_weights = TRUE) {
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
  if (!isTRUE(require_standard_weights)) return(invisible(TRUE))
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

#' @keywords internal
.dsvert_ipw_weighted_rate_box <- function(
    rate_lower, rate_upper, mass_lower, mass_upper) {
  valid <- function(value, upper = Inf) {
    is.numeric(value) && length(value) > 0L && !anyNA(value) &&
      all(is.finite(value)) && all(value >= 0) && all(value <= upper)
  }
  if (!valid(rate_lower, 1) || !valid(rate_upper, 1) ||
      !valid(mass_lower) || !valid(mass_upper) ||
      length(rate_lower) != length(rate_upper) ||
      length(rate_lower) != length(mass_lower) ||
      any(rate_lower > rate_upper) || any(mass_lower > mass_upper)) {
    stop("Invalid release-bound target-arm rate box", call. = FALSE)
  }
  extreme <- function(rates, upper_to_low) {
    order_index <- order(rates, method = "radix")
    values <- vapply(0:length(rates), function(count) {
      mass <- mass_lower
      if (count > 0L) {
        selected <- if (isTRUE(upper_to_low)) {
          order_index[seq_len(count)]
        } else {
          tail(order_index, count)
        }
        mass[selected] <- mass_upper[selected]
      }
      total <- sum(mass)
      if (total > 0) sum(mass * rates) / total else NA_real_
    }, numeric(1L))
    values[is.finite(values)]
  }
  lower <- extreme(rate_lower, upper_to_low = TRUE)
  upper <- extreme(rate_upper, upper_to_low = FALSE)
  if (!length(lower) || !length(upper)) return(c(lower = 0, upper = 1))
  c(lower = min(lower), upper = max(upper))
}

#' @keywords internal
.dsvert_ipw_target_arm_postprocess <- function(
    x, strata, treatment, treated, estimand, event, level) {
  x <- .dsvert_dp_table_contract(x)
  estimand <- match.arg(estimand, c("ATT", "ATC"))
  design <- .dsvert_dp_causal_design(
    x, strata, treatment, treated,
    stats::setNames(rep(1, length(unique(strata))), unique(strata)), event)
  target_column <- if (identical(estimand, "ATT")) "treated" else "control"
  target_rows <- design$row_index[, target_column]
  other_column <- setdiff(c("treated", "control"), target_column)
  other_rows <- design$row_index[, other_column]
  target_counts <- rowSums(x$table)[target_rows]
  target_total <- sum(target_counts)
  radius <- .dsvert_dp_table_simultaneous_radius(x, level)
  if (target_total <= 0 || !is.finite(target_total)) {
    broad <- .dsvert_dp_causal_effect_regions(
      c(lower = 0, upper = 1), c(lower = 0, upper = 1))
    return(structure(list(
      status = "dp_target_population_non_estimable",
      point_estimates = list(risk_treated = NULL, risk_control = NULL,
                             risk_difference = NULL, risk_ratio = NULL,
                             odds_ratio = NULL),
      point_status = stats::setNames(rep("non_estimable_zero_target_arm", 5L),
        c("risk_treated", "risk_control", "risk_difference", "risk_ratio",
          "odds_ratio")),
      mechanism_regions = broad,
      standard_weights = NULL,
      target_weight_source = "target_arm_weights_derived_from_signed_sticky_dp_table",
      target_population = target_column,
      level = level, simultaneous_radius = radius,
      coverage_method = .dsvert_dp_table_coverage_method(x),
      uncertainty_scope = paste(
        "DP mechanism noise only; target arm is non-estimable and sampling",
        "uncertainty is unavailable"),
      additional_privacy_cost = c(epsilon = 0, delta = 0),
      additional_server_calls = 0L, epsilon = x$epsilon, delta = x$delta,
      server = x$server), class = c("ds.vertIPWTargetArm", "list")))
  }
  box <- .dsvert_dp_count_box(x$table, radius)
  target_rate <- function(table) {
    event_total <- sum(table[target_rows, design$event])
    total <- sum(table[target_rows, , drop = FALSE])
    if (total > 0) event_total / total else NA_real_
  }
  row_risk <- x$table[, design$event] / rowSums(x$table)
  other_risk <- row_risk[other_rows]
  weighted_other <- sum(target_counts * other_risk) / target_total
  target_observed <- target_rate(x$table)
  target_region <- .dsvert_dp_risk_bounds(
    sum(box$lower[target_rows, design$event]),
    sum(box$upper[target_rows, design$event]),
    sum(box$lower[target_rows, design$nonevent]),
    sum(box$upper[target_rows, design$nonevent]))
  other_bounds <- t(vapply(other_rows, function(row) {
    .dsvert_dp_risk_bounds(
      box$lower[row, design$event], box$upper[row, design$event],
      box$lower[row, design$nonevent], box$upper[row, design$nonevent])
  }, numeric(2L)))
  other_region <- .dsvert_ipw_weighted_rate_box(
    other_bounds[, "lower"], other_bounds[, "upper"],
    rowSums(box$lower[target_rows, , drop = FALSE]),
    rowSums(box$upper[target_rows, , drop = FALSE]))
  if (identical(target_column, "treated")) {
    point <- .dsvert_dp_causal_point(target_observed, weighted_other)
    regions <- .dsvert_dp_causal_effect_regions(target_region, other_region)
  } else {
    point <- .dsvert_dp_causal_point(weighted_other, target_observed)
    regions <- .dsvert_dp_causal_effect_regions(other_region, target_region)
  }
  structure(list(
    status = if (all(point$status == "ok")) "ok" else
      if (any(point$status == "boundary_infinite")) "boundary" else
        "dp_point_non_estimable",
    point_estimates = point$values, point_status = point$status,
    mechanism_regions = regions,
    standard_weights = stats::setNames(target_counts / target_total, design$strata),
    target_weight_source = "target_arm_weights_derived_from_signed_sticky_dp_table",
    target_population = target_column, level = level,
    simultaneous_radius = radius,
    coverage_method = paste(.dsvert_dp_table_coverage_method(x),
      "The target-arm mass and opposite-arm risks are jointly optimized over",
      "that same simultaneous count box."),
    uncertainty_scope = paste(
      "DP mechanism noise only; ATT/ATC target weights are derived from the",
      "signed sticky DP table, so sampling uncertainty is unavailable"),
    identification_assumptions = c(
      "consistency", "conditional_exchangeability_within_public_strata",
      "positivity", "no_interference", "correct_public_row_mapping",
      "target_arm_distribution_from_signed_sticky_dp_table"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L, epsilon = x$epsilon, delta = x$delta,
    server = x$server), class = c("ds.vertIPWTargetArm", "list"))
}

#' @export
print.ds.vertIPW <- function(x, ...) {
  title <- if (identical(x$propensity_model,
                         "saturated_categorical_treatment_given_stratum")) {
    paste0("dsVert saturated categorical-stratum IPW-equivalent ",
           x$estimand, " from a sticky DP table")
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
  if (!isTRUE(x$sampling_inference_available)) {
    cat("  ATT/ATC regions cover DP mechanism noise only; sampling inference is unavailable.\n")
  }
  invisible(x)
}
