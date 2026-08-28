test_that("intercept-only binary IPW is post-processing of one sticky table", {
  observed <- list()
  table <- structure(list(
    marker = "sticky-table", row_levels = c("control", "yes"),
    artifact_key = "artifact", final_vector_root = "root"),
                     class = "ds.vertDPContingency")
  mechanism <- list(
    risk_treated = c(lower = 0.50, upper = 0.70),
    risk_control = c(lower = 0.20, upper = 0.40),
    risk_difference = c(lower = 0.10, upper = 0.50))
  causal <- list(
    status = "ok",
    point_estimates = list(
      risk_treated = 0.60, risk_control = 0.30,
      risk_difference = 0.30),
    mechanism_regions = mechanism,
    uncertainty_scope = "DP mechanism noise only")
  combined <- list(
    status = "ok",
    combined_regions = mechanism,
    uncertainty_scope = "joint DP and sampling uncertainty")
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(data_name, row_var, col_var,
                                    server = NULL, datasources = NULL) {
      observed$table <<- list(
        data_name = data_name, row_var = row_var, col_var = col_var,
        server = server, datasources = datasources)
      table
    },
    ds.vertDPCausalStandardization = function(
        x, strata, treatment, treated, standard_weights, event, level) {
      observed$causal <<- list(
        x = x, strata = strata, treatment = treatment, treated = treated,
        standard_weights = standard_weights, event = event, level = level)
      causal
    },
    ds.vertDPCausalStandardizationInference = function(
        x, strata, treatment, treated, standard_weights, event, level) {
      observed$inference <<- list(
        x = x, strata = strata, treatment = treatment, treated = treated,
        standard_weights = standard_weights, event = event, level = level)
      combined
    },
    .package = "dsVertClient")

  fit <- ds.vertIPW(
    outcome_formula = outcome ~ treatment,
    propensity_formula = treatment ~ 1,
    data = "cohort", outcome_family = "binomial", treated = "yes",
    event = "event", level = 0.90, server = "site_a",
    datasources = list(site_a = structure(list(), class = "mock")))

  expect_s3_class(fit, "ds.vertIPW")
  expect_identical(fit$estimate, 0.30)
  expect_identical(fit$estimand, "ATE")
  expect_identical(fit$propensity_model, "intercept_only")
  expect_identical(fit$weights_released, FALSE)
  expect_identical(fit$server_calls_for_artifact, 1L)
  expect_identical(fit$additional_server_calls_after_artifact, 0L)
  expect_identical(fit$additional_privacy_cost_after_artifact,
                   c(epsilon = 0, delta = 0))
  expect_identical(observed$table[c("data_name", "row_var", "col_var", "server")],
                   list(data_name = "cohort", row_var = "treatment",
                        col_var = "outcome", server = "site_a"))
  expect_identical(observed$causal$x, table)
  expect_identical(observed$causal$strata, c("overall", "overall"))
  expect_identical(observed$causal$treatment, c("control", "yes"))
  expect_identical(observed$causal$treated, "yes")
  expect_identical(observed$causal$standard_weights, c(overall = 1))
  expect_identical(observed$causal$event, "event")
  expect_identical(observed$causal$level, 0.90)
  expect_identical(observed$inference, observed$causal)
})

test_that("intercept-only IPW rejects unsupported models before DSI", {
  reached_table <- FALSE
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(...) {
      reached_table <<- TRUE
      stop("test observed a DP table call", call. = FALSE)
    },
    .package = "dsVertClient")
  base <- list(
    outcome_formula = outcome ~ treatment,
    propensity_formula = treatment ~ 1,
    data = "cohort", outcome_family = "binomial", treated = "yes")
  expect_error(do.call(ds.vertIPW, utils::modifyList(base, list(
    outcome_formula = outcome ~ treatment + age))), "outcome_formula")
  expect_error(do.call(ds.vertIPW, utils::modifyList(base, list(
    propensity_formula = treatment ~ age))), "arm_column")
  expect_error(do.call(ds.vertIPW, utils::modifyList(base, list(outcome_family = "gaussian"))),
               "binomial")
  expect_error(do.call(ds.vertIPW, utils::modifyList(base, list(weights_column = "custom"))),
               "weights_column")
  without_treated <- base
  without_treated$treated <- NULL
  expect_error(do.call(ds.vertIPW, without_treated), "treated")
  expect_false(reached_table)
})

test_that("saturated categorical-stratum IPW consumes one matched signed arm table", {
  observed <- list()
  arm_levels <- c("s1_control", "s1_treated", "s2_control", "s2_treated")
  arm_strata <- stats::setNames(c("s1", "s1", "s2", "s2"), arm_levels)
  arm_treatment <- stats::setNames(
    c("control", "treated", "control", "treated"), arm_levels)
  weights <- c(s1 = 0.25, s2 = 0.75)
  table <- structure(list(
    marker = "sticky-stratified-table", row_levels = arm_levels,
    artifact_key = "artifact", final_vector_root = "root"),
    class = "ds.vertDPContingency")
  mechanism <- list(
    risk_treated = c(lower = 0.50, upper = 0.70),
    risk_control = c(lower = 0.20, upper = 0.40),
    risk_difference = c(lower = 0.10, upper = 0.50))
  causal <- list(
    status = "ok",
    point_estimates = list(
      risk_treated = 0.60, risk_control = 0.30,
      risk_difference = 0.30),
    mechanism_regions = mechanism,
    uncertainty_scope = "DP mechanism noise only")
  combined <- list(
    status = "ok", combined_regions = mechanism,
    uncertainty_scope = "joint DP and sampling uncertainty")
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(data_name, row_var, col_var,
                                    server = NULL, datasources = NULL) {
      observed$table <<- list(
        data_name = data_name, row_var = row_var, col_var = col_var,
        server = server, datasources = datasources)
      table
    },
    ds.vertDPCausalStandardization = function(
        x, strata, treatment, treated, standard_weights, event, level) {
      observed$causal <<- list(
        x = x, strata = strata, treatment = treatment, treated = treated,
        standard_weights = standard_weights, event = event, level = level)
      causal
    },
    ds.vertDPCausalStandardizationInference = function(
        x, strata, treatment, treated, standard_weights, event, level) {
      observed$inference <<- list(
        x = x, strata = strata, treatment = treatment, treated = treated,
        standard_weights = standard_weights, event = event, level = level)
      combined
    },
    .package = "dsVertClient")

  fit <- ds.vertIPW(
    outcome ~ treatment, treatment ~ stratum, data = "cohort",
    outcome_family = "binomial", treated = "treated", event = "event",
    arm_column = "arm", arm_strata = arm_strata,
    arm_treatment = arm_treatment, standard_weights = weights,
    server = "site_a", datasources = list(site_a = structure(list(), class = "mock")),
    verbose = FALSE)

  expect_s3_class(fit, "ds.vertIPW")
  expect_identical(fit$estimate, 0.30)
  expect_identical(fit$propensity_model,
                   "saturated_categorical_treatment_given_stratum")
  expect_identical(fit$strata_variable, "stratum")
  expect_identical(fit$standard_weights, weights)
  expect_identical(fit$result_contract,
    "saturated_categorical_ipw_g_formula_equivalence_from_sticky_dp_table_v1")
  expect_identical(observed$table[c("data_name", "row_var", "col_var", "server")],
                   list(data_name = "cohort", row_var = "arm",
                        col_var = "outcome", server = "site_a"))
  expect_identical(observed$causal$strata, unname(arm_strata))
  expect_identical(observed$causal$treatment, unname(arm_treatment))
  expect_identical(observed$causal$standard_weights, weights)
  expect_identical(observed$inference, observed$causal)
  expect_output(print(fit), "IPW-equivalent ATE")
})

test_that("saturated categorical-stratum IPW rejects malformed bindings before DSI", {
  reached_table <- FALSE
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(...) {
      reached_table <<- TRUE
      stop("test observed a DP table call", call. = FALSE)
    },
    .package = "dsVertClient")
  arm_levels <- c("s1_control", "s1_treated", "s2_control", "s2_treated")
  base <- list(
    outcome_formula = outcome ~ treatment,
    propensity_formula = treatment ~ stratum,
    data = "cohort", outcome_family = "binomial", treated = "treated",
    arm_column = "arm",
    arm_strata = stats::setNames(c("s1", "s1", "s2", "s2"), arm_levels),
    arm_treatment = stats::setNames(
      c("control", "treated", "control", "treated"), arm_levels),
    standard_weights = c(s1 = 0.5, s2 = 0.5))
  bad_treatment <- base
  names(bad_treatment$arm_treatment)[[4L]] <- "unexpected"
  expect_error(do.call(ds.vertIPW, bad_treatment), "complete binary treatment")
  bad_weights <- base
  bad_weights$standard_weights <- c(s1 = 1, other = 1)
  expect_error(do.call(ds.vertIPW, bad_weights), "standard_weights")
  bad_arm <- base
  bad_arm$arm_column <- "stratum"
  expect_error(do.call(ds.vertIPW, bad_arm), "arm_column must differ")
  expect_false(reached_table)
})

test_that("saturated categorical-stratum IPW rejects mismatched signed arm levels", {
  arm_levels <- c("s1_control", "s1_treated", "s2_control", "s2_treated")
  released <- structure(list(
    row_levels = rev(arm_levels), artifact_key = "artifact",
    final_vector_root = "root"), class = "ds.vertDPContingency")
  causal_called <- FALSE
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(...) released,
    ds.vertDPCausalStandardization = function(...) {
      causal_called <<- TRUE
      stop("test should not reach causal post-processing", call. = FALSE)
    },
    .package = "dsVertClient")
  expect_error(ds.vertIPW(
    outcome ~ treatment, treatment ~ stratum, data = "cohort",
    treated = "treated", arm_column = "arm",
    arm_strata = stats::setNames(c("s1", "s1", "s2", "s2"), arm_levels),
    arm_treatment = stats::setNames(
      c("control", "treated", "control", "treated"), arm_levels),
    standard_weights = c(s1 = 0.5, s2 = 0.5), verbose = FALSE),
    "signed arm levels")
  expect_false(causal_called)
})

test_that("the IPW alias forwards both admitted routes", {
  seen <- NULL
  testthat::local_mocked_bindings(
    ds.vertIPW = function(...) {
      seen <<- list(...)
      structure(list(ok = TRUE), class = "ipw_sentinel")
    },
    .dsvert_datasources = function(datasources) datasources,
    .dsvert_set_frontdoor = function(out, ...) out,
    .dsvert_add_policy = function(out, ...) out,
    .package = "dsVertClient")
  result <- ds.vert.ipw(
    outcome ~ treatment, treatment ~ 1, data = "cohort", treated = "yes",
    datasources = list(site_a = structure(list(), class = "mock")))
  expect_s3_class(result, "ipw_sentinel")
  expect_identical(seen$propensity_formula, treatment ~ 1)

  arm_levels <- c("s1_control", "s1_treated", "s2_control", "s2_treated")
  result <- ds.vert.ipw(
    outcome ~ treatment, treatment ~ stratum, data = "cohort",
    arm_column = "arm",
    arm_strata = stats::setNames(c("s1", "s1", "s2", "s2"), arm_levels),
    arm_treatment = stats::setNames(
      c("control", "treated", "control", "treated"), arm_levels),
    standard_weights = c(s1 = 0.5, s2 = 0.5),
    datasources = list(site_a = structure(list(), class = "mock")))
  expect_s3_class(result, "ipw_sentinel")
  expect_identical(seen$propensity_formula, treatment ~ stratum)
  expect_identical(seen$arm_column, "arm")
})
