
test_that("Cox client accepts only a literal additive Surv formula", {
  value <- .dsvert_dp_cox_grid_cross_formula(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z)
  expect_identical(value$predictors, c("site_a$x", "site_b$z"))
  for (formula in list(y ~ x, Surv(t, e) ~ log(x), Surv(t, e) ~ x * z,
      Surv(t, e) ~ x + strata(z), Surv(t, e) ~ x - z,
      Surv(start, t, e) ~ x, Surv(t, e) ~ x + x)) {
    expect_error(.dsvert_dp_cox_grid_cross_formula(formula),
                 class = "dsvert_dp_public_failure")
  }
})

test_that("Cox cross-owner contracts require both real custodian signatures", {
  fixture <- .cox_cross_client_fixture()
  validate <- function(value) .dsvert_dp_cox_grid_cross_contract_validate(
    value, fixture$policy, fixture$schema_manifest)
  expect_identical(validate(fixture$contract),
                   .dsvert_joint_dp_client_canonical(fixture$contract))
  for (peer in names(fixture$contract$signatures)) {
    changed <- fixture$contract
    changed$signatures[[peer]] <- NULL
    expect_error(validate(changed), class = "dsvert_dp_public_failure")
    changed <- fixture$contract
    changed$signatures[[peer]] <- strrep("A", 86)
    expect_error(validate(changed), class = "dsvert_dp_public_failure")
  }
  for (change in list(
      function(x) { x$spec$ties <- "efron"; x },
      function(x) { x$spec$intercept <- TRUE; x },
      function(x) { x$spec$time_semantics <- "binned"; x },
      function(x) { x$spec$numeric_contract <- list(); x },
      function(x) { x$spec$sensitivity$raw_l1_sensitivity <- 0; x },
      function(x) { x$artifact$implementation_state <- "same_owner_materialized"; x },
      function(x) { x$source_contract$purpose <- "generic_nonlinear_rpc"; x })) {
    expect_error(validate(fixture$sign(change(fixture$contract))),
                 class = "dsvert_dp_public_failure")
  }
})

test_that("Cox client canonical argmin rescales slopes and omits inference", {
  fixture <- .cox_cross_client_fixture()
  spec <- fixture$contract$spec
  coordinates <- rep(0, length(spec$beta_grid))
  moment <- .dsvert_dp_cox_grid_cross_moment(coordinates, spec)
  expect_identical(moment$selected_candidate, 1L)
  expect_equal(unname(moment$coefficients), unlist(spec$beta_grid[[1L]]))
  for (bad in list(-coordinates - 1, coordinates + 0.5,
                   rep(Inf, length(coordinates)), coordinates[-1L],
                   unlist(spec$sensitivity$maximum_coordinates) + 1)) {
    expect_error(.dsvert_dp_cox_grid_cross_moment(bad, spec),
                 class = "dsvert_dp_public_failure")
  }
  release <- function(...) c(fixture[c("contract", "policy", "schema_manifest")],
                             list(coordinates = coordinates))
  result <- .dsvert_dp_cox_grid_cross_impl(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL), release)
  expect_identical(result$selected_candidate, 1L)
  expect_null(result$std_errors)
  expect_null(result$baseline_hazard)
  expect_false(result$production_ready)
  expect_false(result$protected_optimizer_called)
  expect_error(.dsvert_dp_cox_grid_cross_impl(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$wrong,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL), release),
    class = "dsvert_dp_public_failure")
})

test_that("production Cox client is closed and exposes no test evaluator", {
  expect_identical(names(formals(dp_cox_grid)),
    c("formula", "data", "analysis_id", "datasources"))
  expect_error(dp_cox_grid(Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL)),
    class = "dsvert_dp_public_failure")
  expect_false(.dsvert_dp_cox_grid_cross_client_register()$runtime_enabled)
  expect_identical(.dsvert_dp_cox_grid_cross_client_register()$entry, dp_cox_grid)
  expect_false("dp_cox_grid" %in% getNamespaceExports("dsVertClient"))
})

test_that("Cox server and client canonical contracts agree", {
  skip_if_not_installed("dsVert")
  server <- asNamespace("dsVert")
  if (!exists(".dsvert_dp_cox_grid_cross_register", server, inherits=FALSE)) {
    skip("companion server lacks additive Cox registration")
  }
  f <- .cox_cross_client_fixture()
  validated <- get(".dsvert_dp_cox_grid_cross_contract_validate", server)(
    f$contract, f$policy, f$schema_manifest)
  expect_identical(.dsvert_joint_dp_client_json(validated),
    .dsvert_joint_dp_client_json(.dsvert_dp_cox_grid_cross_contract_validate(
      f$contract,f$policy,f$schema_manifest)))
})
