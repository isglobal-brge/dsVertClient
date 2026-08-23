.formal_cox_discrete_public_release <- function(request) {
  list(
    version = "dsvert-formal-cox-discrete-public-result-v1",
    analysis_id = request$analysis_id,
    artifact_id = digest::digest("discrete-artifact", algo = "sha256",
                                 serialize = FALSE),
    certificate_sha256 = digest::digest("discrete-certificate", algo = "sha256",
                                        serialize = FALSE),
    target = "discrete_logit",
    source_formula_sha256 = request$source_formula_sha256,
    model_formula_sha256 = digest::digest("event ~ time_bin + x", algo = "sha256",
                                          serialize = FALSE),
    time_grid_sha256 = digest::digest("fixed-grid", algo = "sha256",
                                      serialize = FALSE),
    coefficients = list(
      list(coefficient = "(Intercept)", signed_steps = "0",
           output_lattice_bits = 20L, value = 0),
      list(coefficient = "time_bin2", signed_steps = "262144",
           output_lattice_bits = 20L, value = 0.25),
      list(coefficient = "x", signed_steps = "-131072",
           output_lattice_bits = 20L, value = -0.125)),
    production_ready = FALSE)
}

test_that("formal discrete-time Cox reads one bound public release at every site", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      calls <<- calls + 1L
      expect_identical(operation, "formal discrete-time public-result retrieval")
      expect_identical(as.character(expr[[1L]]),
                       "dsvertFormalCoxDiscretePublicResultDS")
      request <- as.list(expr[-1L])
      expected <- .formal_cox_discrete_public_release(request)
      stats::setNames(list(expected, expected), names(conns))
    },
    .package = "dsVertClient")
  fit <- ds.vertCoxDiscreteNonDisclosive(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    formal_analysis_id = "primary_discrete", verbose = FALSE,
    datasources = conns)
  expect_identical(calls, 1L)
  expect_s3_class(fit, "dsvert_formal_dp_cox_discrete")
  expect_s3_class(fit, "ds.vertCoxDiscreteNonDisclosive")
  expect_equal(fit$coefficients, c(`(Intercept)` = 0, time_bin2 = 0.25,
                                    x = -0.125))
  expect_equal(fit$hazard_odds_ratio[["x"]], exp(-0.125))
  expect_identical(fit$coefficient_lattice_steps[["x"]], "-131072")
  expect_false(fit$production_ready)
  expect_false(fit$source_values_exposed)
  expect_null(fit$covariance)
  expect_null(fit$std_errors)
})

test_that("formal discrete-time Cox rejects cross-site and legacy-selector mismatches", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      request <- as.list(expr[-1L])
      left <- .formal_cox_discrete_public_release(request)
      right <- left
      right$time_grid_sha256 <- digest::digest("other-grid", algo = "sha256",
                                               serialize = FALSE)
      stats::setNames(list(left, right), names(conns))
    },
    .package = "dsVertClient")
  expect_error(ds.vertCoxDiscreteNonDisclosive(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    formal_analysis_id = "primary_discrete", verbose = FALSE,
    datasources = conns), "different formal discrete-time public releases")

  expect_error(ds.vertCoxDiscreteNonDisclosive(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    formal_analysis_id = "primary_discrete", J = 3L), "J")
  expect_error(ds.vertCoxDiscreteNonDisclosive(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    formal_analysis_id = "primary_discrete", target = "cox_profile"),
    "target")
})

test_that("formal discrete-time Cox rejects a decimal that disagrees with signed lattice steps", {
  request <- .dsvert_formal_cox_discrete_frontdoor_request(
    "primary_discrete", "study", stats::as.formula("Surv(time, status) ~ x"))$value
  release <- .formal_cox_discrete_public_release(request)
  release$coefficients[[2L]]$value <- 0.5
  expect_error(
    .dsvert_formal_cox_discrete_frontdoor_response(release, request),
    "does not match its signed lattice steps")
})

test_that("formal discrete-time Cox requests are deterministic and minimal", {
  first <- .dsvert_formal_cox_discrete_frontdoor_request(
    "primary_discrete", "study", stats::as.formula("Surv(time, status) ~ x"))
  retry <- .dsvert_formal_cox_discrete_frontdoor_request(
    "primary_discrete", "study", stats::as.formula("Surv(time, status) ~ x"))
  expect_identical(first, retry)
  expect_setequal(names(first$value), c(
    "version", "analysis_id", "data_name", "source_formula_sha256",
    "public_selectors", "privacy_controls", "role_selection"))
  expect_identical(unlist(first$value$public_selectors, use.names = FALSE),
                   c("analysis_id", "data_name", "source_formula_sha256"))
})

test_that("the discrete ds.vert.cox selector preserves its distinct estimand", {
  conns <- list(site_a = structure(list(), class = "mock"))
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      request <- as.list(expr[-1L])
      stats::setNames(list(.formal_cox_discrete_public_release(request)),
                      names(conns))
    },
    .package = "dsVertClient")
  fit <- ds.vert.cox(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    method = "discrete", formal_analysis_id = "primary_discrete",
    verbose = FALSE, datasources = conns)
  expect_s3_class(fit, "dsvert_formal_dp_cox_discrete")
  expect_identical(fit$frontdoor, "ds.vert.cox")
  expect_identical(fit$called_via,
                   "ds.vertCoxDiscreteNonDisclosive_formal_analysis_id")
})
