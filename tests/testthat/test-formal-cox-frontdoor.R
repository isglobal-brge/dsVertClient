.formal_cox_public_release <- function(request) {
  list(
    version = "dsvert-formal-cox-public-result-v1",
    analysis_id = request$analysis_id,
    artifact_id = digest::digest("cox-artifact", algo = "sha256", serialize = FALSE),
    certificate_sha256 = digest::digest(
      "cox-certificate", algo = "sha256", serialize = FALSE),
    formula_sha256 = request$formula_sha256,
    coefficients = list(
      list(coefficient = "x", beta_steps = "64", fraction_bits = 8L,
           beta = 0.25, hazard_ratio_lower = 1.28,
           hazard_ratio_upper = 1.29, hazard_ratio_midpoint = 1.285),
      list(coefficient = "z", beta_steps = "-32", fraction_bits = 8L,
           beta = -0.125, hazard_ratio_lower = 0.88,
           hazard_ratio_upper = 0.89, hazard_ratio_midpoint = 0.885)),
    production_ready = FALSE)
}

test_that("formal Cox reads one certified public release at K=2/3/5", {
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      calls <<- calls + 1L
      expect_identical(operation, "formal Cox public-result retrieval")
      expect_identical(as.character(expr[[1L]]),
                       "dsvertFormalCoxPublicResultDS")
      request <- as.list(expr[-1L])
      expected <- .formal_cox_public_release(request)
      stats::setNames(rep(list(expected), length(conns)), names(conns))
    },
    .package = "dsVertClient")
  for (index in seq_along(c(2L, 3L, 5L))) {
    k <- c(2L, 3L, 5L)[[index]]
    peers <- c("site_a", "site_b", "witness_1", "witness_2",
               "witness_3")[seq_len(k)]
    conns <- stats::setNames(lapply(peers, function(...) {
      structure(list(), class = "mock")
    }), peers)
    fit <- ds.vertCox(
      stats::as.formula("Surv(time, status) ~ x + z"), data = "study",
      formal_analysis_id = "primary_cox", verbose = FALSE,
      datasources = conns)
    expect_identical(calls, as.integer(index))
    expect_s3_class(fit, "dsvert_formal_dp_cox")
    expect_s3_class(fit, "ds.vertCox")
    expect_equal(fit$coefficients, c(x = 0.25, z = -0.125))
    expect_equal(fit$hazard_ratio, c(x = 1.285, z = 0.885))
    expect_identical(fit$coefficient_lattice_steps[["z"]], "-32")
    expect_false(fit$production_ready)
    expect_false(fit$source_values_exposed)
    expect_null(fit$covariance)
    expect_null(fit$std_errors)
  }
})

test_that("fresh formal Cox projects only one final committed public release", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_run = function(conns, selector, .aggregate) {
      expect_identical(names(conns), c("site_a", "site_b"))
      list(
        analysis_id = selector$analysis_id, schema_sha256 = strrep("a", 64L),
        total_blocks = 1L, state = "finalizer_committed",
        public_result = list(
          version = "dsvert-formal-cox-public-result-v1",
          artifact_id = strrep("b", 64L), certificate_sha256 = strrep("c", 64L),
          valid = TRUE, coefficients = list(list(
            index = 0L, beta_steps = "64", fraction_bits = 8L, beta = 0.25,
            hazard_ratio_lower = 1.2, hazard_ratio_upper = 1.3,
            hazard_ratio_midpoint = 1.25)), production_ready = FALSE),
        production_ready = FALSE)
    }, .package = "dsVertClient")
  fit <- ds.vertCox(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    fresh_formal_analysis_id = "fresh_cox", verbose = FALSE,
    datasources = conns)
  expect_s3_class(fit, "ds.vertCox")
  expect_identical(fit$called_via, "ds.vertCox_fresh_formal_analysis_id")
  expect_equal(fit$coefficients, c(x = 0.25))
  expect_false(fit$production_ready)
  expect_error(ds.vertCox(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    formal_analysis_id = "stored", fresh_formal_analysis_id = "fresh_cox",
    datasources = conns), "mutually exclusive")
})

test_that("Cox analysis_id is the standard fresh-analysis selector", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_run = function(conns, selector, .aggregate) {
      expect_identical(names(conns), c("site_a", "site_b"))
      expect_identical(selector$analysis_id, "cox_standard")
      list(
        analysis_id = selector$analysis_id, schema_sha256 = strrep("a", 64L),
        total_blocks = 1L, state = "finalizer_committed",
        public_result = list(
          version = "dsvert-formal-cox-public-result-v1",
          artifact_id = strrep("b", 64L), certificate_sha256 = strrep("c", 64L),
          valid = TRUE, coefficients = list(list(
            index = 0L, beta_steps = "64", fraction_bits = 8L, beta = 0.25,
            hazard_ratio_lower = 1.2, hazard_ratio_upper = 1.3,
            hazard_ratio_midpoint = 1.25)), production_ready = FALSE),
        production_ready = FALSE)
    }, .package = "dsVertClient")
  fit <- ds.vertCox(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    analysis_id = "cox_standard", verbose = FALSE, datasources = conns)
  expect_s3_class(fit, "ds.vertCox")
  expect_identical(fit$called_via, "ds.vertCox_analysis_id")
  expect_equal(fit$coefficients, c(x = 0.25))
  cox <- ds.vert.cox(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    analysis_id = "cox_standard", verbose = FALSE, datasources = conns)
  coxph <- ds.vert.coxph(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    analysis_id = "cox_standard", verbose = FALSE, datasources = conns)
  profile <- ds.vertCoxProfileNonDisclosive(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    analysis_id = "cox_standard", verbose = FALSE, datasources = conns)
  expect_identical(cox$frontdoor, "ds.vert.cox")
  expect_identical(coxph$frontdoor, "ds.vert.coxph")
  expect_equal(cox$coefficients, c(x = 0.25))
  expect_equal(coxph$coefficients, c(x = 0.25))
  expect_identical(profile$called_via,
                   "ds.vertCoxProfileNonDisclosive_analysis_id")
  expect_error(ds.vertCox(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    analysis_id = "cox_standard", fresh_formal_analysis_id = "old",
    datasources = conns), "mutually exclusive")
})

test_that("Cox aliases advertise the configured analysis selector", {
  status <- ds.vertMethodStatus(c(
    "ds.vertCox", "ds.vert.cox", "ds.vert.coxph",
    "ds.vertCoxProfileNonDisclosive"))
  expect_true(all(status$status == "promoted"))
  expect_true(all(grepl("analysis_id", status$safe_scope, fixed = TRUE)))
  aliases <- status$method %in% c(
    "ds.vert.coxph", "ds.vertCoxProfileNonDisclosive")
  expect_true(all(grepl("cannot choose source work",
                        status$principal_limitation[aliases], fixed = TRUE)))
})

test_that("fresh formal Cox aliases retain the profile-only committed route", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_run = function(conns, selector, .aggregate) {
      list(
        analysis_id = selector$analysis_id, schema_sha256 = strrep("a", 64L),
        total_blocks = 1L, state = "finalizer_already_public",
        public_result = list(
          version = "dsvert-formal-cox-public-result-v1",
          artifact_id = strrep("b", 64L), certificate_sha256 = strrep("c", 64L),
          valid = TRUE, coefficients = list(list(
            index = 0L, beta_steps = "64", fraction_bits = 8L, beta = 0.25,
            hazard_ratio_lower = 1.2, hazard_ratio_upper = 1.3,
            hazard_ratio_midpoint = 1.25)), production_ready = FALSE),
        production_ready = FALSE)
    }, .package = "dsVertClient")
  formula <- stats::as.formula("Surv(time, status) ~ x")
  cox <- ds.vert.cox(
    formula, data = "study", fresh_formal_analysis_id = "fresh_cox",
    verbose = FALSE, datasources = conns)
  coxph <- ds.vert.coxph(
    formula, data = "study", fresh_formal_analysis_id = "fresh_cox",
    verbose = FALSE, datasources = conns)
  expect_identical(cox$frontdoor, "ds.vert.cox")
  expect_identical(coxph$frontdoor, "ds.vert.coxph")
  expect_equal(cox$coefficients, c(x = 0.25))
  expect_equal(coxph$coefficients, c(x = 0.25))
  expect_error(ds.vert.cox(
    formula, data = "study", method = "discrete",
    fresh_formal_analysis_id = "fresh_cox", datasources = conns),
    "does not accept method='discrete'")
})

test_that("formal Cox rejects a cross-site certificate mismatch", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      request <- as.list(expr[-1L])
      left <- .formal_cox_public_release(request)
      right <- left
      right$certificate_sha256 <- digest::digest(
        "other Cox certificate", algo = "sha256", serialize = FALSE)
      stats::setNames(list(left, right), names(conns))
    },
    .package = "dsVertClient")
  expect_error(ds.vertCox(
    stats::as.formula("Surv(time, status) ~ x + z"), data = "study",
    formal_analysis_id = "primary_cox", verbose = FALSE,
    datasources = conns), "different formal Cox public releases")
})

test_that("formal Cox rejects a decimal that disagrees with signed lattice steps", {
  request <- .dsvert_formal_cox_frontdoor_request(
    "primary_cox", "study", stats::as.formula("Surv(time, status) ~ x + z"))$value
  release <- .formal_cox_public_release(request)
  release$coefficients[[1L]]$beta <- 0.5
  expect_error(
    .dsvert_formal_cox_frontdoor_public_response(release, request),
    "does not match its signed lattice steps")
})

test_that("formal Cox permits only registered public selectors", {
  first <- .dsvert_formal_cox_frontdoor_request(
    "primary_cox", "study", stats::as.formula("Surv(time, status) ~ x + z"))
  retry <- .dsvert_formal_cox_frontdoor_request(
    "primary_cox", "study", stats::as.formula("Surv(time, status) ~ x + z"))
  expect_identical(first, retry)
  expect_setequal(names(first$value), c(
    "version", "analysis_id", "data_name", "formula_sha256",
    "public_selectors", "privacy_controls", "role_selection"))
  expect_identical(unlist(first$value$public_selectors, use.names = FALSE),
                   c("analysis_id", "data_name", "formula_sha256"))
  expect_error(ds.vertCox(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    formal_analysis_id = "primary_cox", max_iter = 10L), "max_iter")
})

test_that("formal Cox aliases preserve the certified release route", {
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      request <- as.list(expr[-1L])
      stats::setNames(list(.formal_cox_public_release(request)), names(conns))
    },
    .package = "dsVertClient")
  fit <- ds.vert.cox(
    stats::as.formula("Surv(time, status) ~ x + z"), data = "study",
    formal_analysis_id = "primary_cox", verbose = FALSE,
    datasources = list(site_a = structure(list(), class = "mock")))
  expect_s3_class(fit, "dsvert_formal_dp_cox")
  profile <- ds.vertCoxProfileNonDisclosive(
    stats::as.formula("Surv(time, status) ~ x + z"), data = "study",
    formal_analysis_id = "primary_cox", verbose = FALSE,
    datasources = list(site_a = structure(list(), class = "mock")))
  expect_s3_class(profile, "dsvert_formal_dp_cox")
  expect_error(ds.vertCoxProfileNonDisclosive(
    stats::as.formula("Surv(time, status) ~ x"), data = "study",
    formal_analysis_id = "primary_cox", max_iter = 2L), "max_iter")
})
