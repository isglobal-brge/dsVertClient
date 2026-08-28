.formal_glm_public_release <- function(request) {
  list(
    version = "dsvert-formal-glm-public-result-v1",
    analysis_id = request$analysis_id,
    artifact_id = digest::digest("artifact", algo = "sha256", serialize = FALSE),
    certificate_sha256 = digest::digest(
      "certificate", algo = "sha256", serialize = FALSE),
    family = request$family, formula_sha256 = request$formula_sha256,
    coefficients = list(
      list(coefficient = "(Intercept)", signed_steps = "524288",
           output_lattice_bits = 20, value = 0.5),
      list(coefficient = "x", signed_steps = "-262144",
           output_lattice_bits = 20, value = -0.25)),
    production_ready = FALSE)
}

test_that("formal GLM reads the same certified public release at every site", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      calls <<- calls + 1L
      expect_identical(operation, "formal GLM public-result retrieval")
      expect_identical(as.character(expr[[1L]]),
                       "dsvertFormalGLMPublicResultDS")
      request <- as.list(expr[-1L])
      expected <- .formal_glm_public_release(request)
      stats::setNames(list(expected, expected), names(conns))
    },
    .package = "dsVertClient")
  fit <- ds.vertGLM(
    y ~ x, data = "study", family = "binomial",
    formal_analysis_id = "primary_logit", verbose = FALSE,
    datasources = conns)
  expect_identical(calls, 1L)
  expect_s3_class(fit, "dsvert_formal_dp_glm")
  expect_s3_class(fit, "ds.glm")
  expect_equal(fit$coefficients, c("(Intercept)" = 0.5, x = -0.25))
  expect_identical(fit$coefficient_lattice_steps[["x"]], "-262144")
  expect_false(fit$production_ready)
  expect_false(fit$source_values_exposed)
  expect_error(ds.vertLR(fit, fit), class = "dsvert_inference_unavailable")
})

test_that("formal GLM rejects a mismatched public certificate before a fit", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      request <- as.list(expr[-1L])
      left <- .formal_glm_public_release(request)
      right <- left
      right$certificate_sha256 <- digest::digest(
        "other certificate", algo = "sha256", serialize = FALSE)
      stats::setNames(list(left, right), names(conns))
    },
    .package = "dsVertClient")
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "binomial",
    formal_analysis_id = "primary_logit", verbose = FALSE,
    datasources = conns),
    "different formal GLM public releases")
})

test_that("formal GLM rejects a decimal that disagrees with signed lattice steps", {
  request <- .dsvert_formal_glm_frontdoor_request(
    "primary_logit", "study", "binomial", y ~ x)$value
  release <- .formal_glm_public_release(request)
  release$coefficients[[1L]]$value <- 0.75
  expect_error(
    .dsvert_formal_glm_frontdoor_public_response(release, request),
    "does not match its signed lattice steps")
})

test_that("formal lattice projection stays within the signed Ring128 domain", {
  expect_equal(.dsvert_formal_lattice_value(
    "170141183460469231731687303715884105727", 127L, 1,
    "formal GLM coefficient"), 1)
  expect_equal(.dsvert_formal_lattice_value(
    "-170141183460469231731687303715884105728", 127L, -1,
    "formal GLM coefficient"), -1)
  expect_error(.dsvert_formal_lattice_value(
    "170141183460469231731687303715884105728", 127L, 1,
    "formal GLM coefficient"), "lattice coordinate")
})

test_that("formal GLM request carries selectors and no privacy controls", {
  first <- .dsvert_formal_glm_frontdoor_request(
    "primary_count", "study", "poisson", y ~ z + x)
  retry <- .dsvert_formal_glm_frontdoor_request(
    "primary_count", "study", "poisson", y ~ x + z)

  expect_identical(first, retry)
  expect_setequal(names(first$value), c(
    "version", "analysis_id", "data_name", "family", "formula_sha256",
    "public_selectors", "privacy_controls", "role_selection"))
  expect_identical(unlist(first$value$public_selectors, use.names = FALSE),
                   c("analysis_id", "data_name", "family", "formula_sha256"))
  forbidden <- c(
    "epsilon", "delta", "bound", "clip", "seed", "noise", "ring",
    "precision", "backend", "role", "peer")
  # The two declarative ownership fields may name privacy/roles; they never
  # contain an analyst-selected value.  The selectors themselves must not.
  expect_false(any(grepl(
    paste(forbidden, collapse = "|"),
    unlist(first$value$public_selectors, use.names = FALSE),
    ignore.case = TRUE)))
  expect_false(any(c(
    "epsilon", "delta", "bounds", "seed", "compute_peers") %in%
    names(first$value)))
  expect_identical(first$value$privacy_controls,
                   "server_owned_not_transmitted_v1")
  expect_identical(first$value$role_selection,
                   "pinned_server_policy_not_transmitted_v1")
})

test_that("formal GLM rejects analyst-owned legacy knobs", {
  legacy_knob_forced <- FALSE
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "binomial",
    lambda = {
      legacy_knob_forced <<- TRUE
      stop("legacy knob was forced", call. = FALSE)
    }, formal_analysis_id = "primary_logit"), "lambda")
  expect_false(legacy_knob_forced)
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "poisson", ring = 127L,
    formal_analysis_id = "primary_count"), "ring")
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "poisson",
    numeric_backend = "multiprecision",
    formal_analysis_id = "primary_count"), "numeric_backend")
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "binomial", missing = "fail",
    formal_analysis_id = "primary_logit"), "missing")
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "gaussian",
    formal_analysis_id = "primary_logit"), "binomial.*poisson")
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "binomial",
    dp_analysis_id = "gaussian", formal_analysis_id = "primary_logit"),
    "mutually exclusive")
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "binomial",
    formal_analysis_id = "primary_logit",
    fresh_formal_analysis_id = "configured_logit"), "mutually exclusive")
})

test_that("fresh formal GLM selectors are sealed before DSI", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_formal_glm_registered_fresh_run = function(...) {
      calls <<- calls + 1L
      stop("must not run", call. = FALSE)
    },
    .dsvert_aggregate_strict = function(...) {
      calls <<- calls + 1L
      stop("must not retrieve", call. = FALSE)
    },
    .package = "dsVertClient")
  expect_error(ds.vertGLM(
    y ~ x, data = "study", family = "binomial",
    fresh_formal_analysis_id = "configured_logit",
    datasources = conns), class = "dsvert_route_unavailable")
  expect_error(ds.vert.glm(
    y ~ x, data = "study", family = "poisson",
    fresh_formal_analysis_id = "configured_count",
    datasources = conns), class = "dsvert_route_unavailable")
  expect_error(dsVertClient:::.dsvert_formal_glm_fresh_frontdoor_adapter(
    character(), y ~ x, "study", "binomial", FALSE, conns,
    "configured_logit"), class = "dsvert_route_unavailable")
  expect_identical(calls, 0L)
})

test_that("formal GLM alias preserves the certified-release route", {
  testthat::local_mocked_bindings(
    .dsvert_quarantine_test_mode = function() FALSE,
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      request <- as.list(expr[-1L])
      stats::setNames(
        list(.formal_glm_public_release(request)), names(conns))
    },
    .package = "dsVertClient")
  fit <- ds.vert.glm(
    y ~ x, data = "study", family = "poisson",
    formal_analysis_id = "primary_count", verbose = FALSE,
    datasources = list(site_a = structure(list(), class = "mock")))
  expect_s3_class(fit, "dsvert_formal_dp_glm")
  expect_identical(fit$family, "poisson")

  expect_error(ds.vert.glm(
    y ~ x, data = "study", family = "binomial", precision = "high",
    formal_analysis_id = "primary_logit",
    datasources = stop("must not be forced", call. = FALSE)),
    "server-owned")

  expect_error(ds.vert.glm(
    y ~ x, data = "study",
    datasources = stop("must not be forced", call. = FALSE)),
    class = "dsvert_route_unavailable")
})

test_that("formal GLM formula contract is narrow and canonical", {
  expect_error(.dsvert_formal_glm_frontdoor_request(
    "a", "study", "binomial", y ~ x:z), "Only additive")
  expect_error(.dsvert_formal_glm_frontdoor_request(
    "a", "study", "binomial", log(y) ~ x), "response")
  expect_error(.dsvert_formal_glm_frontdoor_request(
    "a", "study", "binomial", y ~ .), "data-dependent")
  expect_error(.dsvert_formal_glm_frontdoor_request(
    "bad id", "study", "binomial", y ~ x), "formal_analysis_id")
  expect_error(.dsvert_formal_glm_frontdoor_request(
    "a", "bad data", "binomial", y ~ x), "data")
})
