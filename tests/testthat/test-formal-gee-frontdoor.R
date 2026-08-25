.formal_gee_certified_fit <- function(...) {
  arguments <- list(...)
  structure(list(
    family = arguments$family,
    coefficients = c("(Intercept)" = 0.5, x = -0.25),
    artifact_id = paste(rep("a", 64L), collapse = ""),
    certificate_sha256 = paste(rep("b", 64L), collapse = ""),
    formal_analysis_id = arguments$formal_analysis_id,
    formula_sha256 = paste(rep("c", 64L), collapse = ""),
    production_ready = FALSE,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    called_via = "ds.vertGLM_formal_analysis_id"),
    class = c("dsvert_formal_dp_glm", "ds.glm", "list"))
}

.formal_gee_fresh_fit <- function(...) {
  arguments <- list(...)
  fit <- .formal_gee_certified_fit(
    family = arguments$family,
    formal_analysis_id = arguments$fresh_formal_analysis_id)
  fit$called_via <- "ds.vertGLM_fresh_formal_analysis_id"
  fit
}

.formal_gee_gaussian_fit <- function(...) {
  arguments <- list(...)
  structure(list(
    family = "gaussian",
    coefficients = c("(Intercept)" = 1.25, x = -0.5),
    analysis_id = arguments$dp_analysis_id,
    certificate_sha256 = paste(rep("d", 64L), collapse = ""),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0)),
    class = c("ds.vertDPGaussian", "list"))
}

test_that("formal independent GEE consumes one certified GLM point release", {
  seen <- NULL
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) {
      seen <<- list(...)
      .formal_gee_certified_fit(...)
    },
    .package = "dsVertClient")

  fit <- ds.vertGEE(
    y ~ x, data = "study", family = "binomial",
    corstr = "independence", formal_analysis_id = "gee-logit",
    verbose = FALSE, datasources = list(site_a = structure(list(), class = "mock")))

  expect_s3_class(fit, "ds.vertGEE")
  expect_s3_class(fit, "dsvert_formal_dp_gee")
  expect_identical(fit$corstr, "independence")
  expect_identical(fit$coefficients, c("(Intercept)" = 0.5, x = -0.25))
  expect_null(fit$robust_covariance)
  expect_null(fit$std_errors)
  expect_false(fit$cluster_correlation_estimated)
  expect_false(fit$source_values_exposed)
  expect_false(fit$production_ready)
  expect_match(fit$inference, "unavailable")
  expect_identical(seen$formal_analysis_id, "gee-logit")
  expect_identical(seen$family, "binomial")
  expect_identical(seen$data, "study")
})

test_that("fresh formal independent GEE runs only the matching durable GLM", {
  seen <- NULL
  glm_calls <- 0L
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) {
      glm_calls <<- glm_calls + 1L
      seen <<- list(...)
      .formal_gee_fresh_fit(...)
    },
    .package = "dsVertClient")

  fit <- ds.vertGEE(
    y ~ x, data = "study", family = "poisson",
    corstr = "independence", fresh_formal_analysis_id = "gee-count",
    verbose = FALSE, datasources = list(site_a = structure(list(), class = "mock")))

  expect_s3_class(fit, "dsvert_formal_dp_gee")
  expect_identical(fit$called_via, "ds.vertGEE_fresh_formal_analysis_id")
  expect_identical(fit$fresh_formal_analysis_id, "gee-count")
  expect_null(fit$formal_analysis_id)
  expect_identical(seen$fresh_formal_analysis_id, "gee-count")
  expect_identical(seen$family, "poisson")
  expect_false(fit$source_values_exposed)
  expect_false(fit$production_ready)

  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "poisson",
    formal_analysis_id = "existing", fresh_formal_analysis_id = "gee-count"),
    "mutually exclusive")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "poisson", lambda = 0,
    fresh_formal_analysis_id = "gee-count"),
    "fresh_formal_analysis_id GEE does not accept legacy controls")
  expect_identical(glm_calls, 1L)
})

test_that("independent Gaussian GEE consumes one signed Gaussian Synopsis", {
  calls <- list()
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) {
      calls <<- c(calls, list(list(...)))
      .formal_gee_gaussian_fit(...)
    },
    .package = "dsVertClient")

  conns <- list(site_a = structure(list(), class = "mock"))
  direct <- ds.vertGEE(
    y ~ x, data = "study", family = "gaussian",
    dp_analysis_id = "gee-gaussian", datasources = conns, verbose = FALSE)
  alias <- ds.vert.gee(
    y ~ x, data = "study", dp_analysis_id = "gee-gaussian",
    datasources = conns)

  expect_s3_class(direct, "dsvert_dp_gaussian_gee")
  expect_s3_class(alias, "ds.vertGEE")
  expect_identical(alias$coefficients, direct$coefficients)
  expect_identical(calls[[1L]]$dp_analysis_id, "gee-gaussian")
  expect_identical(calls[[1L]]$family, "gaussian")
  expect_identical(calls[[1L]]$data, "study")
  expect_null(direct$robust_covariance)
  expect_null(direct$std_errors)
  expect_false(direct$cluster_correlation_estimated)
  expect_false(direct$source_values_exposed)
  expect_identical(direct$additional_privacy_cost, c(epsilon = 0, delta = 0))
  expect_length(calls, 2L)

  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "binomial",
    dp_analysis_id = "gee-gaussian", datasources = conns), "gaussian")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", corstr = "exchangeable",
    dp_analysis_id = "gee-gaussian", datasources = conns), "independence")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", id_col = "patient",
    dp_analysis_id = "gee-gaussian", datasources = conns), "cluster")
  expect_length(calls, 2L)
})

test_that("formal independent GEE rejects cluster and legacy controls before GLM", {
  glm_calls <- 0L
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) {
      glm_calls <<- glm_calls + 1L
      .formal_gee_certified_fit(...)
    },
    .package = "dsVertClient")

  base <- list(
    formula = y ~ x, data = "study", family = "poisson",
    formal_analysis_id = "gee-count", datasources = list())
  expect_error(do.call(ds.vertGEE, c(base, list(corstr = "exchangeable"))),
               "independence")
  expect_error(do.call(ds.vertGEE, c(base, list(id_col = "patient"))),
               "does not accept cluster")
  expect_error(do.call(ds.vertGEE, c(base, list(order_col = "visit"))),
               "does not accept cluster")
  expect_error(do.call(ds.vertGEE, c(base, list(lambda = 0))),
               "does not accept legacy controls")
  gaussian <- base
  gaussian$family <- "gaussian"
  expect_error(do.call(ds.vertGEE, gaussian),
               "binomial.*poisson")
  expect_identical(glm_calls, 0L)
})

test_that("the compatibility GEE alias preserves the formal point boundary", {
  testthat::local_mocked_bindings(
    ds.vertGLM = .formal_gee_certified_fit,
    .package = "dsVertClient")

  fit <- ds.vert.gee(
    y ~ x, data = "study", family = "poisson",
    formal_analysis_id = "gee-count", datasources = list())
  expect_s3_class(fit, "ds.vertGEE")
  expect_identical(fit$corstr, "independence")
  expect_false(fit$production_ready)

  expect_error(
    ds.vert.gee(
      y ~ x, data = "study", family = "poisson", precision = "high",
      formal_analysis_id = "gee-count", datasources = list()),
    "does not accept legacy controls")
})

test_that("the compatibility GEE alias preserves the fresh formal boundary", {
  testthat::local_mocked_bindings(
    ds.vertGLM = .formal_gee_fresh_fit,
    .package = "dsVertClient")

  fit <- ds.vert.gee(
    y ~ x, data = "study", family = "binomial",
    fresh_formal_analysis_id = "gee-logit", datasources = list())
  expect_s3_class(fit, "ds.vertGEE")
  expect_identical(fit$fresh_formal_analysis_id, "gee-logit")
  expect_identical(fit$route, "ds.vertGEE")
  expect_error(ds.vert.gee(
    y ~ x, data = "study", family = "binomial", precision = "high",
    fresh_formal_analysis_id = "gee-logit", datasources = list()),
    "does not accept legacy controls")
})
