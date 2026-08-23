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
    intermediate_values_exposed = FALSE),
    class = c("dsvert_formal_dp_glm", "ds.glm", "list"))
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
