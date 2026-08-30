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

.formal_gee_grid_fit <- function(...) {
  arguments <- list(...)
  structure(list(
    family = paste0(arguments$family, "_finite_grid"),
    coefficients = c("(Intercept)" = 0.25, x = 0.5),
    analysis_id = arguments$analysis_id,
    signed_artifact = list(
      version = dsVertClient:::.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS[[
        arguments$family]]),
    certificate_sha256 = paste(rep("e", 64L), collapse = ""),
    covariance = NULL, std_errors = NULL,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0)),
    class = c("dsvert_dp_glm_grid", "ds.glm", "list"))
}

.formal_gee_gaussian_exchangeable_fit <- function(...) {
  structure(list(
    status = "ok",
    signed_artifact = list(
      version =
        "bounded-normalized-random-intercept-fixed-sufficient-statistics-v2",
      outcome = list(column = "y"), cluster = list(column = "patient"),
      predictor_order = "x",
      estimation_scope =
        "bounded_random_intercept_GLS_fixed_effects_finite_signed_variance_ratio_grid_ML_profile_v1"),
    sigma2 = 2, sigma_b2 = 1,
    coefficients = c(`(Intercept)` = 1.25, x = 0.75),
    certificate_sha256 = paste(rep("a", 64L), collapse = ""),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    additional_server_calls_after_synopsis = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0)),
    class = c("ds.vertDPLMM", "list"))
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

test_that("fresh formal independent GEE is sealed before DSI", {
  glm_calls <- 0L
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) {
      glm_calls <<- glm_calls + 1L
      stop("must not run", call. = FALSE)
    },
    .package = "dsVertClient")

  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "poisson",
    corstr = "independence", fresh_formal_analysis_id = "gee-count",
    verbose = FALSE, datasources = list(site_a = structure(list(), class = "mock"))),
    class = "dsvert_route_unavailable")

  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "poisson",
    formal_analysis_id = "existing", fresh_formal_analysis_id = "gee-count"),
    "mutually exclusive")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "poisson", lambda = 0,
    fresh_formal_analysis_id = "gee-count"),
    class = "dsvert_route_unavailable")
  expect_error(ds.vert.gee(
    y ~ x, data = "study", family = "binomial",
    fresh_formal_analysis_id = "gee-logit", datasources = list()),
    class = "dsvert_route_unavailable")
  expect_identical(glm_calls, 0L)
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

test_that("Gaussian exchangeable GEE consumes only a matching signed random-intercept GLS artifact", {
  seen <- NULL
  testthat::local_mocked_bindings(
    ds.vertDPLMM = function(...) {
      seen <<- list(...)
      .formal_gee_gaussian_exchangeable_fit(...)
    },
    .package = "dsVertClient")
  conns <- list(site_a = list(), site_b = list())
  direct <- ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", id_col = "patient",
    corstr = "exchangeable", analysis_id = "gee-gaussian-clustered",
    datasources = conns, verbose = FALSE)
  alias <- ds.vert.gee(
    y ~ x, data = "study", family = "gaussian", id_col = "patient",
    corstr = "exchangeable", analysis_id = "gee-gaussian-clustered",
    datasources = conns)
  expect_s3_class(direct, "dsvert_dp_gaussian_exchangeable_gee")
  expect_s3_class(alias, "ds.vertGEE")
  expect_identical(seen$analysis_id, "gee-gaussian-clustered")
  expect_identical(direct$coefficients,
                   c(`(Intercept)` = 1.25, x = 0.75))
  expect_equal(direct$working_correlation, 1 / 3)
  expect_identical(direct$correlation_estimation,
                   "signed_random_intercept_variance_components")
  expect_true(direct$cluster_correlation_estimated)
  expect_identical(direct$cluster_columns, "patient")
  expect_null(direct$robust_covariance)
  expect_null(direct$std_errors)
  expect_false(direct$source_values_exposed)
  expect_false(direct$intermediate_values_exposed)
  expect_identical(direct$additional_server_calls_after_synopsis, 0L)
  expect_identical(direct$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_identical(alias$frontdoor, "ds.vert.gee")
  expect_output(print(direct), "exchangeable-working GEE")

  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", id_col = "patient",
    corstr = "ar1", analysis_id = "gee-gaussian-clustered",
    datasources = conns), "requires distinct id_col and order_col")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", corstr = "exchangeable",
    analysis_id = "gee-gaussian-clustered", datasources = conns), "id_col")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", id_col = "patient",
    corstr = "exchangeable", analysis_id = "gee-gaussian-clustered",
    lambda = 0, datasources = conns), "legacy controls")
})

test_that("Gaussian exchangeable GEE rejects a substituted or malformed LMM result", {
  bad_artifact <- .formal_gee_gaussian_exchangeable_fit()
  bad_artifact$signed_artifact$version <-
    "bounded-gaussian-random-slope-likelihood-grid-v1"
  testthat::with_mocked_bindings(
    ds.vertDPLMM = function(...) bad_artifact,
    expect_error(ds.vertGEE(
      y ~ x, data = "study", family = "gaussian", id_col = "patient",
      corstr = "exchangeable", analysis_id = "gee-gaussian-clustered",
      datasources = list()), "does not match"),
    .package = "dsVertClient")

  bad_variance <- .formal_gee_gaussian_exchangeable_fit()
  bad_variance$sigma2 <- -1
  testthat::with_mocked_bindings(
    ds.vertDPLMM = function(...) bad_variance,
    expect_error(ds.vertGEE(
      y ~ x, data = "study", family = "gaussian", id_col = "patient",
      corstr = "exchangeable", analysis_id = "gee-gaussian-clustered",
      datasources = list()), "cannot support"),
    .package = "dsVertClient")

  singular_correlation <- .formal_gee_gaussian_exchangeable_fit()
  singular_correlation$sigma2 <- 0
  testthat::with_mocked_bindings(
    ds.vertDPLMM = function(...) singular_correlation,
    {
      clamped <- ds.vertGEE(
        y ~ x, data = "study", family = "gaussian", id_col = "patient",
        corstr = "exchangeable", analysis_id = "gee-gaussian-clustered",
        datasources = list())
      expect_equal(clamped$working_correlation, 1 - 2^-16)
      expect_equal(clamped$working_correlation_raw, 1)
      expect_true(clamped$working_correlation_clamped)
    },
    .package = "dsVertClient")
})

test_that("independent binomial and Poisson GEE consumes one signed finite grid", {
  calls <- list()
  testthat::local_mocked_bindings(
    ds.vertGLM = function(...) {
      calls <<- c(calls, list(list(...)))
      .formal_gee_grid_fit(...)
    },
    .package = "dsVertClient")

  conns <- list(site_a = structure(list(), class = "mock"))
  direct <- ds.vertGEE(
    y ~ x, data = "study", family = "binomial",
    analysis_id = "gee-grid", datasources = conns, verbose = FALSE)
  alias <- ds.vert.gee(
    y ~ x, data = "study", family = "poisson",
    analysis_id = "gee-grid-poisson", datasources = conns)

  expect_s3_class(direct, "dsvert_dp_glm_grid_gee")
  expect_s3_class(alias, "ds.vertGEE")
  expect_identical(direct$coefficients, c("(Intercept)" = 0.25, x = 0.5))
  expect_identical(alias$coefficients, direct$coefficients)
  expect_identical(direct$analysis_id, "gee-grid")
  expect_identical(alias$analysis_id, "gee-grid-poisson")
  expect_identical(calls[[1L]]$analysis_id, "gee-grid")
  expect_identical(calls[[2L]]$analysis_id, "gee-grid-poisson")
  expect_null(direct$robust_covariance)
  expect_null(direct$std_errors)
  expect_false(direct$cluster_correlation_estimated)
  expect_false(direct$source_values_exposed)
  expect_identical(direct$additional_privacy_cost, c(epsilon = 0, delta = 0))

  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", analysis_id = "gee-grid",
    datasources = conns), "Gaussian analysis_id GEE supports only")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "binomial", analysis_id = "gee-grid",
    id_col = "patient", order_col = "visit", datasources = conns),
    "requires one id_col and no order_col")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "binomial", analysis_id = "gee-grid",
    lambda = 0, datasources = conns), "does not accept legacy controls")
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
               "formal_analysis_id GEE supports only")
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
