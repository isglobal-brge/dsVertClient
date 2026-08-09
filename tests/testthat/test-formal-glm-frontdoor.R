.formal_glm_condition <- function(expr) {
  tryCatch(expr, error = function(error) error)
}

test_that("formal GLM is sealed before datasource forcing or DSI", {
  datasource_forced <- FALSE
  condition <- .formal_glm_condition(ds.vertGLM(
    y ~ x, data = "study", family = "binomial",
    formal_analysis_id = "primary_logit", verbose = FALSE,
    datasources = {
      datasource_forced <<- TRUE
      stop("datasource promise was forced", call. = FALSE)
    }))

  expect_s3_class(condition, "dsvert_formal_glm_frontdoor_unavailable")
  expect_identical(
    condition$code,
    "formal_glm_phase19_durable_r_dsi_release_bridge_not_promoted")
  expect_identical(condition$dsi_calls, 0L)
  expect_identical(condition$openings_performed, 0L)
  expect_false(condition$operation_limit)
  expect_false(condition$request_limit)
  expect_false(condition$history_can_deny_operation)
  expect_false(condition$production_ready)
  expect_false(datasource_forced)
  expect_identical(length(condition$missing), 5L)
  expect_true(
    "registered_r_dsi_lifecycle_for_phase18_source_materialization_v1" %in%
      condition$missing)
  expect_false(any(grepl(
    "phase18_registry_to_materializer", condition$missing, fixed = TRUE)))
  expect_match(condition$request_sha256, "^[0-9a-f]{64}$")
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
})

test_that("formal GLM alias preserves the zero-DSI gate", {
  testthat::local_mocked_bindings(
    .dsvert_quarantine_test_mode = function() FALSE,
    .package = "dsVertClient")
  datasource_forced <- FALSE
  condition <- .formal_glm_condition(ds.vert.glm(
    y ~ x, data = "study", family = "poisson",
    formal_analysis_id = "primary_count", verbose = FALSE,
    datasources = {
      datasource_forced <<- TRUE
      stop("datasource promise was forced", call. = FALSE)
    }))
  expect_s3_class(condition, "dsvert_formal_glm_frontdoor_unavailable")
  expect_false(datasource_forced)

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
