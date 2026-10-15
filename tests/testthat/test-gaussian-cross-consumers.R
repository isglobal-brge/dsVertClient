.gaussian_lasso_consumers <- function(fit, pins = NULL) list(
  path = function() ds.vertLASSO(fit, .02, alpha_grid = c(1, .5), trusted_pinset = pins),
  one_step = function() ds.vertLASSO1Step(fit, c(.02, .01), trusted_pinset = pins),
  proximal = function() ds.vertLASSOProximal(fit, .02, trusted_pinset = pins),
  selection = function() ds.vertLASSOCV(fit, c(.02, .01), trusted_pinset = pins),
  alias_path = function() ds.vert.lasso(fit, .02, alpha_grid = c(1, .5), trusted_pinset = pins),
  alias_one_step = function() ds.vert.lasso_1step(fit, c(.02, .01), trusted_pinset = pins),
  alias_proximal = function() ds.vert.lasso_proximal(fit, .02, trusted_pinset = pins),
  alias_selection = function() ds.vert.lasso_cv(fit, c(.02, .01), trusted_pinset = pins))

.gaussian_region_consumers <- function(fit) list(
  confint = function() ds.vertConfint(fit, type = "mechanism"),
  wald = function() ds.vertWald(fit, "x", type = "mechanism"),
  contrast = function() ds.vertContrast(fit, c("x"), type = "mechanism"),
  lr = function() ds.vertLR("(Intercept)", fit, type = "mechanism"),
  alias_confint = function() ds.vert.confint(fit, type = "mechanism"),
  alias_wald = function() ds.vert.wald(fit, "x", type = "mechanism"),
  alias_contrast = function() ds.vert.contrast(fit, c("x"), type = "mechanism"),
  alias_lr = function() ds.vert.lr("(Intercept)", fit, type = "mechanism"))

.gaussian_consumers_no_calls <- function(code) {
  fail <- function(...) stop("unexpected private operation", call. = FALSE)
  testthat::local_mocked_bindings(
    .dsvert_dp_synopsis_vector_run = fail,
    .dsvert_dp_capsule_vector_run = fail,
    .dsvert_aggregate_strict = fail, .dsvert_fanout_by_site = fail,
    .dsvert_exact_gc_run = fail, .package = "dsVertClient")
  testthat::local_mocked_bindings(datashield.aggregate = fail,
    datashield.assign = fail, datashield.connections_find = fail, .package = "DSI")
  force(code)
}

test_that("current cross Gaussian certificates serve every LASSO consumer at K=2,3,5", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .gaussian_current_fixture(k)
    .with_gaussian_current_fixture(fixture, {
      fit <- .gaussian_current_fit(fixture)
      source_bytes <- serialize(fit, NULL, version = 3)
      .gaussian_consumers_no_calls({
        for (call in .gaussian_lasso_consumers(fit, fixture$pins)) {
          result <- call()
          expect_identical(serialize(call(), NULL, version = 3), serialize(result, NULL, version = 3))
          expect_identical(result$additional_privacy_cost, c(epsilon = 0, delta = 0))
          expect_identical(result$additional_server_calls_after_synopsis, 0L)
          expect_identical(result$input_provenance, "signed_dp_gaussian_synopsis")
          expect_identical(result$source_certificate_validation$authenticity, "caller_anchored")
          expect_false(result$inference$sampling_inference_available)
        }
        proximal <- ds.vertLASSOProximal(fit, .02)
        expect_equal(proximal$coefficients, c(`(Intercept)` = .5, x = .12), tolerance = 1e-7)
        expect_true(proximal$kkt$satisfied)
        expect_false(proximal$client_psd_projection_applied)
        expect_false(proximal$implicit_ridge)
        expect_equal(proximal$XtX_over_n, fit$sufficient_statistics_dp$gram_projected / fit$n_obs)
        expect_equal(proximal$cross_over_n, fit$sufficient_statistics_dp$cross_projected / fit$n_obs)
        selected <- ds.vertLASSOCV(fit, c(.02, .01))
        expect_false(selected$cross_validation)
        expect_identical(selected$selection_method, "DP_projected_pseudo_information_criterion")
      })
      expect_identical(serialize(fit, NULL, version = 3), source_bytes)
    })
  }
})

test_that("Gaussian iterative LASSO and its alias consume current cross fits without another release", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .gaussian_current_fixture(k)
    .with_gaussian_current_fixture(fixture, {
      fit <- .gaussian_current_fit(fixture)
      calls <- 0L
      testthat::local_mocked_bindings(ds.vertGLM = function(...) {
        calls <<- calls + 1L
        fit
      }, .package = "dsVertClient")
      .gaussian_consumers_no_calls({
        for (frontdoor in list(ds.vertLASSOIter, ds.vert.lasso_iter)) {
          result <- frontdoor(y ~ x, data = "outcome_data", family = "gaussian",
            lambda = c(.02, .01), dp_analysis_id = "cross_model", verbose = FALSE,
            datasources = fixture$conns)
          expect_s3_class(result, "ds.vertLASSOIter")
          expect_identical(result$additional_server_calls_after_synopsis, 0L)
          expect_identical(result$additional_privacy_cost, c(epsilon = 0, delta = 0))
          expect_true(all(vapply(result$path_certificates, function(cert) cert$kkt$satisfied, logical(1L))))
          expect_identical(result$final_fit, fit)
        }
      })
      expect_identical(calls, 2L)
    })
  }
})

test_that("current cross Gaussian mechanism regions and aliases work at K=2,3,5", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .gaussian_current_fixture(k)
    .with_gaussian_current_fixture(fixture, {
      fit <- .gaussian_current_fit(fixture)
      .gaussian_consumers_no_calls({
        for (call in .gaussian_region_consumers(fit)) {
          result <- call()
          expect_identical(serialize(call(), NULL, version = 3), serialize(result, NULL, version = 3))
          if (is.data.frame(result)) {
            expect_false(attr(result, "sampling_inference"))
            expect_identical(attr(result, "additional_server_calls"), 0L)
            expect_true(all(is.finite(result$lower) & is.finite(result$upper)))
          } else {
            expect_false(result$sampling_inference)
            expect_identical(result$additional_server_calls, 0L)
            expect_null(result$p_value)
            expect_true(all(is.finite(result$lower) & is.finite(result$upper)))
          }
        }
        confint <- ds.vertConfint(fit, type = "mechanism")
        wald <- ds.vertWald(fit, "x", type = "mechanism")
        contrast <- ds.vertContrast(fit, "x", type = "mechanism")
        expect_equal(wald$lower, confint["x", "lower"])
        expect_equal(unname(contrast$lower), confint["x", "lower"])
        lr <- ds.vertLR("(Intercept)", fit, type = "mechanism")
        expect_equal(lr$statistic, 1.6, tolerance = 1e-10)
        expect_error(ds.vertConfint(fit), "sampling|covariance|inference")
        expect_error(ds.vertWald(fit, "x"), "sampling|covariance|inference")
        expect_error(ds.vertContrast(fit, "x"), "sampling|covariance|inference")
        expect_error(ds.vertLR(fit, fit), "sampling|covariance|inference")
      })
    })
  }
})

.gaussian_consumer_reseal <- function(fit) {
  certificate <- fit$provenance_certificate
  certificate$certificate_sha256 <- .dsvert_dp_gaussian_certificate_hash(
    certificate[setdiff(names(certificate), "certificate_sha256")])
  fit$provenance_certificate <- certificate
  fit
}

test_that("Gaussian consumers revalidate moments, bounds, roots and authority attestations", {
  fixture <- .gaussian_current_fixture(3L)
  .with_gaussian_current_fixture(fixture, {
    fit <- .gaussian_current_fit(fixture)
    changed <- list()
    value <- fit; value$n_obs <- value$n_obs + 1; changed$count <- value
    value <- fit; value$sufficient_statistics_dp$gram_projected[1, 1] <- 0; changed$gram <- value
    value <- fit; value$sufficient_statistics_dp$cross_projected[[1L]] <- 0; changed$cross <- value
    value <- fit; value$signed_artifact$predictors$x$lower <- -2; changed$bounds <- value
    value <- fit; value$signed_artifact$implementation_state <- "same_owner_materialized"; changed$state <- value
    value <- fit; value$contract_sha256 <- strrep("0", 64L); changed$root <- value
    for (field in c("signature", "exact_transcript_sha256", "private_layout_sha256",
        "numeric_certificate_sha256", "private_result_exposed")) {
      value <- fit
      evidence <- .dsvert_joint_dp_client_decode(
        value$provenance_certificate$signed_evidence$cross_gaussian_evidence_json,
        "test evidence", 2^24)
      evidence[[1L]][[field]] <- if (identical(field, "private_result_exposed")) TRUE else strrep("0", 64L)
      value$provenance_certificate$signed_evidence$cross_gaussian_evidence_json <-
        .dsvert_joint_dp_client_json(evidence)
      changed[[field]] <- .gaussian_consumer_reseal(value)
    }
    .gaussian_consumers_no_calls({
      for (value in changed) {
        consumers <- c(.gaussian_lasso_consumers(value, fixture$pins), .gaussian_region_consumers(value))
        for (call in consumers) expect_error(call())
      }
      wrong_pins <- fixture$pins
      wrong_pins[[1L]] <- fixture$pins[[2L]]
      for (call in .gaussian_lasso_consumers(fit, wrong_pins)) expect_error(call(), "pinset|trusted")
      hash <- .dsvert_dp_pinset_hash(fixture$pins)
      rm(list = hash, envir = .dsvert_dp_gaussian_trust_cache)
      for (call in .gaussian_lasso_consumers(fit)) expect_error(call(), "anchored")
      for (call in .gaussian_region_consumers(fit)) expect_error(call(), "anchored")
      for (call in .gaussian_lasso_consumers(fit, fixture$pins)) expect_s3_class(call(), "list")
    })
  })
})

test_that("cross Gaussian LASSO rejects historical capsules despite valid old signatures", {
  legacy <- new.env(parent = environment())
  for (node in parse(testthat::test_path("test-dp-gaussian.R"))) {
    if (is.call(node) && identical(node[[1L]], as.name("<-"))) eval(node, legacy)
  }
  fixture <- legacy$.dp_gaussian_cross_frontdoor_fixture(2L)
  fit <- testthat::with_mocked_bindings(
    ds.vertDPGaussian("outcome_data", "cross_model", datasources = fixture$conns),
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(legacy$.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    }, .package = "dsVertClient")
  expect_true(ds.validateDPGaussianCertificate(fit, trusted_pinset = fixture$pins)$integrity_valid)
  .gaussian_consumers_no_calls({
    for (call in .gaussian_lasso_consumers(fit, fixture$pins)) expect_error(call(), "current authenticated Synopsis")
  })
})
