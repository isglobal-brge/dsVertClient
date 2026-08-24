.formal_ordinal_frequency <- function() {
  structure(list(
    released = TRUE,
    status = "ok",
    variable = "severity",
    levels = c("none", "mild", "severe"),
    counts = c(none = 24, mild = 9, severe = 3),
    coordinate_descriptor = list(dataset = "study"),
    release_sha256 = paste(rep("b", 64L), collapse = ""),
    sticky_noise = TRUE,
    sticky_replay = TRUE,
    intermediate_values_exposed = FALSE,
    additional_privacy_cost = c(epsilon = 0, delta = 0)),
    class = c("ds.vertDPFrequency", "list"))
}

test_that("intercept-only ordinal thresholds post-process a validated Frequency", {
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) {
      seen <<- x
      x
    },
    .package = "dsVertClient")

  frequency <- .formal_ordinal_frequency()
  fit <- ds.vertOrdinal(
    severity ~ 1, data = "study",
    levels_ordered = c("none", "mild", "severe"),
    frequency = frequency, verbose = FALSE)

  probabilities <- c(none = 24.5, mild = 9.5, severe = 3.5) / 37.5
  cumulative <- c(none = probabilities[["none"]],
                  mild = probabilities[["none"]] + probabilities[["mild"]])
  expect_s3_class(fit, "ds.vertOrdinal")
  expect_s3_class(fit, "dsvert_dp_frequency_ordinal")
  expect_identical(seen, frequency)
  expect_identical(fit$levels, c("none", "mild", "severe"))
  expect_equal(fit$probabilities, probabilities)
  expect_equal(fit$cumulative_probabilities, cumulative)
  expect_equal(fit$thresholds, stats::qlogis(cumulative))
  expect_null(fit$beta_po)
  expect_null(fit$covariance_po)
  expect_false(fit$source_values_exposed)
  expect_false(fit$production_ready)
  expect_identical(fit$additional_privacy_cost, c(epsilon = 0, delta = 0))
})

test_that("intercept-only ordinal resolves its signed Frequency by source owner", {
  calls <- list()
  frequency <- .formal_ordinal_frequency()
  testthat::local_mocked_bindings(
    ds.vertDPFrequency = function(data_name, variable, server = NULL,
                                  datasources = NULL) {
      calls <<- c(calls, list(list(
        data_name = data_name, variable = variable, server = server,
        datasources = datasources)))
      frequency
    },
    .dsvert_dp_frequency_contract = function(x) x,
    .package = "dsVertClient")

  conns <- list(peer_a = structure(list(), class = "mock_connection"))
  levels <- c("none", "mild", "severe")
  direct <- ds.vertOrdinal(
    severity ~ 1, data = "study", levels_ordered = levels,
    server = "peer_a", datasources = conns)
  alias <- ds.vert.ordinal(
    severity ~ 1, data = "study", levels_ordered = levels,
    server = "peer_a", datasources = conns)

  expect_s3_class(direct, "dsvert_dp_frequency_ordinal")
  expect_identical(alias$thresholds, direct$thresholds)
  expect_identical(calls, rep(list(list(
    data_name = "study", variable = "severity", server = "peer_a",
    datasources = conns)), 2L))
  expect_error(ds.vertOrdinal(
    severity ~ x, data = "study", levels_ordered = levels,
    server = "peer_a", datasources = conns), "intercept-only")
  expect_length(calls, 2L)
  expect_error(ds.vertOrdinal(severity ~ 1, data = "study",
                              levels_ordered = levels, datasources = conns),
               "requires an explicit source owner")
  expect_length(calls, 2L)
  expect_error(ds.vertOrdinal(
    severity ~ 1, levels_ordered = levels, frequency = frequency,
    server = "peer_a"), "legacy controls")
  expect_length(calls, 2L)
})

test_that("formal ordinal requires an explicit signed ordered domain", {
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) {
      calls <<- calls + 1L
      x
    },
    .package = "dsVertClient")

  frequency <- .formal_ordinal_frequency()
  expect_error(ds.vertOrdinal(severity ~ 1, frequency = frequency),
               "complete permutation")
  expect_error(ds.vertOrdinal(severity ~ x,
                              levels_ordered = c("none", "mild", "severe"),
                              frequency = frequency), "intercept-only")
  expect_error(ds.vertOrdinal(outcome ~ 1,
                              levels_ordered = c("none", "mild", "severe"),
                              frequency = frequency), "outcome")
  expect_error(ds.vertOrdinal(severity ~ 1,
                              levels_ordered = c("none", "mild", "other"),
                              frequency = frequency), "permutation")
  expect_error(ds.vertOrdinal(severity ~ 1,
                              levels_ordered = c("none", "mild", "severe"),
                              max_iter = 10, frequency = frequency),
               "legacy controls")
  expect_identical(calls, 3L)
})

test_that("formal ordinal fails closed on an invalid Frequency and alias", {
  frequency <- .formal_ordinal_frequency()
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) {
      stop("x must be a released, validated ds.vertDPFrequency object",
           call. = FALSE)
    },
    .package = "dsVertClient")
  expect_error(ds.vertOrdinal(severity ~ 1,
                              levels_ordered = c("none", "mild", "severe"),
                              frequency = frequency), "released, validated")

  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) x,
    .package = "dsVertClient")
  fit <- ds.vert.ordinal(severity ~ 1,
                         levels_ordered = c("none", "mild", "severe"),
                         frequency = frequency)
  expect_s3_class(fit, "dsvert_dp_frequency_ordinal")
  expect_false(fit$production_ready)
  expect_error(ds.vert.ordinal(severity ~ 1,
                               levels_ordered = c("none", "mild", "severe"),
                               frequency = frequency, datasources = list()),
               "does not accept datasources")
})

test_that("historical ordinal Newton name retains only the Frequency y ~ 1 route", {
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) x,
    .package = "dsVertClient")

  frequency <- .formal_ordinal_frequency()
  levels <- c("none", "mild", "severe")
  direct <- ds.vertOrdinal(
    severity ~ 1, data = "study", levels_ordered = levels,
    frequency = frequency)
  joint <- ds.vertOrdinalJointNewton(
    severity ~ 1, data = "study", levels_ordered = levels,
    frequency = frequency)

  expect_s3_class(joint, "dsvert_dp_frequency_ordinal")
  expect_identical(joint$thresholds, direct$thresholds)
  expect_identical(joint$probabilities, direct$probabilities)
  expect_identical(joint$additional_privacy_cost, c(epsilon = 0, delta = 0))
  expect_false(joint$source_values_exposed)
  expect_false(joint$production_ready)
  expect_identical(joint$called_via, "ds.vertOrdinalJointNewton_frequency")
  expect_error(ds.vertOrdinalJointNewton(
    severity ~ x, levels_ordered = levels, frequency = frequency),
    "intercept-only")
  expect_error(ds.vertOrdinalJointNewton(
    severity ~ 1, levels_ordered = levels, max_outer = 2L,
    frequency = frequency), "legacy controls")
  expect_error(ds.vertOrdinalJointNewton(
    severity ~ 1, levels_ordered = levels, frequency = frequency,
    datasources = list()), "legacy controls")
})
