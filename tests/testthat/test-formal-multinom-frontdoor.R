.formal_multinom_frequency <- function() {
  structure(list(
    released = TRUE,
    status = "ok",
    variable = "status",
    levels = c("control", "case", "other"),
    counts = c(control = 20, case = 10, other = 5),
    coordinate_descriptor = list(dataset = "study"),
    release_sha256 = paste(rep("a", 64L), collapse = ""),
    sticky_noise = TRUE,
    sticky_replay = TRUE,
    intermediate_values_exposed = FALSE,
    additional_privacy_cost = c(epsilon = 0, delta = 0)),
    class = c("ds.vertDPFrequency", "list"))
}

test_that("intercept-only multinomial post-processes a validated Frequency", {
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) {
      seen <<- x
      x
    },
    .package = "dsVertClient")

  frequency <- .formal_multinom_frequency()
  fit <- ds.vertMultinom(
    status ~ 1, data = "study", reference = "control",
    frequency = frequency, verbose = FALSE)

  expected <- c(control = 20.5, case = 10.5, other = 5.5) / 36.5
  expect_s3_class(fit, "ds.vertMultinom")
  expect_s3_class(fit, "dsvert_dp_frequency_multinom")
  expect_identical(seen, frequency)
  expect_identical(fit$classes, c("case", "other"))
  expect_identical(fit$reference, "control")
  expect_equal(fit$probabilities, expected)
  expect_equal(drop(fit$coefficients),
               log(expected[c("case", "other")] / expected[["control"]]))
  expect_null(fit$std_errors)
  expect_false(fit$source_values_exposed)
  expect_false(fit$production_ready)
  expect_identical(fit$additional_privacy_cost, c(epsilon = 0, delta = 0))
})

test_that("intercept-only multinomial resolves its signed Frequency by source owner", {
  calls <- list()
  frequency <- .formal_multinom_frequency()
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
  direct <- ds.vertMultinom(
    status ~ 1, data = "study", reference = "control",
    server = "peer_a", datasources = conns)
  alias <- ds.vert.multinom(
    status ~ 1, data = "study", reference = "control",
    server = "peer_a", datasources = conns)

  expect_s3_class(direct, "dsvert_dp_frequency_multinom")
  expect_identical(alias$coefficients, direct$coefficients)
  expect_identical(calls, rep(list(list(
    data_name = "study", variable = "status", server = "peer_a",
    datasources = conns)), 2L))
  expect_error(ds.vertMultinom(
    status ~ x, data = "study", server = "peer_a", datasources = conns),
    "intercept-only")
  expect_length(calls, 2L)
  expect_error(ds.vertMultinom(status ~ 1, data = "study", datasources = conns),
               "requires an explicit source owner")
  expect_length(calls, 2L)
  expect_error(ds.vertMultinom(
    status ~ 1, frequency = frequency, server = "peer_a"),
    "legacy controls")
  expect_length(calls, 2L)
})

test_that("formal multinomial rejects unsupported designs before Frequency", {
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) {
      calls <<- calls + 1L
      x
    },
    .package = "dsVertClient")

  frequency <- .formal_multinom_frequency()
  expect_error(ds.vertMultinom(status ~ x, frequency = frequency),
               "intercept-only")
  expect_error(ds.vertMultinom(outcome ~ 1, frequency = frequency),
               "outcome")
  expect_error(ds.vertMultinom(
    status ~ 1, classes = c("case", "control", "other"),
    frequency = frequency), "classes")
  expect_error(ds.vertMultinom(status ~ 1, max_iter = 10, frequency = frequency),
               "legacy controls")
  # Formula and legacy controls fail before validation; outcome/category
  # mismatches are checked only against the authenticated Frequency metadata.
  expect_identical(calls, 2L)
})

test_that("formal multinomial fails closed on an invalid Frequency and alias", {
  frequency <- .formal_multinom_frequency()
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) {
      stop("x must be a released, validated ds.vertDPFrequency object",
           call. = FALSE)
    },
    .package = "dsVertClient")
  expect_error(ds.vertMultinom(status ~ 1, frequency = frequency),
               "released, validated")

  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) x,
    .package = "dsVertClient")
  fit <- ds.vert.multinom(status ~ 1, frequency = frequency)
  expect_s3_class(fit, "dsvert_dp_frequency_multinom")
  expect_false(fit$production_ready)
  expect_error(ds.vert.multinom(status ~ 1, frequency = frequency,
                                datasources = list()),
               "does not accept datasources")
})

test_that("historical joint multinomial names retain only the Frequency y ~ 1 route", {
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) x,
    .package = "dsVertClient")

  frequency <- .formal_multinom_frequency()
  direct <- ds.vertMultinom(
    status ~ 1, data = "study", classes = frequency$levels,
    reference = "control", frequency = frequency)
  joint <- ds.vertMultinomJoint(
    status ~ 1, data = "study", levels = frequency$levels,
    frequency = frequency)
  newton <- ds.vertMultinomJointNewton(
    status ~ 1, data = "study", levels = frequency$levels,
    frequency = frequency)

  for (fit in list(joint, newton)) {
    expect_s3_class(fit, "dsvert_dp_frequency_multinom")
    expect_identical(fit$coefficients, direct$coefficients)
    expect_identical(fit$probabilities, direct$probabilities)
    expect_identical(fit$additional_privacy_cost, c(epsilon = 0, delta = 0))
    expect_false(fit$source_values_exposed)
    expect_false(fit$production_ready)
  }
  expect_identical(joint$called_via, "ds.vertMultinomJoint_frequency")
  expect_identical(newton$called_via, "ds.vertMultinomJointNewton_frequency")
  expect_error(ds.vertMultinomJoint(status ~ x, frequency = frequency),
               "intercept-only")
  expect_error(ds.vertMultinomJointNewton(
    status ~ 1, levels = frequency$levels, max_outer = 2L,
    frequency = frequency), "legacy controls")
  expect_error(ds.vertMultinomJoint(
    status ~ 1, levels = frequency$levels, frequency = frequency,
    datasources = list()), "legacy controls")
})
