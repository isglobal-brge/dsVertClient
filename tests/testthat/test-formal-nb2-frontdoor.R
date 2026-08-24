.formal_nb2_frequency <- function() {
  structure(list(
    released = TRUE,
    status = "ok",
    variable = "admissions",
    levels = c("0", "1", "2", "3"),
    counts = c(`0` = 50, `1` = 10, `2` = 5, `3` = 15),
    coordinate_descriptor = list(dataset = "study"),
    release_sha256 = paste(rep("c", 64L), collapse = ""),
    sticky_noise = TRUE,
    sticky_replay = TRUE,
    intermediate_values_exposed = FALSE,
    additional_privacy_cost = c(epsilon = 0, delta = 0)),
    class = c("ds.vertDPFrequency", "list"))
}

test_that("intercept-only NB2 post-processes a validated count Frequency", {
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) {
      seen <<- x
      x
    },
    .package = "dsVertClient")

  frequency <- .formal_nb2_frequency()
  fit <- ds.vertNBFullRegTheta(
    admissions ~ 1, data = "study", frequency = frequency, verbose = FALSE)

  mu <- (10 + 2 * 5 + 3 * 15) / 80
  variance <- (50 * (0 - mu)^2 + 10 * (1 - mu)^2 +
    5 * (2 - mu)^2 + 15 * (3 - mu)^2) / 80
  expect_s3_class(fit, "ds.vertNBFullRegTheta")
  expect_s3_class(fit, "dsvert_dp_frequency_nb2")
  expect_identical(seen, frequency)
  expect_equal(fit$coefficients[["(Intercept)"]], log(mu))
  expect_equal(fit$mean, mu)
  expect_equal(fit$variance, variance)
  expect_equal(fit$theta, mu^2 / (variance - mu))
  expect_output(print(fit), "sticky-DP NB2 intercept-only")
  expect_false(fit$source_values_exposed)
  expect_false(fit$production_ready)
  expect_identical(fit$additional_privacy_cost, c(epsilon = 0, delta = 0))
})

test_that("intercept-only NB2 resolves its signed count Frequency by source owner", {
  calls <- list()
  frequency <- .formal_nb2_frequency()
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
  direct <- ds.vertNBFullRegTheta(
    admissions ~ 1, data = "study", server = "peer_a", datasources = conns)
  alias <- ds.vert.nb(
    admissions ~ 1, server = "peer_a", data = "study", datasources = conns)

  expect_s3_class(direct, "dsvert_dp_frequency_nb2")
  expect_identical(alias$theta, direct$theta)
  expect_identical(calls, rep(list(list(
    data_name = "study", variable = "admissions", server = "peer_a",
    datasources = conns)), 2L))
  expect_error(ds.vertNBFullRegTheta(
    admissions ~ x, data = "study", server = "peer_a", datasources = conns),
    "intercept-only")
  expect_length(calls, 2L)
  expect_error(ds.vertNBFullRegTheta(
    admissions ~ 1, data = "study", datasources = conns),
    "requires an explicit source owner")
  expect_length(calls, 2L)
  expect_error(ds.vertNBFullRegTheta(
    admissions ~ 1, frequency = frequency, server = "peer_a"),
    "legacy controls")
  expect_length(calls, 2L)
})

test_that("frequency-backed NB2 fails closed before validation for unsupported inputs", {
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) {
      calls <<- calls + 1L
      x
    },
    .package = "dsVertClient")

  frequency <- .formal_nb2_frequency()
  expect_error(ds.vertNBFullRegTheta(admissions ~ x, frequency = frequency),
               "intercept-only")
  expect_error(ds.vertNBFullRegTheta(outcome ~ 1, frequency = frequency),
               "outcome")
  expect_error(ds.vertNBFullRegTheta(admissions ~ 1, max_iter = 2,
                                     frequency = frequency),
               "legacy controls")
  expect_identical(calls, 1L)
})

test_that("frequency-backed NB2 rejects non-count support and aliases safely", {
  frequency <- .formal_nb2_frequency()
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_contract = function(x) x,
    .package = "dsVertClient")

  bad <- frequency
  bad$levels <- c("0", "one", "2", "3")
  names(bad$counts) <- bad$levels
  expect_error(ds.vertNBFullRegTheta(admissions ~ 1, frequency = bad),
               "non-negative integer")

  poisson <- frequency
  poisson$counts <- c(`0` = 10, `1` = 20, `2` = 10, `3` = 0)
  fit <- ds.vert.nb(admissions ~ 1, frequency = poisson)
  expect_s3_class(fit, "dsvert_dp_frequency_nb2")
  expect_true(is.infinite(fit$theta))
  expect_identical(fit$dispersion_status, "poisson_limit")
  expect_false(fit$production_ready)
  expect_error(ds.vert.nb(admissions ~ 1, frequency = frequency,
                          datasources = list()),
               "does not accept datasources")
})
