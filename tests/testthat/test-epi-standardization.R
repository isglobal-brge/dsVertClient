test_that("direct standardization matches the Fay-Feuer Poisson-gamma formulas", {
  cases <- c("70-79" = 8, "50-59" = 12, "60-69" = 20)
  person_time <- c("50-59" = 1200, "60-69" = 800, "70-79" = 400)
  standard_population <- c("60-69" = 3000, "70-79" = 1000,
                           "50-59" = 6000)

  fit <- ds.vertDirectStandardization(
    cases, person_time, standard_population, scale = 1e5)

  pt <- person_time[names(cases)]
  std <- standard_population[names(cases)]
  weights <- std / sum(std)
  rate <- sum(weights * cases / pt)
  variance <- sum(weights^2 * cases / pt^2)
  one_event <- max(weights / pt)
  alpha <- 0.05
  expected_lower <- stats::qgamma(
    alpha / 2, shape = rate^2 / variance, scale = variance / rate)
  expected_upper <- stats::qgamma(
    1 - alpha / 2,
    shape = (rate + one_event)^2 / (variance + one_event^2),
    scale = (variance + one_event^2) / (rate + one_event))

  expect_equal(fit$estimate, rate * 1e5)
  expect_equal(fit$std_error, sqrt(variance) * 1e5)
  expect_equal(fit$conf_int, c(lower = expected_lower,
                               upper = expected_upper) * 1e5)
  expect_equal(fit$standard_weights, weights)
  expect_identical(fit$ci_method, "poisson_gamma")
  expect_identical(fit$input_provenance, "caller_supplied_aggregate")
  expect_identical(fit$method_status, "promoted")
  expect_match(fit$correction$name, "one-event")
  expect_true(length(fit$assumptions) >= 3L)
})

test_that("direct standardization offers a bounded normal interval", {
  fit <- ds.vertDirectStandardization(
    cases = c(0, 1), person_time = c(100, 100),
    standard_population = c(1, 1), scale = 1,
    ci_method = "normal")

  estimate <- 0.5 * 0 + 0.5 * 0.01
  std_error <- sqrt(0.5^2 * 0 / 100^2 + 0.5^2 * 1 / 100^2)
  expected_upper <- estimate + stats::qnorm(0.975) * std_error

  expect_equal(fit$estimate, estimate)
  expect_equal(fit$conf_int,
               c(lower = 0, upper = expected_upper))
  expect_true(fit$correction$applied)
  expect_identical(fit$correction$name, "non-negative lower bound")
})

test_that("direct Poisson-gamma interval handles zero events", {
  fit <- ds.vertDirectStandardization(
    cases = c(0, 0), person_time = c(100, 200),
    standard_population = c(1, 1), scale = 1)

  one_event <- max(c(0.5 / 100, 0.5 / 200))
  expect_equal(fit$estimate, 0)
  expect_equal(fit$conf_int[["lower"]], 0)
  expect_equal(fit$conf_int[["upper"]],
               stats::qgamma(0.975, shape = 1, scale = one_event))
})

test_that("indirect standardization matches an exact Garwood interval", {
  observed <- c("50-59" = 4, "60-69" = 7, "70-79" = 9)
  expected <- c("70-79" = 5.5, "50-59" = 3.5, "60-69" = 6)

  fit <- ds.vertIndirectStandardization(
    observed, expected, measure = "SMR", scale = 100)

  total_observed <- sum(observed)
  total_expected <- sum(expected)
  expected_ci <- c(
    lower = 0.5 * stats::qchisq(0.025, 2 * total_observed) /
      total_expected,
    upper = 0.5 * stats::qchisq(0.975, 2 * (total_observed + 1)) /
      total_expected) * 100

  expect_equal(fit$ratio, total_observed / total_expected)
  expect_equal(fit$estimate, 100 * total_observed / total_expected)
  expect_equal(fit$conf_int, expected_ci)
  expect_identical(fit$measure, "SMR")
  expect_identical(fit$ci_method, "exact_poisson_garwood")
  expect_identical(fit$input_provenance, "caller_supplied_aggregate")
  expect_identical(fit$method_status, "promoted")
  expect_identical(fit$correction$applied, FALSE)
  expect_true(length(fit$assumptions) >= 2L)
})

test_that("indirect standardization handles zero events and SIR labeling", {
  fit <- ds.vertIndirectStandardization(
    observed = c(0, 0), expected = c(1.5, 2.5), measure = "SIR",
    scale = 1, level = 0.90)

  expect_identical(fit$measure, "SIR")
  expect_equal(fit$estimate, 0)
  expect_equal(fit$conf_int[["lower"]], 0)
  expect_equal(fit$conf_int[["upper"]],
               0.5 * stats::qchisq(0.95, 2) / 4)
})

test_that("standardization methods fail closed on malformed aggregates", {
  expect_error(
    ds.vertDirectStandardization(c(1, -1), c(10, 10), c(1, 1)),
    "non-negative")
  expect_error(
    ds.vertDirectStandardization(c(1, 1.2), c(10, 10), c(1, 1)),
    "whole-number")
  expect_error(
    ds.vertDirectStandardization(c(a = 1, b = 2),
                                 c(a = 10, c = 20),
                                 c(a = 1, b = 1)),
    "same stratum names")
  expect_error(
    ds.vertDirectStandardization(c(1, 2), c(10, 0), c(1, 1)),
    "strictly positive")
  expect_error(
    ds.vertDirectStandardization(c(1, 2), c(10, 10), c(0, 0)),
    "positive total")
  expect_error(
    ds.vertDirectStandardization(c(1, 2), c(10, 10),
                                 c(.Machine$double.xmax,
                                   .Machine$double.xmax)),
    "finite positive total")
  expect_error(
    ds.vertDirectStandardization(c(1, 2), c(10, 10), c(1, 1), scale = 0),
    "strictly positive")
  expect_error(
    ds.vertDirectStandardization(2, 1, 1, scale = .Machine$double.xmax),
    "reported scale")

  expect_error(
    ds.vertIndirectStandardization(c(1, 2), c(1), measure = "SMR"),
    "same length")
  expect_error(
    ds.vertIndirectStandardization(c(1, 2), c(0, 0), measure = "SMR"),
    "positive total")
  expect_error(
    ds.vertIndirectStandardization(c(1, NA), c(1, 2), measure = "SMR"),
    "finite")
  expect_error(
    ds.vertIndirectStandardization(c(1, 2),
                                   c(.Machine$double.xmax,
                                     .Machine$double.xmax)),
    "finite positive total")
  expect_error(
    ds.vertIndirectStandardization(c(1, 2), c(1, 2), level = 1),
    "in \\(0, 1\\)")
  expect_error(
    ds.vertIndirectStandardization(2, 1, scale = .Machine$double.xmax),
    "reported scale")
})

test_that("standardization print methods return their objects invisibly", {
  direct <- ds.vertDirectStandardization(5, 100, 1, scale = 1000)
  indirect <- ds.vertIndirectStandardization(5, 4, measure = "SIR")

  expect_output(returned_direct <- print(direct),
                "Directly standardized rate")
  expect_output(returned_indirect <- print(indirect),
                "SIR")
  expect_identical(returned_direct, direct)
  expect_identical(returned_indirect, indirect)
})
