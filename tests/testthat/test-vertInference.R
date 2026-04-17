# Tests for ds.vertLR, ds.vertConfint, ds.vertWald.
# These helpers are client-side transformations over a ds.glm object;
# tests construct a minimal fake ds.glm and verify the statistical logic.

library(testthat)

make_fake_fit <- function(coefs, ses, deviance, family = "gaussian",
                          n_obs = 100L) {
  fit <- list(
    coefficients = coefs,
    std_errors = ses,
    z_values = coefs / ses,
    p_values = 2 * pnorm(abs(coefs / ses), lower.tail = FALSE),
    iterations = 5L,
    converged = TRUE,
    family = family,
    n_obs = n_obs,
    n_vars = length(coefs),
    lambda = 1e-4,
    deviance = deviance,
    null_deviance = deviance * 1.5,
    pseudo_r2 = 1 - deviance / (deviance * 1.5))
  class(fit) <- c("ds.glm", "list")
  fit
}

# =============================================================================
# ds.vertLR
# =============================================================================

test_that("LR test matches chi-square on simulated nested fits", {
  reduced <- make_fake_fit(
    coefs = c(`(Intercept)` = 1.0, x1 = 0.5),
    ses = c(`(Intercept)` = 0.1, x1 = 0.1),
    deviance = 120.0)
  full <- make_fake_fit(
    coefs = c(`(Intercept)` = 1.0, x1 = 0.5, x2 = 0.3),
    ses = c(`(Intercept)` = 0.1, x1 = 0.1, x2 = 0.1),
    deviance = 100.0)

  lr <- ds.vertLR(reduced, full)
  expect_equal(lr$statistic, 20.0)
  expect_equal(lr$df, 1L)
  expect_equal(lr$p_value, pchisq(20.0, df = 1L, lower.tail = FALSE))
})

test_that("LR test floors slightly-negative statistic at 0", {
  reduced <- make_fake_fit(
    coefs = c(`(Intercept)` = 1.0, x1 = 0.5),
    ses = c(`(Intercept)` = 0.1, x1 = 0.1),
    deviance = 99.9999)
  full <- make_fake_fit(
    coefs = c(`(Intercept)` = 1.0, x1 = 0.5, x2 = 0.3),
    ses = c(`(Intercept)` = 0.1, x1 = 0.1, x2 = 0.1),
    deviance = 100.0)  # slightly higher than reduced due to numerical noise

  lr <- ds.vertLR(reduced, full)
  expect_equal(lr$statistic, 0)
  expect_equal(lr$p_value, 1)
})

test_that("LR test errors on non-nested or mismatched inputs", {
  r <- make_fake_fit(c(`(Intercept)` = 1), c(`(Intercept)` = 0.1), 100)
  f <- make_fake_fit(c(x1 = 0.5), c(x1 = 0.1), 80)
  expect_error(ds.vertLR(r, f), "subset of `full`")

  r2 <- make_fake_fit(c(`(Intercept)` = 1), c(`(Intercept)` = 0.1), 100,
                      family = "gaussian", n_obs = 50L)
  f2 <- make_fake_fit(c(`(Intercept)` = 1, x = 0.5),
                      c(`(Intercept)` = 0.1, x = 0.1),
                      80, family = "gaussian", n_obs = 100L)
  expect_error(ds.vertLR(r2, f2), "same cohort")

  r3 <- make_fake_fit(c(`(Intercept)` = 1), c(`(Intercept)` = 0.1), 100,
                      family = "gaussian")
  f3 <- make_fake_fit(c(`(Intercept)` = 1, x = 0.5),
                      c(`(Intercept)` = 0.1, x = 0.1),
                      80, family = "binomial")
  expect_error(ds.vertLR(r3, f3), "same family")
})

# =============================================================================
# ds.vertConfint
# =============================================================================

test_that("confint produces symmetric Wald CIs around the estimate", {
  fit <- make_fake_fit(
    coefs = c(`(Intercept)` = 2.0, x = 0.5),
    ses = c(`(Intercept)` = 0.5, x = 0.1),
    deviance = 50)
  ci <- ds.vertConfint(fit, level = 0.95)
  z <- qnorm(0.975)
  expect_equal(ci$lower[1], 2.0 - z * 0.5)
  expect_equal(ci$upper[1], 2.0 + z * 0.5)
  expect_equal(ci$lower[2], 0.5 - z * 0.1)
  expect_equal(ci$upper[2], 0.5 + z * 0.1)
})

test_that("confint level affects interval width", {
  fit <- make_fake_fit(c(`(Intercept)` = 0), c(`(Intercept)` = 1), 50)
  ci90 <- ds.vertConfint(fit, level = 0.90)
  ci99 <- ds.vertConfint(fit, level = 0.99)
  expect_true((ci99$upper[1] - ci99$lower[1]) >
              (ci90$upper[1] - ci90$lower[1]))
})

test_that("confint errors on unknown parm", {
  fit <- make_fake_fit(c(`(Intercept)` = 0), c(`(Intercept)` = 1), 50)
  expect_error(ds.vertConfint(fit, parm = "nope"),
               "Unknown coefficient")
})

# =============================================================================
# ds.vertWald
# =============================================================================

test_that("wald reduces to z-test against null = 0", {
  fit <- make_fake_fit(c(`(Intercept)` = 1.5, x = 0.5),
                      c(`(Intercept)` = 0.5, x = 0.1), 50)
  w <- ds.vertWald(fit, parm = "x")
  expect_equal(w$estimate, 0.5)
  expect_equal(w$std_error, 0.1)
  expect_equal(w$z, 5)
  expect_equal(w$p_value, 2 * pnorm(5, lower.tail = FALSE))
})

test_that("wald supports non-zero null", {
  fit <- make_fake_fit(c(x = 0.7), c(x = 0.1), 50)
  w <- ds.vertWald(fit, parm = "x", null = 0.5)
  expect_equal(w$z, (0.7 - 0.5) / 0.1)
})

test_that("wald errors on missing parm", {
  fit <- make_fake_fit(c(x = 0.7), c(x = 0.1), 50)
  expect_error(ds.vertWald(fit, parm = "nope"), "not in fit")
})
