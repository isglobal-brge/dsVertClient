# Tests for ds.vertContrast.

library(testthat)

make_fit_with_cov <- function() {
  coefs <- c(`(Intercept)` = 1.0, x1 = 0.5, x2 = -0.3, x3 = 0.2)
  cov <- diag(c(0.04, 0.01, 0.01, 0.01))
  dimnames(cov) <- list(names(coefs), names(coefs))
  fit <- list(
    coefficients = coefs,
    std_errors = sqrt(diag(cov)),
    covariance = cov,
    family = "gaussian", n_obs = 200L, n_vars = 4L,
    deviance = 100, null_deviance = 200)
  class(fit) <- c("ds.glm", "list")
  fit
}

test_that("contrast matches single-coef Wald when K is one-hot row", {
  fit <- make_fit_with_cov()
  K <- matrix(c(0, 1, 0, 0), nrow = 1)  # beta_x1
  res <- ds.vertContrast(fit, K = K)

  w <- ds.vertWald(fit, parm = "x1")
  # chi^2 on 1 df equals z^2
  expect_equal(res$statistic, w$z^2, tolerance = 1e-10)
  expect_equal(res$df, 1L)
  expect_equal(res$p_value, w$p_value, tolerance = 1e-10)
})

test_that("contrast accepts character vector of coefs (conjunction)", {
  fit <- make_fit_with_cov()
  res <- ds.vertContrast(fit, K = c("x1", "x2"))
  expect_equal(res$df, 2L)
  # Estimate should be (x1, x2) coefs
  expect_equal(unname(res$estimate), c(0.5, -0.3), tolerance = 1e-10)
})

test_that("contrast supports non-zero null vector", {
  fit <- make_fit_with_cov()
  K <- rbind(c(0, 1, 0, 0), c(0, 0, 1, 0))
  m <- c(0.4, -0.5)  # H0: x1=0.4, x2=-0.5
  res <- ds.vertContrast(fit, K = K, m = m)
  expect_equal(unname(res$estimate), c(0.5 - 0.4, -0.3 - (-0.5)),
               tolerance = 1e-10)
  # With small SE = 0.1, both estimates (0.1 and 0.2) -> W = (0.1/0.1)^2 + (0.2/0.1)^2 = 1+4=5
  expect_equal(res$statistic, 5.0, tolerance = 1e-8)
})

test_that("contrast errors when covariance is missing", {
  fit <- make_fit_with_cov()
  fit$covariance <- NULL
  expect_error(ds.vertContrast(fit, K = "x1"),
               "does not expose the full covariance matrix")
})

test_that("contrast errors on singular K*Cov*Kt", {
  fit <- make_fit_with_cov()
  # Perfectly redundant K rows
  K <- rbind(c(0, 1, 0, 0), c(0, 2, 0, 0))  # linearly dependent
  expect_error(ds.vertContrast(fit, K = K), "singular")
})

test_that("contrast matches lmtest::waldtest-style joint on linear model", {
  # Build a known covariance and verify the joint Wald stat
  coefs <- c(a = 1, b = 2, c = 3)
  cov <- matrix(c(
    1, 0.5, 0,
    0.5, 2, 0,
    0, 0, 0.25), nrow = 3, byrow = TRUE)
  dimnames(cov) <- list(names(coefs), names(coefs))
  fit <- list(coefficients = coefs, std_errors = sqrt(diag(cov)),
              covariance = cov, family = "gaussian", n_obs = 100L,
              n_vars = 3L)
  class(fit) <- c("ds.glm", "list")

  # Test joint hypothesis b = 0, c = 0
  K <- rbind(c(0, 1, 0), c(0, 0, 1))
  res <- ds.vertContrast(fit, K = K)
  # Manually: V = K cov K^T = [[2,0],[0,0.25]]; V^-1 = [[0.5,0],[0,4]]
  # W = [2, 3] V^-1 [2, 3]^T = 2*0.5*2 + 3*4*3 = 2 + 36 = 38
  expect_equal(res$statistic, 38, tolerance = 1e-10)
  expect_equal(res$df, 2L)
})
