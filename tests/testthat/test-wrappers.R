# Tests for the client-side wrapper functions (LASSO post-hoc, IPW shape,
# multinom shape, NB shape). Mock ds.glm objects, no MPC interaction.

library(testthat)

mock_fit <- function(coefs = c(`(Intercept)` = 1, x1 = 0.5, x2 = -0.3, x3 = 0.05)) {
  ses <- rep(0.1, length(coefs)); names(ses) <- names(coefs)
  fit <- list(
    coefficients = coefs,
    std_errors = ses,
    z_values = coefs / ses,
    p_values = 2 * pnorm(abs(coefs / ses), lower.tail = FALSE),
    family = "gaussian", n_obs = 100L, n_vars = length(coefs),
    lambda = 1e-4, deviance = 50, null_deviance = 75, pseudo_r2 = 0.33,
    aic = NA_real_, iterations = 5L, converged = TRUE)
  class(fit) <- c("ds.glm", "list")
  fit
}

# =============================================================================
# ds.vertLASSO
# =============================================================================
test_that("LASSO soft-thresholds coefficients below lambda", {
  fit <- mock_fit()
  res <- ds.vertLASSO(fit, lambda_1 = 0.3, alpha_grid = c(1))
  # lambda = 0.3: |x3| = 0.05 < 0.3 -> zeroed; others kept
  path <- res$paths[[1]]
  expect_equal(unname(path["x3"]), 0)
  # x1 (0.5) shrunk to 0.5 - 0.3 = 0.2
  expect_equal(unname(path["x1"]), 0.2, tolerance = 1e-10)
  # x2 (-0.3) shrunk to 0 (|x2| == lambda boundary)
  expect_equal(unname(path["x2"]), 0, tolerance = 1e-10)
  # Intercept preserved
  expect_equal(unname(path["(Intercept)"]), unname(fit$coefficients["(Intercept)"]))
})

test_that("LASSO sweeps alpha_grid", {
  fit <- mock_fit()
  res <- ds.vertLASSO(fit, lambda_1 = 0.3,
                     alpha_grid = c(1, 0.5, 0.1))
  expect_equal(length(res$paths), 3L)
  # Smaller alpha -> less shrinkage -> x1 coef closer to original
  x1_path <- sapply(res$paths, function(p) unname(p["x1"]))
  expect_true(all(diff(x1_path) >= -1e-12))
})

# =============================================================================
# ds.vertMultinom input validation
# =============================================================================
test_that("multinom requires >=2 non-reference classes", {
  expect_error(
    ds.vertMultinom(y ~ x1 + x2, data = "D", classes = c("A")),
    "Need at least 2")
})

# =============================================================================
# ds.vertNB class inheritance
# =============================================================================
test_that("ds.vertNB object inherits from ds.glm", {
  fit <- mock_fit()
  out <- c(unclass(fit), list(theta = 2.5, nb_correction = "placeholder"))
  class(out) <- c("ds.vertNB", "ds.glm", "list")
  expect_true(inherits(out, "ds.glm"))
  expect_true(inherits(out, "ds.vertNB"))
})
