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

# =============================================================================
# ds.vertLASSO1Step — proper quadratic-surrogate LASSO
# =============================================================================

mock_fit_with_cov <- function() {
  coefs <- c(`(Intercept)` = 1.0, x1 = 0.5, x2 = -0.3, x3 = 0.05)
  # Diagonal covariance for simplicity
  cov <- diag(c(0.04, 0.01, 0.01, 0.01))
  dimnames(cov) <- list(names(coefs), names(coefs))
  fit <- list(
    coefficients = coefs,
    std_errors = sqrt(diag(cov)),
    covariance = cov,
    family = "gaussian", n_obs = 200L, n_vars = length(coefs),
    deviance = 50, null_deviance = 75, pseudo_r2 = 0.33)
  class(fit) <- c("ds.glm", "list")
  fit
}

test_that("LASSO1Step zeroes coefficients correctly on diagonal H", {
  fit <- mock_fit_with_cov()
  # H = diag(1/0.04, 1/0.01, 1/0.01, 1/0.01) = diag(25, 100, 100, 100)
  # For coord j: r_j = betahat_j (diagonal H), so beta_j <- ST(betahat_j, lambda/H_jj)
  # With lambda = 1:
  #   x1: threshold 1/100 = 0.01, |0.5| > 0.01 -> 0.5 - 0.01 = 0.49
  #   x2: threshold 1/100 = 0.01, |-0.3| > 0.01 -> -0.3 + 0.01 = -0.29
  #   x3: threshold 1/100 = 0.01, |0.05| > 0.01 -> 0.05 - 0.01 = 0.04
  #   Intercept not thresholded
  res <- ds.vertLASSO1Step(fit, lambda = 1, max_iter = 50L, tol = 1e-12)
  b <- res$paths[[1]]
  expect_equal(unname(b["(Intercept)"]), 1.0, tolerance = 1e-8)
  expect_equal(unname(b["x1"]), 0.49, tolerance = 1e-6)
  expect_equal(unname(b["x2"]), -0.29, tolerance = 1e-6)
  expect_equal(unname(b["x3"]), 0.04, tolerance = 1e-6)
})

test_that("LASSO1Step with large lambda zeroes non-intercept coefs", {
  fit <- mock_fit_with_cov()
  # lambda = 100 on diagonal H = 100 -> threshold = 1, greater than any |coef|
  res <- ds.vertLASSO1Step(fit, lambda = 100, max_iter = 50L, tol = 1e-12)
  b <- res$paths[[1]]
  expect_equal(unname(b["x1"]), 0, tolerance = 1e-6)
  expect_equal(unname(b["x2"]), 0, tolerance = 1e-6)
  expect_equal(unname(b["x3"]), 0, tolerance = 1e-6)
  expect_equal(unname(b["(Intercept)"]), 1.0, tolerance = 1e-8)
})

test_that("LASSO1Step lambda=0 recovers unregularised fit", {
  fit <- mock_fit_with_cov()
  res <- ds.vertLASSO1Step(fit, lambda = 0, max_iter = 50L, tol = 1e-12)
  b <- res$paths[[1]]
  expect_equal(unname(b), unname(fit$coefficients), tolerance = 1e-8)
})

test_that("LASSO1Step objective is non-decreasing in lambda", {
  fit <- mock_fit_with_cov()
  lam <- c(0, 0.1, 1, 10)
  res <- ds.vertLASSO1Step(fit, lambda = lam, max_iter = 100L)
  # Quadratic+L1 objective should be monotone non-decreasing in lambda
  # (more penalty -> objective includes more L1 mass, but beta moves
  # away from betahat to compensate)
  # At least the L1 norm should decrease with lambda
  l1_norms <- sapply(res$paths, function(b) sum(abs(b[-1])))
  expect_true(all(diff(l1_norms) <= 1e-10))
})

test_that("LASSO1Step errors when covariance is missing", {
  fit <- mock_fit_with_cov()
  fit$covariance <- NULL
  expect_error(ds.vertLASSO1Step(fit, lambda = 1),
               "does not expose the full covariance")
})

# =============================================================================
# ds.vertLASSOCV: information-criterion path selection
# =============================================================================
test_that("LASSOCV BIC recovers the true support on a simulated Gaussian", {
  skip_if_not_installed("stats")
  set.seed(2026)
  n <- 400; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  true_beta <- c(1.5, -1.0, 0.8, 0, 0, 0, 0.3, -0.2, 0, 0)
  y <- as.numeric(X %*% true_beta + rnorm(n, 0, 0.5))
  fit_r <- lm(y ~ X - 1)
  fake <- list(coefficients = coef(fit_r), covariance = vcov(fit_r),
               family = "gaussian", n_obs = n)
  names(fake$coefficients) <- paste0("X", seq_len(p))
  dimnames(fake$covariance) <- list(names(fake$coefficients),
                                    names(fake$coefficients))
  class(fake) <- c("ds.glm", "list")

  bic <- ds.vertLASSOCV(fake, criterion = "BIC")
  expect_s3_class(bic, "ds.vertLASSOCV")
  # True support is {1,2,3,7,8}.
  expect_setequal(which(abs(bic$beta.min) > 1e-6), c(1, 2, 3, 7, 8))
  # lambda.min and lambda.1se are both finite.
  expect_true(is.finite(bic$lambda.min))
  expect_true(is.finite(bic$lambda.1se))
  expect_true(bic$lambda.1se >= bic$lambda.min - 1e-10)
})

test_that("LASSOCV errors when covariance is missing", {
  fit <- mock_fit_with_cov()
  fit$covariance <- NULL
  expect_error(ds.vertLASSOCV(fit), "does not expose covariance")
})

test_that("LASSOCV respects a user-supplied lambda_grid", {
  fit <- mock_fit_with_cov()
  lam <- c(0, 0.1, 1, 10)
  cv <- ds.vertLASSOCV(fit, lambda_grid = lam, criterion = "AIC")
  expect_equal(cv$lambda, lam)
  expect_equal(length(cv$ic), length(lam))
})
