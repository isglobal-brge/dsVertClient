# Tests for print/summary methods with mock objects

# =============================================================================
# print.ds.cor
# =============================================================================

test_that("print.ds.cor displays correctly", {
  mock_cor <- list(
    correlation = matrix(c(1, 0.5, 0.5, 1), 2, 2,
                         dimnames = list(c("x", "y"), c("x", "y"))),
    var_names = c("x", "y"),
    n_obs = 100,
    method = "MHE-CKKS-Threshold",
    servers = c("server1", "server2")
  )
  class(mock_cor) <- c("ds.cor", "list")

  output <- capture.output(print(mock_cor))
  expect_true(any(grepl("MHE-CKKS-Threshold", output)))
  expect_true(any(grepl("100", output)))
  expect_true(any(grepl("server1", output)))
})

test_that("print.ds.cor respects digits parameter", {
  mock_cor <- list(
    correlation = matrix(c(1, 0.123456789, 0.123456789, 1), 2, 2,
                         dimnames = list(c("x", "y"), c("x", "y"))),
    var_names = c("x", "y"),
    n_obs = 50,
    method = "MHE-CKKS-Threshold",
    servers = c("s1", "s2")
  )
  class(mock_cor) <- c("ds.cor", "list")

  output2 <- capture.output(print(mock_cor, digits = 2))
  expect_true(any(grepl("0.12", output2)))
})

# =============================================================================
# print.ds.pca
# =============================================================================

test_that("print.ds.pca displays correctly", {
  mock_pca <- list(
    loadings = matrix(c(0.7, 0.3, -0.3, 0.7), 2, 2,
                      dimnames = list(c("x", "y"), c("PC1", "PC2"))),
    eigenvalues = c(PC1 = 1.5, PC2 = 0.5),
    variance_pct = c(PC1 = 75, PC2 = 25),
    cumulative_pct = c(PC1 = 75, PC2 = 100),
    var_names = c("x", "y"),
    n_obs = 100,
    correlation = matrix(c(1, 0.5, 0.5, 1), 2, 2)
  )
  class(mock_pca) <- c("ds.pca", "list")

  output <- capture.output(print(mock_pca))
  expect_true(any(grepl("Principal Component", output)))
  expect_true(any(grepl("100", output)))
  expect_true(any(grepl("PC1", output)))
})

# =============================================================================
# print.ds.glm / summary.ds.glm / coef.ds.glm
# =============================================================================

test_that("print.ds.glm displays correctly", {
  mock_glm <- list(
    coefficients = c("(Intercept)" = 0.5, "x" = 1.2),
    family = "gaussian",
    n_obs = 100,
    n_iter = 5,
    converged = TRUE,
    deviance = 50,
    formula = "y ~ x"
  )
  class(mock_glm) <- c("ds.glm", "list")

  output <- capture.output(print(mock_glm))
  expect_true(any(grepl("gaussian", output)))
  expect_true(any(grepl("converged", output, ignore.case = TRUE)))
})

test_that("summary.ds.glm returns object invisibly", {
  mock_glm <- list(
    coefficients = c("(Intercept)" = 0.5, "x" = 1.2),
    family = "gaussian",
    n_obs = 100,
    n_iter = 5,
    converged = TRUE,
    deviance = 50,
    formula = "y ~ x"
  )
  class(mock_glm) <- c("ds.glm", "list")

  output <- capture.output(result <- summary(mock_glm))
  expect_s3_class(result, "ds.glm")
})

test_that("coef.ds.glm extracts coefficients", {
  mock_glm <- list(
    coefficients = c("(Intercept)" = 0.5, "x" = 1.2),
    family = "gaussian",
    n_obs = 100,
    n_iter = 5,
    converged = TRUE,
    deviance = 50,
    formula = "y ~ x"
  )
  class(mock_glm) <- c("ds.glm", "list")

  coefs <- coef(mock_glm)
  expect_equal(coefs[["(Intercept)"]], 0.5)
  expect_equal(coefs[["x"]], 1.2)
})
