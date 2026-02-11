# Integration tests for GLM functions

skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")

library(DSLite)
library(dsVert)

# =============================================================================
# Test Gaussian GLM (Linear Regression)
# =============================================================================

test_that("ds.vertGLM fits gaussian family across servers", {
  env <- setup_dslite_env(seed = 300)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Variables from different servers
  x_vars <- list(
    server1 = c("age", "weight"),
    server3 = c("glucose")
  )

  # Fit Gaussian GLM with loose tolerance for test speed
  model <- ds.vertGLM("D_aligned", "outcome_cont", x_vars,
                      family = "gaussian",
                      tol = 1e-4,
                      verbose = FALSE,
                      datasources = conns)

  # Check structure
  expect_s3_class(model, "ds.glm")
  expect_equal(model$family, "gaussian")
  expect_true(model$converged || model$iterations >= 50)  # Either converged or ran many iterations
  expect_equal(model$n_obs, env$n)
  expect_equal(model$n_vars, 3)
  expect_length(model$coefficients, 3)

  # Check deviance metrics
expect_true(!is.null(model$deviance))
  expect_true(!is.null(model$null_deviance))
  expect_true(!is.null(model$pseudo_r2))
  expect_true(!is.null(model$aic))
  expect_true(model$deviance >= 0)
  expect_true(model$null_deviance >= 0)
})

# =============================================================================
# Test Binomial GLM (Logistic Regression)
# =============================================================================

test_that("ds.vertGLM fits binomial family across servers", {
  env <- setup_dslite_env(seed = 301)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(
    server1 = c("age"),
    server2 = c("bmi"),
    server3 = c("glucose")
  )

  # Fit Binomial GLM
  model <- ds.vertGLM("D_aligned", "outcome_binary", x_vars,
                      family = "binomial",
                      verbose = FALSE,
                      datasources = conns)

  expect_s3_class(model, "ds.glm")
  expect_equal(model$family, "binomial")
  expect_true(model$converged || model$iterations >= 20)
  expect_equal(model$n_vars, 3)

  # Deviance should be positive
  expect_true(model$deviance > 0)
  expect_true(model$null_deviance > 0)
  # Pseudo R2 should be between 0 and 1 for reasonable models
  expect_true(model$pseudo_r2 >= -0.5 && model$pseudo_r2 <= 1)
})

# =============================================================================
# Test Poisson GLM
# =============================================================================

test_that("ds.vertGLM fits poisson family across servers", {
  env <- setup_dslite_env(seed = 302)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(
    server1 = c("age", "weight"),
    server2 = c("height")
  )

  # Fit Poisson GLM
  model <- ds.vertGLM("D_aligned", "outcome_count", x_vars,
                      family = "poisson",
                      verbose = FALSE,
                      datasources = conns)

  expect_s3_class(model, "ds.glm")
  expect_equal(model$family, "poisson")
  expect_true(model$converged || model$iterations >= 20)
  expect_true(model$deviance >= 0)
})

# =============================================================================
# Test Gamma GLM
# =============================================================================

test_that("ds.vertGLM fits Gamma family across servers", {
  env <- setup_dslite_env(seed = 303)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(
    server1 = c("age"),
    server2 = c("bmi")
  )

  # Fit Gamma GLM (for positive continuous data)
  model <- ds.vertGLM("D_aligned", "outcome_positive", x_vars,
                      family = "Gamma",
                      verbose = FALSE,
                      datasources = conns)

  expect_s3_class(model, "ds.glm")
  expect_equal(model$family, "Gamma")
  expect_true(model$converged || model$iterations >= 20)
  expect_equal(model$n_vars, 2)
  expect_true(model$deviance >= 0)
  expect_true(model$null_deviance >= 0)
})

# =============================================================================
# Test Inverse Gaussian GLM
# =============================================================================

test_that("ds.vertGLM fits inverse.gaussian family across servers", {
  env <- setup_dslite_env(seed = 304)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(
    server1 = c("weight"),
    server2 = c("height")
  )

  # Fit Inverse Gaussian GLM
  model <- ds.vertGLM("D_aligned", "outcome_positive", x_vars,
                      family = "inverse.gaussian",
                      verbose = FALSE,
                      datasources = conns)

  expect_s3_class(model, "ds.glm")
  expect_equal(model$family, "inverse.gaussian")
  expect_true(model$converged || model$iterations >= 20)
  expect_true(model$deviance >= 0)
})

# =============================================================================
# Test GLM with all three servers
# =============================================================================

test_that("ds.vertGLM works with three partitions", {
  env <- setup_dslite_env(seed = 305)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Variables from ALL three servers
  x_vars <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi"),
    server3 = c("glucose", "cholesterol")
  )

  model <- ds.vertGLM("D_aligned", "outcome_cont", x_vars,
                      family = "gaussian",
                      verbose = FALSE,
                      datasources = conns)

  expect_equal(model$n_vars, 6)
  expect_length(model$coefficients, 6)
  expect_named(model$coefficients,
               c("age", "weight", "height", "bmi", "glucose", "cholesterol"))
})

# =============================================================================
# Test Model Output Methods
# =============================================================================

test_that("print.ds.glm works correctly", {
  env <- setup_dslite_env(seed = 306)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(server1 = c("age"), server2 = c("bmi"))

  model <- ds.vertGLM("D_aligned", "outcome_cont", x_vars,
                      family = "gaussian", verbose = FALSE,
                      datasources = conns)

  expect_output(print(model), "Vertically Partitioned GLM")
  expect_output(print(model), "gaussian")
})

test_that("summary.ds.glm shows deviance metrics", {
  env <- setup_dslite_env(seed = 307)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(server1 = c("age"), server2 = c("bmi"))

  model <- ds.vertGLM("D_aligned", "outcome_cont", x_vars,
                      family = "gaussian", verbose = FALSE,
                      datasources = conns)

  # Summary should include deviance info
  expect_output(summary(model), "Deviance")
  expect_output(summary(model), "Null deviance")
  expect_output(summary(model), "Residual deviance")
  expect_output(summary(model), "Pseudo R-squared")
  expect_output(summary(model), "AIC")
})

test_that("coef.ds.glm extracts coefficients", {
  env <- setup_dslite_env(seed = 308)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(server1 = c("age", "weight"))

  model <- ds.vertGLM("D_aligned", "outcome_cont", x_vars,
                      family = "gaussian", verbose = FALSE,
                      datasources = conns)

  coefs <- coef(model)
  expect_named(coefs, c("age", "weight"))
  expect_length(coefs, 2)
})

# =============================================================================
# Test GLM Convergence Parameters
# =============================================================================

test_that("ds.vertGLM respects max_iter and tol parameters", {
  env <- setup_dslite_env(seed = 309)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(server1 = c("age"))

  # With very few iterations, might not converge
  model_few_iter <- suppressWarnings(
    ds.vertGLM("D_aligned", "outcome_cont", x_vars,
               family = "gaussian",
               max_iter = 2,
               verbose = FALSE,
               datasources = conns)
  )

  expect_true(model_few_iter$iterations <= 2)

  # With loose tolerance, should converge faster
  model_loose <- ds.vertGLM("D_aligned", "outcome_cont", x_vars,
                            family = "gaussian",
                            tol = 0.1,
                            verbose = FALSE,
                            datasources = conns)

  expect_true(model_loose$converged)
})

# =============================================================================
# Test Error Handling
# =============================================================================

test_that("ds.vertGLM validates inputs", {
  env <- setup_dslite_env(seed = 310)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Invalid family
  expect_error(
    ds.vertGLM("D_aligned", "outcome_cont",
               list(server1 = c("age")),
               family = "invalid_family",
               datasources = conns),
    "family must be"
  )

  # x_vars not a list
  expect_error(
    ds.vertGLM("D_aligned", "outcome_cont",
               c("age", "weight"),
               datasources = conns),
    "x_vars must be a named list"
  )

  # Unknown server in x_vars
  expect_error(
    ds.vertGLM("D_aligned", "outcome_cont",
               list(unknown_server = c("age")),
               datasources = conns),
    "Unknown server"
  )
})
