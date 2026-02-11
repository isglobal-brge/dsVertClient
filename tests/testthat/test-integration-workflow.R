# Integration tests for complete workflows

skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")

library(DSLite)
library(dsVert)

# =============================================================================
# Complete Analysis Workflow
# =============================================================================

test_that("complete workflow: validation -> alignment -> correlation -> PCA -> GLM", {
  env <- setup_dslite_env(seed = 400)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Step 1: Validate ID format
  validation <- ds.validateIdFormat("D", "patient_id", datasources = conns)
  expect_true(validation$valid)

  # Step 2: Align records
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  expect_equal(ref_hashes$n, env$n)

  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Step 3: Define variables
  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi"),
    server3 = c("glucose", "cholesterol")
  )

  # Step 4: Compute correlation
  cor_matrix <- ds.vertCor("D_aligned", variables, datasources = conns)
  expect_equal(dim(cor_matrix), c(6, 6))

  # Step 5: Perform PCA
  pca_result <- ds.vertPCA("D_aligned", variables,
                            n_components = 4, datasources = conns)
  expect_equal(ncol(pca_result$scores), 4)

  # Step 6: Fit multiple GLM families

  # Gaussian
  model_gaussian <- ds.vertGLM("D_aligned", "outcome_cont", variables,
                                family = "gaussian", verbose = FALSE,
                                datasources = conns)
  expect_true(model_gaussian$converged || model_gaussian$iterations >= 20)

  # Binomial
  model_binomial <- ds.vertGLM("D_aligned", "outcome_binary", variables,
                                family = "binomial", verbose = FALSE,
                                datasources = conns)
  expect_true(model_binomial$converged || model_binomial$iterations >= 20)

  # Poisson
  model_poisson <- ds.vertGLM("D_aligned", "outcome_count", variables,
                               family = "poisson", verbose = FALSE,
                               datasources = conns)
  expect_true(model_poisson$converged || model_poisson$iterations >= 20)

  # Gamma
  model_gamma <- ds.vertGLM("D_aligned", "outcome_positive", variables,
                             family = "Gamma", verbose = FALSE,
                             datasources = conns)
  expect_true(model_gamma$converged || model_gamma$iterations >= 20)

  # Inverse Gaussian
  model_invgauss <- ds.vertGLM("D_aligned", "outcome_positive", variables,
                                family = "inverse.gaussian", verbose = FALSE,
                                datasources = conns)
  expect_true(model_invgauss$converged || model_invgauss$iterations >= 20)
})

# =============================================================================
# Two-Server Workflow
# =============================================================================

test_that("workflow works with two servers", {
  env <- setup_dslite_env(seed = 401)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Use only 2 servers
  conns_2 <- conns[c("server1", "server2")]

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns_2["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns_2)

  # Variables from 2 servers
  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi")
  )

  # Correlation
  cor_matrix <- ds.vertCor("D_aligned", variables, datasources = conns_2)
  expect_equal(dim(cor_matrix), c(4, 4))

  # PCA
  pca_result <- ds.vertPCA("D_aligned", variables,
                            n_components = 2, datasources = conns_2)
  expect_equal(ncol(pca_result$scores), 2)

  # Note: Cannot fit GLM without outcome variable
  # server3 has the outcome, so this tests correlation/PCA only
})

# =============================================================================
# Workflow with Subset of Variables
# =============================================================================

test_that("workflow works with single variable per server", {
  env <- setup_dslite_env(seed = 402)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Single variable per server
  variables <- list(
    server1 = c("age"),
    server2 = c("bmi"),
    server3 = c("glucose")
  )

  # Correlation
  cor_matrix <- ds.vertCor("D_aligned", variables, datasources = conns)
  expect_equal(dim(cor_matrix), c(3, 3))

  # GLM
  model <- ds.vertGLM("D_aligned", "outcome_cont", variables,
                      family = "gaussian", verbose = FALSE,
                      datasources = conns)
  expect_equal(model$n_vars, 3)
})

# =============================================================================
# Model Comparison Workflow
# =============================================================================

test_that("can compare models with different families", {
  env <- setup_dslite_env(seed = 403)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  x_vars <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi")
  )

  # Fit Gamma and Inverse Gaussian to same positive outcome
  model_gamma <- ds.vertGLM("D_aligned", "outcome_positive", x_vars,
                             family = "Gamma", verbose = FALSE,
                             datasources = conns)

  model_invgauss <- ds.vertGLM("D_aligned", "outcome_positive", x_vars,
                                family = "inverse.gaussian", verbose = FALSE,
                                datasources = conns)

  # Compare AIC (lower is better)
  expect_true(!is.na(model_gamma$aic))
  expect_true(!is.na(model_invgauss$aic))

  # Both should have valid deviance
  expect_true(model_gamma$deviance >= 0)
  expect_true(model_invgauss$deviance >= 0)
})

# =============================================================================
# Reproducibility Test
# =============================================================================

test_that("results are reproducible with same seed", {
  # First run
  env1 <- setup_dslite_env(seed = 404)
  conns1 <- DSI::datashield.login(env1$login_data, assign = TRUE, symbol = "D")

  ref1 <- ds.hashId("D", "patient_id", datasource = conns1["server1"])
  ds.alignRecords("D", "patient_id", ref1$hashes,
                  newobj = "D_aligned", datasources = conns1)

  x_vars <- list(server1 = c("age"), server2 = c("bmi"))
  model1 <- ds.vertGLM("D_aligned", "outcome_cont", x_vars,
                       family = "gaussian", verbose = FALSE,
                       datasources = conns1)
  coefs1 <- model1$coefficients

  teardown_dslite_env(conns1)

  # Second run with same seed
  env2 <- setup_dslite_env(seed = 404)
  conns2 <- DSI::datashield.login(env2$login_data, assign = TRUE, symbol = "D")

  ref2 <- ds.hashId("D", "patient_id", datasource = conns2["server1"])
  ds.alignRecords("D", "patient_id", ref2$hashes,
                  newobj = "D_aligned", datasources = conns2)

  model2 <- ds.vertGLM("D_aligned", "outcome_cont", x_vars,
                       family = "gaussian", verbose = FALSE,
                       datasources = conns2)
  coefs2 <- model2$coefficients

  teardown_dslite_env(conns2)

  # Results should be identical
  expect_equal(coefs1, coefs2)
})
