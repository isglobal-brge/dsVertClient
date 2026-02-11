# Integration tests for complete workflows

skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")

library(DSLite)
library(dsVert)

# Skip if MHE binary not available
skip_if_not(dsVert::mheAvailable(), "MHE tool binary not available")

# NOTE: DSLite has limitations with very large function arguments (~1MB base64 keys)
# Full MHE integration tests require actual DataSHIELD servers.
# Tests that involve ds.vertCor/ds.vertPCA (which use MHE) are skipped for DSLite.
# GLM tests don't use MHE and can still run.

# =============================================================================
# Complete Analysis Workflow
# =============================================================================

test_that("complete workflow: validation -> alignment -> GLM (without MHE)", {
  # NOTE: ds.vertCor and ds.vertPCA use MHE which requires large key transfers
  # DSLite cannot handle these due to argument size limitations
  # This test focuses on the GLM workflow which doesn't use MHE

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

  # Step 4-5: SKIP ds.vertCor and ds.vertPCA (MHE not testable with DSLite)
  # These require actual DataSHIELD servers due to large key parameter sizes

  # Step 6: Fit multiple GLM families (these work without MHE)

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

test_that("workflow works with two servers (alignment only)", {
  # NOTE: ds.vertCor and ds.vertPCA use MHE which can't be tested with DSLite
  # This test verifies the alignment workflow with 2 servers

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

  # SKIP: Correlation and PCA use MHE (DSLite limitation)
  # cor_result <- ds.vertCor("D_aligned", variables, datasources = conns_2)
  # pca_result <- ds.vertPCA("D_aligned", variables, n_components = 2, datasources = conns_2)

  # Verify alignment worked by checking observation count
  obs_count <- DSI::datashield.aggregate(conns_2["server1"],
    as.symbol("getObsCountDS('D_aligned')"))
  expect_equal(obs_count[[1]], env$n)
})

# =============================================================================
# Workflow with Subset of Variables
# =============================================================================

test_that("workflow works with single variable per server (GLM only)", {
  # NOTE: ds.vertCor uses MHE which can't be tested with DSLite
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

  # SKIP: Correlation uses MHE (DSLite limitation)
  # cor_result <- ds.vertCor("D_aligned", variables, datasources = conns)

  # GLM (doesn't use MHE)
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
