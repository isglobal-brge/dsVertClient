# Integration tests for correlation and PCA functions

skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")

library(DSLite)
library(dsVert)

# =============================================================================
# Test Correlation
# =============================================================================

test_that("ds.vertCor computes correlation matrix across servers", {
  env <- setup_dslite_env(seed = 200)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align data first
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Define variables per server
  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi"),
    server3 = c("glucose", "cholesterol")
  )

  # Compute correlation
  cor_matrix <- ds.vertCor("D_aligned", variables, datasources = conns)

  # Check structure
  expect_true(is.matrix(cor_matrix))
  expect_equal(nrow(cor_matrix), 6)
  expect_equal(ncol(cor_matrix), 6)

  # Check diagonal is 1
  expect_equal(unname(diag(cor_matrix)), rep(1, 6), tolerance = 1e-6)

  # Check symmetry
  expect_equal(cor_matrix, t(cor_matrix), tolerance = 1e-6)

  # Check values are in [-1, 1]
  expect_true(all(cor_matrix >= -1 & cor_matrix <= 1))
})

test_that("ds.vertCor works with two servers", {
  env <- setup_dslite_env(seed = 201)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Use only 2 servers
  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi")
  )

  cor_matrix <- ds.vertCor("D_aligned", variables, datasources = conns)

  expect_equal(nrow(cor_matrix), 4)
  expect_equal(ncol(cor_matrix), 4)
})

# =============================================================================
# Test PCA
# =============================================================================

test_that("ds.vertPCA performs PCA across servers", {
  env <- setup_dslite_env(seed = 202)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi"),
    server3 = c("glucose", "cholesterol")
  )

  # Perform PCA
  pca_result <- ds.vertPCA("D_aligned", variables,
                            n_components = 3, datasources = conns)

  # Check structure
  expect_s3_class(pca_result, "ds.pca")
  expect_true("scores" %in% names(pca_result))
  expect_true("loadings" %in% names(pca_result))
  expect_true("variance" %in% names(pca_result))
  expect_true("variance_pct" %in% names(pca_result))
  expect_true("cumulative_pct" %in% names(pca_result))

  # Check dimensions
  expect_equal(nrow(pca_result$scores), env$n)
  expect_equal(ncol(pca_result$scores), 3)
  expect_equal(nrow(pca_result$loadings), 6)  # 6 variables
  expect_equal(ncol(pca_result$loadings), 3)

  # Check variance explained sums to ~100%
  expect_equal(unname(pca_result$cumulative_pct[3]),
               unname(sum(pca_result$variance_pct[1:3])),
               tolerance = 1e-6)

  # Variance should be decreasing
  expect_true(all(diff(pca_result$variance) <= 0))
})

test_that("ds.vertPCA respects n_components parameter", {
  env <- setup_dslite_env(seed = 203)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi")
  )

  # Request 2 components
  pca2 <- ds.vertPCA("D_aligned", variables, n_components = 2, datasources = conns)
  expect_equal(ncol(pca2$scores), 2)

  # Request all components (4 variables = max 4 components)
  pca4 <- ds.vertPCA("D_aligned", variables, n_components = 4, datasources = conns)
  expect_equal(ncol(pca4$scores), 4)

  # Request more than available (should cap at max)
  pca_max <- ds.vertPCA("D_aligned", variables, n_components = 10, datasources = conns)
  expect_equal(ncol(pca_max$scores), 4)
})

test_that("ds.vertPCA print method works", {
  env <- setup_dslite_env(seed = 204)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi")
  )

  pca_result <- ds.vertPCA("D_aligned", variables, n_components = 2, datasources = conns)

  # Print should not error
  expect_output(print(pca_result), "Principal Component Analysis")
})
