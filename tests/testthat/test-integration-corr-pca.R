# Integration tests for correlation and PCA functions (MHE-based)

skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")

library(DSLite)
library(dsVert)

# Skip if MHE binary not available
skip_if_not(dsVert::mheAvailable(), "MHE tool binary not available")

# NOTE: DSLite has limitations with very large function arguments (~1MB base64 keys)
# Full MHE integration tests require actual DataSHIELD servers.
# These tests are skipped for DSLite until a workaround is implemented.
skip("DSLite cannot handle large MHE key parameters - use real DataSHIELD servers")

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

  # Define variables per server (use 2 servers for faster test)
  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi")
  )

  # Compute correlation
  result <- ds.vertCor("D_aligned", variables, datasources = conns)

  # Check structure
  expect_s3_class(result, "ds.cor")
  expect_true(is.matrix(result$correlation))
  expect_equal(nrow(result$correlation), 4)
  expect_equal(ncol(result$correlation), 4)
  expect_equal(result$method, "MHE-CKKS")

  # Check diagonal is approximately 1
  expect_equal(unname(diag(result$correlation)), rep(1, 4), tolerance = 1e-4)

  # Check symmetry
  expect_equal(result$correlation, t(result$correlation), tolerance = 1e-6)

  # Check values are in [-1, 1]
  expect_true(all(result$correlation >= -1 - 1e-4 & result$correlation <= 1 + 1e-4))
})

test_that("ds.vertCor works with three servers", {
  env <- setup_dslite_env(seed = 201)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Align
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Use all 3 servers
  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi"),
    server3 = c("glucose", "cholesterol")
  )

  result <- ds.vertCor("D_aligned", variables, datasources = conns)

  expect_equal(nrow(result$correlation), 6)
  expect_equal(ncol(result$correlation), 6)
  expect_equal(length(result$servers), 3)
})

test_that("ds.vertCor print method works", {
  env <- setup_dslite_env(seed = 202)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi")
  )

  result <- ds.vertCor("D_aligned", variables, datasources = conns)

  expect_output(print(result), "Correlation Matrix")
  expect_output(print(result), "MHE-CKKS")
})

# =============================================================================
# Test PCA
# =============================================================================

test_that("ds.vertPCA performs PCA across servers", {
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

  # Perform PCA
  pca_result <- ds.vertPCA("D_aligned", variables,
                            n_components = 3, datasources = conns)

  # Check structure
  expect_s3_class(pca_result, "ds.pca")
  expect_true("loadings" %in% names(pca_result))
  expect_true("eigenvalues" %in% names(pca_result))
  expect_true("variance_pct" %in% names(pca_result))
  expect_true("cumulative_pct" %in% names(pca_result))
  expect_true("correlation" %in% names(pca_result))

  # Check dimensions
  expect_equal(nrow(pca_result$loadings), 4)  # 4 variables
  expect_equal(ncol(pca_result$loadings), 3)  # 3 components requested

  # Check variance explained sums correctly
  expect_equal(unname(pca_result$cumulative_pct[3]),
               unname(sum(pca_result$variance_pct[1:3])),
               tolerance = 1e-6)

  # Eigenvalues should be non-negative and decreasing
  expect_true(all(pca_result$eigenvalues >= 0))
  expect_true(all(diff(pca_result$eigenvalues) <= 1e-10))  # Allow tiny numerical error
})

test_that("ds.vertPCA respects n_components parameter", {
  env <- setup_dslite_env(seed = 204)
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
  expect_equal(ncol(pca2$loadings), 2)
  expect_equal(length(pca2$eigenvalues), 2)

  # Request all components (4 variables = max 4 components)
  pca4 <- ds.vertPCA("D_aligned", variables, n_components = 4, datasources = conns)
  expect_equal(ncol(pca4$loadings), 4)

  # Request more than available (should cap at max)
  pca_max <- ds.vertPCA("D_aligned", variables, n_components = 10, datasources = conns)
  expect_equal(ncol(pca_max$loadings), 4)
})

test_that("ds.vertPCA print method works", {
  env <- setup_dslite_env(seed = 205)
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

# =============================================================================
# Validation: Compare with ground truth
# =============================================================================

test_that("ds.vertCor matches ground truth correlation", {
  set.seed(400)
  n <- 100

  # Create correlated data with known structure
  ids <- paste0("PAT", sprintf("%04d", 1:n))

  age <- round(rnorm(n, 45, 15))
  weight <- round(rnorm(n, 70, 12), 1)
  height <- round(0.3 * weight + rnorm(n, 120, 8), 1)  # Correlated with weight
  bmi <- round(weight / (height/100)^2, 1)

  # Ground truth
  full_data <- data.frame(age = age, weight = weight, height = height, bmi = bmi)
  expected_cor <- cor(full_data)

  # Create shuffled server data
  order1 <- sample(n)
  order2 <- sample(n)

  data_server1 <- data.frame(
    patient_id = ids[order1],
    age = age[order1],
    weight = weight[order1],
    stringsAsFactors = FALSE
  )

  data_server2 <- data.frame(
    patient_id = ids[order2],
    height = height[order2],
    bmi = bmi[order2],
    stringsAsFactors = FALSE
  )

  # Setup DSLite
  dslite_server <- DSLite::newDSLiteServer(
    tables = list(server1 = data_server1, server2 = data_server2)
  )
  dslite_server$config(DSLite::defaultDSConfiguration(include = c("dsVert")))

  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "server1", url = "dslite_server",
                 table = "server1", driver = "DSLiteDriver")
  builder$append(server = "server2", url = "dslite_server",
                 table = "server2", driver = "DSLiteDriver")

  assign("dslite_server", dslite_server, envir = globalenv())
  conns <- DSI::datashield.login(builder$build(), assign = TRUE, symbol = "D")

  on.exit({
    try(DSI::datashield.logout(conns), silent = TRUE)
    if (exists("dslite_server", envir = globalenv())) {
      rm("dslite_server", envir = globalenv())
    }
  })

  # Align and compute
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(
    server1 = c("age", "weight"),
    server2 = c("height", "bmi")
  )

  result <- ds.vertCor("D_aligned", variables, datasources = conns)

  # Reorder expected to match variable order
  expected_reordered <- expected_cor[c("age", "weight", "height", "bmi"),
                                      c("age", "weight", "height", "bmi")]

  max_error <- max(abs(result$correlation - expected_reordered))
  expect_lt(max_error, 1e-3)
})
