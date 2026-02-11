# Validation tests: Compare local computation with distributed dsVert computation

skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")

library(DSLite)
library(dsVert)

# Skip if MHE binary not available
skip_if_not(dsVert::mheAvailable(), "MHE tool binary not available")

# =============================================================================
# Correlation: Local vs Distributed
# =============================================================================

test_that("ds.vertCor matches local cor() on same data", {
  # Create a known dataset
  set.seed(500)
  n <- 100

  # Common IDs
  ids <- paste0("PAT", sprintf("%04d", 1:n))

  # Generate correlated variables with known relationships
  x1 <- rnorm(n)
  x2 <- 0.5 * x1 + sqrt(0.75) * rnorm(n)  # Correlation ~0.5 with x1
  y1 <- rnorm(n)
  y2 <- -0.3 * y1 + sqrt(0.91) * rnorm(n) # Correlation ~-0.3 with y1

  # Create full dataset for local computation
  full_data <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)

  # LOCAL COMPUTATION: Ground truth
  local_cor <- cor(full_data)

  # DISTRIBUTED COMPUTATION: Split data across servers
  order1 <- sample(n)
  order2 <- sample(n)

  data_server1 <- data.frame(
    patient_id = ids[order1],
    x1 = x1[order1],
    x2 = x2[order1],
    stringsAsFactors = FALSE
  )

  data_server2 <- data.frame(
    patient_id = ids[order2],
    y1 = y1[order2],
    y2 = y2[order2],
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

  # Align data
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(
    server1 = c("x1", "x2"),
    server2 = c("y1", "y2")
  )

  # Compute distributed correlation via MHE
  distributed_result <- ds.vertCor("D_aligned", variables, datasources = conns)
  distributed_cor <- distributed_result$correlation

  # Reorder to match variable order
  local_cor_reordered <- local_cor[c("x1", "x2", "y1", "y2"),
                                    c("x1", "x2", "y1", "y2")]

  # VALIDATION: Results should be very close
  max_error <- max(abs(distributed_cor - local_cor_reordered))

  # MHE has approximation error ~10^-4 to 10^-3
  expect_lt(max_error, 0.01)  # Within 1%

  # Specific cross-server correlations should match
  expect_equal(distributed_cor["x1", "y1"], local_cor["x1", "y1"], tolerance = 0.01)
  expect_equal(distributed_cor["x2", "y2"], local_cor["x2", "y2"], tolerance = 0.01)
})

# =============================================================================
# PCA: Local vs Distributed
# =============================================================================

test_that("ds.vertPCA loadings match local PCA on same data", {
  set.seed(501)
  n <- 100

  ids <- paste0("PAT", sprintf("%04d", 1:n))

  # Generate data
  a1 <- rnorm(n)
  a2 <- 0.8 * a1 + 0.2 * rnorm(n)  # Strong correlation
  b1 <- rnorm(n)
  b2 <- rnorm(n)

  full_data <- data.frame(a1 = a1, a2 = a2, b1 = b1, b2 = b2)

  # LOCAL PCA using correlation matrix
  local_cor <- cor(full_data)
  local_eigen <- eigen(local_cor, symmetric = TRUE)
  local_loadings <- local_eigen$vectors
  local_variance <- 100 * local_eigen$values / sum(local_eigen$values)

  # DISTRIBUTED: Split across servers
  order1 <- sample(n)
  order2 <- sample(n)

  data_server1 <- data.frame(
    patient_id = ids[order1],
    a1 = a1[order1],
    a2 = a2[order1],
    stringsAsFactors = FALSE
  )

  data_server2 <- data.frame(
    patient_id = ids[order2],
    b1 = b1[order2],
    b2 = b2[order2],
    stringsAsFactors = FALSE
  )

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

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(server1 = c("a1", "a2"), server2 = c("b1", "b2"))

  distributed_pca <- ds.vertPCA("D_aligned", variables, datasources = conns)

  # VALIDATION: Variance explained should match closely
  expect_equal(unname(distributed_pca$variance_pct[1:4]), local_variance, tolerance = 1)

  # Note: Loadings signs can flip (eigenvector sign is arbitrary)
  # Compare absolute values
  for (i in 1:4) {
    # Take absolute because sign can flip
    dist_abs <- abs(distributed_pca$loadings[, i])
    local_abs <- abs(local_loadings[, i])

    # Should be close (within ~10% due to HE error)
    max_diff <- max(abs(dist_abs - local_abs))
    expect_lt(max_diff, 0.15)
  }
})

# =============================================================================
# Cross-Server Correlation is Accurate
# =============================================================================

test_that("cross-server correlations are accurately computed", {
  # This test specifically validates that variables on DIFFERENT servers
  # have their correlations correctly computed

  set.seed(502)
  n <- 100

  ids <- paste0("PAT", sprintf("%04d", 1:n))

  # Create data with KNOWN cross-server correlations
  base <- rnorm(n)
  x_on_server1 <- base + 0.1 * rnorm(n)
  y_on_server2 <- 0.7 * base + sqrt(1 - 0.7^2) * rnorm(n)

  # Ground truth correlation between x and y
  true_cor_xy <- cor(x_on_server1, y_on_server2)

  # Split data
  order1 <- sample(n)
  order2 <- sample(n)

  data_server1 <- data.frame(
    patient_id = ids[order1],
    x = x_on_server1[order1],
    stringsAsFactors = FALSE
  )

  data_server2 <- data.frame(
    patient_id = ids[order2],
    y = y_on_server2[order2],
    stringsAsFactors = FALSE
  )

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

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(server1 = c("x"), server2 = c("y"))

  result <- ds.vertCor("D_aligned", variables, datasources = conns)
  computed_cor_xy <- result$correlation["x", "y"]

  # The cross-server correlation should match the true value closely
  expect_equal(computed_cor_xy, true_cor_xy, tolerance = 0.02)

  # Also verify the correlation is significant (should be ~0.7 given our setup)
  expect_gt(abs(computed_cor_xy), 0.5)
})

# =============================================================================
# Three-Server Correlation Matches Local
# =============================================================================

test_that("correlation across 3 servers matches local", {
  set.seed(503)
  n <- 80

  ids <- paste0("PAT", sprintf("%04d", 1:n))

  # Generate correlated data across 3 groups
  v1 <- rnorm(n)
  v2 <- rnorm(n)
  v3 <- 0.4 * v1 + 0.3 * v2 + sqrt(0.75) * rnorm(n)

  full_data <- data.frame(v1 = v1, v2 = v2, v3 = v3)
  local_cor <- cor(full_data)

  # Distribute across 3 servers
  orders <- list(sample(n), sample(n), sample(n))

  data_servers <- list(
    server1 = data.frame(patient_id = ids[orders[[1]]],
                         v1 = v1[orders[[1]]], stringsAsFactors = FALSE),
    server2 = data.frame(patient_id = ids[orders[[2]]],
                         v2 = v2[orders[[2]]], stringsAsFactors = FALSE),
    server3 = data.frame(patient_id = ids[orders[[3]]],
                         v3 = v3[orders[[3]]], stringsAsFactors = FALSE)
  )

  dslite_server <- DSLite::newDSLiteServer(tables = data_servers)
  dslite_server$config(DSLite::defaultDSConfiguration(include = c("dsVert")))

  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "server1", url = "dslite_server",
                 table = "server1", driver = "DSLiteDriver")
  builder$append(server = "server2", url = "dslite_server",
                 table = "server2", driver = "DSLiteDriver")
  builder$append(server = "server3", url = "dslite_server",
                 table = "server3", driver = "DSLiteDriver")

  assign("dslite_server", dslite_server, envir = globalenv())
  conns <- DSI::datashield.login(builder$build(), assign = TRUE, symbol = "D")

  on.exit({
    try(DSI::datashield.logout(conns), silent = TRUE)
    if (exists("dslite_server", envir = globalenv())) {
      rm("dslite_server", envir = globalenv())
    }
  })

  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  variables <- list(server1 = "v1", server2 = "v2", server3 = "v3")

  result <- ds.vertCor("D_aligned", variables, datasources = conns)

  # All cross-server pairs should match local computation
  expect_equal(result$correlation["v1", "v2"], local_cor["v1", "v2"], tolerance = 0.02)
  expect_equal(result$correlation["v1", "v3"], local_cor["v1", "v3"], tolerance = 0.02)
  expect_equal(result$correlation["v2", "v3"], local_cor["v2", "v3"], tolerance = 0.02)
})
