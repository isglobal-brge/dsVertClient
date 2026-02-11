# Integration tests for record alignment functions

skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")

library(DSLite)
library(dsVert)

# =============================================================================
# Test ID Validation (requires dev version of dsVert)
# =============================================================================

test_that("ds.validateIdFormat works across multiple servers", {
  env <- setup_dslite_env(seed = 100)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Validate ID format
  validation <- ds.validateIdFormat("D", "patient_id", datasources = conns)

  expect_s3_class(validation, "ds.id.validation")
  expect_true(validation$valid)
  expect_true(validation$format_match)
  expect_equal(nrow(validation$servers), 3)
  expect_equal(validation$servers$n_obs, rep(env$n, 3))
  expect_equal(validation$servers$n_unique, rep(env$n, 3))
  expect_equal(validation$servers$n_missing, rep(0, 3))
})

test_that("ds.validateIdFormat detects pattern mismatches", {
  # Note: DSLite has issues parsing regex patterns with special characters
  # Pattern matching is tested via unit tests, here we just verify basic functionality
  env <- setup_dslite_env(seed = 101)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Pattern that matches - using simple literal pattern
  validation1 <- ds.validateIdFormat("D", "patient_id",
                                      pattern = "PAT",
                                      datasources = conns)
  expect_true(validation1$pattern_match)

  # Pattern that doesn't match
  validation2 <- ds.validateIdFormat("D", "patient_id",
                                      pattern = "XYZ",
                                      datasources = conns)
  expect_false(validation2$pattern_match)
})

# =============================================================================
# Test Hashing and Alignment
# =============================================================================

test_that("ds.hashId returns hashes from a server", {
  env <- setup_dslite_env(seed = 102)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Get hashes from server1
  hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])

  expect_type(hashes, "list")
  expect_named(hashes, c("hashes", "n"))
  expect_equal(hashes$n, env$n)
  expect_length(hashes$hashes, env$n)
  # SHA-256 produces 64 character hex strings
  expect_true(all(nchar(hashes$hashes) == 64))
})

test_that("ds.alignRecords aligns data across servers", {
  env <- setup_dslite_env(seed = 103)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Get reference hashes from server1
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])

  # Align all servers
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Verify alignment by checking observation counts
  for (server in names(conns)) {
    count_call <- call("getObsCountDS", "D_aligned")
    result <- DSI::datashield.aggregate(conns = conns[server], expr = count_call)
    expect_equal(result[[1]]$n_obs, env$n)
  }
})

test_that("complete alignment workflow with validation", {
  env <- setup_dslite_env(seed = 104)
  conns <- DSI::datashield.login(env$login_data, assign = TRUE, symbol = "D")

  on.exit(teardown_dslite_env(conns))

  # Step 1: Validate IDs
  validation <- ds.validateIdFormat("D", "patient_id", datasources = conns)
  expect_true(validation$valid)

  # Step 2: Get reference hashes
  ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])

  # Step 3: Align
  ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                  newobj = "D_aligned", datasources = conns)

  # Step 4: Verify counts match
  counts <- sapply(names(conns), function(server) {
    result <- DSI::datashield.aggregate(
      conns = conns[server],
      expr = call("getObsCountDS", "D_aligned")
    )
    result[[1]]$n_obs
  })

  expect_true(all(counts == env$n))
})
