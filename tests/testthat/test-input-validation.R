# Tests for input validation in client functions (no server needed)

# =============================================================================
# ds.vertCor input validation
# =============================================================================

test_that("ds.vertCor rejects non-character data_name", {
  expect_error(ds.vertCor(123, list(a = "x", b = "y")),
               "data_name must be a single character string")
  expect_error(ds.vertCor(c("a", "b"), list(a = "x", b = "y")),
               "data_name must be a single character string")
})

test_that("ds.vertCor rejects non-list variables", {
  expect_error(ds.vertCor("D", c("x", "y")),
               "variables must be a named list")
})

test_that("ds.vertCor rejects unnamed list variables", {
  expect_error(ds.vertCor("D", list(c("x"), c("y"))),
               "variables must be a named list")
})

test_that("ds.vertCor rejects single server", {
  expect_error(ds.vertCor("D", list(server1 = c("x", "y"))),
               "At least 2 servers required")
})

# =============================================================================
# ds.vertPCA input validation
# =============================================================================

test_that("ds.vertPCA rejects non-character data_name", {
  expect_error(ds.vertPCA(123, list(a = "x", b = "y")),
               "data_name must be a single character string")
})

test_that("ds.vertPCA rejects non-list variables", {
  expect_error(ds.vertPCA("D", "x"),
               "variables must be a named list")
})

# =============================================================================
# ds.vertGLM input validation
# =============================================================================

test_that("ds.vertGLM rejects non-character data_name", {
  expect_error(ds.vertGLM(123, "y", list(a = "x", b = "z")),
               "data_name must be a single character string")
})

test_that("ds.vertGLM rejects non-list x_vars", {
  expect_error(ds.vertGLM("D", "y", c("x", "z")),
               "x_vars must be a named list")
})

# =============================================================================
# ds.hashId input validation
# =============================================================================

test_that("ds.hashId rejects non-character data_name", {
  expect_error(ds.hashId(123, "id"),
               "data_name must be a single character string")
})

test_that("ds.hashId rejects non-character id_variable", {
  expect_error(ds.hashId("D", 123),
               "id_col must be a single character string")
})

# =============================================================================
# ds.alignRecords input validation
# =============================================================================

test_that("ds.alignRecords rejects non-character data_name", {
  expect_error(ds.alignRecords(123, "id", c("a", "b")),
               "data_name must be a single character string")
})

test_that("ds.alignRecords rejects non-character id_variable", {
  expect_error(ds.alignRecords("D", 123, c("a", "b")),
               "id_col must be a single character string")
})

test_that("ds.alignRecords rejects non-character reference_hashes", {
  expect_error(ds.alignRecords("D", "id", 123),
               "reference_hashes must be a non-empty character vector")
})
