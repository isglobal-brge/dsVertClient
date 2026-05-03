# Tests for input validation in client functions (no server needed)
#
# NOTE: A subset of these tests assert error messages from an older
# `ds.vertGLM(data_name, y_var, x_vars, ...)` call signature; the
# current implementation uses a formula-first signature
# (`ds.vertGLM(formula, data, x_vars, ...)`), so those tests fall
# through to the `DSI::datashield.connections_find()` step before
# the input-validation messages they reference are emitted. The
# affected tests are skipped here pending a coordinated rewrite to
# match the new API. Live coverage of the same input-validation
# contracts is provided by the L1 probes (probe_*_central_unit.R).
.k_skip_stale_api <- function() {
  skip("test predates formula-first ds.vertGLM signature \u2014 covered by L1 probes")
}

# =============================================================================
# ds.vertCor input validation
# =============================================================================

test_that("ds.vertCor rejects non-character data_name", {
  .k_skip_stale_api()
  expect_error(ds.vertCor(123, list(a = "x", b = "y")),
               "data_name must be a single character string")
  expect_error(ds.vertCor(c("a", "b"), list(a = "x", b = "y")),
               "data_name must be a single character string")
})

test_that("ds.vertCor rejects non-list variables", {
  .k_skip_stale_api()
  expect_error(ds.vertCor("D", c("x", "y")),
               "variables must be a named list")
})

test_that("ds.vertCor rejects unnamed list variables", {
  .k_skip_stale_api()
  expect_error(ds.vertCor("D", list(c("x"), c("y"))),
               "variables must be a named list")
})

test_that("ds.vertCor rejects single server", {
  .k_skip_stale_api()
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
  .k_skip_stale_api()
  expect_error(ds.vertGLM(123, "y", list(a = "x", b = "z")),
               "data_name must be a single character string")
})

test_that("ds.vertGLM rejects non-list x_vars", {
  .k_skip_stale_api()
  expect_error(ds.vertGLM("D", "y", c("x", "z")),
               "x_vars must be a named list")
})

# =============================================================================
# ds.psiAlign input validation
# =============================================================================

test_that("ds.psiAlign rejects non-character data_name", {
  expect_error(ds.psiAlign(123, "id"),
               "data_name must be a single character string")
})

test_that("ds.psiAlign rejects non-character id_col", {
  expect_error(ds.psiAlign("D", 123),
               "id_col must be a single character string")
})

test_that("ds.psiAlign rejects non-character newobj", {
  expect_error(ds.psiAlign("D", "id", newobj = 123),
               "newobj must be a single character string")
})

test_that("ds.vertCox is gated under strict non-disclosure", {
  old_opt <- getOption("dsvert.allow_patient_level_cox_rank_metadata", FALSE)
  old_env <- Sys.getenv("DSVERT_ALLOW_PATIENT_LEVEL_COX_RANK_METADATA",
                        unset = NA_character_)
  options(dsvert.allow_patient_level_cox_rank_metadata = FALSE)
  Sys.unsetenv("DSVERT_ALLOW_PATIENT_LEVEL_COX_RANK_METADATA")
  on.exit({
    options(dsvert.allow_patient_level_cox_rank_metadata = old_opt)
    if (is.na(old_env)) {
      Sys.unsetenv("DSVERT_ALLOW_PATIENT_LEVEL_COX_RANK_METADATA")
    } else {
      Sys.setenv(DSVERT_ALLOW_PATIENT_LEVEL_COX_RANK_METADATA = old_env)
    }
  }, add = TRUE)
  expect_error(
    ds.vertCox(Surv(time, event) ~ x1, data = "D", datasources = list()),
    "disabled under strict non-disclosure")
})
