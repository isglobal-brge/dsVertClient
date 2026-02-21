# Tests for secure_agg mode in ds.vertGLM (client-side validation)

# =============================================================================
# eta_privacy="auto" mode selection
# =============================================================================

test_that("ds.vertGLM auto selects secure_agg for K>=3", {
  # With 3 servers (2 non-label), auto should select secure_agg.
  # Since there's no DataSHIELD connection, this will fail after validation
  # but before Phase 0. We just need to verify it gets past the validation.
  expect_error(
    ds.vertGLM("D", "y",
               list(a = "x1", b = "x2", c = "x3"),
               y_server = "a", family = "gaussian"),
    "No DataSHIELD connections"
  )
})

test_that("ds.vertGLM auto selects he_link for K=2 binomial", {
  # With 2 servers (1 non-label) and binomial, auto selects he_link.
  # This should fail at HE-Link validation only if we force eta_privacy,
  # otherwise it should get past validation to connection setup.
  expect_error(
    ds.vertGLM("D", "y",
               list(a = "x1", b = "x2"),
               y_server = "a", family = "binomial"),
    "No DataSHIELD connections"
  )
})

test_that("ds.vertGLM auto selects transport for K=2 gaussian", {
  # With 2 servers and gaussian, auto selects transport.
  expect_error(
    ds.vertGLM("D", "y",
               list(a = "x1", b = "x2"),
               y_server = "a", family = "gaussian"),
    "No DataSHIELD connections"
  )
})

# =============================================================================
# eta_privacy="secure_agg" policy enforcement
# =============================================================================

test_that("ds.vertGLM rejects secure_agg with K=2", {
  expect_error(
    ds.vertGLM("D", "y",
               list(a = "x1", b = "x2"),
               y_server = "a", eta_privacy = "secure_agg"),
    "secure_agg requires >= 3 servers"
  )
})

# =============================================================================
# eta_privacy validation
# =============================================================================

test_that("ds.vertGLM rejects invalid eta_privacy values", {
  expect_error(
    ds.vertGLM("D", "y", list(a = "x1", b = "x2"),
               y_server = "a", eta_privacy = "invalid"),
    "eta_privacy must be"
  )
})

test_that("ds.vertGLM accepts valid eta_privacy values", {
  # These should all get past validation to the connection check
  for (mode in c("auto", "transport", "secure_agg")) {
    # secure_agg needs K>=3
    if (mode == "secure_agg") {
      expect_error(
        ds.vertGLM("D", "y",
                   list(a = "x1", b = "x2", c = "x3"),
                   y_server = "a", eta_privacy = mode),
        "No DataSHIELD connections"
      )
    } else {
      expect_error(
        ds.vertGLM("D", "y",
                   list(a = "x1", b = "x2"),
                   y_server = "a", eta_privacy = mode),
        "No DataSHIELD connections"
      )
    }
  }
})

# =============================================================================
# Existing he_link tests still pass
# =============================================================================

test_that("ds.vertGLM rejects he_link with >2 servers", {
  expect_error(
    ds.vertGLM("D", "y",
               list(a = "x1", b = "x2", c = "x3"),
               y_server = "a", family = "binomial",
               eta_privacy = "he_link", log_n = 14L),
    "exactly 2 servers"
  )
})

test_that("ds.vertGLM rejects he_link with gaussian family", {
  expect_error(
    ds.vertGLM("D", "y",
               list(a = "x1", b = "x2"),
               y_server = "a", family = "gaussian",
               eta_privacy = "he_link", log_n = 14L),
    "binomial family only"
  )
})

test_that("ds.vertGLM rejects he_link with log_n < 14", {
  expect_error(
    ds.vertGLM("D", "y",
               list(a = "x1", b = "x2"),
               y_server = "a", family = "binomial",
               eta_privacy = "he_link", log_n = 12L),
    "log_n >= 14"
  )
})
