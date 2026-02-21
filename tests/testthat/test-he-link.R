# Tests for HE-Link eta_privacy parameter in ds.vertGLM

# =============================================================================
# eta_privacy input validation
# =============================================================================

test_that("ds.vertGLM rejects invalid eta_privacy", {
  expect_error(
    ds.vertGLM("D", "y", list(a = "x1", b = "x2"),
               y_server = "a", eta_privacy = "invalid"),
    "eta_privacy must be"
  )
})

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
