test_that("ordinal DP postprocessing preserves signed order and original units", {
  fixture <- .categorical_cross_client_fixture("ordinal")
  fit <- .dsvert_dp_ordinal_grid_cross_postprocess(fixture$contract, c(100, 90))
  expect_identical(fit$selected_candidate, 2L)
  expect_equal(fit$coefficients, c("site_a$x" = 0.5, "site_b$z" = 0.1))
  expect_equal(fit$thresholds, c(B = -1, A = 0))
  expect_identical(fit$ordered_levels, c("B", "A", "C"))
  expect_equal(fit$selected_dp_negative_log_likelihood, 90 / 256)
  expect_null(fit$standard_errors)
  expect_identical(.dsvert_dp_ordinal_grid_cross_postprocess(
    fixture$contract, c(90, 90))$selected_candidate, 1L)
})

test_that("ordinal signed candidates enforce the certified finite-log domain", {
  fixture <- .categorical_cross_client_fixture("ordinal")
  mutations <- list(
    function(raw) { raw$candidate_grid[[1L]]$beta[[1L]] <- 1; raw },
    function(raw) { raw$candidate_grid[[1L]]$beta <- c(0, 8, 1); raw },
    function(raw) { raw$candidate_grid[[1L]]$thresholds <- c(0, 1/32); raw },
    function(raw) { raw$candidate_grid[[1L]]$thresholds <- c(-2^-58, 1/16-2^-57); raw },
    function(raw) { raw$candidate_grid[[1L]]$thresholds <- c(-8, 8.001); raw },
    function(raw) { raw$candidate_grid[[1L]]$thresholds <- c(1, -1); raw },
    function(raw) { raw$ordered_levels <- c("A", "A", "C"); raw },
    function(raw) { raw$ordered_levels <- c("A", "B", "D"); raw })
  for (mutate in mutations) {
    .categorical_cross_client_reject(function() {
      fixture$registration$spec(mutate(fixture$raw), fixture$policy, fixture$schema)
    })
  }
})
