test_that("LMM ML contracts sign covariance-major candidates and patient caps", {
  f <- .grouped_cross_client_fixture("lmm")
  raw <- f$raw
  raw$parameters <- list(objective = "ml", variance_grid = list(
    list(residual_variance = .25, random_intercept_variance = 0),
    list(residual_variance = 1, random_intercept_variance = .25)))
  spec <- .dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema)
  count <- length(spec$beta_grid)
  expect_length(spec$candidate_grid, 2 * count)
  expect_equal(vapply(spec$candidate_grid, `[[`, numeric(1), "variance_index"),
    rep(1:2, each = count))
  expect_equal(vapply(spec$candidate_grid, `[[`, numeric(1), "beta_index"), rep(seq_len(count), 2))
  caps <- unlist(spec$sensitivity$per_cluster_caps)
  expect_length(caps, 2 * count)
  expect_equal(spec$sensitivity$raw_l1_sensitivity, sum(caps))
  expect_equal(unlist(spec$sensitivity$maximum_coordinates), spec$grouping$cluster_capacity * caps)
  expect_match(spec$numeric_contract$certificate_sha256, "^[0-9a-f]{64}$")
  expect_identical(spec$numeric_contract$profile, "grouped-lmm-ml-variance-f264-q64-log-up-v1")
  expect_false(spec$numeric_contract$full_likelihood_value)
  artifact <- .dsvert_dp_grouped_grid_cross_artifact(spec)
  expect_equal(artifact$coordinate_count, 2 * count)
  unsigned <- list(version = f$contract$version, spec = spec, artifact = artifact,
    source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, artifact))
  validate <- function(x) .dsvert_dp_grouped_grid_cross_contract_validate(f$sign(x), f$policy, f$schema_manifest)
  expect_identical(validate(unsigned)$spec$parameters$objective, "ml")
  for (mutate in list(
    function(x) { x$spec$candidate_order <- rev(x$spec$candidate_order); x },
    function(x) { x$spec$candidate_grid[[1]]$variance_index <- 2; x },
    function(x) { x$spec$numeric_contract$certificate_sha256 <- strrep("0", 64); x },
    function(x) { x$spec$numeric_contract$logdet_addition <- "after_clamp"; x },
    function(x) { x$spec$sensitivity$per_cluster_caps[[1]] <- 1; x })) {
    expect_error(validate(mutate(unsigned)), class = "dsvert_dp_public_failure")
  }
  values <- rep(10, 2 * count)
  values[count+1] <- 0
  fit <- .dsvert_dp_grouped_grid_cross_moment(values, spec, artifact)
  expect_equal(fit$selected_candidate, count+1)
  expect_equal(unname(fit$normalized_coefficients), unlist(spec$beta_grid[[1]], use.names = FALSE))
  expect_identical(fit$parameters, c(list(objective = "ml"), spec$parameters$variance_grid[[2]]))
  expect_identical(fit$loss_objective, spec$numeric_contract$objective)

})

test_that("LMM ML rejects ambiguous covariance order and uncertified objectives", {
  f <- .grouped_cross_client_fixture("lmm")
  pair <- list(residual_variance = 1, random_intercept_variance = .25)
  other <- list(residual_variance = .25, random_intercept_variance = 0)
  for (parameters in list(
    list(objective = "reml", variance_grid = list(pair)),
    list(objective = "ml", variance_grid = list(pair, pair)),
    list(objective = "ml", variance_grid = list(pair, other)),
    list(objective = "ml", variance_grid = list(list(objective = "ml", variance_grid = list(pair)))),
    list(objective = "ml", variance_grid = list(list(residual_variance = .3, random_intercept_variance = 0))))) {
    raw <- f$raw; raw$parameters <- parameters
    expect_error(.dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema),
      class = "dsvert_dp_public_failure")
  }
})
