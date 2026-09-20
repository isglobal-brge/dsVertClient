test_that("poisson GH5 grids bind variance-major candidates and full patient sensitivity", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .grouped_cross_client_fixture("poisson_glmm", owners = owners)
    raw <- f$raw
    raw$beta_grid <- raw$beta_grid[1:2]
    raw$parameters <- list(variance_grid = list(0, .25), quadrature = "gh5_fixed_v1")
    spec <- .dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema)
    count <- length(spec$beta_grid)
    expect_length(spec$candidate_grid, 2 * count)
    expect_equal(vapply(spec$candidate_grid, `[[`, numeric(1L), "variance_index"),
      rep(1:2, each = count))
    expect_equal(vapply(spec$candidate_grid, `[[`, numeric(1L), "beta_index"),
      rep(seq_len(count), 2))
    expected_caps <- unlist(lapply(c(0, .25), function(tau) {
      fixed <- raw
      fixed$parameters <- list(random_intercept_variance = tau, quadrature = "gh5_fixed_v1")
      .dsvert_dp_grouped_grid_cross_spec(fixed, f$policy, f$schema)$sensitivity$per_cluster_caps
    }), use.names = FALSE)
    caps <- unlist(spec$sensitivity$per_cluster_caps, use.names = FALSE)
    expect_equal(caps, expected_caps)
    expect_equal(spec$sensitivity$raw_l1_sensitivity, sum(caps))
    expect_equal(spec$sensitivity$raw_l2_sensitivity, sqrt(sum(caps^2)))
    expect_equal(unlist(spec$sensitivity$maximum_coordinates), spec$grouping$cluster_capacity*caps)
    replacement <- .dsvert_dp_grouped_grid_cross_sensitivity(spec$beta_grid, spec$family,
      spec$max_outcome, spec$numeric_grid_bits, spec$grouping, spec$parameters,
      "replace_one_fixed_cohort", spec$numeric_contract)
    expect_equal(replacement$raw_l1_sensitivity, 2*sum(caps))
    expect_equal(replacement$raw_l2_sensitivity, 2*sqrt(sum(caps^2)))
    expect_identical(spec$numeric_contract$profile_sha256,
      "f72e66abaf2e503a809f23d4563418d2889843174109398ae48b02f0ec7edb84")
    expect_identical(spec$numeric_contract$certificate_sha256,
      "1fe5a19c74969134731d5d3210b9646e52783192570256806e0ac7a1e4ae9c94")
    expect_lt(spec$numeric_contract$composition_error_bound, 1)
    expect_false(spec$numeric_contract$quadrature_error_included)
    artifact <- .dsvert_dp_grouped_grid_cross_artifact(spec)
    expect_equal(artifact$coordinate_count, length(caps))
    source <- .dsvert_dp_grouped_grid_cross_source_contract(spec, artifact)
    outcome <- source$private_layout$blocks[[length(spec$predictor_order)+1L]]
    expect_equal(outcome$value_fraction_bits, 0)
    expect_false(identical(spec$outcome$owner_peer, spec$grouping$owner_peer))
    unsigned <- list(version = f$contract$version, spec = spec, artifact = artifact,
      source_contract = source)
    validate <- function(x) .dsvert_dp_grouped_grid_cross_contract_validate(f$sign(x), f$policy, f$schema_manifest)
    expect_equal(vapply(validate(unsigned)$spec$candidate_grid, `[[`, numeric(1L),
      "variance_index"), rep(1:2, each = count))
    for (mutate in list(
      function(x) { x$spec$candidate_order <- rev(x$spec$candidate_order); x },
      function(x) { x$spec$candidate_grid[[1]]$variance_index <- 2; x },
      function(x) { x$spec$candidate_grid <- x$spec$candidate_grid[-1]; x },
      function(x) { x$artifact$coordinate_count <- count; x },
      function(x) { x$spec$numeric_contract$certificate_sha256 <- strrep("0", 64); x },
      function(x) { x$spec$numeric_contract$composition_error_bound <- 0; x },
      function(x) { x$spec$sensitivity$per_cluster_caps[[1]] <- 1; x },
      function(x) { x$spec$sensitivity$maximum_coordinates[[1]] <- 1; x },
      function(x) { x$spec$sensitivity$raw_l1_sensitivity <- 1; x },
      function(x) { x$spec$sensitivity$raw_l2_sensitivity <- 1; x })) {
      expect_error(validate(mutate(unsigned)), class = "dsvert_dp_public_failure")
    }
  }
})

test_that("poisson GH5 variance grids reject ambiguous or uncertified values", {
  f <- .grouped_cross_client_fixture("poisson_glmm")
  for (grid in list(list(), list(.25, 0), list(0, 0), list(.125), c(0, .25),
      list(a = 0), list(NA_real_), list(Inf), list("0"), list(list(0)))) {
    raw <- f$raw
    raw$parameters <- list(variance_grid = grid, quadrature = "gh5_fixed_v1")
    expect_error(.dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema),
      class = "dsvert_dp_public_failure")
  }
  for (tau in c(0, .25)) {
    raw <- f$raw
    raw$parameters <- list(variance_grid = list(tau), quadrature = "gh5_fixed_v1")
    spec <- .dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema)
    expect_length(spec$candidate_grid, length(spec$beta_grid))
  }
  for (parameters in list(
      list(variance_grid = list(0), random_intercept_variance = 0, quadrature = "gh5_fixed_v1"),
      list(variance_grid = list(0), quadrature = "adaptive"))) {
    raw <- f$raw; raw$parameters <- parameters
    expect_error(.dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema),
      class = "dsvert_dp_public_failure")
  }
  numeric <- .dsvert_dp_grouped_grid_cross_numeric("poisson_glmm",
    list(max_patients_per_cluster = 16), list(variance_grid = list(0, .25)))
  expect_lt(numeric$composition_error_bound, .096)
})

test_that("poisson GH5 C500 admission stays inside the named measurement domain", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .grouped_cross_client_fixture("poisson_glmm", owners = owners)
    raw <- f$raw
    raw$parameters <- list(variance_grid = list(0, .25), quadrature = "gh5_fixed_v1")
    raw$beta_grid <- raw$beta_grid[1:2]
    raw$grouping$cluster_capacity <- 500
    policy <- f$policy; policy$unit_capacity <- 2000
    spec <- .dsvert_dp_grouped_grid_cross_spec(raw, policy, f$schema)
    expect_identical(spec$grouping$capacity_profile, "poisson-glmm-gh5-n2000-c500-b4-p3-j4-v1")
    expect_length(spec$candidate_grid, 4)
    expect_length(spec$participating_peers, owners)
    for (mutate in list(
      function(x) { x$grouping$cluster_capacity <- 499; x },
      function(x) { x$grouping$cluster_capacity <- 501; x },
      function(x) { x$grouping$max_patients_per_cluster <- 5; x },
      function(x) { x$parameters <- f$raw$parameters; x })) {
      expect_error(.dsvert_dp_grouped_grid_cross_spec(mutate(raw), policy, f$schema),
        class = "dsvert_dp_public_failure")
    }
    policy$unit_capacity <- 1999
    expect_error(.dsvert_dp_grouped_grid_cross_spec(raw, policy, f$schema),
      class = "dsvert_dp_public_failure")
  }
})


test_that("poisson GH5 postprocessing selects variance and keeps first ties", {
  f <- .grouped_cross_client_fixture("poisson_glmm")
  raw <- f$raw
  raw$parameters <- list(variance_grid = list(0, .25), quadrature = "gh5_fixed_v1")
  spec <- .dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema)
  artifact <- .dsvert_dp_grouped_grid_cross_artifact(spec)
  values <- rep(10, length(spec$candidate_grid))
  selected <- length(spec$beta_grid)+1L
  values[selected] <- 0
  fit <- .dsvert_dp_grouped_grid_cross_moment(values, spec, artifact)
  expect_equal(fit$selected_candidate, selected)
  expect_equal(unname(fit$normalized_coefficients), unlist(spec$beta_grid[[1]]))
  expect_identical(fit$parameters,
    list(random_intercept_variance = .25, quadrature = "gh5_fixed_v1"))
  expect_identical(fit$loss_objective, "factorial_free_gh5_selection_plus_2_per_live_row_v1")
  values[1] <- 0
  fit <- .dsvert_dp_grouped_grid_cross_moment(values, spec, artifact)
  expect_identical(fit$selected_candidate, 1L)
  expect_identical(fit$parameters$random_intercept_variance, 0)
  contract <- f$sign(list(version = f$contract$version, spec = spec, artifact = artifact,
    source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, artifact)))
  expect_error(dp_glmm_grid(spec$outcome$reference, unlist(spec$predictor_order),
    contract, f$policy, f$schema_manifest, datasources = list()),
    class = "dsvert_dp_public_failure")
})

test_that("poisson GH5 public artifacts rebuild numeric caps and lattice identity", {
  f <- .grouped_cross_client_fixture("poisson_glmm")
  raw <- f$raw
  raw$parameters <- list(variance_grid = list(0, .25), quadrature = "gh5_fixed_v1")
  spec <- .dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema)
  signed <- .dsvert_dp_grouped_grid_cross_artifact(spec)
  contract <- f$sign(list(version = f$contract$version, spec = spec, artifact = signed,
    source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, signed)))
  artifact <- .dsvert_dp_grouped_cross_workload_artifact(contract)
  validate <- function(value) .dsvert_dp_grouped_cross_client_artifact(value,
    spec$dataset, spec$analysis_id, NULL, spec$adjacency, 2^spec$numeric_grid_bits,
    spec$observation_capacity)
  admitted <- validate(artifact)
  expect_equal(admitted$statistic_maximum,
    unlist(spec$sensitivity$maximum_coordinates, use.names = FALSE))
  expect_identical(admitted$source_coordinate_scaling,
    "all_coordinates_already_on_common_numeric_lattice_v1")
  for (mutate in list(
    function(x) { x$spec$candidate_grid <- rev(x$spec$candidate_grid); x },
    function(x) { x$spec$candidate_order <- rev(x$spec$candidate_order); x },
    function(x) { x$spec$numeric_contract$certificate_sha256 <- strrep("0", 64); x },
    function(x) { x$spec$sensitivity$maximum_coordinates[[1]] <- 1; x },
    function(x) { x$source_contract$private_layout$blocks[[1]]$length <- 15; x })) {
    expect_error(validate(.dsvert_dp_grouped_cross_workload_artifact(f$sign(mutate(contract)))),
      class = "dsvert_dp_public_failure")
  }
  artifact$source_coordinate_scaling <- NULL
  expect_error(validate(artifact), class = "dsvert_dp_public_failure")
})


test_that("new poisson variance-grid contracts reject p4 and J5 before execution", {
  f <- .grouped_cross_client_fixture("poisson_glmm", owners = 5L)
  raw <- f$raw
  raw$parameters <- list(variance_grid = list(0), quadrature = "gh5_fixed_v1")
  raw$beta_grid <- lapply(c(-.5, -.25, 0, .25, .5), function(intercept) c(intercept, rep(0, 3)))
  raw$beta_grid <- raw$beta_grid[order(vapply(raw$beta_grid,
    .dsvert_joint_dp_client_json, character(1L)), method = "radix")]
  admitted <- raw; admitted$beta_grid <- raw$beta_grid[1:4]
  expect_length(.dsvert_dp_grouped_grid_cross_spec(admitted, f$policy, f$schema)$candidate_grid, 4)
  expect_error(.dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema),
    class = "dsvert_dp_public_failure")
  admitted$parameters$variance_grid <- list(0, .25)
  admitted$beta_grid <- raw$beta_grid[1:2]
  expect_length(.dsvert_dp_grouped_grid_cross_spec(admitted, f$policy, f$schema)$candidate_grid, 4)
  admitted$beta_grid <- raw$beta_grid[1:3]
  expect_error(.dsvert_dp_grouped_grid_cross_spec(admitted, f$policy, f$schema),
    class = "dsvert_dp_public_failure")
  rounded <- raw; rounded$beta_grid <- list(c(.01, .03, .15, .81))
  expect_error(.dsvert_dp_grouped_grid_cross_spec(rounded, f$policy, f$schema),
    class = "dsvert_dp_public_failure")
  rounded$beta_grid[[1]][4] <- .81 - 2^-50
  expect_length(.dsvert_dp_grouped_grid_cross_spec(rounded, f$policy, f$schema)$candidate_grid, 1)
  schema <- f$schema
  columns <- schema$unsigned$datasets[[raw$dataset]]$columns
  descriptor <- columns[[sub("^[^$]+\\$", "", raw$predictor_order[[1L]])]]
  schema$unsigned$datasets[[raw$dataset]]$columns$extra <- descriptor
  raw$predictor_order <- sort(c(raw$predictor_order, paste0(descriptor$owner_peer, "$extra")),
    method = "radix")
  raw$beta_grid <- list(rep(0, 5))
  fixed <- raw; fixed$parameters <- f$raw$parameters
  expect_length(.dsvert_dp_grouped_grid_cross_spec(fixed, f$policy, schema)$predictors, 4)
  expect_error(.dsvert_dp_grouped_grid_cross_spec(raw, f$policy, schema),
    class = "dsvert_dp_public_failure")
})
