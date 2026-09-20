.lmm_client_admission_fixture <- function(owners = 2L, large = FALSE) {
  f <- .grouped_cross_client_fixture("lmm", owners = owners)
  f$raw$parameters <- list(objective = "ml", variance_grid = list(
    list(residual_variance = .25, random_intercept_variance = 0),
    list(residual_variance = 1, random_intercept_variance = .25)))
  if (large) {
    f$policy$unit_capacity <- 2000
    f$raw$grouping$cluster_capacity <- 500
  }
  spec <- .dsvert_dp_grouped_grid_cross_spec(f$raw, f$policy, f$schema)
  artifact <- .dsvert_dp_grouped_grid_cross_artifact(spec)
  f$contract <- f$sign(list(version = .DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_VERSION,
    spec = spec, artifact = artifact,
    source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, artifact)))
  f
}

test_that("client production admission reconstructs authenticated LMM ML", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .lmm_client_admission_fixture(owners)
    value <- .dsvert_dp_glm_grid_profile_admit(f$contract, f$policy, f$schema_manifest)
    expect_identical(value$spec$family, "lmm")
    expect_equal(value$artifact$coordinate_count, 4)
    expect_identical(value$spec$numeric_contract$profile,
      "grouped-lmm-ml-variance-f264-q64-log-up-v1")
    expect_equal(vapply(value$spec$candidate_grid, `[[`, numeric(1), "variance_index"),
                 c(1, 1, 2, 2))
    for (mutate in list(
      function(x) { x$spec$numeric_contract$certificate_sha256 <- strrep("0", 64); x },
      function(x) { x$spec$candidate_order <- rev(x$spec$candidate_order); x },
      function(x) { x$spec$sensitivity$per_cluster_caps[[1]] <- 1; x })) {
      expect_error(.dsvert_dp_glm_grid_profile_admit(
        f$sign(mutate(f$contract)), f$policy, f$schema_manifest),
        class = "dsvert_dp_public_failure")
    }
    raw <- list(version = "lmm_grid_cross_v1", dataset = "aligned",
      contract = .dsvert_joint_dp_client_json(f$contract))
    fragments <- list(describe = list(), survival = list(),
      gaussian = list(grouped = raw), vertical_cross = list())
    admitted <- .dsvert_dp_capsule_manifest_fragments(fragments)
    expect_identical(admitted$gaussian$grouped, .dsvert_joint_dp_client_canonical(raw))
  }
  for (family in .DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_FAMILIES) {
    f <- .grouped_cross_client_fixture(family)
    expect_error(.dsvert_dp_glm_grid_profile_admit(f$contract, f$policy, f$schema_manifest),
                 class = "dsvert_dp_public_failure")
  }
})

test_that("client n2000 signed bounds are exact and preserve small contracts", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .lmm_client_admission_fixture(owners, large = TRUE)
    value <- .dsvert_dp_glm_grid_profile_admit(f$contract, f$policy, f$schema_manifest)
    expect_identical(value$spec$grouping$capacity_profile,
                     "lmm-ml-n2000-c500-b4-p3-j4-v1")
    expect_equal(value$source_contract$private_layout$padded_units, 2000)
    expect_equal(unlist(value$spec$sensitivity$maximum_coordinates),
                 500 * unlist(value$spec$sensitivity$per_cluster_caps))
    changed <- f$contract
    changed$spec$grouping$capacity_profile <- "capacity-proven"
    expect_error(.dsvert_dp_glm_grid_profile_admit(
      f$sign(changed), f$policy, f$schema_manifest), class = "dsvert_dp_public_failure")
  }
  f <- .lmm_client_admission_fixture(large = TRUE)
  for (mutate in list(
    function(x) { x$raw$grouping$cluster_capacity <- 65; x },
    function(x) { x$raw$grouping$cluster_capacity <- 499; x },
    function(x) { x$raw$grouping$cluster_capacity <- 501; x },
    function(x) { x$raw$grouping$max_patients_per_cluster <- 5; x },
    function(x) { x$policy$unit_capacity <- 1999; x },
    function(x) { x$raw$parameters <- x$raw$parameters$variance_grid[[1]]; x },
    function(x) { x$raw$parameters$variance_grid[[3]] <-
      list(residual_variance = 2, random_intercept_variance = .25); x })) {
    changed <- mutate(f)
    expect_error(.dsvert_dp_grouped_grid_cross_spec(
      changed$raw, changed$policy, changed$schema), class = "dsvert_dp_public_failure")
  }
  small <- .lmm_client_admission_fixture()
  expect_null(small$contract$spec$grouping$capacity_profile)
  other <- .grouped_cross_client_fixture("binomial_glmm")
  other$raw$grouping$cluster_capacity <- 500
  expect_error(.dsvert_dp_grouped_grid_cross_spec(other$raw, other$policy, other$schema),
               class = "dsvert_dp_public_failure")
})
