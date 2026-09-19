test_that("all five grouped contracts require both actual authorities", {
  for (family in .DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_FAMILIES) {
    f <- .grouped_cross_client_fixture(family)
    valid <- .dsvert_dp_grouped_grid_cross_contract_validate(
      f$contract, f$policy, f$schema_manifest)
    expect_identical(valid, .dsvert_joint_dp_client_canonical(f$contract))
    expect_identical(valid$artifact$implementation_state, "cross_owner_exact_gc_materialized")
    expect_identical(valid$artifact$cross_owner_state, "exact_gc_to_joint_dp_vector_v1")
    expect_identical(valid$spec$sensitivity$dp_unit, "patient")
    expect_identical(unlist(valid$spec$participating_peers), c("site_a", "site_b"))
    for (peer in c("site_a", "site_b")) {
      bad <- f$contract; bad$signatures[[peer]] <- NULL
      expect_error(.dsvert_dp_grouped_grid_cross_contract_validate(
        bad, f$policy, f$schema_manifest), class = "dsvert_dp_public_failure")
      bad <- f$schema_manifest; bad$signatures[[peer]] <- NULL
      expect_error(.dsvert_dp_grouped_grid_cross_contract_validate(
        f$contract, f$policy, bad), class = "dsvert_dp_public_failure")
    }
  }
})

test_that("resigned grouped tampering cannot weaken public bounds or semantics", {
  f <- .grouped_cross_client_fixture("poisson_gee", "ar1")
  mutations <- list(
    function(x) { x$spec$sensitivity$raw_l1_sensitivity <- 1; x },
    function(x) { x$spec$sensitivity$raw_l2_sensitivity <- 1; x },
    function(x) { x$spec$sensitivity$maximum_coordinates[[1]] <- 1; x },
    function(x) { x$spec$numeric_contract$per_cluster_error_bound <- 0; x },
    function(x) { x$spec$grouping$patient_rule <- "every_record"; x },
    function(x) { x$spec$grouping$ordering <- "reorder_on_deletion"; x },
    function(x) { x$spec$grouping$routing <- "public_indices"; x },
    function(x) { x$artifact$cross_owner_state <- "reserved_not_materialized"; x },
    function(x) { x$artifact$implementation_state <- "same_owner_materialized"; x },
    function(x) { x$artifact$result_evidence_required <- FALSE; x },
    function(x) { x$artifact$coordinate_order <- "estimating_function_norm"; x },
    function(x) { x$source_contract$purpose <- "plaintext"; x },
    function(x) { x$spec$parameters$rho <- .9; x },
    function(x) { x$spec$predictor_order <- rev(x$spec$predictor_order); x },
    function(x) { x$spec$beta_grid <- rev(x$spec$beta_grid); x })
  for (mutate in mutations) {
    error <- tryCatch(.dsvert_dp_grouped_grid_cross_contract_validate(
      f$sign(mutate(f$contract)), f$policy, f$schema_manifest), error = identity)
    expect_s3_class(error, "dsvert_dp_public_failure")
    expect_identical(conditionMessage(error),
      "[dsvert_dp_public_failure:v1] Protected capsule operation failed.")
  }
})

test_that("grouped postprocessing selects likelihood and omits inference", {
  for (family in .DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_FAMILIES) {
    f <- .grouped_cross_client_fixture(family)
    spec <- f$contract$spec; artifact <- f$contract$artifact
    width <- artifact$coordinate_count / length(spec$beta_grid)
    values <- rep(0, artifact$coordinate_count)
    values[1] <- 2; values[1+width] <- 1
    # Deliberately huge second-candidate bread/meat must not affect selection.
    if (width > 1) values[seq.int(width+2, length(values))] <-
      unlist(spec$sensitivity$maximum_coordinates)[seq.int(width+2, length(values))]
    result <- .dsvert_dp_grouped_grid_cross_moment(values, spec, artifact)
    expect_identical(result$selected_candidate, 2L)
    expect_equal(unname(result$coefficients), c(.25, .25, .25))
    expect_null(result$standard_errors)
    expect_null(result$p_values)
    values[1] <- 1
    expect_identical(.dsvert_dp_grouped_grid_cross_moment(
      values, spec, artifact)$selected_candidate, 1L)
    for (bad in list(values[-1], c(NA_real_, values[-1]), values + .5,
                    c(-1, values[-1]), unlist(spec$sensitivity$maximum_coordinates)+1)) {
      expect_error(.dsvert_dp_grouped_grid_cross_moment(bad, spec, artifact),
                   class = "dsvert_dp_public_failure")
    }
  }
})

test_that("public grouped functions fail closed before fused release promotion", {
  registration <- .dsvert_dp_grouped_grid_cross_register()
  expect_length(registration$versions, 5)
  for (family in .DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_FAMILIES) {
    f <- .grouped_cross_client_fixture(family)
    entry <- if (family == "lmm") dp_lmm_grid else if (grepl("_glmm$", family))
      dp_glmm_grid else dp_gee_grid
    expect_error(entry("site_a$y", c("site_a$x", "site_b$z"), f$contract,
      f$policy, f$schema_manifest), class = "dsvert_dp_public_failure")
    expect_error(entry("site_a$y", "site_a$x", f$contract, f$policy,
      f$schema_manifest), class = "dsvert_dp_public_failure")
  }
})

test_that("patient replacement and movement use two complete cluster ranges", {
  for (family in .DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_FAMILIES) {
    f <- .grouped_cross_client_fixture(family)
    spec <- f$contract$spec
    replacement <- .dsvert_dp_grouped_grid_cross_sensitivity(spec$beta_grid,
      family, spec$max_outcome, spec$numeric_grid_bits, spec$grouping,
      spec$parameters, "replace_one_fixed_cohort", spec$numeric_contract)
    expect_equal(replacement$raw_l1_sensitivity, 2*spec$sensitivity$raw_l1_sensitivity)
    expect_equal(replacement$raw_l2_sensitivity, 2*spec$sensitivity$raw_l2_sensitivity)
    expect_identical(replacement$maximum_coordinates, spec$sensitivity$maximum_coordinates)
  }
})
