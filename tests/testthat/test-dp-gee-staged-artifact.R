test_that("internal staged GEE reconstruction checks the entire beta-only vector", {
  for (family in c("binomial_gee", "poisson_gee")) {
    f <- .grouped_cross_client_fixture(family, "ar1")
    raw <- f$raw
    raw$parameters$composition <- "staged_fixed_rho_v1"
    spec <- .dsvert_dp_grouped_grid_cross_spec(raw, f$policy, f$schema)
    signed <- .dsvert_dp_grouped_grid_cross_artifact(spec)
    contract <- f$sign(list(version = f$contract$version, spec = spec, artifact = signed,
      source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, signed)))
    artifact <- .dsvert_dp_grouped_cross_workload_artifact(contract)
    read <- function(a) .dsvert_dp_grouped_cross_client_artifact(a, spec$dataset,
      spec$analysis_id, spec$owner_peer, spec$adjacency, 2^spec$numeric_grid_bits,
      spec$observation_capacity, family)
    result <- read(artifact)
    expect_equal(length(result$statistic_maximum), signed$coordinate_count)
    expect_equal(result$statistic_maximum, unlist(spec$sensitivity$maximum_coordinates))
    expect_equal(result$candidate_loss_bounds, unlist(spec$sensitivity$per_cluster_caps))
    expect_identical(result$source_coordinate_scaling,
      "all_coordinates_already_on_common_numeric_lattice_v1")
    expect_null(spec$candidate_grid)
    # An internal validator is not public producer/reader admission.
    manifest <- list(workload = list(families = list(gaussian_models = list(
      artifacts = list(gee = artifact)))))
    expect_length(.dsvert_dp_glm_grid_cross_artifacts(manifest), 0)
    for (field in c("coordinate_count", "statistic_maximum", "candidate_loss_bounds")) {
      changed <- artifact; changed[[field]] <- 0
      expect_error(read(changed), class = "dsvert_dp_public_failure")
    }
    for (field in c("beta_encoded", "staged_numeric", "candidate_order", "sensitivity")) {
      changed <- contract
      if (field == "beta_encoded") changed$spec$beta_encoded[[1L]][[1L]] <- "1"
      if (field == "staged_numeric") changed$spec$staged_numeric$BreadCap <- "1"
      if (field == "candidate_order") changed$spec$candidate_order <- rev(spec$candidate_order)
      if (field == "sensitivity") changed$spec$sensitivity$maximum_coordinates[[2L]] <- 0
      expect_error(read(.dsvert_dp_grouped_cross_workload_artifact(changed)),
        class = "dsvert_dp_public_failure")
    }
  }
})
