test_that("staged GEE reconstruction checks the entire beta-only vector", {
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
    manifest <- list(workload = list(families = list(gaussian_models = list(
      artifacts = list(gee = artifact)))))
    expect_named(.dsvert_dp_glm_grid_cross_artifacts(manifest), "gee")
    legacy <- manifest
    legacy$workload$families$gaussian_models$artifacts$gee <-
      .dsvert_dp_grouped_cross_workload_artifact(f$contract)
    expect_length(.dsvert_dp_glm_grid_cross_artifacts(legacy), 0)
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


test_that("fixed-rho GEE public admission and preflight preserve separated owner sources", {
  for (family in c("binomial_gee", "poisson_gee")) for (owners in c(2L, 3L, 5L)) {
    f <- .grouped_cross_client_fixture(family, "exchangeable", owners)
    expect_error(.dsvert_dp_glm_grid_profile_admit(f$contract, f$policy,
      f$schema_manifest), class = "dsvert_dp_public_failure")
    f$raw$parameters$composition <- "staged_fixed_rho_v1"
    spec <- .dsvert_dp_grouped_grid_cross_spec(f$raw, f$policy, f$schema)
    signed <- .dsvert_dp_grouped_grid_cross_artifact(spec)
    contract <- f$sign(list(version = f$contract$version, spec = spec, artifact = signed,
      source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, signed)))
    admitted <- .dsvert_dp_glm_grid_profile_admit(contract, f$policy, f$schema_manifest)
    expect_identical(admitted$spec$family, family)
    expect_identical(admitted$artifact$correlation_contract, "signed-analyst-fixed-rho-v1")
    artifact <- .dsvert_dp_grouped_cross_workload_artifact(contract)
    manifest <- list(admission = list(unit_capacity = 16, adjacency = f$policy$adjacency),
      bounds = list(numeric_grid_bits = 8),
      capsule_identity = list(capsule_id = strrep("a", 64)),
      workload = list(coordinate_count = signed$coordinate_count + 1,
        capsule_mechanism = list(mechanism = "discrete-laplace", source_context_hash = strrep("b", 64)),
        families = list(admitted_count = list(owner_peer = "site_a", dataset = "aligned"),
          numeric_moments = list(artifacts = list()), numeric_pair_moments = list(artifacts = list()),
          gaussian_models = list(artifacts = list(grouped = artifact)),
          fixed_numeric_histograms = list(artifacts = list()), categorical_marginals = list(artifacts = list()),
          categorical_pairs = list(sets = list()), correlation_artifacts = list(),
          describe_artifacts = list(), survival_artifacts = list())))
    preflight <- function(value) .dsvert_dp_glm_grid_cross_preflight(value,
      list(pinset = f$policy$peer_pinset, designated = f$policy$designated_noise_peers),
      .dsvert_joint_dp_client_json(f$schema_manifest))
    expect_silent(preflight(manifest))
    layout <- .dsvert_dp_gaussian_cross_layout_client(manifest)
    expect_equal(unlist(layout$source_peers), names(f$policy$peer_pinset))
    expect_equal(unlist(layout$computation_peers), c("site_a", "site_b"))
    outcome <- layout$blocks[["grouped::site_b$y::value"]]
    route <- layout$blocks[["grouped::site_a$cluster::value"]]
    expect_equal(outcome$fraction_bits, 0)
    expect_equal(outcome$maximum, if (family == "binomial_gee") 1 else 4)
    expect_identical(route$owner_peer, "site_a")
    expect_identical(outcome$owner_peer, "site_b")
    expect_true(route$private_routing_input)
    expect_identical(.dsvert_dp_staged_grouped_tag(artifact, "-staged-v1"),
      "gee-fixed-rho-staged-v1")
    for (mutate in list(
        function(x) { x$version <- "bounded-lmm-cross-grid-v1"; x },
        function(x) { x$composition <- NULL; x },
        function(x) { x$correlation_contract <- "estimated-alpha-v1"; x },
        function(x) { x$statistic_maximum[[1L]] <- 1; x })) {
      changed <- manifest
      changed$workload$families$gaussian_models$artifacts$grouped <- mutate(artifact)
      expect_error(.dsvert_dp_grouped_cross_client_artifact(
        changed$workload$families$gaussian_models$artifacts$grouped,
        "aligned", "grouped", NULL, f$policy$adjacency, 256, 16, family),
        class = "dsvert_dp_public_failure")
    }
    changed <- manifest; changed$workload$coordinate_count <- 52
    expect_error(preflight(changed), class = "dsvert_dp_public_failure")
  }
})
