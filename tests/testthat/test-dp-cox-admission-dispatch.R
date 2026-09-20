test_that("Cox discovery and generic layout preserve private padded source lanes", {
  skip_if_not_installed("dsVert")
  server <- get(".dsvert_dp_cox_cross_transport_layout", asNamespace("dsVert"))
  for (owners in c(2L, 3L, 5L)) {
    f <- .cox_cross_client_fixture(capacity = 5, owners = owners)
    artifact <- .dsvert_dp_cox_cross_workload_artifact(f$contract)
    families <- setNames(rep(list(list(artifacts = list())), 10L), c(
      "admitted_count", "numeric_moments", "numeric_pair_moments", "gaussian_models",
      "fixed_numeric_histograms", "categorical_marginals", "categorical_pairs",
      "correlation_artifacts", "describe_artifacts", "survival_artifacts"))
    families$admitted_count <- list(owner_peer = "site_a", dataset = "aligned")
    families$survival_artifacts <- list()
    families$gaussian_models$artifacts <- list(cox_grid = artifact)
    manifest <- list(capsule_identity = list(capsule_id = strrep("6", 64)),
      admission = list(unit_capacity = 5, adjacency = f$policy$adjacency),
      bounds = list(numeric_grid_bits = 8),
      workload = list(coordinate_count = artifact$coordinate_count + 1L,
        families = families, capsule_mechanism = list(mechanism = "discrete-laplace")))
    expect_named(.dsvert_dp_glm_grid_cross_artifacts(manifest), "cox_grid")
    layout <- .dsvert_dp_gaussian_cross_layout_client(manifest)
    expect_true(layout$enabled)
    expect_identical(layout[setdiff(names(layout), "enabled")], server(manifest, artifact))
    time <- Filter(function(block) isTRUE(block$private_time_validity), layout$blocks)
    expect_length(time, 1L)
    expect_identical(time[[1L]]$kind, "validity")
    expect_identical(time[[1L]]$length, 8L)
    expect_identical(unlist(layout$source_peers), names(f$keys))
    # Signed Cox artifacts use the authenticated staged runtime.
    expect_true(.dsvert_dp_synopsis_supported_glm_grid_cross_v1(manifest))
    context <- list(pinset = f$policy$peer_pinset, designated = f$policy$designated_noise_peers)
    preflight <- function(value = manifest, schema = f$schema_manifest) {
      .dsvert_dp_glm_grid_cross_preflight(value, context, .dsvert_joint_dp_client_json(schema))
    }
    expect_true(preflight())
    changed <- manifest
    changed$workload$families$gaussian_models$artifacts$cox_grid$time$owner_peer <- "site_b"
    expect_error(preflight(changed))
    changed <- f$schema_manifest
    changed$datasets$aligned$columns$time$upper <- 21
    expect_error(preflight(schema = changed))
    changed <- f$contract
    changed$spec$complete_case <- "ignore_missing"
    bad <- manifest
    bad$workload$families$gaussian_models$artifacts$cox_grid <-
      .dsvert_dp_cox_cross_workload_artifact(changed)
    expect_error(preflight(bad))
    actual <- .dsvert_dp_gaussian_artifact(manifest, "aligned", "cox_grid", "site_a",
      f$policy$adjacency, 256, 5)
    expect_identical(actual$statistic_maximum,
      unlist(f$contract$spec$sensitivity$maximum_coordinates, use.names = FALSE))
    expect_identical(actual$time, artifact$time)
    expect_identical(actual$event, artifact$event)
    expect_null(actual$outcome)
  }
})

test_that("Cox rejects mixed-family runtime and transport manifests", {
  f <- .cox_cross_client_fixture()
  cox <- .dsvert_dp_cox_cross_workload_artifact(f$contract)
  grouped <- .grouped_cross_client_fixture()
  lmm <- .dsvert_dp_grouped_cross_workload_artifact(grouped$contract)
  manifest <- list(workload = list(families = list(gaussian_models = list(
    artifacts = list(cox_grid = cox, lmm_grid = lmm)))))
  expect_length(.dsvert_dp_glm_grid_cross_artifacts(manifest), 2L)
  expect_false(.dsvert_dp_synopsis_supported_glm_grid_cross_v1(manifest))
  expect_error(.dsvert_dp_gaussian_cross_layout_client(manifest))
  expect_error(.dsvert_dp_cox_cross_orchestrate(
    .dsvert_joint_dp_client_json(manifest), manifest, list(), list(),
    f$policy, f$schema_manifest,
    function(...) stop("mixed-family transport must not run"), list()),
    class = "dsvert_dp_public_failure")
})

test_that("Cox reconstruction rejects changed arithmetic, source and request bindings", {
  f <- .cox_cross_client_fixture(capacity = 5)
  artifact <- .dsvert_dp_cox_cross_workload_artifact(f$contract)
  args <- list(artifact = artifact, data_name = "aligned", analysis_id = "cox_grid",
    owner_peer = "site_a", adjacency = f$policy$adjacency, scale = 256, capacity = 5)
  for (field in c("data_name", "analysis_id", "owner_peer", "adjacency", "scale", "capacity")) {
    bad <- args
    bad[[field]] <- if (is.numeric(bad[[field]])) bad[[field]] + 1 else "wrong"
    expect_error(do.call(.dsvert_dp_cox_cross_client_artifact, bad))
  }
  for (field in c("numeric_contract", "sensitivity", "source_contract", "artifact")) {
    bad <- f$contract
    if (field %in% c("numeric_contract", "sensitivity")) bad$spec[[field]] <- list() else bad[[field]] <- list()
    changed <- args
    changed$artifact$signed_contract <- .dsvert_joint_dp_client_json(bad)
    expect_error(do.call(.dsvert_dp_cox_cross_client_artifact, changed))
  }
  bad <- artifact; bad$statistic_maximum <- list(0, 0, 0)
  args$artifact <- bad
  expect_error(do.call(.dsvert_dp_cox_cross_client_artifact, args))
  expect_error(.dsvert_dp_glm_grid_cross_client_artifact(artifact, "aligned", "cox_grid",
    "site_a", f$policy$adjacency, 256, 5, "lmm"))
  expect_error(.cox_cross_client_fixture(capacity = 401),
    class = "dsvert_dp_public_failure")
})
