.gee_ar1_grid_fixture <- function() {
  capacity <- 20
  scale <- 256
  candidate_grid <- list(
    list(beta = c(0, 0), rho = 0),
    list(beta = c(0.2, 0.5), rho = 0.5))
  candidates <- dsVertClient:::.dsvert_dp_gee_ar1_grid_candidates(
    candidate_grid, 4L)
  loss_bounds <- vapply(candidates, `[[`, numeric(1L), "loss_bound")
  raw <- ceiling(loss_bounds * scale)
  artifact <- list(
    version = "bounded-gaussian-ar1-working-gls-grid-v1",
    spec_version = "gaussian_ar1_working_gls_grid_v1", analysis_id = "gee_ar1",
    dataset = "protected", owner_peer = "server_a",
    outcome = list(column = "y", lower = -2, upper = 6),
    cluster = list(column = "site", levels = c("a", "b", "c")),
    order = list(column = "visit", lower = 0, upper = 10),
    predictors = list(x = list(column = "x", lower = 0, upper = 10)),
    predictor_order = "x", intercept = TRUE,
    design_terms = c("(Intercept)", "x"), observation_capacity = capacity,
    max_patients_per_cluster = 4L,
    candidate_grid = lapply(candidates, function(candidate) list(
      beta = candidate$beta, rho = candidate$rho)),
    candidate_order = "canonical_beta_rho_grid_v1",
    candidate_loss_bounds = as.list(loss_bounds), numeric_grid_bits = 8L,
    coordinate_count = 2L,
    coordinate_order =
      "signed_candidate_grid_cluster_gaussian_ar1_working_gls_loss_v1",
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    repeated_record_policy = paste(
      "require_one_complete_bounded_row_per_admitted_patient_with_one",
      "consistent_public_cluster_level_and_strict_within_cluster_order_v1"),
    missingness_policy = paste(
      "missing_or_nonfinite_outcome_predictor_or_order_or_missing_or",
      "inconsistent_cluster_excludes_patient_and_order_ties_reject_v1"),
    contribution_domain = paste(
      "one_bounded_patient_can_change_one_clipped_cluster_gaussian_ar1",
      "working_gls_loss_per_signed_candidate_v1"),
    statistic_maximum = as.list(capacity * raw),
    source_raw_l1_sensitivity = sum(raw),
    source_raw_l2_sensitivity = sqrt(sum(raw^2)),
    natural_l1_sensitivity = sum(raw) / scale,
    natural_l2_sensitivity = sqrt(sum(raw^2)) / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste(
      "one_patient_can_change_one_entire_clipped_cluster_ar1_working",
      "gls_loss_by_at_most_its_signed_bound_v1"),
    estimation_scope =
      "bounded_gaussian_ar1_working_gls_finite_signed_beta_rho_grid_v1",
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  list(artifact = artifact, manifest = list(workload = list(families = list(
    gaussian_models = list(artifacts = list(gee_ar1 = artifact))))))
}

test_that("Gaussian AR1 working-GLS artifact validates and selects only a signed point", {
  fixture <- .gee_ar1_grid_fixture()
  artifact <- dsVertClient:::.dsvert_dp_gee_ar1_grid_artifact(
    fixture$manifest, "protected", "gee_ar1", "server_a",
    "add_remove_patient", 256, 20)
  fit <- dsVertClient:::.dsvert_dp_gee_ar1_grid_moment(c(100, 90), artifact)

  expect_identical(artifact$order$column, "visit")
  expect_identical(fit$status, "ok")
  expect_identical(fit$selected_candidate, 2L)
  expect_equal(fit$coefficients, c(`(Intercept)` = 0.2, x = 0.05))
  expect_equal(fit$working_correlation, 0.5)
  tampered <- fixture$manifest
  tampered$workload$families$gaussian_models$artifacts$gee_ar1$
    candidate_loss_bounds[[1L]] <- 0
  expect_error(dsVertClient:::.dsvert_dp_gee_ar1_grid_artifact(
    tampered, "protected", "gee_ar1", "server_a",
    "add_remove_patient", 256, 20), "descriptor is invalid")
  expect_error(dsVertClient:::.dsvert_dp_gee_ar1_grid_moment(
    c(100, max(unlist(artifact$statistic_maximum)) + 1), artifact),
    "violates its signed bounds")
})

test_that("Gaussian AR1 working-GLS Synopsis uses one validated signed block", {
  fixture <- .gee_ar1_grid_fixture()
  artifact <- fixture$artifact
  block <- list(key = "gee_ar1", descriptor = artifact)
  context <- list(
    lattice = list(output_lattice_scale = 256), layout = list(),
    manifest = fixture$manifest, release = list(), adjacency = "add_remove_patient")
  coordinates <- c(100, 90)
  verification <- list(
    integrity_valid = TRUE, authenticity = "session_transport_anchored",
    artifact = artifact, coordinates = coordinates,
    cohort_id = "cohort", logical_snapshot = list(version = "v1"))
  testthat::local_mocked_bindings(
    .dsvert_dp_datasources = function(value) value,
    .dsvert_dp_synopsis_vector_run = function(...) list(run = TRUE),
    .dsvert_dp_vector_context = function(...) context,
    .dsvert_dp_vector_public_metadata = function(...) list(),
    .dsvert_dp_capsule_single_block = function(...) {
      list(start = 1L, end = 1L, length = 1L)
    },
    .dsvert_dp_vector_block_capacity = function(...) 20,
    .dsvert_dp_capsule_vector_blocks = function(...) list(block),
    .dsvert_dp_capsule_vector_values = function(...) coordinates,
    .dsvert_dp_gaussian_synopsis_certificate_build = function(...) {
      list(certificate_sha256 = strrep("a", 64L))
    },
    ds.validateDPGaussianCertificate = function(...) verification,
    .package = "dsVertClient")

  released <- dsVertClient:::.dsvert_dp_gee_ar1_grid_synopsis_release(
    "protected", "gee_ar1", datasources = list(server_a = list()),
    .aggregate = function(...) stop("unexpected aggregate", call. = FALSE))
  expect_identical(released$coordinates, coordinates)
  expect_identical(released$moment$selected_candidate, 2L)
  expect_equal(released$moment$working_correlation, 0.5)
})

test_that("GEE AR1 frontdoor admits only the signed finite-grid contract", {
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_dp_gee_ar1_grid_impl = function(...) {
      seen <<- list(...)
      structure(list(
        status = "public_certified_gaussian_ar1_working_gls_finite_grid",
        family = "gaussian", corstr = "ar1",
        coefficients = c(`(Intercept)` = 0.2, x = 0.05),
        working_correlation = 0.5),
        class = c("dsvert_dp_gaussian_ar1_gee", "ds.vertGEE", "list"))
    },
    .package = "dsVertClient")

  fit <- ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", id_col = "site",
    order_col = "visit", corstr = "ar1", analysis_id = "gee-ar1",
    datasources = list(site_a = structure(list(), class = "mock")))
  expect_s3_class(fit, "dsvert_dp_gaussian_ar1_gee")
  expect_identical(seen$data_name, "study")
  expect_identical(seen$analysis_id, "gee-ar1")
  expect_identical(seen$id_col, "site")
  expect_identical(seen$order_col, "visit")
  expect_output(print(fit), "AR\\(1\\) working-GLS")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", id_col = "site",
    corstr = "ar1", analysis_id = "gee-ar1", datasources = list()),
    "requires distinct id_col and order_col")
  expect_error(ds.vertGEE(
    y ~ x, data = "study", family = "gaussian", id_col = "site",
    order_col = "visit", corstr = "ar1", analysis_id = "gee-ar1",
    lambda = 0, datasources = list()), "does not accept legacy controls")
})
