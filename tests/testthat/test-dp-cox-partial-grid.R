.cox_partial_grid_fixture <- function() {
  capacity <- 20
  scale <- 256
  candidates <- dsVertClient:::.dsvert_dp_cox_partial_grid_candidates(
    list(c(0.5), c(0)), 1L, capacity)
  bounds <- vapply(candidates, `[[`, numeric(1L), "loss_bound")
  raw <- ceiling(bounds * scale)
  artifact <- list(
    version = "bounded-cox-partial-likelihood-grid-v1",
    spec_version = "cox_partial_likelihood_grid_v1", analysis_id = "cox_grid",
    dataset = "study", owner_peer = "server_a",
    time = list(column = "time", lower = 0, upper = 10),
    event = list(column = "event", levels = c("censor", "event"),
      censor = "censor", event_level = "event"),
    time_grid = c(2, 5, 10),
    predictors = list(x = list(column = "x", lower = 0, upper = 10)),
    predictor_order = "x", intercept = FALSE, design_terms = "x",
    observation_capacity = capacity,
    candidate_grid = lapply(candidates, `[[`, "beta"),
    candidate_order = "canonical_beta_grid_v1",
    candidate_loss_bounds = as.list(bounds), numeric_grid_bits = 8L,
    coordinate_count = 2L,
    coordinate_order =
      "signed_candidate_grid_breslow_cox_partial_likelihood_loss_v1",
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    repeated_record_policy =
      "require_one_complete_bounded_row_per_admitted_patient_v1",
    missingness_policy =
      "missing_or_nonfinite_time_or_predictor_or_event_excludes_patient_v1",
    contribution_domain =
      "one_bounded_patient_can_change_the_clipped_cox_cohort_loss_per_signed_candidate_v1",
    statistic_maximum = as.list(raw),
    source_raw_l1_sensitivity = sum(raw),
    source_raw_l2_sensitivity = sqrt(sum(raw^2)),
    natural_l1_sensitivity = sum(raw) / scale,
    natural_l2_sensitivity = sqrt(sum(raw^2)) / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis =
      "one_patient_can_change_one_entire_clipped_breslow_cox_loss_by_at_most_its_signed_bound_v1",
    estimation_scope =
      "bounded_same_owner_breslow_cox_partial_likelihood_finite_signed_beta_grid_v1",
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  list(artifact = artifact, manifest = list(workload = list(families = list(
    survival_artifacts = list(cox_grid = artifact)))))
}

test_that("Cox partial-likelihood grid validates only its signed candidates", {
  fixture <- .cox_partial_grid_fixture()
  artifact <- dsVertClient:::.dsvert_dp_cox_partial_grid_artifact(
    fixture$manifest, "study", "cox_grid", "server_a",
    "add_remove_patient", 256, 20)
  fit <- dsVertClient:::.dsvert_dp_cox_partial_grid_moment(c(90, 100), artifact)
  expect_identical(fit$selected_candidate, 1L)
  expect_equal(fit$coefficients, c(x = 0.05))
  expect_equal(fit$hazard_ratio, c(x = exp(0.05)))
  tampered <- fixture$manifest
  tampered$workload$families$survival_artifacts$cox_grid$event$levels <-
    c("censor", "other")
  expect_error(dsVertClient:::.dsvert_dp_cox_partial_grid_artifact(
    tampered, "study", "cox_grid", "server_a",
    "add_remove_patient", 256, 20), "descriptor is invalid")
})

test_that("Cox partial-likelihood Synopsis uses one certified signed block", {
  fixture <- .cox_partial_grid_fixture()
  artifact <- fixture$artifact
  block <- list(key = "cox_grid", descriptor = artifact)
  context <- list(
    lattice = list(output_lattice_scale = 256), layout = list(),
    manifest = fixture$manifest, release = list(), adjacency = "add_remove_patient")
  coordinates <- c(90, 100)
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

  released <- dsVertClient:::.dsvert_dp_cox_partial_grid_synopsis_release(
    "study", "cox_grid", datasources = list(server_a = list()),
    .aggregate = function(...) stop("unexpected aggregate", call. = FALSE))
  expect_identical(released$coordinates, coordinates)
  expect_identical(released$moment$selected_candidate, 1L)
  expect_equal(released$moment$coefficients, c(x = 0.05))
})

test_that("Cox partial-grid frontdoor rejects legacy controls", {
  called <- NULL
  testthat::local_mocked_bindings(
    .dsvert_dp_cox_partial_grid_impl = function(...) {
      called <<- list(...)
      structure(list(coefficients = c(x = 0.05), hazard_ratio = c(x = exp(0.05)),
        status = "public_certified_breslow_cox_partial_likelihood_finite_grid"),
        class = c("dsvert_dp_cox_partial_grid", "ds.vertCox", "list"))
    }, .package = "dsVertClient")
  fit <- ds.vertCox(Surv(time, event) ~ x, data = "study",
    analysis_id = "cox-grid", datasources = list(site_a = structure(list(), class = "mock")))
  expect_s3_class(fit, "dsvert_dp_cox_partial_grid")
  expect_identical(called$data_name, "study")
  expect_identical(called$analysis_id, "cox-grid")
  expect_error(ds.vertCox(Surv(time, event) ~ x, data = "study",
    analysis_id = "cox-grid", max_iter = 1L, datasources = list()),
    "no legacy Cox controls")
  profile <- ds.vertCoxProfileNonDisclosive(
    Surv(time, event) ~ x, data = "study", analysis_id = "cox-grid",
    datasources = list(site_a = structure(list(), class = "mock")))
  expect_s3_class(profile, "dsvert_dp_cox_partial_grid")
  expect_error(ds.vertCoxProfileNonDisclosive(
    Surv(time, event) ~ x, data = "study", analysis_id = "cox-grid",
    max_iter = 1L, datasources = list()), "does not accept legacy Cox controls")
})
