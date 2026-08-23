.lmm_artifact_fixture <- function() {
  capacity <- 20
  scale <- 256
  cluster_capacity <- 4
  raw_l1 <- 2 * cluster_capacity + 2 + 3 * scale
  raw_l2 <- sqrt(
    2 + (2 * cluster_capacity - 1)^2 + 2 * scale^2 + (scale + 1)^2)
  artifact <- list(
    version = "bounded-normalized-random-intercept-moments-v1",
    spec_version = "random_intercept_v1", analysis_id = "lmm_a",
    dataset = "protected", owner_peer = "server_a",
    outcome = list(column = "y", lower = -2, upper = 6),
    cluster = list(column = "site", levels = c("a", "b", "c")),
    observation_capacity = capacity,
    max_patients_per_cluster = cluster_capacity,
    numeric_grid_bits = 8L, coordinate_count = 6L,
    coordinate_order = paste(
      "n_then_cluster_count_then_cluster_size_sq_then_quantized_sum_y",
      "then_quantized_sum_y_sq_then_quantized_cluster_mean_sq_v1",
      sep = "_"),
    source_coordinate_scaling =
      "three_counts_then_three_common_lattice_moments_v1",
    repeated_record_policy = paste(
      "clip_finite_outcome_then_mean_once_per_admitted_patient_and",
      "require_one_consistent_public_cluster_level_v1", sep = "_"),
    missingness_policy = paste(
      "missing_or_nonfinite_outcome_or_missing_or_inconsistent_cluster",
      "excludes_the_patient_from_all_six_coordinates_v1", sep = "_"),
    contribution_domain = paste(
      "one_bounded_patient_outcome_and_one_consistent_cluster_level",
      "with_public_cluster_size_cap_v1", sep = "_"),
    quantization_contract = list(
      version = "random-intercept-unit-moment-quantization-v1",
      input_rounding = "nearest_integer_ties_to_even_r_v1",
      sum_y_max_abs_error_normalized = capacity / (2 * scale),
      sum_y_sq_max_abs_error_normalized = capacity / (2 * scale),
      cluster_mean_sq_max_abs_error_normalized =
        3 * capacity / (2 * scale) + capacity / (4 * scale^2)),
    statistic_maximum = c(
      capacity, capacity, capacity * cluster_capacity,
      rep(capacity * scale, 3L)),
    source_raw_l1_sensitivity = raw_l1,
    source_raw_l2_sensitivity = raw_l2,
    natural_l1_sensitivity = raw_l1 / scale,
    natural_l2_sensitivity = raw_l2 / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste(
      "one_patient_changes_three_counts_by_1_1_and_at_most",
      "2C_minus_1_and_three_quantized_moments_by_S_S_and_S_plus_1",
      "with_replace_one_as_two_add_remove_changes_v1", sep = "_"),
    estimation_scope =
      "bounded_random_intercept_method_of_moments_no_fixed_covariates_v1",
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  list(
    artifact = artifact,
    manifest = list(workload = list(families = list(
      gaussian_models = list(artifacts = list(lmm_a = artifact))))))
}

test_that("random-intercept LMM artifacts validate their full signed contract", {
  fixture <- .lmm_artifact_fixture()
  artifact <- dsVertClient:::.dsvert_dp_lmm_artifact(
    fixture$manifest, "protected", "lmm_a", "server_a",
    "add_remove_patient", 256, 20)

  expect_identical(artifact$analysis_id, "lmm_a")
  expect_identical(artifact$coordinate_count, 6L)
  expect_identical(artifact$cluster$levels, c("a", "b", "c"))
  tampered <- fixture$manifest
  tampered$workload$families$gaussian_models$artifacts$lmm_a$
    source_raw_l1_sensitivity <- 1
  expect_error(
    dsVertClient:::.dsvert_dp_lmm_artifact(
      tampered,
      "protected", "lmm_a", "server_a", "add_remove_patient", 256, 20),
    "descriptor is invalid")
})

test_that("random-intercept LMM post-processing only consumes DP coordinates", {
  fixture <- .lmm_artifact_fixture()
  artifact <- dsVertClient:::.dsvert_dp_lmm_artifact(
    fixture$manifest, "protected", "lmm_a", "server_a",
    "add_remove_patient", 256, 20)
  result <- dsVertClient:::.dsvert_dp_lmm_moments(
    c(8, 3, 24, 4.2, 3.1, 2.8), artifact)

  expect_identical(result$status, "ok")
  expect_true(all(is.finite(c(result$coefficient, result$sigma2,
                              result$sigma_b2, result$icc))))
  expect_true(result$icc >= 0 && result$icc <= 1)
})

test_that("random-intercept LMM Synopsis release validates its signed block", {
  fixture <- .lmm_artifact_fixture()
  artifact <- fixture$artifact
  block <- list(key = "lmm_a", descriptor = artifact)
  context <- list(
    lattice = list(output_lattice_scale = 256), layout = list(),
    manifest = fixture$manifest, release = list(), adjacency =
      "add_remove_patient")
  coordinates <- c(8, 3, 24, 4.2, 3.1, 2.8)
  verification <- list(
    integrity_valid = TRUE, authenticity = "session_transport_anchored",
    artifact = artifact, coordinates = coordinates,
    validated_moment = dsVertClient:::.dsvert_dp_lmm_moments(
      coordinates, artifact), cohort_id = "cohort",
    logical_snapshot = list(version = "v1"))
  testthat::local_mocked_bindings(
    .dsvert_dp_datasources = function(value) value,
    .dsvert_dp_synopsis_vector_run = function(...) list(run = TRUE),
    .dsvert_dp_vector_context = function(...) context,
    .dsvert_dp_vector_public_metadata = function(...) list(),
    .dsvert_dp_capsule_single_block = function(...) {
      list(start = 1L, end = 1L, length = 1L)
    },
    .dsvert_dp_vector_block_capacity = function(...) 20,
    .dsvert_dp_lmm_artifact = function(...) artifact,
    .dsvert_dp_capsule_vector_blocks = function(...) list(block),
    .dsvert_dp_capsule_vector_values = function(...) coordinates,
    .dsvert_dp_gaussian_synopsis_certificate_build = function(...) {
      list(certificate_sha256 = strrep("a", 64L))
    },
    ds.validateDPGaussianCertificate = function(...) verification,
    .package = "dsVertClient")

  released <- dsVertClient:::.dsvert_dp_lmm_synopsis_release(
    "protected", "lmm_a", datasources = list(server_a = list()),
    .aggregate = function(...) stop("unexpected aggregate", call. = FALSE))
  expect_identical(released$coordinates, coordinates)
  expect_identical(released$moment$status, "ok")

  coordinates <- c(8, 3, 81, 4.2, 3.1, 2.8)
  expect_error(
    dsVertClient:::.dsvert_dp_lmm_synopsis_release(
      "protected", "lmm_a", datasources = list(server_a = list()),
      .aggregate = function(...) stop("unexpected aggregate", call. = FALSE)),
    "violates its signed bounds")
})
