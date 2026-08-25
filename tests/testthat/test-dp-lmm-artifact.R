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

.lmm_fixed_artifact_fixture <- function() {
  capacity <- 20
  scale <- 256
  cluster_capacity <- 4
  dimension <- 2
  summary_count <- dimension * (dimension + 1) / 2 + dimension + 1
  raw_l1 <- (3 + summary_count * (1 + 2 * cluster_capacity^2)) * scale
  raw_l2 <- sqrt(
    (3 + summary_count * (1 + 4 * cluster_capacity^4)) * scale^2)
  artifact <- list(
    version = "bounded-normalized-random-intercept-fixed-sufficient-statistics-v2",
    spec_version = "random_intercept_fixed_v2", analysis_id = "lmm_fixed",
    dataset = "protected", owner_peer = "server_a",
    outcome = list(column = "y", lower = -2, upper = 6),
    cluster = list(column = "site", levels = c("a", "b", "c")),
    predictors = list(x = list(column = "x", lower = 0, upper = 10)),
    predictor_order = "x", intercept = TRUE,
    design_terms = c("(Intercept)", "x"),
    observation_capacity = capacity,
    max_patients_per_cluster = cluster_capacity,
    variance_ratio_grid = c(0, 0.5, 2),
    numeric_grid_bits = 8L,
    coordinate_count = as.integer((cluster_capacity + 1L) *
                                  (summary_count + 1L)),
    coordinate_order = paste(
      "n_then_global_xtx_upper_xty_yty_then_each_cluster_size_from_1",
      "through_C_as_count_xtx_upper_xty_yty_v2", sep = "_"),
    source_coordinate_scaling =
      "counts_left_shifted_to_common_numeric_lattice_v1",
    repeated_record_policy = paste(
      "clip_finite_complete_outcome_predictor_rows_then_mean_once_per",
      "admitted_patient_and_require_one_consistent_public_cluster_level_v2",
      sep = "_"),
    missingness_policy = paste(
      "missing_or_nonfinite_outcome_predictor_or_missing_or_inconsistent",
      "cluster_excludes_the_patient_from_every_LMM_coordinate_v2",
      sep = "_"),
    contribution_domain = paste(
      "one_bounded_patient_vector_and_one_consistent_cluster_with",
      "public_cluster_size_cap_v2", sep = "_"),
    quantization_contract = list(
      version = "random-intercept-fixed-common-lattice-quantization-v1",
      input_rounding = "nearest_integer_ties_to_even_r_v1",
      common_lattice = "numeric_grid_v1"),
    statistic_maximum = c(
      capacity, rep(capacity * scale, summary_count),
      rep(c(capacity, rep(capacity * cluster_capacity * scale,
                         summary_count)), cluster_capacity)),
    source_raw_l1_sensitivity = raw_l1,
    source_raw_l2_sensitivity = raw_l2,
    natural_l1_sensitivity = raw_l1 / scale,
    natural_l2_sensitivity = raw_l2 / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste(
      "one_patient_changes_n_and_at_most_two_cluster_size_bins_and",
      "bounded_squared_grid_cluster_summaries_with_replace_one_as",
      "two_add_remove_changes_v2", sep = "_"),
    estimation_scope = paste(
      "bounded_random_intercept_GLS_fixed_effects_finite_signed",
      "variance_ratio_grid_ML_profile_v1", sep = "_"),
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  list(
    artifact = artifact,
    manifest = list(workload = list(families = list(
      gaussian_models = list(artifacts = list(lmm_fixed = artifact))))))
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

test_that("fixed-effect random-intercept LMM validates and solves GLS coordinates", {
  fixture <- .lmm_fixed_artifact_fixture()
  artifact <- dsVertClient:::.dsvert_dp_lmm_fixed_artifact(
    fixture$manifest, "protected", "lmm_fixed", "server_a",
    "add_remove_patient", 256, 20)
  global <- c(4, 4, 2, 1.68, 2.1, 1.52, 1.43)
  size_one <- c(0, rep(0, 6L))
  size_two <- c(2, 8, 4, 3.28, 4.2, 2.98, 2.81)
  coordinates <- c(global, size_one, size_two, rep(0, 14L))
  fit <- dsVertClient:::.dsvert_dp_lmm_fixed_moments(coordinates, artifact)

  expect_identical(artifact$coordinate_count, 35L)
  expect_identical(artifact$design_terms, c("(Intercept)", "x"))
  expect_identical(fit$status, "ok")
  expect_identical(names(fit$coefficients), c("(Intercept)", "x"))
  expect_true(all(is.finite(c(fit$coefficients, fit$sigma2,
                              fit$sigma_b2, fit$icc))))
  expect_true(fit$icc >= 0 && fit$icc <= 1)
  wire <- fixture$manifest
  wire$workload$families$gaussian_models$artifacts$lmm_fixed$
    design_terms <- list("(Intercept)", "x")
  wire$workload$families$gaussian_models$artifacts$lmm_fixed$
    variance_ratio_grid <- as.list(c(0, 0.5, 2))
  normalized_wire <- dsVertClient:::.dsvert_dp_lmm_fixed_artifact(
    wire, "protected", "lmm_fixed", "server_a", "add_remove_patient", 256,
    20)
  expect_identical(normalized_wire$design_terms, c("(Intercept)", "x"))
  expect_identical(normalized_wire$variance_ratio_grid, c(0, 0.5, 2))
  tampered <- fixture$manifest
  tampered$workload$families$gaussian_models$artifacts$lmm_fixed$
    variance_ratio_grid <- c(0.1, 1)
  expect_error(
    dsVertClient:::.dsvert_dp_lmm_fixed_artifact(
      tampered, "protected", "lmm_fixed", "server_a",
      "add_remove_patient", 256, 20),
    "descriptor is invalid")
})

test_that("LMM frontdoor admits signed fixed-effect artifacts only", {
  artifact <- .lmm_fixed_artifact_fixture()$artifact
  release <- structure(list(
    status = "ok", coefficient = c("(Intercept)" = 1, x = 0.5),
    coefficients = c("(Intercept)" = 1, x = 0.5), sigma2 = 2,
    sigma_b2 = 1, icc = 1 / 3, cluster_count = 2, n_obs = 4,
    signed_artifact = artifact), class = c("ds.vertDPLMM", "list"))
  testthat::local_mocked_bindings(
    ds.vertDPLMM = function(...) release,
    .package = "dsVertClient")
  fit <- ds.vertLMM(
    y ~ x, data = "protected", cluster_col = "site",
    analysis_id = "lmm_fixed", reml = FALSE)

  expect_s3_class(fit, "ds.vertLMM")
  expect_identical(fit$coefficients, c("(Intercept)" = 1, x = 0.5))
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

test_that("historical LMM names only admit the signed moment estimand", {
  release <- structure(list(
    status = "ok", coefficient = c("(Intercept)" = 4.5),
    coefficients = c("(Intercept)" = 4.5), sigma2 = 1,
    sigma_b2 = 0.5, icc = 1 / 3, cluster_count = 3,
    n_obs = 24, signed_artifact = list(
      outcome = list(column = "y"),
      cluster = list(column = "site"),
      estimation_scope =
        "bounded_random_intercept_method_of_moments_no_fixed_covariates_v1")),
    class = c("ds.vertDPLMM", "list"))
  testthat::local_mocked_bindings(
    ds.vertDPLMM = function(...) release,
    .package = "dsVertClient")

  fit <- ds.vertLMM(
    y ~ 1, data = "protected", cluster_col = "site",
    analysis_id = "lmm_a", reml = FALSE)
  expect_s3_class(fit, "ds.vertLMM")
  expect_identical(fit$coefficients, c("(Intercept)" = 4.5))
  expect_identical(fit$reml, FALSE)
  expect_null(fit$cluster_sizes)

  alias <- ds.vert.lmm(
    y ~ 1, data = "protected", cluster_col = "site",
    analysis_id = "lmm_a")
  expect_s3_class(alias, "ds.vertLMM")
  expect_identical(alias$frontdoor, "ds.vert.lmm")
  k3_conns <- stats::setNames(vector("list", 3L), paste0("site", 1:3))
  k3 <- ds.vertLMM.k3(
    y ~ 1, data = "protected", cluster_col = "site",
    analysis_id = "lmm_a", datasources = k3_conns)
  expect_s3_class(k3, "ds.vertLMM")
  expect_identical(k3$frontdoor, "ds.vertLMM.k3")
  expect_error(
    ds.vertLMM.k3(
      y ~ 1, "protected", "site", "lmm_a", datasources = k3_conns[1:2]),
    "at least three DataSHIELD connections")
  expect_error(
    ds.vertLMM(y ~ x, "protected", "site", "lmm_a", reml = FALSE),
    "only an outcome ~ 1 formula")
  expect_error(
    ds.vertLMM(y ~ 1, "protected", "site", "lmm_a", reml = TRUE),
    "reml=FALSE")

  release$signed_artifact$cluster$column <- "other_site"
  expect_error(
    ds.vertLMM(
      y ~ 1, "protected", "site", "lmm_a", reml = FALSE),
    "does not match the signed LMM artifact")
})

test_that("historical GLMM names admit only signed binary moment postprocessing", {
  release <- structure(list(
    status = "ok", coefficient = c("(Intercept)" = 0.65),
    coefficients = c("(Intercept)" = 0.65), sigma2 = 0.18,
    sigma_b2 = 0.04, icc = 0.2, n_obs = 100,
    signed_artifact = list(
      outcome = list(column = "y", lower = 0, upper = 1),
      cluster = list(column = "site"),
      estimation_scope =
        "bounded_random_intercept_method_of_moments_no_fixed_covariates_v1")),
    class = c("ds.vertDPLMM", "list"))
  testthat::local_mocked_bindings(
    ds.vertDPLMM = function(...) release,
    .package = "dsVertClient")

  fit <- ds.vertGLMM(
    y ~ 1, data = "protected", cluster_col = "site",
    analysis_id = "binary_random_intercept")
  expect_s3_class(fit, "ds.vertGLMM")
  expect_identical(fit$family, "binomial")
  expect_identical(
    fit$estimand,
    "binary_random_intercept_population_average_moment_approximation")
  expect_equal(fit$marginal_probability, 0.65)
  expect_equal(fit$coefficients, c("(Intercept)" = stats::qlogis(0.65)))
  expect_equal(fit$icc_observed, 0.2)
  expect_equal(fit$sigma_b2, pi^2 * 0.2 / (3 * 0.8))
  expect_null(fit$standard_errors)
  expect_null(fit$cluster_sizes)
  expect_false(fit$legacy_fallback_called)

  alias <- ds.vert.glmm(
    y ~ 1, data = "protected", cluster_col = "site",
    analysis_id = "binary_random_intercept",
    datasources = list(site_a = list(), site_b = list()))
  expect_s3_class(alias, "ds.vertGLMM")
  expect_identical(alias$frontdoor, "ds.vert.glmm")
  expect_identical(alias$method_frontdoor, "moment")

  expect_error(
    ds.vertGLMM(y ~ x, "protected", "site", "binary_random_intercept"),
    "only an outcome ~ 1 formula")
  expect_error(
    ds.vertGLMM(y ~ 1, "protected", "site", "binary_random_intercept",
                compute_se = TRUE),
    "compute_se=FALSE")

  release$signed_artifact$outcome$upper <- 2
  expect_error(
    ds.vertGLMM(y ~ 1, "protected", "site", "binary_random_intercept"),
    "binary \\[0, 1\\] outcome bounds")

  release$signed_artifact$outcome$upper <- 1
  release$coefficient <- release$coefficients <- c("(Intercept)" = 0)
  projected <- ds.vertGLMM(
    y ~ 1, "protected", "site", "binary_random_intercept")
  expect_equal(projected$marginal_probability, 1 / 200)
  expect_true(projected$probability_projection_applied)
})
