test_that("NB2 signed loss caps match the complete negative-binomial likelihood", {
  beta_grid <- list(c(0, 0), c(0, 1), c(-8, 0), c(8, -8, 8))
  theta_grid <- c(.Machine$double.xmin, 1e-6, 0.25, 1, 2, 64)
  for (max_outcome in c(1L, 2L, 8L, 1024L)) {
    expected <- unlist(lapply(theta_grid, function(theta) {
      vapply(beta_grid, function(beta) {
        eta_bound <- sum(abs(beta))
        max(vapply(c(-eta_bound, eta_bound), function(eta) {
          max(-stats::dnbinom(0:max_outcome, size = theta,
                             mu = exp(eta), log = TRUE))
        }, numeric(1L)))
      }, numeric(1L))
    }), use.names = FALSE)
    actual <- dsVertClient:::.dsvert_dp_nb_grid_loss_bounds(
      beta_grid, theta_grid, max_outcome)
    expect_equal(actual, expected, tolerance = 1e-11)
  }
  expect_equal(dsVertClient:::.dsvert_dp_nb_grid_loss_bounds(
    list(c(0, 0)), 2, 1L),
    -stats::dnbinom(1, size = 2, mu = 1, log = TRUE), tolerance = 1e-12)
})

.nb_grid_artifact_fixture <- function() {
  scale <- 256
  capacity <- 20L
  beta_grid <- list(c(0, 0), c(0, 1))
  theta_grid <- c(1, 2)
  loss_bounds <- dsVertClient:::.dsvert_dp_nb_grid_loss_bounds(
    beta_grid, theta_grid, 8L)
  raw <- ceiling(loss_bounds * scale)
  artifact <- list(
    version = "bounded-negative-binomial-likelihood-grid-v2",
    spec_version = "negative_binomial_grid_v2", analysis_id = "nb_grid",
    dataset = "protected", owner_peer = "server_a",
    outcome = list(column = "y", lower = 0, upper = 8),
    predictors = list(x = list(column = "x", lower = 0, upper = 10)),
    predictor_order = "x", intercept = TRUE,
    design_terms = c("(Intercept)", "x"), observation_capacity = capacity,
    max_outcome = 8L, beta_grid = beta_grid, theta_grid = theta_grid,
    candidate_order = "theta_grid_then_beta_grid_v1",
    candidate_loss_bounds = as.list(loss_bounds), numeric_grid_bits = 8L,
    coordinate_count = 4L,
    coordinate_order = paste(
      "theta_grid_then_beta_grid_negative_binomial_log_likelihood_v2",
      sep = "_"),
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    repeated_record_policy = paste(
      "require_one_bounded_count_outcome_and_mean_once_per_admitted",
      "patient_v1", sep = "_"),
    missingness_policy = paste(
      "noninteger_or_out_of_range_or_missing_outcome_or_missing_or",
      "nonfinite_predictor_excludes_patient_v1", sep = "_"),
    contribution_domain = paste(
      "one_bounded_patient_negative_binomial_log_likelihood",
      "contribution_for_every_signed_candidate_v2", sep = "_"),
    statistic_maximum = as.list(capacity * raw),
    source_raw_l1_sensitivity = sum(raw),
    source_raw_l2_sensitivity = sqrt(sum(raw^2)),
    natural_l1_sensitivity = sum(raw) / scale,
    natural_l2_sensitivity = sqrt(sum(raw^2)) / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste(
      "one_patient_changes_one_candidate_loss_by_at_most_its_signed",
      "negative_binomial_loss_bound_v2", sep = "_"),
    estimation_scope = paste(
      "bounded_negative_binomial_fixed_covariates_finite_signed",
      "beta_theta_grid_v2", sep = "_"),
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  list(artifact = artifact, manifest = list(workload = list(families = list(
    gaussian_models = list(artifacts = list(nb_grid = artifact))))))
}

test_that("NB2 finite grid validates its signed contract and selects a candidate", {
  fixture <- .nb_grid_artifact_fixture()
  artifact <- dsVertClient:::.dsvert_dp_nb_grid_artifact(
    fixture$manifest, "protected", "nb_grid", "server_a",
    "add_remove_patient", 256, 20)
  fit <- dsVertClient:::.dsvert_dp_nb_grid_moment(
    c(100, 90, 95, 80), artifact)

  expect_identical(artifact$coordinate_count, 4L)
  expect_identical(fit$status, "ok")
  expect_identical(fit$selected_candidate, 4L)
  expect_equal(fit$coefficients, c(`(Intercept)` = 0, x = 0.1))
  expect_equal(fit$theta, 2)
  printable <- c(fit, list(selected_dp_negative_log_likelihood =
    fit$selected_dp_negative_log_likelihood))
  class(printable) <- c("dsvert_dp_nb2_grid", "ds.vertNBFullRegTheta", "list")
  expect_match(paste(capture.output(print(printable)), collapse = "\n"),
               "finite-grid", fixed = TRUE)
  tampered <- fixture$manifest
  tampered$workload$families$gaussian_models$artifacts$nb_grid$
    candidate_loss_bounds[[1L]] <- 0
  expect_error(dsVertClient:::.dsvert_dp_nb_grid_artifact(
    tampered, "protected", "nb_grid", "server_a",
    "add_remove_patient", 256, 20), "descriptor is invalid")
})

test_that("NB2 v2 caps bind the corrected lattice sensitivities", {
  fixture <- .nb_grid_artifact_fixture()
  artifact <- dsVertClient:::.dsvert_dp_nb_grid_artifact(
    fixture$manifest, "protected", "nb_grid", "server_a",
    "add_remove_patient", 256, 20)
  raw <- c(1598, 2770, 1896, 3338)
  expect_equal(ceiling(artifact$candidate_loss_bounds * 256), raw)
  expect_equal(artifact$statistic_maximum, 20 * raw)
  expect_equal(artifact$source_raw_l1_sensitivity, sum(raw))
  expect_equal(artifact$source_raw_l2_sensitivity, sqrt(sum(raw^2)))
  expect_equal(artifact$natural_l1_sensitivity, sum(raw) / 256)
  expect_equal(artifact$natural_l2_sensitivity, sqrt(sum(raw^2)) / 256)

  replacement <- fixture$manifest
  signed <- replacement$workload$families$gaussian_models$artifacts$nb_grid
  signed$adjacency <- "replace_one_fixed_cohort"
  for (field in c("source_raw_l1_sensitivity", "source_raw_l2_sensitivity",
                  "natural_l1_sensitivity", "natural_l2_sensitivity")) {
    signed[[field]] <- 2 * signed[[field]]
  }
  replacement$workload$families$gaussian_models$artifacts$nb_grid <- signed
  expect_equal(dsVertClient:::.dsvert_dp_nb_grid_artifact(
    replacement, "protected", "nb_grid", "server_a",
    "replace_one_fixed_cohort", 256, 20)$source_raw_l1_sensitivity,
    2 * sum(raw))
})

test_that("NB2 sealed v1 specifications and mixed descriptors fail closed", {
  fixture <- .nb_grid_artifact_fixture()
  for (field in c("version", "spec_version", "coordinate_order",
                  "contribution_domain", "adjacency_sensitivity_basis",
                  "estimation_scope")) {
    legacy <- fixture$manifest
    signed <- legacy$workload$families$gaussian_models$artifacts$nb_grid
    signed[[field]] <- sub("v2$", "v1", signed[[field]])
    legacy$workload$families$gaussian_models$artifacts$nb_grid <- signed
    expect_error(dsVertClient:::.dsvert_dp_nb_grid_artifact(
      legacy, "protected", "nb_grid", "server_a",
      "add_remove_patient", 256, 20), "valid|invalid")
  }

  fragments <- list(describe = list(), survival = list(), gaussian = list(
    nb_grid = list(version = "negative_binomial_grid_v2",
      dataset = "protected", outcome = "y", predictors = "x",
      intercept = TRUE, max_outcome = 8L,
      beta_grid = list(c(0, 0), c(0, 1)), theta_grid = c(1, 2))),
    vertical_cross = list())
  expect_identical(dsVertClient:::.dsvert_dp_capsule_manifest_fragments(
    fragments)$gaussian$nb_grid$version, "negative_binomial_grid_v2")
  fragments$gaussian$nb_grid$version <- "negative_binomial_grid_v1"
  expect_error(dsVertClient:::.dsvert_dp_capsule_manifest_fragments(fragments),
               "invalid custodian workload specification")
  fragments$gaussian$nb_grid[c("max_outcome", "beta_grid", "theta_grid")] <- NULL
  expect_error(dsVertClient:::.dsvert_dp_capsule_manifest_fragments(fragments),
               "invalid custodian workload specification")
})
