.nb_grid_artifact_fixture <- function() {
  scale <- 256
  capacity <- 20L
  beta_grid <- list(c(0, 0), c(0, 1))
  theta_grid <- c(1, 2)
  loss_bounds <- dsVertClient:::.dsvert_dp_nb_grid_loss_bounds(
    beta_grid, theta_grid, 8L)
  raw <- ceiling(loss_bounds * scale)
  artifact <- list(
    version = "bounded-negative-binomial-likelihood-grid-v1",
    spec_version = "negative_binomial_grid_v1", analysis_id = "nb_grid",
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
      "theta_grid_then_beta_grid_negative_binomial_log_likelihood_v1",
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
      "contribution_for_every_signed_candidate_v1", sep = "_"),
    statistic_maximum = as.list(capacity * raw),
    source_raw_l1_sensitivity = sum(raw),
    source_raw_l2_sensitivity = sqrt(sum(raw^2)),
    natural_l1_sensitivity = sum(raw) / scale,
    natural_l2_sensitivity = sqrt(sum(raw^2)) / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste(
      "one_patient_changes_one_candidate_loss_by_at_most_its_signed",
      "negative_binomial_loss_bound_v1", sep = "_"),
    estimation_scope = paste(
      "bounded_negative_binomial_fixed_covariates_finite_signed",
      "beta_theta_grid_v1", sep = "_"),
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
