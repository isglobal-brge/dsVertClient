.ordinal_grid_artifact_fixture <- function() {
  scale <- 256
  capacity <- 20L
  levels <- c("A", "B", "C")
  candidates <- list(
    list(thresholds = c(-1, 1), beta = c(0, 0)),
    list(thresholds = c(-0.5, 0.5), beta = c(1, 1)))
  loss_bounds <- dsVertClient:::.dsvert_dp_ordinal_grid_loss_bounds(candidates)
  raw <- ceiling(loss_bounds * scale)
  artifact <- list(
    version = "bounded-ordinal-likelihood-grid-v1",
    spec_version = "ordinal_grid_v1", analysis_id = "ordinal_grid",
    dataset = "protected", owner_peer = "server_a",
    outcome = list(column = "y", levels = levels,
                   ordered_levels = levels),
    predictors = list(x = list(column = "x", lower = 0, upper = 10)),
    predictor_order = "x", intercept = TRUE,
    design_terms = c("(Intercept)", "x"), observation_capacity = capacity,
    candidate_grid = candidates,
    candidate_order = "canonical_ordinal_cumulative_logit_grid_v1",
    candidate_loss_bounds = as.list(loss_bounds), numeric_grid_bits = 8L,
    coordinate_count = 2L,
    coordinate_order = paste(
      "canonical_ordinal_candidate_grid_cumulative_logit_negative_log_likelihood_v1",
      sep = "_"),
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    repeated_record_policy = paste(
      "require_one_categorical_outcome_and_mean_once_per_admitted",
      "patient_v1", sep = "_"),
    missingness_policy = paste(
      "missing_outcome_or_missing_or_nonfinite_predictor_excludes_patient",
      "and_unknown_or_conflicting_nonmissing_outcome_rejects_v1", sep = "_"),
    contribution_domain = paste(
      "one_bounded_patient_ordinal_cumulative_logit_negative_log_likelihood",
      "contribution_for_every_signed_candidate_v1", sep = "_"),
    statistic_maximum = as.list(capacity * raw),
    source_raw_l1_sensitivity = sum(raw),
    source_raw_l2_sensitivity = sqrt(sum(raw^2)),
    natural_l1_sensitivity = sum(raw) / scale,
    natural_l2_sensitivity = sqrt(sum(raw^2)) / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste(
      "one_patient_changes_one_candidate_loss_by_at_most_its_signed",
      "ordinal_cumulative_logit_loss_bound_v1", sep = "_"),
    estimation_scope = paste(
      "bounded_ordinal_cumulative_logit_fixed_covariates_finite_signed",
      "candidate_grid_v1", sep = "_"),
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  list(artifact = artifact, manifest = list(workload = list(families = list(
    gaussian_models = list(artifacts = list(ordinal_grid = artifact))))))
}

test_that("ordinal finite grid validates its signed contract and selects a candidate", {
  fixture <- .ordinal_grid_artifact_fixture()
  artifact <- dsVertClient:::.dsvert_dp_ordinal_grid_artifact(
    fixture$manifest, "protected", "ordinal_grid", "server_a",
    "add_remove_patient", 256, 20)
  fit <- dsVertClient:::.dsvert_dp_ordinal_grid_moment(c(100, 90), artifact)

  expect_identical(artifact$coordinate_count, 2L)
  expect_identical(fit$status, "ok")
  expect_identical(fit$selected_candidate, 2L)
  expect_identical(fit$ordered_levels, c("A", "B", "C"))
  expect_equal(fit$coefficients, c(x = 0.1))
  expect_equal(unname(fit$thresholds), c(-1.5, -0.5))
  serialized <- fixture$manifest
  serialized$workload$families$gaussian_models$artifacts$ordinal_grid$
    outcome$levels <- as.list(c("A", "B", "C"))
  expect_identical(dsVertClient:::.dsvert_dp_ordinal_grid_artifact(
    serialized, "protected", "ordinal_grid", "server_a",
    "add_remove_patient", 256, 20)$outcome$ordered_levels,
    c("A", "B", "C"))
  tampered <- fixture$manifest
  tampered$workload$families$gaussian_models$artifacts$ordinal_grid$
    candidate_grid[[2L]]$thresholds <- c(1, 1)
  expect_error(dsVertClient:::.dsvert_dp_ordinal_grid_artifact(
    tampered, "protected", "ordinal_grid", "server_a",
    "add_remove_patient", 256, 20), "ordinal grid descriptor is invalid")
})
