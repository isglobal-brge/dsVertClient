.multinom_grid_artifact_fixture <- function() {
  scale <- 256
  capacity <- 20L
  levels <- c("A", "B", "C")
  reference <- "A"
  beta_grid <- list(c(0, 0, 0, 0), c(0, 1, 0, -1))
  loss_bounds <- dsVertClient:::.dsvert_dp_multinom_grid_loss_bounds(
    beta_grid, levels, reference)
  raw <- ceiling(loss_bounds * scale)
  artifact <- list(
    version = "bounded-multinomial-likelihood-grid-v1",
    spec_version = "multinomial_grid_v1", analysis_id = "multinom_grid",
    dataset = "protected", owner_peer = "server_a",
    outcome = list(column = "y", levels = levels, reference = reference),
    predictors = list(x = list(column = "x", lower = 0, upper = 10)),
    predictor_order = "x", intercept = TRUE,
    design_terms = c("(Intercept)", "x"), observation_capacity = capacity,
    beta_grid = beta_grid, candidate_order = "canonical_beta_grid_softmax_v1",
    candidate_loss_bounds = as.list(loss_bounds), numeric_grid_bits = 8L,
    coordinate_count = 2L,
    coordinate_order = paste(
      "canonical_beta_grid_multinomial_softmax_negative_log_likelihood_v1",
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
      "one_bounded_patient_multinomial_softmax_negative_log_likelihood",
      "contribution_for_every_signed_candidate_v1", sep = "_"),
    statistic_maximum = as.list(capacity * raw),
    source_raw_l1_sensitivity = sum(raw),
    source_raw_l2_sensitivity = sqrt(sum(raw^2)),
    natural_l1_sensitivity = sum(raw) / scale,
    natural_l2_sensitivity = sqrt(sum(raw^2)) / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste(
      "one_patient_changes_one_candidate_loss_by_at_most_its_signed",
      "multinomial_softmax_loss_bound_v1", sep = "_"),
    estimation_scope = paste(
      "bounded_multinomial_fixed_covariates_finite_signed_beta_grid_v1",
      sep = "_"),
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  list(artifact = artifact, manifest = list(workload = list(families = list(
    gaussian_models = list(artifacts = list(multinom_grid = artifact))))))
}

test_that("multinomial finite grid validates its signed contract and selects a candidate", {
  fixture <- .multinom_grid_artifact_fixture()
  artifact <- dsVertClient:::.dsvert_dp_multinom_grid_artifact(
    fixture$manifest, "protected", "multinom_grid", "server_a",
    "add_remove_patient", 256, 20)
  fit <- dsVertClient:::.dsvert_dp_multinom_grid_moment(c(100, 90), artifact)

  expect_identical(artifact$coordinate_count, 2L)
  expect_identical(fit$status, "ok")
  expect_identical(fit$selected_candidate, 2L)
  expect_identical(fit$reference, "A")
  expect_equal(fit$coefficients,
               structure(c(0, 0.1, 0, -0.1), dim = c(2L, 2L),
                         dimnames = list(c("(Intercept)", "x"), c("B", "C"))))
  printable <- c(fit, list(selected_candidate = fit$selected_candidate))
  class(printable) <- c("dsvert_dp_multinom_grid", "ds.vertMultinom", "list")
  expect_match(paste(capture.output(print(printable)), collapse = "\n"),
               "finite-grid", fixed = TRUE)
  serialized <- fixture$manifest
  serialized$workload$families$gaussian_models$artifacts$multinom_grid$
    outcome$levels <- as.list(c("A", "B", "C"))
  expect_identical(dsVertClient:::.dsvert_dp_multinom_grid_artifact(
    serialized, "protected", "multinom_grid", "server_a",
    "add_remove_patient", 256, 20)$outcome$levels, c("A", "B", "C"))
  tampered <- fixture$manifest
  tampered$workload$families$gaussian_models$artifacts$multinom_grid$
    outcome$reference <- "not-a-level"
  expect_error(dsVertClient:::.dsvert_dp_multinom_grid_artifact(
    tampered, "protected", "multinom_grid", "server_a",
    "add_remove_patient", 256, 20), "outcome domain is invalid")
})
