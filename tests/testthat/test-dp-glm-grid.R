.glm_grid_artifact_fixture <- function(family = c("binomial", "poisson")) {
  family <- match.arg(family)
  scale <- 256
  capacity <- 20L
  beta_grid <- list(c(0, 0), c(0, 1))
  maximum <- if (identical(family, "poisson")) 8L else NULL
  loss_bounds <- dsVertClient:::.dsvert_dp_glm_grid_loss_bounds(
    family, beta_grid, maximum)
  raw <- ceiling(loss_bounds * scale)
  artifact <- list(
    version = dsVertClient:::.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS[[family]],
    spec_version = paste0(family, "_grid_v1"),
    analysis_id = paste0(family, "_grid"), dataset = "protected",
    owner_peer = "server_a",
    outcome = list(column = "y", lower = 0,
                   upper = if (identical(family, "poisson")) 8 else 1),
    predictors = list(x = list(column = "x", lower = 0, upper = 10)),
    predictor_order = "x", intercept = TRUE,
    design_terms = c("(Intercept)", "x"), observation_capacity = capacity,
    max_outcome = maximum, beta_grid = beta_grid,
    candidate_order = "canonical_beta_grid_glm_v1",
    candidate_loss_bounds = as.list(loss_bounds), numeric_grid_bits = 8L,
    coordinate_count = 2L,
    coordinate_order = paste(
      "canonical_beta_grid", family, "negative_log_likelihood_v1", sep = "_"),
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    repeated_record_policy = paste(
      "require_one_bounded", family, "outcome_and_mean_once_per_admitted",
      "patient_v1", sep = "_"),
    missingness_policy = paste(
      "noninteger_or_out_of_range_or_missing_outcome_or_missing_or",
      "nonfinite_predictor_excludes_patient_v1", sep = "_"),
    contribution_domain = paste(
      "one_bounded_patient", family, "negative_log_likelihood",
      "contribution_for_every_signed_candidate_v1", sep = "_"),
    statistic_maximum = as.list(capacity * raw),
    source_raw_l1_sensitivity = sum(raw),
    source_raw_l2_sensitivity = sqrt(sum(raw^2)),
    natural_l1_sensitivity = sum(raw) / scale,
    natural_l2_sensitivity = sqrt(sum(raw^2)) / scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste(
      "one_patient_changes_one_candidate_loss_by_at_most_its_signed",
      family, "loss_bound_v1", sep = "_"),
    estimation_scope = paste(
      "bounded", family, "fixed_covariates_finite_signed_beta_grid_v1",
      sep = "_"),
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  list(artifact = artifact, manifest = list(workload = list(families = list(
    gaussian_models = list(artifacts = stats::setNames(
      list(artifact), artifact$analysis_id))))))
}

test_that("finite binomial and Poisson grids validate their signed contracts", {
  for (family in c("binomial", "poisson")) {
    fixture <- .glm_grid_artifact_fixture(family)
    analysis_id <- fixture$artifact$analysis_id
    artifact <- dsVertClient:::.dsvert_dp_glm_grid_artifact(
      fixture$manifest, "protected", analysis_id, "server_a",
      "add_remove_patient", 256, 20, family)
    fit <- dsVertClient:::.dsvert_dp_glm_grid_moment(c(100, 90), artifact,
                                                      family)
    expect_identical(artifact$coordinate_count, 2L)
    expect_identical(fit$status, "ok")
    expect_identical(fit$selected_candidate, 2L)
    expect_equal(fit$coefficients, c("(Intercept)" = 0, x = 0.1))
    tampered <- fixture$manifest
    tampered$workload$families$gaussian_models$artifacts[[analysis_id]]$
      candidate_loss_bounds[[1L]] <- 0
    expect_error(dsVertClient:::.dsvert_dp_glm_grid_artifact(
      tampered, "protected", analysis_id, "server_a",
      "add_remove_patient", 256, 20, family),
      "finite GLM grid descriptor is invalid")
  }
})

test_that("ds.vertGLM admits only a signed finite binomial or Poisson grid", {
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_dp_glm_grid_impl = function(formula, data_name, analysis_id,
                                        family, server = NULL,
                                        datasources = NULL, .aggregate) {
      seen <<- list(formula = formula, data_name = data_name,
                    analysis_id = analysis_id, family = family)
      structure(list(coefficients = c("(Intercept)" = 0, x = 0.1),
                     production_ready = FALSE,
                     source_values_exposed = FALSE),
                class = c("dsvert_dp_glm_grid", "ds.glm", "list"))
    },
    .dsvert_federation_argument = function(data, datasources) {
      list(value = data, datasources = datasources)
    },
    .package = "dsVertClient")
  fit <- ds.vertGLM(
    y ~ x, data = "protected", family = "binomial",
    analysis_id = "binomial_grid", datasources = list(site_a = list()))
  expect_s3_class(fit, "dsvert_dp_glm_grid")
  expect_identical(seen$family, "binomial")
  expect_identical(seen$analysis_id, "binomial_grid")
  expect_error(ds.vertGLM(
    y ~ x, data = "protected", family = "binomial",
    analysis_id = "binomial_grid", lambda = 1),
    "does not accept legacy controls")
  expect_error(ds.vertGLM(
    y ~ x, data = "protected", family = "gaussian",
    analysis_id = "binomial_grid"),
    "binomial or Poisson")
})
