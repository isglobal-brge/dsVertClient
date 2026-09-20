# Internal Cox projection for authenticated source/lifecycle wiring.
# Public workload discovery stays closed until the complete staged path is proved.
.dsvert_dp_cox_cross_workload_artifact <- function(contract) {
  spec <- contract$spec
  artifact <- contract$artifact
  artifact$transcript$producer <- artifact$transcript$operation
  artifact$transcript$operation <- NULL
  c(artifact, list(dataset = spec$dataset, family = spec$family,
    time = spec$time[c("column", "dataset", "owner_peer", "lower", "upper")],
    event = spec$event[c("column", "dataset", "owner_peer", "lower", "upper")],
    predictors = lapply(spec$predictors, function(x)
      x[c("column", "dataset", "owner_peer", "lower", "upper")]),
    predictor_order = unlist(spec$predictor_order, use.names = FALSE),
    input_variable_order = unlist(spec$input_variable_order, use.names = FALSE),
    design_terms = unlist(spec$design_terms, use.names = FALSE),
    beta_grid = spec$beta_grid, intercept = FALSE,
    observation_capacity = spec$observation_capacity, padded_capacity = spec$padded_capacity,
    ties = spec$ties, time_semantics = spec$time_semantics,
    predictor_normalization = spec$predictor_normalization,
    complete_case = spec$complete_case,
    statistic_maximum = spec$sensitivity$maximum_coordinates,
    source_raw_l1_sensitivity = spec$sensitivity$raw_l1_sensitivity,
    source_raw_l2_sensitivity = spec$sensitivity$raw_l2_sensitivity,
    natural_l1_sensitivity = spec$sensitivity$natural_l1_sensitivity,
    natural_l2_sensitivity = spec$sensitivity$natural_l2_sensitivity,
    numeric_certificate = spec$numeric_contract, adjacency = spec$adjacency,
    signed_contract = .dsvert_joint_dp_client_json(.dsvert_joint_dp_client_canonical(contract))))
}

# Internal verifier for the two server publication receipts. It consumes an
# already authenticated DP release/compilation, including on a cold read; it
# neither authenticates the DP vector itself nor enables public Cox dispatch.
.dsvert_dp_cox_cross_public_evidence_set <- function(responses, context, manifest,
    analysis_id, release, compiled, policy, schema_manifest) {
  artifact <- manifest$workload$families$gaussian_models$artifacts[[analysis_id]]
  if (!is.list(artifact) || !identical(artifact$family, "cox") ||
      !identical(artifact$version, .DSVERT_CLIENT_DP_COX_GRID_CROSS_ARTIFACT_VERSION) ||
      !identical(artifact$spec_version, .DSVERT_CLIENT_DP_COX_GRID_CROSS_SPEC_VERSION) ||
      !identical(artifact$analysis_id, analysis_id)) .dsvert_dp_cox_grid_cross_fail()
  contract <- .dsvert_dp_cox_grid_cross_contract_validate(
    .dsvert_dp_glm_grid_cross_embedded_contract(artifact), policy, schema_manifest)
  if (contract$spec$observation_capacity > 400) .dsvert_dp_cox_grid_cross_fail()
  .dsvert_dp_glm_grid_cross_equal(artifact,
    .dsvert_dp_cox_cross_workload_artifact(contract))
  .dsvert_dp_glm_grid_cross_equal(as.list(context$pinset), as.list(policy$peer_pinset))
  if (!identical(sort(context$designated, method = "radix"),
      sort(unlist(contract$spec$computation_peers, use.names = FALSE), method = "radix"))) {
    .dsvert_dp_cox_grid_cross_fail()
  }
  .dsvert_dp_staged_cross_public_evidence_set(
    responses, context, manifest, artifact, release, compiled)
}
