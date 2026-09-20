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

# Internal cold vector reader. As in the shared Synopsis reader, trusted and
# compiled must come from authenticated bundle/compilation validation; neither
# is a public caller override. Source admission and public dispatch stay closed.
.dsvert_dp_cox_cross_read_vector <- function(release_receipts, replay_responses,
    publication_receipts, trusted, compiled, policy, schema_manifest, analysis_id) {
  manifest <- trusted$manifest
  artifact <- manifest$workload$families$gaussian_models$artifacts[[analysis_id]]
  contract <- .dsvert_dp_cox_grid_cross_contract_validate(
    .dsvert_dp_glm_grid_cross_embedded_contract(artifact), policy, schema_manifest)
  spec <- contract$spec
  if (spec$observation_capacity > 400 || !identical(spec$analysis_id, analysis_id)) {
    .dsvert_dp_cox_grid_cross_fail()
  }
  .dsvert_dp_glm_grid_cross_equal(artifact, .dsvert_dp_cox_cross_workload_artifact(contract))
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  .dsvert_dp_glm_grid_cross_equal(compiled$layout, layout)
  .dsvert_dp_glm_grid_cross_equal(compiled$lattice,
    .dsvert_dp_synopsis_client_lattice(manifest, layout))
  block <- .dsvert_dp_capsule_single_block(layout, "gaussian_models",
    dataset = spec$dataset, owner_peer = spec$owner_peer,
    predicate = function(value) identical(value$key, analysis_id))
  positions <- seq.int(block$start, block$end)
  if (length(positions) != artifact$coordinate_count ||
      any(compiled$lattice$scale_shifts[positions] != 0) ||
      compiled$lattice$output_lattice_bits != spec$numeric_grid_bits) {
    .dsvert_dp_cox_grid_cross_fail()
  }
  execution <- .dsvert_dp_synopsis_client_execution(compiled)
  releases <- .dsvert_dp_synopsis_client_release_set(
    release_receipts, compiled, execution, trusted)
  replay <- .dsvert_dp_synopsis_client_replay(
    replay_responses, releases, compiled, execution, trusted)
  evidence <- .dsvert_dp_cox_cross_public_evidence_set(publication_receipts,
    trusted$context, manifest, analysis_id, releases$reference, compiled, policy, schema_manifest)
  # REPLAY has checked exact nonnegative integers <= 2^53-1 and signed caps.
  # Cox consumes that integer lattice directly, with no second scaling/rounding.
  coordinates <- as.numeric(replay$scaled[positions])
  .dsvert_dp_cox_grid_cross_moment(coordinates, spec)
  list(contract = contract, policy = policy, schema_manifest = schema_manifest,
    coordinates = coordinates, publication_evidence = evidence,
    release_receipts = releases$receipts, replay_responses = replay$replay)
}
