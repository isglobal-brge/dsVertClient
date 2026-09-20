# Authenticated Step-2 release adapter shared by NB2 and categorical grids.
.dsvert_dp_likelihood_grid_cross_release <- function(contract, policy, schema, datasources, .aggregate) {
  contract <- .dsvert_dp_glm_grid_profile_admit(contract, policy, schema)
  datasources <- .dsvert_dp_datasources(datasources)
  bootstrap <- .dsvert_dp_synopsis_bootstrap_build_v1(datasources, .aggregate = .aggregate)
  trusted <- .dsvert_dp_synopsis_client_bundle(bootstrap$manifest_bundle, bootstrap$status)
  .dsvert_dp_glm_grid_cross_equal(as.list(policy$peer_pinset), as.list(trusted$context$pinset))
  if (!setequal(unname(policy$designated_noise_peers), trusted$context$designated)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  .dsvert_dp_glm_grid_cross_equal(schema,
    .dsvert_joint_dp_client_decode(bootstrap$manifest_bundle$schema_json,
      "signed likelihood-grid schema", .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES))
  spec <- contract$spec
  check <- function(manifest) {
    artifact <- manifest$workload$families$gaussian_models$artifacts[[spec$analysis_id]]
    if (!is.list(artifact)) .dsvert_dp_glm_grid_cross_fail()
    .dsvert_dp_glm_grid_cross_equal(.dsvert_dp_glm_grid_cross_embedded_contract(artifact), contract)
    .dsvert_dp_glm_grid_cross_client_artifact(artifact, spec$dataset,
      spec$analysis_id, NULL, manifest$admission$adjacency,
      2^manifest$bounds$numeric_grid_bits, manifest$admission$unit_capacity, spec$family)
  }
  check(trusted$manifest)
  run <- .dsvert_dp_synopsis_vector_run(datasources, .aggregate = .aggregate,
    .request_check = check)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  artifact <- check(context$manifest)
  blocks <- .dsvert_dp_capsule_vector_blocks(context$layout, "gaussian_models",
    dataset = spec$dataset, owner_peer = artifact$owner_peer)
  blocks <- blocks[vapply(blocks, function(x) identical(x$key, spec$analysis_id), logical(1L))]
  if (length(blocks) != 1L) .dsvert_dp_glm_grid_cross_fail()
  .dsvert_dp_glm_grid_cross_equal(blocks[[1L]]$descriptor,
    context$manifest$workload$families$gaussian_models$artifacts[[spec$analysis_id]])
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, blocks[[1L]])
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(context, artifact, blocks[[1L]], coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  if (!identical(verification$integrity_valid, TRUE) ||
      !identical(verification$authenticity, "session_transport_anchored") ||
      !identical(as.numeric(verification$output_lattice_scale), 2^spec$numeric_grid_bits)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  coordinates <- unname(verification$coordinates * verification$output_lattice_scale)
  result <- if (identical(spec$family, "nb")) {
    .dsvert_dp_nb_grid_cross_moment(coordinates, spec)
  } else {
    .dsvert_dp_categorical_grid_cross_postprocess(contract, coordinates, spec$family)
  }
  result$implementation_state <- artifact$implementation_state
  c(.dsvert_dp_vector_public_metadata(context), result, list(
    provenance_certificate = certificate, certificate_sha256 = certificate$certificate_sha256,
    signed_contract = contract))
}

