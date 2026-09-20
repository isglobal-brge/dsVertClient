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

# Internal client source layout. Time values and permutation controls stay on
# the time owner; only time-presence validity enters the encrypted source lanes.
.dsvert_dp_cox_cross_source_blocks <- function(artifact, cursor) {
  spec <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)$spec
  if (!identical(spec$family, "cox")) .dsvert_dp_cox_grid_cross_fail()
  blocks <- list()
  descriptors <- c(spec$predictors, setNames(list(spec$event, spec$time),
    c(spec$event$reference, spec$time$reference)))
  for (reference in names(descriptors)) {
    descriptor <- descriptors[[reference]]
    time <- identical(reference, spec$time$reference)
    event <- identical(reference, spec$event$reference)
    for (kind in if (time) "validity" else c("value", "validity")) {
      size <- spec$padded_capacity
      end <- cursor + size - 1
      if (!is.numeric(cursor) || length(cursor) != 1L || !is.finite(cursor) ||
          cursor < 1 || cursor != round(cursor) ||
          end > .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_TRANSPORT_COORDINATES) {
        .dsvert_dp_cox_grid_cross_fail()
      }
      key <- paste(artifact$analysis_id, reference, kind, sep = "::")
      blocks[[key]] <- list(input_family = "cox_grid", analysis_id = artifact$analysis_id,
        reference = reference, variable = descriptor$column, dataset = descriptor$dataset,
        owner_peer = descriptor$owner_peer, kind = kind, outcome = event,
        private_time_validity = time, lower = descriptor$lower, upper = descriptor$upper,
        start = as.integer(cursor), end = as.integer(end), length = as.integer(size),
        fraction_bits = if (kind == "validity" || event) 0L else 50L,
        maximum = if (kind == "validity" || event) 1 else 2^50)
      cursor <- end + 1
    }
  }
  list(blocks = blocks, cursor = cursor)
}


# Owner-first binding relays only the authenticated public route commitment.
# The caller must already have authenticated the artifact and source receipt.
.dsvert_dp_cox_cross_bind <- function(manifest, context, artifact, source_receipt,
    session_id, .remote_context, .aggregate) {
  spec <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)$spec
  peers <- context$designated
  owner <- spec$time$owner_peer
  if (!identical(spec$family, "cox") || spec$observation_capacity > 400 ||
      length(peers) != 2L || anyDuplicated(peers) || !owner %in% peers ||
      !identical(owner, spec$event$owner_peer) || is.null(.remote_context)) {
    .dsvert_dp_cox_grid_cross_fail()
  }
  invoke <- function(peer, routing = "") {
    request <- as.call(c(list(as.name("dsvertDPSynopsisGLMGridCrossDS")),
      .remote_context, list(analysis_id = artifact$analysis_id,
        session_id = session_id, action = "bind", batch = 0,
        routing_receipt_json = routing)))
    response <- .dsvert_fanout_by_site(context$conns[peer], setNames(list(request), peer),
      operation = "Cox owner-first bind", .aggregate = .aggregate)[[peer]]
    value <- .dsvert_joint_dp_client_decode(response, "Cox bind response",
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_RECEIPT_BYTES)
    .dsvert_dp_glm_grid_cross_fields(value, c("bound", "routing_receipt"))
    value
  }
  first <- invoke(owner)
  route <- first$routing_receipt
  expected <- list(version = "dsvert-cox-staged-routing-receipt-v1",
    capsule_id = source_receipt$capsule_id, analysis_id = artifact$analysis_id,
    peer_name = owner, peer_identity_pk = unname(context$pinset[[owner]]),
    semantic_key = .dsvert_dp_capsule_source_hash(list(
      version = "cross-grid-semantic-release-key-v2",
      capsule_id = source_receipt$capsule_id,
      source_contract_sha256 = source_receipt$contract_hash,
      signed_contract = artifact$signed_contract, family = artifact$family,
      version_family = artifact$spec_version,
      mechanism = manifest$workload$capsule_mechanism,
      alignment = spec$alignment, caps = artifact$sensitivity)),
    artifact_sha256 = .dsvert_dp_capsule_source_hash(artifact),
    source_contract_sha256 = source_receipt$contract_hash,
    profile_sha256 = artifact$numeric_certificate$profile_sha256,
    certificate_sha256 = artifact$numeric_certificate$certificate_sha256,
    private_result_exposed = FALSE)
  .dsvert_dp_glm_grid_cross_fields(route,
    c(names(expected), "routing_digest", "stage_plan_digest", "signature"))
  .dsvert_dp_glm_grid_cross_equal(route[names(expected)], expected)
  .dsvert_dp_capsule_source_verify(route, "cross-grid-result", owner, context)
  for (field in c("routing_digest", "stage_plan_digest")) {
    if (!.dsvert_dp_capsule_source_hex(route[[field]]) ||
        identical(route[[field]], strrep("0", 64))) .dsvert_dp_cox_grid_cross_fail()
  }
  evaluator <- setdiff(peers, owner)
  second <- invoke(evaluator, .dsvert_joint_dp_client_json(route))
  .dsvert_dp_glm_grid_cross_equal(second$routing_receipt, route)
  responses <- setNames(lapply(list(first$bound, second$bound),
    .dsvert_joint_dp_client_json), c(owner, evaluator))
  bound <- .dsvert_dp_lmm_cross_receipts(responses, context, artifact, "bound")
  for (field in c("capsule_id", "semantic_key", "source_contract_sha256", "stage_plan_digest")) {
    if (!identical(bound[[field]], route[[field]])) .dsvert_dp_cox_grid_cross_fail()
  }
  bound
}


# Exact projection used by the authenticated single-Cox source contract.
.dsvert_dp_cox_cross_transport_layout <- function(manifest, artifact) {
  release <- .dsvert_dp_capsule_vector_layout(manifest)
  artifacts <- manifest$workload$families$gaussian_models$artifacts
  if (length(artifacts) != 1L || !identical(names(artifacts), artifact$analysis_id) ||
      release$coordinate_count != artifact$coordinate_count + 1L) {
    .dsvert_dp_cox_grid_cross_fail()
  }
  start <- ceiling(release$coordinate_count /
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES) *
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES + 1
  projection <- .dsvert_dp_cox_cross_source_blocks(artifact, start)
  shape <- list(version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_LAYOUT_VERSION,
    capsule_id = manifest$capsule_identity$capsule_id,
    release_coordinate_count = as.integer(release$coordinate_count),
    release_coordinate_order_sha256 = release$sha256,
    private_start = as.integer(start),
    padding_coordinates = as.integer(start - release$coordinate_count - 1),
    transport_coordinate_count = as.integer(projection$cursor - 1),
    blocks = projection$blocks, source_peers = artifact$participating_peers,
    computation_peers = artifact$computation_peers,
    padding_rule = "zero_to_next_source_chunk_boundary_v1",
    payload_rule = "manifest_order_capacity_padded_ring128_value_then_validity_no_exact_release_v1")
  shape$transport_coordinate_order_sha256 <- .dsvert_dp_capsule_source_hash(shape)
  shape
}

# Internal signed entry into the existing durable executor. Public discovery
# remains closed until workload admission and certificate dispatch are wired.
.dsvert_dp_cox_cross_orchestrate <- function(manifest_json, manifest, context,
    source_receipt, policy, schema_manifest, .aggregate, .remote_context) {
  artifacts <- manifest$workload$families$gaussian_models$artifacts
  if (length(artifacts) != 1L) .dsvert_dp_cox_grid_cross_fail()
  artifact <- artifacts[[1L]]
  contract <- .dsvert_dp_cox_grid_cross_contract_validate(
    .dsvert_dp_glm_grid_cross_embedded_contract(artifact), policy, schema_manifest)
  if (contract$spec$observation_capacity > 400 ||
      !identical(names(artifacts), contract$spec$analysis_id) ||
      !setequal(context$designated, unlist(contract$spec$computation_peers))) {
    .dsvert_dp_cox_grid_cross_fail()
  }
  .dsvert_dp_glm_grid_cross_equal(as.list(context$pinset), as.list(policy$peer_pinset))
  .dsvert_dp_glm_grid_cross_equal(artifact, .dsvert_dp_cox_cross_workload_artifact(contract))
  .dsvert_dp_glm_grid_cross_equal(manifest,
    .dsvert_joint_dp_client_decode(manifest_json, "Cox manifest",
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES))
  layout <- .dsvert_dp_cox_cross_transport_layout(manifest, artifact)
  .dsvert_dp_staged_cross_orchestrate(manifest_json, manifest, context,
    source_receipt, artifact, layout, .aggregate, .remote_context)
}
