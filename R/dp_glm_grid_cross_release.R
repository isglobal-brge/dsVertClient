.dsvert_dp_synopsis_supported_glm_grid_cross_v1 <- function(manifest) {
  artifacts <- .dsvert_dp_glm_grid_cross_artifacts(manifest)
  length(artifacts) == 1L && length(manifest$workload$families$gaussian_models$artifacts) == 1L && all(vapply(artifacts, function(artifact) {
    identical(artifact$implementation_state, "cross_owner_exact_gc_materialized") &&
      identical(artifact$cross_owner_state, "exact_gc_to_joint_dp_vector_v1")
  }, logical(1L)))
}

.dsvert_dp_glm_grid_cross_preflight <- function(manifest, context, schema_json) {
  artifacts <- .dsvert_dp_glm_grid_cross_artifacts(manifest)
  if (!length(artifacts)) return(invisible(TRUE))
  .dsvert_dp_glm_grid_cross_noise_policy(manifest)
  .dsvert_dp_glm_grid_cross_integer(manifest$workload$coordinate_count, 1, 51)
  schema <- .dsvert_joint_dp_client_decode(schema_json, "signed grid schema",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES)
  policy <- list(peer_pinset = context$pinset,
    peer_pinset_sha256 = .dsvert_dp_capsule_source_hash(as.list(context$pinset)),
    designated_noise_peers = context$designated,
    unit_capacity = manifest$admission$unit_capacity,
    numeric_grid_bits = manifest$bounds$numeric_grid_bits,
    adjacency = manifest$admission$adjacency)
  for (artifact in artifacts) {
    contract <- .dsvert_dp_glm_grid_cross_embedded_contract(artifact)
    .dsvert_dp_glm_grid_profile_admit(contract, policy, schema)
    expected <- if (identical(contract$spec$family, "lmm")) {
      .dsvert_dp_grouped_cross_workload_artifact(contract)
    } else .dsvert_dp_glm_grid_cross_workload_artifact(contract)
    .dsvert_dp_glm_grid_cross_equal(artifact, expected)
  }
  invisible(TRUE)
}

.dsvert_dp_glm_grid_cross_receipts <- function(responses, context, artifact, phase) {
  peers <- context$designated
  if (!setequal(names(responses), peers)) .dsvert_dp_glm_grid_cross_fail()
  results <- lapply(peers, function(peer) {
    value <- .dsvert_joint_dp_client_decode(responses[[peer]], "cross-grid receipt",
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_RECEIPT_BYTES)
    fields <- c("version", "phase", "capsule_id", "analysis_id", "peer_name",
      "peer_identity_pk", "semantic_key", "artifact_sha256", "source_contract_sha256",
      "profile_sha256", "certificate_sha256", "private_result_exposed", "signature")
    extra <- switch(phase, bound = "batch_count",
      prepared = c("operation_id", "source_key", "output_key", "purpose", "operation", "vector_len", "batch", "persisted"),
      batch_persisted = "batch", complete = c("coordinate_count", "implementation_state", "cross_owner_state"),
      .dsvert_dp_glm_grid_cross_fail())
    .dsvert_dp_glm_grid_cross_fields(value, c(fields, extra))
    .dsvert_dp_capsule_source_verify(value, "cross-grid-result", peer, context)
    if (identical(phase, "prepared") &&
        (!is.logical(value$persisted) || length(value$persisted) != 1L ||
         is.na(value$persisted))) .dsvert_dp_glm_grid_cross_fail()
    if (!identical(value$version, "dsvert-cross-grid-receipt-v2") ||
        !identical(value$phase, phase) ||
        !identical(value$artifact_sha256, .dsvert_dp_capsule_source_hash(artifact)) ||
        !identical(value$profile_sha256, artifact$numeric_certificate$profile_sha256) ||
        !identical(value$certificate_sha256, artifact$numeric_certificate$certificate_sha256) ||
        !identical(value$analysis_id, artifact$analysis_id) ||
        !identical(value$private_result_exposed, FALSE)) .dsvert_dp_glm_grid_cross_fail()
    value
  })
  common <- function(value) value[setdiff(names(value),
    c("peer_name", "peer_identity_pk", "signature",
      if (identical(phase, "prepared")) "persisted"))]
  .dsvert_dp_glm_grid_cross_equal(common(results[[1]]), common(results[[2]]))
  result <- results[[1]]
  # This local decision is derived only after both signatures and all common
  # bindings were checked. A unilateral durable write still reruns the batch.
  if (identical(phase, "prepared")) attr(result, "both_persisted") <-
    all(vapply(results, function(x) isTRUE(x$persisted), logical(1L)))
  result
}

.dsvert_dp_glm_grid_cross_orchestrate <- function(manifest_json, manifest, context,
    source_receipt, .aggregate, .remote_context) {
  artifacts <- .dsvert_dp_glm_grid_cross_artifacts(manifest)
  if (!length(artifacts)) return(NULL)
  if (any(vapply(artifacts, function(artifact) identical(artifact$family, "lmm"), logical(1L)))) {
    return(.dsvert_dp_lmm_cross_orchestrate(manifest_json, manifest, context,
      source_receipt, .aggregate, .remote_context))
  }
  if (is.null(.remote_context)) .dsvert_dp_glm_grid_cross_fail()
  layout <- .dsvert_dp_gaussian_cross_layout_client(manifest)
  peers <- context$designated
  if (!identical(source_receipt$purpose, .DSVERT_CLIENT_DP_GLM_GRID_CROSS_SOURCE_PURPOSE) ||
      !identical(as.numeric(source_receipt$coordinate_count), as.numeric(layout$transport_coordinate_count)) ||
      !identical(as.numeric(source_receipt$release_coordinate_count), as.numeric(layout$release_coordinate_count)) ||
      !identical(source_receipt$private_layout_sha256, layout$transport_coordinate_order_sha256) ||
      !identical(source_receipt$sampler_handoff_ready, FALSE) ||
      !identical(source_receipt$payload_exposed, FALSE)) .dsvert_dp_glm_grid_cross_fail()
  session_id <- .dsvert_uuid4()
  # A dedicated source-computation session needs the same authenticated
  # Synopsis authority before its exact transport can bind the pinned pair.
  prepare_calls <- stats::setNames(lapply(peers, function(peer) {
    as.call(c(list(as.name("dsvertDPSynopsisPrepareDS")),
      list(session_id = session_id), .remote_context))
  }), peers)
  .dsvert_dp_synopsis_runner_json_set(.dsvert_fanout_by_site(context$conns,
    prepare_calls, operation = "cross-grid Synopsis authority", .aggregate = .aggregate),
    peers, "PREPARE", .DSVERT_CLIENT_SYNOPSIS_PREPARE_MAX_OBJECT_BYTES)
  setup <- .dsvert_dp_cross_exact_setup(.dsvert_setup_exact_gc_transport,
    context$all_conns, context$servers, match(peers, context$servers),
    session_id, .aggregate = .aggregate)
  on.exit(.dsvert_dp_cross_exact_cleanup(context$conns, session_id, setup,
    .aggregate, .dsvert_setup_exact_gc_transport), add = TRUE)
  .dsvert_dp_alignment_mask_run(manifest_json, context, layout, source_receipt,
    session_id, .aggregate, .remote_context = .remote_context)
  completed <- list()
  for (id in names(artifacts)) {
    artifact <- artifacts[[id]]
    invoke <- function(action, batch = 0) {
      calls <- stats::setNames(lapply(peers, function(peer) {
        as.call(c(list(as.name("dsvertDPSynopsisGLMGridCrossDS")),
          .remote_context, list(analysis_id = id, session_id = session_id,
            action = action, batch = batch)))
      }), peers)
      .dsvert_fanout_by_site(context$conns, calls,
        operation = paste("cross-grid", action), .aggregate = .aggregate)
    }
    binding <- .dsvert_dp_glm_grid_cross_receipts(invoke("bind"), context, artifact, "bound")
    if (!identical(binding$source_contract_sha256, source_receipt$contract_hash) ||
        !identical(binding$capsule_id, source_receipt$capsule_id) ||
        !isTRUE(all.equal(binding$batch_count, ceiling(artifact$observation_capacity / artifact$transcript$row_batch_size) *
          ceiling(artifact$coordinate_count / 8)))) .dsvert_dp_glm_grid_cross_fail()
    check_binding <- function(receipt) {
      for (field in c("source_contract_sha256", "capsule_id", "semantic_key")) {
        if (!identical(receipt[[field]], binding[[field]])) .dsvert_dp_glm_grid_cross_fail()
      }
      receipt
    }
    for (batch in seq_len(binding$batch_count)) {
      stage <- .dsvert_dp_glm_grid_cross_receipts(invoke("prepare", batch),
        context, artifact, "prepared")
      check_binding(stage)
      if (isTRUE(attr(stage, "both_persisted"))) next
      initialized <- invoke("start", batch)
      .dsvert_exact_gc_run(context$all_conns, server_names = context$servers,
        servers = match(peers, context$servers), session_id = session_id,
        operation_id = stage$operation_id, source_key = stage$source_key,
        output_key = stage$output_key, operation = stage$operation,
        ring = 128L, frac_bits = 0L, vector_len = stage$vector_len,
        purpose = stage$purpose, transport_ready = TRUE, initialized = initialized, .aggregate = .aggregate)
      check_binding(.dsvert_dp_glm_grid_cross_receipts(invoke("store", batch), context,
        artifact, "batch_persisted"))
    }
    completed[[id]] <- check_binding(.dsvert_dp_glm_grid_cross_receipts(invoke("finalize"),
      context, artifact, "complete"))
  }
  list(enabled = TRUE, sampler_handoff_ready = TRUE,
    private_result_exposed = FALSE, receipts = completed)
}
