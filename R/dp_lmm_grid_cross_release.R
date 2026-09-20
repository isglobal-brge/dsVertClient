# Internal orchestration, reached only after the signed workload admission gate.
.dsvert_dp_lmm_cross_public_evidence_set <- function(
    responses, context, manifest, analysis_id, release, compiled) {
  artifact <- .dsvert_dp_glm_grid_cross_artifacts(manifest)[[analysis_id]]
  if (!is.list(artifact) || !.dsvert_dp_staged_grouped_artifact(artifact)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  .dsvert_dp_staged_cross_public_evidence_set(
    responses, context, manifest, artifact, release, compiled)
}

# Shared verifier for an already admitted typed artifact. The family-specific
# caller validates admission; this function never discovers or starts a release.
.dsvert_dp_staged_cross_public_evidence_set <- function(
    responses, context, manifest, artifact, release, compiled) {
  reference <- .dsvert_dp_lmm_cross_receipts(responses, context, artifact, "published")
  for (field in c("artifact_key", "execution_id", "final_vector_root", "result_set_sha256")) {
    if (!.dsvert_dp_capsule_source_hex(reference[[field]]) ||
        !identical(reference[[field]], release[[field]])) .dsvert_dp_glm_grid_cross_fail()
  }
  if (!identical(reference$artifact_key, compiled$artifact$artifact_key)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  namespace <- list(version = "dsvert-stateless-catalog-synopsis-source-contract-v1",
    manifest_capsule_id = manifest$capsule_identity$capsule_id,
    artifact_key = compiled$artifact$artifact_key,
    source_claim_set_sha256 = compiled$artifact$semantic$source_claim_set_sha256)
  capsule_id <- .dsvert_dp_synopsis_client_hash(
    "dsVert/stateless-catalog-synopsis/source-namespace/v1|", namespace)
  semantic_key <- .dsvert_dp_capsule_source_hash(list(
    version = "cross-grid-semantic-release-key-v2", capsule_id = capsule_id,
    source_contract_sha256 = release$source_contract_sha256,
    signed_contract = artifact$signed_contract,
    family = artifact$family, version_family = artifact$spec_version,
    mechanism = manifest$workload$capsule_mechanism,
    alignment = .dsvert_dp_glm_grid_cross_embedded_contract(artifact)$spec$alignment,
    caps = artifact$sensitivity))
  if (!identical(reference$capsule_id, capsule_id) ||
      !identical(reference$source_contract_sha256, release$source_contract_sha256) ||
      !identical(reference$semantic_key, semantic_key)) .dsvert_dp_glm_grid_cross_fail()
  # Keep both actual signatures so an offline/cold certificate never relies on
  # the orchestration's in-memory reference receipt or a synthesized authority.
  stats::setNames(lapply(context$designated, function(peer) {
    .dsvert_joint_dp_client_decode(responses[[peer]], "LMM public evidence",
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_RECEIPT_BYTES)
  }), context$designated)
}

.dsvert_dp_lmm_cross_receipts <- function(responses, context, artifact, phase) {
  peers <- context$designated
  if (!is.list(responses) || length(responses) != 2L ||
      anyDuplicated(names(responses)) || !setequal(names(responses), peers)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  extra <- switch(phase, bound = character(),
    prepared = c("operation_id", "source_key", "output_key", "operation", "purpose", "vector_len", "persisted"),
    persisted = "stage_receipt", complete = c("coordinate_count", "stage_receipt"),
    published = c("coordinate_count", "stage_receipt", "artifact_key", "execution_id",
      "final_vector_root", "result_set_sha256"),
    .dsvert_dp_glm_grid_cross_fail())
  fields <- c("version", "phase", "capsule_id", "analysis_id", "peer_name",
    "peer_identity_pk", "semantic_key", "artifact_sha256", "source_contract_sha256",
    "profile_sha256", "certificate_sha256", "stage_plan_digest", "private_result_exposed", "signature")
  results <- lapply(peers, function(peer) {
    value <- .dsvert_joint_dp_client_decode(responses[[peer]], "staged LMM receipt",
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_RECEIPT_BYTES)
    .dsvert_dp_glm_grid_cross_fields(value, c(fields, extra))
    .dsvert_dp_capsule_source_verify(value, "cross-grid-result", peer, context)
    if (!identical(value$version, .dsvert_dp_staged_grouped_tag(artifact, "-staged-receipt-v1", "dsvert-")) ||
        !identical(value$phase, phase) || !identical(value$analysis_id, artifact$analysis_id) ||
        !identical(value$peer_name, peer) ||
        !identical(value$peer_identity_pk, unname(context$pinset[[peer]])) ||
        !identical(value$artifact_sha256, .dsvert_dp_capsule_source_hash(artifact)) ||
        !identical(value$profile_sha256, artifact$numeric_certificate$profile_sha256) ||
        !identical(value$certificate_sha256, artifact$numeric_certificate$certificate_sha256) ||
        !identical(value$private_result_exposed, FALSE)) .dsvert_dp_glm_grid_cross_fail()
    for (field in c("stage_plan_digest", if (phase %in% c("persisted", "complete", "published")) "stage_receipt")) {
      if (!is.character(value[[field]]) || length(value[[field]]) != 1L || is.na(value[[field]]) ||
          !grepl("^[0-9a-f]{64}$", value[[field]]) || identical(value[[field]], strrep("0", 64))) {
        .dsvert_dp_glm_grid_cross_fail()
      }
    }
    operation <- if (identical(artifact$family, "cox")) "cox-loss-staged-v1" else
      .dsvert_dp_staged_grouped_tag(artifact, "-staged-v1", "grouped-")
    if (phase == "prepared" && (!identical(value$operation, operation) ||
        !is.character(value$purpose) || length(value$purpose) != 1L || is.na(value$purpose) ||
        !grepl(paste0("^", operation, "/[0-9a-f]{64}$"), value$purpose) ||
        !identical(as.numeric(value$vector_len), as.numeric(artifact$coordinate_count)) ||
        !is.logical(value$persisted) || length(value$persisted) != 1L || is.na(value$persisted))) {
      .dsvert_dp_glm_grid_cross_fail()
    }
    if (phase %in% c("complete", "published") && !identical(as.numeric(value$coordinate_count),
        as.numeric(artifact$coordinate_count))) .dsvert_dp_glm_grid_cross_fail()
    value
  })
  common <- function(value) value[setdiff(names(value),
    c("peer_name", "peer_identity_pk", "signature", if (phase == "prepared") "persisted"))]
  .dsvert_dp_glm_grid_cross_equal(common(results[[1L]]), common(results[[2L]]))
  result <- results[[1L]]
  if (phase == "prepared") attr(result, "both_persisted") <-
    all(vapply(results, function(value) isTRUE(value$persisted), logical(1L)))
  result
}

.dsvert_dp_lmm_cross_orchestrate <- function(manifest_json, manifest, context,
    source_receipt, .aggregate, .remote_context) {
  artifacts <- .dsvert_dp_glm_grid_cross_artifacts(manifest)
  if (length(artifacts) != 1L || !.dsvert_dp_staged_grouped_artifact(artifacts[[1L]]) ||
      is.null(.remote_context)) .dsvert_dp_glm_grid_cross_fail()
  artifact <- artifacts[[1L]]
  layout <- .dsvert_dp_gaussian_cross_layout_client(manifest)
  .dsvert_dp_staged_cross_orchestrate(manifest_json, manifest, context,
    source_receipt, artifact, layout, .aggregate, .remote_context)
}

# Shared authenticated executor; callers admit a typed artifact and its layout.
.dsvert_dp_staged_cross_orchestrate <- function(manifest_json, manifest, context,
    source_receipt, artifact, layout, .aggregate, .remote_context) {
  if (is.null(.remote_context)) .dsvert_dp_glm_grid_cross_fail()
  if (!identical(source_receipt$purpose, .DSVERT_CLIENT_DP_GLM_GRID_CROSS_SOURCE_PURPOSE) ||
      !identical(as.numeric(source_receipt$coordinate_count), as.numeric(layout$transport_coordinate_count)) ||
      !identical(as.numeric(source_receipt$release_coordinate_count), as.numeric(layout$release_coordinate_count)) ||
      !identical(source_receipt$private_layout_sha256, layout$transport_coordinate_order_sha256) ||
      !identical(source_receipt$sampler_handoff_ready, FALSE) ||
      !identical(source_receipt$payload_exposed, FALSE)) .dsvert_dp_glm_grid_cross_fail()
  peers <- context$designated
  session_id <- .dsvert_uuid4()
  prepare_calls <- stats::setNames(lapply(peers, function(peer) {
    as.call(c(list(as.name("dsvertDPSynopsisPrepareDS")), list(session_id = session_id), .remote_context))
  }), peers)
  .dsvert_dp_synopsis_runner_json_set(.dsvert_fanout_by_site(context$conns,
    prepare_calls, operation = "staged LMM Synopsis authority", .aggregate = .aggregate),
    peers, "PREPARE", .DSVERT_CLIENT_SYNOPSIS_PREPARE_MAX_OBJECT_BYTES)
  setup <- .dsvert_dp_cross_exact_setup(.dsvert_setup_exact_gc_transport,
    context$all_conns, context$servers, match(peers, context$servers), session_id, .aggregate = .aggregate)
  on.exit(.dsvert_dp_cross_exact_cleanup(context$conns, session_id, setup,
    .aggregate, .dsvert_setup_exact_gc_transport), add = TRUE)
  .dsvert_dp_alignment_mask_run(manifest_json, context, layout, source_receipt,
    session_id, .aggregate, .remote_context = .remote_context)
  invoke <- function(action) {
    calls <- stats::setNames(lapply(peers, function(peer) {
      as.call(c(list(as.name("dsvertDPSynopsisGLMGridCrossDS")), .remote_context,
        list(analysis_id = artifact$analysis_id, session_id = session_id, action = action, batch = 0)))
    }), peers)
    .dsvert_fanout_by_site(context$conns, calls, operation = paste("staged LMM", action), .aggregate = .aggregate)
  }
  bound <- if (identical(artifact$family, "cox")) {
    .dsvert_dp_cox_cross_bind(manifest, context, artifact, source_receipt,
      session_id, .remote_context, .aggregate)
  } else .dsvert_dp_lmm_cross_receipts(invoke("bind"), context, artifact, "bound")
  if (!identical(bound$source_contract_sha256, source_receipt$contract_hash) ||
      !identical(bound$capsule_id, source_receipt$capsule_id)) .dsvert_dp_glm_grid_cross_fail()
  checked <- function(responses, phase) {
    receipt <- .dsvert_dp_lmm_cross_receipts(responses, context, artifact, phase)
    for (field in c("source_contract_sha256", "capsule_id", "semantic_key", "stage_plan_digest")) {
      if (!identical(receipt[[field]], bound[[field]])) .dsvert_dp_glm_grid_cross_fail()
    }
    receipt
  }
  stage <- checked(invoke("prepare"), "prepared")
  persisted <- NULL
  if (!isTRUE(attr(stage, "both_persisted"))) {
    initialized <- invoke("start")
    .dsvert_exact_gc_run(context$all_conns, server_names = context$servers,
      servers = match(peers, context$servers), session_id = session_id,
      operation_id = stage$operation_id, source_key = stage$source_key,
      output_key = stage$output_key, operation = stage$operation, ring = 128L,
      frac_bits = 0L, vector_len = stage$vector_len, purpose = stage$purpose,
      transport_ready = TRUE, initialized = initialized, .aggregate = .aggregate)
    persisted <- checked(invoke("store"), "persisted")
  }
  complete <- checked(invoke("finalize"), "complete")
  if (!is.null(persisted) && !identical(persisted$stage_receipt, complete$stage_receipt)) {
    .dsvert_dp_glm_grid_cross_fail()
  }
  list(enabled = TRUE, sampler_handoff_ready = TRUE, private_result_exposed = FALSE,
    receipts = setNames(list(complete), artifact$analysis_id))
}
