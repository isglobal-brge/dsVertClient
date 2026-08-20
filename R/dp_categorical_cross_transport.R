# Client orchestration for fixed-domain cross-owner categorical tables.
#
# The relay sees only signed public receipts and opaque exact-GC capability
# handles. One-hot inputs, products, cell totals, and additive result shares
# remain inside the two pinned computation peers until the joint-DP vector is
# opened.

.DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION <-
  "fixed-domain-categorical-cross-contingency-v1"
.DSVERT_CLIENT_DP_CATEGORICAL_CROSS_UNIT_POLICY <- paste0(
  "per_owner_consistent_level_else_zero_then_one_joint_cell_per_",
  "admitted_unit_v1")
.DSVERT_CLIENT_DP_CATEGORICAL_CROSS_BIND_VERSION <-
  "dsvert-cross-categorical-exact-binding-v1"
.DSVERT_CLIENT_DP_CATEGORICAL_CROSS_STAGE_VERSION <-
  "dsvert-cross-categorical-exact-stage-v1"
.DSVERT_CLIENT_DP_CATEGORICAL_CROSS_RECEIPT_VERSION <-
  "dsvert-cross-categorical-result-receipt-v1"
.DSVERT_CLIENT_DP_CATEGORICAL_CROSS_PRODUCER <-
  "dp.categorical-cross.v1"
.DSVERT_CLIENT_DP_CATEGORICAL_CROSS_MAX_RECEIPT_BYTES <- 128L * 1024L

.dsvert_dp_categorical_cross_artifacts_client <- function(manifest) {
  artifacts <- tryCatch(
    manifest$workload$families$categorical_pairs$cross_artifacts,
    error = function(error) NULL)
  if (!is.list(artifacts)) return(list())
  artifacts <- artifacts[vapply(artifacts, function(artifact) {
    is.list(artifact) && identical(
      artifact$version,
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION)
  }, logical(1L))]
  if (!length(artifacts)) return(list())
  if (is.null(names(artifacts)) || anyNA(names(artifacts)) ||
      any(!nzchar(names(artifacts))) || anyDuplicated(names(artifacts))) {
    stop("The signed cross-owner categorical artifact set is invalid",
         call. = FALSE)
  }
  artifacts[order(names(artifacts), method = "radix")]
}

.dsvert_dp_categorical_cross_tag_client <- function(
    capsule_id, analysis_id) {
  substr(.dsvert_dp_capsule_source_hash(list(
    protocol = "dsvert-cross-categorical-slot-namespace-v1",
    capsule_id = capsule_id, analysis_id = analysis_id)), 1L, 20L)
}

.dsvert_dp_categorical_cross_verify_signed <- function(
    value, domain, peer, context, expected_fields) {
  if (!.dsvert_dp_has_exact_names(value, expected_fields)) {
    stop("A pinned peer returned a malformed cross-owner categorical receipt",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(value, domain, peer, context)
  value
}

.dsvert_dp_categorical_cross_bind_set <- function(
    responses, context, manifest, layout, analysis_id, source_receipt) {
  peers <- context$designated
  artifact <- .dsvert_dp_categorical_cross_artifacts_client(
    manifest)[[analysis_id]]
  fields <- c(
    "version", "phase", "capsule_id", "analysis_id", "artifact_sha256",
    "source_contract_hash", "private_layout_sha256", "transcript_sha256",
    "numeric_certificate_sha256", "peer_name", "peer_identity_pk",
    "padded_units", "row_level_count", "column_level_count", "ring_bits",
    "frac_bits", "state", "source_values_exposed",
    "alignment_hash_exposed", "alignment_hash_exposed_to_relay",
    "alignment_hash_exposed_to_computation_peers", "exact_intermediates_exposed",
    "fixed_transcript", "signature")
  if (!is.list(responses) || !setequal(names(responses), peers) ||
      is.null(artifact)) {
    stop("Categorical binding did not cover both computation peers",
         call. = FALSE)
  }
  receipts <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_joint_dp_client_decode(
      responses[[peer]], "cross-owner categorical binding receipt",
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_MAX_RECEIPT_BYTES)
    value <- .dsvert_dp_categorical_cross_verify_signed(
      value, "cross-categorical-bind", peer, context, fields)
    valid <- identical(
      value$version, .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_BIND_VERSION) &&
      identical(value$phase, "cross_categorical_private_inputs_bound") &&
      identical(value$capsule_id, source_receipt$capsule_id) &&
      identical(value$analysis_id, analysis_id) &&
      identical(value$artifact_sha256,
                .dsvert_dp_capsule_source_hash(artifact)) &&
      .dsvert_dp_capsule_source_hex(value$source_contract_hash) &&
      identical(value$private_layout_sha256,
                layout$transport_coordinate_order_sha256) &&
      identical(value$transcript_sha256,
                .dsvert_dp_capsule_source_hash(artifact$transcript)) &&
      identical(value$numeric_certificate_sha256,
                .dsvert_dp_capsule_source_hash(artifact$numeric_certificate)) &&
      identical(value$peer_name, peer) &&
      identical(value$peer_identity_pk, unname(context$pinset[[peer]])) &&
      identical(as.numeric(value$padded_units),
                as.numeric(artifact$transcript$padded_units)) &&
      identical(as.numeric(value$row_level_count),
                as.numeric(length(artifact$left$levels))) &&
      identical(as.numeric(value$column_level_count),
                as.numeric(length(artifact$right$levels))) &&
      identical(as.numeric(value$ring_bits), 128) &&
      identical(as.numeric(value$frac_bits),
                as.numeric(artifact$numeric_certificate$frac_bits)) &&
      value$state %in% c("bound", "complete") &&
      identical(value$source_values_exposed, FALSE) &&
      identical(value$alignment_hash_exposed, FALSE) &&
      identical(value$alignment_hash_exposed_to_relay, FALSE) &&
      identical(value$alignment_hash_exposed_to_computation_peers, FALSE) &&
      identical(value$exact_intermediates_exposed, FALSE) &&
      identical(value$fixed_transcript, TRUE)
    if (!isTRUE(valid)) {
      stop("A pinned peer returned a misbound categorical receipt",
           call. = FALSE)
    }
    value
  }), peers)
  stable <- setdiff(fields, c(
    "peer_name", "peer_identity_pk", "state", "signature"))
  if (length(unique(vapply(receipts, function(value) {
        .dsvert_joint_dp_client_json(value[stable])
      }, character(1L)))) != 1L ||
      length(unique(vapply(receipts, `[[`, character(1L), "state"))) != 1L) {
    stop("The peers disagree on the cross-owner categorical binding",
         call. = FALSE)
  }
  list(artifact = artifact, receipts = receipts,
        state = receipts[[1L]]$state,
        capsule_id = receipts[[1L]]$capsule_id,
        source_contract_hash = receipts[[1L]]$source_contract_hash)
}

.dsvert_dp_categorical_cross_stage_contract <- function(
    artifact, capsule_id, analysis_id) {
  capacity <- as.numeric(artifact$transcript$padded_units)
  total_n <- capacity * length(artifact$left$levels) *
    length(artifact$right$levels)
  if (!is.finite(total_n) || total_n != floor(total_n) ||
      total_n < 1 || total_n > 2^31 - 1) {
    stop("The cross-owner categorical exact stage is not representable",
         call. = FALSE)
  }
  tag <- .dsvert_dp_categorical_cross_tag_client(capsule_id, analysis_id)
  list(
    stage = "cell-products", stage_index = 1L,
    total_n = as.integer(total_n),
    purpose = paste0(
      "dp.categorical-cross.", tag, ".cell-products"))
}

.dsvert_dp_categorical_cross_stage_set <- function(
    responses, context, manifest, binding, analysis_id) {
  peers <- context$designated
  artifact <- binding$artifact
  expected <- .dsvert_dp_categorical_cross_stage_contract(
    artifact, binding$capsule_id, analysis_id)
  fields <- c(
    "capability_id", "manifest_handle", "context_hash", "plan_id",
    "ring_bits", "frac_bits", "backend", "bound_x", "bound_y",
    "max_chunk", "total_n", "numeric_policy_id", "version", "state",
    "producer", "purpose", "capsule_id", "analysis_id", "stage",
    "stage_index", "artifact_sha256", "source_contract_hash",
    "transcript_sha256", "numeric_certificate_sha256",
    "exact_intermediates_exposed", "source_values_exposed")
  if (!is.list(responses) || !setequal(names(responses), peers)) {
    stop("Categorical exact preparation did not cover both peers",
         call. = FALSE)
  }
  manifests <- stats::setNames(lapply(peers, function(peer) {
    value <- responses[[peer]]
    bound <- format(
      2^artifact$numeric_certificate$frac_bits,
      scientific = FALSE, trim = TRUE)
    valid <- .dsvert_dp_has_exact_names(value, fields) &&
      identical(value$capability_id, .DSVERT_CLIENT_EXACT_GC_CAPABILITY) &&
      identical(value$version,
                .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_STAGE_VERSION) &&
      identical(value$state, "prepared") &&
      identical(value$producer,
                .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_PRODUCER) &&
      identical(value$purpose, expected$purpose) &&
      identical(value$capsule_id, binding$capsule_id) &&
      identical(value$analysis_id, analysis_id) &&
      identical(value$stage, expected$stage) &&
      identical(as.numeric(value$stage_index), 1) &&
      identical(value$artifact_sha256,
                .dsvert_dp_capsule_source_hash(artifact)) &&
      identical(value$source_contract_hash,
                binding$source_contract_hash) &&
      identical(value$transcript_sha256,
                .dsvert_dp_capsule_source_hash(artifact$transcript)) &&
      identical(value$numeric_certificate_sha256,
                .dsvert_dp_capsule_source_hash(artifact$numeric_certificate)) &&
      .dsvert_dp_is_string(value$manifest_handle) &&
      grepl("^[A-Za-z0-9_-]{43}$", value$manifest_handle) &&
      all(vapply(value[c(
        "context_hash", "plan_id", "numeric_policy_id")],
        .dsvert_dp_capsule_source_hex, logical(1L))) &&
      identical(as.numeric(value$ring_bits), 128) &&
      identical(as.numeric(value$frac_bits),
                as.numeric(artifact$numeric_certificate$frac_bits)) &&
      identical(value$backend, "direct-wide") &&
      identical(value$bound_x, bound) && identical(value$bound_y, bound) &&
      identical(as.numeric(value$total_n), as.numeric(expected$total_n)) &&
      .dsvert_dp_is_integer(value$max_chunk, 1, 256) &&
      identical(value$exact_intermediates_exposed, FALSE) &&
      identical(value$source_values_exposed, FALSE)
    if (!isTRUE(valid)) {
      stop("A computation peer changed the categorical exact stage",
           call. = FALSE)
    }
    value
  }), peers)
  consensus <- setdiff(fields, "manifest_handle")
  if (length(unique(vapply(manifests, function(value) {
        .dsvert_joint_dp_client_json(value[consensus])
      }, character(1L)))) != 1L) {
    stop("The computation peers disagree on the categorical exact stage",
         call. = FALSE)
  }
  manifests
}

.dsvert_dp_categorical_cross_result_set <- function(
    responses, context, manifest, layout, binding, analysis_id) {
  peers <- context$designated
  artifact <- binding$artifact
  release_layout <- .dsvert_dp_capsule_vector_layout(manifest)
  block <- release_layout$blocks[[paste(
    "categorical_pairs", "cross", analysis_id, sep = "::")]]
  fields <- c(
    "version", "phase", "capsule_id", "analysis_id", "peer_name",
    "peer_identity_pk", "artifact_sha256", "source_contract_hash",
    "private_layout_sha256", "transcript_sha256",
    "numeric_certificate_sha256", "exact_transcript_sha256",
    "coordinate_count", "release_start", "release_end",
    "release_coordinate_order_sha256", "ring_bits", "frac_bits", "state",
    "fixed_transcript", "result_share_exposed",
    "exact_intermediates_exposed", "alignment_hash_exposed",
    "alignment_hash_exposed_to_relay",
    "alignment_hash_exposed_to_computation_peers", "signature")
  if (!is.list(responses) || !setequal(names(responses), peers) ||
      !is.list(block)) {
    stop("Categorical finalization did not cover both computation peers",
         call. = FALSE)
  }
  receipts <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_joint_dp_client_decode(
      responses[[peer]], "cross-owner categorical result receipt",
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_MAX_RECEIPT_BYTES)
    value <- .dsvert_dp_categorical_cross_verify_signed(
      value, "cross-categorical-result", peer, context, fields)
    valid <- identical(
      value$version,
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_RECEIPT_VERSION) &&
      identical(value$phase, "cross_categorical_result_share_persisted") &&
      identical(value$capsule_id, binding$capsule_id) &&
      identical(value$analysis_id, analysis_id) &&
      identical(value$peer_name, peer) &&
      identical(value$peer_identity_pk, unname(context$pinset[[peer]])) &&
      identical(value$artifact_sha256,
                .dsvert_dp_capsule_source_hash(artifact)) &&
      identical(value$source_contract_hash,
                binding$source_contract_hash) &&
      identical(value$private_layout_sha256,
                layout$transport_coordinate_order_sha256) &&
      identical(value$transcript_sha256,
                .dsvert_dp_capsule_source_hash(artifact$transcript)) &&
      identical(value$numeric_certificate_sha256,
                .dsvert_dp_capsule_source_hash(artifact$numeric_certificate)) &&
      .dsvert_dp_capsule_source_hex(value$exact_transcript_sha256) &&
      identical(as.numeric(value$coordinate_count),
                as.numeric(block$length)) &&
      identical(as.numeric(value$release_start), as.numeric(block$start)) &&
      identical(as.numeric(value$release_end), as.numeric(block$end)) &&
      identical(value$release_coordinate_order_sha256,
                release_layout$sha256) &&
      identical(as.numeric(value$ring_bits), 128) &&
      identical(as.numeric(value$frac_bits),
                as.numeric(artifact$numeric_certificate$frac_bits)) &&
      identical(value$state, "complete") &&
      identical(value$fixed_transcript, TRUE) &&
      identical(value$result_share_exposed, FALSE) &&
      identical(value$exact_intermediates_exposed, FALSE) &&
      identical(value$alignment_hash_exposed, FALSE) &&
      identical(value$alignment_hash_exposed_to_relay, FALSE) &&
      identical(value$alignment_hash_exposed_to_computation_peers, FALSE)
    if (!isTRUE(valid)) {
      stop("A computation peer returned an invalid categorical result receipt",
           call. = FALSE)
    }
    value
  }), peers)
  comparable <- setdiff(fields, c(
    "peer_name", "peer_identity_pk", "signature"))
  if (length(unique(vapply(receipts, function(value) {
        .dsvert_joint_dp_client_json(value[comparable])
      }, character(1L)))) != 1L) {
    stop("The peers disagree on the exact categorical result",
         call. = FALSE)
  }
  receipts
}

.dsvert_dp_categorical_cross_orchestrate <- function(
    manifest_json, manifest, context, source_receipt, .aggregate,
    .setup_exact = .dsvert_setup_exact_gc_transport,
    .vecmul = .dsvert_exact_gc_vecmul_run,
    .alignment_mask = .dsvert_dp_alignment_mask_run,
    .shared_exact = NULL, .remote_context = NULL) {
  artifacts <- .dsvert_dp_categorical_cross_artifacts_client(manifest)
  if (!length(artifacts)) {
    return(list(enabled = FALSE, sampler_handoff_ready = TRUE,
                exact_intermediates_exposed = FALSE,
                source_values_exposed = FALSE))
  }
  if (!is.function(.setup_exact) || !is.function(.vecmul) ||
      !is.function(.alignment_mask)) {
    stop("Invalid cross-owner categorical exact transport implementation",
         call. = FALSE)
  }
  layout <- .dsvert_dp_gaussian_cross_layout_client(manifest)
  peers <- context$designated
  synopsis <- !is.null(.remote_context)
  openings <- if (isTRUE(synopsis)) NULL else
    .dsvert_dp_capsule_allocation_openings(context)
  if (isTRUE(synopsis)) {
    .remote_context <-
      .dsvert_dp_synopsis_categorical_cross_remote_context_v1(
        .remote_context)
  }
  valid_source <- identical(
      sort(peers, method = "radix"),
      unname(unlist(layout$computation_peers))) &&
    identical(source_receipt$purpose,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CATEGORICAL_CROSS_PURPOSE) &&
    identical(as.numeric(source_receipt$coordinate_count),
              as.numeric(layout$transport_coordinate_count)) &&
    identical(as.numeric(source_receipt$release_coordinate_count),
              as.numeric(layout$release_coordinate_count)) &&
    identical(source_receipt$private_layout_sha256,
              layout$transport_coordinate_order_sha256) &&
    identical(source_receipt$sampler_handoff_ready, FALSE) &&
    identical(source_receipt$payload_exposed, FALSE)
  if (!isTRUE(valid_source)) {
    stop("The cross-owner categorical source transport is misbound",
         call. = FALSE)
  }
  if (is.null(.shared_exact)) {
    session_id <- .dsvert_uuid4()
    setup_result <- .dsvert_dp_cross_exact_setup(
      .setup_exact,
      context$all_conns, context$servers,
      match(peers, context$servers), session_id,
      .aggregate = .aggregate)
    on.exit(.dsvert_dp_cross_exact_cleanup(
      context$conns, session_id, setup_result, .aggregate, .setup_exact),
      add = TRUE)
    alignment <- if (isTRUE(synopsis) && identical(
        .alignment_mask, .dsvert_dp_alignment_mask_run)) {
      .alignment_mask(
        manifest_json, context, layout, source_receipt,
        session_id, .aggregate, .remote_context = .remote_context)
    } else {
      .alignment_mask(
        manifest_json, context, layout, source_receipt,
        session_id, .aggregate)
    }
    .shared_exact <- .dsvert_dp_cross_shared_exact_build(
      manifest_json, manifest, context, layout, source_receipt,
      session_id, alignment)
  } else {
    .shared_exact <- .dsvert_dp_cross_shared_exact_validate(
      .shared_exact, manifest_json, manifest, context, layout,
      source_receipt)
    session_id <- .shared_exact$session_id
  }
  completed <- list()
  for (analysis_id in names(artifacts)) {
    bind_calls <- stats::setNames(lapply(peers, function(peer) {
      if (isTRUE(synopsis)) {
        .dsvert_dp_synopsis_categorical_cross_bind_call_v1(
          .remote_context, analysis_id, session_id)
      } else call(
        name = "dsvertDPCategoricalCrossBindDS",
        manifest_json = manifest_json, analysis_id = analysis_id,
        session_id = session_id,
        first_opening_json = openings[[peers[[1L]]]],
        second_opening_json = openings[[peers[[2L]]]])
    }), peers)
    binding <- .dsvert_dp_categorical_cross_bind_set(
      .dsvert_fanout_by_site(
        context$conns, bind_calls,
        operation = "cross-owner categorical private-input binding",
        .aggregate = .aggregate),
      context, manifest, layout, analysis_id, source_receipt)
    if (!identical(binding$source_contract_hash,
                   source_receipt$contract_hash)) {
      stop("The exact categorical binding changed the source contract",
           call. = FALSE)
    }
    if (!identical(binding$state, "complete")) {
      contract <- .dsvert_dp_categorical_cross_stage_contract(
        binding$artifact, binding$capsule_id, analysis_id)
      prepare_calls <- stats::setNames(lapply(peers, function(peer) call(
        name = "dsvertDPCategoricalCrossPrepareDS",
        analysis_id = analysis_id, session_id = session_id)), peers)
      producer_manifests <- .dsvert_dp_categorical_cross_stage_set(
        .dsvert_fanout_by_site(
          context$conns, prepare_calls,
          operation = "cross-owner categorical exact preparation",
          .aggregate = .aggregate),
        context, manifest, binding, analysis_id)
      .vecmul(
        context$all_conns, server_names = context$servers,
        servers = match(peers, context$servers), session_id = session_id,
        total_n = contract$total_n,
        input_manifests = producer_manifests,
        transport_ready = TRUE, .aggregate = .aggregate)
    }
    finalize_calls <- stats::setNames(lapply(peers, function(peer) {
      if (isTRUE(synopsis)) {
        .dsvert_dp_synopsis_categorical_cross_finalize_call_v1(
          .remote_context, analysis_id, session_id)
      } else call(
        name = "dsvertDPCategoricalCrossFinalizeDS",
        manifest_json = manifest_json, analysis_id = analysis_id,
        session_id = session_id)
    }), peers)
    completed[[analysis_id]] <- .dsvert_dp_categorical_cross_result_set(
      .dsvert_fanout_by_site(
        context$conns, finalize_calls,
        operation = "cross-owner categorical exact finalization",
        .aggregate = .aggregate),
      context, manifest, layout, binding, analysis_id)
  }
  list(
    enabled = TRUE, sampler_handoff_ready = TRUE,
    exact_intermediates_exposed = FALSE, source_values_exposed = FALSE,
    analyses = names(completed), receipts = completed,
    private_layout_sha256 = layout$transport_coordinate_order_sha256)
}

.dsvert_dp_cross_orchestrate <- function(
    manifest_json, manifest, context, source_receipt, .aggregate,
    .setup_exact = .dsvert_setup_exact_gc_transport,
    .alignment_mask = .dsvert_dp_alignment_mask_run) {
  categorical <- .dsvert_dp_categorical_cross_artifacts_client(manifest)
  gaussian <- .dsvert_dp_gaussian_cross_artifacts_client(manifest)
  completed <- list()
  shared_exact <- NULL
  share_session <- length(categorical) > 0L && length(gaussian) > 0L
  if (isTRUE(share_session)) {
    if (!is.function(.setup_exact) || !is.function(.alignment_mask)) {
      stop("Invalid shared cross-owner exact transport implementation",
           call. = FALSE)
    }
    layout <- .dsvert_dp_gaussian_cross_layout_client(manifest)
    peers <- context$designated
    .dsvert_dp_capsule_allocation_openings(context)
    valid_source <- identical(
        sort(peers, method = "radix"),
        unname(unlist(layout$computation_peers))) &&
      identical(
        source_receipt$purpose,
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CATEGORICAL_CROSS_PURPOSE) &&
      identical(as.numeric(source_receipt$coordinate_count),
                as.numeric(layout$transport_coordinate_count)) &&
      identical(as.numeric(source_receipt$release_coordinate_count),
                as.numeric(layout$release_coordinate_count)) &&
      identical(source_receipt$private_layout_sha256,
                layout$transport_coordinate_order_sha256) &&
      identical(source_receipt$sampler_handoff_ready, FALSE) &&
      identical(source_receipt$payload_exposed, FALSE)
    if (!isTRUE(valid_source)) {
      stop("The shared cross-owner source transport is misbound",
           call. = FALSE)
    }
    session_id <- .dsvert_uuid4()
    setup_result <- .dsvert_dp_cross_exact_setup(
      .setup_exact,
      context$all_conns, context$servers,
      match(peers, context$servers), session_id,
      .aggregate = .aggregate)
    on.exit(.dsvert_dp_cross_exact_cleanup(
      context$conns, session_id, setup_result, .aggregate, .setup_exact),
      add = TRUE)
    alignment <- .alignment_mask(
      manifest_json, context, layout, source_receipt,
      session_id, .aggregate)
    shared_exact <- .dsvert_dp_cross_shared_exact_build(
      manifest_json, manifest, context, layout, source_receipt,
      session_id, alignment)
  }
  if (length(categorical)) {
    completed$categorical <- .dsvert_dp_categorical_cross_orchestrate(
      manifest_json, manifest, context, source_receipt, .aggregate,
      .setup_exact = .setup_exact, .alignment_mask = .alignment_mask,
      .shared_exact = shared_exact)
  }
  if (length(gaussian)) {
    completed$gaussian <- .dsvert_dp_gaussian_cross_orchestrate(
      manifest_json, manifest, context, source_receipt, .aggregate,
      .setup_exact = .setup_exact, .alignment_mask = .alignment_mask,
      .shared_exact = shared_exact)
  }
  if (!length(completed)) {
    if (!identical(source_receipt$sampler_handoff_ready, TRUE)) {
      stop("The ordinary capsule source is not ready for joint sampling",
           call. = FALSE)
    }
    return(list(enabled = FALSE, sampler_handoff_ready = TRUE,
                exact_intermediates_exposed = FALSE,
                source_values_exposed = FALSE))
  }
  list(
    enabled = TRUE, sampler_handoff_ready = TRUE,
    exact_intermediates_exposed = FALSE, source_values_exposed = FALSE,
    families = names(completed), results = completed,
    shared_exact_session = isTRUE(share_session),
    exact_session_count = if (isTRUE(share_session)) 1L else length(completed),
    alignment_gate_count = if (isTRUE(share_session)) 1L else length(completed))
}
