# Shared exact-MPC transport for one immutable cross-owner capsule.
#
# The alignment predicate is capsule/source-contract scoped rather than
# analysis-family scoped.  Categorical and Gaussian stages may therefore use
# one authenticated peer session and one completed alignment batch.  Every
# later exact operation remains bound to its own family, analysis, stage,
# producer and purpose by the existing server-owned contracts.

.DSVERT_CLIENT_DP_CROSS_SHARED_EXACT_VERSION <-
  "dsvert-cross-shared-exact-session-v1"
.DSVERT_CLIENT_DP_CROSS_SHARED_EXACT_DOMAINS <-
  "family-analysis-stage-producer-purpose-separated-v1"

.dsvert_dp_cross_exact_setup <- function(
    .setup_exact, datasources, server_names, servers, session_id,
    .aggregate) {
  production <- identical(
    .setup_exact, .dsvert_setup_exact_gc_transport)
  result <- if (isTRUE(production)) {
    .setup_exact(
      datasources, server_names, servers, session_id,
      cleanup_purpose = .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
      .aggregate = .aggregate)
  } else {
    # The alternate implementation is an internal test seam only.
    .setup_exact(
      datasources, server_names, servers, session_id,
      .aggregate = .aggregate)
  }
  if (isTRUE(production) &&
      (is.null(attr(result, "exact_gc_cleanup_capabilities", exact = TRUE)) ||
       !identical(attr(result, "exact_gc_cleanup_purpose", exact = TRUE),
                  .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE))) {
    stop("The exact cross-owner session has no cleanup capability.",
         call. = FALSE)
  }
  result
}

.dsvert_dp_cross_exact_cleanup <- function(
    conns, session_id, setup_result, .aggregate, .setup_exact) {
  capabilities <- attr(
    setup_result, "exact_gc_cleanup_capabilities", exact = TRUE)
  if (!is.null(capabilities)) {
    return(.dsvert_exact_gc_cleanup_best_effort(
      conns, session_id, setup_result, .aggregate = .aggregate))
  }
  # Test seams without a minted capability have no remote cleanup authority.
  invisible(FALSE)
}

.dsvert_dp_cross_shared_session_id <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl(
        paste0("^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-",
               "[89ab][0-9a-f]{3}-[0-9a-f]{12}$"), value)) {
    stop("Invalid shared exact-MPC session identifier", call. = FALSE)
  }
  value
}

.dsvert_dp_cross_shared_bindings <- function(context, layout) {
  peers <- unname(context$designated)
  sources <- unname(unlist(layout$source_peers))
  pins <- context$pinset[peers]
  valid <- length(peers) == 2L && !anyDuplicated(peers) &&
    length(sources) >= 2L && length(sources) <= 64L &&
    !anyDuplicated(sources) && is.character(pins) && length(pins) == 2L &&
    !anyNA(pins) && all(nzchar(pins)) &&
    all(peers %in% context$servers) &&
    identical(
      sort(peers, method = "radix"),
      sort(unname(unlist(layout$computation_peers)), method = "radix"))
  if (!isTRUE(valid)) {
    stop("Invalid shared exact-MPC peer binding", call. = FALSE)
  }
  list(
    peers = peers, sources = sources,
    peer_binding_sha256 = .dsvert_dp_capsule_source_hash(list(
      designated_peers = as.list(peers),
      designated_identity_pins = as.list(unname(pins)),
      ordered_server_pool = as.list(unname(context$servers)))),
    source_binding_sha256 = .dsvert_dp_capsule_source_hash(list(
      ordered_source_peers = as.list(sources))))
}

.dsvert_dp_cross_alignment_validate <- function(
    alignment, layout, source_receipt) {
  fields <- c(
    "capability_id", "state", "batch_operation_id", "contract_hash",
    "source_count", "coordinate_count", "chunk_count",
    "alignment_digest_exposed", "mismatch_source_exposed",
    "gate_share_exposed", "fixed_transcript")
  sources <- unname(unlist(layout$source_peers))
  source_count <- length(sources)
  total <- .dsvert_dp_alignment_mask_private_projection_client(
    layout)$coordinate_count
  chunk_size <- .dsvert_dp_alignment_mask_chunk_size_client(source_count)
  chunk_count <- as.integer(ceiling(total / chunk_size))
  valid <- .dsvert_dp_has_exact_names(alignment, fields) &&
    identical(alignment$capability_id,
              .DSVERT_CLIENT_EXACT_GC_CAPABILITY) &&
    identical(alignment$state, "complete") &&
    is.character(alignment$batch_operation_id) &&
    length(alignment$batch_operation_id) == 1L &&
    !is.na(alignment$batch_operation_id) &&
    grepl("^op_[0-9a-f]{32}$", alignment$batch_operation_id) &&
    identical(alignment$contract_hash, source_receipt$contract_hash) &&
    identical(as.numeric(alignment$source_count),
              as.numeric(source_count)) &&
    identical(as.numeric(alignment$coordinate_count), total) &&
    identical(as.numeric(alignment$chunk_count),
              as.numeric(chunk_count)) &&
    identical(alignment$alignment_digest_exposed, FALSE) &&
    identical(alignment$mismatch_source_exposed, FALSE) &&
    identical(alignment$gate_share_exposed, FALSE) &&
    identical(alignment$fixed_transcript, TRUE)
  if (!isTRUE(valid)) {
    stop("The exact private alignment gate is not capsule-bound",
         call. = FALSE)
  }
  alignment
}

.dsvert_dp_cross_shared_exact_build <- function(
    manifest_json, manifest, context, layout, source_receipt,
    session_id, alignment) {
  session_id <- .dsvert_dp_cross_shared_session_id(session_id)
  bindings <- .dsvert_dp_cross_shared_bindings(context, layout)
  alignment <- .dsvert_dp_cross_alignment_validate(
    alignment, layout, source_receipt)
  capsule_id <- source_receipt$capsule_id
  valid_source <- .dsvert_dp_capsule_source_hex(capsule_id) &&
    .dsvert_dp_capsule_source_hex(source_receipt$contract_hash) &&
    identical(source_receipt$private_layout_sha256,
              layout$transport_coordinate_order_sha256) &&
    identical(as.numeric(source_receipt$coordinate_count),
              as.numeric(layout$transport_coordinate_count)) &&
    is.character(source_receipt$purpose) &&
    length(source_receipt$purpose) == 1L && !is.na(source_receipt$purpose)
  if (!isTRUE(valid_source)) {
    stop("The shared exact-MPC source contract is misbound", call. = FALSE)
  }
  list(
    version = .DSVERT_CLIENT_DP_CROSS_SHARED_EXACT_VERSION,
    session_id = session_id,
    manifest_sha256 = digest::digest(
      charToRaw(manifest_json), algo = "sha256", serialize = FALSE),
    capsule_id = capsule_id,
    source_contract_hash = source_receipt$contract_hash,
    source_purpose = source_receipt$purpose,
    private_layout_sha256 = layout$transport_coordinate_order_sha256,
    peer_binding_sha256 = bindings$peer_binding_sha256,
    source_binding_sha256 = bindings$source_binding_sha256,
    source_count = as.integer(length(bindings$sources)),
    coordinate_count = as.numeric(layout$transport_coordinate_count),
    alignment_batch_operation_id = alignment$batch_operation_id,
    alignment_chunk_count = as.integer(alignment$chunk_count),
    alignment_gate_complete = TRUE,
    alignment_digest_exposed = FALSE,
    mismatch_source_exposed = FALSE,
    gate_share_exposed = FALSE,
    exact_intermediates_exposed = FALSE,
    family_operation_domains =
      .DSVERT_CLIENT_DP_CROSS_SHARED_EXACT_DOMAINS)
}

.dsvert_dp_cross_shared_exact_validate <- function(
    shared, manifest_json, manifest, context, layout, source_receipt) {
  fields <- c(
    "version", "session_id", "manifest_sha256", "capsule_id",
    "source_contract_hash", "source_purpose", "private_layout_sha256",
    "peer_binding_sha256", "source_binding_sha256", "source_count",
    "coordinate_count", "alignment_batch_operation_id",
    "alignment_chunk_count", "alignment_gate_complete",
    "alignment_digest_exposed", "mismatch_source_exposed",
    "gate_share_exposed", "exact_intermediates_exposed",
    "family_operation_domains")
  bindings <- .dsvert_dp_cross_shared_bindings(context, layout)
  expected_chunks <- ceiling(
    .dsvert_dp_alignment_mask_private_projection_client(
      layout)$coordinate_count /
      .dsvert_dp_alignment_mask_chunk_size_client(length(bindings$sources)))
  valid <- .dsvert_dp_has_exact_names(shared, fields) &&
    identical(shared$version,
              .DSVERT_CLIENT_DP_CROSS_SHARED_EXACT_VERSION) &&
    identical(
      tryCatch(.dsvert_dp_cross_shared_session_id(shared$session_id),
               error = function(error) ""),
      shared$session_id) &&
    identical(shared$manifest_sha256, digest::digest(
      charToRaw(manifest_json), algo = "sha256", serialize = FALSE)) &&
    identical(shared$capsule_id, source_receipt$capsule_id) &&
    identical(shared$source_contract_hash,
              source_receipt$contract_hash) &&
    identical(shared$source_purpose, source_receipt$purpose) &&
    identical(shared$private_layout_sha256,
              layout$transport_coordinate_order_sha256) &&
    identical(shared$private_layout_sha256,
              source_receipt$private_layout_sha256) &&
    identical(shared$peer_binding_sha256,
              bindings$peer_binding_sha256) &&
    identical(shared$source_binding_sha256,
              bindings$source_binding_sha256) &&
    identical(as.numeric(shared$source_count),
              as.numeric(length(bindings$sources))) &&
    identical(as.numeric(shared$coordinate_count),
              as.numeric(layout$transport_coordinate_count)) &&
    is.character(shared$alignment_batch_operation_id) &&
    length(shared$alignment_batch_operation_id) == 1L &&
    !is.na(shared$alignment_batch_operation_id) &&
    grepl("^op_[0-9a-f]{32}$", shared$alignment_batch_operation_id) &&
    identical(as.numeric(shared$alignment_chunk_count),
              as.numeric(expected_chunks)) &&
    identical(shared$alignment_gate_complete, TRUE) &&
    identical(shared$alignment_digest_exposed, FALSE) &&
    identical(shared$mismatch_source_exposed, FALSE) &&
    identical(shared$gate_share_exposed, FALSE) &&
    identical(shared$exact_intermediates_exposed, FALSE) &&
    identical(shared$family_operation_domains,
              .DSVERT_CLIENT_DP_CROSS_SHARED_EXACT_DOMAINS)
  if (!isTRUE(valid)) {
    stop("The shared exact-MPC session is not bound to this capsule",
         call. = FALSE)
  }
  shared
}
