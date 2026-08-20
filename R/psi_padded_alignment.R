# Internal client orchestrator for the fixed-capacity pinned PSI
# protocol. The relay handles only authenticated ciphertext envelopes and a
# public, server-signed contract; it never receives identifiers, row maps,
# membership bits, exact counts or GC shares.

.DSVERT_CLIENT_PSI_PADDED_PROTOCOL <- "dsvert-pinned-padded-psi-v5"

.dsvert_validate_psi_padded_attestation <- function(
    value, expected_contract = NULL) {
  invalid <- function() stop(
    "Pinned padded PSI returned an invalid alignment attestation.",
    call. = FALSE)
  required <- c(
    "attestation_version", "alignment_attested", "alignment_protocol",
    "attestation_id", "contract_hash", "policy_id", "alignment_purpose",
    "dataset_id", "dataset_version", "id_column", "source_binding_id",
    "pinset_id",
    "capacity_bucket", "relay_frame_bytes", "inline_max_bytes",
    "peer_count", "reference_peer", "compute_peers")
  scalar <- function(candidate, pattern) {
    is.character(candidate) && length(candidate) == 1L &&
      !is.na(candidate) && grepl(pattern, candidate, perl = TRUE)
  }
  integer <- function(candidate, lower, upper) {
    candidate <- suppressWarnings(as.numeric(candidate))
    length(candidate) == 1L && !is.na(candidate) && is.finite(candidate) &&
      candidate == floor(candidate) && candidate >= lower && candidate <= upper
  }
  if (!is.list(value) || !identical(names(value), required) ||
      !integer(value$attestation_version, 3L, 3L) ||
      !identical(value$alignment_attested, TRUE) ||
      !identical(value$alignment_protocol,
                 .DSVERT_CLIENT_PSI_PADDED_PROTOCOL) ||
      !scalar(value$attestation_id, "^attest_[0-9a-f]{64}$") ||
      !scalar(value$contract_hash, "^[0-9a-f]{64}$") ||
      !scalar(value$policy_id, "^policy_[0-9a-f]{64}$") ||
      !scalar(value$alignment_purpose,
              "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$") ||
      !scalar(value$dataset_id,
              "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$") ||
      !scalar(value$dataset_version,
              "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$") ||
      !scalar(value$id_column,
              "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$") ||
      !scalar(value$source_binding_id, "^source_[0-9a-f]{64}$") ||
      !scalar(value$pinset_id, "^pinset_[0-9a-f]{64}$") ||
      !integer(value$capacity_bucket, 64L, 1048576L) ||
      bitwAnd(as.integer(value$capacity_bucket),
              as.integer(value$capacity_bucket) - 1L) != 0L ||
      !integer(value$relay_frame_bytes, 16L * 1024L, 64L * 1024L^2) ||
      !integer(value$inline_max_bytes, 16L * 1024L, 64L * 1024L^2) ||
      !integer(value$peer_count, 2L, 1000000L) ||
      !scalar(value$reference_peer,
              "^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$") ||
      !is.character(value$compute_peers) ||
      length(value$compute_peers) != 2L || anyNA(value$compute_peers) ||
      anyDuplicated(value$compute_peers) ||
      any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$",
                 value$compute_peers, perl = TRUE)) ||
      !value$reference_peer %in% value$compute_peers ||
      value$peer_count < length(value$compute_peers)) invalid()
  normalized <- value
  for (field in c(
      "attestation_version", "capacity_bucket", "relay_frame_bytes",
      "inline_max_bytes", "peer_count")) {
    normalized[[field]] <- as.integer(value[[field]])
  }
  if (!is.null(expected_contract)) {
    expected <- list(
      attestation_version = 3L,
      alignment_attested = TRUE,
      alignment_protocol = .DSVERT_CLIENT_PSI_PADDED_PROTOCOL,
      attestation_id = expected_contract$attestation_id,
      contract_hash = expected_contract$contract_hash,
      policy_id = expected_contract$policy_id,
      alignment_purpose = expected_contract$alignment_purpose,
      dataset_id = expected_contract$dataset_id,
      dataset_version = expected_contract$dataset_version,
      id_column = expected_contract$id_column,
      source_binding_id = expected_contract$source_binding_id,
      pinset_id = expected_contract$pinset_id,
      capacity_bucket = as.integer(expected_contract$capacity),
      relay_frame_bytes = as.integer(expected_contract$relay_frame_bytes),
      inline_max_bytes = as.integer(expected_contract$inline_max_bytes),
      peer_count = length(expected_contract$peer_names),
      reference_peer = expected_contract$reference_peer,
      compute_peers = expected_contract$compute_peers)
    if (!identical(normalized, expected)) invalid()
  }
  normalized
}

.dsvert_psi_padded_json_b64url <- function(value) {
  .dsvert_exact_gc_b64url_encode(charToRaw(as.character(jsonlite::toJSON(
    value, auto_unbox = TRUE, null = "null", na = "null", digits = 17,
    pretty = FALSE))))
}

.dsvert_psi_padded_contract <- function(bound, server_names, session_id,
                                        operation_id) {
  invalid <- function() stop(
    "Pinned padded PSI peers did not agree on one complete contract.",
    call. = FALSE)
  if (!is.list(bound) || !identical(names(bound), server_names)) invalid()
  contracts <- lapply(bound, `[[`, "contract")
  if (any(vapply(contracts, is.null, logical(1L))) ||
      !all(vapply(contracts[-1L], identical, logical(1L), contracts[[1L]]))) {
    invalid()
  }
  contract <- contracts[[1L]]
  required <- c(
    "protocol", "session_id", "operation_id", "policy_id",
    "alignment_purpose", "dataset_id", "dataset_version", "id_column",
    "source_binding_id", "pinset_id",
    "capacity", "relay_frame_bytes", "inline_max_bytes", "peer_names",
    "peer_ids", "reference_peer",
    "compute_peers", "snapshot_ids", "attestation_id", "contract_hash")
  scalar <- function(value, pattern) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      grepl(pattern, value)
  }
  capacity <- suppressWarnings(as.numeric(contract$capacity))
  if (!is.list(contract) || !identical(names(contract), required) ||
      !identical(contract$protocol, .DSVERT_CLIENT_PSI_PADDED_PROTOCOL) ||
      !identical(contract$session_id, session_id) ||
      !identical(contract$operation_id, operation_id) ||
      !scalar(contract$policy_id, "^policy_[0-9a-f]{64}$") ||
      !scalar(contract$alignment_purpose,
              "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$") ||
      !scalar(contract$dataset_id,
              "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$") ||
      !scalar(contract$dataset_version,
              "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$") ||
      !scalar(contract$id_column,
              "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$") ||
      !scalar(contract$source_binding_id, "^source_[0-9a-f]{64}$") ||
      !scalar(contract$pinset_id, "^pinset_[0-9a-f]{64}$") ||
      !scalar(contract$attestation_id, "^attest_[0-9a-f]{64}$") ||
      !scalar(contract$contract_hash, "^[0-9a-f]{64}$") ||
      length(capacity) != 1L || is.na(capacity) || !is.finite(capacity) ||
      capacity != floor(capacity) || capacity < 64 ||
      capacity > 1048576 || bitwAnd(as.integer(capacity),
                                   as.integer(capacity) - 1L) != 0L ||
      !is.numeric(contract$relay_frame_bytes) ||
      length(contract$relay_frame_bytes) != 1L ||
      is.na(contract$relay_frame_bytes) ||
      contract$relay_frame_bytes < 16L * 1024L ||
      contract$relay_frame_bytes > 64L * 1024L^2 ||
      contract$relay_frame_bytes != floor(contract$relay_frame_bytes) ||
      !is.numeric(contract$inline_max_bytes) ||
      length(contract$inline_max_bytes) != 1L ||
      is.na(contract$inline_max_bytes) ||
      contract$inline_max_bytes < 16L * 1024L ||
      contract$inline_max_bytes > 64L * 1024L^2 ||
      contract$inline_max_bytes != floor(contract$inline_max_bytes) ||
      !is.character(contract$peer_names) ||
      length(contract$peer_names) != length(server_names) ||
      anyNA(contract$peer_names) || anyDuplicated(contract$peer_names) ||
      !setequal(contract$peer_names, server_names) ||
      !is.character(contract$peer_ids) ||
      length(contract$peer_ids) != length(server_names) ||
      anyNA(contract$peer_ids) || anyDuplicated(contract$peer_ids) ||
      any(!grepl("^dsv1_[0-9a-f]{64}$", contract$peer_ids)) ||
      !identical(contract$reference_peer, contract$peer_names[[1L]]) ||
      !identical(contract$compute_peers, contract$peer_names[1:2]) ||
      !is.character(contract$snapshot_ids) ||
      length(contract$snapshot_ids) != length(server_names) ||
      anyNA(contract$snapshot_ids) || anyDuplicated(contract$snapshot_ids) ||
      any(!grepl("^snap_[0-9a-f]{64}$", contract$snapshot_ids))) invalid()
  receipts <- lapply(bound, `[[`, "receipt")
  names(receipts) <- server_names
  list(contract = contract, receipts = receipts)
}

.dsvert_psi_padded_validate_ack <- function(values, field, expected = TRUE) {
  if (!is.list(values) || !length(values) || any(!vapply(
      values, function(value) is.list(value) &&
        identical(value[[field]], expected), logical(1L)))) {
    stop("A pinned padded PSI peer rejected a fixed protocol phase.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_psi_padded_align <- function(
    data_name, id_col, newobj, datasources,
    session_id = .dsvert_uuid4(),
    operation_id = .dsvert_exact_gc_new_context()$operation_id,
    timeout_seconds = getOption("dsvert.exact_gc.timeout_seconds", 900),
    .aggregate = DSI::datashield.aggregate,
    .assign = DSI::datashield.assign.expr,
    .exact_run = .dsvert_exact_gc_run) {
  server_names <- names(datasources)
  if (!is.list(datasources) || length(datasources) < 2L ||
      is.null(server_names) || anyNA(server_names) ||
      any(!nzchar(server_names)) || anyDuplicated(server_names)) {
    stop("Pinned padded PSI requires at least two uniquely named connections.",
         call. = FALSE)
  }
  data_name <- .dsvert_site_character(data_name, datasources, "data_name")
  id_col <- .dsvert_site_character(id_col, datasources, "id_col")
  if (!is.character(newobj) || length(newobj) != 1L || is.na(newobj) ||
      !nzchar(newobj)) {
    stop("Pinned padded PSI destination name must be one non-empty string.",
         call. = FALSE)
  }
  if (!is.function(.aggregate) || !is.function(.assign) ||
      !is.function(.exact_run)) {
    stop("Pinned padded PSI transport functions are unavailable.",
         call. = FALSE)
  }
  relay_controller <- NULL

  transport_call <- function(expression, envelope,
                             relay_descriptor_b64url) {
    if (nzchar(envelope)) expression$envelope <- envelope
    if (nzchar(relay_descriptor_b64url)) {
      expression$relay_descriptor_b64url <- relay_descriptor_b64url
    }
    expression
  }

  offers <- .dsvert_fanout_by_site(
    datasources, stats::setNames(lapply(server_names, function(server) {
      call(
        name = "psiPaddedInitDS", data_name = data_name[[server]],
        id_col = id_col[[server]], session_id = session_id,
        operation_id = operation_id)
    }), server_names),
    operation = "pinned padded PSI initialization", .aggregate = .aggregate)
  encoded_offers <- .dsvert_psi_padded_json_b64url(unname(offers))
  encoded_names <- .dsvert_psi_padded_json_b64url(unname(server_names))
  bound <- .dsvert_aggregate_strict(
    datasources, call(
      name = "psiPaddedBindDS", offers_b64url = encoded_offers,
      connection_names_b64url = encoded_names, session_id = session_id),
    operation = "pinned padded PSI contract binding", .aggregate = .aggregate)
  agreement <- .dsvert_psi_padded_contract(
    bound, server_names, session_id, operation_id)
  contract <- agreement$contract
  relay_controller <- .dsvert_psi_relay_controller(session_id, contract)
  force_relay <- .dsvert_dsi_max_expression_bytes() <
    as.numeric(contract$inline_max_bytes) + 16L * 1024L

  confirmed <- .dsvert_aggregate_strict(
    datasources, call(
      name = "psiPaddedConfirmDS",
      receipts_b64url = .dsvert_psi_padded_json_b64url(agreement$receipts),
      session_id = session_id),
    operation = "pinned padded PSI all-peer confirmation",
    .aggregate = .aggregate)
  .dsvert_psi_padded_validate_ack(confirmed, "confirmed")
  prepared <- .dsvert_fanout_by_site(
    datasources, stats::setNames(lapply(server_names, function(server) {
      call(
        name = "psiPaddedPrepareDS", data_name = data_name[[server]],
        session_id = session_id)
    }), server_names),
    operation = "pinned padded PSI fixed-capacity preparation",
    .aggregate = .aggregate)
  .dsvert_psi_padded_validate_ack(prepared, "prepared")

  reference <- contract$reference_peer
  targets <- setdiff(contract$peer_names, reference)
  for (target in targets) {
    reference_export <- .dsvert_aggregate_strict(
      datasources[reference], call(
        name = "psiPaddedReferenceExportDS", target = target,
        session_id = session_id, force_relay = force_relay),
      operation = "pinned padded PSI reference exchange",
      .aggregate = .aggregate)[[1L]]
    reference_delivery <- .dsvert_psi_relay_deliver(
      reference_export, reference, target, contract, datasources,
      session_id, relay_controller, .aggregate)
    target_export <- .dsvert_aggregate_strict(
      datasources[target], transport_call(call(
        name = "psiPaddedTargetProcessDS",
        session_id = session_id, force_relay = force_relay),
        reference_delivery$envelope,
        reference_delivery$relay_descriptor_b64url),
      operation = "pinned padded PSI target exchange",
      .aggregate = .aggregate)[[1L]]
    target_delivery <- .dsvert_psi_relay_deliver(
      target_export, target, reference, contract, datasources,
      session_id, relay_controller, .aggregate)
    doubled <- .dsvert_aggregate_strict(
      datasources[reference], transport_call(call(
        name = "psiPaddedReferenceDoubleDS", target = target,
        session_id = session_id, force_relay = force_relay),
        target_delivery$envelope,
        target_delivery$relay_descriptor_b64url),
      operation = "pinned padded PSI double masking",
      .aggregate = .aggregate)[[1L]]
    doubled_delivery <- .dsvert_psi_relay_deliver(
      doubled, reference, target, contract, datasources,
      session_id, relay_controller, .aggregate)
    membership <- .dsvert_aggregate_strict(
      datasources[target], transport_call(call(
        name = "psiPaddedTargetMatchDS",
        session_id = session_id, force_relay = force_relay),
        doubled_delivery$envelope,
        doubled_delivery$relay_descriptor_b64url),
      operation = "pinned padded PSI private membership sharing",
      .aggregate = .aggregate)[[1L]]
    membership_deliveries <- lapply(
      contract$compute_peers, function(peer) {
        .dsvert_psi_relay_deliver(
          membership$transports[[peer]], target, peer, contract,
          datasources, session_id, relay_controller, .aggregate)
      })
    names(membership_deliveries) <- contract$compute_peers
    requests <- stats::setNames(lapply(contract$compute_peers, function(peer) {
      transport_call(call(
        name = "psiPaddedMembershipAcceptDS", target = target,
        session_id = session_id),
        membership_deliveries[[peer]]$envelope,
        membership_deliveries[[peer]]$relay_descriptor_b64url)
    }), contract$compute_peers)
    accepted <- .dsvert_fanout_by_site(
      datasources[contract$compute_peers], requests,
      operation = "pinned padded PSI membership delivery",
      .aggregate = .aggregate)
    .dsvert_psi_padded_validate_ack(accepted, "accepted")
  }

  compute_indices <- match(contract$compute_peers, server_names)
  chunk_count <- as.integer(ceiling(as.numeric(contract$capacity) / 4096L))
  for (chunk_index in seq_len(chunk_count)) {
    initialized <- .dsvert_aggregate_strict(
      datasources[contract$compute_peers], call(
        name = "psiPaddedANDStartDS", chunk_index = chunk_index,
        session_id = session_id),
      operation = "pinned padded PSI exact AND start",
      .aggregate = .aggregate)
    exemplar <- initialized[[1L]]
    suffix <- sub("^op_", "", exemplar$operation_id)
    .exact_run(
      datasources = datasources, server_names = server_names,
      servers = compute_indices, session_id = session_id,
      operation_id = exemplar$operation_id,
      source_key = paste0("exact_gc_in_", suffix),
      output_key = paste0("exact_gc_out_", suffix),
      operation = "compare-signed", ring = 63L, frac_bits = 0L,
      vector_len = exemplar$vector_len, purpose = exemplar$purpose,
      transport_ready = TRUE, timeout_seconds = timeout_seconds,
      initialized = initialized, .aggregate = .aggregate)
    outputs <- .dsvert_aggregate_strict(
      datasources[contract$compute_peers], call(
        name = "psiPaddedANDFinalizeDS", chunk_index = chunk_index,
        session_id = session_id, force_relay = force_relay),
      operation = "pinned padded PSI exact AND finalization",
      .aggregate = .aggregate)
    for (peer in contract$compute_peers) {
      output_delivery <- .dsvert_psi_relay_deliver(
        outputs[[peer]], peer, reference, contract, datasources,
        session_id, relay_controller, .aggregate)
      accepted <- .dsvert_aggregate_strict(
        datasources[reference], transport_call(call(
          name = "psiPaddedANDAcceptDS", chunk_index = chunk_index,
          sender = peer, session_id = session_id),
          output_delivery$envelope,
          output_delivery$relay_descriptor_b64url),
        operation = "pinned padded PSI exact AND delivery",
        .aggregate = .aggregate)
      .dsvert_psi_padded_validate_ack(accepted, "accepted")
    }
  }

  final <- .dsvert_aggregate_strict(
    datasources[reference], call(
      name = "psiPaddedFinalPrepareDS", session_id = session_id,
      force_relay = force_relay),
    operation = "pinned padded PSI final selection",
    .aggregate = .aggregate)[[1L]]
  final_deliveries <- lapply(targets, function(target) {
    .dsvert_psi_relay_deliver(
      final$transports[[target]], reference, target, contract, datasources,
      session_id, relay_controller, .aggregate)
  })
  names(final_deliveries) <- targets
  assignments <- stats::setNames(lapply(server_names, function(peer) {
    delivery <- if (identical(peer, reference)) {
      list(envelope = "", relay_descriptor_b64url = "")
    } else {
      final_deliveries[[peer]]
    }
    transport_call(call(
      name = "psiPaddedFilterDS", data_name = data_name[[peer]],
      session_id = session_id),
      delivery$envelope, delivery$relay_descriptor_b64url)
  }), server_names)
  .dsvert_assign_strict(
    datasources, newobj, assignments,
    operation = "pinned padded PSI aligned-object commit", .assign = .assign)
  attestations <- .dsvert_aggregate_strict(
    datasources, call(
      name = "psiPaddedAttestationDS", data_name = newobj,
      session_id = session_id),
    operation = "pinned padded PSI attestation", .aggregate = .aggregate)
  attestations <- lapply(
    attestations, .dsvert_validate_psi_padded_attestation,
    expected_contract = contract)
  if (!all(vapply(attestations[-1L], identical, logical(1L),
                  attestations[[1L]]))) {
    stop("Pinned padded PSI peers returned inconsistent attestations.",
         call. = FALSE)
  }
  .dsvert_psi_relay_journal_delete(relay_controller)
  invisible(attestations[[1L]])
}
