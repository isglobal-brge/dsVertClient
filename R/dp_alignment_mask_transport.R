# Exact private-alignment gate for cross-owner capsule inputs.
#
# The analyst relay transports only authenticated exact-GC frames.  Each
# compute peer keeps one XOR share of every source custodian's alignment
# digest; the digests and mismatch location are reconstructed only inside the
# fixed-shape circuit.  One uniform terminal outcome is opened after all
# public chunks have completed.

.DSVERT_CLIENT_DP_ALIGNMENT_MASK_OPERATION <- "alignment-mask-ring128"
.DSVERT_CLIENT_DP_ALIGNMENT_MASK_OUTPUT <-
  "alignment-masked-ring128-share-v1"
.DSVERT_CLIENT_DP_ALIGNMENT_MASK_TERMINAL_VERSION <-
  "dsvert-alignment-mask-terminal-v1"

.dsvert_dp_alignment_mask_chunk_size_client <- function(source_count) {
  source_count <- as.integer(.dsvert_exact_gc_number(
    source_count, "alignment source count", 2, 64))
  size <- floor((
    .DSVERT_CLIENT_EXACT_GC_MAX_CIRCUIT_TYPE_BITS / 128L -
      4L * source_count - 1L) / 3L)
  if (!is.finite(size) || size < 1L) {
    stop("The private alignment-mask circuit shape is not representable",
         call. = FALSE)
  }
  as.integer(min(size, 4096L))
}

.dsvert_dp_alignment_mask_operation_client <- function(
    batch_operation_id, contract_hash, chunk_index, chunk_count) {
  if (!is.character(batch_operation_id) ||
      length(batch_operation_id) != 1L || is.na(batch_operation_id) ||
      !grepl("^op_[0-9a-f]{32}$", batch_operation_id) ||
      !.dsvert_dp_capsule_source_hex(contract_hash) ||
      !.dsvert_dp_is_integer(chunk_index, 1, 2^31 - 1) ||
      !.dsvert_dp_is_integer(chunk_count, 1, 2^31 - 1)) {
    stop("Invalid private alignment-mask chunk contract", call. = FALSE)
  }
  digest <- digest::digest(charToRaw(paste0(
    "dsVert/alignment-mask/ring128/chunk/v1|", batch_operation_id, "|",
    contract_hash, "|", as.integer(chunk_index), "|",
    as.integer(chunk_count))), algo = "sha256", serialize = FALSE)
  paste0("op_", substr(digest, 1L, 32L))
}

.dsvert_dp_alignment_mask_keys_client <- function(operation_id) {
  if (!is.character(operation_id) || length(operation_id) != 1L ||
      is.na(operation_id) || !grepl("^op_[0-9a-f]{32}$", operation_id)) {
    stop("Invalid private alignment-mask operation", call. = FALSE)
  }
  suffix <- sub("^op_", "", operation_id)
  list(source = paste0("exact_gc_in_", suffix),
       output = paste0("exact_gc_out_", suffix))
}

.dsvert_dp_alignment_mask_terminal_error <- function(source_count,
                                                      coordinate_count) {
  structure(list(
    message = paste0(
      "The private cross-owner alignment contract is invalid; no source ",
      "digest, predicate, or offending custodian was released."),
    call = NULL, code = "alignment_contract_invalid", retryable = FALSE,
    source_count = as.integer(source_count),
    coordinate_count = as.numeric(coordinate_count),
    mismatch_source_exposed = FALSE, alignment_digest_exposed = FALSE),
    class = c("dsvert_alignment_contract_invalid", "error", "condition"))
}

.dsvert_dp_alignment_mask_terminal_set <- function(
    responses, peers, source_count, coordinate_count, chunk_count) {
  fields <- c(
    "capability_id", "version", "state", "terminal_outcome",
    "fixed_transcript", "source_count", "coordinate_count", "chunk_count",
    "alignment_digest_exposed", "mismatch_source_exposed",
    "gate_share_exposed")
  if (!is.list(responses) || !setequal(names(responses), peers)) {
    stop("The private alignment gate did not cover both computation peers",
         call. = FALSE)
  }
  values <- responses[peers]
  valid <- vapply(values, function(value) {
    .dsvert_dp_has_exact_names(value, fields) &&
      identical(value$capability_id,
                .DSVERT_CLIENT_EXACT_GC_CAPABILITY) &&
      identical(value$version,
                .DSVERT_CLIENT_DP_ALIGNMENT_MASK_TERMINAL_VERSION) &&
      value$state %in% c("complete", "alignment_contract_invalid") &&
      identical(value$terminal_outcome, value$state) &&
      identical(value$fixed_transcript, TRUE) &&
      identical(as.numeric(value$source_count), as.numeric(source_count)) &&
      identical(as.numeric(value$coordinate_count),
                as.numeric(coordinate_count)) &&
      identical(as.numeric(value$chunk_count), as.numeric(chunk_count)) &&
      identical(value$alignment_digest_exposed, FALSE) &&
      identical(value$mismatch_source_exposed, FALSE) &&
      identical(value$gate_share_exposed, FALSE)
  }, logical(1L))
  states <- vapply(values, function(value) {
    if (is.list(value) && is.character(value$state) &&
        length(value$state) == 1L) value$state else ""
  }, character(1L))
  if (!all(valid) || length(unique(states)) != 1L) {
    stop("The computation peers disagreed on the private alignment gate",
         call. = FALSE)
  }
  states[[1L]]
}

.dsvert_dp_alignment_mask_run <- function(
    manifest_json, context, layout, source_receipt, session_id, .aggregate,
    .run = .dsvert_exact_gc_run, .remote_context = NULL) {
  peers <- context$designated
  sources <- unname(unlist(layout$source_peers))
  source_count <- length(sources)
  total <- suppressWarnings(as.numeric(layout$transport_coordinate_count))
  contract_hash <- source_receipt$contract_hash
  if (length(peers) != 2L || anyDuplicated(peers) ||
      source_count < 2L || source_count > 64L || anyDuplicated(sources) ||
      length(total) != 1L || is.na(total) || !is.finite(total) ||
      total != floor(total) || total < 1 || total > 2^53 ||
      !.dsvert_dp_capsule_source_hex(contract_hash) ||
      !is.function(.run)) {
    stop("Invalid private alignment-mask orchestration contract",
         call. = FALSE)
  }
  chunk_size <- .dsvert_dp_alignment_mask_chunk_size_client(source_count)
  chunk_count <- as.integer(ceiling(total / chunk_size))
  batch <- .dsvert_exact_gc_new_context()$operation_id
  conns <- context$conns
  synopsis <- !is.null(.remote_context)
  if (isTRUE(synopsis)) {
    .remote_context <-
      .dsvert_dp_synopsis_categorical_cross_remote_context_v1(
        .remote_context)
  }
  start_method <- if (isTRUE(synopsis)) {
    as.name("dsvertDPSynopsisAlignmentMaskStartDS")
  } else as.name("dsvertDPAlignmentMaskStartDS")
  store_method <- if (isTRUE(synopsis)) {
    as.name("dsvertDPSynopsisAlignmentMaskStoreDS")
  } else as.name("dsvertDPAlignmentMaskStoreDS")
  seal_method <- if (isTRUE(synopsis)) {
    as.name("dsvertDPSynopsisAlignmentMaskSealDS")
  } else as.name("dsvertDPAlignmentMaskSealDS")
  receive_method <- if (isTRUE(synopsis)) {
    as.name("dsvertDPSynopsisAlignmentMaskReceiveDS")
  } else as.name("dsvertDPAlignmentMaskReceiveDS")
  endpoint_context <- if (isTRUE(synopsis)) .remote_context else
    list(manifest_json = manifest_json)
  cleanup_operations <- character()
  committed <- FALSE
  on.exit(if (!committed && length(cleanup_operations)) {
    for (operation_id in rev(cleanup_operations)) {
      .dsvert_exact_gc_abort(
        conns, session_id, operation_id, .aggregate)
    }
  }, add = TRUE)

  for (index in seq_len(chunk_count)) {
    operation_id <- .dsvert_dp_alignment_mask_operation_client(
      batch, contract_hash, index, chunk_count)
    cleanup_operations <- c(cleanup_operations, operation_id)
    keys <- .dsvert_dp_alignment_mask_keys_client(operation_id)
    n <- as.integer(min(chunk_size, total - (index - 1) * chunk_size))
    purpose <- paste0(
      "dp.alignment-mask.", substr(contract_hash, 1L, 20L),
      ".c-", index, "-", chunk_count)
    common <- c(endpoint_context, list(
      batch_operation_id = batch,
      operation_id = operation_id, chunk_index = as.integer(index),
      chunk_count = chunk_count, session_id = session_id))
    started <- .dsvert_fanout_by_site(
      conns, stats::setNames(lapply(peers, function(peer) {
        as.call(c(list(start_method), common))
      }), peers), operation = "private alignment-mask start",
      .aggregate = .aggregate)
    .run(
      context$all_conns, server_names = context$servers,
      servers = match(peers, context$servers), session_id = session_id,
      operation_id = operation_id, source_key = keys$source,
      output_key = keys$output,
      operation = .DSVERT_CLIENT_DP_ALIGNMENT_MASK_OPERATION,
      ring = 128L, frac_bits = 0L, vector_len = n, purpose = purpose,
      transport_ready = TRUE, initialized = started,
      alignment_source_count = source_count, .aggregate = .aggregate)
    stored <- .dsvert_fanout_by_site(
      conns, stats::setNames(lapply(peers, function(peer) {
        as.call(c(list(store_method), common))
      }), peers), operation = "private alignment-mask persistence",
      .aggregate = .aggregate)
    for (peer in peers) {
      value <- stored[[peer]]
      if (!is.list(value) ||
          !identical(value$capability_id,
                     .DSVERT_CLIENT_EXACT_GC_CAPABILITY) ||
          !identical(value$state, "stored") || !isTRUE(value$stored) ||
          !identical(as.numeric(value$chunk_index), as.numeric(index)) ||
          !identical(as.numeric(value$chunk_count),
                     as.numeric(chunk_count))) {
        stop("A computation peer failed to persist the private alignment mask",
             call. = FALSE)
      }
    }
  }

  sealed <- .dsvert_fanout_by_site(
    conns, stats::setNames(lapply(peers, function(peer) as.call(c(
      list(seal_method),
      endpoint_context,
      list(batch_operation_id = batch, session_id = session_id)))), peers),
    operation = "private alignment terminal sealing",
    .aggregate = .aggregate)
  for (peer in peers) {
    value <- sealed[[peer]]
    if (!is.list(value) ||
        !identical(value$capability_id,
                   .DSVERT_CLIENT_EXACT_GC_CAPABILITY) ||
        !identical(value$state, "sealed") ||
        !is.character(value$peer_blob) || length(value$peer_blob) != 1L ||
        is.na(value$peer_blob) ||
        !grepl("^[A-Za-z0-9_-]+$", value$peer_blob)) {
      stop("A computation peer failed to seal the private alignment outcome",
           call. = FALSE)
    }
  }
  receive_calls <- stats::setNames(lapply(seq_along(peers), function(i) {
    other <- setdiff(seq_along(peers), i)
    as.call(c(list(receive_method),
      endpoint_context,
      list(peer_blob = sealed[[peers[[other]]]]$peer_blob,
           batch_operation_id = batch, session_id = session_id)))
  }), peers)
  terminal <- .dsvert_fanout_by_site(
    conns, receive_calls, operation = "private alignment terminal opening",
    .aggregate = .aggregate)
  state <- .dsvert_dp_alignment_mask_terminal_set(
    terminal, peers, source_count, total, chunk_count)
  if (!identical(state, "complete")) {
    stop(.dsvert_dp_alignment_mask_terminal_error(source_count, total))
  }
  committed <- TRUE
  invisible(list(
    capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
    state = "complete", batch_operation_id = batch,
    contract_hash = contract_hash, source_count = as.integer(source_count),
    coordinate_count = total, chunk_count = chunk_count,
    alignment_digest_exposed = FALSE, mismatch_source_exposed = FALSE,
    gate_share_exposed = FALSE, fixed_transcript = TRUE))
}
