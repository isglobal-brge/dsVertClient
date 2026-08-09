# Client pump for producer-minted, purpose-bound server blob tickets.

.DSVERT_TYPED_BLOB_FRAME_OVERHEAD_BYTES <- 1024L

.dsvert_typed_blob_resource_model <- function() {
  list(
    version = "dsvert-typed-blob-resource-model-v4",
    producer_source_streaming = FALSE,
    client_source_streaming = FALSE,
    producer_capabilities = list(
      mpcTypedSourceProbeDS = list(
        producer_source_streaming = TRUE,
        client_source_streaming = TRUE,
        statistical_producer = FALSE,
        recipient_consumer_streaming = FALSE),
      biomedical_capsule_source_v2 = list(
        endpoints = c(
          "dsvertDPCapsuleSourcePrepareDS",
          "dsvertDPCapsuleSourceChunkDS",
          "dsvertDPCapsuleSourceAcceptDS"),
        producer_source_streaming = TRUE,
        client_source_streaming = TRUE,
        statistical_producer = TRUE,
        recipient_consumer_streaming = TRUE,
        migrated_scope =
          "cross-owner Gaussian/categorical private tail coordinates",
        producer_incremental_peak_memory = "O(C + chunk) in migrated scope",
        release_prefix_peak_memory = "O(release + C*p_local)",
        full_statistical_materializer_streaming = FALSE,
        client_incremental_peak_memory = "O(chunk)",
        prepare_builds_all_ciphertexts = FALSE,
        emitted_chunk_durable_replay = TRUE,
        unmaterialized_requires_same_snapshot = TRUE,
        bounded_full_transfer_disk_reservation = TRUE)),
    recipient_disk_spool = TRUE,
    recipient_spool_aggregate_bounded = TRUE,
    producer_peak_memory = "O(payload)",
    client_peak_memory = "O(payload)",
    recipient_transport_peak_memory = "O(frame)",
    recipient_consumer_peak_memory = "O(payload)",
    max_payload_chars = 512 * 1024^2,
    blocker = paste(
      "Statistical producer methods and .callMpcTool return one complete Base64 value;",
      "the DSI aggregate response and client hash contract materialize it",
      "before immutable frames can be forwarded, and current consumers read",
      "the committed value in full. True end-to-end streaming needs a private",
      "producer spool and streaming MPC runtime input/output per producer.",
      "Only the data-free source probe currently proves the spool/reader/client",
      "layers for generic blobs. The separately declared biomedical-capsule",
      "source v2 is the sole statistical chunk producer/consumer exception."))
}

.dsvert_dsi_retry_deadline_seconds <- function() {
  value <- getOption(
    "dsvert.dsi.retry_deadline_seconds",
    getOption("dsvert.dsi.timeout_seconds", 300))
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value <= 0 || value > 24 * 60 * 60) {
    stop("DSI retry deadline must be between zero and 24 hours",
         call. = FALSE)
  }
  value
}

.dsvert_monotonic_seconds <- function() {
  unname(proc.time()[["elapsed"]])
}

.dsvert_retry_sleep <- function(seconds) Sys.sleep(seconds)

.dsvert_retry_jitter <- function() stats::runif(1L, 0.75, 1.25)

.dsvert_retry_deadline_condition <- function(operation) {
  structure(
    list(
      message = paste0(
        "DataSHIELD transport remained unavailable until the retry ",
        "deadline during '", operation,
        "'; the idempotent request remains resumable."),
      call = NULL, code = "retry_deadline_exceeded",
      operation = operation),
    class = c(
      "dsvert_retry_deadline_exceeded", "retry_deadline_exceeded",
      "dsvert_transport_unavailable", "transport_unavailable",
      "error", "condition"))
}

# Retry only an unavailable/ambiguous response.  There is deliberately no
# attempt counter: availability is governed by a monotonic deadline, while the
# immutable request bytes and server-side idempotency key never change.
.dsvert_retry_idempotent <- function(
    attempt, classify, operation,
    timeout_seconds = .dsvert_dsi_retry_deadline_seconds(),
    .clock = .dsvert_monotonic_seconds, .sleep = .dsvert_retry_sleep,
    .jitter = .dsvert_retry_jitter) {
  if (!is.function(attempt) || !is.function(classify) ||
      !is.character(operation) || length(operation) != 1L ||
      is.na(operation) || !nzchar(operation) || !is.function(.clock) ||
      !is.function(.sleep) || !is.function(.jitter)) {
    stop("Invalid idempotent DSI retry contract", call. = FALSE)
  }
  timeout_seconds <- suppressWarnings(as.numeric(timeout_seconds))
  if (length(timeout_seconds) != 1L || is.na(timeout_seconds) ||
      !is.finite(timeout_seconds) || timeout_seconds <= 0 ||
      timeout_seconds > 24 * 60 * 60) {
    stop("Invalid idempotent DSI retry deadline", call. = FALSE)
  }
  started <- .clock()
  if (!is.numeric(started) || length(started) != 1L || is.na(started) ||
      !is.finite(started)) {
    stop("Monotonic DSI retry clock is unavailable", call. = FALSE)
  }
  deadline <- started + timeout_seconds
  unavailable <- 0L
  repeat {
    response <- attempt()
    classified <- classify(response)
    if (!is.list(classified) ||
        !is.character(classified$state) ||
        length(classified$state) != 1L || is.na(classified$state)) {
      stop("Invalid idempotent DSI retry classifier", call. = FALSE)
    }
    if (!identical(classified$state, "missing")) return(classified)
    unavailable <- unavailable + 1L
    now <- .clock()
    if (!is.numeric(now) || length(now) != 1L || is.na(now) ||
        !is.finite(now) || now < started) {
      stop("Monotonic DSI retry clock moved backwards", call. = FALSE)
    }
    if (now >= deadline) stop(.dsvert_retry_deadline_condition(operation))
    jitter <- .jitter()
    if (!is.numeric(jitter) || length(jitter) != 1L || is.na(jitter) ||
        !is.finite(jitter) || jitter < 0.5 || jitter > 1.5) {
      stop("Invalid data-independent DSI retry jitter", call. = FALSE)
    }
    exponent <- min(unavailable - 1L, 6L)
    delay <- min(1, 0.025 * 2^exponent) * jitter
    .sleep(min(delay, deadline - now))
  }
}

# Run one single-connection typed-blob DSI attempt.  Availability failures are
# returned as NULL so the immutable request can be replayed. A present callback
# is terminal unless it carries the fixed public capacity-backpressure tag;
# identity, policy, signature, and other contract failures therefore never
# enter the retry loop or become a generic transport timeout.
.dsvert_typed_blob_dsi_attempt <- function(
    conn, expression, .aggregate = DSI::datashield.aggregate) {
  callback_failed <- FALSE
  callback_backpressure <- FALSE
  peer_rejection <- NULL
  resource_oversize <- NULL
  response <- tryCatch(
    .dsvert_transport_aggregate(
      .aggregate = .aggregate,
      conns = conn, expr = expression, async = TRUE,
      error = function(site, message) {
        callback_failed <<- TRUE
        parsed <- .dsvert_client_parse_peer_not_recognized(message)
        if (is.null(peer_rejection) && !is.null(parsed)) {
          peer_rejection <<- parsed
        }
        oversized <- .dsvert_client_parse_resource_oversize(message)
        if (is.null(resource_oversize) && !is.null(oversized)) {
          resource_oversize <<- oversized
        }
        if (!is.null(.dsvert_client_parse_resource_backpressure(message))) {
          callback_backpressure <<- TRUE
        }
        invisible(NULL)
      }, errors.print = FALSE),
    interrupt = function(e) stop(e),
    error = function(e) {
      if (inherits(e, "dsvert_dsi_poisoned_session")) stop(e)
      candidate <- if (inherits(e, "dsvert_peer_not_recognized")) {
        e
      } else {
        .dsvert_client_parse_peer_not_recognized(conditionMessage(e))
      }
      if (!inherits(peer_rejection, "dsvert_peer_not_recognized") &&
          inherits(candidate, "dsvert_peer_not_recognized")) {
        peer_rejection <<- candidate
      }
      oversized <- if (inherits(e, "dsvert_resource_oversize")) {
        e
      } else {
        .dsvert_client_parse_resource_oversize(conditionMessage(e))
      }
      if (is.null(resource_oversize) &&
          inherits(oversized, "dsvert_resource_oversize")) {
        resource_oversize <<- oversized
      }
      NULL
    })
  if (inherits(peer_rejection, "dsvert_peer_not_recognized")) {
    stop(peer_rejection)
  }
  if (inherits(resource_oversize, "dsvert_resource_oversize")) {
    stop(resource_oversize)
  }
  if (isTRUE(callback_failed)) {
    if (isTRUE(callback_backpressure)) return(NULL)
    stop(.dsvert_client_remote_contract_failure(
      "typed-blob single-peer transport"))
  }
  response
}

# Base64url frames and canonical UUIDs need no R-string escaping.  For this
# fixed call schema, the serialized expression is therefore the scalar byte
# lengths plus a small constant syntax overhead.  Keeping a 1 KiB reserve is
# deliberately conservative (the current expression uses fewer than 128
# bytes) and avoids deparsing every large frame solely to count it.
.dsvert_validate_typed_blob_frame_size <- function(
    ticket, chunk, offset, session_id) {
  scalar_ascii <- function(value, pattern, max_bytes) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value) && nchar(value, type = "bytes") <= max_bytes &&
      identical(nchar(value, type = "chars"),
                nchar(value, type = "bytes")) &&
      grepl(pattern, value, perl = TRUE, useBytes = TRUE)
  }
  scalar_bytes <- function(value) {
    if (is.character(value) && length(value) == 1L && !is.na(value)) {
      nchar(value, type = "bytes")
    } else {
      NA_real_
    }
  }
  ticket_bytes <- scalar_bytes(ticket)
  chunk_bytes <- scalar_bytes(chunk)
  if (!is.na(ticket_bytes) && ticket_bytes > 32L * 1024L) {
    stop(.dsvert_client_resource_oversize(
      requested_bytes = ticket_bytes, capacity_bytes = 32L * 1024L,
      scope = "typed-blob ticket"))
  }
  if (!is.na(chunk_bytes) &&
      chunk_bytes > .DSVERT_DSI_PROBE_ABSOLUTE_MAX) {
    stop(.dsvert_client_resource_oversize(
      requested_bytes = chunk_bytes,
      capacity_bytes = .DSVERT_DSI_PROBE_ABSOLUTE_MAX,
      scope = "typed-blob immutable frame"))
  }
  if (!scalar_ascii(ticket, "^[A-Za-z0-9_-]+$", 32L * 1024L) ||
      !scalar_ascii(chunk, "^[A-Za-z0-9_-]+$",
                    .DSVERT_DSI_PROBE_ABSOLUTE_MAX) ||
      !scalar_ascii(
        session_id,
        "^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$",
        36L) ||
      !is.numeric(offset) || length(offset) != 1L || is.na(offset) ||
      !is.finite(offset) || offset < 0 || offset != floor(offset) ||
      offset > 512 * 1024^2) {
    stop("Invalid typed-blob frame expression", call. = FALSE)
  }
  offset_text <- format(offset, scientific = FALSE, trim = TRUE)
  size_bound <- sum(
    nchar(c(ticket, chunk, session_id, offset_text), type = "bytes")) +
    .DSVERT_TYPED_BLOB_FRAME_OVERHEAD_BYTES
  capacity <- .dsvert_dsi_max_expression_bytes()
  if (size_bound > capacity) {
    stop(.dsvert_client_resource_oversize(
      requested_bytes = size_bound, capacity_bytes = capacity,
      scope = "typed-blob DataSHIELD expression"))
  }
  invisible(as.numeric(size_bound))
}

.dsvert_validate_typed_blob_contract <- function(transfer, target) {
  required <- c(
    "ticket", "transfer_id", "capability_id", "sender_name", "recipient_name",
    "payload_chars", "payload_sha256")
  if (!is.list(transfer) || is.null(names(transfer)) ||
      anyDuplicated(names(transfer)) ||
      !identical(sort(names(transfer)), sort(required))) {
    stop("Producer returned an invalid typed-blob transfer contract",
         call. = FALSE)
  }
  scalar <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value)
  }
  if (!scalar(transfer$ticket) ||
      nchar(transfer$ticket, type = "bytes") > 32L * 1024L ||
      !grepl("^[A-Za-z0-9_-]+$", transfer$ticket,
             perl = TRUE, useBytes = TRUE) ||
      !scalar(transfer$transfer_id) ||
      !grepl("^tb_[0-9a-f]{32}$", transfer$transfer_id) ||
      !scalar(transfer$capability_id) ||
      !grepl("^[a-z][a-z0-9]*(?:[.-][a-z0-9]+)*$",
             transfer$capability_id, perl = TRUE) ||
      !scalar(transfer$sender_name) ||
      identical(transfer$sender_name, transfer$recipient_name) ||
      !scalar(transfer$recipient_name) ||
      !identical(transfer$recipient_name, target) ||
      !scalar(transfer$payload_sha256) ||
      !grepl("^[0-9a-f]{64}$", transfer$payload_sha256)) {
    stop("Producer returned an invalid or misrouted typed-blob transfer",
         call. = FALSE)
  }
  expected_chars <- suppressWarnings(as.numeric(transfer$payload_chars))
  if (length(expected_chars) != 1L || is.na(expected_chars) ||
      !is.finite(expected_chars) || expected_chars != floor(expected_chars) ||
      expected_chars < 1 || expected_chars > 512 * 1024^2) {
    stop("Producer returned an invalid typed-blob payload descriptor",
         call. = FALSE)
  }
  transfer$payload_chars <- as.numeric(expected_chars)
  transfer
}

.dsvert_validate_typed_blob_transfer <- function(transfer, blob, target) {
  transfer <- .dsvert_validate_typed_blob_contract(transfer, target)
  if (!is.character(blob) || length(blob) != 1L || is.na(blob) ||
      !nzchar(blob) ||
      !grepl("^[A-Za-z0-9_-]+$", blob, perl = TRUE, useBytes = TRUE) ||
      nchar(blob) %% 4L == 1L ||
      !identical(nchar(blob, type = "chars"),
                 nchar(blob, type = "bytes"))) {
    stop("Typed-blob payload must be non-empty canonical Base64url text",
         call. = FALSE)
  }
  expected_chars <- transfer$payload_chars
  actual_chars <- nchar(blob, type = "bytes")
  if (length(expected_chars) != 1L || is.na(expected_chars) ||
      !is.finite(expected_chars) || expected_chars != floor(expected_chars) ||
      expected_chars < 1 || expected_chars != actual_chars ||
      !identical(paste0(openssl::sha256(charToRaw(blob))),
                 transfer$payload_sha256)) {
    stop("Typed-blob payload does not match its producer contract",
         call. = FALSE)
  }
  transfer$payload_chars <- as.numeric(expected_chars)
  transfer
}

.dsvert_typed_blob_ack_state <- function(response, target, transfer,
                                          expected_committed,
                                          expected_sealed) {
  missing <- list(state = "missing", receipt = NULL,
                  committed_chars = NULL, sealed = NULL)
  invalid <- list(state = "invalid", receipt = NULL,
                  committed_chars = NULL, sealed = NULL)
  if (is.null(response)) return(missing)
  if (!is.list(response) || length(response) != 1L ||
      is.null(names(response)) || !identical(names(response), target)) {
    return(invalid)
  }
  ack <- response[[1L]]
  rejection_fields <- c("version", "operation", "rejected")
  if (is.list(ack) && !is.null(names(ack)) &&
      !anyDuplicated(names(ack)) &&
      identical(sort(names(ack)), sort(rejection_fields)) &&
      identical(ack$version, "dsvert-typed-blob-rejection-v1") &&
      identical(ack$operation, "store") && identical(ack$rejected, TRUE)) {
    return(list(state = "rejected", receipt = NULL,
                committed_chars = NULL, sealed = NULL))
  }
  required <- c(
    "version", "transfer_id", "committed_chars", "total_chars", "sealed",
    "receipt")
  if (is.null(ack)) return(missing)
  if (!is.list(ack) || is.null(names(ack)) || anyDuplicated(names(ack)) ||
      !identical(sort(names(ack)), sort(required)) ||
      !identical(ack$version, "dsvert-typed-blob-v1") ||
      !identical(ack$transfer_id, transfer$transfer_id) ||
      !is.numeric(ack$committed_chars) ||
      length(ack$committed_chars) != 1L || is.na(ack$committed_chars) ||
      !is.finite(ack$committed_chars) ||
      !is.numeric(ack$total_chars) || length(ack$total_chars) != 1L ||
      is.na(ack$total_chars) || !is.finite(ack$total_chars) ||
      as.numeric(ack$total_chars) != transfer$payload_chars ||
      !is.logical(ack$sealed) || length(ack$sealed) != 1L ||
      is.na(ack$sealed)) {
    return(invalid)
  }
  committed <- as.numeric(ack$committed_chars)
  receipt_valid <- is.character(ack$receipt) &&
    length(ack$receipt) == 1L && !is.na(ack$receipt) &&
    nzchar(ack$receipt) &&
    nchar(ack$receipt, type = "bytes") <= 32L * 1024L &&
    grepl("^[A-Za-z0-9_-]+$", ack$receipt,
          perl = TRUE, useBytes = TRUE)
  completed <- identical(ack$sealed, TRUE) &&
    committed == transfer$payload_chars && receipt_valid
  in_progress <- identical(ack$sealed, FALSE) && is.null(ack$receipt) &&
    committed >= expected_committed &&
    committed < transfer$payload_chars
  if (!completed && !in_progress) return(invalid)
  list(
    state = if (completed && !isTRUE(expected_sealed)) "complete" else
      if (in_progress && committed > expected_committed) "ahead" else "ack",
    receipt = ack$receipt, committed_chars = committed,
    sealed = isTRUE(ack$sealed))
}

.dsvert_typed_blob_receipt_ack_state <- function(response, producer,
                                                  transfer) {
  if (is.null(response)) return(list(state = "missing"))
  if (!is.list(response) || length(response) != 1L ||
      is.null(names(response)) || !identical(names(response), producer)) {
    return(list(state = "invalid"))
  }
  ack <- response[[1L]]
  rejection_fields <- c("version", "operation", "rejected")
  if (is.list(ack) && !is.null(names(ack)) &&
      !anyDuplicated(names(ack)) &&
      identical(sort(names(ack)), sort(rejection_fields)) &&
      identical(ack$version, "dsvert-typed-blob-rejection-v1") &&
      identical(ack$operation, "receipt") && identical(ack$rejected, TRUE)) {
    return(list(state = "rejected"))
  }
  required <- c("version", "transfer_id", "confirmed")
  if (is.null(ack)) return(list(state = "missing"))
  if (!is.list(ack) || is.null(names(ack)) || anyDuplicated(names(ack)) ||
      !identical(sort(names(ack)), sort(required)) ||
      !identical(ack$version, "dsvert-typed-blob-receipt-v1") ||
      !identical(ack$transfer_id, transfer$transfer_id) ||
      !identical(ack$confirmed, TRUE)) return(list(state = "invalid"))
  list(state = "ack")
}

.dsvert_confirm_typed_blob_receipt <- function(
    receipt, transfer, producer_conn, session_id, .aggregate) {
  if (length(producer_conn) != 1L || is.null(names(producer_conn)) ||
      !identical(names(producer_conn), transfer$sender_name)) {
    stop("Typed-blob producer connection does not match its signed contract",
         call. = FALSE)
  }
  .dsvert_validate_real_dsi_transport(producer_conn, .aggregate)
  expression <- call(
    name = "mpcTypedBlobReceiptDS", receipt = receipt,
    session_id = session_id)
  .dsvert_validate_dsi_expression_sizes(expression)
  attempt <- function() {
    .dsvert_typed_blob_dsi_attempt(producer_conn, expression, .aggregate)
  }
  outcome <- .dsvert_retry_idempotent(
    attempt = attempt,
    classify = function(response) .dsvert_typed_blob_receipt_ack_state(
      response, transfer$sender_name, transfer),
    operation = "typed-blob signed receipt confirmation")
  if (identical(outcome$state, "invalid")) {
    stop("Typed-blob producer returned a malformed receipt acknowledgement",
         call. = FALSE)
  }
  if (identical(outcome$state, "rejected")) {
    stop("Typed-blob producer rejected the signed peer receipt",
         call. = FALSE)
  }
  if (!identical(outcome$state, "ack")) {
    stop("Typed-blob producer did not acknowledge the signed peer receipt",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_validate_typed_blob_read_size <- function(
    ticket, offset, max_chars, session_id) {
  if (!is.character(ticket) || length(ticket) != 1L || is.na(ticket) ||
      !nzchar(ticket) || nchar(ticket, type = "bytes") > 32L * 1024L ||
      !grepl("^[A-Za-z0-9_-]+$", ticket, perl = TRUE, useBytes = TRUE) ||
      !is.numeric(offset) || length(offset) != 1L || is.na(offset) ||
      !is.finite(offset) || offset < 0 || offset != floor(offset) ||
      !is.numeric(max_chars) || length(max_chars) != 1L ||
      is.na(max_chars) || !is.finite(max_chars) || max_chars < 1 ||
      max_chars != floor(max_chars) || max_chars > 8 * 1024^2 ||
      !is.character(session_id) || length(session_id) != 1L ||
      is.na(session_id) ||
      !grepl("^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$",
             session_id)) {
    stop("Invalid typed-source read expression", call. = FALSE)
  }
  expression <- call(
    name = "mpcTypedBlobReadDS", ticket = ticket, offset = offset,
    max_chars = max_chars, session_id = session_id)
  .dsvert_validate_dsi_expression_sizes(expression)
  invisible(expression)
}

.dsvert_typed_blob_source_state <- function(
    response, producer, transfer, expected_offset, max_chars) {
  if (is.null(response)) return(list(state = "missing"))
  if (!is.list(response) || length(response) != 1L ||
      is.null(names(response)) || !identical(names(response), producer)) {
    return(list(state = "invalid"))
  }
  value <- response[[1L]]
  rejection_fields <- c("version", "operation", "rejected")
  if (is.list(value) && !is.null(names(value)) &&
      !anyDuplicated(names(value)) &&
      identical(sort(names(value)), sort(rejection_fields)) &&
      identical(value$version, "dsvert-typed-blob-rejection-v1") &&
      identical(value$operation, "read") &&
      identical(value$rejected, TRUE)) {
    return(list(state = "rejected"))
  }
  required <- c(
    "version", "transfer_id", "offset", "chunk", "chunk_chars",
    "chunk_sha256", "total_chars", "payload_sha256", "final")
  if (is.null(value)) return(list(state = "missing"))
  if (!is.list(value) || is.null(names(value)) ||
      anyDuplicated(names(value)) ||
      !identical(sort(names(value)), sort(required)) ||
      !identical(value$version, "dsvert-typed-blob-source-v1") ||
      !identical(value$transfer_id, transfer$transfer_id) ||
      !identical(as.numeric(value$offset), as.numeric(expected_offset)) ||
      !is.character(value$chunk) || length(value$chunk) != 1L ||
      is.na(value$chunk) || !nzchar(value$chunk) ||
      !grepl("^[A-Za-z0-9_-]+$", value$chunk,
             perl = TRUE, useBytes = TRUE) ||
      !identical(nchar(value$chunk, type = "chars"),
                 nchar(value$chunk, type = "bytes")) ||
      !is.character(value$chunk_sha256) ||
      length(value$chunk_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", value$chunk_sha256) ||
      !identical(value$payload_sha256, transfer$payload_sha256) ||
      !identical(as.numeric(value$total_chars), transfer$payload_chars) ||
      !is.logical(value$final) || length(value$final) != 1L ||
      is.na(value$final)) {
    return(list(state = "invalid"))
  }
  chunk_chars <- nchar(value$chunk, type = "bytes")
  declared_chars <- suppressWarnings(as.numeric(value$chunk_chars))
  expected_chars <- min(max_chars, transfer$payload_chars - expected_offset)
  expected_final <- identical(
    as.numeric(expected_offset + expected_chars), transfer$payload_chars)
  if (length(declared_chars) != 1L || is.na(declared_chars) ||
      !is.finite(declared_chars) || declared_chars != chunk_chars ||
      chunk_chars != expected_chars ||
      !identical(value$final, expected_final) ||
      !identical(paste0(openssl::sha256(charToRaw(value$chunk))),
                 value$chunk_sha256)) {
    return(list(state = "invalid"))
  }
  list(state = "ack", chunk = value$chunk,
       chunk_chars = as.numeric(chunk_chars), final = value$final)
}

.dsvert_read_typed_blob_source_frame <- function(
    transfer, producer_conn, offset, max_chars, session_id, .aggregate) {
  producer <- transfer$sender_name
  expression <- .dsvert_validate_typed_blob_read_size(
    transfer$ticket, offset, max_chars, session_id)
  attempt <- function() {
    .dsvert_typed_blob_dsi_attempt(producer_conn, expression, .aggregate)
  }
  outcome <- .dsvert_retry_idempotent(
    attempt = attempt,
    classify = function(response) .dsvert_typed_blob_source_state(
      response, producer, transfer, offset, max_chars),
    operation = "typed-blob immutable source frame read")
  if (identical(outcome$state, "rejected")) {
    stop("Typed-blob producer rejected the purpose-bound source read",
         call. = FALSE)
  }
  if (!identical(outcome$state, "ack")) {
    stop("Typed-blob producer returned a malformed source frame",
         call. = FALSE)
  }
  outcome
}

# Commit one frame while fetching the next immutable producer frame in the
# same named DSI fan-out. A partial outer response replays the byte-identical
# pair and advances neither offset; a present malformed response is fatal.
.dsvert_typed_blob_pipeline_step <- function(
    transfer, producer_conn, recipient_conn, store_expression,
    next_offset, max_chars, session_id, .aggregate) {
  producer <- transfer$sender_name
  target <- transfer$recipient_name
  paths <- c(producer_conn, recipient_conn)
  next_expression <- .dsvert_validate_typed_blob_read_size(
    transfer$ticket, next_offset, max_chars, session_id)
  expressions <- stats::setNames(
    list(next_expression, store_expression), c(producer, target))
  expected_committed <- as.numeric(next_offset)

  attempt <- function() {
    callback_failed <- character()
    callback_backpressure <- character()
    peer_rejection <- NULL
    resource_oversize <- NULL
    response <- tryCatch(
      .dsvert_transport_aggregate(
        .aggregate = .aggregate, conns = paths, expr = expressions,
        async = TRUE,
        error = function(site, message) {
          site <- if (length(site)) as.character(site[[1L]]) else "unknown"
          callback_failed <<- union(callback_failed, site)
          parsed <- .dsvert_client_parse_peer_not_recognized(message)
          if (is.null(peer_rejection) && !is.null(parsed)) {
            peer_rejection <<- parsed
          }
          oversized <- .dsvert_client_parse_resource_oversize(message)
          if (is.null(resource_oversize) && !is.null(oversized)) {
            resource_oversize <<- oversized
          }
          if (!is.null(.dsvert_client_parse_resource_backpressure(message))) {
            callback_backpressure <<- union(callback_backpressure, site)
          }
          invisible(NULL)
        }, errors.print = FALSE),
      interrupt = function(e) stop(e),
      error = function(e) {
        if (inherits(e, "dsvert_dsi_poisoned_session")) stop(e)
        candidate <- if (inherits(e, "dsvert_peer_not_recognized")) {
          e
        } else {
          .dsvert_client_parse_peer_not_recognized(conditionMessage(e))
        }
        if (!inherits(peer_rejection, "dsvert_peer_not_recognized") &&
            inherits(candidate, "dsvert_peer_not_recognized")) {
          peer_rejection <<- candidate
        }
        oversized <- if (inherits(e, "dsvert_resource_oversize")) {
          e
        } else {
          .dsvert_client_parse_resource_oversize(conditionMessage(e))
        }
        if (is.null(resource_oversize) &&
            inherits(oversized, "dsvert_resource_oversize")) {
          resource_oversize <<- oversized
        }
        NULL
      })
    if (inherits(peer_rejection, "dsvert_peer_not_recognized")) {
      stop(peer_rejection)
    }
    if (inherits(resource_oversize, "dsvert_resource_oversize")) {
      stop(resource_oversize)
    }
    if (length(callback_failed)) {
      if (setequal(callback_failed, callback_backpressure)) return(NULL)
      stop(.dsvert_client_remote_contract_failure(
        "typed-blob pipelined transport"))
    }
    response
  }
  classify <- function(response) {
    if (is.null(response)) return(list(state = "missing"))
    named <- is.list(response) && length(response) == 2L &&
      !is.null(names(response)) && !anyNA(names(response)) &&
      !anyDuplicated(names(response)) && setequal(names(response), names(paths))
    if (!named) return(list(state = "invalid"))
    response <- response[names(paths)]
    if (any(vapply(response, is.null, logical(1L)))) {
      return(list(state = "missing"))
    }
    sink <- .dsvert_typed_blob_ack_state(
      response[target], target, transfer, expected_committed,
      expected_sealed = FALSE)
    source <- .dsvert_typed_blob_source_state(
      response[producer], producer, transfer, next_offset, max_chars)
    states <- c(sink$state, source$state)
    if (any(states == "missing")) return(list(state = "missing"))
    if (any(states == "rejected")) return(list(state = "rejected"))
    if (!sink$state %in% c("ack", "ahead", "complete") ||
        !identical(source$state, "ack")) {
      return(list(state = "invalid"))
    }
    list(state = "ack", sink = sink, source = source)
  }
  outcome <- .dsvert_retry_idempotent(
    attempt = attempt, classify = classify,
    operation = "typed-blob pipelined frame commit and source read")
  if (identical(outcome$state, "rejected")) {
    stop("Typed-blob peer rejected the pipelined transfer step",
         call. = FALSE)
  }
  if (!identical(outcome$state, "ack")) {
    stop("Typed-blob peer returned a malformed pipelined transfer step",
         call. = FALSE)
  }
  outcome
}

#' Run the data-free typed source-stream diagnostic (internal)
#'
#' This is the only client constructor for the source-stream pilot. It asks one
#' pinned producer to create random diagnostic bytes for one pinned recipient,
#' then pumps the purpose-bound ticket with the same immutable frame protocol
#' used by typed payloads. It never names a server path and is not a statistical
#' producer.
#'
#' @keywords internal
.dsvert_run_typed_source_probe <- function(
    producer_conn, recipient_conn, recipient_pk, payload_bytes, session_id,
    .aggregate = DSI::datashield.aggregate) {
  named_peer <- function(conn) {
    is.list(conn) && length(conn) == 1L && !is.null(names(conn)) &&
      is.character(names(conn)) && length(names(conn)) == 1L &&
      !is.na(names(conn)) && nzchar(names(conn))
  }
  if (!named_peer(producer_conn) || !named_peer(recipient_conn) ||
      identical(names(producer_conn), names(recipient_conn))) {
    stop("Source probe requires two distinct named pinned peers",
         call. = FALSE)
  }
  if (!is.character(recipient_pk) || length(recipient_pk) != 1L ||
      is.na(recipient_pk) || !nzchar(recipient_pk) ||
      nchar(recipient_pk, type = "bytes") > 32L * 1024L ||
      !grepl("^[A-Za-z0-9_-]+$", recipient_pk,
             perl = TRUE, useBytes = TRUE)) {
    stop("Source probe requires one canonical pinned recipient key",
         call. = FALSE)
  }
  payload_bytes <- suppressWarnings(as.numeric(payload_bytes))
  if (length(payload_bytes) != 1L || is.na(payload_bytes) ||
      !is.finite(payload_bytes) || payload_bytes != floor(payload_bytes) ||
      payload_bytes < 1 || payload_bytes > 384 * 1024^2) {
    stop("Source probe byte count must be one integer from 1 to 384 MiB",
         call. = FALSE)
  }
  if (!is.character(session_id) || length(session_id) != 1L ||
      is.na(session_id) ||
      !grepl("^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$",
             session_id)) {
    stop("Source probe requires one canonical session identifier",
         call. = FALSE)
  }

  produced <- .dsvert_aggregate_strict(
    producer_conn,
    call(
      name = "mpcTypedSourceProbeDS", recipient_pk = recipient_pk,
      payload_bytes = payload_bytes, session_id = session_id),
    operation = "data-free typed source-stream probe",
    .aggregate = .aggregate)
  value <- produced[[1L]]
  if (!is.list(value) || is.null(names(value)) || anyDuplicated(names(value)) ||
      !identical(names(value), "source_transfer")) {
    stop("Source probe producer returned an invalid transfer envelope",
         call. = FALSE)
  }
  transfer <- .dsvert_validate_typed_blob_contract(
    value$source_transfer, names(recipient_conn))
  if (!identical(transfer$sender_name, names(producer_conn)) ||
      !identical(transfer$capability_id,
                 "blob.transport.source-probe.v1")) {
    stop("Source probe producer returned a misbound transfer contract",
         call. = FALSE)
  }
  frames <- .dsvert_store_typed_blob_stream(
    transfer = transfer, conn = recipient_conn, session_id = session_id,
    producer_conn = producer_conn, .aggregate = .aggregate)
  invisible(list(transfer = transfer, frames = as.integer(frames)))
}

#' Pump one producer-spooled typed payload without materialising it (internal)
#'
#' Source and destination use immutable absolute offsets. The client verifies
#' every frame immediately. The recipient verifies the complete ticket-bound
#' SHA-256 before atomic commit and signs a receipt binding that hash and length;
#' the producer verifies that receipt before releasing its source spool. Thus
#' the untrusted relay needs no second O(payload) integrity spool and its peak
#' memory and disk use are both O(frame).
#'
#' @keywords internal
.dsvert_store_typed_blob_stream <- function(
    transfer, conn, session_id, producer_conn,
    .aggregate = DSI::datashield.aggregate) {
  if (length(conn) != 1L || is.null(names(conn)) ||
      !is.character(names(conn)) || is.na(names(conn)) ||
      !nzchar(names(conn))) {
    stop("conn must contain exactly one named DataSHIELD target",
         call. = FALSE)
  }
  target <- names(conn)
  transfer <- .dsvert_validate_typed_blob_contract(transfer, target)
  if (length(producer_conn) != 1L || is.null(names(producer_conn)) ||
      !identical(names(producer_conn), transfer$sender_name)) {
    stop("producer_conn must name the signed typed-blob producer",
         call. = FALSE)
  }
  paths <- c(producer_conn, conn)
  if (anyDuplicated(names(paths))) {
    stop("Typed source and recipient must be distinct pinned peers",
         call. = FALSE)
  }
  .dsvert_maybe_negotiate_dsi_chunk_size(paths, .aggregate)
  .dsvert_validate_real_dsi_transport(paths, .aggregate)
  chunk_size <- .dsvert_get_chunk_size()
  frame_count <- .dsvert_chunk_count(transfer$payload_chars, chunk_size)
  if (frame_count > 4096L) {
    stop(.dsvert_client_resource_oversize(
      requested_bytes = frame_count, capacity_bytes = 4096L,
      scope = "typed-blob frame metadata"))
  }

  .dsvert_chunk_env$geometry_locked <- TRUE
  on.exit({
    .dsvert_chunk_env$geometry_locked <- FALSE
  }, add = TRUE)

  final_receipt <- NULL
  source_offset <- 0
  sink_committed <- 0
  source <- .dsvert_read_typed_blob_source_frame(
    transfer, producer_conn, source_offset, chunk_size,
    session_id, .aggregate)
  while (source_offset < transfer$payload_chars) {
    chunk <- source$chunk
    chunk_end <- source_offset + source$chunk_chars
    next_source <- NULL

    if (source_offset >= sink_committed) {
      if (!identical(source_offset, sink_committed)) {
        stop("Typed-blob recipient acknowledgement skipped source geometry",
             call. = FALSE)
      }
      expected_sealed <- identical(chunk_end, transfer$payload_chars)
      expr <- call(
        name = "mpcTypedBlobStoreDS", ticket = transfer$ticket,
        chunk = chunk, offset = source_offset, session_id = session_id)
      .dsvert_validate_typed_blob_frame_size(
        transfer$ticket, chunk, source_offset, session_id)
      if (!expected_sealed) {
        step <- .dsvert_typed_blob_pipeline_step(
          transfer, producer_conn, conn, expr, chunk_end, chunk_size,
          session_id, .aggregate)
        acknowledgement <- step$sink
        next_source <- step$source
      } else {
        attempt <- function() {
          .dsvert_typed_blob_dsi_attempt(conn, expr, .aggregate)
        }
        acknowledgement <- .dsvert_retry_idempotent(
          attempt = attempt,
          classify = function(response) .dsvert_typed_blob_ack_state(
            response, target, transfer, chunk_end, expected_sealed),
          operation = "typed-blob immutable streamed frame commit")
      }
      if (!acknowledgement$state %in% c("ack", "ahead", "complete")) {
        stop("Typed-blob recipient rejected or malformed a streamed frame",
             call. = FALSE)
      }
      sink_committed <- acknowledgement$committed_chars
      if (!is.numeric(sink_committed) || length(sink_committed) != 1L ||
          is.na(sink_committed) || !is.finite(sink_committed) ||
          sink_committed < chunk_end || sink_committed > transfer$payload_chars ||
          (sink_committed < transfer$payload_chars &&
           sink_committed %% chunk_size != 0)) {
        stop("Typed-blob recipient returned an impossible streamed offset",
             call. = FALSE)
      }
      if (isTRUE(acknowledgement$sealed)) {
        final_receipt <- acknowledgement$receipt
        source_offset <- transfer$payload_chars
        break
      }
    }
    source_offset <- chunk_end
    if (source_offset < transfer$payload_chars) {
      if (!is.null(next_source)) {
        source <- next_source
      } else {
        source <- .dsvert_read_typed_blob_source_frame(
          transfer, producer_conn, source_offset, chunk_size,
          session_id, .aggregate)
      }
    }
  }
  if (!identical(sink_committed, transfer$payload_chars) ||
      !is.character(final_receipt) || length(final_receipt) != 1L ||
      is.na(final_receipt) || !nzchar(final_receipt)) {
    stop("Typed-blob recipient did not commit the complete streamed source",
         call. = FALSE)
  }
  .dsvert_confirm_typed_blob_receipt(
    final_receipt, transfer, producer_conn, session_id, .aggregate)
  .dsvert_set_chunk_size(chunk_size)
  invisible(frame_count)
}

#' Send one producer-ticketed opaque payload to its pinned recipient (internal)
#'
#' The call intentionally has no key, filename, recipient or purpose override.
#' Ambiguous DSI responses replay the exact same request until a monotonic
#' availability deadline. The server commits by signed transfer ID, absolute
#' offset and payload hash; no request-count or rate quota is consulted.
#'
#' @keywords internal
.dsvert_store_typed_blob <- function(
    blob, transfer, conn, session_id, producer_conn,
    .aggregate = DSI::datashield.aggregate) {
  if (length(conn) != 1L || is.null(names(conn)) ||
      !is.character(names(conn)) || is.na(names(conn)) ||
      !nzchar(names(conn))) {
    stop("conn must contain exactly one named DataSHIELD target",
         call. = FALSE)
  }
  target <- names(conn)
  .dsvert_maybe_negotiate_dsi_chunk_size(conn, .aggregate)
  transfer <- .dsvert_validate_typed_blob_transfer(transfer, blob, target)
  if (missing(producer_conn) || length(producer_conn) != 1L ||
      is.null(names(producer_conn)) ||
      !identical(names(producer_conn), transfer$sender_name)) {
    stop("producer_conn must name the signed typed-blob producer",
         call. = FALSE)
  }
  .dsvert_validate_real_dsi_transport(conn, .aggregate)
  chunk_size <- .dsvert_get_chunk_size()
  frame_count <- .dsvert_chunk_count(transfer$payload_chars, chunk_size)
  if (frame_count > 4096L) {
    stop(.dsvert_client_resource_oversize(
      requested_bytes = frame_count, capacity_bytes = 4096L,
      scope = "typed-blob frame metadata"))
  }
  .dsvert_chunk_env$geometry_locked <- TRUE
  on.exit({
    .dsvert_chunk_env$geometry_locked <- FALSE
  }, add = TRUE)

  final_receipt <- NULL
  offset <- 0
  while (offset < transfer$payload_chars) {
    end <- min(offset + chunk_size, transfer$payload_chars)
    chunk <- substr(blob, offset + 1, end)
    expected_committed <- as.numeric(end)
    expected_sealed <- identical(expected_committed, transfer$payload_chars)
    expr <- call(
      name = "mpcTypedBlobStoreDS", ticket = transfer$ticket,
      chunk = chunk, offset = offset, session_id = session_id)
    .dsvert_validate_typed_blob_frame_size(
      transfer$ticket, chunk, offset, session_id)

    attempt <- function() {
      .dsvert_typed_blob_dsi_attempt(conn, expr, .aggregate)
    }
    acknowledgement <- .dsvert_retry_idempotent(
      attempt = attempt,
      classify = function(response) .dsvert_typed_blob_ack_state(
        response, target, transfer, expected_committed, expected_sealed),
      operation = "typed-blob immutable frame commit")
    if (identical(acknowledgement$state, "invalid")) {
      stop("Typed-blob store returned a malformed acknowledgement for target '",
           target, "'", call. = FALSE)
    }
    if (identical(acknowledgement$state, "rejected")) {
      stop("Typed-blob store rejected the immutable frame for target '",
           target, "'", call. = FALSE)
    }
    if (!acknowledgement$state %in% c("ack", "ahead", "complete")) {
      stop("Typed-blob frame was not acknowledged by target '", target,
           "'", call. = FALSE)
    }
    committed <- acknowledgement$committed_chars
    if (!is.numeric(committed) || length(committed) != 1L ||
        is.na(committed) || !is.finite(committed) ||
        committed != floor(committed) || committed < expected_committed ||
        committed > transfer$payload_chars ||
        (committed < transfer$payload_chars && committed %% chunk_size != 0)) {
      stop("Typed-blob store returned an impossible committed offset for target '",
           target, "'", call. = FALSE)
    }
    if (isTRUE(acknowledgement$sealed)) {
      final_receipt <- acknowledgement$receipt
      offset <- transfer$payload_chars
    } else {
      offset <- committed
    }
  }
  if (!is.character(final_receipt) || length(final_receipt) != 1L ||
      is.na(final_receipt) || !nzchar(final_receipt)) {
    stop("Typed-blob recipient did not return its signed commit receipt",
         call. = FALSE)
  }
  .dsvert_confirm_typed_blob_receipt(
    final_receipt, transfer, producer_conn, session_id, .aggregate)
  .dsvert_set_chunk_size(chunk_size)
  invisible(frame_count)
}

# Transitional dispatcher. Active-first callers pass a producer contract and
# therefore have no key/slot argument. Non-migrated second-wave and quarantine
# callers may still pass a legacy key until their own protocol migration.
.dsvert_store_transfer_or_legacy <- function(
    blob, contract, conn, session_id, producer_conns = NULL,
    .aggregate = DSI::datashield.aggregate) {
  if (is.list(contract)) {
    if (is.null(producer_conns) ||
        is.null(producer_conns[[contract$sender_name]])) {
      stop("Typed-blob relay requires the producer connection set",
           call. = FALSE)
    }
    return(.dsvert_store_typed_blob(
      blob, contract, conn, session_id,
      producer_conn = producer_conns[contract$sender_name],
      .aggregate = .aggregate))
  }
  if (!is.character(contract) || length(contract) != 1L ||
      is.na(contract) || !nzchar(contract)) {
    stop("A producer transfer contract or legacy key is required",
         call. = FALSE)
  }
  .dsvert_store_blob(
    blob, contract, conn, session_id, .aggregate = .aggregate)
}
