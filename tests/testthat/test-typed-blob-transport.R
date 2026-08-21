.typed_client_transfer <- function(blob, recipient = "peer", sender = "self") {
  list(
    ticket = paste(rep("A", 86L), collapse = ""),
    transfer_id = "tb_00112233445566778899aabbccddeeff",
    capability_id = "blob.input.peer-x.v1",
    sender_name = sender,
    recipient_name = recipient,
    payload_chars = as.numeric(nchar(blob, type = "bytes")),
    payload_sha256 = paste0(openssl::sha256(charToRaw(blob))))
}

.typed_client_ack <- function(expr, total) {
  offset <- as.numeric(expr$offset)
  committed <- offset + nchar(expr$chunk, type = "bytes")
  list(
    version = "dsvert-typed-blob-v1",
    transfer_id = "tb_00112233445566778899aabbccddeeff",
    committed_chars = committed,
    total_chars = as.numeric(total),
    sealed = identical(committed, as.numeric(total)),
    receipt = if (identical(committed, as.numeric(total))) "cmVjZWlwdA" else NULL)
}

.typed_client_receipt_ack <- function() {
  list(
    version = "dsvert-typed-blob-receipt-v1",
    transfer_id = "tb_00112233445566778899aabbccddeeff",
    confirmed = TRUE)
}

.typed_peer_error_token <- function(peer = "peer") {
  paste0(
    "DSVERT_PEER_NOT_RECOGNIZED_V1|peer=", peer,
    "|expected=", strrep("1", 64L),
    "|observed=", strrep("2", 64L))
}

test_that("typed transport does not advertise unimplemented source streaming", {
  model <- .dsvert_typed_blob_resource_model()
  expect_identical(
    model$version, "dsvert-typed-blob-resource-model-v5")
  expect_false(model$producer_source_streaming)
  expect_false(model$client_source_streaming)
  expect_true(model$producer_capabilities$mpcTypedSourceProbeDS$
                producer_source_streaming)
  expect_true(model$producer_capabilities$mpcTypedSourceProbeDS$
                client_source_streaming)
  expect_false(model$producer_capabilities$mpcTypedSourceProbeDS$
                 statistical_producer)
  expect_false(model$producer_capabilities$mpcTypedSourceProbeDS$
                 recipient_consumer_streaming)
  expect_false("formal_cox_blockwise_control_v1" %in%
                 names(model$producer_capabilities))
  capsule <- model$producer_capabilities$biomedical_capsule_source_v2
  expect_true(capsule$producer_source_streaming)
  expect_true(capsule$client_source_streaming)
  expect_true(capsule$statistical_producer)
  expect_true(capsule$recipient_consumer_streaming)
  expect_identical(capsule$producer_incremental_peak_memory,
                   "O(C + chunk) in migrated scope")
  expect_identical(capsule$release_prefix_peak_memory,
                   "O(release + C*p_local)")
  expect_false(capsule$full_statistical_materializer_streaming)
  expect_identical(capsule$client_incremental_peak_memory, "O(chunk)")
  expect_false(capsule$prepare_builds_all_ciphertexts)
  expect_true(capsule$emitted_chunk_durable_replay)
  expect_true(capsule$unmaterialized_requires_same_snapshot)
  expect_true(capsule$bounded_full_transfer_disk_reservation)
  expect_true(model$recipient_disk_spool)
  expect_true(model$recipient_spool_aggregate_bounded)
  expect_identical(model$producer_peak_memory, "O(payload)")
  expect_identical(model$client_peak_memory, "O(payload)")
  expect_identical(model$recipient_transport_peak_memory, "O(frame)")
  expect_identical(model$recipient_consumer_peak_memory, "O(payload)")
  expect_identical(model$max_payload_chars, 512 * 1024^2)
  expect_match(model$blocker, "complete Base64 value")
  expect_match(model$blocker, "streaming MPC runtime")
  expect_match(model$blocker, "biomedical-capsule")

  producer <- paste(deparse(body(.dsvert_validate_typed_blob_transfer)),
                    collapse = "\n")
  pump <- paste(deparse(body(.dsvert_store_typed_blob)), collapse = "\n")
  expect_match(producer, "charToRaw(blob)", fixed = TRUE)
  expect_match(pump, "substr(blob", fixed = TRUE)
})

test_that("source-stream client pumps fixed frames with O-frame state", {
  blob <- paste(rep(c("A", "b", "1", "_"), 41L), collapse = "")
  transfer <- .typed_client_transfer(blob, recipient = "recipient",
                                     sender = "producer")
  calls <- list()
  fanout_calls <- list()
  source_loss <- TRUE
  pipelined_source_loss <- TRUE
  sink_loss <- TRUE
  dispatch <- function(site, expr) {
    calls[[length(calls) + 1L]] <<- expr
    method <- as.character(expr[[1L]])
    if (identical(method, "mpcTypedBlobReadDS")) {
      offset <- as.numeric(expr$offset)
      max_chars <- as.numeric(expr$max_chars)
      end <- min(offset + max_chars, nchar(blob, type = "bytes"))
      chunk <- substr(blob, offset + 1, end)
      value <- list(
        version = "dsvert-typed-blob-source-v1",
        transfer_id = transfer$transfer_id, offset = offset, chunk = chunk,
        chunk_chars = as.numeric(nchar(chunk, type = "bytes")),
        chunk_sha256 = paste0(openssl::sha256(charToRaw(chunk))),
        total_chars = transfer$payload_chars,
        payload_sha256 = transfer$payload_sha256,
        final = identical(as.numeric(end), transfer$payload_chars))
      if (source_loss) {
        source_loss <<- FALSE
        return(NULL)
      }
      if (offset == 46 && pipelined_source_loss) {
        pipelined_source_loss <<- FALSE
        return(NULL)
      }
      return(value)
    }
    if (identical(method, "mpcTypedBlobStoreDS")) {
      value <- .typed_client_ack(expr, nchar(blob, type = "bytes"))
      if (sink_loss) {
        sink_loss <<- FALSE
        return(NULL)
      }
      return(value)
    }
    if (identical(method, "mpcTypedBlobReceiptDS")) {
      return(.typed_client_receipt_ack())
    }
    stop("Unexpected streamed transport method", call. = FALSE)
  }
  aggregate <- function(conns, expr, async, error, errors.print) {
    expressions <- if (is.list(expr) && !is.call(expr)) {
      methods <- vapply(expr, function(value) {
        as.character(value[[1L]])
      }, character(1L))
      if (setequal(methods,
                   c("mpcTypedBlobReadDS", "mpcTypedBlobStoreDS"))) {
        fanout_calls[[length(fanout_calls) + 1L]] <<- expr
      }
      expr[names(conns)]
    } else {
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    }
    stats::setNames(lapply(names(conns), function(site) {
      dispatch(site, expressions[[site]])
    }), names(conns))
  }
  withr::local_options(list(dsvert.chunk_size = 23L))
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  frames <- .dsvert_store_typed_blob_stream(
    transfer,
    conn = list(recipient = structure(list(), class = "test")),
    session_id = "12345678-1234-4234-8234-123456789abc",
    producer_conn = list(producer = structure(list(), class = "test")),
    .aggregate = aggregate)
  expect_identical(as.integer(frames),
                   as.integer(ceiling(nchar(blob) / 23)))
  source_calls <- calls[vapply(calls, function(expr) identical(
    as.character(expr[[1L]]), "mpcTypedBlobReadDS"), logical(1L))]
  expect_identical(source_calls[[1L]], source_calls[[2L]])
  expect_identical(
    unique(vapply(source_calls, function(expr)
      as.numeric(expr$max_chars), numeric(1L))), 23)
  expect_identical(sort(unique(vapply(source_calls, function(expr)
    as.numeric(expr$offset), numeric(1L)))),
    seq(0, nchar(blob) - 1, by = 23))
  sink_calls <- calls[vapply(calls, function(expr) identical(
    as.character(expr[[1L]]), "mpcTypedBlobStoreDS"), logical(1L))]
  expect_identical(sink_calls[[1L]], sink_calls[[2L]])
  expect_gte(length(fanout_calls), 2L)
  expect_identical(fanout_calls[[1L]], fanout_calls[[2L]])
  expect_identical(fanout_calls[[3L]], fanout_calls[[4L]])
  pair_hashes <- vapply(fanout_calls, digest::digest, character(1L),
                        algo = "sha256", serialize = TRUE)
  expect_identical(length(unique(pair_hashes)), as.integer(frames - 1L))
  expect_identical(length(fanout_calls), as.integer(frames + 1L))
  expect_true(all(vapply(fanout_calls, function(expressions) {
    identical(sort(unname(vapply(expressions, function(expr) {
      as.character(expr[[1L]])
    }, character(1L)))), sort(c(
      "mpcTypedBlobReadDS", "mpcTypedBlobStoreDS")))
  }, logical(1L))))
  expect_true(all(vapply(c(source_calls, sink_calls), function(expr) {
    is.null(expr$filename) && is.null(expr$key) && is.null(expr$purpose)
  }, logical(1L))))
  pump <- paste(deparse(body(.dsvert_store_typed_blob_stream)), collapse = "\n")
  expect_false(grepl("tempfile", pump, fixed = TRUE))
  expect_false(grepl("writeBin", pump, fixed = TRUE))
  expect_false(grepl("digest::digest", pump, fixed = TRUE))
})

test_that("source stream cannot confirm a receiver-side terminal hash failure", {
  blob <- strrep("A", 256 * 1024L)
  tampered <- paste0(substr(blob, 1L, nchar(blob) - 1L), "B")
  transfer <- .typed_client_transfer(
    blob, recipient = "recipient", sender = "producer")
  received <- ""
  receipt_calls <- 0L
  dispatch <- function(expr) {
    method <- as.character(expr[[1L]])
    if (identical(method, "mpcTypedBlobReadDS")) {
      offset <- as.numeric(expr$offset)
      end <- min(offset + as.numeric(expr$max_chars), nchar(tampered))
      chunk <- substr(tampered, offset + 1L, end)
      return(list(
        version = "dsvert-typed-blob-source-v1",
        transfer_id = transfer$transfer_id, offset = offset, chunk = chunk,
        chunk_chars = as.numeric(nchar(chunk)),
        chunk_sha256 = paste0(openssl::sha256(charToRaw(chunk))),
        total_chars = transfer$payload_chars,
        payload_sha256 = transfer$payload_sha256,
        final = identical(as.numeric(end), transfer$payload_chars)))
    }
    if (identical(method, "mpcTypedBlobStoreDS")) {
      expect_identical(as.numeric(expr$offset), as.numeric(nchar(received)))
      received <<- paste0(received, expr$chunk)
      if (nchar(received) == transfer$payload_chars &&
          !identical(paste0(openssl::sha256(charToRaw(received))),
                     transfer$payload_sha256)) {
        return(list(
          version = "dsvert-typed-blob-rejection-v1",
          operation = "store", rejected = TRUE))
      }
      return(.typed_client_ack(expr, nchar(tampered)))
    }
    if (identical(method, "mpcTypedBlobReceiptDS")) {
      receipt_calls <<- receipt_calls + 1L
      return(.typed_client_receipt_ack())
    }
    stop("Unexpected streamed transport method", call. = FALSE)
  }
  aggregate <- function(conns, expr, async, error, errors.print) {
    expressions <- if (is.list(expr) && !is.call(expr)) expr[names(conns)] else
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    stats::setNames(lapply(names(conns), function(site) {
      dispatch(expressions[[site]])
    }), names(conns))
  }
  withr::local_options(list(dsvert.chunk_size = 32 * 1024L))
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)

  expect_error(.dsvert_store_typed_blob_stream(
    transfer,
    conn = list(recipient = structure(list(), class = "test")),
    session_id = "12345678-1234-4234-8234-123456789abc",
    producer_conn = list(producer = structure(list(), class = "test")),
    .aggregate = aggregate), "rejected or malformed")
  expect_identical(receipt_calls, 0L)
})

test_that("source-probe constructor binds one producer to one recipient", {
  blob <- paste(rep(c("A", "b", "1", "_"), 13L), collapse = "")
  transfer <- .typed_client_transfer(
    blob, recipient = "recipient", sender = "producer")
  transfer$capability_id <- "blob.transport.source-probe.v1"
  calls <- list()
  dispatch <- function(conns, expr) {
    calls[[length(calls) + 1L]] <<- expr
    method <- as.character(expr[[1L]])
    if (identical(method, "mpcTypedSourceProbeDS")) {
      expect_identical(names(conns), "producer")
      expect_identical(expr$recipient_pk, "cGlubmVkLXJlY2lwaWVudA")
      return(list(source_transfer = transfer))
    }
    if (identical(method, "mpcTypedBlobReadDS")) {
      offset <- as.numeric(expr$offset)
      end <- min(offset + as.numeric(expr$max_chars), nchar(blob))
      chunk <- substr(blob, offset + 1, end)
      return(list(
        version = "dsvert-typed-blob-source-v1",
        transfer_id = transfer$transfer_id, offset = offset, chunk = chunk,
        chunk_chars = as.numeric(nchar(chunk)),
        chunk_sha256 = paste0(openssl::sha256(charToRaw(chunk))),
        total_chars = transfer$payload_chars,
        payload_sha256 = transfer$payload_sha256,
        final = identical(as.numeric(end), transfer$payload_chars)))
    }
    if (identical(method, "mpcTypedBlobStoreDS")) {
      return(.typed_client_ack(expr, nchar(blob)))
    }
    if (identical(method, "mpcTypedBlobReceiptDS")) {
      return(.typed_client_receipt_ack())
    }
    stop("Unexpected source-probe method", call. = FALSE)
  }
  aggregate <- function(conns, expr, async, error, errors.print) {
    expressions <- if (is.list(expr) && !is.call(expr)) expr[names(conns)] else
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    stats::setNames(lapply(names(conns), function(site) {
      dispatch(conns[site], expressions[[site]])
    }), names(conns))
  }
  withr::local_options(list(dsvert.chunk_size = 19L))
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  result <- .dsvert_run_typed_source_probe(
    producer_conn = list(producer = structure(list(), class = "test")),
    recipient_conn = list(recipient = structure(list(), class = "test")),
    recipient_pk = "cGlubmVkLXJlY2lwaWVudA", payload_bytes = 39,
    session_id = "12345678-1234-4234-8234-123456789abc",
    .aggregate = aggregate)
  expect_identical(result$transfer, transfer)
  expect_identical(result$frames, 3L)
  expect_identical(
    as.character(calls[[1L]][[1L]]), "mpcTypedSourceProbeDS")
  expect_true(all(vapply(calls, function(expr) {
    is.null(expr$filename) && is.null(expr$key) && is.null(expr$purpose)
  }, logical(1L))))
})

test_that("source-stream client rejects a bad frame hash before forwarding", {
  blob <- "AbCdEf0123_-"
  transfer <- .typed_client_transfer(
    blob, recipient = "recipient", sender = "producer")
  sink_calls <- 0L
  aggregate <- function(conns, expr, async, error, errors.print) {
    method <- as.character(expr[[1L]])
    if (identical(method, "mpcTypedBlobReadDS")) {
      return(list(producer = list(
        version = "dsvert-typed-blob-source-v1",
        transfer_id = transfer$transfer_id, offset = 0,
        chunk = blob, chunk_chars = as.numeric(nchar(blob)),
        chunk_sha256 = strrep("0", 64L),
        total_chars = transfer$payload_chars,
        payload_sha256 = transfer$payload_sha256, final = TRUE)))
    }
    sink_calls <<- sink_calls + 1L
    NULL
  }
  withr::local_options(list(dsvert.chunk_size = 64L))
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  expect_error(.dsvert_store_typed_blob_stream(
    transfer,
    conn = list(recipient = structure(list(), class = "test")),
    session_id = "12345678-1234-4234-8234-123456789abc",
    producer_conn = list(producer = structure(list(), class = "test")),
    .aggregate = aggregate), "malformed source frame")
  expect_identical(sink_calls, 0L)
})

test_that("typed frame size bound is conservative without deparsing payload", {
  ticket <- strrep("A", 32L * 1024L)
  chunk <- strrep("B", 640L * 1024L)
  sid <- "12345678-1234-4234-8234-123456789abc"
  offset <- 123456789
  expr <- call(
    name = "mpcTypedBlobStoreDS", ticket = ticket, chunk = chunk,
    offset = as.numeric(offset), session_id = sid)
  actual <- nchar(
    paste(deparse(expr, width.cutoff = 500L), collapse = "\n"),
    type = "bytes")

  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  withr::local_options(list(dsvert.dsi_max_expression_bytes = actual + 4096L))
  bound <- .dsvert_validate_typed_blob_frame_size(
    ticket, chunk, offset, sid)
  expect_gte(bound, actual)
  expect_lte(bound - actual,
             .DSVERT_TYPED_BLOB_FRAME_OVERHEAD_BYTES)

  withr::local_options(list(dsvert.dsi_max_expression_bytes = bound - 1L))
  oversized <- tryCatch(.dsvert_validate_typed_blob_frame_size(
    ticket, chunk, offset, sid), error = identity)
  expect_s3_class(oversized, "dsvert_resource_oversize")
  expect_identical(oversized$code, "resource_oversize")
  expect_false(oversized$retryable)
  expect_error(.dsvert_validate_typed_blob_frame_size(
    ticket, chunk, offset, paste0("method_", sid)),
    "Invalid typed-blob frame")
})

test_that("typed blob client uses immutable absolute-offset frames", {
  blob <- paste0(rep(c("A", "b", "1", "_"), 8L), collapse = "")
  transfer <- .typed_client_transfer(blob)
  calls <- list()
  aggregate <- function(conns, expr, async, error, errors.print) {
    calls[[length(calls) + 1L]] <<- expr
    expect_true(async)
    expect_false(errors.print)
    if (identical(as.character(expr[[1L]]), "mpcTypedBlobReceiptDS")) {
      expect_identical(names(conns), "self")
      return(list(self = .typed_client_receipt_ack()))
    }
    expect_identical(names(conns), "peer")
    list(peer = .typed_client_ack(expr, nchar(blob)))
  }
  withr::local_options(list(dsvert.chunk_size = 7L))
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  frames <- .dsvert_store_typed_blob(
    blob, transfer, list(peer = structure(list(), class = "test")),
    "12345678-1234-4234-8234-123456789abc",
    producer_conn = list(self = structure(list(), class = "test")),
    .aggregate = aggregate)
  expect_identical(as.integer(frames), 5L)
  frame_calls <- calls[vapply(calls, function(expr) {
    identical(as.character(expr[[1L]]), "mpcTypedBlobStoreDS")
  }, logical(1L))]
  expect_identical(
    vapply(frame_calls, function(expr) as.numeric(expr$offset), numeric(1L)),
    c(0, 7, 14, 21, 28))
  expect_identical(paste0(vapply(
    frame_calls, function(expr) as.character(expr$chunk), character(1L)),
    collapse = ""), blob)
  expect_true(all(vapply(frame_calls, function(expr) {
    identical(as.character(expr[[1L]]), "mpcTypedBlobStoreDS") &&
      is.null(expr$key) && is.null(expr$purpose) && is.null(expr$filename)
  }, logical(1L))))
})

test_that("typed blob client replays identical ambiguous frames and receipts", {
  blob <- "AbCdEf0123_-"
  transfer <- .typed_client_transfer(blob)
  calls <- list()
  receipt_attempts <- 0L
  aggregate <- function(conns, expr, async, error, errors.print) {
    calls[[length(calls) + 1L]] <<- expr
    if (identical(as.character(expr[[1L]]), "mpcTypedBlobReceiptDS")) {
      receipt_attempts <<- receipt_attempts + 1L
      if (receipt_attempts == 1L) return(NULL)
      return(list(self = .typed_client_receipt_ack()))
    }
    if (sum(vapply(calls, function(value) {
      identical(as.character(value[[1L]]), "mpcTypedBlobStoreDS")
    }, logical(1L))) == 1L) return(NULL)
    list(peer = .typed_client_ack(expr, nchar(blob)))
  }
  withr::local_options(list(dsvert.chunk_size = 64L))
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  expect_invisible(.dsvert_store_typed_blob(
    blob, transfer, list(peer = structure(list(), class = "test")),
    "12345678-1234-4234-8234-123456789abc",
    producer_conn = list(self = structure(list(), class = "test")),
    .aggregate = aggregate))
  frame_calls <- calls[vapply(calls, function(expr) {
    identical(as.character(expr[[1L]]), "mpcTypedBlobStoreDS")
  }, logical(1L))]
  expect_length(frame_calls, 2L)
  expect_identical(frame_calls[[1L]], frame_calls[[2L]])
  receipt_calls <- calls[vapply(calls, function(expr) {
    identical(as.character(expr[[1L]]), "mpcTypedBlobReceiptDS")
  }, logical(1L))]
  expect_length(receipt_calls, 2L)
  expect_identical(receipt_calls[[1L]], receipt_calls[[2L]])
})

test_that("idempotent retries have a deadline but no attempt quota", {
  observed <- new.env(parent = emptyenv())
  observed$attempts <- 0L
  observed$sleeps <- numeric()
  clock_values <- c(0, 0, 0.01, 0.03, 0.06)
  observed$clock_index <- 0L
  clock <- function() {
    observed$clock_index <- observed$clock_index + 1L
    clock_values[[min(observed$clock_index, length(clock_values))]]
  }
  condition <- tryCatch(
    .dsvert_retry_idempotent(
      attempt = function() {
        observed$attempts <- observed$attempts + 1L
        NULL
      },
      classify = function(response) list(state = "missing"),
      operation = "test unavailable frame", timeout_seconds = 0.05,
      .clock = clock,
      .sleep = function(seconds) {
        observed$sleeps <- c(observed$sleeps, seconds)
      },
      .jitter = function() 1),
    error = identity)
  expect_s3_class(condition, "retry_deadline_exceeded")
  expect_s3_class(condition, "transport_unavailable")
  expect_identical(condition$code, "retry_deadline_exceeded")
  expect_gt(observed$attempts, 2L)
  expect_length(observed$sleeps, observed$attempts - 1L)

  observed$attempts <- 0L
  result <- .dsvert_retry_idempotent(
    attempt = function() {
      observed$attempts <- observed$attempts + 1L
      if (observed$attempts < 6L) NULL else TRUE
    },
    classify = function(response) list(
      state = if (isTRUE(response)) "ack" else "missing"),
    operation = "test recovered frame", timeout_seconds = 10,
    .clock = local({ now <- 0; function() { now <<- now + 0.01; now } }),
    .sleep = function(seconds) invisible(NULL),
    .jitter = function() 1)
  expect_identical(result$state, "ack")
  expect_identical(observed$attempts, 6L)
})

test_that("typed blob phases propagate an unrecognized peer without retry", {
  blob <- "AbCdEf0123_-"
  transfer <- .typed_client_transfer(blob)
  conn <- list(peer = structure(list(), class = "test"))
  producer <- list(self = structure(list(), class = "test"))
  sid <- "12345678-1234-4234-8234-123456789abc"

  calls <- 0L
  rejected_store <- function(conns, expr, error, ...) {
    calls <<- calls + 1L
    error("peer", .typed_peer_error_token("peer"))
    list(peer = NULL)
  }
  condition <- tryCatch(
    .dsvert_store_typed_blob(
      blob, transfer, conn, sid, producer_conn = producer,
      .aggregate = rejected_store),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$peer_name, "peer")
  expect_identical(calls, 1L)

  calls <- 0L
  rejected_receipt <- function(conns, expr, error, ...) {
    calls <<- calls + 1L
    if (identical(as.character(expr[[1L]]), "mpcTypedBlobStoreDS")) {
      return(list(peer = .typed_client_ack(expr, nchar(blob))))
    }
    error("self", .typed_peer_error_token("self"))
    list(self = NULL)
  }
  condition <- tryCatch(
    .dsvert_store_typed_blob(
      blob, transfer, conn, sid, producer_conn = producer,
      .aggregate = rejected_receipt),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$peer_name, "self")
  expect_identical(calls, 2L)

  source_transfer <- .typed_client_transfer(
    blob, recipient = "recipient", sender = "producer")
  calls <- 0L
  rejected_read <- function(conns, expr, error, ...) {
    calls <<- calls + 1L
    error("producer", .typed_peer_error_token("producer"))
    list(producer = NULL)
  }
  condition <- tryCatch(
    .dsvert_read_typed_blob_source_frame(
      source_transfer,
      producer_conn = list(producer = structure(list(), class = "test")),
      offset = 0, max_chars = nchar(blob), session_id = sid,
      .aggregate = rejected_read),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$peer_name, "producer")
  expect_identical(calls, 1L)
  expect_match(conditionMessage(condition), "Server administrator")
})

test_that("typed blob propagates a top-level typed peer rejection", {
  rejection <- .dsvert_client_peer_not_recognized_condition(
    "peer", strrep("2", 64L), strrep("1", 64L))
  calls <- 0L
  condition <- tryCatch(
    .dsvert_typed_blob_dsi_attempt(
      list(peer = structure(list(), class = "test")),
      call(name = "mpcTypedBlobStoreDS"),
      .aggregate = function(...) {
        calls <<- calls + 1L
        stop(rejection)
      }),
    dsvert_peer_not_recognized = identity)
  expect_identical(condition, rejection)
  expect_identical(calls, 1L)
})

test_that("typed blob immediately propagates a poisoned DSI session", {
  poisoned <- .dsvert_dsi_poisoned_session_condition("peer")
  calls <- 0L
  condition <- tryCatch(
    .dsvert_typed_blob_dsi_attempt(
      list(peer = structure(list(), class = "test")),
      call(name = "mpcTypedBlobStoreDS"),
      .aggregate = function(...) {
        calls <<- calls + 1L
        stop(poisoned)
      }),
    dsvert_dsi_poisoned_session = identity)
  expect_identical(condition, poisoned)
  expect_identical(calls, 1L)
})

test_that("typed blob preserves a callback peer rejection over a generic error", {
  calls <- 0L
  condition <- tryCatch(
    .dsvert_typed_blob_dsi_attempt(
      list(peer = structure(list(), class = "test")),
      call(name = "mpcTypedBlobStoreDS"),
      .aggregate = function(conns, expr, error, ...) {
        calls <<- calls + 1L
        error("peer", .typed_peer_error_token("peer"))
        stop("generic connector failure", call. = FALSE)
      }),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$peer_name, "peer")
  expect_identical(calls, 1L)
})

test_that("typed blob treats oversize and contract callbacks as terminal", {
  conn <- list(peer = structure(list(), class = "test"))
  expression <- call(name = "mpcTypedBlobStoreDS")
  tag <- paste0(
    "[dsvert_resource_oversize:v1] resource_oversize: ",
    "typed payload cannot fit fixed policy")

  calls <- 0L
  callback_oversize <- function(conns, expr, error, ...) {
    calls <<- calls + 1L
    error("peer", tag)
    list(peer = NULL)
  }
  condition <- tryCatch(.dsvert_typed_blob_dsi_attempt(
    conn, expression, .aggregate = callback_oversize), error = identity)
  expect_s3_class(condition, "dsvert_resource_oversize")
  expect_identical(condition$code, "resource_oversize")
  expect_false(condition$retryable)
  expect_identical(calls, 1L)

  calls <- 0L
  top_oversize <- function(...) {
    calls <<- calls + 1L
    stop(tag, call. = FALSE)
  }
  condition <- tryCatch(.dsvert_typed_blob_dsi_attempt(
    conn, expression, .aggregate = top_oversize), error = identity)
  expect_s3_class(condition, "dsvert_resource_oversize")
  expect_false(condition$retryable)
  expect_identical(calls, 1L)

  calls <- 0L
  signature_failure <- function(conns, expr, error, ...) {
    calls <<- calls + 1L
    error("peer", "signature or pinned contract mismatch")
    list(peer = NULL)
  }
  condition <- tryCatch(.dsvert_typed_blob_dsi_attempt(
    conn, expression, .aggregate = signature_failure), error = identity)
  expect_s3_class(condition, "dsvert_remote_contract_failure")
  expect_false(condition$retryable)
  expect_identical(calls, 1L)
})

test_that("pipelined typed transport never retries terminal oversize", {
  blob <- "AbCdEf0123_-"
  transfer <- .typed_client_transfer(
    blob, recipient = "recipient", sender = "producer")
  conns <- list(
    producer = structure(list(), class = "test"),
    recipient = structure(list(), class = "test"))
  tag <- paste0(
    "[dsvert_resource_oversize:v1] resource_oversize: ",
    "pipeline cannot fit fixed policy")
  calls <- 0L
  aggregate <- function(conns, expr, error, ...) {
    calls <<- calls + 1L
    error("recipient", tag)
    list(producer = NULL, recipient = NULL)
  }
  condition <- tryCatch(.dsvert_typed_blob_pipeline_step(
    transfer, conns["producer"], conns["recipient"],
    call(name = "mpcTypedBlobStoreDS"), next_offset = 1,
    max_chars = 4, session_id =
      "12345678-1234-4234-8234-123456789abc",
    .aggregate = aggregate), error = identity)
  expect_s3_class(condition, "dsvert_resource_oversize")
  expect_false(condition$retryable)
  expect_identical(calls, 1L)
})

test_that("typed blob client resumes recipient progress without retransmission", {
  blob <- paste(rep(c("A", "b", "1", "_"), 8L), collapse = "")
  transfer <- .typed_client_transfer(blob)
  frame_offsets <- numeric()
  aggregate <- function(conns, expr, async, error, errors.print) {
    if (identical(as.character(expr[[1L]]), "mpcTypedBlobReceiptDS")) {
      return(list(self = .typed_client_receipt_ack()))
    }
    frame_offsets <<- c(frame_offsets, as.numeric(expr$offset))
    if (length(frame_offsets) == 1L) {
      return(list(peer = list(
        version = "dsvert-typed-blob-v1",
        transfer_id = transfer$transfer_id,
        committed_chars = 14, total_chars = transfer$payload_chars,
        sealed = FALSE, receipt = NULL)))
    }
    list(peer = .typed_client_ack(expr, nchar(blob)))
  }
  withr::local_options(list(dsvert.chunk_size = 7L))
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  expect_invisible(.dsvert_store_typed_blob(
    blob, transfer, list(peer = structure(list(), class = "test")),
    "12345678-1234-4234-8234-123456789abc",
    producer_conn = list(self = structure(list(), class = "test")),
    .aggregate = aggregate))
  expect_identical(frame_offsets, c(0, 14, 21, 28))
})

test_that("typed blob client rejects misrouting, mutation and malformed ACKs", {
  blob <- "AbCdEf0123_-"
  transfer <- .typed_client_transfer(blob)
  conn <- list(peer = structure(list(), class = "test"))
  producer <- list(self = structure(list(), class = "test"))
  sid <- "12345678-1234-4234-8234-123456789abc"
  expect_error(.dsvert_store_typed_blob(
    blob, .typed_client_transfer(blob, "other"), conn, sid,
    producer_conn = producer,
    .aggregate = function(...) NULL), "misrouted")
  changed <- transfer
  changed$payload_sha256 <- strrep("0", 64L)
  expect_error(.dsvert_store_typed_blob(
    blob, changed, conn, sid, producer_conn = producer,
    .aggregate = function(...) NULL),
    "does not match")
  extra <- transfer
  extra$key <- "attacker_key"
  expect_error(.dsvert_store_typed_blob(
    blob, extra, conn, sid, producer_conn = producer,
    .aggregate = function(...) NULL),
    "invalid typed-blob transfer contract")

  calls <- 0L
  malformed <- function(conns, expr, ...) {
    calls <<- calls + 1L
    list(peer = list(committed_chars = nchar(blob)))
  }
  expect_error(.dsvert_store_typed_blob(
    blob, transfer, conn, sid, producer_conn = producer,
    .aggregate = malformed),
    "malformed acknowledgement")
  expect_identical(calls, 1L)

  calls <- 0L
  rejected <- function(conns, expr, ...) {
    calls <<- calls + 1L
    list(peer = list(
      version = "dsvert-typed-blob-rejection-v1",
      operation = "store", rejected = TRUE))
  }
  expect_error(.dsvert_store_typed_blob(
    blob, transfer, conn, sid, producer_conn = producer,
    .aggregate = rejected), "rejected the immutable frame")
  expect_identical(calls, 1L)
})
