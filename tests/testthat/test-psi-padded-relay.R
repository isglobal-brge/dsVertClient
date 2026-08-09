.psi_client_relay_frame <- function(descriptor, sequence, offset, bytes,
                                    final = FALSE) {
  list(
    version = "dsvert-relay-v1",
    session_id = "52345678-1234-4234-9234-123456789abc",
    operation_id = descriptor$operation_id,
    sender_peer_id = descriptor$sender_peer_id,
    recipient_peer_id = descriptor$recipient_peer_id,
    capability_id = descriptor$capability_id,
    sequence = sequence, offset = offset, chunk_bytes = bytes,
    total_bytes = descriptor$total_bytes, final = final,
    payload_hash = descriptor$payload_hash,
    chunk_hash = strrep(if (sequence == 0) "6" else "7", 64L),
    payload = if (sequence == 0) "QUFBQUFBQUFBQQ" else "QkJCQkJCQkJCQg",
    signature = "signature")
}

test_that("padded PSI client journal resumes frozen frames after a lost delivery ACK", {
  root <- tempfile("dsvert-client-psi-relay-")
  dir.create(root, mode = "0700")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  old <- options(dsvert.client.state_dir = root)
  on.exit(options(old), add = TRUE)
  peers <- c("alpha", "beta")
  peer_ids <- paste0("dsv1_", c(strrep("a", 64L), strrep("b", 64L)))
  contract <- list(
    peer_names = peers, peer_ids = peer_ids,
    contract_hash = strrep("c", 64L), relay_frame_bytes = 16L * 1024L)
  descriptor <- list(
    version = "dsvert-psi-padded-relay-v4",
    protocol = "dsvert-pinned-padded-psi-v4",
    contract_hash = contract$contract_hash,
    operation_id = paste0("op_", strrep("d", 32L)),
    capability_id = "psi.padded.v4",
    sender = "alpha", recipient = "beta",
    sender_peer_id = peer_ids[[1L]], recipient_peer_id = peer_ids[[2L]],
    total_bytes = 20, payload_hash = strrep("e", 64L),
    frame_bytes = 16L * 1024L)
  frames <- list(
    .psi_client_relay_frame(descriptor, 0, 0, 10, FALSE),
    .psi_client_relay_frame(descriptor, 1, 10, 10, TRUE))
  transport <- list(transport = "relay", envelope = NULL,
                    relay = descriptor)
  conns <- list(alpha = list(), beta = list())
  recipient_attempts <- 0L
  source_attempts <- 0L
  fail_recipient <- TRUE
  aggregate <- function(conns, expr, error, async, errors.print, ...) {
    peer <- names(conns)[[1L]]
    route <- jsonlite::fromJSON(
      .dsvert_dsi_text_decode(as.character(expr[["request_json"]])),
      simplifyVector = FALSE)
    if (identical(peer, "alpha")) {
      source_attempts <<- source_attempts + 1L
      cursor <- as.numeric(route[[peer_ids[[1L]]]]$outbox_cursor)
      value <- if (cursor == 0) list(
        peer_id = peer_ids[[1L]], accepted = list(),
        outbox_cursor = 20, outbox_eof = 20, outbound = frames) else list(
        peer_id = peer_ids[[1L]], accepted = list(),
        outbox_cursor = 20, outbox_eof = 20, outbound = list())
    } else {
      recipient_attempts <<- recipient_attempts + 1L
      if (isTRUE(fail_recipient)) stop("simulated lost DSI response")
      accepted <- list(
        # The first response was lost after the server committed both frames.
        # On the byte-frozen retry every duplicate carries the recipient's
        # current cumulative absolute ACK.
        list(status = "duplicate", operation_id = descriptor$operation_id,
             offset = 0, ack_offset = 20, terminal = TRUE,
             receipt = list(
               version = "dsvert-relay-receipt-v1",
               session_id = "52345678-1234-4234-9234-123456789abc",
               operation_id = descriptor$operation_id,
               sender_peer_id = descriptor$sender_peer_id,
               recipient_peer_id = descriptor$recipient_peer_id,
               capability_id = descriptor$capability_id,
               total_bytes = 20, payload_hash = descriptor$payload_hash,
               ack_offset = 20, terminal = TRUE,
               signature = strrep("S", 86L))),
        list(status = "duplicate", operation_id = descriptor$operation_id,
             offset = 10, ack_offset = 20, terminal = TRUE,
             receipt = list(
               version = "dsvert-relay-receipt-v1",
               session_id = "52345678-1234-4234-9234-123456789abc",
               operation_id = descriptor$operation_id,
               sender_peer_id = descriptor$sender_peer_id,
               recipient_peer_id = descriptor$recipient_peer_id,
               capability_id = descriptor$capability_id,
               total_bytes = 20, payload_hash = descriptor$payload_hash,
               ack_offset = 20, terminal = TRUE,
               signature = strrep("S", 86L))))
      value <- list(
        peer_id = peer_ids[[2L]], accepted = accepted,
        outbox_cursor = 0, outbox_eof = 0, outbound = list())
    }
    stats::setNames(list(value), peer)
  }
  session_id <- "52345678-1234-4234-9234-123456789abc"
  controller <- .dsvert_psi_relay_controller(session_id, contract)
  # Bound this first invocation to one transport attempt so the test exercises
  # durable cross-process resumption rather than the independent deadline-based
  # retry loop (which intentionally has no request-count limit).
  expect_error(testthat::with_mocked_bindings(
    .dsvert_psi_relay_deliver(
      transport, "alpha", "beta", contract, conns, session_id,
      controller, aggregate),
    .dsvert_retry_idempotent = function(attempt, classify, ...) {
      classify(attempt())
    }, .package = "dsVertClient"), "transport failed")
  expect_false(is.null(controller$value$pending))
  expect_identical(controller$value$pending$frames, frames)

  # Re-open from disk as a fresh client process would. The exact pending bytes
  # are delivered before the source cursor is advanced.
  fail_recipient <- FALSE
  resumed <- .dsvert_psi_relay_controller(session_id, contract)
  expect_equal(resumed$value$pending$frames, frames)
  result <- .dsvert_psi_relay_deliver(
    transport, "alpha", "beta", contract, conns, session_id,
    resumed, aggregate)
  expect_identical(result$envelope, "")
  expect_true(nzchar(result$relay_descriptor_b64url))
  expect_null(resumed$value$pending)
  expect_identical(as.numeric(resumed$value$cursors[[peer_ids[[1L]]]]), 20)
  expect_identical(source_attempts, 2L)
  expect_identical(recipient_attempts, 2L)
})
