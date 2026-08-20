.DSVERT_CLIENT_PSI_RELAY_VERSION <- "dsvert-psi-padded-relay-v4"
.DSVERT_CLIENT_PSI_RELAY_CAPABILITY <- "psi.padded.v4"
.DSVERT_CLIENT_PSI_RELAY_JOURNAL <- "dsvert-client-psi-relay-journal-v1"

.dsvert_psi_relay_path_is_link <- function(path) {
  target <- Sys.readlink(path)
  length(target) == 1L && !is.na(target) && nzchar(target)
}

.dsvert_psi_relay_state_root <- function() {
  configured <- getOption("dsvert.client.state_dir")
  root <- if (!is.null(configured)) configured else file.path(
    Sys.getenv("HOME", unset = ""), ".dsvertclient")
  if (!is.character(root) || length(root) != 1L || is.na(root) ||
      !nzchar(root)) stop("The dsVertClient persistent state directory is unavailable.",
                          call. = FALSE)
  root <- path.expand(root)
  if (!grepl("^/", root)) stop(
    "The dsVertClient persistent state directory must be absolute.",
    call. = FALSE)
  directory <- file.path(root, "psi-padded-relay-v1")
  if (!dir.exists(directory) && !dir.create(
      directory, recursive = TRUE, mode = "0700", showWarnings = FALSE)) {
    stop("Could not create the private PSI relay journal directory.",
         call. = FALSE)
  }
  if (.dsvert_psi_relay_path_is_link(directory)) stop(
    "The PSI relay journal directory must not be a symbolic link.",
    call. = FALSE)
  Sys.chmod(directory, mode = "0700")
  directory
}

.dsvert_psi_relay_journal_path <- function(session_id) {
  if (!is.character(session_id) || length(session_id) != 1L ||
      is.na(session_id) || !grepl(
        "^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$",
        session_id)) stop("Invalid PSI relay journal session.", call. = FALSE)
  file.path(.dsvert_psi_relay_state_root(), paste0(session_id, ".json"))
}

.dsvert_psi_relay_journal_hash <- function(value) {
  digest::digest(as.character(jsonlite::toJSON(
    value, auto_unbox = TRUE, null = "null", digits = NA)),
    algo = "sha256", serialize = FALSE)
}

.dsvert_psi_relay_journal_write <- function(controller) {
  if (!is.environment(controller) || !is.list(controller$value)) stop(
    "Invalid PSI relay journal controller.", call. = FALSE)
  value <- controller$value
  unsigned <- value[setdiff(names(value), "sha256")]
  value$sha256 <- .dsvert_psi_relay_journal_hash(unsigned)
  path <- controller$path
  temporary <- tempfile(".psi-relay-journal-", tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary, force = TRUE),
          add = TRUE)
  connection <- file(temporary, open = "wb")
  on.exit(try(if (isOpen(connection)) close(connection), silent = TRUE),
          add = TRUE)
  writeBin(charToRaw(as.character(jsonlite::toJSON(
    value, auto_unbox = TRUE, null = "null", digits = NA))), connection)
  flush(connection)
  close(connection)
  Sys.chmod(temporary, mode = "0600")
  if (!file.rename(temporary, path)) stop(
    "Could not atomically commit the PSI relay journal.", call. = FALSE)
  Sys.chmod(path, mode = "0600")
  controller$value <- value
  invisible(TRUE)
}

.dsvert_psi_relay_controller <- function(session_id, contract) {
  path <- .dsvert_psi_relay_journal_path(session_id)
  peer_ids <- stats::setNames(as.list(rep(0, length(contract$peer_names))),
                              contract$peer_ids)
  initial <- list(
    version = .DSVERT_CLIENT_PSI_RELAY_JOURNAL,
    session_id = session_id, contract_hash = contract$contract_hash,
    cursors = peer_ids, pending = NULL, sha256 = "")
  value <- if (file.exists(path)) {
    info <- file.info(path)
    if (nrow(info) != 1L || is.na(info$size) || info$size < 1L ||
        info$size > 4L * 1024L^2 ||
        .dsvert_psi_relay_path_is_link(path)) stop(
      "The PSI relay journal has an invalid representation.", call. = FALSE)
    parsed <- tryCatch(jsonlite::fromJSON(
      rawToChar(readBin(path, "raw", n = info$size + 1L)),
      simplifyVector = FALSE), error = function(e) NULL)
    if (!is.list(parsed) || !identical(
        parsed$version, .DSVERT_CLIENT_PSI_RELAY_JOURNAL) ||
        !identical(parsed$session_id, session_id) ||
        !identical(parsed$contract_hash, contract$contract_hash) ||
        !is.list(parsed$cursors) ||
        !setequal(names(parsed$cursors), contract$peer_ids) ||
        !is.character(parsed$sha256) ||
        !identical(parsed$sha256, .dsvert_psi_relay_journal_hash(
          parsed[setdiff(names(parsed), "sha256")]))) stop(
      "The PSI relay journal cannot be authenticated structurally.",
      call. = FALSE)
    parsed
  } else initial
  controller <- new.env(parent = emptyenv())
  controller$path <- path
  controller$value <- value
  if (!file.exists(path)) .dsvert_psi_relay_journal_write(controller)
  controller
}

.dsvert_psi_relay_journal_delete <- function(controller) {
  if (!is.environment(controller) || !is.character(controller$path) ||
      length(controller$path) != 1L) return(invisible(FALSE))
  if (file.exists(controller$path) &&
      !.dsvert_psi_relay_path_is_link(controller$path)) {
    unlink(controller$path, force = TRUE)
  }
  invisible(TRUE)
}

.dsvert_psi_relay_request_json <- function(routes) {
  as.character(jsonlite::toJSON(
    routes, auto_unbox = TRUE, null = "null", digits = NA))
}

.dsvert_psi_relay_exchange <- function(
    datasources, peer, route, session_id, operation_id, .aggregate,
    terminal_receipt_b64url = "") {
  expression <- call(
    name = "psiPaddedRelayExchangeDS",
    request_json = .dsvert_psi_relay_request_json(route),
    session_id = session_id, outbound_operation_id = operation_id)
  if (nzchar(terminal_receipt_b64url)) {
    expression$terminal_receipt_b64url <- terminal_receipt_b64url
  }
  result <- .dsvert_aggregate_strict(
    datasources[peer], expression,
    operation = "pinned padded PSI framed relay",
    .aggregate = .aggregate)[[1L]]
  required <- c("peer_id", "accepted", "outbox_cursor", "outbox_eof",
                "outbound")
  if (!is.list(result) || !identical(names(result), required) ||
      !is.list(result$accepted) || !is.list(result$outbound)) stop(
    "A padded PSI relay peer returned an invalid exchange response.",
    call. = FALSE)
  result
}

.dsvert_psi_relay_validate_descriptor <- function(
    descriptor, contract, sender, recipient) {
  required <- c(
    "version", "protocol", "contract_hash", "operation_id",
    "capability_id", "sender", "recipient", "sender_peer_id",
    "recipient_peer_id", "total_bytes", "payload_hash", "frame_bytes")
  fail <- function() stop("Invalid padded PSI relay descriptor returned by a peer.",
                          call. = FALSE)
  if (!is.list(descriptor) || !identical(names(descriptor), required) ||
      !identical(descriptor$version, .DSVERT_CLIENT_PSI_RELAY_VERSION) ||
      !identical(descriptor$protocol,
                 .DSVERT_CLIENT_PSI_PADDED_PROTOCOL) ||
      !identical(descriptor$contract_hash, contract$contract_hash) ||
      !identical(descriptor$capability_id,
                 .DSVERT_CLIENT_PSI_RELAY_CAPABILITY) ||
      !identical(descriptor$sender, sender) ||
      !identical(descriptor$recipient, recipient) ||
      !identical(descriptor$sender_peer_id,
                 contract$peer_ids[[match(sender, contract$peer_names)]]) ||
      !identical(descriptor$recipient_peer_id,
                 contract$peer_ids[[match(recipient, contract$peer_names)]]) ||
      !identical(as.numeric(descriptor$frame_bytes),
                 as.numeric(contract$relay_frame_bytes)) ||
      !is.character(descriptor$operation_id) ||
      !grepl("^op_[0-9a-f]{32}$", descriptor$operation_id) ||
      !is.character(descriptor$payload_hash) ||
      !grepl("^[0-9a-f]{64}$", descriptor$payload_hash)) fail()
  descriptor$total_bytes <- suppressWarnings(as.numeric(descriptor$total_bytes))
  if (length(descriptor$total_bytes) != 1L ||
      is.na(descriptor$total_bytes) || !is.finite(descriptor$total_bytes) ||
      descriptor$total_bytes < 1 || descriptor$total_bytes !=
        floor(descriptor$total_bytes) || descriptor$total_bytes >
        64L * 1024L^2) fail()
  descriptor
}

.dsvert_psi_relay_frame_end <- function(frame, descriptor) {
  required <- c(
    "version", "session_id", "operation_id", "sender_peer_id",
    "recipient_peer_id", "capability_id", "sequence", "offset",
    "chunk_bytes", "total_bytes", "final", "payload_hash", "chunk_hash",
    "payload", "signature")
  if (!is.list(frame) || !identical(names(frame), required) ||
      !identical(frame$operation_id, descriptor$operation_id) ||
      !identical(frame$sender_peer_id, descriptor$sender_peer_id) ||
      !identical(frame$recipient_peer_id, descriptor$recipient_peer_id) ||
      !identical(frame$capability_id, descriptor$capability_id) ||
      !identical(frame$payload_hash, descriptor$payload_hash) ||
      !identical(as.numeric(frame$total_bytes), descriptor$total_bytes)) stop(
    "A padded PSI relay peer returned a misbound frame.", call. = FALSE)
  offset <- suppressWarnings(as.numeric(frame$offset))
  chunk <- suppressWarnings(as.numeric(frame$chunk_bytes))
  if (length(offset) != 1L || is.na(offset) || !is.finite(offset) ||
      offset < 0 || offset != floor(offset) || length(chunk) != 1L ||
      is.na(chunk) || !is.finite(chunk) || chunk < 1 || chunk != floor(chunk) ||
      chunk > descriptor$frame_bytes || offset + chunk >
        descriptor$total_bytes) stop("Invalid padded PSI relay frame bounds.",
                                     call. = FALSE)
  offset + chunk
}

.dsvert_psi_relay_deliver <- function(
    transport, sender, recipient, contract, datasources, session_id,
    controller, .aggregate) {
  if (!is.list(transport) ||
      !transport$transport %in% c("inline", "local", "relay")) stop(
    "A padded PSI producer returned an invalid transport selector.",
    call. = FALSE)
  if (identical(transport$transport, "inline")) {
    if (!is.character(transport$envelope) ||
        length(transport$envelope) != 1L || is.na(transport$envelope) ||
        !nzchar(transport$envelope) || !is.null(transport$relay)) stop(
      "A padded PSI producer returned an invalid inline envelope.",
      call. = FALSE)
    return(list(envelope = transport$envelope,
                relay_descriptor_b64url = ""))
  }
  if (!is.null(transport$envelope)) stop(
    "A padded PSI relayed transport exposed an inline envelope.",
    call. = FALSE)
  descriptor <- .dsvert_psi_relay_validate_descriptor(
    transport$relay, contract, sender, recipient)
  if (identical(transport$transport, "local")) {
    if (!identical(sender, recipient)) stop(
      "A padded PSI local transport crossed peer boundaries.",
      call. = FALSE)
    return(list(
      envelope = "",
      relay_descriptor_b64url =
        .dsvert_psi_padded_json_b64url(descriptor)))
  }
  if (identical(sender, recipient)) stop(
    "A padded PSI self-envelope was sent through the external relay.",
    call. = FALSE)
  source_id <- descriptor$sender_peer_id
  recipient_id <- descriptor$recipient_peer_id
  cursor <- suppressWarnings(as.numeric(controller$value$cursors[[source_id]]))
  if (length(cursor) != 1L || is.na(cursor) || !is.finite(cursor) ||
      cursor < 0 || cursor != floor(cursor)) stop(
    "The PSI relay journal contains an invalid source cursor.",
    call. = FALSE)

  pending <- controller$value$pending
  if (!is.null(pending)) {
    same <- identical(pending$operation_id, descriptor$operation_id) &&
      identical(pending$sender, sender) &&
      identical(pending$recipient, recipient) &&
      identical(pending$payload_hash, descriptor$payload_hash)
    if (!isTRUE(same)) stop(
      "The PSI relay journal contains a different unfinished transfer.",
      call. = FALSE)
  }
  repeat {
    if (is.null(pending)) {
      source_route <- stats::setNames(list(list(
        outbox_cursor = cursor, deliveries = list())), source_id)
      source <- .dsvert_psi_relay_exchange(
        datasources, sender, source_route, session_id,
        descriptor$operation_id, .aggregate)
      if (!identical(source$peer_id, source_id)) stop(
        "A padded PSI frame was returned by the wrong peer.", call. = FALSE)
      frames <- source$outbound
      if (!length(frames)) {
        if (!identical(as.numeric(source$outbox_cursor),
                       as.numeric(source$outbox_eof))) stop(
          "A padded PSI relay returned an unexplained byte gap.",
          call. = FALSE)
        controller$value$cursors[[source_id]] <-
          as.numeric(source$outbox_cursor)
        controller$value$pending <- NULL
        .dsvert_psi_relay_journal_write(controller)
        break
      }
      ends <- vapply(frames, .dsvert_psi_relay_frame_end, numeric(1L),
                     descriptor = descriptor)
      offsets <- vapply(frames, function(frame) {
        as.numeric(frame$offset)
      }, numeric(1L))
      chunks <- vapply(frames, function(frame) {
        as.numeric(frame$chunk_bytes)
      }, numeric(1L))
      if (any(diff(offsets) < 0) ||
          (length(offsets) > 1L &&
           any(offsets[-1L] != ends[-length(ends)])) ||
          !identical(as.numeric(source$outbox_cursor),
                     cursor + sum(chunks))) stop(
        "A padded PSI relay returned reordered or discontinuous frames.",
        call. = FALSE)
      pending <- list(
        operation_id = descriptor$operation_id,
        sender = sender, recipient = recipient,
        payload_hash = descriptor$payload_hash,
        next_cursor = as.numeric(source$outbox_cursor), frames = frames)
      controller$value$pending <- pending
      .dsvert_psi_relay_journal_write(controller)
    }
    recipient_cursor <- suppressWarnings(as.numeric(
      controller$value$cursors[[recipient_id]]))
    delivery_route <- stats::setNames(list(list(
      outbox_cursor = recipient_cursor,
      deliveries = pending$frames)), recipient_id)
    accepted <- .dsvert_psi_relay_exchange(
      datasources, recipient, delivery_route, session_id,
      descriptor$operation_id, .aggregate)
    if (!identical(accepted$peer_id, recipient_id) ||
        length(accepted$accepted) != length(pending$frames)) stop(
      "A padded PSI relay recipient returned an invalid acknowledgment.",
      call. = FALSE)
    expected_acks <- vapply(
      pending$frames, .dsvert_psi_relay_frame_end, numeric(1L),
      descriptor = descriptor)
    actual_acks <- vapply(accepted$accepted, function(value) {
      suppressWarnings(as.numeric(value$ack_offset))
    }, numeric(1L))
    # ACKs are cumulative absolute offsets.  A byte-frozen retry after the
    # recipient committed the whole batch therefore returns the already
    # advanced offset for every duplicated frame, rather than replaying the
    # historical per-frame offsets.  Both forms are valid, but no ACK may
    # move backwards, precede the frame it authenticates, or exceed the
    # authenticated envelope boundary.
    if (any(!is.finite(actual_acks)) || any(actual_acks < expected_acks) ||
        any(actual_acks > descriptor$total_bytes) ||
        any(diff(actual_acks) < 0)) stop(
      "A padded PSI relay recipient returned a conflicting absolute ACK.",
      call. = FALSE)
    finals <- which(vapply(pending$frames, function(frame) {
      isTRUE(frame$final)
    }, logical(1L)))
    terminal_receipt <- NULL
    if (length(finals)) {
      if (length(finals) != 1L) stop(
        "A padded PSI relay returned multiple terminal frames.",
        call. = FALSE)
      terminal <- accepted$accepted[[finals]]
      receipt <- terminal$receipt
      receipt_fields <- c(
        "version", "session_id", "operation_id", "sender_peer_id",
        "recipient_peer_id", "capability_id", "total_bytes",
        "payload_hash", "ack_offset", "terminal", "signature")
      if (!isTRUE(terminal$terminal) || !is.list(receipt) ||
          !identical(names(receipt), receipt_fields) ||
          !identical(receipt$version, "dsvert-relay-receipt-v1") ||
          !identical(receipt$session_id, session_id) ||
          !identical(receipt$operation_id, descriptor$operation_id) ||
          !identical(receipt$sender_peer_id, descriptor$sender_peer_id) ||
          !identical(receipt$recipient_peer_id,
                     descriptor$recipient_peer_id) ||
          !identical(receipt$capability_id, descriptor$capability_id) ||
          !identical(receipt$payload_hash, descriptor$payload_hash) ||
          !identical(as.numeric(receipt$total_bytes),
                     descriptor$total_bytes) ||
          !identical(as.numeric(receipt$ack_offset),
                     descriptor$total_bytes) ||
          !isTRUE(receipt$terminal) || !is.character(receipt$signature) ||
          length(receipt$signature) != 1L || is.na(receipt$signature) ||
          !grepl("^[A-Za-z0-9_-]{86}$", receipt$signature)) stop(
        "A padded PSI relay recipient omitted the signed terminal receipt.",
        call. = FALSE)
      terminal_receipt <- receipt
    }
    cursor <- pending$next_cursor
    if (!is.null(terminal_receipt)) {
      # The source verifies the receiving peer's pinned signature and final
      # payload hash before the terminal cursor is allowed to compact its
      # spool.  The analyst/relay can delay this call, but cannot forge it.
      source_route <- stats::setNames(list(list(
        outbox_cursor = cursor, deliveries = list())), source_id)
      finalized <- .dsvert_psi_relay_exchange(
        datasources, sender, source_route, session_id,
        descriptor$operation_id, .aggregate,
        terminal_receipt_b64url =
          .dsvert_psi_padded_json_b64url(terminal_receipt))
      if (!identical(finalized$peer_id, source_id) ||
          length(finalized$outbound) ||
          !identical(as.numeric(finalized$outbox_cursor),
                     as.numeric(cursor))) stop(
        "A padded PSI source did not finalize its signed terminal ACK.",
        call. = FALSE)
    }
    controller$value$cursors[[source_id]] <- cursor
    controller$value$pending <- NULL
    .dsvert_psi_relay_journal_write(controller)
    pending <- NULL
    if (!is.null(terminal_receipt)) break
  }
  list(
    envelope = "",
    relay_descriptor_b64url = .dsvert_psi_padded_json_b64url(descriptor))
}
