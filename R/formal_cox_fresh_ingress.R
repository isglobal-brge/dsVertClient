# Private client ingress for configured fresh formal-Cox sources.

.dsvert_formal_cox_fresh_ingress_sha256 <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[0-9a-f]{64}$", value)
}

.dsvert_formal_cox_fresh_ingress_flag <- function(value) {
  is.logical(value) && length(value) == 1L && !is.na(value) && identical(value, FALSE)
}

.dsvert_formal_cox_fresh_ingress_integer <- function(value, field) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value >= 0 && value <= .Machine$integer.max &&
    value == floor(value)
  if (!isTRUE(valid)) {
    stop(paste0("A configured fresh Cox host returned an invalid ", field, "."),
         call. = FALSE)
  }
  as.integer(value)
}

.dsvert_formal_cox_fresh_ingress_peers <- function(value, count = NULL) {
  valid <- is.character(value) && length(value) >= 2L && !anyNA(value) &&
    !anyDuplicated(value) && all(grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)) &&
    (is.null(count) || length(value) == count)
  if (!isTRUE(valid)) {
    stop("A configured fresh Cox host returned an invalid peer set.",
         call. = FALSE)
  }
  enc2utf8(value)
}

.dsvert_formal_cox_fresh_ingress_payload <- function(value, action, fields) {
  value <- .dsvert_formal_cox_fresh_source_reply(value, action)$payload
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields)
  if (!isTRUE(valid)) {
    stop("A configured fresh Cox host returned an invalid ingress record.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_cox_fresh_ingress_shape <- function(value, peer) {
  fields <- c("version", "analysis_id", "schema_sha256", "source",
              "custodian_peers", "designated_compute_peers", "total_blocks",
              "production_ready")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version, "dsvert-formal-cox-fresh-source-shape-v1") &&
    identical(value$source, peer) &&
    .dsvert_formal_cox_fresh_ingress_sha256(value$schema_sha256) &&
    .dsvert_formal_cox_fresh_ingress_flag(value$production_ready)
  if (!isTRUE(valid)) {
    stop("A configured fresh Cox host returned an invalid source shape.",
         call. = FALSE)
  }
  value$analysis_id <- .dsvert_formal_cox_frontdoor_id(
    value$analysis_id, "source analysis id")
  value$custodian_peers <- .dsvert_formal_cox_fresh_ingress_peers(
    value$custodian_peers)
  value$designated_compute_peers <- .dsvert_formal_cox_fresh_ingress_peers(
    value$designated_compute_peers, count = 2L)
  if (!all(value$designated_compute_peers %in% value$custodian_peers)) {
    stop("A configured fresh Cox host returned an invalid compute peer set.",
         call. = FALSE)
  }
  value$total_blocks <- .dsvert_formal_cox_fresh_ingress_integer(
    value$total_blocks, "block count")
  if (value$total_blocks < 1L) {
    stop("A configured fresh Cox host returned an invalid block count.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_cox_fresh_ingress_delivery <- function(value, recipient) {
  fields <- c("version", "purpose", "receipt", "receipt_sha256",
              "recipient_peer_name", "envelope", "binding")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version, "dsvert-formal-cox-blockwise-source-delivery-v1") &&
    identical(value$purpose,
              "formal-cox-recipient-encrypted-source-delivery-v1") &&
    identical(value$recipient_peer_name, recipient) && is.list(value$receipt) &&
    is.list(value$envelope) && is.list(value$binding) &&
    .dsvert_formal_cox_fresh_ingress_sha256(value$receipt_sha256)
  if (!isTRUE(valid)) {
    stop("A configured fresh Cox host returned an invalid encrypted delivery.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_cox_fresh_ingress_replayed <- function(value) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("A configured fresh Cox host returned an invalid replay flag.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_cox_fresh_ingress_worker <- function(value, peer) {
  fields <- c("version", "peer_name", "plan_sha256", "attempt_id",
              "replayed", "production_ready")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version, "dsvert-formal-cox-blockwise-worker-provision-v1") &&
    identical(value$peer_name, peer) &&
    .dsvert_formal_cox_fresh_ingress_sha256(value$plan_sha256) &&
    .dsvert_formal_cox_fresh_ingress_sha256(value$attempt_id) &&
    is.logical(value$replayed) && length(value$replayed) == 1L &&
    !is.na(value$replayed) &&
    .dsvert_formal_cox_fresh_ingress_flag(value$production_ready)
  if (!isTRUE(valid)) {
    stop("A configured fresh Cox host returned an invalid worker selector.",
         call. = FALSE)
  }
  .dsvert_formal_cox_fresh_ingress_replayed(value$replayed)
  value
}

#' Build private source ingress for one configured fresh formal-Cox analysis
#'
#' The result is intentionally internal: it contains only the two server-side
#' worker selectors and never a Cox fit, a source value or a private envelope.
.dsvert_formal_cox_fresh_ingress <- function(
    conns, selector, .aggregate = DSI::datashield.aggregate) {
  valid_conns <- is.list(conns) && length(conns) >= 2L && !is.null(names(conns)) &&
    !anyNA(names(conns)) && !anyDuplicated(names(conns)) && all(nzchar(names(conns)))
  if (!isTRUE(valid_conns)) {
    stop("Configured fresh Cox ingress requires named custodian connections.",
         call. = FALSE)
  }
  if (!is.list(selector) || !identical(
      names(selector), c("analysis_id", "data_name", "formula_sha256"))) {
    stop("Configured fresh Cox ingress requires one fixed source selector.",
         call. = FALSE)
  }
  selector <- .dsvert_formal_cox_fresh_source_selector(
    selector$analysis_id, selector$data_name, selector$formula_sha256)
  peers <- names(conns)
  call <- function(peer, action, payload) {
    .dsvert_formal_cox_fresh_source_call(
      conns[peer], selector, action, payload, .aggregate = .aggregate)
  }
  shapes <- lapply(peers, function(peer) {
    .dsvert_formal_cox_fresh_ingress_shape(
      .dsvert_formal_cox_fresh_ingress_payload(
        call(peer, "shape", structure(list(), names = character())), "shape",
        c("version", "analysis_id", "schema_sha256", "source",
          "custodian_peers", "designated_compute_peers", "total_blocks",
          "production_ready")), peer)
  })
  names(shapes) <- peers
  reference <- shapes[[1L]]
  normalized <- function(shape) {
    shape$source <- NULL
    shape
  }
  if (!identical(reference$analysis_id, selector$analysis_id) ||
      !identical(reference$custodian_peers, peers) ||
      !all(vapply(shapes, function(shape)
        identical(normalized(shape), normalized(reference)), logical(1L)))) {
    stop("Configured fresh Cox custodians returned incompatible source shapes.",
         call. = FALSE)
  }
  compute_peers <- reference$designated_compute_peers
  tickets <- lapply(compute_peers, function(peer) {
    value <- .dsvert_formal_cox_fresh_ingress_payload(
      call(peer, "ticket", structure(list(), names = character())), "ticket", "ticket")
    if (!is.list(value$ticket)) {
      stop("A configured fresh Cox host returned an invalid recipient ticket.",
           call. = FALSE)
    }
    value$ticket
  })
  recipient_tickets <- unname(tickets)
  deliveries <- stats::setNames(vector("list", length(compute_peers)), compute_peers)
  for (source in peers) for (block in seq.int(0L, reference$total_blocks - 1L)) {
    produced <- .dsvert_formal_cox_fresh_ingress_payload(
      call(source, "produce", list(recipient_tickets = recipient_tickets,
                                    block_index = block)), "produce",
      c("receipt", "receipt_sha256", "replayed"))
    if (!is.list(produced$receipt) ||
        !.dsvert_formal_cox_fresh_ingress_sha256(produced$receipt_sha256)) {
      stop("A configured fresh Cox source did not produce an authenticated block.",
           call. = FALSE)
    }
    .dsvert_formal_cox_fresh_ingress_replayed(produced$replayed)
    for (recipient in compute_peers) {
      delivery <- .dsvert_formal_cox_fresh_ingress_delivery(
        .dsvert_formal_cox_fresh_ingress_payload(
          call(source, "delivery", list(recipient_tickets = recipient_tickets,
                                          block_index = block,
                                          recipient_peer_name = recipient)),
          "delivery", c("version", "purpose", "receipt", "receipt_sha256",
                        "recipient_peer_name", "envelope", "binding")), recipient)
      imported <- .dsvert_formal_cox_fresh_ingress_payload(
        call(recipient, "import", list(recipient_tickets = recipient_tickets,
                                         delivery = delivery)), "import",
        c("version", "purpose", "receipt_sha256", "recipient_peer_name", "replayed"))
      if (!identical(imported$version,
                     "dsvert-formal-cox-blockwise-source-import-receipt-v1") ||
          !identical(imported$purpose,
                     "formal-cox-recipient-encrypted-source-delivery-v1") ||
          !identical(imported$receipt_sha256, delivery$receipt_sha256) ||
          !identical(imported$recipient_peer_name, recipient)) {
        stop("A configured fresh Cox recipient rejected an encrypted block.",
             call. = FALSE)
      }
      .dsvert_formal_cox_fresh_ingress_replayed(imported$replayed)
      deliveries[[recipient]] <- delivery
    }
  }
  workers <- lapply(compute_peers, function(peer) {
    .dsvert_formal_cox_fresh_ingress_worker(
      .dsvert_formal_cox_fresh_ingress_payload(
        call(peer, "provision", list(recipient_tickets = recipient_tickets,
                                      delivery = deliveries[[peer]])), "provision",
        c("version", "peer_name", "plan_sha256", "attempt_id",
          "replayed", "production_ready")), peer)
  })
  names(workers) <- compute_peers
  list(analysis_id = reference$analysis_id, schema_sha256 = reference$schema_sha256,
       total_blocks = reference$total_blocks, compute_peers = compute_peers,
       workers = workers, production_ready = FALSE)
}
