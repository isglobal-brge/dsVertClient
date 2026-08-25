# Private Phase18 coordinator for a configured fresh binomial/Poisson GLM.
#
# It treats every ticket, source receipt, binding and encrypted pair as an
# opaque server record.  The only output is the two host selectors required by
# the private Phase20 coordinator; no source value or fitted result crosses R.

.dsvert_formal_glm_registered_fresh_ingress_sha256 <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[0-9a-f]{64}$", value)
}

.dsvert_formal_glm_registered_fresh_ingress_flag <- function(value) {
  is.logical(value) && length(value) == 1L && !is.na(value) && identical(value, FALSE)
}

.dsvert_formal_glm_registered_fresh_ingress_integer <- function(value, field) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value >= 0 && value <= .Machine$integer.max &&
    value == floor(value)
  if (!isTRUE(valid)) {
    stop(paste0("A registered fresh GLM host returned an invalid ", field, "."),
         call. = FALSE)
  }
  as.integer(value)
}

.dsvert_formal_glm_registered_fresh_ingress_json <- function(value, field) {
  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    nchar(value, type = "bytes") >= 2L && nchar(value, type = "bytes") <= 4L * 1024L^2
  if (!isTRUE(valid)) {
    stop(paste0("A registered fresh GLM host returned an invalid ", field, "."),
         call. = FALSE)
  }
  enc2utf8(value)
}

.dsvert_formal_glm_registered_fresh_ingress_shape <- function(value, peer) {
  fields <- c(
    "version", "artifact_id", "source_contract_sha256", "source",
    "custodian_peers", "designated_compute_peers", "total_blocks",
    "production_ready")
  valid_peers <- function(x, count = NULL) {
    is.character(x) && length(x) >= 2L && !anyNA(x) && !anyDuplicated(x) &&
      (is.null(count) || length(x) == count) &&
      all(grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", x))
  }
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version, "dsvert-formal-glm-registered-fresh-source-shape-v1") &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(value$artifact_id) &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(value$source_contract_sha256) &&
    identical(value$source, peer) && valid_peers(value$custodian_peers) &&
    valid_peers(value$designated_compute_peers, count = 2L) &&
    all(value$designated_compute_peers %in% value$custodian_peers) &&
    .dsvert_formal_glm_registered_fresh_ingress_flag(value$production_ready)
  if (!isTRUE(valid)) {
    stop("A registered fresh GLM host returned an invalid source shape.",
         call. = FALSE)
  }
  value$total_blocks <- .dsvert_formal_glm_registered_fresh_ingress_integer(
    value$total_blocks, "block count")
  if (value$total_blocks < 1L) {
    stop("A registered fresh GLM host returned an invalid block count.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_glm_registered_fresh_ingress_payload <- function(
    value, action, fields) {
  value <- .dsvert_formal_glm_registered_source_reply(value, action)$payload
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields)
  if (!isTRUE(valid)) {
    stop("A registered fresh GLM host returned an invalid ingress record.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_glm_registered_fresh_ingress_replayed <- function(value) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("A registered fresh GLM host returned an invalid replay flag.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_glm_registered_fresh_ingress_chunk <- function(
    value, source, block, offset) {
  fields <- c("version", "purpose", "handle", "artifact_id",
              "source_contract_sha256", "authorization_sha256", "source",
              "block_index", "pair_sha256", "pair_bytes", "offset",
              "chunk_sha256", "chunk_bytes", "complete", "production_ready")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$source, source) &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(value$artifact_id) &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(value$source_contract_sha256) &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(value$authorization_sha256) &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(value$pair_sha256) &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(value$chunk_sha256) &&
    .dsvert_formal_glm_registered_fresh_ingress_flag(value$production_ready) &&
    is.logical(value$complete) && length(value$complete) == 1L && !is.na(value$complete)
  if (!isTRUE(valid)) {
    stop("A registered fresh GLM host returned an invalid encrypted block chunk.",
         call. = FALSE)
  }
  observed_block <- .dsvert_formal_glm_registered_fresh_ingress_integer(
    value$block_index, "block index")
  observed_offset <- .dsvert_formal_glm_registered_fresh_ingress_integer(
    value$offset, "chunk offset")
  bytes <- .dsvert_formal_glm_registered_fresh_ingress_integer(
    value$chunk_bytes, "chunk size")
  pair_bytes <- .dsvert_formal_glm_registered_fresh_ingress_integer(
    value$pair_bytes, "pair size")
  if (!identical(observed_block, block) || !identical(observed_offset, offset) ||
      bytes < 1L || pair_bytes < offset + bytes ||
      identical(value$complete, TRUE) != identical(pair_bytes, offset + bytes)) {
    stop("A registered fresh GLM host returned a non-contiguous encrypted block chunk.",
         call. = FALSE)
  }
  list(bytes = bytes, complete = value$complete)
}

#' Build the durable Phase18 ingress for one configured fresh formal GLM
#'
#' This internal function is intentionally unavailable as an R analytical API.
#' It moves only server-generated tickets, receipts and encrypted chunks across
#' custodians, then returns the two provisioned Phase20 host selectors.
.dsvert_formal_glm_registered_fresh_ingress <- function(
    conns, selector, .aggregate = DSI::datashield.aggregate) {
  valid_conns <- is.list(conns) && length(conns) >= 2L && !is.null(names(conns)) &&
    !anyNA(names(conns)) && !anyDuplicated(names(conns)) && all(nzchar(names(conns)))
  if (!isTRUE(valid_conns)) {
    stop("Registered fresh GLM ingress requires named custodian connections.",
         call. = FALSE)
  }
  selector <- .dsvert_formal_glm_registered_fresh_source_selector(
    selector$analysis_id, selector$data_name, selector$family, selector$formula_sha256)
  peers <- names(conns)
  call <- function(peer, action, payload) {
    .dsvert_formal_glm_registered_fresh_source_call(
      conns[peer], selector, action, payload, .aggregate = .aggregate)
  }
  shapes <- lapply(peers, function(peer) {
    .dsvert_formal_glm_registered_fresh_ingress_shape(
      .dsvert_formal_glm_registered_fresh_ingress_payload(
        call(peer, "shape", structure(list(), names = character())), "shape",
        c("version", "artifact_id", "source_contract_sha256", "source",
          "custodian_peers", "designated_compute_peers", "total_blocks",
          "production_ready")), peer)
  })
  names(shapes) <- peers
  reference_shape <- shapes[[1L]]
  normalized <- function(shape) {
    shape$source <- NULL
    shape
  }
  if (!identical(reference_shape$custodian_peers, peers) ||
      !all(vapply(shapes, function(shape)
        identical(normalized(shape), normalized(reference_shape)), logical(1L)))) {
    stop("Registered fresh GLM custodians returned incompatible source shapes.",
         call. = FALSE)
  }
  compute_peers <- reference_shape$designated_compute_peers
  tickets <- lapply(compute_peers, function(peer) {
    value <- .dsvert_formal_glm_registered_fresh_ingress_payload(
      call(peer, "ticket", structure(list(), names = character())), "ticket",
      c("ticket", "replayed"))
    if (!is.list(value$ticket)) {
      stop("A registered fresh GLM host returned an invalid recipient ticket.",
           call. = FALSE)
    }
    .dsvert_formal_glm_registered_fresh_ingress_replayed(value$replayed)
    value$ticket
  })
  recipient_tickets <- unname(tickets)
  for (peer in peers) {
    value <- .dsvert_formal_glm_registered_fresh_ingress_payload(
      call(peer, "ticket_set", list(recipient_tickets = recipient_tickets)), "ticket_set",
      c("ticket_receipts", "replayed"))
    if (!is.list(value$ticket_receipts) || length(value$ticket_receipts) != 2L) {
      stop("A registered fresh GLM host rejected the recipient ticket set.",
           call. = FALSE)
    }
    .dsvert_formal_glm_registered_fresh_ingress_replayed(value$replayed)
  }
  for (source in peers) {
    for (block in seq.int(0L, reference_shape$total_blocks - 1L)) {
      sealed <- .dsvert_formal_glm_registered_fresh_ingress_payload(
        call(source, "seal_block", list(
          recipient_tickets = recipient_tickets, block_index = block)), "seal_block",
        c("source_receipt", "replayed"))
      if (!is.list(sealed$source_receipt)) {
        stop("A registered fresh GLM source did not seal its block.", call. = FALSE)
      }
      .dsvert_formal_glm_registered_fresh_ingress_replayed(sealed$replayed)
      offset <- 0L
      repeat {
        chunk <- .dsvert_formal_glm_registered_fresh_ingress_payload(
          call(source, "chunk", list(
            recipient_tickets = recipient_tickets, block_index = block, offset = offset)),
          "chunk", c("chunk_receipt", "pair_chunk_base64", "replayed"))
        shape <- .dsvert_formal_glm_registered_fresh_ingress_chunk(
          chunk$chunk_receipt, source, block, offset)
        .dsvert_formal_glm_registered_fresh_ingress_replayed(chunk$replayed)
        for (recipient in compute_peers) {
          delivery <- .dsvert_formal_glm_registered_fresh_ingress_payload(
            call(recipient, "import_chunk", list(
              recipient_tickets = recipient_tickets, chunk_receipt = chunk$chunk_receipt,
              pair_chunk_base64 = chunk$pair_chunk_base64)), "import_chunk",
            c("chunk_delivery", "replayed"))
          if (!is.list(delivery$chunk_delivery)) {
            stop("A registered fresh GLM recipient rejected an encrypted block chunk.",
                 call. = FALSE)
          }
          .dsvert_formal_glm_registered_fresh_ingress_replayed(delivery$replayed)
        }
        if (isTRUE(shape$complete)) break
        offset <- offset + shape$bytes
      }
    }
  }
  local_receipts <- lapply(peers, function(source) {
    value <- .dsvert_formal_glm_registered_fresh_ingress_payload(
      call(source, "local_receipt", list(recipient_tickets = recipient_tickets)),
      "local_receipt", c("local_receipt_json", "replayed"))
    .dsvert_formal_glm_registered_fresh_ingress_replayed(value$replayed)
    .dsvert_formal_glm_registered_fresh_ingress_json(value$local_receipt_json,
                                                      "source receipt")
  })
  for (receipt in local_receipts) for (peer in peers) {
    value <- .dsvert_formal_glm_registered_fresh_ingress_payload(
      call(peer, "receipt_commit", list(local_receipt_json = receipt)), "receipt_commit",
      c("local_receipt_json", "replayed"))
    if (!identical(.dsvert_formal_glm_registered_fresh_ingress_json(
      value$local_receipt_json, "committed source receipt"), receipt)) {
      stop("A registered fresh GLM host changed a source receipt.", call. = FALSE)
    }
    .dsvert_formal_glm_registered_fresh_ingress_replayed(value$replayed)
  }
  receipt_sets <- lapply(peers, function(peer) {
    value <- .dsvert_formal_glm_registered_fresh_ingress_payload(
      call(peer, "receipt_set", structure(list(), names = character())), "receipt_set",
      c("receipt_set_json", "replayed"))
    .dsvert_formal_glm_registered_fresh_ingress_replayed(value$replayed)
    .dsvert_formal_glm_registered_fresh_ingress_json(value$receipt_set_json,
                                                      "receipt set")
  })
  if (!all(vapply(receipt_sets, identical, logical(1L), receipt_sets[[1L]]))) {
    stop("Registered fresh GLM custodians sealed different receipt sets.", call. = FALSE)
  }
  for (peer in peers) {
    value <- .dsvert_formal_glm_registered_fresh_ingress_payload(
      call(peer, "binding", list(recipient_tickets = recipient_tickets)), "binding",
      c("binding_record_json", "replayed"))
    .dsvert_formal_glm_registered_fresh_ingress_replayed(value$replayed)
    .dsvert_formal_glm_registered_fresh_ingress_json(value$binding_record_json,
                                                      "attempt binding")
  }
  hosts <- lapply(compute_peers, function(peer) {
    .dsvert_formal_glm_registered_job_receipt(
      .dsvert_formal_glm_registered_fresh_ingress_payload(
        call(peer, "host_provision", structure(list(), names = character())),
        "host_provision", c("version", "peer", "artifact_id", "receipt_set_sha256",
                            "config_sha256", "replayed", "production_ready")))
  })
  names(hosts) <- compute_peers
  list(artifact_id = reference_shape$artifact_id,
       source_contract_sha256 = reference_shape$source_contract_sha256,
       total_blocks = reference_shape$total_blocks,
       compute_peers = compute_peers, hosts = hosts, production_ready = FALSE)
}

#' Run configured fresh formal-GLM ingress through its durable Phase20 handoff
#'
#' This remains internal until the ordinary frontdoor invokes it for a
#' server-configured source. It completes the two-authority publication but
#' never carries a model result or a private lifecycle record.
.dsvert_formal_glm_registered_fresh_run <- function(
    conns, selector, .aggregate = DSI::datashield.aggregate, max_cycles = NULL) {
  ingress <- .dsvert_formal_glm_registered_fresh_ingress(
    conns, selector, .aggregate = .aggregate)
  fields <- c("artifact_id", "source_contract_sha256", "total_blocks",
              "compute_peers", "hosts", "production_ready")
  valid <- is.list(ingress) && !is.null(names(ingress)) && !anyNA(names(ingress)) &&
    !anyDuplicated(names(ingress)) && setequal(names(ingress), fields) &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(ingress$artifact_id) &&
    .dsvert_formal_glm_registered_fresh_ingress_sha256(ingress$source_contract_sha256) &&
    is.character(ingress$compute_peers) && length(ingress$compute_peers) == 2L &&
    !anyNA(ingress$compute_peers) && !anyDuplicated(ingress$compute_peers) &&
    identical(names(ingress$hosts), ingress$compute_peers) &&
    all(vapply(ingress$hosts, function(host)
      !is.null(.dsvert_formal_glm_registered_job_receipt(host)), logical(1L))) &&
    .dsvert_formal_glm_registered_fresh_ingress_flag(ingress$production_ready)
  if (!isTRUE(valid)) {
    stop("Registered fresh GLM ingress returned an invalid host set.", call. = FALSE)
  }
  blocks <- .dsvert_formal_glm_registered_fresh_ingress_integer(
    ingress$total_blocks, "block count")
  if (blocks < 1L || any(!ingress$compute_peers %in% names(conns))) {
    stop("Registered fresh GLM ingress returned an invalid host set.", call. = FALSE)
  }
  phase20 <- .dsvert_formal_glm_registered_job_run(
    conns[ingress$compute_peers], ingress$hosts,
    .aggregate = .aggregate, max_cycles = max_cycles)
  if (!is.list(phase20) || !identical(names(phase20), c("state", "production_ready")) ||
      !identical(phase20$state, "terminal_complete") ||
      !identical(phase20$production_ready, FALSE)) {
    stop("Registered fresh GLM hosts did not complete their durable handoff.",
         call. = FALSE)
  }
  phase21 <- .dsvert_formal_glm_registered_phase21_run(
    conns[ingress$compute_peers], ingress$hosts,
    .aggregate = .aggregate, max_cycles = max_cycles)
  if (!is.list(phase21) || !identical(phase21, list(
      state = "public_terminal_complete", production_ready = FALSE))) {
    stop("Registered fresh GLM hosts did not complete their durable publication.",
         call. = FALSE)
  }
  list(artifact_id = ingress$artifact_id, total_blocks = blocks,
       state = "public_terminal_complete", production_ready = FALSE)
}
