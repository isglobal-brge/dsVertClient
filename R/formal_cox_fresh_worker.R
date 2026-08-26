# Private client relay for a provisioned fresh formal-Cox worker.

.DSVERT_FORMAL_COX_WORKER_CONTROL_REPLY_VERSION <-
  "dsvert-formal-cox-worker-control-response-v1"
.DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES <- 2L * 1024L * 1024L
.DSVERT_FORMAL_COX_WORKER_CONTROL_ACTIONS <- c(
  "host_start", "bind", "offer", "accept", "confirm", "poll", "relay", "result",
  "completion", "opening", "finalizer_ticket", "finalizer_seal", "finalizer_prepare",
  "finalizer_stage", "commit")

.dsvert_formal_cox_fresh_worker_sha256 <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[0-9a-f]{64}$", value)
}

.dsvert_formal_cox_fresh_worker_selector <- function(worker) {
  fields <- c("version", "peer_name", "plan_sha256", "attempt_id",
              "replayed", "production_ready")
  valid <- is.list(worker) && !is.null(names(worker)) && !anyNA(names(worker)) &&
    !anyDuplicated(names(worker)) && setequal(names(worker), fields) &&
    identical(worker$version, "dsvert-formal-cox-blockwise-worker-provision-v1") &&
    is.character(worker$peer_name) && length(worker$peer_name) == 1L &&
    !is.na(worker$peer_name) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", worker$peer_name) &&
    .dsvert_formal_cox_fresh_worker_sha256(worker$plan_sha256) &&
    .dsvert_formal_cox_fresh_worker_sha256(worker$attempt_id) &&
    is.logical(worker$replayed) && length(worker$replayed) == 1L &&
    !is.na(worker$replayed) && identical(worker$production_ready, FALSE)
  if (!isTRUE(valid)) {
    stop("Configured fresh Cox ingress returned an invalid worker selector.",
         call. = FALSE)
  }
  list(peer_name = enc2utf8(worker$peer_name), plan_sha256 = worker$plan_sha256,
       attempt_id = worker$attempt_id)
}

.dsvert_formal_cox_fresh_worker_payload <- function(action, payload) {
  if (!is.list(payload) || (!is.null(attributes(payload)) &&
      length(setdiff(names(attributes(payload)), "names")))) {
    stop("Configured fresh Cox worker requires a protocol payload.",
         call. = FALSE)
  }
  fields <- names(payload)
  if (is.null(fields)) fields <- character()
  if (anyNA(fields) || anyDuplicated(fields) ||
      (action %in% c("host_start", "offer", "completion", "opening") && length(fields))) {
    stop("Configured fresh Cox worker requires a closed action payload.",
         call. = FALSE)
  }
  if (action %in% c("accept", "confirm") &&
      (!identical(fields, "frame") || !is.character(payload$frame) ||
       length(payload$frame) != 1L || is.na(payload$frame) ||
       nchar(payload$frame, type = "bytes") < 4L ||
       nchar(payload$frame, type = "bytes") >
         .DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES ||
       !grepl("^[A-Za-z0-9+/]+={0,2}$", payload$frame))) {
    stop("Configured fresh Cox worker requires a bounded opaque frame.",
         call. = FALSE)
  }
  if (identical(action, "finalizer_ticket") &&
      (!identical(fields, "headers") || !is.list(payload$headers) ||
       length(payload$headers) != 2L || any(vapply(payload$headers, is.null,
                                                    logical(1L))))) {
    stop("Configured fresh Cox worker requires two finalizer headers.",
         call. = FALSE)
  }
  if (identical(action, "finalizer_seal") &&
      (!identical(fields, c("ticket", "headers")) || !is.list(payload$ticket) ||
       !is.list(payload$headers) || length(payload$headers) != 2L ||
       any(vapply(payload$headers, is.null, logical(1L))))) {
    stop("Configured fresh Cox worker requires one ticket and two finalizer headers.",
         call. = FALSE)
  }
  if (action %in% c("finalizer_prepare", "finalizer_stage") &&
      (!identical(fields, c("ticket", "headers", "envelopes")) ||
       !is.list(payload$ticket) || !is.list(payload$headers) ||
       !is.list(payload$envelopes) || length(payload$headers) != 2L ||
       length(payload$envelopes) != 2L ||
       any(vapply(payload$headers, is.null, logical(1L))) ||
       any(vapply(payload$envelopes, is.null, logical(1L))))) {
    stop("Configured fresh Cox worker requires one ticket, two headers and two envelopes.",
         call. = FALSE)
  }
  if (!identical(action, "host_start")) {
    encoded <- tryCatch(.dsvert_joint_dp_client_json(payload),
                        error = function(error) NULL)
    forbidden <- paste0(
      "\\\"[^\\\"]*(?:private|signing_key|secret|storage|path|config|",
      "source|rows|values|master|pid)[^\\\"]*\\\"\\s*:")
    if (!is.character(encoded) || length(encoded) != 1L || is.na(encoded) ||
        nchar(encoded, type = "bytes") < 2L ||
        nchar(encoded, type = "bytes") > .DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES ||
        grepl(forbidden, encoded, perl = TRUE, ignore.case = TRUE)) {
      stop("Configured fresh Cox worker requires a bounded opaque frame.",
           call. = FALSE)
    }
  }
  payload
}

.dsvert_formal_cox_fresh_worker_reply <- function(value, action) {
  fields <- c("version", "action", "payload", "production_ready")
  payload <- if (is.list(value)) value$payload else NULL
  encoded <- tryCatch(.dsvert_joint_dp_client_json(payload),
                      error = function(error) NULL)
  forbidden <- paste0(
    "\\\"[^\\\"]*(?:private|signing_key|secret|storage|path|config|",
    "source|rows|values|master|pid)[^\\\"]*\\\"\\s*:")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version, .DSVERT_FORMAL_COX_WORKER_CONTROL_REPLY_VERSION) &&
    identical(value$action, action) && is.list(value$payload) &&
    identical(value$production_ready, FALSE) && is.character(encoded) &&
    length(encoded) == 1L && !is.na(encoded) &&
    nchar(encoded, type = "bytes") <= .DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES &&
    !grepl(forbidden, encoded, perl = TRUE, ignore.case = TRUE)
  if (!isTRUE(valid)) {
    stop("A server returned an unsafe configured fresh Cox worker reply.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_cox_fresh_worker_call <- function(
    conn, worker, action, payload,
    .aggregate = DSI::datashield.aggregate) {
  if (!is.list(conn) || length(conn) != 1L || is.null(names(conn)) ||
      !is.character(names(conn)) || length(names(conn)) != 1L ||
      is.na(names(conn)) || !nzchar(names(conn))) {
    stop("Configured fresh Cox worker requires one named server.",
         call. = FALSE)
  }
  worker <- .dsvert_formal_cox_fresh_worker_selector(worker)
  if (!identical(names(conn)[[1L]], worker$peer_name)) {
    stop("Configured fresh Cox worker selector belongs to another server.",
         call. = FALSE)
  }
  if (!is.character(action) || length(action) != 1L || is.na(action) ||
      !action %in% .DSVERT_FORMAL_COX_WORKER_CONTROL_ACTIONS) {
    stop("Configured fresh Cox worker requires one closed action.",
         call. = FALSE)
  }
  payload <- .dsvert_formal_cox_fresh_worker_payload(action, payload)
  replies <- .dsvert_aggregate_strict(
    conns = conn,
    expr = call(
      name = "dsvertFormalCoxWorkerControlDS",
      plan_sha256 = worker$plan_sha256, attempt_id = worker$attempt_id,
      action = action, payload = payload),
    operation = "configured fresh formal Cox worker relay",
    .aggregate = .aggregate)
  .dsvert_formal_cox_fresh_worker_reply(replies[[1L]], action)
}

# The completion marker is an internal barrier before the two-authority
# opening.  It intentionally returns no sealed output, coefficient or share.
.dsvert_formal_cox_fresh_worker_completion <- function(
    conns, workers, .aggregate = DSI::datashield.aggregate) {
  valid <- is.list(conns) && is.list(workers) && length(conns) == 2L &&
    length(workers) == 2L && !is.null(names(conns)) && !is.null(names(workers)) &&
    !anyNA(names(conns)) && !anyNA(names(workers)) && !anyDuplicated(names(conns)) &&
    !anyDuplicated(names(workers)) && identical(names(conns), names(workers))
  if (!isTRUE(valid)) {
    stop("Configured fresh Cox completion requires both named compute peers.",
         call. = FALSE)
  }
  parse <- function(reply) {
    payload <- if (is.list(reply)) reply$payload else NULL
    fields <- if (is.list(payload)) names(payload) else NULL
    if (!is.list(payload) || is.null(fields) || anyNA(fields) ||
        anyDuplicated(fields) || !is.logical(payload$complete) ||
        length(payload$complete) != 1L || is.na(payload$complete)) {
      stop("A configured fresh Cox worker returned an invalid completion.",
           call. = FALSE)
    }
    if (!isTRUE(payload$complete)) {
      if (!identical(fields, "complete")) {
        stop("A configured fresh Cox worker returned an invalid completion.",
             call. = FALSE)
      }
      return(NULL)
    }
    completion <- payload$completion
    expected <- c(
      "version", "plan_sha256", "transcript_sha256", "final_commit_sha256",
      "schedule_steps", "fixed_schedule_complete", "output_kind",
      "production_ready", "completion_sha256")
    hashes <- c("plan_sha256", "transcript_sha256", "final_commit_sha256",
                "completion_sha256")
    if (!identical(fields, c("complete", "completion")) ||
        !is.list(completion) || is.null(names(completion)) ||
        anyNA(names(completion)) || anyDuplicated(names(completion)) ||
        !identical(names(completion), expected) ||
        !identical(completion$version,
                   "dsvert-formal-cox-blockwise-completion-v1") ||
        !all(vapply(hashes, function(field)
          .dsvert_formal_cox_fresh_worker_sha256(completion[[field]]),
          logical(1L))) ||
        !is.numeric(completion$schedule_steps) ||
        length(completion$schedule_steps) != 1L ||
        is.na(completion$schedule_steps) || !is.finite(completion$schedule_steps) ||
        completion$schedule_steps < 1L ||
        completion$schedule_steps > .Machine$integer.max ||
        completion$schedule_steps != floor(completion$schedule_steps) ||
        !identical(completion$fixed_schedule_complete, TRUE) ||
        !identical(completion$output_kind, "sealed_private_result_v1") ||
        !identical(completion$production_ready, FALSE)) {
      stop("A configured fresh Cox worker returned an invalid completion.",
           call. = FALSE)
    }
    completion$schedule_steps <- as.integer(completion$schedule_steps)
    completion
  }
  values <- Map(function(conn, worker) {
    parse(.dsvert_formal_cox_fresh_worker_call(
      conn, worker, "completion", structure(list(), names = character()),
      .aggregate = .aggregate))
  }, conns, workers)
  if (xor(is.null(values[[1L]]), is.null(values[[2L]]))) {
    stop("Configured fresh Cox workers disagree about completion.", call. = FALSE)
  }
  if (is.null(values[[1L]])) return(NULL)
  left <- .dsvert_joint_dp_client_json(values[[1L]])
  right <- .dsvert_joint_dp_client_json(values[[2L]])
  if (!identical(left, right)) {
    stop("Configured fresh Cox workers returned different completion records.",
         call. = FALSE)
  }
  values[[1L]]
}

# Read the two share-free opening headers after both workers have committed
# their fixed schedule.  The headers are an internal handoff to the configured
# finalizer: they contain commitments and signatures, never coefficient or
# validity shares.  The finalizer remains the authority which verifies those
# signatures and may perform the single public opening.
.dsvert_formal_cox_fresh_worker_openings <- function(
    conns, workers, completion, .aggregate = DSI::datashield.aggregate) {
  valid_shape <- is.list(conns) && is.list(workers) && length(conns) == 2L &&
    length(workers) == 2L && !is.null(names(conns)) && !is.null(names(workers)) &&
    !anyNA(names(conns)) && !anyNA(names(workers)) && !anyDuplicated(names(conns)) &&
    !anyDuplicated(names(workers)) && identical(names(conns), names(workers))
  if (!isTRUE(valid_shape)) {
    stop("Configured fresh Cox openings require both named compute peers.",
         call. = FALSE)
  }
  completion_json <- tryCatch(.dsvert_joint_dp_client_json(completion),
                              error = function(error) NULL)
  if (!is.character(completion_json) || length(completion_json) != 1L ||
      is.na(completion_json)) {
    stop("Configured fresh Cox openings require one canonical completion.",
         call. = FALSE)
  }
  header_fields <- c(
    "version", "purpose", "artifact_id", "plan_sha256", "run_id",
    "pinset_sha256", "final_cursor", "completion", "final_receipt",
    "peer_name", "peer_id", "role", "coefficient_count", "ring_bits",
    "fraction_bits", "local_beta_validity_sha256", "signature")
  valid_header <- function(value, worker) {
    encoded <- tryCatch(.dsvert_joint_dp_client_json(value),
                        error = function(error) NULL)
    hashes <- c("artifact_id", "plan_sha256", "run_id", "pinset_sha256",
                "local_beta_validity_sha256")
    valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
      !anyDuplicated(names(value)) && setequal(names(value), header_fields) &&
      identical(value$version, "dsvert-formal-cox-blockwise-sticky-opening-v1") &&
      identical(value$purpose, "formal_cox_one_public_beta_validity_opening_v1") &&
      all(vapply(hashes, function(field)
        .dsvert_formal_cox_fresh_worker_sha256(value[[field]]), logical(1L))) &&
      identical(value$peer_name, worker$peer_name) &&
      is.character(value$peer_id) && length(value$peer_id) == 1L &&
      !is.na(value$peer_id) && nzchar(value$peer_id) &&
      is.character(value$role) && length(value$role) == 1L &&
      !is.na(value$role) && value$role %in% c("garbler", "evaluator") &&
      is.numeric(value$coefficient_count) && length(value$coefficient_count) == 1L &&
      is.finite(value$coefficient_count) && value$coefficient_count >= 1L &&
      value$coefficient_count == floor(value$coefficient_count) &&
      is.numeric(value$ring_bits) && length(value$ring_bits) == 1L &&
      value$ring_bits %in% c(63L, 127L) &&
      is.numeric(value$fraction_bits) && length(value$fraction_bits) == 1L &&
      is.finite(value$fraction_bits) && value$fraction_bits >= 0L &&
      value$fraction_bits < value$ring_bits &&
      is.character(value$signature) && length(value$signature) == 1L &&
      !is.na(value$signature) && nchar(value$signature, type = "bytes") >= 32L &&
      nchar(value$signature, type = "bytes") <= 128L &&
      grepl("^[A-Za-z0-9+/]+={0,2}$", value$signature) &&
      !is.null(encoded) && !grepl(
        "\\\"(?:[^\\\"]*(?:share|secret|storage|path|source|value)[^\\\"]*|coefficient_shares)\\\"\\s*:",
        encoded, perl = TRUE, ignore.case = TRUE) &&
      identical(tryCatch(.dsvert_joint_dp_client_json(value$completion),
                         error = function(error) NULL), completion_json)
    if (!isTRUE(valid)) {
      stop("A configured fresh Cox worker returned an unsafe opening header.",
           call. = FALSE)
    }
    value
  }
  openings <- Map(function(conn, worker) {
    reply <- .dsvert_formal_cox_fresh_worker_call(
      conn, worker, "opening", structure(list(), names = character()),
      .aggregate = .aggregate)$payload
    if (!is.list(reply) || !identical(names(reply), c("header", "replayed")) ||
        !is.list(reply$header) || !is.logical(reply$replayed) ||
        length(reply$replayed) != 1L || is.na(reply$replayed)) {
      stop("A configured fresh Cox worker returned an invalid opening reply.",
           call. = FALSE)
    }
    valid_header(reply$header, worker)
  }, conns, workers)
  shared <- c("version", "purpose", "artifact_id", "plan_sha256", "run_id",
              "pinset_sha256", "coefficient_count", "ring_bits", "fraction_bits")
  shared_equal <- vapply(shared, function(field)
    identical(openings[[1L]][[field]], openings[[2L]][[field]]), logical(1L))
  roles <- unname(sort(vapply(openings, `[[`, character(1L), "role")))
  if (!all(shared_equal) || !identical(roles, c("evaluator", "garbler"))) {
    stop("Configured fresh Cox workers returned incompatible opening headers.",
         call. = FALSE)
  }
  names(openings) <- names(workers)
  openings
}

.dsvert_formal_cox_fresh_worker_base64 <- function(value, min_bytes, max_bytes) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    nchar(value, type = "bytes") >= min_bytes &&
    nchar(value, type = "bytes") <= max_bytes &&
    grepl("^[A-Za-z0-9+/]+={0,2}$", value)
}

# Continue a completed fresh worker only as far as the encrypted finalizer
# handoff.  The return value has signed headers, a public ticket and opaque
# ciphertext envelopes; it is deliberately not a Cox estimate or certificate.
.dsvert_formal_cox_fresh_worker_finalizer_handoff <- function(
    conns, workers, completion, .aggregate = DSI::datashield.aggregate) {
  openings <- .dsvert_formal_cox_fresh_worker_openings(
    conns, workers, completion, .aggregate = .aggregate)
  roles <- vapply(openings, `[[`, character(1L), "role")
  finalizer_index <- which(roles == "garbler")
  if (length(finalizer_index) != 1L) {
    stop("Configured fresh Cox headers have no designated finalizer.", call. = FALSE)
  }
  headers <- unname(openings)
  ticket_reply <- .dsvert_formal_cox_fresh_worker_call(
    conns[[finalizer_index]], workers[[finalizer_index]], "finalizer_ticket",
    list(headers = headers), .aggregate = .aggregate)$payload
  ticket_fields <- c(
    "version", "family", "purpose", "artifact_id", "final_pair_root_sha256",
    "plan_sha256", "pinset_sha256", "finalizer_peer_name", "finalizer_peer_id",
    "finalizer_role", "recipient_x25519_public_key", "transport_key_sha256",
    "issuer_peer_name", "issuer_peer_id", "signature")
  finalizer <- openings[[finalizer_index]]
  valid_ticket <- is.list(ticket_reply) && identical(names(ticket_reply),
                                                       c("ticket", "replayed")) &&
    is.list(ticket_reply$ticket) && !is.null(names(ticket_reply$ticket)) &&
    identical(names(ticket_reply$ticket), ticket_fields) &&
    is.logical(ticket_reply$replayed) && length(ticket_reply$replayed) == 1L &&
    !is.na(ticket_reply$replayed) && identical(ticket_reply$ticket$version,
                                                "dsvert-typed-finalizer-handoff-ticket-v1") &&
    identical(ticket_reply$ticket$family, "formal_cox") &&
    identical(ticket_reply$ticket$purpose, "formal_cox_blockwise_sticky_opening_v1") &&
    identical(ticket_reply$ticket$artifact_id, finalizer$artifact_id) &&
    identical(ticket_reply$ticket$plan_sha256, finalizer$plan_sha256) &&
    identical(ticket_reply$ticket$pinset_sha256, finalizer$pinset_sha256) &&
    identical(ticket_reply$ticket$finalizer_peer_name, finalizer$peer_name) &&
    identical(ticket_reply$ticket$finalizer_peer_id, finalizer$peer_id) &&
    identical(ticket_reply$ticket$finalizer_role, "garbler") &&
    identical(ticket_reply$ticket$issuer_peer_name, finalizer$peer_name) &&
    identical(ticket_reply$ticket$issuer_peer_id, finalizer$peer_id) &&
    all(vapply(c("final_pair_root_sha256", "transport_key_sha256"), function(field)
      .dsvert_formal_cox_fresh_worker_sha256(ticket_reply$ticket[[field]]), logical(1L))) &&
    .dsvert_formal_cox_fresh_worker_base64(
      ticket_reply$ticket$recipient_x25519_public_key, 40L, 48L) &&
    .dsvert_formal_cox_fresh_worker_base64(ticket_reply$ticket$signature, 80L, 96L)
  if (!isTRUE(valid_ticket)) {
    stop("Configured fresh Cox worker returned an invalid finalizer ticket.",
         call. = FALSE)
  }
  ticket <- ticket_reply$ticket
  envelope_fields <- c(
    "version", "family", "purpose", "artifact_id", "final_pair_root_sha256",
    "plan_sha256", "pinset_sha256", "ticket_sha256", "finalizer_peer_name",
    "finalizer_peer_id", "recipient_transport_key_sha256", "sender_peer_name",
    "sender_peer_id", "sender_role", "payload_kind", "payload_sha256",
    "ciphertext_sha256", "ciphertext", "signature")
  envelopes <- Map(function(conn, worker, header) {
    reply <- .dsvert_formal_cox_fresh_worker_call(
      conn, worker, "finalizer_seal", list(ticket = ticket, headers = headers),
      .aggregate = .aggregate)$payload
    valid <- is.list(reply) && identical(names(reply), c("envelope", "replayed")) &&
      is.list(reply$envelope) && !is.null(names(reply$envelope)) &&
      identical(names(reply$envelope), envelope_fields) &&
      is.logical(reply$replayed) && length(reply$replayed) == 1L && !is.na(reply$replayed) &&
      identical(reply$envelope$version, "dsvert-typed-finalizer-handoff-envelope-v1") &&
      identical(reply$envelope$family, ticket$family) &&
      identical(reply$envelope$purpose, ticket$purpose) &&
      identical(reply$envelope$artifact_id, ticket$artifact_id) &&
      identical(reply$envelope$final_pair_root_sha256, ticket$final_pair_root_sha256) &&
      identical(reply$envelope$plan_sha256, ticket$plan_sha256) &&
      identical(reply$envelope$pinset_sha256, ticket$pinset_sha256) &&
      identical(reply$envelope$finalizer_peer_name, ticket$finalizer_peer_name) &&
      identical(reply$envelope$finalizer_peer_id, ticket$finalizer_peer_id) &&
      identical(reply$envelope$recipient_transport_key_sha256, ticket$transport_key_sha256) &&
      identical(reply$envelope$sender_peer_name, header$peer_name) &&
      identical(reply$envelope$sender_peer_id, header$peer_id) &&
      identical(reply$envelope$sender_role, header$role) &&
      identical(reply$envelope$payload_kind, "formal_cox_opening_local_output_v1") &&
      all(vapply(c("ticket_sha256", "payload_sha256", "ciphertext_sha256"), function(field)
        .dsvert_formal_cox_fresh_worker_sha256(reply$envelope[[field]]), logical(1L))) &&
      .dsvert_formal_cox_fresh_worker_base64(reply$envelope$ciphertext, 80L,
                                              .DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES) &&
      .dsvert_formal_cox_fresh_worker_base64(reply$envelope$signature, 80L, 96L)
    if (!isTRUE(valid)) {
      stop("Configured fresh Cox worker returned an invalid finalizer envelope.",
           call. = FALSE)
    }
    reply$envelope
  }, conns, workers, openings)
  names(envelopes) <- names(workers)
  list(headers = openings, ticket = ticket, envelopes = envelopes,
       production_ready = FALSE)
}

# Persist the finalizer's candidate intent from two opaque handoffs. This is
# still not a Cox estimate or public certificate: the sole return is a
# share-free intent bound to the ticket, or a digest of an already-public
# certificate on durable replay.
.dsvert_formal_cox_fresh_worker_prepare_finalizer <- function(
    conns, workers, handoff, .aggregate = DSI::datashield.aggregate) {
  expected_handoff <- c("headers", "ticket", "envelopes", "production_ready")
  valid_handoff <- is.list(handoff) && !is.null(names(handoff)) &&
    identical(names(handoff), expected_handoff) && is.list(handoff$headers) &&
    is.list(handoff$ticket) && is.list(handoff$envelopes) &&
    length(handoff$headers) == 2L && length(handoff$envelopes) == 2L &&
    identical(handoff$production_ready, FALSE)
  if (!isTRUE(valid_handoff) || !is.list(conns) || !is.list(workers) ||
      length(conns) != 2L || length(workers) != 2L ||
      !identical(names(conns), names(workers))) {
    stop("Configured fresh Cox finalizer preparation requires one validated handoff.",
         call. = FALSE)
  }
  roles <- vapply(handoff$headers, function(header) {
    if (is.list(header)) header$role else NA_character_
  }, character(1L))
  finalizer_index <- which(roles == "garbler")
  if (length(finalizer_index) != 1L ||
      !identical(handoff$ticket$finalizer_peer_name,
                 workers[[finalizer_index]]$peer_name) ||
      !identical(handoff$ticket$finalizer_role, "garbler")) {
    stop("Configured fresh Cox handoff has no designated finalizer.", call. = FALSE)
  }
  payload <- .dsvert_formal_cox_fresh_worker_call(
    conns[[finalizer_index]], workers[[finalizer_index]], "finalizer_prepare",
    list(ticket = handoff$ticket, headers = unname(handoff$headers),
         envelopes = unname(handoff$envelopes)), .aggregate = .aggregate)$payload
  fields <- c("intent", "finalized", "certificate_sha256", "replayed")
  valid_common <- is.list(payload) && identical(names(payload), fields) &&
    is.logical(payload$finalized) && length(payload$finalized) == 1L &&
    !is.na(payload$finalized) && is.logical(payload$replayed) &&
    length(payload$replayed) == 1L && !is.na(payload$replayed) &&
    is.character(payload$certificate_sha256) &&
    length(payload$certificate_sha256) == 1L && !is.na(payload$certificate_sha256)
  if (!isTRUE(valid_common)) {
    stop("Configured fresh Cox worker returned an invalid finalizer preparation.",
         call. = FALSE)
  }
  if (isTRUE(payload$finalized)) {
    if (!is.null(payload$intent) ||
        !.dsvert_formal_cox_fresh_worker_sha256(payload$certificate_sha256)) {
      stop("Configured fresh Cox worker returned an invalid finalizer replay.",
           call. = FALSE)
    }
  } else {
    intent_fields <- c("version", "purpose", "artifact_id", "candidate_sha256",
                       "final_pair_root_sha256", "opening_mode", "exp_postprocess_mode")
    intent <- payload$intent
    valid_intent <- is.list(intent) && !is.null(names(intent)) &&
      identical(names(intent), intent_fields) &&
      identical(intent$version, "dsvert-formal-cox-blockwise-sticky-opening-v1") &&
      identical(intent$purpose, "formal_cox_one_public_beta_validity_opening_v1") &&
      identical(intent$artifact_id, handoff$ticket$artifact_id) &&
      identical(intent$final_pair_root_sha256,
                handoff$ticket$final_pair_root_sha256) &&
      .dsvert_formal_cox_fresh_worker_sha256(intent$candidate_sha256) &&
      identical(intent$opening_mode, "dual_authority_additive_ring_and_xor_validity_v1") &&
      identical(intent$exp_postprocess_mode,
                "certified_dyadic_interval_midpoint_v1") &&
      identical(payload$certificate_sha256, "")
    if (!isTRUE(valid_intent)) {
      stop("Configured fresh Cox worker returned an invalid finalizer intent.",
           call. = FALSE)
    }
  }
  list(intent = payload$intent, finalized = payload$finalized,
       certificate_sha256 = payload$certificate_sha256,
       replayed = payload$replayed, production_ready = FALSE)
}

# Stage the already prepared finalizer records at both compute authorities.
# This is still not an opening: it verifies only the role-local, share-free
# record each worker persisted in Rock and returns no candidate to the caller.
.dsvert_formal_cox_fresh_worker_stage_finalizer <- function(
    conns, workers, handoff, prepared, .aggregate = DSI::datashield.aggregate) {
  expected_handoff <- c("headers", "ticket", "envelopes", "production_ready")
  expected_prepared <- c("intent", "finalized", "certificate_sha256", "replayed",
                         "production_ready")
  valid <- is.list(conns) && is.list(workers) && length(conns) == 2L &&
    length(workers) == 2L && identical(names(conns), names(workers)) &&
    is.list(handoff) && identical(names(handoff), expected_handoff) &&
    is.list(handoff$headers) && length(handoff$headers) == 2L &&
    is.list(handoff$ticket) && is.list(handoff$envelopes) &&
    length(handoff$envelopes) == 2L &&
    identical(handoff$production_ready, FALSE) &&
    is.list(prepared) && identical(names(prepared), expected_prepared) &&
    identical(prepared$finalized, FALSE) && is.list(prepared$intent) &&
    .dsvert_formal_cox_fresh_worker_sha256(prepared$intent$artifact_id) &&
    .dsvert_formal_cox_fresh_worker_sha256(prepared$intent$candidate_sha256) &&
    identical(prepared$intent$artifact_id, handoff$ticket$artifact_id)
  if (!isTRUE(valid)) {
    stop("Configured fresh Cox finalizer staging requires one prepared handoff.",
         call. = FALSE)
  }
  stages <- Map(function(conn, worker, header) {
    payload <- .dsvert_formal_cox_fresh_worker_call(
      conn, worker, "finalizer_stage",
      list(ticket = handoff$ticket, headers = unname(handoff$headers),
           envelopes = unname(handoff$envelopes)), .aggregate = .aggregate)$payload
    fields <- c("artifact_id", "candidate_sha256", "local_role", "production_ready")
    valid_stage <- is.list(payload) && identical(names(payload), fields) &&
      identical(payload$artifact_id, prepared$intent$artifact_id) &&
      identical(payload$local_role, header$role) &&
      identical(payload$production_ready, FALSE) &&
      is.character(payload$candidate_sha256) && length(payload$candidate_sha256) == 1L &&
      !is.na(payload$candidate_sha256) &&
      if (identical(header$role, "garbler")) {
        identical(payload$candidate_sha256, prepared$intent$candidate_sha256)
      } else {
        identical(header$role, "evaluator") && identical(payload$candidate_sha256, "")
      }
    if (!isTRUE(valid_stage)) {
      stop("Configured fresh Cox worker returned an invalid finalizer stage.",
           call. = FALSE)
    }
    payload
  }, conns, workers, handoff$headers)
  names(stages) <- names(workers)
  list(artifact_id = prepared$intent$artifact_id, production_ready = FALSE)
}

# Run the fixed two-peer schedule by forwarding only authenticated opaque
# frames.  The client never decodes MPC payloads or turns the completion
# marker into a public result; the opening lifecycle owns that later boundary.
.dsvert_formal_cox_fresh_worker_run <- function(
    conns, workers, .aggregate = DSI::datashield.aggregate) {
  valid <- is.list(conns) && is.list(workers) && length(conns) == 2L &&
    length(workers) == 2L && !is.null(names(conns)) && !is.null(names(workers)) &&
    !anyNA(names(conns)) && !anyNA(names(workers)) && !anyDuplicated(names(conns)) &&
    !anyDuplicated(names(workers)) && identical(names(conns), names(workers))
  if (!isTRUE(valid)) {
    stop("Configured fresh Cox execution requires both named compute peers.",
         call. = FALSE)
  }
  empty <- structure(list(), names = character())
  call <- function(index, action, payload) {
    .dsvert_formal_cox_fresh_worker_call(
      conns[index], workers[[index]], action, payload, .aggregate = .aggregate)
  }
  payload <- function(reply, expected, action) {
    value <- if (is.list(reply)) reply$payload else NULL
    fields <- if (is.list(value)) names(value) else NULL
    if (is.null(fields)) fields <- character()
    if (!is.list(value) || is.null(fields) || anyNA(fields) ||
        anyDuplicated(fields) || !identical(fields, expected)) {
      stop(paste0("A configured fresh Cox worker returned an invalid ", action,
                  " reply."), call. = FALSE)
    }
    value
  }
  offset <- function(value, field) {
    valid_offset <- is.character(value) && length(value) == 1L && !is.na(value) &&
      grepl("^(0|[1-9][0-9]{0,18})$", value) &&
      (nchar(value, type = "bytes") < 19L || value <= "9223372036854775807")
    if (!isTRUE(valid_offset)) {
      stop(paste0("A configured fresh Cox worker returned an invalid ", field,
                  "."), call. = FALSE)
    }
    value
  }
  frame <- function(reply, action) {
    value <- payload(reply, "frame", action)$frame
    valid_frame <- is.character(value) && length(value) == 1L && !is.na(value) &&
      nchar(value, type = "bytes") >= 4L &&
      nchar(value, type = "bytes") <= .DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES &&
      grepl("^[A-Za-z0-9+/]+={0,2}$", value)
    if (!isTRUE(valid_frame)) {
      stop("A configured fresh Cox worker returned an invalid root frame.",
           call. = FALSE)
    }
    value
  }
  poll <- function(reply) {
    value <- payload(reply, c("chunk", "accepted"), "poll")
    value$accepted <- offset(value$accepted, "poll acknowledgement")
    if (is.null(value$chunk)) return(value)
    chunk <- value$chunk
    expected <- c("sender", "offset", "payload_sha256", "payload")
    valid_chunk <- is.list(chunk) && !is.null(names(chunk)) &&
      !anyNA(names(chunk)) && !anyDuplicated(names(chunk)) &&
      identical(names(chunk), expected) && is.character(chunk$sender) &&
      length(chunk$sender) == 1L && !is.na(chunk$sender) &&
      grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", chunk$sender) &&
      .dsvert_formal_cox_fresh_worker_sha256(chunk$payload_sha256) &&
      is.character(chunk$payload) && length(chunk$payload) == 1L &&
      !is.na(chunk$payload) && nchar(chunk$payload, type = "bytes") >= 4L &&
      nchar(chunk$payload, type = "bytes") <=
        .DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES &&
      grepl("^[A-Za-z0-9+/]+={0,2}$", chunk$payload)
    if (!isTRUE(valid_chunk)) {
      stop("A configured fresh Cox worker returned an invalid relay chunk.",
           call. = FALSE)
    }
    chunk$offset <- offset(chunk$offset, "relay offset")
    value$chunk <- chunk
    value
  }
  relay <- function(reply) {
    offset(payload(reply, "accepted", "relay")$accepted, "relay acknowledgement")
  }
  result <- function(reply) {
    value <- payload(reply, c("receipt", "done"), "result")
    if (!is.logical(value$done) || length(value$done) != 1L || is.na(value$done) ||
        (!isTRUE(value$done) && !isFALSE(value$done)) ||
        !is.list(value$receipt)) {
      stop("A configured fresh Cox worker returned an invalid worker result.",
           call. = FALSE)
    }
    value
  }
  for (index in seq_along(workers)) call(index, "host_start", empty)
  for (index in seq_along(workers)) {
    other <- names(workers)[[3L - index]]
    payload(call(index, "bind", list(peer = other)), character(), "bind")
  }
  acknowledgements <- c("0", "0")
  repeat {
    completion <- .dsvert_formal_cox_fresh_worker_completion(
      conns, workers, .aggregate = .aggregate)
    if (!is.null(completion)) return(completion)
    offer <- frame(call(1L, "offer", empty), "offer")
    accept <- frame(call(2L, "accept", list(frame = offer)), "accept")
    payload(call(1L, "confirm", list(frame = accept)), character(), "confirm")
    done <- c(FALSE, FALSE)
    receipts <- vector("list", 2L)
    repeat {
      for (direction in list(c(1L, 2L), c(2L, 1L))) {
        sent <- poll(call(direction[[1L]], "poll",
                          list(acknowledged = acknowledgements[[direction[[1L]]]])))
        if (!is.null(sent$chunk)) {
          acknowledgements[[direction[[1L]]]] <- relay(call(
            direction[[2L]], "relay", list(chunk = sent$chunk)))
        }
      }
      for (index in seq_along(workers)) {
        observed <- result(call(index, "result", empty))
        if (isTRUE(observed$done)) {
          receipts[[index]] <- observed$receipt
          done[[index]] <- TRUE
        }
      }
      if (all(done)) break
      Sys.sleep(0.01)
    }
    for (index in seq_along(workers)) {
      payload(call(index, "commit", list(receipts = unname(receipts))),
              character(), "commit")
    }
  }
}
