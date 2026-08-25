# Private client relay for a provisioned fresh formal-Cox worker.

.DSVERT_FORMAL_COX_WORKER_CONTROL_REPLY_VERSION <-
  "dsvert-formal-cox-worker-control-response-v1"
.DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES <- 2L * 1024L * 1024L
.DSVERT_FORMAL_COX_WORKER_CONTROL_ACTIONS <- c(
  "host_start", "bind", "offer", "accept", "confirm", "poll", "relay", "result",
  "completion", "opening", "commit")

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
