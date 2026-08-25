# Private client seam for the closed registered formal-GLM host relay.
#
# It is intentionally not an analytical frontdoor.  A later fresh-GLM driver
# uses it to move only authenticated peer frames between provisioned hosts.

.DSVERT_FORMAL_GLM_REGISTERED_JOB_CONTROL_REPLY_VERSION <-
  "dsvert-formal-glm-registered-phase20-job-control-response-v1"

.dsvert_formal_glm_registered_job_receipt <- function(value) {
  fields <- c(
    "version", "peer", "artifact_id", "receipt_set_sha256", "config_sha256",
    "replayed", "production_ready")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version,
              "dsvert-formal-glm-registered-phase20-job-host-provision-v1") &&
    is.character(value$peer) && length(value$peer) == 1L && !is.na(value$peer) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value$peer) &&
    all(vapply(c("artifact_id", "receipt_set_sha256", "config_sha256"),
               function(field) is.character(value[[field]]) &&
                 length(value[[field]]) == 1L && !is.na(value[[field]]) &&
                 grepl("^[0-9a-f]{64}$", value[[field]]), logical(1L))) &&
    is.logical(value$replayed) && length(value$replayed) == 1L &&
    !is.na(value$replayed) && identical(value$production_ready, FALSE)
  if (!isTRUE(valid)) {
    stop("Registered formal GLM job control requires one provision receipt.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_glm_registered_job_control_reply <- function(value, action) {
  fields <- c("version", "action", "payload", "production_ready")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version,
              .DSVERT_FORMAL_GLM_REGISTERED_JOB_CONTROL_REPLY_VERSION) &&
    identical(value$action, action) && is.list(value$payload) &&
    identical(value$production_ready, FALSE)
  if (!isTRUE(valid)) {
    stop("A server returned an invalid registered formal GLM control reply.",
         call. = FALSE)
  }
  encoded <- tryCatch(.dsvert_joint_dp_client_json(value$payload),
                      error = function(error) NULL)
  forbidden <- paste0(
    "\\\"[^\\\"]*(?:config|path|key|secret|socket|storage|backend|",
    "share|plan|context|token|capsule|run_id|opening)[^\\\"]*\\\"\\s*:")
  if (!is.character(encoded) || length(encoded) != 1L || is.na(encoded) ||
      nchar(encoded, type = "bytes") > 2L * 1024L * 1024L ||
      grepl(forbidden, encoded, perl = TRUE, ignore.case = TRUE)) {
    stop("A server returned an unsafe registered formal GLM control reply.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_glm_registered_job_control_call <- function(
    conn, receipt, action, payload,
    .aggregate = DSI::datashield.aggregate) {
  if (!is.list(conn) || length(conn) != 1L || is.null(names(conn)) ||
      !is.character(names(conn)) || length(names(conn)) != 1L ||
      is.na(names(conn)) || !nzchar(names(conn))) {
    stop("Registered formal GLM job control requires one named server.",
         call. = FALSE)
  }
  receipt <- .dsvert_formal_glm_registered_job_receipt(receipt)
  actions <- c("negotiate", "start", "health", "job_ref", "bind",
               "heartbeat", "poll", "relay", "compute", "terminal",
               "compute_start", "compute_status", "terminal_start",
               "terminal_status",
               "phase21_preflight", "phase21_preflight_bind",
               "phase21_stage_start", "phase21_stage_status",
               "phase21_stage_record", "phase21_stage_poll",
               "phase21_stage_relay", "phase21_stage_import",
               "phase21_ticket", "phase21_ticket_import", "phase21_seal",
               "phase21_seal_import", "phase21_candidate",
               "phase21_candidate_import", "phase21_candidate_verify",
               "phase21_local_release_import", "phase21_base_certificate",
               "phase21_base_certificate_import", "phase21_authorization",
               "phase21_authorization_import", "phase21_publication",
               "phase21_commit", "phase21_commit_import", "phase21_ack",
               "phase21_ack_import", "phase21_cleanup",
               "phase21_cleanup_import")
  if (!is.character(action) || length(action) != 1L || is.na(action) ||
      !action %in% actions || !is.list(payload)) {
    stop("Registered formal GLM job control requires one closed action and payload.",
         call. = FALSE)
  }
  if (startsWith(action, "phase21_")) {
    if (!identical(names(payload), "frame")) {
      stop("Registered formal GLM Phase21 control requires one opaque frame.",
           call. = FALSE)
    }
    payload$frame <- .dsvert_formal_glm_registered_job_base64(payload$frame)
  }
  replies <- .dsvert_aggregate_strict(
    conns = conn,
    expr = call(
      name = "dsvertFormalGLMRegisteredJobControlDS",
      receipt = receipt, action = action, payload = payload),
    operation = "registered formal GLM job control relay",
    .aggregate = .aggregate)
  .dsvert_formal_glm_registered_job_control_reply(replies[[1L]], action)
}

# This private coordinator deliberately has no analytical result.  It only
# drives the two pre-provisioned Phase20 hosts until both encrypted handoffs
# are durable; Phase21 remains the sole public-release authority.
.dsvert_formal_glm_registered_job_empty <- function() {
  structure(list(), names = character())
}

.dsvert_formal_glm_registered_job_sha256 <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[0-9a-f]{64}$", value)
}

.dsvert_formal_glm_registered_job_base64 <- function(value, allow_empty = FALSE) {
  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    nchar(value, type = "bytes") <= 2L * 1024L * 1024L &&
    (isTRUE(allow_empty) || nzchar(value)) &&
    grepl("^(?:[A-Za-z0-9+/]{4})*(?:[A-Za-z0-9+/]{2}==|[A-Za-z0-9+/]{3}=)?$",
          value, perl = TRUE)
  if (!isTRUE(valid)) {
    stop("A registered formal GLM host returned an invalid opaque frame.",
         call. = FALSE)
  }
  enc2utf8(value)
}

.dsvert_formal_glm_registered_job_flag <- function(value, field) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(paste0("A registered formal GLM host returned an invalid ", field, "."),
         call. = FALSE)
  }
  value
}

.dsvert_formal_glm_registered_job_offset <- function(value, field) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value >= 0 && value <= (2^53 - 1) &&
    value == floor(value)
  if (!isTRUE(valid)) {
    stop(paste0("A registered formal GLM host returned an invalid ", field, "."),
         call. = FALSE)
  }
  value
}

.dsvert_formal_glm_registered_job_ref <- function(value) {
  fields <- c("artifact_id", "receipt_set_sha256", "attempt_id", "job_sha256",
              "transport_sha256", "production_ready")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    all(vapply(fields[seq_len(5L)], function(field)
      .dsvert_formal_glm_registered_job_sha256(value[[field]]), logical(1L))) &&
    identical(.dsvert_formal_glm_registered_job_flag(
      value$production_ready, "production flag"), FALSE)
  if (!isTRUE(valid)) {
    stop("A registered formal GLM host returned an invalid job reference.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_glm_registered_job_chunk <- function(value, ref) {
  fields <- c("job_sha256", "transport_sha256", "offset", "payload_sha256", "payload")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$job_sha256, ref$job_sha256) &&
    identical(value$transport_sha256, ref$transport_sha256) &&
    .dsvert_formal_glm_registered_job_sha256(value$payload_sha256)
  if (!isTRUE(valid)) {
    stop("A registered formal GLM host returned an invalid relay chunk.",
         call. = FALSE)
  }
  value$offset <- .dsvert_formal_glm_registered_job_offset(value$offset, "relay offset")
  value$payload <- .dsvert_formal_glm_registered_job_base64(value$payload)
  value
}

.dsvert_formal_glm_registered_job_result <- function(value, expected_state) {
  fields <- c("state", "outbound", "inspect_only", "production_ready")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$state, expected_state) &&
    identical(.dsvert_formal_glm_registered_job_flag(
      value$inspect_only, "inspect-only flag"), FALSE) &&
    identical(.dsvert_formal_glm_registered_job_flag(
      value$production_ready, "production flag"), FALSE)
  if (!isTRUE(valid)) {
    stop("A registered formal GLM host is not in the expected durable state.",
         call. = FALSE)
  }
  value$outbound <- .dsvert_formal_glm_registered_job_base64(
    value$outbound, allow_empty = TRUE)
  value
}

.dsvert_formal_glm_registered_job_ref_claim <- function(value) {
  fields <- c("ref", "claim")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields)
  if (!isTRUE(valid)) {
    stop("A registered formal GLM host returned an invalid peer claim.",
         call. = FALSE)
  }
  list(ref = .dsvert_formal_glm_registered_job_ref(value$ref),
       claim = .dsvert_formal_glm_registered_job_base64(value$claim))
}

.dsvert_formal_glm_registered_job_poll <- function(value, ref) {
  fields <- names(value)
  valid_fields <- !is.null(fields) && !anyNA(fields) && !anyDuplicated(fields) &&
    setequal(fields, c("state", "accepted_through", "production_ready")) ||
    setequal(fields, c("state", "accepted_through", "relay_chunk", "production_ready"))
  valid <- is.list(value) && isTRUE(valid_fields) &&
    identical(value$state, "running") &&
    identical(.dsvert_formal_glm_registered_job_flag(
      value$production_ready, "production flag"), FALSE)
  if (!isTRUE(valid)) {
    stop("A registered formal GLM host returned an invalid poll result.",
         call. = FALSE)
  }
  value$accepted_through <- .dsvert_formal_glm_registered_job_offset(
    value$accepted_through, "accepted offset")
  if (!is.null(value$relay_chunk)) {
    value$relay_chunk <- .dsvert_formal_glm_registered_job_chunk(value$relay_chunk, ref)
  }
  value
}

.dsvert_formal_glm_registered_job_task <- function(value, terminal = FALSE) {
  fields <- names(value)
  valid_fields <- !is.null(fields) && !anyNA(fields) && !anyDuplicated(fields) &&
    (setequal(fields, c("state", "production_ready")) ||
       (isTRUE(terminal) && setequal(fields, c("state", "commit", "production_ready"))))
  valid <- is.list(value) && isTRUE(valid_fields) &&
    is.character(value$state) && length(value$state) == 1L &&
    value$state %in% c("running", "complete") &&
    identical(.dsvert_formal_glm_registered_job_flag(
      value$production_ready, "production flag"), FALSE)
  if (!isTRUE(valid)) {
    stop("A registered formal GLM host returned an invalid task status.",
         call. = FALSE)
  }
  if (identical(value$state, "complete") && isTRUE(terminal)) {
    commit <- value$commit
    commit_fields <- tolower(names(commit))
    commit_valid <- is.list(commit) && !is.null(names(commit)) &&
      !anyNA(names(commit)) && !anyDuplicated(commit_fields) &&
      setequal(commit_fields, c("sha256", "bytes", "replayed")) &&
      .dsvert_formal_glm_registered_job_sha256(
        commit[[which(commit_fields == "sha256")]]) &&
      .dsvert_formal_glm_registered_job_offset(
        commit[[which(commit_fields == "bytes")]], "handoff bytes") >= 64 &&
      is.logical(commit[[which(commit_fields == "replayed")]]) &&
      length(commit[[which(commit_fields == "replayed")]]) == 1L &&
      !is.na(commit[[which(commit_fields == "replayed")]])
    if (!isTRUE(commit_valid)) {
      stop("A registered formal GLM host returned an invalid handoff commitment.",
           call. = FALSE)
    }
  } else if (!is.null(value$commit)) {
    stop("A registered formal GLM host returned an early handoff commitment.",
         call. = FALSE)
  }
  value
}

#' Drive two provisioned registered-formal-GLM hosts to their durable handoff
#'
#' This is private orchestration for the fresh formal-GLM route. It exchanges
#' only server-authenticated frames and opaque relay chunks. It does not read
#' a result, publish a certificate, or make any production-ready assertion.
.dsvert_formal_glm_registered_job_run <- function(
    conns, receipts, .aggregate = DSI::datashield.aggregate, max_cycles = NULL) {
  valid_conns <- is.list(conns) && length(conns) == 2L && !is.null(names(conns)) &&
    !anyNA(names(conns)) && !anyDuplicated(names(conns)) && all(nzchar(names(conns)))
  valid_receipts <- is.list(receipts) && identical(names(receipts), names(conns))
  if (!isTRUE(valid_conns) || !isTRUE(valid_receipts)) {
    stop("Registered formal GLM job orchestration requires exactly two named hosts.",
         call. = FALSE)
  }
  if (!is.null(max_cycles) && (!is.numeric(max_cycles) || length(max_cycles) != 1L ||
      is.na(max_cycles) || !is.finite(max_cycles) || max_cycles < 1 ||
      max_cycles != floor(max_cycles))) {
    stop("Registered formal GLM job orchestration received an invalid test cycle limit.",
         call. = FALSE)
  }
  lapply(receipts, .dsvert_formal_glm_registered_job_receipt)
  peers <- names(conns)
  call <- function(peer, action, payload) {
    reply <- .dsvert_formal_glm_registered_job_control_call(
      conns[peer], receipts[peer], action, payload, .aggregate = .aggregate)
    .dsvert_formal_glm_registered_job_control_reply(reply, action)$payload
  }
  proposal <- .dsvert_formal_glm_registered_job_result(
    call(peers[[1L]], "negotiate", list(inbound = "")), "negotiating")
  if (!nzchar(proposal$outbound)) {
    stop("The registered formal GLM garbler did not emit a proposal.", call. = FALSE)
  }
  accept <- .dsvert_formal_glm_registered_job_result(
    call(peers[[2L]], "negotiate", list(inbound = proposal$outbound)), "negotiating")
  if (!nzchar(accept$outbound)) {
    stop("The registered formal GLM evaluator did not emit an acceptance.", call. = FALSE)
  }
  .dsvert_formal_glm_registered_job_result(
    call(peers[[1L]], "negotiate", list(inbound = accept$outbound)), "negotiating")
  refs <- lapply(peers, function(peer) {
    .dsvert_formal_glm_registered_job_result(
      call(peer, "start", .dsvert_formal_glm_registered_job_empty()), "running")
    .dsvert_formal_glm_registered_job_ref_claim(
      call(peer, "job_ref", .dsvert_formal_glm_registered_job_empty()))
  })
  names(refs) <- peers
  if (!identical(refs[[1L]]$ref, refs[[2L]]$ref)) {
    stop("Registered formal GLM hosts derived different job references.", call. = FALSE)
  }
  ref <- refs[[1L]]$ref
  for (index in seq_along(peers)) {
    peer <- peers[[index]]
    other <- peers[[3L - index]]
    call(peer, "bind", list(frame = refs[[other]]$claim))
    call(peer, "heartbeat", .dsvert_formal_glm_registered_job_empty())
  }
  acknowledged <- stats::setNames(c(0, 0), peers)
  cycle <- 0
  relay_once <- function() {
    nonlocal_progress <- FALSE
    for (index in seq_along(peers)) {
      peer <- peers[[index]]
      other <- peers[[3L - index]]
      poll <- .dsvert_formal_glm_registered_job_poll(
        call(peer, "poll", list(ref = ref, acknowledged = acknowledged[[peer]])), ref)
      if (!is.null(poll$relay_chunk)) {
        accepted <- call(other, "relay", list(ref = ref, chunk = poll$relay_chunk))
        valid_accepted <- is.list(accepted) && identical(names(accepted), "accepted")
        if (!isTRUE(valid_accepted)) {
          stop("A registered formal GLM host returned an invalid relay acknowledgement.",
               call. = FALSE)
        }
        acknowledged[[peer]] <<- .dsvert_formal_glm_registered_job_offset(
          accepted$accepted, "relay acknowledgement")
        nonlocal_progress <- TRUE
      }
    }
    nonlocal_progress
  }
  repeat {
    cycle <- cycle + 1
    progress <- relay_once()
    if (!progress) break
    if (!is.null(max_cycles) && cycle >= max_cycles) {
      stop("Registered formal GLM job relay did not become idle in the requested test cycles.",
           call. = FALSE)
    }
  }
  for (peer in peers) {
    .dsvert_formal_glm_registered_job_task(
      call(peer, "compute_start", .dsvert_formal_glm_registered_job_empty()))
  }
  repeat {
    complete <- vapply(peers, function(peer) {
      call(peer, "heartbeat", .dsvert_formal_glm_registered_job_empty())
      status <- .dsvert_formal_glm_registered_job_task(
        call(peer, "compute_status", .dsvert_formal_glm_registered_job_empty()))
      identical(status$state, "complete")
    }, logical(1L))
    progress <- relay_once()
    if (all(complete) && !progress) break
    if (!progress) Sys.sleep(0.01)
  }
  for (peer in peers) {
    .dsvert_formal_glm_registered_job_task(
      call(peer, "terminal_start", .dsvert_formal_glm_registered_job_empty()), terminal = TRUE)
  }
  repeat {
    complete <- vapply(peers, function(peer) {
      call(peer, "heartbeat", .dsvert_formal_glm_registered_job_empty())
      status <- .dsvert_formal_glm_registered_job_task(
        call(peer, "terminal_status", .dsvert_formal_glm_registered_job_empty()),
        terminal = TRUE)
      identical(status$state, "complete")
    }, logical(1L))
    progress <- relay_once()
    if (all(complete) && !progress) break
    if (!progress) Sys.sleep(0.01)
  }
  list(state = "terminal_complete", production_ready = FALSE)
}
