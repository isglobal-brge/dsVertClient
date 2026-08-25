# Private client relay for a provisioned fresh formal-Cox worker.

.DSVERT_FORMAL_COX_WORKER_CONTROL_REPLY_VERSION <-
  "dsvert-formal-cox-worker-control-response-v1"
.DSVERT_FORMAL_COX_WORKER_CONTROL_MAX_BYTES <- 2L * 1024L * 1024L
.DSVERT_FORMAL_COX_WORKER_CONTROL_ACTIONS <- c(
  "host_start", "bind", "offer", "accept", "confirm", "poll", "relay", "result",
  "completion", "commit")

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
      (action %in% c("host_start", "offer", "completion") && length(fields))) {
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
