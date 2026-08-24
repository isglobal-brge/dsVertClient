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
               "heartbeat", "poll", "relay", "compute", "terminal")
  if (!is.character(action) || length(action) != 1L || is.na(action) ||
      !action %in% actions || !is.list(payload)) {
    stop("Registered formal GLM job control requires one closed action and payload.",
         call. = FALSE)
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
