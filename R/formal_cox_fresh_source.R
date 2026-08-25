# Private client relay for a configured fresh formal-Cox source.
#
# This moves only selectors, signed tickets and encrypted envelopes.  It is
# deliberately not a Cox frontdoor and never carries a row, model control or
# public result.

.DSVERT_FORMAL_COX_FRESH_SOURCE_REPLY_VERSION <-
  "dsvert-formal-cox-fresh-source-response-v1"
.DSVERT_FORMAL_COX_FRESH_SOURCE_MAX_BYTES <- 4L * 1024L * 1024L

.dsvert_formal_cox_fresh_source_selector <- function(
    analysis_id, data_name, formula_sha256) {
  analysis_id <- .dsvert_formal_cox_frontdoor_id(analysis_id, "analysis_id")
  data_name <- .dsvert_formal_cox_frontdoor_id(data_name, "data_name")
  if (!is.character(formula_sha256) || length(formula_sha256) != 1L ||
      is.na(formula_sha256) || !grepl("^[0-9a-f]{64}$", formula_sha256)) {
    stop("Configured fresh Cox ingress requires one formula digest.",
         call. = FALSE)
  }
  list(analysis_id = analysis_id, data_name = data_name,
       formula_sha256 = formula_sha256)
}

.dsvert_formal_cox_fresh_source_payload <- function(action, payload) {
  if (!is.list(payload) || (!is.null(attributes(payload)) &&
      length(setdiff(names(attributes(payload)), "names")))) {
    stop("Configured fresh Cox ingress requires a protocol payload.",
         call. = FALSE)
  }
  fields <- names(payload)
  if (is.null(fields)) fields <- character()
  expected <- switch(action,
    shape = character(), ticket = character(),
    produce = c("recipient_tickets", "block_index"),
    delivery = c("recipient_tickets", "block_index", "recipient_peer_name"),
    import = c("recipient_tickets", "delivery"),
    provision = c("recipient_tickets", "delivery"))
  if (anyNA(fields) || anyDuplicated(fields) || !identical(fields, expected)) {
    stop("Configured fresh Cox ingress requires a closed action payload.",
         call. = FALSE)
  }
  if ("recipient_tickets" %in% fields &&
      (!is.list(payload$recipient_tickets) ||
       length(payload$recipient_tickets) != 2L ||
       !is.null(names(payload$recipient_tickets)))) {
    stop("Configured fresh Cox ingress requires exactly two recipient tickets.",
         call. = FALSE)
  }
  if ("block_index" %in% fields &&
      (!is.numeric(payload$block_index) || length(payload$block_index) != 1L ||
       is.na(payload$block_index) || !is.finite(payload$block_index) ||
       payload$block_index < 0 || payload$block_index != floor(payload$block_index))) {
    stop("Configured fresh Cox ingress requires one non-negative block index.",
         call. = FALSE)
  }
  if ("recipient_peer_name" %in% fields) {
    payload$recipient_peer_name <- .dsvert_formal_cox_frontdoor_id(
      payload$recipient_peer_name, "recipient peer")
  }
  if ("delivery" %in% fields) {
    encoded <- tryCatch(.dsvert_joint_dp_client_json(payload$delivery),
                        error = function(error) NULL)
    if (!is.list(payload$delivery) || !is.character(encoded) ||
        length(encoded) != 1L || is.na(encoded) ||
        nchar(encoded, type = "bytes") < 2L ||
        nchar(encoded, type = "bytes") > .DSVERT_FORMAL_COX_FRESH_SOURCE_MAX_BYTES) {
      stop("Configured fresh Cox ingress requires one bounded delivery.",
           call. = FALSE)
    }
  }
  payload
}

.dsvert_formal_cox_fresh_source_reply <- function(value, action) {
  fields <- c("version", "action", "payload", "production_ready")
  encoded <- tryCatch(.dsvert_joint_dp_client_json(value$payload),
                      error = function(error) NULL)
  forbidden <- paste0(
    "\\\"[^\\\"]*(?:recipient_signing_key|source_signing_key|",
    "canonical_input_base64|rows|private_key|storage_root|path|secret|config)",
    "[^\\\"]*\\\"\\s*:")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version, .DSVERT_FORMAL_COX_FRESH_SOURCE_REPLY_VERSION) &&
    identical(value$action, action) && is.list(value$payload) &&
    identical(value$production_ready, FALSE) && is.character(encoded) &&
    length(encoded) == 1L && !is.na(encoded) &&
    nchar(encoded, type = "bytes") <= .DSVERT_FORMAL_COX_FRESH_SOURCE_MAX_BYTES &&
    !grepl(forbidden, encoded, perl = TRUE, ignore.case = TRUE)
  if (!isTRUE(valid)) {
    stop("A server returned an unsafe configured fresh Cox source reply.",
         call. = FALSE)
  }
  value
}

.dsvert_formal_cox_fresh_source_call <- function(
    conn, selector, action, payload,
    .aggregate = DSI::datashield.aggregate) {
  if (!is.list(conn) || length(conn) != 1L || is.null(names(conn)) ||
      !is.character(names(conn)) || length(names(conn)) != 1L ||
      is.na(names(conn)) || !nzchar(names(conn))) {
    stop("Configured fresh Cox ingress requires one named server.",
         call. = FALSE)
  }
  if (!is.list(selector) || !identical(
      names(selector), c("analysis_id", "data_name", "formula_sha256"))) {
    stop("Configured fresh Cox ingress requires one fixed source selector.",
         call. = FALSE)
  }
  selector <- .dsvert_formal_cox_fresh_source_selector(
    selector$analysis_id, selector$data_name, selector$formula_sha256)
  actions <- c("shape", "ticket", "produce", "delivery", "import", "provision")
  if (!is.character(action) || length(action) != 1L || is.na(action) ||
      !action %in% actions) {
    stop("Configured fresh Cox ingress requires one closed action.",
         call. = FALSE)
  }
  payload <- .dsvert_formal_cox_fresh_source_payload(action, payload)
  replies <- .dsvert_aggregate_strict(
    conns = conn,
    expr = call(
      name = "dsvertFormalCoxFreshSourceDS",
      analysis_id = selector$analysis_id, data_name = selector$data_name,
      formula_sha256 = selector$formula_sha256,
      action = action, payload = payload),
    operation = "configured fresh formal Cox source relay",
    .aggregate = .aggregate)
  .dsvert_formal_cox_fresh_source_reply(replies[[1L]], action)
}
