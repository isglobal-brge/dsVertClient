# Private client seam for the closed registered formal-GLM source relay.
#
# This is transport plumbing for the fresh binomial/Poisson driver, not a
# model frontdoor.  It relays only signed protocol records and one bounded
# opaque pair chunk at a time.

.DSVERT_FORMAL_GLM_REGISTERED_SOURCE_REPLY_VERSION <-
  "dsvert-formal-glm-registered-phase18-source-response-v1"
.DSVERT_FORMAL_GLM_REGISTERED_FRESH_SOURCE_REPLY_VERSION <-
  "dsvert-formal-glm-registered-fresh-source-response-v1"
.DSVERT_FORMAL_GLM_REGISTERED_SOURCE_MAX_BYTES <- 4L * 1024L * 1024L

.dsvert_formal_glm_registered_source_payload <- function(action, payload) {
  attributes <- attributes(payload)
  if (!is.list(payload) || (!is.null(attributes) &&
      length(setdiff(names(attributes), "names")))) {
    stop("Registered formal GLM source requires a protocol payload.",
         call. = FALSE)
  }
  fields <- names(payload)
  if (is.null(fields)) fields <- character()
  expected <- switch(action,
    ticket = character(), receipt_set = character(), host_provision = character(),
    ticket_set = "recipient_tickets", seal_block = c("recipient_tickets", "block_index"),
    chunk = c("recipient_tickets", "block_index", "offset"),
    import_chunk = c("recipient_tickets", "chunk_receipt", "pair_chunk_base64"),
    local_receipt = "recipient_tickets", receipt_commit = "local_receipt_json",
    binding = "recipient_tickets")
  if (anyNA(fields) || anyDuplicated(fields) || !identical(fields, expected)) {
    stop("Registered formal GLM source requires a closed action payload.",
         call. = FALSE)
  }
  if ("recipient_tickets" %in% fields &&
      (!is.list(payload$recipient_tickets) ||
       length(payload$recipient_tickets) != 2L ||
       !is.null(names(payload$recipient_tickets)))) {
    stop("Registered formal GLM source requires exactly two recipient tickets.",
         call. = FALSE)
  }
  if ("block_index" %in% fields &&
      (!is.numeric(payload$block_index) || length(payload$block_index) != 1L ||
       is.na(payload$block_index) || !is.finite(payload$block_index) ||
       payload$block_index < 0 || payload$block_index != floor(payload$block_index))) {
    stop("Registered formal GLM source requires one non-negative block index.",
         call. = FALSE)
  }
  if ("offset" %in% fields &&
      (!is.numeric(payload$offset) || length(payload$offset) != 1L ||
       is.na(payload$offset) || !is.finite(payload$offset) ||
       payload$offset < 0 || payload$offset != floor(payload$offset))) {
    stop("Registered formal GLM source requires one non-negative chunk offset.",
         call. = FALSE)
  }
  payload
}

.dsvert_formal_glm_registered_source_reply <- function(value, action) {
  fields <- c("version", "action", "payload", "production_ready")
  valid <- is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), fields) &&
    identical(value$version, .DSVERT_FORMAL_GLM_REGISTERED_SOURCE_REPLY_VERSION) &&
    identical(value$action, action) && is.list(value$payload) &&
    identical(value$production_ready, FALSE)
  if (!isTRUE(valid)) {
    stop("A server returned an invalid registered formal GLM source reply.",
         call. = FALSE)
  }
  encoded <- tryCatch(.dsvert_joint_dp_client_json(value$payload),
                      error = function(error) NULL)
  forbidden <- paste0(
    "\\\"[^\\\"]*(?:rows|values|validity|private_consensus|alignment_consensus|",
    "authorization_json|local_signing_key|pair_json|storage|path|secret|",
    "backend|context|token|capsule|run_id|opening)[^\\\"]*\\\"\\s*:")
  if (!is.character(encoded) || length(encoded) != 1L || is.na(encoded) ||
      nchar(encoded, type = "bytes") > .DSVERT_FORMAL_GLM_REGISTERED_SOURCE_MAX_BYTES ||
      grepl(forbidden, encoded, perl = TRUE, ignore.case = TRUE)) {
    stop("A server returned an unsafe registered formal GLM source reply.",
         call. = FALSE)
  }
  if (identical(action, "chunk")) {
    chunk_fields <- c("chunk_receipt", "pair_chunk_base64", "replayed")
    chunk <- tryCatch(jsonlite::base64_dec(value$payload$pair_chunk_base64),
                      error = function(error) raw())
    canonical <- gsub("[\\r\\n]", "", jsonlite::base64_enc(chunk))
    if (!identical(names(value$payload), chunk_fields) ||
        !is.list(value$payload$chunk_receipt) ||
        !is.character(value$payload$pair_chunk_base64) ||
        length(value$payload$pair_chunk_base64) != 1L ||
        !identical(value$payload$pair_chunk_base64, canonical) ||
        !length(chunk) || length(chunk) > 1024L^2 ||
        !is.logical(value$payload$replayed) ||
        length(value$payload$replayed) != 1L || is.na(value$payload$replayed)) {
      if (length(chunk)) chunk[] <- as.raw(0L)
      stop("A server returned an invalid registered formal GLM source chunk.",
           call. = FALSE)
    }
    if (length(chunk)) chunk[] <- as.raw(0L)
  }
  value
}

.dsvert_formal_glm_registered_source_call <- function(
    conn, source_contract_json, action, payload,
    .aggregate = DSI::datashield.aggregate) {
  if (!is.list(conn) || length(conn) != 1L || is.null(names(conn)) ||
      !is.character(names(conn)) || length(names(conn)) != 1L ||
      is.na(names(conn)) || !nzchar(names(conn))) {
    stop("Registered formal GLM source requires one named server.",
         call. = FALSE)
  }
  if (!is.character(source_contract_json) || length(source_contract_json) != 1L ||
      is.na(source_contract_json) || nchar(source_contract_json, type = "bytes") < 2L ||
      nchar(source_contract_json, type = "bytes") > 1024L * 1024L) {
    stop("Registered formal GLM source requires one bounded contract.",
         call. = FALSE)
  }
  actions <- c("ticket", "ticket_set", "seal_block", "chunk", "import_chunk",
               "local_receipt", "receipt_commit", "receipt_set", "binding",
               "host_provision")
  if (!is.character(action) || length(action) != 1L || is.na(action) ||
      !action %in% actions) {
    stop("Registered formal GLM source requires one closed action.",
         call. = FALSE)
  }
  payload <- .dsvert_formal_glm_registered_source_payload(action, payload)
  replies <- .dsvert_aggregate_strict(
    conns = conn,
    expr = call(
      name = "dsvertFormalGLMRegisteredSourceDS",
      source_contract_json = source_contract_json, action = action,
      payload = payload),
    operation = "registered formal GLM source relay",
    .aggregate = .aggregate)
  .dsvert_formal_glm_registered_source_reply(replies[[1L]], action)
}

.dsvert_formal_glm_registered_fresh_source_selector <- function(
    analysis_id, data_name, family, formula_sha256) {
  valid_label <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)
  }
  if (!valid_label(analysis_id) || !valid_label(data_name) ||
      !is.character(family) || length(family) != 1L || is.na(family) ||
      !family %in% c("binomial", "poisson") ||
      !is.character(formula_sha256) || length(formula_sha256) != 1L ||
      is.na(formula_sha256) || !grepl("^[0-9a-f]{64}$", formula_sha256)) {
    stop("Registered fresh GLM source requires one fixed analysis selector.",
         call. = FALSE)
  }
  list(
    analysis_id = enc2utf8(analysis_id), data_name = enc2utf8(data_name),
    family = family, formula_sha256 = formula_sha256)
}

.dsvert_formal_glm_registered_fresh_source_call <- function(
    conn, selector, action, payload,
    .aggregate = DSI::datashield.aggregate) {
  if (!is.list(conn) || length(conn) != 1L || is.null(names(conn)) ||
      !is.character(names(conn)) || length(names(conn)) != 1L ||
      is.na(names(conn)) || !nzchar(names(conn))) {
    stop("Registered fresh GLM source requires one named server.",
         call. = FALSE)
  }
  if (!is.list(selector) || !identical(
      names(selector), c("analysis_id", "data_name", "family", "formula_sha256"))) {
    stop("Registered fresh GLM source requires one fixed analysis selector.",
         call. = FALSE)
  }
  selector <- .dsvert_formal_glm_registered_fresh_source_selector(
    selector$analysis_id, selector$data_name, selector$family,
    selector$formula_sha256)
  actions <- c("ticket", "ticket_set", "seal_block", "chunk", "import_chunk",
               "local_receipt", "receipt_commit", "receipt_set", "binding",
               "host_provision")
  if (!is.character(action) || length(action) != 1L || is.na(action) ||
      !action %in% actions) {
    stop("Registered fresh GLM source requires one closed action.",
         call. = FALSE)
  }
  payload <- .dsvert_formal_glm_registered_source_payload(action, payload)
  replies <- .dsvert_aggregate_strict(
    conns = conn,
    expr = call(
      name = "dsvertFormalGLMRegisteredFreshSourceDS",
      analysis_id = selector$analysis_id, data_name = selector$data_name,
      family = selector$family, formula_sha256 = selector$formula_sha256,
      action = action, payload = payload),
    operation = "registered fresh formal GLM source relay",
    .aggregate = .aggregate)
  value <- replies[[1L]]
  if (!is.list(value) || !identical(names(value),
      c("version", "action", "payload", "production_ready")) ||
      !identical(value$version,
                 .DSVERT_FORMAL_GLM_REGISTERED_FRESH_SOURCE_REPLY_VERSION)) {
    stop("A server returned an invalid registered fresh GLM source reply.",
         call. = FALSE)
  }
  value$version <- .DSVERT_FORMAL_GLM_REGISTERED_SOURCE_REPLY_VERSION
  .dsvert_formal_glm_registered_source_reply(value, action)
}
