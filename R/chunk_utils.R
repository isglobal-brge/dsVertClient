#' @title Acknowledged Chunking for DataSHIELD Transport
#' @description Utilities for sending large payloads through DataSHIELD's
#'   R expression parser with fixed transfer geometry and typed acknowledgements.
#'
#' @details
#' DataSHIELD passes function arguments as inline R expressions via HTTP.
#' Large string arguments (base64-encoded ciphertexts, EC points, encrypted
#' vectors) can exceed the parser or HTTP body size limit. This module
#' implements acknowledged chunking:
#'
#' \enumerate{
#'   \item Use a conservative unprobed default of 640 KiB of Base64 text, then
#'     let production DSI negotiate 688 KiB (enough for 512 KiB raw) or a
#'     smaller common value before any payload is sent. Test harnesses may set
#'     \code{options(dsvert.chunk_size = N)} explicitly.
#'   \item Freeze the chunk size, chunk count, indices, and bytes before the
#'     first request
#'   \item Require an exact logical \code{TRUE} acknowledgement from the
#'     intended DataSHIELD target before advancing
#'   \item If an idempotent store loses its response, replay only the exact
#'     same request until a monotonic availability deadline; never silently
#'     change transfer geometry
#' }
#'
#' Automatic fallback to a smaller geometry is deliberately not attempted
#' after transmission starts: the server may already have committed the first
#' chunk even when its response is lost. Changing \code{n_chunks} at that point
#' would turn a recoverable response loss into a conflicting transfer.
#'
#' @name chunk-utils
NULL

# Package-level environment for caching the effective chunk size across calls.
# Persists for the R session (reset on package reload or explicit reset).
.dsvert_chunk_env <- new.env(parent = emptyenv())
.dsvert_chunk_env$effective_chunk_size <- NULL
.dsvert_chunk_env$effective_expression_bytes <- NULL
.dsvert_chunk_env$probed <- FALSE
.dsvert_chunk_env$geometry_locked <- FALSE
.dsvert_chunk_env$active_probe_cache_key <- NULL
.dsvert_chunk_env$probe_cache <- list()
.dsvert_chunk_env$probe_hint_cache <- list()
.dsvert_chunk_env$probe_clock <- 0

.DSVERT_LEGACY_BLOB_FRAME_OVERHEAD_BYTES <- 1024L

# The normal legacy payload is Base64/Base64url, and its key/session fields use
# the server's storage alphabet. For that fixed schema, adding scalar byte
# lengths to a conservative syntax reserve avoids materialising a second full
# deparse of every frame. Unusual compatibility input falls back to the generic
# exact deparse guard instead of being rejected solely by this optimization.
.dsvert_validate_legacy_blob_frame_size <- function(
    expression, store_function, key, chunk, chunk_index, n_chunks,
    session_id) {
  safe_ascii <- function(value, pattern, max_bytes) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value) && nchar(value, type = "bytes") <= max_bytes &&
      identical(nchar(value, type = "chars"),
                nchar(value, type = "bytes")) &&
      grepl(pattern, value, perl = TRUE, useBytes = TRUE)
  }
  integer_scalar <- function(value, lower, upper) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == floor(value) &&
      value >= lower && value <= upper
  }
  fixed_schema <- is.call(expression) &&
    safe_ascii(store_function, "^mpcStoreBlobDS$", 32L) &&
    safe_ascii(key, "^[A-Za-z0-9._-]+$", 255L) &&
    safe_ascii(chunk, "^[A-Za-z0-9+/_=-]+$",
               .DSVERT_DSI_PROBE_ABSOLUTE_MAX) &&
    safe_ascii(session_id, "^[A-Za-z0-9._-]+$", 128L) &&
    integer_scalar(n_chunks, 1L, 4096L) &&
    integer_scalar(chunk_index, 1L, n_chunks)
  if (!isTRUE(fixed_schema)) {
    return(.dsvert_validate_dsi_expression_sizes(expression))
  }
  numbers <- c(
    format(chunk_index, scientific = FALSE, trim = TRUE),
    format(n_chunks, scientific = FALSE, trim = TRUE))
  size_bound <- sum(nchar(
    c(store_function, key, chunk, session_id, numbers), type = "bytes")) +
    .DSVERT_LEGACY_BLOB_FRAME_OVERHEAD_BYTES
  capacity <- .dsvert_dsi_max_expression_bytes()
  if (size_bound > capacity) {
    stop(.dsvert_client_resource_oversize(
      requested_bytes = size_bound, capacity_bytes = capacity,
      scope = "legacy blob DataSHIELD expression"))
  }
  invisible(as.numeric(size_bound))
}

#' Get the current effective chunk size (internal)
#'
#' Returns the cached effective chunk size if available (from a previous
#' successful transfer), otherwise returns the configured chunk size from
#' \code{getOption("dsvert.chunk_size")} or the conservative unprobed default
#' (640 KiB of Base64 text, representing at most 480 KiB raw). Production DSI
#' negotiation can replace it with a verified 688 KiB text geometry, enough
#' for a 512 KiB raw frame.
#'
#' @return Integer. Chunk size in characters.
#' @keywords internal
.dsvert_get_chunk_size <- function() {
  size <- if (!is.null(.dsvert_chunk_env$effective_chunk_size)) {
    .dsvert_chunk_env$effective_chunk_size
  } else {
    getOption("dsvert.chunk_size", 640L * 1024L)
  }
  if (!is.numeric(size) || length(size) != 1L || is.na(size) ||
      !is.finite(size) || size != floor(size) || size < 1 ||
      size > .DSVERT_DSI_PROBE_ABSOLUTE_MAX) {
    stop("option 'dsvert.chunk_size' must be one positive integer no larger than 8 MiB",
         call. = FALSE)
  }
  as.integer(size)
}

#' Set the effective chunk size (internal)
#' @param size Integer. Chunk size in characters.
#' @keywords internal
.dsvert_set_chunk_size <- function(size) {
  .dsvert_chunk_env$effective_chunk_size <- as.integer(size)
  .dsvert_chunk_env$probed <- TRUE
}

#' Reset the cached chunk size (internal)
#'
#' Forces the configured chunk size to be read again on the next send.
#'
#' @keywords internal
.dsvert_reset_chunk_size <- function() {
  .dsvert_chunk_env$effective_chunk_size <- NULL
  .dsvert_chunk_env$effective_expression_bytes <- NULL
  .dsvert_chunk_env$probed <- FALSE
  .dsvert_chunk_env$geometry_locked <- FALSE
  .dsvert_chunk_env$active_probe_cache_key <- NULL
}

#' Classify a one-target DataSHIELD acknowledgement (internal)
#'
#' @param response Result returned by \code{DSI::datashield.aggregate}.
#' @param target Expected logical DataSHIELD connection name.
#' @return One of \code{"ack"}, \code{"missing"}, or \code{"invalid"}.
#' @keywords internal
.dsvert_chunk_ack_state <- function(response, target) {
  if (is.null(response)) return("missing")
  if (!is.list(response) || length(response) != 1L ||
      is.null(names(response)) || !identical(names(response), target)) {
    return("invalid")
  }
  ack <- response[[1L]]
  if (identical(ack, TRUE)) return("ack")
  if (is.null(ack)) return("missing")
  "invalid"
}

#' Compute a representable immutable chunk count (internal)
#' @keywords internal
.dsvert_chunk_count <- function(n_chars, chunk_size) {
  count <- max(1, ceiling(n_chars / chunk_size))
  if (!is.finite(count) || count > .Machine$integer.max) {
    stop("payload requires more chunks than the protocol can represent",
         call. = FALSE)
  }
  as.integer(count)
}

#' Send chunked data via DataSHIELD with immutable geometry (internal)
#'
#' The complete transfer geometry is fixed before the first request. A
#' callback exception or missing acknowledgement may be replayed
#' when \code{idempotent = TRUE}; every replay receives byte-for-byte identical
#' chunk data and identical indices until an availability deadline. There is
#' no request-attempt quota. Explicit negative or malformed replies are
#' protocol failures and are never retried.
#'
#' @param data Character. The full payload string to send.
#' @param send_one_chunk Function(chunk_str, chunk_index, n_chunks).
#'   Callback that sends one chunk and returns the unmodified result from
#'   \code{DSI::datashield.aggregate}.
#' @param target Expected logical DataSHIELD connection name.
#' @param idempotent Whether exact replay is permitted after an ambiguous
#'   response. This must only be enabled for a server operation whose duplicate
#'   request is explicitly idempotent.
#'
#' @return Integer. Number of chunks sent (invisible).
#' @keywords internal
.dsvert_adaptive_send <- function(data, send_one_chunk,
                                  target, idempotent = FALSE) {
  if (!is.character(data) || length(data) != 1L || is.na(data) ||
      !nzchar(data)) {
    stop("data must be one non-empty character string", call. = FALSE)
  }
  if (!is.function(send_one_chunk)) {
    stop("send_one_chunk must be a function", call. = FALSE)
  }
  if (!is.character(target) || length(target) != 1L || is.na(target) ||
      !nzchar(target)) {
    stop("target must be one non-empty DataSHIELD connection name",
         call. = FALSE)
  }
  if (!is.logical(idempotent) || length(idempotent) != 1L ||
      is.na(idempotent)) {
    stop("idempotent must be TRUE or FALSE", call. = FALSE)
  }

  chunk_size <- .dsvert_get_chunk_size()
  n_chars <- nchar(data, type = "chars")
  n_bytes <- nchar(data, type = "bytes")
  if (!identical(n_chars, n_bytes)) {
    stop("chunk transport accepts ASCII/Base64 payloads only", call. = FALSE)
  }
  n_chunks <- .dsvert_chunk_count(n_chars, chunk_size)

  # From this point onward geometry is immutable.  A failed or ambiguous first
  # payload frame may be replayed exactly, but transport probing/downshifting is
  # forbidden until the enclosing method resets its transport state.
  .dsvert_chunk_env$geometry_locked <- TRUE
  on.exit({
    .dsvert_chunk_env$geometry_locked <- FALSE
  }, add = TRUE)

  send_locked_chunk <- function(chunk, chunk_index) {
    attempt <- function() {
      tryCatch(
        list(ok = TRUE,
             response = send_one_chunk(chunk, chunk_index, n_chunks)),
        interrupt = function(e) stop(e),
        error = function(e) {
          rejection <- if (inherits(e, "dsvert_peer_not_recognized")) {
            e
          } else {
            .dsvert_client_parse_peer_not_recognized(conditionMessage(e))
          }
          if (inherits(rejection, "dsvert_peer_not_recognized")) {
            stop(rejection)
          }
          oversized <- if (inherits(e, "dsvert_resource_oversize")) {
            e
          } else {
            .dsvert_client_parse_resource_oversize(conditionMessage(e))
          }
          if (inherits(oversized, "dsvert_resource_oversize")) {
            stop(oversized)
          }
          list(ok = FALSE, error = e)
        })
    }
    classify <- function(outcome) list(
      state = if (isTRUE(outcome$ok)) {
        .dsvert_chunk_ack_state(outcome$response, target)
      } else {
        "missing"
      })
    outcome <- if (isTRUE(idempotent)) {
      .dsvert_retry_idempotent(
        attempt = attempt, classify = classify,
        operation = "legacy immutable blob frame commit")
    } else {
      classify(attempt())
    }
    if (identical(outcome$state, "ack")) return(invisible(TRUE))
    if (identical(outcome$state, "invalid")) {
      stop("Chunk store returned a malformed acknowledgement for target '",
           target, "'; expected exactly one named logical TRUE.",
           call. = FALSE)
    }
    stop("Chunk store was not acknowledged by target '", target,
         "'; no partial result was accepted.", call. = FALSE)
  }

  for (ch in seq_len(n_chunks)) {
    start <- (ch - 1L) * chunk_size + 1L
    end <- min(ch * chunk_size, n_chars)
    chunk <- substr(data, start, end)
    send_locked_chunk(chunk, as.integer(ch))
  }

  .dsvert_set_chunk_size(chunk_size)
  invisible(n_chunks)
}

#' Store one acknowledged opaque blob on one DataSHIELD target (internal)
#'
#' This is the only legacy client helper that enables automatic chunk replay.
#' Both supported server methods reject conflicting duplicates and return an
#' exact logical \code{TRUE} after either the initial commit or an identical
#' replay.
#'
#' @keywords internal
.dsvert_store_blob <- function(blob, key, conn, session_id,
                               store_function = "mpcStoreBlobDS",
                               .aggregate = DSI::datashield.aggregate) {
  if (!is.character(store_function) || length(store_function) != 1L ||
      is.na(store_function) ||
      !identical(store_function, "mpcStoreBlobDS")) {
    stop("store_function must name a supported idempotent blob endpoint",
         call. = FALSE)
  }
  if (length(conn) != 1L || is.null(names(conn)) ||
      !is.character(names(conn)) || is.na(names(conn)) ||
      !nzchar(names(conn))) {
    stop("conn must contain exactly one named DataSHIELD target",
         call. = FALSE)
  }
  target <- names(conn)
  .dsvert_maybe_negotiate_dsi_chunk_size(conn, .aggregate)
  .dsvert_validate_real_dsi_transport(conn, .aggregate)
  .dsvert_adaptive_send(
    blob,
    function(chunk_str, chunk_idx, n_chunks) {
      expr <- call(name = store_function, key = key, chunk = chunk_str,
                   chunk_index = chunk_idx, n_chunks = n_chunks,
                   session_id = session_id)
      .dsvert_validate_legacy_blob_frame_size(
        expr, store_function, key, chunk_str, chunk_idx, n_chunks,
        session_id)
      callback_failed <- FALSE
      peer_rejection <- NULL
      resource_oversize <- NULL
      result <- tryCatch(
        .dsvert_transport_aggregate(
          .aggregate = .aggregate,
          conns = conn, expr = expr, async = TRUE,
          error = function(site, message) {
            callback_failed <<- TRUE
            parsed <- .dsvert_client_parse_peer_not_recognized(message)
            if (is.null(peer_rejection) && !is.null(parsed)) {
              peer_rejection <<- parsed
            }
            oversized <- .dsvert_client_parse_resource_oversize(message)
            if (is.null(resource_oversize) && !is.null(oversized)) {
              resource_oversize <<- oversized
            }
            invisible(NULL)
          }, errors.print = FALSE),
        interrupt = function(e) stop(e),
        error = function(e) {
          if (inherits(e, "dsvert_dsi_poisoned_session")) stop(e)
          candidate <- if (inherits(e, "dsvert_peer_not_recognized")) {
            e
          } else {
            .dsvert_client_parse_peer_not_recognized(conditionMessage(e))
          }
          if (!inherits(peer_rejection, "dsvert_peer_not_recognized") &&
              inherits(candidate, "dsvert_peer_not_recognized")) {
            peer_rejection <<- candidate
          }
          oversized <- if (inherits(e, "dsvert_resource_oversize")) {
            e
          } else {
            .dsvert_client_parse_resource_oversize(conditionMessage(e))
          }
          if (is.null(resource_oversize) &&
              inherits(oversized, "dsvert_resource_oversize")) {
            resource_oversize <<- oversized
          }
          NULL
        })
      if (inherits(peer_rejection, "dsvert_peer_not_recognized")) {
        stop(peer_rejection)
      }
      if (inherits(resource_oversize, "dsvert_resource_oversize")) {
        stop(resource_oversize)
      }
      # A per-site callback may represent a deterministic server conflict or
      # contract failure.  Keep it distinct from an absent/top-level transport
      # response so the legacy path cannot retry a conflicting mutation.
      if (isTRUE(callback_failed)) return(FALSE)
      result
    },
    target = target,
    idempotent = TRUE)
}
