#' @title Adaptive Chunking for DataSHIELD Transport
#' @description Utilities for sending large payloads through DataSHIELD's
#'   R expression parser with automatic chunk size detection and fallback.
#'
#' @details
#' DataSHIELD passes function arguments as inline R expressions via HTTP.
#' Large string arguments (base64-encoded ciphertexts, EC points, encrypted
#' vectors) can exceed the parser or HTTP body size limit. This module
#' implements adaptive chunking:
#'
#' \enumerate{
#'   \item Start with a large initial chunk size (default 200KB, configurable
#'     via \code{options(dsvert.chunk_size = N)})
#'   \item On first large send, probe whether the chunk size is accepted
#'   \item If the probe fails, reduce the chunk size by 25\% and retry
#'   \item Cache the successful chunk size for all subsequent sends
#' }
#'
#' The probe is only conclusive when a chunk of at least half the target
#' chunk size is successfully sent. Small payloads that succeed trivially
#' do not confirm that the full chunk size is accepted.
#'
#' @name chunk-utils
NULL

# Package-level environment for caching the effective chunk size across calls.
# Persists for the R session (reset on package reload or explicit reset).
.dsvert_chunk_env <- new.env(parent = emptyenv())
.dsvert_chunk_env$effective_chunk_size <- NULL
.dsvert_chunk_env$probed <- FALSE

#' Get the current effective chunk size (internal)
#'
#' Returns the cached effective chunk size if available (from a previous
#' successful probe), otherwise returns the initial chunk size from
#' \code{getOption("dsvert.chunk_size")} or the default (200KB).
#'
#' @return Integer. Chunk size in characters.
#' @keywords internal
.dsvert_get_chunk_size <- function() {
  if (!is.null(.dsvert_chunk_env$effective_chunk_size)) {
    return(.dsvert_chunk_env$effective_chunk_size)
  }
  as.integer(getOption("dsvert.chunk_size", 200000L))
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
#' Forces re-probing on the next send. Useful if connecting to a different
#' Opal instance with different limits.
#'
#' @keywords internal
.dsvert_reset_chunk_size <- function() {
  .dsvert_chunk_env$effective_chunk_size <- NULL
  .dsvert_chunk_env$probed <- FALSE
}

#' Send chunked data via DataSHIELD with adaptive chunk sizing (internal)
#'
#' Wraps the standard chunk-and-send pattern with automatic fallback.
#' On the first large send in a session, probes the server with the initial
#' chunk size. If the probe fails (expression too large), reduces the
#' chunk size by 25\% and retries until a working size is found (or the
#' minimum is reached). The successful chunk size is cached for all subsequent
#' sends.
#'
#' Small payloads (less than half the chunk size) are sent directly without
#' counting as a conclusive probe, since they don't test the actual limit.
#'
#' @param data Character. The full payload string to send.
#' @param send_one_chunk Function(chunk_str, chunk_index, n_chunks).
#'   Callback that sends a single chunk via \code{DSI::datashield.aggregate}.
#'   Must throw an error on failure.
#' @param min_chunk_size Integer. Minimum chunk size before giving up.
#'   Default 10000 (10KB).
#'
#' @return Integer. Number of chunks sent (invisible).
#' @keywords internal
.dsvert_adaptive_send <- function(data, send_one_chunk,
                                  min_chunk_size = 10000L) {
  chunk_size <- .dsvert_get_chunk_size()
  n_chars <- nchar(data)

  # Fast path: already probed, chunk size is known-good.
  # Safety net: if the first chunk fails (e.g., different server with stricter
  # limits), reset the probe and fall through to the probe path. This ensures
  # the cached size converges to the minimum across all servers.
  if (.dsvert_chunk_env$probed) {
    n_chunks <- max(1L, ceiling(n_chars / chunk_size))
    first_end <- min(chunk_size, n_chars)
    first_ok <- tryCatch({
      send_one_chunk(substr(data, 1L, first_end), 1L, as.integer(n_chunks))
      TRUE
    }, error = function(e) {
      if (first_end < chunk_size %/% 2L) stop(e)  # small payload — error is not size-related
      FALSE
    })
    if (first_ok) {
      if (n_chunks > 1L) {
        for (ch in 2:n_chunks) {
          start <- (ch - 1L) * chunk_size + 1L
          end <- min(ch * chunk_size, n_chars)
          send_one_chunk(substr(data, start, end), as.integer(ch), as.integer(n_chunks))
        }
      }
      return(invisible(n_chunks))
    }
    # Fast path failed — different server with stricter limits. Re-probe.
    .dsvert_chunk_env$probed <- FALSE
    chunk_size <- (chunk_size * 3L) %/% 4L  # reduce by 25%
    message(sprintf("[dsVert] Cached chunk size too large for this server, reducing to %dKB",
                    chunk_size %/% 1000L))
  }

  # Small payload: send directly. Don't count as conclusive probe since
  # the data is too small to test the actual parser limit.
  if (n_chars < chunk_size %/% 2L) {
    send_one_chunk(data, 1L, 1L)
    return(invisible(1L))
  }

  # Probe path: first large send — discover the working chunk size
  repeat {
    n_chunks <- max(1L, ceiling(n_chars / chunk_size))

    # Try sending the first chunk as a probe
    first_end <- min(chunk_size, n_chars)
    first_chunk <- substr(data, 1L, first_end)

    probe_ok <- tryCatch({
      send_one_chunk(first_chunk, 1L, as.integer(n_chunks))
      TRUE
    }, error = function(e) {
      # If chunk_size is already at minimum, the error is not about size
      if (chunk_size <= min_chunk_size) stop(e)
      FALSE
    })

    if (probe_ok) {
      # Probe succeeded — send remaining chunks
      if (n_chunks > 1L) {
        for (ch in 2:n_chunks) {
          start <- (ch - 1L) * chunk_size + 1L
          end <- min(ch * chunk_size, n_chars)
          send_one_chunk(substr(data, start, end), as.integer(ch), as.integer(n_chunks))
        }
      }
      # Only cache as conclusive if we actually sent a chunk near chunk_size
      if (first_end >= chunk_size %/% 2L) {
        .dsvert_set_chunk_size(chunk_size)
      }
      return(invisible(n_chunks))
    }

    # Probe failed — reduce chunk size by 25% and retry
    chunk_size <- (chunk_size * 3L) %/% 4L
    if (chunk_size < min_chunk_size) {
      stop("DataSHIELD chunk transfer failed even at ", min_chunk_size,
           " byte chunks. Check Opal/Rock server configuration.",
           call. = FALSE)
    }
    message(sprintf("[dsVert] Expression size limit detected, reducing chunk size to %dKB",
                    chunk_size %/% 1000L))
  }
}
