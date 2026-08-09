local_chunk_size <- function(size) {
  old <- options(dsvert.chunk_size = size)
  dsVertClient:::.dsvert_reset_chunk_size()
  withr::defer({
    options(old)
    dsVertClient:::.dsvert_reset_chunk_size()
  }, envir = parent.frame())
}

legacy_peer_rejection_token <- function(peer = "site_a") {
  paste0(
    "DSVERT_PEER_NOT_RECOGNIZED_V1|peer=", peer,
    "|expected=", strrep("1", 64L),
    "|observed=", strrep("2", 64L),
    " SECRET server trace")
}

test_that("blob transport requires typed TRUE from the exact target", {
  local_chunk_size(100L)
  conn <- list(site_a = structure(list(), class = "fake_conn"))

  malformed <- list(
    setNames(list(1L), "site_a"),
    setNames(list(TRUE), "site_b"),
    setNames(list(TRUE, NULL), c("site_a", "site_b")))

  for (reply in malformed) {
    calls <- 0L
    aggregate <- function(conns, expr, ...) {
      calls <<- calls + 1L
      reply
    }
    expect_error(
      dsVertClient:::.dsvert_store_blob(
        "opaque", "k2_peer_x_share", conn, "session-a",
        .aggregate = aggregate),
      "malformed acknowledgement")
    expect_identical(calls, 1L)
  }
})

test_that("NULL acknowledgement is never treated as a successful chunk", {
  local_chunk_size(4L)
  withr::local_options(list(dsvert.dsi.retry_deadline_seconds = 0.001))
  conn <- list(site_a = structure(list(), class = "fake_conn"))
  calls <- list()
  aggregate <- function(conns, expr, ...) {
    calls[[length(calls) + 1L]] <<- expr
    setNames(list(NULL), "site_a")
  }

  condition <- expect_error(
    dsVertClient:::.dsvert_store_blob(
      "abcdefgh", "k2_peer_x_share", conn, "session-a",
        .aggregate = aggregate),
    class = "retry_deadline_exceeded")

  expect_s3_class(condition, "transport_unavailable")
  expect_gte(length(calls), 1L)
  expect_true(all(vapply(calls, identical, logical(1L), calls[[1L]])))
  expect_identical(calls[[1L]][["chunk_index"]], 1L)
  expect_identical(calls[[1L]][["n_chunks"]], 2L)
})

test_that("lost response replays identical bytes and immutable geometry", {
  local_chunk_size(4L)
  conn <- list(site_a = structure(list(), class = "fake_conn"))
  calls <- list()
  lost_once <- FALSE
  stored <- list()

  aggregate <- function(conns, expr, ...) {
    calls[[length(calls) + 1L]] <<- expr
    index <- expr[["chunk_index"]]
    chunk <- expr[["chunk"]]
    prior <- stored[[as.character(index)]]
    if (!is.null(prior)) expect_identical(chunk, prior)
    stored[[as.character(index)]] <<- chunk

    if (identical(index, 2L) && !lost_once) {
      lost_once <<- TRUE
      return(setNames(list(NULL), "site_a"))
    }
    setNames(list(TRUE), "site_a")
  }

  sent <- dsVertClient:::.dsvert_store_blob(
    "abcdefghij", "k2_peer_x_share", conn, "session-a",
    .aggregate = aggregate)

  expect_identical(sent, 3L)
  indices <- vapply(calls, function(x) x[["chunk_index"]], integer(1L))
  expect_identical(indices, c(1L, 2L, 2L, 3L))
  expect_identical(calls[[2L]], calls[[3L]])
  expect_true(all(vapply(calls, function(x)
    identical(x[["n_chunks"]], 3L), logical(1L))))
  expect_identical(paste0(unlist(stored, use.names = FALSE), collapse = ""),
                   "abcdefghij")
})

test_that("automatic replay is disabled for callbacks not declared idempotent", {
  local_chunk_size(100L)
  calls <- 0L
  expect_error(
    dsVertClient:::.dsvert_adaptive_send(
      "opaque",
      function(chunk, index, total) {
        calls <<- calls + 1L
        stop("connection lost")
      },
      target = "site_a", idempotent = FALSE),
    "not acknowledged.*no partial result")
  expect_identical(calls, 1L)
})

test_that("legacy store callbacks are terminal even on its idempotent path", {
  local_chunk_size(100L)
  conn <- list(site_a = structure(list(), class = "fake_conn"))
  calls <- 0L
  aggregate <- function(conns, expr, async, error, errors.print) {
    calls <<- calls + 1L
    error("site_a", "conflicting server mutation")
    list(site_a = NULL)
  }
  expect_error(dsVertClient:::.dsvert_store_blob(
    "opaque", "k2_peer_x_share", conn, "session-a",
    .aggregate = aggregate), "malformed acknowledgement")
  expect_identical(calls, 1L)
})

test_that("legacy store propagates callback peer rejection without retry", {
  local_chunk_size(100L)
  withr::local_options(list(dsvert.dsi.retry_deadline_seconds = 0.001))
  conn <- list(site_a = structure(list(), class = "fake_conn"))
  calls <- 0L
  aggregate <- function(conns, expr, async, error, errors.print) {
    calls <<- calls + 1L
    error("site_a", legacy_peer_rejection_token())
    list(site_a = NULL)
  }

  condition <- tryCatch(
    dsVertClient:::.dsvert_store_blob(
      "opaque", "k2_peer_x_share", conn, "session-a",
      .aggregate = aggregate),
    dsvert_peer_not_recognized = identity)

  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(calls, 1L)
  expect_identical(condition$peer_name, "site_a")
  expect_false(grepl("SECRET", conditionMessage(condition), fixed = TRUE))
})

test_that("legacy store propagates top-level peer rejection without retry", {
  local_chunk_size(100L)
  withr::local_options(list(dsvert.dsi.retry_deadline_seconds = 0.001))
  conn <- list(site_a = structure(list(), class = "fake_conn"))
  calls <- 0L
  aggregate <- function(...) {
    calls <<- calls + 1L
    stop(legacy_peer_rejection_token(), call. = FALSE)
  }

  condition <- tryCatch(
    dsVertClient:::.dsvert_store_blob(
      "opaque", "k2_peer_x_share", conn, "session-a",
      .aggregate = aggregate),
    dsvert_peer_not_recognized = identity)

  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(calls, 1L)
  expect_identical(condition$peer_name, "site_a")
  expect_false(grepl("SECRET", conditionMessage(condition), fixed = TRUE))
})

test_that("legacy store propagates terminal oversize without retry", {
  local_chunk_size(100L)
  conn <- list(site_a = structure(list(), class = "fake_conn"))
  tag <- paste0(
    "[dsvert_resource_oversize:v1] resource_oversize: ",
    "legacy payload cannot fit fixed policy")

  for (mode in c("callback", "top")) {
    calls <- 0L
    aggregate <- function(conns, expr, async, error, errors.print) {
      calls <<- calls + 1L
      if (identical(mode, "callback")) {
        error("site_a", tag)
        return(list(site_a = NULL))
      }
      stop(tag, call. = FALSE)
    }
    condition <- tryCatch(dsVertClient:::.dsvert_store_blob(
      "opaque", "k2_peer_x_share", conn, "session-a",
      .aggregate = aggregate), error = identity)
    expect_s3_class(condition, "dsvert_resource_oversize")
    expect_identical(condition$code, "resource_oversize")
    expect_false(condition$retryable)
    expect_identical(calls, 1L)
  }
})

test_that("chunk transport rejects non-ASCII payloads and unrepresentable geometry", {
  local_chunk_size(100L)
  calls <- 0L
  expect_error(
    dsVertClient:::.dsvert_adaptive_send(
      "caf\u00e9",
      function(...) {
        calls <<- calls + 1L
        list(site_a = TRUE)
      },
      target = "site_a", idempotent = TRUE),
    "ASCII/Base64")
  expect_identical(calls, 0L)

  expect_error(dsVertClient:::.dsvert_chunk_count(
    .Machine$integer.max + 1, 1L), "more chunks")
})

test_that("legacy Base64 frame bound is conservative without full deparse", {
  chunk <- strrep("A", 640 * 1024L)
  key <- "k2_peer_x_share"
  sid <- "12345678-1234-4234-8234-123456789abc"
  expr <- call(
    name = "mpcStoreBlobDS", key = key, chunk = chunk,
    chunk_index = 17L, n_chunks = 101L, session_id = sid)
  actual <- nchar(
    paste(deparse(expr, width.cutoff = 500L), collapse = "\n"),
    type = "bytes")
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  withr::local_options(list(dsvert.dsi_max_expression_bytes = actual + 4096L))
  bound <- .dsvert_validate_legacy_blob_frame_size(
    expr, "mpcStoreBlobDS", key, chunk, 17L, 101L, sid)
  expect_gte(bound, actual)
  expect_lte(bound - actual, .DSVERT_LEGACY_BLOB_FRAME_OVERHEAD_BYTES)

  # Compatibility strings outside the fast Base64 schema still use the exact
  # generic expression guard rather than acquiring a new rejection rule.
  unusual <- "still opaque: {legacy payload}"
  unusual_expr <- call(
    name = "mpcStoreBlobDS", key = key, chunk = unusual,
    chunk_index = 1L, n_chunks = 1L, session_id = sid)
  expect_identical(
    unname(.dsvert_validate_legacy_blob_frame_size(
      unusual_expr, "mpcStoreBlobDS", key, unusual, 1L, 1L, sid)),
    unname(.dsvert_validate_dsi_expression_sizes(unusual_expr)))
})
