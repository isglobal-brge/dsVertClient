test_that("exact GC setup selects exactly two compute peers from K >= 3", {
  federation <- c("site_a", "site_b", "site_c")
  conns <- stats::setNames(lapply(federation, function(...) {
    structure(list(), class = "mock")
  }), federation)
  calls <- list()
  bind_maps <- list()
  encode <- function(byte, n) {
    value <- gsub("[\r\n]", "", jsonlite::base64_enc(as.raw(rep(byte, n))))
    gsub("=+$", "", gsub("/", "_", gsub("+", "-", value, fixed = TRUE),
                            fixed = TRUE))
  }
  aggregate <- function(conns, expr, ...) {
    calls[[length(calls) + 1L]] <<- names(conns)
    command <- as.character(expr[[1L]])
    if (identical(command, "exactGCTransportInitDS")) {
      return(stats::setNames(lapply(seq_along(conns), function(index) list(
        capability_id = "exact_gc_v1",
        transport_pk = encode(index, 32L),
        identity_pk = encode(index + 3L, 32L),
        signature = encode(index + 6L, 64L))), names(conns)))
    }
    if (identical(command, "exactGCBindPeersDS")) {
      decode_map <- function(value) {
        padded <- gsub("-", "+", gsub("_", "/", value, fixed = TRUE),
                       fixed = TRUE)
        remainder <- nchar(padded) %% 4L
        if (remainder == 2L) padded <- paste0(padded, "==")
        if (remainder == 3L) padded <- paste0(padded, "=")
        jsonlite::fromJSON(
          rawToChar(jsonlite::base64_dec(padded)), simplifyVector = FALSE)
      }
      arguments <- as.list(expr)
      bind_maps[[length(bind_maps) + 1L]] <<- list(
        transport = decode_map(arguments[["transport_keys_b64"]]),
        identity = decode_map(arguments[["identity_info_b64"]]))
      return(stats::setNames(lapply(conns, function(...) list(
        capability_id = "exact_gc_v1", bound = TRUE)), names(conns)))
    }
    stop("unexpected endpoint")
  }

  selected <- .dsvert_setup_exact_gc_transport(
    conns, federation, c(1L, 3L),
    "12345678-1234-4234-9234-123456789abc", .aggregate = aggregate)
  expect_identical(names(selected), c("site_a", "site_c"))
  expect_identical(calls, rep(list(c("site_a", "site_c")), 2L))
  expect_length(bind_maps, 1L)
  expect_identical(sort(names(bind_maps[[1L]]$transport)),
                   c("site_a", "site_c"))
  expect_identical(sort(names(bind_maps[[1L]]$identity)),
                   c("site_a", "site_c"))
  expect_false("site_b" %in% names(bind_maps[[1L]]$transport))

  before <- length(calls)
  expect_error(.dsvert_setup_exact_gc_transport(
    conns, federation, 1:3,
    "12345678-1234-4234-9234-123456789abc", .aggregate = aggregate),
    "exactly two")
  expect_identical(length(calls), before)
})
