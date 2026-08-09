test_that("the public peer error parser accepts only the fixed safe contract", {
  expected <- paste(rep("a", 64L), collapse = "")
  observed <- paste(rep("b", 64L), collapse = "")
  token <- paste0(
    "remote prefix DSVERT_PEER_NOT_RECOGNIZED_V1|peer=site_b",
    "|expected=", expected, "|observed=", observed,
    " server-local detail that must not be copied")

  condition <- .dsvert_client_parse_peer_not_recognized(token)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$code, "peer_not_recognized")
  expect_identical(condition$peer_name, "site_b")
  expect_identical(condition$expected_fingerprint_sha256, expected)
  expect_identical(condition$observed_fingerprint_sha256, observed)
  expect_match(conditionMessage(condition), "Server administrator")
  expect_match(conditionMessage(condition), "verify the observed fingerprint")
  expect_match(conditionMessage(condition), "restart", ignore.case = TRUE)
  expect_match(conditionMessage(condition), "each other participating server")
  expect_match(conditionMessage(condition), "must not pin its own identity")
  expect_false(grepl("server-local detail", conditionMessage(condition),
                     fixed = TRUE))

  expect_null(.dsvert_client_parse_peer_not_recognized(paste0(
    "DSVERT_PEER_NOT_RECOGNIZED_V1|peer=site_b|expected=", expected,
    "|observed=", observed, "f")))
  expect_null(.dsvert_client_parse_peer_not_recognized(
    "ordinary remote error containing no authenticated peer token"))

  inconsistent <- .dsvert_client_parse_peer_not_recognized(paste0(
    "DSVERT_PEER_NOT_RECOGNIZED_V1|peer=site_b|expected=inconsistent",
    "|observed=", observed))
  expect_s3_class(inconsistent, "dsvert_peer_not_recognized")
  expect_identical(inconsistent$expected_fingerprint_sha256, "inconsistent")
  expect_match(conditionMessage(inconsistent), "inconsistent name-bound")
})

test_that("strict aggregate propagates peer rejection without retry or leakage", {
  expected <- paste(rep("1", 64L), collapse = "")
  observed <- paste(rep("2", 64L), collapse = "")
  token <- paste0(
    "DSVERT_PEER_NOT_RECOGNIZED_V1|peer=site_a|expected=", expected,
    "|observed=", observed, " SECRET server trace")
  conns <- list(site_a = structure(1, class = "fake"))
  calls <- 0L
  aggregate <- function(conns, expr, error, ...) {
    calls <<- calls + 1L
    error("site_a", token)
    list(site_a = NULL)
  }

  condition <- tryCatch(
    .dsvert_aggregate_strict(
      conns, call(name = "dsvertJointDPCountStartDS"),
      operation = "pinned phase", .aggregate = aggregate),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(calls, 1L)
  expect_identical(condition$peer_name, "site_a")
  expect_match(conditionMessage(condition), "dsvert.trusted_peers", fixed = TRUE)
  expect_false(grepl("SECRET", conditionMessage(condition), fixed = TRUE))
})

test_that("availability fan-out never hides an unrecognized peer", {
  fingerprint <- paste(rep("c", 64L), collapse = "")
  token <- paste0(
    "DSVERT_PEER_NOT_RECOGNIZED_V1|peer=site_a|expected=unconfigured",
    "|observed=", fingerprint)
  conns <- list(site_a = structure(1, class = "fake"))
  expressions <- list(site_a = call(name = "phaseDS"))
  aggregate <- function(conns, expr, error, ...) {
    error("site_a", token)
    list(site_a = NULL)
  }

  condition <- tryCatch(
    .dsvert_fanout_cycle(
      conns, expressions, operation = "offset phase",
      .aggregate = aggregate),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$expected_fingerprint_sha256, "unconfigured")
})

test_that("a generic top-level DSI failure cannot overwrite a peer rejection", {
  fingerprint <- strrep("d", 64L)
  token <- paste0(
    "DSVERT_PEER_NOT_RECOGNIZED_V1|peer=site_a|expected=unconfigured",
    "|observed=", fingerprint)
  conns <- list(site_a = structure(1, class = "fake"))
  expressions <- list(site_a = call(name = "phaseDS"))
  aggregate <- function(conns, expr, error, ...) {
    error("site_a", token)
    stop("generic connector failure", call. = FALSE)
  }

  condition <- tryCatch(
    .dsvert_fanout_cycle(
      conns, expressions, operation = "pinned phase",
      .aggregate = aggregate),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$peer_name, "site_a")

  assign <- function(conns, symbol, value, error, ...) {
    error("site_a", token)
    stop("generic connector failure", call. = FALSE)
  }
  condition <- tryCatch(
    .dsvert_assign_strict(
      conns, "payload", expressions, operation = "pinned assignment",
      .assign = assign),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$peer_name, "site_a")
})
