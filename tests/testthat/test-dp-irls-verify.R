.irls_arithmetic_fixtures <- function() {
  lines <- readLines(testthat::test_path("..", "..", "inst", "testdata",
                                       "dp_irls_update_v1.jsonl"), warn = FALSE)
  lapply(lines, jsonlite::fromJSON, simplifyVector = FALSE)
}

test_that("independent exact verifier reproduces all shared canonical answers", {
  fixtures <- .irls_arithmetic_fixtures()
  expect_length(fixtures, 269L)
  expect_equal(sum(grepl("^random-", vapply(fixtures, `[[`, "", "name"))), 240L)
  for (fixture in fixtures) {
    request <- .dsvert_dp_irls_arithmetic_json(fixture$request)
    expect_identical(.dsvert_dp_irls_verify_update_json(request),
      .dsvert_dp_irls_arithmetic_json(fixture$expected), info = fixture$name)
  }
})

test_that("malformed or inconsistent public numeric inputs fail closed", {
  fixture <- .irls_arithmetic_fixtures()[[1L]]
  x <- fixture$request
  cases <- list(
    missing = function(v) { v$floor <- NULL; v },
    extra = function(v) { v$extra <- "1"; v },
    version = function(v) { v$version <- "dp_irls_v1"; v },
    number = function(v) { v$floor$n <- 1; v },
    leading_zero = function(v) { v$floor$n <- "01"; v },
    signed_zero = function(v) { v$floor$n <- "-0"; v },
    positive_sign = function(v) { v$floor$n <- "+1"; v },
    unreduced = function(v) { v$floor <- list(d = "2", n = "2"); v },
    zero_denominator = function(v) { v$floor$d <- "0"; v },
    negative_denominator = function(v) { v$floor$d <- "-1"; v },
    zero_floor = function(v) { v$floor$n <- "0"; v },
    zero_quantum = function(v) { v$floor_quantum$n <- "0"; v },
    zero_lattice = function(v) { v$lattice$n <- "0"; v },
    zero_damping = function(v) { v$damping$n <- "0"; v },
    excess_damping = function(v) { v$damping$n <- "2"; v },
    negative_ridge = function(v) { v$ridge$n <- "-1"; v },
    negative_trust = function(v) { v$trust$n <- "-1"; v },
    outside_box = function(v) { v$beta[[1]]$n <- "1001"; v },
    off_lattice = function(v) { v$beta[[1]] <- list(d = "2", n = "1"); v },
    reversed_box = function(v) { v$box_lower[[1]]$n <- "1001"; v },
    triangle_shape = function(v) { v$h[[1]] <- NULL; v },
    coordinate_object = function(v) { names(v$g) <- c("a", "b"); v },
    negative_payload = function(v) { v$g[[1]] <- "-1"; v },
    overflow_payload = function(v) { v$g[[1]] <- "999"; v },
    zero_translation = function(v) { v$g_translation[[1]] <- "0"; v },
    denominator_certificate = function(v) { v$limits$common_denominator <- "2"; v },
    low_kmax = function(v) { v$limits$k_max <- "0"; v },
    high_kmax = function(v) { v$limits$k_max <- "999"; v },
    low_input_bits = function(v) { v$limits$input_bits <- "1"; v },
    low_temporary_bits = function(v) { v$limits$temporary_bits <- "1"; v },
    high_temporary_bits = function(v) { v$limits$temporary_bits <- "5000000"; v },
    oversized_integer = function(v) { v$g[[1]] <- strrep("9", 1501L); v })
  for (name in names(cases)) {
    expect_error(.dsvert_dp_irls_verify_update(cases[[name]](x)),
                 "Invalid dp_irls_update_v1", info = name)
  }
})

test_that("canonical wire bytes reject whitespace duplicate and reordered keys", {
  request <- .irls_arithmetic_fixtures()[[1L]]$request
  wire <- .dsvert_dp_irls_arithmetic_json(request)
  expect_error(.dsvert_dp_irls_verify_update_json(paste0(wire, "\n")), "noncanonical JSON")
  expect_error(.dsvert_dp_irls_verify_update_json(sub("^\\{", "{ ", wire)), "noncanonical JSON")
  duplicate <- sub('"version":"dp_irls_update_v1"',
    '"version":"dp_irls_update_v1","version":"dp_irls_update_v1"', wire, fixed = TRUE)
  expect_error(.dsvert_dp_irls_verify_update_json(duplicate), "Invalid dp_irls_update_v1")
  reordered <- as.character(jsonlite::toJSON(rev(request), auto_unbox = TRUE))
  expect_error(.dsvert_dp_irls_verify_update_json(reordered), "noncanonical JSON")
  expect_error(.dsvert_dp_irls_verify_update_json("null"), "Invalid dp_irls_update_v1")
  expect_error(.dsvert_dp_irls_verify_update_json("{"), "JSON decoding")
})

test_that("one signed width covers every admitted lattice beta and payload", {
  request <- .irls_arithmetic_fixtures()[[1L]]$request
  original <- request$limits
  request$beta <- request$box_lower
  request$g <- lapply(request$g_translation, function(a) as.character(2L * gmp::as.bigz(a)))
  expect_identical(.dsvert_dp_irls_numeric_limits(request), original)
  expect_silent(.dsvert_dp_irls_verify_update(request))
  request <- .irls_arithmetic_fixtures()[[1L]]$request
  request$limits$input_bits <- as.character(as.integer(original$input_bits) + 1L)
  request$limits$temporary_bits <- as.character(as.integer(original$temporary_bits) +
                                                8L * length(request$beta) + 16L)
  expect_identical(.dsvert_dp_irls_verify_update(request),
                    .irls_arithmetic_fixtures()[[1L]]$expected)
})
