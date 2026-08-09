test_that("UUIDv4 session IDs use exactly 16 random bytes", {
  requested <- NULL
  source_bytes <- as.raw(0:15)

  id <- dsVertClient:::.dsvert_uuid4(.random_bytes = function(n) {
    requested <<- n
    source_bytes
  })

  expect_identical(requested, 16L)
  expect_match(
    id,
    "^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$"
  )

  encoded <- as.raw(strtoi(
    substring(gsub("-", "", id, fixed = TRUE),
              seq.int(1L, 31L, by = 2L), seq.int(2L, 32L, by = 2L)),
    base = 16L
  ))
  expect_identical(bitwShiftR(as.integer(encoded[[7L]]), 4L), 4L)
  expect_identical(bitwShiftR(as.integer(encoded[[9L]]), 6L), 2L)
})

test_that("UUIDv4 session IDs are canonically formatted and distinct", {
  ids <- replicate(512L, dsVertClient:::.dsvert_uuid4())

  expect_true(all(grepl(
    "^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$",
    ids
  )))
  expect_length(unique(ids), length(ids))
})

test_that("secure random byte acquisition validates length and fails closed", {
  missing_device <- tempfile("dsvert-no-csprng-")
  bytes <- dsVertClient:::.dsvert_random_bytes(16L)

  expect_type(bytes, "raw")
  expect_length(bytes, 16L)

  expect_error(
    dsVertClient:::.dsvert_random_bytes(
      16L, .urandom_path = missing_device, .openssl_rand_bytes = NULL
    ),
    "cryptographically secure random source is unavailable"
  )
  expect_error(
    dsVertClient:::.dsvert_uuid4(.random_bytes = function(n) raw(15L)),
    "exactly 16 random bytes"
  )
  expect_error(
    dsVertClient:::.dsvert_uuid4(.random_bytes = function(n) rep(0L, 16L)),
    "exactly 16 random bytes"
  )
  expect_error(dsVertClient:::.dsvert_random_bytes(Inf),
               "positive integer")
  expect_error(dsVertClient:::.dsvert_random_bytes(1.5),
               "positive integer")
})
