test_that("public DP replay assembly is exact across many chunks", {
  coordinate_count <- 100003L
  contract <- list(
    coordinate_count = coordinate_count,
    chunk_coordinates = 128L,
    chunk_count = as.integer(ceiling(coordinate_count / 128L)))
  reads <- 0L

  observed <- .dsvert_vector_collect_replay(contract, function(index) {
    reads <<- reads + 1L
    geometry <- .dsvert_vector_chunk_geometry(contract, index)
    as.character(seq.int(
      geometry$offset + 1L, length.out = geometry$count))
  })

  expect_identical(observed, as.character(seq_len(coordinate_count)))
  expect_identical(reads, contract$chunk_count)
})

test_that("public DP replay assembly rejects malformed chunk shapes", {
  contract <- list(
    coordinate_count = 257L, chunk_coordinates = 128L,
    chunk_count = 3L)

  expect_error(
    .dsvert_vector_collect_replay(contract, function(index) {
      geometry <- .dsvert_vector_chunk_geometry(contract, index)
      character(geometry$count - identical(index, 1L))
    }),
    "wrong shape")
  expect_error(
    .dsvert_vector_collect_replay(contract, NULL),
    "Invalid final DP vector replay reader")
})
