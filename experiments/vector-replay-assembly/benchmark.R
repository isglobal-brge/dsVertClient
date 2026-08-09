# Reproducible benchmark for public DP vector replay assembly.
# Run from the dsVertClient repository root:
#   Rscript experiments/vector-replay-assembly/benchmark.R

cases <- data.frame(
  coordinate_count = c(50003L, 100003L, 200003L, 1000000L),
  chunk_coordinates = c(128L, 128L, 128L, 8192L))
repetitions <- 3L

chunks_for <- function(coordinate_count, chunk_coordinates) {
  starts <- seq.int(1L, coordinate_count, by = chunk_coordinates)
  lapply(starts, function(start) {
    as.character(seq.int(
      start, length.out = min(chunk_coordinates, coordinate_count - start + 1L)))
  })
}

append_collect <- function(chunks) {
  result <- character()
  for (chunk in chunks) result <- c(result, chunk)
  result
}

preallocated_collect <- function(chunks, coordinate_count) {
  result <- character(coordinate_count)
  offset <- 0L
  for (chunk in chunks) {
    destination <- seq.int(offset + 1L, length.out = length(chunk))
    result[destination] <- chunk
    offset <- offset + length(chunk)
  }
  result
}

elapsed_median <- function(code) {
  median(replicate(repetitions, system.time(force(code()))[["elapsed"]]))
}

rows <- lapply(seq_len(nrow(cases)), function(case_index) {
  coordinate_count <- cases$coordinate_count[[case_index]]
  chunk_coordinates <- cases$chunk_coordinates[[case_index]]
  chunks <- chunks_for(coordinate_count, chunk_coordinates)
  old <- append_collect(chunks)
  new <- preallocated_collect(chunks, coordinate_count)
  stopifnot(identical(old, new))
  chunk_count <- length(chunks)
  copied_references <- sum(vapply(
    chunks, length, integer(1L)) * rev(seq_len(chunk_count)))
  data.frame(
    coordinate_count = coordinate_count,
    chunk_coordinates = chunk_coordinates,
    chunk_count = chunk_count,
    append_seconds = elapsed_median(function() append_collect(chunks)),
    preallocated_seconds = elapsed_median(
      function() preallocated_collect(chunks, coordinate_count)),
    estimated_append_reference_bytes = copied_references * 8,
    preallocated_reference_bytes = coordinate_count * 8,
    byte_identical = TRUE)
})

print(do.call(rbind, rows), row.names = FALSE)
