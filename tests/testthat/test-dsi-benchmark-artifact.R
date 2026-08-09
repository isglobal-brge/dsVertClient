test_that("committed DSI scaling benchmark is complete and self-consistent", {
  path <- system.file(
    "docs", "benchmarks", "dsi_transport_dslite_20260731.csv",
    package = "dsVertClient")
  if (!nzchar(path)) {
    path <- testthat::test_path(
      "..", "..", "inst", "docs", "benchmarks",
      "dsi_transport_dslite_20260731.csv")
  }
  results <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c(
    "peers", "payload_mib_per_peer", "unique_raw_mib", "chunk_raw_kib",
    "fanout_cycles", "replay_cycles", "reconnects", "wall_seconds",
    "throughput_unique_mib_s", "raw_submitted_mib",
    "base64_submitted_mib", "expression_submitted_mib",
    "base64_expansion", "max_expression_kib", "rss_start_mib",
    "rss_peak_mib", "rss_delta_mib", "spool_peak_mib",
    "source_disk_mib", "integrity")
  expect_true(all(required %in% names(results)))
  expect_true(all(names(results) %in% c("dispatch_engine", required)))
  expect_equal(nrow(results), 9L)
  expect_setequal(
    paste(results$peers, results$payload_mib_per_peer, sep = "x"),
    as.vector(outer(c(2L, 4L, 8L), c(1L, 16L, 64L), paste, sep = "x")))

  expected_cycles <- ceiling(
    results$payload_mib_per_peer * 1024 / results$chunk_raw_kib) + 1L
  expected_replay_mib <-
    results$peers * pmin(
      results$payload_mib_per_peer,
      results$chunk_raw_kib / 1024)
  expect_equal(results$unique_raw_mib,
               results$peers * results$payload_mib_per_peer)
  expect_equal(results$fanout_cycles, expected_cycles)
  expect_equal(results$raw_submitted_mib,
               results$unique_raw_mib + expected_replay_mib)
  expect_equal(results$replay_cycles, rep(1L, 9L))
  expect_equal(results$reconnects, rep(1L, 9L))
  expect_true(all(results$wall_seconds > 0))
  expect_true(all(results$throughput_unique_mib_s > 0))
  expect_true(all(results$base64_expansion >= 4 / 3))
  expect_true(all(results$base64_expansion < 1.333335))
  expect_true(all(results$max_expression_kib < 641))
  expect_true(all(results$rss_peak_mib >= results$rss_start_mib))
  expect_equal(results$rss_delta_mib,
               results$rss_peak_mib - results$rss_start_mib)
  expect_equal(results$spool_peak_mib, results$unique_raw_mib)
  expect_equal(results$source_disk_mib, results$payload_mib_per_peer)
  expect_identical(results$integrity, rep("sha256_identical", 9L))
})

test_that("committed DSLite frame sweep records the parser ceiling", {
  path <- system.file(
    "docs", "benchmarks", "dsi_transport_frame_sweep_dslite_20260731.csv",
    package = "dsVertClient")
  if (!nzchar(path)) {
    path <- testthat::test_path(
      "..", "..", "inst", "docs", "benchmarks",
      "dsi_transport_frame_sweep_dslite_20260731.csv")
  }
  sweep <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_equal(nrow(sweep), 4L)
  expect_identical(sweep$peers, rep(1L, 4L))
  expect_identical(sweep$payload_mib_per_peer, rep(16L, 4L))
  expect_identical(sweep$chunk_raw_kib, c(512L, 1024L, 2048L, 4096L))
  expect_identical(
    sweep$integrity,
    c("sha256_identical", rep("transport_rejected", 3L)))
  expect_identical(sweep$fanout_cycles[[1L]], 33L)
  expect_identical(sweep$replay_cycles[[1L]], 1L)
  expect_identical(sweep$reconnects[[1L]], 1L)
  expect_true(sweep$throughput_unique_mib_s[[1L]] > 0)
  expect_true(sweep$max_expression_kib[[1L]] < 768)
  expect_true(all(is.na(sweep$throughput_unique_mib_s[-1L])))
})

test_that("adaptive capsule-source benchmark records 1e6-coordinate calls", {
  path <- system.file(
    "docs", "benchmarks", "capsule_source_adaptive_window_20260808.csv",
    package = "dsVertClient")
  if (!nzchar(path)) {
    path <- testthat::test_path(
      "..", "..", "inst", "docs", "benchmarks",
      "capsule_source_adaptive_window_20260808.csv")
  }
  results <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_identical(nrow(results), 9L)
  expect_identical(unique(results$coordinate_count), 1000000L)
  expect_identical(unique(results$chunks_per_source), 123L)
  expect_setequal(unique(results$source_count), c(2L, 3L, 5L))
  expect_setequal(unique(results$window_chunks), c(1L, 2L, 8L))
  windows <- ceiling(results$chunks_per_source / results$window_chunks)
  expect_equal(
    results$data_phases, 2L * results$source_count * windows)
  expect_equal(
    results$data_aggregate_invocations,
    3L * results$source_count * windows)
  expect_equal(results$data_plus_incremental_v2_negotiation_phases,
               results$data_phases + 2L)
  expect_equal(
    results$data_plus_incremental_v2_negotiation_aggregate_invocations,
    results$data_aggregate_invocations + 2L + results$source_count)
  expect_equal(
    results$phase_reduction_vs_scalar,
    1 - results$data_phases /
      (2L * results$source_count * results$chunks_per_source),
    tolerance = 1e-7)
  adaptive <- results[results$window_chunks == 8L, ]
  expect_identical(adaptive$data_phases, c(64L, 96L, 160L))
  expect_true(all(adaptive$phase_reduction_vs_scalar > 0.8699))
})

test_that("balanced direct-dispatch A/B keeps bytes equal and records local timing", {
  path <- system.file(
    "docs", "benchmarks", "dsi_dispatch_ab_dslite_20260801.csv",
    package = "dsVertClient")
  if (!nzchar(path)) {
    path <- testthat::test_path(
      "..", "..", "inst", "docs", "benchmarks",
      "dsi_dispatch_ab_dslite_20260801.csv")
  }
  results <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_equal(nrow(results), 8L)
  counts <- table(results$dispatch_engine)
  expect_identical(names(counts), c("direct", "dsi-standard"))
  expect_identical(unname(as.integer(counts)), c(4L, 4L))
  expect_setequal(results$replicate, 1:4)
  order_balance <- table(
    results$order_position, results$dispatch_engine)
  expect_identical(dim(order_balance), c(2L, 2L))
  expect_identical(as.integer(order_balance), rep(2L, 4L))
  expect_true(all(table(results$replicate, results$dispatch_engine) == 1L))
  invariant <- c(
    "peers", "payload_mib_per_peer", "unique_raw_mib", "chunk_raw_kib",
    "fanout_cycles", "replay_cycles", "reconnects", "raw_submitted_mib",
    "base64_submitted_mib", "expression_submitted_mib", "base64_expansion",
    "max_expression_kib", "spool_peak_mib", "source_disk_mib", "integrity")
  direct <- results[results$dispatch_engine == "direct", invariant]
  standard <- results[results$dispatch_engine == "dsi-standard", invariant]
  rownames(direct) <- NULL
  rownames(standard) <- NULL
  expect_identical(direct, standard)
  direct_wall <- median(results$wall_seconds[
    results$dispatch_engine == "direct"])
  standard_wall <- median(results$wall_seconds[
    results$dispatch_engine == "dsi-standard"])
  # This assertion describes only the committed, warmed, order-balanced local
  # DSLite evidence. It is not a portable Opal/Armadillo speed guarantee.
  expect_lt(direct_wall, standard_wall)
  paired <- merge(
    results[results$dispatch_engine == "direct",
            c("replicate", "wall_seconds")],
    results[results$dispatch_engine == "dsi-standard",
            c("replicate", "wall_seconds")],
    by = "replicate", suffixes = c("_direct", "_standard"))
  expect_true(all(
    paired$wall_seconds_direct < paired$wall_seconds_standard))
  expect_true(all(results$integrity == "sha256_identical"))
})

test_that("typed pump memory benchmark proves its non-streaming source model", {
  path <- system.file(
    "docs", "benchmarks", "typed_blob_memory_20260801.csv",
    package = "dsVertClient")
  if (!nzchar(path)) {
    path <- testthat::test_path(
      "..", "..", "inst", "docs", "benchmarks",
      "typed_blob_memory_20260801.csv")
  }
  results <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_identical(names(results), c(
    "payload_mib", "payload_chars", "frame_chars", "frames", "dsi_calls",
    "wall_seconds", "throughput_mib_s", "max_rss_bytes",
    "baseline_max_rss_bytes", "rss_delta_bytes", "source_memory_model",
    "integrity"))
  expect_identical(results$payload_mib, c(1L, 16L, 64L, 128L))
  expect_equal(results$payload_chars, results$payload_mib * 1024^2)
  expect_identical(results$frame_chars, rep(640L * 1024L, 4L))
  expect_equal(results$frames,
               ceiling(results$payload_chars / results$frame_chars))
  expect_equal(results$dsi_calls, results$frames + 1L)
  expect_true(all(results$wall_seconds > 0))
  expect_true(all(results$throughput_mib_s > 0))
  expect_true(all(results$max_rss_bytes >= results$baseline_max_rss_bytes))
  expect_equal(results$rss_delta_bytes,
               results$max_rss_bytes - results$baseline_max_rss_bytes)
  expect_true(results$rss_delta_bytes[[4L]] >
              results$rss_delta_bytes[[1L]])
  expect_identical(results$source_memory_model, rep("O(payload)", 4L))
  expect_identical(results$integrity, rep("sha256_identical", 4L))
})

test_that("typed source-stream benchmark has no relay payload spool", {
  path <- system.file(
    "docs", "benchmarks", "typed_source_streaming_multiprocess_20260801.csv",
    package = "dsVertClient")
  if (!nzchar(path)) {
    path <- testthat::test_path(
      "..", "..", "inst", "docs", "benchmarks",
      "typed_source_streaming_multiprocess_20260801.csv")
  }
  results <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_identical(results$payload_mib, c(64L, 128L, 256L))
  expect_equal(results$frames,
               ceiling(results$payload_chars / (512 * 1024)))
  expect_identical(results$controller_payload_spool_bytes, rep(0L, 3L))
  expect_equal(results$redundant_integrity_io_bytes_avoided,
               2 * results$payload_chars)
  expect_equal(results$controller_rss_delta_bytes,
               results$controller_peak_rss_bytes -
                 results$controller_baseline_rss_bytes)
  expect_identical(results$producer_frame_replays, rep(1L, 3L))
  expect_identical(results$recipient_frame_replays, rep(1L, 3L))
  expect_identical(results$reconnect_count, rep(1L, 3L))
  expect_true(all(tolower(as.character(
    results$source_removed_after_receipt)) == "true"))
  expect_true(all(results$integrity == "sha256_identical"))
})
