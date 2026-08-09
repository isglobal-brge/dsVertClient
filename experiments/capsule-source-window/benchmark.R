#!/usr/bin/env Rscript

# Data-free transport geometry benchmark for the negotiated capsule-source
# byte window. Run from the dsvert-paper repository root.

devtools::load_all("dsVertClient", quiet = TRUE)

make_envelope <- function(recipient, chunk_index, ciphertext_bytes = 132000L) {
  list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_VERSION,
    phase = "encrypted_source_chunk_committed",
    purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE,
    capsule_id = strrep("1", 64L), contract_hash = strrep("2", 64L),
    source_transfer_id = paste0("csrc_", strrep("3", 64L)),
    source_name = "peer_a", source_identity_pk = strrep("A", 43L),
    recipient_name = recipient,
    recipient_identity_pk = if (recipient == "peer_a") {
      strrep("A", 43L)
    } else {
      strrep("B", 43L)
    },
    recipient_ticket_hash = strrep("4", 64L),
    chunk_index = as.numeric(chunk_index), chunk_count = 123,
    coordinate_offset = as.numeric(chunk_index * 8192L),
    coordinates_in_chunk = 8192, chunk_coordinates = 8192,
    ring_bits = 128,
    record_encoding = "little_endian_unsigned_fixed_16_bytes",
    ciphertext_bytes = as.numeric(ciphertext_bytes),
    ciphertext_sha256 = strrep("5", 64L),
    ciphertext = strrep("C", ceiling(4 * ciphertext_bytes / 3)),
    ready_for_sampling = FALSE, signature = strrep("S", 86L))
}

make_bundle <- function(chunk_index) {
  envelopes <- lapply(c("peer_a", "peer_b"), make_envelope,
                      chunk_index = chunk_index)
  first <- envelopes[[1L]]
  list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_BUNDLE_VERSION,
    phase = "encrypted_source_chunk_bundle_committed",
    purpose = first$purpose, capsule_id = first$capsule_id,
    contract_hash = first$contract_hash,
    source_transfer_id = first$source_transfer_id,
    source_name = first$source_name,
    source_identity_pk = first$source_identity_pk,
    recipients = list("peer_a", "peer_b"),
    chunk_index = first$chunk_index, chunk_count = first$chunk_count,
    coordinate_offset = first$coordinate_offset,
    coordinates_in_chunk = first$coordinates_in_chunk,
    chunk_coordinates = first$chunk_coordinates,
    ring_bits = 128,
    record_encoding = first$record_encoding,
    envelopes = envelopes, ready_for_sampling = FALSE)
}

pack_prefix <- function(
    first_index = 0L, maximum_chunks = 8L,
    response_cap = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES) {
  bundles <- list()
  bundle_bytes <- numeric()
  empty <- .dsvert_joint_dp_client_json(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION,
    phase = "encrypted_source_chunk_window_committed",
    bundles = list(), ready_for_sampling = FALSE))
  for (index in seq.int(first_index, length.out = maximum_chunks)) {
    bundle <- make_bundle(index)
    bundle_json <- .dsvert_joint_dp_client_json(bundle)
    candidate_bytes <- nchar(empty, type = "bytes") +
      sum(bundle_bytes) + nchar(bundle_json, type = "bytes") +
      length(bundle_bytes)
    if (candidate_bytes > response_cap) break
    bundles <- c(bundles, list(bundle))
    bundle_bytes <- c(bundle_bytes, nchar(bundle_json, type = "bytes"))
  }
  encoded <- .dsvert_joint_dp_client_json(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION,
    phase = "encrypted_source_chunk_window_committed",
    bundles = bundles, ready_for_sampling = FALSE))
  list(bundles = bundles, encoded = encoded)
}

sample_bundle <- make_bundle(0L)
sample_json <- .dsvert_joint_dp_client_json(sample_bundle)
packed <- pack_prefix()
window_width <- length(packed$bundles)
stopifnot(window_width >= 1L)
stopifnot(identical(
  vapply(packed$bundles, function(value) {
    digest::digest(.dsvert_joint_dp_client_json(value), "sha256",
                   serialize = FALSE)
  }, character(1L)),
  vapply(seq_len(window_width) - 1L, function(index) {
    digest::digest(
      .dsvert_joint_dp_client_json(make_bundle(index)), "sha256",
      serialize = FALSE)
  }, character(1L))))

cases <- expand.grid(
  sources = c(2L, 3L, 5L), chunks_per_source = 123L,
  window_chunks = c(1L, 2L, window_width))
cases$legacy_phases <- with(cases, 2L * sources * chunks_per_source)
cases$window_data_phases <- with(
  cases, 2L * sources * ceiling(chunks_per_source / window_chunks))
cases$data_plus_incremental_v2_negotiation_phases <-
  cases$window_data_phases + 2L
cases$legacy_dsaggregate_invocations <- with(
  cases, 3L * sources * chunks_per_source)
cases$data_plus_incremental_v2_negotiation_aggregate_invocations <- with(
  cases, 3L * sources * ceiling(chunks_per_source / window_chunks) +
    2L + sources)
cases$phase_reduction_vs_scalar <- with(
  cases, 1 - window_data_phases / legacy_phases)

elapsed <- function(code, repetitions = 7L, inner = 5L) {
  median(replicate(repetitions, {
    system.time(for (iteration in seq_len(inner)) code())[["elapsed"]]
  })) / inner
}
legacy_encode_seconds <- elapsed(function() {
  lapply(0:15, function(index) {
    .dsvert_joint_dp_client_json(make_bundle(index))
  })
})
window_encode_seconds <- elapsed(function() {
  start <- seq.int(0L, 15L, by = window_width)
  lapply(start, function(index) pack_prefix(
    index, min(window_width, 16L - index))$encoded)
})

cat("production_like_ciphertext_bytes_per_recipient:", 132000L, "\n")
cat("canonical_pair_bundle_bytes:", nchar(sample_json, type = "bytes"), "\n")
cat("response_probe_padding_bytes:", 8L * 1024L^2, "\n")
cat("response_usable_bytes:",
    .dsvert_dsi_response_probe_usable_bytes(8L * 1024L^2), "\n")
cat("application_window_hard_cap_bytes:",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES, "\n")
cat("effective_window_chunks:", window_width, "\n")
cat("legacy_encode_16_chunks_median_seconds:", legacy_encode_seconds, "\n")
cat("window_encode_16_chunks_median_seconds:", window_encode_seconds, "\n")
print(cases, row.names = FALSE)
