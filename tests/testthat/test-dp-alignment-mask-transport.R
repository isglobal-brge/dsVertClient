.alignment_mask_client_context <- function(k) {
  servers <- paste0("site_", letters[seq_len(max(k, 2L))])
  list(
    designated = servers[1:2], servers = servers,
    conns = stats::setNames(as.list(seq_len(2L)), servers[1:2]),
    all_conns = stats::setNames(as.list(seq_along(servers)), servers))
}

test_that("private alignment gate uses one fixed K-aware transcript", {
  for (k in c(2L, 3L, 5L)) {
    context <- .alignment_mask_client_context(k)
    chunk_size <- .dsvert_dp_alignment_mask_chunk_size_client(k)
    total <- chunk_size + 3L
    layout <- list(
      source_peers = as.list(context$servers[seq_len(k)]),
      transport_coordinate_count = total)
    source <- list(contract_hash = strrep(as.character(k), 64L))
    commands <- character()
    runs <- list()
    aborted <- character()
    terminal_state <- "complete"

    testthat::local_mocked_bindings(
      .dsvert_exact_gc_new_context = function(...) list(
        operation_id = "op_11111111111111111111111111111111"),
      .dsvert_exact_gc_abort = function(conns, session_id, operation_id,
                                        .aggregate) {
        aborted <<- c(aborted, operation_id)
        invisible(TRUE)
      },
      .dsvert_fanout_by_site = function(
          conns, expressions, operation, .aggregate, ...) {
        command <- unique(vapply(expressions, function(value) {
          as.character(value[[1L]])
        }, character(1L)))
        expect_length(command, 1L)
        commands <<- c(commands, command)
        peers <- names(expressions)
        if (identical(command, "dsvertDPAlignmentMaskStartDS")) {
          return(stats::setNames(lapply(peers, function(peer) {
            list(state = "running")
          }), peers))
        }
        if (identical(command, "dsvertDPAlignmentMaskStoreDS")) {
          arguments <- as.list(expressions[[1L]])[-1L]
          return(stats::setNames(lapply(peers, function(peer) list(
            capability_id = "exact_gc_v1", state = "stored", stored = TRUE,
            chunk_index = arguments$chunk_index,
            chunk_count = arguments$chunk_count)), peers))
        }
        if (identical(command, "dsvertDPAlignmentMaskSealDS")) {
          return(stats::setNames(list(
            list(capability_id = "exact_gc_v1", state = "sealed",
                 peer_blob = "QQ"),
            list(capability_id = "exact_gc_v1", state = "sealed",
                 peer_blob = "Qg")), peers))
        }
        if (identical(command, "dsvertDPAlignmentMaskReceiveDS")) {
          return(stats::setNames(lapply(peers, function(peer) list(
            capability_id = "exact_gc_v1",
            version = "dsvert-alignment-mask-terminal-v1",
            state = terminal_state, terminal_outcome = terminal_state,
            fixed_transcript = TRUE, source_count = k,
            coordinate_count = total, chunk_count = 2L,
            alignment_digest_exposed = FALSE,
            mismatch_source_exposed = FALSE,
            gate_share_exposed = FALSE)), peers))
        }
        stop("unexpected test command")
      },
      .package = "dsVertClient")

    result <- .dsvert_dp_alignment_mask_run(
      "canonical-manifest", context, layout, source,
      "11111111-1111-4111-8111-111111111111", function(...) NULL,
      .run = function(..., operation_id, operation, ring, frac_bits,
                      vector_len, alignment_source_count, initialized) {
        runs[[length(runs) + 1L]] <<- list(
          operation_id = operation_id, operation = operation,
          ring = ring, frac_bits = frac_bits, vector_len = vector_len,
          source_count = alignment_source_count,
          initialized = initialized)
        invisible(list())
      })
    expect_identical(result$state, "complete")
    expect_false(result$alignment_digest_exposed)
    expect_false(result$mismatch_source_exposed)
    expect_false(result$gate_share_exposed)
    expect_identical(result$chunk_count, 2L)
    expect_length(runs, 2L)
    expect_identical(vapply(runs, `[[`, integer(1L), "vector_len"),
                     c(chunk_size, 3L))
    expect_true(all(vapply(runs, `[[`, integer(1L), "source_count") == k))
    expect_true(all(vapply(runs, `[[`, character(1L), "operation") ==
                      "alignment-mask-ring128"))
    expect_identical(commands, c(
      "dsvertDPAlignmentMaskStartDS", "dsvertDPAlignmentMaskStoreDS",
      "dsvertDPAlignmentMaskStartDS", "dsvertDPAlignmentMaskStoreDS",
      "dsvertDPAlignmentMaskSealDS", "dsvertDPAlignmentMaskReceiveDS"))
    expect_length(aborted, 0L)
  }
})

test_that("all mismatch locations collapse to one terminal outcome", {
  fields <- function(state) stats::setNames(lapply(c("site_a", "site_b"),
    function(peer) list(
      capability_id = "exact_gc_v1",
      version = "dsvert-alignment-mask-terminal-v1",
      state = state, terminal_outcome = state, fixed_transcript = TRUE,
      source_count = 5L, coordinate_count = 4099,
      chunk_count = 2L, alignment_digest_exposed = FALSE,
      mismatch_source_exposed = FALSE, gate_share_exposed = FALSE)),
    c("site_a", "site_b"))
  for (mismatch in seq_len(5L)) {
    # `mismatch` deliberately never enters the public response. Every source,
    # content, shape and phase failure must end in this same semantic token.
    expect_identical(.dsvert_dp_alignment_mask_terminal_set(
      fields("alignment_contract_invalid"), c("site_a", "site_b"),
      5L, 4099, 2L), "alignment_contract_invalid")
  }
  condition <- .dsvert_dp_alignment_mask_terminal_error(5L, 4099)
  expect_s3_class(condition, "dsvert_alignment_contract_invalid")
  expect_identical(condition$code, "alignment_contract_invalid")
  expect_false(condition$mismatch_source_exposed)
  expect_false(condition$alignment_digest_exposed)
})

test_that("alignment operation identifiers bind chunk geometry", {
  batch <- "op_11111111111111111111111111111111"
  contract <- strrep("a", 64L)
  ids <- vapply(1:5, function(index) {
    .dsvert_dp_alignment_mask_operation_client(
      batch, contract, index, 5L)
  }, character(1L))
  expect_true(all(grepl("^op_[0-9a-f]{32}$", ids)))
  expect_identical(anyDuplicated(ids), 0L)
  expect_false(identical(ids[[1L]],
                         .dsvert_dp_alignment_mask_operation_client(
                           batch, contract, 1L, 6L)))
})
