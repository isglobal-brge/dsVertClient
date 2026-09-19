test_that("staged LMM receipts bind both peers to plan terminal and candidate shape", {
  artifact <- list(analysis_id = "lmm", coordinate_count = 10L,
    numeric_certificate = list(profile_sha256 = strrep("a", 64), certificate_sha256 = strrep("b", 64)))
  peers <- c("peer_a", "peer_b")
  context <- list(designated = peers)
  checked <- character()
  local_mocked_bindings(.dsvert_dp_capsule_source_verify = function(value, purpose, peer, context) {
    checked <<- c(checked, peer)
    expect_identical(purpose, "cross-grid-result")
  }, .package = "dsVertClient")
  receipt <- function(peer, phase, extra = list()) c(list(
    version = "dsvert-lmm-staged-receipt-v1", phase = phase, capsule_id = strrep("c", 64),
    analysis_id = "lmm", peer_name = peer, peer_identity_pk = paste0("pk-", peer),
    semantic_key = strrep("d", 64), artifact_sha256 = .dsvert_dp_capsule_source_hash(artifact),
    source_contract_sha256 = strrep("e", 64), profile_sha256 = strrep("a", 64),
    certificate_sha256 = strrep("b", 64), stage_plan_digest = strrep("f", 64),
    private_result_exposed = FALSE, signature = "signature"), extra)
  wire <- function(values) lapply(values, .dsvert_joint_dp_client_json)
  extra <- list(operation_id = paste0("op_", strrep("a", 32)), source_key = paste0("exact_gc_in_", strrep("a", 32)),
    output_key = paste0("exact_gc_out_", strrep("a", 32)), operation = "grouped-lmm-staged-v1",
    purpose = paste0("grouped-lmm-staged-v1/", strrep("1", 64)), vector_len = 10L, persisted = TRUE)
  values <- setNames(lapply(peers, receipt, phase = "prepared", extra = extra), peers)
  result <- .dsvert_dp_lmm_cross_receipts(wire(values), context, artifact, "prepared")
  expect_true(attr(result, "both_persisted"))
  expect_setequal(checked, peers)
  values[[2L]]$persisted <- FALSE
  expect_false(attr(.dsvert_dp_lmm_cross_receipts(wire(values), context, artifact, "prepared"), "both_persisted"))
  values[[2L]]$stage_plan_digest <- strrep("2", 64)
  expect_error(.dsvert_dp_lmm_cross_receipts(wire(values), context, artifact, "prepared"))
  values[[2L]]$stage_plan_digest <- values[[1L]]$stage_plan_digest
  values[[2L]]$vector_len <- 1L
  expect_error(.dsvert_dp_lmm_cross_receipts(wire(values), context, artifact, "prepared"))
  values <- setNames(lapply(peers, receipt, phase = "complete",
    extra = list(coordinate_count = 10L, stage_receipt = strrep("3", 64))), peers)
  result <- .dsvert_dp_lmm_cross_receipts(wire(values), context, artifact, "complete")
  expect_identical(result$stage_receipt, strrep("3", 64))
  values[[2L]]$stage_receipt <- strrep("4", 64)
  expect_error(.dsvert_dp_lmm_cross_receipts(wire(values), context, artifact, "complete"))
  values[[2L]]$stage_receipt <- strrep("3", 64)
  values[[1L]]$validity_share <- "AA=="
  expect_error(.dsvert_dp_lmm_cross_receipts(wire(values), context, artifact, "complete"))
})
