.psi_padded_client_contract <- function(peers, session_id, operation_id,
                                        capacity = 64L) {
  list(
    protocol = "dsvert-pinned-padded-psi-v4",
    session_id = session_id, operation_id = operation_id,
    policy_id = paste0("policy_", strrep("1", 64L)),
    alignment_purpose = "patient-record-alignment-v1",
    dataset_id = "test-logical-cohort",
    dataset_version = "v1",
    id_column = "patient_id",
    source_binding_id = paste0("source_", strrep("6", 64L)),
    pinset_id = paste0("pinset_", strrep("2", 64L)),
    capacity = capacity, relay_frame_bytes = 64L * 1024L,
    inline_max_bytes = 448L * 1024L, peer_names = peers,
    peer_ids = paste0("dsv1_", vapply(seq_along(peers), function(index) {
      paste(rep(as.character(index), 64L), collapse = "")
    }, character(1L))),
    reference_peer = peers[[1L]], compute_peers = peers[1:2],
    snapshot_ids = paste0("snap_", vapply(seq_along(peers), function(index) {
      paste(rep(letters[[index]], 64L), collapse = "")
    }, character(1L))),
    attestation_id = paste0("attest_", strrep("3", 64L)),
    contract_hash = strrep("4", 64L))
}

test_that("padded PSI client follows one fixed K/B-dependent DSI schedule", {
  client_state <- tempfile("dsvert-client-psi-state-")
  dir.create(client_state, mode = "0700")
  on.exit(unlink(client_state, recursive = TRUE, force = TRUE), add = TRUE)
  old <- options(dsvert.client.state_dir = client_state)
  on.exit(options(old), add = TRUE)
  peers <- c("alpha", "beta", "gamma")
  datasources <- stats::setNames(lapply(peers, function(peer) {
    list(peer = peer)
  }), peers)
  session_id <- "12345678-1234-4234-9234-123456789abc"
  operation_id <- paste0("op_", strrep("a", 32L))
  contract <- .psi_padded_client_contract(
    peers, session_id, operation_id, capacity = 64L)
  transcript <- character()
  assign_calls <- list()
  exact_calls <- list()
  aggregate <- function(conns, expr, error, async, errors.print, ...) {
    expressions <- if (is.call(expr)) {
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    } else expr[names(conns)]
    stats::setNames(lapply(names(conns), function(peer) {
      request <- expressions[[peer]]
      method <- as.character(request[[1L]])
      transcript <<- c(transcript, paste(peer, method, sep = ":"))
      switch(method,
        psiPaddedInitDS = list(peer = peer),
        psiPaddedBindDS = list(
          contract = contract,
          receipt = list(unsigned = list(peer = peer), signature = "sig")),
        psiPaddedConfirmDS = list(confirmed = TRUE),
        psiPaddedPrepareDS = list(prepared = TRUE),
        psiPaddedReferenceExportDS = list(
          transport = "inline", envelope = paste0("ref-", peer),
          relay = NULL),
        psiPaddedTargetProcessDS = list(
          transport = "inline", envelope = paste0("target-", peer),
          relay = NULL),
        psiPaddedReferenceDoubleDS = list(
          transport = "inline", envelope = paste0("double-", peer),
          relay = NULL),
        psiPaddedTargetMatchDS = list(transports = stats::setNames(lapply(
          contract$compute_peers, function(compute) list(
            transport = "inline",
            envelope = paste0("membership-", compute), relay = NULL)),
          contract$compute_peers)),
        psiPaddedMembershipAcceptDS = list(accepted = TRUE),
        psiPaddedANDStartDS = list(
          operation_id = paste0("op_", strrep("b", 32L)),
          vector_len = 64L, purpose = "psi.padded.and.test"),
        psiPaddedANDFinalizeDS = list(
          transport = "inline", envelope = paste0("and-", peer),
          relay = NULL),
        psiPaddedANDAcceptDS = list(accepted = TRUE),
        psiPaddedFinalPrepareDS = list(transports = stats::setNames(lapply(
          setdiff(peers, contract$reference_peer), function(target) list(
            transport = "inline", envelope = paste0("final-", target),
            relay = NULL)), setdiff(peers, contract$reference_peer))),
        psiPaddedAttestationDS = list(
          attestation_version = 2L, alignment_attested = TRUE,
          alignment_protocol = contract$protocol,
          attestation_id = contract$attestation_id,
          contract_hash = contract$contract_hash,
          policy_id = contract$policy_id,
          alignment_purpose = contract$alignment_purpose,
          dataset_id = contract$dataset_id,
          dataset_version = contract$dataset_version,
          id_column = contract$id_column,
          source_binding_id = contract$source_binding_id,
          pinset_id = contract$pinset_id,
          capacity_bucket = contract$capacity,
          relay_frame_bytes = contract$relay_frame_bytes,
          inline_max_bytes = contract$inline_max_bytes,
          peer_count = length(peers),
          reference_peer = contract$reference_peer,
          compute_peers = contract$compute_peers),
        mpcCleanupDS = TRUE,
        stop("Unexpected padded PSI test method: ", method))
    }), names(conns))
  }
  assign <- function(conns, symbol, expr, error, success, async,
                     errors.print, ...) {
    assign_calls[[length(assign_calls) + 1L]] <<- list(
      sites = names(conns), symbol = symbol, expr = expr)
    for (peer in names(conns)) success(peer)
    invisible(NULL)
  }
  exact <- function(...) {
    exact_calls[[length(exact_calls) + 1L]] <<- list(...)
    invisible(TRUE)
  }

  attestation <- .dsvert_psi_padded_align(
    "D", "patient_id", "DA", datasources,
    session_id = session_id, operation_id = operation_id,
    .aggregate = aggregate, .assign = assign, .exact_run = exact)

  expect_true(isTRUE(attestation$alignment_attested))
  expect_length(exact_calls, 1L)
  expect_identical(exact_calls[[1L]]$servers, 1:2)
  expect_identical(exact_calls[[1L]]$vector_len, 64L)
  expect_identical(exact_calls[[1L]]$operation, "compare-signed")
  expect_length(assign_calls, 1L)
  expect_identical(assign_calls[[1L]]$sites, peers)
  expect_identical(assign_calls[[1L]]$symbol, "DA")
  expect_identical(
    vapply(assign_calls[[1L]]$expr, function(value) {
      as.character(value[[1L]])
    }, character(1L)),
    stats::setNames(rep("psiPaddedFilterDS", length(peers)), peers))
  # Two public targets imply two repetitions of each pairwise phase and four
  # membership accepts; no call count depends on zero/one/full intersection.
  methods <- sub("^[^:]+:", "", transcript)
  expect_identical(sum(methods == "psiPaddedReferenceExportDS"), 2L)
  expect_identical(sum(methods == "psiPaddedTargetProcessDS"), 2L)
  expect_identical(sum(methods == "psiPaddedReferenceDoubleDS"), 2L)
  expect_identical(sum(methods == "psiPaddedTargetMatchDS"), 2L)
  expect_identical(sum(methods == "psiPaddedMembershipAcceptDS"), 4L)
  expect_identical(sum(methods == "psiPaddedANDAcceptDS"), 2L)
  expect_false("psiPaddedExactTransportDS" %in% methods)
})

test_that("padded PSI client rejects one substituted contract before matching", {
  peers <- c("alpha", "beta", "gamma")
  session_id <- "22345678-1234-4234-9234-123456789abc"
  operation_id <- paste0("op_", strrep("c", 32L))
  contract <- .psi_padded_client_contract(peers, session_id, operation_id)
  bound <- stats::setNames(lapply(peers, function(peer) list(
    contract = contract, receipt = list(peer = peer))), peers)
  substituted <- contract
  substituted$reference_peer <- "beta"
  bound$gamma$contract <- substituted
  expect_error(
    .dsvert_psi_padded_contract(
      bound, peers, session_id, operation_id, "patient_id"),
    "did not agree")
})
