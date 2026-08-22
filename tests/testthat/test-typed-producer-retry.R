test_that("typed producers replay one identical request without attempt quota", {
  conns <- list(site = structure(list(), class = "test"))
  expr <- call(
    name = "k2GradientR1DS", peer_pk = "pinned-peer-key",
    session_id = "12345678-1234-4234-8234-123456789abc")
  observed <- new.env(parent = emptyenv())
  observed$calls <- 0L
  observed$requests <- list()
  aggregate <- function(conns, expr, error, async, errors.print) {
    observed$calls <- observed$calls + 1L
    observed$requests[[observed$calls]] <- expr
    if (observed$calls < 5L) stop("ambiguous top-level transport loss")
    list(site = list(
      encrypted_r1 = "opaque", encrypted_r1_transfer = "contract"))
  }

  result <- testthat::with_mocked_bindings(
    .dsvert_aggregate_strict(
      conns, expr, operation = "typed producer test",
      .aggregate = aggregate),
    .dsvert_retry_sleep = function(seconds) invisible(NULL),
    .dsvert_retry_jitter = function() 1)
  expect_identical(observed$calls, 5L)
  expect_true(all(vapply(
    observed$requests, identical, logical(1L), expr)))
  request_hashes <- vapply(observed$requests, digest::digest, character(1L),
                           algo = "sha256", serialize = TRUE)
  expect_length(unique(request_hashes), 1L)
  expect_identical(result[["site"]]$encrypted_r1, "opaque")
})

test_that("typed producer callback errors and rejections are never retried", {
  conns <- list(site = structure(list(), class = "test"))
  expr <- call(
    name = "k2ShareInputDS", data_name = "D", x_vars = "x",
    peer_pk = "pinned-peer-key", ring = 127L,
    session_id = "12345678-1234-4234-8234-123456789abc")

  callback_calls <- 0L
  callback_failure <- function(conns, expr, error, async, errors.print) {
    callback_calls <<- callback_calls + 1L
    error("site", "server-side contract failure")
    list(site = NULL)
  }
  expect_error(.dsvert_aggregate_strict(
    conns, expr, operation = "typed producer callback test",
    .aggregate = callback_failure), "partial or invalid")
  expect_identical(callback_calls, 1L)

  rejection_calls <- 0L
  rejected <- function(conns, expr, error, async, errors.print) {
    rejection_calls <<- rejection_calls + 1L
    list(site = list(
      version = "dsvert-typed-blob-rejection-v1",
      operation = "producer", rejected = TRUE))
  }
  expect_error(.dsvert_aggregate_strict(
    conns, expr, operation = "typed producer rejection test",
    .aggregate = rejected), "partial or invalid")
  expect_identical(rejection_calls, 1L)
})

test_that("only the audited producer allowlist receives automatic replay", {
  audited <- list(
    call(name = "k2ShareInputDS"),
    call(name = "glmRing63ShareExtraInputDS"),
    call(name = "glmRing63ExportOwnShareDS"),
    call(name = "k2GradientR1DS"),
    call(name = "k2BeaverShareVectorDS", variable = "x"),
    call(name = "k2BeaverVecmulR1DS"),
    call(name = "k2ShareWeightsDS"),
    call(name = "k2IknpBaseSenderChoicesDS"),
    call(name = "k2IknpBaseReceiverEncryptDS"),
    call(name = "k2IknpReceiverExtendDS"),
    call(name = "k2IknpSenderEncryptDS"),
    call(name = "mpcTypedSourceProbeDS"),
    call(name = "psiPaddedRelayExchangeDS"),
    call(name = "exactGCTransportInitDS"),
    call(name = "exactGCBindPeersDS"),
    call(name = "exactGCChisqProductPrepareDS"),
    call(name = "exactGCVecmulClaimInputsDS"),
    call(name = "exactGCVecmulStartDS"),
    call(name = "exactGCVecmulValidityDS"),
    call(name = "exactGCVecmulValidityReceiveDS"),
    call(name = "exactGCVecmulCommitDS"),
    call(name = "k2ChisqCrossAccumulateCountDS"),
    call(name = "dsvertDPCountCompileDS"),
    call(name = "dsvertDPCountAuthorizeDS"),
    call(name = "dsvertDPCountStartDS"),
    call(name = "dsvertDPCountFinalShareDS"),
    call(name = "dsvertDPCountReleaseDS"),
    call(name = "dsvertDPFrequencyClaimDS"),
    call(name = "dsvertDPFrequencyCompileDS"),
    call(name = "dsvertDPFrequencyAuthorizeDS"),
    call(name = "dsvertDPFrequencyCleanupDS"),
    call(name = "dsvertDPFrequencySourceWindowDS"),
    call(name = "dsvertDPFrequencyFinalizeWindowDS"),
    call(name = "dsvertDPFrequencyReplayDS"),
    call(name = "dsvertDPCapsuleManifestDraftDS"),
    call(name = "dsvertDPCapsuleManifestSignDS"),
    call(name = "dsvertDPCapsuleManifestBuildDS"),
    call(name = "dsvertDPCapsuleSourceTicketDS"),
    call(name = "dsvertDPCapsuleSourcePrepareDS"),
    call(name = "dsvertDPCapsuleSourceChunkDS"),
    call(name = "dsvertDPCapsuleSourceAcceptDS"),
    call(name = "dsvertDPAlignmentMaskStartDS"),
    call(name = "dsvertDPAlignmentMaskStoreDS"),
    call(name = "dsvertDPAlignmentMaskSealDS"),
    call(name = "dsvertDPAlignmentMaskReceiveDS"),
    call(name = "dsvertDPCategoricalCrossBindDS"),
    call(name = "dsvertDPCategoricalCrossPrepareDS"),
    call(name = "dsvertDPCategoricalCrossFinalizeDS"),
    call(name = "dsvertDPGaussianCrossBindDS"),
    call(name = "dsvertDPGaussianCrossPrepareDS"),
    call(name = "dsvertDPGaussianCrossFinalizeDS"),
    call(name = "dsvertJointDPVectorPrepareDS"),
    call(name = "dsvertJointDPVectorStartDS"),
    call(name = "dsvertJointDPVectorResultDS"),
    call(name = "dsvertJointDPVectorFinalShareDS"),
    call(name = "dsvertJointDPVectorReleaseDS"),
    call(name = "dsvertJointDPVectorReplayDS"),
    call(name = "dsvertJointDPVectorFinalizeAckDS"),
    call(name = "dsvertDPSynopsisClaimDS"),
    call(name = "dsvertDPSynopsisCompileDS"),
    call(name = "dsvertDPSynopsisPrepareDS"),
    call(name = "dsvertDPSynopsisStartDS"),
    call(name = "dsvertDPSynopsisResultDS"),
    call(name = "dsvertDPSynopsisFinalShareDS"),
    call(name = "dsvertDPSynopsisReleaseDS"),
    call(name = "dsvertDPSynopsisSourceTicketDS"),
    call(name = "dsvertDPSynopsisSourcePrepareDS"),
    call(name = "dsvertDPSynopsisSourceChunkDS"),
    call(name = "dsvertDPSynopsisSourceAcceptDS"),
    call(name = "dsvertDPSynopsisCategoricalCrossBindDS"),
    call(name = "dsvertDPSynopsisCategoricalCrossFinalizeDS"),
    call(name = "dsvertDPSynopsisAlignmentMaskStartDS"),
    call(name = "dsvertDPSynopsisAlignmentMaskStoreDS"),
    call(name = "dsvertDPSynopsisAlignmentMaskSealDS"),
    call(name = "dsvertDPSynopsisAlignmentMaskReceiveDS"),
    call(name = "dsvertDPSynopsisBootstrapDS"),
    call(name = "dsvertDPSynopsisBindDS"),
    call(name = "dsvertDPSynopsisPublicationDS"),
    call(name = "dsvertDPSynopsisPublishedReplayDS"),
    call(name = "dsvertDPSynopsisFinalizeAckDS"))
  expect_identical(length(audited),
                   length(.DSVERT_IDEMPOTENT_TYPED_PRODUCERS))
  expect_true(all(vapply(
    audited, .dsvert_is_idempotent_typed_producer, logical(1L))))
  expect_false(.dsvert_is_idempotent_typed_producer(call(
    name = "k2BeaverReceiveVectorDS", variable = "x")))
  expect_false(.dsvert_is_idempotent_typed_producer(list(
    site_a = call(name = "k2GradientR1DS"),
    site_b = call(name = "k2GradientR2DS"))))

  calls <- 0L
  expect_error(.dsvert_aggregate_strict(
    list(site = structure(list(), class = "test")),
    call(name = "k2GradientR2DS"), operation = "non-producer test",
    .aggregate = function(...) {
      calls <<- calls + 1L
      stop("transport failure")
    }), "partial or invalid")
  expect_identical(calls, 1L)
})

test_that("exact Chisq phases replay an identical request after a lost ACK", {
  phases <- c(
    "exactGCTransportInitDS", "exactGCBindPeersDS",
    "exactGCChisqProductPrepareDS", "exactGCVecmulClaimInputsDS",
    "exactGCVecmulStartDS", "exactGCVecmulValidityDS",
    "exactGCVecmulValidityReceiveDS", "exactGCVecmulCommitDS",
    "k2ChisqCrossAccumulateCountDS")
  conns <- list(site = structure(list(), class = "test"))

  for (phase in phases) {
    request <- call(name = phase, immutable_request_id = "fixed")
    observed <- list()
    aggregate <- function(conns, expr, error, async, errors.print) {
      observed[[length(observed) + 1L]] <<- expr
      if (length(observed) == 1L) stop("simulated lost response")
      list(site = list(stored = TRUE))
    }
    result <- testthat::with_mocked_bindings(
      .dsvert_aggregate_strict(
        conns, request, operation = paste(phase, "lost-ACK test"),
        .aggregate = aggregate),
      .dsvert_retry_sleep = function(seconds) invisible(NULL),
      .dsvert_retry_jitter = function() 1)
    expect_length(observed, 2L)
    expect_identical(observed[[1L]], observed[[2L]])
    expect_true(is.list(result[["site"]]))
  }
})

test_that("canonical Count phases replay an identical request after a lost ACK", {
  phases <- c(
    "dsvertDPCountCompileDS", "dsvertDPCountAuthorizeDS",
    "dsvertDPCountStartDS", "dsvertDPCountFinalShareDS",
    "dsvertDPCountReleaseDS")
  conns <- list(site = structure(list(), class = "test"))

  expect_false(any(phases %in% names(.DSVERT_DSI_TEXT_REMOTE_FORMALS)))
  for (phase in phases) {
    request <- call(name = phase, immutable_request_id = "fixed")
    observed <- list()
    aggregate <- function(conns, expr, error, async, errors.print) {
      observed[[length(observed) + 1L]] <<- expr
      if (length(observed) == 1L) stop("simulated lost Count response")
      list(site = list(stored = TRUE))
    }
    result <- testthat::with_mocked_bindings(
      .dsvert_aggregate_strict(
        conns, request, operation = paste(phase, "lost-ACK test"),
        .aggregate = aggregate),
      .dsvert_retry_sleep = function(seconds) invisible(NULL),
      .dsvert_retry_jitter = function() 1)
    expect_length(observed, 2L)
    expect_identical(observed[[1L]], observed[[2L]])
    expect_true(is.list(result[["site"]]))
  }
})

test_that("Gaussian cross phases replay an identical request after a lost ACK", {
  phases <- c(
    "dsvertDPGaussianCrossBindDS",
    "dsvertDPGaussianCrossPrepareDS",
    "dsvertDPGaussianCrossFinalizeDS")
  conns <- list(site = structure(list(), class = "test"))

  for (phase in phases) {
    request <- call(
      name = phase,
      capsule_id = "cap_fixed",
      immutable_request_id = "fixed")
    observed <- list()
    aggregate <- function(conns, expr, error, async, errors.print) {
      observed[[length(observed) + 1L]] <<- expr
      if (length(observed) == 1L) stop("simulated lost response")
      list(site = list(stored = TRUE))
    }
    result <- testthat::with_mocked_bindings(
      .dsvert_aggregate_strict(
        conns, request, operation = paste(phase, "lost-ACK test"),
        .aggregate = aggregate),
      .dsvert_retry_sleep = function(seconds) invisible(NULL),
      .dsvert_retry_jitter = function() 1)
    expect_length(observed, 2L)
    expect_identical(observed[[1L]], observed[[2L]])
    expect_true(is.list(result[["site"]]))
  }
})

test_that("legacy producer branches never receive automatic replay", {
  legacy <- list(
    call(name = "k2IknpBaseReceiverEncryptDS",
         ciphertexts_blob_key = "legacy-output"),
    call(name = "k2IknpBaseReceiverEncryptDS",
         points_blob_key = "legacy-input"),
    call(name = "k2IknpReceiverExtendDS",
         u_matrix_blob_key = "legacy-output"),
    call(name = "k2IknpSenderEncryptDS",
         u_matrix_blob_key = "legacy-input"),
    call(name = "k2IknpSenderEncryptDS",
         ciphertexts_blob_key = "legacy-output"))
  legacy[[length(legacy) + 1L]] <- as.call(list(
    as.name("k2IknpReceiverExtendDS"), "y", "operation", 10L))
  duplicate <- list(NULL, NULL)
  names(duplicate) <- c("u_matrix_blob_key", "u_matrix_blob_key")
  legacy[[length(legacy) + 1L]] <- as.call(c(
    list(as.name("k2IknpSenderEncryptDS")), duplicate))
  expect_false(any(vapply(
    legacy, .dsvert_is_idempotent_typed_producer, logical(1L))))

  explicit_typed <- list(
    call(name = "k2IknpBaseReceiverEncryptDS",
         points_blob_key = NULL, ciphertexts_blob_key = ""),
    call(name = "k2IknpReceiverExtendDS", u_matrix_blob_key = NULL),
    call(name = "k2IknpSenderEncryptDS", u_matrix_blob_key = NULL,
         ciphertexts_blob_key = ""))
  expect_true(all(vapply(
    explicit_typed, .dsvert_is_idempotent_typed_producer, logical(1L))))
})
