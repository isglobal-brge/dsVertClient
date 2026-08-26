.formal_cox_fresh_ingress_selector <- function() {
  list(analysis_id = "fresh_cox", data_name = "study",
       formula_sha256 = strrep("a", 64L))
}

.formal_cox_fresh_ingress_reply <- function(action, payload) {
  list(version = "dsvert-formal-cox-fresh-source-response-v1",
       action = action, payload = payload, production_ready = FALSE)
}

.formal_cox_fresh_ingress_shape <- function(source) {
  list(
    version = "dsvert-formal-cox-fresh-source-shape-v1",
    analysis_id = "fresh_cox", schema_sha256 = strrep("b", 64L),
    source = source, custodian_peers = c("site_a", "site_b", "site_c"),
    designated_compute_peers = c("site_a", "site_c"), total_blocks = 2L,
    production_ready = FALSE)
}

.formal_cox_fresh_ingress_delivery <- function(recipient) {
  list(
    version = "dsvert-formal-cox-blockwise-source-delivery-v1",
    purpose = "formal-cox-recipient-encrypted-source-delivery-v1",
    receipt = list(version = "receipt"), receipt_sha256 = strrep("c", 64L),
    recipient_peer_name = recipient, envelope = list(ciphertext = "opaque"),
    binding = list(receipt = "opaque"))
}

.formal_cox_fresh_ingress_worker <- function(peer) {
  list(version = "dsvert-formal-cox-blockwise-worker-provision-v1",
       peer_name = peer, plan_sha256 = strrep("d", 64L),
       attempt_id = strrep("e", 64L), replayed = FALSE,
       production_ready = FALSE)
}

test_that("fresh Cox ingress provisions both compute workers from opaque source blocks", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"),
                site_c = structure(list(), class = "mock"))
  calls <- list()
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_source_call = function(
        conn, selector, action, payload, .aggregate) {
      peer <- names(conn)[[1L]]
      calls[[length(calls) + 1L]] <<- list(peer = peer, action = action,
                                            payload = payload)
      value <- switch(action,
        shape = .formal_cox_fresh_ingress_shape(peer),
        ticket = list(ticket = list(version = "ticket", recipient = peer)),
        produce = list(receipt = list(version = "receipt"),
                       receipt_sha256 = strrep("c", 64L), replayed = FALSE),
        delivery = .formal_cox_fresh_ingress_delivery(payload$recipient_peer_name),
        import = list(
          version = "dsvert-formal-cox-blockwise-source-import-receipt-v1",
          purpose = "formal-cox-recipient-encrypted-source-delivery-v1",
          receipt_sha256 = strrep("c", 64L), recipient_peer_name = peer,
          replayed = FALSE),
        provision = .formal_cox_fresh_ingress_worker(peer),
        stop("unexpected action", call. = FALSE))
      .formal_cox_fresh_ingress_reply(action, value)
    },
    .package = "dsVertClient")

  ingress <- .dsvert_formal_cox_fresh_ingress(
    conns, .formal_cox_fresh_ingress_selector(), .aggregate = identity)
  expect_identical(ingress$schema_sha256, strrep("b", 64L))
  expect_identical(ingress$total_blocks, 2L)
  expect_identical(ingress$compute_peers, c("site_a", "site_c"))
  expect_identical(names(ingress$workers), c("site_a", "site_c"))
  expect_false(ingress$production_ready)
  expect_identical(sum(vapply(calls, function(call) call$action == "produce",
                            logical(1L))), 6L)
  expect_identical(sum(vapply(calls, function(call) call$action == "delivery",
                            logical(1L))), 12L)
  expect_identical(sum(vapply(calls, function(call) call$action == "import",
                            logical(1L))), 12L)
  expect_identical(sum(vapply(calls, function(call) call$action == "provision",
                            logical(1L))), 2L)
})

test_that("fresh Cox ingress fails before tickets on incompatible source shapes", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  calls <- character()
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_source_call = function(
        conn, selector, action, payload, .aggregate) {
      source <- names(conn)[[1L]]
      calls <<- c(calls, action)
      if (!identical(action, "shape")) stop("ticket emitted", call. = FALSE)
      shape <- .formal_cox_fresh_ingress_shape(source)
      if (identical(source, "site_b")) shape$schema_sha256 <- strrep("f", 64L)
      .formal_cox_fresh_ingress_reply(action, shape)
    },
    .package = "dsVertClient")

  expect_error(.dsvert_formal_cox_fresh_ingress(
    conns, .formal_cox_fresh_ingress_selector(), .aggregate = identity),
    "incompatible source shapes")
  expect_identical(calls, c("shape", "shape"))
})

test_that("fresh Cox ingress binds the source analysis before tickets", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  calls <- character()
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_source_call = function(
        conn, selector, action, payload, .aggregate) {
      calls <<- c(calls, action)
      shape <- .formal_cox_fresh_ingress_shape(names(conn)[[1L]])
      shape$analysis_id <- "other_cox"
      .formal_cox_fresh_ingress_reply(action, shape)
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_cox_fresh_ingress(
    conns, .formal_cox_fresh_ingress_selector(), .aggregate = identity),
    "incompatible source shapes")
  expect_identical(calls, c("shape", "shape"))
})

test_that("fresh Cox run stages its finalizer without returning a private result", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"),
                site_c = structure(list(), class = "mock"))
  ingress <- list(
    analysis_id = "fresh_cox", schema_sha256 = strrep("b", 64L),
    total_blocks = 2L, compute_peers = c("site_a", "site_c"),
    workers = list(site_a = .formal_cox_fresh_ingress_worker("site_a"),
                   site_c = .formal_cox_fresh_ingress_worker("site_c")),
    production_ready = FALSE)
  completion <- list(
    version = "dsvert-formal-cox-blockwise-completion-v1",
    plan_sha256 = strrep("d", 64L), transcript_sha256 = strrep("e", 64L),
    final_commit_sha256 = strrep("f", 64L), schedule_steps = 4L,
    fixed_schedule_complete = TRUE, output_kind = "sealed_private_result_v1",
    production_ready = FALSE, completion_sha256 = strrep("1", 64L))
  handoff <- list(headers = list(list(role = "garbler"), list(role = "evaluator")),
                  ticket = list(version = "ticket"),
                  envelopes = list(list(version = "envelope"), list(version = "envelope")),
                  production_ready = FALSE)
  expected_handoff <- handoff
  prepared <- list(
    intent = list(
      version = "dsvert-formal-cox-blockwise-sticky-opening-v1",
      purpose = "formal_cox_one_public_beta_validity_opening_v1",
      artifact_id = strrep("2", 64L), candidate_sha256 = strrep("3", 64L),
      final_pair_root_sha256 = strrep("4", 64L),
      opening_mode = "dual_authority_additive_ring_and_xor_validity_v1",
      exp_postprocess_mode = "certified_dyadic_interval_midpoint_v1"),
    finalized = FALSE, certificate_sha256 = "", replayed = FALSE,
    production_ready = FALSE)
  staged <- list(artifact_id = strrep("2", 64L), production_ready = FALSE)
  calls <- character()
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_ingress = function(conns, selector, .aggregate) {
      calls <<- c(calls, "ingress")
      expect_identical(names(conns), c("site_a", "site_b", "site_c"))
      expect_identical(selector, .formal_cox_fresh_ingress_selector())
      ingress
    },
    .dsvert_formal_cox_fresh_worker_run = function(conns, workers, .aggregate) {
      calls <<- c(calls, "worker")
      expect_identical(names(conns), c("site_a", "site_c"))
      expect_identical(names(workers), c("site_a", "site_c"))
      completion
    },
    .dsvert_formal_cox_fresh_worker_finalizer_handoff = function(
        conns, workers, completion, .aggregate) {
      calls <<- c(calls, "handoff")
      expect_identical(names(conns), c("site_a", "site_c"))
      expect_identical(completion$completion_sha256, strrep("1", 64L))
      handoff
    },
    .dsvert_formal_cox_fresh_worker_prepare_finalizer = function(
        conns, workers, handoff, .aggregate) {
      calls <<- c(calls, "prepare")
      expect_identical(handoff, expected_handoff)
      prepared
    },
    .dsvert_formal_cox_fresh_worker_stage_finalizer = function(
        conns, workers, handoff, prepared, .aggregate) {
      calls <<- c(calls, "stage")
      expect_identical(names(conns), c("site_a", "site_c"))
      expect_identical(handoff, expected_handoff)
      expect_identical(prepared$intent$candidate_sha256, strrep("3", 64L))
      staged
    },
    .dsvert_formal_cox_fresh_worker_finish_finalizer = function(
        conns, workers, handoff, .aggregate) {
      calls <<- c(calls, "finish")
      expect_identical(names(conns), c("site_a", "site_c"))
      expect_identical(handoff, expected_handoff)
      list(certificate_sha256 = strrep("4", 64L), production_ready = FALSE)
    },
    .package = "dsVertClient")

  result <- .dsvert_formal_cox_fresh_run(
    conns, .formal_cox_fresh_ingress_selector(), .aggregate = identity)
  expect_identical(calls, c("ingress", "worker", "handoff", "prepare", "stage", "finish"))
  expect_identical(result, list(
    analysis_id = "fresh_cox", schema_sha256 = strrep("b", 64L),
    total_blocks = 2L, state = "finalizer_committed",
    production_ready = FALSE))
  expect_false(any(grepl("intent|candidate|certificate|envelope|header|ticket",
                         names(result), ignore.case = TRUE)))
})

test_that("fresh Cox run fails closed on an unsafe finalizer state", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  ingress <- list(
    analysis_id = "fresh_cox", schema_sha256 = strrep("b", 64L),
    total_blocks = 1L, compute_peers = c("site_a", "site_b"),
    workers = list(site_a = .formal_cox_fresh_ingress_worker("site_a"),
                   site_b = .formal_cox_fresh_ingress_worker("site_b")),
    production_ready = FALSE)
  called_prepare <- FALSE
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_ingress = function(...) ingress,
    .dsvert_formal_cox_fresh_worker_run = function(...) list(
      version = "dsvert-formal-cox-blockwise-completion-v1",
      plan_sha256 = strrep("d", 64L), transcript_sha256 = strrep("e", 64L),
      final_commit_sha256 = strrep("f", 64L), schedule_steps = 1L,
      fixed_schedule_complete = TRUE, output_kind = "sealed_private_result_v1",
      production_ready = FALSE, completion_sha256 = strrep("1", 64L)),
    .dsvert_formal_cox_fresh_worker_finalizer_handoff = function(...) list(
      headers = list(), ticket = list(), envelopes = list(), production_ready = FALSE),
    .dsvert_formal_cox_fresh_worker_prepare_finalizer = function(...) {
      called_prepare <<- TRUE
      list(intent = list(), finalized = FALSE, certificate_sha256 = "unsafe",
           replayed = FALSE, production_ready = FALSE)
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_cox_fresh_run(
    conns, .formal_cox_fresh_ingress_selector(), .aggregate = identity),
    "invalid finalizer state")
  expect_true(called_prepare)
})
