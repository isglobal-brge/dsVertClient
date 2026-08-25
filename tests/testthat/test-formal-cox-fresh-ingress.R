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
