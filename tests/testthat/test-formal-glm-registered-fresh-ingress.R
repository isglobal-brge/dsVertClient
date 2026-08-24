.formal_glm_registered_fresh_ingress_selector <- function() {
  list(
    analysis_id = "fresh_binomial", data_name = "D", family = "binomial",
    formula_sha256 = strrep("a", 64L))
}

.formal_glm_registered_fresh_ingress_host <- function(peer) {
  list(
    version = "dsvert-formal-glm-registered-phase20-job-host-provision-v1",
    peer = peer, artifact_id = strrep("b", 64L),
    receipt_set_sha256 = strrep("c", 64L), config_sha256 = strrep("d", 64L),
    replayed = FALSE, production_ready = FALSE)
}

.formal_glm_registered_fresh_ingress_reply <- function(action, payload) {
  list(
    version = "dsvert-formal-glm-registered-phase18-source-response-v1",
    action = action, payload = payload, production_ready = FALSE)
}

.formal_glm_registered_fresh_ingress_shape <- function(source) {
  list(
    version = "dsvert-formal-glm-registered-fresh-source-shape-v1",
    artifact_id = strrep("b", 64L), source_contract_sha256 = strrep("c", 64L),
    source = source, custodian_peers = c("site_a", "site_b"),
    designated_compute_peers = c("site_a", "site_b"), total_blocks = 1L,
    production_ready = FALSE)
}

.formal_glm_registered_fresh_ingress_chunk <- function(source) {
  list(
    version = "chunk", purpose = "opaque", handle = strrep("e", 64L),
    artifact_id = strrep("b", 64L), source_contract_sha256 = strrep("c", 64L),
    authorization_sha256 = strrep("d", 64L), source = source, block_index = 0L,
    pair_sha256 = strrep("f", 64L), pair_bytes = 1L, offset = 0L,
    chunk_sha256 = strrep("1", 64L), chunk_bytes = 1L, complete = TRUE,
    production_ready = FALSE)
}

test_that("registered fresh GLM ingress seals K source blocks and provisions two hosts", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  selector <- .formal_glm_registered_fresh_ingress_selector()
  calls <- list()
  testthat::local_mocked_bindings(
    .dsvert_formal_glm_registered_fresh_source_call = function(
        conn, selector, action, payload, .aggregate) {
      source <- names(conn)[[1L]]
      calls[[length(calls) + 1L]] <<- list(source = source, action = action,
                                            payload = payload)
      result <- switch(action,
        shape = .formal_glm_registered_fresh_ingress_shape(source),
        ticket = list(ticket = list(version = "ticket", recipient = source),
                      replayed = FALSE),
        ticket_set = list(ticket_receipts = list(list(version = "a"), list(version = "b")),
                          replayed = FALSE),
        seal_block = list(source_receipt = list(version = "source"), replayed = FALSE),
        chunk = list(chunk_receipt = .formal_glm_registered_fresh_ingress_chunk(source),
                     pair_chunk_base64 = "QQ==", replayed = FALSE),
        import_chunk = list(chunk_delivery = list(version = "delivery"), replayed = FALSE),
        local_receipt = list(local_receipt_json = paste0("{\"source\":\"", source, "\"}"),
                             replayed = FALSE),
        receipt_commit = list(local_receipt_json = payload$local_receipt_json,
                              replayed = FALSE),
        receipt_set = list(receipt_set_json = "{\"version\":\"set\"}", replayed = FALSE),
        binding = list(binding_record_json = "{\"version\":\"binding\"}", replayed = FALSE),
        host_provision = .formal_glm_registered_fresh_ingress_host(source),
        stop("unexpected action", call. = FALSE))
      .formal_glm_registered_fresh_ingress_reply(action, result)
    },
    .package = "dsVertClient")
  ingress <- .dsvert_formal_glm_registered_fresh_ingress(
    conns, selector, .aggregate = identity)
  expect_identical(names(ingress$hosts), c("site_a", "site_b"))
  expect_identical(ingress$artifact_id, strrep("b", 64L))
  expect_identical(ingress$total_blocks, 1L)
  expect_false(ingress$production_ready)
  expect_identical(sum(vapply(calls, function(value)
    identical(value$action, "import_chunk"), logical(1L))), 4L)
  expect_true(all(vapply(calls, function(value)
    !identical(value$action, "ticket") || value$source %in% c("site_a", "site_b"),
    logical(1L))))
})

test_that("registered fresh GLM ingress fails before tickets on incompatible source shape", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  calls <- character()
  testthat::local_mocked_bindings(
    .dsvert_formal_glm_registered_fresh_source_call = function(
        conn, selector, action, payload, .aggregate) {
      source <- names(conn)[[1L]]
      calls <<- c(calls, action)
      if (!identical(action, "shape")) stop("ticket emitted", call. = FALSE)
      shape <- .formal_glm_registered_fresh_ingress_shape(source)
      if (identical(source, "site_b")) shape$artifact_id <- strrep("9", 64L)
      .formal_glm_registered_fresh_ingress_reply(action, shape)
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_glm_registered_fresh_ingress(
    conns, .formal_glm_registered_fresh_ingress_selector(), .aggregate = identity),
    "incompatible source shapes")
  expect_identical(calls, c("shape", "shape"))
})

test_that("registered fresh GLM run composes ingress with only the designated hosts", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"),
                site_c = structure(list(), class = "mock"))
  ingress <- list(
    artifact_id = strrep("b", 64L), source_contract_sha256 = strrep("c", 64L),
    total_blocks = 3L, compute_peers = c("site_a", "site_c"),
    hosts = list(site_a = .formal_glm_registered_fresh_ingress_host("site_a"),
                 site_c = .formal_glm_registered_fresh_ingress_host("site_c")),
    production_ready = FALSE)
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_formal_glm_registered_fresh_ingress = function(
        conns, selector, .aggregate) ingress,
    .dsvert_formal_glm_registered_job_run = function(
        conns, receipts, .aggregate, max_cycles) {
      seen <<- list(conns = conns, receipts = receipts, max_cycles = max_cycles)
      list(state = "terminal_complete", production_ready = FALSE)
    },
    .package = "dsVertClient")
  result <- .dsvert_formal_glm_registered_fresh_run(
    conns, .formal_glm_registered_fresh_ingress_selector(), .aggregate = identity,
    max_cycles = 4L)
  expect_identical(names(seen$conns), c("site_a", "site_c"))
  expect_identical(names(seen$receipts), c("site_a", "site_c"))
  expect_identical(seen$max_cycles, 4L)
  expect_identical(result, list(
    artifact_id = strrep("b", 64L), total_blocks = 3L,
    state = "phase21_publication_pending", production_ready = FALSE))
})
