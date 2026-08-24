.formal_glm_registered_source_chunk_receipt <- function() {
  list(
    version = "dsvert-formal-glm-registered-phase18-source-outbox-chunk-v3",
    purpose = "formal_glm_registered_source_outbox_chunk_v3",
    handle = strrep("a", 64L), artifact_id = strrep("b", 64L),
    source_contract_sha256 = strrep("c", 64L),
    authorization_sha256 = strrep("d", 64L), source = "site_a",
    block_index = 0L, pair_sha256 = strrep("e", 64L), pair_bytes = 6L,
    offset = 0L, chunk_sha256 = strrep("f", 64L), chunk_bytes = 6L,
    complete = TRUE, production_ready = FALSE)
}

test_that("registered formal GLM source client calls the closed relay", {
  conn <- list(site_a = structure(list(), class = "mock"))
  receipt <- .formal_glm_registered_source_chunk_receipt()
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      seen <<- list(conns = conns, expr = expr, operation = operation,
                    aggregate = .aggregate)
      list(site_a = list(
        version = "dsvert-formal-glm-registered-phase18-source-response-v1",
        action = "chunk", payload = list(
          chunk_receipt = receipt, pair_chunk_base64 = "b3BhcXVl",
          replayed = TRUE), production_ready = FALSE))
    },
    .package = "dsVertClient")
  reply <- .dsvert_formal_glm_registered_source_call(
    conn, "{\"contract\":\"registered\"}", "chunk", list(
      recipient_tickets = unname(list(list(ticket = "a"), list(ticket = "b"))),
      block_index = 0L, offset = 0L), .aggregate = identity)
  expect_identical(seen$conns, conn)
  expect_identical(seen$operation, "registered formal GLM source relay")
  expect_identical(as.character(seen$expr[[1L]]),
                   "dsvertFormalGLMRegisteredSourceDS")
  expect_identical(as.list(seen$expr[-1L]), list(
    source_contract_json = "{\"contract\":\"registered\"}", action = "chunk",
    payload = list(
      recipient_tickets = unname(list(list(ticket = "a"), list(ticket = "b"))),
      block_index = 0L, offset = 0L)))
  expect_identical(reply$payload$pair_chunk_base64, "b3BhcXVl")
  expect_false(reply$production_ready)
})

test_that("registered formal GLM source client fails closed", {
  conn <- list(site_a = structure(list(), class = "mock"))
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(...) {
      calls <<- calls + 1L
      list(site_a = list(
        version = "dsvert-formal-glm-registered-phase18-source-response-v1",
        action = "ticket", payload = list(storage_key = "forbidden"),
        production_ready = FALSE))
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_glm_registered_source_call(
    conn, "{\"contract\":\"registered\"}", "unknown",
    structure(list(), names = character())), "closed action")
  expect_identical(calls, 0L)
  expect_error(.dsvert_formal_glm_registered_source_call(
    conn, "{\"contract\":\"registered\"}", "ticket",
    structure(list(), names = character())),
    "unsafe registered formal GLM source reply")
  expect_identical(calls, 1L)
})
