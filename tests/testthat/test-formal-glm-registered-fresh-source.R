test_that("registered fresh GLM source sends selectors but no source contract", {
  selector <- list(
    analysis_id = "fresh_binomial", data_name = "D", family = "binomial",
    formula_sha256 = paste(rep("a", 64L), collapse = ""))
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      seen <<- list(conns = conns, expr = expr, operation = operation)
      list(server_a = list(
        version = "dsvert-formal-glm-registered-fresh-source-response-v1",
        action = "ticket", payload = list(ticket = list(version = "ticket"),
                                            replayed = FALSE),
        production_ready = FALSE))
    },
    .package = "dsVertClient")

  value <- .dsvert_formal_glm_registered_fresh_source_call(
    list(server_a = "connection"), selector, "ticket",
    structure(list(), names = character()), .aggregate = identity)

  expect_identical(value$action, "ticket")
  expect_identical(as.character(seen$conns), "connection")
  expect_identical(as.character(seen$expr[[1L]]),
                   "dsvertFormalGLMRegisteredFreshSourceDS")
  expect_identical(seen$expr$analysis_id, selector$analysis_id)
  expect_identical(seen$expr$data_name, selector$data_name)
  expect_identical(seen$expr$family, selector$family)
  expect_identical(seen$expr$formula_sha256, selector$formula_sha256)
  expect_false("source_contract_json" %in% names(as.list(seen$expr)))
})

test_that("registered fresh GLM source rejects widened selectors and replies", {
  selector <- list(
    analysis_id = "fresh_binomial", data_name = "D", family = "binomial",
    formula_sha256 = paste(rep("a", 64L), collapse = ""))
  expect_error(.dsvert_formal_glm_registered_fresh_source_call(
    list(server_a = "connection"), selector, "ticket", list(extra = TRUE)),
    "closed action payload")
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(...) list(server_a = list(
      version = "wrong", action = "ticket", payload = list(),
      production_ready = FALSE)),
    .package = "dsVertClient")
  expect_error(.dsvert_formal_glm_registered_fresh_source_call(
    list(server_a = "connection"), selector, "ticket",
    structure(list(), names = character())),
    "invalid registered fresh GLM source reply")
})
