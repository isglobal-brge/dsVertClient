test_that("configured fresh Cox source sends only closed selectors and payloads", {
  selector <- list(analysis_id = "fresh_cox", data_name = "study",
                   formula_sha256 = strrep("a", 64L))
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      seen <<- list(conns = conns, expr = expr, operation = operation)
      list(site_a = list(
        version = "dsvert-formal-cox-fresh-source-response-v1", action = "ticket",
        payload = list(ticket = list(version = "ticket")), production_ready = FALSE))
    },
    .package = "dsVertClient")
  reply <- .dsvert_formal_cox_fresh_source_call(
    list(site_a = "connection"), selector, "ticket",
    structure(list(), names = character()), .aggregate = identity)
  expect_identical(reply$action, "ticket")
  expect_identical(as.character(seen$expr[[1L]]), "dsvertFormalCoxFreshSourceDS")
  expect_identical(seen$expr$analysis_id, selector$analysis_id)
  expect_identical(seen$expr$data_name, selector$data_name)
  expect_identical(seen$expr$formula_sha256, selector$formula_sha256)
  expect_false(any(c("path", "key", "source") %in% names(as.list(seen$expr))))
})

test_that("configured fresh Cox source rejects widened payloads and unsafe replies", {
  selector <- list(analysis_id = "fresh_cox", data_name = "study",
                   formula_sha256 = strrep("a", 64L))
  expect_error(.dsvert_formal_cox_fresh_source_call(
    list(site_a = "connection"), selector, "ticket", list(extra = TRUE)),
    "closed action payload")
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(...) list(site_a = list(
      version = "dsvert-formal-cox-fresh-source-response-v1", action = "ticket",
      payload = list(storage_root = "unsafe"), production_ready = FALSE)),
    .package = "dsVertClient")
  expect_error(.dsvert_formal_cox_fresh_source_call(
    list(site_a = "connection"), selector, "ticket",
    structure(list(), names = character())), "unsafe")
})
