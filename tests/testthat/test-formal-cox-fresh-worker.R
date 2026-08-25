.formal_cox_fresh_worker_selector <- function() {
  list(
    version = "dsvert-formal-cox-blockwise-worker-provision-v1",
    peer_name = "site_a", plan_sha256 = strrep("a", 64L),
    attempt_id = strrep("b", 64L), replayed = FALSE,
    production_ready = FALSE)
}

.formal_cox_fresh_worker_reply <- function(action, payload) {
  list(version = "dsvert-formal-cox-worker-control-response-v1", action = action,
       payload = payload, production_ready = FALSE)
}

test_that("configured fresh Cox worker sends only its provisioned selector", {
  seen <- NULL
  worker <- .formal_cox_fresh_worker_selector()
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      seen <<- list(conns = conns, expr = expr, operation = operation)
      list(site_a = .formal_cox_fresh_worker_reply(
        "host_start", list(version = "status", replayed = FALSE)))
    },
    .package = "dsVertClient")
  result <- .dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "host_start",
    structure(list(), names = character()), .aggregate = identity)
  expect_identical(as.character(seen$expr[[1L]]),
                   "dsvertFormalCoxWorkerControlDS")
  expect_identical(seen$expr$plan_sha256, worker$plan_sha256)
  expect_identical(seen$expr$attempt_id, worker$attempt_id)
  expect_false(any(grepl("key|secret|path|source|config|pid",
                         names(as.list(seen$expr)), ignore.case = TRUE)))
  expect_identical(result$payload$version, "status")
})

test_that("configured fresh Cox worker rejects cross-server, widened and unsafe frames", {
  worker <- .formal_cox_fresh_worker_selector()
  expect_error(.dsvert_formal_cox_fresh_worker_call(
    list(site_b = "connection"), worker, "host_start",
    structure(list(), names = character())), "another server")
  expect_error(.dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "host_start", list(extra = TRUE)),
    "closed action payload")
  expect_error(.dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "relay", list(private_key = "x")),
    "bounded opaque frame")
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(...) list(site_a =
      .formal_cox_fresh_worker_reply("result", list(storage_root = "unsafe"))),
    .package = "dsVertClient")
  expect_error(.dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "result", list()), "unsafe")

  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(...) list(site_a = "not-a-reply"),
    .package = "dsVertClient")
  expect_error(.dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "result", list()), "unsafe")
})
