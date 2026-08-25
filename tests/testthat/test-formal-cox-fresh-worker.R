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

test_that("configured fresh Cox worker permits only an empty completion query", {
  worker <- .formal_cox_fresh_worker_selector()
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      seen <<- expr
      list(site_a = .formal_cox_fresh_worker_reply(
        "completion", list(complete = FALSE)))
    },
    .package = "dsVertClient")
  result <- .dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "completion",
    structure(list(), names = character()), .aggregate = identity)
  expect_identical(as.character(seen[[1L]]),
                   "dsvertFormalCoxWorkerControlDS")
  expect_identical(seen$action, "completion")
  expect_identical(result$payload, list(complete = FALSE))
  expect_error(.dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "completion", list(extra = TRUE)),
    "closed action payload")
})

test_that("configured fresh Cox worker accepts only one matching K2 completion", {
  workers <- list(site_a = .formal_cox_fresh_worker_selector(),
                  site_b = utils::modifyList(.formal_cox_fresh_worker_selector(),
                                              list(peer_name = "site_b")))
  completion <- list(
    version = "dsvert-formal-cox-blockwise-completion-v1",
    plan_sha256 = strrep("a", 64L), transcript_sha256 = strrep("b", 64L),
    final_commit_sha256 = strrep("c", 64L), schedule_steps = 12L,
    fixed_schedule_complete = TRUE, output_kind = "sealed_private_result_v1",
    production_ready = FALSE, completion_sha256 = strrep("d", 64L))
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(conn, worker, action, payload,
                                                     .aggregate) {
      expect_identical(action, "completion")
      expect_identical(payload, structure(list(), names = character()))
      list(payload = list(complete = TRUE, completion = completion))
    },
    .package = "dsVertClient")
  observed <- .dsvert_formal_cox_fresh_worker_completion(
    list(site_a = "connection", site_b = "connection"), workers,
    .aggregate = identity)
  expect_identical(observed, completion)

  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(...) {
      list(payload = list(complete = FALSE))
    },
    .package = "dsVertClient")
  expect_null(.dsvert_formal_cox_fresh_worker_completion(
    list(site_a = "connection", site_b = "connection"), workers,
    .aggregate = identity))

  bad <- completion
  bad$final_commit_sha256 <- strrep("e", 64L)
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(conn, worker, action, payload,
                                                     .aggregate) {
      list(payload = list(complete = TRUE,
                          completion = if (identical(worker$peer_name, "site_a")) {
                            completion
                          } else {
                            bad
                          }))
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_cox_fresh_worker_completion(
    list(site_a = "connection", site_b = "connection"), workers,
    .aggregate = identity), "different completion")
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
    list(site_a = "connection"), worker, "root_claim", list()),
    "closed action")
  expect_error(.dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "accept", list(step = 0L)),
    "bounded opaque frame")
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
