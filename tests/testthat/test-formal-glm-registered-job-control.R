.formal_glm_registered_job_control_receipt <- function() {
  list(
    version = "dsvert-formal-glm-registered-phase20-job-host-provision-v1",
    peer = "site_a", artifact_id = strrep("a", 64L),
    receipt_set_sha256 = strrep("b", 64L), config_sha256 = strrep("c", 64L),
    replayed = TRUE, production_ready = FALSE)
}

test_that("registered formal GLM control uses the closed server relay", {
  conn <- list(site_a = structure(list(), class = "mock"))
  receipt <- .formal_glm_registered_job_control_receipt()
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      seen <<- list(conns = conns, expr = expr, operation = operation,
                    aggregate = .aggregate)
      list(site_a = list(
        version = "dsvert-formal-glm-registered-phase20-job-control-response-v1",
        action = "start",
        payload = list(state = "running", outbound = "AQ==",
                       inspect_only = FALSE, production_ready = FALSE),
        production_ready = FALSE))
    },
    .package = "dsVertClient")
  reply <- .dsvert_formal_glm_registered_job_control_call(
    conn, receipt, "start", structure(list(), names = character()),
    .aggregate = identity)
  expect_identical(seen$conns, conn)
  expect_identical(seen$operation, "registered formal GLM job control relay")
  expect_identical(as.character(seen$expr[[1L]]),
                   "dsvertFormalGLMRegisteredJobControlDS")
  expect_identical(as.list(seen$expr[-1L]), list(
    receipt = receipt, action = "start",
    payload = structure(list(), names = character())))
  expect_identical(reply$payload$outbound, "AQ==")
  expect_false(reply$production_ready)
})

test_that("registered formal GLM control fails closed before and after DSI", {
  conn <- list(site_a = structure(list(), class = "mock"))
  receipt <- .formal_glm_registered_job_control_receipt()
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(...) {
      calls <<- calls + 1L
      list(site_a = list(
        version = "dsvert-formal-glm-registered-phase20-job-control-response-v1",
        action = "health", payload = list(storage_key = "forbidden"),
        production_ready = FALSE))
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_glm_registered_job_control_call(
    conn, receipt, "unknown", structure(list(), names = character())),
    "closed action")
  expect_identical(calls, 0L)
  expect_error(.dsvert_formal_glm_registered_job_control_call(
    conn, receipt, "health", structure(list(), names = character())),
    "unsafe registered formal GLM control reply")
  expect_identical(calls, 1L)
})
