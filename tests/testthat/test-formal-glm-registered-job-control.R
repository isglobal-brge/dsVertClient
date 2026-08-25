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

test_that("registered formal GLM control admits opaque task actions", {
  conn <- list(site_a = structure(list(), class = "mock"))
  receipt <- .formal_glm_registered_job_control_receipt()
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      seen <<- list(conns = conns, expr = expr, operation = operation)
      list(site_a = list(
        version = "dsvert-formal-glm-registered-phase20-job-control-response-v1",
        action = "terminal_status",
        payload = list(state = "running", production_ready = FALSE),
        production_ready = FALSE))
    },
    .package = "dsVertClient")
  reply <- .dsvert_formal_glm_registered_job_control_call(
    conn, receipt, "terminal_status", structure(list(), names = character()),
    .aggregate = identity)
  expect_identical(as.character(seen$expr[[1L]]),
                   "dsvertFormalGLMRegisteredJobControlDS")
  expect_identical(as.list(seen$expr[-1L]), list(
    receipt = receipt, action = "terminal_status",
    payload = structure(list(), names = character())))
  expect_identical(reply$payload, list(state = "running", production_ready = FALSE))
  expect_false(reply$production_ready)
})

test_that("registered formal GLM control relays only opaque Phase21 frames", {
  conn <- list(site_a = structure(list(), class = "mock"))
  receipt <- .formal_glm_registered_job_control_receipt()
  seen <- NULL
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      seen <<- list(expr = expr, operation = operation)
      list(site_a = list(
        version = "dsvert-formal-glm-registered-phase20-job-control-response-v1",
        action = "phase21_ticket",
        payload = list(frame = "eyJyZWNvcmQiOnt9fQ=="),
        production_ready = FALSE))
    },
    .package = "dsVertClient")
  reply <- .dsvert_formal_glm_registered_job_control_call(
    conn, receipt, "phase21_ticket", list(frame = "e30="),
    .aggregate = identity)
  expect_identical(seen$operation, "registered formal GLM job control relay")
  expect_identical(as.list(seen$expr[-1L]), list(
    receipt = receipt, action = "phase21_ticket", payload = list(frame = "e30=")))
  expect_identical(reply$payload, list(frame = "eyJyZWNvcmQiOnt9fQ=="))
  expect_error(.dsvert_formal_glm_registered_job_control_call(
    conn, receipt, "phase21_ticket", structure(list(), names = character())),
    "one opaque frame")
})

test_that("registered formal GLM Phase21 driver preserves the signed lifecycle order", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  receipts <- list(site_a = .formal_glm_registered_job_control_receipt(),
                   site_b = utils::modifyList(
                     .formal_glm_registered_job_control_receipt(),
                     list(peer = "site_b")))
  calls <- list()
  complete_status <- "eyJzdGF0ZSI6ImNvbXBsZXRlIiwicHJvZHVjdGlvbl9yZWFkeSI6ZmFsc2V9"
  testthat::local_mocked_bindings(
    .dsvert_formal_glm_registered_job_control_call = function(
        conn, receipt, action, payload, .aggregate) {
      calls[[length(calls) + 1L]] <<- list(
        peer = names(conn)[[1L]], action = action, payload = payload)
      frame <- if (identical(action, "phase21_stage_status")) {
        complete_status
      } else {
        "eyJyZWNvcmQiOnt9fQ=="
      }
      list(
        version = "dsvert-formal-glm-registered-phase20-job-control-response-v1",
        action = action, payload = list(frame = frame), production_ready = FALSE)
    },
    .package = "dsVertClient")
  result <- .dsvert_formal_glm_registered_phase21_run(
    conns, receipts, .aggregate = identity, max_cycles = 2L)
  expect_identical(result, list(
    state = "public_terminal_complete", production_ready = FALSE))
  expect_identical(vapply(calls, `[[`, character(1L), "action"), c(
    "phase21_preflight", "phase21_preflight",
    "phase21_preflight_bind", "phase21_preflight_bind",
    "phase21_stage_start", "phase21_stage_start",
    "phase21_stage_status", "phase21_stage_status",
    "phase21_stage_record", "phase21_stage_record",
    "phase21_stage_import", "phase21_stage_import",
    "phase21_ticket", "phase21_ticket_import",
    "phase21_seal", "phase21_seal",
    "phase21_seal_import", "phase21_seal_import",
    "phase21_candidate", "phase21_candidate_import",
    "phase21_candidate_verify", "phase21_candidate_verify",
    "phase21_local_release_import", "phase21_local_release_import",
    "phase21_base_certificate", "phase21_base_certificate_import",
    "phase21_authorization", "phase21_authorization_import",
    "phase21_authorization", "phase21_authorization_import",
    "phase21_publication", "phase21_commit", "phase21_commit",
    "phase21_commit_import", "phase21_ack", "phase21_ack_import",
    "phase21_cleanup", "phase21_cleanup", "phase21_cleanup_import"))
  expect_true(all(vapply(calls, function(call) {
    identical(names(call$payload), "frame")
  }, logical(1L))))
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

.formal_glm_registered_job_control_reply <- function(action, payload) {
  list(
    version = "dsvert-formal-glm-registered-phase20-job-control-response-v1",
    action = action, payload = payload, production_ready = FALSE)
}

.formal_glm_registered_job_control_ref <- function() {
  list(
    artifact_id = strrep("a", 64L), receipt_set_sha256 = strrep("b", 64L),
    attempt_id = strrep("c", 64L), job_sha256 = strrep("d", 64L),
    transport_sha256 = strrep("e", 64L), production_ready = FALSE)
}

test_that("registered formal GLM job driver completes the closed two-host lifecycle", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  receipts <- list(site_a = .formal_glm_registered_job_control_receipt(),
                   site_b = utils::modifyList(
                     .formal_glm_registered_job_control_receipt(),
                     list(peer = "site_b")))
  ref <- .formal_glm_registered_job_control_ref()
  polls <- c(site_a = 0L, site_b = 0L)
  relays <- list()
  testthat::local_mocked_bindings(
    .dsvert_formal_glm_registered_job_control_call = function(
        conn, receipt, action, payload, .aggregate) {
      site <- names(conn)[[1L]]
      if (identical(action, "negotiate")) {
        outbound <- if (identical(site, "site_a") && !nzchar(payload$inbound)) {
          "cHJvcG9zYWw="
        } else if (identical(site, "site_b")) {
          "YWNjZXB0"
        } else {
          ""
        }
        return(.formal_glm_registered_job_control_reply(action, list(
          state = "negotiating", outbound = outbound, inspect_only = FALSE,
          production_ready = FALSE)))
      }
      if (identical(action, "start")) {
        return(.formal_glm_registered_job_control_reply(action, list(
          state = "running", outbound = "", inspect_only = FALSE,
          production_ready = FALSE)))
      }
      if (identical(action, "job_ref")) {
        return(.formal_glm_registered_job_control_reply(action, list(
          ref = ref, claim = if (identical(site, "site_a")) "YQ==" else "Yg==")))
      }
      if (action %in% c("bind", "heartbeat")) {
        return(.formal_glm_registered_job_control_reply(
          action, structure(list(), names = character())))
      }
      if (identical(action, "poll")) {
        polls[[site]] <<- polls[[site]] + 1L
        chunk <- if (polls[[site]] == 1L) list(
          job_sha256 = ref$job_sha256, transport_sha256 = ref$transport_sha256,
          offset = 0, payload_sha256 = strrep("f", 64L),
          payload = if (identical(site, "site_a")) "QQ==" else "Qg==") else NULL
        payload <- list(state = "running", accepted_through = payload$acknowledged,
                        production_ready = FALSE)
        if (!is.null(chunk)) payload$relay_chunk <- chunk
        return(.formal_glm_registered_job_control_reply(action, payload))
      }
      if (identical(action, "relay")) {
        relays[[length(relays) + 1L]] <<- list(site = site, payload = payload)
        return(.formal_glm_registered_job_control_reply(action, list(accepted = 1)))
      }
      if (action %in% c("compute_start", "terminal_start")) {
        return(.formal_glm_registered_job_control_reply(action, list(
          state = "running", production_ready = FALSE)))
      }
      if (identical(action, "compute_status")) {
        return(.formal_glm_registered_job_control_reply(action, list(
          state = "complete", production_ready = FALSE)))
      }
      if (identical(action, "terminal_status")) {
        return(.formal_glm_registered_job_control_reply(action, list(
          state = "complete", commit = list(
            SHA256 = strrep("1", 64L), Bytes = 64, Replayed = FALSE),
          production_ready = FALSE)))
      }
      stop("unexpected action", call. = FALSE)
    },
    .package = "dsVertClient")
  result <- .dsvert_formal_glm_registered_job_run(
    conns, receipts, .aggregate = identity, max_cycles = 4L)
  expect_identical(result$state, "terminal_complete")
  expect_false(result$production_ready)
  expect_length(relays, 2L)
  expect_true(all(vapply(relays, function(value)
    identical(value$payload$ref, ref), logical(1L))))
})

test_that("registered formal GLM job driver rejects mismatched peer identities before bind", {
  conns <- list(site_a = structure(list(), class = "mock"),
                site_b = structure(list(), class = "mock"))
  receipts <- list(site_a = .formal_glm_registered_job_control_receipt(),
                   site_b = utils::modifyList(
                     .formal_glm_registered_job_control_receipt(),
                     list(peer = "site_b")))
  calls <- character()
  testthat::local_mocked_bindings(
    .dsvert_formal_glm_registered_job_control_call = function(
        conn, receipt, action, payload, .aggregate) {
      site <- names(conn)[[1L]]
      calls <<- c(calls, action)
      if (identical(action, "negotiate")) {
        return(.formal_glm_registered_job_control_reply(action, list(
          state = "negotiating",
          outbound = if (identical(site, "site_a") && !nzchar(payload$inbound)) {
            "cHJvcG9zYWw="
          } else if (identical(site, "site_b")) "YWNjZXB0" else "",
          inspect_only = FALSE, production_ready = FALSE)))
      }
      if (identical(action, "start")) {
        return(.formal_glm_registered_job_control_reply(action, list(
          state = "running", outbound = "", inspect_only = FALSE,
          production_ready = FALSE)))
      }
      if (identical(action, "job_ref")) {
        ref <- .formal_glm_registered_job_control_ref()
        if (identical(site, "site_b")) ref$job_sha256 <- strrep("9", 64L)
        return(.formal_glm_registered_job_control_reply(action,
          list(ref = ref, claim = "YQ==")))
      }
      stop("unexpected action", call. = FALSE)
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_glm_registered_job_run(
    conns, receipts, .aggregate = identity), "different job references")
  expect_false("bind" %in% calls)
})
