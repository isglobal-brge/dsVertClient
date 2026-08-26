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

test_that("configured fresh Cox worker admits the closed finalizer-stage action", {
  worker <- .formal_cox_fresh_worker_selector()
  testthat::local_mocked_bindings(
    .dsvert_aggregate_strict = function(conns, expr, operation, .aggregate) {
      expect_identical(expr$action, "finalizer_stage")
      list(site_a = .formal_cox_fresh_worker_reply(
        "finalizer_stage", list(artifact_id = strrep("a", 64L))))
    },
    .package = "dsVertClient")
  result <- .dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "finalizer_stage",
    list(ticket = list(version = "ticket"), headers = list(list(), list()),
         envelopes = list(list(), list())), .aggregate = identity)
  expect_identical(result$action, "finalizer_stage")
  expect_error(.dsvert_formal_cox_fresh_worker_call(
    list(site_a = "connection"), worker, "finalizer_stage",
    list(ticket = list(version = "ticket"), headers = list(list(), list())),
    .aggregate = identity), "one ticket, two headers and two envelopes")
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

test_that("configured fresh Cox worker returns only matching share-free opening headers", {
  workers <- list(site_a = .formal_cox_fresh_worker_selector(),
                  site_b = utils::modifyList(.formal_cox_fresh_worker_selector(),
                                              list(peer_name = "site_b")))
  completion <- list(
    version = "dsvert-formal-cox-blockwise-completion-v1",
    plan_sha256 = strrep("a", 64L), transcript_sha256 = strrep("b", 64L),
    final_commit_sha256 = strrep("c", 64L), schedule_steps = 12L,
    fixed_schedule_complete = TRUE, output_kind = "sealed_private_result_v1",
    production_ready = FALSE, completion_sha256 = strrep("d", 64L))
  header <- function(peer, role, local_sha) list(
    version = "dsvert-formal-cox-blockwise-sticky-opening-v1",
    purpose = "formal_cox_one_public_beta_validity_opening_v1",
    artifact_id = strrep("e", 64L), plan_sha256 = strrep("a", 64L),
    run_id = strrep("f", 64L), pinset_sha256 = strrep("0", 64L),
    final_cursor = list(schedule_index = 11L), completion = completion,
    final_receipt = list(version = "receipt"), peer_name = peer,
    peer_id = paste0(peer, "-id"), role = role, coefficient_count = 2L,
    ring_bits = 127L, fraction_bits = 40L,
    local_beta_validity_sha256 = local_sha,
    signature = "AQIDBAUGBwgJCgsMDQ4PEBESExQVFhcYGRobHB0eHyA=")
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(conn, worker, action, payload,
                                                     .aggregate) {
      expect_identical(action, "opening")
      expect_identical(payload, structure(list(), names = character()))
      list(payload = list(
        header = header(worker$peer_name,
                        if (identical(worker$peer_name, "site_a")) "garbler" else "evaluator",
                        if (identical(worker$peer_name, "site_a")) strrep("1", 64L) else strrep("2", 64L)),
        replayed = FALSE))
    },
    .package = "dsVertClient")
  openings <- .dsvert_formal_cox_fresh_worker_openings(
    list(site_a = "connection", site_b = "connection"), workers, completion,
    .aggregate = identity)
  expect_identical(names(openings), names(workers))
  expect_false(any(grepl("share|secret|source|storage|path", names(openings[[1L]]),
                         ignore.case = TRUE)))

  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(conn, worker, action, payload,
                                                     .aggregate) {
      value <- header(worker$peer_name, "garbler", strrep("1", 64L))
      list(payload = list(header = value, replayed = FALSE))
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_cox_fresh_worker_openings(
    list(site_a = "connection", site_b = "connection"), workers, completion,
    .aggregate = identity), "incompatible opening headers")

  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(conn, worker, action, payload,
                                                     .aggregate) {
      value <- header(worker$peer_name,
                      if (identical(worker$peer_name, "site_a")) "garbler" else "evaluator",
                      strrep("1", 64L))
      value$final_receipt$coefficient_shares <- "unsafe"
      list(payload = list(header = value, replayed = FALSE))
    },
    .package = "dsVertClient")
  expect_error(.dsvert_formal_cox_fresh_worker_openings(
    list(site_a = "connection", site_b = "connection"), workers, completion,
    .aggregate = identity), "unsafe opening header")
})

test_that("configured fresh Cox worker relays only a finalizer ticket and ciphertext envelopes", {
  workers <- list(site_a = .formal_cox_fresh_worker_selector(),
                  site_b = utils::modifyList(.formal_cox_fresh_worker_selector(),
                                              list(peer_name = "site_b")))
  completion <- list(
    version = "dsvert-formal-cox-blockwise-completion-v1",
    plan_sha256 = strrep("a", 64L), transcript_sha256 = strrep("b", 64L),
    final_commit_sha256 = strrep("c", 64L), schedule_steps = 12L,
    fixed_schedule_complete = TRUE, output_kind = "sealed_private_result_v1",
    production_ready = FALSE, completion_sha256 = strrep("d", 64L))
  header <- function(peer, role, local_sha) list(
    version = "dsvert-formal-cox-blockwise-sticky-opening-v1",
    purpose = "formal_cox_one_public_beta_validity_opening_v1",
    artifact_id = strrep("e", 64L), plan_sha256 = strrep("a", 64L),
    run_id = strrep("f", 64L), pinset_sha256 = strrep("0", 64L),
    final_cursor = list(schedule_index = 11L), completion = completion,
    final_receipt = list(version = "receipt"), peer_name = peer,
    peer_id = paste0(peer, "-id"), role = role, coefficient_count = 2L,
    ring_bits = 127L, fraction_bits = 40L,
    local_beta_validity_sha256 = local_sha,
    signature = "AQIDBAUGBwgJCgsMDQ4PEBESExQVFhcYGRobHB0eHyA=")
  ticket <- list(
    version = "dsvert-typed-finalizer-handoff-ticket-v1", family = "formal_cox",
    purpose = "formal_cox_blockwise_sticky_opening_v1", artifact_id = strrep("e", 64L),
    final_pair_root_sha256 = strrep("1", 64L), plan_sha256 = strrep("a", 64L),
    pinset_sha256 = strrep("0", 64L), finalizer_peer_name = "site_a",
    finalizer_peer_id = "site_a-id", finalizer_role = "garbler",
    recipient_x25519_public_key = strrep("A", 44L), transport_key_sha256 = strrep("2", 64L),
    issuer_peer_name = "site_a", issuer_peer_id = "site_a-id", signature = strrep("A", 88L))
  intent <- list(
    version = "dsvert-formal-cox-blockwise-sticky-opening-v1",
    purpose = "formal_cox_one_public_beta_validity_opening_v1",
    artifact_id = strrep("e", 64L), candidate_sha256 = strrep("8", 64L),
    final_pair_root_sha256 = strrep("1", 64L),
    opening_mode = "dual_authority_additive_ring_and_xor_validity_v1",
    exp_postprocess_mode = "certified_dyadic_interval_midpoint_v1")
  calls <- character()
  bad_stage <- FALSE
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(conn, worker, action, payload,
                                                     .aggregate) {
      calls <<- c(calls, paste(worker$peer_name, action, sep = ":"))
      role <- if (identical(worker$peer_name, "site_a")) "garbler" else "evaluator"
      if (identical(action, "opening")) {
        return(list(payload = list(header = header(worker$peer_name, role,
                                                    if (role == "garbler") strrep("3", 64L) else strrep("4", 64L)),
                                   replayed = FALSE)))
      }
      if (identical(action, "finalizer_ticket")) {
        expect_identical(worker$peer_name, "site_a")
        expect_null(names(payload$headers))
        return(list(payload = list(ticket = ticket, replayed = FALSE)))
      }
      if (identical(action, "finalizer_seal")) {
        expect_identical(payload$ticket, ticket)
        sender <- header(worker$peer_name, role,
                         if (role == "garbler") strrep("3", 64L) else strrep("4", 64L))
        envelope <- list(
          version = "dsvert-typed-finalizer-handoff-envelope-v1", family = "formal_cox",
          purpose = ticket$purpose, artifact_id = ticket$artifact_id,
          final_pair_root_sha256 = ticket$final_pair_root_sha256,
          plan_sha256 = ticket$plan_sha256, pinset_sha256 = ticket$pinset_sha256,
          ticket_sha256 = strrep("5", 64L), finalizer_peer_name = ticket$finalizer_peer_name,
          finalizer_peer_id = ticket$finalizer_peer_id,
          recipient_transport_key_sha256 = ticket$transport_key_sha256,
          sender_peer_name = sender$peer_name, sender_peer_id = sender$peer_id,
          sender_role = sender$role, payload_kind = "formal_cox_opening_local_output_v1",
          payload_sha256 = strrep("6", 64L), ciphertext_sha256 = strrep("7", 64L),
          ciphertext = strrep("A", 80L), signature = strrep("A", 88L))
        return(list(payload = list(envelope = envelope, replayed = FALSE)))
      }
      if (identical(action, "finalizer_prepare")) {
        expect_identical(worker$peer_name, "site_a")
        expect_identical(payload$ticket, ticket)
        expect_null(names(payload$headers))
        expect_null(names(payload$envelopes))
        return(list(payload = list(intent = intent, finalized = FALSE,
                                   certificate_sha256 = "", replayed = FALSE)))
      }
      if (identical(action, "finalizer_stage")) {
        expect_identical(payload$ticket, ticket)
        expect_null(names(payload$headers))
        expect_null(names(payload$envelopes))
        return(list(payload = list(
          artifact_id = intent$artifact_id,
          candidate_sha256 = if (identical(role, "garbler")) {
            if (isTRUE(bad_stage)) strrep("9", 64L) else intent$candidate_sha256
          } else {
            ""
          },
          local_role = role, production_ready = FALSE)))
      }
      stop("unexpected action", call. = FALSE)
    },
    .package = "dsVertClient")
  handoff <- .dsvert_formal_cox_fresh_worker_finalizer_handoff(
    list(site_a = "connection", site_b = "connection"), workers, completion,
    .aggregate = identity)
  expect_false(handoff$production_ready)
  expect_identical(names(handoff$envelopes), names(workers))
  prepared <- .dsvert_formal_cox_fresh_worker_prepare_finalizer(
    list(site_a = "connection", site_b = "connection"), workers, handoff,
    .aggregate = identity)
  expect_false(prepared$production_ready)
  expect_false(prepared$finalized)
  expect_identical(prepared$intent, intent)
  staged <- .dsvert_formal_cox_fresh_worker_stage_finalizer(
    list(site_a = "connection", site_b = "connection"), workers, handoff,
    prepared, .aggregate = identity)
  expect_identical(staged, list(
    artifact_id = intent$artifact_id, production_ready = FALSE))
  bad_stage <- TRUE
  expect_error(.dsvert_formal_cox_fresh_worker_stage_finalizer(
    list(site_a = "connection", site_b = "connection"), workers, handoff,
    prepared, .aggregate = identity), "invalid finalizer stage")
  expect_identical(calls, c("site_a:opening", "site_b:opening",
                            "site_a:finalizer_ticket", "site_a:finalizer_seal",
                            "site_b:finalizer_seal", "site_a:finalizer_prepare",
                            "site_a:finalizer_stage", "site_b:finalizer_stage",
                            "site_a:finalizer_stage"))
  expect_false(any(grepl("share|secret|storage|path|source", names(handoff),
                         ignore.case = TRUE)))
})

test_that("configured fresh Cox worker relays one K2 schedule without outputs", {
  workers <- list(site_a = .formal_cox_fresh_worker_selector(),
                  site_b = utils::modifyList(.formal_cox_fresh_worker_selector(),
                                              list(peer_name = "site_b")))
  completion <- list(
    version = "dsvert-formal-cox-blockwise-completion-v1",
    plan_sha256 = strrep("a", 64L), transcript_sha256 = strrep("b", 64L),
    final_commit_sha256 = strrep("c", 64L), schedule_steps = 1L,
    fixed_schedule_complete = TRUE, output_kind = "sealed_private_result_v1",
    production_ready = FALSE, completion_sha256 = strrep("d", 64L))
  calls <- character()
  committed <- FALSE
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(conn, worker, action, payload,
                                                     .aggregate) {
      calls <<- c(calls, paste(worker$peer_name, action, sep = ":"))
      reply <- switch(action,
        host_start = list(version = "status", replayed = FALSE),
        bind = list(),
        completion = if (committed) {
          list(complete = TRUE, completion = completion)
        } else {
          list(complete = FALSE)
        },
        offer = list(frame = "AQ=="),
        accept = list(frame = "Ag=="),
        confirm = list(),
        poll = list(chunk = NULL, accepted = "0"),
        result = list(receipt = list(version = "receipt", peer = worker$peer_name),
                      done = TRUE),
        commit = {
          committed <<- TRUE
          list()
        },
        stop("unexpected action", call. = FALSE))
      list(payload = reply)
    },
    .package = "dsVertClient")
  observed <- .dsvert_formal_cox_fresh_worker_run(
    list(site_a = "connection", site_b = "connection"), workers,
    .aggregate = identity)
  expect_identical(observed, completion)
  expect_identical(calls, c(
    "site_a:host_start", "site_b:host_start", "site_a:bind", "site_b:bind",
    "site_a:completion", "site_b:completion", "site_a:offer", "site_b:accept",
    "site_a:confirm", "site_a:poll", "site_b:poll", "site_a:result",
    "site_b:result", "site_a:commit", "site_b:commit", "site_a:completion",
    "site_b:completion"))
  expect_false(any(grepl("coefficient|share|secret|source", names(observed),
                         ignore.case = TRUE)))
})

test_that("configured fresh Cox worker relays offsets above 2^53 as strings", {
  workers <- list(site_a = .formal_cox_fresh_worker_selector(),
                  site_b = utils::modifyList(.formal_cox_fresh_worker_selector(),
                                              list(peer_name = "site_b")))
  completion <- list(
    version = "dsvert-formal-cox-blockwise-completion-v1",
    plan_sha256 = strrep("a", 64L), transcript_sha256 = strrep("b", 64L),
    final_commit_sha256 = strrep("c", 64L), schedule_steps = 1L,
    fixed_schedule_complete = TRUE, output_kind = "sealed_private_result_v1",
    production_ready = FALSE, completion_sha256 = strrep("d", 64L))
  committed <- FALSE
  poll_a <- 0L
  relayed <- NULL
  testthat::local_mocked_bindings(
    .dsvert_formal_cox_fresh_worker_call = function(conn, worker, action, payload,
                                                     .aggregate) {
      reply <- switch(action,
        host_start = list(version = "status", replayed = FALSE),
        bind = list(),
        completion = if (committed) {
          list(complete = TRUE, completion = completion)
        } else {
          list(complete = FALSE)
        },
        offer = list(frame = "AQ=="),
        accept = list(frame = "Ag=="),
        confirm = list(),
        poll = {
          if (identical(worker$peer_name, "site_a") && poll_a == 0L) {
            poll_a <<- poll_a + 1L
            list(chunk = list(sender = "site_a", offset = "9007199254740993",
                              payload_sha256 = strrep("e", 64L), payload = "AQID"),
                 accepted = "0")
          } else {
            list(chunk = NULL, accepted = "0")
          }
        },
        relay = {
          relayed <<- payload$chunk
          list(accepted = "9007199254740996")
        },
        result = list(receipt = list(version = "receipt", peer = worker$peer_name),
                      done = TRUE),
        commit = {
          committed <<- TRUE
          list()
        },
        stop("unexpected action", call. = FALSE))
      list(payload = reply)
    },
    .package = "dsVertClient")
  .dsvert_formal_cox_fresh_worker_run(
    list(site_a = "connection", site_b = "connection"), workers,
    .aggregate = identity)
  expect_identical(relayed$offset, "9007199254740993")
  expect_false(is.numeric(relayed$offset))
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
