test_that("exact client pump retries identical fan-out requests without a quota", {
  peer_a <- paste0("dsv1_", strrep("a", 64L))
  peer_b <- paste0("dsv1_", strrep("b", 64L))
  context <- strrep("c", 64L)
  session_id <- "12345678-1234-4234-9234-123456789abc"
  operation_id <- "op_33333333333333333333333333333333"
  servers <- c("site_a", "site_b")
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  peers <- c(site_a = peer_a, site_b = peer_b)
  other <- c(site_a = peer_b, site_b = peer_a)
  roles <- c(site_a = "garbler", site_b = "evaluator")
  output_bytes <- c(site_a = 2, site_b = 1)
  inbound <- c(site_a = 0, site_b = 0)
  exchange_calls <- 0L
  failed_request <- NULL
  retried_identically <- FALSE
  partial_request <- NULL
  partial_retried_identically <- FALSE
  abort_calls <- 0L
  encoded_envelopes <- 0L
  clock <- 0
  initialized <- stats::setNames(lapply(servers, function(server) list(
    capability_id = "exact_gc_v1", peer_id = peers[[server]],
    peer_peer_id = other[[server]], role = roles[[server]],
    context_hash = context, operation = "truncate-floor",
    output_kind = "ring-share", purpose = "test.client-pump",
    source_producer = "test.vecmul",
    ring_bits = 127L, frac_bits = 2L, vector_len = 2L,
    threshold = "", chunk_bytes = 65536L,
    ttl_seconds = if (server == "site_a") 10L else 12L,
    max_runtime_seconds = 120L, worker_heartbeat = 1,
    state = "running", stored = FALSE)), servers)

  negotiated <- .dsvert_exact_gc_validate_init(
    initialized, servers, "truncate-floor", 127L, 2L, 2L,
    "test.client-pump")
  expect_identical(negotiated$ttl_seconds, 10L)

  envelope <- function(server, offset) {
    source <- peers[[server]]
    target <- other[[server]]
    list(
      version = "dsvert-exact-gc-envelope-v1",
      capability_id = "exact_gc_v1", session_id = session_id,
      operation_id = operation_id, context_hash = context,
      sender_peer_id = source, recipient_peer_id = target,
      offset = as.numeric(offset), chunk_bytes = 1,
      payload_hash = strrep(if (server == "site_a") "d" else "e", 64L),
      payload = if (server == "site_a") "AQ" else "Ag",
      signature = "opaque-signature")
  }
  aggregate <- function(conns, expr, ...) {
    if (is.list(expr) && length(expr) == length(conns) &&
        all(vapply(expr, is.call, logical(1L)))) {
      exchange_calls <<- exchange_calls + 1L
      clock <<- clock + 0.2
      if (exchange_calls == 1L) {
        failed_request <<- expr
        stop("simulated transient DSI loss")
      }
      if (exchange_calls == 2L) {
        retried_identically <<- identical(expr, failed_request)
      }
      if (exchange_calls == 3L) {
        partial_request <<- expr
        return(stats::setNames(
          list(list(state = "ignored-partial-response"), NULL), servers))
      }
      if (exchange_calls == 4L) {
        partial_retried_identically <<- identical(expr, partial_request)
      }
      response <- list()
      for (server in names(conns)) {
        site_call <- expr[[server]]
        expect_identical(as.character(site_call[[1L]]), "exactGCExchangeDS")
        self <- peers[[server]]
        expect_true(is.logical(site_call$long_poll))
        expect_identical(site_call$peer_id, self)
        expect_false(other[[server]] %in% names(site_call))
        empty_delivery <- !nzchar(site_call$delivery_payload)
        expect_identical(empty_delivery, all(!nzchar(c(
          site_call$delivery_payload_hash, site_call$delivery_payload,
          site_call$delivery_signature))))
        if (empty_delivery) {
          expect_equal(site_call$delivery_offset, 0)
          expect_equal(site_call$delivery_chunk_bytes, 0)
        }
        delivery <- if (empty_delivery) NULL else list(
          offset = site_call$delivery_offset,
          chunk_bytes = site_call$delivery_chunk_bytes)
        if (!is.null(delivery)) {
          inbound[[server]] <<- delivery$offset + delivery$chunk_bytes
        }
        offset <- as.numeric(site_call$read_offset)
        outgoing <- if (offset < output_bytes[[server]])
          envelope(server, offset) else NULL
        peer_done <- inbound[[server]] >= output_bytes[[setdiff(servers, server)]]
        own_done <- offset >= output_bytes[[server]]
        done <- peer_done && own_done
        response[[server]] <- list(
          capability_id = "exact_gc_v1", peer_id = self,
          state = if (done) "complete" else "running", stored = done,
          inbound_size = inbound[[server]], outbound = outgoing,
          worker_heartbeat = 1 + exchange_calls)
      }
      return(response)
    }
    command <- as.character(expr[[1L]])
    if (identical(command, "exactGCAbortDS")) {
      abort_calls <<- abort_calls + 1L
      return(stats::setNames(as.list(rep(TRUE, length(conns))), names(conns)))
    }
    stop("unexpected mock command: ", command)
  }

  result <- testthat::with_mocked_bindings(
    .dsvert_exact_gc_run(
      conns, server_names = servers, servers = 1:2,
      session_id = session_id, operation_id = operation_id,
      source_key = "exact_gc_in_33333333333333333333333333333333",
      output_key = "exact_gc_out_33333333333333333333333333333333",
      operation = "truncate-floor", ring = 127L, frac_bits = 2L,
      vector_len = 2L, purpose = "test.client-pump",
      transport_ready = TRUE, timeout_seconds = 1, initialized = initialized,
      .aggregate = aggregate),
    .dsvert_exact_gc_monotonic_seconds = function() clock,
    .dsvert_exact_gc_delivery_fields = function(envelope) {
      if (is.null(envelope)) {
        return(list(
          delivery_offset = 0,
          delivery_chunk_bytes = 0,
          delivery_payload_hash = "",
          delivery_payload = "",
          delivery_signature = ""))
      }
      encoded_envelopes <<- encoded_envelopes + 1L
      list(
        delivery_offset = envelope$offset,
        delivery_chunk_bytes = envelope$chunk_bytes,
        delivery_payload_hash = envelope$payload_hash,
        delivery_payload = envelope$payload,
        delivery_signature = envelope$signature)
    },
    .dsvert_exact_gc_sleep = function(seconds) {
      clock <<- clock + min(seconds, 0.01)
    }, .package = "dsVertClient")
  expect_true(retried_identically)
  expect_true(partial_retried_identically)
  expect_gt(exchange_calls, 1L)
  expect_gt(clock, 1)
  expect_identical(abort_calls, 0L)
  expect_gte(encoded_envelopes, as.integer(sum(output_bytes)))
  expect_identical(result$capability_id, "exact_gc_v1")
  expect_false(any(grepl("share|result", names(result), ignore.case = TRUE)))
})

.client_exact_gc_analysis_identity_pk <- function(index) {
  encoded <- jsonlite::base64_enc(as.raw(rep(as.integer(index), 32L)))
  chartr("+/", "-_", sub("=+$", "", encoded))
}

.client_exact_gc_analysis_contract <- function(k = 3L) {
  owners <- paste0("site_", seq_len(k))
  pins <- setNames(vapply(
    seq_along(owners), .client_exact_gc_analysis_identity_pk, character(1L)),
    owners)
  contract <- list(
    version = "dsvert-analysis-contract-v1",
    artifact_key = strrep("0", 64L),
    semantic = list(
      version = "dsvert-analysis-semantic-v1",
      domain = "study-domain",
      cohort_id = "cohort-v1",
      owner_snapshots = setNames(lapply(seq_along(owners), function(index) {
        list(
          version = "dsvert-analysis-snapshot-v1",
          dataset_id = "cohort_table",
          dataset_version = "v1",
          snapshot_commitment = strrep(sprintf("%x", index), 64L))
      }), unname(pins)),
      noise_authorities = unname(pins[seq_len(2L)]),
      analysis = list(
        primitive = "joint-dp-laplace-v2",
        formula = NULL,
        effective_arguments = list(
          statistic = "admitted_privacy_unit_count")),
      privacy = list(
        version = "dsvert-per-analysis-dp-v1",
        adjacency = "add_remove_patient",
        privacy_unit = "patient",
        contribution = list(
          version = "dsvert-contribution-policy-v1",
          max_records_per_unit = 1,
          overflow_policy = "reject_operation",
          constraints = list(
            version = "dsvert-contribution-constraints-v1",
            policy_sha256 = strrep("c", 64L))),
        mechanism = list(
          family = "discrete_laplace",
          version = "discrete-laplace-output-perturbation-tv-v2",
          sensitivity = list(
            version = "dsvert-sensitivity-v1", norm = "l1", value = 1),
          calibration = list(
            version = "dsvert-calibration-v1", noise_scale = 1,
            sampler = "hkdf-sha256-aes128ctr-two-geometric-tv-v2",
            implementation_delta = 1e-9),
          randomness = list(
            version = "dsvert-randomness-plan-v1",
            lanes = list(final_noise = list(
              version = "dsvert-randomness-lane-v1",
              purpose = "privatize_final_vector",
              primitive = "hkdf-sha256-aes128ctr-two-geometric-tv-v2",
              coordinates = 1)))),
        epsilon = 1,
        delta = 1e-6),
      numeric = list(
        version = "dsvert-numeric-semantics-v1",
        value_bits = 127,
        fractional_bits = 0,
        rounding = "toward_zero",
        overflow = "reject",
        sampler_encoding = "aes128ctr_integer_coordinate_v2",
        output_encoding = "twos_complement_integer_v1"),
      public_shape = list(count = 1)),
    execution = list(
      version = "dsvert-analysis-execution-v1",
      peer_pins = as.list(pins),
      backend = list(
        kernel = "joint-dp-laplace-v2",
        ring = "ring127",
        build_sha256 = strrep("a", 64L)),
      transport = list(chunk_coordinates = 4096)))
  contract$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(contract$semantic)
  contract
}

test_that("client derives the same execution-free Count binding for K=2,3,5", {
  bindings <- lapply(c(2L, 3L, 5L), function(k) {
    .dsvert_exact_gc_analysis_binding(
      .client_exact_gc_analysis_contract(k))
  })
  expect_true(all(vapply(bindings, function(value) {
    identical(value$binding$artifact_key, value$contract$artifact_key) &&
      identical(sort(unname(unlist(value$binding$authority_roles)),
                       method = "radix"),
                sort(unlist(
                  value$contract$semantic$noise_authorities),
                  method = "radix"))
  }, logical(1L))))
  expect_identical(
    bindings[[2L]]$sha256,
    "a06bf0f0e116adeac41990142cf11cf95aff778fbc00983fd6ab612b9cbd0221")

  changed_execution <- bindings[[2L]]$contract
  changed_execution$execution$backend$build_sha256 <- strrep("b", 64L)
  changed_execution$execution$transport$chunk_coordinates <- 8192
  expect_identical(
    .dsvert_exact_gc_analysis_binding(changed_execution)$binding,
    bindings[[2L]]$binding)
  expect_false(any(vapply(c(
    "session", "operation", "transport", "build_sha256", "ring",
    "chunk_coordinates"), grepl, logical(1L),
    x = .dsvert_joint_dp_client_json(bindings[[2L]]$binding), fixed = TRUE)))
})

test_that("exact client emits one complete direct scalar delivery shape", {
  expect_identical(.dsvert_exact_gc_delivery_fields(), list(
    delivery_offset = 0,
    delivery_chunk_bytes = 0,
    delivery_payload_hash = "",
    delivery_payload = "",
    delivery_signature = ""))
  envelope <- list(
    version = "dsvert-exact-gc-envelope-v1",
    capability_id = "exact_gc_v1",
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = "op_33333333333333333333333333333333",
    context_hash = strrep("c", 64L),
    sender_peer_id = paste0("dsv1_", strrep("a", 64L)),
    recipient_peer_id = paste0("dsv1_", strrep("b", 64L)),
    offset = 17,
    chunk_bytes = 1,
    payload_hash = strrep("d", 64L),
    payload = "AQ",
    signature = "opaque-signature")
  expect_identical(.dsvert_exact_gc_delivery_fields(envelope), list(
    delivery_offset = 17,
    delivery_chunk_bytes = 1,
    delivery_payload_hash = strrep("d", 64L),
    delivery_payload = "AQ",
    delivery_signature = "opaque-signature"))
  envelope$context_hash <- NULL
  expect_error(.dsvert_exact_gc_delivery_fields(envelope), "malformed")
})

test_that("worker heartbeats renew inactivity but never the total deadline", {
  servers <- c("site_a", "site_b")
  peers <- c(
    site_a = paste0("dsv1_", strrep("a", 64L)),
    site_b = paste0("dsv1_", strrep("b", 64L)))
  initialized <- stats::setNames(lapply(servers, function(server) list(
    capability_id = "exact_gc_v1", peer_id = peers[[server]],
    peer_peer_id = peers[[setdiff(servers, server)]],
    role = if (server == "site_a") "garbler" else "evaluator",
    context_hash = strrep("c", 64L), operation = "truncate-floor",
    output_kind = "ring-share", purpose = "test.total-runtime",
    source_producer = "test.producer", ring_bits = 127L, frac_bits = 2L,
    vector_len = 1L, threshold = "", chunk_bytes = 65536L,
    ttl_seconds = 10L, max_runtime_seconds = 10L,
    worker_heartbeat = 1, state = "running", stored = FALSE)), servers)
  conns <- stats::setNames(list(
    structure(list(), class = "mock"),
    structure(list(), class = "mock")), servers)
  clock <- 0
  exchanges <- 0L
  aborts <- 0L
  aggregate <- function(conns, expr, ...) {
    if (is.list(expr) && all(vapply(expr, is.call, logical(1L)))) {
      exchanges <<- exchanges + 1L
      clock <<- clock + 1
      return(stats::setNames(lapply(servers, function(server) list(
        capability_id = "exact_gc_v1", peer_id = peers[[server]],
        state = "running", stored = FALSE, inbound_size = 0,
        outbound = NULL, worker_heartbeat = 1 + exchanges)), servers))
    }
    aborts <<- aborts + 1L
    stats::setNames(as.list(rep(TRUE, length(conns))), names(conns))
  }

  condition <- tryCatch(testthat::with_mocked_bindings(
    .dsvert_exact_gc_run(
      conns, server_names = servers, servers = 1:2,
      session_id = "12345678-1234-4234-9234-123456789abc",
      operation_id = "op_33333333333333333333333333333333",
      source_key = "exact_gc_in_33333333333333333333333333333333",
      output_key = "exact_gc_out_33333333333333333333333333333333",
      operation = "truncate-floor", ring = 127L, frac_bits = 2L,
      vector_len = 1L, purpose = "test.total-runtime",
      transport_ready = TRUE, timeout_seconds = 1,
      initialized = initialized, .aggregate = aggregate),
    .dsvert_exact_gc_monotonic_seconds = function() clock,
    .dsvert_exact_gc_sleep = function(...) NULL,
    .package = "dsVertClient"), error = identity)
  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "total runtime lease")
  expect_gt(exchanges, 10L)
  expect_gt(aborts, 0L)
})

test_that("present malformed exact response is fatal and raw errors do not leak", {
  servers <- c("site_a", "site_b")
  peer_a <- paste0("dsv1_", strrep("a", 64L))
  peer_b <- paste0("dsv1_", strrep("b", 64L))
  peers <- c(site_a = peer_a, site_b = peer_b)
  initialized <- stats::setNames(lapply(servers, function(server) list(
    capability_id = "exact_gc_v1", peer_id = peers[[server]],
    peer_peer_id = peers[[setdiff(servers, server)]],
    role = if (server == "site_a") "garbler" else "evaluator",
    context_hash = strrep("c", 64L), operation = "truncate-floor",
    output_kind = "ring-share", purpose = "test.malformed",
    source_producer = "test.producer", ring_bits = 127L, frac_bits = 2L,
    vector_len = 1L, threshold = "", chunk_bytes = 65536L,
    ttl_seconds = 10L, max_runtime_seconds = 120L,
    worker_heartbeat = 1,
    state = "running", stored = FALSE)), servers)
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  exchanges <- 0L
  aggregate <- function(conns, expr, ...) {
    if (is.list(expr) && all(vapply(expr, is.call, logical(1L)))) {
      exchanges <<- exchanges + 1L
      return(stats::setNames(rep(list(list(
        private_server_detail = "patient-row-should-never-leak")), 2L),
        servers))
    }
    stats::setNames(as.list(rep(TRUE, length(conns))), names(conns))
  }
  caught <- tryCatch({
    .dsvert_exact_gc_run(
      conns, server_names = servers, servers = 1:2,
      session_id = "12345678-1234-4234-9234-123456789abc",
      operation_id = "op_33333333333333333333333333333333",
      source_key = "exact_gc_in_33333333333333333333333333333333",
      output_key = "exact_gc_out_33333333333333333333333333333333",
      operation = "truncate-floor", ring = 127L, frac_bits = 2L,
      vector_len = 1L, purpose = "test.malformed", transport_ready = TRUE,
      timeout_seconds = 1, initialized = initialized,
      .aggregate = aggregate)
    NULL
  }, error = identity)
  expect_s3_class(caught, "error")
  expect_identical(exchanges, 1L)
  expect_match(conditionMessage(caught), "invalid exchange state")
  expect_false(grepl("patient-row|private_server", conditionMessage(caught)))
})

test_that("exact capability contract advertises only proved rings and floor truncation", {
  contract <- .dsvert_exact_gc_capability_contract()
  expect_identical(contract$supported_ring_bits, 63L:4096L)
  expect_identical(contract$wire_container_bits,
                   c(64L, 128L, 256L, 512L, 1024L, 2048L, 4096L))
  expect_identical(contract$truncation_semantics, "floor")
  expect_true(all(unlist(contract[c(
    "available", "allowed", "canonical_encoding",
    "canonical_input_encoding", "shape_bounds_enforced", "exact_truncation",
    "count_guard_e2e_verified", "clamp_count_e2e_verified",
    "joint_dp_count_e2e_verified", "joint_dp_vector_e2e_verified",
    "alignment_mask_e2e_verified",
    "multiprecision_truncation_e2e_verified")],
    use.names = FALSE)))
  expect_false(contract$exact_comparison)
  expect_false(contract$fail_closed_overflow)
  expect_false(contract$runtime_bounds_enforced)
  expect_false(contract$raw_product_overflow_guard)
  expect_false(contract$checked_mul_truncate)
  expect_false(contract$vecmul_truncation_e2e_verified)
  expect_true(contract$dynamic_ring_fallback)
  expect_match(contract$vecmul_numeric_precondition, "producer-minted")
  expect_true(contract$core_exact_comparison)
  expect_false(contract$comparison_e2e_verified)
  expect_false(contract$workload_glm_e2e_verified)
  expect_identical(contract$operations,
                   c("truncate-floor", "count-guard", "clamp-count",
                     "joint-dp-laplace-v2",
                     "joint-dp-vector-laplace-v3",
                     "alignment-mask-ring128"))
  expect_identical(contract$max_ring_bits, 4096L)
  expect_identical(contract$max_frac_bits, 4095L)
})

test_that("formal GLM schedule is transport-only and remains unadvertised", {
  servers <- c("site_a", "site_b")
  peers <- c(
    site_a = paste0("dsv1_", strrep("a", 64L)),
    site_b = paste0("dsv1_", strrep("b", 64L)))
  purpose <- paste0("formal-glm/phase19/", strrep("c", 64L))
  initialized <- stats::setNames(lapply(servers, function(server) list(
    capability_id = "exact_gc_v1", peer_id = peers[[server]],
    peer_peer_id = peers[[setdiff(servers, server)]],
    role = if (server == "site_a") "garbler" else "evaluator",
    context_hash = strrep("d", 64L),
    operation = "formal-glm-phase19-schedule-v1",
    output_kind = "formal-glm-phase19-ring128-dp-bridge-share-v1",
    purpose = purpose, source_producer = "formal-glm.phase18-materializer",
    ring_bits = 128L, frac_bits = 0L, vector_len = 4L,
    threshold = "", chunk_bytes = 65536L,
    ttl_seconds = 10L, max_runtime_seconds = 120L,
    worker_heartbeat = 1, state = "running", stored = FALSE)), servers)
  conns <- stats::setNames(list(
    structure(list(), class = "mock"),
    structure(list(), class = "mock")), servers)
  aggregate <- function(conns, expr, ...) {
    expect_true(is.list(expr) && !is.call(expr))
    stats::setNames(lapply(servers, function(server) list(
      capability_id = "exact_gc_v1", peer_id = peers[[server]],
      state = "complete", stored = TRUE, inbound_size = 0,
      outbound = NULL, worker_heartbeat = 2)), servers)
  }

  result <- .dsvert_exact_gc_run(
    conns, server_names = servers, servers = 1:2,
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = "op_33333333333333333333333333333333",
    source_key = "exact_gc_in_33333333333333333333333333333333",
    output_key = "exact_gc_out_33333333333333333333333333333333",
    operation = "formal-glm-phase19-schedule-v1", ring = 128L,
    frac_bits = 0L, vector_len = 4L, purpose = purpose,
    transport_ready = TRUE, initialized = initialized,
    timeout_seconds = 1, .aggregate = aggregate)

  expect_identical(result$operation, "formal-glm-phase19-schedule-v1")
  expect_false(any(grepl(
    "share|payload|release|result", names(result), ignore.case = TRUE)))
  contract <- .dsvert_exact_gc_capability_contract()
  expect_false("formal-glm-phase19-schedule-v1" %in% contract$operations)
  expect_false("formal-glm-phase19-schedule-v1" %in% contract$core_operations)

  expect_error(.dsvert_exact_gc_run(
    conns, server_names = servers, servers = 1:2,
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = "op_44444444444444444444444444444444",
    source_key = "exact_gc_in_44444444444444444444444444444444",
    output_key = "exact_gc_out_44444444444444444444444444444444",
    operation = "formal-glm-phase19-schedule-v1", ring = 127L,
    frac_bits = 0L, vector_len = 4L, purpose = purpose,
    transport_ready = TRUE, initialized = initialized,
    timeout_seconds = 1, .aggregate = aggregate),
    "formal-GLM schedule shape")
  expect_error(.dsvert_exact_gc_run(
    conns, server_names = servers, servers = 1:2,
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = "op_55555555555555555555555555555555",
    source_key = "exact_gc_in_55555555555555555555555555555555",
    output_key = "exact_gc_out_55555555555555555555555555555555",
    operation = "formal-glm-phase19-schedule-v1", ring = 128L,
    frac_bits = 0L, vector_len = 5L, purpose = purpose,
    transport_ready = TRUE, initialized = initialized,
    timeout_seconds = 1, .aggregate = aggregate),
    "formal-GLM schedule shape")
})

test_that("client multiprecision shape policy matches the server through Ring4096", {
  expect_identical(vapply(
    c(63L, 64L, 65L, 128L, 129L, 256L, 257L, 512L, 513L,
      1024L, 1025L, 2048L, 2049L, 4096L),
    .dsvert_exact_gc_container_bits, integer(1L)),
    c(64L, 64L, 128L, 128L, 256L, 256L, 512L, 512L, 1024L,
      1024L, 2048L, 2048L, 4096L, 4096L))
  expect_identical(vapply(
    c(63L, 512L, 513L, 1024L, 1025L, 2048L, 2049L, 4096L),
    .dsvert_exact_gc_direct_mul_max_chunk, integer(1L)),
    c(64L, 64L, 16L, 16L, 4L, 4L, 1L, 1L))
  expect_error(.dsvert_exact_gc_container_bits(4097L), "Invalid ring")
})

test_that("Ring4096 exact truncation reaches the pump only within its proved shape", {
  servers <- c("site_a", "site_b")
  peer_a <- paste0("dsv1_", strrep("a", 64L))
  peer_b <- paste0("dsv1_", strrep("b", 64L))
  peers <- c(site_a = peer_a, site_b = peer_b)
  conns <- stats::setNames(list(
    structure(list(), class = "mock"),
    structure(list(), class = "mock")), servers)
  initialized <- stats::setNames(lapply(servers, function(server) list(
    capability_id = "exact_gc_v1", peer_id = peers[[server]],
    peer_peer_id = peers[[setdiff(servers, server)]],
    role = if (server == "site_a") "garbler" else "evaluator",
    context_hash = strrep("c", 64L), operation = "truncate-floor",
    output_kind = "ring-share", purpose = "test.ring4096",
    source_producer = "test.producer", ring_bits = 4096L, frac_bits = 1L,
    vector_len = 42L, threshold = "", chunk_bytes = 65536L,
    ttl_seconds = 10L, max_runtime_seconds = 120L,
    worker_heartbeat = 1, state = "running", stored = FALSE)), servers)
  aggregate <- function(conns, expr, ...) {
    expect_true(is.list(expr) && !is.call(expr))
    stats::setNames(lapply(servers, function(server) list(
      capability_id = "exact_gc_v1", peer_id = peers[[server]],
      state = "complete", stored = TRUE, inbound_size = 0,
      outbound = NULL, worker_heartbeat = 2)), servers)
  }
  result <- .dsvert_exact_gc_run(
    conns, server_names = servers, servers = 1:2,
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = "op_33333333333333333333333333333333",
    source_key = "exact_gc_in_33333333333333333333333333333333",
    output_key = "exact_gc_out_33333333333333333333333333333333",
    operation = "truncate-floor", ring = 4096L, frac_bits = 1L,
    vector_len = 42L, purpose = "test.ring4096", transport_ready = TRUE,
    timeout_seconds = 1, initialized = initialized, .aggregate = aggregate)
  expect_identical(result$ring_bits, 4096L)
  oversized <- tryCatch(.dsvert_exact_gc_run(
    conns, server_names = servers, servers = 1:2,
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = "op_44444444444444444444444444444444",
    source_key = "exact_gc_in_44444444444444444444444444444444",
    output_key = "exact_gc_out_44444444444444444444444444444444",
    operation = "truncate-floor", ring = 4096L, frac_bits = 1L,
    vector_len = 43L, purpose = "test.ring4096", transport_ready = TRUE,
    timeout_seconds = 1, initialized = initialized, .aggregate = aggregate),
    error = identity)
  expect_s3_class(oversized, "dsvert_resource_oversize")
  expect_identical(oversized$code, "resource_oversize")
  expect_false(oversized$retryable)
})

test_that("exact operation contexts have canonical CSPRNG-derived capabilities", {
  context <- .dsvert_exact_gc_new_context(function(n) as.raw(seq_len(n)))
  expect_match(context$operation_id, "^op_[0-9a-f]{32}$")
  expect_match(context$source_key, "^exact_gc_in_[0-9a-f]{32}$")
  expect_match(context$output_key, "^exact_gc_out_[0-9a-f]{32}$")
  expect_error(.dsvert_exact_gc_new_context(function(n) raw(15L)),
               "16 secure random bytes")
})

test_that("joint-DP vector client has no output-share reconstruction surface", {
  expect_false(exists(
    ".dsvert_joint_dp_count_reconstruct", mode = "function",
    inherits = TRUE))
  expect_false(exists(
    ".dsvert_joint_dp_count_decode_payload", mode = "function",
    inherits = TRUE))
  body_text <- paste(c(
    deparse(body(.dsvert_joint_dp_vector_capsule)),
    deparse(body(.dsvert_joint_dp_vector_capsule_once))), collapse = "\n")
  expect_false(grepl("share_raw|payload_b64|base64_dec", body_text))
  expect_match(
    paste(deparse(body(.dsvert_joint_dp_vector_capsule)), collapse = "\n"),
    ".dsvert_joint_dp_vector_capsule_once", fixed = TRUE)
  expect_match(body_text, "dsvertJointDPVectorReleaseDS")
})

test_that("exact vecmul chunk identities cover datasets beyond one circuit", {
  batch <- "op_99999999999999999999999999999999"
  policy <- strrep("a", 64L)
  plan <- strrep("b", 64L)
  operations <- vapply(1:2, function(index) {
    .dsvert_exact_gc_vecmul_chunk_operation(
      batch, index, 2L, policy, plan)
  }, character(1L))
  expect_length(unique(operations), 2L)
  expect_true(all(grepl("^op_[0-9a-f]{32}$", operations)))
  expect_identical(
    .dsvert_exact_gc_vecmul_chunk_operation(
      batch, 1L, 2L, policy, plan),
    operations[[1L]])
  expect_match(
    .dsvert_exact_gc_vecmul_keys(batch)$destination,
    "^k2_exact_vecmul_z_[0-9a-f]{32}$")
  total_n <- 4097L
  chunk_size <- 256L
  lengths <- vapply(seq_len(ceiling(total_n / chunk_size)), function(index) {
    min(chunk_size, total_n - (index - 1L) * chunk_size)
  }, numeric(1L))
  expect_length(lengths, 17L)
  expect_true(all(lengths[1:16] == 256L))
  expect_identical(lengths[[17L]], 1)
})

test_that("direct exact vecmul uses checked chunks and no Beaver relay", {
  servers <- c("site_a", "site_b")
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  commands <- character(0)
  exact_lengths <- integer(0)
  validity_receipts <- list()
  policy <- strrep("a", 64L)
  plan <- strrep("c", 64L)
  aggregate <- function(conns, expr, ...) {
    per_site <- is.list(expr) && !is.call(expr)
    expressions <- if (per_site) expr[names(conns)] else
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    command_names <- vapply(
      expressions, function(value) as.character(value[[1L]]), character(1L))
    if (per_site) {
      expect_length(unique(command_names), 1L)
      validity_receipts[[length(validity_receipts) + 1L]] <<- expressions
    }
    command <- command_names[[1L]]
    commands <<- c(commands, command)
    stats::setNames(lapply(names(conns), function(server) {
      server_expr <- expressions[[server]]
      switch(command_names[[server]],
        exactGCVecmulBindInputsDS = list(
          capability_id = "exact_gc_v1", state = "bound", stored = TRUE,
          context_hash = strrep("b", 64L), policy_id = policy,
          plan_id = plan, backend = "ring127-ot",
          bound_x = "9223372036854775807",
          bound_y = "9223372036854775807",
          ring_bits = 127L, frac_bits = 50L, max_chunk = 256L),
        exactGCVecmulStartDS = list(
          capability_id = "exact_gc_v1", state = "running"),
        exactGCVecmulValidityDS = list(
          capability_id = "exact_gc_v1", state = "sealed",
          peer_blob = if (identical(server, "site_a")) "QQ" else "Qg"),
        exactGCVecmulValidityReceiveDS = list(
          capability_id = "exact_gc_v1", state = "checked", stored = TRUE),
        exactGCVecmulCommitDS = list(
          capability_id = "exact_gc_v1",
          state = if (as.integer(server_expr$chunk_index) ==
                       as.integer(server_expr$chunk_count)) {
            "committed"
          } else {
            "partial"
          },
          stored = TRUE),
        stop("unexpected exact vecmul endpoint: ", command_names[[server]]))
    }), names(conns))
  }
  result <- testthat::with_mocked_bindings(
    .dsvert_exact_gc_vecmul_run(
      conns, server_names = servers, servers = 1:2,
      session_id = "12345678-1234-4234-9234-123456789abc",
      total_n = 4097L, x_key = "arbitrary_x", y_key = "arbitrary_y",
      output_key = "arbitrary_z", transport_ready = TRUE,
      .aggregate = aggregate),
    .dsvert_exact_gc_new_context = function(...) list(
      operation_id = "op_99999999999999999999999999999999"),
    .dsvert_exact_gc_run = function(..., vector_len) {
      exact_lengths <<- c(exact_lengths, as.integer(vector_len))
      invisible(list(capability_id = "exact_gc_v1"))
    },
    .package = "dsVertClient")
  expect_identical(result$destination_key, "arbitrary_z")
  expect_length(exact_lengths, 17L)
  expect_true(all(exact_lengths[1:16] == 256L))
  expect_identical(exact_lengths[[17L]], 1L)
  expect_identical(commands[[1L]], "exactGCVecmulBindInputsDS")
  expect_identical(sum(commands == "exactGCVecmulBindInputsDS"), 1L)
  expect_identical(sum(commands == "exactGCVecmulStartDS"), 17L)
  expect_identical(sum(commands == "exactGCVecmulValidityReceiveDS"), 17L)
  expect_length(validity_receipts, 17L)
  expect_true(all(vapply(validity_receipts, function(requests) {
    identical(requests$site_a$peer_blob, "Qg") &&
      identical(requests$site_b$peer_blob, "QQ")
  }, logical(1L))))
  expect_identical(sum(commands == "exactGCAbortDS"), 0L)
  expect_false(any(commands %in% c(
    "exactGCVecmulPrepDS", "exactGCVecmulR1DS", "exactGCVecmulR2StageDS",
    "exactGCVecmulConsumeDS", "exactGCInitDS")))
  expect_false(any(commands %in% c(
    "k2BeaverVecmulR1DS", "k2BeaverVecmulR2DS", "mpcStoreBlobDS")))
})

test_that("exact vecmul rejects a partial validity fan-out before commit", {
  withr::local_options(list(dsvert.dsi.retry_deadline_seconds = 0.001))
  servers <- c("site_a", "site_b")
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  policy <- strrep("a", 64L)
  plan <- strrep("c", 64L)
  receive_calls <- 0L
  commit_seen <- FALSE
  aggregate <- function(conns, expr, ...) {
    if (is.list(expr) && !is.call(expr)) {
      receive_calls <<- receive_calls + 1L
      commands <- vapply(expr[names(conns)], function(value) {
        as.character(value[[1L]])
      }, character(1L))
      expect_true(all(commands == "exactGCVecmulValidityReceiveDS"))
      return(list(
        site_a = list(capability_id = "exact_gc_v1", state = "checked",
                      stored = TRUE),
        site_b = NULL))
    }
    command <- as.character(expr[[1L]])
    values <- lapply(names(conns), function(server) switch(command,
      exactGCVecmulBindInputsDS = list(
        capability_id = "exact_gc_v1", state = "bound", stored = TRUE,
        context_hash = strrep("b", 64L), policy_id = policy,
        plan_id = plan, backend = "ring127-ot",
        bound_x = "9223372036854775807",
        bound_y = "9223372036854775807",
        ring_bits = 127L, frac_bits = 50L, max_chunk = 256L),
      exactGCVecmulStartDS = list(
        capability_id = "exact_gc_v1", state = "running"),
      exactGCVecmulValidityDS = list(
        capability_id = "exact_gc_v1", state = "sealed",
        peer_blob = if (identical(server, "site_a")) "QQ" else "Qg"),
      exactGCVecmulCommitDS = {
        commit_seen <<- TRUE
        list(capability_id = "exact_gc_v1", state = "committed",
             stored = TRUE)
      },
      exactGCAbortDS = TRUE,
      stop("unexpected exact vecmul endpoint: ", command)))
    stats::setNames(values, names(conns))
  }
  expect_error(testthat::with_mocked_bindings(
    .dsvert_exact_gc_vecmul_run(
      conns, server_names = servers, servers = 1:2,
      session_id = "12345678-1234-4234-9234-123456789abc",
      total_n = 1L, x_key = "arbitrary_x", y_key = "arbitrary_y",
      output_key = "arbitrary_z", transport_ready = TRUE,
      .aggregate = aggregate),
    .dsvert_exact_gc_new_context = function(...) list(
      operation_id = "op_99999999999999999999999999999999"),
    .dsvert_exact_gc_run = function(...) invisible(list()),
    .package = "dsVertClient"),
    "remained unavailable until the retry deadline")
  expect_gte(receive_calls, 1L)
  expect_false(commit_seen)
})

test_that("producer manifests are claimed per peer without remote slot binding", {
  servers <- c("site_a", "site_b")
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  policy <- strrep("a", 64L)
  plan <- strrep("c", 64L)
  handles <- c(site_a = paste0("A", strrep("a", 42L)),
               site_b = paste0("B", strrep("b", 42L)))
  manifests <- lapply(servers, function(server) list(
    manifest_handle = unname(handles[[server]]), total_n = 1L))
  names(manifests) <- servers
  commands <- character(0)
  claimed <- NULL
  aggregate <- function(conns, expr, ...) {
    if (is.list(expr) && !is.call(expr)) {
      command_names <- vapply(
        expr[names(conns)], function(value) as.character(value[[1L]]),
        character(1L))
      expect_true(all(command_names == "exactGCVecmulValidityReceiveDS"))
      commands <<- c(commands, "exactGCVecmulValidityReceiveDS")
      return(stats::setNames(lapply(names(conns), function(...) list(
        capability_id = "exact_gc_v1", state = "checked", stored = TRUE)),
        names(conns)))
    }
    command <- as.character(expr[[1L]])
    commands <<- c(commands, command)
    value <- switch(command,
      exactGCVecmulStartDS = list(
        capability_id = "exact_gc_v1", state = "running"),
      exactGCVecmulValidityDS = list(
        capability_id = "exact_gc_v1", state = "sealed", peer_blob = "QQ"),
      exactGCVecmulValidityReceiveDS = list(
        capability_id = "exact_gc_v1", state = "checked", stored = TRUE),
      exactGCVecmulCommitDS = list(
        capability_id = "exact_gc_v1", state = "committed", stored = TRUE),
      stop("unexpected producer-manifest endpoint: ", command))
    stats::setNames(rep(list(value), length(conns)), names(conns))
  }
  fanout <- function(conns, requests, ...) {
    claimed <<- requests
    stats::setNames(lapply(names(conns), function(server) list(
      capability_id = "exact_gc_v1", state = "claimed", stored = TRUE,
      context_hash = strrep("b", 64L), policy_id = policy,
      plan_id = plan, backend = "ring127-ot",
      bound_x = "4503599627370496", bound_y = "2251799813685248",
      ring_bits = 127L, frac_bits = 50L, max_chunk = 256L)), names(conns))
  }
  result <- testthat::with_mocked_bindings(
    .dsvert_exact_gc_vecmul_run(
      conns, server_names = servers, servers = 1:2,
      session_id = "12345678-1234-4234-9234-123456789abc",
      total_n = 1L, output_key = "producer_output",
      input_manifests = manifests, transport_ready = TRUE,
      .aggregate = aggregate),
    .dsvert_exact_gc_new_context = function(...) list(
      operation_id = "op_99999999999999999999999999999999"),
    .dsvert_exact_gc_run = function(...) invisible(list()),
    .dsvert_fanout_by_site = fanout,
    .package = "dsVertClient")
  expect_identical(result$destination_key, "producer_output")
  expect_identical(names(claimed), servers)
  expect_identical(vapply(claimed, function(request) {
    as.character(request$manifest_handle)
  }, character(1L)), handles)
  expect_true(all(vapply(claimed, function(request) {
    identical(as.character(request[[1L]]), "exactGCVecmulClaimInputsDS")
  }, logical(1L))))
  expect_false("exactGCVecmulBindInputsDS" %in% commands)
})

test_that("Ring127 GLM initializes pinned peers before producer minting", {
  servers <- c("site_a", "site_b")
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  events <- character(0)
  exact_call <- NULL
  manifest <- function(prefix) list(
    manifest_handle = paste0(prefix, strrep(tolower(prefix), 42L)),
    total_n = 3L)
  aggregate <- function(conns, expr, ...) {
    command <- as.character(expr[[1L]])
    events <<- c(events, command)
    if (identical(command, "k2PrepareWeightedResidualShareDS")) {
      values <- list(
        site_a = list(stored = TRUE,
                      exact_vecmul_manifest = manifest("A")),
        site_b = list(stored = TRUE,
                      exact_vecmul_manifest = manifest("B")))
      return(values[names(conns)])
    }
    if (identical(command, "k2FinalizeWeightedResidualShareDS")) {
      return(stats::setNames(as.list(rep(TRUE, length(conns))),
                             names(conns)))
    }
    stop("unexpected GLM exact producer endpoint: ", command)
  }
  testthat::with_mocked_bindings(
    .glm_apply_shared_weight_residual(
      datasources = conns, dcf_parties = servers, dcf_conns = 1:2,
      transport_pks = list(),
      session_id = "12345678-1234-4234-9234-123456789abc",
      n_obs = 3L, .dsAgg = aggregate,
      .sendBlob = function(...) stop("blob relay reached"), ring = 127L),
    .dsvert_setup_exact_gc_transport = function(...) {
      events <<- c(events, "setup-exact")
      invisible(list())
    },
    .dsvert_exact_gc_vecmul_run = function(...) {
      exact_call <<- list(...)
      events <<- c(events, "run-exact")
      invisible(list())
    },
    .package = "dsVertClient")
  expect_lt(match("setup-exact", events),
            match("k2PrepareWeightedResidualShareDS", events))
  expect_true(isTRUE(exact_call$transport_ready))
  expect_identical(names(exact_call$input_manifests), servers)
  expect_identical(
    vapply(exact_call$input_manifests, `[[`, character(1L),
           "manifest_handle"),
    c(site_a = paste0("A", strrep("a", 42L)),
      site_b = paste0("B", strrep("b", 42L))))
})

test_that("failed checked vecmul aborts both its chunk and retained batch", {
  servers <- c("site_a", "site_b")
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  policy <- strrep("a", 64L)
  plan <- strrep("c", 64L)
  batch <- "op_99999999999999999999999999999999"
  aborted <- character(0)
  aggregate <- function(conns, expr, ...) {
    selected <- conns
    command <- as.character(expr[[1L]])
    if (identical(command, "exactGCAbortDS")) {
      aborted <<- c(aborted, as.character(expr$operation_id))
      return(stats::setNames(as.list(rep(TRUE, length(selected))),
                             names(selected)))
    }
    value <- switch(command,
      exactGCVecmulBindInputsDS = list(
        capability_id = "exact_gc_v1", state = "bound", stored = TRUE,
        context_hash = strrep("b", 64L), policy_id = policy,
        plan_id = plan, backend = "ring127-ot",
        bound_x = "9223372036854775807",
        bound_y = "9223372036854775807",
        ring_bits = 127L, frac_bits = 50L, max_chunk = 256L),
      exactGCVecmulStartDS = list(
        capability_id = "exact_gc_v1", state = "running"),
      exactGCVecmulValidityDS = list(
        capability_id = "exact_gc_v1", state = "sealed",
        peer_blob = "not*base64url"),
      stop("unexpected command: ", command))
    stats::setNames(rep(list(value), length(selected)), names(selected))
  }
  expect_error(testthat::with_mocked_bindings(
    .dsvert_exact_gc_vecmul_run(
      conns, server_names = servers, servers = 1:2,
      session_id = "12345678-1234-4234-9234-123456789abc",
      total_n = 1L, transport_ready = TRUE, .aggregate = aggregate),
    .dsvert_exact_gc_new_context = function(...) list(operation_id = batch),
    .dsvert_exact_gc_run = function(...) invisible(list()),
    .package = "dsVertClient"), "Exact MPC multiplication failed")
  chunk <- .dsvert_exact_gc_vecmul_chunk_operation(
    batch, 1L, 1L, policy, plan)
  expect_setequal(unique(aborted), c(batch, chunk))
})

test_that("Ring127 common helper cannot fall back to legacy truncation", {
  called <- NULL
  copied <- NULL
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")),
                           c("site_a", "site_b"))
  testthat::with_mocked_bindings(
    .ring127_vecmul(
      "input_x", "input_y", "output_z", 17L,
      datasources = conns, dealer_ci = 2L,
      server_list = c("site_a", "site_b"),
      server_names = c("site_a", "site_b"), y_server = "site_a",
      nl = "site_b", transport_pks = list(),
      session_id = "12345678-1234-4234-9234-123456789abc",
      .dsAgg = function(...) NULL,
      .sendBlob = function(...) stop("legacy blob relay reached")),
    .ring127_invocation_tag = function() strrep("a", 32L),
    .ring127_local_fanout = function(datasources, server_list, server_names,
                                     y_server, .dsAgg, make_call) {
      copied <<- make_call("site_a", TRUE)
      invisible(NULL)
    },
    .dsvert_exact_gc_vecmul_run = function(...) {
      called <<- list(...)
      invisible(list())
    },
    .package = "dsVertClient")
  expect_identical(called$x_key, "input_x")
  expect_identical(called$y_key, "input_y")
  expect_identical(called$output_key,
                   paste0("__r127_mul_", strrep("a", 32L)))
  expect_identical(called$total_n, 17L)
  expect_false(called$transport_ready)
  expect_identical(copied$a_key, called$output_key)
  expect_identical(copied$output_key, "output_z")
  expect_identical(as.character(copied[[1L]]),
                   "k2Ring127AffineCombineDS")
})

test_that("fresh Ring127 destinations avoid the overwrite-copy round", {
  called <- NULL
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")),
                           c("site_a", "site_b"))
  testthat::with_mocked_bindings(
    .ring127_vecmul(
      "input_x", "input_y", "fresh_output", 17L,
      datasources = conns, dealer_ci = 2L,
      server_list = c("site_a", "site_b"),
      server_names = c("site_a", "site_b"), y_server = "site_a",
      nl = "site_b", transport_pks = list(),
      session_id = "12345678-1234-4234-9234-123456789abc",
      .dsAgg = function(...) stop("fresh output must not be copied"),
      .sendBlob = function(...) NULL, destination_fresh = TRUE),
    .dsvert_exact_gc_vecmul_run = function(...) {
      called <<- list(...)
      invisible(list())
    },
    .package = "dsVertClient")
  expect_identical(called$output_key, "fresh_output")
})

test_that("Ring127 exact transport readiness is explicit and operation-local", {
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")),
                           c("site_a", "site_b"))
  setup <- NULL
  phase_calls <- 0L
  phase <- function(conns, expr) {
    phase_calls <<- phase_calls + 1L
    stats::setNames(rep(list(TRUE), length(conns)), names(conns))
  }
  ready <- testthat::with_mocked_bindings(
    .ring127_transport_once(
      FALSE, conns, names(conns), names(conns),
      "12345678-1234-4234-9234-123456789abc", phase),
    .dsvert_setup_exact_gc_transport = function(...) {
      setup <<- list(...)
      invisible(list())
    },
    .package = "dsVertClient")
  expect_true(ready)
  expect_identical(setup$servers, 1:2)
  expect_identical(
    setup$.aggregate(
      conns, call("phaseDS"), async = TRUE, errors.print = FALSE,
      error = function(...) NULL),
    list(site_a = TRUE, site_b = TRUE))
  expect_identical(phase_calls, 1L)

  expect_true(testthat::with_mocked_bindings(
    .ring127_transport_once(
      TRUE, conns, names(conns), names(conns),
      "12345678-1234-4234-9234-123456789abc", function(...) NULL),
    .dsvert_setup_exact_gc_transport = function(...) {
      stop("ready transport was initialized again")
    },
    .package = "dsVertClient"))
})

test_that("nested Ring127 primitive propagates one ready transport", {
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")),
                           c("site_a", "site_b"))
  setup_calls <- 0L
  ready_flags <- logical(0)
  testthat::with_mocked_bindings(
    .ring127_exp_round_keyed_extended(
      "input", "output", 5L, conns, 1L, names(conns), names(conns),
      "site_a", "site_b", list(),
      "12345678-1234-4234-9234-123456789abc",
      function(...) NULL, function(...) NULL),
    .ring127_transport_once = function(transport_ready, ...) {
      setup_calls <<- setup_calls + 1L
      TRUE
    },
    .ring127_get_half_fp = function() "AQ==",
    .ring127_exact_public_scale = function(..., transport_ready = FALSE) {
      ready_flags <<- c(ready_flags, transport_ready)
      invisible(NULL)
    },
    .ring127_exp_round_keyed = function(..., transport_ready = FALSE) {
      ready_flags <<- c(ready_flags, transport_ready)
      invisible(NULL)
    },
    .ring127_vecmul = function(..., transport_ready = FALSE) {
      ready_flags <<- c(ready_flags, transport_ready)
      invisible(NULL)
    },
    .package = "dsVertClient")
  expect_identical(setup_calls, 1L)
  expect_identical(ready_flags, rep(TRUE, 3L))
})

test_that("Ring127 Clenshaw products are one-shot across stages and calls", {
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")),
                           c("site_a", "site_b"))
  previous <- .ring127_exp_coef_cache$coef_res
  on.exit({
    if (is.null(previous)) {
      rm(list = "coef_res", envir = .ring127_exp_coef_cache)
    } else {
      .ring127_exp_coef_cache$coef_res <- previous
    }
  }, add = TRUE)
  .ring127_exp_coef_cache$coef_res <- list(
    one_over_a = jsonlite::base64_enc(raw(16L)),
    degree = 3L,
    coeffs = jsonlite::base64_enc(raw(4L * 16L)))
  exact_destinations <- character(0)
  fresh_flags <- logical(0)
  tag_index <- 0L
  run_once <- function() .ring127_exp_round_keyed(
    "input", "reusable_output", 5L, conns, 2L, names(conns),
    names(conns), "site_a", "site_b", list(),
    "12345678-1234-4234-9234-123456789abc",
    function(...) NULL, function(...) NULL)
  testthat::with_mocked_bindings({
    run_once()
    run_once()
  },
    .ring127_invocation_tag = function() {
      tag_index <<- tag_index + 1L
      sprintf("%032x", tag_index)
    },
    .ring127_transport_once = function(...) TRUE,
    .ring127_local_fanout = function(...) invisible(NULL),
    .ring127_exact_public_scale = function(
        in_key, scalar_fp, output_key, ...,
        destination_fresh = FALSE) {
      exact_destinations <<- c(exact_destinations, output_key)
      fresh_flags <<- c(fresh_flags, destination_fresh)
      invisible(output_key)
    },
    .ring127_vecmul = function(
        x_key, y_key, output_key, ...,
        destination_fresh = FALSE) {
      exact_destinations <<- c(exact_destinations, output_key)
      fresh_flags <<- c(fresh_flags, destination_fresh)
      invisible(output_key)
    },
    .package = "dsVertClient")
  expect_length(exact_destinations, 8L)
  expect_identical(anyDuplicated(exact_destinations), 0L)
  expect_true(all(fresh_flags))
  expect_false("reusable_output" %in% exact_destinations)
})

test_that("GLM softplus stages claim producer manifests without legacy Bind", {
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")),
                           c("site_a", "site_b"))
  prepared_calls <- NULL
  exact_call <- NULL
  testthat::with_mocked_bindings(
    .ring127_glm_softplus_manifest_mul(
      stage = 9L, invocation_id = strrep("a", 32L),
      output_key = paste0("__r127_spmul09_", strrep("a", 32L)),
      n = 23L, datasources = conns, server_list = names(conns),
      server_names = names(conns),
      session_id = "12345678-1234-4234-9234-123456789abc",
      .dsAgg = function(...) stop("raw aggregate should be wrapped")),
    .dsvert_fanout_by_site = function(conns, expressions, ...) {
      prepared_calls <<- expressions
      stats::setNames(lapply(names(conns), function(server) list(
        manifest_handle = strrep(if (server == "site_a") "A" else "B",
                                 43L),
        total_n = 23L)), names(conns))
    },
    .dsvert_exact_gc_vecmul_run = function(...) {
      exact_call <<- list(...)
      invisible(list())
    },
    .package = "dsVertClient")
  expect_identical(unname(vapply(prepared_calls, function(expr) {
    as.character(expr[[1L]])
  }, character(1L))), rep("exactGCGLMSoftplusPrepareDS", 2L))
  expect_true(all(vapply(prepared_calls, function(expr) {
    identical(expr$stage, 9L)
  }, logical(1L))))
  expect_true(all(vapply(prepared_calls, function(expr) {
    identical(expr$invocation_id, strrep("a", 32L))
  }, logical(1L))))
  expect_identical(exact_call$output_key,
                   paste0("__r127_spmul09_", strrep("a", 32L)))
  expect_identical(exact_call$total_n, 23L)
  expect_true(exact_call$transport_ready)
  expect_length(exact_call$input_manifests, 2L)
  expect_false(any(vapply(prepared_calls, function(expr) {
    identical(as.character(expr[[1L]]), "exactGCVecmulBindInputsDS")
  }, logical(1L))))
})

test_that("GLM softplus uses all fixed producer stages and no slot Bind", {
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")),
                           c("site_a", "site_b"))
  previous <- .ring127_softplus_coef_cache$coef_res
  on.exit({
    if (is.null(previous)) {
      rm(list = "coef_res", envir = .ring127_softplus_coef_cache)
    } else {
      .ring127_softplus_coef_cache$coef_res <- previous
    }
  }, add = TRUE)
  .ring127_softplus_coef_cache$coef_res <- list(
    one_over_a = jsonlite::base64_enc(raw(16L)),
    degree = 36L,
    coeffs = jsonlite::base64_enc(raw(37L * 16L)))
  stages <- integer(0)
  outputs <- character(0)
  invocations <- character(0)
  run_once <- function() .ring127_softplus_round_keyed(
      in_key = "k2_eta_share_fp", out_key = "softplus_share_fp", n = 5L,
      datasources = conns, dealer_ci = 2L, server_list = names(conns),
      server_names = names(conns), y_server = "site_a", nl = "site_b",
      transport_pks = list(),
      session_id = "12345678-1234-4234-9234-123456789abc",
      .dsAgg = function(...) NULL, .sendBlob = function(...) NULL,
      producer_bound_glm = TRUE)
  tag_index <- 0L
  testthat::with_mocked_bindings({
    run_once()
    run_once()
  },
    .ring127_invocation_tag = function() {
      tag_index <<- tag_index + 1L
      strrep(if (tag_index == 1L) "a" else "b", 32L)
    },
    .ring127_transport_once = function(...) TRUE,
    .ring127_local_fanout = function(...) invisible(NULL),
    .ring127_glm_softplus_manifest_mul = function(
        stage, invocation_id, output_key, ...) {
      stages <<- c(stages, stage)
      invocations <<- c(invocations, invocation_id)
      outputs <<- c(outputs, output_key)
      invisible(output_key)
    },
    .ring127_vecmul = function(...) {
      stop("legacy Ring127 slot Bind was reached")
    },
    .package = "dsVertClient")
  expect_identical(stages, rep(0:36, 2L))
  expect_identical(invocations,
                   c(rep(strrep("a", 32L), 37L),
                     rep(strrep("b", 32L), 37L)))
  expect_identical(anyDuplicated(outputs), 0L)
  expect_identical(outputs[[1L]],
                   paste0("__r127_spy_", strrep("a", 32L)))
  expect_identical(outputs[[37L]],
                   paste0("__r127_spmul36_", strrep("a", 32L)))
  expect_identical(outputs[[38L]],
                   paste0("__r127_spy_", strrep("b", 32L)))
  expect_match(paste(deparse(body(.k2_strict_loop)), collapse = "\n"),
               "producer_bound_glm = TRUE", fixed = TRUE)
  expect_match(paste(deparse(body(.k3_ring63_gradient_loop)), collapse = "\n"),
               "producer_bound_glm = TRUE", fixed = TRUE)
})

test_that("public Ring127 scaling uses additive shares and exact truncation", {
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")),
                           c("site_a", "site_b"))
  commands <- character(0)
  constant_keys <- character(0)
  constant_owners <- logical(0)
  aggregate_calls <- 0L
  selected_names <- NULL
  exact_call <- NULL
  aggregate <- function(selected, expr) {
    aggregate_calls <<- aggregate_calls + 1L
    selected_names <<- names(selected)
    commands <<- c(commands, unname(vapply(
      expr, function(site_expr) as.character(site_expr[[1L]]), character(1L))))
    constant_keys <<- c(constant_keys, unname(vapply(
      expr, function(site_expr) as.character(site_expr$output_key), character(1L))))
    constant_owners <<- c(constant_owners, unname(vapply(
      expr, function(site_expr) isTRUE(site_expr$is_party0), logical(1L))))
    invisible(stats::setNames(rep(list(TRUE), length(expr)), names(expr)))
  }

  testthat::with_mocked_bindings(
    .ring127_exact_public_scale(
      in_key = "private_x", scalar_fp = "AQ", output_key = "scaled_x",
      n = 11L, datasources = conns, dealer_ci = 2L,
      server_list = c("site_a", "site_b"),
      server_names = c("site_a", "site_b"), y_server = "site_a",
      nl = "site_b", transport_pks = list(),
      session_id = "12345678-1234-4234-9234-123456789abc",
      .dsAgg = aggregate, .sendBlob = function(...) NULL),
    .ring127_vecmul = function(x_key, y_key, output_key, n, ...) {
      exact_call <<- list(x_key = x_key, y_key = y_key,
                         output_key = output_key, n = n)
      invisible(NULL)
    },
    .package = "dsVertClient")

  expect_identical(commands,
                   rep("k2Ring127AffineCombineDS", length(conns)))
  expect_identical(aggregate_calls, 1L)
  expect_identical(selected_names, names(conns))
  expect_false(any(commands == "k2Ring127LocalScaleDS"))
  expect_length(unique(constant_keys), 1L)
  expect_identical(constant_owners, c(TRUE, FALSE))
  expect_identical(exact_call$x_key, "private_x")
  expect_identical(exact_call$y_key, constant_keys[[1L]])
  expect_identical(exact_call$output_key, "scaled_x")
  expect_identical(exact_call$n, 11L)
})

test_that("exact transport setup uses only exact-specific pinned endpoints", {
  servers <- c("site_a", "site_b")
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  commands <- character(0)
  encode <- function(byte, n) {
    raw <- as.raw(rep(byte, n))
    b64 <- gsub("[\r\n]", "", jsonlite::base64_enc(raw))
    gsub("=+$", "", gsub("/", "_", gsub("+", "-", b64, fixed = TRUE),
                            fixed = TRUE))
  }
  aggregate <- function(conns, expr, ...) {
    command <- as.character(expr[[1L]])
    commands <<- c(commands, command)
    if (identical(command, "exactGCTransportInitDS")) {
      return(stats::setNames(lapply(seq_along(conns), function(index) list(
        capability_id = "exact_gc_v1", transport_pk = encode(index, 32L),
        identity_pk = encode(index + 2L, 32L),
        signature = encode(index + 4L, 64L))), names(conns)))
    }
    if (identical(command, "exactGCBindPeersDS")) {
      expect_false("analysis_contract_b64" %in% names(expr))
      return(stats::setNames(lapply(seq_along(conns), function(index) list(
        capability_id = "exact_gc_v1", bound = TRUE)), names(conns)))
    }
    stop("unexpected endpoint")
  }
  result <- .dsvert_setup_exact_gc_transport(
    conns, servers, 1:2, "12345678-1234-4234-9234-123456789abc",
    .aggregate = aggregate)
  expect_identical(names(result), servers)
  expect_identical(commands,
                   c("exactGCTransportInitDS", "exactGCBindPeersDS"))
  expect_false(any(commands %in% c(
    "glmRing63TransportInitDS", "mpcStoreTransportKeysDS")))
})

test_that("analysis-bound setup selects authorities by full K identity pins", {
  contract <- .client_exact_gc_analysis_contract(3L)
  analysis <- .dsvert_exact_gc_analysis_binding(contract)
  server_names <- names(contract$execution$peer_pins)
  pins <- unlist(contract$execution$peer_pins, use.names = TRUE)
  authorities <- unlist(
    contract$semantic$noise_authorities, use.names = FALSE)
  authority_names <- names(pins)[match(authorities, unname(pins))]
  datasources <- stats::setNames(lapply(
    server_names, function(...) structure(list(), class = "mock")),
    server_names)
  selected <- match(authority_names, server_names)
  bound_call <- NULL
  encode <- function(byte, n) {
    .dsvert_exact_gc_b64url_encode(as.raw(rep(byte, n)))
  }
  aggregate <- function(conns, expr, ...) {
    command <- as.character(expr[[1L]])
    if (identical(command, "exactGCTransportInitDS")) {
      return(stats::setNames(lapply(names(conns), function(server) list(
        capability_id = "exact_gc_v1",
        transport_pk = encode(match(server, server_names) + 10L, 32L),
        identity_pk = unname(pins[[server]]),
        signature = encode(match(server, server_names) + 20L, 64L))),
        names(conns)))
    }
    if (identical(command, "exactGCBindPeersDS")) {
      bound_call <<- expr
      return(stats::setNames(lapply(names(conns), function(...) list(
        capability_id = "exact_gc_v1", bound = TRUE,
        analysis_binding = analysis$binding,
        analysis_binding_sha256 = analysis$sha256)), names(conns)))
    }
    stop("unexpected endpoint")
  }

  result <- .dsvert_setup_exact_gc_transport(
    datasources, server_names, selected,
    "12345678-1234-4234-9234-123456789abc",
    analysis_contract = contract, .aggregate = aggregate)
  expect_setequal(names(result), authority_names)
  expect_false("analysis_contract_b64" %in% names(bound_call))
  expect_identical(bound_call$artifact_key, contract$artifact_key)
  expect_identical(
    attr(result, "exact_gc_analysis_binding", exact = TRUE),
    analysis$binding)

  wrong <- contract
  names(wrong$execution$peer_pins)[[3L]] <- "other_site"
  expect_error(.dsvert_setup_exact_gc_transport(
    datasources, server_names, selected,
    "12345678-1234-4234-9234-123456789abc",
    analysis_contract = wrong, .aggregate = aggregate),
    "full K|peer names")
})

test_that("client accepts only analysis-bound scalar Count initialization", {
  contract <- .client_exact_gc_analysis_contract(3L)
  analysis <- .dsvert_exact_gc_analysis_binding(contract)
  servers <- names(contract$execution$peer_pins)[1:2]
  peer_ids <- vapply(
    unlist(contract$semantic$noise_authorities, use.names = FALSE),
    .dsvert_exact_gc_identity_peer_id, character(1L))
  purpose <- paste0("joint-dp-laplace-v2/", strrep("c", 64L))
  states <- stats::setNames(lapply(seq_along(servers), function(index) list(
    capability_id = "exact_gc_v1",
    peer_id = peer_ids[[index]],
    peer_peer_id = peer_ids[[3L - index]],
    role = if (identical(peer_ids[[index]], sort(peer_ids)[[1L]])) {
      "garbler"
    } else {
      "evaluator"
    },
    context_hash = strrep("d", 64L),
    operation = "joint-dp-laplace-v2",
    output_kind = "joint-dp-ring-share-v2",
    purpose = purpose,
    source_producer = "count.scalar.v1",
    ring_bits = 127L,
    frac_bits = 0L,
    vector_len = 1L,
    threshold = "",
    chunk_bytes = 65536L,
    ttl_seconds = 10L,
    max_runtime_seconds = 120L,
    worker_heartbeat = 1,
    analysis_binding_sha256 = analysis$sha256,
    state = "running",
    stored = FALSE)), servers)
  validated <- .dsvert_exact_gc_validate_init(
    states, servers, "joint-dp-laplace-v2", 127L, 0L, 1L, purpose,
    analysis_binding = analysis)
  expect_identical(validated$analysis_binding_sha256, analysis$sha256)

  conns <- stats::setNames(lapply(
    servers, function(...) structure(list(), class = "mock")), servers)
  aggregate <- function(conns, expr, ...) {
    expect_true(is.list(expr) && !is.call(expr))
    stats::setNames(lapply(names(conns), function(server) list(
      capability_id = "exact_gc_v1",
      peer_id = states[[server]]$peer_id,
      state = "complete", stored = TRUE,
      inbound_size = 0, outbound = NULL, worker_heartbeat = 2)),
      names(conns))
  }
  result <- .dsvert_exact_gc_run(
    conns, server_names = servers, servers = 1:2,
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = "op_33333333333333333333333333333333",
    source_key = "exact_gc_in_33333333333333333333333333333333",
    output_key = "exact_gc_out_33333333333333333333333333333333",
    operation = "joint-dp-laplace-v2", ring = 127L,
    frac_bits = 0L, vector_len = 1L, purpose = purpose,
    transport_ready = TRUE, initialized = states,
    analysis_contract = contract, timeout_seconds = 1,
    .aggregate = aggregate)
  expect_identical(result$analysis_binding_sha256, analysis$sha256)

  bad <- states
  bad[[1L]]$analysis_binding_sha256 <- strrep("e", 64L)
  expect_error(.dsvert_exact_gc_validate_init(
    bad, servers, "joint-dp-laplace-v2", 127L, 0L, 1L, purpose,
    analysis_binding = analysis), "analysis binding")
  expect_error(.dsvert_exact_gc_validate_init(
    states, servers, "joint-dp-laplace-v2", 127L, 0L, 1L, purpose),
    "analysis binding")
  expect_error(.dsvert_exact_gc_validate_init(
    states, servers, "joint-dp-laplace-v2", 128L, 0L, 1L, purpose,
    analysis_binding = analysis))
})

test_that("cross transport mints and consumes peer-specific cleanup capabilities", {
  servers <- c("site_a", "site_b")
  conns <- stats::setNames(list(structure(list(), class = "mock"),
                                structure(list(), class = "mock")), servers)
  session_id <- "12345678-1234-4234-9234-123456789abc"
  commands <- character()
  cleanup_calls <- list()
  encode <- function(byte, n) {
    raw <- as.raw(rep(byte, n))
    b64 <- gsub("[\r\n]", "", jsonlite::base64_enc(raw))
    gsub("=+$", "", gsub("/", "_", gsub("+", "-", b64, fixed = TRUE),
                            fixed = TRUE))
  }
  cleanup_capability <- function(index) {
    as.character(jsonlite::toJSON(list(
      version = .DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION,
      contract = list(
        version = .DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION,
        session_id = session_id,
        cleanup_purpose =
          .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
        operation_scope =
          "all_and_only_operations_in_bound_exact_session_v1",
        peer_binding_digest = strrep(as.character(index), 64L)),
      signature = strrep(if (index == 1L) "A" else "B", 86L)),
      auto_unbox = TRUE, null = "null", na = "null", digits = 17,
      pretty = FALSE))
  }
  capabilities <- stats::setNames(
    lapply(seq_along(servers), cleanup_capability), servers)
  aggregate <- function(conns, expr, ...) {
    expressions <- if (is.list(expr) && !is.call(expr)) {
      expr[names(conns)]
    } else {
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    }
    stats::setNames(lapply(names(conns), function(server) {
      expression <- expressions[[server]]
      command <- as.character(expression[[1L]])
      commands <<- c(commands, command)
      if (identical(command, "exactGCTransportInitDS")) {
        index <- match(server, servers)
        return(list(
          capability_id = "exact_gc_v1",
          transport_pk = encode(index, 32L),
          identity_pk = encode(index + 2L, 32L),
          signature = encode(index + 4L, 64L)))
      }
      if (identical(command, "exactGCBindPeersDS")) {
        expect_identical(
          expression$cleanup_purpose,
          .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE)
        return(list(
          capability_id = "exact_gc_v1", bound = TRUE,
          cleanup_purpose =
            .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
          cleanup_capability_json = capabilities[[server]]))
      }
      if (identical(command, "exactGCCleanupDS")) {
        cleanup_calls[[server]] <<- expression
        return(TRUE)
      }
      stop("unexpected endpoint: ", command)
    }), names(conns))
  }

  setup <- .dsvert_setup_exact_gc_transport(
    conns, servers, 1:2, session_id,
    cleanup_purpose = .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
    .aggregate = aggregate)
  expect_identical(
    attr(setup, "exact_gc_cleanup_capabilities", exact = TRUE),
    capabilities)
  expect_true(.dsvert_exact_gc_cleanup_best_effort(
    conns, session_id, setup, .aggregate = aggregate))
  expect_setequal(names(cleanup_calls), servers)
  for (server in servers) {
    expect_identical(
      cleanup_calls[[server]]$cleanup_capability_json,
      .dsvert_dsi_text_encode(capabilities[[server]]))
  }
  expect_identical(sum(commands == "exactGCCleanupDS"), 2L)
  expect_false(any(commands == "mpcCleanupDS"))

  malformed <- jsonlite::fromJSON(capabilities[[1L]], simplifyVector = FALSE)
  malformed$unexpected <- TRUE
  capabilities[[1L]] <- as.character(jsonlite::toJSON(
    malformed, auto_unbox = TRUE, null = "null", na = "null", digits = 17,
    pretty = FALSE))
  expect_error(.dsvert_setup_exact_gc_transport(
    conns, servers, 1:2, session_id,
    cleanup_purpose = .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
    .aggregate = aggregate), "invalid cleanup capability", fixed = TRUE)
})
