test_that("Count execution keeps analyst inputs out of the cryptographic contract", {
  expect_true(exists(".dsvert_dp_count_execute_v1", mode = "function"))
  expect_identical(
    names(formals(.dsvert_dp_count_execute_v1)),
    c(
      "data_name", "datasources", ".aggregate", ".prepare",
      ".new_context", ".run_exact", ".store_typed", ".cleanup",
      ".abort"))
  expect_false(any(c(
    "epsilon", "delta", "bound", "seed", "budget", "catalog",
    "authorities", "contract") %in%
      names(formals(.dsvert_dp_count_execute_v1))))
})

.count_exec_b64url <- function(value) {
  chartr("+/", "-_", sub(
    "=+$", "", gsub("[\r\n]", "", jsonlite::base64_enc(value)),
    perl = TRUE))
}

.count_exec_plan <- function() {
  list(
    version = "dsvert-joint-dp-laplace-plan-v2",
    sampler = "hkdf-sha256-aes128ctr-two-geometric-tv-v2",
    bernoulli_bits = 8L,
    stop_numerator = "51",
    max_geometric_steps = 76L,
    sensitivity_steps = "1",
    coordinate_count = 1L,
    epsilon_effective_upper_numerator = "1",
    epsilon_effective_upper_denominator = "1",
    implementation_delta_numerator = "1",
    implementation_delta_denominator = "1000000000",
    implementation_delta_bound = "1/1000000000",
    accounting = "exact-rational finite sampler certificate",
    bernoulli_trials = 608L,
    aes_blocks = 38L)
}

.count_exec_fixture_cache <- new.env(parent = emptyenv())

.count_exec_fixture <- function(k) {
  cache_key <- as.character(as.integer(k))
  if (exists(cache_key, envir = .count_exec_fixture_cache,
             inherits = FALSE)) {
    return(get(cache_key, envir = .count_exec_fixture_cache,
               inherits = FALSE))
  }
  private <- lapply(seq_len(k), function(index) {
    openssl::read_ed25519_key(as.raw((seq.int(0L, 31L) + index) %% 256L))
  })
  peers <- paste0("site_", seq_len(k))
  names(private) <- peers
  pins <- stats::setNames(vapply(private, function(key) {
    .count_exec_b64url(as.list(key)$pubkey$data)
  }, character(1L)), peers)
  config <- list(
    version = "dsvert-dp-count-config-v1",
    domain = "study-domain",
    cohort_id = "cohort-v1",
    dataset_id = "cohort-table",
    dataset_version = "v1",
    privacy_unit_column = "patient_id",
    alignment_purpose = "patient-record-alignment-v1",
    count_upper_bound = 1000L,
    max_records_per_unit = 1L,
    overflow_policy = "reject_operation",
    privacy = list(epsilon = 1, delta = 1e-6),
    calibration = list(implementation_delta = 1e-9),
    peer_pins = pins,
    backend_build_sha256 = strrep("a", 64L),
    transport_chunk_coordinates = 4096L)
  receipts <- stats::setNames(lapply(seq_along(peers), function(index) {
    list(
      peer_identity_pk = unname(pins[[peers[[index]]]]),
      snapshot_commitment = sprintf("%064x", index),
      sampler_plan = .count_exec_plan())
  }), peers)
  contract <- .dsvert_dp_count_client_contract_v1(config, receipts)
  roles <- .dsvert_dp_count_client_execution_roles_v1(contract)
  authorities <- unname(roles$peers[c("garbler", "evaluator")])
  session_id <- "12345678-1234-4234-9234-123456789abc"
  authorizations <- stats::setNames(lapply(authorities, function(peer) {
    role <- names(roles$peers)[match(peer, unname(roles$peers))]
    list(
      version = .DSVERT_CLIENT_DP_COUNT_PUBLIC_AUTH_VERSION,
      session_id = session_id,
      artifact_key = contract$artifact_key,
      analysis_binding_sha256 = roles$binding$sha256,
      local_authority = list(
        peer_name = peer,
        identity_pk = unname(roles$identities[[role]]),
        role = unname(role)),
      seed_commitment = strrep(if (role == "garbler") "b" else "c", 64L))
  }), authorities)
  transport <- stats::setNames(
    as.list(paste0("transport_pk_", seq_along(authorities))), authorities)
  attr(transport, "exact_gc_cleanup_capabilities") <- stats::setNames(
    as.list(paste0("cleanup_", seq_along(authorities))), authorities)
  attr(transport, "exact_gc_cleanup_purpose") <-
    .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE
  prepared <- list(
    version = "dsvert-dp-count-prepared-v1",
    mode = "add_remove_dp",
    payload = list(
      session_id = session_id,
      contract = contract,
      authorities = authorities,
      authorizations = authorizations,
      transport = transport))
  result <- list(
    private = private, config = config, contract = contract,
    roles = roles, authorities = authorities, prepared = prepared)
  result$start <- .count_exec_start(result)
  result$final <- .count_exec_final(result)
  result$release <- .count_exec_release(result)
  assign(cache_key, result, envir = .count_exec_fixture_cache)
  result
}

.count_exec_start <- function(fixture) {
  roles <- c("garbler", "evaluator")
  peers <- fixture$authorities
  stats::setNames(lapply(roles, function(role) {
    other <- setdiff(roles, role)
    list(
      capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
      peer_id = .dsvert_exact_gc_identity_peer_id(
        fixture$roles$identities[[role]]),
      peer_peer_id = .dsvert_exact_gc_identity_peer_id(
        fixture$roles$identities[[other]]),
      role = role,
      context_hash = strrep("d", 64L),
      operation = "joint-dp-laplace-v2",
      output_kind = "joint-dp-ring-share-v2",
      purpose = paste0("joint-dp-laplace-v2/", strrep("e", 64L)),
      source_producer = "count.scalar.v1",
      ring_bits = 127L,
      frac_bits = 0L,
      vector_len = 1L,
      threshold = "1",
      chunk_bytes = 16384L,
      ttl_seconds = 180L,
      max_runtime_seconds = 21600L,
      worker_heartbeat = 1L,
      state = "running",
      stored = FALSE,
      analysis_binding_sha256 = fixture$roles$binding$sha256)
  }), peers)
}

.count_exec_final <- function(fixture) {
  ciphertext <- "QUJD"
  garbler <- fixture$authorities[[1L]]
  evaluator <- fixture$authorities[[2L]]
  list(
    version = .DSVERT_CLIENT_DP_COUNT_FINAL_TRANSFER_VERSION,
    state = "final_share_sealed",
    artifact_key = fixture$contract$artifact_key,
    contract_sha256 = .dsvert_dp_count_client_hash_v1(
      .DSVERT_CLIENT_DP_COUNT_CONTRACT_DOMAIN, fixture$contract),
    analysis_binding_sha256 = fixture$roles$binding$sha256,
    circuit = paste0("joint-dp-laplace-v2/", strrep("e", 64L)),
    ciphertext = ciphertext,
    transfer = list(
      ticket = "opaque_ticket",
      transfer_id = paste0("tb_", strrep("1", 32L)),
      capability_id = .DSVERT_CLIENT_DP_COUNT_FINAL_CAPABILITY,
      sender_name = garbler,
      recipient_name = evaluator,
      payload_chars = nchar(ciphertext, type = "bytes"),
      payload_sha256 = paste0(openssl::sha256(charToRaw(ciphertext)))),
    intermediate_values_exposed = FALSE,
    capability_available = TRUE)
}

.count_exec_release <- function(fixture, value = "17") {
  contract <- fixture$contract
  privacy <- contract$semantic$privacy
  mechanism <- privacy$mechanism
  core <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_COUNT_RELEASE_VERSION,
    artifact_key = contract$artifact_key,
    contract_sha256 = .dsvert_dp_count_client_hash_v1(
      .DSVERT_CLIENT_DP_COUNT_CONTRACT_DOMAIN, contract),
    analysis_binding_sha256 = fixture$roles$binding$sha256,
    worker_static_sha256 = .dsvert_dp_count_client_hash_v1(
      .DSVERT_CLIENT_DP_COUNT_WORKER_DOMAIN,
      .dsvert_dp_count_client_worker_static_v1(contract)),
    circuit = paste0("joint-dp-laplace-v2/", strrep("e", 64L)),
    mechanism = list(
      family = mechanism$family,
      version = mechanism$version,
      sampler = mechanism$calibration$sampler,
      epsilon = privacy$epsilon,
      delta = privacy$delta,
      implementation_delta = mechanism$calibration$implementation_delta,
      sensitivity_l1 = 1),
    bounds = list(lower = "0", upper = "1000"),
    value = value,
    source_identity_pk = unname(fixture$roles$identities[["garbler"]]),
    finalizer_identity_pk =
      unname(fixture$roles$identities[["evaluator"]]),
    backend = "exact-gc-joint-dp-laplace-ring127-v2",
    postprocessing =
      "one-joint-noise-draw-and-one-clamp-inside-exact-gc",
    intermediate_values_exposed = FALSE,
    public_openings = 1))
  release_sha256 <- .dsvert_dp_count_client_hash_v1(
    .DSVERT_CLIENT_DP_COUNT_RELEASE_DOMAIN, core)
  signed <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    core, list(release_sha256 = release_sha256)))
  evaluator <- fixture$authorities[[2L]]
  signature <- openssl::ed25519_sign(charToRaw(paste0(
    .DSVERT_CLIENT_DP_COUNT_RELEASE_SIGNATURE_DOMAIN,
    .dsvert_joint_dp_client_json(signed))), fixture$private[[evaluator]])
  .dsvert_dp_analysis_client_canonical_value_v1(c(
    signed, list(signature = .count_exec_b64url(signature))))
}

.count_exec_mock <- function(fixture, mutate = NULL, fail = NULL) {
  state <- new.env(parent = emptyenv())
  state$calls <- list()
  state$run <- list()
  state$stores <- list()
  state$cleanups <- 0L
  state$aborts <- 0L
  state$cleanup_peers <- list()
  state$abort_peers <- list()
  state$authorization_roles <- character()
  state$start <- fixture$start
  state$final <- fixture$final
  state$release <- fixture$release
  aggregate <- function(conns, expr, ...) {
    command <- as.character(expr[[1L]])
    state$calls[[length(state$calls) + 1L]] <- list(
      command = command, peers = names(conns), expression = expr)
    if (identical(fail, command)) {
      callback <- list(...)[["error"]]
      if (is.function(callback)) {
        callback(names(conns)[[1L]], "injected phase failure")
        return(stats::setNames(vector("list", length(conns)), names(conns)))
      }
      stop("injected phase failure")
    }
    if (command == "dsvertDPCountStartDS") {
      decoded <- jsonlite::fromJSON(
        .dsvert_dsi_text_decode(expr$authorizations_json),
        simplifyVector = FALSE)
      state$authorization_roles <- vapply(decoded, function(value) {
        value$local_authority$role
      }, character(1L))
      values <- state$start[names(conns)]
    } else if (command == "dsvertDPCountFinalShareDS") {
      values <- stats::setNames(list(state$final), names(conns))
    } else if (command == "dsvertDPCountReleaseDS") {
      values <- stats::setNames(list(state$release), names(conns))
    } else {
      stop("unexpected Count endpoint: ", command)
    }
    if (is.function(mutate)) values <- mutate(command, values, state)
    values
  }
  prepare <- function(data_name, datasources, .aggregate) fixture$prepared
  run_exact <- function(...) {
    state$run[[length(state$run) + 1L]] <- list(...)
    if (identical(fail, "run_exact")) stop("injected worker failure")
    invisible(list())
  }
  store <- function(...) {
    state$stores[[length(state$stores) + 1L]] <- list(...)
    if (identical(fail, "store_typed")) stop("injected store failure")
    invisible(TRUE)
  }
  cleanup <- function(...) {
    state$cleanups <- state$cleanups + 1L
    state$cleanup_peers[[state$cleanups]] <- names(list(...)[[1L]])
    invisible(TRUE)
  }
  abort <- function(...) {
    state$aborts <- state$aborts + 1L
    state$abort_peers[[state$aborts]] <- names(list(...)[[1L]])
    invisible(NULL)
  }
  list(
    state = state, aggregate = aggregate, prepare = prepare,
    run_exact = run_exact, store = store, cleanup = cleanup, abort = abort)
}

.count_exec_run <- function(
    fixture, mock,
    .new_context = function() list(
      operation_id = paste0("op_", strrep("2", 32L)),
      source_key = paste0("exact_gc_in_", strrep("2", 32L)),
      output_key = paste0("exact_gc_out_", strrep("2", 32L)))) {
  order <- rev(names(fixture$config$peer_pins))
  datasources <- stats::setNames(lapply(
    order, function(...) structure(list(), class = "mock")), order)
  .dsvert_dp_count_execute_v1(
    "D", datasources, .aggregate = mock$aggregate,
    .prepare = mock$prepare, .new_context = .new_context,
    .run_exact = mock$run_exact, .store_typed = mock$store,
    .cleanup = mock$cleanup, .abort = mock$abort)
}

test_that("Count execution is K-generic but computes on exactly two roles", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .count_exec_fixture(k)
    mock <- .count_exec_mock(fixture)
    result <- .count_exec_run(fixture, mock)

    expect_identical(
      vapply(mock$state$calls, `[[`, character(1L), "command"),
      c("dsvertDPCountStartDS", "dsvertDPCountFinalShareDS",
        "dsvertDPCountReleaseDS"))
    expect_identical(mock$state$calls[[1L]]$peers,
                     fixture$authorities)
    expect_identical(mock$state$calls[[2L]]$peers,
                     fixture$authorities[[1L]])
    expect_identical(mock$state$calls[[3L]]$peers,
                     fixture$authorities[[2L]])
    expect_identical(mock$state$authorization_roles,
                     c("garbler", "evaluator"))
    expect_length(mock$state$run, 1L)
    expect_identical(mock$state$run[[1L]]$servers, match(
      fixture$authorities, rev(names(fixture$config$peer_pins))))
    expect_identical(mock$state$run[[1L]]$operation,
                     "joint-dp-laplace-v2")
    expect_identical(mock$state$run[[1L]]$ring, 127L)
    expect_identical(mock$state$run[[1L]]$frac_bits, 0L)
    expect_identical(mock$state$run[[1L]]$vector_len, 1L)
    expect_identical(mock$state$run[[1L]]$transport_ready, TRUE)
    expect_identical(mock$state$run[[1L]]$analysis_contract,
                     fixture$contract)
    expect_length(mock$state$stores, 1L)
    expect_identical(names(mock$state$stores[[1L]]$conn),
                     fixture$authorities[[2L]])
    expect_identical(names(mock$state$stores[[1L]]$producer_conn),
                     fixture$authorities[[1L]])
    expect_identical(names(result), c("version", "mode", "payload"))
    expect_identical(result$version,
                     "dsvert-dp-count-execution-result-v1")
    expect_identical(result$mode, "add_remove_dp")
    expect_identical(names(result$payload), c(
      "release", "finalizer_peer", "accuracy_95_abs",
      "accuracy_95_confidence", "accuracy_95_method"))
    expect_identical(result$payload$release$value, "17")
    expect_identical(result$payload$finalizer_peer,
                     fixture$authorities[[2L]])
    expect_identical(result$payload$accuracy_95_abs, 13L)
    expect_identical(result$payload$accuracy_95_confidence, 0.95)
    expect_identical(result$payload$accuracy_95_method,
      "conservative_truncated_dyadic_two_geometric_tail_bound_v1")
    expect_identical(mock$state$cleanups, 1L)
    expect_identical(mock$state$cleanup_peers[[1L]],
                     fixture$authorities)
    expect_identical(mock$state$aborts, 0L)
  }
})

test_that("Count execution rejects cardinality, order and closed-state tamper", {
  fixture <- .count_exec_fixture(3L)

  bad <- fixture
  bad$prepared$payload$authorities <-
    rev(bad$prepared$payload$authorities)
  mock <- .count_exec_mock(bad)
  expect_error(.count_exec_run(bad, mock), "authority|binding")
  expect_length(mock$state$calls, 0L)
  expect_identical(mock$state$cleanups, 1L)
  expect_identical(mock$state$aborts, 0L)

  mock <- .count_exec_mock(fixture)
  expect_error(.count_exec_run(
    fixture, mock, .new_context = function() {
      stop("injected context failure")
    }), "context failure")
  expect_identical(mock$state$cleanups, 1L)
  expect_identical(mock$state$aborts, 0L)

  mutations <- list(
    start_leak = function(command, values, state) {
      if (command == "dsvertDPCountStartDS") values[[1L]]$raw_count <- 17L
      values
    },
    start_missing = function(command, values, state) {
      if (command == "dsvertDPCountStartDS") values <- values[-1L]
      values
    },
    final_leak = function(command, values, state) {
      if (command == "dsvertDPCountFinalShareDS") {
        values[[1L]]$output_share <- "forbidden"
      }
      values
    },
    final_route = function(command, values, state) {
      if (command == "dsvertDPCountFinalShareDS") {
        values[[1L]]$transfer$recipient_name <-
          fixture$authorities[[1L]]
      }
      values
    },
    release_leak = function(command, values, state) {
      if (command == "dsvertDPCountReleaseDS") {
        values[[1L]]$preclamp_value <- "18"
      }
      values
    },
    release_bound = function(command, values, state) {
      if (command == "dsvertDPCountReleaseDS") values[[1L]]$value <- "1001"
      values
    },
    release_signature = function(command, values, state) {
      if (command == "dsvertDPCountReleaseDS") {
        values[[1L]]$signature <- paste0(
          if (substr(values[[1L]]$signature, 1L, 1L) == "A") "B" else "A",
          substring(values[[1L]]$signature, 2L))
      }
      values
    })
  for (name in names(mutations)) {
    mock <- .count_exec_mock(fixture, mutate = mutations[[name]])
    expect_error(.count_exec_run(fixture, mock), info = name)
    expect_identical(mock$state$cleanups, 1L, info = name)
    expect_identical(mock$state$aborts, 1L, info = name)
  }
})

test_that("Count execution aborts failures and always consumes cleanup authority", {
  fixture <- .count_exec_fixture(3L)
  for (phase in c(
      "dsvertDPCountStartDS", "run_exact", "dsvertDPCountFinalShareDS",
      "store_typed", "dsvertDPCountReleaseDS")) {
    mock <- .count_exec_mock(fixture, fail = phase)
    expect_error(
      .count_exec_run(fixture, mock),
      "injected|DataSHIELD transport failed", info = phase)
    expect_identical(mock$state$aborts, 1L, info = phase)
    expect_identical(mock$state$cleanups, 1L, info = phase)
  }
})

test_that("Count release is deterministic, bounded and contains no execution state", {
  fixture <- .count_exec_fixture(5L)
  first <- .count_exec_mock(fixture)
  second <- .count_exec_mock(fixture)
  result_a <- .count_exec_run(fixture, first)
  result_b <- .count_exec_run(fixture, second)
  release_a <- result_a$payload$release
  release_b <- result_b$payload$release

  expect_identical(result_a, result_b)
  expect_identical(as.numeric(release_a$value), 17)
  expect_false(any(c(
    "session_id", "operation_id", "request_id", "timestamp",
    "raw_count", "exact_count", "preclamp_value", "source_share",
    "output_share", "validity_share", "seed", "ciphertext") %in%
      names(release_a)))
  expect_identical(release_a$intermediate_values_exposed, FALSE)
  expect_identical(as.numeric(release_a$public_openings), 1)
  exposed <- paste(unlist(result_a, use.names = TRUE), collapse = "|")
  expect_false(grepl(
    "session_id|operation_id|source_key|output_key|ciphertext|ticket|seed",
    exposed, ignore.case = TRUE))
})

test_that("Count 95% accuracy radius is an exact conservative dyadic bound", {
  expect_identical(.dsvert_dp_count_client_accuracy_95_v1(
    .count_exec_plan()), 13L)

  plan <- .count_exec_plan()
  plan$stop_numerator <- "1"
  plan$max_geometric_steps <- 1L
  expect_identical(.dsvert_dp_count_client_accuracy_95_v1(plan), 1L)
})

test_that("fixed Count execution returns public consensus without session work", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .count_exec_fixture(k)
    peers <- names(fixture$config$peer_pins)
    declaration <- .dsvert_dp_analysis_client_canonical_value_v1(list(
      version = "dsvert-fixed-cohort-count-declaration-v1",
      domain = "study-domain",
      cohort_id = "cohort-v1",
      dataset_id = "cohort-table",
      dataset_version = "v1",
      privacy_unit_column = "patient_id",
      alignment_purpose = "patient-record-alignment-v1",
      adjacency = "replace_one_fixed_cohort",
      fixed_cohort_size = 1000L,
      peer_pins = as.list(fixture$config$peer_pins)))
    prepared <- list(
      version = "dsvert-dp-count-prepared-v1",
      mode = "fixed_cohort_public",
      payload = list(
        declaration = declaration,
        receipt_set_sha256 = strrep("f", 64L),
        peer_count = k))
    datasources <- stats::setNames(lapply(
      rev(peers), function(...) structure(list(), class = "mock")),
      rev(peers))
    forbidden <- function(...) stop("fixed Count entered session/MPC")
    result <- .dsvert_dp_count_execute_v1(
      "D", datasources,
      .aggregate = forbidden,
      .prepare = function(...) prepared,
      .new_context = forbidden,
      .run_exact = forbidden,
      .store_typed = forbidden,
      .cleanup = forbidden,
      .abort = forbidden)

    expect_identical(result, list(
      version = "dsvert-dp-count-execution-result-v1",
      mode = "fixed_cohort_public",
      payload = list(
        declaration = declaration,
        receipt_set_sha256 = strrep("f", 64L),
        peer_count = k)))
  }
})
