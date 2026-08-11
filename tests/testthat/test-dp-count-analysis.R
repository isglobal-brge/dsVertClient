.count_client_b64url <- function(value) {
  encoded <- gsub("[\r\n]", "", jsonlite::base64_enc(value))
  chartr("+/", "-_", sub("=+$", "", encoded))
}

.count_client_keys <- function(k) {
  private <- lapply(seq_len(k), function(index) {
    openssl::read_ed25519_key(as.raw((seq.int(0L, 31L) + index) %% 256L))
  })
  peers <- paste0("site_", seq_len(k))
  names(private) <- peers
  public <- stats::setNames(vapply(private, function(key) {
    .count_client_b64url(as.list(key)$pubkey$data)
  }, character(1L)), peers)
  list(private = private, public = public)
}

.count_client_hash <- function(domain, value) {
  digest::digest(
    charToRaw(paste0(domain, .dsvert_joint_dp_client_json(value))),
    algo = "sha256", serialize = FALSE)
}

.count_client_plan <- function() {
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

.count_client_config <- function(keys) {
  list(
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
    peer_pins = keys$public,
    backend_build_sha256 = strrep("a", 64L),
    transport_chunk_coordinates = 4096L)
}

.count_client_sign <- function(message, key) {
  .count_client_b64url(openssl::ed25519_sign(message, key))
}

.count_client_receipts <- function(config, keys) {
  wire_config <- config
  wire_config$peer_pins <- as.list(config$peer_pins)
  config_sha256 <- .count_client_hash(
    "dsVert/dp-count/config/v1|", wire_config)
  peers <- sort(names(config$peer_pins), method = "radix")
  receipts <- stats::setNames(lapply(seq_along(peers), function(index) {
    peer <- peers[[index]]
    unsigned <- list(
      version = "dsvert-dp-count-receipt-v1",
      peer_name = peer,
      peer_identity_pk = unname(config$peer_pins[[peer]]),
      config_sha256 = config_sha256,
      psi_run_sha256 = strrep("b", 64L),
      snapshot_commitment = sprintf("%064x", index),
      sampler_plan = .count_client_plan())
    message <- charToRaw(paste0(
      "dsVert/dp-count/receipt/v1|",
      .dsvert_joint_dp_client_json(unsigned)))
    c(unsigned, list(signature = .count_client_sign(
      message, keys$private[[peer]])))
  }), peers)
  receipts
}

.count_client_resign_receipt <- function(receipt, key) {
  unsigned <- receipt[setdiff(names(receipt), "signature")]
  receipt$signature <- .count_client_sign(charToRaw(paste0(
    "dsVert/dp-count/receipt/v1|",
    .dsvert_joint_dp_client_json(unsigned))), key)
  receipt
}

.count_client_fixture <- function(k) {
  keys <- .count_client_keys(k)
  config <- .count_client_config(keys)
  receipts <- .count_client_receipts(config, keys)
  envelopes <- stats::setNames(lapply(names(config$peer_pins), function(peer) {
    list(config = config, receipt = receipts[[peer]])
  }), names(config$peer_pins))
  list(keys = keys, config = config, receipts = receipts,
       envelopes = envelopes)
}

.count_client_public_authorization <- function(
    peer, session_id, compiled, keys) {
  contract <- compiled$contract
  binding <- .dsvert_exact_gc_analysis_binding(contract)
  worker <- .dsvert_dp_count_client_worker_static_v1(contract)
  role <- names(binding$binding$authority_roles)[match(
    unname(compiled$config$peer_pins[[peer]]),
    unlist(binding$binding$authority_roles, use.names = FALSE))]
  context <- worker[[paste0(role, "_commitment_context")]]
  core <- list(
    version = "dsvert-dp-count-public-authorization-v1",
    session_id = session_id,
    artifact_key = contract$artifact_key,
    config_sha256 = compiled$config_sha256,
    receipt_set_sha256 = compiled$receipt_set_sha256,
    psi_run_sha256 = compiled$psi_run_sha256,
    contract_sha256 = compiled$contract_sha256,
    analysis_binding_sha256 = binding$sha256,
    worker_static_sha256 = .count_client_hash(
      "dsVert/dp-count/worker-static/v1|", worker),
    local_authority = list(
      peer_name = peer,
      identity_pk = unname(compiled$config$peer_pins[[peer]]),
      role = unname(role)),
    commitment_context = context,
    seed_commitment = digest::digest(
      c(as.raw(strtoi(substring(context, seq(1L, 63L, 2L),
                                seq(2L, 64L, 2L)), 16L)),
        as.raw(rep(match(peer, names(keys$private)), 32L))),
      algo = "sha256", serialize = FALSE))
  authorization_sha256 <- .count_client_hash(
    "dsVert/dp-count/public-authorization/v1|", core)
  signed <- c(core, list(authorization_sha256 = authorization_sha256))
  c(signed, list(signature = .count_client_sign(
    charToRaw(paste0(
      "dsVert/dp-count/public-authorization-signature/v1|",
      .dsvert_joint_dp_client_json(signed))),
    keys$private[[peer]])))
}

.count_client_resign_authorization <- function(value, key) {
  core <- value[setdiff(names(value), c("authorization_sha256", "signature"))]
  value$authorization_sha256 <- .count_client_hash(
    "dsVert/dp-count/public-authorization/v1|", core)
  signed <- c(core, list(
    authorization_sha256 = value$authorization_sha256))
  value$signature <- .count_client_sign(charToRaw(paste0(
    "dsVert/dp-count/public-authorization-signature/v1|",
    .dsvert_joint_dp_client_json(signed))), key)
  value
}

.count_client_mock <- function(fixture, mutate = NULL) {
  state <- new.env(parent = emptyenv())
  state$commands <- character()
  state$calls <- list()
  state$authorize_wire <- NULL
  state$compiled <- .dsvert_dp_count_client_compile_v1(
    fixture$envelopes, names(fixture$envelopes))
  analysis <- .dsvert_exact_gc_analysis_binding(state$compiled$contract)
  roles <- unlist(
    analysis$binding$authority_roles[c("garbler", "evaluator")],
    use.names = TRUE)
  state$authorities <- names(fixture$config$peer_pins)[match(
    unname(roles), unname(fixture$config$peer_pins))]

  aggregate <- function(conns, expr, ...) {
    command <- as.character(expr[[1L]])
    state$commands <- c(state$commands, command)
    state$calls[[length(state$calls) + 1L]] <- list(
      command = command, peers = names(conns), expression = expr)
    if (identical(command, "dsvertDPCountCompileDS")) {
      values <- fixture$envelopes[names(conns)]
    } else if (identical(command, "dsvertDPCountAuthorizeDS")) {
      config_json <- .dsvert_dsi_text_decode(expr$config_json)
      receipts_json <- .dsvert_dsi_text_decode(expr$receipts_json)
      expect_identical(
        jsonlite::fromJSON(config_json, simplifyVector = FALSE)$peer_pins,
        as.list(state$compiled$config$peer_pins))
      decoded_receipts <- jsonlite::fromJSON(
        receipts_json, simplifyVector = FALSE)
      expect_null(names(decoded_receipts))
      expect_identical(
        vapply(decoded_receipts, `[[`, character(1L), "peer_name"),
        sort(names(fixture$receipts), method = "radix"))
      state$authorize_wire <- list(
        config_json = config_json, receipts_json = receipts_json)
      session_id <- as.character(expr$session_id)
      values <- stats::setNames(lapply(names(conns), function(peer) {
        .count_client_public_authorization(
          peer, session_id, state$compiled, fixture$keys)
      }), names(conns))
    } else if (identical(command, "exactGCTransportInitDS")) {
      values <- stats::setNames(lapply(names(conns), function(peer) {
        index <- match(peer, names(fixture$config$peer_pins))
        list(
          capability_id = "exact_gc_v1",
          transport_pk = .count_client_b64url(as.raw(rep(index + 30L, 32L))),
          identity_pk = unname(fixture$config$peer_pins[[peer]]),
          signature = .count_client_b64url(as.raw(rep(index + 60L, 64L))))
      }), names(conns))
    } else if (identical(command, "exactGCBindPeersDS")) {
      expect_false("analysis_contract_b64" %in% names(expr))
      expect_identical(expr$artifact_key,
                       state$compiled$contract$artifact_key)
      values <- stats::setNames(lapply(names(conns), function(peer) {
        cleanup <- as.character(jsonlite::toJSON(list(
          version = .DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION,
          contract = list(
            version = .DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION,
            session_id = as.character(expr$session_id),
            cleanup_purpose =
              .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
            operation_scope =
              "all_and_only_operations_in_bound_exact_session_v1",
            peer_binding_digest = digest::digest(
              paste0("binding|", peer), algo = "sha256",
              serialize = FALSE)),
          signature = .count_client_b64url(as.raw(rep(
            match(peer, names(fixture$config$peer_pins)) + 80L,
            64L)))), auto_unbox = TRUE, null = "null", na = "null",
          digits = 17, pretty = FALSE))
        list(
          capability_id = "exact_gc_v1", bound = TRUE,
          analysis_binding = analysis$binding,
          analysis_binding_sha256 = analysis$sha256,
          cleanup_purpose =
            .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
          cleanup_capability_json = cleanup)
      }), names(conns))
    } else {
      stop("unexpected endpoint: ", command)
    }
    if (is.function(mutate)) values <- mutate(command, values, state)
    values
  }
  list(state = state, aggregate = aggregate)
}

test_that("Count compile-authorize-bind is K-generic and role stable", {
  expect_false(any(c(
    "config", "receipts", "contract", "epsilon", "delta", "artifact_key") %in%
      names(formals(.dsvert_dp_count_compile_authorize_bind_v1))))
  session_id <- "12345678-1234-4234-9234-123456789abc"
  for (k in c(2L, 3L, 5L)) {
    fixture <- .count_client_fixture(k)
    order <- rev(names(fixture$envelopes))
    datasources <- stats::setNames(lapply(
      order, function(...) structure(list(), class = "mock")), order)
    mock <- .count_client_mock(fixture)

    result <- .dsvert_dp_count_compile_authorize_bind_v1(
      "D", datasources, .aggregate = mock$aggregate,
      .session_id = function() session_id)

    expect_identical(
      mock$state$commands,
      c("dsvertDPCountCompileDS", "dsvertDPCountAuthorizeDS",
        "exactGCTransportInitDS", "exactGCBindPeersDS"))
    expect_identical(mock$state$calls[[1L]]$peers, order)
    expect_identical(
      mock$state$calls[[2L]]$peers, mock$state$authorities)
    expect_identical(
      mock$state$calls[[3L]]$peers, mock$state$authorities)
    expect_identical(
      mock$state$calls[[4L]]$peers, mock$state$authorities)
    expect_identical(names(result), c(
      "session_id", "contract", "authorities", "authorizations",
      "transport"))
    expect_identical(result$session_id, session_id)
    expect_identical(result$contract, mock$state$compiled$contract)
    expect_identical(result$authorities, mock$state$authorities)
    expect_identical(names(result$authorizations), mock$state$authorities)
    expect_identical(names(result$transport), mock$state$authorities)
    expect_identical(
      names(attr(result$transport,
                 "exact_gc_cleanup_capabilities", exact = TRUE)),
      mock$state$authorities)
    expect_identical(
      attr(result$transport, "exact_gc_cleanup_purpose", exact = TRUE),
      .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE)
    bind <- mock$state$calls[[4L]]$expression
    expect_setequal(names(bind)[nzchar(names(bind))], c(
      "transport_keys_b64", "identity_info_b64", "session_id",
      "cleanup_purpose", "artifact_key"))
    expect_identical(
      bind$cleanup_purpose,
      .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE)
    expect_false(any(grepl(
      "analysis_contract|contract_json|receipt|config_json",
      names(bind), ignore.case = TRUE)))
    expect_false(grepl(
      "lifetime|budget|ledger|sqlite|catalog|data_name",
      paste(unlist(mock$state$authorize_wire, use.names = FALSE),
            collapse = ""), ignore.case = TRUE))
  }
})

test_that("Count client compiler has exact K, signature and server hash parity", {
  expected_artifacts <- c(
    `2` = "43699087a778670dc9560e2d08d5cc6ee9a8fe88ba95d14f5cca1cbfa24eb413",
    `3` = "dd7f9d853be96a30d9b5bf6adce94b61bcdeb9f0deb671387a9d81cc78b9050a",
    `5` = "202de715ec4969d33b4379535317e5491ddb1f3ad241a1f14698f66217e3465b")
  for (k in c(2L, 3L, 5L)) {
    fixture <- .count_client_fixture(k)
    compiled <- .dsvert_dp_count_client_compile_v1(
      fixture$envelopes[rev(names(fixture$envelopes))],
      rev(names(fixture$envelopes)))
    expect_identical(
      compiled$contract,
      .dsvert_dp_analysis_contract_validate_v1(compiled$contract))
    expect_identical(compiled$contract$artifact_key,
                     unname(expected_artifacts[[as.character(k)]]))
    expect_identical(length(compiled$config$peer_pins), k)
    expect_identical(names(compiled$receipts),
                     sort(names(fixture$receipts), method = "radix"))
  }

  fixture <- .count_client_fixture(3L)
  compiled <- .dsvert_dp_count_client_compile_v1(
    fixture$envelopes, names(fixture$envelopes))
  expect_identical(
    compiled$config_sha256,
    "080b424640dd19d6ba66486eb55d059db47db77c42536746ef6cb82f74ec846b")
  expect_identical(
    compiled$receipt_set_sha256,
    "2cc5814af6db7281125de27d6ed2dd29cacff7d5471efb96976ab9b625d10986")
  expect_identical(
    compiled$contract_sha256,
    "bd9c848fef424c9c45d89e8dee6402d918f95be205436ce84eec740c3fc364fa")
  expect_identical(
    compiled$binding$sha256,
    "6b6b67dcebfafcaa3cd3656bec73a99364b4f0d82411906f6ad27d8114637817")
  expect_identical(
    compiled$worker_static_sha256,
    "d9fe49a1a8eefe982205019fe77b9a4187edada13a69a8e03a8feb681bd8f141")
})

test_that("Count compile consensus rejects closed-schema and receipt tampering", {
  fixture <- .count_client_fixture(3L)
  compile <- function(envelopes = fixture$envelopes,
                      sites = names(fixture$envelopes)) {
    .dsvert_dp_count_client_compile_v1(envelopes, sites)
  }

  missing <- fixture$envelopes[-1L]
  expect_error(compile(missing), "complete federation")

  extra <- fixture$envelopes
  extra[[1L]]$unexpected <- TRUE
  expect_error(compile(extra), "compile envelope")

  config_extra <- fixture$envelopes
  config_extra[[1L]]$config$request_id <- "analyst-controlled"
  expect_error(compile(config_extra), "configuration")

  disagreement <- fixture$envelopes
  disagreement[[2L]]$config$count_upper_bound <- 999L
  expect_error(compile(disagreement), "disagree")

  high_epsilon <- fixture$envelopes
  high_epsilon <- lapply(high_epsilon, function(envelope) {
    envelope$config$privacy$epsilon <- 8 + 1e-12
    envelope
  })
  names(high_epsilon) <- names(fixture$envelopes)
  expect_error(compile(high_epsilon), "privacy parameters")

  wrong_sites <- paste0("alias_", seq_along(fixture$envelopes))
  expect_error(compile(fixture$envelopes, wrong_sites),
               "complete federation|full-K")

  receipt_extra <- fixture$envelopes
  receipt_extra[[1L]]$receipt$nonce <- "not-allowed"
  expect_error(compile(receipt_extra), "receipt fields")

  duplicate <- fixture$envelopes
  duplicate[[2L]]$receipt <- duplicate[[1L]]$receipt
  expect_error(compile(duplicate), "connected peers|coverage")

  bad_signature <- fixture$envelopes
  bad_signature[[1L]]$receipt$signature <- paste0(
    if (substr(bad_signature[[1L]]$receipt$signature, 1L, 1L) == "A")
      "B" else "A",
    substring(bad_signature[[1L]]$receipt$signature, 2L))
  expect_error(compile(bad_signature), "signature verification")

  different_run <- fixture$envelopes
  peer <- names(different_run)[[3L]]
  different_run[[peer]]$receipt$psi_run_sha256 <- strrep("c", 64L)
  different_run[[peer]]$receipt <- .count_client_resign_receipt(
    different_run[[peer]]$receipt, fixture$keys$private[[peer]])
  expect_error(compile(different_run), "PSI run")

  excessive <- fixture$envelopes
  peer <- names(excessive)[[1L]]
  excessive[[peer]]$receipt$sampler_plan$implementation_delta_numerator <- "2"
  excessive[[peer]]$receipt$sampler_plan$implementation_delta_bound <-
    "2/1000000000"
  excessive[[peer]]$receipt <- .count_client_resign_receipt(
    excessive[[peer]]$receipt, fixture$keys$private[[peer]])
  expect_error(compile(excessive), "exceeds calibration")
})

test_that("Count authorization fails before transport on every public tamper", {
  session_id <- "12345678-1234-4234-9234-123456789abc"
  fixture <- .count_client_fixture(3L)
  compiled <- .dsvert_dp_count_client_compile_v1(
    fixture$envelopes, names(fixture$envelopes))
  authorizations <- stats::setNames(lapply(
    compiled$authorities, function(peer) {
      .count_client_public_authorization(
        peer, session_id, compiled, fixture$keys)
    }), compiled$authorities)
  datasources <- stats::setNames(lapply(
    rev(names(fixture$envelopes)),
    function(...) structure(list(), class = "mock")),
    rev(names(fixture$envelopes)))
  mutations <- list(
    extra = function(command, values, state) {
      if (command == "dsvertDPCountAuthorizeDS") values[[1L]]$extra <- TRUE
      values
    },
    missing = function(command, values, state) {
      if (command == "dsvertDPCountAuthorizeDS") values <- values[-1L]
      values
    },
    contract = function(command, values, state) {
      if (command == "dsvertDPCountAuthorizeDS") {
        values[[1L]]$contract_sha256 <- strrep("f", 64L)
      }
      values
    },
    role = function(command, values, state) {
      if (command == "dsvertDPCountAuthorizeDS") {
        values[[1L]]$local_authority$role <- "evaluator"
      }
      values
    },
    authorization_hash = function(command, values, state) {
      if (command == "dsvertDPCountAuthorizeDS") {
        values[[1L]]$authorization_sha256 <- strrep("f", 64L)
      }
      values
    },
    signature = function(command, values, state) {
      if (command == "dsvertDPCountAuthorizeDS") {
        values[[1L]]$signature <- paste0(
          if (substr(values[[1L]]$signature, 1L, 1L) == "A") "B" else "A",
          substring(values[[1L]]$signature, 2L))
      }
      values
    },
    duplicate_seed = function(command, values, state) {
      if (command == "dsvertDPCountAuthorizeDS") {
        peers <- names(values)
        values[[2L]]$seed_commitment <- values[[1L]]$seed_commitment
        values[[2L]] <- .count_client_resign_authorization(
          values[[2L]], fixture$keys$private[[peers[[2L]]]])
      }
      values
    })
  for (name in names(mutations)) {
    changed <- mutations[[name]](
      "dsvertDPCountAuthorizeDS", authorizations, new.env())
    expect_error(.dsvert_dp_count_client_authorization_set_v1(
      changed, session_id, compiled), info = name)
  }

  mock <- .count_client_mock(fixture, mutations$signature)
  expect_error(.dsvert_dp_count_compile_authorize_bind_v1(
    "D", datasources, .aggregate = mock$aggregate,
    .session_id = function() session_id), "signature")
  expect_identical(
    mock$state$commands,
    c("dsvertDPCountCompileDS", "dsvertDPCountAuthorizeDS"))
})

test_that("Count control-plane result contains no protected statistic", {
  fixture <- .count_client_fixture(3L)
  datasources <- stats::setNames(lapply(
    names(fixture$envelopes),
    function(...) structure(list(), class = "mock")),
    names(fixture$envelopes))
  mock <- .count_client_mock(fixture)
  result <- .dsvert_dp_count_compile_authorize_bind_v1(
    "D", datasources, .aggregate = mock$aggregate,
    .session_id = function() "12345678-1234-4234-9234-123456789abc")
  collect_names <- function(value) {
    if (!is.list(value)) return(character())
    c(names(value), unlist(lapply(value, collect_names), use.names = FALSE))
  }
  expect_false(any(collect_names(result) %in% c(
    "raw_count", "exact_count", "local_count", "source_share",
    "output_share", "validity_share", "seed", "patient_ids", "members",
    "identifiers")))
  expect_false(any(vapply(result$authorizations, function(value) {
    any(c("config", "receipts", "contract", "worker_static") %in%
          names(value))
  }, logical(1L))))
})
