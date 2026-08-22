test_that("Frequency preparation rejects an implicit source before DSI", {
  calls <- 0L
  aggregate <- function(...) {
    calls <<- calls + 1L
    stop("unexpected DSI call")
  }
  datasources <- list(
    site_a = structure(list(), class = "mock"),
    site_b = structure(list(), class = "mock"))

  expect_error(
    .dsvert_dp_frequency_prepare_v1(
      data_name = "aligned", variable_name = "status",
      source_owner = NULL, datasources = datasources,
      .aggregate = aggregate),
    "explicit|source owner")
  expect_identical(calls, 0L)
})

test_that("Frequency preparation phases are classified as idempotent", {
  endpoints <- c(
    "dsvertDPFrequencyClaimDS", "dsvertDPFrequencyCompileDS",
    "dsvertDPFrequencyAuthorizeDS", "dsvertDPFrequencyCleanupDS")
  expect_true(all(endpoints %in% .DSVERT_IDEMPOTENT_TYPED_PRODUCERS))
})

.frequency_client_convolution_plan <- function(request) {
  plan <- .client_fresh_go_laplace_plan_v3(request)
  plan$version <- paste0(
    "dsvert-joint-dp-vector-independent-full-draw-convolution-plan-v3")
  plan$sampler <- .DSVERT_CLIENT_VECTOR_SAMPLER
  plan$maximum_chunk_coordinates <- min(8192L, request$total_coordinate_count)
  plan$independent_noise_peer_count <- 2
  plan$complete_epsilon_per_peer <- TRUE
  plan$epsilon_divided_by_peer_count <- FALSE
  plan$geometric_variables_per_peer_per_coordinate <- 2
  plan$geometric_variables_total_per_coordinate <- 4
  plan$per_peer_implementation_delta_numerator <-
    plan$implementation_delta_numerator
  plan$per_peer_implementation_delta_denominator <-
    plan$implementation_delta_denominator
  plan$per_peer_implementation_delta_bound <- plan$implementation_delta_bound
  plan$release_implementation_delta_aggregation <- "max_per_peer_not_sum"
  plan$two_peer_ideal_transfer_delta_numerator <- "0"
  plan$two_peer_ideal_transfer_delta_denominator <- "1"
  plan$two_peer_ideal_transfer_delta_bound <- "0"
  plan$threat_model <- "one pinned noise peer remains non-colluding"
  plan$privacy_argument <- paste(
    "each complete-epsilon draw is DP; adding the other independent",
    "draw is post-processing")
  plan
}

test_that("Frequency accepts the real convolution planner rational bounds", {
  request <- .dsvert_dp_analysis_frequency_candidate_requests_v2(
    list(adjacency = "add_remove_patient", epsilon = 1, delta = 0.01),
    list(implementation_delta = 0.001), 3L)$convolution
  plan <- .frequency_client_convolution_plan(request)
  profile <- .dsvert_dp_analysis_frequency_profile_v1(
    "independent_full_global_draw_convolution_ring128_v3")
  expect_gt(nchar(plan$per_peer_implementation_delta_numerator), 128L)
  expect_silent(.dsvert_vector_plan_validate(
    plan, .dsvert_vector_hash(plan), list(
      exact_gc = FALSE, gaussian = FALSE, plan_version = profile$plan,
      sampler = profile$sampler, mechanism = NULL,
      manifest_selection = NULL), 3L, request$sensitivity_steps))
  plan$per_peer_implementation_delta_numerator <- strrep("1", 1025L)
  expect_error(.dsvert_vector_plan_validate(
    plan, .dsvert_vector_hash(plan), list(
      exact_gc = FALSE, gaussian = FALSE, plan_version = profile$plan,
      sampler = profile$sampler, mechanism = NULL,
      manifest_selection = NULL), 3L, request$sensitivity_steps),
    "invalid")
})

.frequency_client_b64url <- function(value) {
  sub("=+$", "", chartr("+/", "-_", gsub(
    "[\r\n]", "", jsonlite::base64_enc(value))), perl = TRUE)
}

.frequency_client_keys <- function(k, peers = paste0(
    "site_", letters[seq_len(k)])) {
  stopifnot(length(peers) == k)
  private <- stats::setNames(lapply(seq_len(k), function(index) {
    openssl::read_ed25519_key(
      as.raw((seq.int(0L, 31L) + 7L * index) %% 256L))
  }), peers)
  public <- vapply(private, function(key) {
    .frequency_client_b64url(as.list(key)$pubkey$data)
  }, character(1L))
  list(private = private, public = public)
}

.frequency_client_sign <- function(message, key) {
  .frequency_client_b64url(openssl::ed25519_sign(message, key))
}

.frequency_client_selection <- function(kind, dimension = 3L) {
  privacy <- list(adjacency = "add_remove_patient", epsilon = 1, delta = 0.01)
  calibration <- list(implementation_delta = 0.001)
  requests <- .dsvert_dp_analysis_frequency_candidate_requests_v2(
    privacy, calibration, dimension)
  selector_request <- list(
    version = "dsvert-joint-dp-frequency-backend-selection-request-v1",
    adjacency = privacy$adjacency, coordinate_upper_bound = "1000",
    convolution_request = requests$convolution,
    gaussian_request = requests$gaussian)
  oracle <- .client_fresh_go_vector_plan(
    selector_request, "joint-dp-frequency-backend-select-v1")
  primitives <- .dsvert_dp_analysis_frequency_primitives_v2()
  plans <- list(
    convolution = oracle$convolution_plan, gaussian = oracle$gaussian_plan)
  certificates <- list(
    convolution = oracle$convolution_certificate,
    gaussian = oracle$gaussian_certificate)
  selected <- primitives[[kind]]
  selected_plan_hash <- .dsvert_dp_analysis_frequency_hash_v1(
    "dsVert/frequency/full-plan/v1|", plans[[kind]])
  selection_certificate <- list(
    version = "dsvert-joint-dp-frequency-backend-selection-certificate-v1",
    policy = "minimum_certified_simultaneous_95_abs_convolution_tie_v1",
    objective = "minimum_certified_simultaneous_95_abs",
    selected_primitive = selected,
    selected_plan_sha256 = selected_plan_hash,
    selected_simultaneous_95_abs =
      certificates[[kind]]$simultaneous_95_abs,
    tie_break = "convolution_laplace_v3_on_equal_certified_radius",
    input_scope = paste0(
      "public_adjacency_planner_requests_and_coordinate_upper_bound_only"),
    source_material_consulted = FALSE, private_randomness_consulted = FALSE,
    runtime_failure_consulted = FALSE, automatic_fallback = FALSE,
    utility_optimality_claimed = FALSE)
  candidates <- stats::setNames(lapply(names(primitives), function(candidate) {
    profile <- .dsvert_dp_analysis_frequency_profile_v1(
      primitives[[candidate]])
    certificate <- certificates[[candidate]]
    radius <- certificate$simultaneous_95_abs
    if (identical(kind, "gaussian") && identical(candidate, "convolution")) {
      radius <- as.character(openssl::bignum(
        certificates$gaussian$simultaneous_95_abs) + 1)
    }
    list(
      planner_request_sha256 = .dsvert_dp_analysis_frequency_hash_v1(
        profile$request_domain, requests[[candidate]]),
      full_plan_sha256 = .dsvert_dp_analysis_frequency_hash_v1(
        "dsVert/frequency/full-plan/v1|", plans[[candidate]]),
      accuracy_certificate_sha256 =
        .dsvert_dp_analysis_frequency_hash_v1(
          "dsVert/frequency/backend-selection/accuracy-certificate/v1|",
          certificate),
      simultaneous_95_abs = radius,
      absolute_support = certificate$absolute_support)
  }), names(primitives))
  summary <- list(
    version = "dsvert-frequency-backend-selection-v2",
    policy_sha256 = .dsvert_dp_analysis_frequency_policy_sha256_v2(),
    selection_certificate_sha256 = .dsvert_dp_analysis_frequency_hash_v1(
      "dsVert/frequency/backend-selection/certificate/v1|",
      selection_certificate),
    objective = "minimum_certified_simultaneous_95_abs",
    tie_break = "convolution_laplace_v3_on_equal_certified_radius",
    candidates = candidates, selected_primitive = selected,
    selected_simultaneous_95_abs =
      certificates[[kind]]$simultaneous_95_abs)
  list(
    summary = summary, selected_request = requests[[kind]],
    selected_plan = plans[[kind]],
    selected_accuracy_certificate = certificates[[kind]],
    selection_certificate = selection_certificate)
}

.frequency_client_factor <- function() {
  variable <- "status"
  list(
    version = "dsvert-psi-padded-factor-entry-v1",
    variable_name = variable,
    variable_id = paste0("var_", digest::digest(paste0(
      "dsVert/psi-padded/factor-variable/v1|",
      .dsvert_joint_dp_client_json(list(variable_name = variable))),
      "sha256", serialize = FALSE)),
    levels = as.list(c("case", "control", "other")), dimension = 3L)
}

.frequency_client_authorization <- function(
    peer, session_id, compiled, keys) {
  roles <- compiled$binding$binding$authority_roles
  role <- names(roles)[match(
    unname(compiled$config$peer_pins[[peer]]), unlist(roles, use.names = FALSE))]
  local <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    peer_name = peer,
    identity_pk = unname(compiled$config$peer_pins[[peer]]),
    role = unname(role)))
  private <- list(
    version = "dsvert-dp-frequency-session-authorization-v1",
    session_id = session_id, artifact_key = compiled$contract$artifact_key,
    config = compiled$config, config_sha256 = compiled$config_sha256,
    source_claim_sha256 = compiled$source_claim_sha256,
    receipt_peers = as.list(names(compiled$receipts)),
    receipt_set_sha256 = compiled$receipt_set_sha256,
    psi_run_sha256 = compiled$psi_run_sha256, contract = compiled$contract,
    contract_sha256 = compiled$contract_sha256,
    analysis_binding = compiled$binding$binding,
    analysis_binding_sha256 = compiled$binding$sha256,
    worker_static = compiled$worker_static,
    worker_static_sha256 = compiled$worker_static_sha256,
    local_authority = local)
  wire_private <- private
  wire_private$config$peer_pins <- as.list(wire_private$config$peer_pins)
  authorization_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
    "dsVert/dp-frequency/session-authorization/v1|", wire_private)
  unsigned <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = "dsvert-dp-frequency-public-authorization-v1",
    session_id = session_id, artifact_key = compiled$contract$artifact_key,
    config_sha256 = compiled$config_sha256,
    source_claim_sha256 = compiled$source_claim_sha256,
    receipt_set_sha256 = compiled$receipt_set_sha256,
    psi_run_sha256 = compiled$psi_run_sha256,
    contract_sha256 = compiled$contract_sha256,
    analysis_binding_sha256 = compiled$binding$sha256,
    worker_static_sha256 = compiled$worker_static_sha256,
    local_authority = local,
    commitment_context = compiled$worker_static$commitment_contexts[[role]],
    seed_commitment = digest::digest(
      paste0("Frequency seed|", peer), "sha256", serialize = FALSE),
    authorization_sha256 = authorization_sha256))
  .dsvert_dp_analysis_client_canonical_value_v1(c(
    unsigned, list(signature = .frequency_client_sign(charToRaw(paste0(
      "dsVert/dp-frequency/public-authorization-signature/v1|",
      .dsvert_joint_dp_client_json(unsigned))), keys$private[[peer]]))))
}

.frequency_client_fixture <- function(
    k = 3L, kind = c("convolution", "gaussian"), peers = paste0(
      "site_", letters[seq_len(k)])) {
  kind <- match.arg(kind)
  keys <- .frequency_client_keys(k, peers)
  pins <- keys$public[order(names(keys$public), method = "radix")]
  source_peer <- tail(names(pins), 1L)
  factor <- .frequency_client_factor()
  source <- list(
    alignment_purpose = "patient-record-alignment-v1",
    dataset_id = "cohort_table", dataset_version = "v1",
    id_column = "patient_id")
  claim_unsigned <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = "dsvert-dp-frequency-factor-claim-v2",
    source_peer_name = source_peer,
    source_identity_pk = unname(pins[[source_peer]]),
    psi_run_sha256 = strrep("1", 64L),
    attestation_id = paste0("attest_", strrep("2", 64L)),
    contract_hash = strrep("3", 64L),
    source_binding_id = paste0("source_", digest::digest(
      .dsvert_dp_frequency_client_wire_json_v1(source),
      "sha256", serialize = FALSE)),
    alignment_hash = strrep("4", 64L),
    alignment_purpose = source$alignment_purpose,
    dataset_id = source$dataset_id, dataset_version = source$dataset_version,
    privacy_unit_column = source$id_column,
    pinset_id = .dsvert_dp_frequency_client_pinset_v1(pins),
    capacity_bucket = 64L, factor_entry = factor,
    factor_entry_sha256 = .dsvert_dp_frequency_client_factor_hash_v1(factor)))
  claim <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    claim_unsigned, list(signature = .frequency_client_sign(charToRaw(paste0(
      "dsVert/dp-frequency/factor-claim/v2|",
      .dsvert_joint_dp_client_json(claim_unsigned))),
      keys$private[[source_peer]]))))
  config <- .dsvert_dp_frequency_client_config_v1(list(
    version = "dsvert-dp-frequency-config-v2", domain = "study-domain",
    cohort_id = "cohort-v1", dataset_id = source$dataset_id,
    dataset_version = source$dataset_version,
    privacy_unit_column = source$id_column,
    alignment_purpose = source$alignment_purpose,
    source_owner = list(
      peer_name = source_peer, identity_pk = unname(pins[[source_peer]])),
    source_binding_id = claim$source_binding_id,
    factor_domain = factor, factor_entry_sha256 = claim$factor_entry_sha256,
    coordinate_upper_bound = 1000L, max_records_per_unit = 1L,
    repeated_record_policy =
      "psi_v5_first_eligible_source_record_per_privacy_unit_v1",
    overflow_policy = "clip_to_psi_v5_first_eligible_source_record_v1",
    missingness_policy = "missing_or_out_of_domain_rows_are_ignored",
    privacy = list(
      adjacency = "add_remove_patient", epsilon = 1, delta = 0.01),
    calibration = list(implementation_delta = 0.001),
    peer_pins = as.list(pins), backend_build_sha256 = strrep("5", 64L),
    transport_chunk_coordinates = 3L,
    backend_selection = .frequency_client_selection(kind)))
  config_sha256 <- .dsvert_dp_frequency_client_config_hash_v1(config)
  claim_sha256 <- .dsvert_dp_frequency_client_claim_hash_v1(claim, config)
  receipts <- stats::setNames(lapply(seq_along(pins), function(index) {
    peer <- names(pins)[[index]]
    unsigned <- .dsvert_dp_analysis_client_canonical_value_v1(list(
      version = "dsvert-dp-frequency-receipt-v1", peer_name = peer,
      peer_identity_pk = unname(pins[[peer]]),
      config_sha256 = config_sha256, source_claim_sha256 = claim_sha256,
      psi_run_sha256 = claim$psi_run_sha256,
      snapshot_commitment = sprintf("%064x", index)))
    .dsvert_dp_analysis_client_canonical_value_v1(c(
      unsigned, list(signature = .frequency_client_sign(charToRaw(paste0(
        "dsVert/dp-frequency/receipt/v1|",
        .dsvert_joint_dp_client_json(unsigned))), keys$private[[peer]]))))
  }), names(pins))
  envelopes <- stats::setNames(lapply(names(pins), function(peer) {
    list(config = config, receipt = receipts[[peer]])
  }), names(pins))
  compiled <- .dsvert_dp_frequency_client_compile_v1(
    envelopes, names(envelopes), claim)
  session_id <- "00000000-0000-4000-8000-000000000001"
  authorizations <- stats::setNames(lapply(compiled$authorities, function(peer) {
    .frequency_client_authorization(peer, session_id, compiled, keys)
  }), compiled$authorities)
  transport_raw <- stats::setNames(lapply(seq_along(pins), function(index) {
    as.raw(rep(100L + index, 32L))
  }), names(pins))
  list(
    keys = keys, pins = pins, source_peer = source_peer, claim = claim,
    config = config, receipts = receipts, envelopes = envelopes,
    compiled = compiled, session_id = session_id,
    authorizations = authorizations, transport_raw = transport_raw)
}

test_that("Frequency independently compiles K2/K3/K5 and both real plans", {
  cases <- list(c(2L, "convolution"), c(3L, "gaussian"),
                c(5L, "convolution"))
  for (case in cases) {
    fixture <- .frequency_client_fixture(as.integer(case[[1L]]), case[[2L]])
    expected_secondary <- sort(
      setdiff(unname(fixture$pins),
              unname(fixture$pins[[fixture$source_peer]])),
      method = "radix")[[1L]]
    expect_identical(fixture$compiled$authorities[[1L]], fixture$source_peer)
    expect_identical(
      unname(fixture$pins[[fixture$compiled$authorities[[2L]]]]),
      expected_secondary)
    expect_identical(
      fixture$compiled$worker_static$selected_plan,
      .dsvert_dp_analysis_client_canonical_value_v1(
        fixture$config$backend_selection$selected_plan))
    expect_identical(
      names(fixture$compiled$worker_static$authority_roles),
      c("source_owner", "secondary_noise_authority"))
  }
})

test_that("Frequency worker is invariant across canonical DSV1", {
  for (kind in c("convolution", "gaussian")) {
    fixture <- .frequency_client_fixture(3L, kind)
    wire <- fixture$config
    wire$peer_pins <- as.list(wire$peer_pins)
    decoded <- jsonlite::fromJSON(
      .dsvert_joint_dp_client_json(wire),
      simplifyVector = FALSE)
    decoded <- .dsvert_dp_frequency_client_config_v1(decoded)
    worker <- fixture$compiled$worker_static
    roundtrip <- .dsvert_dp_frequency_client_worker_v1(
      fixture$compiled$contract, decoded, fixture$compiled$binding)

    expect_identical(roundtrip, worker)
    expect_identical(
      worker$selected_request,
      .dsvert_dp_analysis_client_canonical_value_v1(
        decoded$backend_selection$selected_request))
    expect_identical(
      worker$selected_plan,
      .dsvert_dp_analysis_client_canonical_value_v1(
        decoded$backend_selection$selected_plan))
    expect_identical(
      worker$selected_request,
      .dsvert_dp_analysis_client_canonical_value_v1(
        fixture$config$backend_selection$selected_request))
    expect_identical(
      worker$selected_plan,
      .dsvert_dp_analysis_client_canonical_value_v1(
        fixture$config$backend_selection$selected_plan))
  }
})

test_that("Frequency metadata cap includes the factor variable prefix", {
  levels <- list(strrep("x", 16L * 1024L * 1024L - 1L))
  expect_identical(
    .dsvert_dp_analysis_frequency_levels_dimension_v1(levels, 1L), 1)
  expect_null(.dsvert_dp_analysis_frequency_levels_dimension_v1(
    levels, 1L, "ab"))
})

.frequency_client_peer_binding <- function(fixture) {
  frequency <- .dsvert_exact_gc_frequency_binding(fixture$compiled)
  peers <- sort(fixture$compiled$authorities, method = "radix")
  transport <- vapply(fixture$transport_raw[peers],
                      .frequency_client_b64url, character(1L))
  contract <- list(
    version = "dsvert-exact-gc-frequency-peer-binding-v1",
    capability_id = "exact_gc_v1", session_id = fixture$session_id,
    consortium_id = fixture$compiled$contract$artifact_key,
    full_peer_pinset_sha256 = digest::digest(
      .dsvert_dp_frequency_client_wire_json_v1(as.list(fixture$pins)),
      "sha256", serialize = FALSE),
    designated_peers = as.list(peers),
    designated_peer_pinset = as.list(fixture$pins[peers]),
    identity_pks = as.list(fixture$pins[peers]),
    transport_pks = as.list(transport),
    frequency_binding = frequency$binding,
    frequency_binding_sha256 = frequency$sha256)
  list(
    contract = contract, frequency = frequency,
    sha256 = digest::digest(
      .dsvert_dp_frequency_client_wire_json_v1(contract),
      "sha256", serialize = FALSE))
}

.frequency_client_cleanup_capability <- function(fixture, peer, digest) {
  contract <- list(
    version = "dsvert-exact-gc-cleanup-capability-v1",
    session_id = fixture$session_id,
    cleanup_purpose = "dp.frequency.exact-session.v1",
    operation_scope = "all_and_only_operations_in_bound_exact_session_v1",
    peer_binding_digest = digest)
  envelope <- list(
    version = "dsvert-exact-gc-cleanup-capability-v1", contract = contract,
    signature = .frequency_client_sign(charToRaw(paste0(
      "dsVert/exact-gc/cleanup-capability/v1|",
      .dsvert_dp_frequency_client_wire_json_v1(contract))),
      fixture$keys$private[[peer]]))
  .dsvert_dp_frequency_client_wire_json_v1(envelope)
}

.frequency_client_flip <- function(value) {
  paste0(if (substr(value, 1L, 1L) == "A") "B" else "A",
         substring(value, 2L))
}

.frequency_client_mock <- function(fixture, mutate = NULL) {
  state <- new.env(parent = emptyenv())
  state$commands <- character()
  state$calls <- list()
  binding <- .frequency_client_peer_binding(fixture)
  aggregate <- function(conns, expr, ...) {
    command <- as.character(expr[[1L]])
    peers <- names(conns)
    state$commands <- c(state$commands, command)
    state$calls[[length(state$calls) + 1L]] <- list(
      command = command, peers = peers, expression = expr)
    values <- switch(command,
      dsvertDPFrequencyClaimDS = stats::setNames(
        list(fixture$claim), peers),
      dsvertDPFrequencyCompileDS = fixture$envelopes[peers],
      dsvertDPFrequencyAuthorizeDS = fixture$authorizations[peers],
      exactGCTransportInitDS = stats::setNames(lapply(peers, function(peer) {
        raw <- fixture$transport_raw[[peer]]
        list(
          capability_id = "exact_gc_v1",
          transport_pk = .frequency_client_b64url(raw),
          identity_pk = unname(fixture$pins[[peer]]),
          signature = .frequency_client_sign(raw,
            fixture$keys$private[[peer]]))
      }), peers),
      exactGCBindPeersDS = {
        expect_false("artifact_key" %in% names(expr))
        expect_identical(
          as.character(expr$cleanup_purpose),
          "dp.frequency.exact-session.v1")
        stats::setNames(lapply(peers, function(peer) list(
          bound = TRUE, capability_id = "exact_gc_v1",
          frequency_binding = binding$frequency$binding,
          frequency_binding_sha256 = binding$frequency$sha256,
          cleanup_purpose = "dp.frequency.exact-session.v1",
          cleanup_capability_json = .frequency_client_cleanup_capability(
            fixture, peer, binding$sha256))), peers)
      },
      dsvertDPFrequencyCleanupDS = stats::setNames(lapply(
        peers, function(...) list(cleaned = TRUE, state = "cleaned")), peers),
      stop("unexpected Frequency endpoint: ", command))
    if (is.function(mutate)) values <- mutate(command, values, state)
    values
  }
  list(state = state, aggregate = aggregate, binding = binding)
}

test_that("Frequency preparation uses one source, full K, then two roles", {
  cases <- list(c(2L, "convolution"), c(3L, "gaussian"),
                c(5L, "convolution"))
  for (case in cases) {
    fixture <- .frequency_client_fixture(as.integer(case[[1L]]), case[[2L]])
    mock <- .frequency_client_mock(fixture)
    order <- rev(names(fixture$pins))
    datasources <- stats::setNames(lapply(
      order, function(...) structure(list(), class = "mock")), order)
    result <- .dsvert_dp_frequency_prepare_v1(
      "D", "status", fixture$source_peer, datasources,
      .aggregate = mock$aggregate,
      .session_id = function() fixture$session_id)
    expect_identical(mock$state$commands, c(
      "dsvertDPFrequencyClaimDS", "dsvertDPFrequencyCompileDS",
      "dsvertDPFrequencyAuthorizeDS", "exactGCTransportInitDS",
      "exactGCBindPeersDS"))
    expect_identical(mock$state$calls[[1L]]$peers, fixture$source_peer)
    expect_identical(mock$state$calls[[2L]]$peers, order)
    for (index in 3:5) {
      expect_identical(
        mock$state$calls[[index]]$peers, fixture$compiled$authorities)
    }
    expect_false(any(grepl("Run|Worker|Exchange", mock$state$commands)))
    expect_identical(result$authorities, fixture$compiled$authorities)
    expect_identical(names(result$authorizations), result$authorities)
    expect_identical(
      attr(result$transport, "exact_gc_frequency_binding", exact = TRUE),
      mock$binding$frequency$binding)
    expect_identical(
      attr(result$transport, "exact_gc_peer_binding_digest", exact = TRUE),
      mock$binding$sha256)
    expect_identical(
      names(attr(result$transport, "exact_gc_cleanup_capabilities",
                 exact = TRUE)), result$authorities)
  }
})

test_that("Frequency transport orders peer names by radix bytes", {
  fixture <- .frequency_client_fixture(
    3L, "convolution", c("site_a", "Site_b", "site-c"))
  mock <- .frequency_client_mock(fixture)
  datasources <- stats::setNames(lapply(
    rev(names(fixture$pins)), function(...) structure(list(), class = "mock")),
    rev(names(fixture$pins)))

  result <- .dsvert_dp_frequency_prepare_v1(
    "D", "status", fixture$source_peer, datasources,
    .aggregate = mock$aggregate,
    .session_id = function() fixture$session_id)

  expect_identical(
    attr(result$transport, "exact_gc_peer_binding_digest", exact = TRUE),
    mock$binding$sha256)
})

test_that("Frequency preparation stops at every signed boundary", {
  fixture <- .frequency_client_fixture(3L, "convolution")
  datasources <- stats::setNames(lapply(
    rev(names(fixture$pins)),
    function(...) structure(list(), class = "mock")),
    rev(names(fixture$pins)))
  mutations <- list(
    claim = function(command, values, state) {
      if (command == "dsvertDPFrequencyClaimDS")
        values[[1L]]$signature <- .frequency_client_flip(values[[1L]]$signature)
      values
    },
    receipt = function(command, values, state) {
      if (command == "dsvertDPFrequencyCompileDS")
        values[[1L]]$receipt$signature <-
          .frequency_client_flip(values[[1L]]$receipt$signature)
      values
    },
    authorization = function(command, values, state) {
      if (command == "dsvertDPFrequencyAuthorizeDS")
        values[[1L]]$authorization_sha256 <- strrep("f", 64L)
      values
    },
    authorization_failure = function(command, values, state) {
      if (command == "dsvertDPFrequencyAuthorizeDS")
        stop("[dsvert_test_terminal:v1] simulated partial authorization failure")
      values
    },
    transport = function(command, values, state) {
      if (command == "exactGCTransportInitDS")
        values[[1L]]$signature <- .frequency_client_flip(values[[1L]]$signature)
      values
    },
    binding = function(command, values, state) {
      if (command == "exactGCBindPeersDS")
        values[[1L]]$frequency_binding_sha256 <- strrep("e", 64L)
      values
    },
    cleanup = function(command, values, state) {
      if (command == "exactGCBindPeersDS") {
        value <- jsonlite::fromJSON(
          values[[1L]]$cleanup_capability_json, simplifyVector = FALSE)
        value$contract$peer_binding_digest <- strrep("d", 64L)
        values[[1L]]$cleanup_capability_json <-
          .dsvert_dp_frequency_client_wire_json_v1(value)
      }
      values
    })
  for (name in names(mutations)) {
    mock <- .frequency_client_mock(fixture, mutations[[name]])
    expect_error(.dsvert_dp_frequency_prepare_v1(
      "D", "status", fixture$source_peer, datasources,
      .aggregate = mock$aggregate,
      .session_id = function() fixture$session_id),
      "Frequency|frequency|exact MPC")
    cleanup <- which(
      mock$state$commands == "dsvertDPFrequencyCleanupDS")
    if (name %in% c(
        "authorization", "authorization_failure", "transport", "binding",
        "cleanup")) {
      expect_length(cleanup, 1L)
      expect_identical(
        mock$state$calls[[cleanup]]$peers, fixture$compiled$authorities)
    } else {
      expect_length(cleanup, 0L)
    }
  }
})
