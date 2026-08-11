# Internal execution of one canonical Count analysis.  The client relays only
# signed public authorizations and opaque MPC ciphertext; it never receives an
# exact count, an output share, or either authority's sticky seed.

.DSVERT_CLIENT_DP_COUNT_FINAL_TRANSFER_VERSION <-
  "dsvert-dp-count-final-share-transfer-v1"
.DSVERT_CLIENT_DP_COUNT_FINAL_CAPABILITY <-
  "blob.analysis-dp.count-final-share.v1"
.DSVERT_CLIENT_DP_COUNT_RELEASE_VERSION <-
  "dsvert-dp-count-release-v1"
.DSVERT_CLIENT_DP_COUNT_RELEASE_DOMAIN <-
  "dsVert/dp-count/release/v1|"
.DSVERT_CLIENT_DP_COUNT_RELEASE_SIGNATURE_DOMAIN <-
  "dsVert/dp-count/release-signature/v1|"
.DSVERT_CLIENT_DP_COUNT_EXECUTION_RESULT_VERSION <-
  "dsvert-dp-count-execution-result-v1"
.DSVERT_CLIENT_DP_COUNT_ACCURACY_95_METHOD <-
  "conservative_truncated_dyadic_two_geometric_tail_bound_v1"
.DSVERT_CLIENT_DP_COUNT_AUTHORIZATIONS_MAX_BYTES <- 1024L * 1024L
.DSVERT_CLIENT_DP_COUNT_FINAL_CIPHERTEXT_MAX_BYTES <- 1024L * 1024L

.dsvert_dp_count_client_accuracy_95_v1 <- function(plan) {
  plan <- .dsvert_dp_count_client_plan_v1(plan)
  denominator <- openssl::bignum(2)^
    as.integer(plan$bernoulli_bits)
  continuation <- denominator - openssl::bignum(plan$stop_numerator)
  left <- 40 * continuation
  right <- denominator + continuation
  maximum <- as.integer(plan$max_geometric_steps)
  for (radius in 0:maximum) {
    if (left <= right) return(as.integer(radius))
    left <- left * continuation
    right <- right * denominator
  }
  maximum
}

.dsvert_dp_count_client_prepared_v1 <- function(prepared, datasources) {
  fields <- c("version", "mode", "payload")
  valid <- is.list(prepared) && !is.null(names(prepared)) &&
    !anyNA(names(prepared)) && !anyDuplicated(names(prepared)) &&
    setequal(names(prepared), fields) &&
    identical(prepared$version, .DSVERT_CLIENT_DP_COUNT_PREPARED_VERSION) &&
    is.character(prepared$mode) && length(prepared$mode) == 1L &&
    !is.na(prepared$mode) &&
    prepared$mode %in% c("add_remove_dp", "fixed_cohort_public") &&
    is.list(prepared$payload) &&
    is.list(datasources) && length(datasources) >= 2L &&
    !is.null(names(datasources)) && !anyNA(names(datasources)) &&
    !any(!nzchar(names(datasources))) && !anyDuplicated(names(datasources))
  if (!isTRUE(valid)) {
    stop("Invalid prepared Count execution mode.", call. = FALSE)
  }
  if (identical(prepared$mode, "add_remove_dp")) return(prepared)

  payload <- prepared$payload
  if (is.null(names(payload)) || anyNA(names(payload)) ||
      anyDuplicated(names(payload)) || !setequal(
        names(payload),
        c("declaration", "receipt_set_sha256", "peer_count"))) {
    stop("Invalid prepared fixed-cohort Count result.", call. = FALSE)
  }
  declaration <- .dsvert_dp_count_client_fixed_declaration_v1(
    payload$declaration, names(datasources))
  receipt_set_sha256 <- .dsvert_dp_count_client_hex_v1(
    payload$receipt_set_sha256, "fixed receipt-set hash")
  peer_count <- .dsvert_dp_count_client_positive_integer_v1(
    payload$peer_count, "fixed peer count", 4096)
  if (!identical(peer_count, as.numeric(length(datasources)))) {
    stop("Invalid prepared fixed-cohort Count peer count.", call. = FALSE)
  }
  list(
    version = .DSVERT_CLIENT_DP_COUNT_PREPARED_VERSION,
    mode = "fixed_cohort_public",
    payload = list(
      declaration = declaration,
      receipt_set_sha256 = receipt_set_sha256,
      peer_count = as.integer(peer_count)))
}

.dsvert_dp_count_client_execution_roles_v1 <- function(contract) {
  binding <- .dsvert_exact_gc_analysis_binding(contract)
  identities <- unlist(
    binding$binding$authority_roles[c("garbler", "evaluator")],
    use.names = TRUE)
  peers <- names(binding$full_pins)[match(
    unname(identities), unname(binding$full_pins))]
  if (length(peers) != 2L || anyNA(peers) || anyDuplicated(peers)) {
    stop("Count execution authorities do not match the full-K pinset.",
         call. = FALSE)
  }
  list(
    binding = binding,
    peers = stats::setNames(unname(peers), names(identities)),
    identities = identities)
}

.dsvert_dp_count_client_execution_context_v1 <- function(
    prepared, datasources) {
  fields <- c(
    "session_id", "contract", "authorities", "authorizations",
    "transport")
  if (!is.list(prepared) || is.null(names(prepared)) ||
      anyNA(names(prepared)) || anyDuplicated(names(prepared)) ||
      !setequal(names(prepared), fields) ||
      !is.list(datasources) || length(datasources) < 2L ||
      is.null(names(datasources)) || anyNA(names(datasources)) ||
      any(!nzchar(names(datasources))) || anyDuplicated(names(datasources))) {
    stop("Invalid prepared Count execution context.", call. = FALSE)
  }
  session_id <- prepared$session_id
  if (!is.character(session_id) || length(session_id) != 1L ||
      is.na(session_id) || !grepl(
        paste0("^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-",
               "[89ab][0-9a-f]{3}-[0-9a-f]{12}$"), session_id)) {
    stop("Invalid prepared Count session identifier.", call. = FALSE)
  }
  contract <- .dsvert_dp_analysis_contract_validate_v1(prepared$contract)
  roles <- .dsvert_dp_count_client_execution_roles_v1(contract)
  authorities <- unname(roles$peers[c("garbler", "evaluator")])
  transport <- prepared$transport
  capabilities <- attr(
    transport, "exact_gc_cleanup_capabilities", exact = TRUE)
  cleanup_purpose <- attr(
    transport, "exact_gc_cleanup_purpose", exact = TRUE)
  valid <- identical(prepared$authorities, authorities) &&
    all(authorities %in% names(datasources)) &&
    is.list(prepared$authorizations) &&
    !is.null(names(prepared$authorizations)) &&
    !anyNA(names(prepared$authorizations)) &&
    !anyDuplicated(names(prepared$authorizations)) &&
    setequal(names(prepared$authorizations), authorities) &&
    is.list(transport) && !is.null(names(transport)) &&
    !anyNA(names(transport)) && !anyDuplicated(names(transport)) &&
    setequal(names(transport), authorities) &&
    all(vapply(transport, function(value) {
      is.character(value) && length(value) == 1L && !is.na(value) &&
        nzchar(value)
    }, logical(1L))) &&
    is.list(capabilities) && !is.null(names(capabilities)) &&
    !anyNA(names(capabilities)) && !anyDuplicated(names(capabilities)) &&
    setequal(names(capabilities), authorities) &&
    identical(cleanup_purpose,
              .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE)
  if (!isTRUE(valid)) {
    stop("Invalid prepared Count authority or cleanup binding.",
         call. = FALSE)
  }
  list(
    session_id = session_id,
    contract = contract,
    roles = roles,
    authorities = authorities,
    conns = datasources[authorities],
    authorizations = prepared$authorizations[authorities],
    transport = transport[authorities])
}

.dsvert_dp_count_client_authorizations_wire_v1 <- function(context) {
  role_order <- c("garbler", "evaluator")
  values <- lapply(role_order, function(role) {
    peer <- unname(context$roles$peers[[role]])
    value <- context$authorizations[[peer]]
    local <- if (is.list(value)) value$local_authority else NULL
    valid <- is.list(value) && is.list(local) &&
      identical(value$version,
                .DSVERT_CLIENT_DP_COUNT_PUBLIC_AUTH_VERSION) &&
      identical(value$session_id, context$session_id) &&
      identical(value$artifact_key, context$contract$artifact_key) &&
      identical(value$analysis_binding_sha256,
                context$roles$binding$sha256) &&
      identical(local$peer_name, peer) &&
      identical(local$role, role) &&
      identical(local$identity_pk,
                unname(context$roles$identities[[role]]))
    if (!isTRUE(valid)) {
      stop("Invalid prepared Count public authorization order or binding.",
           call. = FALSE)
    }
    .dsvert_dp_analysis_client_canonical_value_v1(value)
  })
  json <- .dsvert_joint_dp_client_json(values)
  if (nchar(json, type = "bytes") >
      .DSVERT_CLIENT_DP_COUNT_AUTHORIZATIONS_MAX_BYTES) {
    stop("Count authorization input exceeds its bounded wire contract.",
         call. = FALSE)
  }
  json
}

.dsvert_dp_count_client_operation_v1 <- function(value) {
  fields <- c("operation_id", "source_key", "output_key")
  if (!is.list(value) || is.null(names(value)) ||
      anyNA(names(value)) || anyDuplicated(names(value)) ||
      !setequal(names(value), fields) ||
      !is.character(value$operation_id) ||
      !grepl("^op_[0-9a-f]{32}$", value$operation_id) ||
      !is.character(value$source_key) ||
      !grepl("^exact_gc_in_[0-9a-f]{32}$", value$source_key) ||
      !is.character(value$output_key) ||
      !grepl("^exact_gc_out_[0-9a-f]{32}$", value$output_key) ||
      !identical(sub("^op_", "", value$operation_id),
                 sub("^exact_gc_in_", "", value$source_key)) ||
      !identical(sub("^op_", "", value$operation_id),
                 sub("^exact_gc_out_", "", value$output_key))) {
    stop("Invalid Count exact-MPC operation context.", call. = FALSE)
  }
  value
}

.dsvert_dp_count_client_start_v1 <- function(values, context) {
  fields <- c(
    "capability_id", "peer_id", "peer_peer_id", "role", "context_hash",
    "operation", "output_kind", "purpose", "source_producer",
    "ring_bits", "frac_bits", "vector_len", "threshold", "chunk_bytes",
    "ttl_seconds", "max_runtime_seconds", "worker_heartbeat", "state",
    "stored", "analysis_binding_sha256")
  peers <- context$authorities
  if (!is.list(values) || is.null(names(values)) || anyNA(names(values)) ||
      anyDuplicated(names(values)) || !setequal(names(values), peers)) {
    stop("Count Start did not cover both bound authorities.", call. = FALSE)
  }
  values <- values[peers]
  for (role in c("garbler", "evaluator")) {
    peer <- unname(context$roles$peers[[role]])
    value <- values[[peer]]
    valid <- is.list(value) && !is.null(names(value)) &&
      !anyNA(names(value)) && !anyDuplicated(names(value)) &&
      setequal(names(value), fields) &&
      identical(value$capability_id,
                .DSVERT_CLIENT_EXACT_GC_CAPABILITY) &&
      identical(value$peer_id, .dsvert_exact_gc_identity_peer_id(
        context$roles$identities[[role]])) &&
      identical(value$peer_peer_id, .dsvert_exact_gc_identity_peer_id(
        context$roles$identities[[setdiff(
          c("garbler", "evaluator"), role)]])) &&
      identical(value$role, role) &&
      identical(value$operation, "joint-dp-laplace-v2") &&
      identical(value$output_kind, "joint-dp-ring-share-v2") &&
      identical(value$source_producer, "count.scalar.v1") &&
      identical(as.numeric(value$ring_bits), 127) &&
      identical(as.numeric(value$frac_bits), 0) &&
      identical(as.numeric(value$vector_len), 1) &&
      identical(value$state, "running") && identical(value$stored, FALSE) &&
      identical(value$analysis_binding_sha256,
                context$roles$binding$sha256) &&
      is.character(value$context_hash) &&
      grepl("^[0-9a-f]{64}$", value$context_hash) &&
      is.character(value$purpose) &&
      grepl("^joint-dp-laplace-v2/[0-9a-f]{64}$", value$purpose)
    if (!isTRUE(valid)) {
      stop("A Count authority returned an invalid closed Start state.",
           call. = FALSE)
    }
  }
  contexts <- vapply(values, `[[`, character(1L), "context_hash")
  purposes <- vapply(values, `[[`, character(1L), "purpose")
  if (length(unique(contexts)) != 1L || length(unique(purposes)) != 1L) {
    stop("Count authorities disagree on the exact-MPC Start contract.",
         call. = FALSE)
  }
  list(values = values, purpose = purposes[[1L]])
}

.dsvert_dp_count_client_final_share_v1 <- function(
    value, context, purpose) {
  fields <- c(
    "version", "state", "artifact_key", "contract_sha256",
    "analysis_binding_sha256", "circuit", "ciphertext", "transfer",
    "intermediate_values_exposed", "capability_available")
  contract_sha256 <- .dsvert_dp_count_client_hash_v1(
    .DSVERT_CLIENT_DP_COUNT_CONTRACT_DOMAIN, context$contract)
  garbler <- unname(context$roles$peers[["garbler"]])
  evaluator <- unname(context$roles$peers[["evaluator"]])
  valid <- is.list(value) && !is.null(names(value)) &&
    !anyNA(names(value)) && !anyDuplicated(names(value)) &&
    setequal(names(value), fields) &&
    identical(value$version,
              .DSVERT_CLIENT_DP_COUNT_FINAL_TRANSFER_VERSION) &&
    identical(value$state, "final_share_sealed") &&
    identical(value$artifact_key, context$contract$artifact_key) &&
    identical(value$contract_sha256, contract_sha256) &&
    identical(value$analysis_binding_sha256,
              context$roles$binding$sha256) &&
    identical(value$circuit, purpose) &&
    is.character(value$ciphertext) && length(value$ciphertext) == 1L &&
    !is.na(value$ciphertext) && nzchar(value$ciphertext) &&
    nchar(value$ciphertext, type = "bytes") <=
      .DSVERT_CLIENT_DP_COUNT_FINAL_CIPHERTEXT_MAX_BYTES &&
    identical(value$intermediate_values_exposed, FALSE) &&
    identical(value$capability_available, TRUE)
  if (!isTRUE(valid)) {
    stop("The Count garbler returned an invalid closed final share.",
         call. = FALSE)
  }
  transfer <- .dsvert_validate_typed_blob_transfer(
    value$transfer, value$ciphertext, evaluator)
  if (!identical(transfer$capability_id,
                 .DSVERT_CLIENT_DP_COUNT_FINAL_CAPABILITY) ||
      !identical(transfer$sender_name, garbler) ||
      !identical(transfer$recipient_name, evaluator)) {
    stop("The Count final share has an invalid typed capability route.",
         call. = FALSE)
  }
  value$transfer <- transfer
  value
}

.dsvert_dp_count_client_release_v1 <- function(
    value, contract, purpose) {
  fields <- c(
    "version", "artifact_key", "contract_sha256",
    "analysis_binding_sha256", "worker_static_sha256", "circuit",
    "mechanism", "bounds", "value", "source_identity_pk",
    "finalizer_identity_pk", "backend", "postprocessing",
    "intermediate_values_exposed", "public_openings", "release_sha256",
    "signature")
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      anyDuplicated(names(value)) || !setequal(names(value), fields)) {
    stop("The Count evaluator returned an invalid closed release.",
         call. = FALSE)
  }
  contract <- .dsvert_dp_analysis_contract_validate_v1(contract)
  roles <- .dsvert_dp_count_client_execution_roles_v1(contract)
  privacy <- contract$semantic$privacy
  mechanism <- privacy$mechanism
  bounds <- contract$semantic$analysis$effective_arguments$count_bounds
  upper <- .dsvert_dp_count_client_positive_integer_v1(
    bounds$upper, "release upper bound", 1000000)
  if (!is.character(value$value) || length(value$value) != 1L ||
      is.na(value$value) || !grepl("^(0|[1-9][0-9]*)$", value$value) ||
      nchar(value$value, type = "bytes") > 7L ||
      as.numeric(value$value) > upper) {
    stop("The Count release value violates its signed bounds.",
         call. = FALSE)
  }
  expected_core <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_COUNT_RELEASE_VERSION,
    artifact_key = contract$artifact_key,
    contract_sha256 = .dsvert_dp_count_client_hash_v1(
      .DSVERT_CLIENT_DP_COUNT_CONTRACT_DOMAIN, contract),
    analysis_binding_sha256 = roles$binding$sha256,
    worker_static_sha256 = .dsvert_dp_count_client_hash_v1(
      .DSVERT_CLIENT_DP_COUNT_WORKER_DOMAIN,
      .dsvert_dp_count_client_worker_static_v1(contract)),
    circuit = purpose,
    mechanism = list(
      family = mechanism$family,
      version = mechanism$version,
      sampler = mechanism$calibration$sampler,
      epsilon = privacy$epsilon,
      delta = privacy$delta,
      implementation_delta = mechanism$calibration$implementation_delta,
      sensitivity_l1 = 1),
    bounds = list(lower = "0", upper = as.character(as.integer(upper))),
    value = value$value,
    source_identity_pk = unname(roles$identities[["garbler"]]),
    finalizer_identity_pk = unname(roles$identities[["evaluator"]]),
    backend = "exact-gc-joint-dp-laplace-ring127-v2",
    postprocessing =
      "one-joint-noise-draw-and-one-clamp-inside-exact-gc",
    intermediate_values_exposed = FALSE,
    public_openings = 1))
  core <- .dsvert_dp_analysis_client_canonical_value_v1(
    value[setdiff(names(value), c("release_sha256", "signature"))])
  if (!identical(core, expected_core)) {
    stop("The Count evaluator returned a misbound release.", call. = FALSE)
  }
  release_sha256 <- .dsvert_dp_count_client_hex_v1(
    value$release_sha256, "release hash")
  expected_hash <- .dsvert_dp_count_client_hash_v1(
    .DSVERT_CLIENT_DP_COUNT_RELEASE_DOMAIN, core)
  if (!identical(release_sha256, expected_hash)) {
    stop("The Count evaluator returned a corrupt release.", call. = FALSE)
  }
  signed <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    core, list(release_sha256 = expected_hash)))
  public <- .dsvert_joint_dp_client_b64url(
    roles$identities[["evaluator"]], 32L,
    "Count release identity public key")
  signature <- .dsvert_joint_dp_client_b64url(
    value$signature, 64L, "Count release signature")
  verified <- tryCatch(openssl::ed25519_verify(
    charToRaw(paste0(
      .DSVERT_CLIENT_DP_COUNT_RELEASE_SIGNATURE_DOMAIN,
      .dsvert_joint_dp_client_json(signed))),
    signature, openssl::read_ed25519_pubkey(public)),
    error = function(error) FALSE)
  if (!isTRUE(verified)) {
    stop("Count release signature verification failed.", call. = FALSE)
  }
  .dsvert_dp_analysis_client_canonical_value_v1(value)
}

.dsvert_dp_count_execute_v1 <- function(
    data_name, datasources,
    .aggregate = DSI::datashield.aggregate,
    .prepare = .dsvert_dp_count_compile_authorize_bind_v1,
    .new_context = .dsvert_exact_gc_new_context,
    .run_exact = .dsvert_exact_gc_run,
    .store_typed = .dsvert_store_typed_blob,
    .cleanup = .dsvert_exact_gc_cleanup_best_effort,
    .abort = .dsvert_exact_gc_abort) {
  dependencies <- list(
    .aggregate, .prepare, .new_context, .run_exact, .store_typed,
    .cleanup, .abort)
  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) ||
      !grepl("^[a-zA-Z._][a-zA-Z0-9._]*$", data_name) ||
      any(!vapply(dependencies, is.function, logical(1L)))) {
    stop("Invalid Count execution request or dependency.", call. = FALSE)
  }
  prepared <- .prepare(
    data_name, datasources, .aggregate = .aggregate)
  fixed_candidate <- is.list(prepared) &&
    identical(prepared$version, .DSVERT_CLIENT_DP_COUNT_PREPARED_VERSION) &&
    identical(prepared$mode, "fixed_cohort_public")
  if (isTRUE(fixed_candidate)) {
    prepared <- .dsvert_dp_count_client_prepared_v1(prepared, datasources)
    return(list(
      version = .DSVERT_CLIENT_DP_COUNT_EXECUTION_RESULT_VERSION,
      mode = "fixed_cohort_public",
      payload = prepared$payload))
  }
  prepared_payload <- if (is.list(prepared)) prepared$payload else NULL
  cleanup_conns <- list()
  cleanup_session <- if (is.list(prepared_payload) &&
      is.character(prepared_payload$session_id) &&
      length(prepared_payload$session_id) == 1L &&
      !is.na(prepared_payload$session_id)) {
    prepared_payload$session_id
  } else {
    ""
  }
  cleanup_transport <- if (is.list(prepared_payload)) {
    prepared_payload$transport
  } else {
    NULL
  }
  cleanup_authorities <- if (is.list(cleanup_transport) &&
      !is.null(names(cleanup_transport)) &&
      length(names(cleanup_transport)) == 2L &&
      !anyNA(names(cleanup_transport)) &&
      !anyDuplicated(names(cleanup_transport))) {
    names(cleanup_transport)
  } else if (is.list(prepared_payload) &&
      is.character(prepared_payload$authorities) &&
      length(prepared_payload$authorities) == 2L &&
      !anyNA(prepared_payload$authorities) &&
      !anyDuplicated(prepared_payload$authorities)) {
    prepared_payload$authorities
  } else {
    character()
  }
  if (is.list(datasources) && !is.null(names(datasources)) &&
      length(cleanup_authorities) == 2L &&
      all(cleanup_authorities %in% names(datasources))) {
    cleanup_conns <- datasources[cleanup_authorities]
  }
  context <- NULL
  operation <- NULL
  released <- FALSE
  start_attempted <- FALSE
  on.exit({
    if (isTRUE(start_attempted) && !isTRUE(released) &&
        !is.null(context) && !is.null(operation)) {
      tryCatch(.abort(
        context$conns, context$session_id, operation$operation_id,
        .aggregate), error = function(error) NULL)
    }
    tryCatch(.cleanup(
      cleanup_conns, cleanup_session, cleanup_transport,
      .aggregate = .aggregate), error = function(error) NULL)
  }, add = TRUE)
  prepared <- .dsvert_dp_count_client_prepared_v1(prepared, datasources)
  context <- .dsvert_dp_count_client_execution_context_v1(
    prepared$payload, datasources)
  cleanup_conns <- context$conns
  cleanup_session <- context$session_id
  cleanup_transport <- prepared$payload$transport
  operation <- .dsvert_dp_count_client_operation_v1(.new_context())

  authorizations_json <-
    .dsvert_dp_count_client_authorizations_wire_v1(context)
  start_call <- call(
    name = "dsvertDPCountStartDS", data_name = data_name,
    session_id = context$session_id,
    operation_id = operation$operation_id,
    source_key = operation$source_key, output_key = operation$output_key,
    authorizations_json = .dsvert_dsi_text_encode(
      authorizations_json, "Count execution authorizations"))
  start_attempted <- TRUE
  initialized <- .dsvert_aggregate_strict(
    context$conns, start_call, operation = "Count exact-MPC Start",
    .aggregate = .aggregate)
  started <- .dsvert_dp_count_client_start_v1(initialized, context)
  selected <- match(context$authorities, names(datasources))
  .run_exact(
    datasources = datasources, server_names = names(datasources),
    servers = selected, session_id = context$session_id,
    operation_id = operation$operation_id,
    source_key = operation$source_key, output_key = operation$output_key,
    operation = "joint-dp-laplace-v2", ring = 127L, frac_bits = 0L,
    vector_len = 1L, purpose = started$purpose,
    transport_ready = TRUE, initialized = started$values,
    analysis_contract = context$contract, .aggregate = .aggregate)

  garbler <- unname(context$roles$peers[["garbler"]])
  evaluator <- unname(context$roles$peers[["evaluator"]])
  final_raw <- .dsvert_aggregate_strict(
    context$conns[garbler], call(
      name = "dsvertDPCountFinalShareDS",
      session_id = context$session_id,
      operation_id = operation$operation_id,
      output_key = operation$output_key,
      recipient_pk = unname(context$transport[[evaluator]])),
    operation = "Count encrypted final-share mint",
    .aggregate = .aggregate)[[garbler]]
  final <- .dsvert_dp_count_client_final_share_v1(
    final_raw, context, started$purpose)
  .store_typed(
    blob = final$ciphertext, transfer = final$transfer,
    conn = context$conns[evaluator], session_id = context$session_id,
    producer_conn = context$conns[garbler], .aggregate = .aggregate)
  release <- .dsvert_aggregate_strict(
    context$conns[evaluator], call(
      name = "dsvertDPCountReleaseDS",
      session_id = context$session_id,
      operation_id = operation$operation_id,
      output_key = operation$output_key),
    operation = "Count single public release",
    .aggregate = .aggregate)[[evaluator]]
  release <- .dsvert_dp_count_client_release_v1(
    release, context$contract, started$purpose)
  accuracy_95_abs <- .dsvert_dp_count_client_accuracy_95_v1(
    context$contract$semantic$analysis$effective_arguments$sampler_plan)
  released <- TRUE
  list(
    version = .DSVERT_CLIENT_DP_COUNT_EXECUTION_RESULT_VERSION,
    mode = "add_remove_dp",
    payload = list(
      release = release,
      finalizer_peer = evaluator,
      accuracy_95_abs = accuracy_95_abs,
      accuracy_95_confidence = 0.95,
      accuracy_95_method = .DSVERT_CLIENT_DP_COUNT_ACCURACY_95_METHOD))
}
