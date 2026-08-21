# Internal DSI pump for the exact two-peer garbled-circuit backend.  The client
# relays only signed, context-bound ciphertext envelopes and never receives a
# source share, output share, wire label, or peer-derived channel key.

.DSVERT_CLIENT_EXACT_GC_CAPABILITY <- "exact_gc_v1"
.DSVERT_CLIENT_EXACT_GC_PEER_RE <- "^dsv1_[0-9a-f]{64}$"
.DSVERT_CLIENT_EXACT_GC_MIN_RING_BITS <- 63L
.DSVERT_CLIENT_EXACT_GC_MAX_RING_BITS <- 4096L
.DSVERT_CLIENT_EXACT_GC_MAX_BOUND_DIGITS <- 1200L
.DSVERT_CLIENT_EXACT_GC_MAX_CIRCUIT_TYPE_BITS <- 512L * 1024L
.DSVERT_CLIENT_EXACT_GC_DIRECT_MUL_BIT_WORK <- 512L * 512L * 64L
.DSVERT_CLIENT_JOINT_DP_VECTOR_MAX_CHUNK <- 128L
.DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE <-
  "dp.cross-owner.exact-session.v1"
.DSVERT_CLIENT_EXACT_GC_FREQUENCY_CLEANUP_PURPOSE <-
  "dp.frequency.exact-session.v1"
.DSVERT_CLIENT_EXACT_GC_ANALYSIS_BINDING_VERSION <-
  "dsvert-exact-gc-analysis-binding-v1"
.DSVERT_CLIENT_EXACT_GC_FREQUENCY_BINDING_VERSION <-
  "dsvert-exact-gc-frequency-binding-v1"
.DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION <-
  "dsvert-exact-gc-cleanup-capability-v1"
.DSVERT_CLIENT_EXACT_GC_WIRE_CONTAINER_BITS <-
  c(64L, 128L, 256L, 512L, 1024L, 2048L, 4096L)

.dsvert_exact_gc_monotonic_seconds <- function() {
  unname(proc.time()[["elapsed"]])
}

.dsvert_exact_gc_sleep <- function(seconds) Sys.sleep(seconds)

.dsvert_exact_gc_new_context <- function(.random_bytes = .dsvert_random_bytes) {
  bytes <- .random_bytes(16L)
  if (!is.raw(bytes) || length(bytes) != 16L) {
    stop("Exact MPC context generation requires 16 secure random bytes.",
         call. = FALSE)
  }
  suffix <- paste0(sprintf("%02x", as.integer(bytes)), collapse = "")
  list(
    operation_id = paste0("op_", suffix),
    source_key = paste0("exact_gc_in_", suffix),
    output_key = paste0("exact_gc_out_", suffix))
}

.dsvert_exact_gc_b64url_encode <- function(value) {
  if (!is.raw(value)) stop("Exact MPC request must be raw.", call. = FALSE)
  encoded <- gsub("[\r\n]", "", jsonlite::base64_enc(value))
  encoded <- gsub("+", "-", encoded, fixed = TRUE)
  encoded <- gsub("/", "_", encoded, fixed = TRUE)
  gsub("=+$", "", encoded)
}

.dsvert_exact_gc_identity_peer_id <- function(identity_pk) {
  identity_pk <- .dsvert_dp_analysis_client_identity_pk(identity_pk)
  raw_pk <- tryCatch(
    jsonlite::base64_dec(chartr("-_", "+/", paste0(identity_pk, "="))),
    error = function(error) raw(0L))
  if (!is.raw(raw_pk) || length(raw_pk) != 32L) {
    stop("Invalid exact MPC analysis identity.", call. = FALSE)
  }
  paste0("dsv1_", digest::digest(
    c(charToRaw("dsVert/peer-capability/v1|"), raw_pk),
    algo = "sha256", serialize = FALSE))
}

.dsvert_exact_gc_analysis_binding <- function(contract) {
  contract <- .dsvert_dp_analysis_contract_validate_v1(contract)
  semantic <- contract$semantic
  if (!identical(semantic$analysis$primitive, "joint-dp-laplace-v2") ||
      !identical(contract$execution$backend$kernel,
                 "joint-dp-laplace-v2") ||
      !identical(contract$execution$backend$ring, "ring127")) {
    stop("The exact MPC analysis contract is not the audited Count shape.",
         call. = FALSE)
  }
  full_pins <- unlist(contract$execution$peer_pins, use.names = TRUE)
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", names(full_pins))) ||
      any(nchar(names(full_pins), type = "bytes") > 128L)) {
    stop("Invalid full K peer names in the exact MPC analysis contract.",
         call. = FALSE)
  }
  full_pins <- full_pins[order(names(full_pins), method = "radix")]
  authorities <- unlist(
    semantic$noise_authorities, use.names = FALSE)
  authority_names <- names(full_pins)[match(authorities, unname(full_pins))]
  authority_peer_ids <- vapply(
    authorities, .dsvert_exact_gc_identity_peer_id, character(1L),
    USE.NAMES = FALSE)
  authority_order <- order(authority_peer_ids, method = "radix")
  binding <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_EXACT_GC_ANALYSIS_BINDING_VERSION,
    artifact_key = contract$artifact_key,
    semantic_contract_sha256 = digest::digest(
      .dsvert_joint_dp_client_json(semantic), algo = "sha256",
      serialize = FALSE),
    authority_roles = as.list(stats::setNames(
      authorities[authority_order], c("garbler", "evaluator")))))
  list(
    contract = contract,
    full_pins = full_pins,
    authority_names = authority_names,
    binding = binding,
    sha256 = digest::digest(
      .dsvert_joint_dp_client_json(binding), algo = "sha256",
      serialize = FALSE))
}

.dsvert_exact_gc_frequency_analysis_binding <- function(contract) {
  contract <- .dsvert_dp_analysis_contract_validate_v1(contract)
  roles <- contract$semantic$noise_authority_roles
  role_names <- c("source_owner", "secondary_noise_authority")
  authority_ids <- unlist(roles$authority_ids, use.names = FALSE)
  pins <- unlist(contract$execution$peer_pins, use.names = TRUE)
  pins <- pins[order(names(pins), method = "radix")]
  authority_names <- names(pins)[match(authority_ids, unname(pins))]
  if (!identical(unlist(roles$role_order, use.names = FALSE), role_names) ||
      length(authority_ids) != 2L || anyNA(authority_names) ||
      anyDuplicated(authority_names) ||
      !identical(contract$execution$backend$ring, "ring128")) {
    stop("Invalid exact MPC Frequency analysis contract.", call. = FALSE)
  }
  binding <- list(
    version = "dsvert-dp-frequency-analysis-binding-v1",
    artifact_key = contract$artifact_key,
    semantic_contract_sha256 = .dsvert_dp_analysis_frequency_hash_v1(
      "dsVert/dp-frequency/semantic-contract/v1|", contract$semantic),
    authority_roles = as.list(stats::setNames(authority_ids, role_names)))
  list(contract = contract, full_pins = pins,
       authority_names = unname(authority_names), binding = binding,
       sha256 = .dsvert_dp_analysis_frequency_hash_v1(
         "dsVert/dp-frequency/analysis-binding/v1|", binding))
}

.dsvert_exact_gc_frequency_binding <- function(frequency_analysis) {
  required <- c("contract", "binding", "worker_static", "config_sha256",
    "source_claim_sha256", "receipt_set_sha256", "psi_run_sha256",
    "contract_sha256", "worker_static_sha256")
  if (!is.list(frequency_analysis) ||
      any(!required %in% names(frequency_analysis))) {
    stop("Invalid exact MPC Frequency preparation.", call. = FALSE)
  }
  analysis <- .dsvert_exact_gc_frequency_analysis_binding(
    frequency_analysis$contract)
  common <- c(
    artifact_key = analysis$contract$artifact_key,
    config_sha256 = frequency_analysis$config_sha256,
    source_claim_sha256 = frequency_analysis$source_claim_sha256,
    receipt_set_sha256 = frequency_analysis$receipt_set_sha256,
    psi_run_sha256 = frequency_analysis$psi_run_sha256,
    contract_sha256 = frequency_analysis$contract_sha256,
    analysis_binding_sha256 = analysis$sha256,
    worker_static_sha256 = frequency_analysis$worker_static_sha256)
  if (!identical(frequency_analysis$binding, analysis) ||
      any(!vapply(common, .dsvert_dp_analysis_frequency_hex_v1, logical(1L))) ||
      !.dsvert_dp_analysis_frequency_hex_v1(
        frequency_analysis$worker_static$release_contract_hash)) {
    stop("Invalid exact MPC Frequency preparation.", call. = FALSE)
  }
  binding <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    list(version = .DSVERT_CLIENT_EXACT_GC_FREQUENCY_BINDING_VERSION),
    as.list(common), list(
      authority_roles = analysis$binding$authority_roles,
      release_contract_hash =
        frequency_analysis$worker_static$release_contract_hash)))
  list(analysis = analysis, binding = binding, sha256 = digest::digest(
    .dsvert_joint_dp_client_json(binding), "sha256", serialize = FALSE))
}

.dsvert_setup_exact_gc_transport <- function(
    datasources, server_names, servers, session_id,
    cleanup_purpose = "",
    analysis_contract = NULL,
    frequency_analysis = NULL,
    .aggregate = DSI::datashield.aggregate) {
  if (!is.null(analysis_contract) && !is.null(frequency_analysis)) {
    stop("Exact MPC transport accepts one analysis binding.", call. = FALSE)
  }
  count_analysis <- if (is.null(analysis_contract)) NULL else
    .dsvert_exact_gc_analysis_binding(analysis_contract)
  frequency <- if (is.null(frequency_analysis)) NULL else
    .dsvert_exact_gc_frequency_binding(frequency_analysis)
  analysis <- if (is.null(frequency)) count_analysis else frequency$analysis
  if (!is.null(analysis) &&
      (length(server_names) != length(datasources) || anyNA(server_names) ||
       any(!nzchar(server_names)) || anyDuplicated(server_names) ||
       !setequal(server_names, names(analysis$full_pins)))) {
    stop("Exact MPC analysis requires the full K peer names and pins.",
         call. = FALSE)
  }
  if (!is.numeric(servers) || length(servers) != 2L || anyNA(servers) ||
      any(!is.finite(servers)) || any(servers != floor(servers)) ||
      anyDuplicated(servers) || any(servers < 1L) ||
      any(servers > length(datasources))) {
    stop("Exact MPC transport requires exactly two distinct compute peers.",
         call. = FALSE)
  }
  servers <- as.integer(servers)
  selected_names <- server_names[servers]
  if (length(selected_names) != 2L || anyNA(selected_names) ||
      any(!nzchar(selected_names)) || anyDuplicated(selected_names)) {
    stop("Exact MPC transport requires two unique federation names.",
         call. = FALSE)
  }
  if (!is.null(analysis)) {
    if (!setequal(selected_names, analysis$authority_names)) {
      stop("Exact MPC compute peers must be the two noise authorities.",
           call. = FALSE)
    }
  }
  conns <- datasources[servers]
  names(conns) <- selected_names
  if (!is.character(cleanup_purpose) || length(cleanup_purpose) != 1L ||
      is.na(cleanup_purpose) ||
      !cleanup_purpose %in% c(
        "", .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
        .DSVERT_CLIENT_EXACT_GC_FREQUENCY_CLEANUP_PURPOSE)) {
    stop("Invalid exact MPC cleanup capability purpose.", call. = FALSE)
  }
  initialized <- .dsvert_aggregate_strict(
    conns, call(name = "exactGCTransportInitDS", session_id = session_id),
    operation = "exact MPC transport initialization",
    .aggregate = .aggregate)
  transport <- list()
  identities <- list()
  for (server in selected_names) {
    value <- initialized[[server]]
    if (!is.list(value) ||
        !identical(value$capability_id,
                   .DSVERT_CLIENT_EXACT_GC_CAPABILITY) ||
        !is.character(value$transport_pk) ||
        !is.character(value$identity_pk) ||
        !is.character(value$signature)) {
      stop("An exact MPC peer returned an invalid signed transport key.",
           call. = FALSE)
    }
    if (!is.null(frequency)) {
      transport_raw <- .dsvert_joint_dp_client_b64url(
        value$transport_pk, 32L, "Frequency transport public key")
      value$transport_pk <- .dsvert_exact_gc_b64url_encode(transport_raw)
      value$identity_pk <- .dsvert_dp_analysis_client_identity_pk(
        value$identity_pk)
      .dsvert_dp_frequency_client_verify_v1(
        transport_raw, value$identity_pk, value$signature, "transport")
    }
    transport[[server]] <- value$transport_pk
    identities[[server]] <- list(
      identity_pk = value$identity_pk, signature = value$signature)
    if (!is.null(analysis) &&
        !identical(
          .dsvert_dp_analysis_client_identity_pk(value$identity_pk),
          unname(analysis$full_pins[[server]]))) {
      stop("An exact MPC peer identity disagrees with the full K pinset.",
           call. = FALSE)
    }
  }
  ordered <- sort(selected_names, method = "radix")
  frequency_peer_digest <- NULL
  if (!is.null(frequency)) {
    contract <- list(
      version = "dsvert-exact-gc-frequency-peer-binding-v1",
      capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
      session_id = session_id,
      consortium_id = analysis$contract$artifact_key,
      full_peer_pinset_sha256 = digest::digest(
        .dsvert_dp_frequency_client_wire_json_v1(
          as.list(analysis$full_pins)), "sha256", serialize = FALSE),
      designated_peers = as.list(ordered),
      designated_peer_pinset = as.list(analysis$full_pins[ordered]),
      identity_pks = lapply(identities[ordered], `[[`, "identity_pk"),
      transport_pks = as.list(unlist(transport[ordered], use.names = TRUE)),
      frequency_binding = frequency$binding,
      frequency_binding_sha256 = frequency$sha256)
    frequency_peer_digest <- digest::digest(
      .dsvert_dp_frequency_client_wire_json_v1(contract), "sha256",
      serialize = FALSE)
  }
  encode_map <- function(value) {
    .dsvert_exact_gc_b64url_encode(charToRaw(as.character(jsonlite::toJSON(
      value[ordered], auto_unbox = TRUE, null = "null", digits = NA))))
  }
  bind_call <- if (is.null(analysis) || !is.null(frequency)) {
    expression <- call(
      name = "exactGCBindPeersDS",
      transport_keys_b64 = encode_map(transport),
      identity_info_b64 = encode_map(identities), session_id = session_id)
    if (nzchar(cleanup_purpose)) {
      expression$cleanup_purpose <- cleanup_purpose
    }
    expression
  } else {
    expression <- call(
      name = "exactGCBindPeersDS",
      transport_keys_b64 = encode_map(transport),
      identity_info_b64 = encode_map(identities), session_id = session_id,
      artifact_key = analysis$contract$artifact_key)
    if (nzchar(cleanup_purpose)) {
      expression$cleanup_purpose <- cleanup_purpose
    }
    expression
  }
  bound <- .dsvert_aggregate_strict(
    conns, bind_call,
    operation = "exact MPC pinned-peer binding",
    .aggregate = .aggregate)
  cleanup_capabilities <- list()
  for (server in selected_names) {
    value <- bound[[server]]
    if (!is.list(value) || !isTRUE(value$bound) ||
        !identical(value$capability_id,
                   .DSVERT_CLIENT_EXACT_GC_CAPABILITY)) {
      stop("An exact MPC peer did not commit its pinned peer binding.",
           call. = FALSE)
    }
    if (!is.null(analysis) &&
        ((is.null(frequency) &&
          (!identical(value$analysis_binding, analysis$binding) ||
           !identical(value$analysis_binding_sha256, analysis$sha256))) ||
         (!is.null(frequency) &&
          (!is.null(value$analysis_binding) ||
           !identical(value$frequency_binding, frequency$binding) ||
           !identical(value$frequency_binding_sha256, frequency$sha256))))) {
      stop("An exact MPC peer returned a conflicting analysis binding.",
           call. = FALSE)
    }
    if (nzchar(cleanup_purpose)) {
      capability <- tryCatch(jsonlite::fromJSON(
        value$cleanup_capability_json, simplifyVector = FALSE),
        error = function(error) NULL)
      contract <- if (is.list(capability)) capability$contract else NULL
      canonical <- if (is.list(capability)) tryCatch(
        as.character(jsonlite::toJSON(
          capability, auto_unbox = TRUE, null = "null", na = "null",
          digits = 17, pretty = FALSE)),
        error = function(error) NULL) else NULL
      if (!identical(value$cleanup_purpose, cleanup_purpose) ||
          !is.character(value$cleanup_capability_json) ||
          length(value$cleanup_capability_json) != 1L ||
          nchar(value$cleanup_capability_json, type = "bytes") > 16384L ||
          !is.list(capability) ||
          !identical(names(capability),
                     c("version", "contract", "signature")) ||
          !identical(canonical, value$cleanup_capability_json) ||
          !identical(capability$version,
                     .DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION) ||
          !is.list(contract) ||
          !identical(names(contract), c(
            "version", "session_id", "cleanup_purpose",
            "operation_scope", "peer_binding_digest")) ||
          !identical(contract$version,
                     .DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION) ||
          !identical(contract$session_id, session_id) ||
          !identical(contract$cleanup_purpose, cleanup_purpose) ||
          !identical(contract$operation_scope,
                     "all_and_only_operations_in_bound_exact_session_v1") ||
          !is.character(contract$peer_binding_digest) ||
          length(contract$peer_binding_digest) != 1L ||
          !grepl("^[0-9a-f]{64}$", contract$peer_binding_digest) ||
          !is.character(capability$signature) ||
          length(capability$signature) != 1L ||
          !grepl("^[A-Za-z0-9_-]{86}$", capability$signature)) {
        stop("An exact MPC peer returned an invalid cleanup capability.",
             call. = FALSE)
      }
      if (!is.null(frequency)) {
        if (!identical(contract$peer_binding_digest,
                       frequency_peer_digest)) {
          stop("An exact MPC peer returned a misbound cleanup capability.",
               call. = FALSE)
        }
        .dsvert_dp_frequency_client_verify_v1(
          charToRaw(paste0("dsVert/exact-gc/cleanup-capability/v1|",
            .dsvert_dp_frequency_client_wire_json_v1(contract))),
          unname(analysis$full_pins[[server]]), capability$signature,
          "cleanup capability")
      }
      cleanup_capabilities[[server]] <- value$cleanup_capability_json
    }
  }
  if (nzchar(cleanup_purpose)) {
    attr(transport, "exact_gc_cleanup_capabilities") <-
      cleanup_capabilities[selected_names]
    attr(transport, "exact_gc_cleanup_purpose") <- cleanup_purpose
  }
  if (!is.null(analysis)) {
    prefix <- if (is.null(frequency)) "exact_gc_analysis_binding" else
      "exact_gc_frequency_binding"
    binding <- if (is.null(frequency)) analysis else frequency
    attr(transport, prefix) <- binding$binding
    attr(transport, paste0(prefix, "_sha256")) <- binding$sha256
  }
  if (!is.null(frequency)) {
    attr(transport, "exact_gc_peer_binding_digest") <- frequency_peer_digest
  }
  transport
}

.dsvert_setup_frequency_transport <- function(
    datasources, server_names, servers, session_id, frequency_analysis,
    .aggregate = DSI::datashield.aggregate) {
  .dsvert_setup_exact_gc_transport(
    datasources, server_names, servers, session_id,
    cleanup_purpose = .DSVERT_CLIENT_EXACT_GC_FREQUENCY_CLEANUP_PURPOSE,
    frequency_analysis = frequency_analysis, .aggregate = .aggregate)
}

.dsvert_exact_gc_cleanup_best_effort <- function(
    conns, session_id, setup_result,
    .aggregate = DSI::datashield.aggregate) {
  capabilities <- attr(
    setup_result, "exact_gc_cleanup_capabilities", exact = TRUE)
  purpose <- attr(setup_result, "exact_gc_cleanup_purpose", exact = TRUE)
  if (!is.list(capabilities) || !length(capabilities) ||
      is.null(names(capabilities)) || anyNA(names(capabilities)) ||
      any(!nzchar(names(capabilities))) || anyDuplicated(names(capabilities)) ||
      !purpose %in% c(
        .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
        .DSVERT_CLIENT_EXACT_GC_FREQUENCY_CLEANUP_PURPOSE) ||
      !is.list(conns) || !setequal(names(conns), names(capabilities))) {
    return(invisible(FALSE))
  }
  selected <- conns[names(capabilities)]
  expressions <- stats::setNames(lapply(names(capabilities), function(site) {
    capability <- capabilities[[site]]
    if (!is.character(capability) || length(capability) != 1L ||
        is.na(capability) || !nzchar(capability)) {
      return(NULL)
    }
    call(
      name = "exactGCCleanupDS", session_id = session_id,
      cleanup_capability_json = capability)
  }), names(capabilities))
  if (any(!vapply(expressions, is.call, logical(1L)))) {
    return(invisible(FALSE))
  }
  .dsvert_cleanup_best_effort(
    selected, expressions, .aggregate = .aggregate)
}

.dsvert_exact_gc_number <- function(value, what, minimum = 0,
                                    maximum = 2^53) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value != floor(value) || value < minimum || value > maximum) {
    stop("Invalid ", what, " returned by an exact MPC peer.", call. = FALSE)
  }
  value
}

.dsvert_exact_gc_container_bits <- function(ring_bits) {
  ring_bits <- as.integer(.dsvert_exact_gc_number(
    ring_bits, "ring", .DSVERT_CLIENT_EXACT_GC_MIN_RING_BITS,
    .DSVERT_CLIENT_EXACT_GC_MAX_RING_BITS))
  .DSVERT_CLIENT_EXACT_GC_WIRE_CONTAINER_BITS[
    which(.DSVERT_CLIENT_EXACT_GC_WIRE_CONTAINER_BITS >= ring_bits)[[1L]]]
}

.dsvert_exact_gc_direct_mul_max_chunk <- function(ring_bits) {
  container_bits <- .dsvert_exact_gc_container_bits(ring_bits)
  as.integer(min(
    64,
    floor((.DSVERT_CLIENT_EXACT_GC_MAX_CIRCUIT_TYPE_BITS /
             container_bits - 1) / 7),
    floor(.DSVERT_CLIENT_EXACT_GC_DIRECT_MUL_BIT_WORK /
            container_bits^2)))
}

.dsvert_exact_gc_peer_id <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl(.DSVERT_CLIENT_EXACT_GC_PEER_RE, value)) {
    stop("An exact MPC peer returned an invalid pinned identity.",
         call. = FALSE)
  }
  value
}

.dsvert_exact_gc_validate_init <- function(results, server_names,
                                           operation, ring, frac_bits,
                                           vector_len, purpose,
                                           mul_contract = NULL,
                                           alignment_source_count = NULL,
                                           analysis_binding = NULL) {
  analysis <- NULL
  if (!is.null(analysis_binding)) {
    if (!is.list(analysis_binding) || is.null(analysis_binding$contract)) {
      stop("Invalid exact MPC analysis binding.", call. = FALSE)
    }
    analysis <- .dsvert_exact_gc_analysis_binding(
      analysis_binding$contract)
    if (!identical(analysis_binding$binding, analysis$binding) ||
        !identical(analysis_binding$sha256, analysis$sha256)) {
      stop("Invalid exact MPC analysis binding.", call. = FALSE)
    }
  }
  if (identical(operation, "joint-dp-laplace-v2")) {
    if (is.null(analysis) || ring != 127L || frac_bits != 0L ||
        vector_len != 1L) {
      stop("Scalar Count initialization requires its exact analysis binding.",
           call. = FALSE)
    }
  } else if (!is.null(analysis)) {
    stop("An exact MPC analysis binding is valid only for scalar Count.",
         call. = FALSE)
  }
  if (!is.null(analysis)) {
    if (!setequal(server_names, analysis$authority_names)) {
      stop("Scalar Count peers do not match the analysis noise authorities.",
           call. = FALSE)
    }
  }
  if (!is.list(results) || !all(server_names %in% names(results))) {
    stop("Exact MPC peers did not return initialization state.", call. = FALSE)
  }
  states <- results[server_names]
  required <- c(
    "capability_id", "peer_id", "peer_peer_id", "role", "context_hash",
    "operation", "output_kind", "purpose", "source_producer", "ring_bits",
    "frac_bits", "vector_len", "threshold", "chunk_bytes", "state",
    "stored", "worker_heartbeat")
  if (!is.null(analysis)) {
    required <- c(required, "analysis_binding_sha256")
  }
  expected_kind <- if (identical(
      operation, "joint-dp-vector-laplace-v3")) {
    "joint-dp-vector-ring128-share-v1"
  } else if (identical(operation, "joint-dp-laplace-v2")) {
    "joint-dp-ring-share-v2"
  } else if (identical(operation, "formal-glm-phase19-schedule-v1")) {
    "formal-glm-phase19-ring128-dp-bridge-share-v1"
  } else if (identical(operation, "alignment-mask-ring128")) {
    "alignment-masked-ring128-share-v1"
  } else if (operation %in% c(
      "compare-signed", "truncate-floor", "clamp-count"))
    "ring-share" else if (identical(operation, "mul-truncate-checked"))
      "checked-ring-share" else "xor-bit-share"
  for (server in server_names) {
    value <- states[[server]]
    if (!is.list(value) || !all(required %in% names(value)) ||
        !identical(value$capability_id,
                   .DSVERT_CLIENT_EXACT_GC_CAPABILITY) ||
        !identical(value$operation, operation) ||
        !identical(value$output_kind, expected_kind) ||
        !identical(value$purpose, purpose) ||
        !is.character(value$source_producer) ||
        length(value$source_producer) != 1L ||
        !grepl("^[a-z][a-z0-9_.:/-]*$", value$source_producer) ||
        .dsvert_exact_gc_number(value$ring_bits, "ring") != ring ||
        .dsvert_exact_gc_number(value$frac_bits, "fractional-bit count") !=
          frac_bits ||
        .dsvert_exact_gc_number(value$vector_len, "vector length", 1, 4096) !=
          vector_len ||
        !identical(value$state, "running") || isTRUE(value$stored)) {
      stop("An exact MPC peer rejected or changed the operation contract.",
           call. = FALSE)
    }
    if (!is.null(analysis) &&
        !identical(value$analysis_binding_sha256, analysis$sha256)) {
      stop("An exact MPC peer returned a conflicting analysis binding.",
           call. = FALSE)
    }
    .dsvert_exact_gc_peer_id(value$peer_id)
    .dsvert_exact_gc_peer_id(value$peer_peer_id)
    if (!is.character(value$context_hash) || length(value$context_hash) != 1L ||
        !grepl("^[0-9a-f]{64}$", value$context_hash) ||
        !value$role %in% c("garbler", "evaluator")) {
      stop("An exact MPC peer returned invalid protocol context.",
           call. = FALSE)
    }
    if (identical(operation, "mul-truncate-checked")) {
      if (!is.list(mul_contract) ||
          !identical(value$mul_plan_id, mul_contract$plan_id) ||
          !identical(value$mul_backend, mul_contract$backend) ||
          !identical(value$bound_x, mul_contract$bound_x) ||
          !identical(value$bound_y, mul_contract$bound_y)) {
        stop("An exact MPC peer changed the multiplication plan.",
             call. = FALSE)
      }
    }
    .dsvert_exact_gc_number(value$chunk_bytes, "chunk size",
                            16 * 1024, 8 * 1024^2)
    .dsvert_exact_gc_number(
      value$worker_heartbeat, "worker heartbeat", 1, 2^53 - 1)
  }
  peer_ids <- vapply(states, `[[`, character(1L), "peer_id")
  peer_peers <- vapply(states, `[[`, character(1L), "peer_peer_id")
  roles <- vapply(states, `[[`, character(1L), "role")
  contexts <- vapply(states, `[[`, character(1L), "context_hash")
  thresholds <- vapply(states, function(value) as.character(value$threshold),
                       character(1L))
  chunks <- vapply(states, function(value) as.numeric(value$chunk_bytes),
                   numeric(1L))
  # Current servers publish their active inactivity lease.  Treat omission by
  # an older compatible server as the current conservative 180-second default;
  # the client always adopts the shortest peer lease.
  ttls <- vapply(states, function(value) {
    if (is.null(value$ttl_seconds)) return(180)
    .dsvert_exact_gc_number(
      value$ttl_seconds, "inactivity lease", 10, 86400)
  }, numeric(1L))
  max_runtimes <- vapply(seq_along(states), function(index) {
    value <- states[[index]]
    runtime <- if (is.null(value$max_runtime_seconds)) 21600 else
      .dsvert_exact_gc_number(
        value$max_runtime_seconds, "total runtime lease", 10, 7 * 86400)
    if (runtime < ttls[[index]]) {
      stop("An exact MPC peer returned an invalid total runtime lease.",
           call. = FALSE)
    }
    runtime
  }, numeric(1L))
  heartbeats <- vapply(states, function(value) {
    .dsvert_exact_gc_number(
      value$worker_heartbeat, "worker heartbeat", 1, 2^53 - 1)
  }, numeric(1L))
  producers <- vapply(states, `[[`, character(1L), "source_producer")
  if (!is.null(analysis)) {
    expected_peer_ids <- vapply(
      analysis$full_pins[server_names],
      .dsvert_exact_gc_identity_peer_id, character(1L), USE.NAMES = TRUE)
    expected_roles <- stats::setNames(vapply(
      unname(analysis$full_pins[server_names]), function(identity_pk) {
        role <- names(analysis$binding$authority_roles)[match(
          identity_pk, unlist(
            analysis$binding$authority_roles, use.names = FALSE))]
        unname(role)
      }, character(1L)), server_names)
    if (!identical(peer_ids, expected_peer_ids) ||
        !identical(roles, expected_roles)) {
      stop("Exact MPC peer roles disagree with the analysis binding.",
           call. = FALSE)
    }
  }
  if (anyDuplicated(peer_ids) ||
      !identical(unname(peer_peers), unname(rev(peer_ids))) ||
      !setequal(roles, c("garbler", "evaluator")) ||
      length(unique(contexts)) != 1L || length(unique(thresholds)) != 1L ||
      length(unique(chunks)) != 1L || length(unique(producers)) != 1L) {
    stop("Exact MPC peers did not agree on one pinned protocol context.",
         call. = FALSE)
  }
  if (identical(operation, "alignment-mask-ring128")) {
    expected_k <- .dsvert_exact_gc_number(
      alignment_source_count, "alignment source count", 2, 64)
    if (!identical(thresholds[[1L]], as.character(expected_k))) {
      stop("Exact MPC peers changed the alignment source count.",
           call. = FALSE)
    }
  } else if (!is.null(alignment_source_count)) {
    stop("Unexpected exact MPC alignment source count.", call. = FALSE)
  }
  result <- list(
       states = states, peer_ids = peer_ids, context_hash = contexts[[1L]],
       threshold = thresholds[[1L]], chunk_bytes = as.integer(chunks[[1L]]),
       ttl_seconds = as.integer(min(ttls)),
       max_runtime_seconds = as.integer(min(max_runtimes)),
       worker_heartbeats = heartbeats)
  if (!is.null(analysis)) {
    result$analysis_binding_sha256 <- analysis$sha256
  }
  result
}

.dsvert_exact_gc_validate_exchange <- function(value, expected_peer_id) {
  required <- c("capability_id", "peer_id", "state", "stored",
                "inbound_size", "outbound", "worker_heartbeat")
  if (!is.list(value) || !all(required %in% names(value)) ||
      !identical(value$capability_id,
                 .DSVERT_CLIENT_EXACT_GC_CAPABILITY) ||
      !identical(.dsvert_exact_gc_peer_id(value$peer_id), expected_peer_id) ||
      !value$state %in% c("running", "complete") ||
      !is.logical(value$stored) || length(value$stored) != 1L ||
      is.na(value$stored) ||
      !identical(isTRUE(value$stored), identical(value$state, "complete"))) {
    stop("An exact MPC peer returned invalid exchange state.", call. = FALSE)
  }
  value$inbound_size <- .dsvert_exact_gc_number(
    value$inbound_size, "inbound acknowledgment")
  value$worker_heartbeat <- .dsvert_exact_gc_number(
    value$worker_heartbeat, "worker heartbeat", 1, 2^53 - 1)
  value
}

.dsvert_exact_gc_validate_envelope <- function(envelope, source_peer,
                                               target_peer, session_id,
                                               operation_id, context_hash,
                                               chunk_bytes) {
  required <- c(
    "version", "capability_id", "session_id", "operation_id",
    "context_hash", "sender_peer_id", "recipient_peer_id", "offset",
    "chunk_bytes", "payload_hash", "payload", "signature")
  if (!is.list(envelope) || !identical(sort(names(envelope)), sort(required)) ||
      !identical(envelope$version, "dsvert-exact-gc-envelope-v1") ||
      !identical(envelope$capability_id,
                 .DSVERT_CLIENT_EXACT_GC_CAPABILITY) ||
      !identical(envelope$session_id, session_id) ||
      !identical(envelope$operation_id, operation_id) ||
      !identical(envelope$context_hash, context_hash) ||
      !identical(envelope$sender_peer_id, source_peer) ||
      !identical(envelope$recipient_peer_id, target_peer) ||
      !is.character(envelope$payload_hash) ||
      !grepl("^[0-9a-f]{64}$", envelope$payload_hash) ||
      !is.character(envelope$payload) || length(envelope$payload) != 1L ||
      !is.character(envelope$signature) || length(envelope$signature) != 1L) {
    stop("An exact MPC peer returned a malformed opaque envelope.",
         call. = FALSE)
  }
  envelope$offset <- .dsvert_exact_gc_number(
    envelope$offset, "envelope offset")
  envelope$chunk_bytes <- .dsvert_exact_gc_number(
    envelope$chunk_bytes, "envelope length", 1, chunk_bytes)
  max_encoded <- 4 * ceiling(envelope$chunk_bytes / 3) + 4
  if (nchar(envelope$payload, type = "bytes") > max_encoded ||
      !grepl("^[A-Za-z0-9_-]+$", envelope$payload)) {
    stop("An exact MPC peer returned an oversized opaque envelope.",
         call. = FALSE)
  }
  envelope
}

.dsvert_exact_gc_abort <- function(datasources, session_id, operation_id,
                                   .aggregate) {
  tryCatch(
    .dsvert_aggregate_strict(
      datasources, call(
        name = "exactGCAbortDS", session_id = session_id,
        operation_id = operation_id),
      operation = "exact MPC cleanup", result_contract = "logical_true",
      .aggregate = .aggregate), error = function(e) NULL)
  invisible(NULL)
}

.dsvert_exact_gc_delivery_fields <- function(envelope = NULL) {
  if (is.null(envelope)) {
    return(list(
      delivery_offset = 0,
      delivery_chunk_bytes = 0,
      delivery_payload_hash = "",
      delivery_payload = "",
      delivery_signature = ""))
  }
  required <- c(
    "version", "capability_id", "session_id", "operation_id",
    "context_hash", "sender_peer_id", "recipient_peer_id", "offset",
    "chunk_bytes", "payload_hash", "payload", "signature")
  if (!is.list(envelope) || !identical(sort(names(envelope)), sort(required))) {
    stop("An exact MPC peer returned a malformed opaque envelope.",
         call. = FALSE)
  }
  list(
    delivery_offset = envelope$offset,
    delivery_chunk_bytes = envelope$chunk_bytes,
    delivery_payload_hash = envelope$payload_hash,
    delivery_payload = envelope$payload,
    delivery_signature = envelope$signature)
}

.dsvert_exact_gc_vecmul_chunk_operation <- function(
    batch_operation_id, chunk_index, chunk_count, policy_id, plan_id) {
  digest <- digest::digest(charToRaw(paste0(
    "dsVert/exact-gc/checked-mul-chunk/v3|", batch_operation_id, "|",
    policy_id, "|", plan_id, "|", as.integer(chunk_index), "|",
    as.integer(chunk_count))),
    algo = "sha256", serialize = FALSE)
  paste0("op_", substr(digest, 1L, 32L))
}

.dsvert_exact_gc_vecmul_keys <- function(operation_id) {
  suffix <- sub("^op_", "", operation_id)
  list(source = paste0("exact_gc_in_", suffix),
       output = paste0("exact_gc_out_", suffix),
       destination = paste0("k2_exact_vecmul_z_", suffix),
       x = paste0("k2_exact_vecmul_x_", suffix),
       y = paste0("k2_exact_vecmul_y_", suffix))
}

# Direct checked Ring127 multiplication with exact floor truncation.
.dsvert_exact_gc_vecmul_run <- function(
    datasources, server_names = names(datasources),
    servers = seq_along(datasources), session_id, total_n,
    x_key = "k2_beaver_x", y_key = "k2_beaver_y", output_key = NULL,
    input_manifests = NULL,
    transport_ready = TRUE,
    timeout_seconds = getOption("dsvert.exact_gc.timeout_seconds", 900),
    .aggregate = DSI::datashield.aggregate) {
  if (length(servers) != 2L || anyDuplicated(servers)) {
    stop("Exact MPC multiplication requires exactly two peers.",
         call. = FALSE)
  }
  total_n <- suppressWarnings(as.numeric(total_n))
  if (length(total_n) != 1L || is.na(total_n) || !is.finite(total_n) ||
      total_n != floor(total_n) || total_n < 1 || total_n > 2^31 - 1) {
    stop("Invalid exact MPC multiplication contract.", call. = FALSE)
  }
  total_n <- as.integer(total_n)
  conns <- datasources[servers]
  selected_names <- server_names[servers]
  names(conns) <- selected_names
  if (!isTRUE(transport_ready)) {
    .dsvert_setup_exact_gc_transport(
      datasources, server_names, servers, session_id,
      .aggregate = .aggregate)
  }
  batch <- .dsvert_exact_gc_new_context()$operation_id
  batch_committed <- FALSE
  cleanup_operations <- batch
  on.exit(if (!batch_committed) {
    for (operation in rev(unique(cleanup_operations))) {
      .dsvert_exact_gc_abort(conns, session_id, operation, .aggregate)
    }
  }, add = TRUE)
  if (is.null(output_key)) {
    output_key <- .dsvert_exact_gc_vecmul_keys(batch)$destination
  }
  if (is.null(input_manifests)) {
    .dsvert_block_retired_remote_route("legacy_glm")
    bound <- .dsvert_aggregate_strict(
      conns, call(
        name = "exactGCVecmulBindInputsDS",
        x_key = x_key, y_key = y_key, output_key = output_key,
        total_n = total_n, batch_operation_id = batch,
        session_id = session_id),
      operation = "exact MPC legacy input binding", .aggregate = .aggregate)
    expected_state <- "bound"
  } else {
    if (!is.list(input_manifests) || is.null(names(input_manifests)) ||
        anyDuplicated(names(input_manifests)) ||
        !setequal(names(input_manifests), selected_names)) {
      stop("Exact MPC multiplication requires one producer manifest per peer.",
           call. = FALSE)
    }
    requests <- stats::setNames(lapply(selected_names, function(server) {
      manifest <- input_manifests[[server]]
      if (!is.list(manifest) ||
          !is.character(manifest$manifest_handle) ||
          length(manifest$manifest_handle) != 1L ||
          !grepl("^[A-Za-z0-9_-]{43}$", manifest$manifest_handle) ||
          !identical(as.integer(manifest$total_n), total_n)) {
        stop("Exact MPC multiplication received an invalid producer manifest.",
             call. = FALSE)
      }
      call(
        name = "exactGCVecmulClaimInputsDS",
        manifest_handle = manifest$manifest_handle,
        batch_operation_id = batch, session_id = session_id)
    }), selected_names)
    bound <- .dsvert_fanout_by_site(
      conns, requests, operation = "exact MPC producer-manifest claim",
      .aggregate = .aggregate)
    expected_state <- "claimed"
  }
  contexts <- policies <- plans <- backends <- bounds_x <- bounds_y <-
    character()
  ring_bits <- frac_bits <- max_chunks <- integer()
  for (server in selected_names) {
    value <- bound[[server]]
    required <- c(
      "stored", "state", "capability_id", "context_hash", "policy_id",
      "plan_id", "backend", "bound_x", "bound_y", "ring_bits",
      "frac_bits", "max_chunk")
    if (!is.list(value) || !all(required %in% names(value)) ||
        !isTRUE(value$stored) ||
        !identical(value$state, expected_state) ||
        !identical(value$capability_id,
                   .DSVERT_CLIENT_EXACT_GC_CAPABILITY) ||
        !is.character(value$context_hash) || length(value$context_hash) != 1L ||
        is.na(value$context_hash) ||
        !grepl("^[0-9a-f]{64}$", value$context_hash) ||
        !is.character(value$policy_id) || length(value$policy_id) != 1L ||
        is.na(value$policy_id) ||
        !grepl("^[0-9a-f]{64}$", value$policy_id) ||
        !is.character(value$plan_id) || length(value$plan_id) != 1L ||
        is.na(value$plan_id) ||
        !grepl("^[0-9a-f]{64}$", value$plan_id) ||
        !is.character(value$backend) || length(value$backend) != 1L ||
        is.na(value$backend) ||
        !value$backend %in% c("ring127-ot", "direct-wide") ||
        !is.character(value$bound_x) || length(value$bound_x) != 1L ||
        is.na(value$bound_x) ||
        nchar(value$bound_x, type = "bytes") >
          .DSVERT_CLIENT_EXACT_GC_MAX_BOUND_DIGITS ||
        !grepl("^[1-9][0-9]*$", value$bound_x) ||
        !is.character(value$bound_y) || length(value$bound_y) != 1L ||
        is.na(value$bound_y) ||
        nchar(value$bound_y, type = "bytes") >
          .DSVERT_CLIENT_EXACT_GC_MAX_BOUND_DIGITS ||
        !grepl("^[1-9][0-9]*$", value$bound_y)) {
      stop("Exact MPC multiplication failed.", call. = FALSE)
    }
    parsed <- tryCatch(list(
      ring_bits = as.integer(.dsvert_exact_gc_number(
        value$ring_bits, "multiplication ring",
        .DSVERT_CLIENT_EXACT_GC_MIN_RING_BITS,
        .DSVERT_CLIENT_EXACT_GC_MAX_RING_BITS)),
      frac_bits = as.integer(.dsvert_exact_gc_number(
        value$frac_bits, "multiplication fractional-bit count", 0,
        .DSVERT_CLIENT_EXACT_GC_MAX_RING_BITS - 1L)),
      max_chunk = as.integer(.dsvert_exact_gc_number(
        value$max_chunk, "multiplication chunk size", 1, 256))),
      error = function(e) NULL)
    if (is.null(parsed) || parsed$frac_bits >= parsed$ring_bits) {
      stop("Exact MPC multiplication failed.", call. = FALSE)
    }
    expected_max_chunk <- if (identical(value$backend, "ring127-ot")) {
      if (parsed$ring_bits != 127L || parsed$frac_bits != 50L) {
        stop("Exact MPC multiplication failed.", call. = FALSE)
      }
      256L
    } else {
      .dsvert_exact_gc_direct_mul_max_chunk(parsed$ring_bits)
    }
    if (!identical(parsed$max_chunk, expected_max_chunk)) {
      stop("Exact MPC multiplication failed.", call. = FALSE)
    }
    contexts <- c(contexts, value$context_hash)
    policies <- c(policies, value$policy_id)
    plans <- c(plans, value$plan_id)
    backends <- c(backends, value$backend)
    bounds_x <- c(bounds_x, value$bound_x)
    bounds_y <- c(bounds_y, value$bound_y)
    ring_bits <- c(ring_bits, parsed$ring_bits)
    frac_bits <- c(frac_bits, parsed$frac_bits)
    max_chunks <- c(max_chunks, parsed$max_chunk)
  }
  if (length(unique(contexts)) != 1L || length(unique(policies)) != 1L ||
      length(unique(plans)) != 1L || length(unique(backends)) != 1L ||
      length(unique(bounds_x)) != 1L || length(unique(bounds_y)) != 1L ||
      length(unique(ring_bits)) != 1L || length(unique(frac_bits)) != 1L ||
      length(unique(max_chunks)) != 1L) {
    stop("Exact MPC multiplication failed.", call. = FALSE)
  }
  policy_id <- policies[[1L]]
  plan <- list(
    plan_id = plans[[1L]], backend = backends[[1L]],
    bound_x = bounds_x[[1L]], bound_y = bounds_y[[1L]],
    ring_bits = ring_bits[[1L]], frac_bits = frac_bits[[1L]],
    max_chunk = max_chunks[[1L]])
  chunk_size <- plan$max_chunk
  chunk_count <- as.integer(ceiling(total_n / chunk_size))
  for (index in seq_len(chunk_count)) {
    offset <- (index - 1L) * chunk_size
    n <- as.integer(min(chunk_size, total_n - offset))
    operation_id <- .dsvert_exact_gc_vecmul_chunk_operation(
      batch, index, chunk_count, policy_id, plan$plan_id)
    cleanup_operations <- c(cleanup_operations, operation_id)
    common <- list(
      n = n, total_n = total_n, chunk_index = as.integer(index),
      chunk_count = chunk_count, batch_operation_id = batch,
      session_id = session_id, operation_id = operation_id,
      policy_id = policy_id, plan_id = plan$plan_id)
    started <- .dsvert_aggregate_strict(
      conns, as.call(c(list(as.name("exactGCVecmulStartDS")), common)),
      operation = "exact MPC checked multiplication start",
      .aggregate = .aggregate)
    keys <- .dsvert_exact_gc_vecmul_keys(operation_id)
    purpose <- paste0(
      "k2.vecmul.mul-truncate.v3.p-", policy_id,
      ".m-", substr(plan$plan_id, 1L, 16L),
      ".c-", index, "-", chunk_count)
    .dsvert_exact_gc_run(
      conns, server_names = selected_names, servers = 1:2,
      session_id = session_id, operation_id = operation_id,
      source_key = keys$source, output_key = keys$output,
      operation = "mul-truncate-checked", ring = plan$ring_bits,
      frac_bits = plan$frac_bits,
      vector_len = n, purpose = purpose, transport_ready = TRUE,
      timeout_seconds = timeout_seconds, initialized = started,
      mul_contract = plan,
      .aggregate = .aggregate)
    validity <- .dsvert_aggregate_strict(
      conns, as.call(c(list(as.name("exactGCVecmulValidityDS")), common)),
      operation = "exact MPC validity sealing", .aggregate = .aggregate)
    for (server in selected_names) {
      value <- validity[[server]]
      if (!is.list(value) || !identical(value$state, "sealed") ||
          !identical(value$capability_id,
                     .DSVERT_CLIENT_EXACT_GC_CAPABILITY) ||
          !is.character(value$peer_blob) ||
          !grepl("^[A-Za-z0-9_-]+$", value$peer_blob)) {
        stop("Exact MPC multiplication failed.", call. = FALSE)
      }
    }
    receive_requests <- stats::setNames(lapply(
      seq_along(selected_names), function(target_index) {
        source_index <- setdiff(seq_along(selected_names), target_index)
        as.call(c(
          list(as.name("exactGCVecmulValidityReceiveDS")),
          list(peer_blob = validity[[selected_names[[source_index]]]]$peer_blob),
          common))
      }), selected_names)
    received <- .dsvert_aggregate_strict(
      conns, receive_requests, operation = "exact MPC validity receipt",
      .aggregate = .aggregate)
    for (target_index in seq_along(selected_names)) {
      value <- received[[selected_names[[target_index]]]]
      if (!is.list(value) || !isTRUE(value$stored) ||
          !identical(value$state, "checked")) {
        stop("Exact MPC multiplication failed.", call. = FALSE)
      }
    }
    committed <- .dsvert_aggregate_strict(
      conns, as.call(c(list(as.name("exactGCVecmulCommitDS")), common)),
      operation = "exact MPC checked multiplication commit",
      .aggregate = .aggregate)
    expected <- if (index == chunk_count) "committed" else "partial"
    for (server in selected_names) {
      value <- committed[[server]]
      if (!is.list(value) || !isTRUE(value$stored) ||
          !identical(value$state, expected)) {
        stop("Exact MPC multiplication failed.", call. = FALSE)
      }
    }
  }
  batch_committed <- TRUE
  invisible(list(
    capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
    operation_id = batch, destination_key = output_key,
    policy_id = policy_id, plan_id = plan$plan_id,
    ring_bits = plan$ring_bits, frac_bits = plan$frac_bits,
    backend = plan$backend))
}

#' Run one exact server-to-server operation over the DSI connections
#'
#' Internal orchestration primitive.  It performs one fan-out exchange per
#' cycle and retries the exact same request after transient DSI failures.  There
#' is deliberately no request counter or rate limit; the optional monotonic
#' timeout is an inactivity lease renewed only by verified bytes, offsets or a
#' terminal state transition.
#' @keywords internal
.dsvert_exact_gc_run <- function(
    datasources, server_names = names(datasources), servers = seq_along(datasources),
    session_id, operation_id, source_key, output_key, operation,
    ring, frac_bits = 0L, vector_len, purpose,
    transport_ready = FALSE,
    timeout_seconds = getOption("dsvert.exact_gc.timeout_seconds", 900),
    initialized = NULL, mul_contract = NULL,
    alignment_source_count = NULL,
    analysis_contract = NULL,
    .aggregate = DSI::datashield.aggregate) {
  analysis <- if (is.null(analysis_contract)) NULL else
    .dsvert_exact_gc_analysis_binding(analysis_contract)
  if (length(servers) != 2L || anyDuplicated(servers) ||
      any(servers < 1L) || any(servers > length(datasources))) {
    stop("Exact MPC requires exactly two distinct peer servers.", call. = FALSE)
  }
  conns <- datasources[servers]
  names_selected <- server_names[servers]
  if (length(names_selected) != 2L || anyNA(names_selected) ||
      any(!nzchar(names_selected)) || anyDuplicated(names_selected)) {
    stop("Exact MPC servers require two unique federation names.",
         call. = FALSE)
  }
  names(conns) <- names_selected
  if (!is.character(session_id) || length(session_id) != 1L ||
      is.na(session_id) ||
      !grepl("^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$",
             session_id) ||
      !is.character(operation_id) || length(operation_id) != 1L ||
      is.na(operation_id) ||
      !grepl("^op_[0-9a-f]{32}$", operation_id)) {
    stop("Exact MPC requires canonical session and operation identifiers.",
         call. = FALSE)
  }
  ring <- suppressWarnings(as.numeric(ring))
  frac_bits <- suppressWarnings(as.numeric(frac_bits))
  vector_len <- suppressWarnings(as.numeric(vector_len))
  allowed_operations <- c(
      "compare-signed", "truncate-floor", "mul-truncate-checked",
      "count-guard", "clamp-count", "joint-dp-vector-laplace-v3",
      "alignment-mask-ring128", "formal-glm-phase19-schedule-v1")
  if (!is.null(analysis)) {
    allowed_operations <- c(allowed_operations, "joint-dp-laplace-v2")
  }
  if (!is.character(operation) || length(operation) != 1L ||
      is.na(operation) || !operation %in% allowed_operations ||
      !is.character(purpose) || length(purpose) != 1L || is.na(purpose) ||
      nchar(purpose, type = "bytes") > 128L ||
      !grepl("^[a-z][a-z0-9_.:/-]*$", purpose) ||
      length(ring) != 1L || is.na(ring) || !is.finite(ring) ||
      ring != floor(ring) ||
      ring < .DSVERT_CLIENT_EXACT_GC_MIN_RING_BITS ||
      ring > .DSVERT_CLIENT_EXACT_GC_MAX_RING_BITS ||
      length(frac_bits) != 1L || is.na(frac_bits) || !is.finite(frac_bits) ||
      frac_bits != floor(frac_bits) ||
      frac_bits < 0 || frac_bits > ring - 1L ||
      (operation %in% c(
         "compare-signed", "count-guard", "clamp-count",
         "joint-dp-laplace-v2",
         "joint-dp-vector-laplace-v3", "alignment-mask-ring128",
         "formal-glm-phase19-schedule-v1") &&
       frac_bits != 0L) ||
      length(vector_len) != 1L || is.na(vector_len) ||
      !is.finite(vector_len) ||
      vector_len != floor(vector_len) || vector_len < 1 || vector_len > 4096) {
    stop("Invalid exact MPC operation contract.", call. = FALSE)
  }
  ring <- as.integer(ring)
  frac_bits <- as.integer(frac_bits)
  vector_len <- as.integer(vector_len)
  if (!is.null(analysis) &&
      (!identical(operation, "joint-dp-laplace-v2") || ring != 127L ||
       frac_bits != 0L || vector_len != 1L)) {
    stop("Invalid analysis-bound exact MPC Count shape.", call. = FALSE)
  }
  if (identical(operation, "joint-dp-vector-laplace-v3") &&
      (ring != 128L || vector_len >
         .DSVERT_CLIENT_JOINT_DP_VECTOR_MAX_CHUNK)) {
    stop("Invalid exact MPC biomedical vector shape.", call. = FALSE)
  }
  if (identical(operation, "formal-glm-phase19-schedule-v1") &&
      (ring != 128L || frac_bits != 0L || vector_len > 4L)) {
    stop("Invalid internal formal-GLM schedule shape.", call. = FALSE)
  }
  alignment_k <- if (identical(operation, "alignment-mask-ring128")) {
    if (ring != 128L) {
      stop("Invalid exact MPC alignment-mask ring.", call. = FALSE)
    }
    as.integer(.dsvert_exact_gc_number(
      alignment_source_count, "alignment source count", 2, 64))
  } else {
    if (!is.null(alignment_source_count)) {
      stop("Unexpected exact MPC alignment source count.", call. = FALSE)
    }
    NULL
  }
  container_bits <- .dsvert_exact_gc_container_bits(ring)
  if (identical(operation, "mul-truncate-checked")) {
    plan_max_chunk <- if (is.list(mul_contract)) {
      suppressWarnings(as.numeric(mul_contract$max_chunk))
    } else {
      numeric()
    }
    expected_max_chunk <- if (is.list(mul_contract) &&
                              identical(mul_contract$backend, "ring127-ot")) {
      if (ring != 127L || frac_bits != 50L) NA_integer_ else 256L
    } else if (is.list(mul_contract) &&
               identical(mul_contract$backend, "direct-wide")) {
      .dsvert_exact_gc_direct_mul_max_chunk(ring)
    } else {
      NA_integer_
    }
    if (length(plan_max_chunk) != 1L || is.na(plan_max_chunk) ||
        !is.finite(plan_max_chunk) || plan_max_chunk != floor(plan_max_chunk) ||
        !identical(as.integer(plan_max_chunk), expected_max_chunk) ||
        vector_len > plan_max_chunk) {
      stop("Invalid exact MPC multiplication shape.", call. = FALSE)
    }
  }
  input_containers <- if (identical(operation, "mul-truncate-checked")) {
    7L * vector_len + 1L
  } else if (identical(operation, "alignment-mask-ring128")) {
    3L * vector_len + 4L * alignment_k + 1L
  } else if (identical(operation, "count-guard")) {
    2L * vector_len + 1L
  } else if (identical(operation, "truncate-floor")) {
    3L * vector_len
  } else {
    3L * vector_len
  }
  if (container_bits * input_containers >
      .DSVERT_CLIENT_EXACT_GC_MAX_CIRCUIT_TYPE_BITS) {
    stop(.dsvert_client_resource_oversize(
      requested_bytes = container_bits * input_containers,
      capacity_bytes = .DSVERT_CLIENT_EXACT_GC_MAX_CIRCUIT_TYPE_BITS,
      scope = "exact MPC circuit shape"))
  }
  timeout_seconds <- suppressWarnings(as.numeric(timeout_seconds))
  if (length(timeout_seconds) != 1L || is.na(timeout_seconds) ||
      timeout_seconds <= 0) {
    stop("Exact MPC timeout must be positive.", call. = FALSE)
  }
  if (!isTRUE(transport_ready)) {
    if (is.null(analysis)) {
      .dsvert_setup_exact_gc_transport(
        datasources, server_names, servers, session_id,
        .aggregate = .aggregate)
    } else {
      .dsvert_setup_exact_gc_transport(
        datasources, server_names, servers, session_id,
        analysis_contract = analysis$contract, .aggregate = .aggregate)
    }
  }
  committed <- FALSE
  on.exit(if (!committed) {
    .dsvert_exact_gc_abort(conns, session_id, operation_id, .aggregate)
  }, add = TRUE)

  if (is.null(initialized)) {
    stop("Exact MPC requires a purpose-specific server-minted start.",
         call. = FALSE)
  }
  init <- initialized
  contract <- .dsvert_exact_gc_validate_init(
    init, names_selected, operation, ring, frac_bits, vector_len, purpose,
    mul_contract = mul_contract,
    alignment_source_count = alignment_k,
    analysis_binding = analysis)
  timeout_seconds <- min(timeout_seconds, contract$ttl_seconds)
  peer_ids <- contract$peer_ids
  names(peer_ids) <- names_selected
  server_for_peer <- stats::setNames(names_selected, unname(peer_ids))
  read_offset <- stats::setNames(rep(0, 2L), unname(peer_ids))
  pending <- stats::setNames(vector("list", 2L), unname(peer_ids))
  inbound_seen <- stats::setNames(rep(0, 2L), unname(peer_ids))
  heartbeat_seen <- stats::setNames(
    as.numeric(contract$worker_heartbeats[names_selected]),
    unname(peer_ids))
  complete <- stats::setNames(rep(FALSE, 2L), unname(peer_ids))
  preferred_peer <- NULL
  now <- .dsvert_exact_gc_monotonic_seconds()
  if (!is.numeric(now) || length(now) != 1L || is.na(now) ||
      !is.finite(now)) {
    stop("Exact MPC monotonic clock is unavailable.", call. = FALSE)
  }
  deadline <- if (is.finite(timeout_seconds)) now + timeout_seconds else Inf
  runtime_deadline <- now + contract$max_runtime_seconds

  repeat {
    now <- .dsvert_exact_gc_monotonic_seconds()
    if (!is.numeric(now) || length(now) != 1L || is.na(now) ||
        !is.finite(now)) {
      stop("Exact MPC monotonic clock is unavailable.", call. = FALSE)
    }
    if (now > runtime_deadline) {
      stop("Exact MPC total runtime lease expired before both peers committed output.",
           call. = FALSE)
    }
    if (now > deadline) {
      stop("Exact MPC transport inactivity lease expired before both peers committed output.",
           call. = FALSE)
    }
    requests <- list()
    pending_sources <- unname(peer_ids)[!vapply(
      pending[unname(peer_ids)], is.null, logical(1L))]
    had_pending <- length(pending_sources) > 0L
    # A recipient's reply is the only event that can acknowledge the one
    # pending source envelope. Polling that source in the same cycle can only
    # return its cached duplicate because its read offset has not advanced.
    # Avoid the redundant authenticated DSI call, but retain the full fan-out
    # when both directions have an outstanding envelope.
    polled_peers <- if (length(pending_sources) == 1L) {
      setdiff(unname(peer_ids), pending_sources)
    } else {
      unname(peer_ids)
    }
    preferred_before <- preferred_peer
    for (target in polled_peers) {
      source <- setdiff(unname(peer_ids), target)
      server <- server_for_peer[[target]]
      fields <- .dsvert_exact_gc_delivery_fields(pending[[source]])
      arguments <- c(list(
        session_id = session_id,
        operation_id = operation_id,
        peer_id = target,
        read_offset = read_offset[[target]]), fields, list(
        long_poll = !had_pending &&
          (is.null(preferred_before) || identical(target, preferred_before))))
      requests[[server]] <- as.call(c(
        list(as.name("exactGCExchangeDS")), arguments))
    }
    cycle <- .dsvert_fanout_cycle(
      conns[unname(server_for_peer[polled_peers])], requests,
      operation = "exact MPC byte exchange",
      .aggregate = .aggregate)
    if (identical(cycle$state, "unavailable")) {
      # State and offsets remain unchanged: the next iteration is an exact
      # retry.  No request count is consumed or consulted.
      .dsvert_exact_gc_sleep(0.01)
      next
    }
    response <- cycle$result
    validated <- list()
    for (peer in polled_peers) {
      server <- server_for_peer[[peer]]
      validated[[peer]] <- .dsvert_exact_gc_validate_exchange(
        response[[server]], peer)
    }
    progressed <- FALSE
    new_sources <- character()

    # A fixed-cadence, server-verified worker heartbeat is valid progress. An
    # unchanged counter (including any number of successful empty relay polls)
    # never renews this client-side inactivity deadline.
    for (peer in polled_peers) {
      heartbeat <- validated[[peer]]$worker_heartbeat
      if (heartbeat < heartbeat_seen[[peer]]) {
        stop("An exact MPC peer rolled back its worker heartbeat.",
             call. = FALSE)
      }
      if (heartbeat > heartbeat_seen[[peer]]) {
        heartbeat_seen[[peer]] <- heartbeat
        progressed <- TRUE
      }
    }

    # A target acknowledges the one pending envelope emitted by its peer.
    for (target in polled_peers) {
      source <- setdiff(unname(peer_ids), target)
      acknowledgment <- validated[[target]]$inbound_size
      if (acknowledgment < inbound_seen[[target]]) {
        stop("An exact MPC peer rolled back an acknowledged byte offset.",
             call. = FALSE)
      }
      inbound_seen[[target]] <- acknowledgment
      delivery <- pending[[source]]
      if (!is.null(delivery)) {
        expected <- delivery$offset + delivery$chunk_bytes
        if (acknowledgment != expected) {
          stop("An exact MPC peer returned a conflicting byte acknowledgment.",
               call. = FALSE)
        }
        read_offset[[source]] <- expected
        pending[[source]] <- NULL
        progressed <- TRUE
      }
    }

    # Capture at most one new source envelope.  If an acknowledgment advanced
    # the offset in this same fan-out cycle, its cached duplicate is ignored.
    for (source in polled_peers) {
      target <- setdiff(unname(peer_ids), source)
      envelope <- validated[[source]]$outbound
      if (is.null(envelope)) next
      envelope <- .dsvert_exact_gc_validate_envelope(
        envelope, source, target, session_id, operation_id,
        contract$context_hash, contract$chunk_bytes)
      if (envelope$offset < read_offset[[source]]) next
      if (envelope$offset != read_offset[[source]]) {
        stop("An exact MPC peer returned an outbound byte gap.", call. = FALSE)
      }
      existing <- pending[[source]]
      if (!is.null(existing)) {
        if (!identical(existing, envelope)) {
          stop("An exact MPC peer changed an unacknowledged envelope.",
               call. = FALSE)
        }
      } else {
        pending[[source]] <- envelope
        new_sources <- c(new_sources, source)
        progressed <- TRUE
      }
    }
    if (length(new_sources) == 1L) {
      preferred_peer <- new_sources[[1L]]
    } else if (length(new_sources) > 1L) {
      preferred_peer <- NULL
    } else if (!had_pending && !is.null(preferred_before)) {
      # The predicted sender reached a public protocol boundary. The next
      # fan-out probes both peers so a direction switch costs at most one fixed
      # long-poll window and cannot stall behind a stale prediction.
      preferred_peer <- NULL
    }
    complete_before <- complete
    for (peer in polled_peers) {
      complete[[peer]] <- identical(validated[[peer]]$state, "complete") &&
        isTRUE(validated[[peer]]$stored)
    }
    if (any(complete & !complete_before)) progressed <- TRUE
    if (isTRUE(progressed) && is.finite(timeout_seconds)) {
      now <- .dsvert_exact_gc_monotonic_seconds()
      if (!is.numeric(now) || length(now) != 1L || is.na(now) ||
          !is.finite(now)) {
        stop("Exact MPC monotonic clock is unavailable.", call. = FALSE)
      }
      deadline <- now + timeout_seconds
    }
    if (all(complete) && all(vapply(pending, is.null, logical(1L)))) break
    if (!progressed) .dsvert_exact_gc_sleep(0.005)
  }

  # The terminating exchange response from each pinned peer has already been
  # validated as complete+stored under this operation's context.  A separate
  # remotely callable status probe would add surface without strengthening the
  # commit condition.
  committed <- TRUE
  result <- list(
    capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
    operation_id = operation_id, context_hash = contract$context_hash,
    output_key = output_key, operation = operation, ring_bits = ring,
    frac_bits = frac_bits, vector_len = vector_len,
    threshold = contract$threshold)
  if (!is.null(analysis)) {
    result$analysis_binding_sha256 <- analysis$sha256
  }
  invisible(result)
}

.dsvert_exact_gc_capability_contract <- function() {
  list(
    capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
    available = TRUE, allowed = TRUE, e2e_verified = TRUE,
    canonical_encoding = TRUE, canonical_input_encoding = TRUE,
    shape_bounds_enforced = TRUE,
    fail_closed_overflow = FALSE, runtime_bounds_enforced = FALSE,
    raw_product_overflow_guard = FALSE, checked_mul_truncate = FALSE,
    dynamic_ring_fallback = TRUE,
    vecmul_numeric_precondition =
      "strict producer-minted input manifest required for promotion",
    exact_truncation = TRUE, exact_comparison = FALSE,
    core_exact_comparison = TRUE, comparison_e2e_verified = FALSE,
    vecmul_truncation_e2e_verified = FALSE,
    count_guard_e2e_verified = TRUE,
    clamp_count_e2e_verified = TRUE,
    joint_dp_count_e2e_verified = TRUE,
    joint_dp_vector_e2e_verified = TRUE,
    alignment_mask_e2e_verified = TRUE,
    multiprecision_truncation_e2e_verified = TRUE,
    workload_glm_e2e_verified = FALSE,
    operations = c(
      "truncate-floor", "count-guard", "clamp-count",
      "joint-dp-laplace-v2", "joint-dp-vector-laplace-v3",
      "alignment-mask-ring128"),
    core_operations = c(
      "compare-signed", "truncate-floor", "mul-truncate-checked",
      "count-guard", "clamp-count", "joint-dp-laplace-v2",
      "joint-dp-vector-laplace-v3", "alignment-mask-ring128"),
    verified_purposes = c(
      "count-guard", "multiprecision-truncate", "joint-dp-count-clamp",
      "joint-dp-count-one-draw",
      "joint-dp-biomedical-vector-one-draw",
      "private-alignment-mask-ring128"),
    truncation_semantics = "floor",
    supported_ring_bits = .DSVERT_CLIENT_EXACT_GC_MIN_RING_BITS:
      .DSVERT_CLIENT_EXACT_GC_MAX_RING_BITS,
    wire_container_bits = .DSVERT_CLIENT_EXACT_GC_WIRE_CONTAINER_BITS,
    max_ring_bits = .DSVERT_CLIENT_EXACT_GC_MAX_RING_BITS,
    max_frac_bits = .DSVERT_CLIENT_EXACT_GC_MAX_RING_BITS - 1L)
}
