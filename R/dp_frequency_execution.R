.DSVERT_CLIENT_DP_FREQUENCY_SOURCE_VERSION <- "dsvert-dp-frequency-source-window-v1"
.DSVERT_CLIENT_DP_FREQUENCY_FINAL_VERSION <- "dsvert-dp-frequency-finalization-v1"
.DSVERT_CLIENT_DP_FREQUENCY_RELEASE_VERSION <- "dsvert-dp-frequency-release-v1"
.DSVERT_CLIENT_DP_FREQUENCY_RESULT_VERSION <- "dsvert-dp-frequency-execution-result-v1"
.DSVERT_CLIENT_DP_FREQUENCY_SOURCE_CAPABILITY <- "blob.analysis-dp.frequency-source-to-finalizer.v1"
.DSVERT_CLIENT_DP_FREQUENCY_CONTEXT_VERSION <- "dsvert-dp-frequency-typed-context-v1"
.DSVERT_CLIENT_DP_FREQUENCY_CONTEXT_PURPOSE <- "analysis-dp.frequency-source-to-finalizer.v1"
.DSVERT_CLIENT_DP_FREQUENCY_AUTH_SET_DOMAIN <- "dsVert/dp-frequency/authorization-set/v1|"
.DSVERT_CLIENT_DP_FREQUENCY_CHUNK_DOMAIN <- "dsVert/dp-frequency/final-binary-chunk/v1|"
.DSVERT_CLIENT_DP_FREQUENCY_WINDOW_DOMAIN <- "dsVert/dp-frequency/final-window/v1|"
.DSVERT_CLIENT_DP_FREQUENCY_RELEASE_DOMAIN <- "dsVert/dp-frequency/release/v1|"
.DSVERT_CLIENT_DP_FREQUENCY_RELEASE_SIGNATURE_DOMAIN <- "dsVert/dp-frequency/release-signature/v1|"
.DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES <- 8192L
.DSVERT_CLIENT_DP_FREQUENCY_WINDOW_COORDINATES <- 65536L
.DSVERT_CLIENT_DP_FREQUENCY_CIPHERTEXT_CHARS <- 1409104L
.DSVERT_CLIENT_DP_FREQUENCY_AUTHORIZATIONS_MAX_BYTES <- 1024L^2
.dsvert_dp_frequency_execution_integer_v1 <- function(value, what, minimum = 0, maximum = 1e6) {
  if (length(value) != 1L || is.list(value) || is.object(value) || is.na(value) ||
      (!is.numeric(value) && !is.character(value)) ||
      (is.character(value) && (nchar(value, type = "bytes") > 7L ||
        !grepl("^(0|[1-9][0-9]*)$", value)))) {
    stop("Invalid Frequency execution ", what, ".", call. = FALSE)
  }
  value <- suppressWarnings(as.numeric(value))
  if (!is.finite(value) || value != floor(value) ||
      value < minimum || value > maximum) {
    stop("Invalid Frequency execution ", what, ".", call. = FALSE)
  }
  as.integer(value)
}
.dsvert_dp_frequency_execution_integer_matches_v1 <- function(value, expected, what) {
  parsed <- tryCatch(.dsvert_dp_frequency_execution_integer_v1(
    value, what), error = function(error) NULL)
  !is.null(parsed) && identical(parsed, as.integer(expected))
}
.dsvert_dp_frequency_execution_geometry_v1 <- function(worker) {
  if (!is.list(worker)) {
    stop("Invalid prepared Frequency execution geometry.", call. = FALSE)
  }
  d <- .dsvert_dp_frequency_execution_integer_v1(
    worker$d, "dimension", 1L, 1e6)
  chunks <- as.integer(ceiling(
    d / .DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES))
  valid <- identical(as.numeric(worker$ring_bits), 128) &&
    identical(as.numeric(worker$frac_bits), 0) &&
    identical(as.numeric(worker$chunk_coordinates),
              as.numeric(min(.DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES,
                             d))) &&
    identical(as.numeric(worker$chunk_count), as.numeric(chunks))
  if (!isTRUE(valid)) {
    stop("Invalid prepared Frequency execution geometry.", call. = FALSE)
  }
  list(d = d, chunk_count = chunks, window_count = as.integer(ceiling(
    d / .DSVERT_CLIENT_DP_FREQUENCY_WINDOW_COORDINATES)))
}
.dsvert_dp_frequency_execution_window_v1 <- function(geometry, index) {
  index <- .dsvert_dp_frequency_execution_integer_v1(
    index, "window index", 0L, geometry$window_count - 1L)
  offset <- index * .DSVERT_CLIENT_DP_FREQUENCY_WINDOW_COORDINATES
  count <- min(.DSVERT_CLIENT_DP_FREQUENCY_WINDOW_COORDINATES,
               geometry$d - offset)
  list(
    index = index, count = as.integer(count), offset = as.integer(offset),
    first_chunk = as.integer(8L * index), chunks = as.integer(ceiling(
      count / .DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES)))
}
.dsvert_dp_frequency_execution_compiled_v1 <- function(prepared, peers) {
  fields <- c(
    "version", "session_id", "source_owner", "claim", "config",
    "receipts", "contract", "worker_static", "authorities",
    "authorizations", "transport")
  if (!.dsvert_dp_frequency_client_object_v1(prepared, fields) ||
      !identical(prepared$version,
                 .DSVERT_CLIENT_DP_FREQUENCY_PREPARED_VERSION) ||
      !is.list(prepared$receipts) || is.null(names(prepared$receipts)) ||
      anyNA(names(prepared$receipts)) || anyDuplicated(names(prepared$receipts)) ||
      !setequal(names(prepared$receipts), peers)) {
    stop("Invalid prepared Frequency execution context.", call. = FALSE)
  }
  envelopes <- stats::setNames(lapply(peers, function(peer) list(
    config = prepared$config, receipt = prepared$receipts[[peer]])), peers)
  compiled <- .dsvert_dp_frequency_client_compile_v1(
    envelopes, peers, prepared$claim)
  same <- identical(prepared$source_owner,
                    compiled$config$source_owner$peer_name) &&
    identical(prepared$claim, compiled$claim) &&
    identical(prepared$config, compiled$config) &&
    identical(prepared$receipts, compiled$receipts) &&
    identical(prepared$contract, compiled$contract) &&
    identical(prepared$worker_static, compiled$worker_static) &&
    identical(prepared$authorities, compiled$authorities)
  if (!isTRUE(same)) stop(
    "Prepared Frequency proof disagrees with its K-wide compilation.",
    call. = FALSE)
  compiled
}
.dsvert_dp_frequency_execution_cleanup_capability_v1 <- function(value, peer, session_id, digest, pins) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || nchar(value, type = "bytes") > 16384L) {
    stop("Invalid prepared Frequency cleanup capability.", call. = FALSE)
  }
  parsed <- tryCatch(jsonlite::fromJSON(
    value, simplifyVector = FALSE), error = function(error) NULL)
  contract <- if (is.list(parsed)) parsed$contract else NULL
  valid <- .dsvert_dp_frequency_client_object_v1(
      parsed, c("version", "contract", "signature")) &&
    identical(parsed$version,
              .DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION) &&
    .dsvert_dp_frequency_client_object_v1(contract, c(
      "version", "session_id", "cleanup_purpose", "operation_scope",
      "peer_binding_digest")) &&
    identical(contract$version,
              .DSVERT_CLIENT_EXACT_GC_CLEANUP_CAPABILITY_VERSION) &&
    identical(contract$session_id, session_id) &&
    identical(contract$cleanup_purpose,
              .DSVERT_CLIENT_EXACT_GC_FREQUENCY_CLEANUP_PURPOSE) &&
    identical(contract$operation_scope,
              "all_and_only_operations_in_bound_exact_session_v1") &&
    identical(contract$peer_binding_digest, digest) &&
    identical(.dsvert_dp_frequency_client_wire_json_v1(parsed), value)
  if (!isTRUE(valid)) {
    stop("Invalid prepared Frequency cleanup capability.", call. = FALSE)
  }
  .dsvert_dp_frequency_client_verify_v1(
    charToRaw(paste0("dsVert/exact-gc/cleanup-capability/v1|",
                     .dsvert_dp_frequency_client_wire_json_v1(contract))),
    unname(pins[[peer]]), parsed$signature, "cleanup capability")
  invisible(TRUE)
}
.dsvert_dp_frequency_execution_context_v1 <- function(prepared, datasources) {
  datasources <- .dsvert_dp_datasources(datasources)
  peers <- names(datasources)
  compiled <- .dsvert_dp_frequency_execution_compiled_v1(prepared, peers)
  session_id <- prepared$session_id
  if (!is.character(session_id) || length(session_id) != 1L ||
      is.na(session_id) || !grepl(paste0(
        "^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-",
        "[89ab][0-9a-f]{3}-[0-9a-f]{12}$"), session_id)) {
    stop("Invalid prepared Frequency session identifier.", call. = FALSE)
  }
  authorizations <- .dsvert_dp_frequency_client_authorizations_v1(
    prepared$authorizations, session_id, compiled)
  authorities <- compiled$authorities
  transport <- prepared$transport
  capabilities <- attr(transport, "exact_gc_cleanup_capabilities", exact = TRUE)
  purpose <- attr(transport, "exact_gc_cleanup_purpose", exact = TRUE)
  peer_digest <- attr(transport, "exact_gc_peer_binding_digest", exact = TRUE)
  frequency <- .dsvert_exact_gc_frequency_binding(compiled)
  valid <- is.list(transport) && !is.null(names(transport)) &&
    !anyNA(names(transport)) && !anyDuplicated(names(transport)) &&
    setequal(names(transport), authorities) &&
    is.list(capabilities) && !is.null(names(capabilities)) &&
    !anyNA(names(capabilities)) && !anyDuplicated(names(capabilities)) &&
    setequal(names(capabilities), authorities) &&
    identical(purpose,
              .DSVERT_CLIENT_EXACT_GC_FREQUENCY_CLEANUP_PURPOSE) &&
    .dsvert_dp_analysis_frequency_hex_v1(peer_digest) &&
    identical(attr(transport, "exact_gc_frequency_binding", exact = TRUE),
              frequency$binding) &&
    identical(attr(
      transport, "exact_gc_frequency_binding_sha256", exact = TRUE),
      frequency$sha256)
  if (!isTRUE(valid)) stop(
    "Invalid prepared Frequency transport or cleanup binding.", call. = FALSE)
  for (peer in authorities) {
    raw <- .dsvert_joint_dp_client_b64url(
      transport[[peer]], 32L, "Frequency transport public key")
    if (!identical(.dsvert_exact_gc_b64url_encode(raw), transport[[peer]])) {
      stop("Invalid prepared Frequency transport key.", call. = FALSE)
    }
    .dsvert_dp_frequency_execution_cleanup_capability_v1(
      capabilities[[peer]], peer, session_id, peer_digest,
      compiled$config$peer_pins)
  }
  geometry <- .dsvert_dp_frequency_execution_geometry_v1(
    compiled$worker_static)
  auth_hash <- .dsvert_dp_frequency_client_hash_v1(
    .DSVERT_CLIENT_DP_FREQUENCY_AUTH_SET_DOMAIN,
    unname(authorizations[authorities]))
  list(
    prepared = prepared, compiled = compiled, session_id = session_id,
    authorities = authorities, source = authorities[[1L]],
    secondary = authorities[[2L]], conns = datasources[authorities],
    authorizations = authorizations[authorities], transport = transport,
    geometry = geometry, authorization_set_sha256 = auth_hash)
}
.dsvert_dp_frequency_execution_authorizations_v1 <- function(context) {
  value <- .dsvert_joint_dp_client_json(
    unname(context$authorizations[context$authorities]))
  if (nchar(value, type = "bytes") >
      .DSVERT_CLIENT_DP_FREQUENCY_AUTHORIZATIONS_MAX_BYTES) {
    stop("Frequency authorizations exceed their bounded wire contract.",
         call. = FALSE)
  }
  .dsvert_dsi_text_encode(value, "Frequency execution authorizations")
}
.dsvert_dp_frequency_execution_expected_context_v1 <- function(context, operation_id, window) {
  compiled <- context$compiled
  roles <- stats::setNames(as.list(context$authorities),
                           c("source_owner", "secondary_noise_authority"))
  .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_CONTEXT_VERSION,
    purpose = .DSVERT_CLIENT_DP_FREQUENCY_CONTEXT_PURPOSE,
    artifact_key = compiled$contract$artifact_key,
    contract_sha256 = compiled$contract_sha256,
    analysis_binding_sha256 = compiled$binding$sha256,
    worker_static_sha256 = compiled$worker_static_sha256,
    authorization_set_sha256 = context$authorization_set_sha256,
    release_contract_hash = compiled$worker_static$release_contract_hash,
    operation_id = operation_id, window_index = as.character(window$index),
    window_count = as.character(context$geometry$window_count),
    first_chunk_index = as.character(window$first_chunk),
    chunks_in_window = as.character(window$chunks),
    coordinate_offset = as.character(window$offset),
    coordinate_count = as.character(window$count),
    padded_coordinate_count = "65536", ring_bits = "128", frac_bits = "0",
    roles = roles, sender = context$source, recipient = context$secondary))
}
.dsvert_dp_frequency_execution_source_v1 <- function(value, context, expected_context, window, issued = TRUE) {
  common <- c(
    "version", "state", "artifact_key", "window_index", "window_count",
    "intermediate_values_exposed")
  fields <- if (issued) c(common, "context", "ciphertext_chars", "transfer")
  else common
  valid <- .dsvert_dp_frequency_client_object_v1(value, fields) &&
    identical(value$version,
              .DSVERT_CLIENT_DP_FREQUENCY_SOURCE_VERSION) &&
    identical(value$state, if (issued) "issued" else "delivered") &&
    identical(value$artifact_key, context$compiled$contract$artifact_key) &&
    .dsvert_dp_frequency_execution_integer_matches_v1(
      value$window_index, window$index, "source window index") &&
    .dsvert_dp_frequency_execution_integer_matches_v1(
      value$window_count, context$geometry$window_count,
      "source window count") &&
    identical(value$intermediate_values_exposed, FALSE)
  if (!isTRUE(valid)) {
    stop("Frequency source returned an invalid closed window state.",
         call. = FALSE)
  }
  if (!issued) return(.dsvert_dp_analysis_client_canonical_value_v1(value))
  ciphertext <- value$ciphertext_chars
  context_fields <- names(expected_context)
  context_valid <- .dsvert_dp_frequency_client_object_v1(
    value$context, context_fields) &&
    .dsvert_dp_frequency_client_object_v1(
      value$context$roles, names(expected_context$roles)) &&
    all(vapply(value$context[setdiff(context_fields, "roles")], function(x) {
      is.character(x) && length(x) == 1L && !is.na(x) &&
        nchar(x, type = "bytes") <= 1024L
    }, logical(1L))) && all(vapply(value$context$roles, function(x) {
      is.character(x) && length(x) == 1L && !is.na(x) &&
        nchar(x, type = "bytes") <= 128L
    }, logical(1L)))
  valid <- context_valid && identical(value$context, expected_context) &&
    is.character(ciphertext) &&
    length(ciphertext) == 1L && !is.na(ciphertext) &&
    nchar(ciphertext, type = "bytes") ==
      .DSVERT_CLIENT_DP_FREQUENCY_CIPHERTEXT_CHARS &&
    grepl("^[A-Za-z0-9_-]+$", ciphertext, perl = TRUE, useBytes = TRUE)
  transfer <- if (isTRUE(valid)) tryCatch(
    .dsvert_validate_typed_blob_transfer(
      value$transfer, ciphertext, context$secondary),
    error = function(error) NULL) else NULL
  valid <- valid && is.list(transfer) &&
    identical(transfer$capability_id,
              .DSVERT_CLIENT_DP_FREQUENCY_SOURCE_CAPABILITY) &&
    identical(transfer$sender_name, context$source) &&
    identical(transfer$recipient_name, context$secondary) &&
    identical(as.numeric(transfer$payload_chars),
              as.numeric(.DSVERT_CLIENT_DP_FREQUENCY_CIPHERTEXT_CHARS))
  if (!isTRUE(valid)) {
    stop("Frequency source returned an invalid ciphertext transfer.",
         call. = FALSE)
  }
  value$transfer <- transfer
  value
}
.dsvert_dp_frequency_execution_final_v1 <- function(value, context, window, final) {
  fields <- if (final) c(
    "version", "state", "artifact_key", "release",
    "intermediate_values_exposed") else c(
      "version", "state", "artifact_key", "window_index", "window_count",
      "intermediate_values_exposed")
  valid <- .dsvert_dp_frequency_client_object_v1(value, fields) &&
    identical(value$version, .DSVERT_CLIENT_DP_FREQUENCY_FINAL_VERSION) &&
    identical(value$state,
              if (final) "release_committed" else "window_committed") &&
    identical(value$artifact_key, context$compiled$contract$artifact_key) &&
    identical(value$intermediate_values_exposed, FALSE)
  if (!final) valid <- valid &&
    .dsvert_dp_frequency_execution_integer_matches_v1(
      value$window_index, window$index, "final window index") &&
    .dsvert_dp_frequency_execution_integer_matches_v1(
      value$window_count, context$geometry$window_count,
      "final window count")
  if (!isTRUE(valid)) {
    stop("Frequency finalizer returned an invalid closed state.",
         call. = FALSE)
  }
  value
}
.dsvert_dp_frequency_execution_release_v1 <- function(value, context) {
  fields <- c(
    "version", "artifact_key", "contract_sha256",
    "analysis_binding_sha256", "worker_static_sha256",
    "authorization_set_sha256", "release_contract_hash", "primitive",
    "mechanism", "sampler", "d", "chunk_coordinates", "chunk_count",
    "window_count", "coordinate_order_sha256", "bounds", "authority_roles",
    "final_chunk_hashes", "final_window_hashes", "final_vector_root",
    "postprocessing", "intermediate_values_exposed", "public_openings",
    "release_sha256", "signature")
  if (!.dsvert_dp_frequency_client_object_v1(value, fields)) {
    stop("Frequency finalizer returned an invalid closed release.",
         call. = FALSE)
  }
  worker <- context$compiled$worker_static
  contract <- context$compiled$contract
  profile <- worker$selected_profile
  hashes_valid <- is.list(value$final_chunk_hashes) &&
    is.null(names(value$final_chunk_hashes)) &&
    length(value$final_chunk_hashes) == context$geometry$chunk_count &&
    is.list(value$final_window_hashes) &&
    is.null(names(value$final_window_hashes)) &&
    length(value$final_window_hashes) == context$geometry$window_count
  hash_values <- if (isTRUE(hashes_valid)) c(
    value$final_chunk_hashes, value$final_window_hashes,
    list(value$final_vector_root)) else list()
  hashes_valid <- hashes_valid && all(vapply(hash_values, function(hash) {
    is.character(hash) && length(hash) == 1L && !is.na(hash) &&
      nchar(hash, type = "bytes") == 64L &&
      .dsvert_dp_analysis_frequency_hex_v1(hash)
  }, logical(1L)))
  expected_static <- list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_RELEASE_VERSION,
    artifact_key = contract$artifact_key,
    contract_sha256 = context$compiled$contract_sha256,
    analysis_binding_sha256 = context$compiled$binding$sha256,
    worker_static_sha256 = context$compiled$worker_static_sha256,
    authorization_set_sha256 = context$authorization_set_sha256,
    release_contract_hash = worker$release_contract_hash,
    primitive = worker$selected_primitive,
    mechanism = profile$mechanism, sampler = profile$sampler,
    d = worker$d, chunk_coordinates = 8192L,
    chunk_count = worker$chunk_count,
    window_count = context$geometry$window_count,
    coordinate_order_sha256 = contract$semantic$analysis$
      effective_arguments$sampler_plan$coordinate_order_sha256,
    bounds = worker$raw_bound, authority_roles = worker$authority_roles,
    postprocessing = profile$output_transform,
    intermediate_values_exposed = FALSE, public_openings = 1L)
  static <- value[names(expected_static)]
  structured <- c("bounds", "authority_roles")
  scalar <- function(x, bytes = 1024L) length(x) == 1L && !is.list(x) &&
    !anyNA(x) && (!is.character(x) || nchar(x, type = "bytes") <= bytes)
  static_valid <- all(vapply(
    static[setdiff(names(static), structured)], scalar, logical(1L))) &&
    .dsvert_dp_frequency_client_object_v1(
      static$bounds, names(expected_static$bounds)) &&
    .dsvert_dp_frequency_client_object_v1(
      static$authority_roles, names(expected_static$authority_roles)) &&
    all(vapply(c(static$bounds, static$authority_roles),
               scalar, logical(1L))) && identical(
    .dsvert_joint_dp_client_json(
      .dsvert_dp_analysis_client_canonical_value_v1(static)),
    .dsvert_joint_dp_client_json(
      .dsvert_dp_analysis_client_canonical_value_v1(expected_static)))
  signature_valid <- is.character(value$signature) &&
    length(value$signature) == 1L && !is.na(value$signature) &&
    nchar(value$signature, type = "bytes") == 86L &&
    grepl("^[A-Za-z0-9_-]{86}$", value$signature)
  release_hash_valid <- is.character(value$release_sha256) &&
    length(value$release_sha256) == 1L && !is.na(value$release_sha256) &&
    nchar(value$release_sha256, type = "bytes") == 64L &&
    .dsvert_dp_analysis_frequency_hex_v1(value$release_sha256)
  if (!isTRUE(hashes_valid) || !isTRUE(static_valid) ||
      !isTRUE(signature_valid) || !isTRUE(release_hash_valid)) {
    stop("Frequency finalizer returned a corrupt or misbound release.",
         call. = FALSE)
  }
  core <- .dsvert_dp_analysis_client_canonical_value_v1(
    value[setdiff(names(value), c("release_sha256", "signature"))])
  expected_hash <- .dsvert_dp_frequency_client_hash_v1(
    .DSVERT_CLIENT_DP_FREQUENCY_RELEASE_DOMAIN, core)
  if (!identical(value$release_sha256, expected_hash)) {
    stop("Frequency finalizer returned a corrupt or misbound release.",
         call. = FALSE)
  }
  signed <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    core, list(release_sha256 = expected_hash)))
  secondary_identity <- unname(
    worker$authority_roles$secondary_noise_authority)
  .dsvert_dp_frequency_client_verify_v1(
    charToRaw(paste0(
      .DSVERT_CLIENT_DP_FREQUENCY_RELEASE_SIGNATURE_DOMAIN,
      .dsvert_joint_dp_client_json(signed))),
    secondary_identity, value$signature, "release")
  .dsvert_dp_analysis_client_canonical_value_v1(value)
}
.dsvert_dp_frequency_execution_ring128_v1 <- function(values) {
  values <- as.integer(values)
  output <- raw(16L * length(values))
  if (length(values)) {
    starts <- seq.int(1L, length(output), by = 16L)
    for (byte in 0:3) {
      output[starts + byte] <- as.raw(bitwAnd(
        bitwShiftR(values, 8L * byte), 255L))
    }
  }
  output
}
.dsvert_dp_frequency_execution_chunk_hash_v1 <- function(values, index, offset) {
  raw <- .dsvert_dp_frequency_execution_ring128_v1(values)
  .dsvert_dp_frequency_client_hash_v1(
    .DSVERT_CLIENT_DP_FREQUENCY_CHUNK_DOMAIN, list(
      version = "dsvert-dp-frequency-final-binary-chunk-v1",
      chunk_index = as.integer(index), coordinate_offset = as.integer(offset),
      coordinate_count = length(values), ring128_b64 = gsub(
        "[\r\n]", "", jsonlite::base64_enc(raw))))
}
.dsvert_dp_frequency_execution_merkle_v1 <- function(hashes) {
  hash <- function(domain, values) digest::digest(c(
    charToRaw(domain), as.raw(0L),
    unlist(lapply(values, .dsvert_vector_hex_raw), recursive = FALSE)),
    "sha256", serialize = FALSE)
  nodes <- vapply(hashes, function(value) hash(
    "dsVert/dp-frequency/merkle-leaf/v1", value), character(1L))
  while (length(nodes) > 1L) {
    if (length(nodes) %% 2L) nodes <- c(nodes, tail(nodes, 1L))
    nodes <- vapply(seq.int(1L, length(nodes), by = 2L), function(index) {
      hash("dsVert/dp-frequency/merkle-parent/v1", nodes[index + 0:1])
    }, character(1L))
  }
  unname(nodes[[1L]])
}
.dsvert_dp_frequency_execution_replay_v1 <- function(value, context, window, release) {
  fields <- c(
    "version", "state", "release", "window", "intermediate_values_exposed")
  valid <- .dsvert_dp_frequency_client_object_v1(value, fields) &&
    identical(value$version, .DSVERT_CLIENT_DP_FREQUENCY_FINAL_VERSION) &&
    identical(value$state, "release_replay") &&
    identical(value$intermediate_values_exposed, FALSE)
  replayed_release <- value$release
  valid <- valid && is.list(replayed_release) &&
    length(replayed_release) == length(release) &&
    identical(names(replayed_release), names(release)) &&
    identical(replayed_release$release_sha256, release$release_sha256) &&
    identical(replayed_release$signature, release$signature) &&
    identical(replayed_release$final_vector_root, release$final_vector_root)
  public <- value$window
  valid <- valid && .dsvert_dp_frequency_client_object_v1(public, c(
    "version", "window_index", "coordinate_offset", "coordinate_count",
    "chunks", "window_sha256")) &&
    identical(public$version, "dsvert-dp-frequency-final-window-v1") &&
    .dsvert_dp_frequency_execution_integer_matches_v1(
      public$window_index, window$index, "replay window index") &&
    .dsvert_dp_frequency_execution_integer_matches_v1(
      public$coordinate_offset, window$offset, "replay coordinate offset") &&
    .dsvert_dp_frequency_execution_integer_matches_v1(
      public$coordinate_count, window$count, "replay coordinate count") &&
    is.list(public$chunks) && is.null(names(public$chunks)) &&
    length(public$chunks) == window$chunks
  if (!isTRUE(valid)) {
    stop("Frequency replay returned an invalid closed window.",
         call. = FALSE)
  }
  upper <- .dsvert_dp_frequency_execution_integer_v1(
    context$compiled$worker_static$raw_bound$upper, "release upper bound",
    0L, 1e6)
  hashes <- character(window$chunks)
  values <- numeric(window$count)
  for (position in seq_len(window$chunks)) {
    chunk <- public$chunks[[position]]
    local <- (position - 1L) *
      .DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES
    count <- min(
      .DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES, window$count - local)
    offset <- window$offset + local
    index <- window$first_chunk + position - 1L
    valid <- .dsvert_dp_frequency_client_object_v1(chunk, c(
      "version", "chunk_index", "coordinate_offset", "coordinate_count",
      "values", "chunk_sha256")) &&
      identical(chunk$version, "dsvert-dp-frequency-final-chunk-v1") &&
      .dsvert_dp_frequency_execution_integer_matches_v1(
        chunk$chunk_index, index, "chunk index") &&
      .dsvert_dp_frequency_execution_integer_matches_v1(
        chunk$coordinate_offset, offset, "chunk coordinate offset") &&
      .dsvert_dp_frequency_execution_integer_matches_v1(
        chunk$coordinate_count, count, "chunk coordinate count") &&
      is.list(chunk$values) && is.null(names(chunk$values)) &&
      length(chunk$values) == count
    scalar <- if (isTRUE(valid)) vapply(chunk$values, function(text) {
      is.character(text) && length(text) == 1L && !is.na(text) &&
        nchar(text, type = "bytes") <= 7L &&
        grepl("^(0|[1-9][0-9]*)$", text)
    }, logical(1L)) else FALSE
    valid <- valid && all(scalar)
    texts <- if (isTRUE(valid)) unlist(chunk$values, use.names = FALSE) else NULL
    parsed <- suppressWarnings(as.numeric(texts))
    valid <- valid && length(parsed) == count && all(is.finite(parsed)) &&
      all(parsed <= upper)
    expected <- if (isTRUE(valid))
      .dsvert_dp_frequency_execution_chunk_hash_v1(parsed, index, offset)
    else NULL
    valid <- valid && identical(chunk$chunk_sha256, expected)
    if (!isTRUE(valid)) {
      stop("Frequency replay returned a corrupt final chunk.", call. = FALSE)
    }
    range <- seq.int(local + 1L, local + count)
    values[range] <- parsed
    hashes[[position]] <- expected
  }
  core <- list(
    version = "dsvert-dp-frequency-final-window-v1",
    window_index = window$index, coordinate_offset = window$offset,
    coordinate_count = window$count, chunk_hashes = as.list(hashes))
  expected_window <- .dsvert_dp_frequency_client_hash_v1(
    .DSVERT_CLIENT_DP_FREQUENCY_WINDOW_DOMAIN, core)
  if (!identical(public$window_sha256, expected_window)) {
    stop("Frequency replay returned a corrupt final window.", call. = FALSE)
  }
  list(values = values, chunk_hashes = hashes,
       window_sha256 = expected_window)
}
.dsvert_dp_frequency_execute_v1 <- function(prepared, data_name, datasources,
    .aggregate = DSI::datashield.aggregate,
    .new_context = .dsvert_exact_gc_new_context,
    .store_typed = .dsvert_store_typed_blob,
    .retry = .dsvert_retry_idempotent,
    .execution_context = .dsvert_dp_frequency_execution_context_v1,
    .frequency_cleanup = .dsvert_cleanup_best_effort,
    .exact_cleanup = .dsvert_exact_gc_cleanup_best_effort) {
  dependencies <- list(
    .aggregate, .new_context, .store_typed, .retry, .execution_context,
    .frequency_cleanup, .exact_cleanup)
  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) || !grepl("^[A-Za-z._][A-Za-z0-9._]*$", data_name) ||
      any(!vapply(dependencies, is.function, logical(1L)))) {
    stop("Invalid Frequency execution request or dependency.", call. = FALSE)
  }
  datasources <- .dsvert_dp_datasources(datasources)
  cleanup_session <- if (is.list(prepared) &&
      is.character(prepared$session_id) && length(prepared$session_id) == 1L &&
      !is.na(prepared$session_id)) prepared$session_id else ""
  cleanup_authorities <- tryCatch(
    .dsvert_exact_gc_frequency_analysis_binding(
      prepared$contract)$authority_names, error = function(error) character())
  cleanup_conns <- if (length(cleanup_authorities) == 2L &&
      all(cleanup_authorities %in% names(datasources)))
    datasources[cleanup_authorities] else list()
  cleanup_transport <- if (is.list(prepared)) prepared$transport else NULL
  on.exit({
    tryCatch(.frequency_cleanup(
      cleanup_conns,
      call(name = "dsvertDPFrequencyCleanupDS",
           session_id = cleanup_session),
      .aggregate = .aggregate), error = function(error) NULL)
    tryCatch(.exact_cleanup(
      cleanup_conns, cleanup_session, cleanup_transport,
      .aggregate = .aggregate), error = function(error) NULL)
  }, add = TRUE)
  context <- .execution_context(prepared, datasources)
  cleanup_conns <- context$conns
  cleanup_session <- context$session_id
  cleanup_transport <- context$transport
  generated <- .new_context()
  operation_id <- if (is.list(generated)) generated$operation_id else NULL
  if (!is.character(operation_id) || length(operation_id) != 1L ||
      is.na(operation_id) || !grepl("^op_[0-9a-f]{32}$", operation_id)) {
    stop("Invalid Frequency execution operation identifier.", call. = FALSE)
  }
  authorizations <- .dsvert_dp_frequency_execution_authorizations_v1(context)
  claim <- .dsvert_dsi_text_encode(
    .dsvert_joint_dp_client_json(context$compiled$claim),
    "Frequency source Claim")
  release <- NULL
  for (index in seq_len(context$geometry$window_count) - 1L) {
    window <- .dsvert_dp_frequency_execution_window_v1(
      context$geometry, index)
    expected_context <- .dsvert_dp_frequency_execution_expected_context_v1(
      context, operation_id, window)
    source_call <- call(
      name = "dsvertDPFrequencySourceWindowDS", data_name = data_name,
      session_id = context$session_id, operation_id = operation_id,
      window_index = index,
      public_authorizations_json = authorizations,
      source_claim_json = claim)
    source <- .dsvert_aggregate_strict(
      context$conns[context$source], source_call,
      operation = "Frequency encrypted source window",
      .aggregate = .aggregate)[[context$source]]
    source <- .dsvert_dp_frequency_execution_source_v1(
      source, context, expected_context, window, issued = TRUE)
    .store_typed(
      blob = source$ciphertext_chars, transfer = source$transfer,
      conn = context$conns[context$secondary],
      producer_conn = context$conns[context$source],
      session_id = context$session_id, .aggregate = .aggregate)
    acknowledged <- .retry(
      attempt = function() .dsvert_aggregate_strict(
        context$conns[context$source], source_call,
        operation = "Frequency source delivery acknowledgement",
        .aggregate = .aggregate)[[context$source]],
      classify = function(value) {
        if (is.list(value) && identical(value$state, "issued")) {
          repeated <- .dsvert_dp_frequency_execution_source_v1(
            value, context, expected_context, window, issued = TRUE)
          if (!identical(repeated, source)) {
            stop("Frequency source ciphertext changed during acknowledgement.",
                 call. = FALSE)
          }
          return(list(state = "missing", result = repeated))
        }
        list(
          state = "ack",
          result = .dsvert_dp_frequency_execution_source_v1(
            value, context, expected_context, window, issued = FALSE))
      }, operation = "Frequency source-window delivery acknowledgement")
    if (!identical(acknowledged$state, "ack")) {
      stop("Frequency source delivery was not acknowledged.", call. = FALSE)
    }
    final <- identical(index + 1L, context$geometry$window_count)
    final_call <- call(
      name = "dsvertDPFrequencyFinalizeWindowDS",
      session_id = context$session_id, operation_id = operation_id,
      window_index = index,
      public_authorizations_json = authorizations)
    finalized <- .dsvert_aggregate_strict(
      context$conns[context$secondary], final_call,
      operation = "Frequency secondary window finalization",
      .aggregate = .aggregate)[[context$secondary]]
    finalized <- .dsvert_dp_frequency_execution_final_v1(
      finalized, context, window, final)
    if (final) release <- .dsvert_dp_frequency_execution_release_v1(
      finalized$release, context)
  }
  if (is.null(release)) {
    stop("Frequency execution did not commit its unique release.",
         call. = FALSE)
  }
  values <- numeric(context$geometry$d)
  chunk_hashes <- character(context$geometry$chunk_count)
  window_hashes <- character(context$geometry$window_count)
  for (index in seq_len(context$geometry$window_count) - 1L) {
    window <- .dsvert_dp_frequency_execution_window_v1(
      context$geometry, index)
    replayed <- .dsvert_aggregate_strict(
      context$conns[context$secondary], call(
        name = "dsvertDPFrequencyReplayDS",
        session_id = context$session_id, operation_id = operation_id,
        window_index = index), operation = "Frequency signed window replay",
      .aggregate = .aggregate)[[context$secondary]]
    replayed <- .dsvert_dp_frequency_execution_replay_v1(
      replayed, context, window, release)
    value_range <- seq.int(window$offset + 1L, window$offset + window$count)
    hash_range <- seq.int(
      window$first_chunk + 1L, window$first_chunk + window$chunks)
    values[value_range] <- replayed$values
    chunk_hashes[hash_range] <- replayed$chunk_hashes
    window_hashes[[index + 1L]] <- replayed$window_sha256
  }
  release_chunks <- unlist(release$final_chunk_hashes, use.names = FALSE)
  release_windows <- unlist(release$final_window_hashes, use.names = FALSE)
  if (!identical(chunk_hashes, release_chunks) ||
      !identical(window_hashes, release_windows) ||
      !identical(.dsvert_dp_frequency_execution_merkle_v1(chunk_hashes),
                 release$final_vector_root)) {
    stop("Frequency replay disagrees with its signed release commitment.",
         call. = FALSE)
  }
  proof <- list(
    session_id = context$session_id, claim = context$compiled$claim,
    config = context$compiled$config, receipts = context$compiled$receipts,
    contract = context$compiled$contract,
    worker_static = context$compiled$worker_static,
    authorities = context$authorities,
    authorizations = context$authorizations, release = release)
  list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_RESULT_VERSION,
    operation_id = operation_id,
    variable_name = context$compiled$config$factor_domain$variable_name,
    levels = unlist(
      context$compiled$config$factor_domain$levels, use.names = FALSE),
    values = unname(values), source_owner = context$source,
    finalizer_peer = context$secondary, proof = proof)
}
