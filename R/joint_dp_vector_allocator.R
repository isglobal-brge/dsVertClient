# Client relay for the server-authoritative biomedical capsule allocator.  The
# relay can only forward and verify pinned signatures; it never supplies a
# proposal, privacy parameter, mechanism, peer set or capsule identity.

.DSVERT_CLIENT_VECTOR_ALLOC_PREPARE_VERSION <-
  "dsvert-joint-dp-prepare-receipt-v1"
.DSVERT_CLIENT_VECTOR_ALLOC_COMMIT_VERSION <-
  "dsvert-joint-dp-commit-receipt-v1"
.DSVERT_CLIENT_VECTOR_ALLOC_AUTHORIZE_VERSION <-
  "dsvert-joint-dp-authorize-receipt-v1"
.DSVERT_CLIENT_VECTOR_ALLOC_OPEN_VERSION <-
  "dsvert-joint-dp-opening-token-v1"
.DSVERT_CLIENT_VECTOR_ALLOC_MAX_BYTES <- 64L * 1024L

.dsvert_vector_allocation_fields <- function(version) {
  prefix <- c(
    "version", "phase", "consortium_id", "peer_name",
    "peer_identity_pk", "capsule_id", "query_id")
  fields <- switch(version,
    "dsvert-joint-dp-prepare-receipt-v1" = c(
      prefix, "privacy_epoch", "noise_key_id", "mechanism_hash",
      "allocation_index", "epsilon", "delta", "previous_chain",
      "snapshot_binding", "seed_commitment", "rollback_mode", "signature"),
    "dsvert-joint-dp-commit-receipt-v1" = c(
      prefix, "allocation_index", "previous_chain", "new_chain",
      "joint_record_hash", "prepare_set_hash", "seed_commitment",
      "external_anchor", "signature"),
    "dsvert-joint-dp-authorize-receipt-v1" = c(
      prefix, "allocation_index", "new_chain", "joint_record_hash",
      "commit_set_hash", "peer_commit_stored", "signature"),
    "dsvert-joint-dp-opening-token-v1" = c(
      prefix, "allocation_index", "new_chain", "joint_record_hash",
      "authorization_set_hash", "seed_commitment", "release_scope",
      "capability_available", "signature"),
    NULL)
  if (is.null(fields)) {
    stop("Unknown biomedical allocation receipt version", call. = FALSE)
  }
  fields
}

.dsvert_vector_allocation_verify <- function(
    receipt_json, version, phase, peer, context, capsule_id) {
  receipt <- .dsvert_joint_dp_client_decode(
    receipt_json, "biomedical allocation receipt",
    .DSVERT_CLIENT_VECTOR_ALLOC_MAX_BYTES)
  hex <- function(value) .dsvert_vector_hex(value)
  chain <- function(value) identical(value, "GENESIS") || hex(value)
  index <- function(value, minimum = 0) {
    .dsvert_vector_string(value, "^(0|[1-9][0-9]*)$", 32L) &&
      is.finite(suppressWarnings(as.numeric(value))) &&
      as.numeric(value) >= minimum && as.numeric(value) <= 2^53 - 1
  }
  decimal <- function(value, minimum, maximum, open = FALSE) {
    number <- suppressWarnings(as.numeric(value))
    .dsvert_vector_string(value, maximum_bytes = 128L) &&
      is.finite(number) && number <= maximum &&
      if (isTRUE(open)) number > minimum else number >= minimum
  }
  valid <- .dsvert_vector_exact_names(
    receipt, .dsvert_vector_allocation_fields(version)) &&
    identical(receipt$version, version) &&
    identical(receipt$phase, phase) &&
    identical(receipt$peer_name, peer) &&
    identical(receipt$peer_identity_pk, unname(context$pinset[[peer]])) &&
    .dsvert_vector_string(
      receipt$consortium_id, "^jdpc1_[0-9a-f]{64}$", 70L) &&
    identical(receipt$capsule_id, capsule_id) &&
    identical(receipt$query_id, capsule_id) &&
    hex(receipt$capsule_id) && index(receipt$allocation_index)
  phase_valid <- if (identical(
      version, .DSVERT_CLIENT_VECTOR_ALLOC_PREPARE_VERSION)) {
    root <- context$status[[peer]]$noise_root
    index(receipt$privacy_epoch, 1) && is.list(root) &&
      identical(as.numeric(receipt$privacy_epoch),
                as.numeric(root$privacy_epoch)) &&
      identical(receipt$noise_key_id, root$key_id) &&
      .dsvert_vector_string(receipt$noise_key_id, maximum_bytes = 256L) &&
      hex(receipt$mechanism_hash) &&
      decimal(receipt$epsilon, 0, 8, TRUE) &&
      decimal(receipt$delta, 0, 1) &&
      chain(receipt$previous_chain) && hex(receipt$snapshot_binding) &&
      hex(receipt$seed_commitment) && receipt$rollback_mode %in% c(
        "cross_signed_two_peer",
        "cross_signed_two_peer_plus_external_cas")
  } else if (identical(
      version, .DSVERT_CLIENT_VECTOR_ALLOC_COMMIT_VERSION)) {
    chain(receipt$previous_chain) && hex(receipt$new_chain) &&
      hex(receipt$joint_record_hash) && hex(receipt$prepare_set_hash) &&
      hex(receipt$seed_commitment) &&
      is.logical(receipt$external_anchor) &&
      length(receipt$external_anchor) == 1L &&
      !is.na(receipt$external_anchor)
  } else if (identical(
      version, .DSVERT_CLIENT_VECTOR_ALLOC_AUTHORIZE_VERSION)) {
    hex(receipt$new_chain) && hex(receipt$joint_record_hash) &&
      hex(receipt$commit_set_hash) &&
      identical(receipt$peer_commit_stored, TRUE)
  } else {
    hex(receipt$new_chain) && hex(receipt$joint_record_hash) &&
      hex(receipt$authorization_set_hash) &&
      hex(receipt$seed_commitment) &&
      identical(receipt$release_scope, "joint_mpc_single_opening") &&
      identical(receipt$capability_available, FALSE)
  }
  if (!isTRUE(valid) || !isTRUE(phase_valid)) {
    stop("A pinned peer returned a misbound biomedical allocation receipt",
         call. = FALSE)
  }
  unsigned <- receipt[setdiff(names(receipt), "signature")]
  message <- charToRaw(paste0(
    .DSVERT_CLIENT_JOINT_DP_RECEIPT_DOMAIN,
    .dsvert_joint_dp_client_json(unsigned)))
  public <- .dsvert_joint_dp_client_b64url(
    unname(context$pinset[[peer]]), 32L,
    "biomedical allocation identity public key")
  signature <- .dsvert_joint_dp_client_b64url(
    receipt$signature, 64L, "biomedical allocation signature")
  verified <- tryCatch(openssl::ed25519_verify(
    message, signature, openssl::read_ed25519_pubkey(public)),
    error = function(error) FALSE)
  if (!isTRUE(verified)) {
    stop("A pinned peer returned an invalid biomedical allocation signature",
         call. = FALSE)
  }
  receipt
}

.dsvert_vector_allocation_set <- function(
    responses, version, phase, context, capsule_id) {
  peers <- context$designated
  if (!is.list(responses) || !identical(
        sort(names(responses), method = "radix"), peers)) {
    stop("Biomedical allocation receipts do not cover both designated peers",
         call. = FALSE)
  }
  receipts <- stats::setNames(lapply(peers, function(peer) {
    .dsvert_vector_allocation_verify(
      responses[[peer]], version, phase, peer, context, capsule_id)
  }), peers)
  stable <- switch(version,
    "dsvert-joint-dp-prepare-receipt-v1" = c(
      "consortium_id", "capsule_id", "query_id", "mechanism_hash",
      "allocation_index", "epsilon", "delta", "previous_chain"),
    "dsvert-joint-dp-commit-receipt-v1" = c(
      "consortium_id", "capsule_id", "query_id", "allocation_index",
      "previous_chain", "new_chain", "joint_record_hash",
      "prepare_set_hash"),
    "dsvert-joint-dp-authorize-receipt-v1" = c(
      "consortium_id", "capsule_id", "query_id", "allocation_index",
      "new_chain", "joint_record_hash", "commit_set_hash"),
    "dsvert-joint-dp-opening-token-v1" = c(
      "consortium_id", "capsule_id", "query_id", "allocation_index",
      "new_chain", "joint_record_hash", "authorization_set_hash",
      "release_scope", "capability_available"))
  consensus <- vapply(receipts, function(value) {
    .dsvert_joint_dp_client_json(value[stable])
  }, character(1L))
  if (length(unique(consensus)) != 1L) {
    stop("The designated peers disagree on the biomedical allocation",
         call. = FALSE)
  }
  list(receipts = receipts, json = responses[peers],
       reference = receipts[[1L]])
}

.dsvert_joint_dp_vector_allocation <- function(
    manifest_json, capsule_id, context,
    .aggregate = DSI::datashield.aggregate) {
  peers <- context$designated
  if (!is.character(peers) || length(peers) != 2L || anyNA(peers) ||
      anyDuplicated(peers) ||
      !identical(peers, sort(peers, method = "radix")) ||
      !identical(names(context$conns), peers) ||
      !.dsvert_vector_hex(capsule_id)) {
    stop("Invalid biomedical allocation relay context", call. = FALSE)
  }

  proof_calls <- stats::setNames(lapply(peers, function(peer) call(
    name = "dsvertJointDPVectorAllocationProofDS",
    manifest_json = manifest_json)), peers)
  proof <- .dsvert_vector_try_phase(.dsvert_fanout_by_site(
    context$conns, proof_calls,
    operation = "joint-DP allocation proof replay",
    .aggregate = .aggregate))
  replayed <- isTRUE(proof$ok)
  if (replayed) {
    openings <- .dsvert_vector_allocation_set(
      proof$value, .DSVERT_CLIENT_VECTOR_ALLOC_OPEN_VERSION,
      "open_authorized", context, capsule_id)
  } else {
    leader <- peers[[1L]]
    follower <- peers[[2L]]
    leader_raw <- .dsvert_fanout_by_site(
      context$conns[leader], stats::setNames(list(call(
        name = "dsvertJointDPVectorAllocationPrepareDS",
        manifest_json = manifest_json)), leader),
      operation = "joint-DP allocation leader prepare",
      .aggregate = .aggregate)
    leader_receipt <- .dsvert_vector_allocation_verify(
      leader_raw[[leader]], .DSVERT_CLIENT_VECTOR_ALLOC_PREPARE_VERSION,
      "prepared", leader, context, capsule_id)
    follower_raw <- .dsvert_fanout_by_site(
      context$conns[follower], stats::setNames(list(call(
        name = "dsvertJointDPVectorAllocationPrepareDS",
        manifest_json = manifest_json,
        leader_prepare_json = leader_raw[[leader]])), follower),
      operation = "joint-DP allocation follower prepare",
      .aggregate = .aggregate)
    prepares <- .dsvert_vector_allocation_set(
      c(leader_raw, follower_raw),
      .DSVERT_CLIENT_VECTOR_ALLOC_PREPARE_VERSION, "prepared",
      context, capsule_id)

    commit_calls <- stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertJointDPVectorAllocationCommitDS",
      first_prepare_json = prepares$json[[peers[[1L]]]],
      second_prepare_json = prepares$json[[peers[[2L]]]])), peers)
    commits <- .dsvert_vector_allocation_set(.dsvert_fanout_by_site(
      context$conns, commit_calls, operation = "joint-DP allocation commit",
      .aggregate = .aggregate),
      .DSVERT_CLIENT_VECTOR_ALLOC_COMMIT_VERSION, "committed",
      context, capsule_id)

    authorize_calls <- stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertJointDPVectorAllocationAuthorizeDS",
      first_commit_json = commits$json[[peers[[1L]]]],
      second_commit_json = commits$json[[peers[[2L]]]])), peers)
    authorizations <- .dsvert_vector_allocation_set(
      .dsvert_fanout_by_site(
        context$conns, authorize_calls,
        operation = "joint-DP allocation authorization",
        .aggregate = .aggregate),
      .DSVERT_CLIENT_VECTOR_ALLOC_AUTHORIZE_VERSION, "authorized",
      context, capsule_id)

    open_calls <- stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertJointDPVectorAllocationOpenDS",
      first_authorization_json = authorizations$json[[peers[[1L]]]],
      second_authorization_json = authorizations$json[[peers[[2L]]]])),
      peers)
    openings <- .dsvert_vector_allocation_set(.dsvert_fanout_by_site(
      context$conns, open_calls, operation = "joint-DP allocation opening",
      .aggregate = .aggregate),
      .DSVERT_CLIENT_VECTOR_ALLOC_OPEN_VERSION, "open_authorized",
      context, capsule_id)
  }

  context$allocation_openings <- openings$json
  context$allocation_certificate <- list(
    version = "dsvert-joint-dp-biomedical-vector-allocation-client-v1",
    capsule_id = capsule_id,
    consortium_id = openings$reference$consortium_id,
    allocation_index = openings$reference$allocation_index,
    opening_set_hash = .dsvert_vector_hash(openings$receipts),
    designated_noise_peers = peers,
    proof_replayed = replayed,
    relay_is_authority = FALSE, capability_available = FALSE,
    data_access = FALSE, history_gate = TRUE,
    request_limit = FALSE, operation_limit = TRUE)
  context
}
