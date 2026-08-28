# Offline-verifiable public provenance for bounded Gaussian capsule views.
# The certificate contains only the released DP coordinates needed by the
# selected Gaussian block, signed public receipts and public commitments.

.DSVERT_DP_GAUSSIAN_CERTIFICATE_VERSION <-
  "dsvert-dp-gaussian-provenance-certificate-v3"
.DSVERT_DP_GAUSSIAN_CERTIFICATE_COMMITMENT_VERSION <-
  "dsvert-dp-gaussian-artifact-membership-v1"
.DSVERT_DP_GAUSSIAN_CERTIFICATE_EVIDENCE_VERSION <-
  "dsvert-dp-gaussian-signed-evidence-v1"
.DSVERT_DP_GAUSSIAN_CERTIFICATE_CHUNK_VERSION <-
  "dsvert-dp-gaussian-public-chunk-evidence-v1"
.DSVERT_DP_GAUSSIAN_SYNOPSIS_CERTIFICATE_VERSION <-
  "dsvert-dp-gaussian-synopsis-provenance-certificate-v1"
.DSVERT_DP_GAUSSIAN_SYNOPSIS_EVIDENCE_VERSION <-
  "dsvert-dp-gaussian-synopsis-signed-evidence-v1"

.dsvert_dp_gaussian_trust_cache <- new.env(parent = emptyenv())

.dsvert_dp_gaussian_certificate_hash <- function(value) {
  .dsvert_dp_capsule_manifest_hash(value)
}

.dsvert_dp_gaussian_synopsis_no_legacy_fields <- function(value) {
  forbidden <- paste(
    c("release_instance", "privacy_epoch", "noise_key", "capsule",
      "lifetime", "ledger", "reservation", "request", "rate", "catalog",
      "quota"), collapse = "|")
  visit <- function(node, path = "certificate") {
    if (!is.list(node)) return(invisible(TRUE))
    fields <- names(node)
    if (!is.null(fields)) {
      invalid <- fields[grepl(forbidden, fields, ignore.case = TRUE)]
      if (length(invalid)) {
        stop(
          "The Gaussian Synopsis certificate contains a forbidden legacy ",
          "field at ", path, ": ", invalid[[1L]], call. = FALSE)
      }
      for (index in seq_along(node)) {
        visit(node[[index]], paste0(path, "$", fields[[index]]))
      }
    } else {
      for (index in seq_along(node)) {
        visit(node[[index]], paste0(path, "[[", index, "]]"))
      }
    }
    invisible(TRUE)
  }
  visit(value)
}

.dsvert_dp_gaussian_synopsis_evidence_json <- function(value, what) {
  json <- .dsvert_joint_dp_client_json(value)
  .dsvert_dp_synopsis_client_json(
    json, what, .DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES)
  json
}

.dsvert_dp_gaussian_synopsis_evidence_decode <- function(value, what) {
  if (!.dsvert_dp_is_string(value)) {
    stop("Invalid Gaussian Synopsis ", what, call. = FALSE)
  }
  .dsvert_dp_synopsis_client_json(
    value, paste("Gaussian Synopsis", what),
    .DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES)
}

.dsvert_dp_gaussian_pinset <- function(value, what = "certificate pinset") {
  if (is.list(value) && !is.null(names(value))) {
    valid <- all(vapply(value, .dsvert_dp_is_string, logical(1L)))
    if (isTRUE(valid)) value <- unlist(value, use.names = TRUE)
  }
  if (!is.character(value) || !length(value) || is.null(names(value)) ||
      anyNA(value) || anyNA(names(value)) || any(!nzchar(names(value))) ||
      anyDuplicated(names(value)) || anyDuplicated(value)) {
    stop("Invalid ", what, call. = FALSE)
  }
  value <- value[order(names(value), method = "radix")]
  normalized <- vapply(value, function(key) {
    result <- .dsvert_dp_normalize_identity_pk(key)
    if (is.null(result)) stop("Invalid ", what, call. = FALSE)
    result
  }, character(1L))
  names(normalized) <- names(value)
  normalized
}

.dsvert_dp_gaussian_cache_pinset <- function(pinset) {
  pinset <- .dsvert_dp_gaussian_pinset(pinset)
  key <- .dsvert_dp_pinset_hash(pinset)
  assign(key, pinset, envir = .dsvert_dp_gaussian_trust_cache)
  invisible(key)
}

.dsvert_dp_gaussian_scaled_text <- function(values, scale) {
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0) || !.dsvert_dp_is_number(scale, 0, 2^62,
                                               lower_open = TRUE)) {
    stop("Invalid public DP Gaussian lattice values", call. = FALSE)
  }
  scaled <- values * scale
  if (anyNA(scaled) || any(!is.finite(scaled)) || any(scaled < 0) ||
      any(scaled > 2^53) || any(scaled != floor(scaled))) {
    stop("A public DP Gaussian coordinate is off its signed lattice",
         call. = FALSE)
  }
  unname(sprintf("%.0f", scaled))
}

.dsvert_dp_gaussian_artifact_proof <- function(index, analysis_id) {
  entries <- tryCatch(index$gaussian_models, error = function(error) NULL)
  if (!is.list(entries) || is.null(names(entries)) ||
      anyNA(names(entries)) || anyDuplicated(names(entries)) ||
      !analysis_id %in% names(entries)) {
    stop("The signed manifest has no Gaussian artifact commitment",
         call. = FALSE)
  }
  ordered <- sort(names(entries), method = "radix")
  if (!identical(names(entries), ordered)) {
    stop("The Gaussian artifact commitment index is not canonical",
         call. = FALSE)
  }
  hashes <- vapply(entries, function(entry) {
    if (!is.list(entry) || !.dsvert_dp_capsule_source_hex(
          entry$node_sha256)) {
      stop("Invalid Gaussian artifact commitment", call. = FALSE)
    }
    entry$node_sha256
  }, character(1L))
  position <- match(analysis_id, ordered)
  list(
    version = .DSVERT_DP_GAUSSIAN_CERTIFICATE_COMMITMENT_VERSION,
    context = index$context,
    entry = entries[[analysis_id]],
    leaf_index = as.numeric(position),
    leaf_count = as.numeric(length(entries)),
    merkle_proof = .dsvert_dp_capsule_artifact_merkle_proof_client(
      hashes, position))
}

.dsvert_dp_gaussian_chunk_evidence <- function(
    release, block, release_receipt) {
  hashes <- unlist(release_receipt$final_chunk_hashes, use.names = FALSE)
  chunk_count <- as.integer(release_receipt$chunk_count)
  if (!is.character(hashes) || length(hashes) != chunk_count ||
      anyNA(hashes) || any(!vapply(
        hashes, .dsvert_dp_capsule_source_hex, logical(1L)))) {
    stop("The signed final DP vector chunk set is invalid", call. = FALSE)
  }
  first <- as.integer((block$start - 1L) %/%
                        .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES)
  last <- as.integer((block$end - 1L) %/%
                       .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES)
  indices <- seq.int(first, last)
  scale <- as.numeric(release_receipt$output_lattice_scale)
  lapply(indices, function(index) {
    geometry <- .dsvert_vector_chunk_geometry(list(
      coordinate_count = as.integer(release$coordinate_count),
      chunk_coordinates = .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES,
      chunk_count = chunk_count), index)
    positions <- seq.int(geometry$offset + 1L,
                         geometry$offset + geometry$count)
    chunk <- .dsvert_joint_dp_client_canonical(list(
      version = .DSVERT_CLIENT_VECTOR_CHUNK_VERSION,
      capsule_id = release_receipt$capsule_id,
      release_contract_hash = release_receipt$release_contract_hash,
      chunk_index = as.numeric(index),
      chunk_count = as.numeric(chunk_count),
      coordinate_offset = as.numeric(geometry$offset),
      coordinates_in_chunk = as.numeric(geometry$count),
      output_lattice_bits =
        as.numeric(release_receipt$output_lattice_bits),
      output_lattice_scale = scale,
      scaled_values = as.list(.dsvert_dp_gaussian_scaled_text(
        release$values[positions], scale)),
      value_encoding = "nonnegative-decimal-integer-common-lattice-v1",
      postprocessing = paste0(
        "signed-Ring128-decode-then-fixed-public-coordinate-clamp-v1"),
      source_values_exposed = FALSE,
      preclamp_values_exposed = FALSE))
    observed <- .dsvert_vector_hash(chunk)
    if (!identical(observed, hashes[[index + 1L]])) {
      stop("The Gaussian block does not match the signed final DP vector",
           call. = FALSE)
    }
    list(
      version = .DSVERT_DP_GAUSSIAN_CERTIFICATE_CHUNK_VERSION,
      chunk_hash = observed, chunk = chunk,
      merkle_proof = .dsvert_vector_merkle_proof(hashes, index))
  })
}

.dsvert_dp_gaussian_certificate_build <- function(
    context, artifact, block, coordinates) {
  bundle <- context$manifest_bundle
  release <- context$release
  provenance <- release$signed_provenance
  required_bundle <- c(
    "schema_json", "schema_sha256", "workload_contract_sha256",
    "manifest_sha256", "capsule_id", "artifact_commitments",
    "artifact_commitment_count", "artifact_commitments_root",
    "manifest_build_receipts")
  required_provenance <- c(
    "version", "ordered_peer_pinset", "peer_pinset_sha256",
    "designated_noise_peers", "release_instance_id", "release_instance",
    "prepare_receipts", "release_receipts",
    "finalization_receipts", "protected_shares_included",
    "preclamp_values_included", "source_values_included")
  if (!is.list(bundle) || !all(required_bundle %in% names(bundle)) ||
      !.dsvert_dp_has_exact_names(provenance, required_provenance) ||
      !identical(provenance$version,
                 "dsvert-joint-dp-vector-public-provenance-v1") ||
      !identical(provenance$protected_shares_included, FALSE) ||
      !identical(provenance$preclamp_values_included, FALSE) ||
      !identical(provenance$source_values_included, FALSE)) {
    stop("The Gaussian result lacks signed public capsule provenance",
         call. = FALSE)
  }
  pinset <- .dsvert_dp_gaussian_pinset(provenance$ordered_peer_pinset)
  if (!identical(.dsvert_dp_pinset_hash(pinset),
                 provenance$peer_pinset_sha256)) {
    stop("The Gaussian provenance pinset is inconsistent", call. = FALSE)
  }
  designated <- unlist(provenance$designated_noise_peers, use.names = FALSE)
  designated <- enc2utf8(unname(as.character(designated)))
  release_receipts <- provenance$release_receipts
  if (!is.list(release_receipts) ||
      !setequal(names(release_receipts), designated)) {
    stop("The Gaussian provenance lacks both vector releases",
         call. = FALSE)
  }
  release_reference <- release_receipts[[designated[[1L]]]]
  if (!identical(release_reference$final_vector_root,
                 release$final_vector_root)) {
    stop("The Gaussian release root changed before certification",
         call. = FALSE)
  }
  commitment <- .dsvert_dp_gaussian_artifact_proof(
    bundle$artifact_commitments, artifact$analysis_id)
  if (!identical(as.numeric(commitment$leaf_count),
                 as.numeric(bundle$artifact_commitment_count)) ||
      !identical(.dsvert_dp_capsule_artifact_merkle_root_client(
        vapply(bundle$artifact_commitments$gaussian_models, `[[`,
               character(1L), "node_sha256"),
        bundle$artifact_commitments$context),
        bundle$artifact_commitments_root)) {
    stop("The Gaussian artifact commitment root is inconsistent",
         call. = FALSE)
  }
  entry <- commitment$entry
  profile <- .dsvert_dp_vector_release_profile(release, context$manifest)
  prepare_reference <- provenance$prepare_receipts[[designated[[1L]]]]
  if (!identical(provenance$release_instance_id,
                 release$release_instance_id) ||
      !identical(
        .dsvert_joint_dp_client_canonical(provenance$release_instance),
        .dsvert_joint_dp_client_canonical(release$release_instance)) ||
      !identical(prepare_reference$release_instance_id,
                 provenance$release_instance_id) ||
      !identical(
        .dsvert_joint_dp_client_canonical(prepare_reference$release_instance),
        .dsvert_joint_dp_client_canonical(provenance$release_instance)) ||
      !identical(.dsvert_vector_hash(provenance$release_instance),
                 provenance$release_instance_id)) {
    stop("The Gaussian release-instance provenance is inconsistent",
         call. = FALSE)
  }
  .dsvert_vector_plan_validate(
    prepare_reference$mechanism_plan, prepare_reference$plan_sha256,
    profile, release$coordinate_count, prepare_reference$sensitivity_steps)
  if (!identical(release_reference$backend, profile$backend) ||
      !identical(release_reference$sampler, profile$sampler) ||
      !identical(release_reference$mechanism,
                 profile$release_mechanism)) {
    stop("The Gaussian release mechanism profile is misbound",
         call. = FALSE)
  }
  # Commit the exact descriptor representation carried by the signed
  # manifest/layout. The validated working copy may normalize JSON arrays to
  # character vectors, which is semantically equivalent but hashes
  # differently in R.
  signed_artifact <- block$descriptor
  artifact_sha256 <- .dsvert_dp_gaussian_certificate_hash(signed_artifact)
  block_public <- block[c(
    "family", "key", "start", "end", "length", "owner_peer", "dataset")]
  names(block_public)[names(block_public) == "key"] <- "analysis_id"
  if (!identical(entry$descriptor_sha256, artifact_sha256) ||
      !identical(.dsvert_joint_dp_client_canonical(
        entry[setdiff(names(entry), c("version", "node_sha256",
                                     "descriptor_sha256"))]),
        .dsvert_joint_dp_client_canonical(block_public))) {
    stop("The Gaussian artifact is not the committed coordinate block",
         call. = FALSE)
  }
  chunk_evidence <- .dsvert_dp_gaussian_chunk_evidence(
    release, block, release_reference)
  included_coordinates <- sum(vapply(chunk_evidence, function(item) {
    as.numeric(item$chunk$coordinates_in_chunk)
  }, numeric(1L)))
  chunk_proof_geometry <- list(
    granularity = "whole_intersecting_public_chunks_v1",
    block_coordinates = as.numeric(block$length),
    included_public_coordinates = included_coordinates,
    overfetch_coordinates = included_coordinates - as.numeric(block$length),
    intersecting_chunk_count = as.numeric(length(chunk_evidence)))
  block_scaled <- .dsvert_dp_gaussian_scaled_text(
    coordinates, as.numeric(release_reference$output_lattice_scale))
  block_values_sha256 <- .dsvert_dp_gaussian_certificate_hash(list(
    protocol = "dsvert-dp-gaussian-released-block-values-v1",
    start = as.numeric(block$start), end = as.numeric(block$end),
    scaled_values = as.list(block_scaled)))
  query_contract_sha256 <- .dsvert_dp_gaussian_certificate_hash(list(
    protocol = "dsvert-dp-gaussian-query-contract-v1",
    capsule_id = bundle$capsule_id, analysis_id = artifact$analysis_id,
    release_instance_id = provenance$release_instance_id,
    release_instance = provenance$release_instance,
    artifact_sha256 = artifact_sha256,
    artifact_node_sha256 = entry$node_sha256,
    release_contract_hash = release_reference$release_contract_hash,
    capsule_mechanism = profile$mechanism,
    mechanism_plan_sha256 = prepare_reference$plan_sha256,
    backend = profile$backend, sampler = profile$sampler,
    release_mechanism = profile$release_mechanism,
    final_vector_root = release_reference$final_vector_root,
    coordinate_order_sha256 =
      commitment$context$coordinate_order_sha256,
    block_values_sha256 = block_values_sha256))
  evidence <- list(
    version = .DSVERT_DP_GAUSSIAN_CERTIFICATE_EVIDENCE_VERSION,
    signed_schema_json = bundle$schema_json,
    manifest_build_receipts = bundle$manifest_build_receipts,
    vector_prepare_receipts = provenance$prepare_receipts,
    vector_release_receipts = provenance$release_receipts,
    vector_finalization_receipts = provenance$finalization_receipts)
  unsigned <- list(
    version = .DSVERT_DP_GAUSSIAN_CERTIFICATE_VERSION,
    analysis_id = artifact$analysis_id,
    dataset = artifact$dataset,
    owner_peer = artifact$owner_peer,
    cross_owner_state = artifact$cross_owner_state,
    query_contract_sha256 = query_contract_sha256,
    artifact_sha256 = artifact_sha256,
    artifact = signed_artifact,
    artifact_commitment = commitment,
    block = block_public,
    block_values_sha256 = block_values_sha256,
    chunk_evidence = chunk_evidence,
    chunk_proof_geometry = chunk_proof_geometry,
    capsule_id = bundle$capsule_id,
    release_instance_id = provenance$release_instance_id,
    release_instance = provenance$release_instance,
    manifest_sha256 = bundle$manifest_sha256,
    schema_sha256 = bundle$schema_sha256,
    workload_contract_sha256 = bundle$workload_contract_sha256,
    release_contract_hash = release_reference$release_contract_hash,
    source_contract_hash =
      provenance$prepare_receipts[[designated[[1L]]]]$source_contract_hash,
    transcript_hash = release_reference$transcript_hash,
    final_vector_root = release_reference$final_vector_root,
    coordinate_order_sha256 =
      commitment$context$coordinate_order_sha256,
    coordinate_count = as.numeric(release$coordinate_count),
    output_lattice_bits =
      as.numeric(release_reference$output_lattice_bits),
    output_lattice_scale =
      as.numeric(release_reference$output_lattice_scale),
    capsule_mechanism = profile$mechanism,
    mechanism_plan_sha256 = prepare_reference$plan_sha256,
    backend = profile$backend,
    sampler = profile$sampler,
    mechanism = release_reference$mechanism,
    epsilon = release_reference$epsilon,
    delta = release_reference$delta,
    implementation_delta = paste0(
      release_reference$implementation_delta_numerator, "/",
      release_reference$implementation_delta_denominator),
    peer_context = list(
      ordered_peer_pinset = as.list(pinset),
      peer_pinset_sha256 = provenance$peer_pinset_sha256,
      designated_noise_peers = as.list(designated),
      privacy_epochs = as.list(stats::setNames(vapply(
        bundle$manifest_build_receipts, function(receipt) {
          as.numeric(receipt$privacy_epoch)
        }, numeric(1L)), names(bundle$manifest_build_receipts))),
      noise_key_ids = as.list(stats::setNames(vapply(
        bundle$manifest_build_receipts, `[[`, character(1L),
        "noise_key_id"), names(bundle$manifest_build_receipts))),
      noise_provider_ids = as.list(stats::setNames(vapply(
        provenance$release_instance$peer_noise_roots,
        `[[`, character(1L), "provider_id"), designated)),
      release_domain_generations = as.list(stats::setNames(vapply(
        provenance$release_instance$peer_noise_roots, function(root) {
          as.numeric(root$release_domain_generation)
        }, numeric(1L)), designated)),
      release_domain_ids = as.list(stats::setNames(vapply(
        provenance$release_instance$peer_noise_roots,
        `[[`, character(1L), "release_domain_id"), designated))),
    signed_evidence = evidence,
    public_dp_coordinates_only = TRUE,
    protected_shares_included = FALSE,
    preclamp_values_included = FALSE,
    patient_derived_identifiers_included = FALSE)
  certificate <- c(unsigned, list(
    certificate_sha256 = .dsvert_dp_gaussian_certificate_hash(unsigned)))
  .dsvert_dp_gaussian_cache_pinset(pinset)
  certificate
}

.dsvert_dp_gaussian_artifact_root_from_proof <- function(commitment) {
  required <- c(
    "version", "context", "entry", "leaf_index", "leaf_count",
    "merkle_proof")
  if (!.dsvert_dp_has_exact_names(commitment, required) ||
      !identical(commitment$version,
                 .DSVERT_DP_GAUSSIAN_CERTIFICATE_COMMITMENT_VERSION) ||
      !.dsvert_dp_is_integer(commitment$leaf_index, 1) ||
      !.dsvert_dp_is_integer(commitment$leaf_count, 1) ||
      commitment$leaf_index > commitment$leaf_count ||
      !is.list(commitment$merkle_proof)) {
    stop("Invalid Gaussian artifact membership proof", call. = FALSE)
  }
  entry <- commitment$entry
  expected_entry <- c(
    "version", "family", "analysis_id", "dataset", "owner_peer",
    "start", "end", "length", "descriptor_sha256", "node_sha256")
  if (!.dsvert_dp_has_exact_names(entry, expected_entry) ||
      !identical(entry$version,
                 .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_ENTRY_VERSION) ||
      !identical(entry$family, "gaussian_models") ||
      !.dsvert_dp_capsule_source_hex(entry$descriptor_sha256) ||
      !.dsvert_dp_capsule_source_hex(entry$node_sha256) ||
      !.dsvert_dp_is_integer(entry$start, 1) ||
      !.dsvert_dp_is_integer(entry$end, entry$start) ||
      !.dsvert_dp_is_integer(entry$length, 1) ||
      entry$end - entry$start + 1 != entry$length) {
    stop("Invalid Gaussian artifact commitment entry", call. = FALSE)
  }
  unsigned_entry <- entry[setdiff(names(entry), "node_sha256")]
  observed <- .dsvert_dp_gaussian_certificate_hash(list(
    protocol = "dsvert-biomedical-capsule-artifact-merkle-leaf-v1",
    context = commitment$context, entry = unsigned_entry))
  if (!identical(observed, entry$node_sha256)) {
    stop("The Gaussian artifact leaf commitment is invalid", call. = FALSE)
  }
  expected_depth <- if (commitment$leaf_count == 1) 0L else
    ceiling(log2(commitment$leaf_count))
  if (length(commitment$merkle_proof) != expected_depth) {
    stop("The Gaussian artifact proof has the wrong depth", call. = FALSE)
  }
  root <- observed
  for (item in commitment$merkle_proof) {
    if (!.dsvert_dp_has_exact_names(item, c("side", "sha256")) ||
        !item$side %in% c("left", "right") ||
        !.dsvert_dp_capsule_source_hex(item$sha256)) {
      stop("Invalid Gaussian artifact Merkle path", call. = FALSE)
    }
    root <- if (identical(item$side, "left")) {
      .dsvert_dp_capsule_artifact_merkle_parent_client(item$sha256, root)
    } else {
      .dsvert_dp_capsule_artifact_merkle_parent_client(root, item$sha256)
    }
  }
  root
}

.dsvert_dp_gaussian_certificate_build_receipts <- function(
    certificate, pinset, context) {
  receipts <- certificate$signed_evidence$manifest_build_receipts
  peers <- names(pinset)
  expected <- c(
    "version", "phase", "peer_name", "peer_identity_pk",
    "peer_pinset_sha256", "schema_sha256", "workload_contract_sha256",
    "manifest_sha256", "manifest_bytes", "capsule_id", "privacy_epoch",
    "noise_key_id", "artifact_commitment_count",
    "artifact_commitments_root", "durable_memoization",
    "deterministic_replay", "data_access", "operation_limit",
    "request_limit", "history_can_deny_operation", "signature")
  if (!is.list(receipts) || !setequal(names(receipts), peers)) {
    stop("The certificate build receipts do not cover the pinset",
         call. = FALSE)
  }
  receipts <- receipts[peers]
  for (peer in peers) {
    receipt <- receipts[[peer]]
    valid <- .dsvert_dp_has_exact_names(receipt, expected) &&
      identical(receipt$version,
                .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_BUILD_VERSION) &&
      identical(receipt$phase,
                "server_authoritative_manifest_memoized") &&
      identical(receipt$peer_name, peer) &&
      identical(receipt$peer_identity_pk, unname(pinset[[peer]])) &&
      identical(receipt$peer_pinset_sha256,
                certificate$peer_context$peer_pinset_sha256) &&
      identical(receipt$schema_sha256, certificate$schema_sha256) &&
      identical(receipt$workload_contract_sha256,
                certificate$workload_contract_sha256) &&
      identical(receipt$manifest_sha256, certificate$manifest_sha256) &&
      identical(receipt$capsule_id, certificate$capsule_id) &&
      identical(as.numeric(receipt$privacy_epoch),
                as.numeric(certificate$peer_context$
                  privacy_epochs[[peer]])) &&
      identical(as.numeric(receipt$artifact_commitment_count),
                as.numeric(certificate$artifact_commitment$leaf_count)) &&
      identical(receipt$artifact_commitments_root,
                context$artifact_commitments_root) &&
      .dsvert_dp_is_string(receipt$noise_key_id) &&
      identical(receipt$durable_memoization, TRUE) &&
      identical(receipt$deterministic_replay, TRUE) &&
      identical(receipt$data_access, FALSE) &&
      identical(receipt$operation_limit, FALSE) &&
      identical(receipt$request_limit, FALSE) &&
      identical(receipt$history_can_deny_operation, FALSE)
    if (!isTRUE(valid)) {
      stop("A Gaussian certificate build receipt is misbound",
           call. = FALSE)
    }
    .dsvert_dp_capsule_manifest_verify(receipt, "build", peer,
                                       list(pinset = pinset))
  }
  receipts
}

.dsvert_dp_gaussian_certificate_schema <- function(
    certificate, pinset, context) {
  signed <- .dsvert_joint_dp_client_decode(
    certificate$signed_evidence$signed_schema_json,
    "Gaussian certificate signed schema", 8L * 1024L^2)
  if (!.dsvert_dp_has_exact_names(signed, c(
        "version", "logical_snapshot", "peer_pinset_sha256",
        "datasets", "signatures")) ||
      !is.list(signed$signatures) ||
      !setequal(names(signed$signatures), names(pinset)) ||
      !identical(signed$peer_pinset_sha256,
                 certificate$peer_context$peer_pinset_sha256) ||
      !identical(.dsvert_joint_dp_client_canonical(
        signed$logical_snapshot),
        .dsvert_joint_dp_client_canonical(context$logical_snapshot))) {
    stop("The Gaussian certificate signed schema is invalid", call. = FALSE)
  }
  unsigned <- signed[setdiff(names(signed), "signatures")]
  if (!identical(.dsvert_dp_capsule_manifest_hash(unsigned),
                 certificate$schema_sha256)) {
    stop("The Gaussian certificate schema hash is invalid", call. = FALSE)
  }
  for (peer in names(pinset)) {
    .dsvert_dp_capsule_schema_verify(
      unsigned, signed$signatures[[peer]], peer, list(pinset = pinset))
  }
  signed
}

.dsvert_dp_gaussian_certificate_vector <- function(
    certificate, pinset, index_context) {
  evidence <- certificate$signed_evidence
  designated <- unlist(
    certificate$peer_context$designated_noise_peers, use.names = FALSE)
  designated <- enc2utf8(unname(as.character(designated)))
  if (length(designated) != 2L || anyDuplicated(designated) ||
      !identical(designated, sort(designated, method = "radix")) ||
      !all(designated %in% names(pinset))) {
    stop("Invalid designated peer set in Gaussian certificate",
         call. = FALSE)
  }
  vector_context <- list(
    pinset = pinset, designated = designated, servers = names(pinset))
  prepare <- evidence$vector_prepare_receipts
  if (!is.list(prepare) || !setequal(names(prepare), designated)) {
    stop("The Gaussian certificate lacks vector prepare receipts",
         call. = FALSE)
  }
  prepare <- prepare[designated]
  backends <- vapply(prepare, function(value) {
    if (!is.list(value) || !.dsvert_vector_string(
          value$backend, maximum_bytes = 256L)) {
      stop("A Gaussian certificate prepare has no signed backend",
           call. = FALSE)
    }
    value$backend
  }, character(1L))
  if (length(unique(backends)) != 1L) {
    stop("The Gaussian certificate peers disagree on the vector backend",
         call. = FALSE)
  }
  profile <- .dsvert_vector_profile(
    index_context$capsule_mechanism,
    index_context$mechanism_selection, backend = backends[[1L]])
  prepare <- stats::setNames(lapply(designated, function(peer) {
    .dsvert_vector_verify_receipt(
      .dsvert_joint_dp_client_json(prepare[[peer]]),
      .DSVERT_CLIENT_VECTOR_PREPARE_VERSION, "vector_prepared", peer,
      vector_context, profile = profile)
  }), designated)
  stable <- c(
    "capsule_id", "release_instance_id", "release_instance",
    "release_contract_hash", "transcript_hash",
    "coordinate_count", "chunk_count", "backend", "sampler",
    "history_gate", "request_limit", "source_contract_hash",
    "coordinate_order_sha256", "lattice_transform_sha256",
    "mechanism_plan", "plan_sha256",
    "epsilon", "allocated_delta", "sensitivity_steps",
    "complete_epsilon_per_peer", "delta_aggregation",
    "capability_available",
    if (isTRUE(profile$selection_bound)) c(
      "backend_assessment", "backend_selection",
      "backend_selection_sha256", "one_joint_draw") else character())
  if (length(unique(vapply(prepare, function(receipt) {
        .dsvert_joint_dp_client_json(receipt[stable])
      }, character(1L)))) != 1L) {
    stop("The Gaussian certificate prepare receipts disagree",
         call. = FALSE)
  }
  reference <- prepare[[1L]]
  chunk_coordinates <- if (isTRUE(profile$exact_gc)) {
    as.integer(reference$mechanism_plan$maximum_chunk_coordinates)
  } else .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES
  expected_chunks <- ceiling(certificate$coordinate_count /
                               chunk_coordinates)
  valid <- identical(reference$capsule_id, certificate$capsule_id) &&
    identical(reference$release_instance_id,
              certificate$release_instance_id) &&
    identical(
      .dsvert_joint_dp_client_canonical(reference$release_instance),
      .dsvert_joint_dp_client_canonical(certificate$release_instance)) &&
    identical(.dsvert_vector_hash(certificate$release_instance),
              certificate$release_instance_id) &&
    identical(reference$release_contract_hash,
              certificate$release_contract_hash) &&
    identical(reference$transcript_hash, certificate$transcript_hash) &&
    identical(reference$source_contract_hash,
              certificate$source_contract_hash) &&
    identical(reference$coordinate_order_sha256,
              certificate$coordinate_order_sha256) &&
    identical(reference$coordinate_order_sha256,
              index_context$coordinate_order_sha256) &&
    identical(as.numeric(reference$coordinate_count),
              as.numeric(certificate$coordinate_count)) &&
    identical(as.numeric(reference$chunk_count),
              as.numeric(expected_chunks)) &&
    identical(reference$epsilon, certificate$epsilon) &&
    identical(reference$allocated_delta, certificate$delta) &&
    identical(reference$plan_sha256,
              certificate$mechanism_plan_sha256) &&
    identical(reference$backend, certificate$backend) &&
    identical(reference$sampler, certificate$sampler) &&
    identical(profile$mechanism, certificate$capsule_mechanism) &&
    identical(profile$release_mechanism, certificate$mechanism) &&
    identical(reference$complete_epsilon_per_peer,
              profile$complete_epsilon_per_peer) &&
    identical(reference$delta_aggregation,
              profile$delta_aggregation) &&
    identical(reference$capability_available, TRUE)
  if (!isTRUE(valid)) {
    stop("The Gaussian certificate vector contract is misbound",
         call. = FALSE)
  }
  .dsvert_vector_plan_validate(
    reference$mechanism_plan, reference$plan_sha256, profile,
    certificate$coordinate_count, reference$sensitivity_steps)
  contract <- list(
    capsule_id = reference$capsule_id,
    release_instance_id = reference$release_instance_id,
    release_instance = reference$release_instance,
    release_contract_hash = reference$release_contract_hash,
    transcript_hash = reference$transcript_hash,
    coordinate_count = as.integer(certificate$coordinate_count),
    chunk_coordinates = chunk_coordinates,
    chunk_count = as.integer(expected_chunks),
    source_contract_hash = reference$source_contract_hash,
    coordinate_order_sha256 = reference$coordinate_order_sha256,
    mechanism_plan = reference$mechanism_plan,
    plan_sha256 = reference$plan_sha256,
    manifest_sha256 = certificate$manifest_sha256,
    backend_assessment = reference$backend_assessment,
    backend_selection = reference$backend_selection,
    backend_selection_sha256 = reference$backend_selection_sha256,
    one_joint_draw = isTRUE(profile$exact_gc),
    capsule_mechanism = profile$mechanism,
    profile = profile,
    epsilon = reference$epsilon, delta = reference$allocated_delta,
    sensitivity_steps = reference$sensitivity_steps)
  releases <- evidence$vector_release_receipts
  if (!is.list(releases) || !setequal(names(releases), designated)) {
    stop("The Gaussian certificate lacks vector release receipts",
         call. = FALSE)
  }
  releases <- releases[designated]
  releases <- stats::setNames(lapply(designated, function(peer) {
    .dsvert_vector_verify_receipt(
      .dsvert_joint_dp_client_json(releases[[peer]]),
      .DSVERT_CLIENT_VECTOR_RELEASE_VERSION, "vector_released", peer,
      vector_context, contract)
  }), designated)
  release_stable <- c(
    "result_set_hash", "final_vector_root", "final_chunk_hashes",
    "output_lattice_bits", "output_lattice_scale", "mechanism", "epsilon",
    "delta", "implementation_delta_numerator",
    "implementation_delta_denominator", "delta_aggregation",
    "postprocessing", "intermediate_payload_exposed", "durable_replay",
    "capability_available")
  if (length(unique(vapply(releases, function(receipt) {
        .dsvert_joint_dp_client_json(receipt[release_stable])
      }, character(1L)))) != 1L) {
    stop("The Gaussian certificate release receipts disagree",
         call. = FALSE)
  }
  release <- releases[[1L]]
  hashes <- unlist(release$final_chunk_hashes, use.names = FALSE)
  accounting <- paste0(release$implementation_delta_numerator, "/",
                       release$implementation_delta_denominator)
  valid <- is.character(hashes) && length(hashes) == expected_chunks &&
    !anyNA(hashes) && all(vapply(
      hashes, .dsvert_dp_capsule_source_hex, logical(1L))) &&
    identical(.dsvert_vector_merkle_root(hashes),
              certificate$final_vector_root) &&
    identical(release$final_vector_root, certificate$final_vector_root) &&
    identical(as.numeric(release$output_lattice_bits),
              as.numeric(certificate$output_lattice_bits)) &&
    identical(as.numeric(release$output_lattice_scale),
              as.numeric(certificate$output_lattice_scale)) &&
    identical(release$mechanism, certificate$mechanism) &&
    identical(release$mechanism, profile$release_mechanism) &&
    identical(release$epsilon, certificate$epsilon) &&
    identical(release$delta, certificate$delta) &&
    identical(accounting, certificate$implementation_delta) &&
    identical(release$delta_aggregation, profile$delta_aggregation) &&
    identical(release$postprocessing,
              profile$postprocessing) &&
    identical(release$intermediate_payload_exposed, FALSE) &&
    identical(release$durable_replay, TRUE) &&
    identical(release$capability_available, TRUE)
  if (!isTRUE(valid)) {
    stop("The Gaussian certificate final vector is invalid", call. = FALSE)
  }
  acks <- evidence$vector_finalization_receipts
  if (!is.list(acks) || !setequal(names(acks), names(pinset))) {
    stop("The Gaussian certificate lacks finalization receipts",
         call. = FALSE)
  }
  acks <- acks[names(pinset)]
  for (peer in names(pinset)) {
    ack <- .dsvert_vector_verify_receipt(
      .dsvert_joint_dp_client_json(acks[[peer]]),
      .DSVERT_CLIENT_VECTOR_ACK_VERSION,
      "vector_finalized_and_compacted", peer, vector_context, contract)
    valid <- identical(ack$final_vector_root,
                       certificate$final_vector_root) &&
      identical(ack$source_intermediates_compacted, TRUE) &&
      identical(ack$sampler_intermediates_compacted, TRUE) &&
      identical(ack$final_chunks_retained, TRUE) &&
      identical(ack$durable_replay_retained, TRUE) &&
      identical(ack$idempotent, TRUE)
    if (!isTRUE(valid)) {
      stop("A Gaussian certificate finalization receipt is invalid",
           call. = FALSE)
    }
  }
  list(contract = contract, release = release, hashes = hashes)
}

.dsvert_dp_gaussian_certificate_coordinates <- function(
    certificate, vector) {
  block <- certificate$block
  expected_chunks <- seq.int(
    as.integer((block$start - 1L) %/%
                 .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES),
    as.integer((block$end - 1L) %/%
                 .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES))
  evidence <- certificate$chunk_evidence
  if (!is.list(evidence) || length(evidence) != length(expected_chunks)) {
    stop("The Gaussian certificate has the wrong intersecting chunk set",
         call. = FALSE)
  }
  geometry_summary <- certificate$chunk_proof_geometry
  geometry_fields <- c(
    "granularity", "block_coordinates", "included_public_coordinates",
    "overfetch_coordinates", "intersecting_chunk_count")
  included_coordinates <- sum(vapply(evidence, function(item) {
    if (!is.list(item) || !is.list(item$chunk)) return(NA_real_)
    as.numeric(item$chunk$coordinates_in_chunk)
  }, numeric(1L)))
  geometry_valid <- .dsvert_dp_has_exact_names(
    geometry_summary, geometry_fields) &&
    identical(geometry_summary$granularity,
              "whole_intersecting_public_chunks_v1") &&
    identical(as.numeric(geometry_summary$block_coordinates),
              as.numeric(block$length)) &&
    identical(as.numeric(geometry_summary$included_public_coordinates),
              included_coordinates) &&
    identical(as.numeric(geometry_summary$overfetch_coordinates),
              included_coordinates - as.numeric(block$length)) &&
    identical(as.numeric(geometry_summary$intersecting_chunk_count),
              as.numeric(length(evidence))) &&
    is.finite(included_coordinates) &&
    included_coordinates >= as.numeric(block$length)
  if (!isTRUE(geometry_valid)) {
    stop("The Gaussian whole-chunk proof geometry is invalid",
         call. = FALSE)
  }
  selected <- character()
  for (offset in seq_along(evidence)) {
    item <- evidence[[offset]]
    index <- expected_chunks[[offset]]
    if (!.dsvert_dp_has_exact_names(
          item, c("version", "chunk_hash", "chunk", "merkle_proof")) ||
        !identical(item$version,
                   .DSVERT_DP_GAUSSIAN_CERTIFICATE_CHUNK_VERSION) ||
        !.dsvert_dp_capsule_source_hex(item$chunk_hash)) {
      stop("Invalid Gaussian public chunk evidence", call. = FALSE)
    }
    chunk <- item$chunk
    fields <- c(
      "version", "capsule_id", "release_contract_hash", "chunk_index",
      "chunk_count", "coordinate_offset", "coordinates_in_chunk",
      "output_lattice_bits", "output_lattice_scale", "scaled_values",
      "value_encoding", "postprocessing", "source_values_exposed",
      "preclamp_values_exposed")
    geometry <- .dsvert_vector_chunk_geometry(vector$contract, index)
    values <- if (is.list(chunk$scaled_values)) {
      unlist(chunk$scaled_values, use.names = FALSE)
    } else NULL
    valid <- .dsvert_dp_has_exact_names(chunk, fields) &&
      identical(chunk$version, .DSVERT_CLIENT_VECTOR_CHUNK_VERSION) &&
      identical(chunk$capsule_id, certificate$capsule_id) &&
      identical(chunk$release_contract_hash,
                certificate$release_contract_hash) &&
      identical(as.numeric(chunk$chunk_index), as.numeric(index)) &&
      identical(as.numeric(chunk$chunk_count),
                as.numeric(vector$contract$chunk_count)) &&
      identical(as.numeric(chunk$coordinate_offset),
                as.numeric(geometry$offset)) &&
      identical(as.numeric(chunk$coordinates_in_chunk),
                as.numeric(geometry$count)) &&
      identical(as.numeric(chunk$output_lattice_bits),
                as.numeric(certificate$output_lattice_bits)) &&
      identical(as.numeric(chunk$output_lattice_scale),
                as.numeric(certificate$output_lattice_scale)) &&
      is.character(values) && length(values) == geometry$count &&
      !anyNA(values) &&
      all(vapply(values, .dsvert_vector_integer_text, logical(1L))) &&
      identical(chunk$value_encoding,
                "nonnegative-decimal-integer-common-lattice-v1") &&
      identical(chunk$postprocessing,
                paste0("signed-Ring128-decode-then-fixed-public-coordinate-",
                       "clamp-v1")) &&
      identical(chunk$source_values_exposed, FALSE) &&
      identical(chunk$preclamp_values_exposed, FALSE) &&
      identical(.dsvert_vector_hash(chunk), item$chunk_hash) &&
      identical(item$chunk_hash, vector$hashes[[index + 1L]]) &&
      identical(.dsvert_joint_dp_client_canonical(item$merkle_proof),
                .dsvert_joint_dp_client_canonical(
                  .dsvert_vector_merkle_proof(vector$hashes, index)))
    if (!isTRUE(valid)) {
      stop("A Gaussian public chunk failed its signed Merkle proof",
           call. = FALSE)
    }
    global <- seq.int(geometry$offset + 1L,
                      geometry$offset + geometry$count)
    keep <- which(global >= block$start & global <= block$end)
    selected <- c(selected, values[keep])
  }
  if (length(selected) != block$length ||
      !identical(.dsvert_dp_gaussian_certificate_hash(list(
        protocol = "dsvert-dp-gaussian-released-block-values-v1",
        start = as.numeric(block$start), end = as.numeric(block$end),
        scaled_values = as.list(selected))),
        certificate$block_values_sha256)) {
    stop("The Gaussian released block commitment is invalid", call. = FALSE)
  }
  .dsvert_vector_scaled_to_double(
    selected, as.numeric(certificate$output_lattice_scale))
}

.dsvert_dp_gaussian_certificate_result_match <- function(
    object, certificate, artifact, moment) {
  if (identical(object, certificate)) return(invisible(TRUE))
  if (!is.list(object) ||
      !identical(object$analysis_id, certificate$analysis_id) ||
      !identical(
        .dsvert_joint_dp_client_json(object$signed_artifact),
        .dsvert_joint_dp_client_json(artifact)) ||
      !isTRUE(all.equal(as.numeric(object$n_obs), as.numeric(moment$n),
                        tolerance = 0)) ||
      !is.list(object$sufficient_statistics_dp) ||
      !isTRUE(all.equal(
        object$sufficient_statistics_dp$gram_projected, moment$gram,
        tolerance = 0)) ||
      !isTRUE(all.equal(
        object$sufficient_statistics_dp$cross_projected, moment$cross,
        tolerance = 0)) ||
      !isTRUE(all.equal(
        object$sufficient_statistics_dp$outcome_square_projected,
        moment$outcome_square, tolerance = 0)) ||
      !identical(object$capsule_id, certificate$capsule_id) ||
      !identical(object$manifest_sha256, certificate$manifest_sha256) ||
      !identical(object$final_vector_root,
                 certificate$final_vector_root) ||
      !identical(object$coordinate_order_sha256,
                 certificate$coordinate_order_sha256)) {
    stop("The Gaussian result does not match its signed certificate",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_gaussian_synopsis_certificate_build <- function(
    context, artifact, block, coordinates) {
  release <- context$release
  provenance <- if (is.list(release)) release$signed_provenance else NULL
  bundle <- context$manifest_bundle
  compilation <- context$verification_compilation
  required_roots <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  required_provenance <- c(
    "version", "ordered_peer_pinset", "peer_pinset_sha256",
    "designated_noise_peers", required_roots, "full_plan_sha256",
    "compile_receipts", "release_receipts", "replay_responses",
    "protected_shares_included", "preclamp_values_included",
    "source_values_included", "intermediate_payload_exposed",
    "durable_replay")
  same_owner_artifact <- (
    identical(artifact$version,
              "bounded-normalized-gaussian-sufficient-statistics-v1") &&
      identical(artifact$spec_version, "v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "random_intercept_v1")) ||
    ((identical(artifact$version,
                .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION) &&
        identical(artifact$spec_version, "random_intercept_fixed_v2")) ||
     (identical(artifact$version,
                .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION) &&
        identical(artifact$spec_version, "random_intercept_fixed_v3") &&
        identical(artifact$estimation_profile, "reml"))) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "gaussian_random_slope_grid_v1")) ||
    (artifact$version %in% c(
       .DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION,
       .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, if (identical(
         artifact$version, .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION)) {
           "gaussian_ar1_robust_working_gls_grid_v1"
         } else {
           "gaussian_ar1_working_gls_grid_v1"
         })) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_GLMM_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "binary_random_intercept_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_POISSON_GLMM_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "poisson_random_intercept_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_POISSON_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "poisson_random_slope_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "binary_random_slope_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS[["binomial"]]) &&
       identical(artifact$spec_version, "binomial_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS[["poisson"]]) &&
       identical(artifact$spec_version, "poisson_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_NB_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "negative_binomial_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "multinomial_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_ORDINAL_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "ordinal_grid_v1")) ||
    (identical(artifact$version,
               .DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION) &&
       identical(artifact$spec_version, "cox_partial_likelihood_grid_v1"))
  cross_owner_artifact <- identical(
    artifact$version, .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "v2") &&
    identical(artifact$implementation_state,
              "cross_owner_exact_gc_materialized") &&
    identical(artifact$cross_owner_state, "exact_gc_to_joint_dp_vector_v1")
  cross_evidence <- context$cross_gaussian_evidence[[artifact$analysis_id]] %||%
    NULL
  if (!isTRUE(context$synopsis) || !is.list(bundle) || !is.list(compilation) ||
      !is.list(provenance) || !all(required_provenance %in% names(provenance)) ||
      !identical(
        provenance$version,
        "dsvert-stateless-synopsis-public-provenance-v1") ||
      !identical(provenance$protected_shares_included, FALSE) ||
      !identical(provenance$preclamp_values_included, FALSE) ||
      !identical(provenance$source_values_included, FALSE) ||
      !identical(provenance$intermediate_payload_exposed, FALSE) ||
      !identical(provenance$durable_replay, TRUE) ||
      !(isTRUE(same_owner_artifact) || isTRUE(cross_owner_artifact)) ||
      (isTRUE(same_owner_artifact) &&
       (!identical(artifact$implementation_state, "same_owner_materialized") ||
        !identical(artifact$cross_owner_state,
                   "reserved_not_materialized"))) ||
      (isTRUE(cross_owner_artifact) &&
       (!is.list(cross_evidence) || !length(cross_evidence)))) {
    stop("The Gaussian result lacks closed Synopsis provenance",
         call. = FALSE)
  }
  if (!all(vapply(required_roots, function(field) {
        .dsvert_vector_hex(release[[field]]) &&
          identical(release[[field]], provenance[[field]])
      }, logical(1L)))) {
    stop("The Gaussian Synopsis roots are detached", call. = FALSE)
  }
  pinset <- .dsvert_dp_gaussian_pinset(provenance$ordered_peer_pinset)
  designated <- unlist(provenance$designated_noise_peers, use.names = FALSE)
  if (!identical(.dsvert_dp_pinset_hash(pinset),
                 provenance$peer_pinset_sha256) ||
      length(designated) != 2L || anyNA(designated) ||
      anyDuplicated(designated) || !all(designated %in% names(pinset)) ||
      !identical(
        .dsvert_joint_dp_client_json(compilation$receipts),
        .dsvert_joint_dp_client_json(provenance$compile_receipts))) {
    stop("The Gaussian Synopsis peer or compilation binding is invalid",
         call. = FALSE)
  }
  if (!is.character(release$scaled_values) ||
      length(release$scaled_values) != release$coordinate_count ||
      anyNA(release$scaled_values) ||
      !all(vapply(release$scaled_values, .dsvert_vector_integer_text,
                  logical(1L))) ||
      !is.numeric(coordinates) || length(coordinates) != block$length) {
    stop("The Gaussian Synopsis released coordinates are invalid",
         call. = FALSE)
  }
  positions <- seq.int(block$start, block$end)
  scaled <- release$scaled_values[positions]
  descriptor <- block$descriptor
  descriptor_sha256 <- .dsvert_dp_gaussian_certificate_hash(descriptor)
  block_public <- block[c(
    "family", "key", "start", "end", "length", "owner_peer", "dataset")]
  names(block_public)[names(block_public) == "key"] <- "analysis_id"
  block_values_sha256 <- .dsvert_dp_gaussian_certificate_hash(list(
    protocol = "dsvert-dp-gaussian-synopsis-block-values-v1",
    start = as.numeric(block$start), end = as.numeric(block$end),
    scaled_values = as.list(scaled)))
  public_bundle <- bundle
  public_bundle$context <- NULL
  evidence <- list(
    version = .DSVERT_DP_GAUSSIAN_SYNOPSIS_EVIDENCE_VERSION,
    federation_state_json =
      .dsvert_dp_gaussian_synopsis_evidence_json(
        unclass(context$status), "federation state"),
    manifest_bundle_json =
      .dsvert_dp_gaussian_synopsis_evidence_json(
        public_bundle, "manifest bundle"),
    compilation_json =
      .dsvert_dp_gaussian_synopsis_evidence_json(
        compilation, "compilation"),
    release_set_json =
      .dsvert_dp_gaussian_synopsis_evidence_json(
        provenance$release_receipts, "release set"),
    replay_set_json =
      .dsvert_dp_gaussian_synopsis_evidence_json(
        provenance$replay_responses, "replay set"))
  if (isTRUE(cross_owner_artifact)) {
    evidence$cross_gaussian_evidence_json <-
      .dsvert_dp_gaussian_synopsis_evidence_json(
        cross_evidence, "cross-owner Gaussian evidence")
  }
  unsigned <- list(
    version = .DSVERT_DP_GAUSSIAN_SYNOPSIS_CERTIFICATE_VERSION,
    analysis_id = artifact$analysis_id,
    dataset = artifact$dataset,
    owner_peer = artifact$owner_peer,
    cross_owner_state = artifact$cross_owner_state,
    descriptor_sha256 = descriptor_sha256,
    descriptor = descriptor,
    block = block_public,
    block_values_sha256 = block_values_sha256,
    manifest_sha256 = release$manifest_sha256,
    schema_sha256 = bundle$schema_sha256,
    workload_contract_sha256 = bundle$workload_contract_sha256,
    artifact_key = release$artifact_key,
    execution_id = release$execution_id,
    contract_sha256 = release$contract_sha256,
    attempt_sha256 = release$attempt_sha256,
    source_contract_sha256 = release$source_contract_sha256,
    result_set_sha256 = release$result_set_sha256,
    final_vector_root = release$final_vector_root,
    coordinate_order_sha256 = release$coordinate_order_sha256,
    coordinate_count = as.numeric(release$coordinate_count),
    output_lattice_bits = as.numeric(release$output_lattice_bits),
    output_lattice_scale = as.numeric(release$output_lattice_scale),
    plan_sha256 = release$plan_sha256,
    backend = release$backend,
    sampler = release$sampler,
    mechanism = release$mechanism,
    epsilon = release$epsilon,
    delta = release$delta,
    implementation_delta = release$implementation_delta,
    peer_context = list(
      ordered_peer_pinset = as.list(pinset),
      peer_pinset_sha256 = provenance$peer_pinset_sha256,
      designated_noise_peers = as.list(designated)),
    signed_evidence = evidence,
    public_dp_coordinates_only = TRUE,
    protected_shares_included = FALSE,
    preclamp_values_included = FALSE,
    patient_derived_identifiers_included = FALSE)
  certificate <- c(unsigned, list(
    certificate_sha256 = .dsvert_dp_gaussian_certificate_hash(unsigned)))
  .dsvert_dp_gaussian_synopsis_no_legacy_fields(certificate)
  .dsvert_dp_gaussian_cache_pinset(pinset)
  certificate
}

.dsvert_dp_gaussian_synopsis_trusted <- function(certificate) {
  evidence <- certificate$signed_evidence
  status <- .dsvert_dp_gaussian_synopsis_evidence_decode(
    evidence$federation_state_json, "federation state")
  if (!is.list(status) || is.null(names(status)) || length(status) < 2L ||
      anyNA(names(status)) || any(!nzchar(names(status))) ||
      anyDuplicated(names(status))) {
    stop("The Gaussian Synopsis federation evidence is invalid",
         call. = FALSE)
  }
  peers <- sort(names(status), method = "radix")
  if (!identical(names(status), peers)) {
    stop("The Gaussian Synopsis peer order is not canonical", call. = FALSE)
  }
  class(status) <- c("ds.vertDPSynopsisStatus", "list")
  policy <- tryCatch(status[[peers[[1L]]]]$policy,
                     error = function(error) NULL)
  policy_info <- .dsvert_dp_synopsis_client_policy_v1(policy, peers)
  placeholder <- stats::setNames(lapply(peers, function(peer) {
    structure(list(peer = peer), class = "dsvert_certificate_peer")
  }), peers)
  context <- list(
    status = status, servers = peers, pinset = policy_info$pins,
    designated = policy_info$designated,
    conns = placeholder[policy_info$designated], all_conns = placeholder,
    policy = policy_info$value)
  bundle <- .dsvert_dp_gaussian_synopsis_evidence_decode(
    evidence$manifest_bundle_json, "manifest bundle")
  if (!is.list(bundle) || "context" %in% names(bundle)) {
    stop("The Gaussian Synopsis manifest evidence is invalid",
         call. = FALSE)
  }
  manifest_preview <- .dsvert_joint_dp_client_decode(
    bundle$manifest_json, "Gaussian Synopsis manifest evidence",
    .DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES)
  bundle$logical_snapshot <- .dsvert_joint_dp_client_canonical(
    manifest_preview$logical_snapshot)
  bundle$artifact_commitments <- .dsvert_joint_dp_client_canonical(
    bundle$artifact_commitments)
  bundle$context <- context
  trusted <- .dsvert_dp_synopsis_client_bundle(bundle, status)
  list(status = status, bundle = bundle, trusted = trusted)
}

.dsvert_dp_gaussian_synopsis_result_match <- function(
    object, certificate, artifact, moment) {
  if (identical(object, certificate)) return(invisible(TRUE))
  roots <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  valid <- is.list(object) &&
    identical(object$analysis_id, certificate$analysis_id) &&
    identical(.dsvert_joint_dp_client_json(object$signed_artifact),
              .dsvert_joint_dp_client_json(artifact)) &&
    isTRUE(all.equal(as.numeric(object$n_obs), as.numeric(moment$n),
                     tolerance = 0)) &&
    is.list(object$sufficient_statistics_dp) &&
    isTRUE(all.equal(object$sufficient_statistics_dp$gram_projected,
                     moment$gram, tolerance = 0)) &&
    isTRUE(all.equal(object$sufficient_statistics_dp$cross_projected,
                     moment$cross, tolerance = 0)) &&
    isTRUE(all.equal(
      object$sufficient_statistics_dp$outcome_square_projected,
      moment$outcome_square, tolerance = 0)) &&
    all(vapply(roots, function(field) {
      identical(object[[field]], certificate[[field]])
    }, logical(1L))) &&
    identical(object$coordinate_order_sha256,
              certificate$coordinate_order_sha256) &&
    is.list(object$provenance_certificate) &&
    identical(object$provenance_certificate$certificate_sha256,
              certificate$certificate_sha256)
  if (!isTRUE(valid)) {
    stop("The Gaussian result does not match its signed Synopsis certificate",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_gaussian_synopsis_certificate_validate <- function(
    object, certificate, trusted_pinset = NULL) {
  required <- c(
    "version", "analysis_id", "dataset", "owner_peer",
    "cross_owner_state", "descriptor_sha256", "descriptor", "block",
    "block_values_sha256", "manifest_sha256", "schema_sha256",
    "workload_contract_sha256", "artifact_key", "execution_id",
    "contract_sha256", "attempt_sha256", "source_contract_sha256",
    "result_set_sha256", "final_vector_root", "coordinate_order_sha256",
    "coordinate_count", "output_lattice_bits", "output_lattice_scale",
    "plan_sha256", "backend", "sampler", "mechanism", "epsilon", "delta",
    "implementation_delta", "peer_context", "signed_evidence",
    "public_dp_coordinates_only", "protected_shares_included",
    "preclamp_values_included", "patient_derived_identifiers_included",
    "certificate_sha256")
  cross_owner <- identical(
    certificate$cross_owner_state, "exact_gc_to_joint_dp_vector_v1")
  if (!.dsvert_dp_has_exact_names(certificate, required) ||
      !identical(certificate$version,
                 .DSVERT_DP_GAUSSIAN_SYNOPSIS_CERTIFICATE_VERSION) ||
      !certificate$cross_owner_state %in% c(
        "reserved_not_materialized", "exact_gc_to_joint_dp_vector_v1") ||
      !identical(certificate$public_dp_coordinates_only, TRUE) ||
      !identical(certificate$protected_shares_included, FALSE) ||
      !identical(certificate$preclamp_values_included, FALSE) ||
      !identical(certificate$patient_derived_identifiers_included, FALSE) ||
      !identical(certificate$certificate_sha256,
        .dsvert_dp_gaussian_certificate_hash(certificate[
          setdiff(names(certificate), "certificate_sha256")]))) {
    stop("Invalid Gaussian Synopsis provenance certificate", call. = FALSE)
  }
  .dsvert_dp_gaussian_synopsis_no_legacy_fields(certificate)
  peer_fields <- c(
    "ordered_peer_pinset", "peer_pinset_sha256", "designated_noise_peers")
  evidence_fields <- c(
    "version", "federation_state_json", "manifest_bundle_json",
    "compilation_json", "release_set_json", "replay_set_json")
  if (isTRUE(cross_owner)) {
    evidence_fields <- c(evidence_fields, "cross_gaussian_evidence_json")
  }
  if (!.dsvert_dp_has_exact_names(certificate$peer_context, peer_fields) ||
      !.dsvert_dp_has_exact_names(
        certificate$signed_evidence, evidence_fields) ||
      !identical(certificate$signed_evidence$version,
                 .DSVERT_DP_GAUSSIAN_SYNOPSIS_EVIDENCE_VERSION)) {
    stop("Invalid Gaussian Synopsis evidence envelope", call. = FALSE)
  }
  pinset <- .dsvert_dp_gaussian_pinset(
    certificate$peer_context$ordered_peer_pinset)
  pin_hash <- .dsvert_dp_pinset_hash(pinset)
  designated <- unlist(
    certificate$peer_context$designated_noise_peers, use.names = FALSE)
  if (!identical(pin_hash, certificate$peer_context$peer_pinset_sha256) ||
      length(designated) != 2L || anyNA(designated) ||
      anyDuplicated(designated) || !all(designated %in% names(pinset))) {
    stop("Invalid Gaussian Synopsis peer context", call. = FALSE)
  }
  authenticity <- "unanchored"
  if (!is.null(trusted_pinset)) {
    trusted_pins <- .dsvert_dp_gaussian_pinset(
      trusted_pinset, "trusted pinset")
    if (!identical(trusted_pins, pinset)) {
      stop("The Gaussian certificate does not match the trusted pinset",
           call. = FALSE)
    }
    authenticity <- "caller_anchored"
  } else if (exists(pin_hash, envir = .dsvert_dp_gaussian_trust_cache,
                    inherits = FALSE) && identical(
                      get(pin_hash, envir =
                            .dsvert_dp_gaussian_trust_cache,
                          inherits = FALSE), pinset)) {
    authenticity <- "session_transport_anchored"
  }
  reconstructed <- .dsvert_dp_gaussian_synopsis_trusted(certificate)
  trusted <- reconstructed$trusted
  bundle <- reconstructed$bundle
  if (!identical(trusted$context$pinset, pinset) ||
      !identical(trusted$context$designated, designated) ||
      !identical(bundle$manifest_sha256, certificate$manifest_sha256) ||
      !identical(bundle$schema_sha256, certificate$schema_sha256) ||
      !identical(bundle$workload_contract_sha256,
                 certificate$workload_contract_sha256)) {
    stop("The Gaussian Synopsis manifest bindings changed", call. = FALSE)
  }
  compilation <- .dsvert_dp_gaussian_synopsis_evidence_decode(
    certificate$signed_evidence$compilation_json, "compilation")
  compiled <- .dsvert_dp_synopsis_client_compile(
    compilation, trusted, bundle)
  execution <- .dsvert_dp_synopsis_client_execution(compiled)
  release_values <- .dsvert_dp_gaussian_synopsis_evidence_decode(
    certificate$signed_evidence$release_set_json, "release set")
  if (!is.list(release_values)) {
    stop("The Gaussian Synopsis release evidence is invalid", call. = FALSE)
  }
  release_json <- lapply(release_values, .dsvert_joint_dp_client_json)
  releases <- .dsvert_dp_synopsis_client_release_set(
    release_json, compiled, execution, trusted)
  replay_values <- .dsvert_dp_gaussian_synopsis_evidence_decode(
    certificate$signed_evidence$replay_set_json, "replay set")
  chunk_count <- as.integer(releases$reference$public_chunk_count)
  if (!is.list(replay_values) || length(replay_values) != chunk_count) {
    stop("The Gaussian Synopsis replay evidence is incomplete",
         call. = FALSE)
  }
  replay_input <- stats::setNames(lapply(seq_len(chunk_count), function(index) {
    json <- .dsvert_joint_dp_client_json(replay_values[[index]])
    stats::setNames(rep(list(json), 2L), designated)
  }), as.character(seq.int(0L, chunk_count - 1L)))
  replay <- .dsvert_dp_synopsis_client_replay(
    replay_input, releases, compiled, execution, trusted)
  reference <- releases$reference
  roots <- c(
    artifact_key = compiled$artifact$artifact_key,
    execution_id = execution$execution_id,
    contract_sha256 = reference$contract_sha256,
    attempt_sha256 = reference$attempt_sha256,
    source_contract_sha256 = reference$source_contract_sha256,
    result_set_sha256 = reference$result_set_sha256,
    final_vector_root = reference$final_vector_root)
  if (!all(vapply(names(roots), function(field) {
        identical(certificate[[field]], unname(roots[[field]]))
      }, logical(1L))) ||
      !identical(certificate$coordinate_order_sha256,
                 compiled$layout$sha256) ||
      !identical(as.numeric(certificate$coordinate_count),
                 as.numeric(compiled$layout$coordinate_count)) ||
      !identical(as.numeric(certificate$output_lattice_bits),
                 as.numeric(compiled$lattice$output_lattice_bits)) ||
      !identical(as.numeric(certificate$output_lattice_scale),
                 as.numeric(compiled$lattice$output_lattice_scale)) ||
      !identical(certificate$plan_sha256,
                 compiled$physical$full_plan_sha256) ||
      !identical(certificate$backend, compiled$profile$backend) ||
      !identical(certificate$sampler, compiled$profile$sampler) ||
      !identical(certificate$mechanism, reference$mechanism) ||
      !identical(as.numeric(certificate$epsilon),
                 as.numeric(reference$epsilon)) ||
      !identical(as.numeric(certificate$delta),
                 as.numeric(reference$delta)) ||
      !identical(certificate$implementation_delta, paste0(
        reference$implementation_delta_numerator, "/",
        reference$implementation_delta_denominator))) {
    stop("The Gaussian Synopsis execution bindings changed", call. = FALSE)
  }
  count_block <- .dsvert_dp_capsule_single_block(
    compiled$layout, "admitted_count",
    description = "signed admitted-count capacity block")
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  cox_partial_grid <- identical(
    certificate$descriptor$version,
    .DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION)
  artifact <- if (isTRUE(cox_partial_grid)) {
    .dsvert_dp_cox_partial_grid_artifact(
      trusted$manifest, certificate$dataset, certificate$analysis_id,
      certificate$owner_peer, trusted$manifest$admission$adjacency,
      compiled$lattice$output_lattice_scale, capacity)
  } else {
    .dsvert_dp_gaussian_artifact(
      trusted$manifest, certificate$dataset, certificate$analysis_id,
      certificate$owner_peer, trusted$manifest$admission$adjacency,
      compiled$lattice$output_lattice_scale, capacity)
  }
  blocks <- .dsvert_dp_capsule_vector_blocks(
    compiled$layout,
    if (isTRUE(cox_partial_grid)) "survival_artifacts" else "gaussian_models",
    dataset = certificate$dataset,
    owner_peer = certificate$owner_peer)
  blocks <- blocks[vapply(blocks, function(block) {
    identical(block$key, certificate$analysis_id)
  }, logical(1L))]
  block <- if (length(blocks) == 1L) blocks[[1L]] else NULL
  block_public <- if (is.list(block)) block[c(
    "family", "key", "start", "end", "length", "owner_peer", "dataset")]
  else NULL
  if (is.list(block_public)) {
    names(block_public)[names(block_public) == "key"] <- "analysis_id"
  }
  if (!is.list(block) ||
      !identical(.dsvert_joint_dp_client_json(block_public),
                 .dsvert_joint_dp_client_json(certificate$block)) ||
      !identical(.dsvert_joint_dp_client_json(block$descriptor),
                 .dsvert_joint_dp_client_json(certificate$descriptor)) ||
      !identical(.dsvert_dp_gaussian_certificate_hash(
        certificate$descriptor), certificate$descriptor_sha256)) {
    stop("The Gaussian Synopsis descriptor or layout changed", call. = FALSE)
  }
  if (isTRUE(cross_owner)) {
    public <- .dsvert_dp_gaussian_synopsis_evidence_decode(
      certificate$signed_evidence$cross_gaussian_evidence_json,
      "cross-owner Gaussian evidence")
    receipts <- public
    fields <- c(
      "version", "phase", "analysis_id", "peer_name", "peer_identity_pk",
      "artifact_sha256", "source_contract_sha256", "private_layout_sha256",
      "transcript_sha256", "numeric_certificate_sha256",
      "exact_transcript_sha256", "coordinate_count", "public_start",
      "public_end", "public_coordinate_order_sha256", "ring_bits", "frac_bits",
      "state", "fixed_transcript", "private_result_exposed",
      "exact_intermediates_exposed", "alignment_hash_exposed", "signature")
    if (!is.list(receipts) || !identical(names(receipts), designated)) {
      stop("The Gaussian certificate lacks both cross-owner attestations",
           call. = FALSE)
    }
    receipts <- lapply(designated, function(peer) {
      value <- receipts[[peer]]
      valid <- .dsvert_dp_has_exact_names(value, fields) &&
        identical(value$version,
                  .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_PUBLIC_EVIDENCE_VERSION) &&
        identical(value$phase, "cross_gaussian_public_result_evidence") &&
        identical(value$analysis_id, certificate$analysis_id) &&
        identical(value$peer_name, peer) &&
        identical(value$peer_identity_pk, unname(pinset[[peer]])) &&
        identical(value$artifact_sha256,
                  .dsvert_dp_capsule_source_hash(artifact)) &&
        identical(value$source_contract_sha256,
                  certificate$source_contract_sha256) &&
        identical(value$private_layout_sha256,
                  .dsvert_dp_gaussian_cross_layout_client(
                    trusted$manifest)$transport_coordinate_order_sha256) &&
        identical(value$transcript_sha256,
                  .dsvert_dp_capsule_source_hash(artifact$transcript)) &&
        identical(value$numeric_certificate_sha256,
                  .dsvert_dp_capsule_source_hash(artifact$numeric_certificate)) &&
        .dsvert_dp_capsule_source_hex(value$exact_transcript_sha256) &&
        identical(as.numeric(value$coordinate_count), as.numeric(block$length)) &&
        identical(as.numeric(value$public_start), as.numeric(block$start)) &&
        identical(as.numeric(value$public_end), as.numeric(block$end)) &&
        identical(value$public_coordinate_order_sha256, compiled$layout$sha256) &&
        identical(as.numeric(value$ring_bits), 128) &&
        identical(as.numeric(value$frac_bits),
                  as.numeric(artifact$numeric_grid_bits)) &&
        identical(value$state, "complete") &&
        identical(value$fixed_transcript, TRUE) &&
        identical(value$private_result_exposed, FALSE) &&
        identical(value$exact_intermediates_exposed, FALSE) &&
        identical(value$alignment_hash_exposed, FALSE) &&
        .dsvert_dp_capsule_source_verify(
          value, "cross-gaussian-synopsis-evidence", peer,
          list(pinset = pinset))
      if (!isTRUE(valid)) {
        stop("A Gaussian cross-owner attestation is invalid", call. = FALSE)
      }
      value
    })
    stable <- setdiff(fields, c("peer_name", "peer_identity_pk", "signature"))
    if (length(unique(vapply(receipts, function(value) {
          .dsvert_joint_dp_client_json(value[stable])
        }, character(1L)))) != 1L) {
      stop("The Gaussian cross-owner attestations disagree", call. = FALSE)
    }
  }
  scaled <- replay$scaled[seq.int(block$start, block$end)]
  block_hash <- .dsvert_dp_gaussian_certificate_hash(list(
    protocol = "dsvert-dp-gaussian-synopsis-block-values-v1",
    start = as.numeric(block$start), end = as.numeric(block$end),
    scaled_values = as.list(scaled)))
  if (!identical(block_hash, certificate$block_values_sha256)) {
    stop("The Gaussian Synopsis released block commitment changed",
         call. = FALSE)
  }
  coordinates <- .dsvert_vector_scaled_to_double(
    scaled, compiled$lattice$output_lattice_scale)
  coordinate_upper <- if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION)) {
    .dsvert_dp_lmm_fixed_coordinate_upper(artifact)
  } else if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_GLMM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_POISSON_GLMM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_POISSON_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else if (identical(artifact$version,
                       .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else if (artifact$version %in% c(
               .DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION,
               .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else if (artifact$version %in%
             unname(.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS)) {
    as.numeric(artifact$statistic_maximum)
  } else if (identical(artifact$version,
                       .DSVERT_CLIENT_DP_NB_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else if (identical(artifact$version,
                       .DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else if (identical(artifact$version,
                       .DSVERT_CLIENT_DP_ORDINAL_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else if (identical(artifact$version,
                       .DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else if (identical(
        artifact$version,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION)) {
    c(capacity, capacity,
      capacity * artifact$max_patients_per_cluster,
      rep(capacity, 3L))
  } else {
    rep(capacity, artifact$coordinate_count)
  }
  if (length(coordinates) != artifact$coordinate_count ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 0) || any(coordinates > coordinate_upper)) {
    stop("The certified Gaussian Synopsis block violates its signed bounds",
         call. = FALSE)
  }
  if (identical(artifact$version,
                .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION)) {
    moment <- .dsvert_dp_lmm_moments(coordinates, artifact)
    accuracy_release <- list(
      manifest_sha256 = certificate$manifest_sha256,
      epsilon = certificate$epsilon,
      mechanism = certificate$mechanism,
      implementation_delta = certificate$implementation_delta,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256,
      backend = certificate$backend)
    accuracy_simultaneous_95 <- .dsvert_dp_vector_accuracy_radius(
      accuracy_release, trusted$manifest,
      coordinate_count = artifact$coordinate_count,
      confidence = 0.95, maximum_error = max(coordinate_upper))
    return(list(
      integrity_valid = TRUE, authenticity = authenticity,
      artifact = artifact,
      bounds = list(outcome = artifact$outcome, cluster = artifact$cluster),
      design_terms = NULL, coordinates = coordinates,
      validated_moment = moment,
      coordinate_capacity = capacity,
      output_lattice_scale = compiled$lattice$output_lattice_scale,
      accuracy_simultaneous_95 = accuracy_simultaneous_95,
      sufficient_statistics_dp = list(
        projected_summary = moment$projected_summary),
      n_obs = moment$projected_summary[["n"]],
      cohort_id = trusted$status[[trusted$context$servers[[1L]]]]$policy$cohort_id,
      logical_snapshot = trusted$manifest$logical_snapshot,
      analysis_id = certificate$analysis_id,
      manifest_sha256 = certificate$manifest_sha256,
      artifact_key = certificate$artifact_key,
      execution_id = certificate$execution_id,
      contract_sha256 = certificate$contract_sha256,
      attempt_sha256 = certificate$attempt_sha256,
      source_contract_sha256 = certificate$source_contract_sha256,
      result_set_sha256 = certificate$result_set_sha256,
      final_vector_root = certificate$final_vector_root,
      coordinate_order_sha256 = certificate$coordinate_order_sha256,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256,
      backend = certificate$backend, sampler = certificate$sampler,
      mechanism = certificate$mechanism,
      epsilon = certificate$epsilon, delta = certificate$delta,
      implementation_delta = certificate$implementation_delta,
      privacy = list(
        version = "dsvert-per-synopsis-dp-v1",
        per_artifact_epsilon = certificate$epsilon,
        per_artifact_delta = certificate$delta,
        unlimited_replay = TRUE,
        replay_is_postprocessing = TRUE,
        finite_global_composition_claim = FALSE),
      certificate = certificate))
  }
  if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION)) {
    moment <- .dsvert_dp_lmm_fixed_moments(coordinates, artifact)
    accuracy_release <- list(
      manifest_sha256 = certificate$manifest_sha256,
      epsilon = certificate$epsilon,
      mechanism = certificate$mechanism,
      implementation_delta = certificate$implementation_delta,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256,
      backend = certificate$backend)
    accuracy_simultaneous_95 <- .dsvert_dp_vector_accuracy_radius(
      accuracy_release, trusted$manifest,
      coordinate_count = artifact$coordinate_count,
      confidence = 0.95, maximum_error = max(coordinate_upper))
    return(list(
      integrity_valid = TRUE, authenticity = authenticity,
      artifact = artifact,
      bounds = list(outcome = artifact$outcome, cluster = artifact$cluster,
                    predictors = artifact$predictors),
      design_terms = artifact$design_terms, coordinates = coordinates,
      validated_moment = moment,
      coordinate_capacity = capacity,
      output_lattice_scale = compiled$lattice$output_lattice_scale,
      accuracy_simultaneous_95 = accuracy_simultaneous_95,
      sufficient_statistics_dp = list(
        projected_summary = moment$projected_summary),
      n_obs = moment$n_obs %||% NULL,
      cohort_id = trusted$status[[trusted$context$servers[[1L]]]]$policy$cohort_id,
      logical_snapshot = trusted$manifest$logical_snapshot,
      analysis_id = certificate$analysis_id,
      manifest_sha256 = certificate$manifest_sha256,
      artifact_key = certificate$artifact_key,
      execution_id = certificate$execution_id,
      contract_sha256 = certificate$contract_sha256,
      attempt_sha256 = certificate$attempt_sha256,
      source_contract_sha256 = certificate$source_contract_sha256,
      result_set_sha256 = certificate$result_set_sha256,
      final_vector_root = certificate$final_vector_root,
      coordinate_order_sha256 = certificate$coordinate_order_sha256,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256,
      backend = certificate$backend, sampler = certificate$sampler,
      mechanism = certificate$mechanism,
      epsilon = certificate$epsilon, delta = certificate$delta,
      implementation_delta = certificate$implementation_delta,
      privacy = list(
        version = "dsvert-per-synopsis-dp-v1",
        per_artifact_epsilon = certificate$epsilon,
        per_artifact_delta = certificate$delta,
        unlimited_replay = TRUE,
        replay_is_postprocessing = TRUE,
        finite_global_composition_claim = FALSE),
      certificate = certificate))
  }
  if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_GLMM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_POISSON_GLMM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_POISSON_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION)) {
    moment <- if (identical(artifact$version,
                            .DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION)) {
      .dsvert_dp_cox_partial_grid_moment(coordinates, artifact)
    } else if (artifact$version %in% c(
                 .DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION,
                 .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION)) {
      .dsvert_dp_gee_ar1_grid_moment(coordinates, artifact)
    } else if (identical(artifact$version,
                          .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
      .dsvert_dp_lmm_random_slope_grid_moment(coordinates, artifact)
    } else if (identical(artifact$version,
                          .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
      .dsvert_dp_glmm_random_slope_grid_moment(coordinates, artifact)
    } else if (identical(artifact$version,
                          .DSVERT_CLIENT_DP_POISSON_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION)) {
      .dsvert_dp_poisson_glmm_random_slope_grid_moment(coordinates, artifact)
    } else .dsvert_dp_glmm_grid_moment(coordinates, artifact)
    accuracy_release <- list(
      manifest_sha256 = certificate$manifest_sha256,
      epsilon = certificate$epsilon,
      mechanism = certificate$mechanism,
      implementation_delta = certificate$implementation_delta,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256,
      backend = certificate$backend)
    accuracy_simultaneous_95 <- .dsvert_dp_vector_accuracy_radius(
      accuracy_release, trusted$manifest,
      coordinate_count = artifact$coordinate_count,
      confidence = 0.95, maximum_error = max(coordinate_upper))
    return(list(
      integrity_valid = TRUE, authenticity = authenticity,
      artifact = artifact,
      bounds = if (identical(artifact$version,
                             .DSVERT_CLIENT_DP_COX_PARTIAL_GRID_ARTIFACT_VERSION)) {
        list(time = artifact$time, event = artifact$event,
             predictors = artifact$predictors)
      } else {
        list(outcome = artifact$outcome, cluster = artifact$cluster,
             predictors = artifact$predictors)
      },
      design_terms = artifact$design_terms, coordinates = coordinates,
      validated_moment = moment,
      coordinate_capacity = capacity,
      output_lattice_scale = compiled$lattice$output_lattice_scale,
      accuracy_simultaneous_95 = accuracy_simultaneous_95,
      sufficient_statistics_dp = list(
        candidate_negative_log_likelihoods = coordinates),
      n_obs = NULL,
      cohort_id = trusted$status[[trusted$context$servers[[1L]]]]$policy$cohort_id,
      logical_snapshot = trusted$manifest$logical_snapshot,
      analysis_id = certificate$analysis_id,
      manifest_sha256 = certificate$manifest_sha256,
      artifact_key = certificate$artifact_key,
      execution_id = certificate$execution_id,
      contract_sha256 = certificate$contract_sha256,
      attempt_sha256 = certificate$attempt_sha256,
      source_contract_sha256 = certificate$source_contract_sha256,
      result_set_sha256 = certificate$result_set_sha256,
      final_vector_root = certificate$final_vector_root,
      coordinate_order_sha256 = certificate$coordinate_order_sha256,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256,
      backend = certificate$backend, sampler = certificate$sampler,
      mechanism = certificate$mechanism,
      epsilon = certificate$epsilon, delta = certificate$delta,
      implementation_delta = certificate$implementation_delta,
      privacy = list(
        version = "dsvert-per-synopsis-dp-v1",
        per_artifact_epsilon = certificate$epsilon,
        per_artifact_delta = certificate$delta,
        unlimited_replay = TRUE,
        replay_is_postprocessing = TRUE,
        finite_global_composition_claim = FALSE),
      certificate = certificate))
  }
  if (artifact$version %in%
      unname(.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS)) {
    family_name <- names(.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS)[[
      match(artifact$version,
            unname(.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS))]]
    moment <- .dsvert_dp_glm_grid_moment(coordinates, artifact, family_name)
    accuracy_release <- list(
      manifest_sha256 = certificate$manifest_sha256,
      epsilon = certificate$epsilon, mechanism = certificate$mechanism,
      implementation_delta = certificate$implementation_delta,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256, backend = certificate$backend)
    accuracy_simultaneous_95 <- .dsvert_dp_vector_accuracy_radius(
      accuracy_release, trusted$manifest,
      coordinate_count = artifact$coordinate_count, confidence = 0.95,
      maximum_error = max(coordinate_upper))
    return(list(
      integrity_valid = TRUE, authenticity = authenticity, artifact = artifact,
      bounds = list(outcome = artifact$outcome, predictors = artifact$predictors),
      design_terms = artifact$design_terms, coordinates = coordinates,
      validated_moment = moment, coordinate_capacity = capacity,
      output_lattice_scale = compiled$lattice$output_lattice_scale,
      accuracy_simultaneous_95 = accuracy_simultaneous_95,
      sufficient_statistics_dp = list(
        candidate_negative_log_likelihoods = coordinates), n_obs = NULL,
      cohort_id = trusted$status[[trusted$context$servers[[1L]]]]$policy$cohort_id,
      logical_snapshot = trusted$manifest$logical_snapshot,
      analysis_id = certificate$analysis_id,
      manifest_sha256 = certificate$manifest_sha256,
      artifact_key = certificate$artifact_key,
      execution_id = certificate$execution_id,
      contract_sha256 = certificate$contract_sha256,
      attempt_sha256 = certificate$attempt_sha256,
      source_contract_sha256 = certificate$source_contract_sha256,
      result_set_sha256 = certificate$result_set_sha256,
      final_vector_root = certificate$final_vector_root,
      coordinate_order_sha256 = certificate$coordinate_order_sha256,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256, backend = certificate$backend,
      sampler = certificate$sampler, mechanism = certificate$mechanism,
      epsilon = certificate$epsilon, delta = certificate$delta,
      implementation_delta = certificate$implementation_delta,
      privacy = list(version = "dsvert-per-synopsis-dp-v1",
        per_artifact_epsilon = certificate$epsilon,
        per_artifact_delta = certificate$delta, unlimited_replay = TRUE,
        replay_is_postprocessing = TRUE,
        finite_global_composition_claim = FALSE), certificate = certificate))
  }
  if (identical(artifact$version,
                .DSVERT_CLIENT_DP_NB_GRID_ARTIFACT_VERSION)) {
    moment <- .dsvert_dp_nb_grid_moment(coordinates, artifact)
    accuracy_release <- list(
      manifest_sha256 = certificate$manifest_sha256,
      epsilon = certificate$epsilon, mechanism = certificate$mechanism,
      implementation_delta = certificate$implementation_delta,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256, backend = certificate$backend)
    accuracy_simultaneous_95 <- .dsvert_dp_vector_accuracy_radius(
      accuracy_release, trusted$manifest,
      coordinate_count = artifact$coordinate_count, confidence = 0.95,
      maximum_error = max(coordinate_upper))
    return(list(
      integrity_valid = TRUE, authenticity = authenticity, artifact = artifact,
      bounds = list(outcome = artifact$outcome, predictors = artifact$predictors),
      design_terms = artifact$design_terms, coordinates = coordinates,
      validated_moment = moment, coordinate_capacity = capacity,
      output_lattice_scale = compiled$lattice$output_lattice_scale,
      accuracy_simultaneous_95 = accuracy_simultaneous_95,
      sufficient_statistics_dp = list(
        candidate_negative_log_likelihoods = coordinates), n_obs = NULL,
      cohort_id = trusted$status[[trusted$context$servers[[1L]]]]$policy$cohort_id,
      logical_snapshot = trusted$manifest$logical_snapshot,
      analysis_id = certificate$analysis_id,
      manifest_sha256 = certificate$manifest_sha256,
      artifact_key = certificate$artifact_key,
      execution_id = certificate$execution_id,
      contract_sha256 = certificate$contract_sha256,
      attempt_sha256 = certificate$attempt_sha256,
      source_contract_sha256 = certificate$source_contract_sha256,
      result_set_sha256 = certificate$result_set_sha256,
      final_vector_root = certificate$final_vector_root,
      coordinate_order_sha256 = certificate$coordinate_order_sha256,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256, backend = certificate$backend,
      sampler = certificate$sampler, mechanism = certificate$mechanism,
      epsilon = certificate$epsilon, delta = certificate$delta,
      implementation_delta = certificate$implementation_delta,
      privacy = list(version = "dsvert-per-synopsis-dp-v1",
        per_artifact_epsilon = certificate$epsilon,
        per_artifact_delta = certificate$delta, unlimited_replay = TRUE,
        replay_is_postprocessing = TRUE,
        finite_global_composition_claim = FALSE), certificate = certificate))
  }
  if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_ORDINAL_GRID_ARTIFACT_VERSION)) {
    moment <- if (identical(artifact$version,
                            .DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION)) {
      .dsvert_dp_multinom_grid_moment(coordinates, artifact)
    } else {
      .dsvert_dp_ordinal_grid_moment(coordinates, artifact)
    }
    accuracy_release <- list(
      manifest_sha256 = certificate$manifest_sha256,
      epsilon = certificate$epsilon, mechanism = certificate$mechanism,
      implementation_delta = certificate$implementation_delta,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256, backend = certificate$backend)
    accuracy_simultaneous_95 <- .dsvert_dp_vector_accuracy_radius(
      accuracy_release, trusted$manifest,
      coordinate_count = artifact$coordinate_count, confidence = 0.95,
      maximum_error = max(coordinate_upper))
    return(list(
      integrity_valid = TRUE, authenticity = authenticity, artifact = artifact,
      bounds = list(outcome = artifact$outcome, predictors = artifact$predictors),
      design_terms = artifact$design_terms, coordinates = coordinates,
      validated_moment = moment, coordinate_capacity = capacity,
      output_lattice_scale = compiled$lattice$output_lattice_scale,
      accuracy_simultaneous_95 = accuracy_simultaneous_95,
      sufficient_statistics_dp = list(
        candidate_negative_log_likelihoods = coordinates), n_obs = NULL,
      cohort_id = trusted$status[[trusted$context$servers[[1L]]]]$policy$cohort_id,
      logical_snapshot = trusted$manifest$logical_snapshot,
      analysis_id = certificate$analysis_id,
      manifest_sha256 = certificate$manifest_sha256,
      artifact_key = certificate$artifact_key,
      execution_id = certificate$execution_id,
      contract_sha256 = certificate$contract_sha256,
      attempt_sha256 = certificate$attempt_sha256,
      source_contract_sha256 = certificate$source_contract_sha256,
      result_set_sha256 = certificate$result_set_sha256,
      final_vector_root = certificate$final_vector_root,
      coordinate_order_sha256 = certificate$coordinate_order_sha256,
      mechanism_plan = compiled$physical$full_plan,
      plan_sha256 = certificate$plan_sha256, backend = certificate$backend,
      sampler = certificate$sampler, mechanism = certificate$mechanism,
      epsilon = certificate$epsilon, delta = certificate$delta,
      implementation_delta = certificate$implementation_delta,
      privacy = list(version = "dsvert-per-synopsis-dp-v1",
        per_artifact_epsilon = certificate$epsilon,
        per_artifact_delta = certificate$delta, unlimited_replay = TRUE,
        replay_is_postprocessing = TRUE,
        finite_global_composition_claim = FALSE), certificate = certificate))
  }
  moment <- .dsvert_dp_gaussian_unpack(coordinates, artifact, capacity)
  .dsvert_dp_gaussian_synopsis_result_match(
    object, certificate, artifact, moment)
  accuracy_release <- list(
    manifest_sha256 = certificate$manifest_sha256,
    epsilon = certificate$epsilon,
    mechanism = certificate$mechanism,
    implementation_delta = certificate$implementation_delta,
    mechanism_plan = compiled$physical$full_plan,
    plan_sha256 = certificate$plan_sha256,
    backend = certificate$backend)
  accuracy_simultaneous_95 <- .dsvert_dp_vector_accuracy_radius(
    accuracy_release, trusted$manifest,
    coordinate_count = artifact$coordinate_count,
    confidence = 0.95, maximum_error = capacity)
  list(
    integrity_valid = TRUE, authenticity = authenticity,
    artifact = artifact,
    bounds = list(outcome = artifact$outcome,
                  predictors = artifact$predictors),
    design_terms = artifact$design_terms,
    coordinates = coordinates,
    validated_moment = moment,
    coordinate_capacity = capacity,
    output_lattice_scale = compiled$lattice$output_lattice_scale,
    accuracy_simultaneous_95 = accuracy_simultaneous_95,
    sufficient_statistics_dp = list(
      gram_projected = moment$gram, cross_projected = moment$cross,
      outcome_square_projected = moment$outcome_square),
    n_obs = moment$n,
    cohort_id = trusted$status[[trusted$context$servers[[1L]]]]$policy$cohort_id,
    logical_snapshot = trusted$manifest$logical_snapshot,
    analysis_id = certificate$analysis_id,
    manifest_sha256 = certificate$manifest_sha256,
    artifact_key = certificate$artifact_key,
    execution_id = certificate$execution_id,
    contract_sha256 = certificate$contract_sha256,
    attempt_sha256 = certificate$attempt_sha256,
    source_contract_sha256 = certificate$source_contract_sha256,
    result_set_sha256 = certificate$result_set_sha256,
    final_vector_root = certificate$final_vector_root,
    coordinate_order_sha256 = certificate$coordinate_order_sha256,
    mechanism_plan = compiled$physical$full_plan,
    plan_sha256 = certificate$plan_sha256,
    backend = certificate$backend, sampler = certificate$sampler,
    mechanism = certificate$mechanism,
    epsilon = certificate$epsilon, delta = certificate$delta,
    implementation_delta = certificate$implementation_delta,
    privacy = list(
      version = "dsvert-per-synopsis-dp-v1",
      per_artifact_epsilon = certificate$epsilon,
      per_artifact_delta = certificate$delta,
      unlimited_replay = TRUE,
      replay_is_postprocessing = TRUE,
      finite_global_composition_claim = FALSE),
    certificate = certificate)
}

#' Verify a bounded Gaussian provenance certificate without DSI
#'
#' Dispatches between the same-owner no-lifetime Synopsis provenance
#' certificate v1 and the byte-compatible legacy capsule certificate v3. It
#' revalidates the signed schema/manifest, pinned-peer receipts, Gaussian
#' descriptor and exact released DP block before algebra. A self-contained
#' certificate proves internal integrity. Strong
#' peer authenticity additionally requires a caller-supplied trusted pinset or
#' the ephemeral pinset cache populated by an online fit in this R session.
#'
#' @param x A `ds.vertDPGaussian` result or its `provenance_certificate`.
#' @param trusted_pinset Optional trusted named peer-to-Ed25519-public-key map.
#' @return A verified public provenance view with separate `integrity_valid`
#'   and `authenticity` fields.
#' @export
ds.validateDPGaussianCertificate <- function(x, trusted_pinset = NULL) {
  certificate <- if (is.list(x) && (
      identical(x$version, .DSVERT_DP_GAUSSIAN_CERTIFICATE_VERSION) ||
      identical(x$version,
                .DSVERT_DP_GAUSSIAN_SYNOPSIS_CERTIFICATE_VERSION))) {
    x
  } else if (is.list(x)) {
    x$provenance_certificate
  } else NULL
  if (is.list(certificate) && identical(
      certificate$version,
      .DSVERT_DP_GAUSSIAN_SYNOPSIS_CERTIFICATE_VERSION)) {
    return(.dsvert_dp_gaussian_synopsis_certificate_validate(
      x, certificate, trusted_pinset = trusted_pinset))
  }
  required <- c(
    "version", "analysis_id", "dataset", "owner_peer",
    "cross_owner_state", "query_contract_sha256", "artifact_sha256",
    "artifact", "artifact_commitment", "block", "block_values_sha256",
    "chunk_evidence", "chunk_proof_geometry", "capsule_id",
    "release_instance_id", "release_instance",
    "manifest_sha256", "schema_sha256",
    "workload_contract_sha256", "release_contract_hash",
    "source_contract_hash", "transcript_hash", "final_vector_root",
    "coordinate_order_sha256", "coordinate_count", "output_lattice_bits",
    "output_lattice_scale", "capsule_mechanism",
    "mechanism_plan_sha256", "backend", "sampler", "mechanism",
    "epsilon", "delta",
    "implementation_delta", "peer_context", "signed_evidence",
    "public_dp_coordinates_only", "protected_shares_included",
    "preclamp_values_included", "patient_derived_identifiers_included",
    "certificate_sha256")
  if (!.dsvert_dp_has_exact_names(certificate, required) ||
      !identical(certificate$version,
                 .DSVERT_DP_GAUSSIAN_CERTIFICATE_VERSION) ||
      !identical(certificate$public_dp_coordinates_only, TRUE) ||
      !identical(certificate$protected_shares_included, FALSE) ||
      !identical(certificate$preclamp_values_included, FALSE) ||
      !identical(certificate$patient_derived_identifiers_included, FALSE) ||
      !identical(certificate$certificate_sha256,
                 .dsvert_dp_gaussian_certificate_hash(certificate[
                   setdiff(names(certificate), "certificate_sha256")]))) {
    stop("Invalid Gaussian provenance certificate", call. = FALSE)
  }
  peer_fields <- c(
    "ordered_peer_pinset", "peer_pinset_sha256",
    "designated_noise_peers", "privacy_epochs", "noise_key_ids",
    "noise_provider_ids", "release_domain_generations",
    "release_domain_ids")
  if (!.dsvert_dp_has_exact_names(certificate$peer_context, peer_fields)) {
    stop("Invalid Gaussian certificate peer context", call. = FALSE)
  }
  pinset <- .dsvert_dp_gaussian_pinset(
    certificate$peer_context$ordered_peer_pinset)
  pin_hash <- .dsvert_dp_pinset_hash(pinset)
  if (!identical(pin_hash, certificate$peer_context$peer_pinset_sha256)) {
    stop("The Gaussian certificate pinset hash is invalid", call. = FALSE)
  }
  instance <- certificate$release_instance
  instance_fields <- c("version", "capsule_id", "peer_noise_roots")
  designated <- unlist(
    certificate$peer_context$designated_noise_peers, use.names = FALSE)
  roots <- if (is.list(instance)) instance$peer_noise_roots else NULL
  expected_roots <- stats::setNames(lapply(designated, function(peer) {
    list(
      privacy_epoch = as.numeric(
        certificate$peer_context$privacy_epochs[[peer]]),
      noise_key_id = certificate$peer_context$noise_key_ids[[peer]],
      provider_id = certificate$peer_context$noise_provider_ids[[peer]],
      release_domain_generation = as.numeric(
        certificate$peer_context$release_domain_generations[[peer]]),
      release_domain_id =
        certificate$peer_context$release_domain_ids[[peer]])
  }), designated)
  if (!.dsvert_dp_has_exact_names(instance, instance_fields) ||
      !identical(instance$version,
                 .DSVERT_CLIENT_VECTOR_RELEASE_INSTANCE_VERSION) ||
      !identical(instance$capsule_id, certificate$capsule_id) ||
      !identical(names(roots), designated) ||
      anyDuplicated(vapply(roots, `[[`, character(1L), "noise_key_id")) ||
      !identical(.dsvert_joint_dp_client_canonical(roots),
                 .dsvert_joint_dp_client_canonical(expected_roots)) ||
      !identical(.dsvert_vector_hash(instance),
                 certificate$release_instance_id)) {
    stop("The Gaussian certificate release instance is invalid",
         call. = FALSE)
  }
  authenticity <- "unanchored"
  if (!is.null(trusted_pinset)) {
    trusted <- .dsvert_dp_gaussian_pinset(trusted_pinset, "trusted pinset")
    if (!identical(trusted, pinset)) {
      stop("The Gaussian certificate does not match the trusted pinset",
           call. = FALSE)
    }
    authenticity <- "caller_anchored"
  } else if (exists(pin_hash, envir = .dsvert_dp_gaussian_trust_cache,
                    inherits = FALSE) && identical(
                      get(pin_hash, envir =
                            .dsvert_dp_gaussian_trust_cache,
                          inherits = FALSE), pinset)) {
    authenticity <- "session_transport_anchored"
  }
  commitment <- certificate$artifact_commitment
  root <- .dsvert_dp_gaussian_artifact_root_from_proof(commitment)
  index_context <- commitment$context
  context_fields <- c(
    "version", "manifest_sha256", "capsule_id", "domain", "cohort_id",
    "peer_pinset_sha256", "designated_noise_peers",
    "privacy_epoch_scope",
    "epsilon", "delta", "adjacency", "unit_capacity",
    "max_records_per_unit", "overflow_policy", "consortium_id",
    "policy_contract_hash", "logical_snapshot", "capsule_schema",
    "admission_sha256", "bounds_sha256", "workload_sha256",
    "release_lattice", "capsule_mechanism", "mechanism_selection",
    "coordinate_layout_version",
    "coordinate_count", "coordinate_order_sha256")
  if (!.dsvert_dp_has_exact_names(index_context, context_fields) ||
      !identical(index_context$version,
                 .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_VERSION) ||
      !identical(index_context$manifest_sha256,
                 certificate$manifest_sha256) ||
      !identical(index_context$capsule_id, certificate$capsule_id) ||
      !identical(index_context$peer_pinset_sha256, pin_hash) ||
      !identical(index_context$privacy_epoch_scope,
                 "per_peer_signed_receipts_v1") ||
      !identical(as.numeric(index_context$coordinate_count),
                 as.numeric(certificate$coordinate_count)) ||
      !identical(index_context$coordinate_order_sha256,
                 certificate$coordinate_order_sha256) ||
      !identical(root, certificate$signed_evidence$
                   manifest_build_receipts[[1L]]$
                   artifact_commitments_root)) {
    stop("The Gaussian artifact commitment context is invalid",
         call. = FALSE)
  }
  evidence_fields <- c(
    "version", "signed_schema_json", "manifest_build_receipts",
    "vector_prepare_receipts", "vector_release_receipts",
    "vector_finalization_receipts")
  if (!.dsvert_dp_has_exact_names(
        certificate$signed_evidence, evidence_fields) ||
      !identical(certificate$signed_evidence$version,
                 .DSVERT_DP_GAUSSIAN_CERTIFICATE_EVIDENCE_VERSION)) {
    stop("Invalid Gaussian signed evidence envelope", call. = FALSE)
  }
  .dsvert_dp_gaussian_certificate_schema(certificate, pinset, index_context)
  build_receipts <- .dsvert_dp_gaussian_certificate_build_receipts(
    certificate, pinset, c(index_context, list(
      artifact_commitments_root = root)))
  vector <- .dsvert_dp_gaussian_certificate_vector(
    certificate, pinset, index_context)
  entry <- commitment$entry
  block <- certificate$block
  block_expected <- entry[c(
    "family", "analysis_id", "start", "end", "length", "owner_peer",
    "dataset")]
  if (!.dsvert_dp_has_exact_names(block, names(block_expected)) ||
      !identical(.dsvert_joint_dp_client_canonical(block),
                 .dsvert_joint_dp_client_canonical(block_expected)) ||
      !identical(entry$analysis_id, certificate$analysis_id) ||
      !identical(entry$dataset, certificate$dataset) ||
      !identical(entry$owner_peer, certificate$owner_peer) ||
      !identical(entry$descriptor_sha256, certificate$artifact_sha256) ||
      !identical(.dsvert_dp_gaussian_certificate_hash(
        certificate$artifact), certificate$artifact_sha256)) {
    stop("The Gaussian artifact descriptor or coordinate range changed",
         call. = FALSE)
  }
  coordinates <- .dsvert_dp_gaussian_certificate_coordinates(
    certificate, vector)
  scale <- as.numeric(certificate$output_lattice_scale)
  capacity <- as.numeric(index_context$unit_capacity)
  pseudo_manifest <- list(workload = list(families = list(
    gaussian_models = list(artifacts = stats::setNames(
      list(certificate$artifact), certificate$analysis_id)))))
  artifact <- .dsvert_dp_gaussian_artifact(
    pseudo_manifest, certificate$dataset, certificate$analysis_id,
    certificate$owner_peer, index_context$adjacency, scale, capacity)
  coordinate_upper <- if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_FIXED_REML_ARTIFACT_VERSION)) {
    .dsvert_dp_lmm_fixed_coordinate_upper(artifact)
  } else if (identical(artifact$version,
                       .DSVERT_CLIENT_DP_RANDOM_INTERCEPT_ARTIFACT_VERSION)) {
    c(capacity, capacity, capacity * artifact$max_patients_per_cluster,
      rep(capacity, 3L))
  } else if (artifact$version %in% c(
        .DSVERT_CLIENT_DP_GLMM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_POISSON_GLMM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_POISSON_GLMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_LMM_RANDOM_SLOPE_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GEE_AR1_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_GEE_AR1_ROBUST_GRID_ARTIFACT_VERSION,
        unname(.DSVERT_CLIENT_DP_GLM_GRID_ARTIFACT_VERSIONS),
        .DSVERT_CLIENT_DP_NB_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_MULTINOM_GRID_ARTIFACT_VERSION,
        .DSVERT_CLIENT_DP_ORDINAL_GRID_ARTIFACT_VERSION)) {
    as.numeric(artifact$statistic_maximum)
  } else rep(capacity, artifact$coordinate_count)
  if (length(coordinates) != artifact$coordinate_count ||
      any(coordinates < 0) || any(coordinates > coordinate_upper)) {
    stop("The certified Gaussian block violates its signed bounds",
         call. = FALSE)
  }
  expected_query <- .dsvert_dp_gaussian_certificate_hash(list(
    protocol = "dsvert-dp-gaussian-query-contract-v1",
    capsule_id = certificate$capsule_id,
    release_instance_id = certificate$release_instance_id,
    release_instance = certificate$release_instance,
    analysis_id = certificate$analysis_id,
    artifact_sha256 = certificate$artifact_sha256,
    artifact_node_sha256 = entry$node_sha256,
    release_contract_hash = certificate$release_contract_hash,
    capsule_mechanism = certificate$capsule_mechanism,
    mechanism_plan_sha256 = certificate$mechanism_plan_sha256,
    backend = certificate$backend, sampler = certificate$sampler,
    release_mechanism = certificate$mechanism,
    final_vector_root = certificate$final_vector_root,
    coordinate_order_sha256 = certificate$coordinate_order_sha256,
    block_values_sha256 = certificate$block_values_sha256))
  if (!identical(expected_query, certificate$query_contract_sha256)) {
    stop("The Gaussian query contract hash is invalid", call. = FALSE)
  }
  moment <- .dsvert_dp_gaussian_unpack(coordinates, artifact, capacity)
  .dsvert_dp_gaussian_certificate_result_match(
    x, certificate, artifact, moment)
  privacy_epochs <- vapply(build_receipts, function(receipt) {
    as.numeric(receipt$privacy_epoch)
  }, numeric(1L))
  noise_key_ids <- vapply(build_receipts, `[[`, character(1L),
                          "noise_key_id")
  certificate_epochs <- tryCatch(vapply(
    certificate$peer_context$privacy_epochs, as.numeric, numeric(1L)),
    error = function(error) numeric())
  certificate_key_ids <- tryCatch(vapply(
    certificate$peer_context$noise_key_ids, as.character, character(1L)),
    error = function(error) character())
  if (!identical(privacy_epochs, certificate_epochs) ||
      !identical(noise_key_ids, certificate_key_ids)) {
    stop("The Gaussian certificate privacy epochs or key ids changed",
         call. = FALSE)
  }
  result <- list(
    integrity_valid = TRUE,
    authenticity = authenticity,
    artifact = artifact,
    bounds = list(outcome = artifact$outcome,
                  predictors = artifact$predictors),
    design_terms = artifact$design_terms,
    sufficient_statistics_dp = list(
      gram_projected = moment$gram,
      cross_projected = moment$cross,
      outcome_square_projected = moment$outcome_square),
    n_obs = moment$n,
    cohort_id = index_context$cohort_id,
    logical_snapshot = index_context$logical_snapshot,
    analysis_id = certificate$analysis_id,
    capsule_id = certificate$capsule_id,
    release_instance_id = certificate$release_instance_id,
    release_instance = certificate$release_instance,
    manifest_sha256 = certificate$manifest_sha256,
    query_contract_sha256 = certificate$query_contract_sha256,
    release_contract_hash = certificate$release_contract_hash,
    final_vector_root = certificate$final_vector_root,
    coordinate_order_sha256 = certificate$coordinate_order_sha256,
    privacy_epochs = privacy_epochs,
    noise_key_ids = noise_key_ids,
    noise_provider_ids = unlist(
      certificate$peer_context$noise_provider_ids, use.names = TRUE),
    release_domain_generations = unlist(
      certificate$peer_context$release_domain_generations,
      use.names = TRUE),
    release_domain_ids = unlist(
      certificate$peer_context$release_domain_ids, use.names = TRUE),
    capsule_mechanism = certificate$capsule_mechanism,
    backend = certificate$backend,
    sampler = certificate$sampler,
    mechanism = certificate$mechanism,
    mechanism_plan = vector$contract$mechanism_plan,
    mechanism_plan_sha256 = certificate$mechanism_plan_sha256,
    epsilon = certificate$epsilon,
    delta = certificate$delta,
    proof_geometry = c(certificate$chunk_proof_geometry, list(
      certificate_bytes = as.numeric(nchar(
        .dsvert_joint_dp_client_json(certificate), type = "bytes")),
      coordinate_overfetch_ratio = if (
          certificate$chunk_proof_geometry$block_coordinates > 0) {
        certificate$chunk_proof_geometry$overfetch_coordinates /
          certificate$chunk_proof_geometry$block_coordinates
      } else 0)),
    certificate = certificate)
  class(result) <- c("ds.vertDPGaussianCertificateValidation", "list")
  result
}
