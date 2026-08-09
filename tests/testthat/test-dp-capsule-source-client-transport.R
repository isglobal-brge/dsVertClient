.capsule_source_client_b64url <- function(value) {
  sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))),
    perl = TRUE)
}

.capsule_source_client_fixture <- function(
    k = 3L, coordinate_count = 8193L, source_count = k,
    ciphertext_bytes = 96L, window_transport = FALSE) {
  stopifnot(k %in% 2:5, source_count >= 1L, source_count <= k)
  peers <- paste0("peer_", letters[seq_len(k)])
  designated <- peers[1:2]
  sources <- peers[seq_len(source_count)]
  keys <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  pins <- vapply(keys, function(key) {
    .capsule_source_client_b64url(tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  pin_hash <- .dsvert_dp_pinset_hash(pins)

  lifetime_max <- 8
  capsules <- 7
  releases <- 0
  lifetime_epsilon <- .dsvert_joint_dp_lifetime_exact_total(1, lifetime_max)
  lifetime_delta <- .dsvert_joint_dp_lifetime_exact_total(
    2^-100, lifetime_max)
  telemetry <- list(
    capsules_created = capsules,
    lifetime_max_distinct_capsules = lifetime_max,
    remaining_distinct_capsules = lifetime_max - capsules,
    capsule_epsilon = 1,
    capsule_delta = 2^-100,
    cumulative_epsilon_upper_bound = capsules,
    cumulative_delta_upper_bound = capsules * 2^-100,
    cumulative_delta_vacuous = FALSE,
    composition_role = "basic_composition_authenticated_lifetime_bound",
    registration_policy =
      "allocator_admitted_distinct_capsules_up_to_lifetime_limit",
    privacy_budget_gate = TRUE, operation_limit = TRUE,
    request_limit = FALSE, history_can_deny_new_release = TRUE,
    admission_role = "allocator_reservation_before_protected_access")
  release_telemetry <- list(
    releases_published = releases,
    lifetime_max_distinct_capsules = lifetime_max,
    remaining_distinct_capsules = lifetime_max - capsules,
    release_epsilon = 1,
    release_delta = 2^-100,
    cumulative_epsilon_upper_bound = 0,
    cumulative_delta_upper_bound = 0,
    cumulative_delta_vacuous = FALSE,
    composition_role = "basic_composition_authenticated_lifetime_bound",
    release_accounting = "one_public_release_instance_per_capsule_id",
    replay_accounting = "same_instance_replay_not_recounted",
    rotation_accounting = paste0(
      "prepublication_rotation_only_postpublication_replay_or_",
      "fail_closed"),
    privacy_budget_gate = TRUE, operation_limit = TRUE,
    request_limit = FALSE, history_can_deny_operation = TRUE,
    admission_role = "authenticated_lifetime_gate_before_sampler")
  statuses <- stats::setNames(lapply(peers, function(peer) {
    is_designated <- peer %in% designated
    list(
      version = "dsvert-joint-dp-capsule-status-v5", enabled = TRUE,
      privacy_contract = list(
        definition = "bounded_lifetime_epsilon_delta_dp",
        scope = paste0(
          "at_most_N_immutable_snapshot_workload_capsules_per_stable_",
          "privacy_accountant_namespace"),
        adversary_model = "authenticated_semi_honest_noncollusion",
        assumptions = paste0(
          "declared_adjacency_bounds_immutable_snapshot_protocol_",
          "compliant_peers_at_least_one_noncolluding_designated_noise_peer_",
          "retains_and_uses_complete_authenticated_monotonic_history_",
          "stable_unique_privacy_accountant_namespace_per_protected_",
          "privacy_universe"),
        simultaneous_designated_history_rollback_protection =
          "not_claimed_without_external_linearizable_cas",
        transcript_security = "computational_mpc_and_csprng",
        malicious_security = FALSE,
        operation_accounting =
          "one_per_distinct_capsule_allocator_commit",
        privacy_budget_gate = TRUE,
        operation_limit = TRUE, request_limit = FALSE,
        history_can_deny_operation = TRUE,
        release_instance_accounting =
          "one_public_release_instance_per_capsule_id",
        accuracy_depends_on_request_history = FALSE,
        reuse = "unlimited_sticky_postprocessing",
        new_capsules = "allowed_until_authenticated_lifetime_bound",
        lifetime_max_distinct_capsules = lifetime_max,
        lifetime_epsilon_upper_bound = lifetime_epsilon,
        lifetime_delta_upper_bound = lifetime_delta),
      policy = list(
        contract = "immutable_reusable_capsule_v1",
        domain = "source-client-study", cohort_id = "source-client-cohort",
        peer_name = peer, own_identity_pk = unname(pins[[peer]]),
        peer_pinset = pins, peer_pinset_sha256 = pin_hash,
        peer_count = as.integer(k), designated_noise_peers = designated,
        capsule_epsilon = 1, capsule_delta = 2^-100,
        lifetime_max_distinct_capsules = lifetime_max,
        lifetime_epsilon_upper_bound = lifetime_epsilon,
        lifetime_delta_upper_bound = lifetime_delta,
        adjacency = "add_remove_patient", patient_column = "patient_id",
        unit_capacity = 1000000, max_records_per_unit = 2,
        overflow_policy = "reject_snapshot",
        sampler = "two_private_hmac_seeds_one_gc_sample_v1"),
      noise_root = list(
        protocol = "dsvert-dp-noise-root-v1",
        provider_id = paste0("provider-", peer),
        key_id = paste0("key-", peer), privacy_epoch = 1,
        external = FALSE, storage = "owner_only_file",
        automatic_generation = TRUE, automatic_recovery = TRUE,
        automatic_rotation = FALSE, rotation_count = 0,
        key_material_exposed = FALSE),
      release_domain = if (is_designated) list(
        version = "dsvert-joint-dp-release-domain-v1",
        generation = 1,
        domain_id = paste0("rd_", digest::digest(
          paste0("source-client/", peer),
          "sha256", serialize = FALSE)),
        rotation_count = 0,
        automatic_generation = TRUE, automatic_rotation = TRUE,
        snapshot_derived = FALSE, key_material_exposed = FALSE) else NULL,
      role = list(
        designated_noise_peer = is_designated,
        allocator = if (is_designated) "authenticated_ready" else
          "not_applicable_policy_attestor"),
      composition_telemetry = if (is_designated) telemetry else NULL,
      release_instance_telemetry = if (is_designated) {
        release_telemetry
      } else NULL)
  }), peers)

  logical_snapshot <- list(
    logical_snapshot_id = "fixed-source-client-snapshot", version = "v1",
    alignment_protocol_version = 1)
  admission <- list(protocol = "fixed-public-admission-v1")
  bounds <- list(protocol = "fixed-public-bounds-v1")
  scope_selection <- .dsvert_joint_dp_client_canonical(list(
    mode = "catalog_v1",
    explicit_catalog = list(
      numeric_moments = list(), categorical_marginals = list(),
      categorical_pairs = list(), correlations = list()),
    referenced_by_signed_specs = list(
      numeric = list(), categorical = list(), describe = list(),
      survival = list(), gaussian = NULL, vertical_cross = NULL),
    included = list(
      numeric_moments = list(), categorical_marginals = list(),
      same_owner_categorical_pairs = list(),
      same_owner_correlations = list())))
  primitive_scope <- list(
    version = "dsvert-biomedical-capsule-primitive-scope-v1",
    mode = "catalog_v1",
    authority = "custodian_policy_and_signed_workload_specs_only",
    analyst_expandable = FALSE,
    client_query_can_add_coordinates = FALSE,
    consensus = paste0(
      "byte_identical_manifest_hash_with_all_pinned_peer_build_",
      "signatures_required_before_source_access"),
    mismatch_behavior = "reject_before_protected_snapshot_resolution",
    compatibility_default = "all_schema",
    recommended_deployment_mode = "catalog_v1",
    selection_sha256 = .dsvert_dp_capsule_source_hash(scope_selection),
    selection = scope_selection,
    projected_cost = list(
      schema_numeric_column_count = 0,
      schema_categorical_column_count = 0,
      possible_same_owner_numeric_pair_count = 0,
      possible_same_owner_categorical_pair_count = 0,
      included_numeric_moment_count = 0,
      included_categorical_marginal_count = 0,
      included_numeric_pair_count = 0,
      included_categorical_pair_count = 0,
      included_cross_categorical_pair_count = 0,
      numeric_moment_coordinate_count = 0,
      numeric_pair_coordinate_count = 0,
      categorical_marginal_coordinate_count = 0,
      categorical_pair_coordinate_count = 0,
      gaussian_model_coordinate_count = 0,
      projected_coordinate_count = as.numeric(coordinate_count),
      projected_integer_l1_sensitivity = 1,
      projected_integer_l2_sensitivity = 1,
      automatic_pair_expansion = "none",
      scaling_contract = paste0(
        "linear_in_declared_univariates_plus_explicit_pairs_and_",
        "declared_model_cross_products")))
  workload <- list(
    workload_version = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION,
    coordinate_count = as.numeric(coordinate_count),
    capsule_mechanism = list(source_context_hash = strrep("3", 64L)),
    primitive_scope = primitive_scope,
    registered_release_lifecycle =
      .dsvert_dp_capsule_client_registered_release_lifecycle(),
    declared_workload_fully_materialized = TRUE,
    package_family_coverage_complete = FALSE,
    execution_state = .DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE)
  identity_context <- list(
    status = statuses, servers = peers, pinset = pins,
    designated = designated)
  policy_identity <- .dsvert_dp_capsule_source_policy_identity(
    identity_context)
  identity_contract <- list(
    protocol = "dsvert-joint-dp-capsule-identity-v3",
    consortium_id = policy_identity$consortium_id,
    policy_contract_hash = policy_identity$policy_contract_hash,
    peer_pinset_sha256 = pin_hash,
    logical_snapshot = logical_snapshot,
    capsule_schema = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION,
    admission = admission, bounds = bounds, workload = workload,
    privacy_epoch_scope = "per_peer_signed_receipts_v1")
  capsule_id <- .dsvert_dp_capsule_source_hash(identity_contract)
  manifest <- list(
    version = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION,
    logical_snapshot = logical_snapshot,
    capsule_schema = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION,
    admission = admission, bounds = bounds, workload = workload,
    capsule_identity = list(
      capsule_id = capsule_id, contract = identity_contract),
    execution_state = .DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE)
  manifest_json <- .dsvert_joint_dp_client_json(manifest)
  contract_hash <- strrep("4", 64L)
  chunk_count <- ceiling(
    coordinate_count / .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES)

  sign_value <- function(unsigned, peer, domain, omit_ciphertext = FALSE) {
    signed <- unsigned
    if (isTRUE(omit_ciphertext)) signed$ciphertext <- NULL
    message <- charToRaw(paste0(
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SIGNATURE_DOMAIN, domain, "|",
      .dsvert_joint_dp_client_json(signed)))
    c(unsigned, list(signature = .capsule_source_client_b64url(
      openssl::ed25519_sign(message, keys[[peer]]))))
  }
  sign_json <- function(unsigned, peer, domain, omit_ciphertext = FALSE) {
    .dsvert_joint_dp_client_json(sign_value(
      unsigned, peer, domain, omit_ciphertext))
  }
  decode <- function(value) .dsvert_joint_dp_client_decode(
    value, "test source object",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES)

  transport_public <- stats::setNames(lapply(seq_along(designated),
    function(index) .capsule_source_client_b64url(as.raw(
      (seq_len(32L) + 71L * index) %% 256L))), designated)
  ticket_values <- ticket_json <- vector("list", length(designated))
  names(ticket_values) <- names(ticket_json) <- designated
  for (recipient in designated) {
    unsigned <- list(
      version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_VERSION,
      phase = "recipient_key_committed",
      purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE,
      capsule_id = capsule_id, contract_hash = contract_hash,
      recipient_name = recipient,
      recipient_identity_pk = unname(pins[[recipient]]),
      transport_key_id = .dsvert_dp_capsule_source_hash(list(
        protocol = "dsvert-biomedical-capsule-source-key-id-v1",
        capsule_id = capsule_id, recipient_name = recipient,
        transport_pk = transport_public[[recipient]])),
      transport_pk = transport_public[[recipient]],
      peer_pinset_sha256 = pin_hash,
      designated_noise_peers = as.list(designated),
      source_peers = as.list(sources),
      coordinate_count = as.numeric(coordinate_count),
      chunk_coordinates = as.numeric(
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES),
      chunk_count = as.numeric(chunk_count), persistent = TRUE,
      ready_for_sampling = FALSE)
    ticket_values[[recipient]] <- sign_value(
      unsigned, recipient, "recipient-ticket")
    ticket_json[[recipient]] <- .dsvert_joint_dp_client_json(
      ticket_values[[recipient]])
  }
  window_capability <- .dsvert_dp_capsule_source_window_capability()
  ticket_negotiation_values <- ticket_negotiation_json <-
    vector("list", length(designated))
  names(ticket_negotiation_values) <- names(ticket_negotiation_json) <-
    designated
  for (recipient in designated) {
    unsigned <- list(
      version =
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_NEGOTIATION_VERSION,
      phase = "recipient_transport_window_attested",
      ticket_json = ticket_json[[recipient]],
      ticket_sha256 = .dsvert_dp_capsule_source_hash(
        ticket_values[[recipient]]),
      capability = window_capability)
    ticket_negotiation_values[[recipient]] <- sign_value(
      unsigned, recipient, "recipient-window-capability")
    ticket_negotiation_json[[recipient]] <- .dsvert_joint_dp_client_json(
      ticket_negotiation_values[[recipient]])
  }
  ticket_hashes <- vapply(ticket_values,
                          .dsvert_dp_capsule_source_hash, character(1L))

  transfer_ids <- stats::setNames(vapply(sources, function(source) {
    .dsvert_dp_capsule_source_transfer_id(
      contract_hash, capsule_id, source)
  }, character(1L)), sources)
  summary_values <- summary_json <- vector("list", length(sources))
  names(summary_values) <- names(summary_json) <- sources
  for (source in sources) {
    unsigned <- list(
      version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SUMMARY_VERSION,
      phase = "source_chunk_stream_ready",
      purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE,
      capsule_id = capsule_id, contract_hash = contract_hash,
      source_transfer_id = transfer_ids[[source]], source_name = source,
      source_identity_pk = unname(pins[[source]]),
      recipients = as.list(designated),
      coordinate_count = as.numeric(coordinate_count),
      chunk_coordinates = as.numeric(
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES),
      chunk_count = as.numeric(chunk_count), ring_bits = 128,
      record_encoding = "little_endian_unsigned_fixed_16_bytes",
      emitted_chunk_durable_replay = TRUE,
      unmaterialized_requires_same_snapshot = TRUE,
      complete_durable_replay = FALSE, history_gate = FALSE,
      ready_for_sampling = FALSE)
    summary_values[[source]] <- sign_value(
      unsigned, source, "source-summary")
    summary_json[[source]] <- .dsvert_joint_dp_client_json(
      summary_values[[source]])
  }
  negotiation_set_sha256 <- .dsvert_dp_capsule_source_hash(list(
    protocol = "dsvert-biomedical-capsule-source-negotiation-set-v1",
    negotiation_hashes = unname(lapply(
      ticket_negotiation_values,
      .dsvert_dp_capsule_source_hash))))
  summary_negotiation_values <- summary_negotiation_json <-
    vector("list", length(sources))
  names(summary_negotiation_values) <- names(summary_negotiation_json) <-
    sources
  for (source in sources) {
    unsigned <- list(
      version =
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SUMMARY_NEGOTIATION_VERSION,
      phase = "source_transport_window_attested",
      summary_json = summary_json[[source]],
      summary_sha256 = .dsvert_dp_capsule_source_hash(
        summary_values[[source]]),
      ticket_negotiation_set_sha256 = negotiation_set_sha256,
      capability = window_capability)
    summary_negotiation_values[[source]] <- sign_value(
      unsigned, source, "source-window-capability")
    summary_negotiation_json[[source]] <- .dsvert_joint_dp_client_json(
      summary_negotiation_values[[source]])
  }
  source_capability_values <- source_capability_json <-
    vector("list", length(sources))
  names(source_capability_values) <- names(source_capability_json) <- sources
  for (source in sources) {
    unsigned <- list(
      version =
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ATTESTATION_VERSION,
      phase = "source_transport_capability_attested",
      capsule_id = capsule_id, contract_hash = contract_hash,
      source_name = source,
      source_identity_pk = unname(pins[[source]]),
      capability = window_capability)
    source_capability_values[[source]] <- sign_value(
      unsigned, source, "source-window-capability-advertisement")
    source_capability_json[[source]] <- .dsvert_joint_dp_client_json(
      source_capability_values[[source]])
  }

  state <- new.env(parent = emptyenv())
  state$calls <- list()
  state$fault <- NULL
  state$ciphertext_epoch <- 0L
  state$accepted <- list()
  state$bundle_sizes <- numeric()
  state$window_transport <- isTRUE(window_transport)
  state$negotiation_error <- FALSE
  state$negotiation_oversize <- FALSE
  state$source_capability_unsupported <- character()

  make_envelope <- function(source, recipient, chunk_index) {
    offset <- chunk_index *
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES
    count <- min(
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES,
      coordinate_count - offset)
    source_index <- match(source, peers)
    recipient_index <- match(recipient, peers)
    raw <- as.raw((seq_len(ciphertext_bytes) + 19L * source_index +
      37L * recipient_index + 53L * chunk_index +
      97L * state$ciphertext_epoch) %% 256L)
    unsigned <- list(
      version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_VERSION,
      phase = "encrypted_source_chunk_committed",
      purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE,
      capsule_id = capsule_id, contract_hash = contract_hash,
      source_transfer_id = transfer_ids[[source]], source_name = source,
      source_identity_pk = unname(pins[[source]]),
      recipient_name = recipient,
      recipient_identity_pk = unname(pins[[recipient]]),
      recipient_ticket_hash = ticket_hashes[[recipient]],
      chunk_index = as.numeric(chunk_index),
      chunk_count = as.numeric(chunk_count),
      coordinate_offset = as.numeric(offset),
      coordinates_in_chunk = as.numeric(count),
      chunk_coordinates = as.numeric(
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES),
      ring_bits = 128,
      record_encoding = "little_endian_unsigned_fixed_16_bytes",
      ciphertext_bytes = as.numeric(length(raw)),
      ciphertext_sha256 = digest::digest(
        raw, algo = "sha256", serialize = FALSE),
      ciphertext = .capsule_source_client_b64url(raw),
      ready_for_sampling = FALSE)
    sign_value(unsigned, source, "encrypted-chunk", omit_ciphertext = TRUE)
  }
  make_bundle <- function(source, chunk_index) {
    envelopes <- lapply(designated, function(recipient) {
      make_envelope(source, recipient, chunk_index)
    })
    first <- envelopes[[1L]]
    bundle <- list(
      version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_BUNDLE_VERSION,
      phase = "encrypted_source_chunk_bundle_committed",
      purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE,
      capsule_id = capsule_id, contract_hash = contract_hash,
      source_transfer_id = transfer_ids[[source]], source_name = source,
      source_identity_pk = unname(pins[[source]]),
      recipients = as.list(designated), chunk_index = first$chunk_index,
      chunk_count = first$chunk_count,
      coordinate_offset = first$coordinate_offset,
      coordinates_in_chunk = first$coordinates_in_chunk,
      chunk_coordinates = first$chunk_coordinates, ring_bits = 128,
      record_encoding = "little_endian_unsigned_fixed_16_bytes",
      envelopes = unname(envelopes), ready_for_sampling = FALSE)
    encoded <- .dsvert_joint_dp_client_json(bundle)
    state$bundle_sizes <- c(state$bundle_sizes,
                            nchar(encoded, type = "bytes"))
    encoded
  }
  make_bundle_window <- function(source, indices) {
    if (length(indices) == 1L) return(make_bundle(source, indices[[1L]]))
    bundles <- list()
    bundle_bytes <- numeric()
    empty <- .dsvert_joint_dp_client_json(list(
      version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION,
      phase = "encrypted_source_chunk_window_committed",
      bundles = list(), ready_for_sampling = FALSE))
    for (index in indices) {
      bundle_json <- make_bundle(source, index)
      bundle <- decode(bundle_json)
      candidate_bytes <- nchar(empty, type = "bytes") +
        sum(bundle_bytes) + nchar(bundle_json, type = "bytes") +
        length(bundle_bytes)
      if (candidate_bytes >
          .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES) {
        if (!length(bundles)) return(bundle_json)
        break
      }
      bundles <- c(bundles, list(bundle))
      bundle_bytes <- c(bundle_bytes, nchar(bundle_json, type = "bytes"))
    }
    .dsvert_joint_dp_client_json(list(
      version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION,
      phase = "encrypted_source_chunk_window_committed",
      bundles = bundles, ready_for_sampling = FALSE))
  }
  make_ack <- function(envelope, recipient) {
    source_index <- match(envelope$source_name, sources)
    source_complete <- envelope$chunk_index == chunk_count - 1L
    aggregation_complete <- source_complete && source_index == length(sources)
    unsigned <- list(
      version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ACK_VERSION,
      phase = "source_chunk_aggregated",
      purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE,
      capsule_id = capsule_id, contract_hash = contract_hash,
      source_transfer_id = envelope$source_transfer_id,
      source_name = envelope$source_name,
      source_identity_pk = envelope$source_identity_pk,
      recipient_name = recipient,
      recipient_identity_pk = unname(pins[[recipient]]),
      recipient_ticket_hash = ticket_hashes[[recipient]],
      chunk_index = envelope$chunk_index,
      chunk_count = envelope$chunk_count,
      ciphertext_sha256 = envelope$ciphertext_sha256,
      source_complete = source_complete,
      capsule_aggregation_complete = aggregation_complete,
      history_gate = FALSE, ready_for_sampling = FALSE)
    sign_json(unsigned, recipient, "aggregate-ack")
  }

  aggregate <- function(conns, expr, async, error, errors.print) {
    wire_expressions <- if (is.list(expr) && !is.call(expr)) expr else
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    expressions <- lapply(wire_expressions, function(value) {
      method <- as.character(value[[1L]])
      for (argument in .DSVERT_DSI_TEXT_REMOTE_FORMALS[[method]]) {
        index <- which(names(value) == argument)
        if (length(index) == 1L) {
          value[[index]] <- .dsvert_dsi_text_decode(value[[index]])
        }
      }
      value
    })
    endpoints <- vapply(expressions, function(value) {
      as.character(value[[1L]])
    }, character(1L))
    stopifnot(length(unique(endpoints)) == 1L)
    endpoint <- endpoints[[1L]]
    state$calls[[length(state$calls) + 1L]] <- list(
      endpoint = endpoint, targets = names(conns), expressions = expressions,
      wire_expressions = wire_expressions, async = async)
    if (identical(endpoint, "dsvertDPCapsuleSourceTicketDS") &&
        identical(async, FALSE) && isTRUE(state$negotiation_error)) {
      for (target in names(conns)) error(target, "unused argument")
      return(stats::setNames(vector("list", length(conns)), names(conns)))
    }
    if (identical(endpoint, "dsvertDPCapsuleSourceTicketDS") &&
        identical(async, FALSE) && isTRUE(state$negotiation_oversize)) {
      return(stats::setNames(lapply(names(conns), function(target) {
        strrep("A", .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_TICKET_BYTES + 1L)
      }), names(conns)))
    }
    result <- switch(endpoint,
      dsvertJointDPCapsuleStatusDS = statuses[names(conns)],
      dsvertDPCapsuleSourceTicketDS = stats::setNames(lapply(
        names(conns), function(target) {
          arguments <- as.list(expressions[[target]])[-1L]
          capability_only <- identical(
            arguments$transport_contract,
            .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ONLY_REQUEST_V2)
          if (isTRUE(capability_only)) {
            if (target %in% state$source_capability_unsupported) {
              error(target, "unused argument")
              return(NULL)
            }
            return(source_capability_json[[target]])
          }
          requested <- identical(
            arguments$transport_contract,
            .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)
          if (isTRUE(state$window_transport) && requested) {
            ticket_negotiation_json[[target]]
          } else {
            ticket_json[[target]]
          }
        }), names(conns)),
      dsvertDPCapsuleSourcePrepareDS = stats::setNames(lapply(
        names(conns), function(target) {
          arguments <- as.list(expressions[[target]])[-1L]
          first <- decode(arguments$first_ticket_json)
          if (identical(
              first$version,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_NEGOTIATION_VERSION)) {
            summary_negotiation_json[[target]]
          } else {
            summary_json[[target]]
          }
        }), names(conns)),
      dsvertDPCapsuleSourceChunkDS = {
        source <- names(conns)[[1L]]
        arguments <- as.list(expressions[[source]])[-1L]
        stopifnot(identical(
          arguments$source_transfer_id, transfer_ids[[source]]))
        stats::setNames(list(make_bundle_window(
          source, as.numeric(arguments$chunk_index))), source)
      },
      dsvertDPCapsuleSourceAcceptDS = {
        responses <- stats::setNames(vector("list", length(conns)),
                                     names(conns))
        for (recipient in names(conns)) {
          arguments <- as.list(expressions[[recipient]])[-1L]
          payload <- decode(arguments$envelope_json)
          envelopes <- if (identical(
              payload$version,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION)) {
            payload$envelopes
          } else {
            list(payload)
          }
          acknowledgements <- list()
          conflict <- FALSE
          for (envelope in envelopes) {
            envelope_hash <- .dsvert_dp_capsule_source_hash(envelope)
            key <- paste(recipient, envelope$source_name,
                         envelope$chunk_index, sep = "|")
            prior <- state$accepted[[key]]
            if (!is.null(prior) && !identical(prior, envelope_hash)) {
              error(recipient, "conflicting retry")
              responses[[recipient]] <- NULL
              conflict <- TRUE
              break
            }
            state$accepted[[key]] <- envelope_hash
            acknowledgements[[length(acknowledgements) + 1L]] <-
              decode(make_ack(envelope, recipient))
          }
          if (isTRUE(conflict)) next
          responses[[recipient]] <- if (length(envelopes) == 1L) {
            .dsvert_joint_dp_client_json(acknowledgements[[1L]])
          } else {
            .dsvert_joint_dp_client_json(list(
              version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ACK_WINDOW_VERSION,
              phase = "source_chunk_window_aggregated",
              acknowledgements = acknowledgements,
              ready_for_sampling = FALSE))
          }
        }
        responses
      },
      stop("unexpected endpoint: ", endpoint))
    if (is.function(state$fault)) {
      result <- state$fault(
        endpoint, result, names(conns), expressions, length(state$calls))
    }
    result
  }
  state$sign_value <- sign_value
  state$sign_json <- sign_json
  state$decode <- decode
  state$make_bundle <- make_bundle
  list(
    peers = peers, designated = designated, sources = sources,
    keys = keys, pins = pins, status = statuses,
    manifest = manifest, manifest_json = manifest_json,
    capsule_id = capsule_id, contract_hash = contract_hash,
    tickets = ticket_values, ticket_json = ticket_json,
    ticket_negotiation_json = ticket_negotiation_json,
    allocation_openings = stats::setNames(lapply(
      designated, function(peer) .dsvert_joint_dp_client_json(list(
        version = "test-allocation-opening-v1", peer_name = peer))),
      designated),
    summaries = summary_values, summary_json = summary_json,
    summary_negotiation_json = summary_negotiation_json,
    source_capability_json = source_capability_json,
    transfer_ids = transfer_ids, chunk_count = chunk_count,
    conns = stats::setNames(as.list(paste0("connection-", peers)), peers),
    state = state, aggregate = aggregate)
}

.capsule_source_client_recursive_names <- function(value) {
  if (!is.list(value)) return(character())
  c(names(value), unlist(lapply(
    value, .capsule_source_client_recursive_names), use.names = FALSE))
}

test_that("source policy identity covers the lifetime-bound server common", {
  fixture <- .capsule_source_client_fixture(k = 3L)
  policy <- fixture$status[[fixture$peers[[1L]]]]$policy
  common <- list(
    protocol = "dsvert-joint-dp-control-v3",
    release_scope = "joint_mpc_single_opening",
    capability_id = "joint_mpc_single_opening_v1",
    domain = policy$domain,
    cohort_id = policy$cohort_id,
    ordered_peer_pinset = as.list(fixture$pins),
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    peer_count = length(fixture$pins),
    designated_noise_peers = fixture$designated,
    designated_noise_peer_pinset = as.list(
      fixture$pins[fixture$designated]),
    epsilon_capsule = as.numeric(policy$capsule_epsilon),
    delta_capsule = as.numeric(policy$capsule_delta),
    lifetime_max_distinct_capsules =
      as.numeric(policy$lifetime_max_distinct_capsules),
    lifetime_epsilon_upper_bound = policy$lifetime_epsilon_upper_bound,
    lifetime_delta_upper_bound = policy$lifetime_delta_upper_bound,
    privacy_accounting =
      "bounded_distinct_capsules_one_public_instance_each_v1",
    adjacency = policy$adjacency,
    patient_column = policy$patient_column,
    unit_capacity = policy$unit_capacity,
    max_records_per_unit = policy$max_records_per_unit,
    overflow_policy = policy$overflow_policy,
    sampler = policy$sampler)
  expected_hash <- .dsvert_dp_capsule_source_hash(common)
  observed <- .dsvert_dp_capsule_source_policy_identity(list(
    status = fixture$status,
    servers = fixture$peers,
    pinset = fixture$pins,
    designated = fixture$designated))

  expect_identical(observed$policy_contract_hash, expected_hash)
  expect_identical(observed$consortium_id, paste0("jdpc1_", expected_hash))
  expect_identical(
    fixture$manifest$capsule_identity$contract$policy_contract_hash,
    expected_hash)
  expect_identical(
    fixture$manifest$capsule_identity$contract$consortium_id,
    paste0("jdpc1_", expected_hash))

  legacy_common <- common
  legacy_common[c(
    "lifetime_max_distinct_capsules", "lifetime_epsilon_upper_bound",
    "lifetime_delta_upper_bound", "privacy_accounting")] <- NULL
  expect_false(identical(
    .dsvert_dp_capsule_source_hash(legacy_common), expected_hash))
})

test_that("source policy identity matches the server canonical golden", {
  pins <- c(
    peer_a = "AQIDBAUGBwgJCgsMDQ4PEBESExQVFhcYGRobHB0eHyA",
    peer_b = "ISIjJCUmJygpKissLS4vMDEyMzQ1Njc4OTo7PD0-P0A")
  policy <- list(
    domain = "joint-test-study", cohort_id = "joint-test-cohort",
    peer_pinset_sha256 =
      "b2a95e01d260e2f2cfe606c808bbb3e405a4e775413bae06d0c5938524a8ed13",
    capsule_epsilon = 1, capsule_delta = 1e-6,
    lifetime_max_distinct_capsules = 8,
    lifetime_epsilon_upper_bound = "8",
    lifetime_delta_upper_bound = "7.9999999999999996e-6",
    adjacency = "add_remove_patient", patient_column = "patient_id",
    unit_capacity = 100L, max_records_per_unit = 2L,
    overflow_policy = "reject_snapshot",
    sampler = "two_private_hmac_seeds_one_gc_sample_v1")
  status <- list(peer_a = list(policy = policy))
  observed <- .dsvert_dp_capsule_source_policy_identity(list(
    status = status, servers = c("peer_a", "peer_b"), pinset = pins,
    designated = c("peer_a", "peer_b")))

  expect_identical(
    observed$consortium_id,
    paste0("jdpc1_",
           "0a4af0da437caf103562897a7fd2aad2deab92f3adc8b72b3cdf8fcf7969576f"))
  expect_identical(
    observed$policy_contract_hash,
    "0a4af0da437caf103562897a7fd2aad2deab92f3adc8b72b3cdf8fcf7969576f")
})

test_that("source orchestration fails closed without both allocation openings", {
  fixture <- .capsule_source_client_fixture(k = 3L)
  expect_error(.dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    .aggregate = fixture$aggregate), "no canonical cross-signed allocation")
  swapped <- rev(fixture$allocation_openings)
  expect_error(.dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = swapped, .aggregate = fixture$aggregate),
    "no canonical cross-signed allocation")
})

test_that("K=2 through K=5 source orchestration is sticky, ordered and redacted", {
  for (k in 2:5) {
    fixture <- .capsule_source_client_fixture(k = k)
    first <- .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate)
    first_call_count <- length(fixture$state$calls)
    second <- .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate)
    expect_identical(second, first)
    expect_equal(length(fixture$state$calls), 2L * first_call_count)
    expect_identical(first$source_peers, fixture$sources)
    expect_identical(first$designated_noise_peers, fixture$designated)
    expect_true(first$sampler_handoff_ready)
    expect_false(first$payload_exposed)
    expect_false(first$operation_or_request_limit)
    expect_false(first$history_can_deny_operation)
    forbidden <- "cipher|share|seed|mask|plaintext|patient|value|digest"
    expect_false(any(grepl(
      forbidden, .capsule_source_client_recursive_names(first),
      ignore.case = TRUE)))
  }
})

test_that("one bundle fetch and one parallel accept serve each owner chunk", {
  fixture <- .capsule_source_client_fixture(k = 3L)
  .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  calls <- fixture$state$calls
  endpoints <- vapply(calls, `[[`, character(1L), "endpoint")
  expect_identical(endpoints[1:4], c(
    "dsvertJointDPCapsuleStatusDS",
    "dsvertDPCapsuleSourceTicketDS",
    "dsvertDPCapsuleSourceTicketDS",
    "dsvertDPCapsuleSourcePrepareDS"))
  expect_identical(calls[[1L]]$targets, fixture$peers)
  expect_identical(calls[[2L]]$targets, fixture$designated)
  expect_identical(calls[[2L]]$async, TRUE)
  expect_identical(calls[[3L]]$targets, fixture$designated)
  expect_identical(calls[[3L]]$async, FALSE)
  expect_identical(calls[[4L]]$targets, fixture$sources)
  prepare_arguments <- as.list(calls[[4L]]$expressions[[1L]])[-1L]
  expect_identical(prepare_arguments$first_opening_json,
                   fixture$allocation_openings[[1L]])
  expect_identical(prepare_arguments$second_opening_json,
                   fixture$allocation_openings[[2L]])

  expected_cycles <- length(fixture$sources) * fixture$chunk_count
  expect_equal(sum(endpoints == "dsvertDPCapsuleSourceChunkDS"),
               expected_cycles)
  expect_equal(sum(endpoints == "dsvertDPCapsuleSourceAcceptDS"),
               expected_cycles)
  tail_calls <- calls[-seq_len(4L)]
  expected_source <- rep(fixture$sources, each = fixture$chunk_count * 2L)
  observed_source <- vapply(seq_along(tail_calls), function(index) {
    item <- tail_calls[[index]]
    if (identical(item$endpoint, "dsvertDPCapsuleSourceChunkDS")) {
      item$targets[[1L]]
    } else {
      envelope <- fixture$state$decode(
        as.list(item$expressions[[1L]])[-1L]$envelope_json)
      envelope$source_name
    }
  }, character(1L))
  expect_identical(observed_source, expected_source)
  for (index in seq.int(2L, length(tail_calls), by = 2L)) {
    expect_identical(tail_calls[[index]]$endpoint,
                     "dsvertDPCapsuleSourceAcceptDS")
    expect_identical(tail_calls[[index]]$targets, fixture$designated)
  }
  expect_false("dsvertDPStatusDS" %in% endpoints)
})

test_that("byte windows of signed artifacts reduce calls without changing envelopes", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .capsule_source_client_fixture(
      k = k, coordinate_count = 32769L,
      ciphertext_bytes = 120000L, window_transport = TRUE)
    window_result <- .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate)
    window_calls <- fixture$state$calls
    window_accepted <- as.list(fixture$state$accepted)
    endpoints <- vapply(window_calls, `[[`, character(1L), "endpoint")
    transfer_calls <- window_calls[endpoints %in% c(
      "dsvertDPCapsuleSourceChunkDS", "dsvertDPCapsuleSourceAcceptDS")]
    observed_invocations <- sum(vapply(
      transfer_calls, function(value) length(value$targets), integer(1L)))
    scalar_invocations <- 3L * length(fixture$sources) *
      fixture$chunk_count
    expect_lt(observed_invocations, scalar_invocations)
    expect_equal(length(window_accepted),
                 2L * length(fixture$sources) * fixture$chunk_count)
    expect_true(any(vapply(transfer_calls, function(value) {
      if (!identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")) {
        return(FALSE)
      }
      length(as.list(value$expressions[[1L]])[-1L]$chunk_index) > 1L
    }, logical(1L))))
    chunk_calls <- Filter(function(value) {
      identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
    }, transfer_calls)
    expect_equal(
      vapply(chunk_calls, function(value) {
        length(as.list(value$expressions[[1L]])[-1L]$chunk_index)
      }, integer(1L)),
      rep(fixture$chunk_count, length(fixture$sources)))
    expect_true(all(vapply(chunk_calls, function(value) {
      identical(
        as.list(value$expressions[[1L]])[-1L]$transport_contract,
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)
    }, logical(1L))))
    accept_calls <- Filter(function(value) {
      identical(value$endpoint, "dsvertDPCapsuleSourceAcceptDS")
    }, transfer_calls)
    expect_true(all(vapply(accept_calls, function(value) {
      identical(
        as.list(value$expressions[[1L]])[-1L]$transport_contract,
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)
    }, logical(1L))))

    fixture$state$calls <- list()
    fixture$state$window_transport <- FALSE
    scalar_result <- .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate)
    scalar_calls <- fixture$state$calls
    scalar_endpoints <- vapply(
      scalar_calls, `[[`, character(1L), "endpoint")
    scalar_transfer <- scalar_calls[scalar_endpoints %in% c(
      "dsvertDPCapsuleSourceChunkDS", "dsvertDPCapsuleSourceAcceptDS")]
    expect_equal(sum(vapply(
      scalar_transfer, function(value) length(value$targets), integer(1L))),
      scalar_invocations)
    expect_identical(as.list(fixture$state$accepted), window_accepted)
    expect_identical(scalar_result, window_result)
  }
})

test_that("heterogeneous public node ceilings select a per-source window", {
  active_limits <- NULL
  local_mocked_bindings(
    .dsvert_dsi_transport_site_limits = function(conns, .aggregate) {
      lapply(active_limits, function(value) value[names(conns)])
    },
    .package = "dsVertClient")
  for (k in c(2L, 3L, 5L)) {
    fixture <- .capsule_source_client_fixture(
      k = k, coordinate_count = 9L * 8192L + 1L,
      ciphertext_bytes = 120000L, window_transport = TRUE)
    targets <- fixture$peers
    response_candidates <- rep(c(8, 4, 2) * 1024^2,
                               length.out = length(targets))
    active_limits <- list(
      request_payload_bytes = stats::setNames(
        rep(8 * 1024^2, length(targets)), targets),
      expression_bytes = stats::setNames(
        rep(8 * 1024^2 + .DSVERT_DSI_PROBE_EXPRESSION_RESERVE,
            length(targets)), targets),
      response_bytes = stats::setNames(vapply(
        response_candidates,
        .dsvert_dsi_response_probe_usable_bytes, numeric(1L)), targets),
      response_probe_supported = stats::setNames(
        rep(TRUE, length(targets)), targets))

    result <- .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate)
    expect_true(result$sampler_handoff_ready)
    chunk_calls <- Filter(function(value) {
      identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
    }, fixture$state$calls)
    for (source in fixture$sources) {
      observed <- vapply(Filter(function(value) {
        identical(value$targets, source)
      }, chunk_calls), function(value) {
        length(as.list(value$expressions[[1L]])[-1L]$chunk_index)
      }, integer(1L))
      window <- .dsvert_dp_capsule_source_effective_window(
        source, fixture$designated,
        .dsvert_dp_capsule_source_window_capability(), active_limits)
      full <- fixture$chunk_count %/% window
      remainder <- fixture$chunk_count %% window
      expected <- c(rep(window, full), if (remainder) remainder else integer())
      expect_equal(observed, expected, info = paste("K", k, source))
    }
    expect_equal(
      length(fixture$state$accepted),
      2L * length(fixture$sources) * fixture$chunk_count)
  }
})

test_that("effective W=1 negotiation preserves legacy artifacts exactly", {
  fixture <- .capsule_source_client_fixture(
    k = 3L, coordinate_count = 16385L, ciphertext_bytes = 250000L,
    window_transport = TRUE)
  local_mocked_bindings(
    .dsvert_dsi_transport_site_limits = function(conns, .aggregate) {
      targets <- names(conns)
      list(
        request_payload_bytes = stats::setNames(
          rep(768L * 1024L, length(targets)), targets),
        expression_bytes = stats::setNames(
          rep(768L * 1024L, length(targets)), targets),
        response_bytes = stats::setNames(
          rep(768L * 1024L, length(targets)), targets),
        response_probe_supported = stats::setNames(
          rep(TRUE, length(targets)), targets))
    },
    .package = "dsVertClient")
  for (peer in fixture$designated) {
    wrapper <- fixture$state$decode(
      fixture$ticket_negotiation_json[[peer]])
    expect_identical(wrapper$ticket_json, fixture$ticket_json[[peer]])
  }
  for (source in fixture$sources) {
    wrapper <- fixture$state$decode(
      fixture$summary_negotiation_json[[source]])
    expect_identical(wrapper$summary_json, fixture$summary_json[[source]])
  }
  first <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  first_accepted <- as.list(fixture$state$accepted)
  first_calls <- fixture$state$calls
  chunk_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
  }, first_calls)
  request_widths <- vapply(chunk_calls, function(value) {
    length(as.list(value$expressions[[1L]])[-1L]$chunk_index)
  }, integer(1L))
  expect_identical(
    request_widths,
    rep(1L, length(fixture$sources) * fixture$chunk_count))
  scalar_transfer <- Filter(function(value) {
    value$endpoint %in% c(
      "dsvertDPCapsuleSourceChunkDS", "dsvertDPCapsuleSourceAcceptDS")
  }, first_calls)
  expect_true(all(vapply(scalar_transfer, function(value) {
    !"transport_contract" %in% names(
      as.list(value$expressions[[1L]])[-1L])
  }, logical(1L))))
  fixture$state$window_transport <- FALSE
  second <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  expect_identical(second, first)
  expect_identical(as.list(fixture$state$accepted), first_accepted)
})

test_that("an Armadillo-style old server declines sync without poisoning", {
  fixture <- .capsule_source_client_fixture(
    k = 3L, coordinate_count = 16385L, window_transport = FALSE)
  fixture$state$negotiation_error <- TRUE
  poisoned <- 0L
  local_mocked_bindings(
    .dsvert_dsi_poison_sessions = function(...) poisoned <<- poisoned + 1L,
    .package = "dsVertClient")
  result <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  expect_true(result$sampler_handoff_ready)
  expect_identical(poisoned, 0L)
  ticket_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourceTicketDS")
  }, fixture$state$calls)
  expect_length(ticket_calls, 2L)
  expect_identical(vapply(ticket_calls, `[[`, logical(1L), "async"),
                   c(TRUE, FALSE))
  expect_identical(
    as.list(ticket_calls[[2L]]$expressions[[1L]])[-1L]$transport_contract,
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)
  chunk_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
  }, fixture$state$calls)
  expect_true(all(vapply(chunk_calls, function(value) {
    length(as.list(value$expressions[[1L]])[-1L]$chunk_index) == 1L
  }, logical(1L))))

  oversized <- .capsule_source_client_fixture(
    k = 3L, coordinate_count = 16385L, window_transport = TRUE)
  oversized$state$negotiation_oversize <- TRUE
  oversized_result <- .dsvert_dp_capsule_source_transport(
    oversized$manifest_json, oversized$conns, status = oversized$status,
    allocation_openings = oversized$allocation_openings,
    .aggregate = oversized$aggregate)
  expect_true(oversized_result$sampler_handoff_ready)
  oversized_chunks <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
  }, oversized$state$calls)
  expect_true(all(vapply(oversized_chunks, function(value) {
    length(as.list(value$expressions[[1L]])[-1L]$chunk_index) == 1L
  }, logical(1L))))
  expect_identical(poisoned, 0L)
})

test_that("a signed legacy 768 KiB capability falls back to scalar", {
  fixture <- .capsule_source_client_fixture(
    k = 3L, coordinate_count = 16385L, window_transport = TRUE)
  fixture$state$fault <- function(
      endpoint, result, targets, expressions, call) {
    if (!identical(endpoint, "dsvertDPCapsuleSourceTicketDS")) {
      return(result)
    }
    arguments <- as.list(expressions[[1L]])[-1L]
    if (!identical(
        arguments$transport_contract,
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)) return(result)
    for (peer in targets) {
      wrapper <- fixture$state$decode(result[[peer]])
      unsigned <- wrapper[setdiff(names(wrapper), "signature")]
      unsigned$capability <-
        .dsvert_dp_capsule_source_legacy_window_capability()
      result[[peer]] <- fixture$state$sign_json(
        unsigned, peer, "recipient-window-capability")
    }
    result
  }

  result <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  expect_true(result$sampler_handoff_ready)
  chunk_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
  }, fixture$state$calls)
  expect_true(all(vapply(chunk_calls, function(value) {
    length(as.list(value$expressions[[1L]])[-1L]$chunk_index) == 1L
  }, logical(1L))))
  capability_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourceTicketDS") &&
      identical(
        as.list(value$expressions[[1L]])[-1L]$transport_contract,
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ONLY_REQUEST_V2)
  }, fixture$state$calls)
  expect_length(capability_calls, 0L)
})

test_that("mixed v1 and v2 recipients or sources stay scalar", {
  for (mixed_role in c("recipient", "source")) {
    fixture <- .capsule_source_client_fixture(
      k = 3L, coordinate_count = 16385L, window_transport = TRUE)
    fixture$state$fault <- local({
      role <- mixed_role
      function(endpoint, result, targets, expressions, call) {
        if (!identical(endpoint, "dsvertDPCapsuleSourceTicketDS")) {
          return(result)
        }
        arguments <- as.list(expressions[[1L]])[-1L]
        recipient_phase <- identical(
          arguments$transport_contract,
          .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)
        source_phase <- identical(
          arguments$transport_contract,
          .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ONLY_REQUEST_V2)
        if ((identical(role, "recipient") && recipient_phase) ||
            (identical(role, "source") && source_phase)) {
          peer <- targets[[length(targets)]]
          wrapper <- fixture$state$decode(result[[peer]])
          unsigned <- wrapper[setdiff(names(wrapper), "signature")]
          unsigned$capability <-
            .dsvert_dp_capsule_source_legacy_window_capability()
          domain <- if (recipient_phase) {
            "recipient-window-capability"
          } else {
            "source-window-capability-advertisement"
          }
          result[[peer]] <- fixture$state$sign_json(unsigned, peer, domain)
        }
        result
      }
    })

    result <- .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate)
    expect_true(result$sampler_handoff_ready, info = mixed_role)
    transfer <- Filter(function(value) {
      value$endpoint %in% c(
        "dsvertDPCapsuleSourceChunkDS", "dsvertDPCapsuleSourceAcceptDS")
    }, fixture$state$calls)
    expect_true(all(vapply(transfer, function(value) {
      arguments <- as.list(value$expressions[[1L]])[-1L]
      index <- arguments$chunk_index
      (is.null(index) || length(index) == 1L) &&
        !"transport_contract" %in% names(arguments)
    }, logical(1L))), info = mixed_role)
  }
})

test_that("optional negotiation expression oversize falls back to scalar", {
  fixture <- .capsule_source_client_fixture(
    k = 3L, coordinate_count = 16385L, window_transport = TRUE)
  original_validate <- .dsvert_validate_dsi_expression_sizes
  poisoned <- 0L
  local_mocked_bindings(
    .dsvert_validate_dsi_expression_sizes = function(
        expressions, capacity_bytes = NULL) {
      calls <- if (is.list(expressions) && !is.call(expressions)) {
        expressions
      } else {
        list(expressions)
      }
      negotiation <- any(vapply(calls, function(value) {
        "transport_contract" %in% names(as.list(value)[-1L])
      }, logical(1L)))
      if (isTRUE(negotiation)) {
        stop(.dsvert_client_resource_oversize(
          requested_bytes = 768L * 1024L,
          capacity_bytes = 768L * 1024L - 1L,
          scope = "test optional negotiation expression"))
      }
      original_validate(expressions, capacity_bytes = capacity_bytes)
    },
    .dsvert_dsi_poison_sessions = function(...) poisoned <<- poisoned + 1L,
    .package = "dsVertClient")
  result <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  expect_true(result$sampler_handoff_ready)
  expect_identical(poisoned, 0L)
  chunk_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
  }, fixture$state$calls)
  expect_true(all(vapply(chunk_calls, function(value) {
    length(as.list(value$expressions[[1L]])[-1L]$chunk_index) == 1L
  }, logical(1L))))
})

test_that("optional negotiation preserves typed peer rejection evidence", {
  fixture <- .capsule_source_client_fixture(
    k = 3L, coordinate_count = 16385L, window_transport = TRUE)
  rejection <- .dsvert_client_peer_not_recognized_condition(
    "peer_a", strrep("1", 64L), strrep("2", 64L))
  aggregate <- function(conns, expr, async, error, errors.print) {
    expressions <- if (is.list(expr) && !is.call(expr)) expr else list(expr)
    arguments <- as.list(expressions[[1L]])[-1L]
    if (identical(async, FALSE) && identical(
        arguments$transport_contract,
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)) {
      stop(rejection)
    }
    fixture$aggregate(conns, expr, async, error, errors.print)
  }

  condition <- tryCatch(
    .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = aggregate),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  prepare_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourcePrepareDS")
  }, fixture$state$calls)
  expect_length(prepare_calls, 0L)
})

test_that("optional negotiation never reuses a poisoned DSI session", {
  calls <- 0L
  local_mocked_bindings(
    .dsvert_validate_real_dsi_transport = function(...) invisible(TRUE),
    .dsvert_dsi_job_session_key = function(connection) "poisoned-session",
    .dsvert_dsi_session_is_poisoned = function(key) TRUE,
    .package = "dsVertClient")
  condition <- tryCatch(
    .dsvert_dp_capsule_source_optional_negotiation(
      list(site = list()), list(site = call("publicCapabilityDS")),
      function(...) {
        calls <<- calls + 1L
        list(site = "unreachable")
      }),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(condition, "dsvert_dsi_poisoned_session")
  expect_identical(calls, 0L)
})

test_that("an ambiguous sync transport failure poisons before Prepare", {
  for (message in c(
      "connection reset after execute",
      "Internal server error: response failed after execute")) {
    fixture <- .capsule_source_client_fixture(
      k = 3L, coordinate_count = 16385L, window_transport = TRUE)
    .dsvert_dsi_clear_poisoned_sessions()
    withr::defer(.dsvert_dsi_clear_poisoned_sessions())
    keys <- stats::setNames(paste0("ambiguous-", fixture$designated),
                            fixture$designated)
    aggregate <- function(conns, expr, async, error, errors.print) {
      expressions <- if (is.list(expr) && !is.call(expr)) expr else list(expr)
      arguments <- as.list(expressions[[1L]])[-1L]
      if (identical(async, FALSE) && identical(
          arguments$transport_contract,
          .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)) {
        stop(simpleError(message))
      }
      fixture$aggregate(conns, expr, async, error, errors.print)
    }
    local_mocked_bindings(
      .dsvert_dsi_job_session_key = function(connection) {
        peer <- sub("^connection-", "", as.character(connection))
        unname(keys[[peer]])
      },
      .package = "dsVertClient")

    condition <- tryCatch(
      .dsvert_dp_capsule_source_transport(
        fixture$manifest_json, fixture$conns, status = fixture$status,
        allocation_openings = fixture$allocation_openings,
        .aggregate = aggregate),
      dsvert_dsi_poisoned_session = identity)
    expect_s3_class(condition, "dsvert_dsi_poisoned_session")
    expect_match(conditionMessage(condition), "fresh DSI login connections",
                 info = message)
    endpoints <- vapply(fixture$state$calls, `[[`, character(1L), "endpoint")
    expect_false(any(endpoints %in% c(
      "dsvertDPCapsuleSourcePrepareDS", "dsvertDPCapsuleSourceChunkDS")),
      info = message)
    expect_true(all(vapply(
      unname(keys), .dsvert_dsi_session_is_poisoned, logical(1L))),
      info = message)
  }
})

test_that("a mixed old source is detected before protected Prepare", {
  for (k in c(3L, 5L)) {
    fixture <- .capsule_source_client_fixture(
      k = k, coordinate_count = 16385L, window_transport = TRUE)
    old_source <- setdiff(fixture$sources, fixture$designated)[[1L]]
    fixture$state$source_capability_unsupported <- old_source
    poisoned <- 0L
    local_mocked_bindings(
      .dsvert_dsi_poison_sessions =
        function(...) poisoned <<- poisoned + 1L,
      .package = "dsVertClient")
    mixed <- .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate)
    mixed_accepted <- as.list(fixture$state$accepted)
    prepare_calls <- Filter(function(value) {
      identical(value$endpoint, "dsvertDPCapsuleSourcePrepareDS")
    }, fixture$state$calls)
    expect_length(prepare_calls, 1L)
    expect_true(all(vapply(prepare_calls[[1L]]$expressions, function(value) {
      arguments <- as.list(value)[-1L]
      first <- fixture$state$decode(arguments$first_ticket_json)
      identical(first$version,
                .DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_VERSION)
    }, logical(1L))))
    chunk_calls <- Filter(function(value) {
      identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
    }, fixture$state$calls)
    expect_true(all(vapply(chunk_calls, function(value) {
      length(as.list(value$expressions[[1L]])[-1L]$chunk_index) == 1L
    }, logical(1L))))
    expect_identical(poisoned, 0L)

    fixture$state$calls <- list()
    fixture$state$window_transport <- FALSE
    fixture$state$source_capability_unsupported <- character()
    scalar <- .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate)
    expect_identical(scalar, mixed)
    expect_identical(as.list(fixture$state$accepted), mixed_accepted)
  }
})

test_that("a forged source capability fails before protected Prepare", {
  fixture <- .capsule_source_client_fixture(
    k = 3L, coordinate_count = 16385L, window_transport = TRUE)
  fixture$state$fault <- function(
      endpoint, result, targets, expressions, call) {
    if (!identical(endpoint, "dsvertDPCapsuleSourceTicketDS")) {
      return(result)
    }
    arguments <- as.list(expressions[[1L]])[-1L]
    if (!identical(
        arguments$transport_contract,
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ONLY_REQUEST_V2)) {
      return(result)
    }
    source <- targets[[1L]]
    forged <- fixture$state$decode(result[[source]])
    forged$capability$maximum_window_chunks <- 1
    result[[source]] <- .dsvert_joint_dp_client_json(forged)
    result
  }

  expect_error(
    .dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate),
    "invalid source capability attestation")
  prepare_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourcePrepareDS")
  }, fixture$state$calls)
  expect_length(prepare_calls, 0L)
})

test_that("window downgrade, replay and malformed windows fail safely", {
  fixture <- .capsule_source_client_fixture(
    k = 3L, coordinate_count = 16385L, window_transport = TRUE)
  fixture$state$fault <- function(endpoint, result, targets, expressions, call) {
    if (identical(endpoint, "dsvertDPCapsuleSourcePrepareDS") &&
        identical(targets[[1L]], fixture$sources[[1L]])) {
      result[[targets[[1L]]]] <- fixture$summary_json[[targets[[1L]]]]
    }
    result
  }
  .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  chunk_calls <- Filter(function(value) {
    identical(value$endpoint, "dsvertDPCapsuleSourceChunkDS")
  }, fixture$state$calls)
  expect_true(all(vapply(chunk_calls, function(value) {
    length(as.list(value$expressions[[1L]])[-1L]$chunk_index) == 1L
  }, logical(1L))))

  for (fault_case in c("reorder", "duplicate", "tamper", "oversize")) {
    fixture <- .capsule_source_client_fixture(
      k = 3L, coordinate_count = 16385L, window_transport = TRUE)
    fixture$state$fault <- local({
      selected <- fault_case
      function(endpoint, result, targets, expressions, call) {
        if (!identical(endpoint, "dsvertDPCapsuleSourceChunkDS")) {
          return(result)
        }
        source <- targets[[1L]]
        value <- fixture$state$decode(result[[source]])
        if (!identical(value$version,
                       .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION)) {
          return(result)
        }
        if (identical(selected, "reorder")) {
          value$bundles <- rev(value$bundles)
        } else if (identical(selected, "duplicate")) {
          value$bundles[[2L]] <- value$bundles[[1L]]
        } else if (identical(selected, "tamper")) {
          value$bundles[[2L]]$envelopes[[1L]]$ciphertext_sha256 <-
            strrep("8", 64L)
        } else {
          result[[source]] <- strrep(
            "A", .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES + 1L)
          return(result)
        }
        result[[source]] <- .dsvert_joint_dp_client_json(value)
        result
      }
    })
    expect_error(.dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate), "invalid|conflicting", ignore.case = TRUE)
  }

  fixture <- .capsule_source_client_fixture(
    k = 5L, source_count = 2L, coordinate_count = 16385L,
    window_transport = TRUE)
  first <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  accepted <- as.list(fixture$state$accepted)
  second <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  expect_identical(second, first)
  expect_identical(as.list(fixture$state$accepted), accepted)
})

test_that("tickets and summaries reject pin, recipient and source substitution", {
  fixture <- .capsule_source_client_fixture(k = 3L, coordinate_count = 7L)
  bad_status <- fixture$status
  bad_status$peer_c$policy$own_identity_pk <- fixture$pins[["peer_a"]]
  expect_error(.dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = bad_status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate), "invalid reusable joint-DP capsule policy")

  fixture$state$fault <- function(endpoint, result, targets, expressions, call) {
    if (!identical(endpoint, "dsvertDPCapsuleSourceTicketDS")) return(result)
    peer <- fixture$designated[[1L]]
    value <- fixture$state$decode(result[[peer]])
    value$recipient_name <- fixture$designated[[2L]]
    unsigned <- value[setdiff(names(value), "signature")]
    result[[peer]] <- fixture$state$sign_json(
      unsigned, peer, "recipient-ticket")
    result
  }
  expect_error(.dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate), "invalid capsule source ticket")

  fixture <- .capsule_source_client_fixture(k = 3L, coordinate_count = 7L)
  fixture$state$fault <- function(endpoint, result, targets, expressions, call) {
    if (!identical(endpoint, "dsvertDPCapsuleSourcePrepareDS")) return(result)
    source <- fixture$sources[[1L]]
    value <- fixture$state$decode(result[[source]])
    value$source_name <- fixture$sources[[2L]]
    unsigned <- value[setdiff(names(value), "signature")]
    result[[source]] <- fixture$state$sign_json(
      unsigned, source, "source-summary")
    result
  }
  expect_error(.dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate), "invalid capsule source summary")
})

test_that("envelope signature, body, hash, recipient and order fail closed", {
  cases <- c("signature", "ciphertext", "hash", "recipient", "order")
  for (case in cases) {
    fixture <- .capsule_source_client_fixture(
      k = 3L, coordinate_count = 7L)
    fixture$state$fault <- local({
      fault_case <- case
      function(endpoint, result, targets, expressions, call) {
        if (!identical(endpoint, "dsvertDPCapsuleSourceChunkDS")) {
          return(result)
        }
        source <- targets[[1L]]
        bundle <- fixture$state$decode(result[[source]])
        if (identical(fault_case, "order")) {
          bundle$envelopes <- rev(bundle$envelopes)
        } else {
          envelope <- bundle$envelopes[[1L]]
          if (identical(fault_case, "signature")) {
            envelope$signature <- paste0(
              substr(envelope$signature, 1L, 85L),
              if (substr(envelope$signature, 86L, 86L) == "A") "B" else "A")
          } else if (identical(fault_case, "ciphertext")) {
            envelope$ciphertext <- paste0(
              if (substr(envelope$ciphertext, 1L, 1L) == "A") "B" else "A",
              substr(envelope$ciphertext, 2L, nchar(envelope$ciphertext)))
          } else if (identical(fault_case, "hash")) {
            envelope$ciphertext_sha256 <- strrep("9", 64L)
            unsigned <- envelope[setdiff(names(envelope), "signature")]
            envelope <- fixture$state$sign_value(
              unsigned, source, "encrypted-chunk", TRUE)
          } else {
            envelope$recipient_name <- fixture$designated[[2L]]
            unsigned <- envelope[setdiff(names(envelope), "signature")]
            envelope <- fixture$state$sign_value(
              unsigned, source, "encrypted-chunk", TRUE)
          }
          bundle$envelopes[[1L]] <- envelope
        }
        result[[source]] <- .dsvert_joint_dp_client_json(bundle)
        result
      }
    })
    expect_error(.dsvert_dp_capsule_source_transport(
      fixture$manifest_json, fixture$conns, status = fixture$status,
      allocation_openings = fixture$allocation_openings,
      .aggregate = fixture$aggregate),
      "invalid|conflicting", ignore.case = TRUE)
  }
})

test_that("ACK binding, completion, partial responses and conflicts fail closed", {
  fixture <- .capsule_source_client_fixture(k = 2L, coordinate_count = 7L)
  fixture$state$fault <- function(endpoint, result, targets, expressions, call) {
    if (!identical(endpoint, "dsvertDPCapsuleSourceAcceptDS")) return(result)
    peer <- fixture$designated[[1L]]
    value <- fixture$state$decode(result[[peer]])
    value$capsule_aggregation_complete <- FALSE
    unsigned <- value[setdiff(names(value), "signature")]
    result[[peer]] <- fixture$state$sign_json(
      unsigned, peer, "aggregate-ack")
    result
  }
  expect_error(.dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate), "invalid capsule source acknowledgement")

  fixture <- .capsule_source_client_fixture(k = 2L, coordinate_count = 7L)
  fixture$state$fault <- function(endpoint, result, targets, expressions, call) {
    if (identical(endpoint, "dsvertDPCapsuleSourceAcceptDS")) {
      result[[targets[[2L]]]] <- NULL
    }
    result
  }
  expect_error(.dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate), "DataSHIELD transport failed")

  fixture <- .capsule_source_client_fixture(k = 2L, coordinate_count = 7L)
  first <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  expect_true(first$sampler_handoff_ready)
  fixture$state$ciphertext_epoch <- 1L
  expect_error(.dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate), "DataSHIELD transport failed")
})

test_that("near-bound bundles stay one response and final output stays redacted", {
  fixture <- .capsule_source_client_fixture(
    k = 2L, coordinate_count = 7L, source_count = 1L,
    ciphertext_bytes = 250000L)
  result <- .dsvert_dp_capsule_source_transport(
    fixture$manifest_json, fixture$conns, status = fixture$status,
    allocation_openings = fixture$allocation_openings,
    .aggregate = fixture$aggregate)
  expect_gt(fixture$state$bundle_sizes[[1L]], 640L * 1024L)
  expect_lt(fixture$state$bundle_sizes[[1L]],
            .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_BUNDLE_BYTES)
  endpoints <- vapply(
    fixture$state$calls, `[[`, character(1L), "endpoint")
  expect_equal(sum(endpoints == "dsvertDPCapsuleSourceChunkDS"), 1L)
  expect_equal(sum(endpoints == "dsvertDPCapsuleSourceAcceptDS"), 1L)
  expect_false(any(grepl(
    "cipher|share|seed|mask|plaintext|patient|value|digest",
    .capsule_source_client_recursive_names(result), ignore.case = TRUE)))
})

test_that("manifest v6 lifecycle and consortium binding are exact", {
  fixture <- .capsule_source_client_fixture(k = 3L, coordinate_count = 7L)
  context <- .dsvert_joint_dp_client_context(
    fixture$conns, status = fixture$status,
    .aggregate = fixture$aggregate)
  tamper <- function(change) {
    value <- fixture$manifest
    value <- change(value)
    .dsvert_joint_dp_client_json(value)
  }
  expect_error(.dsvert_dp_capsule_source_manifest(tamper(function(value) {
    value$version <- "dsvert-biomedical-capsule-workload-v1"
    value
  }), context), "manifest is invalid or misbound")
  expect_error(.dsvert_dp_capsule_source_manifest(tamper(function(value) {
    value$capsule_schema <- "dsvert-biomedical-capsule-workload-v1"
    value
  }), context), "manifest is invalid or misbound")
  expect_error(.dsvert_dp_capsule_source_manifest(tamper(function(value) {
    value$workload$registered_release_lifecycle$sampler_integration <- FALSE
    value$capsule_identity$contract$workload <- value$workload
    value$capsule_identity$capsule_id <- .dsvert_dp_capsule_source_hash(
      value$capsule_identity$contract)
    value
  }), context), "manifest is invalid or misbound")
  expect_error(.dsvert_dp_capsule_source_manifest(tamper(function(value) {
    value$capsule_identity$contract$consortium_id <-
      paste0("jdpc1_", strrep("f", 64L))
    value$capsule_identity$capsule_id <- .dsvert_dp_capsule_source_hash(
      value$capsule_identity$contract)
    value
  }), context), "manifest is invalid or misbound")

  invalid_coverage <- fixture$manifest
  invalid_coverage$workload$package_family_coverage_complete <- TRUE
  invalid_coverage$capsule_identity$contract$workload <-
    invalid_coverage$workload
  invalid_coverage$capsule_identity$capsule_id <-
    .dsvert_dp_capsule_source_hash(
      invalid_coverage$capsule_identity$contract)
  expect_error(.dsvert_dp_capsule_source_manifest(
    .dsvert_joint_dp_client_json(invalid_coverage), context),
    "manifest is invalid or misbound")
  validated <- .dsvert_dp_capsule_source_manifest(
    fixture$manifest_json, context)
  expect_identical(
    validated$value$workload$registered_release_lifecycle,
    .dsvert_dp_capsule_client_registered_release_lifecycle())
  lifecycle <- validated$value$workload$registered_release_lifecycle
  for (name in names(lifecycle)) {
    downgraded <- fixture$manifest
    value <- lifecycle[[name]]
    downgraded$workload$registered_release_lifecycle[[name]] <-
      if (is.logical(value)) !value else paste0(value, "-tampered")
    downgraded$capsule_identity$contract$workload <- downgraded$workload
    downgraded$capsule_identity$capsule_id <-
      .dsvert_dp_capsule_source_hash(
        downgraded$capsule_identity$contract)
    expect_error(.dsvert_dp_capsule_source_manifest(
      .dsvert_joint_dp_client_json(downgraded), context),
      "manifest is invalid or misbound",
      info = paste("downgraded lifecycle field", name))
  }
})

test_that("client validates the signed server-authoritative primitive scope", {
  fixture <- .capsule_source_client_fixture(k = 3L, coordinate_count = 7L)
  scope <- .dsvert_dp_capsule_primitive_scope_validate(fixture$manifest)
  expect_identical(scope$mode, "catalog_v1")
  expect_false(scope$analyst_expandable)
  expect_false(scope$client_query_can_add_coordinates)
  expect_identical(scope$projected_cost$projected_coordinate_count, 7)

  all_schema <- fixture$manifest
  all_schema$workload$primitive_scope$mode <- "all_schema"
  all_schema$workload$primitive_scope$selection$mode <- "all_schema"
  all_schema$workload$primitive_scope$selection_sha256 <-
    .dsvert_dp_capsule_source_hash(
      all_schema$workload$primitive_scope$selection)
  all_schema$workload$primitive_scope$projected_cost$
    automatic_pair_expansion <- "explicit_all_schema"
  all_schema$workload$primitive_scope$projected_cost$scaling_contract <-
    "explicit_schema_wide_pair_families_may_be_quadratic"
  expect_identical(
    .dsvert_dp_capsule_primitive_scope_validate(all_schema)$mode,
    "all_schema")

  forged <- fixture$manifest
  forged$workload$primitive_scope$analyst_expandable <- TRUE
  expect_error(
    .dsvert_dp_capsule_primitive_scope_validate(forged),
    "primitive scope")

  forged <- fixture$manifest
  forged$workload$primitive_scope$selection$included$
    numeric_moments <- list("undeclared")
  expect_error(
    .dsvert_dp_capsule_primitive_scope_validate(forged),
    "primitive scope")

  forged <- fixture$manifest
  forged$workload$primitive_scope$projected_cost$
    projected_coordinate_count <- 8
  expect_error(
    .dsvert_dp_capsule_primitive_scope_validate(forged),
    "primitive scope")

  extra <- fixture$manifest
  extra$workload$primitive_scope$selection$explicit_catalog$unexpected <-
    list()
  extra$workload$primitive_scope$selection_sha256 <-
    .dsvert_dp_capsule_source_hash(
      extra$workload$primitive_scope$selection)
  expect_error(
    .dsvert_dp_capsule_primitive_scope_validate(extra),
    "primitive scope")

  missing <- fixture$manifest
  missing$workload$primitive_scope$selection$
    referenced_by_signed_specs$gaussian <- NULL
  missing$workload$primitive_scope$selection_sha256 <-
    .dsvert_dp_capsule_source_hash(
      missing$workload$primitive_scope$selection)
  expect_error(
    .dsvert_dp_capsule_primitive_scope_validate(missing),
    "primitive scope")

  noncanonical_empty <- fixture$manifest
  noncanonical_empty$workload$primitive_scope$selection$included[
    "numeric_moments"] <- list(NULL)
  noncanonical_empty$workload$primitive_scope$selection_sha256 <-
    .dsvert_dp_capsule_source_hash(
      noncanonical_empty$workload$primitive_scope$selection)
  expect_error(
    .dsvert_dp_capsule_primitive_scope_validate(noncanonical_empty),
    "primitive scope")

  reordered <- fixture$manifest
  catalog <- reordered$workload$primitive_scope$selection$explicit_catalog
  reordered$workload$primitive_scope$selection$explicit_catalog <-
    catalog[rev(names(catalog))]
  expect_error(
    .dsvert_dp_capsule_primitive_scope_validate(reordered),
    "primitive scope")

  reordered <- fixture$manifest
  references <- reordered$workload$primitive_scope$selection$
    referenced_by_signed_specs
  reordered$workload$primitive_scope$selection$
    referenced_by_signed_specs <- references[rev(names(references))]
  expect_error(
    .dsvert_dp_capsule_primitive_scope_validate(reordered),
    "primitive scope")
})

test_that("source bridge has no legacy status, quota or request-history path", {
  code <- paste(readLines(.dsvert_client_source_file(
    "dp_capsule_source_transport.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("dsvertDPStatusDS", code, fixed = TRUE))
  expect_false(grepl("remaining_epsilon|remaining_delta|quota", code))
  for (endpoint in c(
      "dsvertDPCapsuleSourceTicketDS",
      "dsvertDPCapsuleSourcePrepareDS",
      "dsvertDPCapsuleSourceChunkDS",
      "dsvertDPCapsuleSourceAcceptDS")) {
    expect_match(code, endpoint, fixed = TRUE)
  }
})
