.vector_client_b64url <- function(value) {
  sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))),
    perl = TRUE)
}

.vector_client_laplace_plan <- function(coordinate_count,
                                        sensitivity_steps) {
  list(
    version = "dsvert-joint-dp-vector-laplace-plan-v3",
    sampler = .DSVERT_CLIENT_VECTOR_EXACT_SAMPLER,
    stop_bits = 128,
    stop_numerator = "1", uniform_bits = 128,
    binary_geometric_bits = 1,
    bernoulli_thresholds = list("1"),
    sensitivity_steps = sensitivity_steps,
    total_coordinate_count = as.numeric(coordinate_count),
    epsilon_effective_upper_numerator = "1",
    epsilon_effective_upper_denominator = "1",
    one_geometric_tv_numerator = "1",
    one_geometric_tv_denominator =
      "1267650600228229401496703205376",
    tail_upper_numerator = "1",
    tail_upper_denominator = "1267650600228229401496703205376",
    rounding_upper_numerator = "1",
    rounding_upper_denominator = "1267650600228229401496703205376",
    implementation_delta_numerator = "1",
    implementation_delta_denominator =
      "1267650600228229401496703205376",
    implementation_delta_bound = "7.888609052210118e-31",
    maximum_noise_magnitude = "1048576",
    maximum_chunk_coordinates = 128,
    private_stream_bytes_per_coordinate = 64,
    accounting = paste(
      "global iid discrete Laplace calibrated once to the workload joint",
      "L1 sensitivity; exact binary-geometric coupling"),
    capability_available = TRUE)
}

.vector_client_convolution_plan <- function(coordinate_count,
                                            sensitivity_steps) {
  plan <- .vector_client_laplace_plan(
    coordinate_count, sensitivity_steps)
  plan$version <- paste0(
    "dsvert-joint-dp-vector-independent-full-draw-",
    "convolution-plan-v3")
  plan$sampler <- .DSVERT_CLIENT_VECTOR_SAMPLER
  plan$maximum_chunk_coordinates <- as.numeric(min(
    .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES, coordinate_count))
  plan$independent_noise_peer_count <- 2
  plan$complete_epsilon_per_peer <- TRUE
  plan$epsilon_divided_by_peer_count <- FALSE
  plan$geometric_variables_per_peer_per_coordinate <- 2
  plan$geometric_variables_total_per_coordinate <- 4
  plan$per_peer_implementation_delta_numerator <-
    plan$implementation_delta_numerator
  plan$per_peer_implementation_delta_denominator <-
    plan$implementation_delta_denominator
  plan$per_peer_implementation_delta_bound <-
    plan$implementation_delta_bound
  plan$release_implementation_delta_aggregation <-
    "max_per_peer_not_sum"
  plan$two_peer_ideal_transfer_delta_numerator <- "0"
  plan$two_peer_ideal_transfer_delta_denominator <- "1"
  plan$two_peer_ideal_transfer_delta_bound <- "0"
  plan$threat_model <- "one pinned noise peer remains non-colluding"
  plan$privacy_argument <- paste(
    "each complete-epsilon draw is DP; adding the other independent",
    "draw is post-processing")
  plan
}

.vector_client_gaussian_plan <- function(coordinate_count) {
  .client_complete_gaussian_plan_v2(list(
    version = .DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION,
    mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
    sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER,
    request_binding_sha256 = digest::digest(
      "gaussian-request", "sha256", serialize = FALSE),
    total_coordinate_count = as.numeric(coordinate_count),
    maximum_chunk_coordinates = as.numeric(min(
      .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES, coordinate_count)),
    independent_noise_peer_count = 2,
    complete_epsilon_per_peer = TRUE,
    epsilon_divided_by_peer_count = FALSE,
    capability_available = TRUE,
    release_delta_aggregation = "max_per_peer_not_sum",
    exact_rational_sampler = FALSE,
    finite_support_transfer_charged = TRUE,
    fixed_work_sampler = TRUE,
    sampler_branches_on_protected_values = FALSE,
    sampler_branches_on_private_randomness = FALSE,
    host_constant_time_claim = FALSE,
    transcript_dp_claim = TRUE,
    sampler_candidate_count = 1,
    sampler_random_bits_per_coordinate = 128,
    sampler_random_bytes_per_coordinate = 17,
    sampler_table_precision_bits = 192,
    sampler_magnitude_count = 33,
    sampler_search_steps = 6,
    vector_tail_tv_upper_numerator = "1",
    vector_tail_tv_upper_denominator =
      "1267650600228229401496703205376",
    vector_sampler_tv_upper_numerator = "1",
    vector_sampler_tv_upper_denominator =
      "1267650600228229401496703205376",
    vector_total_tv_upper_numerator = "2",
    vector_total_tv_upper_denominator =
      "1267650600228229401496703205376",
    observable_worker_shape = "fixed dyadic CDF fixture",
    per_peer_implementation_delta_numerator = "1",
    per_peer_implementation_delta_denominator =
      "1267650600228229401496703205376"))
}

.vector_client_fixture <- function(k = 3L, coordinate_count = 3L,
                                   gaussian = FALSE,
                                   epochs = seq_len(k), keys = NULL,
                                   laplace_plan = NULL) {
  peers <- paste0("site_", letters[seq_len(k)])
  designated <- peers[1:2]
  if (!is.numeric(epochs) || length(epochs) != k || anyNA(epochs) ||
      any(!is.finite(epochs)) || any(epochs < 1) ||
      any(epochs != floor(epochs))) stop("invalid fixture epochs")
  if (is.null(keys)) {
    keys <- stats::setNames(lapply(peers, function(peer) {
      openssl::ed25519_keygen()
    }), peers)
  }
  pins <- vapply(keys, function(key) {
    .vector_client_b64url(tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  status <- stats::setNames(lapply(seq_along(peers), function(index) list(
    policy = list(capsule_epsilon = 1, capsule_delta = 2^-100),
    noise_root = list(
      privacy_epoch = as.numeric(epochs[[index]]),
      provider_id = paste0("provider-", peers[[index]]),
      key_id = paste0(
        "noise-root-", peers[[index]], "-epoch-", epochs[[index]])),
    release_domain = list(
      version = "dsvert-joint-dp-release-domain-v1",
      generation = as.numeric(epochs[[index]]),
      domain_id = paste0("rd_", digest::digest(paste0(
        "vector-client/", peers[[index]], "/", epochs[[index]]),
        "sha256", serialize = FALSE)),
      rotation_count = as.numeric(epochs[[index]] - 1),
      automatic_generation = TRUE, automatic_rotation = TRUE,
      snapshot_derived = FALSE, key_material_exposed = FALSE))),
    peers)
  conns <- stats::setNames(lapply(peers, function(peer) {
    structure(list(peer = peer), class = "fake")
  }), peers)
  context <- list(
    status = status, servers = peers, pinset = pins,
    designated = designated, conns = conns[designated],
    all_conns = conns)
  capsule_id <- digest::digest("vector-client-capsule", "sha256",
                               serialize = FALSE)
  marginal_artifacts <- list()
  if (coordinate_count > 1L) {
    marginal_artifacts$fixture_levels <- list(
      owner_peer = peers[[1L]], dataset = "cohort", column = "fixture",
      levels = paste0("level_", seq_len(coordinate_count - 1L)))
  }
  sensitivity_steps <- if (isTRUE(gaussian)) {
    "3.1415926535897940e+00"
  } else "1048576"
  exact_plan <- if (!isTRUE(gaussian)) {
    if (is.null(laplace_plan)) {
      .vector_client_laplace_plan(coordinate_count, sensitivity_steps)
    } else laplace_plan
  } else NULL
  selected_backend <- if (!isTRUE(gaussian) && coordinate_count > 1L) {
    .DSVERT_CLIENT_VECTOR_BACKEND
  } else NULL
  mechanism_plan <- if (isTRUE(gaussian)) {
    .vector_client_gaussian_plan(coordinate_count)
  } else if (!is.null(selected_backend)) {
    .vector_client_convolution_plan(coordinate_count, sensitivity_steps)
  } else {
    exact_plan
  }
  capsule_mechanism <- if (isTRUE(gaussian)) {
    list(
      mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
      sensitivity_norm = "l2",
      certificate = list(
        gaussian_calibration_request = list(
          epsilon = "9.9999999999999989e-01",
          delta = format(2^-100, digits = 17L,
                         scientific = TRUE, trim = TRUE),
          l2_sensitivity_steps = sensitivity_steps,
          total_coordinate_count = as.numeric(coordinate_count)),
        gaussian_plan = mechanism_plan,
        gaussian_plan_sha256 = .dsvert_vector_hash(mechanism_plan)))
  } else {
    list(mechanism = "discrete-laplace", sensitivity_norm = "l1")
  }
  profile <- .dsvert_vector_profile(
    capsule_mechanism, backend = selected_backend)
  manifest <- list(
    workload = list(
      coordinate_count = as.numeric(coordinate_count),
      capsule_mechanism = capsule_mechanism,
      families = list(
        admitted_count = list(
          owner_peer = peers[[1L]], dataset = "cohort"),
        numeric_moments = list(artifacts = list()),
        numeric_pair_moments = list(artifacts = list()),
        gaussian_models = list(artifacts = list()),
        fixed_numeric_histograms = list(artifacts = list()),
        categorical_marginals = list(artifacts = marginal_artifacts),
        categorical_pairs = list(sets = list()), correlation_artifacts = list(),
        describe_artifacts = list(), survival_artifacts = list())),
    capsule_identity = list(capsule_id = capsule_id))
  coordinate_order_sha256 <-
    .dsvert_dp_capsule_vector_layout(manifest)$sha256
  manifest_bundle <- list(
    manifest_json = "canonical-server-manifest",
    manifest_sha256 = digest::digest(
      "canonical-server-manifest", "sha256", serialize = FALSE),
    capsule_id = capsule_id, context = context)
  if (!isTRUE(gaussian)) {
    promoted <- coordinate_count <= 1L
    selection_reason <- if (promoted) {
      "within_public_exact_gc_cost_ceiling"
    } else {
      "above_public_exact_gc_cost_ceiling"
    }
    assessment_unsigned <- list(
      version = "dsvert-joint-dp-vector-exact-gc-assessment-v2",
      manifest_sha256 = manifest_bundle$manifest_sha256,
      representable = TRUE,
      exact_gc_capability_id =
        "joint_dp_biomedical_vector_exact_gc_v1",
      plan_sha256 = .dsvert_vector_hash(exact_plan),
      maximum_chunk_coordinates = 128L,
      cost_policy_version =
        "dsvert-joint-dp-vector-exact-gc-cost-policy-v1",
      total_coordinate_count = as.integer(coordinate_count),
      maximum_promoted_coordinates = 1L,
      promoted = promoted, selection_reason = selection_reason,
      private_material_accessed = FALSE,
      runtime_failure_consulted = FALSE)
    backend_assessment <- c(assessment_unsigned, list(
      assessment_sha256 =
        .dsvert_joint_dp_vector_exact_gc_client_hash(
          assessment_unsigned)))
    selection_unsigned <- list(
      version = "dsvert-joint-dp-vector-backend-selection-v2",
      manifest_sha256 = manifest_bundle$manifest_sha256,
      backend = profile$backend, one_draw = promoted,
      cost_policy_version =
        "dsvert-joint-dp-vector-exact-gc-cost-policy-v1",
      total_coordinate_count = as.integer(coordinate_count),
      maximum_promoted_coordinates = 1L,
      selection_reason = selection_reason,
      assessment_sha256 = backend_assessment$assessment_sha256,
      exact_gc_plan_sha256 = backend_assessment$plan_sha256,
      exact_gc_maximum_chunk_coordinates = 128L,
      selected_before_private_material = TRUE,
      retry_may_change_backend = FALSE)
    backend_selection <- c(selection_unsigned, list(
      selection_sha256 =
        .dsvert_joint_dp_vector_exact_gc_client_hash(
          selection_unsigned)))
  } else {
    backend_assessment <- backend_selection <- NULL
  }
  release_instance <- .dsvert_vector_release_instance(context, capsule_id)
  contract <- list(
    capsule_id = capsule_id,
    release_instance_id = release_instance$id,
    release_instance = release_instance$value,
    release_contract_hash = digest::digest(
      "release-contract", "sha256", serialize = FALSE),
    transcript_hash = digest::digest(
      "vector-transcript", "sha256", serialize = FALSE),
    manifest_sha256 = manifest_bundle$manifest_sha256,
    coordinate_count = as.integer(coordinate_count),
    chunk_coordinates = if (isTRUE(profile$exact_gc)) 128L else
      .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES,
    chunk_count = as.integer(ceiling(
      coordinate_count / if (isTRUE(profile$exact_gc)) 128L else
        .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES)),
    source_contract_hash = digest::digest(
      "source-contract", "sha256", serialize = FALSE),
    epsilon = "1e+00", delta = format(2^-100, digits = 17L,
                                      scientific = TRUE, trim = TRUE),
    sensitivity_steps = sensitivity_steps)
  contract$mechanism_plan <- mechanism_plan
  contract$plan_sha256 <- .dsvert_vector_hash(mechanism_plan)
  contract$profile <- profile
  contract$backend_assessment <- backend_assessment
  contract$backend_selection <- backend_selection
  contract$backend_selection_sha256 <- if (is.list(backend_selection)) {
    backend_selection$selection_sha256
  } else NULL
  contract$one_joint_draw <- isTRUE(profile$exact_gc)

  sign_receipt <- function(unsigned, peer) {
    message <- charToRaw(paste0(
      .DSVERT_CLIENT_JOINT_DP_RECEIPT_DOMAIN,
      .dsvert_joint_dp_client_json(unsigned)))
    .dsvert_joint_dp_client_json(c(unsigned, list(
      signature = .vector_client_b64url(
        openssl::ed25519_sign(message, keys[[peer]])))))
  }
  common <- function(peer, version, phase) list(
    version = version, phase = phase, capsule_id = capsule_id,
    release_instance_id = release_instance$id,
    release_instance = release_instance$value,
    release_contract_hash = contract$release_contract_hash,
    transcript_hash = contract$transcript_hash, peer_name = peer,
    peer_identity_pk = unname(pins[[peer]]),
    coordinate_count = as.numeric(coordinate_count),
    chunk_count = as.numeric(contract$chunk_count),
    backend = profile$backend,
    sampler = profile$sampler,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE)

  allocation_consortium <- paste0(
    "jdpc1_", digest::digest("vector-client-consortium", "sha256",
                              serialize = FALSE))
  allocation_prefix <- function(peer, version, phase) list(
    version = version, phase = phase,
    consortium_id = allocation_consortium,
    peer_name = peer, peer_identity_pk = unname(pins[[peer]]),
    capsule_id = capsule_id, query_id = capsule_id)
  allocation_mechanism <- digest::digest(
    "vector-client-allocation-mechanism", "sha256", serialize = FALSE)
  allocation_joint <- digest::digest(
    "vector-client-allocation-joint", "sha256", serialize = FALSE)
  allocation_prepare_set <- digest::digest(
    "vector-client-allocation-prepare-set", "sha256", serialize = FALSE)
  allocation_chain <- digest::digest(
    "vector-client-allocation-chain", "sha256", serialize = FALSE)
  allocation_commit_set <- digest::digest(
    "vector-client-allocation-commit-set", "sha256", serialize = FALSE)
  allocation_authorization_set <- digest::digest(
    "vector-client-allocation-authorization-set", "sha256",
    serialize = FALSE)
  allocation_prepare <- stats::setNames(lapply(designated, function(peer) {
    root <- status[[peer]]$noise_root
    sign_receipt(c(allocation_prefix(
      peer, .DSVERT_CLIENT_VECTOR_ALLOC_PREPARE_VERSION, "prepared"), list(
        privacy_epoch = as.character(root$privacy_epoch),
        noise_key_id = root$key_id,
        mechanism_hash = allocation_mechanism,
        allocation_index = "0", epsilon = "1e+00",
        delta = format(2^-100, digits = 17L,
                       scientific = TRUE, trim = TRUE),
        previous_chain = "GENESIS",
        snapshot_binding = digest::digest(
          paste0("allocation-snapshot-", peer), "sha256",
          serialize = FALSE),
        seed_commitment = digest::digest(
          paste0("allocation-seed-", peer), "sha256", serialize = FALSE),
        rollback_mode = "cross_signed_two_peer")), peer)
  }), designated)
  allocation_commit <- stats::setNames(lapply(designated, function(peer) {
    sign_receipt(c(allocation_prefix(
      peer, .DSVERT_CLIENT_VECTOR_ALLOC_COMMIT_VERSION, "committed"), list(
        allocation_index = "0", previous_chain = "GENESIS",
        new_chain = allocation_chain,
        joint_record_hash = allocation_joint,
        prepare_set_hash = allocation_prepare_set,
        seed_commitment = digest::digest(
          paste0("allocation-seed-", peer), "sha256", serialize = FALSE),
        external_anchor = FALSE)), peer)
  }), designated)
  allocation_authorization <- stats::setNames(lapply(
    designated, function(peer) sign_receipt(c(allocation_prefix(
      peer, .DSVERT_CLIENT_VECTOR_ALLOC_AUTHORIZE_VERSION, "authorized"),
      list(
        allocation_index = "0", new_chain = allocation_chain,
        joint_record_hash = allocation_joint,
        commit_set_hash = allocation_commit_set,
        peer_commit_stored = TRUE)), peer)), designated)
  allocation_opening <- stats::setNames(lapply(designated, function(peer) {
    sign_receipt(c(allocation_prefix(
      peer, .DSVERT_CLIENT_VECTOR_ALLOC_OPEN_VERSION, "open_authorized"),
      list(
        allocation_index = "0", new_chain = allocation_chain,
        joint_record_hash = allocation_joint,
        authorization_set_hash = allocation_authorization_set,
        seed_commitment = digest::digest(
          paste0("allocation-seed-", peer), "sha256", serialize = FALSE),
        release_scope = "joint_mpc_single_opening",
        capability_available = FALSE)), peer)
  }), designated)

  prepare <- stats::setNames(lapply(designated, function(peer) {
    sign_receipt(c(common(
      peer, .DSVERT_CLIENT_VECTOR_PREPARE_VERSION, "vector_prepared"), list(
        source_contract_hash = contract$source_contract_hash,
        coordinate_order_sha256 = coordinate_order_sha256,
        lattice_transform_sha256 = digest::digest(
          "lattice", "sha256", serialize = FALSE),
        mechanism_plan = mechanism_plan,
        plan_sha256 = contract$plan_sha256,
        epsilon = contract$epsilon, allocated_delta = contract$delta,
        sensitivity_steps = contract$sensitivity_steps,
        commitment_context = digest::digest(
          paste0("context-", peer), "sha256", serialize = FALSE),
        seed_commitment = digest::digest(
          paste0("seed-", peer), "sha256", serialize = FALSE),
        complete_epsilon_per_peer = profile$complete_epsilon_per_peer,
        delta_aggregation = profile$delta_aggregation,
        capability_available = TRUE), if (isTRUE(profile$selection_bound)) {
          list(
            backend_assessment = backend_assessment,
            backend_selection = backend_selection,
            backend_selection_sha256 =
              backend_selection$selection_sha256,
            one_joint_draw = isTRUE(profile$exact_gc))
        } else list()), peer)
  }), designated)

  implementation_delta_fields <- if (isTRUE(profile$exact_gc)) {
    c("implementation_delta_numerator", "implementation_delta_denominator")
  } else {
    c("per_peer_implementation_delta_numerator",
      "per_peer_implementation_delta_denominator")
  }
  implementation_delta_values <- mechanism_plan[
    implementation_delta_fields]

  start <- stats::setNames(lapply(designated, function(peer) {
    circuit <- digest::digest(
      "fixture-vector-circuit", "sha256", serialize = FALSE)
    identity <- list(
      domain = "dsVert/joint-dp/vector/exact-gc-operation/v1",
      manifest_sha256 = manifest_bundle$manifest_sha256,
      release_contract_hash = contract$release_contract_hash,
      selection_sha256 = contract$backend_selection_sha256,
      transcript_hash = contract$transcript_hash,
      chunk_index = 0L, coordinate_count = as.integer(coordinate_count),
      circuit_digest = circuit,
      purpose = paste0("joint-dp-vector-laplace-v3/", circuit))
    suffix <- substr(
      .dsvert_joint_dp_vector_exact_gc_client_hash(identity), 1L, 32L)
    binding_unsigned <- c(list(
      version = "dsvert-joint-dp-vector-exact-gc-binding-v1"),
      identity, list(
        operation_id = paste0("op_", suffix),
        source_key = paste0("exact_gc_in_", suffix),
        output_key = paste0("exact_gc_out_", suffix),
        operation = "joint-dp-vector-laplace-v3",
        output_kind = "joint-dp-vector-ring128-share-v1",
        source_producer = "joint.dp.vector.one-draw.v1"))
    binding <- c(binding_unsigned, list(
      binding_sha256 =
        .dsvert_joint_dp_vector_exact_gc_client_hash(binding_unsigned)))
    sign_receipt(c(common(
      peer, .DSVERT_CLIENT_VECTOR_START_VERSION, "vector_chunk_noised"),
      list(
        chunk_index = 0, coordinate_offset = 0,
        coordinates_in_chunk = as.numeric(coordinate_count),
        implementation_delta_numerator =
          implementation_delta_values[[1L]],
        implementation_delta_denominator =
          implementation_delta_values[[2L]],
        intermediate_payload_exposed = FALSE,
        capability_available = TRUE), if (isTRUE(profile$exact_gc)) list(
          backend_selection = backend_selection,
          backend_selection_sha256 = backend_selection$selection_sha256,
          binding = binding, binding_sha256 = binding$binding_sha256,
          operation_id = binding$operation_id, purpose = binding$purpose,
          initialization = list(state = "running"),
          source_share_exposed = FALSE, private_seed_exposed = FALSE,
          preclamp_values_exposed = FALSE) else list(
          noised_share_sha256 = digest::digest(
            paste0("share-", peer), "sha256", serialize = FALSE),
          sampler_contract_hash = digest::digest(
            paste0("sampler-", peer), "sha256", serialize = FALSE))), peer)
  }), designated)

  result_values <- stats::setNames(lapply(designated, function(peer) {
    commitment <- digest::digest(
      paste0("share-", peer), "sha256", serialize = FALSE)
    unsigned <- c(common(
      peer, .DSVERT_CLIENT_VECTOR_RESULT_VERSION,
      "vector_local_result_committed"), list(
        local_chunk_commitments = list(commitment),
        local_chunk_set_root = .dsvert_vector_merkle_root(commitment),
        local_chunk_set_sha256 = .dsvert_vector_hash(list(
          protocol = "dsvert-joint-dp-vector-local-chunk-set-v3",
          peer_name = peer, commitments = list(commitment))),
        all_chunks_durable = TRUE, intermediate_payload_exposed = FALSE,
        capability_available = TRUE))
    json <- sign_receipt(unsigned, peer)
    list(json = json, value = .dsvert_joint_dp_client_decode(
      json, "test result", 2L * 1024L^2))
  }), designated)
  result_json <- lapply(result_values, `[[`, "json")
  result_receipts <- lapply(result_values, `[[`, "value")
  result_set_hash <- .dsvert_vector_hash(list(
    protocol = "dsvert-joint-dp-vector-result-set-v3",
    ordered_results = unname(result_receipts)))

  scale <- 2^20
  scaled_values <- as.list(as.character(
    seq_len(coordinate_count) * scale))
  public_chunk <- list(
    version = .DSVERT_CLIENT_VECTOR_CHUNK_VERSION,
    capsule_id = capsule_id,
    release_contract_hash = contract$release_contract_hash,
    chunk_index = 0, chunk_count = 1,
    coordinate_offset = 0,
    coordinates_in_chunk = as.numeric(coordinate_count),
    output_lattice_bits = 20, output_lattice_scale = scale,
    scaled_values = scaled_values,
    value_encoding = "nonnegative-decimal-integer-common-lattice-v1",
    postprocessing = profile$postprocessing,
    source_values_exposed = FALSE, preclamp_values_exposed = FALSE)
  chunk_hash <- .dsvert_vector_hash(public_chunk)
  root <- .dsvert_vector_merkle_root(chunk_hash)
  release <- stats::setNames(lapply(designated, function(peer) {
    sign_receipt(c(common(
      peer, .DSVERT_CLIENT_VECTOR_RELEASE_VERSION, "vector_released"), list(
        result_set_hash = result_set_hash, final_vector_root = root,
        final_chunk_hashes = list(chunk_hash), output_lattice_bits = 20,
        output_lattice_scale = scale,
        mechanism = profile$release_mechanism,
        epsilon = contract$epsilon, delta = contract$delta,
        implementation_delta_numerator =
          implementation_delta_values[[1L]],
        implementation_delta_denominator =
          implementation_delta_values[[2L]],
        delta_aggregation = profile$delta_aggregation,
        postprocessing = profile$postprocessing,
        intermediate_payload_exposed = FALSE, durable_replay = TRUE,
        capability_available = TRUE)), peer)
  }), designated)
  replay <- .dsvert_joint_dp_client_json(list(
    version = "dsvert-joint-dp-vector-replay-v4", capsule_id = capsule_id,
    release_contract_hash = contract$release_contract_hash,
    result_set_hash = result_set_hash, final_vector_root = root,
    chunk_hash = chunk_hash, chunk = public_chunk, merkle_proof = list(),
    durable_replay = TRUE, source_store_read = FALSE,
    sampler_invoked = FALSE, history_gate = TRUE, request_limit = FALSE,
    operation_limit = TRUE))
  ack <- stats::setNames(lapply(peers, function(peer) sign_receipt(c(common(
    peer, .DSVERT_CLIENT_VECTOR_ACK_VERSION,
    "vector_finalized_and_compacted"), list(
      final_vector_root = root, source_intermediates_compacted = TRUE,
      sampler_intermediates_compacted = TRUE, final_chunks_retained = TRUE,
      durable_replay_retained = TRUE, idempotent = TRUE)), peer)), peers)

  make_transfer <- function(sender) {
    recipient <- setdiff(designated, sender)
    ciphertext <- .vector_client_b64url(charToRaw(paste0(
      "opaque-vector-ciphertext-", sender)))
    transfer <- list(
      ticket = .vector_client_b64url(charToRaw(paste0("ticket-", sender))),
      transfer_id = paste0("tb_", substr(digest::digest(
        paste0("transfer-", sender), "sha256", serialize = FALSE), 1L, 32L)),
      capability_id = .DSVERT_CLIENT_VECTOR_TYPED_CAPABILITY,
      sender_name = sender, recipient_name = recipient,
      payload_chars = as.numeric(nchar(ciphertext, type = "bytes")),
      payload_sha256 = digest::digest(
        charToRaw(ciphertext), "sha256", serialize = FALSE))
    list(
      ciphertext = ciphertext, transfer = transfer, capsule_id = capsule_id,
      release_contract_hash = contract$release_contract_hash,
      result_set_hash = result_set_hash, chunk_index = 0,
      intermediate_payload_exposed = FALSE, capability_available = TRUE)
  }

  state <- new.env(parent = emptyenv())
  state$allocation_ready <- FALSE
  state$result_ready <- FALSE
  state$release_ready <- FALSE
  state$source_calls <- 0L
  state$setup_calls <- 0L
  state$store_calls <- 0L
  state$cleanup_calls <- 0L
  state$aggregate_calls <- list()
  aggregate <- function(conns, expr, error, async = TRUE,
                        errors.print = FALSE, ...) {
    calls <- if (is.list(expr) && !is.call(expr)) expr else
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    name <- as.character(calls[[1L]][[1L]])
    state$aggregate_calls[[length(state$aggregate_calls) + 1L]] <- list(
      name = name, targets = names(conns), expressions = calls)
    values <- switch(name,
      dsvertJointDPVectorAllocationProofDS = {
        if (!isTRUE(state$allocation_ready)) {
          for (peer in names(conns)) {
            error(peer, "[dsvert_phase_not_ready:v1]")
          }
          stats::setNames(rep(list(NULL), length(conns)), names(conns))
        } else allocation_opening[names(conns)]
      },
      dsvertJointDPVectorAllocationPrepareDS =
        allocation_prepare[names(conns)],
      dsvertJointDPVectorAllocationCommitDS =
        allocation_commit[names(conns)],
      dsvertJointDPVectorAllocationAuthorizeDS =
        allocation_authorization[names(conns)],
      dsvertJointDPVectorAllocationOpenDS = {
        state$allocation_ready <- TRUE
        allocation_opening[names(conns)]
      },
      dsvertJointDPVectorPrepareDS = prepare[names(conns)],
      dsvertJointDPVectorStartDS = {
        if (!isTRUE(profile$exact_gc)) state$result_ready <- TRUE
        start[names(conns)]
      },
      dsvertJointDPVectorResultDS = {
        if (!isTRUE(state$result_ready)) {
          for (peer in names(conns)) {
            error(peer, "[dsvert_phase_not_ready:v1]")
          }
          stats::setNames(rep(list(NULL), length(conns)), names(conns))
        } else result_json[names(conns)]
      },
      dsvertJointDPVectorFinalShareDS = stats::setNames(
        lapply(names(conns), make_transfer), names(conns)),
      dsvertJointDPVectorReleaseDS = {
        if (!isTRUE(state$release_ready)) {
          for (peer in names(conns)) {
            error(peer, "[dsvert_phase_not_ready:v1]")
          }
          stats::setNames(rep(list(NULL), length(conns)), names(conns))
        } else release[names(conns)]
      },
      dsvertJointDPVectorReplayDS = stats::setNames(
        rep(list(replay), length(conns)), names(conns)),
      dsvertJointDPVectorFinalizeAckDS = ack[names(conns)],
      mpcCleanupDS = {
        state$cleanup_calls <- state$cleanup_calls + 1L
        stats::setNames(rep(list(TRUE), length(conns)), names(conns))
      },
      stop("unexpected client vector call: ", name))
    values
  }
  source_transport <- function(manifest_json, context, aggregate) {
    state$source_calls <- state$source_calls + 1L
    expect_identical(names(context$allocation_openings), designated)
    expect_identical(context$allocation_openings, allocation_opening)
    list(
      capsule_id = capsule_id,
      contract_hash = contract$source_contract_hash,
      coordinate_count = as.numeric(coordinate_count),
      chunk_count = as.numeric(contract$chunk_count),
      sampler_handoff_ready = TRUE, payload_exposed = FALSE)
  }
  setup_transport <- function(datasources, server_names, servers, session_id,
                              .aggregate) {
    state$setup_calls <- state$setup_calls + 1L
    expect_identical(servers, peers)
    expect_identical(names(datasources), peers)
    invisible(stats::setNames(rep("pk", length(peers)), peers))
  }
  setup_exact <- function(datasources, server_names, servers, session_id,
                          .aggregate) {
    state$setup_calls <- state$setup_calls + 1L
    expect_identical(as.integer(servers), match(designated, peers))
    expect_identical(names(datasources), peers)
    invisible(stats::setNames(rep("pk", length(designated)), designated))
  }
  run_exact <- function(
      datasources, server_names, servers, session_id, binding,
      manifest_sha256, release_contract_hash, selection_sha256,
      transcript_hash, chunk_index, coordinate_count, initialized,
      transport_ready, .aggregate) {
    expect_identical(as.integer(servers), match(designated, peers))
    expect_identical(manifest_sha256, manifest_bundle$manifest_sha256)
    expect_identical(release_contract_hash,
                     contract$release_contract_hash)
    expect_identical(selection_sha256,
                     backend_selection$selection_sha256)
    expect_true(isTRUE(transport_ready))
    expect_length(initialized, 2L)
    state$result_ready <- TRUE
    invisible(list(complete = TRUE))
  }
  store_typed <- function(blob, transfer, conn, session_id, producer_conn,
                          .aggregate) {
    state$store_calls <- state$store_calls + 1L
    expect_identical(names(producer_conn), transfer$sender_name)
    expect_identical(names(conn), transfer$recipient_name)
    if (state$store_calls >= 2L) state$release_ready <- TRUE
    invisible(1L)
  }
  list(
    context = context, manifest = manifest,
    manifest_bundle = manifest_bundle, release_instance = release_instance,
    contract = contract,
    state = state, keys = keys,
    aggregate = aggregate, source_transport = source_transport,
    setup_transport = setup_transport, setup_exact = setup_exact,
    run_exact = run_exact, store_typed = store_typed,
    allocation_opening = allocation_opening,
    start = start, result_set_hash = result_set_hash,
    release = release, replay = replay, sign_receipt = sign_receipt)
}

test_that("K=3 through K=5 vector orchestration releases once then replays", {
  for (k in 3:5) {
    fixture <- .vector_client_fixture(k = k)
    run <- function() testthat::with_mocked_bindings(
      .dsvert_joint_dp_vector_capsule(
        fixture$context$all_conns,
        manifest_bundle = fixture$manifest_bundle,
        .aggregate = fixture$aggregate,
        .source_transport = fixture$source_transport,
        .setup_transport = fixture$setup_transport,
        .setup_exact = fixture$setup_exact,
        .run_exact = fixture$run_exact,
        .store_typed = fixture$store_typed),
      .dsvert_dp_capsule_source_manifest = function(manifest_json, context) {
        list(value = fixture$manifest,
             capsule_id = fixture$manifest$capsule_identity$capsule_id)
      },
      .package = "dsVertClient")

    first <- run()
    expect_s3_class(first, "dsvert_joint_dp_vector")
    expect_equal(first$values, c(1, 2, 3), tolerance = 0)
    expect_identical(
      first$coordinate_order_sha256,
      .dsvert_dp_capsule_vector_layout(first$manifest)$sha256)
    expect_identical(fixture$state$source_calls, 1L)
    expect_identical(fixture$state$setup_calls, 1L)
    expect_identical(fixture$state$store_calls, 2L)
    first_methods <- vapply(
      fixture$state$aggregate_calls, `[[`, character(1L), "name")
    expect_identical(sum(
      first_methods == "dsvertJointDPVectorAllocationProofDS"), 1L)
    expect_identical(sum(
      first_methods == "dsvertJointDPVectorResultDS"), 2L)
    expect_identical(sum(
      first_methods == "dsvertJointDPVectorReleaseDS"), 2L)
    expect_identical(first$history_gate, TRUE)
    expect_identical(first$request_limit, FALSE)
    expect_identical(first$operation_limit, TRUE)
    expect_false(any(grepl(
      "seed|share|ciphertext|prepare|receipt|scaled", names(first),
      ignore.case = TRUE)))

    source_before <- fixture$state$source_calls
    setup_before <- fixture$state$setup_calls
    store_before <- fixture$state$store_calls
    second <- run()
    expect_identical(second$values, first$values)
    expect_identical(second$final_vector_root, first$final_vector_root)
    expect_identical(fixture$state$source_calls, source_before)
    expect_identical(fixture$state$setup_calls, setup_before)
    expect_identical(fixture$state$store_calls, store_before)
  }
})

test_that("vector phase lookup catches only phase-not-ready", {
  phase <- .dsvert_client_parse_phase_not_ready(
    "[dsvert_phase_not_ready:v1]")
  expect_identical(
    .dsvert_vector_try_phase(stop(phase)),
    list(ok = FALSE, value = NULL))

  propagated <- list(
    poison = .dsvert_dsi_poisoned_session_condition("site_a"),
    lifetime = .dsvert_client_parse_dp_lifetime_budget_exhausted(
      "[dsvert_dp_lifetime_budget_exhausted:v1]"),
    retry = .dsvert_client_parse_release_instance_retry(
      "[dsvert_retry_current_release_instance:new_release_instance]"),
    oversize = .dsvert_client_resource_oversize(),
    backpressure = .dsvert_client_parse_resource_backpressure(
      "[dsvert_resource_backpressure:v1]"),
    deadline = .dsvert_retry_deadline_condition("vector phase lookup"),
    generic = simpleError("generic vector lookup failure"))
  for (kind in names(propagated)) {
    condition <- propagated[[kind]]
    observed <- tryCatch(
      .dsvert_vector_try_phase(stop(condition)), error = identity)
    expect_identical(observed, condition, info = kind)
  }
})

test_that("invalid allocation blocks source access and the exact cross gate", {
  fixture <- .vector_client_fixture(k = 2L)
  touched <- new.env(parent = emptyenv())
  touched$allocation <- 0L
  touched$source <- 0L
  touched$aggregate <- 0L
  touched$cross <- 0L
  condition <- testthat::with_mocked_bindings(
    tryCatch(.dsvert_joint_dp_vector_capsule_once(
      fixture$context$all_conns,
      manifest_bundle = fixture$manifest_bundle,
      .aggregate = function(...) {
        touched$aggregate <- touched$aggregate + 1L
        stop("DSI must not be reached", call. = FALSE)
      },
      .allocation = function(...) {
        touched$allocation <- touched$allocation + 1L
        # Deliberately lacks the two cross-signed openings and certificate.
        fixture$context
      },
      .source_transport = function(...) {
        touched$source <- touched$source + 1L
        stop("source material must not be reached", call. = FALSE)
      },
      .cross_orchestrate = function(...) {
        touched$cross <- touched$cross + 1L
        stop("the exact cross gate must not be reached", call. = FALSE)
      }), error = identity),
    .dsvert_dp_capsule_source_manifest = function(manifest_json, context) {
      list(value = fixture$manifest,
           capsule_id = fixture$manifest$capsule_identity$capsule_id)
    },
    .dsvert_dp_gaussian_cross_artifacts_client = function(manifest) {
      list(cross_model = list())
    },
    .dsvert_dp_categorical_cross_artifacts_client = function(manifest) list(),
    .package = "dsVertClient")

  expect_s3_class(condition, "error")
  expect_match(condition$message, "verified cross-signed allocation")
  expect_identical(touched$allocation, 1L)
  expect_identical(touched$source, 0L)
  expect_identical(touched$aggregate, 0L)
  expect_identical(touched$cross, 0L)
})

test_that("unilateral rotation changes only the pre-claim release candidate", {
  first <- .vector_client_fixture(k = 2L, epochs = c(1, 2))
  rotated <- .vector_client_fixture(
    k = 2L, epochs = c(1, 3), keys = first$keys)
  first_candidate <- .dsvert_vector_release_instance(
    first$context, first$manifest$capsule_identity$capsule_id)
  replay_candidate <- .dsvert_vector_release_instance(
    first$context, first$manifest$capsule_identity$capsule_id)
  rotated_candidate <- .dsvert_vector_release_instance(
    rotated$context, rotated$manifest$capsule_identity$capsule_id)

  expect_identical(replay_candidate, first_candidate)
  expect_identical(rotated_candidate$value$capsule_id,
                   first_candidate$value$capsule_id)
  expect_false(identical(rotated_candidate$id, first_candidate$id))
  expect_equal(
    rotated_candidate$value$peer_noise_roots$site_b$privacy_epoch,
    3, tolerance = 0)
  expect_identical(first$state$source_calls, 0L)
  expect_identical(rotated$state$source_calls, 0L)
})

test_that("release-domain rotation changes the pre-claim candidate", {
  fixture <- .vector_client_fixture(k = 3L, epochs = c(1, 1, 1))
  first <- .dsvert_vector_release_instance(
    fixture$context, fixture$manifest$capsule_identity$capsule_id)
  rotated <- fixture$context
  peer <- rotated$designated[[1L]]
  rotated$status[[peer]]$release_domain$generation <- 2
  rotated$status[[peer]]$release_domain$rotation_count <- 1
  rotated$status[[peer]]$release_domain$domain_id <- paste0(
    "rd_", digest::digest("domain-only-rotation", "sha256",
                           serialize = FALSE))
  second <- .dsvert_vector_release_instance(
    rotated, fixture$manifest$capsule_identity$capsule_id)

  expect_identical(first$value$capsule_id, second$value$capsule_id)
  expect_false(identical(first$id, second$id))
  expect_identical(
    first$value$peer_noise_roots[[peer]]$noise_key_id,
    second$value$peer_noise_roots[[peer]]$noise_key_id)
  expect_identical(
    second$value$peer_noise_roots[[peer]]$release_domain_generation, 2)
})

test_that("a pre-claim stale prepare refreshes roots once and completes", {
  stale <- .vector_client_fixture(k = 2L, epochs = c(1, 1))
  current <- .vector_client_fixture(
    k = 2L, epochs = c(1, 2), keys = stale$keys)
  active <- stale
  injected <- FALSE
  refreshes <- 0L
  successful_instances <- character()
  aggregate <- function(conns, expr, error, async = TRUE,
                        errors.print = FALSE, ...) {
    calls <- if (is.list(expr) && !is.call(expr)) expr else
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    name <- as.character(calls[[1L]][[1L]])
    if (identical(name, "dsvertJointDPVectorPrepareDS") &&
        identical(active, stale) && !isTRUE(injected)) {
      injected <<- TRUE
      for (peer in names(conns)) error(
        peer, paste0("[",
          .DSVERT_CLIENT_VECTOR_RETRY_CURRENT_INSTANCE_TOKEN,
          ":new_release_instance] rotated root"))
      return(stats::setNames(rep(list(NULL), length(conns)), names(conns)))
    }
    value <- active$aggregate(
      conns, expr, error, async = async, errors.print = errors.print, ...)
    if (identical(name, "dsvertJointDPVectorReleaseDS") &&
        isTRUE(active$state$release_ready) &&
        all(vapply(value, is.character, logical(1L)))) {
      successful_instances <<- c(
        successful_instances, active$release_instance$id)
    }
    value
  }
  source_transport <- function(...) active$source_transport(...)
  setup_transport <- function(...) active$setup_transport(...)
  setup_exact <- function(...) active$setup_exact(...)
  run_exact <- function(...) active$run_exact(...)
  store_typed <- function(...) active$store_typed(...)
  refresh <- function(datasources, aggregate) {
    refreshes <<- refreshes + 1L
    active <<- current
    list(status = NULL, manifest_bundle = current$manifest_bundle)
  }
  run <- function() testthat::with_mocked_bindings(
    .dsvert_joint_dp_vector_capsule(
      active$context$all_conns,
      manifest_bundle = active$manifest_bundle,
      .aggregate = aggregate,
      .source_transport = source_transport,
      .setup_transport = setup_transport,
      .setup_exact = setup_exact,
      .run_exact = run_exact,
      .store_typed = store_typed,
      .retry_refresh = refresh),
    .dsvert_dp_capsule_source_manifest = function(manifest_json, context) {
      list(value = active$manifest,
           capsule_id = active$manifest$capsule_identity$capsule_id)
    },
    .package = "dsVertClient")

  released <- run()
  replay <- run()
  expect_true(injected)
  expect_identical(refreshes, 1L)
  expect_identical(released$release_instance_id,
                   current$release_instance$id)
  expect_identical(replay$release_instance_id,
                   released$release_instance_id)
  expect_identical(replay$final_vector_root,
                   released$final_vector_root)
  expect_length(unique(successful_instances), 1L)
})

test_that("a post-claim lifetime refusal never refreshes or re-enters source", {
  fixture <- .vector_client_fixture(k = 2L)
  refreshes <- 0L
  release_calls <- 0L
  aggregate <- function(conns, expr, error, async = TRUE,
                        errors.print = FALSE, ...) {
    calls <- if (is.list(expr) && !is.call(expr)) expr else
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    name <- as.character(calls[[1L]][[1L]])
    if (identical(name, "dsvertJointDPVectorReleaseDS")) {
      release_calls <<- release_calls + 1L
      for (peer in names(conns)) {
        error(peer, "[dsvert_dp_lifetime_budget_exhausted:v1]")
      }
      return(stats::setNames(rep(list(NULL), length(conns)), names(conns)))
    }
    fixture$aggregate(
      conns, expr, error, async = async, errors.print = errors.print, ...)
  }
  condition <- testthat::with_mocked_bindings(
    tryCatch(.dsvert_joint_dp_vector_capsule(
      fixture$context$all_conns,
      manifest_bundle = fixture$manifest_bundle,
      .aggregate = aggregate,
      .source_transport = fixture$source_transport,
      .setup_transport = fixture$setup_transport,
      .setup_exact = fixture$setup_exact,
      .run_exact = fixture$run_exact,
      .store_typed = fixture$store_typed,
      .retry_refresh = function(...) {
        refreshes <<- refreshes + 1L
        stop("refresh must not run", call. = FALSE)
      }), error = identity),
    .dsvert_dp_capsule_source_manifest = function(manifest_json, context) {
      list(value = fixture$manifest,
           capsule_id = fixture$manifest$capsule_identity$capsule_id)
    },
    .package = "dsVertClient")

  expect_s3_class(condition, "dsvert_dp_lifetime_budget_exhausted")
  expect_identical(refreshes, 0L)
  expect_identical(release_calls, 1L)
  expect_identical(fixture$state$source_calls, 1L)
})

test_that("duplicate designated noise roots fail with an actionable typed error", {
  check_duplicate <- function(k, same_provider) {
    fixture <- .vector_client_fixture(k = k)
    context <- fixture$context
    first <- context$designated[[1L]]
    second <- context$designated[[2L]]
    context$status[[second]]$noise_root$key_id <-
      context$status[[first]]$noise_root$key_id
    if (isTRUE(same_provider)) {
      context$status[[second]]$noise_root$provider_id <-
        context$status[[first]]$noise_root$provider_id
    }
    condition <- tryCatch(
      .dsvert_vector_release_instance(
        context, fixture$manifest$capsule_identity$capsule_id),
      error = identity)
    expect_s3_class(condition, "dsvert_noise_root_not_independent")
    expect_identical(condition$code, "noise_root_not_independent")
    expect_match(condition$message, "regenerate only")
  }
  check_duplicate(2L, same_provider = TRUE)
  check_duplicate(2L, same_provider = FALSE)
  check_duplicate(3L, same_provider = FALSE)

  fixture <- .vector_client_fixture(k = 3L)
  context <- fixture$context
  context$status$site_c$noise_root$key_id <-
    context$status$site_a$noise_root$key_id
  expect_silent(.dsvert_vector_release_instance(
    context, fixture$manifest$capsule_identity$capsule_id))
})

test_that("fixed-work Gaussian v2 vector profile is dispatched end to end", {
  fixture <- .vector_client_fixture(k = 2L, gaussian = TRUE)
  result <- testthat::with_mocked_bindings(
    .dsvert_joint_dp_vector_capsule(
      fixture$context$all_conns,
      manifest_bundle = fixture$manifest_bundle,
      .aggregate = fixture$aggregate,
      .source_transport = fixture$source_transport,
      .setup_transport = fixture$setup_transport,
      .setup_exact = fixture$setup_exact,
      .run_exact = fixture$run_exact,
      .store_typed = fixture$store_typed),
    .dsvert_dp_capsule_source_manifest = function(manifest_json, context) {
      list(value = fixture$manifest,
           capsule_id = fixture$manifest$capsule_identity$capsule_id)
    },
    .package = "dsVertClient")
  expect_identical(result$mechanism,
                   .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
  prepare <- result$signed_provenance$prepare_receipts[[1L]]
  expect_identical(prepare$backend,
                   .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND)
  expect_identical(prepare$sampler,
                   .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER)
  expect_identical(prepare$mechanism_plan$mechanism,
                   .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM)
  expect_identical(.dsvert_vector_hash(prepare$mechanism_plan),
                   prepare$plan_sha256)
  expect_match(prepare$sensitivity_steps, "[.eE]")
})

test_that("Gaussian v2 threat-model and transcript fields fail closed", {
  fixture <- .vector_client_fixture(k = 2L, gaussian = TRUE)
  capsule <- fixture$manifest$workload$capsule_mechanism
  selection <- capsule$certificate
  request <- selection$gaussian_calibration_request
  plan <- selection$gaussian_plan
  validate <- function(candidate) {
    changed <- capsule
    changed$certificate$gaussian_plan <- candidate
    changed$certificate$gaussian_plan_sha256 <-
      .dsvert_vector_hash(candidate)
    profile <- .dsvert_vector_profile(changed)
    .dsvert_vector_plan_validate(
      candidate, changed$certificate$gaussian_plan_sha256, profile,
      request$total_coordinate_count, request$l2_sensitivity_steps)
  }
  expect_silent(validate(plan))
  expect_setequal(names(plan), .dsvert_vector_gaussian_plan_fields())

  semantic_tamper <- list(
    sampler_full_scan_steps = plan$sampler_full_scan_steps - 1L,
    sampler_cdf_table_bytes = plan$sampler_cdf_table_bytes + 1L,
    nominal_variance_multiplier = 1L,
    nominal_standard_deviation_factor = "unknown",
    at_least_one_honest_noise_peer = FALSE,
    maximum_colluding_noise_peers = 2L,
    adversary_view = "analyst_only",
    adversary_view_privacy_argument = "unsupported",
    source_share_hiding_precondition = "none",
    sampler_branches_on_private_randomness = TRUE,
    transcript_dp_claim = FALSE,
    logical_transcript_fixed_shape = FALSE,
    physical_timing_dp_claim = TRUE)
  for (field in names(semantic_tamper)) {
    changed <- plan
    changed[[field]] <- semantic_tamper[[field]]
    expect_error(validate(changed), "misbound", fixed = TRUE,
                 info = paste("tampered client Gaussian field", field))
  }
  missing <- plan
  missing$physical_timing_dp_claim <- NULL
  expect_error(validate(missing), "misbound", fixed = TRUE)
  extra <- plan
  extra$future_unreviewed_claim <- TRUE
  expect_error(validate(extra), "misbound", fixed = TRUE)
})

test_that("fresh server Go Gaussian plan validates in the client", {
  request <- list(
    epsilon = "1", delta = "0.000001", l2_sensitivity_steps = "3",
    total_coordinate_count = 3L)
  plan <- .client_fresh_go_gaussian_plan_v2(request)
  selection <- list(
    gaussian_calibration_request = request,
    gaussian_plan = plan,
    gaussian_plan_sha256 = .dsvert_vector_hash(plan))
  capsule <- list(
    mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
    sensitivity_norm = "l2", certificate = selection)
  profile <- .dsvert_vector_profile(capsule)
  expect_silent(.dsvert_vector_plan_validate(
    plan, selection$gaussian_plan_sha256, profile,
    request$total_coordinate_count, request$l2_sensitivity_steps))
  expect_setequal(names(plan), .dsvert_vector_gaussian_plan_fields())
  expect_false(plan$sampler_branches_on_private_randomness)
  expect_true(plan$logical_transcript_fixed_shape)
  expect_true(plan$transcript_dp_claim)
  expect_false(plan$physical_timing_dp_claim)
})

test_that("fresh server Go exact-Laplace plans validate after JSON roundtrip", {
  profile <- .dsvert_vector_profile(
    "discrete-laplace", backend = .DSVERT_CLIENT_VECTOR_EXACT_BACKEND)
  validate <- function(request) {
    plan <- .client_fresh_go_laplace_plan_v3(request)
    plan_sha256 <- .dsvert_vector_hash(plan)
    encoded <- .dsvert_joint_dp_client_json(plan)
    roundtrip <- .dsvert_joint_dp_client_decode(
      encoded, "exact-Laplace plan",
      .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES)
    expect_identical(.dsvert_vector_hash(roundtrip), plan_sha256)
    expect_silent(.dsvert_vector_plan_validate(
      roundtrip, plan_sha256, profile,
      request$total_coordinate_count, request$sensitivity_steps))
    roundtrip
  }

  default <- validate(list(
    epsilon = "4.0000000000000000e+00",
    delta = format(2^-100, digits = 17L, scientific = TRUE, trim = TRUE),
    sensitivity_steps = "65536", total_coordinate_count = 1L))
  expect_gt(nchar(default$implementation_delta_denominator,
                  type = "bytes"), 128L)
  expect_gt(nchar(default$implementation_delta_bound,
                  type = "bytes"), 256L)

  extreme <- validate(list(
    epsilon = "8.0000000000000000e+00",
    delta = format(2^-100, digits = 17L, scientific = TRUE, trim = TRUE),
    sensitivity_steps = "9007199254740991",
    total_coordinate_count = 1000000L))
  expect_gt(nchar(extreme$implementation_delta_denominator,
                  type = "bytes"),
            nchar(default$implementation_delta_denominator,
                  type = "bytes"))

  receipt_plan <- validate(list(
    epsilon = "1e+00",
    delta = format(2^-100, digits = 17L, scientific = TRUE, trim = TRUE),
    sensitivity_steps = "1048576", total_coordinate_count = 1L))
  fixture <- .vector_client_fixture(
    k = 2L, coordinate_count = 1L, laplace_plan = receipt_plan)
  expected_delta <- paste0(
    receipt_plan$implementation_delta_numerator, "/",
    receipt_plan$implementation_delta_denominator)
  expect_identical(
    .dsvert_vector_implementation_delta(fixture$contract), expected_delta)
  expect_silent(.dsvert_vector_start_set(
    fixture$start, fixture$context, fixture$contract, 0L))
  expect_silent(.dsvert_vector_release_set(
    fixture$release, fixture$context, fixture$contract,
    fixture$result_set_hash, expected_delta))

  oversized <- default
  oversized$implementation_delta_denominator <- paste0(
    "1", strrep("0", .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES))
  oversized_sha256 <- .dsvert_vector_hash(oversized)
  expect_error(.dsvert_vector_plan_validate(
    oversized, oversized_sha256, profile, 1L, "65536"),
    "exact-GC vector plan is invalid", fixed = TRUE)
  oversized_json <- .dsvert_joint_dp_client_json(oversized)
  expect_gt(nchar(oversized_json, type = "bytes"),
            .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES)
  expect_error(.dsvert_joint_dp_client_decode(
    oversized_json, "exact-Laplace plan",
    .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES),
    "Invalid joint-DP exact-Laplace plan", fixed = TRUE)
  oversized_contract <- fixture$contract
  oversized_contract$mechanism_plan$implementation_delta_denominator <-
    oversized$implementation_delta_denominator
  expect_error(.dsvert_vector_implementation_delta(oversized_contract),
               "implementation-delta certificate is invalid", fixed = TRUE)
})

test_that("vector client rejects signature and bilateral replay tampering", {
  fixture <- .vector_client_fixture(k = 2L)
  prepare <- fixture$aggregate(
    fixture$context$conns,
    stats::setNames(lapply(fixture$context$designated, function(peer) call(
      name = "dsvertJointDPVectorPrepareDS")), fixture$context$designated),
    error = function(...) NULL)
  tampered_manifest <- fixture$manifest
  tampered_manifest$workload$families$admitted_count$dataset <-
    "tampered-cohort"
  expect_error(.dsvert_vector_prepare_set(
    prepare, fixture$context, tampered_manifest,
    fixture$release_instance, fixture$manifest_bundle$manifest_sha256),
    "coordinate order")

  tampered_prepare <- stats::setNames(lapply(
    fixture$context$designated, function(peer) {
      value <- .dsvert_joint_dp_client_decode(
        prepare[[peer]], "test prepare", 2L * 1024L^2)
      unsigned <- value[setdiff(names(value), "signature")]
      unsigned$coordinate_order_sha256 <- paste(rep("0", 64L),
                                                 collapse = "")
      fixture$sign_receipt(unsigned, peer)
    }), fixture$context$designated)
  expect_error(.dsvert_vector_prepare_set(
    tampered_prepare, fixture$context, fixture$manifest,
    fixture$release_instance, fixture$manifest_bundle$manifest_sha256),
    "coordinate order")

  bad_prepare <- fixture$aggregate(
    fixture$context$conns,
    stats::setNames(lapply(fixture$context$designated, function(peer) call(
      name = "dsvertJointDPVectorPrepareDS")), fixture$context$designated),
    error = function(...) NULL)
  decoded <- .dsvert_joint_dp_client_decode(
    bad_prepare[[1L]], "test prepare", 2L * 1024L^2)
  first <- substr(decoded$signature, 1L, 1L)
  decoded$signature <- paste0(if (identical(first, "A")) "B" else "A",
                              substr(decoded$signature, 2L,
                                     nchar(decoded$signature)))
  bad_prepare[[1L]] <- .dsvert_joint_dp_client_json(decoded)
  expect_error(.dsvert_vector_prepare_set(
    bad_prepare, fixture$context, fixture$manifest,
    fixture$release_instance, fixture$manifest_bundle$manifest_sha256),
    "invalid vector signature")

  canonical_signature <- decoded$signature
  last <- substr(canonical_signature, nchar(canonical_signature),
                 nchar(canonical_signature))
  noncanonical_last <- c(A = "B", Q = "R", g = "h", w = "x")[[last]]
  expect_false(is.null(noncanonical_last))
  noncanonical_signature <- paste0(substr(
    canonical_signature, 1L, nchar(canonical_signature) - 1L),
    noncanonical_last)
  expect_error(.dsvert_joint_dp_client_b64url(
    noncanonical_signature, 64L, "signature"), "Invalid signature")

  prepared <- .dsvert_vector_prepare_set(
    fixture$aggregate(
      fixture$context$conns,
      stats::setNames(lapply(fixture$context$designated, function(peer) call(
        name = "dsvertJointDPVectorPrepareDS")), fixture$context$designated),
      error = function(...) NULL), fixture$context, fixture$manifest,
    fixture$release_instance, fixture$manifest_bundle$manifest_sha256)
  releases <- lapply(fixture$release, function(value) {
    .dsvert_joint_dp_client_decode(value, "test release", 2L * 1024L^2)
  })
  replay <- stats::setNames(rep(list(fixture$replay), 2L),
                            fixture$context$designated)
  replay[[2L]] <- paste0(replay[[2L]], " ")
  expect_error(.dsvert_vector_replay_chunk(
    replay, fixture$context, prepared$contract, releases[[1L]], 0L),
    "different vector replay")
})
