.dp_gaussian_artifact <- function(capacity = 100, scale = 256) {
  list(
    version = "bounded-normalized-gaussian-sufficient-statistics-v1",
    spec_version = "v1", analysis_id = "gaussian_primary",
    dataset = "cohort", owner_peer = "site_a",
    outcome = list(column = "y", lower = 10, upper = 20),
    predictors = list(
      x = list(column = "x", lower = -2, upper = 2)),
    predictor_order = "x", intercept = TRUE,
    design_terms = c("(Intercept)", "x"), numeric_grid_bits = 8L,
    coordinate_count = 7L,
    coordinate_order = paste0(
      "n_then_xtx_upper_column_major_then_xty_design_order_then_yty_",
      "v1"),
    repeated_record_policy = paste0(
      "clip_finite_rows_then_mean_each_variable_once_per_admitted_unit_",
      "v1"),
    missingness_policy =
      "complete_case_across_outcome_and_all_predictors_v1",
    contribution_domain =
      "one_vector_of_normalized_monomials_in_closed_unit_interval_v1",
    count_gram_intercept_policy = paste0(
      "n_is_complete_case_count_and_moment_upper_bound_gram11_governs_",
      "the_solve_no_averaging_v1"),
    statistic_maximum = c(capacity, rep(capacity * scale, 6L)),
    source_raw_l1_sensitivity = 1 + 6 * scale,
    source_raw_l2_sensitivity = sqrt(1 + 6 * scale^2) *
      (1 + 64 * .Machine$double.eps),
    natural_l1_sensitivity = 7,
    natural_l2_sensitivity = sqrt(7) *
      (1 + 64 * .Machine$double.eps),
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste0(
      "zero_missing_complete_case_vs_all_one_complete_unit_is_worst_",
      "case_for_add_remove_and_replace_one"),
    regularization_policy =
      "none_in_release_explicit_client_postprocessing_only_v1",
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
}

.dp_gaussian_laplace_plan <- function(coordinate_count,
                                      sensitivity_steps) {
  list(
    version =
      "dsvert-joint-dp-vector-independent-full-draw-convolution-plan-v3",
    sampler = .DSVERT_CLIENT_VECTOR_SAMPLER,
    stop_bits = 128,
    stop_numerator = "166112941255448185114116409827120804",
    one_geometric_tv_numerator =
      "6125082604576892348297572878053259906739",
    one_geometric_tv_denominator = paste0(
      "231584178474632390847141970017375815706539969331281128078915",
      "168015826259279872"),
    sensitivity_steps = sensitivity_steps,
    total_coordinate_count = as.numeric(coordinate_count),
    maximum_chunk_coordinates = as.numeric(min(
      .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES, coordinate_count)),
    independent_noise_peer_count = 2,
    complete_epsilon_per_peer = TRUE,
    epsilon_divided_by_peer_count = FALSE,
    capability_available = TRUE,
    release_implementation_delta_aggregation = "max_per_peer_not_sum",
    per_peer_implementation_delta_numerator = "1",
    per_peer_implementation_delta_denominator =
      "1267650600228229401496703205376")
}

.dp_gaussian_certificate_fixture <- function(
    manifest, layout, release, status, capacity, scale,
    peers = c("site_a", "site_b"), designated = peers[1:2], keys = NULL) {
  stopifnot(
    is.character(peers), length(peers) >= 2L, !anyNA(peers),
    !anyDuplicated(peers), identical(peers, sort(peers, method = "radix")),
    is.character(designated), length(designated) == 2L,
    !anyNA(designated), !anyDuplicated(designated),
    identical(designated, sort(designated, method = "radix")),
    all(designated %in% peers))
  if (is.null(keys)) {
    keys <- stats::setNames(lapply(peers, function(peer) {
      openssl::ed25519_keygen()
    }), peers)
  }
  stopifnot(is.list(keys), identical(names(keys), peers))
  b64url <- function(value) sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))),
    perl = TRUE)
  pins <- vapply(keys, function(key) {
    b64url(tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  pin_hash <- .dsvert_dp_pinset_hash(pins)
  capsule_id <- release$capsule_id
  logical_snapshot <- list(
    logical_snapshot_id = "cohort", version = "schema-v1-fixture",
    alignment_protocol_version = 1)
  manifest$logical_snapshot <- logical_snapshot
  manifest$capsule_schema <- "dsvert-biomedical-capsule-workload-v2"
  manifest$admission <- list(version = "fixture-admission-v1")
  manifest$bounds <- list(version = "fixture-bounds-v1")
  manifest$workload$capsule_mechanism <- list(
    mechanism = "discrete-laplace", sensitivity_norm = "l1")
  manifest$capsule_identity <- list(
    capsule_id = capsule_id,
    contract = list(
      consortium_id = "fixture-consortium",
      policy_contract_hash = digest::digest(
        "fixture-policy", "sha256", serialize = FALSE)))
  release$manifest <- manifest
  release$manifest_sha256 <- digest::digest(
    "fixture-manifest", "sha256", serialize = FALSE)
  policy <- list(
    domain = "fixture-domain", cohort_id = "cohort",
    peer_pinset_sha256 = pin_hash,
    designated_noise_peers = designated,
    capsule_epsilon = release$epsilon, capsule_delta = release$delta,
    adjacency = "add_remove_patient", unit_capacity = capacity,
    max_records_per_unit = 1, overflow_policy = "reject_snapshot")
  status <- stats::setNames(lapply(seq_along(peers), function(index) {
    peer <- peers[[index]]
    list(
    policy = policy,
    noise_root = list(
      privacy_epoch = as.numeric(index), epoch = as.numeric(index),
      provider_id = paste0("provider-", peer),
      key_id = paste0("root-", substr(peer, 6L, 6L))),
    release_domain = list(
      generation = as.numeric(index),
      domain_id = paste0("rd_", digest::digest(
        paste0("gaussian-fixture/", peer, "/", index),
        "sha256", serialize = FALSE))))
  }), peers)
  release_instance <- .dsvert_vector_release_instance(
    list(designated = designated, status = status), capsule_id)
  release$release_instance_id <- release_instance$id
  release$release_instance <- release_instance$value
  index <- .dsvert_dp_capsule_artifact_commitment_index_client(
    manifest, policy, release$manifest_sha256)

  sign_manifest <- function(unsigned, peer, domain) {
    message <- .dsvert_dp_capsule_manifest_message(domain, unsigned)
    c(unsigned, list(signature = b64url(
      openssl::ed25519_sign(message, keys[[peer]]))))
  }
  sign_vector <- function(unsigned, peer) {
    message <- charToRaw(paste0(
      .DSVERT_CLIENT_JOINT_DP_RECEIPT_DOMAIN,
      .dsvert_joint_dp_client_json(unsigned)))
    c(unsigned, list(signature = b64url(
      openssl::ed25519_sign(message, keys[[peer]]))))
  }
  schema_unsigned <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = logical_snapshot,
    peer_pinset_sha256 = pin_hash,
    datasets = list()))
  schema_signatures <- stats::setNames(lapply(peers, function(peer) {
    message <- charToRaw(paste0(
      .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
      .dsvert_joint_dp_client_json(schema_unsigned)))
    b64url(openssl::ed25519_sign(message, keys[[peer]]))
  }), peers)
  schema_json <- .dsvert_joint_dp_client_json(c(
    schema_unsigned, list(signatures = schema_signatures)))
  schema_sha256 <- .dsvert_dp_capsule_manifest_hash(schema_unsigned)
  workload_sha256 <- digest::digest(
    "fixture-workload", "sha256", serialize = FALSE)
  build_receipts <- stats::setNames(lapply(peers, function(peer) {
    sign_manifest(list(
      version = .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_BUILD_VERSION,
      phase = "server_authoritative_manifest_memoized",
      peer_name = peer, peer_identity_pk = unname(pins[[peer]]),
      peer_pinset_sha256 = pin_hash, schema_sha256 = schema_sha256,
      workload_contract_sha256 = workload_sha256,
      manifest_sha256 = release$manifest_sha256,
      manifest_bytes = 1024, capsule_id = capsule_id,
      privacy_epoch = status[[peer]]$noise_root$privacy_epoch,
      noise_key_id = status[[peer]]$noise_root$key_id,
      artifact_commitment_count = index$count,
      artifact_commitments_root = index$root,
      durable_memoization = TRUE, deterministic_replay = TRUE,
      data_access = FALSE, operation_limit = FALSE,
      request_limit = FALSE, history_can_deny_operation = FALSE),
      peer, "build")
  }), peers)

  release_contract_hash <- digest::digest(
    "fixture-release-contract", "sha256", serialize = FALSE)
  transcript_hash <- digest::digest(
    "fixture-transcript", "sha256", serialize = FALSE)
  source_contract_hash <- digest::digest(
    "fixture-source-contract", "sha256", serialize = FALSE)
  sensitivity_steps <- as.character(8 * scale)
  mechanism_plan <- .dp_gaussian_laplace_plan(
    release$coordinate_count, sensitivity_steps)
  plan_sha256 <- .dsvert_vector_hash(mechanism_plan)
  assessment_unsigned <- list(
    version = "dsvert-joint-dp-vector-exact-gc-assessment-v2",
    manifest_sha256 = release$manifest_sha256,
    representable = TRUE,
    exact_gc_capability_id = "joint_dp_biomedical_vector_exact_gc_v1",
    plan_sha256 = digest::digest(
      "fixture-exact-plan-request", "sha256", serialize = FALSE),
    maximum_chunk_coordinates = 128L,
    cost_policy_version =
      "dsvert-joint-dp-vector-exact-gc-cost-policy-v1",
    total_coordinate_count = as.integer(release$coordinate_count),
    maximum_promoted_coordinates = 1L, promoted = FALSE,
    selection_reason = "above_public_exact_gc_cost_ceiling",
    private_material_accessed = FALSE,
    runtime_failure_consulted = FALSE)
  backend_assessment <- c(assessment_unsigned, list(
    assessment_sha256 =
      .dsvert_joint_dp_vector_exact_gc_client_hash(assessment_unsigned)))
  selection_unsigned <- list(
    version = "dsvert-joint-dp-vector-backend-selection-v2",
    manifest_sha256 = release$manifest_sha256,
    backend = .DSVERT_CLIENT_VECTOR_BACKEND,
    one_draw = FALSE,
    cost_policy_version =
      "dsvert-joint-dp-vector-exact-gc-cost-policy-v1",
    total_coordinate_count = as.integer(release$coordinate_count),
    maximum_promoted_coordinates = 1L,
    selection_reason = "above_public_exact_gc_cost_ceiling",
    assessment_sha256 = backend_assessment$assessment_sha256,
    exact_gc_plan_sha256 = backend_assessment$plan_sha256,
    exact_gc_maximum_chunk_coordinates = 128L,
    selected_before_private_material = TRUE,
    retry_may_change_backend = FALSE)
  backend_selection <- c(selection_unsigned, list(
    selection_sha256 =
      .dsvert_joint_dp_vector_exact_gc_client_hash(selection_unsigned)))
  release$backend <- .DSVERT_CLIENT_VECTOR_BACKEND
  release$backend_assessment <- backend_assessment
  release$backend_selection <- backend_selection
  release$mechanism_plan <- mechanism_plan
  release$plan_sha256 <- plan_sha256
  chunk <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_VECTOR_CHUNK_VERSION,
    capsule_id = capsule_id,
    release_contract_hash = release_contract_hash,
    chunk_index = 0, chunk_count = 1,
    coordinate_offset = 0,
    coordinates_in_chunk = as.numeric(release$coordinate_count),
    output_lattice_bits = as.numeric(log2(scale)),
    output_lattice_scale = scale,
    scaled_values = as.list(.dsvert_dp_gaussian_scaled_text(
      release$values, scale)),
    value_encoding = "nonnegative-decimal-integer-common-lattice-v1",
    postprocessing = paste0(
      "signed-Ring128-decode-then-fixed-public-coordinate-clamp-v1"),
    source_values_exposed = FALSE, preclamp_values_exposed = FALSE))
  chunk_hash <- .dsvert_vector_hash(chunk)
  final_root <- .dsvert_vector_merkle_root(chunk_hash)
  release$final_vector_root <- final_root
  common <- function(peer, version, phase) list(
    version = version, phase = phase, capsule_id = capsule_id,
    release_instance_id = release_instance$id,
    release_instance = release_instance$value,
    release_contract_hash = release_contract_hash,
    transcript_hash = transcript_hash, peer_name = peer,
    peer_identity_pk = unname(pins[[peer]]),
    coordinate_count = as.numeric(release$coordinate_count),
    chunk_count = 1, backend = .DSVERT_CLIENT_VECTOR_BACKEND,
    sampler = .DSVERT_CLIENT_VECTOR_SAMPLER,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE)
  prepare <- stats::setNames(lapply(designated, function(peer) {
    sign_vector(c(common(
      peer, .DSVERT_CLIENT_VECTOR_PREPARE_VERSION, "vector_prepared"), list(
        source_contract_hash = source_contract_hash,
        coordinate_order_sha256 = layout$sha256,
        lattice_transform_sha256 = digest::digest(
          "fixture-lattice", "sha256", serialize = FALSE),
        mechanism_plan = mechanism_plan,
        plan_sha256 = plan_sha256,
        epsilon = "1e+00",
        allocated_delta = format(
          release$delta, digits = 17, scientific = TRUE, trim = TRUE),
        sensitivity_steps = sensitivity_steps,
        commitment_context = digest::digest(
          paste0("context-", peer), "sha256", serialize = FALSE),
        seed_commitment = digest::digest(
          paste0("seed-", peer), "sha256", serialize = FALSE),
        complete_epsilon_per_peer = TRUE,
        delta_aggregation = "max_per_peer_not_sum",
        capability_available = TRUE,
        backend_assessment = backend_assessment,
        backend_selection = backend_selection,
        backend_selection_sha256 = backend_selection$selection_sha256,
        one_joint_draw = FALSE)), peer)
  }), designated)
  result_set_hash <- digest::digest(
    "fixture-result-set", "sha256", serialize = FALSE)
  releases <- stats::setNames(lapply(designated, function(peer) {
    sign_vector(c(common(
      peer, .DSVERT_CLIENT_VECTOR_RELEASE_VERSION, "vector_released"), list(
        result_set_hash = result_set_hash,
        final_vector_root = final_root,
        final_chunk_hashes = list(chunk_hash),
        output_lattice_bits = as.numeric(log2(scale)),
        output_lattice_scale = scale,
        mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
        epsilon = "1e+00",
        delta = format(
          release$delta, digits = 17, scientific = TRUE, trim = TRUE),
        implementation_delta_numerator = "1",
        implementation_delta_denominator =
          "1267650600228229401496703205376",
        delta_aggregation = "max_per_peer_not_sum",
        postprocessing = paste0(
          "signed-Ring128-decode-then-fixed-public-coordinate-clamp-v1"),
        intermediate_payload_exposed = FALSE,
        durable_replay = TRUE, capability_available = TRUE)), peer)
  }), designated)
  acks <- stats::setNames(lapply(peers, function(peer) {
    sign_vector(c(common(
      peer, .DSVERT_CLIENT_VECTOR_ACK_VERSION,
      "vector_finalized_and_compacted"), list(
        final_vector_root = final_root,
        source_intermediates_compacted = TRUE,
        sampler_intermediates_compacted = TRUE,
        final_chunks_retained = TRUE,
        durable_replay_retained = TRUE, idempotent = TRUE)), peer)
  }), peers)
  release$signed_provenance <- list(
    version = "dsvert-joint-dp-vector-public-provenance-v1",
    ordered_peer_pinset = as.list(pins), peer_pinset_sha256 = pin_hash,
    designated_noise_peers = as.list(designated),
    release_instance_id = release_instance$id,
    release_instance = release_instance$value,
    prepare_receipts = prepare, release_receipts = releases,
    finalization_receipts = acks, protected_shares_included = FALSE,
    preclamp_values_included = FALSE, source_values_included = FALSE)
  bundle <- list(
    schema_json = schema_json, schema_sha256 = schema_sha256,
    workload_contract_sha256 = workload_sha256,
    manifest_sha256 = release$manifest_sha256, capsule_id = capsule_id,
    artifact_commitments = index$value,
    artifact_commitment_count = index$count,
    artifact_commitments_root = index$root,
    manifest_build_receipts = build_receipts)
  list(
    manifest = manifest, release = release, status = status,
    bundle = bundle, pins = pins)
}

.dp_gaussian_fixture <- function(singular = FALSE) {
  capacity <- 100
  scale <- 256
  artifact <- .dp_gaussian_artifact(capacity, scale)
  count <- list(
    owner_peer = "site_a", dataset = "cohort",
    statistic_maximum = capacity, l1_sensitivity = 1)
  manifest <- list(workload = list(
    coordinate_count = 8,
    release_lattice = list(
      version = "biomedical-capsule-common-lattice-v1",
      output_lattice_bits = 8L, output_lattice_scale = scale,
      natural_l1_sensitivity = 8,
      integer_l1_sensitivity_steps = 8 * scale,
      natural_l2_sensitivity = sqrt(8),
      integer_l2_sensitivity_steps = sqrt(8) * scale),
    families = list(
      admitted_count = count,
      numeric_moments = list(artifacts = list()),
      numeric_pair_moments = list(
        artifacts = list(), natural_l1_sensitivity = 0),
      gaussian_models = list(
        artifacts = list(gaussian_primary = artifact),
        natural_l1_sensitivity = 7),
      fixed_numeric_histograms = list(artifacts = list()),
      categorical_marginals = list(artifacts = list()),
      categorical_pairs = list(sets = list()),
      correlation_artifacts = list(), describe_artifacts = list(),
      survival_artifacts = list())))
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  # Normalized exact oracle: x=(0,.5,1), y=(.25,.5,.75), so
  # y=.25+.5*x. Original bounds imply y=15+1.25*x_original.
  gaussian_values <- if (isTRUE(singular)) {
    c(3, 3, 1.5, .75, 1.5, .75, .875)
  } else {
    c(3, 3, 1.5, 1.25, 1.5, 1, .875)
  }
  release <- list(
    capsule_id = strrep("a", 64L), manifest_sha256 = strrep("b", 64L),
    final_vector_root = strrep("c", 64L),
    coordinate_order_sha256 = layout$sha256,
    coordinate_count = 8L, values = c(3, gaussian_values),
    mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
    epsilon = 1, delta = 2^-100,
    implementation_delta = "1/1267650600228229401496703205376",
    delta_aggregation = "max_per_peer_not_sum", sticky_replay = TRUE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    manifest = manifest)
  class(release) <- c("dsvert_joint_dp_vector", "list")
  status <- list(
    site_a = list(
      policy = list(adjacency = "add_remove_patient"),
      noise_root = list(
        epoch = 1, provider_id = "provider-site_a", key_id = "root-a")),
    site_b = list(
      policy = list(adjacency = "add_remove_patient"),
      noise_root = list(
        epoch = 1, provider_id = "provider-site_b", key_id = "root-b")))
  certified <- .dp_gaussian_certificate_fixture(
    manifest, layout, release, status, capacity, scale)
  manifest <- certified$manifest
  release <- certified$release
  status <- certified$status
  conns <- list(
    site_a = structure(list(), class = "fake"),
    site_b = structure(list(), class = "fake"))
  list(
    artifact = artifact, manifest = manifest, layout = layout,
    release = release, status = status, conns = conns,
    run = list(
      release = release, layout = layout, status = status,
      manifest_bundle = certified$bundle),
    pins = certified$pins)
}

.dp_gaussian_cross_frontdoor_fixture <- function(k = 3L) {
  stopifnot(k %in% c(2L, 3L, 5L))
  peers <- c("site_a", "site_b", "site_c", "site_d", "site_e")[seq_len(k)]
  designated <- peers[1:2]
  predictor_owner <- if (k == 2L) "site_a" else "site_c"
  capacity <- 100
  scale <- 256
  coordinate_count <- 7L
  product_error <- 2 + 1 / (4 * scale)
  artifact <- list(
    version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION,
    spec_version = "v2", analysis_id = "cross_model",
    dataset = "outcome_data", owner_peer = "site_b",
    outcome = list(
      column = "y", dataset = "outcome_data", owner_peer = "site_b",
      lower = 0, upper = 1),
    predictors = list(x = list(
      column = "x", dataset = "predictor_data", owner_peer = predictor_owner,
      lower = -1, upper = 1)),
    predictor_order = list("x"), input_variable_order = list("x", "y"),
    participating_peers = as.list(sort(
      c("site_b", predictor_owner), method = "radix")),
    computation_peers = as.list(designated), intercept = TRUE,
    design_terms = list("(Intercept)", "x"), numeric_grid_bits = 8,
    coordinate_count = coordinate_count,
    coordinate_order = paste0(
      "n_then_xtx_upper_column_major_then_xty_design_order_then_yty_v2"),
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    private_input_layout = paste0(
      "capacity_padded_value_then_validity_per_signed_variable_",
      "manifest_order_v1"),
    repeated_record_policy = paste0(
      "clip_finite_rows_then_mean_each_variable_once_per_admitted_unit_v1"),
    missingness_policy = paste0(
      "complete_case_mask_remains_secret_shared_through_joint_noise_v1"),
    contribution_domain = paste0(
      "round_normalized_inputs_then_exact_floor_ring128_products_on_",
      "closed_unit_interval_v1"),
    count_gram_intercept_policy = paste0(
      "n_and_all_moments_share_one_secret_complete_case_mask_and_are_",
      "released_only_after_joint_dp_v1"),
    statistic_maximum = as.list(rep(
      capacity * scale, coordinate_count)),
    source_raw_l1_sensitivity = coordinate_count * scale,
    source_raw_l2_sensitivity = sqrt(coordinate_count) * scale,
    natural_l1_sensitivity = coordinate_count,
    natural_l2_sensitivity = sqrt(coordinate_count),
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste0(
      "zero_missing_complete_case_vs_all_one_complete_unit_is_worst_",
      "case_for_add_remove_and_replace_one"),
    quantization_contract = list(
      input_rounding = "nearest_integer_ties_to_even_r_v1",
      product_rounding = "exact_signed_floor_after_division_by_scale_v1",
      per_product_max_abs_error_lattice_steps = product_error,
      per_product_max_abs_error_normalized = product_error / scale,
      per_sum_max_abs_error_lattice_steps = capacity * product_error,
      same_owner_v1_numerically_identical = FALSE),
    numeric_certificate = list(
      version = "dsvert-cross-gaussian-numeric-certificate-v1",
      ring_bits = 128, frac_bits = 8, required_signed_bits = 26,
      operand_maximum = scale, raw_product_maximum = scale^2,
      accumulated_coordinate_maximum = capacity * scale,
      truncation = "exact_signed_floor_gc_ot_or_direct_wide_v1",
      comparison = "not_used_after_custodian_bound_clipping",
      modular_wrap_proved_absent = TRUE,
      overflow_behavior = "typed_abort_before_commit"),
    transcript = list(
      version = "dsvert-cross-gaussian-fixed-transcript-v1",
      padded_units = capacity, variable_count = 2,
      validity_product_rounds = 1, masked_value_rounds = 1,
      moment_product_rounds = 1, data_dependent_branches = 0,
      exact_intermediate_release_count = 0),
    alignment_contract = list(
      version = "private-psi-ordered-manifest-consensus-v1",
      public_patient_dependent_hash = FALSE,
      mismatch_behavior = "typed_non_prealigned_cohort_failure"),
    regularization_policy =
      "none_in_release_explicit_client_postprocessing_only_v1",
    implementation_state = "cross_owner_exact_gc_materialized",
    cross_owner_state = "exact_gc_to_joint_dp_vector_v1")
  count <- list(
    owner_peer = "site_a", dataset = "count_data",
    statistic_maximum = capacity, l1_sensitivity = 1)
  manifest <- list(workload = list(
    coordinate_count = 8,
    release_lattice = list(
      version = "biomedical-capsule-common-lattice-v1",
      output_lattice_bits = 8L, output_lattice_scale = scale,
      natural_l1_sensitivity = 8,
      integer_l1_sensitivity_steps = 8 * scale,
      natural_l2_sensitivity = sqrt(8),
      integer_l2_sensitivity_steps = sqrt(8) * scale),
    families = list(
      admitted_count = count,
      numeric_moments = list(artifacts = list()),
      numeric_pair_moments = list(
        artifacts = list(), natural_l1_sensitivity = 0),
      gaussian_models = list(
        artifacts = list(cross_model = artifact),
        natural_l1_sensitivity = coordinate_count),
      fixed_numeric_histograms = list(artifacts = list()),
      categorical_marginals = list(artifacts = list()),
      categorical_pairs = list(sets = list()),
      correlation_artifacts = list(), describe_artifacts = list(),
      survival_artifacts = list())))
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  gaussian_values <- c(3, 3, 1.5, 1.25, 1.5, 1, 0.875)
  release <- list(
    capsule_id = strrep("a", 64L), manifest_sha256 = strrep("b", 64L),
    final_vector_root = strrep("c", 64L),
    coordinate_order_sha256 = layout$sha256,
    coordinate_count = 8L, values = c(3, gaussian_values),
    mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
    epsilon = 1, delta = 2^-100,
    implementation_delta = "1/1267650600228229401496703205376",
    delta_aggregation = "max_per_peer_not_sum", sticky_replay = TRUE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    manifest = manifest)
  class(release) <- c("dsvert_joint_dp_vector", "list")
  certified <- .dp_gaussian_certificate_fixture(
    manifest, layout, release, status = list(), capacity, scale,
    peers = peers, designated = designated)
  conns <- stats::setNames(lapply(peers, function(peer) {
    structure(list(peer = peer), class = "fake")
  }), peers)
  list(
    artifact = artifact, manifest = certified$manifest, layout = layout,
    release = certified$release, status = certified$status, conns = conns,
    run = list(
      release = certified$release, layout = layout,
      status = certified$status, manifest_bundle = certified$bundle),
    pins = certified$pins)
}

.dp_gaussian_recertify <- function(fixture) {
  scale <- fixture$release$manifest$workload$release_lattice$
    output_lattice_scale
  capacity <- fixture$artifact$statistic_maximum[[1L]]
  certified <- .dp_gaussian_certificate_fixture(
    fixture$manifest, fixture$layout, fixture$release, fixture$status,
    capacity, scale)
  fixture$manifest <- certified$manifest
  fixture$release <- certified$release
  fixture$status <- certified$status
  fixture$run$release <- certified$release
  fixture$run$status <- certified$status
  fixture$run$manifest_bundle <- certified$bundle
  fixture$pins <- certified$pins
  fixture
}

# The algebraic and legacy-v3 certificate cases below deliberately start after
# the current Synopsis reader.  The Synopsis reader itself is tested in
# test-dp-gaussian-synopsis-release.R and the K=2/K=3/K=5 Rock E2E test; this
# helper keeps the older certificate compatibility cases from pretending to
# exercise a retired production runner.
.dp_gaussian_legacy_released <- function(
    fixture, data_name, analysis_id, server = NULL, datasources = NULL,
    .aggregate = NULL) {
  context <- .dsvert_dp_vector_context(fixture$run, allow_synopsis = FALSE)
  metadata <- .dsvert_dp_vector_public_metadata(context)
  count_block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count",
    description = "signed admitted-count capacity block")
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  artifact <- .dsvert_dp_gaussian_artifact(
    context$manifest, data_name, analysis_id, server, context$adjacency,
    as.numeric(context$lattice$output_lattice_scale), capacity)
  blocks <- .dsvert_dp_capsule_vector_blocks(
    context$layout, "gaussian_models", dataset = data_name,
    owner_peer = artifact$owner_peer)
  expect_length(blocks, 1L)
  block <- blocks[[1L]]
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, block)
  certificate <- .dsvert_dp_gaussian_certificate_build(
    context, artifact, block, coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate)
  list(
    context = context, metadata = metadata, artifact = verification$artifact,
    block = block, coordinates = coordinates,
    moment = .dsvert_dp_gaussian_unpack(
      coordinates, verification$artifact, capacity),
    certificate = certificate, verification = verification,
    scale = as.numeric(context$lattice$output_lattice_scale),
    capacity = capacity)
}

test_that("DP Gaussian matches normalized and original-scale oracle", {
  fixture <- .dp_gaussian_fixture()
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      calls <<- calls + 1L
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .dsvert_aggregate_strict = function(...) {
      stop("legacy aggregate route called", call. = FALSE)
    },
    .dsvert_fanout_by_site = function(...) {
      stop("legacy fanout route called", call. = FALSE)
    },
    .package = "dsVertClient")
  fit <- .dsvert_dp_gaussian_impl(
    "cohort", "gaussian_primary", 0, "site_a", fixture$conns,
    function(...) stop("raw DSI call", call. = FALSE))
  expect_identical(calls, 1L)
  expect_equal(
    fit$coefficients_normalized,
    c(`(Intercept)` = .25, x = .5), tolerance = 1e-10)
  expect_equal(
    fit$coefficients_original_scale,
    c(`(Intercept)` = 15, x = 1.25), tolerance = 1e-10)
  expect_false(fit$regularized)
  expect_true(fit$identifiability$full_rank)
  expect_identical(
    fit$count_gram_intercept_reconciliation$policy,
    fixture$artifact$count_gram_intercept_policy)
  expect_false(fit$count_gram_intercept_reconciliation$silently_averaged)
  expect_equal(fit$count_gram_intercept_reconciliation$discrepancy, 0)
  expect_false(fit$augmented_moment_projection$
                 exact_nearest_with_all_moment_constraints)
  expect_identical(fit$source_values_exposed, FALSE)
  expect_identical(fit$intermediate_values_exposed, FALSE)
  expect_identical(fit$legacy_fallback_called, FALSE)
  expect_identical(fit$additional_server_calls_after_synopsis, 0L)
  expect_identical(
    fit$additional_privacy_cost, c(epsilon = 0, delta = 0))
  expect_false(fit$inference$sampling_inference_available)
  expect_identical(fit$accuracy$coordinate_count, 7L)
  expect_false(fit$accuracy$coefficient_regions_available)
  expect_true(fit$accuracy$max_abs_quantization_per_non_count_sum > 0)
  expect_false(any(c(
    "data", "rows", "fitted_values", "residuals") %in% names(fit)))
})

test_that("Gaussian post-processing is deterministic for one sticky release", {
  fixture <- .dp_gaussian_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .package = "dsVertClient")
  first <- .dsvert_dp_gaussian_impl(
    "cohort", "gaussian_primary", 0, NULL, fixture$conns,
    function(...) NULL)
  second <- .dsvert_dp_gaussian_impl(
    "cohort", "gaussian_primary", 0, NULL, fixture$conns,
    function(...) NULL)
  expect_identical(second, first)
})

test_that("Gaussian mechanism regions are certificate-bound post-processing", {
  fixture <- .dp_gaussian_fixture()
  artifact <- fixture$artifact
  capacity <- artifact$statistic_maximum[[1L]]
  coordinates <- fixture$release$values[-1L]
  moment <- .dsvert_dp_gaussian_unpack(coordinates, artifact, capacity)
  ridge <- 1000
  normalized <- .dsvert_dp_gaussian_solve(moment, artifact, ridge)$coefficients
  original <- .dsvert_dp_gaussian_original_coefficients(normalized, artifact)
  verified <- list(
    integrity_valid = TRUE,
    authenticity = "session_transport_anchored",
    artifact = artifact,
    validated_moment = moment,
    accuracy_simultaneous_95 = list(confidence = 0.95, radius = 0.01),
    output_lattice_scale = 256,
    coordinate_capacity = capacity,
    coordinates = coordinates)
  fit <- structure(list(
    ridge = ridge,
    coefficients_normalized = normalized,
    coefficients_original_scale = original,
    coefficients = original,
    accuracy = list(simultaneous_abs_mechanism_radius = 100)),
    class = c("ds.vertDPGaussian", "list"))
  testthat::local_mocked_bindings(
    ds.validateDPGaussianCertificate = function(...) verified,
    .package = "dsVertClient")

  region <- ds.vertConfint(fit, type = "mechanism")
  expect_s3_class(region, "data.frame")
  expect_identical(attr(region, "interval_scope"),
                   "simultaneous_dp_mechanism")
  expect_identical(attr(region, "sampling_inference"), FALSE)
  expect_true(all(is.finite(as.matrix(region))))
  expect_true(all(region$lower <= region$estimate))
  expect_true(all(region$estimate <= region$upper))
  expect_true(all(region$mechanism_radius >= 0))
  selected <- ds.vertConfint(fit, parm = "x", type = "mechanism")
  expect_identical(rownames(selected), "x")
  aliased <- ds.vert.confint(fit, type = "mechanism")
  expect_identical(aliased$estimate, region$estimate)
  expect_identical(aliased$mechanism_radius, region$mechanism_radius)

  wald <- ds.vertWald(fit, parm = "x", null = 0, type = "mechanism")
  expect_identical(wald$distribution, "simultaneous_dp_mechanism_region")
  expect_false(wald$sampling_inference)
  expect_null(wald$p_value)
  expect_equal(wald$lower, region["x", "lower"])
  expect_equal(wald$upper, region["x", "upper"])
  expect_identical(wald$null_excluded,
                   isTRUE(0 < wald$lower || 0 > wald$upper))
  expect_identical(
    ds.vert.wald(fit, parm = "x", null = 0, type = "mechanism")$lower,
    wald$lower)

  contrast_matrix <- matrix(c(1, -2), nrow = 1L)
  contrast <- ds.vertContrast(
    fit, K = contrast_matrix, m = 0, type = "mechanism")
  expected_estimate <- region["(Intercept)", "estimate"] -
    2 * region["x", "estimate"]
  expected_radius <- region["(Intercept)", "mechanism_radius"] +
    2 * region["x", "mechanism_radius"]
  expect_identical(contrast$distribution,
                   "simultaneous_dp_mechanism_region")
  expect_false(contrast$sampling_inference)
  expect_null(contrast$p_value)
  expect_equal(contrast$estimate, expected_estimate)
  expect_equal(contrast$mechanism_radius, expected_radius)
  expect_equal(contrast$lower, expected_estimate - expected_radius)
  expect_equal(contrast$upper, expected_estimate + expected_radius)
  expect_identical(contrast$null_excluded,
                   contrast$lower > 0 || contrast$upper < 0)
  expect_equal(
    ds.vert.contrast(fit, K = contrast_matrix, type = "mechanism")$upper,
    contrast$upper)
  expect_error(ds.vertWald(fit, parm = "missing", type = "mechanism"),
               "Unknown coefficient")
  expect_error(ds.vertContrast(fit, K = matrix(1, nrow = 1L, ncol = 1L),
                               type = "mechanism"),
               "ncol\\(K\\)")

  zero_solution <- .dsvert_dp_gaussian_solve(moment, artifact, ridge = 0)
  zero_original <- .dsvert_dp_gaussian_original_coefficients(
    zero_solution$coefficients, artifact)
  zero_fit <- fit
  zero_fit$ridge <- 0
  zero_fit$coefficients_normalized <- zero_solution$coefficients
  zero_fit$coefficients_original_scale <- zero_original
  zero_fit$coefficients <- zero_original
  reduction <- ds.vertLR(
    reduced = "(Intercept)", full = zero_fit, type = "mechanism")
  expected_reduced <- sum(moment$cross["(Intercept)"]^2 /
                            moment$gram["(Intercept)", "(Intercept)"])
  expected_full <- sum(zero_solution$coefficients * moment$cross)
  expect_identical(reduction$distribution,
                   "simultaneous_dp_mechanism_rss_reduction")
  expect_false(reduction$sampling_inference)
  expect_null(reduction$p_value)
  expect_equal(reduction$statistic, expected_full - expected_reduced)
  expect_true(is.finite(reduction$lower))
  expect_true(is.finite(reduction$upper))
  expect_gte(reduction$lower, 0)
  expect_lte(reduction$lower, reduction$statistic)
  expect_gte(reduction$upper, reduction$statistic)
  expect_equal(
    ds.vert.lr(reduced = "(Intercept)", full = zero_fit,
               type = "mechanism")$upper,
    reduction$upper)
  expect_error(ds.vertLR("(Intercept)", fit, type = "mechanism"),
               "unpenalized")
  expect_error(ds.vertLR("unknown", zero_fit, type = "mechanism"),
               "Unknown reduced-model")

  tampered_accuracy <- fit
  tampered_accuracy$accuracy$simultaneous_abs_mechanism_radius <- 0
  expect_identical(
    ds.vertConfint(tampered_accuracy, type = "mechanism"), region)
  expect_error(ds.vertConfint(fit, type = "mechanism", level = 0.9),
               "current signed Synopsis certificate")
  tampered_coefficients <- fit
  tampered_coefficients$coefficients[[1L]] <-
    tampered_coefficients$coefficients[[1L]] + 1
  expect_error(ds.vertConfint(tampered_coefficients, type = "mechanism"),
               "does not match its signed sufficient statistics")
  expect_error(ds.vertConfint(fit, parm = "missing", type = "mechanism"),
               "Unknown coefficient")

  too_wide <- verified
  too_wide$accuracy_simultaneous_95$radius <- capacity
  condition <- tryCatch(
    testthat::with_mocked_bindings(
      ds.validateDPGaussianCertificate = function(...) too_wide,
      ds.vertConfint(fit, type = "mechanism"),
      .package = "dsVertClient"),
    non_identifiable = function(error) error)
  expect_s3_class(condition, "non_identifiable")
  expect_identical(condition$reason, "dp_mechanism_region_not_identifiable")
})

test_that("Gaussian provenance is offline-verifiable with calibrated trust", {
  fixture <- .dp_gaussian_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .package = "dsVertClient")
  fit <- .dsvert_dp_gaussian_impl(
    "cohort", "gaussian_primary", 0, NULL, fixture$conns,
    function(...) NULL)
  online <- ds.validateDPGaussianCertificate(fit)
  expect_true(online$integrity_valid)
  expect_identical(online$authenticity, "session_transport_anchored")
  expect_identical(fit$provenance_authenticity,
                   "session_transport_anchored")
  expect_identical(
    .dsvert_vector_hash(online$mechanism_plan),
    online$mechanism_plan_sha256)
  expect_identical(online$capsule_mechanism, "discrete-laplace")
  expect_identical(online$backend, .DSVERT_CLIENT_VECTOR_BACKEND)
  expect_identical(online$sampler, .DSVERT_CLIENT_VECTOR_SAMPLER)
  expect_identical(online$mechanism,
                   .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM)
  expect_identical(
    online$proof_geometry$granularity,
    "whole_intersecting_public_chunks_v1")
  expect_gte(online$proof_geometry$included_public_coordinates,
             online$proof_geometry$block_coordinates)
  expect_identical(
    online$proof_geometry$overfetch_coordinates,
    online$proof_geometry$included_public_coordinates -
      online$proof_geometry$block_coordinates)
  expect_gt(online$proof_geometry$certificate_bytes, 0)

  certificate <- unserialize(serialize(
    fit$provenance_certificate, NULL, version = 3L))
  pin_hash <- certificate$peer_context$peer_pinset_sha256
  if (exists(pin_hash, envir = .dsvert_dp_gaussian_trust_cache,
             inherits = FALSE)) {
    rm(list = pin_hash, envir = .dsvert_dp_gaussian_trust_cache)
  }
  offline <- ds.validateDPGaussianCertificate(certificate)
  expect_true(offline$integrity_valid)
  expect_identical(offline$authenticity, "unanchored")
  caller <- ds.validateDPGaussianCertificate(
    certificate, trusted_pinset = fixture$pins)
  expect_identical(caller$authenticity, "caller_anchored")

  wrong <- fixture$pins
  wrong[[1L]] <- wrong[[2L]]
  expect_error(ds.validateDPGaussianCertificate(
    certificate, trusted_pinset = wrong), "trusted pinset")
})

test_that("one authentic Gaussian release supports LASSO and honest pseudo-IC", {
  fixture <- .dp_gaussian_fixture()
  # Preserve the fitted normal equations while leaving a strictly positive
  # projected residual moment, so the explicitly heuristic pseudo-IC is
  # defined. The augmented moment matrix remains PSD.
  fixture$release$values[[8L]] <- 1.125
  fixture <- .dp_gaussian_recertify(fixture)
  capsule_calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      capsule_calls <<- capsule_calls + 1L
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .dsvert_aggregate_strict = function(...) {
      stop("unexpected additional DSI call", call. = FALSE)
    },
    .dsvert_fanout_by_site = function(...) {
      stop("unexpected additional DSI fan-out", call. = FALSE)
    },
    .package = "dsVertClient")

  gaussian <- ds.vertDPGaussian(
    "cohort", "gaussian_primary", datasources = fixture$conns)
  certificate_hash <- gaussian$provenance_certificate$certificate_sha256
  verified <- ds.validateDPGaussianCertificate(gaussian)
  proximal <- ds.vertLASSOProximal(
    gaussian, lambda = 0.05, max_iter = 5000L, tol = 1e-12)
  selected <- ds.vertLASSOCV(
    gaussian, lambda_grid = c(0.1, 0.05, 0), criterion = "BIC")
  caller_anchored_proximal <- ds.vertLASSOProximal(
    gaussian, lambda = 0.05, max_iter = 5000L, tol = 1e-12,
    trusted_pinset = fixture$pins)

  expect_identical(capsule_calls, 1L)
  expect_true(verified$integrity_valid)
  expect_true(proximal$source_certificate_validation$integrity_valid)
  expect_true(selected$source_certificate_validation$integrity_valid)
  expect_identical(
    proximal$source_certificate_validation$authenticity,
    "session_transport_anchored")
  expect_identical(
    selected$source_certificate_validation$authenticity,
    "session_transport_anchored")
  expect_identical(
    caller_anchored_proximal$source_certificate_validation$authenticity,
    "caller_anchored")
  expect_identical(
    gaussian$provenance_certificate$certificate_sha256, certificate_hash)
  expect_true(proximal$kkt$satisfied)
  expect_lte(proximal$kkt$max_violation, proximal$kkt$tolerance)
  expect_true(all(vapply(
    selected$path_certificates,
    function(value) isTRUE(value$kkt$satisfied), logical(1L))))
  expect_identical(proximal$additional_server_calls_after_capsule, 0L)
  expect_identical(selected$additional_server_calls_after_capsule, 0L)
  expect_equal(proximal$additional_privacy_cost,
               c(epsilon = 0, delta = 0))
  expect_equal(selected$additional_privacy_cost,
               c(epsilon = 0, delta = 0))
  expect_true(selected$selection_available)
  expect_identical(
    selected$selection_method,
    "DP_projected_pseudo_information_criterion")
  expect_false(selected$classical_information_criterion)
  expect_false(selected$cross_validation)
  expect_false(selected$one_standard_error_rule)
  expect_false(selected$inference$sampling_inference_available)
})

test_that("Gaussian provenance rejects tampering and forged signed evidence", {
  fixture <- .dp_gaussian_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .package = "dsVertClient")
  fit <- .dsvert_dp_gaussian_impl(
    "cohort", "gaussian_primary", 0, NULL, fixture$conns,
    function(...) NULL)
  rehash <- function(certificate) {
    certificate$certificate_sha256 <-
      .dsvert_dp_gaussian_certificate_hash(certificate[
        setdiff(names(certificate), "certificate_sha256")])
    certificate
  }

  chunk_tamper <- fit$provenance_certificate
  chunk_tamper$chunk_evidence[[1L]]$chunk$scaled_values[[2L]] <- "1"
  chunk_tamper <- rehash(chunk_tamper)
  expect_error(ds.validateDPGaussianCertificate(chunk_tamper),
               "chunk|commitment|Merkle")

  release_forgery <- fit$provenance_certificate
  peer <- names(release_forgery$signed_evidence$
                  vector_release_receipts)[[1L]]
  release_forgery$signed_evidence$vector_release_receipts[[peer]]$
    mechanism <- .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM
  release_forgery <- rehash(release_forgery)
  expect_error(ds.validateDPGaussianCertificate(release_forgery),
               "signature|receipt")

  plan_forgery <- fit$provenance_certificate
  peer <- names(plan_forgery$signed_evidence$
                  vector_prepare_receipts)[[1L]]
  plan_forgery$signed_evidence$vector_prepare_receipts[[peer]]$
    mechanism_plan$total_coordinate_count <- 999
  plan_forgery <- rehash(plan_forgery)
  expect_error(ds.validateDPGaussianCertificate(plan_forgery),
               "signature|receipt")

  profile_tamper <- fit$provenance_certificate
  profile_tamper$backend <- .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND
  profile_tamper <- rehash(profile_tamper)
  expect_error(ds.validateDPGaussianCertificate(profile_tamper),
               "contract|profile|query|vector")
})

test_that("Gaussian count and intercept Gram coordinates are never averaged", {
  fixture <- .dp_gaussian_fixture()
  # Coordinate 2 is Gaussian n; coordinate 3 is Gram[(Intercept),(Intercept)].
  fixture$release$values[[3L]] <- 2.5
  fixture <- .dp_gaussian_recertify(fixture)
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .package = "dsVertClient")
  fit <- .dsvert_dp_gaussian_impl(
    "cohort", "gaussian_primary", 0, NULL, fixture$conns,
    function(...) NULL)
  expect_equal(fit$n_obs, 3)
  expect_identical(
    fit$n_obs_definition,
    paste0(
      "DP_noisy_complete_case_count_coordinate_not_averaged_with_",
      "Gram11"))
  expect_equal(
    fit$count_gram_intercept_reconciliation$
      gram_intercept_intercept_dp,
    2.5)
  expect_equal(
    fit$count_gram_intercept_reconciliation$discrepancy,
    -0.5)
  expect_false(fit$count_gram_intercept_reconciliation$silently_averaged)
})

test_that("singular Gaussian design needs explicit regularization", {
  fixture <- .dp_gaussian_fixture(singular = TRUE)
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .package = "dsVertClient")
  condition <- tryCatch(
    .dsvert_dp_gaussian_impl(
      "cohort", "gaussian_primary", 0, NULL, fixture$conns,
      function(...) NULL),
    error = identity)
  expect_s3_class(condition, "non_identifiable")
  expect_identical(condition$reason, "singular_dp_gaussian_design")
  expect_match(condition$message, "singular")
  ridge <- .dsvert_dp_gaussian_impl(
    "cohort", "gaussian_primary", 0.1, NULL, fixture$conns,
    function(...) NULL)
  expect_true(ridge$regularized)
  expect_identical(
    ridge$regularization_estimand,
    "explicit_predictor_only_ridge_on_normalized_design")
  expect_equal(ridge$solve_diagnostics$ridge, 0.1)
})

test_that("Gaussian non-estimability failures are typed", {
  fixture <- .dp_gaussian_fixture()
  count_condition <- tryCatch(
    .dsvert_dp_gaussian_unpack(
      c(0, fixture$release$values[-c(1L, 2L)]), fixture$artifact,
      fixture$artifact$statistic_maximum[[1L]]),
    error = identity)
  expect_s3_class(count_condition, "non_identifiable")
  expect_identical(
    count_condition$reason, "non_positive_dp_complete_case_count")

  moment <- list(
    gram = matrix(c(0, 0, 0, 1), 2L), cross = c(0, 0),
    outcome_square = 0)
  artifact <- fixture$artifact
  invert_condition <- tryCatch(
    .dsvert_dp_gaussian_solve(moment, artifact, ridge = 0.1),
    error = identity)
  expect_s3_class(invert_condition, "non_identifiable")
  expect_identical(
    invert_condition$reason, "noninvertible_dp_gaussian_system")
})

test_that("Gaussian rejects descriptor tampering and wrong owner", {
  fixture <- .dp_gaussian_fixture()
  wrong_owner <- fixture
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = wrong_owner), list(...)))
    },
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_gaussian_impl(
      "cohort", "gaussian_primary", 0, "site_b", fixture$conns,
      function(...) NULL),
    "reserved_not_materialized")

  tampered <- fixture
  tampered$manifest$workload$families$gaussian_models$artifacts$
    gaussian_primary$missingness_policy <- "available_case"
  tampered$release$manifest <- tampered$manifest
  tampered$run$release <- tampered$release
  tampered$run$layout <- .dsvert_dp_capsule_vector_layout(tampered$manifest)
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = tampered), list(...)))
    },
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_gaussian_impl(
      "cohort", "gaussian_primary", 0, NULL, fixture$conns,
      function(...) NULL),
    "descriptor is invalid")

  layout_tampered <- fixture
  gaussian_block <- which(vapply(
    layout_tampered$run$layout$blocks,
    function(block) identical(block$family, "gaussian_models"),
    logical(1L)))
  expect_length(gaussian_block, 1L)
  layout_tampered$run$layout$blocks[[gaussian_block]]$descriptor$
    outcome$lower <- 11
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = layout_tampered), list(...)))
    },
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_gaussian_impl(
      "cohort", "gaussian_primary", 0, NULL, fixture$conns,
      function(...) NULL),
    "committed coordinate block")

  release_tampered <- fixture
  release_tampered$release$values[[2L]] <- 101
  release_tampered$run$release <- release_tampered$release
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = release_tampered), list(...)))
    },
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_gaussian_impl(
      "cohort", "gaussian_primary", 0, NULL, fixture$conns,
      function(...) NULL),
    "signed final DP vector")
})

test_that("ds.vertGLM Gaussian capsule adapter is explicit and legacy-free", {
  fixture <- .dp_gaussian_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .dsvert_aggregate_strict = function(...) {
      stop("column discovery or legacy aggregate called", call. = FALSE)
    },
    .package = "dsVertClient")
  fit <- ds.vertGLM(
    y ~ x, data = "cohort", family = "gaussian",
    dp_analysis_id = "gaussian_primary",
    missing = "complete_case_capsule", verbose = FALSE,
    datasources = fixture$conns)
  expect_s3_class(fit, "ds.vertDPGaussian")
  expect_identical(fit$called_via, "ds.vertGLM_explicit_dp_analysis_id")
  expect_identical(fit$legacy_glm_estimand, FALSE)
  expect_error(
    ds.vertGLM(
      y ~ x, data = "cohort", family = "gaussian", max_iter = 2,
      dp_analysis_id = "gaussian_primary", verbose = FALSE,
      datasources = fixture$conns),
    "legacy iterative estimand")
  expect_error(
    ds.vertGLM(
      y ~ x, data = "cohort", family = "binomial",
      dp_analysis_id = "gaussian_primary", verbose = FALSE,
      datasources = fixture$conns),
    "only for family='gaussian'")
})

test_that("Gaussian front door cannot reach legacy model endpoints", {
  namespace <- asNamespace("dsVertClient")
  reachable <- character()
  queue <- c("ds.vertDPGaussian", ".dsvert_dp_gaussian_impl")
  while (length(queue)) {
    name <- queue[[1L]]
    queue <- queue[-1L]
    if (name %in% reachable ||
        !exists(name, namespace, mode = "function", inherits = FALSE)) next
    value <- get(name, namespace, inherits = FALSE)
    reachable <- c(reachable, name)
    globals <- tryCatch(
      unique(unlist(codetools::findGlobals(value, merge = FALSE),
                    use.names = FALSE)),
      error = function(error) character())
    queue <- unique(c(queue, intersect(
      globals, ls(namespace, all.names = TRUE))))
  }
  forbidden <- c(
    "dsvertColNamesDS", "glmStandardizeDS", "glmPartialFitDS",
    "k2GradientR1DS", "k2GradientR2DS", "glmSecureGradientDS")
  bodies <- paste(vapply(reachable, function(name) {
    paste(deparse(body(get(name, namespace, inherits = FALSE))),
          collapse = "\n")
  }, character(1L)), collapse = "\n")
  expect_contains(reachable, ".dsvert_dp_synopsis_vector_run")
  expect_false(".dsvert_dp_capsule_vector_run" %in% reachable)
  expect_length(intersect(reachable, forbidden), 0L)
  for (endpoint in forbidden) {
    expect_false(grepl(endpoint, bodies, fixed = TRUE), info = endpoint)
  }
})

test_that("cross-owner Gaussian K=3 and K=5 are served by both public front doors", {
  for (k in c(3L, 5L)) {
    fixture <- .dp_gaussian_cross_frontdoor_fixture(k)
    capsule_calls <- 0L
    testthat::local_mocked_bindings(
      .dsvert_dp_gaussian_synopsis_release = function(
          data_name, analysis_id, server = NULL, datasources = NULL, .aggregate) {
        capsule_calls <<- capsule_calls + 1L
        expect_identical(names(datasources), names(fixture$conns))
        .dp_gaussian_legacy_released(
          fixture, data_name, analysis_id, server, datasources, .aggregate)
      },
      .dsvert_aggregate_strict = function(...) {
        stop("legacy aggregate route called", call. = FALSE)
      },
      .dsvert_fanout_by_site = function(...) {
        stop("unexpected post-capsule fan-out", call. = FALSE)
      },
      .package = "dsVertClient")

    direct <- ds.vertDPGaussian(
      "outcome_data", "cross_model", server = "site_b",
      datasources = fixture$conns)
    adapted <- ds.vertGLM(
      y ~ x, data = "outcome_data", x_vars = list(site_c = "x"),
      y_server = "site_b", family = "gaussian", lambda = 0,
      dp_analysis_id = "cross_model", missing = "complete_case_capsule",
      verbose = FALSE, datasources = fixture$conns)

    expect_identical(capsule_calls, 2L)
    for (fit in list(direct, adapted)) {
      expect_s3_class(fit, "ds.vertDPGaussian")
      expect_equal(fit$coefficients_normalized,
                   c(`(Intercept)` = 0.25, x = 0.5), tolerance = 1e-12)
      expect_equal(fit$coefficients,
                   c(`(Intercept)` = 0.5, x = 0.25), tolerance = 1e-12)
      expect_identical(fit$server, "site_b")
      expect_identical(fit$participating_servers, c("site_b", "site_c"))
      expect_identical(fit$computation_servers, c("site_a", "site_b"))
      expect_identical(fit$provenance_certificate$artifact$version,
                       .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION)
      expect_identical(fit$signed_artifact$spec_version, "v2")
      expect_identical(
        fit$provenance_certificate$artifact_commitment$context$capsule_schema,
        "dsvert-biomedical-capsule-workload-v2")
      expect_identical(
        fit$cross_owner_state, "exact_gc_to_joint_dp_vector_v1")
      expect_identical(
        fit$signed_artifact$transcript$exact_intermediate_release_count, 0)
      expect_identical(fit$source_values_exposed, FALSE)
      expect_identical(fit$intermediate_values_exposed, FALSE)
      expect_identical(fit$legacy_fallback_called, FALSE)
      expect_identical(fit$additional_server_calls_after_synopsis, 0L)
      expect_identical(fit$history_gate, TRUE)
      expect_identical(fit$request_limit, FALSE)
      expect_identical(fit$operation_limit, TRUE)
      expect_true(fit$provenance_integrity)
      expect_identical(fit$provenance_authenticity,
                       "session_transport_anchored")
      verified <- ds.validateDPGaussianCertificate(fit)
      expect_true(verified$integrity_valid)
      caller_anchored <- ds.validateDPGaussianCertificate(
        fit, trusted_pinset = fixture$pins)
      expect_identical(caller_anchored$authenticity, "caller_anchored")
      expect_length(
        verified$certificate$peer_context$ordered_peer_pinset, k)
      expect_identical(
        names(verified$certificate$signed_evidence$manifest_build_receipts),
        names(fixture$conns))
      expect_identical(
        names(verified$certificate$signed_evidence$vector_release_receipts),
        c("site_a", "site_b"))
      expect_true(all(vapply(
        verified$certificate$signed_evidence$vector_release_receipts,
        function(receipt) identical(receipt$history_gate, TRUE) &&
          identical(receipt$request_limit, FALSE) &&
          identical(receipt$operation_limit, TRUE), logical(1L))))
      expect_true(all(vapply(
        verified$certificate$signed_evidence$manifest_build_receipts,
        function(receipt) identical(receipt$request_limit, FALSE) &&
          identical(receipt$operation_limit, FALSE), logical(1L))))
    }
    expect_identical(adapted$called_via,
                     "ds.vertGLM_explicit_dp_analysis_id")
    expect_identical(adapted$legacy_glm_estimand, FALSE)
    expect_identical(
      adapted$provenance_certificate$certificate_sha256,
      direct$provenance_certificate$certificate_sha256)
  }
})

test_that("cross-owner Gaussian K=2 uses both owner-computation peers", {
  fixture <- .dp_gaussian_cross_frontdoor_fixture(k = 2L)
  capsule_calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(
        data_name, analysis_id, server = NULL, datasources = NULL, .aggregate) {
      capsule_calls <<- capsule_calls + 1L
      expect_identical(names(datasources), c("site_a", "site_b"))
      .dp_gaussian_legacy_released(
        fixture, data_name, analysis_id, server, datasources, .aggregate)
    },
    .dsvert_aggregate_strict = function(...) {
      stop("legacy aggregate route called", call. = FALSE)
    },
    .dsvert_fanout_by_site = function(...) {
      stop("unexpected post-capsule fan-out", call. = FALSE)
    },
    .package = "dsVertClient")

  direct <- ds.vertDPGaussian(
    "outcome_data", "cross_model", server = "site_b",
    datasources = fixture$conns)
  adapted <- ds.vertGLM(
    y ~ x, data = "outcome_data", x_vars = list(site_a = "x"),
    y_server = "site_b", family = "gaussian", lambda = 0,
    dp_analysis_id = "cross_model", missing = "complete_case_capsule",
    verbose = FALSE, datasources = fixture$conns)

  expect_identical(capsule_calls, 2L)
  for (fit in list(direct, adapted)) {
    expect_s3_class(fit, "ds.vertDPGaussian")
    expect_equal(fit$coefficients_normalized,
                 c(`(Intercept)` = 0.25, x = 0.5), tolerance = 1e-12)
    expect_equal(fit$coefficients,
                 c(`(Intercept)` = 0.5, x = 0.25), tolerance = 1e-12)
    expect_identical(fit$server, "site_b")
    expect_identical(fit$participating_servers, c("site_a", "site_b"))
    expect_identical(fit$computation_servers, c("site_a", "site_b"))
    expect_identical(
      fit$signed_artifact$predictors$x$owner_peer, "site_a")
    expect_identical(fit$signed_artifact$outcome$owner_peer, "site_b")
    expect_identical(fit$signed_artifact$spec_version, "v2")
    expect_identical(
      fit$cross_owner_state, "exact_gc_to_joint_dp_vector_v1")
    expect_identical(fit$legacy_fallback_called, FALSE)
    expect_identical(fit$additional_server_calls_after_synopsis, 0L)
    expect_identical(fit$history_gate, TRUE)
    expect_identical(fit$request_limit, FALSE)
    expect_identical(fit$operation_limit, TRUE)
    verified <- ds.validateDPGaussianCertificate(
      fit, trusted_pinset = fixture$pins)
    expect_true(verified$integrity_valid)
    expect_identical(verified$authenticity, "caller_anchored")
    expect_length(
      verified$certificate$peer_context$ordered_peer_pinset, 2L)
    expect_identical(
      names(verified$certificate$signed_evidence$manifest_build_receipts),
      c("site_a", "site_b"))
    expect_identical(
      names(verified$certificate$signed_evidence$vector_release_receipts),
      c("site_a", "site_b"))
  }
  expect_identical(adapted$called_via,
                   "ds.vertGLM_explicit_dp_analysis_id")
  expect_identical(adapted$legacy_glm_estimand, FALSE)
  expect_identical(
    adapted$provenance_certificate$certificate_sha256,
    direct$provenance_certificate$certificate_sha256)
})

test_that("cross-owner Gaussian front doors reject owner and signature tampering", {
  fixture <- .dp_gaussian_cross_frontdoor_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = fixture), list(...)))
    },
    .dsvert_aggregate_strict = function(...) {
      stop("legacy fallback called", call. = FALSE)
    },
    .package = "dsVertClient")
  expect_error(ds.vertGLM(
    y ~ x, data = "outcome_data", x_vars = list(site_a = "x"),
    y_server = "site_b", family = "gaussian",
    dp_analysis_id = "cross_model", missing = "complete_case_capsule",
    verbose = FALSE, datasources = fixture$conns),
    "formula does not match.*no legacy fallback")

  forged <- fixture
  forged$run$release$signed_provenance$release_receipts$site_a$
    request_limit <- TRUE
  testthat::local_mocked_bindings(
    .dsvert_dp_gaussian_synopsis_release = function(...) {
      do.call(.dp_gaussian_legacy_released, c(list(fixture = forged), list(...)))
    },
    .package = "dsVertClient")
  expect_error(ds.vertDPGaussian(
    "outcome_data", "cross_model", server = "site_b",
    datasources = forged$conns), "signature|receipt")
})
