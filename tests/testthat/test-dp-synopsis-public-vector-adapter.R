.synopsis_public_helpers <- local({
  environment <- new.env(parent = asNamespace("dsVertClient"))
  for (file in c(
      "test-dp-capsule-source-client-transport.R",
      "test-joint-dp-vector-capsule-client.R")) {
    for (expression in parse(testthat::test_path(file))) {
      if (is.call(expression) &&
          identical(as.character(expression[[1L]]), "test_that")) break
      eval(expression, envir = environment)
    }
  }
  environment
})

.synopsis_public_domain_hash <- function(domain, value) {
  digest::digest(charToRaw(paste0(
    domain, .dsvert_joint_dp_client_json(value))),
  algo = "sha256", serialize = FALSE)
}

.synopsis_public_fixture <- function(
    k = 3L, scaled_values = c("256", "512"), gaussian_artifact = NULL,
    .source = NULL, parent_authorization = NULL) {
  dimension <- if (is.null(gaussian_artifact)) 2L else
    1L + as.integer(gaussian_artifact$coordinate_count)
  if (!is.character(scaled_values) || length(scaled_values) != dimension) {
    stop("scaled_values must cover the complete Synopsis vector")
  }
  source <- if (is.null(.source)) {
    .synopsis_public_helpers$.capsule_source_client_fixture(
      k = k, coordinate_count = dimension, source_count = 1L)
  } else {
    .source
  }
  legacy_context <- .dsvert_joint_dp_client_context(
    source$conns, status = source$status, .aggregate = source$aggregate)
  peers <- legacy_context$servers
  authorities <- legacy_context$designated
  owner <- peers[[1L]]
  pins <- legacy_context$pinset
  legacy_policy <- source$status[[peers[[1L]]]]$policy
  policy <- list(
    version = .DSVERT_CLIENT_SYNOPSIS_POLICY_VERSION,
    privacy_scope = "per_canonical_artifact_v1",
    global_composition_claim = FALSE,
    domain = legacy_policy$domain, cohort_id = legacy_policy$cohort_id,
    ordered_peer_pinset = as.list(pins),
    peer_pinset_sha256 = .dsvert_vector_hash(as.list(pins)),
    peer_count = as.integer(k),
    designated_noise_peers = as.list(authorities),
    artifact_epsilon = 1, artifact_delta = 2^-100,
    adjacency = legacy_policy$adjacency,
    patient_column = legacy_policy$patient_column,
    unit_capacity = as.integer(legacy_policy$unit_capacity),
    fixed_cohort_size = NULL,
    max_records_per_unit = as.integer(legacy_policy$max_records_per_unit),
    overflow_policy = legacy_policy$overflow_policy,
    contingency_unit_aggregation_policy =
      "consistent_cell_else_exclude_v1",
    numeric_grid_bits = 8L,
    primitive_scope = .dsvert_dp_synopsis_manifest_policy_scope_v1(
      source$manifest),
    mechanism_version = "fixture-synopsis-mechanism-v1",
    sticky_release = TRUE, distinct_artifact_gate = FALSE,
    request_limit = FALSE, rate_limit = FALSE, catalog_limit = FALSE)
  sign_outer <- function(value, peer, domain) {
    unsigned <- value
    if (identical(domain, .DSVERT_CLIENT_SYNOPSIS_BIND_DOMAIN)) {
      unsigned$manifest_json <- NULL
      unsigned$artifact_commitments <- NULL
    }
    message <- charToRaw(paste0(
      domain, .dsvert_joint_dp_client_json(unsigned)))
    .synopsis_public_helpers$.capsule_source_client_b64url(
      openssl::ed25519_sign(message, source$keys[[peer]]))
  }
  bootstraps <- stats::setNames(lapply(peers, function(peer) {
    draft <- list(version = "fixture-signed-draft-v1", peer_name = peer)
    if (!is.null(parent_authorization)) {
      draft$parent_authorization <- parent_authorization
    }
    unsigned <- list(
      version = .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_VERSION,
      phase = "custodian_policy_bootstrap", peer_name = peer,
      peer_identity_pk = unname(pins[[peer]]), policy = policy,
      draft = draft,
      data_access = FALSE, patient_derived_metadata = FALSE,
      request_limit = FALSE, rate_limit = FALSE, catalog_limit = FALSE)
    c(unsigned, list(signature = sign_outer(
      unsigned, peer, .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_DOMAIN)))
  }), peers)
  status <- structure(bootstraps, class = c(
    "ds.vertDPSynopsisStatus", "list"))
  context <- list(
    status = status, servers = peers, pinset = pins,
    designated = authorities, conns = source$conns[authorities],
    all_conns = source$conns, policy = policy)
  manifest <- source$manifest

  schema_unsigned <- list(
    version = "synopsis-client-test-schema-v1",
    logical_snapshot = .dsvert_joint_dp_client_canonical(
      manifest$logical_snapshot),
    peer_pinset_sha256 = policy$peer_pinset_sha256)
  schema_message <- charToRaw(paste0(
    .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
    .dsvert_joint_dp_client_json(schema_unsigned)))
  schema_signatures <- stats::setNames(lapply(peers, function(peer) {
    .synopsis_public_helpers$.capsule_source_client_b64url(
      openssl::ed25519_sign(schema_message, source$keys[[peer]]))
  }), peers)
  signed_schema <- c(schema_unsigned, list(signatures = schema_signatures))
  schema_json <- .dsvert_joint_dp_client_json(signed_schema)
  schema_sha256 <- .dsvert_vector_hash(schema_unsigned)

  families <- list(
    admitted_count = list(
      owner_peer = owner, dataset = paste0("data_", owner),
      statistic_maximum = 100, l1_sensitivity = 1),
    numeric_moments = list(artifacts = list()),
    numeric_pair_moments = list(
      artifacts = list(), natural_l1_sensitivity = 0),
    gaussian_models = list(
      artifacts = if (is.null(gaussian_artifact)) list() else
        stats::setNames(list(gaussian_artifact),
                        gaussian_artifact$analysis_id),
      natural_l1_sensitivity = if (is.null(gaussian_artifact)) 0 else
        gaussian_artifact$natural_l1_sensitivity),
    fixed_numeric_histograms = list(artifacts = list(), l1_sensitivity = 0),
    categorical_marginals = list(
      artifacts = if (is.null(gaussian_artifact)) list(flag = list(
        owner_peer = owner, dataset = paste0("data_", owner),
        column = "flag", levels = "yes", statistic_maximum = 100)) else
          list(),
      l1_sensitivity = if (is.null(gaussian_artifact)) 1 else 0),
    categorical_pairs = list(
      sets = list(), cross_artifacts = list(), l1_sensitivity = 0),
    correlation_artifacts = list(), describe_artifacts = list(),
    survival_artifacts = list())
  bits <- 8L
  scale <- 2^bits
  manifest$workload$coordinate_count <- dimension
  manifest$admission$adjacency <- policy$adjacency
  manifest$admission$unit_capacity <- policy$unit_capacity
  manifest$workload$families <- families
  manifest$workload$vertical_crosses <- list()
  manifest$workload$schema_attestation <- list(
    manifest_sha256 = schema_sha256, signers = as.list(peers),
    signatures = schema_signatures)
  manifest$workload$release_lattice <- list(
    version = "biomedical-capsule-common-lattice-v1",
    transform_rule = "raw_coordinate_left_shift_to_common_numeric_grid_v1",
    output_lattice_bits = bits, output_lattice_scale = scale,
    natural_l1_sensitivity = dimension,
    integer_l1_sensitivity_steps = dimension * scale,
    natural_l2_sensitivity = sqrt(dimension),
    integer_l2_sensitivity_steps = sqrt(dimension) * scale)
  manifest$workload$sensitivity <- list(
    l1 = dimension * scale, l2 = sqrt(dimension) * scale)
  manifest$workload$mechanism_selection <- list(version = "test-selection")
  manifest$workload$capsule_mechanism <- c(
    manifest$workload$capsule_mechanism,
    list(
      mechanism = "discrete-laplace", sensitivity_norm = "l1",
      sensitivity = dimension * scale, coordinate_count = dimension,
      clipping_hash = strrep("c", 64L),
      ring_bits = 128L, frac_bits = 0L))
  policy_hash <- .dsvert_vector_hash(policy)
  manifest$capsule_identity$contract <- list(
    protocol = "dsvert-stateless-catalog-synopsis-capsule-identity-v1",
    consortium_id = paste0("dpsc1_", policy_hash),
    policy_contract_hash = policy_hash,
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    logical_snapshot = manifest$logical_snapshot,
    capsule_schema = manifest$capsule_schema,
    admission = manifest$admission, bounds = manifest$bounds,
    workload = manifest$workload,
    privacy_epoch_scope = "per_canonical_artifact_sticky_v1")
  manifest$capsule_identity$capsule_id <-
    .dsvert_vector_hash(manifest$capsule_identity$contract)
  manifest_json <- .dsvert_joint_dp_client_json(manifest)
  manifest_sha256 <- digest::digest(
    manifest_json, algo = "sha256", serialize = FALSE)
  layout <- .dsvert_dp_capsule_vector_layout(manifest)

  catalog <- list(
    version = "dsvert-stateless-catalog-synopsis-catalog-v1",
    domain = policy$domain, cohort_id = policy$cohort_id,
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    logical_snapshot = .dsvert_joint_dp_client_canonical(
      manifest$logical_snapshot),
    capsule_schema = manifest$capsule_schema,
    schema_manifest_sha256 = schema_sha256,
    admission = manifest$admission, bounds = manifest$bounds,
    families = manifest$workload$families,
    vertical_crosses = manifest$workload$vertical_crosses,
    primitive_scope = manifest$workload$primitive_scope,
    release_lattice = manifest$workload$release_lattice,
    sensitivity = manifest$workload$sensitivity,
    coordinate_count = dimension, coordinate_order_sha256 = layout$sha256,
    clipping_sha256 = manifest$workload$capsule_mechanism$clipping_hash)
  rename_ids <- function(value) {
    if (!is.list(value)) return(value)
    if (!is.null(names(value))) {
      names(value)[names(value) == "analysis_id"] <- "catalog_entry_id"
    }
    lapply(value, rename_ids)
  }
  catalog <- rename_ids(catalog)
  projection <- list(
    version = "dsvert-stateless-catalog-synopsis-projection-v1",
    sha256 = .synopsis_public_domain_hash(
      "dsVert/stateless-catalog-synopsis/catalog/v1|", catalog),
    catalog = catalog)

  sensitivity_steps <- as.character(dimension * scale)
  plan <- .synopsis_public_helpers$.vector_client_convolution_plan(
    dimension, sensitivity_steps)
  # Use the same certified dyadic geometric law as the Gaussian fixtures.
  # The generic structural fixture's 1/2^128 stop probability is deliberately
  # tiny and cannot yield a representable utility radius.
  plan$stop_numerator <- "166112941255448185114116409827120804"
  plan$one_geometric_tv_numerator <-
    "6125082604576892348297572878053259906739"
  plan$one_geometric_tv_denominator <- paste0(
    "231584178474632390847141970017375815706539969331281128078915",
    "168015826259279872")
  # Keep the finite sampler-transfer certificate strictly below the declared
  # decimal delta.  The legacy vector fixture uses equality at binary 2^-100,
  # which is fractionally above the synopsis decimal canonicalization.
  safe_delta_denominator <- "2535301200456458802993406410752"
  plan$implementation_delta_denominator <- safe_delta_denominator
  plan$implementation_delta_bound <- "3.944304526105059e-31"
  plan$per_peer_implementation_delta_denominator <- safe_delta_denominator
  plan$per_peer_implementation_delta_bound <- "3.944304526105059e-31"
  full_plan_sha256 <- .dsvert_vector_hash(plan)
  profile <- .dsvert_vector_profile(
    manifest$workload$capsule_mechanism,
    backend = .DSVERT_CLIENT_VECTOR_BACKEND)
  profile_value <- list(
    version = "dsvert-stateless-catalog-synopsis-backend-profile-v1",
    mechanism = "discrete-laplace", backend = profile$backend,
    plan_version = profile$plan_version, sampler = profile$sampler,
    release_mechanism = profile$release_mechanism, draw_count = 2L,
    complete_epsilon_per_draw = TRUE,
    delta_aggregation = profile$delta_aggregation,
    adversary_model = paste0(
      "analyst_plus_at_most_one_of_two_noncolluding_",
      "noise_authorities_v1"),
    output_transform = profile$postprocessing,
    commitment_purpose = "convolution")
  selection <- list(
    version = "dsvert-stateless-catalog-synopsis-backend-selection-v1",
    rule = "public_coordinate_ceiling_v1",
    selected_before_private_material = TRUE,
    retry_may_change_backend = FALSE,
    policy_version =
      .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_COST_POLICY_VERSION,
    total_coordinate_count = dimension, maximum_promoted_coordinates = 1L,
    promoted = FALSE, backend = profile$backend,
    selection_reason = "above_public_exact_gc_cost_ceiling")
  lattice_info <- .dsvert_dp_synopsis_client_lattice(manifest, layout)
  lattice_transform <- list(
    version = manifest$workload$release_lattice$version,
    transform_rule = manifest$workload$release_lattice$transform_rule,
    output_lattice_bits = bits,
    scale_shifts = as.list(lattice_info$scale_shifts),
    raw_upper_bounds = as.list(lattice_info$raw_upper_bounds),
    sensitivity_norm = "l1", sensitivity_steps = sensitivity_steps)
  lattice <- list(
    version = "dsvert-stateless-catalog-synopsis-lattice-v1",
    coordinate_count = dimension, coordinate_order_sha256 = layout$sha256,
    clipping_sha256 = strrep("c", 64L),
    transform_sha256 = .dsvert_vector_hash(lattice_transform),
    output_lattice_bits = bits, output_lattice_scale = scale,
    sensitivity_norm = "l1", sensitivity_steps = sensitivity_steps,
    ring_bits = 128L, fractional_bits = 0L)
  epsilon <- formatC(1, digits = 18L, format = "e", decimal.mark = ".")
  delta <- formatC(2^-100, digits = 18L, format = "e", decimal.mark = ".")
  request <- list(
    epsilon = epsilon, delta = delta, sensitivity_steps = sensitivity_steps,
    total_coordinate_count = dimension)
  convolution_fields <- c(
    "version", "sampler", "stop_bits", "stop_numerator", "uniform_bits",
    "binary_geometric_bits", "bernoulli_thresholds", "sensitivity_steps",
    "total_coordinate_count", "epsilon_effective_upper_numerator",
    "epsilon_effective_upper_denominator", "one_geometric_tv_numerator",
    "one_geometric_tv_denominator", "tail_upper_numerator",
    "tail_upper_denominator", "rounding_upper_numerator",
    "rounding_upper_denominator", "implementation_delta_numerator",
    "implementation_delta_denominator", "maximum_noise_magnitude",
    "independent_noise_peer_count", "complete_epsilon_per_peer",
    "epsilon_divided_by_peer_count",
    "geometric_variables_per_peer_per_coordinate",
    "geometric_variables_total_per_coordinate",
    "per_peer_implementation_delta_numerator",
    "per_peer_implementation_delta_denominator",
    "release_implementation_delta_aggregation",
    "two_peer_ideal_transfer_delta_numerator",
    "two_peer_ideal_transfer_delta_denominator")
  draw_law <- plan[convolution_fields]
  physical_identity <- list(
    version = "dsvert-stateless-catalog-synopsis-physical-plan-v1",
    request = request, profile = profile_value, lattice = lattice,
    backend_selection = selection, draw_law = draw_law,
    draw_law_sha256 = .synopsis_public_domain_hash(
      "dsVert/stateless-catalog-synopsis/draw-law/v1|", draw_law))
  semantic <- list(
    version = "dsvert-analysis-semantic-stateless-catalog-synopsis-v1",
    catalog_projection = projection,
    source_claim_set_sha256 = strrep("a", 64L),
    noise_authority_roles = list(
      version = "dsvert-synopsis-noise-authority-roles-v1",
      role_order = list(
        "primary_noise_authority", "secondary_noise_authority"),
      authority_ids = as.list(unname(context$pinset[authorities]))),
    privacy = list(
      version = "dsvert-per-synopsis-dp-v1",
      adjacency = policy$adjacency, privacy_unit = "patient",
      epsilon = epsilon, delta = delta,
      mechanism = list(
        version = profile$release_mechanism, family = "discrete_laplace",
        sensitivity = list(
          version = "dsvert-sensitivity-v1", norm = "l1",
          steps = sensitivity_steps),
        randomness = list(
          version = "dsvert-randomness-plan-v1",
          lanes = list(final_noise = list(
            version = "dsvert-randomness-lane-v1",
            purpose = "privatize_final_vector", primitive = profile$sampler,
            coordinates = dimension))))),
    release = physical_identity,
    public_shape = list(
      version = "dsvert-stateless-catalog-synopsis-shape-v1",
      coordinates = dimension))
  artifact_value <- list(
    semantic = semantic,
    artifact_key = .synopsis_public_domain_hash(
      "dsVert/analysis-artifact-key/v1|", semantic),
    physical_plan = c(physical_identity, list(
      full_plan = plan, full_plan_sha256 = full_plan_sha256)))

  compile_receipts <- stats::setNames(lapply(peers, function(peer) {
    unsigned <- list(
      version = "dsvert-stateless-catalog-synopsis-compile-receipt-v1",
      peer_name = peer, peer_identity_pk = unname(context$pinset[[peer]]),
      manifest_sha256 = manifest_sha256,
      artifact_key = artifact_value$artifact_key,
      source_claim_set_sha256 = semantic$source_claim_set_sha256,
      full_plan_sha256 = full_plan_sha256)
    message <- charToRaw(paste0(
      "dsVert/stateless-catalog-synopsis/compile-receipt/v1|",
      .dsvert_joint_dp_client_json(unsigned)))
    c(unsigned, list(signature =
      .synopsis_public_helpers$.capsule_source_client_b64url(
        openssl::ed25519_sign(message, source$keys[[peer]]))))
  }), peers)
  compile_unsigned <- unname(lapply(compile_receipts, function(value) {
    value[setdiff(names(value), "signature")]
  }))
  compilation <- list(
    version = "dsvert-stateless-catalog-synopsis-compile-receipt-set-v1",
    artifact = artifact_value, receipts = compile_receipts,
    receipt_set_sha256 = .synopsis_public_domain_hash(
      "dsVert/stateless-catalog-synopsis/compile-receipt-set/v1|",
      list(
        version = "dsvert-stateless-catalog-synopsis-compile-receipt-set-v1",
        receipts = compile_unsigned)))

  artifact_index <- .dsvert_dp_capsule_artifact_commitment_index_client(
    manifest, .dsvert_dp_synopsis_policy_view_v1(policy), manifest_sha256)
  workload_contract_json <- .dsvert_joint_dp_client_json(list(
    version = "synopsis-client-test-workload-v1"))
  workload_contract_sha256 <- digest::digest(
    workload_contract_json, algo = "sha256", serialize = FALSE)
  bootstrap_set_sha256 <- .dsvert_vector_hash(list(
    version = .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_VERSION,
    bootstraps = unname(unclass(status))))
  bound_manifests <- stats::setNames(lapply(peers, function(peer) {
    unsigned <- list(
      version = .DSVERT_CLIENT_SYNOPSIS_BOUND_MANIFEST_VERSION,
      phase = "server_authoritative_synopsis_manifest_memoized",
      peer_name = peer, peer_identity_pk = unname(context$pinset[[peer]]),
      peer_pinset_sha256 = policy$peer_pinset_sha256,
      bootstrap_set_sha256 = bootstrap_set_sha256,
      schema_sha256 = schema_sha256,
      workload_contract_sha256 = workload_contract_sha256,
      manifest_sha256 = manifest_sha256,
      manifest_bytes = nchar(manifest_json, type = "bytes"),
      capsule_id = manifest$capsule_identity$capsule_id,
      artifact_commitment_count = artifact_index$count,
      artifact_commitments_root = artifact_index$root,
      artifact_commitments = artifact_index$value,
      manifest_json = manifest_json,
      privacy_scope = "per_canonical_artifact_v1",
      global_composition_claim = FALSE,
      durable_memoization = TRUE, deterministic_replay = TRUE,
      data_access = FALSE, request_limit = FALSE, rate_limit = FALSE,
      catalog_limit = FALSE)
    c(unsigned, list(signature = sign_outer(
      unsigned, peer, .DSVERT_CLIENT_SYNOPSIS_BIND_DOMAIN)))
  }), peers)
  build_receipts <- stats::setNames(lapply(bound_manifests, function(value) {
    value[setdiff(names(value), c("manifest_json", "artifact_commitments"))]
  }), peers)
  manifest_bundle <- list(
    schema_json = schema_json, schema_sha256 = schema_sha256,
    workload_contract_json = workload_contract_json,
    workload_contract_sha256 = workload_contract_sha256,
    logical_snapshot = .dsvert_joint_dp_client_canonical(
      manifest$logical_snapshot),
    manifest_json = manifest_json, manifest_sha256 = manifest_sha256,
    capsule_id = manifest$capsule_identity$capsule_id,
    artifact_commitments = artifact_index$value,
    artifact_commitment_count = artifact_index$count,
    artifact_commitments_root = artifact_index$root,
    manifest_build_receipts = build_receipts,
    manifest_signatures = lapply(bound_manifests, `[[`, "signature"),
    deterministic_replay = TRUE, operation_or_request_limit = FALSE,
    history_can_deny_operation = FALSE, context = context)

  base_contract <- list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CONTRACT_VERSION,
    purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE,
    capsule_id = manifest$capsule_identity$capsule_id,
    logical_snapshot_sha256 = .dsvert_vector_hash(manifest$logical_snapshot),
    workload_sha256 = .dsvert_vector_hash(manifest$workload),
    source_context_hash = manifest$workload$capsule_mechanism$source_context_hash,
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    source_peers = list(owner), designated_noise_peers = as.list(authorities),
    coordinate_count = dimension, coordinate_order_sha256 = layout$sha256,
    ring_bits = 128L, record_bytes = 16L,
    record_encoding = "little_endian_unsigned_fixed_16_bytes",
    chunk_coordinates = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES,
    chunk_count = 1L, chunk_shape = "fixed_public_coordinate_slices",
    history_gate = FALSE, ready_for_sampling = FALSE)
  binding <- list(
    version = "dsvert-stateless-catalog-synopsis-source-contract-v1",
    manifest_capsule_id = base_contract$capsule_id,
    artifact_key = artifact_value$artifact_key,
    source_claim_set_sha256 = semantic$source_claim_set_sha256)
  base_contract$capsule_id <- .synopsis_public_domain_hash(
    "dsVert/stateless-catalog-synopsis/source-namespace/v1|", binding)
  base_contract$synopsis_binding <- binding
  source_contract_sha256 <- .dsvert_vector_hash(base_contract)
  execution_id <- digest::digest(charToRaw(paste0(
    "dsVert/stateless-catalog-synopsis/execution-id/v1|",
    artifact_value$artifact_key)), algo = "sha256", serialize = FALSE)
  support <- openssl::bignum(plan$maximum_noise_magnitude) * 2L
  positive_limit <- (openssl::bignum(2) ^ 127L) - 1L
  maximum_scaled_source <- openssl::bignum(0)
  for (index in seq_len(dimension)) {
    candidate <- openssl::bignum(lattice_info$raw_upper_bounds[[index]]) *
      (openssl::bignum(2) ^ lattice_info$scale_shifts[[index]])
    if (candidate > maximum_scaled_source) maximum_scaled_source <- candidate
  }
  ring_unsigned <- list(
    version = "dsvert-stateless-catalog-synopsis-ring128-certificate-v1",
    ring_bits = 128L, fractional_bits = 0L,
    maximum_scaled_source_coordinate = as.character(maximum_scaled_source),
    maximum_release_noise_magnitude = as.character(support),
    positive_limit = as.character(positive_limit), no_wrap_certified = TRUE)
  ring <- c(ring_unsigned, list(sha256 = .synopsis_public_domain_hash(
    "dsVert/stateless-catalog-synopsis/ring128-certificate/v1|",
    ring_unsigned)))
  contract_value <- list(
    version = "dsvert-stateless-catalog-synopsis-execution-v1",
    execution_id = execution_id, artifact_key = artifact_value$artifact_key,
    source_claim_set_sha256 = semantic$source_claim_set_sha256,
    semantic_sha256 = .dsvert_vector_hash(semantic),
    draw_law_sha256 = physical_identity$draw_law_sha256,
    authority_roles = semantic$noise_authority_roles,
    authority_peers = as.list(authorities),
    geometry = list(
      coordinate_count = dimension, public_chunk_coordinates = dimension,
      public_chunk_count = 1L, coordinate_order_sha256 = layout$sha256,
      lattice_transform_sha256 = lattice$transform_sha256),
    mechanism = list(
      mechanism = "discrete-laplace", backend = profile$backend,
      sampler = profile$sampler,
      release_mechanism = profile$release_mechanism), ring = ring)
  contract_sha256 <- .synopsis_public_domain_hash(
    "dsVert/stateless-catalog-synopsis/execution-contract/v1|",
    contract_value)
  attempt_value <- list(
    version = "dsvert-stateless-catalog-synopsis-execution-attempt-v1",
    execution_id = execution_id, artifact_key = artifact_value$artifact_key,
    manifest_sha256 = manifest_sha256,
    manifest_capsule_id = manifest$capsule_identity$capsule_id,
    source_contract_sha256 = source_contract_sha256,
    full_plan_sha256 = full_plan_sha256,
    execution_geometry = list(
      chunk_coordinates = dimension, chunk_count = 1L))
  attempt_sha256 <- .synopsis_public_domain_hash(
    "dsVert/stateless-catalog-synopsis/execution-attempt/v1|",
    attempt_value)
  result_set_sha256 <- strrep("b", 64L)
  public <- list(
    version = "dsvert-stateless-catalog-synopsis-public-chunk-v1",
    artifact_key = artifact_value$artifact_key, execution_id = execution_id,
    contract_sha256 = contract_sha256, attempt_sha256 = attempt_sha256,
    result_set_sha256 = result_set_sha256, public_chunk_index = 0L,
    public_chunk_count = 1L, coordinate_offset = 0L,
    coordinate_count = dimension,
    output_lattice_bits = bits, output_lattice_scale = scale,
    scaled_values = as.list(scaled_values),
    value_encoding = "nonnegative-decimal-integer-common-lattice-v1",
    postprocessing = profile$postprocessing,
    source_values_exposed = FALSE, preclamp_values_exposed = FALSE)
  chunk_sha256 <- .synopsis_public_domain_hash(
    "dsVert/stateless-catalog-synopsis/public-chunk/v1|", public)
  root <- .dsvert_vector_merkle_root(chunk_sha256)
  release_values <- stats::setNames(lapply(seq_along(authorities), function(i) {
    peer <- authorities[[i]]
    unsigned <- list(
      version = "dsvert-stateless-catalog-synopsis-release-v1",
      phase = "synopsis_released", execution_id = execution_id,
      artifact_key = artifact_value$artifact_key,
      contract_sha256 = contract_sha256, attempt_sha256 = attempt_sha256,
      source_contract_sha256 = source_contract_sha256,
      result_set_sha256 = result_set_sha256,
      local_authority = list(
        peer_name = peer, identity_pk = unname(context$pinset[[peer]]),
        role = c("primary_noise_authority",
                 "secondary_noise_authority")[[i]]),
      public_chunk_count = 1L, final_chunk_hashes = list(chunk_sha256),
      final_vector_root = root, output_lattice_bits = bits,
      output_lattice_scale = as.character(scale),
      mechanism = profile$release_mechanism,
      epsilon = epsilon, delta = delta,
      implementation_delta_numerator =
        plan$per_peer_implementation_delta_numerator,
      implementation_delta_denominator =
        plan$per_peer_implementation_delta_denominator,
      delta_aggregation = profile$delta_aggregation,
      postprocessing = profile$postprocessing,
      all_public_chunks_durable = TRUE,
      intermediate_payload_exposed = FALSE, durable_replay = TRUE,
      capability_available = TRUE)
    message <- charToRaw(paste0(
      "dsVert/stateless-catalog-synopsis/release/v1|",
      .dsvert_joint_dp_client_json(unsigned)))
    c(unsigned, list(signature =
      .synopsis_public_helpers$.capsule_source_client_b64url(
        openssl::ed25519_sign(message, source$keys[[peer]]))))
  }), authorities)
  release_receipts <- lapply(release_values, .dsvert_joint_dp_client_json)
  replay <- list(
    version = "dsvert-stateless-catalog-synopsis-replay-v1",
    phase = "synopsis_public_chunk_replayed", execution_id = execution_id,
    artifact_key = artifact_value$artifact_key,
    contract_sha256 = contract_sha256, attempt_sha256 = attempt_sha256,
    source_contract_sha256 = source_contract_sha256,
    result_set_sha256 = result_set_sha256, final_vector_root = root,
    public_chunk_index = 0L, public_chunk_count = 1L,
    chunk_sha256 = chunk_sha256, chunk = public, merkle_proof = list(),
    durable_replay = TRUE, source_store_read = FALSE,
    sampler_invoked = FALSE, finalizer_invoked = FALSE,
    transport_read = FALSE)
  replay_json <- .dsvert_joint_dp_client_json(replay)
  replay_responses <- list(`0` = stats::setNames(
    rep(list(replay_json), 2L), authorities))
  list(
    release_receipts = release_receipts,
    replay_responses = replay_responses,
    manifest_bundle = manifest_bundle, status = context$status,
    artifact = compilation, values = scaled_values,
    authorities = authorities, source = source)
}

test_that("synopsis public vector has a closed trusted-input ABI", {
  expect_identical(
    names(formals(.dsvert_dp_synopsis_public_vector_v1)),
    c("release_receipts", "replay_responses", "manifest_bundle",
      "status", "artifact"))
  expect_error(.dsvert_dp_synopsis_public_vector_v1(
    list(), list(), list(manifest_json = "{}"), list(), list()),
  "manifest|trusted|context|bundle", ignore.case = TRUE)
  expect_identical(
    .dsvert_dp_synopsis_client_lattice_scale(62L),
    "4611686018427387904")
})

test_that("synopsis evidence uses an application bound above portable DSI", {
  expect_identical(
    .DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES,
    32L * 1024L^2)
  expect_identical(
    .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES,
    2L * 1024L^2)
  expect_identical(
    .DSVERT_CLIENT_SYNOPSIS_PREPARE_MAX_OBJECT_BYTES,
    64L * 1024L)

  padding_bytes <- 8L * 1024L^2 + 4096L
  wire <- paste0('{"padding":"', strrep("A", padding_bytes), '"}')
  parsed <- .dsvert_dp_synopsis_client_json(wire, "large synopsis evidence")
  expect_equal(nchar(parsed$padding, type = "bytes"), padding_bytes)
  expect_error(
    .dsvert_dp_synopsis_client_json(
      wire, "oversized synopsis receipt",
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES),
    "Invalid joint-DP oversized synopsis receipt")

  expression <- .dsvert_dsi_text_frame_call(call(
    name = "dsvertDPSynopsisCompileDS",
    manifest_sha256 = strrep("1", 64L), claims_json = wire))
  condition <- tryCatch(
    .dsvert_validate_dsi_expression_sizes(
      expression, capacity_bytes = 8L * 1024L^2),
    error = identity)
  expect_s3_class(condition, "dsvert_resource_oversize")
  expect_identical(condition$scope, "DataSHIELD negotiated expression")
})

test_that("synopsis RELEASE and bilateral REPLAY become one exact public vector", {
  fixture <- .synopsis_public_fixture()
  value <- do.call(.dsvert_dp_synopsis_public_vector_v1, fixture[c(
    "release_receipts", "replay_responses", "manifest_bundle",
    "status", "artifact")])
  expect_s3_class(value, "dsvert_joint_dp_vector")
  expect_identical(value$values, c(1, 2))
  expect_identical(value$scaled_values, fixture$values)
  expect_identical(
    value$signed_provenance$version,
    "dsvert-stateless-synopsis-public-provenance-v1")
  expect_false(
    value$signed_provenance$execution_hashes_client_reconstructed)
  expect_false(any(c(
    "allocation_certificate", "release_instance", "prepare_receipts",
    "finalization_receipts") %in% names(value$signed_provenance)))
  expect_false(any(c(
    "capsule_id", "release_instance_id", "privacy_epoch", "noise_key_id") %in%
      names(value)))

  reversed <- fixture
  reversed$release_receipts <- rev(reversed$release_receipts)
  reversed$replay_responses[[1L]] <- rev(reversed$replay_responses[[1L]])
  expect_identical(
    do.call(.dsvert_dp_synopsis_public_vector_v1, reversed[c(
      "release_receipts", "replay_responses", "manifest_bundle",
      "status", "artifact")]), value)
})

test_that(paste(
  "unrelated signed parent metadata preserves manifest artifact and release",
  "for K=2/3/5"), {
  arguments <- c(
    "release_receipts", "replay_responses", "manifest_bundle",
    "status", "artifact")
  parent_authorization <- list(
    schema_sha256 = strrep("7", 64L),
    workload_contract_sha256 = strrep("8", 64L),
    catalog = list(
      unrelated_column = "stage",
      unrelated_pair = c("disease", "stage"),
      unrelated_gaussian = list(
        outcome = "gaussian_y", predictors = "gaussian_x")))
  for (k in c(2L, 3L, 5L)) {
    selected_only <- .synopsis_public_fixture(k = k)
    expanded_parent <- .synopsis_public_fixture(
      k = k, .source = selected_only$source,
      parent_authorization = parent_authorization)
    selected_release <- do.call(
      .dsvert_dp_synopsis_public_vector_v1,
      selected_only[arguments])
    expanded_release <- do.call(
      .dsvert_dp_synopsis_public_vector_v1,
      expanded_parent[arguments])

    expect_identical(
      expanded_parent$manifest_bundle$manifest_sha256,
      selected_only$manifest_bundle$manifest_sha256,
      info = paste("manifest K =", k))
    expect_identical(
      expanded_release$artifact_key, selected_release$artifact_key,
      info = paste("artifact K =", k))
    expect_identical(
      expanded_release$final_vector_root,
      selected_release$final_vector_root,
      info = paste("release K =", k))
    expect_identical(
      expanded_release$scaled_values, selected_release$scaled_values,
      info = paste("draw K =", k))
    expect_identical(
      expanded_parent$release_receipts,
      selected_only$release_receipts,
      info = paste("signed RELEASE K =", k))
  }
})

test_that("synopsis public vector keeps K witnesses out of the release pair", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .synopsis_public_fixture(k = k)
    value <- do.call(.dsvert_dp_synopsis_public_vector_v1, fixture[c(
      "release_receipts", "replay_responses", "manifest_bundle",
      "status", "artifact")])
    expect_identical(value$values, c(1, 2))
    expect_length(value$signed_provenance$ordered_peer_pinset, k)
    expect_identical(
      value$signed_provenance$designated_noise_peers,
      as.list(fixture$authorities))
  }
})

test_that("synopsis public vector fails closed on consensus and wire tampering", {
  base <- .synopsis_public_fixture()
  reject <- function(change, pattern) {
    candidate <- change(base)
    expect_error(do.call(.dsvert_dp_synopsis_public_vector_v1, candidate[c(
      "release_receipts", "replay_responses", "manifest_bundle",
      "status", "artifact")]), pattern, ignore.case = TRUE)
  }
  reject(function(x) {
    x$artifact$receipts <- x$artifact$receipts[-1L]
    x
  }, "compile|receipt|coverage|peer")
  reject(function(x) {
    receipt <- x$manifest_bundle$manifest_build_receipts[[1L]]
    prefix <- if (startsWith(receipt$signature, "A")) "B" else "A"
    receipt$signature <- paste0(prefix, substring(receipt$signature, 2L))
    x$manifest_bundle$manifest_build_receipts[[1L]] <- receipt
    x$manifest_bundle$manifest_signatures[[1L]] <- receipt$signature
    x
  }, "manifest|signature")
  reject(function(x) {
    x$artifact$artifact$semantic$public_shape$coordinates <- 3L
    x
  }, "artifact|shape|key")
  reject(function(x) {
    x$artifact$artifact$physical_plan$backend_selection$backend <-
      .DSVERT_CLIENT_VECTOR_EXACT_BACKEND
    x
  }, "artifact|backend|selection|physical plan")
  reject(function(x) {
    x$artifact$artifact$semantic$privacy$mechanism$randomness$lanes$
      final_noise$purpose <- "adaptive_retry_noise"
    x
  }, "artifact|privacy|key")
  reject(function(x) {
    x$artifact$artifact$physical_plan$draw_law$maximum_noise_magnitude <-
      "1"
    x$artifact$artifact$physical_plan$draw_law_sha256 <-
      .synopsis_public_domain_hash(
        "dsVert/stateless-catalog-synopsis/draw-law/v1|",
        x$artifact$artifact$physical_plan$draw_law)
    x$artifact$artifact$semantic$release <-
      x$artifact$artifact$physical_plan[c(
        "version", "request", "profile", "lattice", "backend_selection",
        "draw_law", "draw_law_sha256")]
    x$artifact$artifact$artifact_key <- .synopsis_public_domain_hash(
      "dsVert/analysis-artifact-key/v1|",
      x$artifact$artifact$semantic)
    x
  }, "draw law|full plan|compile|artifact")
  reject(function(x) {
    receipt <- .dsvert_joint_dp_client_decode(
      x$release_receipts[[1L]], "test release", 2L * 1024L^2)
    prefix <- if (startsWith(receipt$signature, "A")) "B" else "A"
    receipt$signature <- paste0(prefix, substring(receipt$signature, 2L))
    x$release_receipts[[1L]] <- .dsvert_joint_dp_client_json(receipt)
    x
  }, "signature|release")
  reject(function(x) {
    x$release_receipts <- lapply(x$release_receipts, function(json) {
      receipt <- .dsvert_joint_dp_client_decode(
        json, "test release", 2L * 1024L^2)
      receipt$final_chunk_hashes <- receipt$final_chunk_hashes[[1L]]
      peer <- receipt$local_authority$peer_name
      unsigned <- receipt[setdiff(names(receipt), "signature")]
      receipt$signature <-
        .synopsis_public_helpers$.capsule_source_client_b64url(
          openssl::ed25519_sign(charToRaw(paste0(
            "dsVert/stateless-catalog-synopsis/release/v1|",
            .dsvert_joint_dp_client_json(unsigned))),
          x$source$keys[[peer]]))
      .dsvert_joint_dp_client_json(receipt)
    })
    x
  }, "release|chunk")
  reject(function(x) {
    x$release_receipts <- lapply(x$release_receipts, function(json) {
      receipt <- .dsvert_joint_dp_client_decode(
        json, "test release", 2L * 1024L^2)
      receipt$final_chunk_hashes <- list(receipt$final_chunk_hashes)
      peer <- receipt$local_authority$peer_name
      unsigned <- receipt[setdiff(names(receipt), "signature")]
      receipt$signature <-
        .synopsis_public_helpers$.capsule_source_client_b64url(
          openssl::ed25519_sign(charToRaw(paste0(
            "dsVert/stateless-catalog-synopsis/release/v1|",
            .dsvert_joint_dp_client_json(unsigned))),
          x$source$keys[[peer]]))
      .dsvert_joint_dp_client_json(receipt)
    })
    x
  }, "release|chunk")
  reject(function(x) {
    replay <- .dsvert_joint_dp_client_decode(
      x$replay_responses[[1L]][[2L]], "test replay", 2L * 1024L^2)
    replay$transport_read <- TRUE
    x$replay_responses[[1L]][[2L]] <-
      .dsvert_joint_dp_client_json(replay)
    x
  }, "replay|different|bilateral")
  reject(function(x) {
    x$status[[1L]]$policy$unit_capacity <-
      x$status[[1L]]$policy$unit_capacity + 1
    x
  }, "status|policy|consensus|trusted")
})

test_that("synopsis integer bounds and exact client domain are enforced", {
  bound <- .synopsis_public_fixture(scaled_values = c("25601", "512"))
  expect_error(do.call(.dsvert_dp_synopsis_public_vector_v1, bound[c(
    "release_receipts", "replay_responses", "manifest_bundle",
    "status", "artifact")]), "bound|coordinate")

  large <- .synopsis_public_fixture(
    scaled_values = c("9007199254740992", "512"))
  expect_error(do.call(.dsvert_dp_synopsis_public_vector_v1, large[c(
    "release_receipts", "replay_responses", "manifest_bundle",
    "status", "artifact")]), "exact client|integer domain|bound")

  noncanonical <- .synopsis_public_fixture(scaled_values = c("0256", "512"))
  expect_error(do.call(.dsvert_dp_synopsis_public_vector_v1, noncanonical[c(
    "release_receipts", "replay_responses", "manifest_bundle",
    "status", "artifact")]), "integer|REPLAY|coordinate")

  scalar <- .synopsis_public_fixture()
  replay <- .dsvert_joint_dp_client_decode(
    scalar$replay_responses[[1L]][[1L]], "test replay", 2L * 1024L^2)
  replay$chunk$scaled_values <- replay$chunk$scaled_values[[1L]]
  replay$chunk_sha256 <- .synopsis_public_domain_hash(
    "dsVert/stateless-catalog-synopsis/public-chunk/v1|", replay$chunk)
  replay$final_vector_root <- .dsvert_vector_merkle_root(
    replay$chunk_sha256)
  replay_json <- .dsvert_joint_dp_client_json(replay)
  scalar$replay_responses[[1L]] <- stats::setNames(
    rep(list(replay_json), 2L), scalar$authorities)
  scalar$release_receipts <- lapply(scalar$release_receipts, function(json) {
    receipt <- .dsvert_joint_dp_client_decode(
      json, "test release", 2L * 1024L^2)
    receipt$final_chunk_hashes <- list(replay$chunk_sha256)
    receipt$final_vector_root <- replay$final_vector_root
    peer <- receipt$local_authority$peer_name
    unsigned <- receipt[setdiff(names(receipt), "signature")]
    receipt$signature <-
      .synopsis_public_helpers$.capsule_source_client_b64url(
        openssl::ed25519_sign(charToRaw(paste0(
          "dsVert/stateless-catalog-synopsis/release/v1|",
          .dsvert_joint_dp_client_json(unsigned))),
        scalar$source$keys[[peer]]))
    .dsvert_joint_dp_client_json(receipt)
  })
  expect_error(do.call(.dsvert_dp_synopsis_public_vector_v1, scalar[c(
    "release_receipts", "replay_responses", "manifest_bundle",
    "status", "artifact")]), "REPLAY|chunk|scaled", ignore.case = TRUE)
})

test_that("synopsis lattice recycles scalar categorical bounds", {
  fixture <- .synopsis_public_fixture()
  manifest <- .dsvert_joint_dp_client_decode(
    fixture$manifest_bundle$manifest_json, "test synopsis manifest",
    8L * 1024L^2)
  manifest$workload$families$categorical_marginals$artifacts$flag$levels <-
    c("yes", "no")
  manifest$workload$coordinate_count <- 3L
  manifest$workload$capsule_mechanism$coordinate_count <- 3L
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  lattice <- .dsvert_dp_synopsis_client_lattice(manifest, layout)
  expect_identical(layout$coordinate_count, 3L)
  expect_identical(lattice$raw_upper_bounds, rep("100", 3L))
})

test_that("synopsis profile reconstruction keeps manifest calibration", {
  fixture <- .synopsis_public_fixture()
  expected <- .dsvert_joint_dp_client_decode(
    fixture$manifest_bundle$manifest_json, "test synopsis manifest",
    8L * 1024L^2)$workload$mechanism_selection
  observed <- NULL
  original <- .dsvert_vector_profile
  testthat::local_mocked_bindings(
    .dsvert_vector_profile = function(
        capsule_mechanism, mechanism_selection = NULL, backend = NULL) {
      observed <<- mechanism_selection
      original(capsule_mechanism, mechanism_selection, backend)
    },
    .package = "dsVertClient")
  do.call(.dsvert_dp_synopsis_public_vector_v1, fixture[c(
    "release_receipts", "replay_responses", "manifest_bundle",
    "status", "artifact")])
  expect_identical(observed, expected)
})

test_that("synopsis accepts canonical server release field ordering", {
  fixture <- .synopsis_public_fixture()
  fixture$artifact$artifact$semantic$release <-
    .dsvert_joint_dp_client_canonical(
      fixture$artifact$artifact$semantic$release)
  value <- do.call(.dsvert_dp_synopsis_public_vector_v1, fixture[c(
    "release_receipts", "replay_responses", "manifest_bundle",
    "status", "artifact")])
  expect_identical(value$values, c(1, 2))
})

test_that("client rejects zero delta certificates from finite Laplace v3", {
  fixture <- .synopsis_public_fixture()
  trusted <- .dsvert_dp_synopsis_client_bundle(
    fixture$manifest_bundle, fixture$status)
  zero <- "0.000000000000000000e+00"
  expect_true(.dsvert_dp_synopsis_client_fraction_leq(
    "0", "1", zero, "implementation delta"))
  expect_false(.dsvert_dp_synopsis_client_fraction_leq(
    "1", "1000000000000000000000000000000", zero,
    "implementation delta"))
  long_denominator <- paste0("1", strrep("0", 255L))
  expect_true(.dsvert_dp_synopsis_client_fraction_leq(
    "1", long_denominator, "1.000000000000000000e+00",
    "implementation delta"))
  expect_error(.dsvert_dp_synopsis_client_fraction_leq(
    "1", paste0("1", strrep("0", 512L)),
    "1.000000000000000000e+00", "implementation delta"),
    "Invalid synopsis implementation delta")
  for (peer in trusted$context$servers) {
    trusted$status[[peer]]$policy$capsule_delta <- 0
  }
  compilation <- fixture$artifact
  physical <- compilation$artifact$physical_plan
  physical$request$delta <- zero
  physical$draw_law$implementation_delta_numerator <- "0"
  physical$draw_law$implementation_delta_denominator <- "1"
  physical$draw_law$per_peer_implementation_delta_numerator <- "0"
  physical$draw_law$per_peer_implementation_delta_denominator <- "1"
  physical$full_plan$implementation_delta_numerator <- "0"
  physical$full_plan$implementation_delta_denominator <- "1"
  physical$full_plan$implementation_delta_bound <- "0"
  physical$full_plan$per_peer_implementation_delta_numerator <- "0"
  physical$full_plan$per_peer_implementation_delta_denominator <- "1"
  physical$full_plan$per_peer_implementation_delta_bound <- "0"
  physical$full_plan_sha256 <- .dsvert_vector_hash(physical$full_plan)
  physical$draw_law_sha256 <- .synopsis_public_domain_hash(
    "dsVert/stateless-catalog-synopsis/draw-law/v1|",
    physical$draw_law)
  compilation$artifact$physical_plan <- physical
  compilation$artifact$semantic$privacy$delta <- zero
  compilation$artifact$semantic$release <- physical[setdiff(
    names(physical), c("full_plan", "full_plan_sha256"))]
  compilation$artifact$artifact_key <- .synopsis_public_domain_hash(
    "dsVert/analysis-artifact-key/v1|",
    compilation$artifact$semantic)

  expect_error(
    .dsvert_dp_synopsis_client_compile(
      compilation, trusted, fixture$manifest_bundle),
    "pure DP|finite Laplace|positive implementation delta")
})
