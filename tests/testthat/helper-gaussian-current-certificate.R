# Current Gaussian certificate fixture. RELEASE/REPLAY signatures, cross-owner
# attestations, descriptors, layouts and result binding use production readers.
# Only authenticated bundle/compilation reconstruction and public accuracy
# calibration are supplied as test seams; this is not a federation test.
.gaussian_current_fixture <- function(k = 2L, artifact = NULL,
                                      coordinates = NULL) {
  legacy <- new.env(parent = environment())
  for (node in parse(testthat::test_path("test-dp-gaussian.R"))) {
    if (is.call(node) && identical(node[[1L]], as.name("<-"))) eval(node, legacy)
  }
  original <- legacy$.dp_gaussian_cross_frontdoor_fixture(k)
  if (is.null(artifact)) artifact <- original$artifact
  if (is.null(coordinates)) coordinates <- c(80, 80, 40, 30, 40, 24, 30)
  peers <- names(original$conns)
  authorities <- peers[1:2]
  keys <- stats::setNames(lapply(peers, function(peer) openssl::ed25519_keygen()), peers)
  b64 <- function(value) sub("=+$", "", chartr("+/", "-_",
    gsub("[\r\n]", "", jsonlite::base64_enc(value))), perl = TRUE)
  pins <- vapply(keys, function(key) b64(tail(as.raw(as.list(key)$pubkey), 32L)), character(1L))
  scale <- 2^as.numeric(artifact$numeric_grid_bits)
  capacity <- as.numeric(artifact$transcript$padded_units)
  if (!length(capacity)) capacity <- as.numeric(artifact$statistic_maximum[[1L]])
  manifest <- original$manifest
  manifest$admission <- list(unit_capacity = capacity, adjacency = artifact$adjacency)
  manifest$bounds <- list(numeric_grid_bits = artifact$numeric_grid_bits)
  manifest$workload$coordinate_count <- 1 + length(coordinates)
  manifest$workload$families$admitted_count$statistic_maximum <- capacity
  manifest$workload$families$gaussian_models$artifacts <- stats::setNames(list(artifact), artifact$analysis_id)
  manifest$workload$release_lattice <- list(
    version = "biomedical-capsule-common-lattice-v1",
    transform_rule = "raw_coordinate_left_shift_to_common_numeric_grid_v1",
    output_lattice_bits = as.numeric(artifact$numeric_grid_bits), output_lattice_scale = scale,
    integer_l1_sensitivity_steps = scale * (length(coordinates) + 1))
  manifest$workload$capsule_mechanism <- list(mechanism = "discrete-laplace", sensitivity_norm = "l1")
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  block <- layout$blocks[[paste0("gaussian_models::", artifact$analysis_id)]]
  normalized <- .dsvert_dp_gaussian_artifact(manifest, artifact$dataset,
    artifact$analysis_id, artifact$owner_peer, artifact$adjacency, scale, capacity)
  compiled <- list(
    artifact = list(artifact_key = strrep("c", 64L),
      semantic = list(source_claim_set_sha256 = strrep("d", 64L))),
    layout = layout, lattice = .dsvert_dp_synopsis_client_lattice(manifest, layout),
    profile = .dsvert_vector_profile("discrete-laplace",
      backend = .DSVERT_CLIENT_VECTOR_EXACT_BACKEND),
    physical = list(request = list(epsilon = "8", delta = "0"),
      full_plan_sha256 = strrep("9", 64L),
      full_plan = list(maximum_noise_magnitude = "100",
        implementation_delta_numerator = "0", implementation_delta_denominator = "1")),
    receipts = list())
  compiled$artifact$artifact_key <- digest::digest(charToRaw(paste0(
    "dsVert/analysis-artifact-key/v1|",
    .dsvert_joint_dp_client_json(.dsvert_joint_dp_client_canonical(
      compiled$artifact$semantic)))), algo = "sha256", serialize = FALSE)
  execution <- .dsvert_dp_synopsis_client_execution(compiled)
  ncoords <- length(coordinates) + 1L
  public <- list(version = .DSVERT_CLIENT_SYNOPSIS_PUBLIC_VERSION,
    artifact_key = compiled$artifact$artifact_key, execution_id = execution$execution_id,
    contract_sha256 = strrep("7", 64L), attempt_sha256 = strrep("8", 64L),
    result_set_sha256 = strrep("6", 64L), public_chunk_index = 0,
    public_chunk_count = 1, coordinate_offset = 0, coordinate_count = ncoords,
    output_lattice_bits = as.numeric(artifact$numeric_grid_bits), output_lattice_scale = scale,
    scaled_values = as.list(format(round(c(coordinates[[1L]], coordinates) * scale),
      scientific = FALSE, trim = TRUE)),
    value_encoding = "nonnegative-decimal-integer-common-lattice-v1",
    postprocessing = compiled$profile$postprocessing,
    source_values_exposed = FALSE, preclamp_values_exposed = FALSE)
  chunk_hash <- .dsvert_dp_synopsis_client_hash(
    "dsVert/stateless-catalog-synopsis/public-chunk/v1|", public)
  root <- .dsvert_vector_merkle_root(chunk_hash)
  signed <- function(value, peer, domain) c(value, list(signature = b64(
    openssl::ed25519_sign(charToRaw(paste0(domain,
      .dsvert_joint_dp_client_json(value))), keys[[peer]]))))
  release_values <- stats::setNames(lapply(seq_along(authorities), function(index) {
    peer <- authorities[[index]]
    signed(list(version = .DSVERT_CLIENT_SYNOPSIS_RELEASE_VERSION,
      phase = "synopsis_released", execution_id = execution$execution_id,
      artifact_key = public$artifact_key, contract_sha256 = public$contract_sha256,
      attempt_sha256 = public$attempt_sha256, source_contract_sha256 = strrep("e", 64L),
      result_set_sha256 = public$result_set_sha256,
      local_authority = list(peer_name = peer, identity_pk = unname(pins[[peer]]),
        role = c("primary_noise_authority", "secondary_noise_authority")[[index]]),
      public_chunk_count = 1, final_chunk_hashes = list(chunk_hash), final_vector_root = root,
      output_lattice_bits = as.numeric(artifact$numeric_grid_bits),
      output_lattice_scale = as.character(scale), mechanism = compiled$profile$release_mechanism,
      epsilon = "8", delta = "0", implementation_delta_numerator = "0",
      implementation_delta_denominator = "1", delta_aggregation = compiled$profile$delta_aggregation,
      postprocessing = compiled$profile$postprocessing, all_public_chunks_durable = TRUE,
      intermediate_payload_exposed = FALSE, durable_replay = TRUE, capability_available = TRUE),
      peer, "dsVert/stateless-catalog-synopsis/release/v1|")
  }), authorities)
  vector_release <- lapply(release_values, .dsvert_joint_dp_client_json)
  replay <- list(version = .DSVERT_CLIENT_SYNOPSIS_REPLAY_VERSION,
    phase = "synopsis_public_chunk_replayed", execution_id = execution$execution_id,
    artifact_key = public$artifact_key, contract_sha256 = public$contract_sha256,
    attempt_sha256 = public$attempt_sha256, source_contract_sha256 = strrep("e", 64L),
    result_set_sha256 = public$result_set_sha256, final_vector_root = root,
    public_chunk_index = 0, public_chunk_count = 1, chunk_sha256 = chunk_hash,
    chunk = public, merkle_proof = list(), durable_replay = TRUE,
    source_store_read = FALSE, sampler_invoked = FALSE, finalizer_invoked = FALSE, transport_read = FALSE)
  vector_replay <- list(`0` = stats::setNames(rep(list(.dsvert_joint_dp_client_json(replay)), 2L), authorities))
  evidence <- if (identical(artifact$version, .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION)) {
    stats::setNames(lapply(authorities, function(peer) signed(list(
      version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_PUBLIC_EVIDENCE_VERSION,
      phase = "cross_gaussian_public_result_evidence", analysis_id = artifact$analysis_id,
      peer_name = peer, peer_identity_pk = unname(pins[[peer]]),
      artifact_sha256 = .dsvert_dp_capsule_source_hash(normalized),
      source_contract_sha256 = strrep("e", 64L),
      private_layout_sha256 = .dsvert_dp_gaussian_cross_layout_client(manifest)$transport_coordinate_order_sha256,
      transcript_sha256 = .dsvert_dp_capsule_source_hash(normalized$transcript),
      numeric_certificate_sha256 = .dsvert_dp_capsule_source_hash(normalized$numeric_certificate),
      exact_transcript_sha256 = strrep("f", 64L), coordinate_count = block$length,
      public_start = block$start, public_end = block$end,
      public_coordinate_order_sha256 = layout$sha256, ring_bits = 128,
      frac_bits = as.numeric(artifact$numeric_grid_bits), state = "complete", fixed_transcript = TRUE,
      private_result_exposed = FALSE, exact_intermediates_exposed = FALSE, alignment_hash_exposed = FALSE),
      peer, paste0(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_SIGNATURE_DOMAIN, "cross-gaussian-synopsis-evidence|"))), authorities)
  } else NULL
  context <- list(servers = peers, designated = authorities, pinset = pins)
  policy <- list(peer_pinset = as.list(pins), designated_noise_peers = as.list(authorities),
    peer_pinset_sha256 = .dsvert_dp_pinset_hash(pins),
    unit_capacity = capacity, numeric_grid_bits = artifact$numeric_grid_bits,
    adjacency = artifact$adjacency, cohort_id = "cohort")
  status <- stats::setNames(lapply(peers, function(peer) list(policy = policy)), peers)
  trusted <- list(context = context, manifest = manifest, status = status)
  bundle <- list(manifest_sha256 = strrep("a", 64L), schema_sha256 = strrep("b", 64L),
    workload_contract_sha256 = strrep("c", 64L), schema_json = "{}")
  list(artifact = normalized, raw_artifact = artifact, coordinates = coordinates,
    compiled = compiled, compilation = list(receipts = list()), trusted = trusted,
    bundle = bundle, status = status, vector_release = vector_release, vector_replay = vector_replay,
    evidence = evidence, pins = pins, keys = keys, conns = original$conns,
    scale = scale, capacity = capacity, block = block)
}

.with_gaussian_current_fixture <- function(fixture, code) {
  testthat::local_mocked_bindings(
    .dsvert_dp_synopsis_client_bundle = function(...) fixture$trusted,
    .dsvert_dp_synopsis_client_compile = function(...) fixture$compiled,
    .dsvert_dp_gaussian_synopsis_trusted = function(...) {
      list(trusted = fixture$trusted, bundle = fixture$bundle)
    },
    .dsvert_dp_vector_accuracy_radius = function(...) list(
      radius = 0.01, confidence = 0.95, method = "fixture_public_calibration",
      implementation_tv_upper_bound = 0, sampler_tv_upper_bound = 0,
      additional_privacy_cost = c(epsilon = 0, delta = 0)),
    .package = "dsVertClient")
  force(code)
}

.gaussian_current_source <- function(fixture) {
  vector <- .dsvert_dp_synopsis_public_vector_v1(fixture$vector_release,
    fixture$vector_replay, fixture$bundle, fixture$status, fixture$compilation)
  context <- list(synopsis = TRUE, release = vector, status = fixture$status,
    manifest_bundle = fixture$bundle, verification_compilation = fixture$compilation,
    cross_gaussian_evidence = stats::setNames(list(fixture$evidence), fixture$artifact$analysis_id),
    manifest = fixture$trusted$manifest, layout = fixture$compiled$layout,
    lattice = fixture$compiled$lattice, adjacency = fixture$artifact$adjacency)
  certificate <- .dsvert_dp_gaussian_synopsis_certificate_build(
    context, fixture$artifact, fixture$block, fixture$coordinates)
  verification <- ds.validateDPGaussianCertificate(certificate, trusted_pinset = fixture$pins)
  list(context = context, metadata = .dsvert_dp_vector_public_metadata(context),
    artifact = verification$artifact, coordinates = verification$coordinates,
    moment = verification$validated_moment, certificate = certificate, verification = verification,
    scale = fixture$scale, capacity = fixture$capacity)
}

.gaussian_current_fit <- function(fixture, ridge = 0) {
  source <- .gaussian_current_source(fixture)
  testthat::with_mocked_bindings(
    ds.vertDPGaussian(fixture$artifact$dataset, fixture$artifact$analysis_id,
      ridge = ridge, datasources = fixture$conns),
    .dsvert_dp_gaussian_synopsis_release = function(...) source,
    .package = "dsVertClient")
}
