test_that("Cox public evidence authenticates cold two-authority source and release bindings", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .cox_cross_client_fixture(capacity = 4, owners = owners)
    f$artifact <- .dsvert_dp_cox_cross_workload_artifact(f$contract)
    f$manifest <- list(capsule_identity = list(capsule_id = strrep("a", 64)),
      workload = list(capsule_mechanism = list(mechanism = "discrete-laplace"),
        families = list(gaussian_models = list(artifacts = list(cox_grid = f$artifact)))))
    peers <- c("site_a", "site_b")
    keys <- f$keys; b64 <- f$b64
    context <- list(designated = peers, pinset = f$policy$peer_pinset)
    compiled <- list(artifact = list(artifact_key = strrep("c", 64),
      semantic = list(source_claim_set_sha256 = strrep("d", 64))))
    release <- list(source_contract_sha256 = strrep("e", 64), artifact_key = strrep("c", 64),
      execution_id = strrep("4", 64), final_vector_root = strrep("5", 64), result_set_sha256 = strrep("6", 64))
    capsule_id <- .dsvert_dp_synopsis_client_hash(
      "dsVert/stateless-catalog-synopsis/source-namespace/v1|",
      list(version = "dsvert-stateless-catalog-synopsis-source-contract-v1",
        manifest_capsule_id = f$manifest$capsule_identity$capsule_id,
        artifact_key = compiled$artifact$artifact_key,
        source_claim_set_sha256 = compiled$artifact$semantic$source_claim_set_sha256))
    semantic_key <- .dsvert_dp_capsule_source_hash(list(
      version = "cross-grid-semantic-release-key-v2", capsule_id = capsule_id,
      source_contract_sha256 = release$source_contract_sha256,
      signed_contract = f$artifact$signed_contract, family = "cox", version_family = "cox_grid_cross_v1",
      mechanism = f$manifest$workload$capsule_mechanism, alignment = f$contract$spec$alignment,
      caps = f$artifact$sensitivity))
    receipts <- setNames(lapply(peers, function(peer) list(
      version = "dsvert-cox-staged-receipt-v1", phase = "published", capsule_id = capsule_id,
      analysis_id = "cox_grid", peer_name = peer, peer_identity_pk = unname(context$pinset[[peer]]),
      semantic_key = semantic_key, artifact_sha256 = .dsvert_dp_capsule_source_hash(f$artifact),
      source_contract_sha256 = release$source_contract_sha256,
      profile_sha256 = f$artifact$numeric_certificate$profile_sha256,
      certificate_sha256 = f$artifact$numeric_certificate$certificate_sha256,
      stage_plan_digest = strrep("f", 64), private_result_exposed = FALSE,
      coordinate_count = f$artifact$coordinate_count, stage_receipt = strrep("1", 64),
      artifact_key = release$artifact_key, execution_id = release$execution_id,
      final_vector_root = release$final_vector_root, result_set_sha256 = release$result_set_sha256)), peers)
    wire <- function(values) setNames(lapply(names(values), function(peer) {
      value <- values[[peer]]
      value$signature <- b64(openssl::ed25519_sign(charToRaw(paste0(
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SIGNATURE_DOMAIN, "cross-grid-result|",
        .dsvert_joint_dp_client_json(value))), keys[[peer]]))
      .dsvert_joint_dp_client_json(value)
    }), names(values))
    validate <- function(values, released = release) .dsvert_dp_cox_cross_public_evidence_set(
      values, context, f$manifest, "cox_grid", released, compiled, f$policy, f$schema_manifest)
    signed <- wire(receipts)
    evidence <- validate(signed)
    expect_named(evidence, peers)
    expect_identical(lapply(evidence, .dsvert_joint_dp_client_json), signed)
    # Persisted JSON remains verifiable with no orchestration result in memory.
    expect_identical(validate(lapply(.dsvert_joint_dp_client_decode(
      .dsvert_joint_dp_client_json(evidence), "evidence", 2^20), .dsvert_joint_dp_client_json)), evidence)
    expect_error(validate(signed[1]))
    expect_error(validate(c(signed, signed[1])))
    swapped <- signed; swapped[[2]] <- signed[[1]]
    expect_error(validate(swapped))
    for (field in c("capsule_id", "source_contract_sha256", "semantic_key", "artifact_sha256",
        "artifact_key", "execution_id", "final_vector_root", "result_set_sha256")) {
      changed <- receipts
      changed <- lapply(changed, function(value) { value[[field]] <- strrep("2", 64); value })
      expect_error(validate(wire(changed)))
    }
    for (field in c("stage_receipt", "stage_plan_digest")) {
      changed <- receipts; changed[[2]][[field]] <- strrep("2", 64)
      expect_error(validate(wire(changed)))
      changed <- lapply(receipts, function(value) { value[[field]] <- strrep("0", 64); value })
      expect_error(validate(wire(changed)))
    }
    expect_error(validate(signed, list(source_contract_sha256 = strrep("2", 64))))

    # A valid signature from a source owner cannot substitute for an authority.
    if (owners > 2L) {
      impostor <- receipts[[2L]]
      impostor$peer_name <- "site_c"
      impostor$peer_identity_pk <- unname(context$pinset[["site_c"]])
      expect_error(validate(wire(list(site_a = receipts[[1L]], site_c = impostor))))
    }
    forged <- signed
    decoded <- .dsvert_joint_dp_client_decode(forged[[1L]], "receipt", 2^20)
    decoded$stage_receipt <- strrep("2", 64)
    forged[[1L]] <- .dsvert_joint_dp_client_json(decoded)
    expect_error(validate(forged))
    for (field in c("coordinate_count", "private_result_exposed", "version", "phase",
        "analysis_id", "profile_sha256", "certificate_sha256")) {
      changed <- lapply(receipts, function(value) { value[[field]] <- if (field == "coordinate_count") -1 else "wrong"; value })
      expect_error(validate(wire(changed)))
    }
    for (field in c("share", "validity_share", "cox_loss")) {
      changed <- lapply(receipts, function(value) { value[[field]] <- "private"; value })
      expect_error(validate(wire(changed)))
    }
    original <- f$manifest
    for (field in c("family", "version", "spec_version", "signed_contract", "statistic_maximum")) {
      f$manifest$workload$families$gaussian_models$artifacts$cox_grid[[field]] <- "wrong"
      expect_error(validate(signed))
      f$manifest <- original
    }
    original_context <- context
    context$designated <- rev(peers)
    expect_named(validate(signed), rev(peers))
    context$designated <- rep(peers[[1L]], 2L)
    expect_error(validate(signed))
    context <- original_context
    context$pinset[[1L]] <- strrep("A", 43)
    expect_error(validate(signed))
    context <- original_context
    # Compose real RELEASE signatures and hashed bilateral REPLAY with the Cox
    # publication signatures. Compilation here is an authenticated-input fixture;
    # this does not claim that the public Cox compiler/runner is enabled.
    manifest <- f$manifest
    manifest$workload$coordinate_count <- 4L
    manifest$workload$families <- c(manifest$workload$families, list(
      admitted_count = list(owner_peer = "site_a", dataset = "aligned", statistic_maximum = 4),
      numeric_moments = list(artifacts = list()), numeric_pair_moments = list(artifacts = list()),
      fixed_numeric_histograms = list(artifacts = list()), categorical_marginals = list(artifacts = list()),
      categorical_pairs = list(sets = list()), correlation_artifacts = list(),
      describe_artifacts = list(), survival_artifacts = list()))
    manifest$workload$release_lattice <- list(version = "biomedical-capsule-common-lattice-v1",
      transform_rule = "raw_coordinate_left_shift_to_common_numeric_grid_v1",
      output_lattice_bits = 8, output_lattice_scale = 256,
      integer_l1_sensitivity_steps = 256 + f$contract$spec$sensitivity$raw_l1_sensitivity)
    compiled$layout <- .dsvert_dp_capsule_vector_layout(manifest)
    compiled$lattice <- .dsvert_dp_synopsis_client_lattice(manifest, compiled$layout)
    compiled$profile <- list(exact_gc = TRUE, gaussian = FALSE,
      release_mechanism = "discrete-laplace", delta_aggregation = "fixture",
      postprocessing = "fixture-clamp")
    compiled$physical <- list(request = list(epsilon = "8", delta = "0"),
      full_plan = list(maximum_noise_magnitude = "100", implementation_delta_numerator = "0",
        implementation_delta_denominator = "1"))
    execution <- .dsvert_dp_synopsis_client_execution(compiled)
    public <- list(version = .DSVERT_CLIENT_SYNOPSIS_PUBLIC_VERSION,
      artifact_key = release$artifact_key, execution_id = execution$execution_id,
      contract_sha256 = strrep("7", 64), attempt_sha256 = strrep("8", 64),
      result_set_sha256 = release$result_set_sha256, public_chunk_index = 0,
      public_chunk_count = 1, coordinate_offset = 0, coordinate_count = 4,
      output_lattice_bits = 8, output_lattice_scale = 256,
      scaled_values = list("1024", "117", "83", "83"),
      value_encoding = "nonnegative-decimal-integer-common-lattice-v1",
      postprocessing = compiled$profile$postprocessing,
      source_values_exposed = FALSE, preclamp_values_exposed = FALSE)
    chunk_hash <- .dsvert_dp_synopsis_client_hash(
      "dsVert/stateless-catalog-synopsis/public-chunk/v1|", public)
    root <- .dsvert_vector_merkle_root(chunk_hash)
    release_values <- setNames(lapply(seq_along(peers), function(i) {
      peer <- peers[[i]]
      list(version = .DSVERT_CLIENT_SYNOPSIS_RELEASE_VERSION, phase = "synopsis_released",
        execution_id = execution$execution_id, artifact_key = release$artifact_key,
        contract_sha256 = public$contract_sha256, attempt_sha256 = public$attempt_sha256,
        source_contract_sha256 = release$source_contract_sha256, result_set_sha256 = release$result_set_sha256,
        local_authority = list(peer_name = peer, identity_pk = unname(context$pinset[[peer]]),
          role = c("primary_noise_authority", "secondary_noise_authority")[[i]]),
        public_chunk_count = 1, final_chunk_hashes = list(chunk_hash), final_vector_root = root,
        output_lattice_bits = 8, output_lattice_scale = "256", mechanism = "discrete-laplace",
        epsilon = "8", delta = "0", implementation_delta_numerator = "0",
        implementation_delta_denominator = "1", delta_aggregation = compiled$profile$delta_aggregation,
        postprocessing = compiled$profile$postprocessing, all_public_chunks_durable = TRUE,
        intermediate_payload_exposed = FALSE, durable_replay = TRUE, capability_available = TRUE)
    }), peers)
    release_wire <- function(values) setNames(lapply(peers, function(peer) {
      value <- values[[peer]]
      value$signature <- b64(openssl::ed25519_sign(charToRaw(paste0(
        "dsVert/stateless-catalog-synopsis/release/v1|", .dsvert_joint_dp_client_json(value))), keys[[peer]]))
      .dsvert_joint_dp_client_json(value)
    }), peers)
    vector_release <- release_wire(release_values)
    replay <- list(version = .DSVERT_CLIENT_SYNOPSIS_REPLAY_VERSION,
      phase = "synopsis_public_chunk_replayed", execution_id = execution$execution_id,
      artifact_key = release$artifact_key, contract_sha256 = public$contract_sha256,
      attempt_sha256 = public$attempt_sha256, source_contract_sha256 = release$source_contract_sha256,
      result_set_sha256 = release$result_set_sha256, final_vector_root = root,
      public_chunk_index = 0, public_chunk_count = 1, chunk_sha256 = chunk_hash,
      chunk = public, merkle_proof = list(), durable_replay = TRUE,
      source_store_read = FALSE, sampler_invoked = FALSE, finalizer_invoked = FALSE, transport_read = FALSE)
    replay_wire <- function(value) list(`0` = setNames(rep(list(.dsvert_joint_dp_client_json(value)), 2), peers))
    vector_replay <- replay_wire(replay)
    published <- wire(lapply(receipts, function(value) {
      value$execution_id <- execution$execution_id; value$final_vector_root <- root; value
    }))
    trusted <- list(context = context, manifest = manifest)
    read <- function(r = vector_release, p = vector_replay, e = published, c = compiled) {
      .dsvert_dp_cox_cross_read_vector(r, p, e, trusted, c, f$policy, f$schema_manifest, "cox_grid")
    }
    result <- read()
    expect_identical(result$coordinates, c(117, 83, 83))
    expect_identical(.dsvert_dp_cox_grid_cross_moment(result$coordinates, result$contract$spec)$selected_candidate, 2L)
    expect_identical(read(p = replay_wire(.dsvert_joint_dp_client_decode(
      .dsvert_joint_dp_client_json(replay), "cold replay", 2^20))), result)
    expect_error(read(r = vector_release[1]))
    expect_error(read(p = list()))
    expect_error(read(e = published[1]))
    expect_error(read(e = signed)) # Authentic publication from a different vector.
    changed <- replay; changed$chunk$scaled_values[[2]] <- "118"
    expect_error(read(p = replay_wire(changed)))
    changed <- vector_replay; changed[[1]][[2]] <- .dsvert_joint_dp_client_json(replay$chunk)
    expect_error(read(p = changed))
    changed <- release_values; changed[[2]]$epsilon <- "9"
    expect_error(read(r = release_wire(changed)))
    changed <- vector_release
    decoded <- .dsvert_joint_dp_client_decode(changed[[1]], "release", 2^20)
    decoded$final_vector_root <- strrep("2", 64)
    changed[[1]] <- .dsvert_joint_dp_client_json(decoded)
    expect_error(read(r = changed))
    changed <- compiled; changed$layout$blocks[[2]]$start <- 1L
    expect_error(read(c = changed))
    changed <- compiled; changed$lattice$scale_shifts[[2]] <- 8L
    expect_error(read(c = changed))
    # Even an otherwise correctly signed vector must obey the signed Cox caps.
    reseal <- function(values) {
      changed <- replay
      changed$chunk$scaled_values <- as.list(values)
      changed$chunk_sha256 <- .dsvert_dp_synopsis_client_hash(
        "dsVert/stateless-catalog-synopsis/public-chunk/v1|", changed$chunk)
      changed$final_vector_root <- .dsvert_vector_merkle_root(changed$chunk_sha256)
      r <- lapply(release_values, function(value) {
        value$final_chunk_hashes <- list(changed$chunk_sha256)
        value$final_vector_root <- changed$final_vector_root; value
      })
      e <- wire(lapply(receipts, function(value) {
        value$execution_id <- execution$execution_id
        value$final_vector_root <- changed$final_vector_root; value
      }))
      read(r = release_wire(r), p = replay_wire(changed), e = e)
    }
    upper <- unlist(f$contract$spec$sensitivity$maximum_coordinates, use.names = FALSE)
    expect_identical(reseal(c("1024", as.character(upper)))$coordinates, as.numeric(upper))
    expect_error(reseal(c("1024", as.character(upper + 1))))
    expect_error(reseal(c("1024", "-1", "83", "83")))
    expect_error(reseal(c("1024", "1.5", "83", "83")))
    formula <- stats::reformulate(unlist(f$contract$spec$predictor_order),
      response = "Surv(site_a$time, site_a$event)")
    point <- .dsvert_dp_cox_grid_cross_impl(formula, "aligned", "cox_grid",
      setNames(rep(list(NULL), owners), names(f$policy$peer_pinset)),
      .release = function(...) read())
    expect_identical(point$selected_candidate, 2L)
    expect_identical(point$selected_dp_partial_loss, 83 / 256)
    expect_false(point$source_values_exposed)
    # Discovery recognizes Cox; contract-only state still prevents public release.
    expect_length(.dsvert_dp_glm_grid_cross_artifacts(f$manifest), 1L)
    expect_false(.dsvert_dp_synopsis_supported_glm_grid_cross_v1(f$manifest))
    expect_error(.dsvert_dp_cox_grid_cross_release("aligned", "cox_grid", list()))
  }
})

test_that("Cox public evidence retains the private-v2 400-row scope", {
  f <- .cox_cross_client_fixture(capacity = 401)
  artifact <- .dsvert_dp_cox_cross_workload_artifact(f$contract)
  manifest <- list(workload = list(families = list(
    gaussian_models = list(artifacts = list(cox_grid = artifact)))))
  expect_error(.dsvert_dp_cox_cross_public_evidence_set(list(), list(), manifest,
    "cox_grid", list(), list(), f$policy, f$schema_manifest),
    class = "dsvert_dp_public_failure")
})
