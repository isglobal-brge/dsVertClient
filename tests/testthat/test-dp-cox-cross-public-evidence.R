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
    # No public source discovery or release route is opened by this verifier.
    expect_length(.dsvert_dp_glm_grid_cross_artifacts(f$manifest), 0L)
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
