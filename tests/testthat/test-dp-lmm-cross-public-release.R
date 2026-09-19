.lmm_release_client_fixture <- function(owners = 2L) {
  f <- .grouped_cross_client_fixture("lmm", owners = owners)
  f$policy$unit_capacity <- 15
  f$raw$parameters <- list(objective = "ml", variance_grid = list(
    list(residual_variance = .25, random_intercept_variance = 0),
    list(residual_variance = 1, random_intercept_variance = .25)))
  spec <- .dsvert_dp_grouped_grid_cross_spec(f$raw, f$policy, f$schema)
  artifact <- .dsvert_dp_grouped_grid_cross_artifact(spec)
  f$contract <- f$sign(list(version = f$contract$version, spec = spec, artifact = artifact,
    source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, artifact)))
  f$artifact <- .dsvert_dp_grouped_cross_workload_artifact(f$contract)
  f$manifest <- list(admission = list(unit_capacity = 15, adjacency = f$policy$adjacency),
    bounds = list(numeric_grid_bits = 8),
    capsule_identity = list(capsule_id = strrep("a", 64)),
    workload = list(coordinate_count = artifact$coordinate_count + 1,
      capsule_mechanism = list(mechanism = "discrete-laplace", source_context_hash = strrep("b", 64)),
      families = list(admitted_count = list(owner_peer = "site_a", dataset = "aligned"),
        numeric_moments = list(artifacts = list()), numeric_pair_moments = list(artifacts = list()),
        gaussian_models = list(artifacts = list(grouped = f$artifact)),
        fixed_numeric_histograms = list(artifacts = list()), categorical_marginals = list(artifacts = list()),
        categorical_pairs = list(sets = list()), correlation_artifacts = list(),
        describe_artifacts = list(), survival_artifacts = list())))
  f
}

test_that("LMM ML discovery and private source projection preserve signed padded owners", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .lmm_release_client_fixture(owners)
    artifact <- f$artifact
    expect_type(artifact$participating_peers, "list")
    expect_type(artifact$computation_peers, "list")
    expect_identical(.dsvert_joint_dp_client_json(artifact),
      .dsvert_joint_dp_client_json(jsonlite::fromJSON(.dsvert_joint_dp_client_json(artifact),
        simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)))
    expect_identical(artifact$signed_contract, .dsvert_joint_dp_client_json(f$contract))
    expect_named(.dsvert_dp_glm_grid_cross_artifacts(f$manifest), "grouped")
    wire_schema <- .dsvert_joint_dp_client_decode(.dsvert_joint_dp_client_json(f$schema_manifest),
      "signed grouping schema", 2^20)
    expect_equal(.dsvert_dp_glm_grid_cross_schema_validate(f$policy,
      wire_schema$logical_snapshot, wire_schema, .dsvert_dp_glm_grid_cross_verify), f$schema)
    expect_silent(.dsvert_dp_glm_grid_cross_preflight(f$manifest,
      list(pinset = f$policy$peer_pinset, designated = f$policy$designated_noise_peers),
      .dsvert_joint_dp_client_json(f$schema_manifest)))
    layout <- .dsvert_dp_gaussian_cross_layout_client(f$manifest)
    blocks <- layout$blocks
    expect_equal(length(blocks), 2 * (length(artifact$predictors) + 2))
    expect_true(all(vapply(blocks, function(block) block$length == 16, logical(1))))
    expect_equal(unlist(layout$source_peers), names(f$policy$peer_pinset))
    expect_equal(unlist(layout$computation_peers), c("site_a", "site_b"))
    for (reference in c(unlist(f$contract$spec$predictor_order), "site_b$y")) {
      value <- blocks[[paste("grouped", reference, "value", sep = "::")]]
      validity <- blocks[[paste("grouped", reference, "validity", sep = "::")]]
      expect_equal(value$fraction_bits, 50)
      expect_equal(value$maximum, 2^50)
      expect_equal(validity$fraction_bits, 0)
      expect_equal(validity$maximum, 1)
      expect_equal(validity$start, value$end + 1)
    }
    labels <- blocks[["grouped::site_a$cluster::value"]]
    expect_true(labels$private_routing_input)
    expect_equal(labels$levels, c("a", "b", "c", "d"))
    expect_identical(labels$owner_peer, "site_a")
    expect_equal(labels$fraction_bits, 0)
    expect_equal(labels$maximum, 3)
    expect_equal(layout$transport_coordinate_count,
      layout$private_start - 1 + 16 * length(blocks))
    source <- .dsvert_dp_synopsis_source_manifest_v1(f$manifest,
      list(status = structure(list(), class = "ds.vertDPSynopsisStatus")))
    expect_identical(source$private_layout_sha256, layout$transport_coordinate_order_sha256)
    expect_identical(source$purpose, .DSVERT_CLIENT_DP_GLM_GRID_CROSS_SOURCE_PURPOSE)
    expect_equal(source$coordinate_count, layout$transport_coordinate_count)
    for (family in setdiff(.DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_FAMILIES, "lmm")) {
      unavailable <- f$manifest
      unavailable$workload$families$gaussian_models$artifacts$grouped$version <-
        paste0("bounded-", gsub("_", "-", family), "-cross-grid-v1")
      expect_length(.dsvert_dp_glm_grid_cross_artifacts(unavailable), 0)
    }
  }
  for (levels in list(list(), list(a = "a"), list("a", 1), list("a", NA_character_))) {
    expect_error(.dsvert_dp_glm_grid_cross_schema_levels(levels))
  }
})

test_that("every LMM candidate enters the release lattice without a second shift", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .lmm_release_client_fixture(owners)
    manifest <- f$manifest
    bits <- f$contract$spec$numeric_grid_bits
    manifest$workload$families$admitted_count$statistic_maximum <-
      manifest$admission$unit_capacity
    manifest$workload$release_lattice <- list(
      version = "biomedical-capsule-common-lattice-v1",
      transform_rule = "raw_coordinate_left_shift_to_common_numeric_grid_v1",
      output_lattice_bits = bits, output_lattice_scale = 2^bits,
      integer_l1_sensitivity_steps = 2^bits + f$contract$spec$sensitivity$raw_l1_sensitivity)
    layout <- .dsvert_dp_capsule_vector_layout(manifest)
    lattice <- .dsvert_dp_synopsis_client_lattice(manifest, layout)
    block <- layout$blocks[["gaussian_models::grouped"]]
    positions <- seq.int(block$start, block$end)
    expect_equal(lattice$scale_shifts, c(bits, rep(0L, length(positions))))
    maxima <- unlist(f$contract$spec$sensitivity$maximum_coordinates, use.names = FALSE)
    expect_equal(as.numeric(lattice$raw_upper_bounds[positions]), maxima)
    expect_equal(as.numeric(lattice$raw_upper_bounds[positions]) *
      2^lattice$scale_shifts[positions] / lattice$output_lattice_scale, maxima / 2^bits)
    detached <- f$artifact
    detached$source_coordinate_scaling <- NULL
    expect_error(.dsvert_dp_grouped_cross_client_artifact(detached, "aligned", "grouped",
      NULL, "add_remove_patient", 2^bits, 15), class = "dsvert_dp_public_failure")
  }
})

test_that("LMM public artifacts rebuild ML bounds certificate and variance-major ordering", {
  f <- .lmm_release_client_fixture()
  validate <- function(artifact) .dsvert_dp_grouped_cross_client_artifact(
    artifact, "aligned", "grouped", NULL, "add_remove_patient", 256, 15)
  admitted <- validate(f$artifact)
  expect_equal(admitted$statistic_maximum,
    unlist(f$contract$spec$sensitivity$maximum_coordinates, use.names = FALSE))
  expect_identical(.dsvert_dp_gaussian_artifact(f$manifest, "aligned", "grouped", NULL,
    "add_remove_patient", 256, 15), admitted)
  for (mutate in list(
    function(x) { x$spec$candidate_grid <- rev(x$spec$candidate_grid); x },
    function(x) { x$spec$candidate_order <- rev(x$spec$candidate_order); x },
    function(x) { x$spec$numeric_contract$certificate_sha256 <- strrep("0", 64); x },
    function(x) { x$spec$sensitivity$maximum_coordinates[[1]] <- 1; x },
    function(x) { x$source_contract$private_layout$blocks[[1]]$length <- 15; x })) {
    expect_error(validate(.dsvert_dp_grouped_cross_workload_artifact(f$sign(mutate(f$contract)))),
      class = "dsvert_dp_public_failure")
  }
  bad <- f$artifact
  bad$statistic_maximum[[1]] <- 1
  expect_error(validate(bad), class = "dsvert_dp_public_failure")
  bad <- f$manifest
  bad$workload$coordinate_count <- 52
  expect_error(.dsvert_dp_glm_grid_cross_preflight(bad,
    list(pinset = f$policy$peer_pinset, designated = f$policy$designated_noise_peers),
    .dsvert_joint_dp_client_json(f$schema_manifest)), class = "dsvert_dp_public_failure")
  f$raw$beta_grid <- f$raw$beta_grid[1L]
  f$raw$parameters$variance_grid <- f$raw$parameters$variance_grid[1L]
  spec <- .dsvert_dp_grouped_grid_cross_spec(f$raw, f$policy, f$schema)
  artifact <- .dsvert_dp_grouped_grid_cross_artifact(spec)
  contract <- f$sign(list(version = f$contract$version, spec = spec, artifact = artifact,
    source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, artifact)))
  projected <- .dsvert_dp_grouped_cross_workload_artifact(contract)
  expect_identical(.dsvert_joint_dp_client_json(projected),
    .dsvert_joint_dp_client_json(jsonlite::fromJSON(.dsvert_joint_dp_client_json(projected),
      simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)))
  expect_silent(validate(projected))
})

test_that("cold LMM public evidence requires both real signatures and source semantic binding", {
  f <- .lmm_release_client_fixture()
  peers <- c("site_a", "site_b")
  keys <- setNames(lapply(peers, function(peer) openssl::ed25519_keygen()), peers)
  b64 <- function(value) sub("=+$", "", chartr("+/", "-_",
    gsub("[\r\n]", "", jsonlite::base64_enc(value))))
  context <- list(designated = peers, pinset = vapply(keys,
    function(key) b64(tail(as.raw(as.list(key)$pubkey), 32L)), character(1)))
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
    signed_contract = f$artifact$signed_contract, family = "lmm", version_family = "lmm_grid_cross_v1",
    mechanism = f$manifest$workload$capsule_mechanism, alignment = f$contract$spec$alignment,
    caps = f$artifact$sensitivity))
  receipts <- setNames(lapply(peers, function(peer) list(
    version = "dsvert-lmm-staged-receipt-v1", phase = "published", capsule_id = capsule_id,
    analysis_id = "grouped", peer_name = peer, peer_identity_pk = unname(context$pinset[[peer]]),
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
  validate <- function(values, released = release) .dsvert_dp_lmm_cross_public_evidence_set(
    values, context, f$manifest, "grouped", released, compiled)
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

  # Exercise the cold orchestration branch with already authenticated public
  # publication boundaries. Any source claim or sampler invocation must fail.
  context$servers <- peers
  context$all_conns <- context$conns <- setNames(rep(list(list()), 2L), peers)
  compilation <- list(version = "test-retained-public-compilation")
  published <- list(release = release, verification_compilation = compilation)
  bootstrap <- list(status = list(), manifest_bundle = list(
    manifest_sha256 = strrep("7", 64), schema_json = .dsvert_joint_dp_client_json(f$schema_manifest)))
  requests <- 0L
  local_mocked_bindings(
    .dsvert_dp_datasources = function(value) value,
    .dsvert_dp_synopsis_bootstrap_build_v1 = function(...) bootstrap,
    .dsvert_dp_synopsis_client_bundle = function(...) list(manifest = f$manifest, context = context),
    .dsvert_dp_glm_grid_cross_preflight = function(...) invisible(TRUE),
    .dsvert_dp_synopsis_publication_resume_v1 = function(...) published,
    .dsvert_vector_profile = function(...) list(),
    .dsvert_dp_synopsis_client_compile = function(value, ...) {
      expect_identical(value, compilation)
      compiled
    },
    .dsvert_dp_synopsis_runner_compile = function(...) stop("cold source resolver forbidden"),
    .dsvert_exact_gc_run = function(...) stop("cold sampler forbidden"),
    .dsvert_fanout_by_site = function(conns, calls, operation, ...) {
      expect_identical(operation, "LMM public evidence")
      expect_named(calls, peers)
      for (call in calls) {
        args <- as.list(call)
        expect_identical(args$action, "evidence")
        expect_identical(args$claim_set_json, "{}")
        expect_identical(args$compilation_json, .dsvert_joint_dp_client_json(compilation))
      }
      signed
    }, .package = "dsVertClient")
  cold <- .dsvert_dp_synopsis_vector_run(context$all_conns, .request_check = function(manifest) {
    expect_identical(manifest, f$manifest)
    requests <<- requests + 1L
  })
  expect_identical(cold$release, release)
  expect_identical(cold$cross_lmm_evidence$grouped, evidence)
  expect_equal(requests, 1L)
})
