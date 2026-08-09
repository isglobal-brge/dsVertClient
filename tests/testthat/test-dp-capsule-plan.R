.dp_capsule_plan_b64url <- function(value) {
  sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))),
    perl = TRUE)
}

.dp_capsule_plan_manifest <- function(peer) {
  scope_selection <- .dsvert_joint_dp_client_canonical(list(
    mode = "catalog_v1",
    explicit_catalog = list(
      numeric_moments = list("x"),
      categorical_marginals = list("group"),
      categorical_pairs = list(c("group", "status")),
      correlations = list(c("x", "z"))),
    referenced_by_signed_specs = list(
      numeric = list(), categorical = list(), describe = list(),
      survival = list(), gaussian = NULL, vertical_cross = NULL),
    included = list(
      numeric_moments = list("x", "z"),
      categorical_marginals = list("group", "status"),
      same_owner_categorical_pairs = list(c("group", "status")),
      same_owner_correlations = list(c("x", "z")))))
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
      schema_numeric_column_count = 2,
      schema_categorical_column_count = 2,
      possible_same_owner_numeric_pair_count = 1,
      possible_same_owner_categorical_pair_count = 1,
      included_numeric_moment_count = 2,
      included_categorical_marginal_count = 2,
      included_numeric_pair_count = 1,
      included_categorical_pair_count = 1,
      included_cross_categorical_pair_count = 0,
      numeric_moment_coordinate_count = 6,
      numeric_pair_coordinate_count = 6,
      categorical_marginal_coordinate_count = 4,
      categorical_pair_coordinate_count = 4,
      gaussian_model_coordinate_count = 0,
      projected_coordinate_count = 26,
      projected_integer_l1_sensitivity = 1792,
      projected_integer_l2_sensitivity = 1280,
      automatic_pair_expansion = "none",
      scaling_contract = paste0(
        "linear_in_declared_univariates_plus_explicit_pairs_and_",
        "declared_model_cross_products")))
  families <- list(
    admitted_count = list(owner_peer = peer, dataset = "cohort"),
    numeric_moments = list(artifacts = list(
      x = list(owner_peer = peer, dataset = "cohort", column = "x"),
      z = list(owner_peer = peer, dataset = "cohort", column = "z"))),
    numeric_pair_moments = list(artifacts = list(x_z = list(
      owner_peer = peer, dataset = "cohort", left = "x", right = "z"))),
    gaussian_models = list(artifacts = list()),
    fixed_numeric_histograms = list(artifacts = list(x_hist = list(
      owner_peer = peer, dataset = "cohort", column = "x",
      coordinate_count = 3L))),
    categorical_marginals = list(artifacts = list(
      group = list(
        owner_peer = peer, dataset = "cohort", column = "group",
        levels = c("control", "case")),
      status = list(
        owner_peer = peer, dataset = "cohort", column = "status",
        levels = c("no", "yes")))),
    categorical_pairs = list(
      sets = list(primary = list(
        owner_peer = peer, dataset = "cohort",
        columns = list(
          list(column = "group", levels = c("control", "case")),
          list(column = "status", levels = c("no", "yes"))),
        included_pairs = list(c("group", "status")),
        repeated_record_policy = "consistent_cell_else_exclude_v1",
        missingness_policy = "missing_or_out_of_domain_excluded_v1",
        coordinate_count = 4L, pair_count = 1L)),
      cross_artifacts = list()),
    correlation_artifacts = list(x_z = list()),
    describe_artifacts = list(primary = list()),
    survival_artifacts = list(primary = list(
      owner_peer = peer, dataset = "cohort", coordinate_count = 2L)))
  mechanism_selection <- list(
    version = "biomedical-capsule-noise-selection-v4",
    selector = "formal_fixed_work_backend_minimum_simultaneous_95_radius_v2",
    objective = "simultaneous_95_abs",
    tie_break = "laplace_unless_fixed_work_gaussian_strictly_improves",
    coordinate_count = 26L, epsilon = 1, allocated_delta = 2^-100,
    gaussian_eligible = TRUE, positive_delta_reserved = TRUE,
    deployed_backends = "joint-discrete-laplace-v3",
    gaussian_backend_available = FALSE,
    gaussian_unavailable_reason =
      "formal_fixed_work_gaussian_plan_unavailable",
    gaussian_calibration_request = list(version = "test-request-v1"),
    gaussian_calibration_rounding = list(
      declared_policy_values = "binary64_policy_values",
      epsilon = "inward_128_machine_epsilon_relative_guard_v1",
      delta = "inward_128_machine_epsilon_relative_guard_v1",
      l2_sensitivity = "outward_128_machine_epsilon_relative_guard_v1"),
    gaussian_plan = NULL, gaussian_plan_sha256 = NULL,
    laplace_simultaneous_95_abs = 12,
    gaussian_simultaneous_95_abs = NULL,
    deployment_rule = "formal_backend_or_explicit_laplace_fallback",
    winner = "laplace", utility_winner = "laplace",
    decision = "fixed_work_gaussian_unavailable_explicit_laplace_fallback",
    canonical_selector_invoked = TRUE,
    selector_certificate_sha256 = strrep("7", 64L))
  workload <- list(
    workload_version = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION,
    schema_attestation = list(), coordinate_count = 26,
    sensitivity = list(l1 = 1792, l2 = 1280), families = families,
    vertical_crosses = list(), primitive_scope = primitive_scope,
    release_lattice = list(
      version = "biomedical-capsule-common-lattice-v1",
      transform_rule =
        "raw_coordinate_left_shift_to_common_numeric_grid_v1",
      output_lattice_bits = 8L, output_lattice_scale = 256,
      natural_l1_sensitivity = 7, integer_l1_sensitivity_steps = 1792,
      natural_l2_sensitivity = 5, integer_l2_sensitivity_steps = 1280),
    mechanism_selection = mechanism_selection,
    registered_release_lifecycle =
      .dsvert_dp_capsule_client_registered_release_lifecycle(),
    declared_workload_fully_materialized = TRUE,
    package_family_coverage_complete = FALSE,
    execution_state = .DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE,
    capsule_mechanism = list(
      release_scope = "joint_mpc_single_opening",
      capability_id = "joint_mpc_single_opening_v1",
      producer = "biomedical.capsule.vector.v2",
      purpose = "biomedical_capsule_full_vector",
      source_context_hash = strrep("3", 64L),
      mechanism = "discrete-laplace",
      mechanism_version = "biomedical-capsule-vector-v2",
      sampler = "two_private_hmac_seeds_one_gc_sample_v1",
      sensitivity_norm = "l1", sensitivity = 1792,
      coordinate_count = 26, uses_delta = TRUE,
      clipping_hash = strrep("4", 64L), ring_bits = 128L,
      frac_bits = 0L))
  list(
    version = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION,
    logical_snapshot = list(
      logical_snapshot_id = "plan-snapshot", version = "v1",
      alignment_protocol_version = 1),
    capsule_schema = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION,
    admission = list(protocol = "test-admission-v1"),
    bounds = list(protocol = "test-bounds-v1"), workload = workload,
    capsule_identity = list(capsule_id = strrep("5", 64L), contract = list()),
    execution_state = .DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE)
}

.dp_capsule_plan_fixture <- function(k) {
  stopifnot(k %in% c(2L, 3L, 5L))
  peers <- paste0("peer_", letters[seq_len(k)])
  keys <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  pins <- vapply(keys, function(key) {
    .dp_capsule_plan_b64url(tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  pin_hash <- .dsvert_dp_pinset_hash(pins)
  statuses <- stats::setNames(lapply(seq_along(peers), function(index) {
    peer <- peers[[index]]
    list(
      policy = list(
        domain = "plan-study", cohort_id = "plan-cohort",
        peer_pinset_sha256 = pin_hash, capsule_epsilon = 1,
        capsule_delta = 2^-100),
      noise_root = list(
        privacy_epoch = as.numeric(index), key_id = paste0("noise-", peer)))
  }), peers)
  context <- list(
    servers = peers, pinset = pins, status = statuses,
    designated = peers[1:2],
    all_conns = stats::setNames(lapply(peers, function(peer) {
      structure(list(), class = "plan-test-connection")
    }), peers))
  sign <- function(message, peer) {
    .dp_capsule_plan_b64url(openssl::ed25519_sign(message, keys[[peer]]))
  }
  sign_object <- function(unsigned, peer, domain) {
    c(unsigned, list(signature = sign(
      .dsvert_dp_capsule_manifest_message(domain, unsigned), peer)))
  }
  drafts <- stats::setNames(lapply(seq_along(peers), function(index) {
    peer <- peers[[index]]
    column <- paste0("x_", peer)
    columns <- stats::setNames(list(list(
      kind = "numeric", owner_peer = peer, lower = 0, upper = 1)), column)
    unsigned <- list(
      version = .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_DRAFT_VERSION,
      phase = "custodian_policy_draft", peer_name = peer,
      peer_identity_pk = unname(pins[[peer]]),
      peer_pinset_sha256 = pin_hash, domain = "plan-study",
      cohort_id = "plan-cohort",
      dataset_mapping_mode = "automatic_single_local_dataset",
      datasets = stats::setNames(list(list(
        dataset_id = paste0("cohort-", peer), dataset_version = "v1",
        schema_version = .DSVERT_CLIENT_DP_CAPSULE_POLICY_SCHEMA_VERSION,
        alignment_group = "plan-cohort", alignment_protocol_version = 1,
        patient_column = "patient_id", columns = columns)),
        paste0("data_", peer)),
      workload_fragments = list(
        describe = list(), survival = list(), gaussian = list(),
        vertical_cross = list()),
      data_access = FALSE, patient_derived_metadata = FALSE,
      operation_limit = FALSE, request_limit = FALSE,
      history_can_deny_operation = FALSE)
    .dsvert_joint_dp_client_json(sign_object(unsigned, peer, "draft"))
  }), peers)
  state <- new.env(parent = emptyenv())
  state$calls <- list()
  state$different_build <- FALSE
  state$invalid_build_signature <- FALSE
  state$tamper_scope <- FALSE
  state$tamper_layout <- FALSE
  state$tamper_cost <- FALSE
  state$tamper_calibration <- FALSE

  artifact_index <- function(manifest, policy, manifest_sha256) {
    value <- .dsvert_joint_dp_client_canonical(list(
      version = .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_VERSION,
      context = list(privacy_epoch_scope =
                       "per_peer_signed_receipts_v1"),
      gaussian_models = list()))
    list(
      value = value, count = 0,
      root = digest::digest(
        .dsvert_joint_dp_client_json(value), "sha256", serialize = FALSE))
  }
  source_manifest <- function(manifest_json, context) {
    value <- .dp_capsule_plan_manifest(peers[[1L]])
    build_call <- state$calls[[length(state$calls)]][[1L]]
    schema <- .dsvert_joint_dp_client_decode(
      as.list(build_call)$schema_json, "test signed schema", 8L * 1024L^2)
    unsigned <- schema[setdiff(names(schema), "signatures")]
    manifest_value <- .dsvert_joint_dp_client_decode(
      manifest_json, "test manifest", 8L * 1024L^2)
    value$logical_snapshot <- unsigned$logical_snapshot
    value$workload$schema_attestation <- list(
      manifest_sha256 = manifest_value$schema_sha256,
      signers = as.list(sort(names(schema$signatures), method = "radix")),
      signatures = schema$signatures)
    if (isTRUE(state$tamper_scope)) {
      value$workload$primitive_scope$selection_sha256 <- strrep("0", 64L)
    }
    if (isTRUE(state$tamper_layout)) {
      value$workload$families$fixed_numeric_histograms$artifacts$x_hist$
        coordinate_count <- 2L
    }
    if (isTRUE(state$tamper_cost)) {
      value$workload$primitive_scope$projected_cost$
        included_cross_categorical_pair_count <- 1
    }
    if (isTRUE(state$tamper_calibration)) {
      value$workload$mechanism_selection$winner <- "gaussian"
      value$workload$mechanism_selection$utility_winner <- "gaussian"
    }
    list(
      value = value,
      capsule_id = digest::digest(
        manifest_json, "sha256", serialize = FALSE),
      release_coordinate_count = value$workload$coordinate_count)
  }
  fanout <- function(conns, calls, operation, .aggregate) {
    state$calls[[length(state$calls) + 1L]] <- calls
    endpoint <- as.character(as.list(calls[[1L]])[[1L]])
    if (identical(endpoint, "dsvertDPCapsuleManifestDraftDS")) {
      return(drafts[names(calls)])
    }
    if (identical(endpoint, "dsvertDPCapsuleManifestSignDS")) {
      return(stats::setNames(lapply(names(calls), function(peer) {
        arguments <- as.list(calls[[peer]])
        schema <- .dsvert_joint_dp_client_decode(
          arguments$schema_json, "test schema", 8L * 1024L^2)
        schema_signature <- sign(charToRaw(paste0(
          .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
          .dsvert_joint_dp_client_json(schema))), peer)
        .dsvert_joint_dp_client_json(list(
          version = .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_SIGN_VERSION,
          phase = "global_schema_policy_verified", peer_name = peer,
          peer_identity_pk = unname(pins[[peer]]),
          peer_pinset_sha256 = pin_hash,
          schema_sha256 = digest::digest(
            arguments$schema_json, "sha256", serialize = FALSE),
          workload_contract_sha256 = digest::digest(
            arguments$workload_contract_json, "sha256", serialize = FALSE),
          logical_snapshot = schema$logical_snapshot,
          schema_signature = schema_signature, data_access = FALSE,
          operation_limit = FALSE, request_limit = FALSE,
          history_can_deny_operation = FALSE))
      }), names(calls)))
    }
    if (!identical(endpoint, "dsvertDPCapsuleManifestBuildDS")) {
      stop("protected or release endpoint called")
    }
    stats::setNames(lapply(names(calls), function(peer) {
      arguments <- as.list(calls[[peer]])
      signed_schema <- .dsvert_joint_dp_client_decode(
        arguments$schema_json, "test signed schema", 8L * 1024L^2)
      unsigned_schema <- signed_schema[
        setdiff(names(signed_schema), "signatures")]
      schema_json <- .dsvert_joint_dp_client_json(unsigned_schema)
      schema_hash <- digest::digest(
        schema_json, "sha256", serialize = FALSE)
      workload_hash <- digest::digest(
        arguments$workload_contract_json, "sha256", serialize = FALSE)
      manifest_json <- .dsvert_joint_dp_client_json(list(
        protocol = "test-plan-manifest-v1", schema_sha256 = schema_hash,
        workload_contract_sha256 = workload_hash,
        variant = if (isTRUE(state$different_build) &&
                      identical(peer, tail(peers, 1L))) 2 else 1))
      manifest_hash <- digest::digest(
        manifest_json, "sha256", serialize = FALSE)
      index <- artifact_index(NULL, NULL, manifest_hash)
      unsigned <- list(
        version = .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_BUILD_VERSION,
        phase = "server_authoritative_manifest_memoized",
        peer_name = peer, peer_identity_pk = unname(pins[[peer]]),
        peer_pinset_sha256 = pin_hash, schema_sha256 = schema_hash,
        workload_contract_sha256 = workload_hash,
        manifest_sha256 = manifest_hash,
        manifest_bytes = nchar(manifest_json, type = "bytes"),
        capsule_id = digest::digest(
          manifest_json, "sha256", serialize = FALSE),
        privacy_epoch = statuses[[peer]]$noise_root$privacy_epoch,
        noise_key_id = statuses[[peer]]$noise_root$key_id,
        artifact_commitment_count = index$count,
        artifact_commitments_root = index$root,
        artifact_commitments = index$value, manifest_json = manifest_json,
        durable_memoization = TRUE, deterministic_replay = TRUE,
        data_access = FALSE, operation_limit = FALSE,
        request_limit = FALSE, history_can_deny_operation = FALSE)
      signed <- sign_object(unsigned, peer, "build")
      if (isTRUE(state$invalid_build_signature) &&
          identical(peer, tail(peers, 1L))) {
        signed$signature <- strrep("A", 86L)
      }
      .dsvert_joint_dp_client_json(signed)
    }), names(calls))
  }
  list(
    peers = peers, context = context, state = state, fanout = fanout,
    source_manifest = source_manifest, artifact_index = artifact_index)
}

.dp_capsule_plan_run <- function(fixture, active = FALSE) {
  blocked <- function(...) stop("release/protected route was called")
  input <- if (isTRUE(active)) NULL else fixture$context$all_conns
  testthat::with_mocked_bindings(
    .dsvert_dp_capsule_plan_impl(
      input, .aggregate = function(...) {
        stop("unexpected aggregate fallback")
      }),
    .dsvert_dp_datasources = function(value) {
      fixture$state$normalizer_input <- list(value)
      if (is.null(value)) fixture$context$all_conns else value
    },
    .dsvert_joint_dp_client_context = function(...) fixture$context,
    .dsvert_fanout_by_site = fixture$fanout,
    .dsvert_dp_capsule_source_manifest = fixture$source_manifest,
    .dsvert_dp_capsule_artifact_commitment_index_client =
      fixture$artifact_index,
    .dsvert_joint_dp_vector_capsule = blocked,
    .dsvert_dp_capsule_source_transport = blocked,
    .package = "dsVertClient")
}

test_that("capsule plan fans out every metadata phase for K=2,3,5", {
  shapes <- list()
  for (k in c(2L, 3L, 5L)) {
    fixture <- .dp_capsule_plan_fixture(k)
    result <- .dp_capsule_plan_run(fixture)
    endpoints <- vapply(fixture$state$calls, function(calls) {
      as.character(as.list(calls[[1L]])[[1L]])
    }, character(1L))

    expect_identical(endpoints, c(
      "dsvertDPCapsuleManifestDraftDS",
      "dsvertDPCapsuleManifestSignDS",
      "dsvertDPCapsuleManifestBuildDS"))
    expect_true(all(vapply(
      fixture$state$calls,
      function(calls) identical(names(calls), fixture$peers), logical(1L))))
    expect_s3_class(result, "ds.vertDPCapsulePlan")
    expect_identical(result$consortium$K, k)
    expect_identical(result$consortium$peers, fixture$peers)
    expect_false(result$guarantees$data_access)
    expect_false(result$guarantees$release_created)
    expect_identical(
      result$guarantees$privacy_cost, c(epsilon = 0, delta = 0))
    expect_false(result$guarantees$operation_limit)
    expect_false(result$guarantees$request_limit)
    expect_false(result$guarantees$history_can_deny_operation)
    expect_identical(result$projected_cost$projected_coordinate_count, 26)
    expect_identical(result$artifacts$coordinate_layout$coordinate_count, 26L)
    expect_identical(result$mechanism$mechanism, "discrete-laplace")
    expect_identical(result$calibration$winner, "laplace")
    shapes[[as.character(k)]] <- lapply(result, names)
  }
  expect_identical(shapes[["2"]], shapes[["3"]])
  expect_identical(shapes[["3"]], shapes[["5"]])
})

test_that("capsule plan fails closed on signatures, consensus, scope and layout", {
  signature <- .dp_capsule_plan_fixture(3L)
  signature$state$invalid_build_signature <- TRUE
  expect_error(
    .dp_capsule_plan_run(signature), "invalid capsule manifest signature")

  consensus <- .dp_capsule_plan_fixture(3L)
  consensus$state$different_build <- TRUE
  expect_error(
    .dp_capsule_plan_run(consensus),
    "built different biomedical capsule manifests")

  scope <- .dp_capsule_plan_fixture(2L)
  scope$state$tamper_scope <- TRUE
  expect_error(.dp_capsule_plan_run(scope), "primitive scope is invalid")

  layout <- .dp_capsule_plan_fixture(5L)
  layout$state$tamper_layout <- TRUE
  expect_error(.dp_capsule_plan_run(layout), "does not match its manifest")

  cost <- .dp_capsule_plan_fixture(3L)
  cost$state$tamper_cost <- TRUE
  expect_error(.dp_capsule_plan_run(cost), "projected workload cost")

  calibration <- .dp_capsule_plan_fixture(2L)
  calibration$state$tamper_calibration <- TRUE
  expect_error(
    .dp_capsule_plan_run(calibration), "mechanism or calibration contract")
})

test_that("capsule plan has a stable compact print contract", {
  result <- .dp_capsule_plan_run(.dp_capsule_plan_fixture(2L))
  expect_named(result, c(
    "version", "capsule", "primitive_scope", "selection",
    "projected_cost", "artifacts", "families", "mechanism",
    "calibration", "consortium", "guarantees"))
  expect_named(result$capsule, c(
    "capsule_id", "manifest_sha256", "logical_snapshot"))
  expect_named(result$selection, c(
    "sha256", "mode", "explicit_catalog_counts",
    "signed_spec_reference_counts", "included_counts"))
  expect_named(result$guarantees, c(
    "data_access", "release_created", "privacy_cost",
    "operation_limit", "request_limit", "history_can_deny_operation"))
  expect_named(result$families, c(
    "admitted_count", "numeric_moments", "numeric_pair_moments",
    "gaussian_models", "fixed_numeric_histograms",
    "categorical_marginals", "categorical_pairs",
    "correlation_artifacts", "describe_artifacts",
    "survival_artifacts"))
  printed <- capture.output(print(result))
  expect_lte(length(printed), 8L)
  expect_match(paste(printed, collapse = "\n"), "dry-run.*no data access")
  expect_false(grepl("schema_json|manifest_json|peer_identity_pk",
                     paste(printed, collapse = "\n")))
  expect_false(any(c("x", "z", "group", "status") %in%
                     unlist(result$selection, use.names = FALSE)))
})

test_that("public capsule plan is plug-and-play and hides transport injection", {
  expect_identical(names(formals(ds.vertDPCapsulePlan)),
                   c("datasources", "status"))
  expect_null(formals(ds.vertDPCapsulePlan)$datasources)
  observed <- NULL
  marker <- structure(list(ok = TRUE), class = "plan-wrapper-marker")
  result <- testthat::with_mocked_bindings(
    ds.vertDPCapsulePlan(status = "reusable-status"),
    .dsvert_dp_capsule_plan_impl = function(
        datasources, status, .aggregate) {
      observed <<- list(
        datasources = list(datasources), status = status,
        aggregate = .aggregate)
      marker
    },
    .package = "dsVertClient")
  expect_identical(result, marker)
  expect_null(observed$datasources[[1L]])
  expect_identical(observed$status, "reusable-status")
  expect_identical(observed$aggregate, DSI::datashield.aggregate)

  fixture <- .dp_capsule_plan_fixture(2L)
  active <- .dp_capsule_plan_run(fixture, active = TRUE)
  expect_s3_class(active, "ds.vertDPCapsulePlan")
  expect_null(fixture$state$normalizer_input[[1L]])
})
