.categorical_alignment_complete <- function(
    manifest_json, context, layout, source_receipt, session_id, .aggregate) {
  source_count <- length(unlist(layout$source_peers))
  total <- as.numeric(layout$transport_coordinate_count)
  list(
    capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
    state = "complete",
    batch_operation_id = "op_99999999999999999999999999999999",
    contract_hash = source_receipt$contract_hash,
    source_count = source_count, coordinate_count = total,
    chunk_count = as.integer(ceiling(
      total / .dsvert_dp_alignment_mask_chunk_size_client(source_count))),
    alignment_digest_exposed = FALSE,
    mismatch_source_exposed = FALSE, gate_share_exposed = FALSE,
    fixed_transcript = TRUE)
}

.categorical_cross_client_b64url <- function(value) {
  sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))),
    perl = TRUE)
}

.categorical_cross_client_fixture <- function(
    bind_state = "bound", k = 3L) {
  stopifnot(k %in% 2:5)
  peers <- paste0("site_", letters[seq_len(k)])
  designated <- peers[1:2]
  participants <- if (k == 2L) peers else tail(peers, 2L)
  keys <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  pins <- vapply(keys, function(key) {
    .categorical_cross_client_b64url(
      tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  conns <- stats::setNames(lapply(peers, function(peer) {
    structure(list(peer = peer), class = "fake")
  }), peers)
  context <- list(
    servers = peers, designated = designated, pinset = pins,
    all_conns = conns, conns = conns[designated],
    allocation_openings = stats::setNames(lapply(
      designated, function(peer) .dsvert_joint_dp_client_json(list(
        version = "test-allocation-opening-v1", peer_name = peer))),
      designated))
  capacity <- 2L
  scale <- 256
  coordinate_count <- 6L
  artifact <- list(
    version = .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION,
    spec_version = "v2", analysis_id = "cross_table",
    alignment_group = "aligned-main",
    left = list(
      dataset = "left_data", column = "left_cat",
      owner_peer = participants[[1L]], levels = list("A", "B")),
    right = list(
      dataset = "right_data", column = "right_cat",
      owner_peer = participants[[2L]], levels = list("X", "Y", "Z")),
    participating_peers = as.list(participants),
    computation_peers = as.list(designated),
    coordinate_count = coordinate_count,
    coordinate_order = paste0(
      "canonical_left_level_rows_then_canonical_right_level_columns_",
      "column_major_v1"),
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    private_input_layout = paste0(
      "capacity_padded_one_hot_by_public_level_then_side_manifest_order_v1"),
    repeated_record_policy = paste0(
      "per_owner_consistent_level_else_zero_then_one_joint_cell_per_",
      "admitted_unit_v1"),
    missingness_policy = paste0(
      "missing_out_of_domain_or_conflicting_side_is_all_zero_and_",
      "contributes_no_joint_cell_v1"),
    statistic_maximum = as.list(rep(capacity * scale, coordinate_count)),
    add_remove_l1_sensitivity = 1,
    add_remove_l2_sensitivity = 1,
    replace_one_l1_sensitivity = 2,
    replace_one_l2_sensitivity = sqrt(2) * (1 + 32 * .Machine$double.eps),
    selected_l1_sensitivity = 1,
    selected_l2_sensitivity = 1,
    source_raw_l1_sensitivity = scale,
    source_raw_l2_sensitivity = scale,
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste0(
      "one_hot_joint_cell_add_remove_and_difference_of_two_cells_",
      "replace_one_v1"),
    numeric_certificate = list(
      version = "dsvert-cross-categorical-numeric-certificate-v1",
      ring_bits = 128, frac_bits = 8, required_signed_bits = 13,
      operand_maximum = scale, product_maximum = scale,
      accumulated_coordinate_maximum = capacity * scale,
      truncation = "exact_signed_floor_gc_ot_v1",
      modular_wrap_proved_absent = TRUE,
      overflow_behavior = "typed_abort_before_commit"),
    transcript = list(
      version = "dsvert-cross-categorical-fixed-transcript-v1",
      padded_units = capacity, row_level_count = 2,
      column_level_count = 3, product_coordinate_count = 12,
      exact_multiplication_rounds = 1,
      data_dependent_branches = 0,
      exact_intermediate_release_count = 0),
    alignment_contract = list(
      version = "private-psi-ordered-manifest-consensus-v1",
      public_patient_dependent_hash = FALSE,
      mismatch_behavior = "typed_non_prealigned_cohort_failure"),
    implementation_state = "cross_owner_exact_gc_materialized",
    cross_owner_state = "exact_gc_to_joint_dp_vector_v1")
  manifest <- list(
    admission = list(unit_capacity = capacity),
    workload = list(
      coordinate_count = 7,
      families = list(
        admitted_count = list(
          owner_peer = "site_a", dataset = "count_data"),
        numeric_moments = list(artifacts = list()),
        numeric_pair_moments = list(artifacts = list()),
        gaussian_models = list(artifacts = list()),
        fixed_numeric_histograms = list(artifacts = list()),
        categorical_marginals = list(artifacts = list()),
        categorical_pairs = list(
          sets = list(), cross_artifacts = list(cross_table = artifact)),
        correlation_artifacts = list(), describe_artifacts = list(),
        survival_artifacts = list())),
    capsule_identity = list(capsule_id = strrep("a", 64L)))
  layout <- .dsvert_dp_gaussian_cross_layout_client(manifest)
  release_layout <- .dsvert_dp_capsule_vector_layout(manifest)
  source_hash <- strrep("b", 64L)
  artifact_hash <- .dsvert_dp_capsule_source_hash(artifact)
  transcript_hash <- .dsvert_dp_capsule_source_hash(artifact$transcript)
  certificate_hash <- .dsvert_dp_capsule_source_hash(
    artifact$numeric_certificate)
  exact_hash <- strrep("c", 64L)
  sign_value <- function(unsigned, peer, domain) {
    message <- charToRaw(paste0(
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SIGNATURE_DOMAIN, domain, "|",
      .dsvert_joint_dp_client_json(unsigned)))
    c(unsigned, list(signature = .categorical_cross_client_b64url(
      openssl::ed25519_sign(message, keys[[peer]]))))
  }
  bind <- stats::setNames(lapply(designated, function(peer) {
    .dsvert_joint_dp_client_json(sign_value(list(
      version = .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_BIND_VERSION,
      phase = "cross_categorical_private_inputs_bound",
      capsule_id = manifest$capsule_identity$capsule_id,
      analysis_id = "cross_table", artifact_sha256 = artifact_hash,
      source_contract_hash = source_hash,
      private_layout_sha256 = layout$transport_coordinate_order_sha256,
      transcript_sha256 = transcript_hash,
      numeric_certificate_sha256 = certificate_hash,
      peer_name = peer, peer_identity_pk = unname(pins[[peer]]),
      padded_units = capacity, row_level_count = 2,
      column_level_count = 3, ring_bits = 128, frac_bits = 8,
      state = bind_state, source_values_exposed = FALSE,
      alignment_hash_exposed = FALSE,
      alignment_hash_exposed_to_relay = FALSE,
      alignment_hash_exposed_to_computation_peers = FALSE,
      exact_intermediates_exposed = FALSE,
      fixed_transcript = TRUE), peer, "cross-categorical-bind"))
  }), designated)
  result <- stats::setNames(lapply(designated, function(peer) {
    .dsvert_joint_dp_client_json(sign_value(list(
      version = .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_RECEIPT_VERSION,
      phase = "cross_categorical_result_share_persisted",
      capsule_id = manifest$capsule_identity$capsule_id,
      analysis_id = "cross_table", peer_name = peer,
      peer_identity_pk = unname(pins[[peer]]),
      artifact_sha256 = artifact_hash,
      source_contract_hash = source_hash,
      private_layout_sha256 = layout$transport_coordinate_order_sha256,
      transcript_sha256 = transcript_hash,
      numeric_certificate_sha256 = certificate_hash,
      exact_transcript_sha256 = exact_hash,
      coordinate_count = coordinate_count, release_start = 2,
      release_end = 7,
      release_coordinate_order_sha256 = release_layout$sha256,
      ring_bits = 128, frac_bits = 8, state = "complete",
      fixed_transcript = TRUE, result_share_exposed = FALSE,
      exact_intermediates_exposed = FALSE,
      alignment_hash_exposed = FALSE,
      alignment_hash_exposed_to_relay = FALSE,
      alignment_hash_exposed_to_computation_peers = FALSE),
      peer, "cross-categorical-result"))
  }), designated)
  stage_response <- function(peer) {
    contract <- .dsvert_dp_categorical_cross_stage_contract(
      artifact, manifest$capsule_identity$capsule_id, "cross_table")
    prefix <- if (identical(peer, designated[[1L]])) "A" else "B"
    handle <- paste0(prefix, substr(digest::digest(
      peer, "sha256", serialize = FALSE), 1L, 42L))
    list(
      capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
      manifest_handle = handle, context_hash = strrep("1", 64L),
      plan_id = strrep("2", 64L), ring_bits = 128, frac_bits = 8,
      backend = "direct-wide", bound_x = "256", bound_y = "256",
      max_chunk = 16, total_n = contract$total_n,
      numeric_policy_id = strrep("3", 64L),
      version = .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_STAGE_VERSION,
      state = "prepared",
      producer = .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_PRODUCER,
      purpose = contract$purpose,
      capsule_id = manifest$capsule_identity$capsule_id,
      analysis_id = "cross_table", stage = "cell-products",
      stage_index = 1, artifact_sha256 = artifact_hash,
      source_contract_hash = source_hash,
      transcript_sha256 = transcript_hash,
      numeric_certificate_sha256 = certificate_hash,
      exact_intermediates_exposed = FALSE, source_values_exposed = FALSE)
  }
  state <- new.env(parent = emptyenv())
  state$commands <- character()
  state$calls <- list()
  aggregate <- function(conns, expr, ...) {
    stats::setNames(lapply(names(conns), function(peer) {
      expression <- if (is.list(expr) && !is.call(expr)) expr[[peer]] else expr
      command <- as.character(expression[[1L]])
      state$commands <- c(state$commands, command)
      state$calls[[length(state$calls) + 1L]] <- list(
        peer = peer, expression = expression)
      switch(command,
        dsvertDPCategoricalCrossBindDS = bind[[peer]],
        dsvertDPCategoricalCrossPrepareDS = stage_response(peer),
        dsvertDPCategoricalCrossFinalizeDS = result[[peer]],
        mpcCleanupDS = TRUE,
        stop("unexpected categorical command: ", command))
    }), names(conns))
  }
  source_receipt <- list(
    purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CATEGORICAL_CROSS_PURPOSE,
    capsule_id = manifest$capsule_identity$capsule_id,
    contract_hash = source_hash,
    coordinate_count = layout$transport_coordinate_count,
    release_coordinate_count = layout$release_coordinate_count,
    private_layout_sha256 = layout$transport_coordinate_order_sha256,
    sampler_handoff_ready = FALSE, payload_exposed = FALSE)
  list(
    context = context, manifest = manifest, layout = layout,
    source_receipt = source_receipt, aggregate = aggregate, state = state,
    artifact = artifact, bind = bind, result = result)
}

test_that("categorical K=2 through K=4 separates source and compute peers", {
  for (k in 2:4) {
    fixture <- .categorical_cross_client_fixture(k = k)
    expect_length(fixture$context$all_conns, k)
    expect_identical(
      unlist(fixture$layout$computation_peers), c("site_a", "site_b"))
    expect_identical(
      unlist(fixture$layout$source_peers),
      if (k == 2L) c("site_a", "site_b") else
        tail(fixture$context$servers, 2L))
    expect_identical(fixture$layout$release_coordinate_count, 7L)
    expect_identical(fixture$layout$private_start, 8193L)
  }
})

test_that("categorical binding fails closed without allocation openings", {
  fixture <- .categorical_cross_client_fixture()
  fixture$context$allocation_openings <- NULL
  expect_error(.dsvert_dp_categorical_cross_orchestrate(
    "canonical-cross-manifest", fixture$manifest, fixture$context,
    fixture$source_receipt, fixture$aggregate,
    .setup_exact = function(...) stop("exact setup must not start"),
    .vecmul = function(...) stop("exact multiplication must not start")),
    "no canonical cross-signed allocation")
  expect_length(fixture$state$commands, 0L)
})

test_that("categorical uses one fixed exact multiplication per artifact", {
  for (k in 2:4) {
    fixture <- .categorical_cross_client_fixture(k = k)
    setup <- 0L
    multiplications <- list()
    receipt <- .dsvert_dp_categorical_cross_orchestrate(
      "canonical-cross-manifest", fixture$manifest, fixture$context,
      fixture$source_receipt, fixture$aggregate,
      .setup_exact = function(datasources, server_names, servers, session_id,
                              .aggregate) {
        setup <<- setup + 1L
        expect_identical(server_names[servers], c("site_a", "site_b"))
        invisible(list())
      },
      .vecmul = function(..., input_manifests, total_n) {
        multiplications[[length(multiplications) + 1L]] <<- list(
          manifests = input_manifests, total_n = total_n)
        invisible(list())
      },
      .alignment_mask = .categorical_alignment_complete)
    expect_true(receipt$sampler_handoff_ready)
    expect_false(receipt$exact_intermediates_exposed)
    expect_false(receipt$source_values_exposed)
    expect_identical(setup, 1L)
    expect_length(multiplications, 1L)
    expect_identical(multiplications[[1L]]$total_n, 12L)
    expect_setequal(names(multiplications[[1L]]$manifests),
                    c("site_a", "site_b"))
    expect_identical(sum(fixture$state$commands ==
                           "dsvertDPCategoricalCrossPrepareDS"), 2L)
    bind_calls <- Filter(function(value) identical(
      as.character(value$expression[[1L]]),
      "dsvertDPCategoricalCrossBindDS"), fixture$state$calls)
    expect_length(bind_calls, 2L)
    for (value in bind_calls) {
      arguments <- as.list(value$expression)[-1L]
      expect_identical(arguments$first_opening_json,
                       .dsvert_dsi_text_encode(
                         fixture$context$allocation_openings[[1L]]))
      expect_identical(arguments$second_opening_json,
                       .dsvert_dsi_text_encode(
                         fixture$context$allocation_openings[[2L]]))
    }
    expect_false(any(grepl("exactGCChisq|dsvertOneHot|dsvertColNames",
                           fixture$state$commands)))
  }
})

test_that("categorical consumes a shared exact session without reopening it", {
  for (k in 2:5) {
    fixture <- .categorical_cross_client_fixture(k = k)
    manifest_json <- "canonical-cross-manifest"
    session_id <- paste0(
      "1000000", k, "-0000-4000-8000-00000000000", k)
    alignment <- .categorical_alignment_complete(
      manifest_json, fixture$context, fixture$layout,
      fixture$source_receipt, session_id, fixture$aggregate)
    shared <- .dsvert_dp_cross_shared_exact_build(
      manifest_json, fixture$manifest, fixture$context, fixture$layout,
      fixture$source_receipt, session_id, alignment)
    multiplications <- 0L
    receipt <- .dsvert_dp_categorical_cross_orchestrate(
      manifest_json, fixture$manifest, fixture$context,
      fixture$source_receipt, fixture$aggregate,
      .setup_exact = function(...) stop("duplicate exact setup"),
      .alignment_mask = function(...) stop("duplicate alignment gate"),
      .vecmul = function(...) {
        multiplications <<- multiplications + 1L
        invisible(list())
      },
      .shared_exact = shared)
    expect_true(receipt$sampler_handoff_ready)
    expect_identical(multiplications, 1L)
    expect_false(any(fixture$state$commands == "mpcCleanupDS"))
  }
})

test_that("completed categorical result replays and tampering fails closed", {
  fixture <- .categorical_cross_client_fixture(bind_state = "complete")
  multiplications <- 0L
  receipt <- .dsvert_dp_categorical_cross_orchestrate(
    "canonical-cross-manifest", fixture$manifest, fixture$context,
    fixture$source_receipt, fixture$aggregate,
    .setup_exact = function(...) invisible(list()),
    .vecmul = function(...) {
      multiplications <<- multiplications + 1L
      invisible(list())
    },
    .alignment_mask = .categorical_alignment_complete)
  expect_true(receipt$sampler_handoff_ready)
  expect_identical(multiplications, 0L)
  expect_false(any(fixture$state$commands ==
                     "dsvertDPCategoricalCrossPrepareDS"))

  decoded <- .dsvert_joint_dp_client_decode(
    fixture$result[["site_b"]], "categorical test result", 128L * 1024L)
  decoded$exact_transcript_sha256 <- strrep("d", 64L)
  bad_result <- fixture$result
  bad_result[["site_b"]] <- .dsvert_joint_dp_client_json(decoded)
  binding <- .dsvert_dp_categorical_cross_bind_set(
    fixture$bind, fixture$context, fixture$manifest, fixture$layout,
    "cross_table")
  expect_error(.dsvert_dp_categorical_cross_result_set(
    bad_result, fixture$context, fixture$manifest, fixture$layout,
    binding, "cross_table"), "invalid biomedical capsule source signature")

  bad_source <- fixture$source_receipt
  bad_source$private_layout_sha256 <- strrep("0", 64L)
  expect_error(.dsvert_dp_categorical_cross_orchestrate(
    "canonical-cross-manifest", fixture$manifest, fixture$context,
    bad_source, fixture$aggregate,
    .setup_exact = function(...) invisible(list()),
    .vecmul = function(...) invisible(list())),
    "source transport is misbound")
})
