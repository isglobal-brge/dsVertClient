.gaussian_alignment_complete <- function(
    manifest_json, context, layout, source_receipt, session_id, .aggregate) {
  source_count <- length(unlist(layout$source_peers))
  total <- .dsvert_dp_alignment_mask_private_projection_client(
    layout)$coordinate_count
  list(
    capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
    state = "complete",
    batch_operation_id = "op_88888888888888888888888888888888",
    contract_hash = source_receipt$contract_hash,
    source_count = source_count, coordinate_count = total,
    chunk_count = as.integer(ceiling(
      total / .dsvert_dp_alignment_mask_chunk_size_client(source_count))),
    alignment_digest_exposed = FALSE,
    mismatch_source_exposed = FALSE, gate_share_exposed = FALSE,
    fixed_transcript = TRUE)
}

.gaussian_cross_client_b64url <- function(value) {
  sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))),
    perl = TRUE)
}

test_that("cross Gaussian references bind an optional public owner", {
  qualified <- .dsvert_dp_gaussian_reference("site_02$age")
  expect_identical(qualified$server, "site_02")
  expect_identical(qualified$column, "age")
  expect_identical(.dsvert_dp_gaussian_reference("age")$server, NULL)
  expect_null(.dsvert_dp_gaussian_reference("site_02$other$age"))
  expect_null(.dsvert_dp_gaussian_reference("$age"))
  descriptor <- list(owner_peer = "site_02", column = "age")
  expect_true(.dsvert_dp_gaussian_reference_matches(
    "site_02$age", descriptor))
  expect_false(.dsvert_dp_gaussian_reference_matches(
    "site_03$age", descriptor))
})

.gaussian_cross_client_fixture <- function(bind_state = "bound", k = 3L) {
  stopifnot(k %in% 3:5)
  peers <- paste0("site_", letters[seq_len(k)])
  designated <- peers[1:2]
  keys <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  pins <- vapply(keys, function(key) {
    .gaussian_cross_client_b64url(tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  conns <- stats::setNames(lapply(peers, function(peer) {
    structure(list(peer = peer), class = "fake")
  }), peers)
  context <- list(
    servers = peers, designated = designated, pinset = pins,
    all_conns = conns, conns = conns[designated])
  capacity <- 2L
  scale <- 256
  q <- 2L
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
      column = "x", dataset = "predictor_data", owner_peer = "site_c",
      lower = -1, upper = 1)),
    predictor_order = list("x"), input_variable_order = list("x", "y"),
    participating_peers = as.list(c("site_b", "site_c")),
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
    statistic_maximum = as.list(rep(capacity * scale, coordinate_count)),
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
      ring_bits = 128, frac_bits = 8, required_signed_bits = 21,
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
  manifest <- list(
    admission = list(unit_capacity = capacity),
    workload = list(
      coordinate_count = 8,
      families = list(
        admitted_count = list(owner_peer = "site_a", dataset = "count_data"),
        numeric_moments = list(artifacts = list()),
        numeric_pair_moments = list(artifacts = list()),
        gaussian_models = list(artifacts = list(cross_model = artifact)),
        fixed_numeric_histograms = list(artifacts = list()),
        categorical_marginals = list(artifacts = list()),
        categorical_pairs = list(sets = list()),
        correlation_artifacts = list(), describe_artifacts = list(),
        survival_artifacts = list())),
    capsule_identity = list(capsule_id = strrep("a", 64L)))
  source_capsule_id <- strrep("d", 64L)
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
    c(unsigned, list(signature = .gaussian_cross_client_b64url(
      openssl::ed25519_sign(message, keys[[peer]]))))
  }
  bind <- stats::setNames(lapply(designated, function(peer) {
    .dsvert_joint_dp_client_json(sign_value(list(
      version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_BIND_VERSION,
      phase = "cross_gaussian_private_inputs_bound",
      capsule_id = source_capsule_id,
      analysis_id = "cross_model", artifact_sha256 = artifact_hash,
      source_contract_hash = source_hash,
      private_layout_sha256 = layout$transport_coordinate_order_sha256,
      transcript_sha256 = transcript_hash,
      numeric_certificate_sha256 = certificate_hash,
      peer_name = peer, peer_identity_pk = unname(pins[[peer]]),
      padded_units = capacity, variable_count = 2,
      ring_bits = 128, frac_bits = 8, state = bind_state,
      source_values_exposed = FALSE, alignment_hash_exposed = FALSE,
      alignment_hash_exposed_to_relay = FALSE,
      alignment_hash_exposed_to_computation_peers = FALSE,
      exact_intermediates_exposed = FALSE, fixed_transcript = TRUE),
      peer, "cross-gaussian-bind"))
  }), designated)
  result <- stats::setNames(lapply(designated, function(peer) {
    .dsvert_joint_dp_client_json(sign_value(list(
      version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_RECEIPT_VERSION,
      phase = "cross_gaussian_result_share_persisted",
      capsule_id = source_capsule_id,
      analysis_id = "cross_model", peer_name = peer,
      peer_identity_pk = unname(pins[[peer]]),
      artifact_sha256 = artifact_hash,
      source_contract_hash = source_hash,
      private_layout_sha256 = layout$transport_coordinate_order_sha256,
      transcript_sha256 = transcript_hash,
      numeric_certificate_sha256 = certificate_hash,
      exact_transcript_sha256 = exact_hash,
      coordinate_count = coordinate_count, release_start = 2,
      release_end = 8,
      release_coordinate_order_sha256 = release_layout$sha256,
      ring_bits = 128, frac_bits = 8, state = "complete",
      fixed_transcript = TRUE, result_share_exposed = FALSE,
      exact_intermediates_exposed = FALSE,
      alignment_hash_exposed = FALSE,
      alignment_hash_exposed_to_relay = FALSE,
      alignment_hash_exposed_to_computation_peers = FALSE),
      peer, "cross-gaussian-result"))
  }), designated)
  stage_response <- function(peer, expression) {
    stage <- as.character(expression$stage)
    index <- as.integer(expression$stage_index)
    contract <- .dsvert_dp_gaussian_cross_stage_contract(
      artifact, source_capsule_id, "cross_model",
      stage, index)
    handle_prefix <- if (identical(peer, designated[[1L]])) "A" else "B"
    handle <- paste0(handle_prefix, substr(digest::digest(
      paste(stage, index, peer), "sha256", serialize = FALSE), 1L, 42L))
    c(list(
      capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
      manifest_handle = handle, context_hash = strrep("1", 64L),
      plan_id = strrep("2", 64L), ring_bits = 128, frac_bits = 8,
      backend = "direct-wide", bound_x = "256", bound_y = "256",
      max_chunk = 16, total_n = contract$total_n,
      numeric_policy_id = strrep("3", 64L),
      version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_STAGE_VERSION,
      state = "prepared", producer = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_PRODUCER,
      purpose = contract$purpose,
      capsule_id = source_capsule_id,
      analysis_id = "cross_model", stage = stage, stage_index = index,
      artifact_sha256 = artifact_hash, source_contract_hash = source_hash,
      transcript_sha256 = transcript_hash,
      numeric_certificate_sha256 = certificate_hash,
      exact_intermediates_exposed = FALSE, source_values_exposed = FALSE))
  }
  state <- new.env(parent = emptyenv())
  state$commands <- character()
  aggregate <- function(conns, expr, ...) {
    values <- stats::setNames(lapply(names(conns), function(peer) {
      expression <- if (is.list(expr) && !is.call(expr)) expr[[peer]] else expr
      command <- as.character(expression[[1L]])
      state$commands <- c(state$commands, command)
      switch(command,
        dsvertDPGaussianCrossBindDS = bind[[peer]],
        dsvertDPGaussianCrossPrepareDS = stage_response(peer, expression),
        dsvertDPGaussianCrossFinalizeDS = result[[peer]],
        mpcCleanupDS = TRUE,
        stop("unexpected cross Gaussian client command: ", command))
    }), names(conns))
    values
  }
  source_receipt <- list(
    purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CROSS_PURPOSE,
    capsule_id = source_capsule_id,
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

test_that("cross Gaussian K=3 through K=5 separates owners and compute peers", {
  for (k in 3:5) {
    fixture <- .gaussian_cross_client_fixture(k = k)
    expect_length(fixture$context$all_conns, k)
    expect_true(fixture$layout$enabled)
    expect_identical(fixture$layout$release_coordinate_count, 8L)
    expect_identical(fixture$layout$private_start, 8193L)
    expect_identical(fixture$layout$transport_coordinate_count, 8200L)
    expect_identical(
      unlist(fixture$layout$source_peers), c("site_b", "site_c"))
    expect_identical(
      unlist(fixture$layout$computation_peers), c("site_a", "site_b"))
    expect_identical(
      fixture$layout$transport_coordinate_order_sha256,
      .dsvert_dp_capsule_source_hash(fixture$layout[
        setdiff(names(fixture$layout), c(
          "enabled", "transport_coordinate_order_sha256"))]))
  }
})

test_that("cross Gaussian descriptor and released moments match a central oracle", {
  fixture <- .gaussian_cross_client_fixture()
  artifact <- .dsvert_dp_gaussian_artifact(
    fixture$manifest, "outcome_data", "cross_model", "site_b",
    "add_remove_patient", 256, 2)
  expect_identical(artifact$participating_peers, c("site_b", "site_c"))
  expect_true(artifact$numeric_certificate$modular_wrap_proved_absent)
  moment <- .dsvert_dp_gaussian_unpack(
    c(2, 2, 1, 1, 1, 1, 1), artifact, 2)
  fit <- .dsvert_dp_gaussian_solve(moment, artifact, ridge = 0)
  expect_equal(fit$coefficients, c(`(Intercept)` = 0, x = 1),
               tolerance = 1e-12)

  tampered <- fixture$manifest
  tampered$workload$families$gaussian_models$artifacts$
    cross_model$numeric_certificate$modular_wrap_proved_absent <- FALSE
  expect_error(.dsvert_dp_gaussian_artifact(
    tampered, "outcome_data", "cross_model", "site_b",
    "add_remove_patient", 256, 2), "descriptor is invalid")
})

test_that("cross Gaussian K=3 through K=5 uses one fixed typed transcript", {
  for (k in 3:5) {
    fixture <- .gaussian_cross_client_fixture(k = k)
    setup <- 0L
    multiplications <- list()
    receipt <- .dsvert_dp_gaussian_cross_orchestrate(
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
      .alignment_mask = .gaussian_alignment_complete)
    expect_true(receipt$sampler_handoff_ready)
    expect_false(receipt$exact_intermediates_exposed)
    expect_false(receipt$source_values_exposed)
    expect_identical(setup, 1L)
    expect_identical(vapply(multiplications, `[[`, integer(1L), "total_n"),
                     c(2L, 4L, 6L))
    expect_true(all(vapply(multiplications, function(value) {
      setequal(names(value$manifests), c("site_a", "site_b"))
    }, logical(1L))))
    expect_false("exactGCVecmulBindInputsDS" %in% fixture$state$commands)
    expect_identical(sum(
      fixture$state$commands == "dsvertDPGaussianCrossPrepareDS"), 6L)
  }
})

test_that("Gaussian consumes a shared exact session without reopening it", {
  for (k in 3:5) {
    fixture <- .gaussian_cross_client_fixture(k = k)
    manifest_json <- "canonical-cross-manifest"
    session_id <- paste0(
      "2000000", k, "-0000-4000-8000-00000000000", k)
    alignment <- .gaussian_alignment_complete(
      manifest_json, fixture$context, fixture$layout,
      fixture$source_receipt, session_id, fixture$aggregate)
    shared <- .dsvert_dp_cross_shared_exact_build(
      manifest_json, fixture$manifest, fixture$context, fixture$layout,
      fixture$source_receipt, session_id, alignment)
    multiplications <- 0L
    receipt <- .dsvert_dp_gaussian_cross_orchestrate(
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
    expect_identical(multiplications, 3L)
    expect_false(any(fixture$state$commands == "mpcCleanupDS"))
  }
})

test_that("completed cross Gaussian results replay without recomputation", {
  fixture <- .gaussian_cross_client_fixture(bind_state = "complete")
  multiplications <- 0L
  receipt <- .dsvert_dp_gaussian_cross_orchestrate(
    "canonical-cross-manifest", fixture$manifest, fixture$context,
    fixture$source_receipt, fixture$aggregate,
    .setup_exact = function(...) invisible(list()),
    .vecmul = function(...) {
      multiplications <<- multiplications + 1L
      invisible(list())
    },
    .alignment_mask = .gaussian_alignment_complete)
  expect_true(receipt$sampler_handoff_ready)
  expect_identical(multiplications, 0L)
  expect_false(any(
    fixture$state$commands == "dsvertDPGaussianCrossPrepareDS"))
})

test_that("cross Gaussian result disagreement and source misbinding fail closed", {
  fixture <- .gaussian_cross_client_fixture()
  bad_source <- fixture$source_receipt
  bad_source$private_layout_sha256 <- strrep("0", 64L)
  expect_error(.dsvert_dp_gaussian_cross_orchestrate(
    "canonical-cross-manifest", fixture$manifest, fixture$context,
    bad_source, fixture$aggregate,
    .setup_exact = function(...) invisible(list()),
    .vecmul = function(...) invisible(list())), "source transport is misbound")

  fixture <- .gaussian_cross_client_fixture(bind_state = "complete")
  decoded <- .dsvert_joint_dp_client_decode(
    fixture$result[["site_b"]], "test Gaussian result", 128L * 1024L)
  decoded$exact_transcript_sha256 <- strrep("d", 64L)
  # A changed value with the old signature must be rejected before consensus.
  bad_result <- fixture$result
  bad_result[["site_b"]] <- .dsvert_joint_dp_client_json(decoded)
  binding <- .dsvert_dp_gaussian_cross_bind_set(
    fixture$bind, fixture$context, fixture$manifest, fixture$layout,
    "cross_model", fixture$source_receipt)
  expect_error(.dsvert_dp_gaussian_cross_result_set(
    bad_result, fixture$context, fixture$manifest, fixture$layout,
    binding, "cross_model"), "invalid biomedical capsule source signature")
})
