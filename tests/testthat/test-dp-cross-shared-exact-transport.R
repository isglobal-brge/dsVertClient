.shared_cross_exact_fixture <- function(k) {
  stopifnot(k %in% 2:5)
  servers <- paste0("site_", letters[seq_len(k)])
  peers <- servers[1:2]
  pins <- stats::setNames(vapply(seq_along(servers), function(index) {
    paste0("pin_", index, "_", strrep(letters[[index]], 40L))
  }, character(1L)), servers)
  context <- list(
    servers = servers, designated = peers, pinset = pins,
    all_conns = stats::setNames(as.list(seq_along(servers)), servers),
    conns = stats::setNames(as.list(seq_along(peers)), peers),
    allocation_openings = stats::setNames(as.list(
      paste0("opening-", seq_along(peers))), peers))
  private_layout <- digest::digest(
    paste0("private-layout-k", k), algo = "sha256", serialize = FALSE)
  source_count <- length(servers)
  chunk_size <- .dsvert_dp_alignment_mask_chunk_size_client(source_count)
  total <- as.numeric(chunk_size + 3L)
  layout <- list(
    computation_peers = as.list(peers),
    source_peers = as.list(servers),
    transport_coordinate_count = total,
    release_coordinate_count = 17,
    transport_coordinate_order_sha256 = private_layout)
  capsule_id <- digest::digest(
    paste0("capsule-k", k), algo = "sha256", serialize = FALSE)
  contract_hash <- digest::digest(
    paste0("source-contract-k", k), algo = "sha256", serialize = FALSE)
  manifest <- list(capsule_identity = list(capsule_id = capsule_id))
  manifest_json <- as.character(jsonlite::toJSON(
    manifest, auto_unbox = TRUE, null = "null"))
  source_receipt <- list(
    purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CATEGORICAL_CROSS_PURPOSE,
    capsule_id = capsule_id, contract_hash = contract_hash,
    coordinate_count = total, release_coordinate_count = 17,
    private_layout_sha256 = private_layout,
    sampler_handoff_ready = FALSE, payload_exposed = FALSE)
  alignment <- list(
    capability_id = .DSVERT_CLIENT_EXACT_GC_CAPABILITY,
    state = "complete",
    batch_operation_id = paste0("op_", strrep(as.character(k), 32L)),
    contract_hash = contract_hash, source_count = source_count,
    coordinate_count = total, chunk_count = 2L,
    alignment_digest_exposed = FALSE, mismatch_source_exposed = FALSE,
    gate_share_exposed = FALSE, fixed_transcript = TRUE)
  list(
    k = k, context = context, layout = layout, manifest = manifest,
    manifest_json = manifest_json, source_receipt = source_receipt,
    alignment = alignment,
    session_id = paste0(
      "0000000", k, "-0000-4000-8000-00000000000", k))
}

test_that("shared exact context is replay-exact and capsule-bound for K=2/3/4/5", {
  for (k in 2:5) {
    fixture <- .shared_cross_exact_fixture(k)
    shared <- .dsvert_dp_cross_shared_exact_build(
      fixture$manifest_json, fixture$manifest, fixture$context,
      fixture$layout, fixture$source_receipt, fixture$session_id,
      fixture$alignment)
    replay <- .dsvert_dp_cross_shared_exact_validate(
      shared, fixture$manifest_json, fixture$manifest, fixture$context,
      fixture$layout, fixture$source_receipt)
    expect_identical(replay, shared)
    expect_true(shared$alignment_gate_complete)
    expect_false(shared$alignment_digest_exposed)
    expect_false(shared$mismatch_source_exposed)
    expect_false(shared$gate_share_exposed)
    expect_false(shared$exact_intermediates_exposed)
    expect_identical(
      shared$family_operation_domains,
      .DSVERT_CLIENT_DP_CROSS_SHARED_EXACT_DOMAINS)
    expect_false(any(c(
      "alignment_hash", "alignment_share", "validity_share",
      "source_values", "exact_intermediates") %in% names(shared)))

    tampered <- shared
    tampered$source_contract_hash <- strrep("f", 64L)
    expect_error(.dsvert_dp_cross_shared_exact_validate(
      tampered, fixture$manifest_json, fixture$manifest, fixture$context,
      fixture$layout, fixture$source_receipt), "not bound")
    tampered <- shared
    tampered$alignment_batch_operation_id <-
      "op_ffffffffffffffffffffffffffffffff"
    # A syntactically valid but different batch cannot grant server authority;
    # bind endpoints still require the completed capsule/contract batch in the
    # authenticated session. The client token remains a liveness hint only.
    expect_identical(.dsvert_dp_cross_shared_exact_validate(
      tampered, fixture$manifest_json, fixture$manifest, fixture$context,
      fixture$layout, fixture$source_receipt), tampered)
    tampered <- shared
    tampered$alignment_batch_operation_id <- "op_bad"
    expect_error(.dsvert_dp_cross_shared_exact_validate(
      tampered, fixture$manifest_json, fixture$manifest, fixture$context,
      fixture$layout, fixture$source_receipt), "not bound")
    tampered <- shared
    tampered$family_operation_domains <- "mixed"
    expect_error(.dsvert_dp_cross_shared_exact_validate(
      tampered, fixture$manifest_json, fixture$manifest, fixture$context,
      fixture$layout, fixture$source_receipt), "not bound")
    tampered <- c(shared, list(private_gate = TRUE))
    expect_error(.dsvert_dp_cross_shared_exact_validate(
      tampered, fixture$manifest_json, fixture$manifest, fixture$context,
      fixture$layout, fixture$source_receipt), "not bound")

    changed_context <- fixture$context
    changed_context$pinset[[fixture$context$designated[[1L]]]] <-
      "different-pinned-identity"
    expect_error(.dsvert_dp_cross_shared_exact_validate(
      shared, fixture$manifest_json, fixture$manifest, changed_context,
      fixture$layout, fixture$source_receipt), "not bound")
    changed_layout <- fixture$layout
    changed_layout$source_peers <- rev(changed_layout$source_peers)
    expect_error(.dsvert_dp_cross_shared_exact_validate(
      shared, fixture$manifest_json, fixture$manifest, fixture$context,
      changed_layout, fixture$source_receipt), "not bound")
    expect_error(.dsvert_dp_cross_shared_exact_validate(
      shared, paste0(fixture$manifest_json, " "), fixture$manifest,
      fixture$context, fixture$layout, fixture$source_receipt), "not bound")
  }
})

test_that("one shared session and gate serve both families for K=2/3/4/5", {
  for (k in 2:5) {
    fixture <- .shared_cross_exact_fixture(k)
    setup_count <- gate_count <- cleanup_count <- 0L
    family_sessions <- character()
    setup <- function(datasources, server_names, servers, session_id,
                      .aggregate) {
      setup_count <<- setup_count + 1L
      expect_identical(server_names[servers], fixture$context$designated)
      expect_identical(session_id, fixture$session_id)
      invisible(list())
    }
    gate <- function(manifest_json, context, layout, source_receipt,
                     session_id, .aggregate) {
      gate_count <<- gate_count + 1L
      expect_identical(session_id, fixture$session_id)
      fixture$alignment
    }
    family <- function(
        manifest_json, manifest, context, source_receipt, .aggregate,
        .setup_exact, .alignment_mask, .shared_exact) {
      validated <- .dsvert_dp_cross_shared_exact_validate(
        .shared_exact, manifest_json, manifest, context, fixture$layout,
        source_receipt)
      family_sessions <<- c(family_sessions, validated$session_id)
      list(
        enabled = TRUE, sampler_handoff_ready = TRUE,
        exact_intermediates_exposed = FALSE,
        source_values_exposed = FALSE)
    }
    aggregate <- function(...) stop("unexpected aggregate call")
    testthat::local_mocked_bindings(
      .dsvert_dp_categorical_cross_artifacts_client = function(manifest) {
        list(categorical = TRUE)
      },
      .dsvert_dp_gaussian_cross_artifacts_client = function(manifest) {
        list(gaussian = TRUE)
      },
      .dsvert_dp_gaussian_cross_layout_client = function(manifest) {
        fixture$layout
      },
      .dsvert_dp_capsule_allocation_openings = function(context) {
        context$allocation_openings
      },
      .dsvert_dp_categorical_cross_orchestrate = family,
      .dsvert_dp_gaussian_cross_orchestrate = family,
      .dsvert_uuid4 = function() fixture$session_id,
      .dsvert_cleanup_best_effort = function(...) {
        cleanup_count <<- cleanup_count + 1L
        invisible(NULL)
      },
      .package = "dsVertClient")
    result <- .dsvert_dp_cross_orchestrate(
      fixture$manifest_json, fixture$manifest, fixture$context,
      fixture$source_receipt, aggregate,
      .setup_exact = setup, .alignment_mask = gate)
    expect_true(result$enabled)
    expect_identical(result$families, c("categorical", "gaussian"))
    expect_true(result$shared_exact_session)
    expect_identical(result$exact_session_count, 1L)
    expect_identical(result$alignment_gate_count, 1L)
    expect_identical(setup_count, 1L)
    expect_identical(gate_count, 1L)
    expect_identical(cleanup_count, 0L)
    expect_identical(family_sessions,
                     rep(fixture$session_id, 2L))
  }
})

test_that("tampered shared gate aborts before either family starts", {
  fixture <- .shared_cross_exact_fixture(5L)
  fixture$alignment$contract_hash <- strrep("e", 64L)
  family_calls <- 0L
  family <- function(...) {
    family_calls <<- family_calls + 1L
    stop("family must not start")
  }
  cleanup_count <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_categorical_cross_artifacts_client = function(manifest) {
      list(categorical = TRUE)
    },
    .dsvert_dp_gaussian_cross_artifacts_client = function(manifest) {
      list(gaussian = TRUE)
    },
    .dsvert_dp_gaussian_cross_layout_client = function(manifest) {
      fixture$layout
    },
    .dsvert_dp_capsule_allocation_openings = function(context) {
      context$allocation_openings
    },
    .dsvert_dp_categorical_cross_orchestrate = family,
    .dsvert_dp_gaussian_cross_orchestrate = family,
    .dsvert_uuid4 = function() fixture$session_id,
    .dsvert_cleanup_best_effort = function(...) {
      cleanup_count <<- cleanup_count + 1L
      invisible(NULL)
    },
    .package = "dsVertClient")
  expect_error(.dsvert_dp_cross_orchestrate(
    fixture$manifest_json, fixture$manifest, fixture$context,
    fixture$source_receipt, function(...) NULL,
    .setup_exact = function(...) invisible(list()),
    .alignment_mask = function(...) fixture$alignment),
    "not capsule-bound")
  expect_identical(family_calls, 0L)
  expect_identical(cleanup_count, 0L)
})

test_that("family exact domains remain distinct inside a shared session", {
  expect_false(identical(
    .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_PRODUCER,
    .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_PRODUCER))
  capsule_id <- strrep("a", 64L)
  analysis_id <- "same_analysis_id"
  categorical <- .dsvert_dp_categorical_cross_stage_contract(list(
    transcript = list(padded_units = 2),
    left = list(levels = list("A", "B")),
    right = list(levels = list("X", "Y"))),
    capsule_id, analysis_id)
  gaussian <- .dsvert_dp_gaussian_cross_stage_contract(list(
    transcript = list(padded_units = 2),
    input_variable_order = list("x", "y")),
    capsule_id, analysis_id, "validity", 1L)
  expect_match(categorical$purpose, "dp.categorical-cross", fixed = TRUE)
  expect_match(gaussian$purpose, "dp.gaussian-cross", fixed = TRUE)
  expect_false(identical(categorical$purpose, gaussian$purpose))
})
