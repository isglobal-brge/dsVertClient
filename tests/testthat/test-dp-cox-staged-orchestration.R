.cox_staged_manifest <- function(f) {
  artifact <- .dsvert_dp_cox_cross_workload_artifact(f$contract)
  families <- setNames(rep(list(list(artifacts = list())), 10L), c(
    "admitted_count", "numeric_moments", "numeric_pair_moments", "gaussian_models",
    "fixed_numeric_histograms", "categorical_marginals", "categorical_pairs",
    "correlation_artifacts", "describe_artifacts", "survival_artifacts"))
  families$admitted_count <- list(owner_peer = "site_a", dataset = "aligned")
  families$survival_artifacts <- list()
  families$gaussian_models$artifacts <- list(cox_grid = artifact)
  list(logical_snapshot = f$contract$spec$logical_snapshot,
    capsule_identity = list(capsule_id = strrep("6", 64)),
    admission = list(unit_capacity = f$policy$unit_capacity, adjacency = f$policy$adjacency),
    bounds = list(numeric_grid_bits = f$policy$numeric_grid_bits),
    workload = list(coordinate_count = artifact$coordinate_count + 1L,
      families = families, capsule_mechanism = list(source_context_hash = strrep("a", 64))))
}

test_that("Cox transport commitment agrees with the server at every topology", {
  skip_if_not_installed("dsVert")
  server <- get(".dsvert_dp_cox_cross_transport_layout", asNamespace("dsVert"))
  for (owners in c(2L, 3L, 5L)) {
    f <- .cox_cross_client_fixture(capacity = 5, owners = owners)
    manifest <- .cox_staged_manifest(f)
    artifact <- manifest$workload$families$gaussian_models$artifacts[[1L]]
    expect_identical(.dsvert_dp_cox_cross_transport_layout(manifest, artifact),
      server(manifest, artifact))
    extra <- manifest
    extra$workload$families$gaussian_models$artifacts$other <- artifact
    expect_error(.dsvert_dp_cox_cross_transport_layout(extra, artifact))
  }
})

test_that("Cox signed orchestration preserves bilateral persistence and terminal binding", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .cox_cross_client_fixture(capacity = 4, owners = owners)
    manifest <- .cox_staged_manifest(f)
    artifact <- manifest$workload$families$gaussian_models$artifacts[[1L]]
    layout <- .dsvert_dp_cox_cross_transport_layout(manifest, artifact)
    peers <- c("site_a", "site_b")
    context <- list(designated = rev(peers), pinset = f$policy$peer_pinset,
      conns = setNames(as.list(peers), peers), all_conns = as.list(names(f$keys)),
      servers = names(f$keys))
    source <- list(purpose = .DSVERT_CLIENT_DP_GLM_GRID_CROSS_SOURCE_PURPOSE,
      capsule_id = strrep("a", 64), contract_hash = strrep("b", 64),
      coordinate_count = layout$transport_coordinate_count,
      release_coordinate_count = layout$release_coordinate_count,
      private_layout_sha256 = layout$transport_coordinate_order_sha256,
      sampler_handoff_ready = FALSE, payload_exposed = FALSE)
    semantic <- .dsvert_dp_capsule_source_hash(list(
      version = "cross-grid-semantic-release-key-v2", capsule_id = source$capsule_id,
      source_contract_sha256 = source$contract_hash, signed_contract = artifact$signed_contract,
      family = artifact$family, version_family = artifact$spec_version,
      mechanism = manifest$workload$capsule_mechanism, alignment = f$contract$spec$alignment,
      caps = artifact$sensitivity))
    sign <- function(value, peer) {
      value$signature <- f$b64(openssl::ed25519_sign(charToRaw(paste0(
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SIGNATURE_DOMAIN, "cross-grid-result|",
        .dsvert_joint_dp_client_json(value))), f$keys[[peer]]))
      value
    }
    persisted <- c(site_a = FALSE, site_b = FALSE)
    change <- NULL
    actions <- character(); runs <- 0L; cleanups <- 0L
    receipt <- function(phase, peer) {
      value <- list(version = "dsvert-cox-staged-receipt-v1", phase = phase,
        capsule_id = source$capsule_id, analysis_id = artifact$analysis_id,
        peer_name = peer, peer_identity_pk = unname(context$pinset[[peer]]),
        semantic_key = semantic, artifact_sha256 = .dsvert_dp_capsule_source_hash(artifact),
        source_contract_sha256 = source$contract_hash,
        profile_sha256 = artifact$numeric_certificate$profile_sha256,
        certificate_sha256 = artifact$numeric_certificate$certificate_sha256,
        stage_plan_digest = strrep("c", 64), private_result_exposed = FALSE)
      if (phase == "prepared") value <- c(value, list(operation_id = "op_cox",
        source_key = "source", output_key = "output", operation = "cox-loss-staged-v1",
        purpose = paste0("cox-loss-staged-v1/", strrep("d", 64)),
        vector_len = artifact$coordinate_count, persisted = unname(persisted[[peer]])))
      if (phase %in% c("persisted", "complete")) value$stage_receipt <- strrep("e", 64)
      if (phase == "complete") value$coordinate_count <- artifact$coordinate_count
      if (!is.null(change) && identical(phase, change$phase)) value[[change$field]] <- change$value
      sign(value, peer)
    }
    route <- function() {
      value <- receipt("bound", "site_a")
      value$signature <- value$phase <- NULL
      value$version <- "dsvert-cox-staged-routing-receipt-v1"
      value$routing_digest <- strrep("f", 64)
      sign(value, "site_a")
    }
    testthat::local_mocked_bindings(
      .dsvert_dp_cross_exact_setup = function(...) list(),
      .dsvert_dp_cross_exact_cleanup = function(...) { cleanups <<- cleanups + 1L },
      .dsvert_dp_synopsis_runner_json_set = function(...) NULL,
      .dsvert_dp_alignment_mask_run = function(manifest_json, context, layout, ...) {
        projection <- .dsvert_dp_alignment_mask_private_projection_client(layout)
        expect_identical(projection$source_offset, as.numeric(layout$private_start - 1L))
        expect_identical(projection$coordinate_count,
          as.numeric(layout$transport_coordinate_count - layout$private_start + 1L))
        actions <<- c(actions, "alignment")
      },
      .dsvert_exact_gc_run = function(...) {
        args <- list(...)
        expect_identical(args$operation, "cox-loss-staged-v1")
        expect_identical(args$ring, 128L)
        expect_identical(args$frac_bits, 0L)
        runs <<- runs + 1L
      },
      .dsvert_fanout_by_site = function(conns, calls, operation, .aggregate) {
        request <- as.list(calls[[1L]])
        action <- request$action
        if (is.null(action)) return(setNames(rep(list("{}"), length(calls)), names(calls)))
        actions <<- c(actions, if (action == "bind") paste0("bind:", names(calls)) else action)
        if (action == "bind") {
          peer <- names(calls)
          if (peer == "site_b") expect_identical(request$routing_receipt_json,
            .dsvert_joint_dp_client_json(route()))
          return(setNames(list(.dsvert_joint_dp_client_json(list(
            bound = receipt("bound", peer), routing_receipt = route()))), peer))
        }
        phase <- switch(action, prepare = "prepared", store = "persisted", finalize = "complete")
        if (action == "start") return(setNames(rep(list("{}"), length(calls)), names(calls)))
        setNames(lapply(names(calls), function(peer)
          .dsvert_joint_dp_client_json(receipt(phase, peer))), names(calls))
      }, .package = "dsVertClient")
    run <- function(schema = f$schema_manifest) .dsvert_dp_glm_grid_cross_orchestrate(
      .dsvert_joint_dp_client_json(manifest), manifest, context, source,
      function(...) stop("unexpected"),
      list(manifest_sha256 = strrep("1", 64), claim_set_json = "{}", compilation_json = "{}"),
      .schema_json = .dsvert_joint_dp_client_json(schema))
    bad_schema <- f$schema_manifest
    bad_schema$datasets$aligned$columns$time$upper <- 21
    expect_error(run(bad_schema))
    expect_length(actions, 0L)
    for (state in list(c(FALSE, FALSE), c(TRUE, FALSE), c(FALSE, TRUE), c(TRUE, TRUE))) {
      persisted[] <- state; actions <- character(); runs <- 0L
      before <- cleanups
      result <- run()
      expect_true(result$sampler_handoff_ready)
      expect_false(result$private_result_exposed)
      expect_identical(result$receipts$cox_grid$stage_receipt, strrep("e", 64))
      expect_identical(actions, c("alignment", "bind:site_a", "bind:site_b", "prepare",
        if (!all(state)) c("start", "store"), "finalize"))
      expect_identical(runs, if (all(state)) 0L else 1L)
      expect_identical(cleanups, before + 1L)
    }
    persisted[] <- FALSE
    for (bad in list(
        list(phase = "prepared", field = "operation", value = "grouped-cox-staged-v1"),
        list(phase = "prepared", field = "purpose", value = paste0("grouped-cox-staged-v1/", strrep("d", 64))),
        list(phase = "prepared", field = "stage_plan_digest", value = strrep("2", 64)),
        list(phase = "complete", field = "stage_receipt", value = strrep("2", 64)))) {
      change <- bad; before <- cleanups; actions <- character()
      expect_error(run())
      expect_identical(cleanups, before + 1L)
      if (bad$phase == "prepared") expect_false("start" %in% actions)
    }
    change <- NULL
    original <- source
    for (field in c("purpose", "private_layout_sha256", "coordinate_count")) {
      source[[field]] <- if (field == "coordinate_count") -1 else "wrong"
      actions <- character()
      expect_error(run()); expect_length(actions, 0)
      source <- original
    }
    manifest$workload$families$gaussian_models$artifacts[[1L]]$statistic_maximum <- -1
    actions <- character(); expect_error(run()); expect_length(actions, 0)
  }
})
