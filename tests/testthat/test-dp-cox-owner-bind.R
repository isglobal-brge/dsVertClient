test_that("Cox binds the owner first and authenticates the public routing relay", {
  skip_if_not_installed("dsVert")
  for (owners in c(2L, 3L, 5L)) {
    f <- .cox_cross_client_fixture(capacity = 4, owners = owners)
    artifact <- .dsvert_dp_cox_cross_workload_artifact(f$contract)
    manifest <- list(workload = list(capsule_mechanism = list(mechanism = "discrete-laplace")))
    source <- list(capsule_id = strrep("a", 64))
    source_receipt <- list(capsule_id = source$capsule_id,
      contract_hash = .dsvert_dp_capsule_source_hash(source))
    semantic <- get(".dsvert_dp_glm_grid_cross_key", asNamespace("dsVert"))(
      manifest, artifact, source)
    context <- list(designated = c("site_b", "site_a"), pinset = f$policy$peer_pinset,
      conns = list(site_a = "a", site_b = "b"))
    route <- list(version = "dsvert-cox-staged-routing-receipt-v1",
      capsule_id = source$capsule_id, analysis_id = artifact$analysis_id,
      peer_name = "site_a", peer_identity_pk = unname(context$pinset[["site_a"]]),
      semantic_key = semantic, artifact_sha256 = .dsvert_dp_capsule_source_hash(artifact),
      source_contract_sha256 = source_receipt$contract_hash,
      profile_sha256 = artifact$numeric_certificate$profile_sha256,
      certificate_sha256 = artifact$numeric_certificate$certificate_sha256,
      private_result_exposed = FALSE, routing_digest = strrep("b", 64),
      stage_plan_digest = strrep("c", 64))
    sign <- function(value, peer) {
      value$signature <- NULL
      value$signature <- f$b64(openssl::ed25519_sign(charToRaw(paste0(
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SIGNATURE_DOMAIN, "cross-grid-result|",
        .dsvert_joint_dp_client_json(value))), f$keys[[peer]]))
      value
    }
    bounds <- setNames(lapply(context$designated, function(peer) {
      value <- route[setdiff(names(route), "routing_digest")]
      value$version <- "dsvert-cox-staged-receipt-v1"
      value$phase <- "bound"
      value$peer_name <- peer
      value$peer_identity_pk <- unname(context$pinset[[peer]])
      sign(value, peer)
    }), context$designated)
    calls <- list()
    supplied_route <- sign(route, "site_a")
    second_route <- NULL
    testthat::local_mocked_bindings(.dsvert_fanout_by_site = function(
        conns, requests, operation, .aggregate) {
      peer <- names(requests)
      calls[[length(calls) + 1L]] <<- list(peer = peer, request = as.list(requests[[1]]))
      setNames(list(.dsvert_joint_dp_client_json(list(bound = bounds[[peer]],
        routing_receipt = if (peer == "site_b" && !is.null(second_route))
          second_route else supplied_route))), peer)
    }, .package = "dsVertClient")
    bind <- function() .dsvert_dp_cox_cross_bind(manifest, context, artifact,
      source_receipt, "cox-test", list(manifest_sha256 = strrep("d", 64),
        claim_set_json = "{}", compilation_json = "{}"), function(...) stop("unexpected"))
    result <- bind()
    expect_identical(vapply(calls, `[[`, "", "peer"), c("site_a", "site_b"))
    expect_identical(calls[[1]]$request$routing_receipt_json, "")
    expect_identical(calls[[2]]$request$routing_receipt_json,
      .dsvert_joint_dp_client_json(supplied_route))
    expect_identical(result$semantic_key, semantic)
    expect_false(result$private_result_exposed)
    for (field in c("capsule_id", "semantic_key", "source_contract_sha256",
        "artifact_sha256", "profile_sha256", "certificate_sha256")) {
      changed <- route; changed[[field]] <- strrep("e", 64)
      supplied_route <- sign(changed, "site_a"); calls <- list()
      expect_error(bind())
      expect_length(calls, 1)
    }
    for (field in c("routing_digest", "stage_plan_digest")) {
      changed <- route; changed[[field]] <- strrep("0", 64)
      supplied_route <- sign(changed, "site_a"); calls <- list()
      expect_error(bind())
      expect_length(calls, 1)
    }
    supplied_route <- sign(route, "site_b"); calls <- list()
    expect_error(bind(), "signature")
    expect_length(calls, 1)
    supplied_route <- sign(route, "site_a")
    second_route <- supplied_route; second_route$routing_digest <- strrep("e", 64)
    expect_error(bind())
    second_route <- NULL
    changed <- bounds$site_b; changed$stage_plan_digest <- strrep("e", 64)
    bounds$site_b <- sign(changed, "site_b")
    expect_error(bind())
  }
})
