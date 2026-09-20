test_that("Cox signed drafts preserve contracts and reject a different advertising owner", {
  for (owners in c(2L, 3L, 5L)) {
    f <- .cox_cross_client_fixture(capacity = 5, owners = owners)
    fragments <- list(describe = list(), survival = list(), vertical_cross = list(),
      gaussian = list(cox_grid = list(version = "cox_grid_cross_v1", dataset = "aligned",
        contract = .dsvert_joint_dp_client_json(f$contract))))
    context <- list(pinset = f$policy$peer_pinset, status = list())
    for (peer in names(context$pinset)) context$status[[peer]] <- list(policy = list(
      domain = "cox-study", cohort_id = "group",
      peer_pinset_sha256 = f$policy$peer_pinset_sha256))
    draft <- function(peer) list(
      version = .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_DRAFT_VERSION,
      phase = "custodian_policy_draft", peer_name = peer,
      peer_identity_pk = unname(context$pinset[[peer]]),
      peer_pinset_sha256 = f$policy$peer_pinset_sha256,
      domain = "cox-study", cohort_id = "group",
      dataset_mapping_mode = "custodian_explicit_dataset_mapping_v1",
      datasets = list(aligned = list(dataset_id = "aligned", dataset_version = "v1",
        schema_version = .DSVERT_CLIENT_DP_CAPSULE_POLICY_SCHEMA_VERSION,
        alignment_group = "group", alignment_protocol_version = 1,
        patient_column = "patient", columns = Filter(function(column) {
          identical(column$owner_peer, peer)
        }, f$schema$unsigned$datasets$aligned$columns))),
      workload_fragments = fragments, data_access = FALSE,
      patient_derived_metadata = FALSE, operation_limit = FALSE,
      request_limit = FALSE, history_can_deny_operation = FALSE)
    wire <- function(value, peer) {
      value$signature <- f$b64(openssl::ed25519_sign(
        .dsvert_dp_capsule_manifest_message("draft", value), f$keys[[peer]]))
      .dsvert_joint_dp_client_json(value)
    }
    value <- draft("site_a")
    parsed <- .dsvert_dp_capsule_manifest_draft(wire(value, "site_a"), "site_a", context)
    expect_identical(.dsvert_joint_dp_client_json(parsed$workload_fragments$gaussian),
      .dsvert_joint_dp_client_json(fragments$gaussian))
    bootstrap <- .dsvert_dp_synopsis_bootstrap_draft_v1(value, "site_a", context)
    expect_identical(.dsvert_joint_dp_client_json(bootstrap$workload_fragments$gaussian),
      .dsvert_joint_dp_client_json(fragments$gaussian))
    value <- draft("site_b")
    expect_error(.dsvert_dp_capsule_manifest_draft(wire(value, "site_b"), "site_b", context))
    expect_error(.dsvert_dp_synopsis_bootstrap_draft_v1(value, "site_b", context))
    value <- draft("site_a")
    signed <- .dsvert_joint_dp_client_decode(wire(value, "site_a"), "draft", 2^20)
    signed$workload_fragments$gaussian$cox_grid$dataset <- "other"
    expect_error(.dsvert_dp_capsule_manifest_draft(
      .dsvert_joint_dp_client_json(signed), "site_a", context), "signature")
  }
})
