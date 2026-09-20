test_that("Cox client source lanes match the authenticated server layout", {
  skip_if_not_installed("dsVert")
  server <- asNamespace("dsVert")
  server_blocks <- get(".dsvert_dp_cox_cross_source_blocks", server)
  server_artifact <- get(".dsvert_dp_cox_cross_workload_artifact", server)
  for (owners in c(2L, 3L, 5L)) {
    f <- .cox_cross_client_fixture(capacity = 5, owners = owners)
    contract <- .dsvert_dp_cox_grid_cross_contract_validate(
      f$contract, f$policy, f$schema_manifest)
    artifact <- .dsvert_dp_cox_cross_workload_artifact(contract)
    expect_identical(.dsvert_joint_dp_client_json(artifact),
      .dsvert_joint_dp_client_json(server_artifact(contract)))
    for (cursor in c(1, 17)) {
      actual <- .dsvert_dp_cox_cross_source_blocks(artifact, cursor)
      expect_identical(actual, server_blocks(artifact, cursor))
      blocks <- actual$blocks
      expect_setequal(vapply(blocks, `[[`, "", "owner_peer"), names(f$keys))
      expect_true(all(vapply(blocks, function(b) b$length == 8L, logical(1))))
      expect_identical(vapply(blocks, `[[`, 0L, "start"),
        setNames(as.integer(seq(cursor, by = 8, length.out = length(blocks))), names(blocks)))
      time <- Filter(function(b) identical(b$reference, "site_a$time"), blocks)
      expect_length(time, 1)
      expect_identical(time[[1]]$kind, "validity")
      expect_true(time[[1]]$private_time_validity)
      expect_identical(time[[1]]$fraction_bits, 0L)
      event <- Filter(function(b) identical(b$reference, "site_a$event"), blocks)
      expect_true(all(vapply(event, function(b) b$maximum == 1 && b$fraction_bits == 0L,
        logical(1))))
      predictors <- Filter(function(b) !b$outcome && !b$private_time_validity &&
        b$kind == "value", blocks)
      expect_true(all(vapply(predictors, function(b) b$maximum == 2^50 &&
        b$fraction_bits == 50L, logical(1))))
      expect_equal(actual$cursor, cursor + length(blocks) * 8)
    }
    for (cursor in c(0, -1, 1.5, Inf,
        .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_TRANSPORT_COORDINATES)) {
      expect_error(.dsvert_dp_cox_cross_source_blocks(artifact, cursor))
    }
    # Public discovery is still closed until the complete owner-first runner lands.
    manifest <- list(workload = list(families = list(gaussian_models = list(
      artifacts = list(cox_grid = artifact)))))
    expect_length(.dsvert_dp_glm_grid_cross_artifacts(manifest), 0)
  }
})
