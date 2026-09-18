test_that("grouped prototypes are inventoried without advertising promotion", {
  methods <- c("dp_lmm_grid", "dp_glmm_grid", "dp_gee_grid")
  status <- ds.vertMethodStatus(methods)
  expect_setequal(status$method, methods)
  expect_true(all(status$status == "quarantine"))
  expect_true(all(status$release_contract == "disclosure_safe_protocol_no_statistic"))
  expect_false(any(status$currently_numerically_certified))
  inventory <- .dsvert_capsule_method_inventory()
  rows <- inventory[inventory$method %in% methods, , drop = FALSE]
  expect_equal(nrow(rows), 3L)
  expect_true(all(rows$current_route_status == "signed_workload_unavailable_quarantine"))
  expect_true(all(rows$artifact_implementation_state == "secure_artifact_not_implemented"))
  expect_true(all(vapply(rows$legacy_remote_call_evidence, nrow, integer(1L)) == 0L))
})
