# These tests exercise the two production schema normalizers. The server policy
# context seam supplies only public pins; schema/workload normalization and all
# hashes are computed by the unchanged server implementation.
.cross_gaussian_schema_fixture <- function() {
  helpers <- new.env(parent = environment())
  for (node in parse(testthat::test_path("test-dp-capsule-manifest-client.R"))) {
    if (is.call(node) && identical(node[[1L]], as.name("<-"))) eval(node, helpers)
  }
  fixture <- helpers$.capsule_manifest_client_fixture()
  fixture$drafts <- lapply(fixture$drafts, function(raw) {
    value <- .dsvert_joint_dp_client_decode(raw, "signed test draft", 8L * 1024L^2)
    value$signature <- NULL
    value$workload_fragments <- list(
      describe = list(), survival = list(), gaussian = list(), vertical_cross = list())
    value
  })
  fixture$drafts$peer_a$workload_fragments$gaussian$correlation <- list(
    version = "v2", dataset = "data_peer_a", outcome = "x_peer_a",
    predictors = "time_peer_b", intercept = TRUE)
  fixture$policy <- c(fixture$context$status$peer_a$policy, list(
    peer_name = "peer_a", peer_pinset = fixture$pins,
    datasets = list(local = list(alignment_manifest_version = 1))))
  fixture
}

.cross_gaussian_schema_server <- function(fixture) {
  testthat::skip_if_not_installed("dsVert")
  testthat::with_mocked_bindings(
    dsVert:::.dsvert_dp_synopsis_manifest_schema_v1(
      lapply(fixture$drafts, function(draft) list(draft = draft)), fixture$policy),
    .dsvert_dp_synopsis_policy_context_v1 = function(...) list(pins = fixture$pins),
    .package = "dsVert")
}

.cross_gaussian_schema_parity <- function(fixture) {
  client <- .dsvert_dp_capsule_manifest_schema(fixture$drafts, fixture$context)
  server <- .cross_gaussian_schema_server(fixture)
  expect_identical(client$json, server$schema_json)
  expect_identical(client$sha256, server$schema_sha256)
  expect_identical(client$workload_contract$json, server$workload_json)
  expect_identical(client$workload_contract$sha256, server$workload_sha256)
  expect_identical(client$logical_snapshot, server$logical_snapshot)
  client
}

test_that("cross Gaussian schema resolves qualified unique references like the server", {
  fixture <- .cross_gaussian_schema_fixture()
  canonical <- .cross_gaussian_schema_parity(fixture)
  original <- .dsvert_dp_capsule_manifest_workload_contract(fixture$drafts, fixture$context)
  expect_identical(canonical$workload_contract, original)
  spec <- fixture$drafts$peer_a$workload_fragments$gaussian$correlation
  spec$outcome <- "peer_a$x_peer_a"
  spec$predictors <- "peer_b$time_peer_b"
  fixture$drafts$peer_a$workload_fragments$gaussian$correlation <- spec
  qualified <- .cross_gaussian_schema_parity(fixture)
  expect_identical(qualified$json, canonical$json)
  expect_identical(qualified$workload_contract, canonical$workload_contract)
})

test_that("cross Gaussian schema preserves valid duplicate owner-qualified keys", {
  for (raw_names in list(c("x", "x"), c("peer_a$x", "peer_b$x"),
                         c("x", "peer_b$x"))) {
    fixture <- .cross_gaussian_schema_fixture()
    a <- fixture$drafts$peer_a$datasets$data_peer_a$columns$x_peer_a
    b <- fixture$drafts$peer_b$datasets$data_peer_b$columns$time_peer_b
    fixture$drafts$peer_a$datasets$data_peer_a$columns <- setNames(list(a), raw_names[[1L]])
    fixture$drafts$peer_b$datasets$data_peer_b$columns <- setNames(list(b), raw_names[[2L]])
    spec <- fixture$drafts$peer_a$workload_fragments$gaussian$correlation
    spec$outcome <- "peer_a$x"
    spec$predictors <- "peer_b$x"
    fixture$drafts$peer_a$workload_fragments$gaussian$correlation <- spec
    result <- .cross_gaussian_schema_parity(fixture)
    expect_identical(names(result$unsigned$datasets$data_peer_a$columns), "peer_a$x")
    expect_identical(names(result$unsigned$datasets$data_peer_b$columns), "peer_b$x")
    expect_identical(result$workload_contract,
      .dsvert_dp_capsule_manifest_workload_contract(fixture$drafts, fixture$context))
  }
})

test_that("cross Gaussian schema rejects ambiguous physical columns across datasets", {
  fixture <- .cross_gaussian_schema_fixture()
  fixture$drafts$peer_b$datasets$other <- fixture$drafts$peer_b$datasets$data_peer_b
  expect_error(.dsvert_dp_capsule_manifest_schema(fixture$drafts, fixture$context),
               "ambiguous across datasets")
  expect_error(.cross_gaussian_schema_server(fixture), "invalid or ambiguous")

  fixture <- .cross_gaussian_schema_fixture()
  fixture$drafts$peer_b$datasets$data_peer_b$columns$x_peer_a <-
    fixture$drafts$peer_b$datasets$data_peer_b$columns$time_peer_b
  fixture$drafts$peer_a$workload_fragments$gaussian$correlation$predictors <- "x_peer_a"
  expect_error(.dsvert_dp_capsule_manifest_schema(fixture$drafts, fixture$context),
               "missing or ambiguous")
  expect_error(.cross_gaussian_schema_server(fixture), "missing or ambiguous")
})

test_that("same-owner schema and immutable grid contracts retain their existing bytes", {
  fixture <- .cross_gaussian_schema_fixture()
  fixture$drafts$peer_a$workload_fragments$gaussian$correlation <- list(
    version = "v1", dataset = "data_peer_a", outcome = "x_peer_a",
    predictors = character(), intercept = TRUE)
  result <- .cross_gaussian_schema_parity(fixture)
  original <- .dsvert_dp_capsule_manifest_workload_contract(fixture$drafts, fixture$context)
  expect_identical(result$workload_contract, original)
  # This helper's scope is v2: signed extension contracts stay opaque and exact.
  original$value$gaussian$grid <- list(owner_peer = "peer_b", spec = list(
    version = "negative_binomial_grid_cross_v1",
    contract = list(predictors = "peer_a$x_peer_a", signature = "immutable")))
  original$value <- .dsvert_joint_dp_client_canonical(original$value)
  original$json <- .dsvert_joint_dp_client_json(original$value)
  original$sha256 <- digest::digest(original$json, algo = "sha256", serialize = FALSE)
  expect_identical(.dsvert_dp_capsule_manifest_cross_gaussian_references(
    original, result$unsigned$datasets), original)
})

test_that("server still rejects noncanonical raw local Gaussian registration", {
  fixture <- .cross_gaussian_schema_fixture()
  local <- fixture$drafts$peer_a$workload_fragments
  canonical <- .cross_gaussian_schema_parity(fixture)
  local$gaussian$correlation$outcome <- "peer_a$x_peer_a"
  local$gaussian$correlation$predictors <- "peer_b$time_peer_b"
  testthat::with_mocked_bindings(
    expect_error(dsVert:::.dsvert_dp_capsule_manifest_workload_contract(
      canonical$workload_contract$json, fixture$policy),
      "changed a custodian-owned fragment"),
    .dsvert_dp_capsule_manifest_local_specs = function(...) local,
    .package = "dsVert")
})
