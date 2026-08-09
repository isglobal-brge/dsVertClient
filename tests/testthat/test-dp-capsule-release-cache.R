.capsule_release_cache_status <- function(
    domain = strrep("1", 64L), namespace = strrep("a", 64L),
    capsules = 0, releases = 0) {
  structure(list(peer_a = list(
    version = "dsvert-joint-dp-capsule-status-v6",
    enabled = TRUE,
    privacy_contract = list(
      scope = paste0(
        "at_most_N_immutable_snapshot_workload_capsules_per_stable_",
        "privacy_accountant_namespace"),
      privacy_accountant_namespace_id = paste0("jdpc1_", namespace),
      privacy_accountant_namespace_enforcement =
        "identity_bound_immutable_receipt_v1"),
    policy = list(peer_pinset_sha256 = strrep("2", 64L)),
    noise_root = list(epoch = 1, key_id = "root-1"),
    release_domain = list(
      generation = 1, domain_id = paste0("rd_", domain)),
    role = list(designated_noise_peer = TRUE),
    composition_telemetry = list(capsules_created = capsules),
    release_instance_telemetry = list(releases_published = releases))),
    class = c("ds.vertDPStatus", "ds.vertJointDPCapsuleStatus", "list"))
}

.capsule_release_cache_value <- function() {
  release <- structure(list(
    values = 1, coordinate_count = 1L,
    coordinate_order_sha256 = strrep("3", 64L),
    manifest = list(workload = list())),
    class = c("dsvert_joint_dp_vector", "list"))
  list(
    release = release,
    layout = list(
      coordinate_count = 1L, sha256 = strrep("3", 64L)),
    status = .capsule_release_cache_status(),
    manifest_bundle = list(
      schema_json = "{}", manifest_sha256 = strrep("4", 64L),
      capsule_id = strrep("7", 64L),
      context = new.env(parent = emptyenv())))
}

test_that("public release cache is bounded, LRU and connection-free", {
  .dsvert_dp_capsule_release_cache_clear()
  on.exit(.dsvert_dp_capsule_release_cache_clear(), add = TRUE)
  datasources <- list(peer_a = list())
  status <- .capsule_release_cache_status()
  bundle <- .capsule_release_cache_value()$manifest_bundle
  key1 <- .dsvert_dp_capsule_release_cache_key(datasources, status, bundle)
  key2 <- .dsvert_dp_capsule_release_cache_key(
    datasources, .capsule_release_cache_status(strrep("5", 64L)), bundle)
  key3 <- .dsvert_dp_capsule_release_cache_key(
    datasources, .capsule_release_cache_status(strrep("6", 64L)), bundle)
  expect_match(key1, "^[0-9a-f]{64}$")
  expect_false(identical(key1, key2))
  expect_false(identical(
    key1,
    .dsvert_dp_capsule_release_cache_key(
      datasources,
      .capsule_release_cache_status(namespace = strrep("b", 64L)),
      bundle)))
  expect_identical(
    key1,
    .dsvert_dp_capsule_release_cache_key(
      datasources, .capsule_release_cache_status(capsules = 9, releases = 7),
      bundle))
  changed_bundle <- bundle
  changed_bundle$manifest_sha256 <- strrep("8", 64L)
  expect_false(identical(
    key1,
    .dsvert_dp_capsule_release_cache_key(
      datasources, status, changed_bundle)))
  two_peer_status <- unclass(status)
  two_peer_status$peer_b <- two_peer_status$peer_a
  class(two_peer_status) <- class(status)
  expect_identical(
    .dsvert_dp_capsule_release_cache_key(
      list(peer_a = list(), peer_b = list()), two_peer_status, bundle),
    .dsvert_dp_capsule_release_cache_key(
      list(peer_b = list(), peer_a = list()), two_peer_status, bundle))
  malformed_status <- status
  malformed_status[[1L]]$role <- NULL
  expect_null(.dsvert_dp_capsule_release_cache_key(
    datasources, malformed_status, bundle))

  value <- .capsule_release_cache_value()
  expect_true(.dsvert_dp_capsule_release_cache_put(
    key1, value, .max_entries = 2L, .max_bytes = 1024^2))
  cached <- .dsvert_dp_capsule_release_cache_get(key1)
  expect_null(cached$manifest_bundle$context)
  expect_identical(cached$release, value$release)
  expect_true(.dsvert_dp_capsule_release_cache_put(
    key2, value, .max_entries = 2L, .max_bytes = 1024^2))
  expect_identical(.dsvert_dp_capsule_release_cache_get(key1)$release,
                   value$release)
  expect_true(.dsvert_dp_capsule_release_cache_put(
    key3, value, .max_entries = 2L, .max_bytes = 1024^2))
  expect_null(.dsvert_dp_capsule_release_cache_get(key2))
  expect_identical(.dsvert_dp_capsule_release_cache_get(key1)$release,
                   value$release)
  expect_identical(.dsvert_dp_capsule_release_cache_get(key3)$release,
                   value$release)

  .dsvert_dp_capsule_release_cache_clear()
  expect_false(.dsvert_dp_capsule_release_cache_put(
    key1, value, .max_entries = 2L, .max_bytes = 1L))
  expect_null(.dsvert_dp_capsule_release_cache_get(key1))
})

test_that("same authenticated control state reuses the validated vector", {
  .dsvert_dp_capsule_release_cache_clear()
  on.exit(.dsvert_dp_capsule_release_cache_clear(), add = TRUE)
  datasources <- list(peer_a = list())
  statuses <- list(
    .capsule_release_cache_status(capsules = 0, releases = 0),
    .capsule_release_cache_status(capsules = 1, releases = 1))
  calls <- new.env(parent = emptyenv())
  calls$status <- calls$manifest <- calls$vector <- 0L
  manifest_bundle <- .capsule_release_cache_value()$manifest_bundle
  release <- .capsule_release_cache_value()$release
  layout <- .capsule_release_cache_value()$layout

  run <- function() testthat::with_mocked_bindings(
    .dsvert_dp_capsule_vector_run(datasources),
    .dsvert_joint_dp_capsule_status_impl = function(...) {
      calls$status <- calls$status + 1L
      statuses[[calls$status]]
    },
    .dsvert_dp_capsule_manifest_build = function(...) {
      calls$manifest <- calls$manifest + 1L
      manifest_bundle
    },
    .dsvert_joint_dp_vector_capsule = function(...) {
      calls$vector <- calls$vector + 1L
      release
    },
    .dsvert_dp_capsule_vector_layout = function(...) layout,
    .package = "dsVertClient")

  first <- run()
  second <- run()
  expect_identical(first$release, second$release)
  expect_identical(first$layout, second$layout)
  expect_identical(second$manifest_bundle, manifest_bundle)
  expect_identical(calls$status, 2L)
  expect_identical(calls$manifest, 2L)
  expect_identical(calls$vector, 1L)
})

test_that("a changed authoritative manifest cannot reuse a stale vector", {
  .dsvert_dp_capsule_release_cache_clear()
  on.exit(.dsvert_dp_capsule_release_cache_clear(), add = TRUE)
  datasources <- list(peer_a = list())
  manifests <- list(
    .capsule_release_cache_value()$manifest_bundle,
    .capsule_release_cache_value()$manifest_bundle)
  manifests[[2L]]$manifest_sha256 <- strrep("8", 64L)
  manifests[[2L]]$capsule_id <- strrep("9", 64L)
  calls <- new.env(parent = emptyenv())
  calls$manifest <- calls$vector <- 0L
  release <- .capsule_release_cache_value()$release
  layout <- .capsule_release_cache_value()$layout

  run <- function() testthat::with_mocked_bindings(
    .dsvert_dp_capsule_vector_run(datasources),
    .dsvert_joint_dp_capsule_status_impl = function(...) {
      .capsule_release_cache_status()
    },
    .dsvert_dp_capsule_manifest_build = function(...) {
      calls$manifest <- calls$manifest + 1L
      manifests[[calls$manifest]]
    },
    .dsvert_joint_dp_vector_capsule = function(...) {
      calls$vector <- calls$vector + 1L
      release
    },
    .dsvert_dp_capsule_vector_layout = function(...) layout,
    .package = "dsVertClient")

  run()
  run()
  expect_identical(calls$manifest, 2L)
  expect_identical(calls$vector, 2L)
  expect_length(.dsvert_dp_capsule_release_cache$entries, 2L)
})

test_that("failed vector attempts are never cached", {
  .dsvert_dp_capsule_release_cache_clear()
  on.exit(.dsvert_dp_capsule_release_cache_clear(), add = TRUE)
  datasources <- list(peer_a = list())
  attempts <- 0L
  run <- function() testthat::with_mocked_bindings(
    .dsvert_dp_capsule_vector_run(datasources),
    .dsvert_joint_dp_capsule_status_impl = function(...) {
      .capsule_release_cache_status()
    },
    .dsvert_dp_capsule_manifest_build = function(...) {
      .capsule_release_cache_value()$manifest_bundle
    },
    .dsvert_joint_dp_vector_capsule = function(...) {
      attempts <<- attempts + 1L
      stop("test vector failure", call. = FALSE)
    },
    .package = "dsVertClient")
  expect_error(run(), "test vector failure")
  expect_error(run(), "test vector failure")
  expect_identical(attempts, 2L)
  expect_length(.dsvert_dp_capsule_release_cache$entries, 0L)
})
