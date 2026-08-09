.client_exact_gc_vector_binding <- function() {
  hex <- function(label) digest::digest(
    label, algo = "sha256", serialize = FALSE)
  identity <- list(
    domain = "dsVert/joint-dp/vector/exact-gc-operation/v1",
    manifest_sha256 = hex("manifest"),
    release_contract_hash = hex("release"),
    selection_sha256 = hex("selection"),
    transcript_hash = hex("transcript"),
    chunk_index = 0L, coordinate_count = 2L,
    circuit_digest = hex("circuit"))
  identity$purpose <- paste0(
    "joint-dp-vector-laplace-v3/", identity$circuit_digest)
  suffix <- substr(
    .dsvert_joint_dp_vector_exact_gc_client_hash(identity), 1L, 32L)
  unsigned <- c(list(
    version = "dsvert-joint-dp-vector-exact-gc-binding-v1"),
    identity, list(
      operation_id = paste0("op_", suffix),
      source_key = paste0("exact_gc_in_", suffix),
      output_key = paste0("exact_gc_out_", suffix),
      operation = "joint-dp-vector-laplace-v3",
      output_kind = "joint-dp-vector-ring128-share-v1",
      source_producer = "joint.dp.vector.one-draw.v1"))
  binding <- c(unsigned, list(
    binding_sha256 =
      .dsvert_joint_dp_vector_exact_gc_client_hash(unsigned)))
  list(hex = hex, identity = identity, binding = binding)
}

.client_exact_gc_vector_initializations <- function(binding) {
  peers <- c(
    site_a = paste0("dsv1_", strrep("a", 64L)),
    site_b = paste0("dsv1_", strrep("b", 64L)))
  roles <- c(site_a = "garbler", site_b = "evaluator")
  stats::setNames(lapply(names(peers), function(server) list(
    capability_id = "exact_gc_v1",
    peer_id = peers[[server]],
    peer_peer_id = peers[[setdiff(names(peers), server)]],
    role = roles[[server]], context_hash = strrep("c", 64L),
    operation = "joint-dp-vector-laplace-v3",
    output_kind = "joint-dp-vector-ring128-share-v1",
    purpose = binding$purpose,
    source_producer = "joint.dp.vector.one-draw.v1",
    ring_bits = 128L, frac_bits = 0L, vector_len = 2L,
    threshold = "", chunk_bytes = 65536L, ttl_seconds = 180L,
    max_runtime_seconds = 21600L, worker_heartbeat = 1,
    state = "running", stored = FALSE)), names(peers))
}

test_that("client validates the deterministic manifest-bound operation", {
  fixture <- .client_exact_gc_vector_binding()
  expect_invisible(.dsvert_joint_dp_vector_exact_gc_client_binding(
    fixture$binding, fixture$identity$manifest_sha256,
    fixture$identity$release_contract_hash,
    fixture$identity$selection_sha256,
    fixture$identity$transcript_hash, 0L, 2L))

  tampered <- fixture$binding
  tampered$operation_id <- paste0("op_", strrep("f", 32L))
  tampered$source_key <- paste0("exact_gc_in_", strrep("f", 32L))
  tampered$output_key <- paste0("exact_gc_out_", strrep("f", 32L))
  tampered$binding_sha256 <-
    .dsvert_joint_dp_vector_exact_gc_client_hash(
      tampered[setdiff(names(tampered), "binding_sha256")])
  expect_error(.dsvert_joint_dp_vector_exact_gc_client_binding(
    tampered, fixture$identity$manifest_sha256,
    fixture$identity$release_contract_hash,
    fixture$identity$selection_sha256,
    fixture$identity$transcript_hash, 0L, 2L), "invalid")
})

test_that("generic exact-GC validation accepts only vector share metadata", {
  fixture <- .client_exact_gc_vector_binding()
  initialized <- .client_exact_gc_vector_initializations(fixture$binding)
  validated <- .dsvert_exact_gc_validate_init(
    initialized, names(initialized), "joint-dp-vector-laplace-v3",
    128L, 0L, 2L, fixture$binding$purpose)
  expect_identical(validated$peer_ids,
                   vapply(initialized, `[[`, character(1L), "peer_id"))

  initialized$site_b$output_kind <- "ring-share"
  expect_error(.dsvert_exact_gc_validate_init(
    initialized, names(initialized), "joint-dp-vector-laplace-v3",
    128L, 0L, 2L, fixture$binding$purpose), "rejected|changed")
})

test_that("one-draw pump uses two peers and returns no protected payload", {
  fixture <- .client_exact_gc_vector_binding()
  initialized <- .client_exact_gc_vector_initializations(fixture$binding)
  starts <- lapply(initialized, function(value) list(
    backend = .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BACKEND,
    binding_sha256 = fixture$binding$binding_sha256,
    operation_id = fixture$binding$operation_id,
    purpose = fixture$binding$purpose,
    initialization = value,
    intermediate_payload_exposed = FALSE,
    source_share_exposed = FALSE,
    private_seed_exposed = FALSE,
    preclamp_values_exposed = FALSE))
  observed <- NULL
  runner <- function(...) {
    observed <<- list(...)
    invisible(list(context_hash = strrep("c", 64L)))
  }
  datasources <- structure(list(list(), list()),
                           names = c("site_a", "site_b"))
  result <- .dsvert_joint_dp_vector_exact_gc_run(
    datasources, names(datasources), 1:2,
    "00000000-0000-4000-8000-000000000001",
    fixture$binding, fixture$identity$manifest_sha256,
    fixture$identity$release_contract_hash,
    fixture$identity$selection_sha256,
    fixture$identity$transcript_hash, 0L, 2L, starts,
    .run = runner, .aggregate = function(...) NULL)
  expect_identical(observed$operation,
                   "joint-dp-vector-laplace-v3")
  expect_identical(observed$ring, 128L)
  expect_identical(observed$vector_len, 2L)
  expect_identical(names(observed$initialized), c("site_a", "site_b"))
  expect_true(result$complete)
  expect_false(result$intermediate_payload_exposed)
  expect_false(any(c("share", "source_share", "private_seed",
                     "validity_share") %in% names(result)))
})

test_that("timeout retry preserves operation and never selects fallback", {
  fixture <- .client_exact_gc_vector_binding()
  initialized <- .client_exact_gc_vector_initializations(fixture$binding)
  datasources <- structure(list(list(), list()),
                           names = c("site_a", "site_b"))
  calls <- list()
  runner <- function(...) {
    calls[[length(calls) + 1L]] <<- list(...)
    if (length(calls) == 1L) stop("transport timeout", call. = FALSE)
    invisible(TRUE)
  }
  invoke <- function() .dsvert_joint_dp_vector_exact_gc_run(
    datasources, names(datasources), 1:2,
    "00000000-0000-4000-8000-000000000001",
    fixture$binding, fixture$identity$manifest_sha256,
    fixture$identity$release_contract_hash,
    fixture$identity$selection_sha256,
    fixture$identity$transcript_hash, 0L, 2L, initialized,
    .run = runner, .aggregate = function(...) NULL)
  expect_error(invoke(), "timeout")
  expect_true(invoke()$complete)
  expect_identical(calls[[1L]]$operation_id, calls[[2L]]$operation_id)
  expect_identical(calls[[1L]]$purpose, calls[[2L]]$purpose)
  expect_false(any(vapply(calls, function(value) {
    any(grepl("fallback|convolution", unlist(value), fixed = FALSE))
  }, logical(1L))))
})

test_that("unsafe start wrappers are rejected before transport", {
  fixture <- .client_exact_gc_vector_binding()
  initialized <- .client_exact_gc_vector_initializations(fixture$binding)
  starts <- lapply(initialized, function(value) list(
    backend = .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BACKEND,
    binding_sha256 = fixture$binding$binding_sha256,
    operation_id = fixture$binding$operation_id,
    purpose = fixture$binding$purpose,
    initialization = value,
    intermediate_payload_exposed = FALSE,
    source_share_exposed = FALSE,
    private_seed_exposed = TRUE,
    preclamp_values_exposed = FALSE))
  expect_error(.dsvert_joint_dp_vector_exact_gc_initializations(
    starts, names(starts), fixture$binding), "unsafe")
})
