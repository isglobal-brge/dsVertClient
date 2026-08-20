test_that("Gaussian Synopsis release validates one bounded signed block", {
  artifact <- list(owner_peer = "peer_a", coordinate_count = 2L)
  block <- list(key = "primary", descriptor = list(key = "primary"))
  context <- list(
    lattice = list(output_lattice_scale = 256), layout = list(),
    manifest = list(workload = list(families = list(
      gaussian_models = list(artifacts = list(primary = list(key = "primary")))))),
    release = list())
  verification <- list(
    integrity_valid = TRUE, authenticity = "session_transport_anchored",
    artifact = artifact, coordinates = c(10, 15),
    validated_moment = list(n = 10))
  saw_synopsis <- FALSE
  coordinates <- verification$coordinates
  testthat::local_mocked_bindings(
    .dsvert_dp_datasources = function(value) value,
    .dsvert_dp_synopsis_vector_run = function(...) list(run = TRUE),
    .dsvert_dp_vector_context = function(run, allow_synopsis = FALSE) {
      saw_synopsis <<- isTRUE(allow_synopsis)
      context
    },
    .dsvert_dp_vector_public_metadata = function(context) list(),
    .dsvert_dp_capsule_single_block = function(...) {
      list(start = 1L, end = 1L, length = 1L)
    },
    .dsvert_dp_vector_block_capacity = function(block) 100,
    .dsvert_dp_gaussian_artifact = function(...) artifact,
    .dsvert_dp_capsule_vector_blocks = function(...) list(block),
    .dsvert_dp_capsule_vector_values = function(...) coordinates,
    .dsvert_dp_gaussian_synopsis_certificate_build = function(...) {
      list(certificate_sha256 = strrep("a", 64L))
    },
    ds.validateDPGaussianCertificate = function(...) verification,
    .package = "dsVertClient")

  released <- .dsvert_dp_gaussian_synopsis_release(
    "cohort", "primary", datasources = list(peer_a = list()),
    .aggregate = function(...) stop("unexpected aggregate", call. = FALSE))

  expect_true(saw_synopsis)
  expect_identical(released$artifact, artifact)
  expect_identical(released$coordinates, verification$coordinates)
  expect_identical(released$moment, verification$validated_moment)
  expect_identical(released$capacity, 100)

  coordinates <- c(10, 101)
  expect_error(
    .dsvert_dp_gaussian_synopsis_release(
      "cohort", "primary", datasources = list(peer_a = list()),
      .aggregate = function(...) stop("unexpected aggregate", call. = FALSE)),
    "violates its signed bounds")
})
