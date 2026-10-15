# Exercise certificate context reconstruction with real signed bootstrap and
# manifest evidence. No closure, signature or trusted-bundle validator is mocked.
.gaussian_synopsis_provenance_fixture <- function(k) {
  helpers <- new.env(parent = asNamespace("dsVertClient"))
  for (node in parse(testthat::test_path("test-dp-synopsis-public-vector-adapter.R"))) {
    if (is.call(node) && identical(node[[1L]], as.name("<-"))) eval(node, helpers)
  }
  fixture <- helpers$.synopsis_public_fixture(k)
  public_bundle <- fixture$manifest_bundle
  public_bundle$context <- NULL
  fixture$certificate <- list(signed_evidence = list(
    federation_state_json = .dsvert_dp_gaussian_synopsis_evidence_json(
      unclass(fixture$status), "federation state"),
    manifest_bundle_json = .dsvert_dp_gaussian_synopsis_evidence_json(
      public_bundle, "manifest bundle")))
  fixture
}

test_that("Gaussian certificate reconstruction retains the signed provenance closure", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .gaussian_synopsis_provenance_fixture(k)
    before <- serialize(fixture$certificate, NULL, version = 3L)
    rebuilt <- .dsvert_dp_gaussian_synopsis_trusted(fixture$certificate)
    expect_identical(rebuilt$trusted$context$authorization_dependencies,
                     fixture$manifest_bundle$context$authorization_dependencies)
    expect_length(rebuilt$trusted$context$authorization_dependencies$dependencies, k)
    expect_identical(rebuilt$trusted$context$pinset,
                     fixture$manifest_bundle$context$pinset)
    expect_identical(rebuilt$bundle$manifest_sha256,
                     fixture$manifest_bundle$manifest_sha256)
    expect_identical(serialize(fixture$certificate, NULL, version = 3L), before)
    expect_true(all(vapply(rebuilt$trusted$context$all_conns, inherits,
                          logical(1L), "dsvert_certificate_peer")))
  }
})

test_that("reconstructed Gaussian provenance still requires signed complete epochs", {
  fixture <- .gaussian_synopsis_provenance_fixture(3L)
  rewrite <- function(status) {
    certificate <- fixture$certificate
    certificate$signed_evidence$federation_state_json <-
      .dsvert_dp_gaussian_synopsis_evidence_json(unclass(status), "federation state")
    certificate
  }
  rotated <- fixture$status
  rotated[[1L]]$authorization_dependencies[[1L]]$authorization_epochs[[3L]]$
    authorization_epoch <- strrep("e", 64L)
  expect_error(.dsvert_dp_gaussian_synopsis_trusted(rewrite(rotated)), "signature")
  missing <- fixture$status
  missing[[1L]]$authorization_dependencies[[1L]]$authorization_epochs <-
    missing[[1L]]$authorization_dependencies[[1L]]$authorization_epochs[-3L]
  expect_error(.dsvert_dp_gaussian_synopsis_trusted(rewrite(missing)),
               "authorization epoch")
})
