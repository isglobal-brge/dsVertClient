test_that("the source custodians are exactly the computation and noise authorities", {
  for (family in c("multinomial", "ordinal")) {
    fixture <- .categorical_cross_client_fixture(family)
    policy <- fixture$policy
    policy$peer_pinset <- c(policy$peer_pinset,
      site_c = unname(policy$peer_pinset[1L]))
    policy$designated_noise_peers <- c("site_a", "site_c")
    .categorical_cross_client_reject(function() {
      fixture$registration$spec(fixture$raw, policy, fixture$schema)
    })
    policy$designated_noise_peers <- c("site_a", "site_b")
    schema <- fixture$schema
    schema$unsigned$datasets$aligned$columns$x$owner_peer <- "site_c"
    raw <- fixture$raw
    raw$predictor_order <- c("site_b$z", "site_c$x")
    .categorical_cross_client_reject(function() {
      fixture$registration$spec(raw, policy, schema)
    })
  }
})

test_that("both family validators bind the final numeric certificate", {
  for (family in c("multinomial", "ordinal")) {
    fixture <- .categorical_cross_client_fixture(family)
    changed <- fixture$contract
    changed$spec$numeric_contract$certificate_sha256 <- strrep("0", 64)
    changed <- fixture$sign(changed)
    .categorical_cross_client_reject(function() {
      fixture$registration$validate(changed, fixture$policy,
                                     fixture$schema_manifest)
    })
  }
})
