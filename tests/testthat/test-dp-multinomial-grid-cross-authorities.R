test_that("categorical compute authorities must belong to the signed source owners", {
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
    schema$unsigned$datasets$aligned$columns$y$owner_peer <- "site_c"
    raw <- fixture$raw
    raw$predictor_order <- c("site_b$z", "site_c$x")
    raw$outcome <- "site_c$y"
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


test_that("unpromoted categorical frontdoors stay outside the public inventory", {
  ns <- asNamespace("dsVertClient")
  for (family in c("multinomial", "ordinal")) {
    frontdoor <- paste0("dp_", family, "_grid")
    expect_false(frontdoor %in% getNamespaceExports("dsVertClient"))
    expect_true(is.function(get(frontdoor, envir = ns, inherits = FALSE)))
    registration <- get(paste0(".dsvert_dp_", family, "_grid_cross_register"), ns)()
    expect_false(registration$release_enabled)
  }
})


test_that("categorical admission authenticates K2/K3/K5 with two compute authorities", {
  for (family in c("multinomial", "ordinal")) for (k in c(2L, 3L, 5L)) {
    fixture <- .categorical_cross_client_fixture(family, peer_count = k)
    admitted <- .dsvert_dp_glm_grid_profile_admit(fixture$contract,
      fixture$policy, fixture$schema_manifest)
    expect_true(.dsvert_dp_glm_grid_cross_equal(admitted, fixture$contract))
    expect_length(admitted$spec$participating_peers, k)
    expect_length(admitted$spec$computation_peers, 2L)
    artifact <- .dsvert_dp_glm_grid_cross_workload_artifact(admitted)
    checked <- .dsvert_dp_glm_grid_cross_client_artifact(artifact,
      admitted$spec$dataset, admitted$spec$analysis_id, NULL,
      admitted$spec$adjacency, 2^admitted$spec$numeric_grid_bits,
      admitted$spec$observation_capacity, family)
    expect_equal(checked$class_order, admitted$spec$class_order)
    changed <- fixture$contract
    changed$signatures[[k]] <- NULL
    expect_error(.dsvert_dp_glm_grid_profile_admit(changed,
      fixture$policy, fixture$schema_manifest), class = "dsvert_dp_public_failure")
  }
})
