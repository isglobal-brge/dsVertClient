.grid_cross_client_fixture <- function(family = "binomial") {
  peers <- c("site_a", "site_b", "site_c")
  keys <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  b64 <- function(value) sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))))
  pins <- vapply(keys, function(key) {
    b64(tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  policy <- list(
    peer_pinset = pins,
    peer_pinset_sha256 = .dsvert_dp_capsule_source_hash(as.list(pins)),
    designated_noise_peers = peers[1:2], unit_capacity = 20,
    numeric_grid_bits = 18, adjacency = "add_remove_patient")
  maximum <- if (identical(family, "binomial")) 1 else 8
  snapshot <- list(logical_snapshot_id = "cohort", version = "v1",
                   alignment_protocol_version = 1)
  schema_unsigned <- list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = snapshot,
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    datasets = list(aligned = list(
      dataset_id = "aligned", dataset_version = "v1", schema_version = "v1",
      alignment_group = "group",
      patient_keys = list(site_b = "patient", site_c = "patient"),
      columns = list(
        x = list(kind = "numeric", owner_peer = "site_b", lower = -1, upper = 1),
        y = list(kind = "numeric", owner_peer = "site_c", lower = 0,
                 upper = maximum),
        z = list(kind = "numeric", owner_peer = "site_b", lower = 0,
                 upper = 10)))))
  sign_schema <- function(value) {
    value$signatures <- NULL
    message <- charToRaw(paste0(.DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
                               .dsvert_joint_dp_client_json(value)))
    c(value, list(signatures = lapply(keys, function(key) {
      b64(openssl::ed25519_sign(message, key))
    })))
  }
  schema_manifest <- sign_schema(schema_unsigned)
  schema <- .dsvert_dp_glm_grid_cross_schema_validate(
    policy, snapshot, schema_manifest, .dsvert_dp_glm_grid_cross_verify)
  raw <- list(version = paste0(family, "_grid_cross_v1"),
    analysis_id = "cross_grid", dataset = "aligned", outcome = "site_c$y",
    predictor_order = c("site_b$x", "site_b$z"),
    beta_grid = list(c(0, 0, 0), c(0, 1, 0)), max_outcome = maximum,
    alignment = list(version = "existing_prealigned_logical_dataset_v1",
      method = "pinned_psi_ordered_manifest_v1", alignment_group = "group",
      public_alignment_contract_sha256 = .dsvert_dp_capsule_source_hash(list(
        logical_snapshot = schema$unsigned$logical_snapshot,
        alignment_group = "group", method = "pinned_psi_ordered_manifest_v1")),
      public_patient_dependent_hash = FALSE))
  spec <- .dsvert_dp_glm_grid_cross_spec(raw, policy, schema)
  artifact <- .dsvert_dp_glm_grid_cross_artifact(spec)
  source <- .dsvert_dp_glm_grid_cross_source_contract(spec, artifact)
  sign <- function(value) {
    value$signatures <- NULL
    message <- .dsvert_dp_glm_grid_cross_message(value)
    c(value, list(signatures = lapply(keys, function(key) {
      b64(openssl::ed25519_sign(message, key))
    })))
  }
  contract <- sign(list(version = .DSVERT_CLIENT_DP_GLM_GRID_CROSS_CONTRACT_VERSION,
                       spec = spec, artifact = artifact,
                       source_contract = source))
  list(contract = contract, policy = policy, schema_manifest = schema_manifest,
       schema = schema, sign = sign, sign_schema = sign_schema, raw = raw)
}

.grid_cross_client_reject <- function(value, fixture, schema = fixture$schema_manifest) {
  error <- tryCatch(.dsvert_dp_glm_grid_cross_contract_validate(
    value, fixture$policy, schema), error = identity)
  expect_s3_class(error, "dsvert_dp_public_failure")
  expect_identical(conditionMessage(error),
    "[dsvert_dp_public_failure:v1] Protected capsule operation failed.")
}

test_that("new cross grid contracts accept only complete authenticated public plans", {
  for (family in c("binomial", "poisson")) {
    fixture <- .grid_cross_client_fixture(family)
    contract <- .dsvert_dp_glm_grid_cross_contract_validate(
      fixture$contract, fixture$policy, fixture$schema_manifest)
    expect_identical(contract, .dsvert_joint_dp_client_canonical(fixture$contract))
    expect_identical(unlist(contract$spec$participating_peers),
                     c("site_b", "site_c"))
    expect_identical(unlist(contract$spec$computation_peers),
                     c("site_a", "site_b"))
    expect_identical(contract$artifact$implementation_state,
                     "cross_owner_exact_gc_materialized")
    expect_identical(contract$artifact$cross_owner_state,
                     "exact_gc_to_joint_dp_vector_v1")
    layout <- contract$source_contract$private_layout
    expect_identical(layout$blocks[[1L]]$value_fraction_bits, 50)
    expect_identical(layout$blocks[[3L]]$value_fraction_bits, 0)
    expect_identical(layout$blocks[[3L]]$validity_maximum, 1)
    expect_identical(layout$padding_validity, 0)
    wire <- jsonlite::fromJSON(.dsvert_joint_dp_client_json(contract),
                               simplifyVector = FALSE)
    expect_identical(.dsvert_dp_glm_grid_cross_contract_validate(
      wire, fixture$policy, fixture$schema_manifest), contract)
  }
})

test_that("unsigned and changed bounds and incomplete owner signatures fail closed", {
  fixture <- .grid_cross_client_fixture()
  bad <- fixture$contract
  bad$signatures <- NULL
  .grid_cross_client_reject(bad, fixture)
  bad <- fixture$contract
  bad$signatures$site_c <- NULL
  .grid_cross_client_reject(bad, fixture)
  schema <- fixture$schema_manifest
  schema$signatures$site_c <- NULL
  .grid_cross_client_reject(fixture$contract, fixture, schema)
  schema <- fixture$schema_manifest
  schema$datasets$aligned$columns$x$upper <- 999
  .grid_cross_client_reject(fixture$contract, fixture, schema)
  .grid_cross_client_reject(fixture$contract, fixture,
                             fixture$sign_schema(schema))
  for (field in c("lower", "upper")) {
    bad <- fixture$contract
    bad$spec$predictors[["site_b$x"]][[field]] <- 999
    .grid_cross_client_reject(fixture$sign(bad), fixture)
  }
})

test_that("beta column association and public source roles are immutable", {
  fixture <- .grid_cross_client_fixture()
  mutations <- list(
    function(x) { x$spec$predictor_order <- rev(x$spec$predictor_order); x },
    function(x) { x$spec$beta_grid[[2L]] <- list(0, 0, 1); x },
    function(x) { x$spec$beta_encoded[[2L]][[2L]] <- "0"; x },
    function(x) { x$spec$design_terms <- rev(x$spec$design_terms); x },
    function(x) { x$spec$computation_peers <- list("site_a", "site_c"); x },
    function(x) { x$source_contract$recipients <- list("site_a", "site_c"); x },
    function(x) { x$spec$alignment$method <- "unjoined_sources"; x },
    function(x) { x$source_contract$alignment$method <- "unjoined_sources"; x },
    function(x) { x$source_contract$private_layout$blocks[[1L]]$owner_peer <- "site_c"; x },
    function(x) { x$source_contract$private_layout$blocks[[3L]]$value_fraction_bits <- 50; x },
    function(x) { x$source_contract$private_layout$padding_validity <- 1; x })
  for (mutate in mutations) {
    .grid_cross_client_reject(fixture$sign(mutate(fixture$contract)), fixture)
  }
})

test_that("profile and artifact versions are closed and transcript safe", {
  fixture <- .grid_cross_client_fixture()
  mutations <- list(
    function(x) { x$spec$version <- "binomial_grid_v1"; x },
    function(x) { x$artifact$version <- "bounded-binomial-likelihood-grid-v1"; x },
    function(x) { x$source_contract$version <- "dsvert-biomedical-capsule-source-contract-v1"; x },
    function(x) { x$spec$numeric_contract$profile_sha256 <- strrep("0", 64L); x },
    function(x) { x$spec$numeric_contract$rounding_rule <- "floor"; x },
    function(x) { x$spec$numeric_contract$input_fraction_bits <- 48; x },
    function(x) { x$spec$numeric_contract$arithmetic_width_bits <- 128; x },
    function(x) { x$artifact$transcript$row_batch_size <- 100; x },
    function(x) { x$artifact$implementation_state <- "same_owner_materialized"; x },
    function(x) { x$artifact$unexpected <- "private-candidate-patient-loss"; x },
    function(x) { x$source_contract$purpose <- "general_nonlinear_rpc"; x })
  for (mutate in mutations) {
    .grid_cross_client_reject(fixture$sign(mutate(fixture$contract)), fixture)
  }
  bad <- fixture$contract
  bad$signatures$site_a <- strrep("A", 86L)
  .grid_cross_client_reject(bad, fixture)
})

test_that("individual new validators reject malformed objects with one public error", {
  fixture <- .grid_cross_client_fixture()
  validate <- list(
    function() .dsvert_dp_glm_grid_cross_numeric_validate(list(), "binomial"),
    function() .dsvert_dp_glm_grid_cross_spec_validate(
      list(), fixture$policy, fixture$schema),
    function() .dsvert_dp_glm_grid_cross_artifact_validate(
      list(), fixture$contract$spec),
    function() .dsvert_dp_glm_grid_cross_source_contract_validate(
      list(), fixture$contract$spec, fixture$contract$artifact))
  for (operation in validate) {
    error <- tryCatch(operation(), error = identity)
    expect_s3_class(error, "dsvert_dp_public_failure")
    expect_identical(conditionMessage(error),
      "[dsvert_dp_public_failure:v1] Protected capsule operation failed.")
  }
})

test_that("signed schema structure cannot add ownership or snapshot ambiguity", {
  fixture <- .grid_cross_client_fixture()
  mutations <- list(
    function(x) { x$version <- "unsupported"; x },
    function(x) { x$logical_snapshot$alignment_protocol_version <- 0; x },
    function(x) { x$logical_snapshot$unknown <- TRUE; x },
    function(x) { x$peer_pinset_sha256 <- strrep("0", 64L); x },
    function(x) { x$datasets$aligned$unknown <- TRUE; x },
    function(x) { x$datasets$aligned$patient_keys$site_b <- NULL; x },
    function(x) { x$datasets$aligned$patient_keys$outsider <- "id"; x },
    function(x) { x$datasets$aligned$columns$x$owner_peer <- "outsider"; x },
    function(x) { x$datasets$aligned$columns$x$lower <- 1; x },
    function(x) { x$datasets$aligned$columns$x$upper <- Inf; x },
    function(x) { x$datasets$aligned$columns$x$kind <- "private_computed"; x },
    function(x) { x$datasets$aligned$columns$x$unknown <- TRUE; x },
    function(x) { names(x$datasets$aligned$columns)[1L] <- "site_c$x"; x },
    function(x) { names(x$datasets$aligned$columns)[1L] <- "site_b$x$"; x },
    function(x) { x$datasets$aligned$columns[["site_b$x"]] <- x$datasets$aligned$columns$x; x },
    function(x) { x$datasets$other <- x$datasets$aligned; x })
  for (mutate in mutations) {
    # Retain original signatures: the parser must reject structure before it
    # could use a changed ownership/bounds claim as authority.
    .grid_cross_client_reject(fixture$contract, fixture,
                               mutate(fixture$schema_manifest))
  }
})

test_that("every numeric proof field is mandatory and public profile contents are pinned", {
  for (family in c("binomial", "poisson")) {
    fixture <- .grid_cross_client_fixture(family)
    numeric <- fixture$contract$spec$numeric_contract
    for (field in names(numeric)) {
      bad <- numeric
      bad[[field]] <- NULL
      expect_error(.dsvert_dp_glm_grid_cross_numeric_validate(bad, family),
                   class = "dsvert_dp_public_failure")
    }
    for (field in names(numeric$per_operation_bounds)) {
      bad <- numeric
      bad$per_operation_bounds[[field]] <- NULL
      expect_error(.dsvert_dp_glm_grid_cross_numeric_validate(bad, family),
                   class = "dsvert_dp_public_failure")
    }
  }
})

test_that("cross grid dimensions coefficients and exact integer caps fail closed", {
  fixture <- .grid_cross_client_fixture()
  mutations <- list(
    function(x) { x$spec$analysis_id <- "invalid id"; x },
    function(x) { x$spec$dataset <- "unregistered"; x },
    function(x) { x$spec$predictor_order <- list("x", "z"); x },
    function(x) { x$spec$predictor_order <- list("site_b$x", "site_b$x"); x },
    function(x) { x$spec$predictor_order <- as.list(paste0("site_b$x", 1:17)); x },
    function(x) { x$spec$outcome$reference <- "site_b$x"; x },
    function(x) { x$spec$beta_grid <- list(); x },
    function(x) { x$spec$beta_grid <- rep(x$spec$beta_grid[1L], 257L); x },
    function(x) { x$spec$beta_grid <- rep(x$spec$beta_grid[1L], 2L); x },
    function(x) { x$spec$beta_grid[[2L]] <- list(0, 1); x },
    function(x) { x$spec$beta_grid[[2L]] <- list(0, 9, 0); x },
    function(x) { x$spec$beta_grid[[2L]] <- list(8, 8, 1); x },
    function(x) { x$spec$beta_grid <- rev(x$spec$beta_grid); x },
    function(x) { x$spec$max_outcome <- 0.5; x },
    function(x) { x$spec$numeric_grid_bits <- 7; x },
    function(x) { x$spec$observation_capacity <- 0; x },
    function(x) { x$spec$sensitivity$maximum_coordinates[[1L]] <- 2^53; x },
    function(x) { x$spec$sensitivity$raw_l1_sensitivity <- 0; x },
    function(x) { x$spec$sensitivity$natural_l2_sensitivity <- 0; x },
    function(x) { x$spec$unknown <- TRUE; x },
    function(x) { x$spec$family <- "gaussian"; x })
  for (mutate in mutations) {
    .grid_cross_client_reject(fixture$sign(mutate(fixture$contract)), fixture)
  }
  for (field in c("numeric_grid_bits", "unit_capacity")) {
    bad <- fixture
    bad$policy[[field]] <- 0
    .grid_cross_client_reject(bad$contract, bad)
  }
  bad <- fixture
  bad$policy$adjacency <- "unsupported"
  .grid_cross_client_reject(bad$contract, bad)
})
