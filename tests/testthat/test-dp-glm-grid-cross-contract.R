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

test_that("cross-grid profile coefficients, proof and numeric template are immutable", {
  fixture <- jsonlite::fromJSON(system.file(
    "cross-grid-v1", "numeric_profile_v1.json", package = "dsVertClient"),
    simplifyVector = FALSE)
  mutations <- list(
    function(x) { x$profile$softplus_coefficients_q64[[1]] <- "0"; x },
    function(x) { x$certificate$eta_error_bound <- "0"; x },
    function(x) { x$numeric_contract$binomial$input_fraction_bits <- 48; x })
  for (mutate in mutations) {
    changed <- mutate(fixture)
    testthat::with_mocked_bindings({
      expect_error(.dsvert_dp_glm_grid_cross_numeric("binomial"),
                   class = "dsvert_dp_public_failure")
    }, fromJSON = function(...) changed, .package = "jsonlite")
  }
})

test_that("signed categorical schema labels share the server canonical representations", {
  fixture <- .grid_cross_client_fixture()
  cases <- list(
    list(value = c("B", "A"), labels = c("A", "B")),
    list(value = c(TRUE, FALSE), labels = c("FALSE", "TRUE")),
    list(value = c(10L, -2L), labels = c("-2", "10")),
    list(value = c(1e10, -0), labels = c("0", "10000000000")),
    list(value = factor(c("B", "A")), labels = c("A", "B")),
    list(value = ordered(c("B", "A")), labels = c("A", "B")))
  for (case in cases) {
    schema <- fixture$schema_manifest
    schema$datasets$aligned$columns$category <- list(
      kind = "categorical", owner_peer = "site_b", levels = case$labels)
    signed <- fixture$sign_schema(schema)
    signed$datasets$aligned$columns$category$levels <- case$value
    authenticated <- .dsvert_dp_glm_grid_cross_schema_validate(
      fixture$policy, signed$logical_snapshot, signed,
      .dsvert_dp_glm_grid_cross_verify)
    expect_identical(authenticated$unsigned$datasets$aligned$columns$category$levels,
                     case$labels)
  }
})

test_that("real signatures cannot authorize unsupported categorical label encodings", {
  fixture <- .grid_cross_client_fixture()
  cases <- list(
    list(value = c(1.5, 2.5), labels = c("1.5", "2.5")),
    list(value = Inf, labels = "Inf"),
    list(value = 2^53, labels = "9007199254740992"),
    list(value = as.raw(1), labels = "01"),
    list(value = 1+1i, labels = "1+1i"),
    list(value = structure("A", class = "custom"), labels = "A"),
    list(value = matrix("A"), labels = "A"),
    list(value = structure(2L, class = "factor", levels = "A"), labels = "A"))
  for (case in cases) {
    schema <- fixture$schema_manifest
    schema$datasets$aligned$columns$category <- list(
      kind = "categorical", owner_peer = "site_b", levels = case$labels)
    signed <- fixture$sign_schema(schema)
    signed$datasets$aligned$columns$category$levels <- case$value
    expect_error(.dsvert_dp_glm_grid_cross_schema_validate(
      fixture$policy, signed$logical_snapshot, signed,
      .dsvert_dp_glm_grid_cross_verify), class = "dsvert_dp_public_failure")
  }
})

.grid_cross_client_validate <- function(fixture, value) {
  .dsvert_dp_glm_grid_cross_contract_validate(
    value, fixture$policy, fixture$schema_manifest)
}

test_that("cross-grid signed zero encoding survives canonical round trips", {
  fixture <- .grid_cross_client_fixture()
  for (zero in c(-0, -2^-52)) {
    raw <- fixture$raw
    raw$beta_grid <- list(c(zero, 0, 0))
    spec <- .dsvert_dp_glm_grid_cross_spec(raw, fixture$policy,
                                         fixture$schema)
    expect_identical(spec$beta_encoded[[1]], as.list(c("0", "0", "0")))
    artifact <- .dsvert_dp_glm_grid_cross_artifact(spec)
    contract <- fixture$sign(list(
      version = .DSVERT_CLIENT_DP_GLM_GRID_CROSS_CONTRACT_VERSION,
      spec = spec, artifact = artifact,
      source_contract = .dsvert_dp_glm_grid_cross_source_contract(spec, artifact)))
    validated <- .grid_cross_client_validate(fixture, contract)
    encoded <- jsonlite::fromJSON(.dsvert_joint_dp_client_json(validated),
                                  simplifyVector = FALSE)
    expect_identical(.grid_cross_client_validate(fixture, encoded), validated)
  }
})

test_that("cross-grid beta L1 comparison is exact at the public boundary", {
  expect_true(.dsvert_dp_glm_grid_cross_beta_l1_valid(c(8, 8)))
  expect_true(.dsvert_dp_glm_grid_cross_beta_l1_valid(rep(0.5, 32)))
  expect_true(.dsvert_dp_glm_grid_cross_beta_l1_valid(
    c(8 - 6 * 2^-50, rep(2 + 3 * 2^-51, 4))))
  expect_false(.dsvert_dp_glm_grid_cross_beta_l1_valid(rep(3.2, 5)))
  expect_false(.dsvert_dp_glm_grid_cross_beta_l1_valid(c(8, 8, 2^-60)))
  expect_false(.dsvert_dp_glm_grid_cross_beta_l1_valid(c(8, 8, 2^-1074)))
  fixture <- .grid_cross_client_fixture()
  raw <- fixture$raw
  for (beta in list(c(8, 8, 0), c(7.5, 7.5, 1))) {
    raw$beta_grid <- list(beta)
    expect_no_error(.dsvert_dp_glm_grid_cross_spec(
      raw, fixture$policy, fixture$schema))
  }
  raw$beta_grid <- list(c(8, 8, 2^-60))
  expect_error(.dsvert_dp_glm_grid_cross_spec(
    raw, fixture$policy, fixture$schema), class = "dsvert_dp_public_failure")
})

test_that("cross-grid scalar identifiers reject JSON array representations", {
  fixture <- .grid_cross_client_fixture()
  for (field in c("analysis_id", "dataset", "outcome")) {
    raw <- fixture$raw
    raw[[field]] <- as.list(raw[[field]])
    expect_error(.dsvert_dp_glm_grid_cross_spec(
      raw, fixture$policy, fixture$schema), class = "dsvert_dp_public_failure")
  }
  forged <- fixture$contract[setdiff(names(fixture$contract), "signatures")]
  forged$spec$analysis_id <- list("grid")
  forged$artifact <- .dsvert_dp_glm_grid_cross_artifact(forged$spec)
  forged$source_contract <- .dsvert_dp_glm_grid_cross_source_contract(
    forged$spec, forged$artifact)
  expect_error(.grid_cross_client_validate(fixture, fixture$sign(forged)),
               class = "dsvert_dp_public_failure")
})

test_that("piecewise admission pins new arithmetic, caps and all signatures", {
  for (family in c("binomial", "poisson")) {
    fixture <- .grid_cross_client_fixture(family)
    spec <- .dsvert_dp_glm_grid_cross_spec(
      fixture$raw, fixture$policy, fixture$schema, "piecewise_v2")
    artifact <- .dsvert_dp_glm_grid_cross_artifact(spec)
    value <- fixture$sign(list(version = fixture$contract$version, spec = spec,
      artifact = artifact,
      source_contract = .dsvert_dp_glm_grid_cross_source_contract(spec, artifact)))
    admitted <- .dsvert_dp_glm_grid_profile_admit(
      value, fixture$policy, fixture$schema_manifest)
    expect_identical(admitted$spec$numeric_contract$profile_identity,
                     .DSVERT_DP_GLM_GRID_PROFILE_V2)
    caps <- vapply(spec$sensitivity$candidate_bounds, `[[`, numeric(1),
                   "per_patient_cap")
    expect_equal(spec$sensitivity$raw_l1_sensitivity, sum(caps))
    expect_gte(spec$sensitivity$raw_l2_sensitivity, sqrt(sum(caps^2)))
    expect_identical(spec$beta_encoded, fixture$contract$spec$beta_encoded)
    expect_identical(.dsvert_dp_glm_grid_cross_layout(spec),
                     fixture$contract$source_contract$private_layout)
    expect_error(.dsvert_dp_glm_grid_profile_admit(
      fixture$contract, fixture$policy, fixture$schema_manifest),
      class = "dsvert_dp_public_failure")
    for (field in c("profile_sha256", "certificate_sha256", "profile_envelope")) {
      bad <- value
      bad$spec$numeric_contract[[field]] <- "tampered"
      expect_error(.dsvert_dp_glm_grid_profile_admit(
        fixture$sign(bad), fixture$policy, fixture$schema_manifest),
        class = "dsvert_dp_public_failure")
    }
    bad <- value
    bad$spec$sensitivity$raw_l1_sensitivity <- 1
    expect_error(.dsvert_dp_glm_grid_profile_admit(
      fixture$sign(bad), fixture$policy, fixture$schema_manifest),
      class = "dsvert_dp_public_failure")
    bad <- value
    bad$signatures$site_b <- NULL
    expect_error(.dsvert_dp_glm_grid_profile_admit(
      bad, fixture$policy, fixture$schema_manifest), class = "dsvert_dp_public_failure")
  }
})

test_that("materialized cross grids select only from authenticated signed candidates", {
  for (family in c("binomial", "poisson")) {
    f <- .grid_cross_client_fixture(family)
    spec <- .dsvert_dp_glm_grid_cross_spec(f$raw, f$policy, f$schema, "piecewise_v2")
    core <- .dsvert_dp_glm_grid_cross_artifact(spec)
    signed <- f$sign(list(version = f$contract$version, spec = spec, artifact = core,
      source_contract = .dsvert_dp_glm_grid_cross_source_contract(spec, core)))
    artifact <- .dsvert_dp_glm_grid_cross_workload_artifact(signed)
    manifest <- list(admission = list(unit_capacity = spec$observation_capacity,
      adjacency = spec$adjacency), bounds = list(numeric_grid_bits = spec$numeric_grid_bits),
      workload = list(coordinate_count = artifact$coordinate_count + 1,
        capsule_mechanism = list(mechanism = "discrete-laplace"),
        families = list(gaussian_models = list(
        artifacts = list(cross_grid = artifact)))))
    context <- list(pinset = f$policy$peer_pinset, designated = f$policy$designated_noise_peers)
    fragments <- list(describe = list(), survival = list(), vertical_cross = list(),
      gaussian = list(cross_grid = list(version = spec$version, dataset = spec$dataset,
        contract = artifact$signed_contract)))
    expect_equal(.dsvert_dp_capsule_manifest_fragments(fragments),
      .dsvert_joint_dp_client_canonical(fragments))
    schema_json <- .dsvert_joint_dp_client_json(f$schema_manifest)
    expect_true(.dsvert_dp_glm_grid_cross_preflight(manifest, context, schema_json))
    unsupported <- manifest
    unsupported$workload$capsule_mechanism$mechanism <- "discrete-gaussian"
    expect_error(.dsvert_dp_glm_grid_cross_preflight(unsupported, context, schema_json),
      class = "dsvert_dp_public_failure")
    validated <- .dsvert_dp_glm_grid_artifact(manifest, spec$dataset, spec$analysis_id,
      NULL, spec$adjacency, 2^spec$numeric_grid_bits, spec$observation_capacity, family)
    selected <- .dsvert_dp_glm_grid_moment(c(20, 10), validated, family)
    expect_equal(selected$selected_candidate, 2L)
    expect_equal(unname(selected$normalized_coefficients), unlist(spec$beta_grid[[2]]))
    expect_identical(validated$implementation_state, "cross_owner_exact_gc_materialized")
    for (mutation in c("signature", "cap", "profile", "wrapper")) {
      bad <- signed
      if (mutation == "signature") bad$signatures[[1]] <- NULL
      if (mutation == "cap") bad$spec$sensitivity$candidate_bounds[[1]]$per_patient_cap <- 1
      if (mutation == "profile") bad$spec$numeric_contract$profile_sha256 <- paste(rep("0", 64), collapse = "")
      altered <- .dsvert_dp_glm_grid_cross_workload_artifact(bad)
      if (mutation == "wrapper") altered$statistic_maximum[[1]] <- 1
      broken <- manifest
      broken$workload$families$gaussian_models$artifacts$cross_grid <- altered
      expect_error(.dsvert_dp_glm_grid_cross_preflight(broken, context, schema_json),
        class = "dsvert_dp_public_failure")
    }
  }
  expect_identical(.dsvert_dp_glm_grid_formula_reference(quote(site_a$x)), "site_a$x")
  expect_null(.dsvert_dp_glm_grid_formula_reference(quote(log(site_a$x))))
  expect_null(.dsvert_dp_glm_grid_formula_reference(quote(site_a$x$y)))
})

test_that("grid Claim coverage includes owners without public moment blocks", {
  peers <- c("site_a", "site_b", "site_c")
  context <- list(servers = peers, all_conns = stats::setNames(as.list(peers), peers))
  artifact <- list(version = unname(.DSVERT_CLIENT_DP_GLM_GRID_CROSS_ARTIFACT_VERSIONS[[1]]),
    participating_peers = as.list(peers))
  trusted <- list(manifest = list(workload = list(families = list(
    gaussian_models = list(artifacts = list(grid = artifact))))))
  observed <- NULL
  testthat::local_mocked_bindings(.dsvert_fanout_by_site = function(conns, calls, ...) {
    observed <<- names(calls)
    stop("claim coverage captured")
  })
  expect_error(.dsvert_dp_synopsis_runner_compile(context,
    list(manifest_sha256 = strrep("a", 64)), trusted,
    list(blocks = list(list(owner_peer = "site_a"))), function(...) NULL),
    "claim coverage captured")
  expect_identical(observed, peers)
})

test_that("both certified profiles admit the full signed ten-predictor envelope", {
  for (family in c("binomial", "poisson")) {
    f <- .grid_cross_client_fixture(family)
    f$policy$unit_capacity <- 10000
    columns <- sprintf("x%02d", 1:10)
    owners <- rep(c("site_a", "site_b", "site_c"), c(4, 3, 3))
    schema <- f$schema_manifest
    schema$datasets$aligned$patient_keys <- list(site_a = "patient",
      site_b = "patient", site_c = "patient")
    schema$datasets$aligned$columns <- c(stats::setNames(lapply(seq_along(columns),
      function(i) list(kind = "numeric", owner_peer = owners[[i]], lower = 0,
                       upper = 1)), columns), list(y = list(kind = "numeric",
        owner_peer = "site_c", lower = 0, upper = f$raw$max_outcome)))
    schema <- f$sign_schema(schema)
    authenticated <- .dsvert_dp_glm_grid_cross_schema_validate(
      f$policy, schema$logical_snapshot, schema, .dsvert_dp_glm_grid_cross_verify)
    raw <- f$raw
    raw$predictor_order <- paste0(owners, "$", columns)
    raw$beta_grid <- lapply(seq(-0.25, 0.25, length.out = 50), function(offset)
      c(-0.5 + offset, rep(c(0.25, -0.25), 5)))
    raw$beta_grid <- raw$beta_grid[order(vapply(raw$beta_grid,
      function(beta) .dsvert_joint_dp_client_json(as.list(beta)), character(1)),
      method = "radix")]
    spec <- .dsvert_dp_glm_grid_cross_spec(raw, f$policy, authenticated, "piecewise_v2")
    artifact <- .dsvert_dp_glm_grid_cross_artifact(spec)
    contract <- f$sign(list(version = f$contract$version, spec = spec,
      artifact = artifact,
      source_contract = .dsvert_dp_glm_grid_cross_source_contract(spec, artifact)))
    admitted <- .dsvert_dp_glm_grid_profile_admit(contract, f$policy, schema)
    expect_equal(admitted$spec$observation_capacity, 10000)
    expect_length(admitted$spec$predictor_order, 10)
    expect_length(admitted$spec$beta_grid, 50)
    expect_identical(unlist(admitted$spec$participating_peers),
                     c("site_a", "site_b", "site_c"))
    expect_identical(admitted$artifact$implementation_state,
                     "cross_owner_exact_gc_materialized")
    malformed <- raw
    malformed$predictor_order[c(8, 9, 10)] <- raw$predictor_order[c(9, 10, 8)]
    expect_error(.dsvert_dp_glm_grid_cross_spec(malformed, f$policy, authenticated,
      "piecewise_v2"), class = "dsvert_dp_public_failure")
  }
})
