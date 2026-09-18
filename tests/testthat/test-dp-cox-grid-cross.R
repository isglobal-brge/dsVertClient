.cox_cross_client_fixture <- function(capacity = 16, beta_grid = NULL) {
  peers <- c("site_a", "site_b")
  keys <- stats::setNames(lapply(peers, function(peer) openssl::ed25519_keygen()), peers)
  b64 <- function(value) sub("=+$", "", chartr("+/", "-_",
    gsub("[\r\n]", "", jsonlite::base64_enc(value))))
  pins <- vapply(keys, function(key) b64(tail(as.raw(as.list(key)$pubkey), 32L)),
                 character(1L))
  policy <- list(peer_pinset = pins,
    peer_pinset_sha256 = .dsvert_dp_capsule_source_hash(as.list(pins)),
    designated_noise_peers = peers, unit_capacity = capacity,
    numeric_grid_bits = 8, adjacency = "add_remove_patient")
  snapshot <- list(logical_snapshot_id = "cohort", version = "v1",
                   alignment_protocol_version = 1)
  schema_unsigned <- list(version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = snapshot, peer_pinset_sha256 = policy$peer_pinset_sha256,
    datasets = list(aligned = list(dataset_id = "aligned", dataset_version = "v1",
      schema_version = "v1", alignment_group = "group",
      patient_keys = list(site_a = "patient", site_b = "patient"),
      columns = list(
        event = list(kind = "numeric", owner_peer = "site_a", lower = 0, upper = 1),
        time = list(kind = "numeric", owner_peer = "site_a", lower = 0, upper = 20),
        x = list(kind = "numeric", owner_peer = "site_a", lower = 0, upper = 1),
        z = list(kind = "numeric", owner_peer = "site_b", lower = 0, upper = 1)))))
  sign_schema <- function(value) {
    value$signatures <- NULL
    message <- charToRaw(paste0(.DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
                               .dsvert_joint_dp_client_json(value)))
    c(value, list(signatures = lapply(keys, function(key) {
      b64(openssl::ed25519_sign(message, key))
    })))
  }
  schema_manifest <- sign_schema(schema_unsigned)
  schema <- .dsvert_dp_glm_grid_cross_schema_validate(policy, snapshot,
    schema_manifest, .dsvert_dp_glm_grid_cross_verify)
  if (is.null(beta_grid)) beta_grid <- list(c(0, 0), c(1, -0.5), c(1, 0))
  beta_grid <- beta_grid[order(vapply(beta_grid, function(beta) {
    .dsvert_joint_dp_client_json(as.list(beta))
  }, character(1L)), method = "radix")]
  raw <- list(version = "cox_grid_cross_v1", analysis_id = "cox_grid",
    dataset = "aligned", time = "site_a$time", event = "site_a$event",
    predictor_order = c("site_a$x", "site_b$z"), beta_grid = beta_grid,
    alignment = list(version = "existing_prealigned_logical_dataset_v1",
      method = "pinned_psi_ordered_manifest_v1", alignment_group = "group",
      public_alignment_contract_sha256 = .dsvert_dp_capsule_source_hash(list(
        logical_snapshot = schema$unsigned$logical_snapshot,
        alignment_group = "group", method = "pinned_psi_ordered_manifest_v1")),
      public_patient_dependent_hash = FALSE))
  spec <- .dsvert_dp_cox_grid_cross_spec(raw, policy, schema)
  artifact <- .dsvert_dp_cox_grid_cross_artifact(spec)
  sign <- function(value) {
    value$signatures <- NULL
    message <- .dsvert_dp_cox_grid_cross_message(value)
    c(value, list(signatures = lapply(keys, function(key) {
      b64(openssl::ed25519_sign(message, key))
    })))
  }
  contract <- sign(list(version = .DSVERT_CLIENT_DP_COX_GRID_CROSS_CONTRACT_VERSION,
    spec = spec, artifact = artifact,
    source_contract = .dsvert_dp_cox_grid_cross_source_contract(spec, artifact)))
  list(contract = contract, policy = policy, schema_manifest = schema_manifest,
       schema = schema, raw = raw, keys = keys, b64 = b64,
       sign = sign, sign_schema = sign_schema)
}

test_that("Cox client accepts only a literal additive Surv formula", {
  value <- .dsvert_dp_cox_grid_cross_formula(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z)
  expect_identical(value$predictors, c("site_a$x", "site_b$z"))
  for (formula in list(y ~ x, Surv(t, e) ~ log(x), Surv(t, e) ~ x * z,
      Surv(t, e) ~ x + strata(z), Surv(t, e) ~ x - z,
      Surv(start, t, e) ~ x, Surv(t, e) ~ x + x)) {
    expect_error(.dsvert_dp_cox_grid_cross_formula(formula),
                 class = "dsvert_dp_public_failure")
  }
})

test_that("Cox cross-owner contracts require both real custodian signatures", {
  fixture <- .cox_cross_client_fixture()
  validate <- function(value) .dsvert_dp_cox_grid_cross_contract_validate(
    value, fixture$policy, fixture$schema_manifest)
  expect_identical(validate(fixture$contract),
                   .dsvert_joint_dp_client_canonical(fixture$contract))
  for (peer in names(fixture$contract$signatures)) {
    changed <- fixture$contract
    changed$signatures[[peer]] <- NULL
    expect_error(validate(changed), class = "dsvert_dp_public_failure")
    changed <- fixture$contract
    changed$signatures[[peer]] <- strrep("A", 86)
    expect_error(validate(changed), class = "dsvert_dp_public_failure")
  }
  for (change in list(
      function(x) { x$spec$ties <- "efron"; x },
      function(x) { x$spec$intercept <- TRUE; x },
      function(x) { x$spec$time_semantics <- "binned"; x },
      function(x) { x$spec$numeric_contract <- list(); x },
      function(x) { x$spec$sensitivity$raw_l1_sensitivity <- 0; x },
      function(x) { x$artifact$implementation_state <- "same_owner_materialized"; x },
      function(x) { x$source_contract$purpose <- "generic_nonlinear_rpc"; x })) {
    expect_error(validate(fixture$sign(change(fixture$contract))),
                 class = "dsvert_dp_public_failure")
  }
})

test_that("Cox client canonical argmin rescales slopes and omits inference", {
  fixture <- .cox_cross_client_fixture()
  spec <- fixture$contract$spec
  coordinates <- rep(0, length(spec$beta_grid))
  moment <- .dsvert_dp_cox_grid_cross_moment(coordinates, spec)
  expect_identical(moment$selected_candidate, 1L)
  expect_equal(unname(moment$coefficients), unlist(spec$beta_grid[[1L]]))
  for (bad in list(-coordinates - 1, coordinates + 0.5,
                   rep(Inf, length(coordinates)), coordinates[-1L],
                   unlist(spec$sensitivity$maximum_coordinates) + 1)) {
    expect_error(.dsvert_dp_cox_grid_cross_moment(bad, spec),
                 class = "dsvert_dp_public_failure")
  }
  release <- function(...) c(fixture[c("contract", "policy", "schema_manifest")],
                             list(coordinates = coordinates))
  result <- .dsvert_dp_cox_grid_cross_impl(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL), release)
  expect_identical(result$selected_candidate, 1L)
  expect_null(result$std_errors)
  expect_null(result$baseline_hazard)
  expect_false(result$production_ready)
  expect_false(result$protected_optimizer_called)
  expect_error(.dsvert_dp_cox_grid_cross_impl(
    Surv(site_a$time, site_a$event) ~ site_a$x + site_b$wrong,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL), release),
    class = "dsvert_dp_public_failure")
})

test_that("production Cox client is closed and exposes no test evaluator", {
  expect_identical(names(formals(dp_cox_grid)),
    c("formula", "data", "analysis_id", "datasources"))
  expect_error(dp_cox_grid(Surv(site_a$time, site_a$event) ~ site_a$x + site_b$z,
    "aligned", "cox_grid", list(site_a = NULL, site_b = NULL)),
    class = "dsvert_dp_public_failure")
  expect_false(.dsvert_dp_cox_grid_cross_client_register()$runtime_enabled)
})
