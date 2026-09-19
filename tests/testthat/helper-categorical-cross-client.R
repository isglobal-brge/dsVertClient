.categorical_cross_client_fixture <- function(family = "multinomial",
                                             bits = 8L,
                                             adjacency = "add_remove_patient") {
  peers <- c("site_a", "site_b")
  keys <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  b64 <- function(value) sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))))
  pins <- vapply(keys, function(key) {
    b64(tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  policy <- list(peer_pinset = pins,
    peer_pinset_sha256 = .dsvert_dp_capsule_source_hash(as.list(pins)),
    designated_noise_peers = peers, unit_capacity = 20,
    numeric_grid_bits = bits, adjacency = adjacency)
  snapshot <- list(logical_snapshot_id = "cohort", version = "v1",
                   alignment_protocol_version = 1)
  schema_unsigned <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = snapshot,
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    datasets = list(aligned = list(
      dataset_id = "aligned", dataset_version = "v1", schema_version = "v1",
      alignment_group = "group",
      patient_keys = list(site_a = "patient", site_b = "patient"),
      columns = list(
        x = list(kind = "numeric", owner_peer = "site_a", lower = -1, upper = 1),
        y = list(kind = "categorical", owner_peer = "site_a",
                 levels = c("A", "B", "C")),
        z = list(kind = "numeric", owner_peer = "site_b", lower = 0,
                 upper = 10))))))
  sign_schema <- function(value) {
    value$signatures <- NULL
    value <- .dsvert_joint_dp_client_canonical(value)
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
    analysis_id = "cross_grid", dataset = "aligned", outcome = "site_a$y",
    predictor_order = c("site_a$x", "site_b$z"),
    alignment = list(version = "existing_prealigned_logical_dataset_v1",
      method = "pinned_psi_ordered_manifest_v1", alignment_group = "group",
      public_alignment_contract_sha256 = .dsvert_dp_capsule_source_hash(list(
        logical_snapshot = schema$unsigned$logical_snapshot,
        alignment_group = "group", method = "pinned_psi_ordered_manifest_v1")),
      public_patient_dependent_hash = FALSE))
  if (identical(family, "multinomial")) {
    raw$levels <- c("A", "B", "C")
    raw$reference <- "A"
    raw$beta_grid <- list(c(0, 0, 0, 0, 0, 0), c(0, 1, 1, 0, -1, -1))
  } else {
    raw$ordered_levels <- c("B", "A", "C")
    raw$candidate_grid <- list(
      list(beta = c(0, 0, 0), thresholds = c(-1, 1)),
      list(beta = c(0, 1, 1), thresholds = c(-0.5, 0.5)))
  }
  registration <- if (identical(family, "multinomial")) {
    .dsvert_dp_multinomial_grid_cross_register()
  } else .dsvert_dp_ordinal_grid_cross_register()
  spec <- registration$spec(raw, policy, schema)
  artifact <- registration$artifact(spec)
  source <- .dsvert_dp_glm_grid_cross_source_contract(spec, artifact)
  sign <- function(value) {
    value$signatures <- NULL
    message <- .dsvert_dp_glm_grid_cross_message(value)
    c(value, list(signatures = lapply(keys, function(key) {
      b64(openssl::ed25519_sign(message, key))
    })))
  }
  contract <- sign(list(version = .DSVERT_CLIENT_DP_GLM_GRID_CROSS_CONTRACT_VERSION,
                       spec = spec, artifact = artifact, source_contract = source))
  list(contract = contract, policy = policy, schema_manifest = schema_manifest,
       schema = schema, sign = sign, sign_schema = sign_schema, raw = raw,
       registration = registration)
}

.categorical_cross_client_reject <- function(operation) {
  error <- tryCatch(operation(), error = identity)
  expect_s3_class(error, "dsvert_dp_public_failure")
  expect_identical(conditionMessage(error),
    "[dsvert_dp_public_failure:v1] Protected capsule operation failed.")
}
