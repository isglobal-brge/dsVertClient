.cox_cross_client_fixture <- function(capacity = 16, beta_grid = NULL, owners = 2L) {
  peers <- paste0("site_", letters[seq_len(owners)])
  keys <- stats::setNames(lapply(peers, function(peer) openssl::ed25519_keygen()), peers)
  b64 <- function(value) sub("=+$", "", chartr("+/", "-_",
    gsub("[\r\n]", "", jsonlite::base64_enc(value))))
  pins <- vapply(keys, function(key) b64(tail(as.raw(as.list(key)$pubkey), 32L)),
                 character(1L))
  policy <- list(peer_pinset = pins,
    peer_pinset_sha256 = .dsvert_dp_capsule_source_hash(as.list(pins)),
    designated_noise_peers = peers[1:2], unit_capacity = capacity,
    numeric_grid_bits = 8, adjacency = "add_remove_patient")
  snapshot <- list(logical_snapshot_id = "cohort", version = "v1",
                   alignment_protocol_version = 1)
  schema_unsigned <- list(version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = snapshot, peer_pinset_sha256 = policy$peer_pinset_sha256,
    datasets = list(aligned = list(dataset_id = "aligned", dataset_version = "v1",
      schema_version = "v1", alignment_group = "group",
      patient_keys = setNames(rep(list("patient"), owners), peers),
      columns = list(
        event = list(kind = "numeric", owner_peer = "site_a", lower = 0, upper = 1),
        time = list(kind = "numeric", owner_peer = "site_a", lower = 0, upper = 20),
        x = list(kind = "numeric", owner_peer = "site_a", lower = 0, upper = 1),
        z = list(kind = "numeric", owner_peer = "site_b", lower = 0, upper = 1)))))
  extra_predictors <- character()
  if (owners > 2L) for (peer in peers[-c(1L, 2L)]) {
    column <- paste0("x_", peer)
    schema_unsigned$datasets$aligned$columns[[column]] <- list(
      kind = "numeric", owner_peer = peer, lower = 0, upper = 1)
    extra_predictors <- c(extra_predictors, paste0(peer, "$", column))
  }
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
  if (is.null(beta_grid)) beta_grid <- lapply(list(c(0, 0), c(1, -0.5), c(1, 0)),
    function(beta) c(beta, rep(0, length(extra_predictors))))
  beta_grid <- beta_grid[order(vapply(beta_grid, function(beta) {
    .dsvert_joint_dp_client_json(as.list(beta))
  }, character(1L)), method = "radix")]
  raw <- list(version = "cox_grid_cross_v1", analysis_id = "cox_grid",
    dataset = "aligned", time = "site_a$time", event = "site_a$event",
    predictor_order = c("site_a$x", "site_b$z", extra_predictors), beta_grid = beta_grid,
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
