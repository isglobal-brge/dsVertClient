.grouped_cross_client_fixture <- function(family = "lmm", correlation = "independence", owners = 2L) {
  peers <- paste0("site_", letters[seq_len(owners)])
  keys <- stats::setNames(lapply(peers, function(peer) openssl::ed25519_keygen()), peers)
  b64 <- function(value) sub("=+$", "", chartr("+/", "-_",
    gsub("[\r\n]", "", jsonlite::base64_enc(value))))
  pins <- vapply(keys, function(key) b64(tail(as.raw(as.list(key)$pubkey), 32L)),
                 character(1L))
  policy <- list(peer_pinset = pins,
    peer_pinset_sha256 = .dsvert_dp_capsule_source_hash(as.list(pins)),
    designated_noise_peers = peers[1:2], unit_capacity = 16,
    numeric_grid_bits = 8, adjacency = "add_remove_patient")
  maximum <- if (startsWith(family, "poisson")) 4 else 1
  covariate_peers <- if (owners == 2L) peers else peers[seq.int(3L, owners)]
  covariate_names <- c("x", "z", "w")[seq_along(covariate_peers)]
  covariates <- setNames(lapply(seq_along(covariate_peers), function(i) {
    list(kind = "numeric", owner_peer = covariate_peers[i],
      lower = if (i == 1L) -1 else 0, upper = if (i == 1L) 1 else 2)
  }), covariate_names)
  predictor_order <- paste0(covariate_peers, "$", covariate_names)
  snapshot <- list(logical_snapshot_id = "cohort", version = "v1",
                   alignment_protocol_version = 1)
  unsigned <- list(version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = snapshot, peer_pinset_sha256 = policy$peer_pinset_sha256,
    datasets = list(aligned = list(dataset_id = "aligned", dataset_version = "v1",
      schema_version = "v1", alignment_group = "group",
      patient_keys = setNames(as.list(rep("patient", owners)), peers),
      columns = c(covariates, list(
        y = list(kind = "numeric", owner_peer = "site_b", lower = 0, upper = maximum),
        cluster = list(kind = "categorical", owner_peer = "site_a",
                       levels = c("a", "b", "c", "d")))))))
  sign_message <- function(message) lapply(keys, function(key) {
    b64(openssl::ed25519_sign(message, key))
  })
  sign_schema <- function(value) {
    value$signatures <- NULL
    value$signatures <- sign_message(charToRaw(paste0(
      .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
      .dsvert_joint_dp_client_json(value))))
    value
  }
  manifest <- sign_schema(unsigned)
  schema <- .dsvert_dp_glm_grid_cross_schema_validate(policy, snapshot,
    manifest, .dsvert_dp_glm_grid_cross_verify)
  parameters <- if (family == "lmm") {
    list(residual_variance = 1, random_intercept_variance = .25)
  } else if (grepl("_glmm$", family)) {
    list(random_intercept_variance = .25, quadrature = "gh5_fixed_v1")
  } else list(correlation = correlation, rho = if (correlation == "independence") 0 else .25,
              score_clip = 2)
  raw <- list(version = paste0(family, "_grid_cross_v1"), analysis_id = "grouped",
    dataset = "aligned", outcome = "site_b$y",
    predictor_order = predictor_order,
    beta_grid = list(rep(0, 1+length(covariates)),
      c(0, rep(if (length(covariates) > 2L) .25 else .5, length(covariates)))),
    max_outcome = maximum,
    alignment = list(version = "existing_prealigned_logical_dataset_v1",
      method = "pinned_psi_ordered_manifest_v1", alignment_group = "group",
      public_alignment_contract_sha256 = .dsvert_dp_capsule_source_hash(list(
        logical_snapshot = snapshot, alignment_group = "group",
        method = "pinned_psi_ordered_manifest_v1")), public_patient_dependent_hash = FALSE),
    grouping = list(reference = "site_a$cluster", cluster_capacity = 4,
      max_patients_per_cluster = 4, patient_rule = "one_analysis_row_per_patient_v1",
      ordering = "stable_signed_slots_preserve_gaps_v1"), parameters = parameters)
  raw$beta_grid <- raw$beta_grid[order(vapply(raw$beta_grid,
    .dsvert_joint_dp_client_json, character(1L)), method = "radix")]
  spec <- .dsvert_dp_grouped_grid_cross_spec(raw, policy, schema)
  artifact <- .dsvert_dp_grouped_grid_cross_artifact(spec)
  sign <- function(value) {
    value$signatures <- NULL
    value$signatures <- sign_message(.dsvert_dp_grouped_grid_cross_message(value))
    value
  }
  contract <- sign(list(version = .DSVERT_CLIENT_DP_GROUPED_GRID_CROSS_VERSION,
    spec = spec, artifact = artifact,
    source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, artifact)))
  list(contract = contract, policy = policy, schema = schema,
       schema_manifest = manifest, raw = raw, sign = sign, sign_schema = sign_schema)
}
