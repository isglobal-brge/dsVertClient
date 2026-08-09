# Internal client orchestration for the purpose-bound Ring128 capsule source
# transport. The analyst relays authenticated ciphertext only. Plaintext
# coordinates and additive shares remain inside the pinned data nodes.

.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CONTRACT_VERSION <-
  "dsvert-biomedical-capsule-source-contract-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_VERSION <-
  "dsvert-biomedical-capsule-source-recipient-ticket-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_SUMMARY_VERSION <-
  "dsvert-biomedical-capsule-source-summary-v2-streaming"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_VERSION <-
  "dsvert-biomedical-capsule-source-encrypted-chunk-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_BUNDLE_VERSION <-
  "dsvert-biomedical-capsule-source-encrypted-bundle-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_ACK_VERSION <-
  "dsvert-biomedical-capsule-source-ack-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION <-
  "dsvert-biomedical-capsule-source-byte-window-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT <-
  "dsvert-biomedical-capsule-source-adaptive-window-v2"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_NEGOTIATION_VERSION <-
  "dsvert-biomedical-capsule-source-ticket-negotiation-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ATTESTATION_VERSION <-
  "dsvert-biomedical-capsule-source-capability-attestation-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ONLY_REQUEST <-
  "dsvert-biomedical-capsule-source-capability-only-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ONLY_REQUEST_V2 <-
  "dsvert-biomedical-capsule-source-capability-only-v2"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_SUMMARY_NEGOTIATION_VERSION <-
  "dsvert-biomedical-capsule-source-summary-negotiation-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_ACK_WINDOW_VERSION <-
  "dsvert-biomedical-capsule-source-ack-window-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_RECEIPT_VERSION <-
  "dsvert-biomedical-capsule-source-client-receipt-v1"
.DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION <-
  "dsvert-biomedical-capsule-workload-v7"
.DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE <-
  "registered_lifecycle_available_requires_runtime_preflight"
.DSVERT_CLIENT_DP_CAPSULE_RELEASE_LIFECYCLE_VERSION <-
  "dsvert-biomedical-capsule-release-lifecycle-v1"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE <-
  "biomedical_capsule_ring128_source_shares_only"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CROSS_PURPOSE <-
  "biomedical_capsule_ring128_cross_gaussian_inputs_and_release_shares_only"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CATEGORICAL_CROSS_PURPOSE <-
  "biomedical_capsule_ring128_cross_categorical_inputs_and_release_shares_only"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_SIGNATURE_DOMAIN <-
  "dsVert/biomedical-capsule/source-transport/v1/"
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES <- 8192L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES <- 32L * 1024L^2
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_TICKET_BYTES <- 64L * 1024L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_SUMMARY_BYTES <- 64L * 1024L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_ENVELOPE_BYTES <- 1024L * 1024L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_BUNDLE_BYTES <- 768L * 1024L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_LEGACY_MAX_WINDOW_BYTES <- 768L * 1024L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_LEGACY_MAX_ACCEPT_WINDOW_BYTES <-
  1024L * 1024L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES <- 8L * 1024L^2
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_ACCEPT_WINDOW_BYTES <- 8L * 1024L^2
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_FRAMING_RESERVE <- 64L * 1024L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_ALLOCATION_BYTES <- 64L * 1024L
.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_CHUNKS <- 8L

.dsvert_dp_capsule_client_registered_release_lifecycle <- function() {
  .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_RELEASE_LIFECYCLE_VERSION,
    state = "registered_productive_joint_dp_vector_lifecycle",
    manifest_authority = TRUE,
    included_coordinate_producers = TRUE,
    source_secret_sharing = TRUE,
    recipient_encrypted_transport = TRUE,
    sampler_integration = TRUE,
    confidential_finalizer = TRUE,
    durable_sticky_replay = TRUE,
    package_integration_verified = TRUE,
    live_connector_execution = "runtime_preflight_required",
    raw_intermediate_releasable = FALSE,
    analyst_can_bypass_lifecycle = FALSE))
}

.dsvert_dp_capsule_source_hex <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[0-9a-f]{64}$", value)
}

.dsvert_dp_capsule_source_window_capability <- function() {
  list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT,
    maximum_window_chunks =
      as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_CHUNKS),
    maximum_response_bytes =
      as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES),
    maximum_accept_bytes =
      as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_ACCEPT_WINDOW_BYTES),
    byte_bounded_prefix = TRUE,
    scalar_legacy_byte_identical = TRUE,
    operation_or_request_limit = FALSE)
}

.dsvert_dp_capsule_source_legacy_window_capability <- function() {
  list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION,
    maximum_window_chunks =
      as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_CHUNKS),
    maximum_response_bytes =
      as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_LEGACY_MAX_WINDOW_BYTES),
    maximum_accept_bytes =
      as.numeric(
        .DSVERT_CLIENT_DP_CAPSULE_SOURCE_LEGACY_MAX_ACCEPT_WINDOW_BYTES),
    byte_bounded_prefix = TRUE,
    scalar_legacy_byte_identical = TRUE,
    operation_or_request_limit = FALSE)
}

.dsvert_dp_capsule_source_hash <- function(value) {
  digest::digest(
    .dsvert_joint_dp_client_json(value), algo = "sha256",
    serialize = FALSE)
}

.dsvert_dp_capsule_allocation_openings <- function(context) {
  openings <- context$allocation_openings
  designated <- context$designated
  valid <- is.list(openings) && !is.null(names(openings)) &&
    !anyNA(names(openings)) && !anyDuplicated(names(openings)) &&
    identical(names(openings), designated) &&
    all(vapply(openings, function(value) {
      is.character(value) && length(value) == 1L && !is.na(value) &&
        nzchar(value) && nchar(value, type = "bytes") <=
          .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_ALLOCATION_BYTES
    }, logical(1L)))
  if (!isTRUE(valid)) {
    stop("The capsule source has no canonical cross-signed allocation openings",
         call. = FALSE)
  }
  invisible(lapply(openings, function(value) {
    .dsvert_joint_dp_client_decode(
      value, "biomedical allocation opening",
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_ALLOCATION_BYTES)
    value
  }))
}

.dsvert_dp_capsule_source_string_list <- function(value, what,
                                                   expected_length = NULL) {
  valid <- is.list(value) && is.null(names(value)) && length(value) > 0L &&
    all(vapply(value, function(item) {
      .dsvert_dp_is_string(item) &&
        grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", item)
    }, logical(1L)))
  result <- if (isTRUE(valid)) unname(unlist(value, use.names = FALSE)) else
    character()
  if (!isTRUE(valid) || anyDuplicated(result) ||
      !identical(result, sort(result, method = "radix")) ||
      (!is.null(expected_length) && length(result) != expected_length)) {
    stop("Invalid biomedical capsule source ", what, " returned by a peer",
         call. = FALSE)
  }
  result
}

.dsvert_dp_capsule_source_decode_b64url <- function(value, maximum, what) {
  if (!.dsvert_dp_is_string(value) ||
      !grepl("^[A-Za-z0-9_-]+$", value) || nchar(value) %% 4L == 1L ||
      nchar(value, type = "bytes") > 4 * ceiling(maximum / 3) + 4L) {
    stop("Invalid biomedical capsule source ", what, call. = FALSE)
  }
  padding <- (4L - nchar(value) %% 4L) %% 4L
  decoded <- tryCatch(
    jsonlite::base64_dec(paste0(
      chartr("-_", "+/", value), strrep("=", padding))),
    error = function(e) raw(0L))
  canonical <- if (is.raw(decoded)) {
    sub("=+$", "", chartr(
      "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(decoded))),
      perl = TRUE)
  } else {
    ""
  }
  if (!is.raw(decoded) || length(decoded) > maximum ||
      !identical(canonical, value)) {
    stop("Invalid biomedical capsule source ", what, call. = FALSE)
  }
  decoded
}

.dsvert_dp_capsule_source_verify <- function(
    value, domain, peer, context, ciphertext = FALSE) {
  if (!is.list(value) || !"signature" %in% names(value) ||
      !peer %in% names(context$pinset)) {
    stop("Invalid biomedical capsule source signed object", call. = FALSE)
  }
  unsigned <- value[setdiff(names(value), "signature")]
  if (isTRUE(ciphertext)) unsigned$ciphertext <- NULL
  message <- charToRaw(paste0(
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SIGNATURE_DOMAIN, domain, "|",
    .dsvert_joint_dp_client_json(unsigned)))
  public <- .dsvert_joint_dp_client_b64url(
    unname(context$pinset[[peer]]), 32L, "capsule source identity key")
  signature <- .dsvert_joint_dp_client_b64url(
    value$signature, 64L, "capsule source signature")
  verified <- tryCatch(openssl::ed25519_verify(
    message, signature, openssl::read_ed25519_pubkey(public)),
    error = function(e) FALSE)
  if (!isTRUE(verified)) {
    stop("A pinned peer returned an invalid biomedical capsule source signature",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_capsule_source_policy_identity <- function(context) {
  policy <- context$status[[context$servers[[1L]]]]$policy
  common <- list(
    protocol = "dsvert-joint-dp-control-v3",
    release_scope = "joint_mpc_single_opening",
    capability_id = "joint_mpc_single_opening_v1",
    domain = policy$domain,
    cohort_id = policy$cohort_id,
    ordered_peer_pinset = as.list(context$pinset),
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    peer_count = length(context$pinset),
    designated_noise_peers = context$designated,
    designated_noise_peer_pinset = as.list(
      context$pinset[context$designated]),
    epsilon_capsule = as.numeric(policy$capsule_epsilon),
    delta_capsule = as.numeric(policy$capsule_delta),
    lifetime_max_distinct_capsules =
      as.numeric(policy$lifetime_max_distinct_capsules),
    lifetime_epsilon_upper_bound = policy$lifetime_epsilon_upper_bound,
    lifetime_delta_upper_bound = policy$lifetime_delta_upper_bound,
    privacy_accounting =
      "bounded_distinct_capsules_one_public_instance_each_v1",
    adjacency = policy$adjacency,
    patient_column = policy$patient_column,
    unit_capacity = as.numeric(policy$unit_capacity),
    max_records_per_unit = as.numeric(policy$max_records_per_unit),
    overflow_policy = policy$overflow_policy,
    sampler = policy$sampler)
  policy_hash <- .dsvert_dp_capsule_source_hash(common)
  list(
    consortium_id = paste0("jdpc1_", policy_hash),
    policy_contract_hash = policy_hash)
}

.dsvert_dp_capsule_scope_strings <- function(
    value, what, allow_null = FALSE) {
  if (is.null(value) && isTRUE(allow_null)) value <- list()
  if (is.character(value) && is.null(names(value))) value <- as.list(value)
  if (!is.list(value) || !is.null(names(value)) ||
      !all(vapply(value, function(item) {
        is.character(item) && length(item) == 1L && !is.na(item) &&
          nzchar(item)
      }, logical(1L)))) {
    stop("The signed biomedical primitive scope has invalid ", what,
         call. = FALSE)
  }
  if (!length(value)) return(character())
  result <- unname(unlist(value, use.names = FALSE))
  if (anyDuplicated(result) ||
      !identical(result, sort(result, method = "radix"))) {
    stop("The signed biomedical primitive scope has invalid ", what,
         call. = FALSE)
  }
  result
}

.dsvert_dp_capsule_scope_pairs <- function(value, what) {
  if (!is.list(value) || !is.null(names(value))) {
    stop("The signed biomedical primitive scope has invalid ", what,
         call. = FALSE)
  }
  pairs <- lapply(value, function(pair) {
    if (is.list(pair) && is.null(names(pair))) {
      valid <- length(pair) == 2L && all(vapply(pair, function(item) {
        is.character(item) && length(item) == 1L && !is.na(item) &&
          nzchar(item)
      }, logical(1L)) )
      if (valid) pair <- unname(unlist(pair, use.names = FALSE))
    }
    if (!is.character(pair) || length(pair) != 2L ||
        !is.null(names(pair)) || anyNA(pair) || any(!nzchar(pair)) ||
        anyDuplicated(pair) ||
        !identical(pair, sort(pair, method = "radix"))) {
      stop("The signed biomedical primitive scope has invalid ", what,
           call. = FALSE)
    }
    unname(pair)
  })
  keys <- vapply(
    pairs, function(pair) .dsvert_joint_dp_client_json(as.list(pair)),
    character(1L))
  if (anyDuplicated(keys) ||
      !identical(keys, sort(keys, method = "radix"))) {
    stop("The signed biomedical primitive scope has invalid ", what,
         call. = FALSE)
  }
  pairs
}

.dsvert_dp_capsule_primitive_scope_validate <- function(manifest) {
  invalid <- function() {
    stop("The signed biomedical primitive scope is invalid", call. = FALSE)
  }
  if (!is.list(manifest)) invalid()
  workload <- tryCatch(manifest$workload, error = function(error) NULL)
  scope <- if (is.list(workload)) workload$primitive_scope else NULL
  scope_fields <- c(
    "version", "mode", "authority", "analyst_expandable",
    "client_query_can_add_coordinates", "consensus", "mismatch_behavior",
    "compatibility_default", "recommended_deployment_mode",
    "selection_sha256", "selection", "projected_cost")
  selection_fields <- c(
    "mode", "explicit_catalog", "referenced_by_signed_specs", "included")
  explicit_fields <- c(
    "numeric_moments", "categorical_marginals", "categorical_pairs",
    "correlations")
  reference_fields <- c(
    "numeric", "categorical", "describe", "survival", "gaussian",
    "vertical_cross")
  included_fields <- c(
    "numeric_moments", "categorical_marginals",
    "same_owner_categorical_pairs", "same_owner_correlations")
  cost_fields <- c(
    "schema_numeric_column_count", "schema_categorical_column_count",
    "possible_same_owner_numeric_pair_count",
    "possible_same_owner_categorical_pair_count",
    "included_numeric_moment_count",
    "included_categorical_marginal_count",
    "included_numeric_pair_count", "included_categorical_pair_count",
    "included_cross_categorical_pair_count",
    "numeric_moment_coordinate_count", "numeric_pair_coordinate_count",
    "categorical_marginal_coordinate_count",
    "categorical_pair_coordinate_count", "gaussian_model_coordinate_count",
    "projected_coordinate_count", "projected_integer_l1_sensitivity",
    "projected_integer_l2_sensitivity", "automatic_pair_expansion",
    "scaling_contract")
  if (!identical(manifest$version,
                 .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION) ||
      !identical(manifest$capsule_schema,
                 .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION) ||
      !.dsvert_dp_is_integer(workload$coordinate_count, 1,
                             .DSVERT_DP_MAX_COORDINATES) ||
      !.dsvert_dp_has_exact_names(scope, scope_fields) ||
      !.dsvert_dp_has_exact_names(scope$selection, selection_fields) ||
      !identical(names(scope$selection),
                 sort(selection_fields, method = "radix")) ||
      !.dsvert_dp_has_exact_names(
        scope$selection$explicit_catalog, explicit_fields) ||
      !identical(names(scope$selection$explicit_catalog),
                 sort(explicit_fields, method = "radix")) ||
      !.dsvert_dp_has_exact_names(
        scope$selection$referenced_by_signed_specs, reference_fields) ||
      !identical(names(scope$selection$referenced_by_signed_specs),
                 sort(reference_fields, method = "radix")) ||
      !.dsvert_dp_has_exact_names(scope$selection$included,
                                  included_fields) ||
      !identical(names(scope$selection$included),
                 sort(included_fields, method = "radix")) ||
      !.dsvert_dp_has_exact_names(scope$projected_cost, cost_fields) ||
      !identical(scope$version,
                 "dsvert-biomedical-capsule-primitive-scope-v1") ||
      !.dsvert_dp_is_string(scope$mode) ||
      !scope$mode %in% c("all_schema", "catalog_v1") ||
      !identical(scope$authority,
                 "custodian_policy_and_signed_workload_specs_only") ||
      !identical(scope$analyst_expandable, FALSE) ||
      !identical(scope$client_query_can_add_coordinates, FALSE) ||
      !identical(scope$consensus, paste0(
        "byte_identical_manifest_hash_with_all_pinned_peer_build_",
        "signatures_required_before_source_access")) ||
      !identical(scope$mismatch_behavior,
                 "reject_before_protected_snapshot_resolution") ||
      !identical(scope$compatibility_default, "all_schema") ||
      !identical(scope$recommended_deployment_mode, "catalog_v1") ||
      !identical(scope$selection$mode, scope$mode) ||
      !.dsvert_dp_capsule_source_hex(scope$selection_sha256) ||
      !identical(scope$selection_sha256,
                 .dsvert_dp_capsule_source_hash(scope$selection))) {
    invalid()
  }
  included <- tryCatch(list(
    numeric = .dsvert_dp_capsule_scope_strings(
      scope$selection$included$numeric_moments, "numeric moments"),
    categorical = .dsvert_dp_capsule_scope_strings(
      scope$selection$included$categorical_marginals,
      "categorical marginals"),
    categorical_pairs = .dsvert_dp_capsule_scope_pairs(
      scope$selection$included$same_owner_categorical_pairs,
      "categorical pairs"),
    correlations = .dsvert_dp_capsule_scope_pairs(
      scope$selection$included$same_owner_correlations,
      "correlation pairs")), error = function(error) NULL)
  explicit <- tryCatch(list(
    numeric = .dsvert_dp_capsule_scope_strings(
      scope$selection$explicit_catalog$numeric_moments,
      "explicit numeric moments"),
    categorical = .dsvert_dp_capsule_scope_strings(
      scope$selection$explicit_catalog$categorical_marginals,
      "explicit categorical marginals"),
    categorical_pairs = .dsvert_dp_capsule_scope_pairs(
      scope$selection$explicit_catalog$categorical_pairs,
      "explicit categorical pairs"),
    correlations = .dsvert_dp_capsule_scope_pairs(
      scope$selection$explicit_catalog$correlations,
      "explicit correlation pairs")), error = function(error) NULL)
  references <- tryCatch(list(
    numeric = .dsvert_dp_capsule_scope_strings(
      scope$selection$referenced_by_signed_specs$numeric,
      "signed-spec numeric references"),
    categorical = .dsvert_dp_capsule_scope_strings(
      scope$selection$referenced_by_signed_specs$categorical,
      "signed-spec categorical references"),
    describe = .dsvert_dp_capsule_scope_strings(
      scope$selection$referenced_by_signed_specs$describe,
      "signed describe references", allow_null = TRUE),
    survival = .dsvert_dp_capsule_scope_strings(
      scope$selection$referenced_by_signed_specs$survival,
      "signed survival references", allow_null = TRUE),
    gaussian = .dsvert_dp_capsule_scope_strings(
      scope$selection$referenced_by_signed_specs$gaussian,
      "signed Gaussian references", allow_null = TRUE),
    vertical_cross = .dsvert_dp_capsule_scope_strings(
      scope$selection$referenced_by_signed_specs$vertical_cross,
      "signed vertical-cross references", allow_null = TRUE)),
    error = function(error) NULL)
  if (is.null(included) || is.null(explicit) || is.null(references)) invalid()
  pair_keys <- function(value) vapply(
    value, function(pair) .dsvert_joint_dp_client_json(as.list(pair)),
    character(1L))
  explicit_pair_keys <- pair_keys(explicit$categorical_pairs)
  explicit_correlation_keys <- pair_keys(explicit$correlations)
  included_pair_keys <- pair_keys(included$categorical_pairs)
  included_correlation_keys <- pair_keys(included$correlations)
  expected_catalog_numeric <- sort(unique(c(
    explicit$numeric, references$numeric,
    unlist(explicit$correlations, use.names = FALSE))), method = "radix")
  expected_catalog_categorical <- sort(unique(c(
    explicit$categorical, references$categorical,
    unlist(explicit$categorical_pairs, use.names = FALSE))),
    method = "radix")
  catalog_valid <- if (identical(scope$mode, "catalog_v1")) {
    identical(included$numeric, expected_catalog_numeric) &&
      identical(included$categorical, expected_catalog_categorical) &&
      identical(included_pair_keys, explicit_pair_keys) &&
      identical(included_correlation_keys, explicit_correlation_keys)
  } else {
    !length(explicit$numeric) && !length(explicit$categorical) &&
      !length(explicit$categorical_pairs) && !length(explicit$correlations) &&
      all(references$numeric %in% included$numeric) &&
      all(references$categorical %in% included$categorical)
  }
  if (!isTRUE(catalog_valid)) invalid()

  cost <- scope$projected_cost
  count_names <- cost_fields[seq_len(15L)]
  counts <- vapply(count_names, function(field) {
    value <- cost[[field]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 0 || value != floor(value) ||
        value > 2^53 - 1) return(NA_real_)
    as.numeric(value)
  }, numeric(1L))
  l1 <- cost$projected_integer_l1_sensitivity
  l2 <- cost$projected_integer_l2_sensitivity
  expected_automatic <- if (identical(scope$mode, "all_schema")) {
    "explicit_all_schema"
  } else {
    "none"
  }
  expected_scaling <- if (identical(scope$mode, "all_schema")) {
    "explicit_schema_wide_pair_families_may_be_quadratic"
  } else {
    paste0(
      "linear_in_declared_univariates_plus_explicit_pairs_and_",
      "declared_model_cross_products")
  }
  if (anyNA(counts) ||
      !is.numeric(l1) || length(l1) != 1L || is.na(l1) ||
      !is.finite(l1) || l1 <= 0 ||
      !is.numeric(l2) || length(l2) != 1L || is.na(l2) ||
      !is.finite(l2) || l2 <= 0 ||
      !identical(cost$automatic_pair_expansion, expected_automatic) ||
      !identical(cost$scaling_contract, expected_scaling) ||
      counts[["projected_coordinate_count"]] !=
        as.numeric(workload$coordinate_count) ||
      counts[["included_numeric_moment_count"]] !=
        length(included$numeric) ||
      counts[["included_categorical_marginal_count"]] !=
        length(included$categorical) ||
      counts[["included_numeric_pair_count"]] !=
        length(included$correlations) ||
      counts[["included_categorical_pair_count"]] !=
        length(included$categorical_pairs) ||
      (identical(scope$mode, "all_schema") &&
       (counts[["included_numeric_moment_count"]] !=
          counts[["schema_numeric_column_count"]] ||
        counts[["included_categorical_marginal_count"]] !=
          counts[["schema_categorical_column_count"]] ||
        counts[["included_numeric_pair_count"]] !=
          counts[["possible_same_owner_numeric_pair_count"]] ||
        counts[["included_categorical_pair_count"]] !=
          counts[["possible_same_owner_categorical_pair_count"]])) ||
      counts[["possible_same_owner_numeric_pair_count"]] <
        counts[["included_numeric_pair_count"]] ||
      counts[["possible_same_owner_categorical_pair_count"]] <
        counts[["included_categorical_pair_count"]]) {
    invalid()
  }
  scope
}

.dsvert_dp_capsule_source_manifest <- function(manifest_json, context) {
  manifest <- .dsvert_joint_dp_client_decode(
    manifest_json, "biomedical capsule manifest",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_MANIFEST_BYTES)
  manifest_fields <- c(
    "version", "logical_snapshot", "capsule_schema", "admission", "bounds",
    "workload", "capsule_identity", "execution_state")
  identity_fields <- c("capsule_id", "contract")
  contract_fields <- c(
    "protocol", "consortium_id", "policy_contract_hash",
    "peer_pinset_sha256", "logical_snapshot", "capsule_schema",
    "admission", "bounds", "workload", "privacy_epoch_scope")
  identity <- manifest$capsule_identity
  contract <- if (is.list(identity)) identity$contract else NULL
  release_coordinate_count <- if (is.list(manifest$workload)) {
    manifest$workload$coordinate_count
  } else NULL
  mechanism <- if (is.list(manifest$workload)) {
    manifest$workload$capsule_mechanism
  } else NULL
  policy_identity <- .dsvert_dp_capsule_source_policy_identity(context)
  reference <- context$status[[context$servers[[1L]]]]
  primitive_scope <- tryCatch(
    .dsvert_dp_capsule_primitive_scope_validate(manifest),
    error = function(error) NULL)
  valid <- .dsvert_dp_has_exact_names(manifest, manifest_fields) &&
    identical(manifest$version,
              .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION) &&
    identical(manifest$capsule_schema,
              .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION) &&
    .dsvert_dp_has_exact_names(identity, identity_fields) &&
    .dsvert_dp_has_exact_names(contract, contract_fields) &&
    identical(contract$protocol, "dsvert-joint-dp-capsule-identity-v3") &&
    .dsvert_dp_capsule_source_hex(identity$capsule_id) &&
    identical(identity$capsule_id,
              .dsvert_dp_capsule_source_hash(contract)) &&
    identical(contract$logical_snapshot, manifest$logical_snapshot) &&
    identical(contract$capsule_schema, manifest$capsule_schema) &&
    identical(contract$admission, manifest$admission) &&
    identical(contract$bounds, manifest$bounds) &&
    identical(contract$workload, manifest$workload) &&
    identical(contract$consortium_id, policy_identity$consortium_id) &&
    identical(contract$policy_contract_hash,
              policy_identity$policy_contract_hash) &&
    identical(contract$peer_pinset_sha256,
              reference$policy$peer_pinset_sha256) &&
    .dsvert_dp_is_integer(
      release_coordinate_count, 1, .DSVERT_DP_MAX_COORDINATES) &&
    identical(manifest$workload$workload_version,
              .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION) &&
    identical(manifest$execution_state,
              .DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE) &&
    identical(manifest$workload$execution_state,
              .DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE) &&
    identical(manifest$workload$declared_workload_fully_materialized, TRUE) &&
    identical(manifest$workload$package_family_coverage_complete, FALSE) &&
    identical(
      manifest$workload$registered_release_lifecycle,
      .dsvert_dp_capsule_client_registered_release_lifecycle()) &&
    is.list(primitive_scope) &&
    is.list(mechanism) &&
    .dsvert_dp_capsule_source_hex(mechanism$source_context_hash) &&
    identical(contract$privacy_epoch_scope,
              "per_peer_signed_receipts_v1")
  if (!isTRUE(valid)) {
    stop("The immutable biomedical capsule manifest is invalid or misbound",
         call. = FALSE)
  }
  cross_artifacts <- .dsvert_dp_gaussian_cross_artifacts_client(manifest)
  categorical_cross_artifacts <-
    .dsvert_dp_categorical_cross_artifacts_client(manifest)
  cross_layout <- if (length(cross_artifacts) ||
                      length(categorical_cross_artifacts)) {
    .dsvert_dp_gaussian_cross_layout_client(manifest)
  } else {
    NULL
  }
  coordinate_count <- if (!is.null(cross_layout)) {
    cross_layout$transport_coordinate_count
  } else {
    release_coordinate_count
  }
  list(
    value = manifest, capsule_id = identity$capsule_id,
    coordinate_count = as.numeric(coordinate_count),
    release_coordinate_count = as.numeric(release_coordinate_count),
    release_coordinate_order_sha256 = if (!is.null(cross_layout)) {
      cross_layout$release_coordinate_order_sha256
    } else NULL,
    private_layout_sha256 = if (!is.null(cross_layout)) {
      cross_layout$transport_coordinate_order_sha256
    } else NULL,
    cross_enabled = !is.null(cross_layout),
    purpose = if (length(categorical_cross_artifacts)) {
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CATEGORICAL_CROSS_PURPOSE
    } else if (!is.null(cross_layout)) {
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CROSS_PURPOSE
    } else {
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE
    },
    source_context_hash = mechanism$source_context_hash)
}

.dsvert_dp_capsule_source_ticket <- function(
    ticket_json, peer, context, manifest) {
  ticket <- .dsvert_joint_dp_client_decode(
    ticket_json, "biomedical capsule source recipient ticket",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_TICKET_BYTES)
  fields <- c(
    "version", "phase", "purpose", "capsule_id", "contract_hash",
    "recipient_name", "recipient_identity_pk", "transport_key_id",
    "transport_pk", "peer_pinset_sha256", "designated_noise_peers",
    "source_peers", "coordinate_count", "chunk_coordinates",
    "chunk_count", "persistent", "ready_for_sampling", "signature")
  sources <- tryCatch(.dsvert_dp_capsule_source_string_list(
    ticket$source_peers, "source peer list"), error = function(e) NULL)
  recipients <- tryCatch(.dsvert_dp_capsule_source_string_list(
    ticket$designated_noise_peers, "recipient list", 2L),
    error = function(e) NULL)
  transport <- tryCatch(.dsvert_dp_capsule_source_decode_b64url(
    ticket$transport_pk, 32L, "transport public key"),
    error = function(e) NULL)
  expected_chunks <- ceiling(
    manifest$coordinate_count /
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES)
  valid <- .dsvert_dp_has_exact_names(ticket, fields) &&
    identical(ticket$version,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_VERSION) &&
    identical(ticket$phase, "recipient_key_committed") &&
    identical(ticket$purpose, manifest$purpose) &&
    identical(ticket$capsule_id, manifest$capsule_id) &&
    .dsvert_dp_capsule_source_hex(ticket$contract_hash) &&
    identical(ticket$recipient_name, peer) &&
    identical(ticket$recipient_identity_pk,
              unname(context$pinset[[peer]])) &&
    identical(ticket$peer_pinset_sha256,
              context$status[[peer]]$policy$peer_pinset_sha256) &&
    !is.null(sources) && all(sources %in% context$servers) &&
    !is.null(recipients) && identical(recipients, context$designated) &&
    .dsvert_dp_is_integer(
      ticket$coordinate_count, 1,
      if (isTRUE(manifest$cross_enabled)) {
        .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_TRANSPORT_COORDINATES
      } else {
        .DSVERT_DP_MAX_COORDINATES
      }) &&
    identical(as.numeric(ticket$coordinate_count),
              manifest$coordinate_count) &&
    .dsvert_dp_is_integer(ticket$chunk_coordinates, 1,
                          .DSVERT_DP_MAX_COORDINATES) &&
    identical(as.numeric(ticket$chunk_coordinates),
              as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES)) &&
    .dsvert_dp_is_integer(ticket$chunk_count, 1,
                          .DSVERT_DP_MAX_COORDINATES) &&
    identical(as.numeric(ticket$chunk_count), as.numeric(expected_chunks)) &&
    identical(ticket$persistent, TRUE) &&
    identical(ticket$ready_for_sampling, FALSE) &&
    !is.null(transport) && length(transport) == 32L &&
    .dsvert_dp_capsule_source_hex(ticket$transport_key_id) &&
    identical(ticket$transport_key_id,
              .dsvert_dp_capsule_source_hash(list(
                protocol = "dsvert-biomedical-capsule-source-key-id-v1",
                capsule_id = ticket$capsule_id,
                recipient_name = peer,
                transport_pk = ticket$transport_pk)))
  if (!isTRUE(valid)) {
    stop("A designated peer returned an invalid capsule source ticket",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(
    ticket, "recipient-ticket", peer, context)
  list(
    value = ticket, json = ticket_json,
    hash = .dsvert_dp_capsule_source_hash(ticket), sources = sources,
    recipients = recipients)
}

.dsvert_dp_capsule_source_ticket_response <- function(
    response_json, peer, context, manifest) {
  preview <- .dsvert_joint_dp_client_decode(
    response_json, "biomedical capsule source ticket negotiation",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_TICKET_BYTES)
  if (identical(preview$version,
                .DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_VERSION)) {
    ticket <- .dsvert_dp_capsule_source_ticket(
      response_json, peer, context, manifest)
    ticket$negotiation <- NULL
    return(ticket)
  }
  fields <- c(
    "version", "phase", "ticket_json", "ticket_sha256", "capability",
    "signature")
  ticket <- tryCatch(
    .dsvert_dp_capsule_source_ticket(
      preview$ticket_json, peer, context, manifest),
    error = function(error) NULL)
  current_capability <- identical(
    .dsvert_joint_dp_client_json(preview$capability),
    .dsvert_joint_dp_client_json(
      .dsvert_dp_capsule_source_window_capability()))
  legacy_capability <- identical(
    .dsvert_joint_dp_client_json(preview$capability),
    .dsvert_joint_dp_client_json(
      .dsvert_dp_capsule_source_legacy_window_capability()))
  valid <- .dsvert_dp_has_exact_names(preview, fields) &&
    identical(
      preview$version,
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_TICKET_NEGOTIATION_VERSION) &&
    identical(preview$phase, "recipient_transport_window_attested") &&
    is.list(ticket) &&
    identical(preview$ticket_sha256, ticket$hash) &&
    (isTRUE(current_capability) || isTRUE(legacy_capability))
  if (!isTRUE(valid)) {
    stop("A designated peer returned an invalid source transport negotiation",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(
    preview, "recipient-window-capability", peer, context)
  if (isTRUE(legacy_capability)) {
    ticket$negotiation <- NULL
    return(ticket)
  }
  ticket$negotiation <- list(
    value = preview, json = response_json,
    hash = .dsvert_dp_capsule_source_hash(preview),
    capability = preview$capability)
  ticket
}

.dsvert_dp_capsule_source_ticket_set <- function(
    responses, context, manifest) {
  peers <- context$designated
  if (!is.list(responses) || is.null(names(responses)) ||
      !setequal(names(responses), peers)) {
    stop("Capsule source tickets did not cover both designated peers",
         call. = FALSE)
  }
  tickets <- stats::setNames(lapply(peers, function(peer) {
    .dsvert_dp_capsule_source_ticket_response(
      responses[[peer]], peer, context, manifest)
  }), peers)
  common <- c(
    "version", "phase", "purpose", "capsule_id", "contract_hash",
    "peer_pinset_sha256", "designated_noise_peers", "source_peers",
    "coordinate_count", "chunk_coordinates", "chunk_count", "persistent",
    "ready_for_sampling")
  if (!identical(tickets[[1L]]$value[common],
                 tickets[[2L]]$value[common]) ||
      identical(tickets[[1L]]$value$transport_pk,
                tickets[[2L]]$value$transport_pk)) {
    stop("The designated peers returned conflicting capsule source tickets",
         call. = FALSE)
  }
  negotiations <- lapply(tickets, `[[`, "negotiation")
  if (all(!vapply(negotiations, is.null, logical(1L))) &&
      identical(negotiations[[1L]]$capability,
                negotiations[[2L]]$capability)) {
    attr(tickets, "window_capability") <- negotiations[[1L]]$capability
  }
  tickets
}

.dsvert_dp_capsule_source_capability_attestation <- function(
    response_json, source, context, manifest, contract_hash) {
  value <- .dsvert_joint_dp_client_decode(
    response_json, "biomedical capsule source capability attestation",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_TICKET_BYTES)
  fields <- c(
    "version", "phase", "capsule_id", "contract_hash", "source_name",
    "source_identity_pk", "capability", "signature")
  current_capability <- identical(
    .dsvert_joint_dp_client_json(value$capability),
    .dsvert_joint_dp_client_json(
      .dsvert_dp_capsule_source_window_capability()))
  legacy_capability <- identical(
    .dsvert_joint_dp_client_json(value$capability),
    .dsvert_joint_dp_client_json(
      .dsvert_dp_capsule_source_legacy_window_capability()))
  valid <- .dsvert_dp_has_exact_names(value, fields) &&
    identical(
      value$version,
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ATTESTATION_VERSION) &&
    identical(value$phase, "source_transport_capability_attested") &&
    identical(value$capsule_id, manifest$capsule_id) &&
    identical(value$contract_hash, contract_hash) &&
    identical(value$source_name, source) &&
    identical(value$source_identity_pk, unname(context$pinset[[source]])) &&
    (isTRUE(current_capability) || isTRUE(legacy_capability))
  if (!isTRUE(valid)) {
    stop("A source owner returned an invalid source capability attestation",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(
    value, "source-window-capability-advertisement", source, context)
  isTRUE(current_capability)
}

.dsvert_dp_capsule_source_optional_negotiation <- function(
    conns, expressions, .aggregate) {
  expressions <- .dsvert_dsi_text_frame_expressions(expressions)
  expression_fits <- tryCatch({
    .dsvert_validate_dsi_expression_sizes(expressions)
    TRUE
  }, dsvert_resource_oversize = function(error) FALSE)
  if (!isTRUE(expression_fits)) return(NULL)
  .dsvert_validate_real_dsi_transport(conns, .aggregate)
  targets <- names(conns)
  session_keys <- lapply(conns, .dsvert_dsi_job_session_key)
  poisoned <- targets[vapply(session_keys, function(key) {
    !is.null(key) && .dsvert_dsi_session_is_poisoned(key)
  }, logical(1L))]
  if (length(poisoned)) {
    stop(.dsvert_dsi_poisoned_session_condition(poisoned))
  }
  failed <- FALSE
  peer_rejection <- NULL
  resource_oversize <- NULL
  callback <- function(site, message) {
    failed <<- TRUE
    rejection <- .dsvert_client_parse_peer_not_recognized(message)
    oversized <- .dsvert_client_parse_resource_oversize(message)
    if (is.null(peer_rejection) && !is.null(rejection)) {
      peer_rejection <<- rejection
    }
    if (is.null(resource_oversize) && !is.null(oversized)) {
      resource_oversize <<- oversized
    }
    invisible(NULL)
  }
  result <- tryCatch(
    .dsvert_transport_aggregate(
      .aggregate = .aggregate, conns = conns, expr = expressions,
      async = FALSE, error = callback, errors.print = FALSE,
      require_settled_sync_failure = TRUE),
    interrupt = function(error) {
      .dsvert_dsi_poison_ambiguous_sessions(conns)
      stop(error)
    },
    error = function(error) {
      if (inherits(error, c(
          "dsvert_peer_not_recognized", "dsvert_dsi_poisoned_session"))) {
        stop(error)
      }
      rejection <- .dsvert_client_parse_peer_not_recognized(
        conditionMessage(error))
      oversized <- if (inherits(error, "dsvert_resource_oversize")) {
        error
      } else {
        .dsvert_client_parse_resource_oversize(conditionMessage(error))
      }
      if (!is.null(rejection)) peer_rejection <<- rejection
      if (!is.null(oversized)) resource_oversize <<- oversized
      if (is.null(rejection) && is.null(oversized)) {
        stop(.dsvert_dsi_poison_ambiguous_sessions(conns))
      }
      failed <<- TRUE
      NULL
    })
  if (inherits(peer_rejection, "dsvert_peer_not_recognized")) {
    stop(peer_rejection)
  }
  if (inherits(resource_oversize, "dsvert_resource_oversize")) return(NULL)
  valid <- !isTRUE(failed) && is.list(result) &&
    length(result) == length(targets) && !is.null(names(result)) &&
    !anyNA(names(result)) && !anyDuplicated(names(result)) &&
    setequal(names(result), targets) &&
    !any(vapply(result, is.null, logical(1L)))
  if (!isTRUE(valid)) return(NULL)
  result[targets]
}

.dsvert_dp_capsule_source_transfer_id <- function(
    contract_hash, capsule_id, source) {
  paste0("csrc_", .dsvert_dp_capsule_source_hash(list(
    protocol = "dsvert-biomedical-capsule-source-transfer-id-v1",
    contract_hash = contract_hash, capsule_id = capsule_id,
    source_name = source)))
}

.dsvert_dp_capsule_source_summary <- function(
    summary_json, source, context, ticket_reference) {
  summary <- .dsvert_joint_dp_client_decode(
    summary_json, "biomedical capsule source summary",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_SUMMARY_BYTES)
  fields <- c(
    "version", "phase", "purpose", "capsule_id", "contract_hash",
    "source_transfer_id", "source_name", "source_identity_pk",
    "recipients", "coordinate_count", "chunk_coordinates", "chunk_count",
    "ring_bits", "record_encoding", "emitted_chunk_durable_replay",
    "unmaterialized_requires_same_snapshot", "complete_durable_replay",
    "history_gate", "ready_for_sampling", "signature")
  recipients <- tryCatch(.dsvert_dp_capsule_source_string_list(
    summary$recipients, "summary recipient list", 2L),
    error = function(e) NULL)
  ticket <- ticket_reference$value
  valid <- .dsvert_dp_has_exact_names(summary, fields) &&
    identical(summary$version,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SUMMARY_VERSION) &&
    identical(summary$phase, "source_chunk_stream_ready") &&
    identical(summary$purpose, ticket$purpose) &&
    identical(summary$capsule_id, ticket$capsule_id) &&
    identical(summary$contract_hash, ticket$contract_hash) &&
    identical(summary$source_name, source) &&
    identical(summary$source_identity_pk,
              unname(context$pinset[[source]])) &&
    identical(summary$source_transfer_id,
              .dsvert_dp_capsule_source_transfer_id(
                ticket$contract_hash, ticket$capsule_id, source)) &&
    !is.null(recipients) && identical(recipients, context$designated) &&
    identical(as.numeric(summary$coordinate_count),
              as.numeric(ticket$coordinate_count)) &&
    identical(as.numeric(summary$chunk_coordinates),
              as.numeric(ticket$chunk_coordinates)) &&
    identical(as.numeric(summary$chunk_count),
              as.numeric(ticket$chunk_count)) &&
    identical(as.numeric(summary$ring_bits), 128) &&
    identical(summary$record_encoding,
              "little_endian_unsigned_fixed_16_bytes") &&
    identical(summary$emitted_chunk_durable_replay, TRUE) &&
    identical(summary$unmaterialized_requires_same_snapshot, TRUE) &&
    identical(summary$complete_durable_replay, FALSE) &&
    identical(summary$history_gate, FALSE) &&
    identical(summary$ready_for_sampling, FALSE)
  if (!isTRUE(valid)) {
    stop("A source owner returned an invalid capsule source summary",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(
    summary, "source-summary", source, context)
  list(value = summary, json = summary_json)
}

.dsvert_dp_capsule_source_summary_response <- function(
    response_json, source, context, ticket_reference, tickets) {
  preview <- .dsvert_joint_dp_client_decode(
    response_json, "biomedical capsule source summary negotiation",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_SUMMARY_BYTES)
  if (identical(preview$version,
                .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SUMMARY_VERSION)) {
    summary <- .dsvert_dp_capsule_source_summary(
      response_json, source, context, ticket_reference)
    summary$negotiation <- NULL
    return(summary)
  }
  fields <- c(
    "version", "phase", "summary_json", "summary_sha256",
    "ticket_negotiation_set_sha256", "capability", "signature")
  capability <- attr(tickets, "window_capability", exact = TRUE)
  negotiations <- lapply(context$designated, function(peer) {
    tickets[[peer]]$negotiation
  })
  negotiation_set_sha256 <- if (!is.null(capability) &&
      all(!vapply(negotiations, is.null, logical(1L)))) {
    .dsvert_dp_capsule_source_hash(list(
      protocol = "dsvert-biomedical-capsule-source-negotiation-set-v1",
      negotiation_hashes = unname(lapply(negotiations, `[[`, "hash"))))
  } else {
    NULL
  }
  summary <- tryCatch(
    .dsvert_dp_capsule_source_summary(
      preview$summary_json, source, context, ticket_reference),
    error = function(error) NULL)
  valid <- .dsvert_dp_has_exact_names(preview, fields) &&
    identical(
      preview$version,
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_SUMMARY_NEGOTIATION_VERSION) &&
    identical(preview$phase, "source_transport_window_attested") &&
    is.list(summary) && !is.null(capability) &&
    identical(preview$summary_sha256,
              .dsvert_dp_capsule_source_hash(summary$value)) &&
    identical(preview$ticket_negotiation_set_sha256,
              negotiation_set_sha256) &&
    identical(preview$capability, capability)
  if (!isTRUE(valid)) {
    stop("A source owner returned an invalid source transport negotiation",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(
    preview, "source-window-capability", source, context)
  summary$negotiation <- list(
    value = preview, json = response_json,
    hash = .dsvert_dp_capsule_source_hash(preview),
    capability = preview$capability)
  summary
}

.dsvert_dp_capsule_source_geometry <- function(ticket, chunk_index) {
  if (!.dsvert_dp_is_integer(chunk_index, 0, ticket$chunk_count - 1)) {
    stop("Invalid biomedical capsule source chunk index", call. = FALSE)
  }
  offset <- as.numeric(chunk_index) * as.numeric(ticket$chunk_coordinates)
  count <- min(as.numeric(ticket$chunk_coordinates),
               as.numeric(ticket$coordinate_count) - offset)
  list(index = as.numeric(chunk_index), offset = offset, count = count)
}

.dsvert_dp_capsule_source_envelope <- function(
    envelope, source, recipient, context, ticket, summary, geometry) {
  fields <- c(
    "version", "phase", "purpose", "capsule_id", "contract_hash",
    "source_transfer_id", "source_name", "source_identity_pk",
    "recipient_name", "recipient_identity_pk", "recipient_ticket_hash",
    "chunk_index", "chunk_count", "coordinate_offset",
    "coordinates_in_chunk", "chunk_coordinates", "ring_bits",
    "record_encoding", "ciphertext_bytes", "ciphertext_sha256",
    "ciphertext", "ready_for_sampling", "signature")
  ciphertext <- tryCatch(.dsvert_dp_capsule_source_decode_b64url(
    envelope$ciphertext,
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_ENVELOPE_BYTES, "ciphertext"),
    error = function(e) NULL)
  expected_ticket_hash <- ticket$hash
  valid <- .dsvert_dp_has_exact_names(envelope, fields) &&
    identical(envelope$version,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_VERSION) &&
    identical(envelope$phase, "encrypted_source_chunk_committed") &&
    identical(envelope$purpose, ticket$value$purpose) &&
    identical(envelope$capsule_id, summary$value$capsule_id) &&
    identical(envelope$contract_hash, summary$value$contract_hash) &&
    identical(envelope$source_transfer_id,
              summary$value$source_transfer_id) &&
    identical(envelope$source_name, source) &&
    identical(envelope$source_identity_pk,
              unname(context$pinset[[source]])) &&
    identical(envelope$recipient_name, recipient) &&
    identical(envelope$recipient_identity_pk,
              unname(context$pinset[[recipient]])) &&
    identical(envelope$recipient_ticket_hash, expected_ticket_hash) &&
    identical(as.numeric(envelope$chunk_index), geometry$index) &&
    identical(as.numeric(envelope$chunk_count),
              as.numeric(summary$value$chunk_count)) &&
    identical(as.numeric(envelope$coordinate_offset), geometry$offset) &&
    identical(as.numeric(envelope$coordinates_in_chunk), geometry$count) &&
    identical(as.numeric(envelope$chunk_coordinates),
              as.numeric(summary$value$chunk_coordinates)) &&
    identical(as.numeric(envelope$ring_bits), 128) &&
    identical(envelope$record_encoding,
              "little_endian_unsigned_fixed_16_bytes") &&
    identical(envelope$ready_for_sampling, FALSE) &&
    !is.null(ciphertext) && length(ciphertext) >= 60L &&
    .dsvert_dp_is_integer(
      envelope$ciphertext_bytes, 60,
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_ENVELOPE_BYTES) &&
    identical(as.numeric(envelope$ciphertext_bytes),
              as.numeric(length(ciphertext))) &&
    .dsvert_dp_capsule_source_hex(envelope$ciphertext_sha256) &&
    identical(envelope$ciphertext_sha256, digest::digest(
      ciphertext, algo = "sha256", serialize = FALSE))
  if (!isTRUE(valid)) {
    stop("A source owner returned an invalid capsule ciphertext envelope",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(
    envelope, "encrypted-chunk", source, context, ciphertext = TRUE)
  list(
    value = envelope,
    json = .dsvert_joint_dp_client_json(envelope),
    hash = .dsvert_dp_capsule_source_hash(envelope))
}

.dsvert_dp_capsule_source_bundle <- function(
    bundle_json, source, context, tickets, summary, chunk_index) {
  bundle <- .dsvert_joint_dp_client_decode(
    bundle_json, "biomedical capsule source ciphertext bundle",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_BUNDLE_BYTES)
  fields <- c(
    "version", "phase", "purpose", "capsule_id", "contract_hash",
    "source_transfer_id", "source_name", "source_identity_pk",
    "recipients", "chunk_index", "chunk_count", "coordinate_offset",
    "coordinates_in_chunk", "chunk_coordinates", "ring_bits",
    "record_encoding", "envelopes", "ready_for_sampling")
  recipients <- tryCatch(.dsvert_dp_capsule_source_string_list(
    bundle$recipients, "bundle recipient list", 2L),
    error = function(e) NULL)
  geometry <- .dsvert_dp_capsule_source_geometry(
    tickets[[1L]]$value, chunk_index)
  valid <- .dsvert_dp_has_exact_names(bundle, fields) &&
    identical(bundle$version,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_BUNDLE_VERSION) &&
    identical(bundle$phase, "encrypted_source_chunk_bundle_committed") &&
    identical(bundle$purpose, summary$value$purpose) &&
    identical(bundle$capsule_id, summary$value$capsule_id) &&
    identical(bundle$contract_hash, summary$value$contract_hash) &&
    identical(bundle$source_transfer_id,
              summary$value$source_transfer_id) &&
    identical(bundle$source_name, source) &&
    identical(bundle$source_identity_pk,
              unname(context$pinset[[source]])) &&
    !is.null(recipients) && identical(recipients, context$designated) &&
    identical(as.numeric(bundle$chunk_index), geometry$index) &&
    identical(as.numeric(bundle$chunk_count),
              as.numeric(summary$value$chunk_count)) &&
    identical(as.numeric(bundle$coordinate_offset), geometry$offset) &&
    identical(as.numeric(bundle$coordinates_in_chunk), geometry$count) &&
    identical(as.numeric(bundle$chunk_coordinates),
              as.numeric(summary$value$chunk_coordinates)) &&
    identical(as.numeric(bundle$ring_bits), 128) &&
    identical(bundle$record_encoding,
              "little_endian_unsigned_fixed_16_bytes") &&
    identical(bundle$ready_for_sampling, FALSE) &&
    is.list(bundle$envelopes) && is.null(names(bundle$envelopes)) &&
    length(bundle$envelopes) == 2L
  if (!isTRUE(valid)) {
    stop("A source owner returned an invalid capsule ciphertext bundle",
         call. = FALSE)
  }
  envelopes <- stats::setNames(lapply(seq_along(context$designated),
    function(index) {
      recipient <- context$designated[[index]]
      .dsvert_dp_capsule_source_envelope(
        bundle$envelopes[[index]], source, recipient, context,
        tickets[[recipient]], summary, geometry)
    }), context$designated)
  common <- c(
    "capsule_id", "contract_hash", "source_transfer_id", "source_name",
    "source_identity_pk", "chunk_index", "chunk_count", "coordinate_offset",
    "coordinates_in_chunk", "chunk_coordinates", "ring_bits",
    "record_encoding")
  if (!identical(envelopes[[1L]]$value[common],
                 envelopes[[2L]]$value[common])) {
    stop("The paired capsule ciphertext envelopes have conflicting geometry",
         call. = FALSE)
  }
  envelopes
}

.dsvert_dp_capsule_source_bundle_window <- function(
    response_json, source, context, tickets, summary, requested_indices) {
  response <- .dsvert_joint_dp_client_decode(
    response_json, "biomedical capsule source ciphertext window",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES)
  if (identical(response$version,
                .DSVERT_CLIENT_DP_CAPSULE_SOURCE_BUNDLE_VERSION)) {
    return(list(.dsvert_dp_capsule_source_bundle(
      response_json, source, context, tickets, summary,
      requested_indices[[1L]])))
  }
  fields <- c("version", "phase", "bundles", "ready_for_sampling")
  valid <- .dsvert_dp_has_exact_names(response, fields) &&
    identical(response$version,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION) &&
    identical(response$phase,
              "encrypted_source_chunk_window_committed") &&
    identical(response$ready_for_sampling, FALSE) &&
    is.list(response$bundles) && is.null(names(response$bundles)) &&
    length(response$bundles) >= 1L &&
    length(response$bundles) <= length(requested_indices) &&
    length(response$bundles) <=
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_CHUNKS
  if (!isTRUE(valid)) {
    stop("A source owner returned an invalid capsule ciphertext window",
         call. = FALSE)
  }
  indices <- requested_indices[seq_along(response$bundles)]
  lapply(seq_along(response$bundles), function(index) {
    bundle_json <- .dsvert_joint_dp_client_json(
      response$bundles[[index]])
    .dsvert_dp_capsule_source_bundle(
      bundle_json, source, context, tickets, summary, indices[[index]])
  })
}

.dsvert_dp_capsule_source_ack <- function(
    ack_json, recipient, source, context, ticket, summary, envelope,
    source_complete, aggregation_complete) {
  ack <- .dsvert_joint_dp_client_decode(
    ack_json, "biomedical capsule source acknowledgement",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_SUMMARY_BYTES)
  fields <- c(
    "version", "phase", "purpose", "capsule_id", "contract_hash",
    "source_transfer_id", "source_name", "source_identity_pk",
    "recipient_name", "recipient_identity_pk", "recipient_ticket_hash",
    "chunk_index", "chunk_count", "ciphertext_sha256", "source_complete",
    "capsule_aggregation_complete", "history_gate", "ready_for_sampling",
    "signature")
  value <- envelope$value
  valid <- .dsvert_dp_has_exact_names(ack, fields) &&
    identical(ack$version, .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ACK_VERSION) &&
    identical(ack$phase, "source_chunk_aggregated") &&
    identical(ack$purpose, ticket$value$purpose) &&
    identical(ack$capsule_id, summary$value$capsule_id) &&
    identical(ack$contract_hash, summary$value$contract_hash) &&
    identical(ack$source_transfer_id, summary$value$source_transfer_id) &&
    identical(ack$source_name, source) &&
    identical(ack$source_identity_pk, unname(context$pinset[[source]])) &&
    identical(ack$recipient_name, recipient) &&
    identical(ack$recipient_identity_pk,
              unname(context$pinset[[recipient]])) &&
    identical(ack$recipient_ticket_hash, ticket$hash) &&
    identical(as.numeric(ack$chunk_index), as.numeric(value$chunk_index)) &&
    identical(as.numeric(ack$chunk_count), as.numeric(value$chunk_count)) &&
    identical(ack$ciphertext_sha256, value$ciphertext_sha256) &&
    identical(ack$source_complete, source_complete) &&
    identical(ack$capsule_aggregation_complete, aggregation_complete) &&
    identical(ack$history_gate, FALSE) &&
    identical(ack$ready_for_sampling, FALSE)
  if (!isTRUE(valid)) {
    stop("A designated peer returned an invalid capsule source acknowledgement",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(
    ack, "aggregate-ack", recipient, context)
  ack
}

.dsvert_dp_capsule_source_accept_expressions <- function(
    peers, envelope_batches) {
  stats::setNames(lapply(peers, function(peer) {
    envelopes <- lapply(envelope_batches, function(batch) {
      batch[[peer]]$value
    })
    payload <- if (length(envelopes) == 1L) {
      envelope_batches[[1L]][[peer]]$json
    } else {
      .dsvert_joint_dp_client_json(list(
        version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION,
        phase = "recipient_encrypted_chunk_window",
        envelopes = unname(envelopes), ready_for_sampling = FALSE))
    }
    if (length(envelopes) == 1L) {
      call(name = "dsvertDPCapsuleSourceAcceptDS", envelope_json = payload)
    } else {
      call(
        name = "dsvertDPCapsuleSourceAcceptDS", envelope_json = payload,
        transport_contract =
          .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)
    }
  }), peers)
}

.dsvert_dp_capsule_source_accept_prefix <- function(
    peers, envelope_batches, capacity_bytes = NULL) {
  for (count in rev(seq_along(envelope_batches))) {
    expressions <- .dsvert_dp_capsule_source_accept_expressions(
      peers, envelope_batches[seq_len(count)])
    wire_expressions <- .dsvert_dsi_text_frame_expressions(expressions)
    fits <- tryCatch({
      .dsvert_validate_dsi_expression_sizes(
        wire_expressions, capacity_bytes = capacity_bytes)
      TRUE
    }, dsvert_resource_oversize = function(error) FALSE)
    if (isTRUE(fits)) {
      return(list(count = count, expressions = expressions))
    }
  }
  stop(.dsvert_client_resource_oversize(
    scope = "biomedical capsule source recipient window"))
}

.dsvert_dp_capsule_source_chunks_for_bytes <- function(bytes, maximum) {
  bytes <- suppressWarnings(as.numeric(bytes))
  maximum <- suppressWarnings(as.numeric(maximum))
  if (length(bytes) != 1L || is.na(bytes) || !is.finite(bytes) ||
      bytes != floor(bytes) || bytes < 0 || length(maximum) != 1L ||
      is.na(maximum) || !is.finite(maximum) || maximum != floor(maximum) ||
      maximum < 1) return(1L)
  usable <- min(
    bytes, maximum,
    as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_BYTES)) -
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_FRAMING_RESERVE
  per_chunk <-
    as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_BUNDLE_BYTES) + 1
  as.integer(max(1, min(
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_WINDOW_CHUNKS,
    floor(max(0, usable) / per_chunk))))
}

.dsvert_dp_capsule_source_effective_window <- function(
    source, peers, capability, transport_limits) {
  if (is.null(capability) ||
      !identical(
        .dsvert_joint_dp_client_json(capability),
        .dsvert_joint_dp_client_json(
          .dsvert_dp_capsule_source_window_capability())) ||
      !is.list(transport_limits) ||
      !source %in% names(transport_limits$response_bytes) ||
      !all(peers %in% names(transport_limits$request_payload_bytes)) ||
      !source %in% names(transport_limits$response_probe_supported) ||
      !isTRUE(transport_limits$response_probe_supported[[source]])) {
    return(1L)
  }
  response_chunks <- .dsvert_dp_capsule_source_chunks_for_bytes(
    transport_limits$response_bytes[[source]],
    capability$maximum_response_bytes)
  request_chunks <- vapply(peers, function(peer) {
    .dsvert_dp_capsule_source_chunks_for_bytes(
      transport_limits$request_payload_bytes[[peer]],
      capability$maximum_accept_bytes)
  }, integer(1L))
  as.integer(min(
    capability$maximum_window_chunks, response_chunks, request_chunks))
}

.dsvert_dp_capsule_source_ack_window <- function(
    response_json, recipient, source, context, ticket, summary,
    envelope_batches, source_index, source_count) {
  if (length(envelope_batches) == 1L) {
    envelope <- envelope_batches[[1L]][[recipient]]
    source_complete <- identical(
      as.numeric(envelope$value$chunk_index),
      as.numeric(envelope$value$chunk_count) - 1)
    return(list(.dsvert_dp_capsule_source_ack(
      response_json, recipient, source, context, ticket, summary, envelope,
      source_complete, source_complete && source_index == source_count)))
  }
  response <- .dsvert_joint_dp_client_decode(
    response_json, "biomedical capsule source acknowledgement window",
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_SUMMARY_BYTES)
  fields <- c(
    "version", "phase", "acknowledgements", "ready_for_sampling")
  valid <- .dsvert_dp_has_exact_names(response, fields) &&
    identical(response$version,
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ACK_WINDOW_VERSION) &&
    identical(response$phase, "source_chunk_window_aggregated") &&
    identical(response$ready_for_sampling, FALSE) &&
    is.list(response$acknowledgements) &&
    is.null(names(response$acknowledgements)) &&
    length(response$acknowledgements) == length(envelope_batches)
  if (!isTRUE(valid)) {
    stop("A designated peer returned an invalid source acknowledgement window",
         call. = FALSE)
  }
  lapply(seq_along(envelope_batches), function(index) {
    envelope <- envelope_batches[[index]][[recipient]]
    source_complete <- identical(
      as.numeric(envelope$value$chunk_index),
      as.numeric(envelope$value$chunk_count) - 1)
    .dsvert_dp_capsule_source_ack(
      .dsvert_joint_dp_client_json(response$acknowledgements[[index]]),
      recipient, source, context, ticket, summary, envelope,
      source_complete, source_complete && source_index == source_count)
  })
}

.dsvert_dp_capsule_source_transport_context <- function(
    manifest_json, context, .aggregate) {
  if (!is.list(context$all_conns) ||
      !identical(names(context$all_conns), context$servers) ||
      !identical(names(context$conns), context$designated)) {
    stop("Capsule source transport requires the complete pinned connection set",
         call. = FALSE)
  }
  manifest <- .dsvert_dp_capsule_source_manifest(manifest_json, context)
  peers <- context$designated
  openings <- .dsvert_dp_capsule_allocation_openings(context)
  .dsvert_maybe_negotiate_dsi_chunk_size(
    context$all_conns, .aggregate)
  transport_limits <- .dsvert_dsi_transport_site_limits(
    context$all_conns, .aggregate)

  ticket_expressions <- stats::setNames(lapply(peers, function(peer) call(
    name = "dsvertDPCapsuleSourceTicketDS",
    manifest_json = manifest_json)), peers)
  ticket_responses <- .dsvert_fanout_by_site(
    context$conns, ticket_expressions,
    operation = "biomedical capsule legacy recipient-ticket fan-out",
    .aggregate = .aggregate)
  legacy_tickets <- .dsvert_dp_capsule_source_ticket_set(
    ticket_responses, context, manifest)
  chunk_count <- as.integer(legacy_tickets[[1L]]$value$chunk_count)
  sources <- legacy_tickets[[1L]]$sources
  tickets <- legacy_tickets
  if (chunk_count >= 2L) {
    negotiation_expressions <- stats::setNames(lapply(peers, function(peer) {
      call(
        name = "dsvertDPCapsuleSourceTicketDS",
        manifest_json = manifest_json,
        transport_contract =
          .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)
    }), peers)
    negotiation_responses <- .dsvert_dp_capsule_source_optional_negotiation(
      context$conns, negotiation_expressions, .aggregate)
    response_fits <- !is.null(negotiation_responses) && all(vapply(
      negotiation_responses, function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
          nchar(value, type = "bytes") <=
            .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_TICKET_BYTES
      }, logical(1L)))
    if (isTRUE(response_fits)) {
      candidate <- .dsvert_dp_capsule_source_ticket_set(
        negotiation_responses, context, manifest)
      capability <- attr(candidate, "window_capability", exact = TRUE)
      if (!is.null(capability)) {
        inner_identical <- all(vapply(peers, function(peer) {
          identical(candidate[[peer]]$json, legacy_tickets[[peer]]$json)
        }, logical(1L)))
        if (!isTRUE(inner_identical)) {
          stop("The source transport negotiation changed its persistent ticket",
               call. = FALSE)
        }
        source_expressions <- stats::setNames(lapply(sources, function(source) {
          call(
            name = "dsvertDPCapsuleSourceTicketDS",
            manifest_json = manifest_json,
            transport_contract =
              .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CAPABILITY_ONLY_REQUEST_V2)
        }), sources)
        source_responses <- .dsvert_dp_capsule_source_optional_negotiation(
          context$all_conns[sources], source_expressions, .aggregate)
        source_response_fits <- !is.null(source_responses) && all(vapply(
          source_responses, function(value) {
            is.character(value) && length(value) == 1L && !is.na(value) &&
              nchar(value, type = "bytes") <=
                .DSVERT_CLIENT_DP_CAPSULE_SOURCE_MAX_TICKET_BYTES
          }, logical(1L)))
        if (isTRUE(source_response_fits)) {
          source_compatible <- vapply(sources, function(source) {
            .dsvert_dp_capsule_source_capability_attestation(
              source_responses[[source]], source, context, manifest,
              legacy_tickets[[1L]]$value$contract_hash)
          }, logical(1L))
          if (all(source_compatible)) tickets <- candidate
        }
      }
    }
  }
  negotiated_capability <- attr(
    tickets, "window_capability", exact = TRUE)
  ticket_payloads <- if (is.null(negotiated_capability)) {
    lapply(tickets, `[[`, "json")
  } else {
    lapply(tickets, function(ticket) ticket$negotiation$json)
  }

  prepare_expressions <- stats::setNames(lapply(sources, function(source) call(
    name = "dsvertDPCapsuleSourcePrepareDS",
    manifest_json = manifest_json,
    first_ticket_json = ticket_payloads[[peers[[1L]]]],
    second_ticket_json = ticket_payloads[[peers[[2L]]]],
    first_opening_json = openings[[peers[[1L]]]],
    second_opening_json = openings[[peers[[2L]]]])), sources)
  prepare_responses <- .dsvert_fanout_by_site(
    context$all_conns[sources], prepare_expressions,
    operation = "biomedical capsule source-owner preparation fan-out",
    .aggregate = .aggregate)
  summaries <- stats::setNames(lapply(sources, function(source) {
    .dsvert_dp_capsule_source_summary_response(
      prepare_responses[[source]], source, context, tickets[[1L]], tickets)
  }), sources)
  source_negotiations <- lapply(summaries, `[[`, "negotiation")
  window_capability <- if (!is.null(negotiated_capability) &&
      all(!vapply(source_negotiations, is.null, logical(1L))) &&
      all(vapply(source_negotiations, function(value) {
        identical(value$capability, negotiated_capability)
      }, logical(1L)))) {
    negotiated_capability
  } else {
    NULL
  }

  final_acks <- NULL
  accept_capacity <- transport_limits$expression_bytes[peers]
  if (!is.numeric(accept_capacity) || anyNA(accept_capacity) ||
      any(!is.finite(accept_capacity))) accept_capacity <- NULL
  for (source_index in seq_along(sources)) {
    source <- sources[[source_index]]
    chunk_index <- 0L
    source_window_chunks <- .dsvert_dp_capsule_source_effective_window(
      source, peers, window_capability, transport_limits)
    while (chunk_index < chunk_count) {
      requested_count <- min(
        source_window_chunks, chunk_count - chunk_index)
      requested_indices <- seq.int(
        chunk_index, length.out = requested_count)
      fetch_expression <- if (length(requested_indices) == 1L) {
        call(
          name = "dsvertDPCapsuleSourceChunkDS",
          source_transfer_id = summaries[[source]]$value$source_transfer_id,
          chunk_index = as.integer(requested_indices))
      } else {
        call(
          name = "dsvertDPCapsuleSourceChunkDS",
          source_transfer_id = summaries[[source]]$value$source_transfer_id,
          chunk_index = as.integer(requested_indices),
          transport_contract =
            .DSVERT_CLIENT_DP_CAPSULE_SOURCE_ADAPTIVE_TRANSPORT)
      }
      fetched <- .dsvert_aggregate_strict(
        context$all_conns[source], fetch_expression,
        operation = "biomedical capsule paired ciphertext fetch",
        .aggregate = .aggregate)
      envelope_batches <- .dsvert_dp_capsule_source_bundle_window(
        fetched[[source]], source, context, tickets, summaries[[source]],
        requested_indices)
      if (length(envelope_batches) < requested_count) {
        source_window_chunks <- length(envelope_batches)
      }
      offset <- 1L
      while (offset <= length(envelope_batches)) {
        remaining <- envelope_batches[seq.int(
          offset, length(envelope_batches))]
        prefix <- .dsvert_dp_capsule_source_accept_prefix(
          peers, remaining, capacity_bytes = accept_capacity)
        current <- remaining[seq_len(prefix$count)]
        accepted <- .dsvert_fanout_by_site(
          context$conns, prefix$expressions,
          operation = "biomedical capsule paired recipient acceptance",
          .aggregate = .aggregate,
          .expression_capacity_bytes = accept_capacity)
        ack_windows <- stats::setNames(lapply(peers, function(peer) {
          .dsvert_dp_capsule_source_ack_window(
            accepted[[peer]], peer, source, context, tickets[[peer]],
            summaries[[source]], current, source_index, length(sources))
        }), peers)
        final_acks <- stats::setNames(lapply(peers, function(peer) {
          tail(ack_windows[[peer]], 1L)[[1L]]
        }), peers)
        offset <- offset + prefix$count
      }
      chunk_index <- chunk_index + length(envelope_batches)
    }
  }
  if (is.null(final_acks) || !all(vapply(
        final_acks, `[[`, logical(1L), "capsule_aggregation_complete"))) {
    stop("The pinned peers did not complete the capsule source aggregation",
         call. = FALSE)
  }
  receipt <- list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_RECEIPT_VERSION,
    phase = "recipient_aggregates_committed",
    purpose = manifest$purpose,
    capsule_id = manifest$capsule_id,
    contract_hash = tickets[[1L]]$value$contract_hash,
    peer_pinset_sha256 = tickets[[1L]]$value$peer_pinset_sha256,
    source_peers = sources,
    designated_noise_peers = peers,
    coordinate_count = as.numeric(tickets[[1L]]$value$coordinate_count),
    chunk_coordinates = as.numeric(tickets[[1L]]$value$chunk_coordinates),
    chunk_count = as.numeric(tickets[[1L]]$value$chunk_count),
    ring_bits = 128,
    record_encoding = "little_endian_unsigned_fixed_16_bytes",
    durable_replay = TRUE,
    operation_or_request_limit = FALSE,
    history_can_deny_operation = FALSE,
    payload_exposed = FALSE,
    sampler_handoff_ready = !isTRUE(manifest$cross_enabled))
  if (isTRUE(manifest$cross_enabled)) {
    receipt <- c(receipt, list(
      release_coordinate_count = manifest$release_coordinate_count,
      release_coordinate_order_sha256 =
        manifest$release_coordinate_order_sha256,
      private_layout_sha256 = manifest$private_layout_sha256))
  }
  receipt
}

#' Relay one immutable biomedical capsule source into its two pinned peers
#'
#' Internal orchestration only. It performs a complete reusable-capsule status
#' handshake, relays paired authenticated Ring128 ciphertexts in canonical
#' owner/chunk order, and returns public completion metadata. It never returns
#' a coordinate, share, mask, seed, plaintext, or patient-derived digest.
#'
#' @param manifest_json Canonical immutable biomedical capsule manifest.
#' @param datasources Complete named DataSHIELD connection set.
#' @param status Optional already validated full capsule-status handshake.
#' @param allocation_openings Canonical cross-signed allocation openings,
#'   named exactly by the two designated peers.
#' @param .aggregate Injectable DSI aggregate implementation for tests.
#' @return A redacted internal sampler-handoff receipt.
#' @keywords internal
.dsvert_dp_capsule_source_transport <- function(
    manifest_json, datasources, status = NULL, allocation_openings = NULL,
    .aggregate = DSI::datashield.aggregate) {
  context <- .dsvert_joint_dp_client_context(
    datasources, status = status, .aggregate = .aggregate)
  context$allocation_openings <- allocation_openings
  .dsvert_dp_capsule_source_transport_context(
    manifest_json, context, .aggregate)
}
