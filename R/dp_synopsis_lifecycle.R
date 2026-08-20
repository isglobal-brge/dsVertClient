# Client trust establishment and durable replay for the lifetime-independent
# synopsis.  This protocol has no request, rate, catalogue, or lifetime
# entitlement.  Its status object contains only signed custodian bootstraps.

.DSVERT_CLIENT_SYNOPSIS_POLICY_VERSION <-
  "dsvert-stateless-catalog-synopsis-policy-v1"
.DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_VERSION <-
  "dsvert-stateless-catalog-synopsis-bootstrap-v1"
.DSVERT_CLIENT_SYNOPSIS_BIND_REQUEST_VERSION <-
  "dsvert-stateless-catalog-synopsis-bind-request-v1"
.DSVERT_CLIENT_SYNOPSIS_BIND_SIGNATURE_VERSION <-
  "dsvert-stateless-catalog-synopsis-bind-signature-v1"
.DSVERT_CLIENT_SYNOPSIS_BOUND_MANIFEST_VERSION <-
  "dsvert-stateless-catalog-synopsis-bound-manifest-v1"
.DSVERT_CLIENT_SYNOPSIS_PUBLICATION_VERSION <-
  "dsvert-stateless-catalog-synopsis-publication-v1"
.DSVERT_CLIENT_SYNOPSIS_FINALIZE_ACK_VERSION <-
  "dsvert-stateless-catalog-synopsis-finalize-ack-v1"
.DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_DOMAIN <-
  "dsVert/stateless-catalog-synopsis/bootstrap/v1|"
.DSVERT_CLIENT_SYNOPSIS_BIND_DOMAIN <-
  "dsVert/stateless-catalog-synopsis/bind/v1|"
.DSVERT_CLIENT_SYNOPSIS_FINALIZE_ACK_DOMAIN <-
  "dsVert/stateless-catalog-synopsis/finalize-ack/v1|"
.DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_REQUEST_VERSION <-
  "dsvert-stateless-catalog-synopsis-local-categorical-pair-request-v1"
.DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_SELECTOR_VERSION <-
  "dsvert-stateless-catalog-synopsis-local-categorical-pair-selector-v1"
.DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_SELECTOR_DOMAIN <-
  "dsVert/stateless-catalog-synopsis/local-categorical-pair-selector/v1|"
.DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION <-
  "dsvert-stateless-catalog-synopsis-cross-categorical-pair-request-v1"
.DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_SELECTOR_VERSION <-
  "dsvert-stateless-catalog-synopsis-cross-categorical-pair-selector-v1"
.DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_SELECTOR_DOMAIN <-
  "dsVert/stateless-catalog-synopsis/cross-categorical-pair-selector/v1|"
.DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION <-
  "dsvert-stateless-catalog-synopsis-categorical-pair-request-v1"

.dsvert_dp_synopsis_client_string_list_v1 <- function(
    value, what, minimum = 1L) {
  result <- if (is.character(value) && is.null(names(value))) {
    unname(value)
  } else if (is.list(value) && is.null(names(value)) &&
             all(vapply(value, .dsvert_dp_is_string, logical(1L)))) {
    unname(unlist(value, use.names = FALSE))
  } else character()
  if (length(result) < minimum || anyNA(result) || any(!nzchar(result)) ||
      anyDuplicated(result)) {
    stop("Invalid synopsis ", what, call. = FALSE)
  }
  enc2utf8(result)
}

.dsvert_dp_synopsis_client_pinset_v1 <- function(value, peers) {
  if (is.character(value) && !is.null(names(value))) {
    pins <- value
  } else if (is.list(value) && !is.null(names(value)) &&
             all(vapply(value, .dsvert_dp_is_string, logical(1L)))) {
    pins <- unlist(value, use.names = TRUE)
  } else {
    pins <- character()
  }
  if (!identical(names(pins), peers) || anyDuplicated(unname(pins)) ||
      any(vapply(unname(pins), function(value) {
        is.null(.dsvert_dp_normalize_identity_pk(value))
      }, logical(1L)))) {
    stop("The synopsis bootstrap has an invalid ordered pinset",
         call. = FALSE)
  }
  stats::setNames(vapply(unname(pins),
    .dsvert_dp_normalize_identity_pk, character(1L)), peers)
}

.dsvert_dp_synopsis_connection_identities_v1 <- function(
    conns, peers, .aggregate) {
  raw <- .dsvert_fanout_by_site(
    conns, stats::setNames(rep(list(call(name = "dsvertIdentityPkDS")),
                                  length(peers)), peers),
    operation = "synopsis connection identity fan-out",
    .aggregate = .aggregate)
  stats::setNames(vapply(peers, function(peer) {
    value <- raw[[peer]]
    normalized <- if (is.list(value) && identical(names(value), "identity_pk")) {
      .dsvert_dp_normalize_identity_pk(value$identity_pk)
    } else NULL
    if (is.null(normalized)) {
      stop("Connection '", peer, "' returned an invalid identity",
           call. = FALSE)
    }
    normalized
  }, character(1L)), peers)
}

.dsvert_dp_synopsis_client_policy_v1 <- function(value, peers) {
  fields <- c(
    "version", "privacy_scope", "global_composition_claim", "domain",
    "cohort_id", "ordered_peer_pinset", "peer_pinset_sha256",
    "peer_count", "designated_noise_peers", "artifact_epsilon",
    "artifact_delta", "adjacency", "patient_column", "unit_capacity",
    "fixed_cohort_size", "max_records_per_unit", "overflow_policy",
    "contingency_unit_aggregation_policy", "numeric_grid_bits",
    "primitive_scope", "mechanism_version", "sticky_release",
    "distinct_artifact_gate", "request_limit", "rate_limit",
    "catalog_limit")
  valid <- .dsvert_dp_has_exact_names(value, fields) &&
    identical(value$version, .DSVERT_CLIENT_SYNOPSIS_POLICY_VERSION) &&
    identical(value$privacy_scope, "per_canonical_artifact_v1") &&
    identical(value$global_composition_claim, FALSE) &&
    .dsvert_dp_is_string(value$domain) &&
    .dsvert_dp_is_string(value$cohort_id) &&
    .dsvert_dp_is_integer(value$peer_count, 2L, 4096L) &&
    identical(as.numeric(value$peer_count), as.numeric(length(peers))) &&
    .dsvert_dp_is_number(value$artifact_epsilon, 0, 8,
                         lower_open = TRUE) &&
    .dsvert_dp_is_number(value$artifact_delta, 0, 1) &&
    value$artifact_delta < 1 &&
    .dsvert_dp_is_integer(value$unit_capacity, 1L) &&
    (is.null(value$fixed_cohort_size) ||
       .dsvert_dp_is_integer(value$fixed_cohort_size, 1L)) &&
    .dsvert_dp_is_integer(value$max_records_per_unit, 1L) &&
    .dsvert_dp_is_integer(value$numeric_grid_bits, 8L, 18L) &&
    is.list(value$primitive_scope) &&
    identical(value$sticky_release, TRUE) &&
    identical(value$distinct_artifact_gate, FALSE) &&
    identical(value$request_limit, FALSE) &&
    identical(value$rate_limit, FALSE) &&
    identical(value$catalog_limit, FALSE)
  if (!isTRUE(valid)) {
    stop("A peer returned an invalid no-lifetime synopsis policy",
         call. = FALSE)
  }
  pins <- .dsvert_dp_synopsis_client_pinset_v1(
    value$ordered_peer_pinset, peers)
  if (!identical(value$peer_pinset_sha256,
                 .dsvert_vector_hash(as.list(pins)))) {
    stop("The synopsis bootstrap pinset hash is invalid", call. = FALSE)
  }
  designated <- .dsvert_dp_synopsis_client_string_list_v1(
    value$designated_noise_peers, "noise authority set", 2L)
  if (length(designated) != 2L ||
      !identical(designated, sort(designated, method = "radix")) ||
      !all(designated %in% peers)) {
    stop("The synopsis bootstrap authorities are invalid", call. = FALSE)
  }
  list(value = value, pins = pins, designated = designated)
}

.dsvert_dp_synopsis_bootstrap_draft_v1 <- function(
    draft, peer, context) {
  fields <- c(
    "version", "phase", "peer_name", "peer_identity_pk",
    "peer_pinset_sha256", "domain", "cohort_id",
    "dataset_mapping_mode", "datasets", "workload_fragments",
    "data_access", "patient_derived_metadata", "operation_limit",
    "request_limit", "history_can_deny_operation")
  reference <- context$status[[peer]]$policy
  valid <- .dsvert_dp_has_exact_names(draft, fields) &&
    identical(draft$version,
              .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_DRAFT_VERSION) &&
    identical(draft$phase, "custodian_policy_draft") &&
    identical(draft$peer_name, peer) &&
    identical(draft$peer_identity_pk, unname(context$pinset[[peer]])) &&
    identical(draft$peer_pinset_sha256,
              reference$peer_pinset_sha256) &&
    identical(draft$domain, reference$domain) &&
    identical(draft$cohort_id, reference$cohort_id) &&
    draft$dataset_mapping_mode %in% c(
      "automatic_single_local_dataset",
      "custodian_explicit_dataset_mapping_v1") &&
    identical(draft$data_access, FALSE) &&
    identical(draft$patient_derived_metadata, FALSE) &&
    identical(draft$operation_limit, FALSE) &&
    identical(draft$request_limit, FALSE) &&
    identical(draft$history_can_deny_operation, FALSE) &&
    is.list(draft$datasets) && length(draft$datasets) > 0L &&
    !is.null(names(draft$datasets)) && !anyNA(names(draft$datasets)) &&
    !anyDuplicated(names(draft$datasets))
  if (!isTRUE(valid)) {
    stop("A signed synopsis bootstrap contains an invalid draft",
         call. = FALSE)
  }
  fragments <- .dsvert_dp_capsule_manifest_fragments(
    draft$workload_fragments)
  datasets <- draft$datasets[order(names(draft$datasets), method = "radix")]
  normalized <- vector("list", length(datasets))
  names(normalized) <- names(datasets)
  seen <- character()
  for (data_name in names(datasets)) {
    dataset <- datasets[[data_name]]
    expected <- c(
      "dataset_id", "dataset_version", "schema_version",
      "alignment_group", "alignment_protocol_version", "patient_column",
      "columns")
    valid_dataset <- grepl(
      "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", data_name) &&
      .dsvert_dp_has_exact_names(dataset, expected) &&
      all(vapply(dataset[c(
        "dataset_id", "dataset_version", "schema_version",
        "alignment_group", "patient_column")], function(value) {
          .dsvert_dp_is_string(value) && grepl(
            "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)
        }, logical(1L))) &&
      identical(dataset$schema_version,
                .DSVERT_CLIENT_DP_CAPSULE_POLICY_SCHEMA_VERSION) &&
      identical(dataset$alignment_group, reference$cohort_id) &&
      .dsvert_dp_is_integer(dataset$alignment_protocol_version, 1L) &&
      is.list(dataset$columns) && length(dataset$columns) > 0L &&
      !is.null(names(dataset$columns)) && !anyNA(names(dataset$columns)) &&
      !anyDuplicated(names(dataset$columns))
    if (!isTRUE(valid_dataset)) {
      stop("A signed synopsis bootstrap contains invalid dataset metadata",
           call. = FALSE)
    }
    columns <- dataset$columns[order(names(dataset$columns), method = "radix")]
    if (any(names(columns) %in% seen)) {
      stop("Synopsis bootstrap columns are globally ambiguous",
           call. = FALSE)
    }
    seen <- c(seen, names(columns))
    columns <- lapply(columns, .dsvert_dp_capsule_manifest_column,
                      peer = peer)
    normalized[[data_name]] <- list(
      dataset_id = dataset$dataset_id,
      dataset_version = dataset$dataset_version,
      schema_version = dataset$schema_version,
      alignment_group = dataset$alignment_group,
      alignment_protocol_version =
        as.numeric(dataset$alignment_protocol_version),
      patient_column = dataset$patient_column,
      columns = columns)
  }
  list(value = draft, datasets = normalized,
       workload_fragments = fragments)
}

.dsvert_dp_synopsis_outer_signature_v1 <- function(
    value, peer, context, domain) {
  unsigned <- value[setdiff(names(value), "signature")]
  if (identical(domain, .DSVERT_CLIENT_SYNOPSIS_BIND_DOMAIN)) {
    unsigned$manifest_json <- NULL
    unsigned$artifact_commitments <- NULL
  }
  .dsvert_dp_synopsis_client_signature(
    unsigned, value$signature, unname(context$pinset[[peer]]), domain,
    "synopsis lifecycle")
}

.dsvert_dp_synopsis_bound_manifest_response_v1 <- function(
    response_json, peer, context, schema, bootstrap_set_sha256) {
  response <- .dsvert_dp_synopsis_client_json(
    response_json, paste("synopsis bound manifest from", peer))
  fields <- c(
    "version", "phase", "peer_name", "peer_identity_pk",
    "peer_pinset_sha256", "bootstrap_set_sha256", "schema_sha256",
    "workload_contract_sha256", "manifest_sha256", "manifest_bytes",
    "capsule_id", "artifact_commitment_count",
    "artifact_commitments_root", "artifact_commitments", "manifest_json",
    "privacy_scope", "global_composition_claim", "durable_memoization",
    "deterministic_replay", "data_access", "request_limit", "rate_limit",
    "catalog_limit", "signature")
  valid <- .dsvert_dp_has_exact_names(response, fields) &&
    identical(response$version,
              .DSVERT_CLIENT_SYNOPSIS_BOUND_MANIFEST_VERSION) &&
    identical(response$phase,
              "server_authoritative_synopsis_manifest_memoized") &&
    identical(response$peer_name, peer) &&
    identical(response$peer_identity_pk, unname(context$pinset[[peer]])) &&
    identical(response$peer_pinset_sha256,
              context$policy$peer_pinset_sha256) &&
    identical(response$bootstrap_set_sha256, bootstrap_set_sha256) &&
    identical(response$schema_sha256, schema$sha256) &&
    identical(response$workload_contract_sha256,
              schema$workload_contract$sha256) &&
    .dsvert_vector_hex(response$manifest_sha256) &&
    .dsvert_vector_hex(response$capsule_id) &&
    .dsvert_dp_is_integer(response$manifest_bytes, 1L,
                          .DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES) &&
    identical(as.numeric(response$manifest_bytes),
              as.numeric(nchar(response$manifest_json, type = "bytes"))) &&
    identical(response$manifest_sha256, digest::digest(
      response$manifest_json, algo = "sha256", serialize = FALSE)) &&
    .dsvert_dp_is_integer(response$artifact_commitment_count, 0L,
                          .DSVERT_DP_MAX_COORDINATES) &&
    .dsvert_vector_hex(response$artifact_commitments_root) &&
    is.list(response$artifact_commitments) &&
    identical(response$privacy_scope, "per_canonical_artifact_v1") &&
    identical(response$global_composition_claim, FALSE) &&
    identical(response$durable_memoization, TRUE) &&
    identical(response$deterministic_replay, TRUE) &&
    identical(response$data_access, FALSE) &&
    identical(response$request_limit, FALSE) &&
    identical(response$rate_limit, FALSE) &&
    identical(response$catalog_limit, FALSE)
  if (!isTRUE(valid)) {
    stop("A peer returned an invalid bound synopsis manifest", call. = FALSE)
  }
  .dsvert_dp_synopsis_outer_signature_v1(
    response, peer, context, .DSVERT_CLIENT_SYNOPSIS_BIND_DOMAIN)
  response
}

.dsvert_dp_synopsis_manifest_policy_scope_v1 <- function(manifest) {
  scope <- tryCatch(
    manifest$workload$primitive_scope, error = function(error) NULL)
  if (is.list(scope) && identical(scope$mode, "all_schema")) {
    return(list(mode = "all_schema"))
  }
  explicit <- tryCatch(
    scope$selection$explicit_catalog, error = function(error) NULL)
  fields <- c(
    "numeric_moments", "categorical_marginals",
    "categorical_pairs", "correlations")
  if (!is.list(scope) || !identical(scope$mode, "catalog_v1") ||
      !is.list(explicit) || !setequal(names(explicit), fields)) {
    stop("The bound synopsis has an invalid projected policy scope",
         call. = FALSE)
  }
  string_array <- function(value) {
    if (is.null(value)) return(character())
    if (is.character(value)) return(unname(value))
    if (is.list(value) && !length(value)) return(character())
    if (is.list(value) && is.null(names(value)) &&
        all(vapply(value, .dsvert_dp_is_string, logical(1L)))) {
      return(unname(unlist(value, use.names = FALSE)))
    }
    value
  }
  explicit$numeric_moments <- string_array(explicit$numeric_moments)
  explicit$categorical_marginals <-
    string_array(explicit$categorical_marginals)
  explicit$categorical_pairs <- unname(lapply(
    unname(explicit$categorical_pairs), string_array))
  explicit$correlations <- unname(lapply(
    unname(explicit$correlations), string_array))
  .dsvert_joint_dp_client_canonical(c(
    list(mode = "catalog_v1"), explicit[fields]))
}

.dsvert_dp_synopsis_manifest_value_v1 <- function(
    manifest_json, context, schema) {
  manifest <- .dsvert_joint_dp_client_decode(
    manifest_json, "no-lifetime synopsis manifest",
    .DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES)
  fields <- c(
    "version", "logical_snapshot", "capsule_schema", "admission", "bounds",
    "workload", "capsule_identity", "execution_state")
  identity <- manifest$capsule_identity
  contract <- if (is.list(identity)) identity$contract else NULL
  projected_policy <- context$policy
  projected_policy$primitive_scope <-
    .dsvert_dp_synopsis_manifest_policy_scope_v1(manifest)
  policy_hash <- .dsvert_vector_hash(projected_policy)
  valid <- .dsvert_dp_has_exact_names(manifest, fields) &&
    identical(manifest$version,
              .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION) &&
    identical(manifest$capsule_schema,
              .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION) &&
    identical(manifest$execution_state,
              .DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE) &&
    .dsvert_dp_has_exact_names(identity, c("capsule_id", "contract")) &&
    is.list(contract) &&
    identical(contract$protocol,
              "dsvert-stateless-catalog-synopsis-capsule-identity-v1") &&
    identical(identity$capsule_id, .dsvert_vector_hash(contract)) &&
    identical(contract$consortium_id, paste0("dpsc1_", policy_hash)) &&
    identical(contract$policy_contract_hash, policy_hash) &&
    identical(contract$peer_pinset_sha256,
              context$policy$peer_pinset_sha256) &&
    identical(contract$logical_snapshot, manifest$logical_snapshot) &&
    identical(contract$capsule_schema, manifest$capsule_schema) &&
    identical(contract$admission, manifest$admission) &&
    identical(contract$bounds, manifest$bounds) &&
    identical(contract$workload, manifest$workload) &&
    identical(contract$privacy_epoch_scope,
              "per_canonical_artifact_sticky_v1") &&
    identical(.dsvert_joint_dp_client_canonical(manifest$logical_snapshot),
              schema$logical_snapshot) &&
    is.list(manifest$workload) &&
    identical(manifest$workload$workload_version,
              .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_VERSION) &&
    identical(manifest$workload$execution_state,
              .DSVERT_CLIENT_DP_CAPSULE_EXECUTION_STATE) &&
    identical(manifest$workload$registered_release_lifecycle,
              .dsvert_dp_capsule_client_registered_release_lifecycle()) &&
    identical(manifest$workload$declared_workload_fully_materialized, TRUE) &&
    identical(manifest$workload$package_family_coverage_complete, FALSE)
  if (!isTRUE(valid)) {
    stop("The bound synopsis manifest is invalid or misbound", call. = FALSE)
  }
  attestation <- manifest$workload$schema_attestation
  signers <- tryCatch(.dsvert_dp_capsule_manifest_string_array(
    attestation$signers, "synopsis schema signers"),
    error = function(error) character())
  if (!is.list(attestation) ||
      !identical(attestation$manifest_sha256, schema$sha256) ||
      !identical(.dsvert_joint_dp_client_canonical(attestation$signatures),
                 .dsvert_joint_dp_client_canonical(schema$signatures)) ||
      !identical(signers, context$servers)) {
    stop("The synopsis manifest schema attestation is invalid", call. = FALSE)
  }
  manifest
}

.dsvert_dp_synopsis_local_pair_id_v1 <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)) {
    stop("Invalid ", what, call. = FALSE)
  }
  enc2utf8(value)
}

.dsvert_dp_synopsis_local_pair_request_v1 <- function(value) {
  fields <- c("version", "family", "dataset", "columns", "owner_peer")
  if (!is.list(value) || !.dsvert_dp_has_exact_names(value, fields) ||
      !identical(
        value$version, .DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_REQUEST_VERSION) ||
      !identical(value$family, "categorical_pairs")) {
    stop("Invalid local categorical Synopsis projection request",
         call. = FALSE)
  }
  columns <- value$columns
  if (is.list(columns) && is.null(names(columns)) &&
      all(vapply(columns, .dsvert_dp_is_string, logical(1L)))) {
    columns <- unname(unlist(columns, use.names = FALSE))
  }
  if (!is.character(columns) || length(columns) != 2L || anyNA(columns) ||
      !is.null(names(columns))) {
    stop("Invalid local categorical Synopsis projection columns",
         call. = FALSE)
  }
  columns <- sort(vapply(
    columns, .dsvert_dp_synopsis_local_pair_id_v1, character(1L),
    what = "local categorical Synopsis column"), method = "radix")
  if (anyDuplicated(columns)) {
    stop("A local categorical Synopsis projection needs two columns",
         call. = FALSE)
  }
  owner <- value$owner_peer
  if (!is.null(owner)) owner <- .dsvert_dp_synopsis_local_pair_id_v1(
    owner, "local categorical Synopsis owner")
  list(
    version = .DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_REQUEST_VERSION,
    family = "categorical_pairs",
    dataset = .dsvert_dp_synopsis_local_pair_id_v1(
      value$dataset, "local categorical Synopsis dataset"),
    columns = unname(columns), owner_peer = owner)
}

.dsvert_dp_synopsis_local_pair_reference_v1 <- function(reference, owner) {
  if (!is.character(reference) || length(reference) != 1L ||
      is.na(reference)) {
    stop("Invalid signed categorical column reference", call. = FALSE)
  }
  reference <- enc2utf8(reference)
  prefix <- paste0(owner, "$")
  if (startsWith(reference, prefix)) {
    reference <- substring(reference, nchar(prefix, type = "chars") + 1L)
  } else if (grepl("$", reference, fixed = TRUE)) {
    stop("The signed categorical column reference has another owner",
         call. = FALSE)
  }
  .dsvert_dp_synopsis_local_pair_id_v1(
    reference, "signed categorical column reference")
}

.dsvert_dp_synopsis_local_pair_selector_hash_v1 <- function(value) {
  .dsvert_dp_synopsis_client_hash(
    .DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_SELECTOR_DOMAIN, value)
}

.dsvert_dp_synopsis_local_pair_project_components_v1 <- function(
    schema, dataset_name, owner_peer, references, context) {
  if (!is.list(schema) || !is.list(schema$unsigned) ||
      !is.list(schema$unsigned$datasets) ||
      !dataset_name %in% names(schema$unsigned$datasets) ||
      !is.character(references) || length(references) != 2L ||
      anyNA(references) || anyDuplicated(references)) {
    stop("Invalid local categorical Synopsis projection source",
         call. = FALSE)
  }
  references <- sort(unname(references), method = "radix")
  source <- schema$unsigned$datasets[[dataset_name]]
  columns <- source$columns[references]
  valid_columns <- length(columns) == 2L &&
    all(vapply(columns, function(column) {
      is.list(column) && identical(column$kind, "categorical") &&
        identical(column$owner_peer, owner_peer) &&
        is.atomic(column$levels) && length(column$levels) > 0L
    }, logical(1L)))
  if (!isTRUE(valid_columns) ||
      !owner_peer %in% names(source$patient_keys)) {
    stop("Invalid local categorical Synopsis projection source",
         call. = FALSE)
  }
  columns <- lapply(columns, function(column) {
    column$levels <- .dsvert_dp_capsule_manifest_string_array(
      as.list(column$levels), "local categorical projection levels")
    .dsvert_joint_dp_client_canonical(column)
  })
  projected_dataset <- source[c(
    "dataset_id", "dataset_version", "schema_version", "alignment_group")]
  projected_dataset$patient_keys <- source$patient_keys[owner_peer]
  projected_dataset$columns <- columns
  datasets <- stats::setNames(list(projected_dataset), dataset_name)
  families <- c("describe", "survival", "gaussian", "vertical_cross")
  workload <- .dsvert_joint_dp_client_canonical(c(
    list(version = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_CONTRACT_VERSION),
    stats::setNames(rep(list(list()), length(families)), families)))
  parent_snapshot <- schema$logical_snapshot
  alignment <- tryCatch(
    as.integer(parent_snapshot$alignment_protocol_version),
    error = function(error) NA_integer_)
  if (!is.list(parent_snapshot) ||
      !identical(parent_snapshot$logical_snapshot_id,
                 context$policy$cohort_id) ||
      length(alignment) != 1L || is.na(alignment) || alignment < 1L ||
      !identical(as.numeric(alignment),
                 as.numeric(parent_snapshot$alignment_protocol_version))) {
    stop("Invalid local categorical Synopsis projection snapshot",
         call. = FALSE)
  }
  fingerprint <- .dsvert_dp_capsule_manifest_hash(list(
    protocol = "dsvert-biomedical-capsule-logical-snapshot-v1",
    domain = context$policy$domain,
    cohort_id = context$policy$cohort_id,
    peer_pinset_sha256 = context$policy$peer_pinset_sha256,
    alignment_protocol_version = alignment,
    datasets = datasets, workload_contract = workload))
  logical_snapshot <- .dsvert_joint_dp_client_canonical(list(
    logical_snapshot_id = context$policy$cohort_id,
    version = paste0("schema-v1-", fingerprint),
    alignment_protocol_version = alignment))
  projected_schema <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = logical_snapshot,
    peer_pinset_sha256 = context$policy$peer_pinset_sha256,
    datasets = datasets))
  scope <- .dsvert_joint_dp_client_canonical(list(
    mode = "catalog_v1", numeric_moments = character(),
    categorical_marginals = character(),
    categorical_pairs = list(references), correlations = list()))
  projected_policy <- context$policy
  projected_policy$primitive_scope <- scope
  schema_json <- .dsvert_joint_dp_client_json(projected_schema)
  workload_json <- .dsvert_joint_dp_client_json(workload)
  list(
    unsigned = projected_schema, json = schema_json,
    sha256 = digest::digest(
      schema_json, algo = "sha256", serialize = FALSE),
    logical_snapshot = logical_snapshot,
    workload_contract = list(
      value = workload, json = workload_json,
      sha256 = digest::digest(
        workload_json, algo = "sha256", serialize = FALSE)),
    primitive_scope = scope,
    policy_sha256 = .dsvert_vector_hash(projected_policy))
}

.dsvert_dp_synopsis_local_pair_resolve_v1 <- function(
    request, schema, context) {
  request <- .dsvert_dp_synopsis_local_pair_request_v1(request)
  datasets <- if (is.list(schema) && is.list(schema$unsigned)) {
    schema$unsigned$datasets
  } else NULL
  if (!is.list(datasets) || !request$dataset %in% names(datasets)) {
    stop("The signed Synopsis schema has no requested categorical dataset",
         call. = FALSE)
  }
  dataset <- datasets[[request$dataset]]
  qualified <- lapply(names(dataset$columns), function(reference) {
    column <- dataset$columns[[reference]]
    owner <- if (is.list(column)) column$owner_peer else NULL
    list(
      reference = reference,
      column = if (.dsvert_dp_is_string(owner)) {
        .dsvert_dp_synopsis_local_pair_reference_v1(reference, owner)
      } else "",
      kind = if (is.list(column)) column$kind else NULL,
      owner_peer = owner,
      levels = if (is.list(column)) column$levels else NULL)
  })
  names(qualified) <- names(dataset$columns)
  requested <- request$columns
  candidates <- lapply(requested, function(column) {
    names(qualified)[vapply(qualified, function(candidate) {
      identical(candidate$kind, "categorical") &&
        identical(candidate$column, column) &&
        (is.null(request$owner_peer) ||
           identical(candidate$owner_peer, request$owner_peer))
    }, logical(1L))]
  })
  if (any(lengths(candidates) == 0L)) {
    stop("The requested local categorical pair is absent from signed metadata",
         call. = FALSE)
  }
  combinations <- expand.grid(
    left = candidates[[1L]], right = candidates[[2L]],
    stringsAsFactors = FALSE)
  valid <- vapply(seq_len(nrow(combinations)), function(index) {
    left <- qualified[[combinations$left[[index]]]]
    right <- qualified[[combinations$right[[index]]]]
    !identical(combinations$left[[index]], combinations$right[[index]]) &&
      identical(left$owner_peer, right$owner_peer)
  }, logical(1L))
  combinations <- combinations[valid, , drop = FALSE]
  if (nrow(combinations) != 1L) {
    stop("The signed local categorical pair is missing or ambiguous",
         call. = FALSE)
  }
  references <- sort(unname(c(
    combinations$left[[1L]], combinations$right[[1L]])), method = "radix")
  owner <- qualified[[references[[1L]]]]$owner_peer
  scope <- context$policy$primitive_scope
  if (!is.list(scope) || !.dsvert_dp_is_string(scope$mode) ||
      !scope$mode %in% c("all_schema", "catalog_v1")) {
    stop("The signed Synopsis primitive scope is invalid", call. = FALSE)
  }
  if (identical(scope$mode, "catalog_v1")) {
    pairs <- scope$categorical_pairs
    if (is.null(pairs)) pairs <- list()
    if (!is.list(pairs) || !is.null(names(pairs))) {
      stop("The signed categorical-pair catalog is invalid", call. = FALSE)
    }
    keys <- vapply(pairs, function(pair) {
      if (is.list(pair) && is.null(names(pair))) {
        pair <- unname(unlist(pair, use.names = FALSE))
      }
      if (!is.character(pair) || length(pair) != 2L || anyNA(pair)) return("")
      .dsvert_joint_dp_client_json(as.list(sort(pair, method = "radix")))
    }, character(1L))
    if (!.dsvert_joint_dp_client_json(as.list(references)) %in% keys) {
      stop("The requested categorical pair is not in the signed catalog",
           call. = FALSE)
    }
  }
  columns <- lapply(references, function(reference) {
    column <- qualified[[reference]]
    levels <- .dsvert_dp_capsule_manifest_string_array(
      as.list(column$levels), "local categorical levels")
    list(
      reference = reference, column = column$column,
      levels_sha256 = .dsvert_vector_hash(as.list(levels)))
  })
  projection <- .dsvert_dp_synopsis_local_pair_project_components_v1(
    schema, request$dataset, owner, references, context)
  parent <- list(
    schema_sha256 = projection$sha256,
    workload_contract_sha256 =
      projection$workload_contract$sha256,
    logical_snapshot = projection$logical_snapshot,
    policy_sha256 = projection$policy_sha256)
  unsigned <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_SELECTOR_VERSION,
    family = "categorical_pairs", dataset = request$dataset,
    owner_peer = owner, columns = columns, parent = parent))
  .dsvert_joint_dp_client_canonical(c(unsigned, list(
    sha256 = .dsvert_dp_synopsis_local_pair_selector_hash_v1(unsigned))))
}

.dsvert_dp_synopsis_local_pair_project_v1 <- function(
    schema, selector, context) {
  physical <- tryCatch(vapply(
    selector$columns, `[[`, character(1L), "column"),
    error = function(error) character())
  expected <- if (length(physical) == 2L) tryCatch(
    .dsvert_dp_synopsis_local_pair_resolve_v1(list(
      version = .DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_REQUEST_VERSION,
      family = "categorical_pairs", dataset = selector$dataset,
      columns = unname(physical), owner_peer = selector$owner_peer),
    schema, context), error = function(error) NULL) else NULL
  if (is.null(expected) || !identical(
      .dsvert_joint_dp_client_json(selector),
      .dsvert_joint_dp_client_json(expected))) {
    stop("The categorical Synopsis selector is detached from signed metadata",
         call. = FALSE)
  }
  references <- vapply(
    expected$columns, `[[`, character(1L), "reference")
  projection <- .dsvert_dp_synopsis_local_pair_project_components_v1(
    schema, expected$dataset, expected$owner_peer, references, context)
  projection$primitive_scope <- NULL
  projection$policy_sha256 <- NULL
  projection
}

.dsvert_dp_synopsis_local_pair_manifest_v1 <- function(manifest, selector) {
  references <- vapply(
    selector$columns, `[[`, character(1L), "reference")
  scope <- tryCatch(manifest$workload$primitive_scope,
                    error = function(error) NULL)
  explicit <- tryCatch(
    scope$selection$explicit_catalog$categorical_pairs,
    error = function(error) NULL)
  pairs <- tryCatch(
    manifest$workload$families$categorical_pairs,
    error = function(error) NULL)
  sets <- if (is.list(pairs)) pairs$sets else NULL
  cross <- if (is.list(pairs)) pairs$cross_artifacts else NULL
  families <- manifest$workload$families
  empty <- function(value) is.list(value) && length(value) == 0L
  non_pair_empty <- is.list(families) &&
    empty(families$numeric_moments$artifacts) &&
    empty(families$numeric_pair_moments$artifacts) &&
    empty(families$gaussian_models$artifacts) &&
    empty(families$fixed_numeric_histograms$artifacts) &&
    empty(families$correlation_artifacts) &&
    empty(families$describe_artifacts) && empty(families$survival_artifacts)
  set <- if (is.list(sets) && length(sets) == 1L) sets[[1L]] else NULL
  included <- if (is.list(set) && is.list(set$included_pairs) &&
      length(set$included_pairs) == 1L) {
    unname(unlist(set$included_pairs[[1L]], use.names = FALSE))
  } else character()
  valid <- identical(
    .dsvert_joint_dp_client_json(manifest$logical_snapshot),
    .dsvert_joint_dp_client_json(selector$parent$logical_snapshot)) &&
    is.list(scope) && identical(scope$mode, "catalog_v1") &&
    is.list(explicit) && length(explicit) == 1L &&
    identical(sort(unname(unlist(explicit[[1L]], use.names = FALSE)),
                   method = "radix"), references) &&
    isTRUE(non_pair_empty) && is.list(sets) && length(sets) == 1L &&
    empty(cross) && is.list(set) &&
    identical(set$dataset, selector$dataset) &&
    identical(set$owner_peer, selector$owner_peer) &&
    identical(sort(included, method = "radix"), references)
  if (!isTRUE(valid)) {
    stop("The signed local categorical Synopsis projection is invalid",
         call. = FALSE)
  }
  invisible(manifest)
}

.dsvert_dp_synopsis_policy_view_v1 <- function(policy) {
  list(
    domain = policy$domain, cohort_id = policy$cohort_id,
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    designated_noise_peers = unlist(
      policy$designated_noise_peers, use.names = FALSE),
    capsule_epsilon = policy$artifact_epsilon,
    capsule_delta = policy$artifact_delta,
    adjacency = policy$adjacency, unit_capacity = policy$unit_capacity,
    max_records_per_unit = policy$max_records_per_unit,
    overflow_policy = policy$overflow_policy)
}

.dsvert_dp_synopsis_bootstrap_build_v1 <- function(
    datasources, status = NULL, local_projection = NULL,
    .aggregate = DSI::datashield.aggregate) {
  if (!is.list(datasources) || length(datasources) < 2L ||
      is.null(names(datasources)) || anyNA(names(datasources)) ||
      any(!nzchar(names(datasources))) || anyDuplicated(names(datasources))) {
    stop("The synopsis requires the complete named datasource set",
         call. = FALSE)
  }
  peers <- sort(names(datasources), method = "radix")
  conns <- datasources[peers]
  if (!is.null(status)) {
    if (inherits(status, "dsvert_synopsis_bootstrap_v1") &&
        identical(status$context$servers, peers)) {
      if (!is.null(local_projection)) {
        if (!isTRUE(.dsvert_dp_synopsis_projection_request_compatible_v1(
              status$local_projection, local_projection))) {
          stop("The saved synopsis bootstrap targets another projection",
               call. = FALSE)
        }
      } else if (!is.null(status$local_projection)) {
        stop("A projected synopsis bootstrap requires its selector",
             call. = FALSE)
      }
      identities <- .dsvert_dp_synopsis_connection_identities_v1(
        conns, peers, .aggregate)
      if (!identical(identities, status$context$pinset)) {
        stop("The saved synopsis bootstrap does not match these connections",
             call. = FALSE)
      }
      status$context$all_conns <- conns
      status$context$conns <- status$context$all_conns[
        status$context$designated]
      status$manifest_bundle$context <- status$context
      .dsvert_dp_synopsis_trusted_bundle_v1(
        status$manifest_bundle, status$status)
      return(status)
    }
    stop("Legacy DP status is not valid for the no-lifetime synopsis",
         call. = FALSE)
  }
  identities <- .dsvert_dp_synopsis_connection_identities_v1(
    conns, peers, .aggregate)

  bootstrap_raw <- .dsvert_fanout_by_site(
    conns, stats::setNames(rep(list(call(
      name = "dsvertDPSynopsisBootstrapDS")), length(peers)), peers),
    operation = "no-lifetime synopsis bootstrap fan-out",
    .aggregate = .aggregate)
  bootstraps <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_dp_synopsis_client_json(
      bootstrap_raw[[peer]], paste("synopsis bootstrap from", peer))
    fields <- c(
      "version", "phase", "peer_name", "peer_identity_pk", "policy",
      "draft", "data_access", "patient_derived_metadata", "request_limit",
      "rate_limit", "catalog_limit", "signature")
    if (!.dsvert_dp_has_exact_names(value, fields) ||
        !identical(value$version,
                   .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_VERSION) ||
        !identical(value$phase, "custodian_policy_bootstrap") ||
        !identical(value$peer_name, peer) ||
        !identical(value$peer_identity_pk, identities[[peer]]) ||
        !identical(value$data_access, FALSE) ||
        !identical(value$patient_derived_metadata, FALSE) ||
        !identical(value$request_limit, FALSE) ||
        !identical(value$rate_limit, FALSE) ||
        !identical(value$catalog_limit, FALSE)) {
      stop("A connection returned an invalid synopsis bootstrap",
           call. = FALSE)
    }
    value
  }), peers)
  policy_info <- lapply(bootstraps, function(value) {
    .dsvert_dp_synopsis_client_policy_v1(value$policy, peers)
  })
  reference_policy <- policy_info[[1L]]
  if (!all(vapply(seq_along(policy_info), function(index) {
        identical(.dsvert_joint_dp_client_json(policy_info[[index]]$value),
                  .dsvert_joint_dp_client_json(reference_policy$value))
      }, logical(1L))) ||
      !identical(identities, reference_policy$pins)) {
    stop("Connected peers disagree on the signed synopsis federation",
         call. = FALSE)
  }
  status_value <- bootstraps
  class(status_value) <- c("ds.vertDPSynopsisStatus", "list")
  context <- list(
    status = status_value, servers = peers, pinset = reference_policy$pins,
    designated = reference_policy$designated,
    conns = conns[reference_policy$designated], all_conns = conns,
    policy = reference_policy$value)
  for (peer in peers) {
    .dsvert_dp_synopsis_outer_signature_v1(
      bootstraps[[peer]], peer, context,
      .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_DOMAIN)
  }
  drafts <- stats::setNames(lapply(peers, function(peer) {
    .dsvert_dp_synopsis_bootstrap_draft_v1(
      bootstraps[[peer]]$draft, peer, context)
  }), peers)
  authorization_schema <- .dsvert_dp_capsule_manifest_schema(
    drafts, context)
  local_projection <- if (is.null(local_projection)) NULL else
    .dsvert_dp_synopsis_projection_resolve_v1(
      local_projection, authorization_schema, context)
  schema <- if (is.null(local_projection)) {
    authorization_schema
  } else {
    .dsvert_dp_synopsis_projection_project_v1(
      authorization_schema, local_projection, context)
  }
  bootstrap_set_sha256 <- .dsvert_vector_hash(list(
    version = .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_VERSION,
    bootstraps = unname(bootstraps)))
  bind_request <- function(phase, signatures) {
    request <- list(
      version = .DSVERT_CLIENT_SYNOPSIS_BIND_REQUEST_VERSION,
      phase = phase, bootstraps = bootstraps,
      schema_signatures = signatures)
    if (!is.null(local_projection)) {
      request$local_projection <- local_projection
    }
    .dsvert_joint_dp_client_json(request)
  }
  schema_request <- bind_request("schema_signature", NULL)
  schema_raw <- .dsvert_fanout_by_site(
    conns, stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertDPSynopsisBindDS",
      bootstrap_set_json = schema_request)), peers),
    operation = "synopsis global schema binding", .aggregate = .aggregate)
  schema_responses <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_dp_synopsis_client_json(
      schema_raw[[peer]], paste("synopsis schema binding from", peer))
    fields <- c(
      "version", "phase", "peer_name", "peer_identity_pk",
      "peer_pinset_sha256", "bootstrap_set_sha256", "schema_sha256",
      "workload_contract_sha256", "logical_snapshot", "schema_signature",
      "data_access", "request_limit", "rate_limit", "catalog_limit",
      "signature")
    valid <- .dsvert_dp_has_exact_names(value, fields) &&
      identical(value$version,
                .DSVERT_CLIENT_SYNOPSIS_BIND_SIGNATURE_VERSION) &&
      identical(value$phase, "global_schema_verified") &&
      identical(value$peer_name, peer) &&
      identical(value$peer_identity_pk, unname(context$pinset[[peer]])) &&
      identical(value$peer_pinset_sha256,
                context$policy$peer_pinset_sha256) &&
      identical(value$bootstrap_set_sha256, bootstrap_set_sha256) &&
      identical(value$schema_sha256, schema$sha256) &&
      identical(value$workload_contract_sha256,
                schema$workload_contract$sha256) &&
      identical(.dsvert_joint_dp_client_canonical(value$logical_snapshot),
                schema$logical_snapshot) &&
      identical(value$data_access, FALSE) &&
      identical(value$request_limit, FALSE) &&
      identical(value$rate_limit, FALSE) &&
      identical(value$catalog_limit, FALSE)
    if (!isTRUE(valid)) {
      stop("A peer returned an invalid synopsis schema binding",
           call. = FALSE)
    }
    .dsvert_dp_capsule_schema_verify(
      schema$unsigned, value$schema_signature, peer, context)
    .dsvert_dp_synopsis_outer_signature_v1(
      value, peer, context, .DSVERT_CLIENT_SYNOPSIS_BIND_DOMAIN)
    value
  }), peers)
  signatures <- stats::setNames(lapply(
    schema_responses, `[[`, "schema_signature"), peers)
  schema$signatures <- signatures
  manifest_request <- bind_request("manifest_build", signatures)
  manifest_raw <- .dsvert_fanout_by_site(
    conns, stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertDPSynopsisBindDS",
      bootstrap_set_json = manifest_request)), peers),
    operation = "synopsis authoritative manifest binding",
    .aggregate = .aggregate)
  builds <- stats::setNames(lapply(peers, function(peer) {
    .dsvert_dp_synopsis_bound_manifest_response_v1(
      manifest_raw[[peer]], peer, context, schema, bootstrap_set_sha256)
  }), peers)
  common <- setdiff(names(builds[[1L]]),
                    c("peer_name", "peer_identity_pk", "signature"))
  if (!all(vapply(builds, function(value) {
        identical(.dsvert_joint_dp_client_json(value[common]),
                  .dsvert_joint_dp_client_json(builds[[1L]][common]))
      }, logical(1L)))) {
    stop("Pinned peers built different synopsis manifests", call. = FALSE)
  }
  signed_schema <- .dsvert_joint_dp_client_canonical(c(
    schema$unsigned, list(signatures = signatures)))
  manifest <- .dsvert_dp_synopsis_manifest_value_v1(
    builds[[1L]]$manifest_json, context, schema)
  if (!is.null(local_projection)) {
    if (.dsvert_dp_synopsis_projection_is_cross_v1(local_projection)) {
      .dsvert_dp_synopsis_cross_pair_manifest_v1(
        manifest, local_projection)
    } else {
      .dsvert_dp_synopsis_local_pair_manifest_v1(
        manifest, local_projection)
    }
  }
  artifact_index <- .dsvert_dp_capsule_artifact_commitment_index_client(
    manifest, .dsvert_dp_synopsis_policy_view_v1(context$policy),
    builds[[1L]]$manifest_sha256)
  if (!identical(builds[[1L]]$capsule_id,
                 manifest$capsule_identity$capsule_id) ||
      !identical(.dsvert_joint_dp_client_canonical(
        builds[[1L]]$artifact_commitments), artifact_index$value) ||
      !identical(as.numeric(builds[[1L]]$artifact_commitment_count),
                 as.numeric(artifact_index$count)) ||
      !identical(builds[[1L]]$artifact_commitments_root,
                 artifact_index$root)) {
    stop("The bound synopsis artifact index is invalid", call. = FALSE)
  }
  manifest_bundle <- list(
    schema_json = .dsvert_joint_dp_client_json(signed_schema),
    schema_sha256 = schema$sha256,
    workload_contract_json = schema$workload_contract$json,
    workload_contract_sha256 = schema$workload_contract$sha256,
    logical_snapshot = schema$logical_snapshot,
    manifest_json = builds[[1L]]$manifest_json,
    manifest_sha256 = builds[[1L]]$manifest_sha256,
    capsule_id = builds[[1L]]$capsule_id,
    artifact_commitments = artifact_index$value,
    artifact_commitment_count = artifact_index$count,
    artifact_commitments_root = artifact_index$root,
    manifest_build_receipts = stats::setNames(lapply(peers, function(peer) {
      builds[[peer]][setdiff(names(builds[[peer]]),
                            c("manifest_json", "artifact_commitments"))]
    }), peers),
    manifest_signatures = stats::setNames(lapply(
      builds, `[[`, "signature"), peers),
    deterministic_replay = TRUE,
    operation_or_request_limit = FALSE,
    history_can_deny_operation = FALSE,
    context = context)
  result <- list(
    version = "dsvert-stateless-catalog-synopsis-client-bootstrap-v1",
    status = status_value, manifest_bundle = manifest_bundle,
    context = context,
    layout = .dsvert_dp_capsule_vector_layout(manifest))
  if (!is.null(local_projection)) {
    result$local_projection <- local_projection
  }
  class(result) <- c("dsvert_synopsis_bootstrap_v1", "list")
  result
}

.dsvert_dp_synopsis_publication_set_v1 <- function(raw, bootstrap) {
  context <- bootstrap$context
  authorities <- context$designated
  fields <- c(
    "version", "published", "manifest_sha256", "artifact", "artifact_key",
    "compile_receipts", "receipt_set_sha256", "public_chunk_count",
    "release_receipt", "durable_publication", "session_required",
    "source_store_read", "request_limit", "rate_limit", "catalog_limit")
  publications <- stats::setNames(lapply(authorities, function(peer) {
    value <- .dsvert_dp_synopsis_client_json(
      raw[[peer]], paste("durable synopsis publication from", peer))
    valid <- .dsvert_dp_has_exact_names(value, fields) &&
      identical(value$version,
                .DSVERT_CLIENT_SYNOPSIS_PUBLICATION_VERSION) &&
      identical(value$published, TRUE) &&
      identical(value$manifest_sha256,
                bootstrap$manifest_bundle$manifest_sha256) &&
      is.list(value$artifact) && .dsvert_vector_hex(value$artifact_key) &&
      identical(value$artifact$artifact_key, value$artifact_key) &&
      is.list(value$compile_receipts) &&
      .dsvert_vector_hex(value$receipt_set_sha256) &&
      .dsvert_dp_is_integer(value$public_chunk_count, 1L,
                            .DSVERT_DP_MAX_COORDINATES) &&
      is.list(value$release_receipt) &&
      identical(value$durable_publication, TRUE) &&
      identical(value$session_required, FALSE) &&
      identical(value$source_store_read, FALSE) &&
      identical(value$request_limit, FALSE) &&
      identical(value$rate_limit, FALSE) &&
      identical(value$catalog_limit, FALSE)
    if (!isTRUE(valid)) {
      stop("A noise authority returned an invalid durable publication",
           call. = FALSE)
    }
    value
  }), authorities)
  common <- setdiff(fields, "release_receipt")
  if (!identical(
        .dsvert_joint_dp_client_json(publications[[1L]][common]),
        .dsvert_joint_dp_client_json(publications[[2L]][common]))) {
    stop("The two durable synopsis publications disagree", call. = FALSE)
  }
  receipts <- publications[[1L]]$compile_receipts
  receipt_peers <- vapply(receipts, function(value) {
    if (is.list(value) && .dsvert_dp_is_string(value$peer_name)) {
      value$peer_name
    } else ""
  }, character(1L))
  if (any(!nzchar(receipt_peers)) || anyDuplicated(receipt_peers) ||
      !setequal(receipt_peers, context$servers)) {
    stop("The durable publication lacks exact-K compilation", call. = FALSE)
  }
  names(receipts) <- receipt_peers
  receipts <- receipts[context$servers]
  compilation <- list(
    version = .DSVERT_CLIENT_SYNOPSIS_COMPILE_SET_VERSION,
    artifact = publications[[1L]]$artifact,
    receipts = receipts,
    receipt_set_sha256 = publications[[1L]]$receipt_set_sha256)
  release_json <- stats::setNames(lapply(publications, function(value) {
    .dsvert_joint_dp_client_json(value$release_receipt)
  }), authorities)
  trusted <- .dsvert_dp_synopsis_client_bundle(
    bootstrap$manifest_bundle, bootstrap$status)
  compiled <- .dsvert_dp_synopsis_client_compile(
    compilation, trusted, bootstrap$manifest_bundle)
  execution <- .dsvert_dp_synopsis_client_execution(compiled)
  releases <- .dsvert_dp_synopsis_client_release_set(
    release_json, compiled, execution, trusted)
  if (!identical(as.numeric(publications[[1L]]$public_chunk_count),
                 as.numeric(execution$geometry$public_chunk_count))) {
    stop("The durable synopsis publication geometry is invalid",
         call. = FALSE)
  }
  list(
    publications = publications,
    json = stats::setNames(lapply(publications,
      .dsvert_joint_dp_client_json), authorities),
    compilation = compilation, compiled = compiled,
    release_json = releases$json,
    public_chunk_count = as.integer(publications[[1L]]$public_chunk_count))
}

.dsvert_dp_synopsis_finalize_ack_set_v1 <- function(
    raw, bootstrap, publication) {
  peers <- bootstrap$context$servers
  fields <- c(
    "version", "artifact_key", "manifest_sha256", "peer_name",
    "peer_identity_pk", "release_receipts_sha256",
    "compaction_authorization_sha256", "source_compaction_receipt_sha256",
    "source_intermediates_compacted", "durable_replay_retained",
    "idempotent", "session_required", "request_limit", "rate_limit",
    "catalog_limit", "signature")
  acks <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_dp_synopsis_client_json(
      raw[[peer]], paste("synopsis FinalizeAck from", peer),
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
    valid <- .dsvert_dp_has_exact_names(value, fields) &&
      identical(value$version,
                .DSVERT_CLIENT_SYNOPSIS_FINALIZE_ACK_VERSION) &&
      identical(value$artifact_key,
                publication$compiled$artifact$artifact_key) &&
      identical(value$manifest_sha256,
                bootstrap$manifest_bundle$manifest_sha256) &&
      identical(value$peer_name, peer) &&
      identical(value$peer_identity_pk,
                unname(bootstrap$context$pinset[[peer]])) &&
      .dsvert_vector_hex(value$release_receipts_sha256) &&
      .dsvert_vector_hex(value$compaction_authorization_sha256) &&
      .dsvert_vector_hex(value$source_compaction_receipt_sha256) &&
      identical(value$source_intermediates_compacted, TRUE) &&
      identical(value$durable_replay_retained, TRUE) &&
      identical(value$idempotent, TRUE) &&
      identical(value$session_required, FALSE) &&
      identical(value$request_limit, FALSE) &&
      identical(value$rate_limit, FALSE) &&
      identical(value$catalog_limit, FALSE)
    if (!isTRUE(valid)) {
      stop("A peer returned an invalid synopsis FinalizeAck", call. = FALSE)
    }
    .dsvert_dp_synopsis_outer_signature_v1(
      value, peer, bootstrap$context,
      .DSVERT_CLIENT_SYNOPSIS_FINALIZE_ACK_DOMAIN)
    value
  }), peers)
  if (length(unique(vapply(acks, `[[`, "", "release_receipts_sha256"))) !=
      1L) {
    stop("The K synopsis FinalizeAck receipts disagree", call. = FALSE)
  }
  invisible(acks)
}

.dsvert_dp_synopsis_cleanup_transient_v1 <- function(error) {
  inherits(error, c(
    "dsvert_resource_backpressure", "dsvert_retry_deadline_exceeded",
    "dsvert_transport_unavailable", "dsvert_phase_not_ready"))
}

.dsvert_dp_synopsis_publication_resume_v1 <- function(
    bootstrap, .aggregate = DSI::datashield.aggregate) {
  if (!inherits(bootstrap, "dsvert_synopsis_bootstrap_v1") ||
      !is.list(bootstrap$context) ||
      !identical(bootstrap$context$servers,
                 sort(names(bootstrap$context$all_conns), method = "radix"))) {
    stop("A validated synopsis bootstrap is required", call. = FALSE)
  }
  context <- bootstrap$context
  authorities <- context$designated
  publication_calls <- stats::setNames(lapply(authorities, function(peer) call(
    name = "dsvertDPSynopsisPublicationDS",
    manifest_sha256 = bootstrap$manifest_bundle$manifest_sha256)),
    authorities)
  available <- .dsvert_vector_try_phase(.dsvert_fanout_by_site(
    context$conns, publication_calls,
    operation = "durable synopsis publication lookup",
    .aggregate = .aggregate))
  if (!isTRUE(available$ok)) return(NULL)
  publication <- .dsvert_dp_synopsis_publication_set_v1(
    available$value, bootstrap)
  replay <- vector("list", publication$public_chunk_count)
  names(replay) <- as.character(seq.int(
    0L, publication$public_chunk_count - 1L))
  for (chunk_index in seq.int(0L, publication$public_chunk_count - 1L)) {
    calls <- stats::setNames(lapply(authorities, function(peer) call(
      name = "dsvertDPSynopsisPublishedReplayDS",
      artifact_key = publication$compiled$artifact$artifact_key,
      first_release_json = publication$json[[authorities[[1L]]]],
      second_release_json = publication$json[[authorities[[2L]]]],
      public_chunk_index = as.integer(chunk_index))), authorities)
    replay[[as.character(chunk_index)]] <-
      .dsvert_dp_synopsis_runner_json_set(.dsvert_fanout_by_site(
        context$conns, calls,
        operation = "durable synopsis published REPLAY",
        .aggregate = .aggregate), authorities, "published REPLAY",
        .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
  }
  release <- .dsvert_dp_synopsis_public_vector_v1(
    release_receipts = publication$release_json,
    replay_responses = replay,
    manifest_bundle = bootstrap$manifest_bundle,
    status = bootstrap$status,
    artifact = publication$compilation)
  finalize_calls <- stats::setNames(lapply(context$servers, function(peer) call(
    name = "dsvertDPSynopsisFinalizeAckDS",
    manifest_sha256 = bootstrap$manifest_bundle$manifest_sha256,
    first_release_json = publication$json[[authorities[[1L]]]],
    second_release_json = publication$json[[authorities[[2L]]]])),
    context$servers)
  finalize_raw <- tryCatch(
    .dsvert_fanout_by_site(
      context$all_conns, finalize_calls,
      operation = "K-wide synopsis FinalizeAck", .aggregate = .aggregate),
    error = function(error) {
      if (!.dsvert_dp_synopsis_cleanup_transient_v1(error)) stop(error)
      NULL
    })
  cleanup_pending <- is.null(finalize_raw)
  if (!cleanup_pending) {
    .dsvert_dp_synopsis_finalize_ack_set_v1(
      finalize_raw, bootstrap, publication)
  }
  list(
    release = release, layout = publication$compiled$layout,
    status = bootstrap$status,
    manifest_bundle = bootstrap$manifest_bundle,
    verification_compilation = publication$compilation,
    cleanup_pending = cleanup_pending)
}

.dsvert_dp_synopsis_trusted_bundle_v1 <- function(
    manifest_bundle, status) {
  context <- manifest_bundle$context
  if (!inherits(status, "ds.vertDPSynopsisStatus") ||
      !is.list(context) || !identical(context$servers, names(status)) ||
      !identical(unclass(context$status), unclass(status)) ||
      !identical(context$servers,
                 sort(names(context$all_conns), method = "radix")) ||
      !identical(context$pinset,
                 .dsvert_dp_synopsis_client_pinset_v1(
                   context$policy$ordered_peer_pinset, context$servers)) ||
      !identical(context$designated,
                 unlist(context$policy$designated_noise_peers,
                        use.names = FALSE))) {
    stop("The trusted synopsis status differs from its connection context",
         call. = FALSE)
  }
  policy <- .dsvert_dp_synopsis_client_policy_v1(
    context$policy, context$servers)
  bootstrap_fields <- c(
    "version", "phase", "peer_name", "peer_identity_pk", "policy",
    "draft", "data_access", "patient_derived_metadata", "request_limit",
    "rate_limit", "catalog_limit", "signature")
  for (peer in context$servers) {
    value <- status[[peer]]
    valid <- .dsvert_dp_has_exact_names(value, bootstrap_fields) &&
      identical(value$version, .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_VERSION) &&
      identical(value$phase, "custodian_policy_bootstrap") &&
      identical(value$peer_name, peer) &&
      identical(value$peer_identity_pk, unname(context$pinset[[peer]])) &&
      identical(.dsvert_joint_dp_client_json(value$policy),
                .dsvert_joint_dp_client_json(policy$value)) &&
      identical(value$data_access, FALSE) &&
      identical(value$patient_derived_metadata, FALSE) &&
      identical(value$request_limit, FALSE) &&
      identical(value$rate_limit, FALSE) &&
      identical(value$catalog_limit, FALSE)
    if (!isTRUE(valid)) {
      stop("The trusted synopsis bootstrap set is invalid", call. = FALSE)
    }
    .dsvert_dp_synopsis_outer_signature_v1(
      value, peer, context, .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_DOMAIN)
  }
  bootstrap_set_sha256 <- .dsvert_vector_hash(list(
    version = .DSVERT_CLIENT_SYNOPSIS_BOOTSTRAP_VERSION,
    bootstraps = unname(unclass(status))))
  if (!.dsvert_vector_hex(manifest_bundle$manifest_sha256) ||
      !identical(manifest_bundle$manifest_sha256, digest::digest(
        manifest_bundle$manifest_json, "sha256", serialize = FALSE)) ||
      !identical(manifest_bundle$deterministic_replay, TRUE) ||
      !identical(manifest_bundle$operation_or_request_limit, FALSE) ||
      !identical(manifest_bundle$history_can_deny_operation, FALSE)) {
    stop("The trusted no-lifetime synopsis bundle is invalid",
         call. = FALSE)
  }
  signed <- .dsvert_joint_dp_client_decode(
    manifest_bundle$schema_json, "signed no-lifetime synopsis schema",
    .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_SCHEMA_BYTES)
  signatures <- signed$signatures
  unsigned <- signed
  unsigned$signatures <- NULL
  if (!is.list(signatures) || is.null(names(signatures)) ||
      !setequal(names(signatures), context$servers) ||
      !identical(manifest_bundle$schema_sha256,
                 .dsvert_vector_hash(unsigned))) {
    stop("The trusted synopsis schema is invalid", call. = FALSE)
  }
  for (peer in context$servers) {
    .dsvert_dp_capsule_schema_verify(
      unsigned, signatures[[peer]], peer, context)
  }
  schema <- list(
    unsigned = unsigned, sha256 = manifest_bundle$schema_sha256,
    logical_snapshot = manifest_bundle$logical_snapshot,
    workload_contract = list(
      json = manifest_bundle$workload_contract_json,
      sha256 = manifest_bundle$workload_contract_sha256),
    signatures = signatures)
  receipts <- manifest_bundle$manifest_build_receipts
  if (!is.list(receipts) || is.null(names(receipts)) ||
      !setequal(names(receipts), context$servers) ||
      !is.list(manifest_bundle$manifest_signatures) ||
      !identical(names(manifest_bundle$manifest_signatures),
                 context$servers)) {
    stop("The synopsis manifest receipt coverage is invalid", call. = FALSE)
  }
  builds <- stats::setNames(lapply(context$servers, function(peer) {
    value <- receipts[[peer]]
    if (!identical(value$signature,
                   manifest_bundle$manifest_signatures[[peer]])) {
      stop("The synopsis manifest signature index is invalid",
           call. = FALSE)
    }
    value$manifest_json <- manifest_bundle$manifest_json
    value$artifact_commitments <- manifest_bundle$artifact_commitments
    .dsvert_dp_synopsis_bound_manifest_response_v1(
      .dsvert_joint_dp_client_json(value), peer, context, schema,
      bootstrap_set_sha256)
  }), context$servers)
  common <- setdiff(names(builds[[1L]]),
                    c("peer_name", "peer_identity_pk", "signature"))
  if (!all(vapply(builds, function(value) identical(
        .dsvert_joint_dp_client_json(value[common]),
        .dsvert_joint_dp_client_json(builds[[1L]][common])), logical(1L)))) {
    stop("The trusted synopsis manifest receipts disagree", call. = FALSE)
  }
  manifest <- .dsvert_dp_synopsis_manifest_value_v1(
    manifest_bundle$manifest_json, context, schema)
  index <- .dsvert_dp_capsule_artifact_commitment_index_client(
    manifest, .dsvert_dp_synopsis_policy_view_v1(context$policy),
    manifest_bundle$manifest_sha256)
  if (!identical(index$value, manifest_bundle$artifact_commitments) ||
      !identical(as.numeric(index$count),
                 as.numeric(manifest_bundle$artifact_commitment_count)) ||
      !identical(index$root, manifest_bundle$artifact_commitments_root)) {
    stop("The trusted synopsis artifact index is invalid", call. = FALSE)
  }
  list(
    context = context, status = status, manifest = manifest,
    validated_manifest = list(value = manifest,
                              capsule_id = manifest_bundle$capsule_id))
}

.dsvert_dp_synopsis_source_manifest_v1 <- function(manifest, context) {
  workload <- if (is.list(manifest)) manifest$workload else NULL
  mechanism <- if (is.list(workload)) workload$capsule_mechanism else NULL
  count <- if (is.list(workload)) workload$coordinate_count else NULL
  if (!is.list(context) || !inherits(context$status,
                                     "ds.vertDPSynopsisStatus") ||
      !is.list(manifest$capsule_identity) ||
      !.dsvert_vector_hex(manifest$capsule_identity$capsule_id) ||
      !.dsvert_dp_is_integer(count, 1L, .DSVERT_DP_MAX_COORDINATES) ||
      !is.list(mechanism) ||
      !.dsvert_vector_hex(mechanism$source_context_hash)) {
    stop("The trusted synopsis source manifest is invalid", call. = FALSE)
  }
  list(
    value = manifest,
    capsule_id = manifest$capsule_identity$capsule_id,
    coordinate_count = as.numeric(count),
    release_coordinate_count = as.numeric(count),
    release_coordinate_order_sha256 = NULL,
    private_layout_sha256 = NULL,
    cross_enabled = FALSE,
    purpose = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_PURPOSE,
    source_context_hash = mechanism$source_context_hash)
}
