# Client orchestration for the purpose-bound one-draw vector exact-GC route.
# Only the four generic exact-GC transport endpoints are used.  All values
# visible here are public contracts, signed initialization state, or opaque
# authenticated envelopes.

.DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_SELECTION_VERSION <-
  "dsvert-joint-dp-vector-backend-selection-v2"
.DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_COST_POLICY_VERSION <-
  "dsvert-joint-dp-vector-exact-gc-cost-policy-v1"
.DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_MAX_PROMOTED_COORDINATES <- 1L
.DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BINDING_VERSION <-
  "dsvert-joint-dp-vector-exact-gc-binding-v1"
.DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BACKEND <-
  "exact_gc_one_joint_discrete_laplace_draw_ring128_v3"
.DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_OPERATION <-
  "joint-dp-vector-laplace-v3"
.DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_OUTPUT_KIND <-
  "joint-dp-vector-ring128-share-v1"

.dsvert_joint_dp_vector_exact_gc_client_hash <- function(value) {
  digest::digest(
    .dsvert_joint_dp_client_json(value),
    algo = "sha256", serialize = FALSE)
}

.dsvert_joint_dp_vector_exact_gc_client_hex <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[0-9a-f]{64}$", value)) {
    stop("Invalid ", what, ".", call. = FALSE)
  }
  value
}

.dsvert_joint_dp_vector_exact_gc_client_assessment <- function(
    assessment, manifest_sha256) {
  manifest_sha256 <- .dsvert_joint_dp_vector_exact_gc_client_hex(
    manifest_sha256, "biomedical manifest hash")
  base <- c(
    "version", "manifest_sha256", "representable",
    "exact_gc_capability_id", "plan_sha256",
    "maximum_chunk_coordinates", "cost_policy_version",
    "total_coordinate_count", "maximum_promoted_coordinates",
    "promoted", "selection_reason", "private_material_accessed",
    "runtime_failure_consulted",
    "assessment_sha256")
  required <- base
  maximum <- if (is.list(assessment)) suppressWarnings(as.numeric(
    assessment$maximum_chunk_coordinates)) else NA_real_
  total <- if (is.list(assessment)) suppressWarnings(as.numeric(
    assessment$total_coordinate_count)) else NA_real_
  promoted <- length(total) == 1L && !is.na(total) && is.finite(total) &&
    total == floor(total) && total >= 1L && total <= 1000000L &&
    total <= .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_MAX_PROMOTED_COORDINATES
  reason <- if (isTRUE(promoted)) {
    "within_public_exact_gc_cost_ceiling"
  } else {
    "above_public_exact_gc_cost_ceiling"
  }
  coherent <- identical(assessment$representable, TRUE) &&
    length(maximum) == 1L && !is.na(maximum) && is.finite(maximum) &&
    maximum == floor(maximum) && maximum >= 1L &&
    maximum <= .DSVERT_CLIENT_JOINT_DP_VECTOR_MAX_CHUNK &&
    identical(assessment$cost_policy_version,
              .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_COST_POLICY_VERSION) &&
    identical(as.numeric(assessment$maximum_promoted_coordinates), 1) &&
    identical(assessment$promoted, promoted) &&
    identical(assessment$selection_reason, reason)
  valid <- is.list(assessment) && setequal(names(assessment), required) &&
    identical(assessment$version,
              "dsvert-joint-dp-vector-exact-gc-assessment-v2") &&
    identical(assessment$manifest_sha256, manifest_sha256) &&
    identical(assessment$exact_gc_capability_id,
              "joint_dp_biomedical_vector_exact_gc_v1") &&
    identical(assessment$private_material_accessed, FALSE) &&
    identical(assessment$runtime_failure_consulted, FALSE) && coherent &&
    identical(
      .dsvert_joint_dp_vector_exact_gc_client_hash(
        assessment[setdiff(names(assessment), "assessment_sha256")]),
      assessment$assessment_sha256)
  if (!isTRUE(valid)) {
    stop("The signed exact-GC vector assessment is invalid.",
         call. = FALSE)
  }
  invisible(assessment)
}

.dsvert_joint_dp_vector_exact_gc_client_selection <- function(
    selection, manifest_sha256, require_exact = NULL) {
  manifest_sha256 <- .dsvert_joint_dp_vector_exact_gc_client_hex(
    manifest_sha256, "biomedical manifest hash")
  required <- c(
    "version", "manifest_sha256", "backend", "one_draw",
    "cost_policy_version", "total_coordinate_count",
    "maximum_promoted_coordinates", "selection_reason",
    "assessment_sha256",
    "exact_gc_plan_sha256", "exact_gc_maximum_chunk_coordinates",
    "selected_before_private_material", "retry_may_change_backend",
    "selection_sha256")
  maximum <- if (is.list(selection)) suppressWarnings(as.numeric(
    selection$exact_gc_maximum_chunk_coordinates)) else NA_real_
  exact <- is.list(selection) && identical(
    selection$backend,
    .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BACKEND)
  total <- if (is.list(selection)) suppressWarnings(as.numeric(
    selection$total_coordinate_count)) else NA_real_
  promoted <- length(total) == 1L && !is.na(total) && is.finite(total) &&
    total == floor(total) && total >= 1L && total <= 1000000L &&
    total <= .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_MAX_PROMOTED_COORDINATES
  expected_backend <- if (isTRUE(promoted)) {
    .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BACKEND
  } else {
    .DSVERT_CLIENT_VECTOR_BACKEND
  }
  expected_reason <- if (isTRUE(promoted)) {
    "within_public_exact_gc_cost_ceiling"
  } else {
    "above_public_exact_gc_cost_ceiling"
  }
  coherent <- identical(selection$backend, expected_backend) &&
    identical(selection$one_draw, promoted) &&
    identical(selection$cost_policy_version,
              .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_COST_POLICY_VERSION) &&
    identical(as.numeric(selection$maximum_promoted_coordinates), 1) &&
    identical(selection$selection_reason, expected_reason) &&
    length(maximum) == 1L && !is.na(maximum) && is.finite(maximum) &&
    maximum == floor(maximum) && maximum >= 1L &&
    maximum <= .DSVERT_CLIENT_JOINT_DP_VECTOR_MAX_CHUNK
  valid <- is.list(selection) && setequal(names(selection), required) &&
    identical(selection$version,
      .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_SELECTION_VERSION) &&
    identical(selection$manifest_sha256, manifest_sha256) &&
    identical(selection$selected_before_private_material, TRUE) &&
    identical(selection$retry_may_change_backend, FALSE) &&
    .dsvert_joint_dp_vector_exact_gc_client_hex(
      selection$assessment_sha256, "exact-GC assessment hash") ==
      selection$assessment_sha256 &&
    .dsvert_joint_dp_vector_exact_gc_client_hex(
      selection$exact_gc_plan_sha256, "exact-GC plan hash") ==
      selection$exact_gc_plan_sha256 &&
    coherent &&
    (is.null(require_exact) || identical(exact, isTRUE(require_exact))) &&
    identical(
      .dsvert_joint_dp_vector_exact_gc_client_hash(
        selection[setdiff(names(selection), "selection_sha256")]),
      selection$selection_sha256)
  if (!isTRUE(valid)) {
    stop("The signed one-draw exact-GC vector selection is invalid.",
         call. = FALSE)
  }
  invisible(selection)
}

.dsvert_joint_dp_vector_exact_gc_client_binding <- function(
    binding, manifest_sha256, release_contract_hash,
    selection_sha256, transcript_hash, chunk_index, coordinate_count) {
  manifest_sha256 <- .dsvert_joint_dp_vector_exact_gc_client_hex(
    manifest_sha256, "biomedical manifest hash")
  release_contract_hash <- .dsvert_joint_dp_vector_exact_gc_client_hex(
    release_contract_hash, "vector release contract hash")
  selection_sha256 <- .dsvert_joint_dp_vector_exact_gc_client_hex(
    selection_sha256, "vector backend selection hash")
  transcript_hash <- .dsvert_joint_dp_vector_exact_gc_client_hex(
    transcript_hash, "vector transcript hash")
  chunk_index <- suppressWarnings(as.numeric(chunk_index))
  coordinate_count <- suppressWarnings(as.numeric(coordinate_count))
  required <- c(
    "version", "domain", "manifest_sha256", "release_contract_hash",
    "selection_sha256", "transcript_hash", "chunk_index",
    "coordinate_count", "circuit_digest", "purpose", "operation_id",
    "source_key", "output_key", "operation", "output_kind",
    "source_producer", "binding_sha256")
  identity_names <- c(
    "domain", "manifest_sha256", "release_contract_hash",
    "selection_sha256", "transcript_hash", "chunk_index",
    "coordinate_count", "circuit_digest", "purpose")
  expected_suffix <- if (is.list(binding) &&
                         all(identity_names %in% names(binding))) {
    substr(.dsvert_joint_dp_vector_exact_gc_client_hash(
      binding[identity_names]), 1L, 32L)
  } else {
    ""
  }
  if (!is.list(binding) || !setequal(names(binding), required) ||
      !identical(binding$version,
        .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BINDING_VERSION) ||
      !identical(binding$domain,
        "dsVert/joint-dp/vector/exact-gc-operation/v1") ||
      !identical(binding$manifest_sha256, manifest_sha256) ||
      !identical(binding$release_contract_hash, release_contract_hash) ||
      !identical(binding$selection_sha256, selection_sha256) ||
      !identical(binding$transcript_hash, transcript_hash) ||
      length(chunk_index) != 1L || is.na(chunk_index) ||
      !is.finite(chunk_index) || chunk_index != floor(chunk_index) ||
      chunk_index < 0 || chunk_index > 1000000 ||
      !identical(as.numeric(binding$chunk_index), chunk_index) ||
      length(coordinate_count) != 1L || is.na(coordinate_count) ||
      !is.finite(coordinate_count) ||
      coordinate_count != floor(coordinate_count) ||
      coordinate_count < 1 ||
      coordinate_count > .DSVERT_CLIENT_JOINT_DP_VECTOR_MAX_CHUNK ||
      !identical(as.numeric(binding$coordinate_count), coordinate_count) ||
      !identical(binding$operation,
        .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_OPERATION) ||
      !identical(binding$output_kind,
        .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_OUTPUT_KIND) ||
      !identical(binding$purpose, paste0(
        .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_OPERATION, "/",
        binding$circuit_digest)) ||
      !grepl("^[0-9a-f]{64}$", binding$circuit_digest) ||
      !grepl("^op_[0-9a-f]{32}$", binding$operation_id) ||
      !grepl("^exact_gc_in_[0-9a-f]{32}$", binding$source_key) ||
      !grepl("^exact_gc_out_[0-9a-f]{32}$", binding$output_key) ||
      !identical(sub("^op_", "", binding$operation_id),
                 expected_suffix) ||
      !identical(sub("^op_", "", binding$operation_id),
                 sub("^exact_gc_in_", "", binding$source_key)) ||
      !identical(sub("^op_", "", binding$operation_id),
                 sub("^exact_gc_out_", "", binding$output_key)) ||
      !identical(
        .dsvert_joint_dp_vector_exact_gc_client_hash(
          binding[setdiff(names(binding), "binding_sha256")]),
        binding$binding_sha256)) {
    stop("The signed one-draw exact-GC vector binding is invalid.",
         call. = FALSE)
  }
  invisible(binding)
}

.dsvert_joint_dp_vector_exact_gc_initializations <- function(
    started, server_names, binding) {
  if (!is.list(started) || !is.character(server_names) ||
      length(server_names) != 2L || anyNA(server_names) ||
      anyDuplicated(server_names) || !all(server_names %in% names(started)) ||
      !is.list(binding) ||
      !grepl("^[0-9a-f]{64}$", binding$binding_sha256) ||
      !grepl("^op_[0-9a-f]{32}$", binding$operation_id)) {
    stop("The one-draw exact-GC start set is incomplete.", call. = FALSE)
  }
  result <- stats::setNames(lapply(server_names, function(server) {
    value <- started[[server]]
    if (!is.list(value) ||
        !identical(value$backend,
          .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BACKEND) ||
        !identical(value$binding_sha256, binding$binding_sha256) ||
        !identical(value$operation_id, binding$operation_id) ||
        !identical(value$purpose, binding$purpose) ||
        !identical(value$intermediate_payload_exposed, FALSE) ||
        !identical(value$source_share_exposed, FALSE) ||
        !identical(value$private_seed_exposed, FALSE) ||
        !identical(value$preclamp_values_exposed, FALSE) ||
        !is.list(value$initialization)) {
      stop("A vector peer returned an unsafe exact-GC start receipt.",
           call. = FALSE)
    }
    value$initialization
  }), server_names)
  result
}

.dsvert_joint_dp_vector_exact_gc_run <- function(
    datasources, server_names = names(datasources), servers,
    session_id, binding, manifest_sha256, release_contract_hash,
    selection_sha256, transcript_hash, chunk_index, coordinate_count,
    initialized, transport_ready = TRUE,
    timeout_seconds = getOption("dsvert.exact_gc.timeout_seconds", 900),
    .run = .dsvert_exact_gc_run,
    .aggregate = DSI::datashield.aggregate) {
  if (!is.function(.run) || length(servers) != 2L ||
      anyNA(servers) || anyDuplicated(servers)) {
    stop("One-draw exact-GC requires exactly two compute peers.",
         call. = FALSE)
  }
  .dsvert_joint_dp_vector_exact_gc_client_binding(
    binding, manifest_sha256, release_contract_hash, selection_sha256,
    transcript_hash, chunk_index, coordinate_count)
  selected_names <- server_names[servers]
  if (length(selected_names) != 2L || anyNA(selected_names) ||
      any(!nzchar(selected_names)) || anyDuplicated(selected_names)) {
    stop("One-draw exact-GC requires two unique pinned peer names.",
         call. = FALSE)
  }
  if (is.list(initialized) && all(selected_names %in% names(initialized)) &&
      all(vapply(initialized[selected_names], function(value) {
        is.list(value) && !is.null(value$initialization)
      }, logical(1L)))) {
    initialized <- .dsvert_joint_dp_vector_exact_gc_initializations(
      initialized, selected_names, binding)
  }
  result <- .run(
    datasources, server_names = server_names, servers = servers,
    session_id = session_id, operation_id = binding$operation_id,
    source_key = binding$source_key, output_key = binding$output_key,
    operation = .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_OPERATION,
    ring = 128L, frac_bits = 0L,
    vector_len = as.integer(coordinate_count), purpose = binding$purpose,
    transport_ready = transport_ready,
    timeout_seconds = timeout_seconds, initialized = initialized,
    .aggregate = .aggregate)
  # The generic pump result is public liveness only.  Restate a deliberately
  # narrow contract so a future transport field cannot accidentally surface.
  invisible(list(
    backend = .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_BACKEND,
    binding_sha256 = binding$binding_sha256,
    operation_id = binding$operation_id,
    purpose = binding$purpose,
    complete = TRUE,
    intermediate_payload_exposed = FALSE,
    source_share_exposed = FALSE,
    private_seed_exposed = FALSE,
    preclamp_values_exposed = FALSE))
}
