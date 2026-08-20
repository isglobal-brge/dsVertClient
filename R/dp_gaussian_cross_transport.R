# Client orchestration for the fixed cross-owner Gaussian exact transcript.
# The relay sees only signed public receipts and server-minted manifest handles;
# values, validity masks, products, sums, and additive shares stay at the two
# pinned computation peers.

.DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION <-
  "bounded-normalized-gaussian-cross-sufficient-statistics-v1"
.DSVERT_CLIENT_DP_GAUSSIAN_CROSS_LAYOUT_VERSION <-
  "dsvert-cross-gaussian-private-source-layout-v1"
.DSVERT_CLIENT_DP_GAUSSIAN_CROSS_BIND_VERSION <-
  "dsvert-cross-gaussian-exact-binding-v1"
.DSVERT_CLIENT_DP_GAUSSIAN_CROSS_STAGE_VERSION <-
  "dsvert-cross-gaussian-exact-stage-v1"
.DSVERT_CLIENT_DP_GAUSSIAN_CROSS_RECEIPT_VERSION <-
  "dsvert-cross-gaussian-result-receipt-v1"
.DSVERT_CLIENT_DP_GAUSSIAN_CROSS_PRODUCER <- "dp.gaussian-cross.v1"
.DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_TRANSPORT_COORDINATES <-
  64L * 1024L^2
.DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_RECEIPT_BYTES <- 128L * 1024L

.dsvert_dp_gaussian_cross_artifacts_client <- function(manifest) {
  artifacts <- tryCatch(
    manifest$workload$families$gaussian_models$artifacts,
    error = function(error) NULL)
  if (!is.list(artifacts)) return(list())
  result <- artifacts[vapply(artifacts, function(artifact) {
    is.list(artifact) && identical(
      artifact$version,
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION)
  }, logical(1L))]
  if (!length(result)) return(list())
  if (is.null(names(result)) || anyNA(names(result)) ||
      any(!nzchar(names(result))) || anyDuplicated(names(result))) {
    stop("The signed cross-owner Gaussian artifact set is invalid",
         call. = FALSE)
  }
  result[order(names(result), method = "radix")]
}

.dsvert_dp_gaussian_cross_names_client <- function(value, what,
                                                   qualified = FALSE) {
  if (is.list(value) && is.null(names(value))) {
    valid <- all(vapply(value, function(item) {
      is.character(item) && length(item) == 1L && !is.na(item)
    }, logical(1L)))
    if (isTRUE(valid)) value <- unname(unlist(value, use.names = FALSE))
  }
  valid_item <- if (isTRUE(qualified)) {
    function(item) !is.null(.dsvert_dp_gaussian_reference(item))
  } else {
    function(item) grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", item)
  }
  if (!is.character(value) || !length(value) || !is.null(names(value)) ||
      anyNA(value) || any(!vapply(value, valid_item, logical(1L)))) {
    stop("Invalid signed cross-owner Gaussian ", what, call. = FALSE)
  }
  unname(enc2utf8(value))
}

.dsvert_dp_gaussian_cross_variable_client <- function(artifact, variable) {
  descriptor <- artifact$predictors[[variable]]
  if (is.list(descriptor)) return(descriptor)
  if (.dsvert_dp_gaussian_reference_matches(variable, artifact$outcome)) {
    return(artifact$outcome)
  }
  stop("The signed cross-owner Gaussian variable order is invalid",
       call. = FALSE)
}

.dsvert_dp_gaussian_cross_layout_client <- function(
    manifest, release_layout = NULL) {
  if (is.null(release_layout)) {
    release_layout <- .dsvert_dp_capsule_vector_layout(manifest)
  }
  artifacts <- .dsvert_dp_gaussian_cross_artifacts_client(manifest)
  categorical_artifacts <-
    .dsvert_dp_categorical_cross_artifacts_client(manifest)
  if (!length(artifacts) && !length(categorical_artifacts)) {
    return(list(
      version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_LAYOUT_VERSION,
      enabled = FALSE,
      release_coordinate_count = release_layout$coordinate_count,
      transport_coordinate_count = release_layout$coordinate_count,
      release_coordinate_order_sha256 = release_layout$sha256,
      transport_coordinate_order_sha256 = release_layout$sha256,
      padding_coordinates = 0L, blocks = list(),
      source_peers = character(), computation_peers = character()))
  }
  capacity <- suppressWarnings(as.numeric(manifest$admission$unit_capacity))
  if (length(capacity) != 1L || is.na(capacity) || !is.finite(capacity) ||
      capacity != floor(capacity) || capacity < 1 ||
      capacity > .DSVERT_DP_MAX_COORDINATES) {
    stop("The signed cross-owner Gaussian padded capacity is invalid",
         call. = FALSE)
  }
  capacity <- as.integer(capacity)
  release_count <- as.numeric(release_layout$coordinate_count)
  chunk <- as.numeric(.DSVERT_CLIENT_DP_CAPSULE_SOURCE_CHUNK_COORDINATES)
  private_start <- ceiling(release_count / chunk) * chunk + 1
  padding <- private_start - release_count - 1
  cursor <- private_start
  blocks <- list()
  source_peers <- computation_peers <- character()
  for (analysis_id in names(artifacts)) {
    artifact <- artifacts[[analysis_id]]
    variables <- .dsvert_dp_gaussian_cross_names_client(
      artifact$input_variable_order, "input-variable order",
      qualified = TRUE)
    participants <- .dsvert_dp_gaussian_cross_names_client(
      artifact$participating_peers, "participant list")
    computation <- .dsvert_dp_gaussian_cross_names_client(
      artifact$computation_peers, "computation-peer list")
    if (!identical(artifact$analysis_id, analysis_id) ||
        anyDuplicated(variables) || length(variables) < 2L ||
        length(computation) != 2L || anyDuplicated(computation) ||
        !identical(computation, sort(computation, method = "radix")) ||
        anyDuplicated(participants) ||
        !identical(participants, sort(participants, method = "radix"))) {
      stop("The signed cross-owner Gaussian private layout is invalid",
           call. = FALSE)
    }
    source_peers <- c(source_peers, participants)
    computation_peers <- c(computation_peers, computation)
    for (variable in variables) {
      descriptor <- .dsvert_dp_gaussian_cross_variable_client(
        artifact, variable)
      required <- c("column", "dataset", "owner_peer", "lower", "upper")
      if (!.dsvert_dp_has_exact_names(descriptor, required) ||
          !.dsvert_dp_gaussian_reference_matches(variable, descriptor) ||
          !descriptor$owner_peer %in% participants ||
          !.dsvert_dp_is_number(descriptor$lower) ||
          !.dsvert_dp_is_number(descriptor$upper) ||
          descriptor$lower >= descriptor$upper) {
        stop("The signed cross-owner Gaussian input descriptor is invalid",
             call. = FALSE)
      }
      for (kind in c("value", "validity")) {
        end <- cursor + capacity - 1
        if (!is.finite(end) ||
            end > .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_TRANSPORT_COORDINATES) {
          stop("The cross-owner Gaussian private transport shape is too large",
               call. = FALSE)
        }
        key <- paste(analysis_id, variable, kind, sep = "::")
        blocks[[key]] <- list(
          analysis_id = analysis_id, variable = descriptor$column,
          kind = kind,
          dataset = descriptor$dataset, owner_peer = descriptor$owner_peer,
          lower = descriptor$lower, upper = descriptor$upper,
          start = as.integer(cursor), end = as.integer(end),
          length = capacity)
        cursor <- end + 1
      }
    }
  }
  for (analysis_id in names(categorical_artifacts)) {
    artifact <- categorical_artifacts[[analysis_id]]
    participants <- .dsvert_dp_gaussian_cross_names_client(
      artifact$participating_peers, "categorical participant list")
    computation <- .dsvert_dp_gaussian_cross_names_client(
      artifact$computation_peers, "categorical computation-peer list")
    if (!identical(artifact$analysis_id, analysis_id) ||
        length(participants) != 2L || anyDuplicated(participants) ||
        length(computation) != 2L || anyDuplicated(computation) ||
        !identical(participants, sort(participants, method = "radix")) ||
        !identical(computation, sort(computation, method = "radix"))) {
      stop("The signed cross-owner categorical private layout is invalid",
           call. = FALSE)
    }
    source_peers <- c(source_peers, participants)
    computation_peers <- c(computation_peers, computation)
    for (side in c("left", "right")) {
      descriptor <- artifact[[side]]
      required <- c("dataset", "column", "owner_peer", "levels")
      levels <- tryCatch(.dsvert_dp_capsule_manifest_strings(
        descriptor$levels, paste(side, "categorical levels"),
        sorted = TRUE), error = function(error) NULL)
      if (!.dsvert_dp_has_exact_names(descriptor, required) ||
          !descriptor$owner_peer %in% participants ||
          is.null(levels)) {
        stop("The signed cross-owner categorical input descriptor is invalid",
             call. = FALSE)
      }
      shapes <- c(
        one_hot = as.double(capacity) * length(levels),
        validity = capacity)
      for (kind in names(shapes)) {
        block_length <- shapes[[kind]]
        end <- cursor + block_length - 1
        if (!is.finite(end) || block_length < 1L ||
            block_length > 2^31 - 1 ||
            end > .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_TRANSPORT_COORDINATES) {
          stop("The cross-owner categorical private transport shape is too large",
               call. = FALSE)
        }
        key <- paste("categorical", analysis_id, side, kind, sep = "::")
        blocks[[key]] <- list(
          input_family = "categorical", analysis_id = analysis_id,
          side = side, variable = descriptor$column, kind = kind,
          dataset = descriptor$dataset, owner_peer = descriptor$owner_peer,
          levels = descriptor$levels,
          start = as.integer(cursor), end = as.integer(end),
          length = as.integer(block_length))
        cursor <- end + 1
      }
    }
  }
  source_peers <- sort(unique(source_peers), method = "radix")
  computation_peers <- sort(unique(computation_peers), method = "radix")
  if (length(computation_peers) != 2L) {
    stop("Cross-owner artifacts disagree on their two computation peers",
         call. = FALSE)
  }
  shape <- list(
    version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_LAYOUT_VERSION,
    capsule_id = manifest$capsule_identity$capsule_id,
    release_coordinate_count = as.integer(release_count),
    release_coordinate_order_sha256 = release_layout$sha256,
    private_start = as.integer(private_start),
    padding_coordinates = as.integer(padding),
    transport_coordinate_count = as.integer(cursor - 1),
    blocks = blocks, source_peers = as.list(source_peers),
    computation_peers = as.list(computation_peers),
    padding_rule = "zero_to_next_source_chunk_boundary_v1",
    payload_rule = if (length(categorical_artifacts)) paste0(
      "manifest_order_capacity_padded_ring128_value_one_hot_and_validity_",
      "no_exact_release_v1") else
      "manifest_order_capacity_padded_ring128_value_then_validity_no_exact_release_v1")
  shape$transport_coordinate_order_sha256 <-
    .dsvert_dp_capsule_source_hash(shape)
  c(shape, list(enabled = TRUE))
}

.dsvert_dp_gaussian_cross_tag_client <- function(capsule_id, analysis_id) {
  substr(.dsvert_dp_capsule_source_hash(list(
    protocol = "dsvert-cross-gaussian-slot-namespace-v1",
    capsule_id = capsule_id, analysis_id = analysis_id)), 1L, 20L)
}

.dsvert_dp_gaussian_cross_verify_signed <- function(
    value, domain, peer, context, expected_fields) {
  if (!.dsvert_dp_has_exact_names(value, expected_fields)) {
    stop("A pinned peer returned a malformed cross-owner Gaussian receipt",
         call. = FALSE)
  }
  .dsvert_dp_capsule_source_verify(value, domain, peer, context)
  value
}

.dsvert_dp_gaussian_cross_bind_set <- function(
    responses, context, manifest, layout, analysis_id) {
  peers <- context$designated
  if (!is.list(responses) || !setequal(names(responses), peers)) {
    stop("Cross-owner Gaussian binding did not cover both computation peers",
         call. = FALSE)
  }
  artifact <- .dsvert_dp_gaussian_cross_artifacts_client(manifest)[[analysis_id]]
  fields <- c(
    "version", "phase", "capsule_id", "analysis_id", "artifact_sha256",
    "source_contract_hash", "private_layout_sha256", "transcript_sha256",
    "numeric_certificate_sha256", "peer_name", "peer_identity_pk",
    "padded_units", "variable_count", "ring_bits", "frac_bits", "state",
    "source_values_exposed", "alignment_hash_exposed",
    "alignment_hash_exposed_to_relay",
    "alignment_hash_exposed_to_computation_peers",
    "exact_intermediates_exposed", "fixed_transcript", "signature")
  receipts <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_joint_dp_client_decode(
      responses[[peer]], "cross-owner Gaussian binding receipt",
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_RECEIPT_BYTES)
    value <- .dsvert_dp_gaussian_cross_verify_signed(
      value, "cross-gaussian-bind", peer, context, fields)
    valid <- identical(value$version,
                       .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_BIND_VERSION) &&
      identical(value$phase, "cross_gaussian_private_inputs_bound") &&
      identical(value$capsule_id, manifest$capsule_identity$capsule_id) &&
      identical(value$analysis_id, analysis_id) &&
      identical(value$artifact_sha256,
                .dsvert_dp_capsule_source_hash(artifact)) &&
      .dsvert_dp_capsule_source_hex(value$source_contract_hash) &&
      identical(value$private_layout_sha256,
                layout$transport_coordinate_order_sha256) &&
      identical(value$transcript_sha256,
                .dsvert_dp_capsule_source_hash(artifact$transcript)) &&
      identical(value$numeric_certificate_sha256,
                .dsvert_dp_capsule_source_hash(artifact$numeric_certificate)) &&
      identical(value$peer_name, peer) &&
      identical(value$peer_identity_pk, unname(context$pinset[[peer]])) &&
      identical(as.numeric(value$padded_units),
                as.numeric(artifact$transcript$padded_units)) &&
      identical(as.numeric(value$variable_count),
                as.numeric(length(artifact$input_variable_order))) &&
      identical(as.numeric(value$ring_bits), 128) &&
      identical(as.numeric(value$frac_bits),
                as.numeric(artifact$numeric_grid_bits)) &&
      value$state %in% c("bound", "complete") &&
      identical(value$source_values_exposed, FALSE) &&
      identical(value$alignment_hash_exposed, FALSE) &&
      identical(value$alignment_hash_exposed_to_relay, FALSE) &&
      identical(value$alignment_hash_exposed_to_computation_peers, FALSE) &&
      identical(value$exact_intermediates_exposed, FALSE) &&
      identical(value$fixed_transcript, TRUE)
    if (!isTRUE(valid)) {
      stop("A pinned peer returned a misbound cross-owner Gaussian receipt",
           call. = FALSE)
    }
    value
  }), peers)
  stable <- setdiff(fields, c(
    "peer_name", "peer_identity_pk", "state", "signature"))
  if (length(unique(vapply(receipts, function(value) {
        .dsvert_joint_dp_client_json(value[stable])
      }, character(1L)))) != 1L ||
      length(unique(vapply(receipts, `[[`, character(1L), "state"))) != 1L) {
    stop("The computation peers disagree on the cross-owner Gaussian binding",
         call. = FALSE)
  }
  list(
    artifact = artifact, receipts = receipts,
    state = receipts[[1L]]$state,
    source_contract_hash = receipts[[1L]]$source_contract_hash)
}

.dsvert_dp_gaussian_cross_stage_contract <- function(
    artifact, capsule_id, analysis_id, stage, stage_index) {
  m <- length(artifact$input_variable_order)
  capacity <- as.numeric(artifact$transcript$padded_units)
  tag <- .dsvert_dp_gaussian_cross_tag_client(capsule_id, analysis_id)
  if (identical(stage, "validity")) {
    if (stage_index < 1L || stage_index > m - 1L) {
      stop("Invalid cross-owner Gaussian validity stage", call. = FALSE)
    }
    total_n <- capacity
    purpose <- paste0(
      "dp.gaussian-cross.", tag, ".validity-",
      sprintf("%04d", as.integer(stage_index)))
  } else if (identical(stage, "masked-values") && stage_index == 1L) {
    total_n <- capacity * m
    purpose <- paste0("dp.gaussian-cross.", tag, ".masked-values")
  } else if (identical(stage, "moments") && stage_index == 1L) {
    total_n <- capacity * m * (m + 1L) / 2L
    purpose <- paste0("dp.gaussian-cross.", tag, ".moments")
  } else {
    stop("Invalid cross-owner Gaussian exact stage", call. = FALSE)
  }
  if (!is.finite(total_n) || total_n != floor(total_n) || total_n < 1 ||
      total_n > 2^31 - 1) {
    stop("The cross-owner Gaussian exact stage is not representable",
         call. = FALSE)
  }
  list(
    stage = stage, stage_index = as.integer(stage_index),
    total_n = as.integer(total_n), purpose = purpose)
}

.dsvert_dp_gaussian_cross_stage_set <- function(
    responses, context, manifest, binding, analysis_id, stage, stage_index,
    session_id) {
  peers <- context$designated
  artifact <- binding$artifact
  expected <- .dsvert_dp_gaussian_cross_stage_contract(
    artifact, manifest$capsule_identity$capsule_id, analysis_id,
    stage, stage_index)
  fields <- c(
    "capability_id", "manifest_handle", "context_hash", "plan_id",
    "ring_bits", "frac_bits", "backend", "bound_x", "bound_y",
    "max_chunk", "total_n", "numeric_policy_id", "version", "state",
    "producer", "purpose", "capsule_id", "analysis_id", "stage",
    "stage_index", "artifact_sha256", "source_contract_hash",
    "transcript_sha256", "numeric_certificate_sha256",
    "exact_intermediates_exposed", "source_values_exposed")
  if (!is.list(responses) || !setequal(names(responses), peers)) {
    stop("Cross-owner Gaussian preparation did not cover both peers",
         call. = FALSE)
  }
  manifests <- stats::setNames(lapply(peers, function(peer) {
    value <- responses[[peer]]
    valid <- .dsvert_dp_has_exact_names(value, fields) &&
      identical(value$capability_id, .DSVERT_CLIENT_EXACT_GC_CAPABILITY) &&
      identical(value$version,
                .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_STAGE_VERSION) &&
      value$state %in% c("prepared", "complete") &&
      identical(value$producer,
                .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_PRODUCER) &&
      identical(value$purpose, expected$purpose) &&
      identical(value$capsule_id, manifest$capsule_identity$capsule_id) &&
      identical(value$analysis_id, analysis_id) &&
      identical(value$stage, expected$stage) &&
      identical(as.numeric(value$stage_index),
                as.numeric(expected$stage_index)) &&
      identical(value$artifact_sha256,
                .dsvert_dp_capsule_source_hash(artifact)) &&
      identical(value$source_contract_hash,
                binding$source_contract_hash) &&
      identical(value$transcript_sha256,
                .dsvert_dp_capsule_source_hash(artifact$transcript)) &&
      identical(value$numeric_certificate_sha256,
                .dsvert_dp_capsule_source_hash(artifact$numeric_certificate)) &&
      is.character(value$manifest_handle) &&
      length(value$manifest_handle) == 1L &&
      grepl("^[A-Za-z0-9_-]{43}$", value$manifest_handle) &&
      all(vapply(value[c(
        "context_hash", "plan_id", "numeric_policy_id")],
        .dsvert_dp_capsule_source_hex, logical(1L))) &&
      identical(as.numeric(value$ring_bits), 128) &&
      identical(as.numeric(value$frac_bits),
                as.numeric(artifact$numeric_grid_bits)) &&
      identical(value$backend, "direct-wide") &&
      identical(value$bound_x,
                format(2^artifact$numeric_grid_bits,
                       scientific = FALSE, trim = TRUE)) &&
      identical(value$bound_y, value$bound_x) &&
      identical(as.numeric(value$total_n), as.numeric(expected$total_n)) &&
      .dsvert_dp_is_integer(value$max_chunk, 1, 256) &&
      identical(value$exact_intermediates_exposed, FALSE) &&
      identical(value$source_values_exposed, FALSE)
    if (!isTRUE(valid)) {
      stop("A computation peer changed the cross-owner Gaussian exact stage",
           call. = FALSE)
    }
    value
  }), peers)
  consensus <- c(
    "context_hash", "plan_id", "ring_bits", "frac_bits", "backend",
    "bound_x", "bound_y", "max_chunk", "total_n", "numeric_policy_id",
    "version", "state", "producer", "purpose", "capsule_id",
    "analysis_id", "stage", "stage_index", "artifact_sha256",
    "source_contract_hash", "transcript_sha256",
    "numeric_certificate_sha256", "exact_intermediates_exposed",
    "source_values_exposed")
  if (length(unique(vapply(manifests, function(value) {
        .dsvert_joint_dp_client_json(value[consensus])
      }, character(1L)))) != 1L) {
    stop("The computation peers disagree on the exact Gaussian stage",
         call. = FALSE)
  }
  if (identical(manifests[[1L]]$state, "complete")) {
    stop("A completed exact Gaussian stage cannot be claimed twice",
         call. = FALSE)
  }
  manifests
}

.dsvert_dp_gaussian_cross_result_set <- function(
    responses, context, manifest, layout, binding, analysis_id) {
  peers <- context$designated
  artifact <- binding$artifact
  release_layout <- .dsvert_dp_capsule_vector_layout(manifest)
  block <- release_layout$blocks[[paste(
    "gaussian_models", analysis_id, sep = "::")]]
  fields <- c(
    "version", "phase", "capsule_id", "analysis_id", "peer_name",
    "peer_identity_pk", "artifact_sha256", "source_contract_hash",
    "private_layout_sha256", "transcript_sha256",
    "numeric_certificate_sha256", "exact_transcript_sha256",
    "coordinate_count", "release_start", "release_end",
    "release_coordinate_order_sha256", "ring_bits", "frac_bits", "state",
    "fixed_transcript", "result_share_exposed",
    "exact_intermediates_exposed", "alignment_hash_exposed",
    "alignment_hash_exposed_to_relay",
    "alignment_hash_exposed_to_computation_peers", "signature")
  receipts <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_joint_dp_client_decode(
      responses[[peer]], "cross-owner Gaussian result receipt",
      .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_MAX_RECEIPT_BYTES)
    value <- .dsvert_dp_gaussian_cross_verify_signed(
      value, "cross-gaussian-result", peer, context, fields)
    valid <- identical(value$version,
                       .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_RECEIPT_VERSION) &&
      identical(value$phase, "cross_gaussian_result_share_persisted") &&
      identical(value$capsule_id, manifest$capsule_identity$capsule_id) &&
      identical(value$analysis_id, analysis_id) &&
      identical(value$peer_name, peer) &&
      identical(value$peer_identity_pk, unname(context$pinset[[peer]])) &&
      identical(value$artifact_sha256,
                .dsvert_dp_capsule_source_hash(artifact)) &&
      identical(value$source_contract_hash, binding$source_contract_hash) &&
      identical(value$private_layout_sha256,
                layout$transport_coordinate_order_sha256) &&
      identical(value$transcript_sha256,
                .dsvert_dp_capsule_source_hash(artifact$transcript)) &&
      identical(value$numeric_certificate_sha256,
                .dsvert_dp_capsule_source_hash(artifact$numeric_certificate)) &&
      .dsvert_dp_capsule_source_hex(value$exact_transcript_sha256) &&
      identical(as.numeric(value$coordinate_count),
                as.numeric(block$length)) &&
      identical(as.numeric(value$release_start), as.numeric(block$start)) &&
      identical(as.numeric(value$release_end), as.numeric(block$end)) &&
      identical(value$release_coordinate_order_sha256,
                release_layout$sha256) &&
      identical(as.numeric(value$ring_bits), 128) &&
      identical(as.numeric(value$frac_bits),
                as.numeric(artifact$numeric_grid_bits)) &&
      identical(value$state, "complete") &&
      identical(value$fixed_transcript, TRUE) &&
      identical(value$result_share_exposed, FALSE) &&
      identical(value$exact_intermediates_exposed, FALSE) &&
      identical(value$alignment_hash_exposed, FALSE) &&
      identical(value$alignment_hash_exposed_to_relay, FALSE) &&
      identical(value$alignment_hash_exposed_to_computation_peers, FALSE)
    if (!isTRUE(valid)) {
      stop("A computation peer returned an invalid Gaussian result receipt",
           call. = FALSE)
    }
    value
  }), peers)
  comparable <- setdiff(fields, c(
    "peer_name", "peer_identity_pk", "signature"))
  if (length(unique(vapply(receipts, function(value) {
        .dsvert_joint_dp_client_json(value[comparable])
      }, character(1L)))) != 1L) {
    stop("The computation peers disagree on the exact Gaussian result",
         call. = FALSE)
  }
  receipts
}

.dsvert_dp_gaussian_cross_orchestrate <- function(
    manifest_json, manifest, context, source_receipt, .aggregate,
    .setup_exact = .dsvert_setup_exact_gc_transport,
    .vecmul = .dsvert_exact_gc_vecmul_run,
    .alignment_mask = .dsvert_dp_alignment_mask_run,
    .shared_exact = NULL) {
  artifacts <- .dsvert_dp_gaussian_cross_artifacts_client(manifest)
  if (!length(artifacts)) {
    if (!identical(source_receipt$sampler_handoff_ready, TRUE)) {
      stop("The ordinary capsule source is not ready for joint sampling",
           call. = FALSE)
    }
    return(list(enabled = FALSE, sampler_handoff_ready = TRUE))
  }
  if (!is.function(.setup_exact) || !is.function(.vecmul) ||
      !is.function(.alignment_mask)) {
    stop("Invalid cross-owner Gaussian exact transport implementation",
         call. = FALSE)
  }
  layout <- .dsvert_dp_gaussian_cross_layout_client(manifest)
  peers <- context$designated
  expected_source_purpose <- if (length(
      .dsvert_dp_categorical_cross_artifacts_client(manifest))) {
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CATEGORICAL_CROSS_PURPOSE
  } else {
    .DSVERT_CLIENT_DP_CAPSULE_SOURCE_CROSS_PURPOSE
  }
  if (!identical(sort(peers, method = "radix"),
                 unname(unlist(layout$computation_peers))) ||
      !identical(source_receipt$purpose,
                 expected_source_purpose) ||
      !identical(as.numeric(source_receipt$coordinate_count),
                 as.numeric(layout$transport_coordinate_count)) ||
      !identical(as.numeric(source_receipt$release_coordinate_count),
                 as.numeric(layout$release_coordinate_count)) ||
      !identical(source_receipt$private_layout_sha256,
                 layout$transport_coordinate_order_sha256) ||
      !identical(source_receipt$sampler_handoff_ready, FALSE) ||
      !identical(source_receipt$payload_exposed, FALSE)) {
    stop("The cross-owner Gaussian source transport is misbound",
         call. = FALSE)
  }
  if (is.null(.shared_exact)) {
    session_id <- .dsvert_uuid4()
    setup_result <- .dsvert_dp_cross_exact_setup(
      .setup_exact,
      context$all_conns, context$servers,
      match(peers, context$servers), session_id,
      .aggregate = .aggregate)
    on.exit(.dsvert_dp_cross_exact_cleanup(
      context$conns, session_id, setup_result, .aggregate, .setup_exact),
      add = TRUE)
    alignment <- .alignment_mask(
      manifest_json, context, layout, source_receipt,
      session_id, .aggregate)
    .shared_exact <- .dsvert_dp_cross_shared_exact_build(
      manifest_json, manifest, context, layout, source_receipt,
      session_id, alignment)
  } else {
    .shared_exact <- .dsvert_dp_cross_shared_exact_validate(
      .shared_exact, manifest_json, manifest, context, layout,
      source_receipt)
    session_id <- .shared_exact$session_id
  }
  completed <- list()
  for (analysis_id in names(artifacts)) {
    bind_calls <- stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertDPGaussianCrossBindDS",
      manifest_json = manifest_json, analysis_id = analysis_id,
      session_id = session_id)), peers)
    binding <- .dsvert_dp_gaussian_cross_bind_set(
      .dsvert_fanout_by_site(
        context$conns, bind_calls,
        operation = "cross-owner Gaussian private-input binding",
        .aggregate = .aggregate),
      context, manifest, layout, analysis_id)
    if (!identical(binding$source_contract_hash,
                   source_receipt$contract_hash)) {
      stop("The exact Gaussian binding changed the source contract",
           call. = FALSE)
    }
    if (!identical(binding$state, "complete")) {
      m <- length(binding$artifact$input_variable_order)
      stages <- c(
        lapply(seq_len(m - 1L), function(index) {
          list(stage = "validity", stage_index = as.integer(index))
        }),
        list(list(stage = "masked-values", stage_index = 1L)),
        list(list(stage = "moments", stage_index = 1L)))
      for (stage in stages) {
        contract <- .dsvert_dp_gaussian_cross_stage_contract(
          binding$artifact, manifest$capsule_identity$capsule_id,
          analysis_id, stage$stage, stage$stage_index)
        prepare_calls <- stats::setNames(lapply(peers, function(peer) call(
          name = "dsvertDPGaussianCrossPrepareDS",
          analysis_id = analysis_id, stage = stage$stage,
          stage_index = stage$stage_index,
          session_id = session_id)), peers)
        producer_manifests <- .dsvert_dp_gaussian_cross_stage_set(
          .dsvert_fanout_by_site(
            context$conns, prepare_calls,
            operation = paste("cross-owner Gaussian", stage$stage,
                              "preparation"),
            .aggregate = .aggregate),
          context, manifest, binding, analysis_id,
          stage$stage, stage$stage_index, session_id)
        .vecmul(
          context$all_conns, server_names = context$servers,
          servers = match(peers, context$servers), session_id = session_id,
          total_n = contract$total_n,
          input_manifests = producer_manifests,
          transport_ready = TRUE, .aggregate = .aggregate)
      }
    }
    finalize_calls <- stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertDPGaussianCrossFinalizeDS",
      manifest_json = manifest_json, analysis_id = analysis_id,
      session_id = session_id)), peers)
    completed[[analysis_id]] <- .dsvert_dp_gaussian_cross_result_set(
      .dsvert_fanout_by_site(
        context$conns, finalize_calls,
        operation = "cross-owner Gaussian exact finalization",
        .aggregate = .aggregate),
      context, manifest, layout, binding, analysis_id)
  }
  list(
    enabled = TRUE, sampler_handoff_ready = TRUE,
    exact_intermediates_exposed = FALSE, source_values_exposed = FALSE,
    analyses = names(completed), receipts = completed,
    private_layout_sha256 = layout$transport_coordinate_order_sha256)
}
