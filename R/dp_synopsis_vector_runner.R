# Internal orchestration for the stateless sticky synopsis.  The runner keeps
# all private material on the pinned servers and relays only canonical,
# artifact-bound protocol evidence.

.DSVERT_CLIENT_SYNOPSIS_EXACT_START_RESPONSE_VERSION <-
  "dsvert-stateless-catalog-synopsis-exact-gc-start-response-v1"
.DSVERT_CLIENT_SYNOPSIS_EXACT_START_VERSION <-
  "dsvert-stateless-catalog-synopsis-exact-gc-start-v1"
.DSVERT_CLIENT_SYNOPSIS_EXACT_START_DOMAIN <-
  "dsVert/stateless-catalog-synopsis/exact-gc-start/v1|"

.dsvert_dp_synopsis_effective_cross_v1 <- function(manifest) {
  workload <- if (is.list(manifest)) manifest$workload else NULL
  vertical <- if (is.list(workload)) workload$vertical_crosses else NULL
  zero_coordinates <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == 0
  }
  reserved_only <- function(value) {
    if (is.null(value)) return(FALSE)
    if (!is.list(value)) return(TRUE)
    any(vapply(value, function(entry) {
      if (!is.list(entry) || !identical(
            entry$implementation_state, "reserved_not_materialized")) {
        return(TRUE)
      }
      included <- entry$included_coordinate_count
      !is.null(included) && !zero_coordinates(included)
    }, logical(1L)))
  }

  if (!is.null(vertical)) {
    if (!is.list(vertical)) return(TRUE)
    state <- vertical$implementation_state
    known_unmaterialized <- is.character(state) && length(state) == 1L &&
      !is.na(state) && state %in% c(
        "not_required_by_signed_schema", "reserved_not_materialized")
    if (!is.null(state) && !known_unmaterialized) {
      return(TRUE)
    }
    included <- vertical$included_coordinate_count
    if (!is.null(included) && !zero_coordinates(included)) return(TRUE)
    configured <- vertical$configured_crosses
    if (!is.null(configured) &&
        (!is.list(configured) || length(configured) > 0L)) return(TRUE)
    if (reserved_only(vertical$cross_owner_sets) ||
        reserved_only(vertical$categorical_pair_sets)) return(TRUE)
  }

  families <- if (is.list(workload)) workload$families else NULL
  if (!is.null(families) && !is.list(families)) return(TRUE)
  categorical <- if (is.list(families)) families$categorical_pairs else NULL
  if (!is.null(categorical) && !is.list(categorical)) return(TRUE)
  pair_cross <- if (is.list(categorical)) {
    categorical$cross_artifacts
  } else NULL
  if (!is.null(pair_cross) &&
      (!is.list(pair_cross) || length(pair_cross) > 0L)) return(TRUE)
  gaussian_family <- if (is.list(families)) {
    families$gaussian_models
  } else NULL
  if (!is.null(gaussian_family) && !is.list(gaussian_family)) return(TRUE)
  gaussian <- if (is.list(gaussian_family)) gaussian_family$artifacts else NULL
  if (!is.null(gaussian) && !is.list(gaussian)) return(TRUE)
  if (is.list(gaussian) && any(vapply(gaussian, function(artifact) {
        !is.list(artifact) || !artifact$version %in% c(
          "bounded-normalized-gaussian-sufficient-statistics-v1",
          "bounded-normalized-random-intercept-moments-v1",
          "bounded-normalized-random-intercept-fixed-sufficient-statistics-v2",
          "bounded-normalized-random-intercept-fixed-sufficient-statistics-v3",
          "bounded-binary-random-intercept-likelihood-grid-v1",
          "bounded-negative-binomial-likelihood-grid-v1",
          "bounded-multinomial-likelihood-grid-v1")
      }, logical(1L)))) return(TRUE)
  FALSE
}

.dsvert_dp_synopsis_runner_cross <- function(manifest) {
  .dsvert_dp_synopsis_effective_cross_v1(manifest)
}

.dsvert_dp_synopsis_runner_json_set <- function(
    responses, peers, what,
    maximum_bytes = .DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES) {
  if (!is.list(responses) || is.null(names(responses)) ||
      anyNA(names(responses)) || anyDuplicated(names(responses)) ||
      !setequal(names(responses), peers)) {
    stop("The synopsis ", what, " did not cover its pinned peers",
         call. = FALSE)
  }
  responses <- responses[peers]
  for (peer in peers) {
    .dsvert_dp_synopsis_client_json(
      responses[[peer]], paste(what, "from", peer), maximum_bytes)
  }
  responses
}

.dsvert_dp_synopsis_runner_exact_start_set <- function(
    responses, authorities, trusted, compiled, execution, chunk_index) {
  response_json <- .dsvert_dp_synopsis_runner_json_set(
    responses, authorities, "exact-GC START",
    .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
  expected_offset <- as.integer(
    chunk_index * execution$geometry$public_chunk_coordinates)
  expected_count <- as.integer(min(
    execution$geometry$public_chunk_coordinates,
    execution$geometry$coordinate_count - expected_offset))
  fields <- c(
    "version", "phase", "execution_id", "artifact_key",
    "contract_sha256", "attempt_sha256", "source_contract_sha256",
    "local_authority", "chunk_index", "coordinate_offset",
    "coordinate_count", "backend_selection_sha256",
    "worker_contract_sha256", "binding_sha256", "operation_id",
    "purpose", "local_chunk_durable", "intermediate_payload_exposed",
    "source_share_exposed", "private_seed_exposed",
    "preclamp_values_exposed", "signature")
  decoded <- stats::setNames(lapply(authorities, function(peer) {
    response <- .dsvert_dp_synopsis_client_json(
      response_json[[peer]], paste("exact-GC START from", peer),
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
    if (!.dsvert_dp_has_exact_names(
          response, c("version", "receipt", "initialization")) ||
        !identical(response$version,
                   .DSVERT_CLIENT_SYNOPSIS_EXACT_START_RESPONSE_VERSION) ||
        !is.list(response$receipt) || !is.list(response$initialization)) {
      stop("Invalid exact-GC START response", call. = FALSE)
    }
    receipt <- response$receipt
    role <- match(peer, authorities)
    authority <- list(
      peer_name = peer,
      identity_pk = unname(trusted$context$pinset[[peer]]),
      role = c("primary_noise_authority", "secondary_noise_authority")[[role]])
    flags <- c("intermediate_payload_exposed", "source_share_exposed",
               "private_seed_exposed", "preclamp_values_exposed")
    valid <- .dsvert_dp_has_exact_names(receipt, fields) &&
      identical(receipt$version, .DSVERT_CLIENT_SYNOPSIS_EXACT_START_VERSION) &&
      identical(receipt$phase, "synopsis_exact_gc_initialized") &&
      identical(receipt$execution_id, execution$execution_id) &&
      identical(receipt$artifact_key, compiled$artifact$artifact_key) &&
      identical(
        .dsvert_joint_dp_client_json(receipt$local_authority),
        .dsvert_joint_dp_client_json(authority)) &&
      identical(as.numeric(receipt$chunk_index), as.numeric(chunk_index)) &&
      identical(as.numeric(receipt$coordinate_offset),
                as.numeric(expected_offset)) &&
      identical(as.numeric(receipt$coordinate_count),
                as.numeric(expected_count)) &&
      all(vapply(receipt[c(
        "contract_sha256", "attempt_sha256", "source_contract_sha256",
        "backend_selection_sha256", "worker_contract_sha256",
        "binding_sha256")], .dsvert_vector_hex, logical(1L))) &&
      is.character(receipt$operation_id) && length(receipt$operation_id) == 1L &&
      !is.na(receipt$operation_id) &&
      grepl("^op_[0-9a-f]{32}$", receipt$operation_id) &&
      is.character(receipt$purpose) && length(receipt$purpose) == 1L &&
      !is.na(receipt$purpose) &&
      grepl("^[a-z][a-z0-9_.:/-]*$", receipt$purpose) &&
      identical(receipt$local_chunk_durable, FALSE) &&
      all(vapply(flags, function(field) identical(receipt[[field]], FALSE),
                 logical(1L)))
    if (!isTRUE(valid)) stop("Invalid or misbound exact-GC START receipt",
                             call. = FALSE)
    unsigned <- receipt[setdiff(names(receipt), "signature")]
    .dsvert_dp_synopsis_client_signature(
      unsigned, receipt$signature, authority$identity_pk,
      .DSVERT_CLIENT_SYNOPSIS_EXACT_START_DOMAIN, "exact-GC START")
    initialization <- response$initialization
    state_ok <- is.character(initialization$state) &&
      length(initialization$state) == 1L &&
      !is.na(initialization$state) &&
      initialization$state %in% c("running", "complete") &&
      is.logical(initialization$stored) && length(initialization$stored) == 1L &&
      !is.na(initialization$stored) &&
      identical(initialization$stored,
                identical(initialization$state, "complete"))
    if (!isTRUE(state_ok)) {
      stop("Invalid exact-GC START initialization", call. = FALSE)
    }
    list(receipt = receipt, initialization = initialization)
  }), authorities)
  common <- setdiff(fields, c("local_authority", "signature"))
  if (!identical(
        .dsvert_joint_dp_client_json(decoded[[1L]]$receipt[common]),
        .dsvert_joint_dp_client_json(decoded[[2L]]$receipt[common]))) {
    stop("The exact-GC START receipts disagree", call. = FALSE)
  }
  states <- vapply(decoded, function(value) value$initialization$state,
                   character(1L))
  if (length(unique(states)) != 1L) {
    stop("The exact-GC START peers disagree on liveness", call. = FALSE)
  }
  list(
    receipts = lapply(decoded, `[[`, "receipt"),
    initialized = lapply(decoded, `[[`, "initialization"),
    complete = identical(unname(states[[1L]]), "complete"))
}

.dsvert_dp_synopsis_runner_count <- function(value, what) {
  if (!.dsvert_dp_is_integer(value, 1L, .DSVERT_DP_MAX_COORDINATES)) {
    stop("Invalid synopsis ", what, call. = FALSE)
  }
  as.integer(value)
}

.dsvert_dp_synopsis_runner_compile <- function(
    context, manifest_bundle, trusted, layout, .aggregate) {
  servers <- context$servers
  owners <- tryCatch(vapply(
    layout$blocks, `[[`, character(1L), "owner_peer"),
    error = function(error) character())
  artifacts <- .dsvert_dp_categorical_cross_artifacts_client(
    trusted$manifest)
  if (length(artifacts)) {
    if (length(artifacts) != 1L) {
      stop("The categorical Synopsis source-owner binding is ambiguous",
           call. = FALSE)
    }
    cross_owners <- tryCatch(unname(unlist(lapply(
      artifacts, function(artifact) c(
        artifact$left$owner_peer, artifact$right$owner_peer)),
      use.names = FALSE)), error = function(error) character())
    if (length(cross_owners) != 2L || anyNA(cross_owners) ||
        any(!nzchar(cross_owners)) || anyDuplicated(cross_owners)) {
      stop("The categorical Synopsis source-owner binding is invalid",
           call. = FALSE)
    }
    owners <- c(owners, cross_owners)
  }
  sources <- sort(unique(owners), method = "radix")
  if (!length(sources) || anyNA(sources) || any(!nzchar(sources)) ||
      !all(sources %in% servers)) {
    stop("The synopsis source-peer assignment is invalid", call. = FALSE)
  }

  claim_calls <- stats::setNames(lapply(sources, function(peer) call(
    name = "dsvertDPSynopsisClaimDS",
    manifest_sha256 = manifest_bundle$manifest_sha256)), sources)
  claim_raw <- .dsvert_fanout_by_site(
    context$all_conns[sources], claim_calls,
    operation = "stateless synopsis source Claim",
    .aggregate = .aggregate)
  claims <- stats::setNames(lapply(sources, function(peer) {
    value <- .dsvert_dp_synopsis_client_json(
      claim_raw[[peer]], paste("synopsis Claim from", peer))
    if (!.dsvert_dp_has_exact_names(
          value, c("version", "projection", "claim")) ||
        !identical(
          value$version,
          "dsvert-stateless-catalog-synopsis-local-source-vector-claim-v1") ||
        !is.list(value$claim) ||
        !identical(value$claim$source_peer_name, peer)) {
      stop("Invalid local synopsis Claim", call. = FALSE)
    }
    value
  }), sources)
  projection_json <- vapply(claims, function(value) {
    .dsvert_joint_dp_client_json(value$projection)
  }, character(1L))
  if (length(unique(projection_json)) != 1L) {
    stop("The local synopsis Claims disagree on their catalog",
         call. = FALSE)
  }
  claims_json <- .dsvert_joint_dp_client_json(unname(lapply(
    claims, `[[`, "claim")))

  compile_calls <- stats::setNames(lapply(servers, function(peer) call(
    name = "dsvertDPSynopsisCompileDS",
    manifest_sha256 = manifest_bundle$manifest_sha256,
    claims_json = claims_json)), servers)
  compile_raw <- .dsvert_fanout_by_site(
    context$all_conns, compile_calls,
    operation = "stateless synopsis K-peer compilation",
    .aggregate = .aggregate)
  envelopes <- stats::setNames(lapply(servers, function(peer) {
    value <- .dsvert_dp_synopsis_client_json(
      compile_raw[[peer]], paste("synopsis compilation from", peer))
    if (!.dsvert_dp_has_exact_names(
          value, c("version", "claim_set", "artifact", "receipt")) ||
        !identical(
          value$version,
          "dsvert-stateless-catalog-synopsis-local-compile-envelope-v1") ||
        !is.list(value$claim_set) || !is.list(value$artifact) ||
        !is.list(value$receipt) ||
        !identical(value$receipt$peer_name, peer)) {
      stop("Invalid local synopsis compilation", call. = FALSE)
    }
    value
  }), servers)
  reference <- envelopes[[servers[[1L]]]]
  if (!all(vapply(envelopes, function(value) {
        identical(
          .dsvert_joint_dp_client_json(value$claim_set),
          .dsvert_joint_dp_client_json(reference$claim_set)) &&
          identical(
            .dsvert_joint_dp_client_json(value$artifact),
            .dsvert_joint_dp_client_json(reference$artifact))
      }, logical(1L))) ||
      !identical(reference$artifact$semantic$source_claim_set_sha256,
                 reference$claim_set$sha256)) {
    stop("The K synopsis compilations disagree", call. = FALSE)
  }
  receipts <- stats::setNames(lapply(
    envelopes[servers], `[[`, "receipt"), servers)
  unsigned <- unname(lapply(receipts, function(receipt) {
    receipt[setdiff(names(receipt), "signature")]
  }))
  compilation <- list(
    version = .DSVERT_CLIENT_SYNOPSIS_COMPILE_SET_VERSION,
    artifact = reference$artifact, receipts = receipts,
    receipt_set_sha256 = .dsvert_dp_synopsis_client_hash(
      "dsVert/stateless-catalog-synopsis/compile-receipt-set/v1|",
      list(
        version = .DSVERT_CLIENT_SYNOPSIS_COMPILE_SET_VERSION,
        receipts = unsigned)))
  list(
    claim_set = reference$claim_set, compilation = compilation,
    compiled = .dsvert_dp_synopsis_client_compile(
      compilation, trusted, manifest_bundle), sources = sources)
}

.dsvert_dp_synopsis_vector_run <- function(
    datasources, status = NULL, local_projection = NULL,
    .aggregate = DSI::datashield.aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  bootstrap <- if (is.null(local_projection)) {
    .dsvert_dp_synopsis_bootstrap_build_v1(
      datasources, status = status, .aggregate = .aggregate)
  } else {
    .dsvert_dp_synopsis_bootstrap_build_v1(
      datasources, status = status, local_projection = local_projection,
      .aggregate = .aggregate)
  }
  status <- bootstrap$status
  manifest_bundle <- bootstrap$manifest_bundle
  trusted <- .dsvert_dp_synopsis_client_bundle(manifest_bundle, status)
  cross <- .dsvert_dp_synopsis_supported_categorical_cross_v1(
    trusted$manifest)
  if (.dsvert_dp_synopsis_runner_cross(trusted$manifest) &&
      !isTRUE(cross)) {
    stop("Cross-owner synopsis catalogs are not supported", call. = FALSE)
  }
  published <- .dsvert_dp_synopsis_publication_resume_v1(
    bootstrap, .aggregate = .aggregate)
  if (!is.null(published)) return(published)
  layout <- .dsvert_dp_capsule_vector_layout(trusted$manifest)
  context <- trusted$context
  valid_context <- is.list(context) && is.character(context$servers) &&
    length(context$servers) >= 2L &&
    identical(context$servers, sort(context$servers, method = "radix")) &&
    is.character(context$designated) && length(context$designated) == 2L &&
    !anyDuplicated(context$designated) &&
    all(context$designated %in% context$servers) &&
    is.list(context$all_conns) &&
    identical(names(context$all_conns), context$servers) &&
    is.list(context$conns) &&
    identical(names(context$conns), context$designated)
  if (!isTRUE(valid_context)) {
    stop("Invalid trusted synopsis peer context", call. = FALSE)
  }
  servers <- context$servers
  authorities <- context$designated
  authority_indices <- match(authorities, servers)
  if (length(authority_indices) != 2L || anyNA(authority_indices) ||
      anyDuplicated(authority_indices)) {
    stop("The synopsis authorities are absent from the pinned federation",
         call. = FALSE)
  }
  mechanism <- trusted$manifest$workload$capsule_mechanism
  backend <- if (identical(
      mechanism$mechanism, .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM)) {
    NULL
  } else if (layout$coordinate_count <=
             .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_MAX_PROMOTED_COORDINATES) {
    .DSVERT_CLIENT_VECTOR_EXACT_BACKEND
  } else {
    .DSVERT_CLIENT_VECTOR_BACKEND
  }
  cold_profile <- .dsvert_vector_profile(
    mechanism,
    mechanism_selection = trusted$manifest$workload$mechanism_selection,
    backend = backend)
  built <- .dsvert_dp_synopsis_runner_compile(
    context, manifest_bundle, trusted, layout, .aggregate)
  compiled <- built$compiled
  execution <- .dsvert_dp_synopsis_client_execution(compiled)
  public_chunk_count <- .dsvert_dp_synopsis_runner_count(
    execution$geometry$public_chunk_count, "public chunk count")
  claim_set_json <- .dsvert_joint_dp_client_json(built$claim_set)
  compilation_json <- .dsvert_joint_dp_client_json(built$compilation)
  session_id <- .dsvert_uuid4()

  prepare_calls <- stats::setNames(lapply(authorities, function(peer) call(
    name = "dsvertDPSynopsisPrepareDS", session_id = session_id,
    manifest_sha256 = manifest_bundle$manifest_sha256,
    claim_set_json = claim_set_json,
    compilation_json = compilation_json)), authorities)
  prepare_json <- .dsvert_dp_synopsis_runner_json_set(
    .dsvert_fanout_by_site(
      context$conns, prepare_calls,
      operation = "stateless synopsis PREPARE",
      .aggregate = .aggregate),
    authorities, "PREPARE",
    .DSVERT_CLIENT_SYNOPSIS_PREPARE_MAX_OBJECT_BYTES)

  result_calls <- stats::setNames(lapply(authorities, function(peer) call(
    name = "dsvertDPSynopsisResultDS", session_id = session_id,
    first_prepare_json = prepare_json[[authorities[[1L]]]],
    second_prepare_json = prepare_json[[authorities[[2L]]]])), authorities)
  # A first RESULT is deliberately not probed over DSI.  The server's
  # phase-not-ready condition is serialized by DataSHIELD as an opaque
  # transport failure, so a probe cannot distinguish the normal pre-START
  # state from a real failure.  START and RESULT are already durable,
  # idempotent transitions; replaying START is the safe recovery path.
  existing_result <- list(ok = FALSE, value = NULL)
  if (isTRUE(existing_result$ok)) {
    result_json <- .dsvert_dp_synopsis_runner_json_set(
      existing_result$value, authorities, "RESULT",
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
  } else {
    exact_transport_ready <- FALSE
    if (isTRUE(cold_profile$exact_gc)) {
      .dsvert_setup_exact_gc_transport(
        context$all_conns, servers, authority_indices, session_id,
        .aggregate = .aggregate)
      exact_transport_ready <- TRUE
    }
    source_manifest <- .dsvert_dp_synopsis_source_manifest_v1(
      trusted$manifest, context)
    source_binding <- list(
      version = "dsvert-stateless-catalog-synopsis-source-contract-v1",
      manifest_capsule_id = manifest_bundle$capsule_id,
      artifact_key = compiled$artifact$artifact_key,
      source_claim_set_sha256 = built$claim_set$sha256)
    source_capsule_id <- .dsvert_dp_synopsis_client_hash(
      "dsVert/stateless-catalog-synopsis/source-namespace/v1|",
      source_binding)
    if (!source_manifest$capsule_id %in%
        c(manifest_bundle$capsule_id, source_capsule_id)) {
      stop("The synopsis source namespace is detached", call. = FALSE)
    }
    source_manifest$capsule_id <- source_capsule_id
    valid_source_geometry <- if (isTRUE(cross)) {
      isTRUE(source_manifest$cross_enabled) &&
        identical(as.numeric(source_manifest$release_coordinate_count),
                  as.numeric(compiled$layout$coordinate_count))
    } else {
      !isTRUE(source_manifest$cross_enabled) &&
        identical(as.numeric(source_manifest$coordinate_count),
                  as.numeric(compiled$layout$coordinate_count))
    }
    if (!isTRUE(valid_source_geometry)) {
      stop("The synopsis source geometry is invalid", call. = FALSE)
    }

    ticket_calls <- stats::setNames(lapply(authorities, function(peer) call(
      name = "dsvertDPSynopsisSourceTicketDS",
      manifest_sha256 = manifest_bundle$manifest_sha256,
      claim_set_json = claim_set_json,
      compilation_json = compilation_json)), authorities)
    ticket_raw <- .dsvert_fanout_by_site(
      context$conns, ticket_calls,
      operation = "stateless synopsis source tickets",
      .aggregate = .aggregate)
    tickets <- .dsvert_dp_capsule_source_ticket_set(
      ticket_raw, context, source_manifest)
    source_peers <- unname(unlist(
      tickets[[authorities[[1L]]]]$value$source_peers,
      use.names = FALSE))
    if (!identical(source_peers, built$sources)) {
      stop("The synopsis source tickets changed their source set",
           call. = FALSE)
    }
    source_prepare_calls <- stats::setNames(lapply(
      source_peers, function(source) call(
        name = "dsvertDPSynopsisSourcePrepareDS",
        manifest_sha256 = manifest_bundle$manifest_sha256,
        claim_set_json = claim_set_json,
        compilation_json = compilation_json,
        first_ticket_json = tickets[[authorities[[1L]]]]$json,
        second_ticket_json = tickets[[authorities[[2L]]]]$json)),
      source_peers)
    source_prepare_raw <- .dsvert_fanout_by_site(
      context$all_conns[source_peers], source_prepare_calls,
      operation = "stateless synopsis source preparation",
      .aggregate = .aggregate)
    summaries <- stats::setNames(lapply(source_peers, function(source) {
      .dsvert_dp_capsule_source_summary_response(
        source_prepare_raw[[source]], source, context,
        tickets[[authorities[[1L]]]], tickets)
    }), source_peers)
    chunk_count <- .dsvert_dp_synopsis_runner_count(
      tickets[[authorities[[1L]]]]$value$chunk_count,
      "source chunk count")
    final_acks <- NULL
    for (source_index in seq_along(source_peers)) {
      source <- source_peers[[source_index]]
      for (chunk_index in seq.int(0L, chunk_count - 1L)) {
        chunk_call <- stats::setNames(list(call(
          name = "dsvertDPSynopsisSourceChunkDS",
          manifest_sha256 = manifest_bundle$manifest_sha256,
          claim_set_json = claim_set_json,
          compilation_json = compilation_json,
          source_transfer_id = summaries[[source]]$value$source_transfer_id,
          chunk_index = as.integer(chunk_index))), source)
        fetched <- .dsvert_fanout_by_site(
          context$all_conns[source], chunk_call,
          operation = "stateless synopsis source chunk",
          .aggregate = .aggregate)
        envelope_batches <- .dsvert_dp_capsule_source_bundle_window(
          fetched[[source]], source, context, tickets, summaries[[source]],
          as.integer(chunk_index))
        if (length(envelope_batches) != 1L) {
          stop("The synopsis source returned an invalid chunk window",
               call. = FALSE)
        }
        accept_calls <- stats::setNames(lapply(authorities, function(peer) {
          call(
            name = "dsvertDPSynopsisSourceAcceptDS",
            manifest_sha256 = manifest_bundle$manifest_sha256,
            envelope_json = envelope_batches[[1L]][[peer]]$json)
        }), authorities)
        accepted <- .dsvert_fanout_by_site(
          context$conns, accept_calls,
          operation = "stateless synopsis source acceptance",
          .aggregate = .aggregate)
        final_acks <- stats::setNames(lapply(authorities, function(peer) {
          tail(.dsvert_dp_capsule_source_ack_window(
            accepted[[peer]], peer, source, context, tickets[[peer]],
            summaries[[source]], envelope_batches, source_index,
            length(source_peers)), 1L)[[1L]]
        }), authorities)
      }
    }
    if (is.null(final_acks) || !all(vapply(
          final_acks, `[[`, logical(1L),
          "capsule_aggregation_complete"))) {
      stop("The synopsis source aggregation did not complete",
           call. = FALSE)
    }

    if (isTRUE(cross)) {
      source_receipt <-
        .dsvert_dp_synopsis_categorical_cross_source_receipt_v1(
          source_manifest, tickets, authorities, source_peers)
      cross_manifest <- trusted$manifest
      remote_context <- list(
        manifest_sha256 = manifest_bundle$manifest_sha256,
        claim_set_json = claim_set_json,
        compilation_json = compilation_json)
      cross_receipt <- .dsvert_dp_categorical_cross_orchestrate(
        manifest_bundle$manifest_json, cross_manifest, context,
        source_receipt, .aggregate, .remote_context = remote_context)
      if (!is.list(cross_receipt) ||
          !identical(cross_receipt$enabled, TRUE) ||
          !identical(cross_receipt$sampler_handoff_ready, TRUE) ||
          !identical(cross_receipt$exact_intermediates_exposed, FALSE) ||
          !identical(cross_receipt$source_values_exposed, FALSE)) {
        stop("The exact categorical Synopsis source is not ready for sampling",
             call. = FALSE)
      }
    }

    for (chunk_index in seq.int(0L, public_chunk_count - 1L)) {
      start_calls <- stats::setNames(lapply(authorities, function(peer) call(
        name = "dsvertDPSynopsisStartDS", session_id = session_id,
        first_prepare_json = prepare_json[[authorities[[1L]]]],
        second_prepare_json = prepare_json[[authorities[[2L]]]],
        chunk_index = as.integer(chunk_index))), authorities)
      started_raw <- .dsvert_fanout_by_site(
        context$conns, start_calls,
        operation = "stateless synopsis START",
        .aggregate = .aggregate)
      if (isTRUE(cold_profile$exact_gc)) {
        started <- .dsvert_dp_synopsis_runner_exact_start_set(
          started_raw, authorities, trusted, compiled, execution,
          chunk_index)
        if (!isTRUE(started$complete)) {
          receipt <- started$receipts[[1L]]
          suffix <- sub("^op_", "", receipt$operation_id)
          .dsvert_exact_gc_run(
            datasources = context$all_conns, server_names = servers,
            servers = authority_indices, session_id = session_id,
            operation_id = receipt$operation_id,
            source_key = paste0("exact_gc_in_", suffix),
            output_key = paste0("exact_gc_out_", suffix),
            operation = "joint-dp-vector-laplace-v3", ring = 128L,
            frac_bits = 0L, vector_len = receipt$coordinate_count,
            purpose = receipt$purpose, transport_ready = exact_transport_ready,
            initialized = started$initialized, .aggregate = .aggregate)
        }
      } else {
        .dsvert_dp_synopsis_runner_json_set(
          started_raw, authorities, "START",
          .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
      }
    }
    result_json <- .dsvert_dp_synopsis_runner_json_set(
      .dsvert_fanout_by_site(
        context$conns, result_calls,
        operation = "stateless synopsis RESULT commit",
        .aggregate = .aggregate),
      authorities, "RESULT",
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
  }

  release_calls <- stats::setNames(lapply(authorities, function(peer) call(
    name = "dsvertDPSynopsisReleaseDS", session_id = session_id,
    first_result_json = result_json[[authorities[[1L]]]],
    second_result_json = result_json[[authorities[[2L]]]])), authorities)
  existing_release <- .dsvert_vector_try_phase(.dsvert_fanout_by_site(
    context$conns, release_calls,
    operation = "stateless synopsis RELEASE lookup",
    .aggregate = .aggregate))
  if (isTRUE(existing_release$ok)) {
    release_json <- .dsvert_dp_synopsis_runner_json_set(
      existing_release$value, authorities, "RELEASE",
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
  } else {
    if (!exists("exact_transport_ready", inherits = FALSE) ||
        !isTRUE(exact_transport_ready)) {
      .dsvert_setup_exact_gc_transport(
        context$all_conns, servers, authority_indices, session_id,
        .aggregate = .aggregate)
    }
    for (chunk_index in seq.int(0L, public_chunk_count - 1L)) {
      share_calls <- stats::setNames(lapply(authorities, function(peer) call(
        name = "dsvertDPSynopsisFinalShareDS", session_id = session_id,
        first_result_json = result_json[[authorities[[1L]]]],
        second_result_json = result_json[[authorities[[2L]]]],
        public_chunk_index = as.integer(chunk_index))), authorities)
      share_raw <- .dsvert_dp_synopsis_runner_json_set(
        .dsvert_fanout_by_site(
          context$conns, share_calls,
          operation = "stateless synopsis FINAL_SHARE",
          .aggregate = .aggregate),
        authorities, "FINAL_SHARE",
        .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
      for (sender in authorities) {
        recipient <- setdiff(authorities, sender)
        share <- .dsvert_dp_synopsis_client_json(
          share_raw[[sender]], paste("synopsis FINAL_SHARE from", sender),
          .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
        if (!is.list(share$transfer) ||
            !identical(share$transfer$sender_name, sender) ||
            !identical(share$transfer$recipient_name, recipient) ||
            !.dsvert_dp_is_string(share$ciphertext)) {
          stop("Invalid synopsis FINAL_SHARE transfer", call. = FALSE)
        }
        .dsvert_store_typed_blob(
          share$ciphertext, share$transfer,
          context$all_conns[recipient], session_id,
          producer_conn = context$all_conns[sender],
          .aggregate = .aggregate)
      }
    }
    release_json <- .dsvert_dp_synopsis_runner_json_set(
      .dsvert_fanout_by_site(
        context$conns, release_calls,
        operation = "stateless synopsis RELEASE",
        .aggregate = .aggregate),
      authorities, "RELEASE",
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
  }

  published <- .dsvert_dp_synopsis_publication_resume_v1(
    bootstrap, .aggregate = .aggregate)
  if (is.null(published)) {
    stop("The durable synopsis publication was not available after RELEASE",
         call. = FALSE)
  }
  published
}
