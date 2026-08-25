.capsule_manifest_client_b64url <- function(value) {
  sub("=+$", "", chartr(
    "+/", "-_", gsub("[\r\n]", "", jsonlite::base64_enc(value))),
    perl = TRUE)
}

.capsule_manifest_client_fixture <- function(k = 2L) {
  stopifnot(k %in% 2:5)
  peers <- paste0("peer_", letters[seq_len(k)])
  keys <- stats::setNames(lapply(peers, function(peer) {
    openssl::ed25519_keygen()
  }), peers)
  pins <- vapply(keys, function(key) {
    .capsule_manifest_client_b64url(
      tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  pin_hash <- .dsvert_dp_pinset_hash(pins)
  statuses <- stats::setNames(lapply(seq_along(peers), function(index) {
    peer <- peers[[index]]
    list(
    policy = list(
      domain = "manifest-client-study",
      cohort_id = "manifest-client-cohort",
      peer_pinset_sha256 = pin_hash),
    noise_root = list(
      privacy_epoch = as.numeric(index), key_id = paste0("noise-", peer)))
  }), peers)
  context <- list(
    servers = peers, pinset = pins, status = statuses,
    all_conns = stats::setNames(lapply(peers, function(peer) {
      structure(list(), class = "manifest-test-connection")
    }), peers))

  sign <- function(message, peer) {
    .capsule_manifest_client_b64url(
      openssl::ed25519_sign(message, keys[[peer]]))
  }
  sign_object <- function(unsigned, peer, domain) {
    c(unsigned, list(signature = sign(
      .dsvert_dp_capsule_manifest_message(domain, unsigned), peer)))
  }
  drafts <- stats::setNames(lapply(seq_along(peers), function(index) {
    peer <- peers[[index]]
    data_name <- paste0("data_", peer)
    columns <- list()
    fragments <- list(
      describe = list(), survival = list(), gaussian = list(),
      vertical_cross = list())
    if (index == 1L) {
      columns$x_peer_a <- list(
        kind = "numeric", owner_peer = peer, lower = 0, upper = 10)
      fragments$describe$owner_a_describe <- list(
        version = "v1", dataset = data_name, variables = "x_peer_a",
        histogram_grids = list(x_peer_a = c(5, 10)),
        allocation = list(
          count = 0.25, histogram = 0.25, sum = 0.25, sumsq = 0.25))
      fragments$vertical_cross$owner_a_cross <- list(
        version = "v1", left_dataset = "data_peer_a",
        right_dataset = "data_peer_b", left = "x_peer_a",
        right = "time_peer_b", family = "numeric_cross_moment")
    } else if (index == 2L) {
      columns$status_peer_b <- list(
        kind = "categorical", owner_peer = peer,
        levels = c("0", "1"))
      columns$time_peer_b <- list(
        kind = "numeric", owner_peer = peer, lower = 0, upper = 10)
      fragments$survival$owner_b_survival <- list(
        version = "v1", dataset = data_name,
        time = "time_peer_b", event = "status_peer_b", censor = "0",
        time_grid = c(5, 10), entry = NULL)
    } else {
      columns[[paste0("x_", peer)]] <- list(
        kind = "numeric", owner_peer = peer, lower = -1, upper = 1)
    }
    unsigned <- list(
      version = .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_DRAFT_VERSION,
      phase = "custodian_policy_draft", peer_name = peer,
      peer_identity_pk = unname(pins[[peer]]),
      peer_pinset_sha256 = pin_hash,
      domain = "manifest-client-study",
      cohort_id = "manifest-client-cohort",
      dataset_mapping_mode = "automatic_single_local_dataset",
      datasets = stats::setNames(list(list(
        dataset_id = paste0("cohort-", peer), dataset_version = "v1",
        schema_version = .DSVERT_CLIENT_DP_CAPSULE_POLICY_SCHEMA_VERSION,
        alignment_group = "manifest-client-cohort",
        alignment_protocol_version = 1,
        patient_column = "patient_id", columns = columns)), data_name),
      workload_fragments = fragments,
      data_access = FALSE, patient_derived_metadata = FALSE,
      operation_limit = FALSE, request_limit = FALSE,
      history_can_deny_operation = FALSE)
    .dsvert_joint_dp_client_json(sign_object(unsigned, peer, "draft"))
  }), peers)

  state <- new.env(parent = emptyenv())
  state$calls <- list()
  state$different_build <- FALSE
  state$drafts <- drafts
  artifact_index <- function(manifest, policy, manifest_sha256) {
    context <- list(
      version = .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_VERSION,
      manifest_sha256 = manifest_sha256,
      privacy_epoch_scope = "per_peer_signed_receipts_v1")
    value <- .dsvert_joint_dp_client_canonical(list(
      version = .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_VERSION,
      context = context, gaussian_models = list()))
    list(
      value = value, count = 0,
      root = .dsvert_dp_capsule_artifact_merkle_root_client(
        character(), context))
  }
  fanout <- function(conns, calls, operation, .aggregate) {
    state$calls[[length(state$calls) + 1L]] <- calls
    endpoint <- as.character(as.list(calls[[1L]])[[1L]])
    if (identical(endpoint, "dsvertDPCapsuleManifestDraftDS")) {
      return(state$drafts[names(calls)])
    }
    if (identical(endpoint, "dsvertDPCapsuleManifestSignDS")) {
      return(stats::setNames(lapply(names(calls), function(peer) {
        arguments <- as.list(calls[[peer]])
        schema <- .dsvert_joint_dp_client_decode(
          arguments$schema_json, "test unsigned schema", 8L * 1024L^2)
        schema_hash <- digest::digest(
          arguments$schema_json, algo = "sha256", serialize = FALSE)
        workload_hash <- digest::digest(
          arguments$workload_contract_json, algo = "sha256",
          serialize = FALSE)
        schema_signature <- sign(charToRaw(paste0(
          .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
          .dsvert_joint_dp_client_json(schema))), peer)
        .dsvert_joint_dp_client_json(list(
          version = .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_SIGN_VERSION,
          phase = "global_schema_policy_verified", peer_name = peer,
          peer_identity_pk = unname(pins[[peer]]),
          peer_pinset_sha256 = pin_hash,
          schema_sha256 = schema_hash,
          workload_contract_sha256 = workload_hash,
          logical_snapshot = schema$logical_snapshot,
          schema_signature = schema_signature, data_access = FALSE,
          operation_limit = FALSE, request_limit = FALSE,
          history_can_deny_operation = FALSE))
      }), names(calls)))
    }
    if (!identical(endpoint, "dsvertDPCapsuleManifestBuildDS")) {
      stop("unexpected manifest endpoint")
    }
    stats::setNames(lapply(names(calls), function(peer) {
      arguments <- as.list(calls[[peer]])
      signed_schema <- .dsvert_joint_dp_client_decode(
        arguments$schema_json, "test signed schema", 8L * 1024L^2)
      unsigned_schema <- signed_schema[
        setdiff(names(signed_schema), "signatures")]
      schema_json <- .dsvert_joint_dp_client_json(unsigned_schema)
      schema_hash <- digest::digest(
        schema_json, algo = "sha256", serialize = FALSE)
      workload_hash <- digest::digest(
        arguments$workload_contract_json, algo = "sha256",
        serialize = FALSE)
      manifest_json <- .dsvert_joint_dp_client_json(list(
        protocol = "test-server-authoritative-manifest-v1",
        schema_sha256 = schema_hash,
        workload_contract_sha256 = workload_hash,
        variant = if (isTRUE(state$different_build) &&
                      identical(peer, tail(peers, 1L))) 2 else 1))
      capsule_id <- digest::digest(
        manifest_json, algo = "sha256", serialize = FALSE)
      index <- artifact_index(
        NULL, NULL,
        digest::digest(manifest_json, algo = "sha256", serialize = FALSE))
      unsigned <- list(
        version = .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_BUILD_VERSION,
        phase = "server_authoritative_manifest_memoized",
        peer_name = peer, peer_identity_pk = unname(pins[[peer]]),
        peer_pinset_sha256 = pin_hash, schema_sha256 = schema_hash,
        workload_contract_sha256 = workload_hash,
        manifest_sha256 = digest::digest(
          manifest_json, algo = "sha256", serialize = FALSE),
        manifest_bytes = nchar(manifest_json, type = "bytes"),
        capsule_id = capsule_id,
        privacy_epoch = statuses[[peer]]$noise_root$privacy_epoch,
        noise_key_id = paste0("noise-", peer),
        artifact_commitment_count = index$count,
        artifact_commitments_root = index$root,
        artifact_commitments = index$value,
        manifest_json = manifest_json,
        durable_memoization = TRUE, deterministic_replay = TRUE,
        data_access = FALSE, operation_limit = FALSE,
        request_limit = FALSE, history_can_deny_operation = FALSE)
      .dsvert_joint_dp_client_json(sign_object(unsigned, peer, "build"))
    }), names(calls))
  }
  source_manifest <- function(manifest_json, context) {
    value <- .dsvert_joint_dp_client_decode(
      manifest_json, "test authoritative manifest", 32L * 1024L^2)
    schema_call <- state$calls[[length(state$calls)]][[1L]]
    schema <- .dsvert_joint_dp_client_decode(
      as.list(schema_call)$schema_json, "test signed schema", 8L * 1024L^2)
    unsigned <- schema[setdiff(names(schema), "signatures")]
    list(
      capsule_id = digest::digest(
        manifest_json, algo = "sha256", serialize = FALSE),
      value = list(
        logical_snapshot = unsigned$logical_snapshot,
        workload = list(schema_attestation = list(
          manifest_sha256 = value$schema_sha256,
          signers = as.list(sort(names(schema$signatures),
                                 method = "radix")),
          signatures = schema$signatures))))
  }
  list(
    peers = peers, keys = keys, pins = pins, context = context,
    drafts = drafts, state = state, fanout = fanout,
    source_manifest = source_manifest, artifact_index = artifact_index)
}

.capsule_manifest_client_run <- function(fixture) {
  testthat::with_mocked_bindings(
    .dsvert_dp_capsule_manifest_build(
      fixture$context$all_conns, .aggregate = function(...) NULL),
    .dsvert_joint_dp_client_context = function(...) fixture$context,
    .dsvert_fanout_by_site = fixture$fanout,
    .dsvert_dp_capsule_source_manifest = fixture$source_manifest,
    .dsvert_dp_capsule_artifact_commitment_index_client =
      fixture$artifact_index,
    .package = "dsVertClient")
}

.capsule_manifest_client_recursive_names <- function(value) {
  if (!is.list(value)) return(character())
  c(names(value), unlist(lapply(
    value, .capsule_manifest_client_recursive_names), use.names = FALSE))
}

test_that("K=2 through K=5 relay the signed owner union without analyst metadata", {
  for (k in 2:5) {
    fixture <- .capsule_manifest_client_fixture(k)
    first <- .capsule_manifest_client_run(fixture)
    first_calls <- fixture$state$calls
    second <- .capsule_manifest_client_run(fixture)
    second_calls <- fixture$state$calls[seq_along(first_calls) +
                                         length(first_calls)]

    expect_identical(first$manifest_json, second$manifest_json)
    expect_identical(first$schema_json, second$schema_json)
    expect_identical(
      first$workload_contract_json, second$workload_contract_json)
    expect_setequal(names(first$manifest_signatures), fixture$peers)
    expect_true(all(vapply(
      first$manifest_signatures, function(value) {
        .dsvert_dp_is_string(value) &&
          grepl("^[A-Za-z0-9_-]{86}$", value)
      }, logical(1L))))
    expect_identical(first_calls, second_calls)
    expect_false(first$operation_or_request_limit)
    expect_false(first$history_can_deny_operation)
    expect_identical(unname(vapply(
      first$manifest_build_receipts, `[[`, numeric(1L), "privacy_epoch")),
      as.numeric(seq_len(k)))
    expect_identical(
      first$artifact_commitments$context$privacy_epoch_scope,
      "per_peer_signed_receipts_v1")
    expect_false("privacy_epoch" %in%
                   names(first$artifact_commitments$context))

    contract <- .dsvert_joint_dp_client_decode(
      first$workload_contract_json, "test workload contract",
      8L * 1024L^2)
    expect_identical(
      contract$describe$owner_a_describe$owner_peer, "peer_a")
    expect_identical(
      contract$survival$owner_b_survival$owner_peer, "peer_b")
    expect_identical(
      contract$vertical_cross$owner_a_cross$owner_peer, "peer_a")

    draft_values <- lapply(fixture$drafts, function(value) {
      .dsvert_joint_dp_client_decode(
        value, "test signed draft", 8L * 1024L^2)
    })
    fields <- unlist(lapply(
      draft_values, .capsule_manifest_client_recursive_names),
      use.names = FALSE)
    expect_false(any(c(
      "snapshot_sha256", "alignment_manifest_hash", "row_count",
      "patient_hash", "patient_digest", "alignment_hash") %in% fields))

    for (calls in first_calls) {
      endpoint <- as.character(as.list(calls[[1L]])[[1L]])
      arguments <- names(as.list(calls[[1L]]))[-1L]
      if (identical(endpoint, "dsvertDPCapsuleManifestDraftDS")) {
        expect_length(arguments, 0L)
      } else {
        expect_setequal(
          arguments, c("schema_json", "workload_contract_json"))
      }
    }
  }
})

test_that("client accepts only the signed Gaussian fragment grammar", {
  fragments <- list(
    describe = list(), survival = list(),
    gaussian = list(gaussian_primary = list(
      version = "v1", dataset = "cohort", outcome = "y",
      predictors = list("x"), intercept = TRUE)),
    vertical_cross = list())
  normalized <- .dsvert_dp_capsule_manifest_fragments(fragments)
  expect_named(normalized$gaussian, "gaussian_primary")
  expect_identical(
    unlist(normalized$gaussian$gaussian_primary$predictors,
           use.names = FALSE),
    "x")
  invalid <- fragments
  invalid$gaussian$gaussian_primary$predictors <- list()
  expect_error(
    .dsvert_dp_capsule_manifest_fragments(invalid),
    "invalid custodian workload specification")

  random_intercept <- list(
    describe = list(), survival = list(), vertical_cross = list(),
    gaussian = list(lmm_primary = list(
      version = "random_intercept_v1", dataset = "cohort", outcome = "y",
      cluster = "site", max_patients_per_cluster = 20L)))
  normalized_random_intercept <- .dsvert_dp_capsule_manifest_fragments(
    random_intercept)
  expect_identical(
    normalized_random_intercept$gaussian$lmm_primary$cluster, "site")
  random_intercept$gaussian$lmm_primary$max_patients_per_cluster <- 1L
  expect_error(
    .dsvert_dp_capsule_manifest_fragments(random_intercept),
    "invalid custodian workload specification")

  fixed_random_intercept <- list(
    describe = list(), survival = list(), vertical_cross = list(),
    gaussian = list(lmm_fixed = list(
      version = "random_intercept_fixed_v2", dataset = "cohort",
      outcome = "y", predictors = "x", intercept = TRUE, cluster = "site",
      max_patients_per_cluster = 20L, variance_ratio_grid = c(0, 0.5, 2))))
  normalized_fixed <- .dsvert_dp_capsule_manifest_fragments(
    fixed_random_intercept)
  expect_identical(
    normalized_fixed$gaussian$lmm_fixed$variance_ratio_grid,
    c(0, 0.5, 2))
  fixed_random_intercept$gaussian$lmm_fixed$predictors <- list("x")
  fixed_random_intercept$gaussian$lmm_fixed$variance_ratio_grid <-
    as.list(c(0, 0.5, 2))
  normalized_fixed <- .dsvert_dp_capsule_manifest_fragments(
    fixed_random_intercept)
  expect_identical(
    normalized_fixed$gaussian$lmm_fixed$variance_ratio_grid,
    c(0, 0.5, 2))
  fixed_random_intercept$gaussian$lmm_fixed$variance_ratio_grid <- c(0.1, 1)
  expect_error(
    .dsvert_dp_capsule_manifest_fragments(fixed_random_intercept),
    "invalid custodian workload specification")
})

test_that("draft tampering, owner conflicts and divergent builds are rejected", {
  fixture <- .capsule_manifest_client_fixture(2L)
  changed <- .dsvert_joint_dp_client_decode(
    fixture$drafts$peer_a, "test draft", 8L * 1024L^2)
  changed$datasets$data_peer_a$columns$x_peer_a$upper <- 11
  fixture$state$drafts$peer_a <- .dsvert_joint_dp_client_json(changed)
  expect_error(
    .capsule_manifest_client_run(fixture),
    "invalid capsule manifest signature")

  duplicate <- .capsule_manifest_client_fixture(2L)
  peer_b <- .dsvert_joint_dp_client_decode(
    duplicate$drafts$peer_b, "test draft", 8L * 1024L^2)
  peer_b$workload_fragments$describe$owner_a_describe <- list(
    version = "v1", dataset = "data_peer_b", variables = "time_peer_b",
    histogram_grids = list(time_peer_b = c(5, 10)),
    allocation = list(
      count = 0.25, histogram = 0.25, sum = 0.25, sumsq = 0.25))
  peer_b$signature <- NULL
  peer_b <- c(peer_b, list(signature =
    .capsule_manifest_client_b64url(openssl::ed25519_sign(
      .dsvert_dp_capsule_manifest_message("draft", peer_b),
      duplicate$keys$peer_b))))
  duplicate$state$drafts$peer_b <- .dsvert_joint_dp_client_json(peer_b)
  expect_error(
    .capsule_manifest_client_run(duplicate),
    "same describe workload id")

  divergent <- .capsule_manifest_client_fixture(3L)
  divergent$state$different_build <- TRUE
  expect_error(
    .capsule_manifest_client_run(divergent),
    "built different biomedical capsule manifests")
})

test_that("manifest bootstrap producers retry exact expressions without quota", {
  endpoints <- list(
    call(name = "dsvertDPCapsuleManifestDraftDS"),
    call(name = "dsvertDPCapsuleManifestSignDS",
         schema_json = "fixed-schema",
         workload_contract_json = "fixed-workload"),
    call(name = "dsvertDPCapsuleManifestBuildDS",
         schema_json = "fixed-signed-schema",
         workload_contract_json = "fixed-workload"))
  for (request in endpoints) {
    observed <- list()
    aggregate <- function(conns, expr, error, async, errors.print) {
      observed[[length(observed) + 1L]] <<- expr
      if (length(observed) < 4L) stop("simulated lost response")
      list(site = list(ok = TRUE, request_limit = FALSE))
    }
    result <- testthat::with_mocked_bindings(
      .dsvert_aggregate_strict(
        list(site = structure(list(), class = "test")), request,
        operation = "manifest retry test", .aggregate = aggregate),
      .dsvert_retry_sleep = function(seconds) invisible(NULL),
      .dsvert_retry_jitter = function() 1)
    expect_length(observed, 4L)
    expect_true(all(vapply(
      observed, identical, logical(1L), observed[[1L]])))
    expect_false(result$site$request_limit)
  }
})
