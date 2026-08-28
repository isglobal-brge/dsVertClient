# Client relay for the server-authoritative biomedical capsule manifest.
# There are deliberately no analyst-provided schema, bound, domain, snapshot,
# or workload arguments in this protocol.

.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_DRAFT_VERSION <-
  "dsvert-biomedical-capsule-manifest-draft-v1"
.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_SIGN_VERSION <-
  "dsvert-biomedical-capsule-manifest-schema-signature-v1"
.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_BUILD_VERSION <-
  "dsvert-biomedical-capsule-manifest-build-v2"
.DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_VERSION <-
  "dsvert-biomedical-capsule-artifact-commitment-index-v3"
.DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_ENTRY_VERSION <-
  "dsvert-biomedical-capsule-artifact-commitment-v1"
.DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_MAX_BYTES <- 8L * 1024L^2
.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_REJECTION_VERSION <-
  "dsvert-biomedical-capsule-manifest-rejection-v1"
.DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION <-
  "dsvert-biomedical-capsule-schema-v1"
.DSVERT_CLIENT_DP_CAPSULE_POLICY_SCHEMA_VERSION <- "custodian-policy-v1"
.DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_CONTRACT_VERSION <-
  "dsvert-biomedical-capsule-workload-contract-v2"
.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_SIGNATURE_DOMAIN <-
  "dsVert/dp/biomedical-capsule-manifest/v1/"
.DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN <-
  "dsVert/dp/biomedical-capsule-schema/v1|"
.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_DRAFT_BYTES <- 8L * 1024L^2
.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_SCHEMA_BYTES <- 8L * 1024L^2
.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_WORKLOAD_BYTES <- 8L * 1024L^2
.DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_BUILD_BYTES <- 48L * 1024L^2

.dsvert_dp_capsule_manifest_rejection <- function(value, phase) {
  expected <- c(
    "version", "rejected", "phase", "reason_code", "retryable",
    "operation_limit", "request_limit", "history_can_deny_operation")
  if (is.list(value) && .dsvert_dp_has_exact_names(value, expected) &&
      identical(value$version,
                .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_REJECTION_VERSION) &&
      identical(value$rejected, TRUE) && identical(value$phase, phase) &&
      .dsvert_dp_is_string(value$reason_code) &&
      grepl("^[a-z][a-z0-9_]{0,63}$", value$reason_code) &&
      identical(value$retryable, FALSE) &&
      identical(value$operation_limit, FALSE) &&
      identical(value$request_limit, FALSE) &&
      identical(value$history_can_deny_operation, FALSE)) {
    stop("A custodian rejected capsule manifest phase '", phase,
         "' (", value$reason_code, ")", call. = FALSE)
  }
  invisible(FALSE)
}

.dsvert_dp_capsule_manifest_hash <- function(value) {
  digest::digest(
    .dsvert_joint_dp_client_json(value), algo = "sha256",
    serialize = FALSE)
}

.dsvert_dp_capsule_manifest_message <- function(domain, unsigned) {
  signed <- unsigned
  if (identical(domain, "build") && is.list(signed)) {
    signed$manifest_json <- NULL
    signed$artifact_commitments <- NULL
  }
  charToRaw(paste0(
    .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_SIGNATURE_DOMAIN, domain, "|",
    .dsvert_joint_dp_client_json(signed)))
}

.dsvert_dp_capsule_artifact_merkle_parent_client <- function(left, right) {
  .dsvert_dp_capsule_manifest_hash(list(
    protocol = "dsvert-biomedical-capsule-artifact-merkle-parent-v1",
    left = left, right = right))
}

.dsvert_dp_capsule_artifact_merkle_root_client <- function(
    hashes, context) {
  hashes <- unname(as.character(hashes))
  if (!length(hashes)) {
    return(.dsvert_dp_capsule_manifest_hash(list(
      protocol = "dsvert-biomedical-capsule-artifact-merkle-empty-v1",
      context = context)))
  }
  if (anyNA(hashes) || any(!grepl("^[0-9a-f]{64}$", hashes))) {
    stop("Invalid biomedical artifact commitment leaf", call. = FALSE)
  }
  nodes <- hashes
  while (length(nodes) > 1L) {
    if (length(nodes) %% 2L) nodes <- c(nodes, tail(nodes, 1L))
    nodes <- vapply(seq.int(1L, length(nodes), by = 2L), function(index) {
      .dsvert_dp_capsule_artifact_merkle_parent_client(
        nodes[[index]], nodes[[index + 1L]])
    }, character(1L))
  }
  nodes[[1L]]
}

.dsvert_dp_capsule_artifact_merkle_proof_client <- function(
    hashes, index) {
  hashes <- unname(as.character(hashes))
  if (!length(hashes) || !.dsvert_dp_is_integer(index, 1, length(hashes))) {
    stop("Invalid biomedical artifact proof index", call. = FALSE)
  }
  position <- as.integer(index)
  nodes <- hashes
  proof <- list()
  while (length(nodes) > 1L) {
    original_length <- length(nodes)
    if (original_length %% 2L) nodes <- c(nodes, tail(nodes, 1L))
    sibling <- if (position %% 2L) position + 1L else position - 1L
    proof[[length(proof) + 1L]] <- list(
      side = if (position %% 2L) "right" else "left",
      sha256 = nodes[[sibling]])
    nodes <- vapply(seq.int(1L, length(nodes), by = 2L), function(node) {
      .dsvert_dp_capsule_artifact_merkle_parent_client(
        nodes[[node]], nodes[[node + 1L]])
    }, character(1L))
    position <- as.integer((position + 1L) %/% 2L)
  }
  proof
}

.dsvert_dp_capsule_artifact_commitment_index_client <- function(
    manifest, policy, manifest_sha256) {
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  identity <- manifest$capsule_identity
  contract <- identity$contract
  context <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_VERSION,
    manifest_sha256 = manifest_sha256,
    capsule_id = identity$capsule_id,
    domain = policy$domain,
    cohort_id = policy$cohort_id,
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    designated_noise_peers = as.list(
      sort(policy$designated_noise_peers, method = "radix")),
    privacy_epoch_scope = if (identical(
      manifest$capsule_identity$contract$protocol,
      "dsvert-stateless-catalog-synopsis-capsule-identity-v1")) {
      "per_canonical_artifact_sticky_v1"
    } else {
      "per_peer_signed_receipts_v1"
    },
    epsilon = as.numeric(policy$capsule_epsilon),
    delta = as.numeric(policy$capsule_delta),
    adjacency = policy$adjacency,
    unit_capacity = as.numeric(policy$unit_capacity),
    max_records_per_unit = as.numeric(policy$max_records_per_unit),
    overflow_policy = policy$overflow_policy,
    consortium_id = contract$consortium_id,
    policy_contract_hash = contract$policy_contract_hash,
    logical_snapshot = manifest$logical_snapshot,
    capsule_schema = manifest$capsule_schema,
    admission_sha256 = .dsvert_dp_capsule_manifest_hash(manifest$admission),
    bounds_sha256 = .dsvert_dp_capsule_manifest_hash(manifest$bounds),
    workload_sha256 = .dsvert_dp_capsule_manifest_hash(manifest$workload),
    release_lattice = manifest$workload$release_lattice,
    capsule_mechanism = manifest$workload$capsule_mechanism,
    mechanism_selection = manifest$workload$mechanism_selection,
    coordinate_layout_version = layout$version,
    coordinate_count = as.numeric(layout$coordinate_count),
    coordinate_order_sha256 = layout$sha256))
  gaussian <- layout$blocks[vapply(layout$blocks, function(block) {
    identical(block$family, "gaussian_models")
  }, logical(1L))]
  if (length(gaussian)) {
    gaussian <- gaussian[order(vapply(
      gaussian, `[[`, character(1L), "key"), method = "radix")]
  }
  entries <- lapply(gaussian, function(block) {
    entry <- .dsvert_joint_dp_client_canonical(list(
      version = .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_ENTRY_VERSION,
      family = block$family, analysis_id = block$key,
      dataset = block$dataset, owner_peer = block$owner_peer,
      start = as.numeric(block$start), end = as.numeric(block$end),
      length = as.numeric(block$length),
      descriptor_sha256 = .dsvert_dp_capsule_manifest_hash(
        block$descriptor)))
    c(entry, list(node_sha256 = .dsvert_dp_capsule_manifest_hash(list(
      protocol = "dsvert-biomedical-capsule-artifact-merkle-leaf-v1",
      context = context, entry = entry))))
  })
  names(entries) <- vapply(gaussian, `[[`, character(1L), "key")
  value <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_VERSION,
    context = context, gaussian_models = entries))
  if (nchar(.dsvert_joint_dp_client_json(value), type = "bytes") >
      .DSVERT_CLIENT_DP_CAPSULE_ARTIFACT_INDEX_MAX_BYTES) {
    stop("The public biomedical artifact commitment index is too large",
         call. = FALSE)
  }
  hashes <- vapply(entries, `[[`, character(1L), "node_sha256")
  list(
    value = value, count = as.numeric(length(entries)),
    root = .dsvert_dp_capsule_artifact_merkle_root_client(hashes, context))
}

.dsvert_dp_capsule_manifest_verify <- function(
    value, domain, peer, context) {
  if (!is.list(value) || !"signature" %in% names(value) ||
      !peer %in% names(context$pinset)) {
    stop("Invalid signed capsule manifest object", call. = FALSE)
  }
  unsigned <- value[setdiff(names(value), "signature")]
  message <- .dsvert_dp_capsule_manifest_message(domain, unsigned)
  public <- .dsvert_joint_dp_client_b64url(
    unname(context$pinset[[peer]]), 32L, "capsule manifest identity key")
  signature <- .dsvert_joint_dp_client_b64url(
    value$signature, 64L, "capsule manifest signature")
  valid <- tryCatch(openssl::ed25519_verify(
    message, signature, openssl::read_ed25519_pubkey(public)),
    error = function(error) FALSE)
  if (!isTRUE(valid)) {
    stop("A pinned peer returned an invalid capsule manifest signature",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_capsule_schema_verify <- function(
    unsigned, signature, peer, context) {
  message <- charToRaw(paste0(
    .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
    .dsvert_joint_dp_client_json(unsigned)))
  public <- .dsvert_joint_dp_client_b64url(
    unname(context$pinset[[peer]]), 32L, "capsule schema identity key")
  signature <- .dsvert_joint_dp_client_b64url(
    signature, 64L, "capsule schema signature")
  valid <- tryCatch(openssl::ed25519_verify(
    message, signature, openssl::read_ed25519_pubkey(public)),
    error = function(error) FALSE)
  if (!isTRUE(valid)) {
    stop("A pinned peer returned an invalid global capsule schema signature",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_capsule_manifest_string_array <- function(value, what) {
  valid <- (is.character(value) && length(value) > 0L &&
              is.null(names(value)) && !anyNA(value)) ||
    (is.list(value) && is.null(names(value)) && length(value) > 0L &&
       all(vapply(value, .dsvert_dp_is_string, logical(1L))))
  result <- if (isTRUE(valid)) {
    enc2utf8(unname(unlist(value, use.names = FALSE)))
  } else {
    character()
  }
  if (!isTRUE(valid) || anyDuplicated(result) ||
      any(!nzchar(trimws(result)))) {
    stop("Invalid capsule manifest ", what, " returned by a peer",
         call. = FALSE)
  }
  sort(result, method = "radix")
}

.dsvert_dp_capsule_manifest_column <- function(value, peer) {
  if (!is.list(value) || !.dsvert_dp_has_exact_names(
        value, if (identical(value$kind, "numeric")) {
          c("kind", "owner_peer", "lower", "upper")
        } else {
          c("kind", "owner_peer", "levels")
        }) ||
      !identical(value$owner_peer, peer) ||
      !value$kind %in% c("numeric", "categorical")) {
    stop("Invalid custodian capsule column metadata", call. = FALSE)
  }
  if (identical(value$kind, "numeric")) {
    if (!.dsvert_dp_is_number(value$lower) ||
        !.dsvert_dp_is_number(value$upper) ||
        value$lower >= value$upper ||
        !is.finite(value$upper - value$lower)) {
      stop("Invalid custodian numeric bounds", call. = FALSE)
    }
    return(list(
      kind = "numeric", owner_peer = peer,
      lower = as.numeric(value$lower), upper = as.numeric(value$upper)))
  }
  list(
    kind = "categorical", owner_peer = peer,
    levels = .dsvert_dp_capsule_manifest_string_array(
      value$levels, "categorical domain"))
}

.dsvert_dp_capsule_manifest_number_array <- function(value, what) {
  result <- if (is.numeric(value) && length(value) > 0L &&
               is.null(names(value))) {
    as.numeric(value)
  } else if (is.list(value) && is.null(names(value)) && length(value) > 0L &&
             all(vapply(value, .dsvert_dp_is_number, logical(1L)))) {
    as.numeric(unlist(value, use.names = FALSE))
  } else {
    numeric()
  }
  if (!length(result) || anyNA(result) || any(!is.finite(result))) {
    stop("Invalid capsule manifest ", what, " returned by a peer",
         call. = FALSE)
  }
  result
}

.dsvert_dp_capsule_manifest_fragments <- function(value) {
  families <- c("describe", "survival", "gaussian", "vertical_cross")
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      anyDuplicated(names(value)) || !setequal(names(value), families)) {
    stop("A peer returned invalid custodian workload fragments",
         call. = FALSE)
  }
  result <- stats::setNames(vector("list", length(families)), families)
  identifier <- function(value) {
    .dsvert_dp_is_string(value) &&
      grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)
  }
  column_reference <- function(value) {
    .dsvert_dp_is_string(value) && grepl(
      paste0(
        "^(?:[A-Za-z0-9][A-Za-z0-9._:-]{0,127}\\$)?",
        "[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$"),
      value)
  }
  for (family in families) {
    entries <- value[[family]]
    if (is.null(entries)) entries <- list()
    if (!is.list(entries) || (length(entries) &&
        (is.null(names(entries)) || anyNA(names(entries)) ||
         any(!grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$",
                    names(entries))) || anyDuplicated(names(entries))))) {
      stop("A peer returned invalid custodian workload fragments",
           call. = FALSE)
    }
    if (length(entries)) {
      entries <- entries[order(names(entries), method = "radix")]
    }
    normalized <- vector("list", length(entries))
    names(normalized) <- names(entries)
    for (analysis_id in names(entries)) {
      spec <- entries[[analysis_id]]
      if (!is.list(spec) || is.null(names(spec)) || anyNA(names(spec)) ||
          anyDuplicated(names(spec))) {
        stop("A peer returned an invalid custodian workload specification",
             call. = FALSE)
      }
      if (identical(family, "describe")) {
        expected <- c(
          "version", "dataset", "variables", "histogram_grids",
          "allocation")
        valid <- setequal(names(spec), expected) &&
          identifier(spec$version) && identifier(spec$dataset)
        variables <- if (isTRUE(valid)) {
          .dsvert_dp_capsule_manifest_string_array(
            spec$variables, "describe variables")
        } else character()
        grids <- spec$histogram_grids
        allocation <- spec$allocation
        valid <- isTRUE(valid) && is.list(grids) &&
          !is.null(names(grids)) && identical(names(grids), variables) &&
          is.list(allocation) && !is.null(names(allocation)) &&
          setequal(names(allocation),
                   c("count", "sum", "sumsq", "histogram")) &&
          all(vapply(allocation, .dsvert_dp_is_number, logical(1L)))
        if (isTRUE(valid)) {
          grid_values <- lapply(grids, function(grid) {
            .dsvert_dp_capsule_manifest_number_array(
              grid, "describe histogram grid")
          })
          weights <- as.numeric(unlist(allocation[
            c("count", "sum", "sumsq", "histogram")],
            use.names = FALSE))
          valid <- all(vapply(grid_values, function(grid) {
            all(diff(grid) > 0)
          }, logical(1L))) && all(weights > 0) &&
            abs(sum(weights) - 1) <= 1024 * .Machine$double.eps
        }
      } else if (identical(family, "survival")) {
        if (identical(spec$version, "cox_partial_likelihood_grid_v1")) {
          expected <- c(
            "version", "dataset", "time", "event", "censor",
            "event_level", "time_grid", "predictors", "intercept",
            "candidate_grid")
          valid <- setequal(names(spec), expected) &&
            all(vapply(spec[c("version", "dataset", "time", "event")],
                       identifier, logical(1L))) &&
            .dsvert_dp_is_string(spec$censor) && nzchar(spec$censor) &&
            .dsvert_dp_is_string(spec$event_level) && nzchar(spec$event_level) &&
            !identical(spec$censor, spec$event_level) &&
            identical(spec$intercept, FALSE)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "Cox fixed predictors"),
            error = function(error) character()) else character()
          grid <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_number_array(
              spec$time_grid, "Cox time grid"),
            error = function(error) numeric()) else numeric()
          candidates <- if (isTRUE(valid) && is.list(spec$candidate_grid)) {
            lapply(spec$candidate_grid, function(beta) tryCatch(
              .dsvert_dp_capsule_manifest_number_array(
                beta, "Cox beta grid row"),
              error = function(error) numeric()))
          } else list()
          valid <- isTRUE(valid) && length(predictors) && length(grid) &&
            all(diff(grid) > 0) && length(candidates) &&
            all(vapply(candidates, function(beta) {
              length(beta) == length(predictors) && all(is.finite(beta)) &&
                all(abs(beta) <= 4) && sum(abs(beta)) <= 8
            }, logical(1L)))
        } else {
          expected <- c(
            "version", "dataset", "time", "event", "censor",
            "time_grid", "entry")
          valid <- setequal(names(spec), expected) &&
            all(vapply(spec[c(
              "version", "dataset", "time", "event")],
              identifier, logical(1L))) &&
            .dsvert_dp_is_string(spec$censor) && nzchar(spec$censor) &&
            (is.null(spec$entry) || identifier(spec$entry))
          if (isTRUE(valid)) {
            grid <- .dsvert_dp_capsule_manifest_number_array(
              spec$time_grid, "survival time grid")
            valid <- all(diff(grid) > 0)
          }
        }
      } else if (identical(family, "gaussian")) {
        if (identical(spec$version, "ordinal_grid_v1")) {
          expected <- c(
            "version", "dataset", "outcome", "predictors", "intercept",
            "ordered_levels", "candidate_grid")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            isTRUE(spec$intercept)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "ordinal fixed predictors"),
            error = function(error) character()) else character()
          ordered_levels <- if (isTRUE(valid) &&
              ((is.character(spec$ordered_levels) &&
                is.null(names(spec$ordered_levels))) ||
               (is.list(spec$ordered_levels) &&
                is.null(names(spec$ordered_levels)) &&
                all(vapply(spec$ordered_levels, .dsvert_dp_is_string,
                           logical(1L)))))) {
            enc2utf8(unname(unlist(spec$ordered_levels, use.names = FALSE)))
          } else character()
          candidate_grid <- spec$candidate_grid
          if (!is.list(candidate_grid) || !is.null(names(candidate_grid))) {
            candidate_grid <- list()
          } else {
            candidate_grid <- lapply(candidate_grid, function(candidate) {
              if (!is.list(candidate) || is.null(names(candidate)) ||
                  !setequal(names(candidate), c("thresholds", "beta"))) {
                return(NULL)
              }
              thresholds <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
                candidate$thresholds, "ordinal threshold row"),
                error = function(error) numeric())
              beta <- tryCatch(.dsvert_dp_capsule_manifest_number_array(
                candidate$beta, "ordinal beta row"),
                error = function(error) numeric())
              list(thresholds = thresholds, beta = beta)
            })
          }
          candidate_keys <- if (length(candidate_grid) && all(vapply(
                candidate_grid, is.list, logical(1L)))) vapply(
            candidate_grid, .dsvert_joint_dp_client_json, character(1L)) else character()
          expected_dimension <- 1L + length(predictors)
          valid <- isTRUE(valid) && length(predictors) &&
            !anyDuplicated(predictors) && !spec$outcome %in% predictors &&
            all(vapply(predictors, column_reference, logical(1L))) &&
            length(ordered_levels) >= 3L && !anyNA(ordered_levels) &&
            anyDuplicated(ordered_levels) == 0L &&
            length(candidate_grid) && length(candidate_grid) <= 256L &&
            all(vapply(candidate_grid, function(candidate) {
              length(candidate$thresholds) == length(ordered_levels) - 1L &&
                length(candidate$beta) == expected_dimension &&
                !anyNA(candidate$thresholds) && !anyNA(candidate$beta) &&
                all(is.finite(candidate$thresholds)) &&
                all(is.finite(candidate$beta)) &&
                all(abs(candidate$thresholds) <= 8) &&
                all(abs(candidate$beta) <= 8) &&
                all(diff(candidate$thresholds) >= 1 / 256)
            }, logical(1L))) && !anyDuplicated(candidate_keys) &&
            identical(candidate_keys, sort(candidate_keys, method = "radix"))
        } else if (identical(spec$version, "multinomial_grid_v1")) {
          expected <- c(
            "version", "dataset", "outcome", "predictors", "intercept",
            "levels", "reference", "beta_grid")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            isTRUE(spec$intercept)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "multinomial fixed predictors"),
            error = function(error) character()) else character()
          levels <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$levels, "multinomial outcome levels"),
            error = function(error) character()) else character()
          reference <- if (isTRUE(valid) && .dsvert_dp_is_string(spec$reference)) {
            enc2utf8(spec$reference)
          } else NA_character_
          beta_grid <- spec$beta_grid
          if (!is.list(beta_grid) || !is.null(names(beta_grid))) {
            beta_grid <- list()
          } else {
            beta_grid <- lapply(beta_grid, function(beta) tryCatch(
              .dsvert_dp_capsule_manifest_number_array(
                beta, "multinomial beta grid row"),
              error = function(error) numeric()))
          }
          expected_dimension <- (length(levels) - 1L) * (1L + length(predictors))
          beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
            .dsvert_joint_dp_client_json(as.list(beta))
          }, character(1L)) else character()
          valid <- isTRUE(valid) && length(predictors) &&
            !anyDuplicated(predictors) && !spec$outcome %in% predictors &&
            all(vapply(predictors, column_reference, logical(1L))) &&
            length(levels) >= 3L && !anyDuplicated(levels) &&
            identical(levels, sort(levels, method = "radix")) &&
            !is.na(reference) && reference %in% levels &&
            length(beta_grid) && length(beta_grid) <= 256L &&
            all(vapply(beta_grid, function(beta) {
              length(beta) == expected_dimension && !anyNA(beta) &&
                all(is.finite(beta)) && all(abs(beta) <= 8)
            }, logical(1L))) && !anyDuplicated(beta_keys)
          if (isTRUE(valid)) beta_grid <- beta_grid[order(beta_keys)]
        } else if (spec$version %in% c("binomial_grid_v1",
                                       "poisson_grid_v1")) {
          poisson <- identical(spec$version, "poisson_grid_v1")
          expected <- c(
            "version", "dataset", "outcome", "predictors", "intercept",
            "beta_grid")
          if (isTRUE(poisson)) expected <- c(expected, "max_outcome")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            isTRUE(spec$intercept)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "finite GLM fixed predictors"),
            error = function(error) character()) else character()
          beta_grid <- spec$beta_grid
          if (!is.list(beta_grid) || !is.null(names(beta_grid))) {
            beta_grid <- list()
          } else {
            beta_grid <- lapply(beta_grid, function(beta) tryCatch(
              .dsvert_dp_capsule_manifest_number_array(
                beta, "finite GLM beta grid row"),
              error = function(error) numeric()))
          }
          expected_dimension <- 1L + length(predictors)
          beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
            .dsvert_joint_dp_client_json(as.list(beta))
          }, character(1L)) else character()
          maximum_valid <- if (isTRUE(poisson)) {
            .dsvert_dp_is_integer(spec$max_outcome, 1L, 1024L)
          } else {
            TRUE
          }
          valid <- isTRUE(valid) && isTRUE(maximum_valid) &&
            length(predictors) && !anyDuplicated(predictors) &&
            !spec$outcome %in% predictors &&
            all(vapply(predictors, column_reference, logical(1L))) &&
            length(beta_grid) && length(beta_grid) <= 256L &&
            all(vapply(beta_grid, function(beta) {
              length(beta) == expected_dimension && !anyNA(beta) &&
                all(is.finite(beta)) && all(abs(beta) <= 8) &&
                sum(abs(beta)) <= 16
            }, logical(1L))) && !anyDuplicated(beta_keys)
          if (isTRUE(valid)) beta_grid <- beta_grid[order(beta_keys)]
        } else if (identical(spec$version, "negative_binomial_grid_v1")) {
          expected <- c(
            "version", "dataset", "outcome", "predictors", "intercept",
            "max_outcome", "beta_grid", "theta_grid")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            isTRUE(spec$intercept) && .dsvert_dp_is_integer(
              spec$max_outcome, 1L, 1024L)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "NB2 fixed predictors"),
            error = function(error) character()) else character()
          theta_grid <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_number_array(
              spec$theta_grid, "NB2 theta grid"),
            error = function(error) numeric()) else numeric()
          beta_grid <- spec$beta_grid
          if (!is.list(beta_grid) || !is.null(names(beta_grid))) {
            beta_grid <- list()
          } else {
            beta_grid <- lapply(beta_grid, function(beta) tryCatch(
              .dsvert_dp_capsule_manifest_number_array(
                beta, "NB2 beta grid row"), error = function(error) numeric()))
          }
          expected_dimension <- 1L + length(predictors)
          beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
            .dsvert_joint_dp_client_json(as.list(beta))
          }, character(1L)) else character()
          valid <- isTRUE(valid) && length(predictors) &&
            !anyDuplicated(predictors) && !spec$outcome %in% predictors &&
            all(vapply(predictors, column_reference, logical(1L))) &&
            length(theta_grid) && !anyNA(theta_grid) &&
            all(is.finite(theta_grid)) && all(theta_grid > 0) &&
            all(theta_grid <= 64) && all(diff(theta_grid) > 0) &&
            length(beta_grid) && length(beta_grid) * length(theta_grid) <= 256L &&
            all(vapply(beta_grid, function(beta) {
              length(beta) == expected_dimension && !anyNA(beta) &&
                all(is.finite(beta)) && all(abs(beta) <= 8)
            }, logical(1L))) && !anyDuplicated(beta_keys)
          if (isTRUE(valid)) beta_grid <- beta_grid[order(beta_keys)]
        } else if (spec$version %in% c("gaussian_ar1_working_gls_grid_v1",
                                       "gaussian_ar1_robust_working_gls_grid_v1")) {
          robust_sandwich <- identical(
            spec$version, "gaussian_ar1_robust_working_gls_grid_v1")
          expected <- c(
            "version", "dataset", "outcome", "cluster", "order",
            "predictors", "intercept", "max_patients_per_cluster",
            "candidate_grid")
          if (isTRUE(robust_sandwich)) expected <- c(expected, "score_clip")
          valid <- setequal(names(spec), expected) && identifier(spec$dataset) &&
            column_reference(spec$outcome) && column_reference(spec$cluster) &&
            column_reference(spec$order) && isTRUE(spec$intercept) &&
            .dsvert_dp_is_integer(spec$max_patients_per_cluster, 2L) &&
            (!isTRUE(robust_sandwich) ||
             .dsvert_dp_is_integer(spec$max_patients_per_cluster, 2L, 32L))
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "GEE AR1 fixed predictors"),
            error = function(error) character()) else character()
          valid <- isTRUE(valid) && length(predictors) &&
            !anyDuplicated(predictors) && all(vapply(
              predictors, column_reference, logical(1L))) &&
            (!isTRUE(robust_sandwich) || length(predictors) <= 3L) &&
            length(unique(c(spec$outcome, spec$cluster, spec$order,
                            predictors))) == 3L + length(predictors)
          candidates <- if (isTRUE(valid)) .dsvert_dp_gee_ar1_grid_candidates(
            spec$candidate_grid, as.numeric(spec$max_patients_per_cluster)) else list()
          score_clip <- suppressWarnings(as.numeric(spec$score_clip))
          valid <- isTRUE(valid) && length(candidates) && length(candidates) <=
            (if (isTRUE(robust_sandwich)) 64L else 128L) &&
            (!isTRUE(robust_sandwich) ||
             (length(score_clip) == 1L && is.finite(score_clip) &&
              score_clip >= 0.25 && score_clip <= 32)) &&
            all(vapply(candidates, function(candidate) {
              length(candidate$beta) == 1L + length(predictors)
            }, logical(1L)))
          if (isTRUE(valid)) {
            candidate_grid <- lapply(candidates, function(candidate) list(
              beta = unname(candidate$beta), rho = candidate$rho))
            if (isTRUE(robust_sandwich)) spec$score_clip <- score_clip
          }
        } else if (identical(spec$version, "gaussian_random_slope_grid_v1")) {
          expected <- c(
            "version", "dataset", "outcome", "cluster", "predictors",
            "random_slopes", "intercept", "max_patients_per_cluster",
            "candidate_grid")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            column_reference(spec$cluster) && !identical(
              spec$outcome, spec$cluster) && isTRUE(spec$intercept) &&
            .dsvert_dp_is_integer(spec$max_patients_per_cluster, 2L)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "LMM fixed predictors"),
            error = function(error) character()) else character()
          random_slopes <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$random_slopes, "LMM random slopes"),
            error = function(error) character()) else character()
          valid <- isTRUE(valid) && length(predictors) &&
            !anyDuplicated(predictors) && !spec$outcome %in% predictors &&
            !spec$cluster %in% predictors && all(vapply(
              predictors, column_reference, logical(1L))) &&
            length(random_slopes) && !anyDuplicated(random_slopes) &&
            all(random_slopes %in% predictors)
          effects <- c("(Intercept)", sort(random_slopes, method = "radix"))
          candidates <- if (isTRUE(valid)) .dsvert_dp_lmm_random_slope_candidates(
            spec$candidate_grid, 1L + length(predictors), effects,
            as.numeric(spec$max_patients_per_cluster)) else list()
          valid <- isTRUE(valid) && length(candidates) && length(candidates) <= 128L
          if (isTRUE(valid)) {
            keys <- vapply(candidates, function(candidate) .dsvert_joint_dp_client_json(
              list(beta = as.list(candidate$beta), sigma2 = candidate$sigma2,
                   covariance = as.list(as.vector(t(candidate$covariance))))),
              character(1L))
            candidate_grid <- lapply(candidates, function(candidate) list(
              beta = unname(candidate$beta), sigma2 = candidate$sigma2,
              covariance = unname(as.vector(t(candidate$covariance)))))
          }
        } else if (identical(spec$version, "binary_random_intercept_grid_v1")) {
          expected <- c(
            "version", "dataset", "outcome", "cluster", "predictors",
            "intercept", "max_patients_per_cluster", "beta_grid",
            "variance_grid")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            column_reference(spec$cluster) && !identical(
              spec$outcome, spec$cluster) && isTRUE(spec$intercept) &&
            .dsvert_dp_is_integer(spec$max_patients_per_cluster, 2L)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "GLMM fixed predictors"),
            error = function(error) character()) else character()
          variance_grid <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_number_array(
              spec$variance_grid, "GLMM variance grid"),
            error = function(error) numeric()) else numeric()
          beta_grid <- spec$beta_grid
          if (!is.list(beta_grid) || !is.null(names(beta_grid))) {
            beta_grid <- list()
          } else {
            beta_grid <- lapply(beta_grid, function(beta) {
              tryCatch(.dsvert_dp_capsule_manifest_number_array(
                beta, "GLMM beta grid row"), error = function(error) numeric())
            })
          }
          expected_dimension <- 1L + length(predictors)
          beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
            .dsvert_joint_dp_client_json(as.list(beta))
          }, character(1L)) else character()
          valid <- isTRUE(valid) && length(predictors) &&
            !anyDuplicated(predictors) && !spec$outcome %in% predictors &&
            !spec$cluster %in% predictors && all(vapply(
              predictors, column_reference, logical(1L))) &&
            length(variance_grid) && !anyNA(variance_grid) &&
            all(is.finite(variance_grid)) && all(variance_grid >= 0) &&
            all(variance_grid <= 16) && all(diff(variance_grid) > 0) &&
            length(beta_grid) && length(beta_grid) * length(variance_grid) <= 256L &&
            all(vapply(beta_grid, function(beta) {
              length(beta) == expected_dimension && !anyNA(beta) &&
                all(is.finite(beta)) && all(abs(beta) <= 8)
            }, logical(1L))) && !anyDuplicated(beta_keys)
          if (isTRUE(valid)) beta_grid <- beta_grid[order(beta_keys)]
        } else if (identical(
              spec$version, "poisson_random_intercept_grid_v1")) {
          expected <- c(
            "version", "dataset", "outcome", "cluster", "predictors",
            "intercept", "max_patients_per_cluster", "max_outcome",
            "beta_grid", "variance_grid")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            column_reference(spec$cluster) && !identical(
              spec$outcome, spec$cluster) && isTRUE(spec$intercept) &&
            .dsvert_dp_is_integer(spec$max_patients_per_cluster, 2L) &&
            .dsvert_dp_is_integer(spec$max_outcome, 1L, 1024L)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "Poisson GLMM fixed predictors"),
            error = function(error) character()) else character()
          variance_grid <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_number_array(
              spec$variance_grid, "Poisson GLMM variance grid"),
            error = function(error) numeric()) else numeric()
          beta_grid <- spec$beta_grid
          if (!is.list(beta_grid) || !is.null(names(beta_grid))) {
            beta_grid <- list()
          } else {
            beta_grid <- lapply(beta_grid, function(beta) {
              tryCatch(.dsvert_dp_capsule_manifest_number_array(
                beta, "Poisson GLMM beta grid row"),
              error = function(error) numeric())
            })
          }
          expected_dimension <- 1L + length(predictors)
          beta_keys <- if (length(beta_grid)) vapply(beta_grid, function(beta) {
            .dsvert_joint_dp_client_json(as.list(beta))
          }, character(1L)) else character()
          valid <- isTRUE(valid) && length(predictors) &&
            !anyDuplicated(predictors) && !spec$outcome %in% predictors &&
            !spec$cluster %in% predictors && all(vapply(
              predictors, column_reference, logical(1L))) &&
            identical(predictors, sort(predictors, method = "radix")) &&
            length(variance_grid) && !anyNA(variance_grid) &&
            all(is.finite(variance_grid)) && all(variance_grid >= 0) &&
            all(variance_grid <= 16) && all(diff(variance_grid) > 0) &&
            length(beta_grid) && length(beta_grid) * length(variance_grid) <= 256L &&
            all(vapply(beta_grid, function(beta) {
              length(beta) == expected_dimension && !anyNA(beta) &&
                all(is.finite(beta)) && all(abs(beta) <= 8) &&
                sum(abs(beta)) <= 16
            }, logical(1L))) && !anyDuplicated(beta_keys)
          if (isTRUE(valid)) beta_grid <- beta_grid[order(beta_keys)]
        } else if (identical(spec$version, "binary_random_slope_grid_v1")) {
          expected <- c(
            "version", "dataset", "outcome", "cluster", "predictors",
            "random_slopes", "intercept", "max_patients_per_cluster",
            "candidate_grid")
          valid <- setequal(names(spec), expected) && identifier(spec$dataset) &&
            column_reference(spec$outcome) && column_reference(spec$cluster) &&
            !identical(spec$outcome, spec$cluster) && isTRUE(spec$intercept) &&
            .dsvert_dp_is_integer(spec$max_patients_per_cluster, 2L)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "GLMM fixed predictors"),
            error = function(error) character()) else character()
          random_slopes <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$random_slopes, "GLMM random slopes"),
            error = function(error) character()) else character()
          valid <- isTRUE(valid) && length(predictors) && !anyDuplicated(predictors) &&
            !spec$outcome %in% predictors && !spec$cluster %in% predictors &&
            all(vapply(predictors, column_reference, logical(1L))) &&
            length(random_slopes) >= 1L && length(random_slopes) <= 3L &&
            all(random_slopes %in% predictors) &&
            identical(predictors, sort(predictors, method = "radix")) &&
            identical(random_slopes, sort(random_slopes, method = "radix"))
          effects <- c("(Intercept)", sort(random_slopes, method = "radix"))
          candidates <- if (isTRUE(valid)) .dsvert_dp_glmm_random_slope_candidates(
            spec$candidate_grid, 1L + length(predictors), effects) else list()
          valid <- isTRUE(valid) && length(candidates) && length(candidates) <= 128L
          if (isTRUE(valid)) {
            candidate_grid <- lapply(candidates, function(candidate) list(
              beta = unname(candidate$beta),
              covariance = unname(as.vector(t(candidate$covariance)))))
          }
        } else if (spec$version %in% c("random_intercept_fixed_v2",
                                "random_intercept_fixed_v3")) {
          expected <- c(
            "version", "dataset", "outcome", "cluster", "predictors",
            "intercept", "max_patients_per_cluster", "variance_ratio_grid")
          if (identical(spec$version, "random_intercept_fixed_v3")) {
            expected <- c(expected, "estimation_profile")
          }
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            column_reference(spec$cluster) && !identical(
              spec$outcome, spec$cluster) && isTRUE(spec$intercept) &&
            .dsvert_dp_is_integer(spec$max_patients_per_cluster, 2L)
          predictors <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "LMM fixed predictors"),
            error = function(error) character()) else character()
          grid <- if (isTRUE(valid)) tryCatch(
            .dsvert_dp_capsule_manifest_number_array(
              spec$variance_ratio_grid, "LMM fixed variance-ratio grid"),
            error = function(error) numeric()) else numeric()
          valid <- isTRUE(valid) && length(predictors) &&
            !anyDuplicated(predictors) && !spec$outcome %in% predictors &&
            !spec$cluster %in% predictors && all(vapply(
              predictors, column_reference, logical(1L))) &&
            length(grid) &&
            !anyNA(grid) && all(is.finite(grid)) && all(grid >= 0) &&
            identical(grid[[1L]], 0) && all(diff(grid) > 0)
          if (identical(spec$version, "random_intercept_fixed_v3")) {
            valid <- isTRUE(valid) && .dsvert_dp_is_string(
              spec$estimation_profile) && spec$estimation_profile %in%
              c("ml", "reml")
          }
        } else if (identical(spec$version, "random_intercept_v1")) {
          expected <- c(
            "version", "dataset", "outcome", "cluster",
            "max_patients_per_cluster")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$dataset) && column_reference(spec$outcome) &&
            column_reference(spec$cluster) &&
            !identical(spec$outcome, spec$cluster) &&
            .dsvert_dp_is_integer(spec$max_patients_per_cluster, 2L)
        } else {
          expected <- c(
            "version", "dataset", "outcome", "predictors", "intercept")
          valid <- setequal(names(spec), expected) &&
            identifier(spec$version) && identifier(spec$dataset) &&
            column_reference(spec$outcome) &&
            is.logical(spec$intercept) && length(spec$intercept) == 1L &&
            !is.na(spec$intercept)
          predictors <- if (isTRUE(valid)) {
            tryCatch(.dsvert_dp_capsule_manifest_string_array(
              spec$predictors, "Gaussian predictors"),
            error = function(error) character())
          } else character()
          valid <- isTRUE(valid) && length(predictors) > 0L &&
            all(vapply(predictors, column_reference, logical(1L))) &&
            !spec$outcome %in% predictors
        }
      } else {
        expected <- c(
          "version", "left_dataset", "right_dataset", "left", "right",
          "family")
        valid <- setequal(names(spec), expected) &&
          all(vapply(spec[c("version", "left_dataset", "right_dataset")],
                     identifier, logical(1L))) &&
          column_reference(spec$left) && column_reference(spec$right) &&
          .dsvert_dp_is_string(spec$family) && spec$family %in% c(
            "categorical_pair", "numeric_cross_moment",
            "numeric_by_category") && !identical(spec$left, spec$right)
      }
      if (!isTRUE(valid)) {
        stop("A peer returned an invalid custodian workload specification",
             call. = FALSE)
      }
      if (identical(family, "gaussian") && spec$version %in% c(
            "gaussian_ar1_working_gls_grid_v1",
            "gaussian_ar1_robust_working_gls_grid_v1")) {
        spec$predictors <- predictors
        spec$candidate_grid <- candidate_grid
      } else if (identical(family, "gaussian") &&
                 identical(spec$version, "gaussian_random_slope_grid_v1")) {
        spec$predictors <- predictors
        spec$random_slopes <- sort(random_slopes, method = "radix")
        spec$candidate_grid <- candidate_grid
      } else if (identical(family, "gaussian") &&
                 identical(spec$version, "binary_random_slope_grid_v1")) {
        spec$predictors <- predictors
        spec$random_slopes <- sort(random_slopes, method = "radix")
        spec$candidate_grid <- candidate_grid
      } else if (identical(family, "gaussian") && spec$version %in% c(
            "random_intercept_fixed_v2", "random_intercept_fixed_v3")) {
        spec$predictors <- predictors
        spec$variance_ratio_grid <- grid
      } else if (identical(family, "gaussian") &&
                 identical(spec$version, "binary_random_intercept_grid_v1")) {
        spec$predictors <- predictors
        spec$beta_grid <- beta_grid
        spec$variance_grid <- variance_grid
      } else if (identical(family, "gaussian") &&
                 identical(spec$version, "poisson_random_intercept_grid_v1")) {
        spec$predictors <- predictors
        spec$beta_grid <- beta_grid
        spec$variance_grid <- variance_grid
      } else if (identical(family, "gaussian") &&
                 spec$version %in% c("binomial_grid_v1",
                                     "poisson_grid_v1")) {
        spec$predictors <- predictors
        spec$beta_grid <- beta_grid
      } else if (identical(family, "gaussian") &&
                 identical(spec$version, "negative_binomial_grid_v1")) {
        spec$predictors <- predictors
        spec$beta_grid <- beta_grid
        spec$theta_grid <- theta_grid
      } else if (identical(family, "gaussian") &&
                 identical(spec$version, "multinomial_grid_v1")) {
        spec$predictors <- predictors
        spec$levels <- levels
        spec$reference <- reference
        spec$beta_grid <- beta_grid
      } else if (identical(family, "gaussian") &&
                 identical(spec$version, "ordinal_grid_v1")) {
        spec$predictors <- predictors
        spec$ordered_levels <- ordered_levels
        spec$candidate_grid <- candidate_grid
      }
      normalized[[analysis_id]] <- spec
    }
    result[[family]] <- normalized
  }
  .dsvert_joint_dp_client_canonical(result)
}

.dsvert_dp_capsule_manifest_draft <- function(
    draft_json, peer, context) {
  .dsvert_dp_capsule_manifest_rejection(draft_json, "draft")
  draft <- .dsvert_joint_dp_client_decode(
    draft_json, "biomedical capsule policy draft",
    .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_DRAFT_BYTES)
  fields <- c(
    "version", "phase", "peer_name", "peer_identity_pk",
    "peer_pinset_sha256", "domain", "cohort_id",
    "dataset_mapping_mode", "datasets", "workload_fragments", "data_access",
    "patient_derived_metadata", "operation_limit", "request_limit",
    "history_can_deny_operation", "signature")
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
    stop("A peer returned an invalid capsule policy draft", call. = FALSE)
  }
  .dsvert_dp_capsule_manifest_verify(draft, "draft", peer, context)
  workload_fragments <- .dsvert_dp_capsule_manifest_fragments(
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
      .dsvert_dp_is_integer(dataset$alignment_protocol_version, 1) &&
      is.list(dataset$columns) && length(dataset$columns) > 0L &&
      !is.null(names(dataset$columns)) && !anyNA(names(dataset$columns)) &&
      !anyDuplicated(names(dataset$columns))
    if (!isTRUE(valid_dataset)) {
      stop("A peer returned invalid capsule dataset metadata", call. = FALSE)
    }
    columns <- dataset$columns[order(
      names(dataset$columns), method = "radix")]
    if (any(!grepl(
          "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", names(columns))) ||
        any(names(columns) %in% seen)) {
      stop("A peer returned ambiguous capsule column metadata", call. = FALSE)
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
  list(
    value = draft, datasets = normalized,
    workload_fragments = workload_fragments)
}

.dsvert_dp_capsule_manifest_workload_contract <- function(
    drafts, context) {
  peers <- context$servers
  families <- c("describe", "survival", "gaussian", "vertical_cross")
  contract <- stats::setNames(vector("list", length(families)), families)
  for (family in families) {
    entries <- list()
    for (peer in peers) {
      local <- drafts[[peer]]$workload_fragments[[family]]
      duplicate <- intersect(names(entries), names(local))
      if (length(duplicate)) {
        stop("Pinned custodians assigned the same ", family,
             " workload id to different owners", call. = FALSE)
      }
      for (analysis_id in names(local)) {
        entries[[analysis_id]] <- list(
          owner_peer = peer, spec = local[[analysis_id]])
      }
    }
    if (length(entries)) {
      entries <- entries[order(names(entries), method = "radix")]
    }
    contract[[family]] <- entries
  }
  value <- .dsvert_joint_dp_client_canonical(c(
    list(version = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_CONTRACT_VERSION),
    contract))
  json <- .dsvert_joint_dp_client_json(value)
  if (nchar(json, type = "bytes") >
      .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_WORKLOAD_BYTES) {
    stop("The custodian workload contract is too large", call. = FALSE)
  }
  list(
    value = value, json = json,
    sha256 = digest::digest(json, algo = "sha256", serialize = FALSE))
}

.dsvert_dp_capsule_manifest_expected_snapshot <- function(
    context, datasets, alignment_protocol_version, workload_contract) {
  reference <- context$status[[context$servers[[1L]]]]$policy
  fingerprint <- .dsvert_dp_capsule_manifest_hash(list(
    protocol = "dsvert-biomedical-capsule-logical-snapshot-v1",
    domain = reference$domain, cohort_id = reference$cohort_id,
    peer_pinset_sha256 = reference$peer_pinset_sha256,
    alignment_protocol_version = as.numeric(alignment_protocol_version),
    datasets = datasets,
    workload_contract = workload_contract))
  .dsvert_joint_dp_client_canonical(list(
    logical_snapshot_id = reference$cohort_id,
    version = paste0("schema-v1-", fingerprint),
    alignment_protocol_version = as.numeric(alignment_protocol_version)))
}

.dsvert_dp_capsule_manifest_schema <- function(drafts, context) {
  peers <- context$servers
  if (!is.list(drafts) || is.null(names(drafts)) ||
      !setequal(names(drafts), peers)) {
    stop("Capsule policy drafts do not cover every pinned peer",
         call. = FALSE)
  }
  column_names <- unlist(lapply(peers, function(peer) {
    unlist(lapply(drafts[[peer]]$datasets, function(dataset) {
      names(dataset$columns)
    }), use.names = FALSE)
  }), use.names = FALSE)
  duplicate_columns <- unique(column_names[
    duplicated(column_names) | duplicated(column_names, fromLast = TRUE)])
  datasets <- list()
  alignment_versions <- numeric()
  for (peer in peers) {
    for (data_name in names(drafts[[peer]]$datasets)) {
      local <- drafts[[peer]]$datasets[[data_name]]
      alignment_versions <- c(
        alignment_versions, local$alignment_protocol_version)
      dataset_common <- local[c(
        "dataset_id", "dataset_version", "schema_version",
        "alignment_group")]
      if (is.null(datasets[[data_name]])) {
        datasets[[data_name]] <- c(dataset_common, list(
          patient_keys = list(), columns = list()))
      } else if (!identical(
          datasets[[data_name]][names(dataset_common)], dataset_common)) {
        stop("Pinned peers assigned conflicting identities to dataset '",
             data_name, "'", call. = FALSE)
      }
      if (peer %in% names(datasets[[data_name]]$patient_keys)) {
        stop("A peer duplicated its capsule dataset ownership", call. = FALSE)
      }
      datasets[[data_name]]$patient_keys[[peer]] <- local$patient_column
      local_names <- names(local$columns)
      signed_names <- ifelse(
        local_names %in% duplicate_columns,
        paste0(peer, "$", local_names), local_names)
      if (any(signed_names %in% names(datasets[[data_name]]$columns))) {
        stop("Biomedical capsule columns do not identify a unique signed ",
             "owner/dataset/column triplet", call. = FALSE)
      }
      names(local$columns) <- signed_names
      datasets[[data_name]]$columns <- c(
        datasets[[data_name]]$columns, local$columns)
    }
  }
  if (length(unique(alignment_versions)) != 1L) {
    stop("Pinned datasets disagree on the alignment protocol version",
         call. = FALSE)
  }
  datasets <- datasets[order(names(datasets), method = "radix")]
  datasets <- lapply(datasets, function(dataset) {
    dataset$patient_keys <- dataset$patient_keys[
      order(names(dataset$patient_keys), method = "radix")]
    dataset$columns <- dataset$columns[
      order(names(dataset$columns), method = "radix")]
    dataset
  })
  workload_contract <- .dsvert_dp_capsule_manifest_workload_contract(
    drafts, context)
  snapshot <- .dsvert_dp_capsule_manifest_expected_snapshot(
    context, datasets, unique(alignment_versions),
    workload_contract$value)
  reference <- context$status[[peers[[1L]]]]$policy
  unsigned <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = snapshot,
    peer_pinset_sha256 = reference$peer_pinset_sha256,
    datasets = datasets))
  list(
    unsigned = unsigned,
    json = .dsvert_joint_dp_client_json(unsigned),
    sha256 = .dsvert_dp_capsule_manifest_hash(unsigned),
    logical_snapshot = snapshot,
    workload_contract = workload_contract)
}

.dsvert_dp_capsule_manifest_schema_signature <- function(
    response_json, peer, context, schema) {
  .dsvert_dp_capsule_manifest_rejection(response_json, "schema_sign")
  response <- .dsvert_joint_dp_client_decode(
    response_json, "biomedical capsule schema signature",
    64L * 1024L)
  fields <- c(
    "version", "phase", "peer_name", "peer_identity_pk",
    "peer_pinset_sha256", "schema_sha256",
    "workload_contract_sha256", "logical_snapshot", "schema_signature",
    "data_access", "operation_limit", "request_limit",
    "history_can_deny_operation")
  reference <- context$status[[peer]]$policy
  valid <- .dsvert_dp_has_exact_names(response, fields) &&
    identical(response$version,
              .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_SIGN_VERSION) &&
    identical(response$phase, "global_schema_policy_verified") &&
    identical(response$peer_name, peer) &&
    identical(response$peer_identity_pk, unname(context$pinset[[peer]])) &&
    identical(response$peer_pinset_sha256,
              reference$peer_pinset_sha256) &&
    identical(response$schema_sha256, schema$sha256) &&
    identical(response$workload_contract_sha256,
              schema$workload_contract$sha256) &&
    identical(
      .dsvert_joint_dp_client_canonical(response$logical_snapshot),
      schema$logical_snapshot) &&
    identical(response$data_access, FALSE) &&
    identical(response$operation_limit, FALSE) &&
    identical(response$request_limit, FALSE) &&
    identical(response$history_can_deny_operation, FALSE)
  if (!isTRUE(valid)) {
    stop("A peer returned an invalid capsule schema-signing response",
         call. = FALSE)
  }
  .dsvert_dp_capsule_schema_verify(
    schema$unsigned, response$schema_signature, peer, context)
  response$schema_signature
}

.dsvert_dp_capsule_manifest_build_response <- function(
    response_json, peer, context, schema) {
  .dsvert_dp_capsule_manifest_rejection(response_json, "build")
  response <- .dsvert_joint_dp_client_decode(
    response_json, "server-authoritative biomedical capsule manifest",
    .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_BUILD_BYTES)
  fields <- c(
    "version", "phase", "peer_name", "peer_identity_pk",
    "peer_pinset_sha256", "schema_sha256",
    "workload_contract_sha256", "manifest_sha256",
    "manifest_bytes", "capsule_id", "privacy_epoch", "noise_key_id",
    "artifact_commitment_count", "artifact_commitments_root",
    "artifact_commitments", "manifest_json",
    "durable_memoization", "deterministic_replay", "data_access",
    "operation_limit", "request_limit", "history_can_deny_operation",
    "signature")
  reference <- context$status[[peer]]$policy
  noise_root <- context$status[[peer]]$noise_root
  valid <- .dsvert_dp_has_exact_names(response, fields) &&
    identical(response$version,
              .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_BUILD_VERSION) &&
    identical(response$phase,
              "server_authoritative_manifest_memoized") &&
    identical(response$peer_name, peer) &&
    identical(response$peer_identity_pk, unname(context$pinset[[peer]])) &&
    identical(response$peer_pinset_sha256,
              reference$peer_pinset_sha256) &&
    identical(response$schema_sha256, schema$sha256) &&
    identical(response$workload_contract_sha256,
              schema$workload_contract$sha256) &&
    .dsvert_dp_capsule_source_hex(response$manifest_sha256) &&
    .dsvert_dp_capsule_source_hex(response$capsule_id) &&
    .dsvert_dp_is_integer(response$privacy_epoch, 1) &&
    identical(as.numeric(response$privacy_epoch),
              as.numeric(noise_root$privacy_epoch)) &&
    .dsvert_dp_is_string(response$noise_key_id) &&
    identical(response$noise_key_id, noise_root$key_id) &&
    .dsvert_dp_is_integer(response$artifact_commitment_count, 0,
                          .DSVERT_DP_MAX_COORDINATES) &&
    .dsvert_dp_capsule_source_hex(response$artifact_commitments_root) &&
    is.list(response$artifact_commitments) &&
    .dsvert_dp_is_string(response$manifest_json) &&
    .dsvert_dp_is_integer(
      response$manifest_bytes, 1,
      .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_BUILD_BYTES) &&
    identical(as.numeric(response$manifest_bytes), as.numeric(nchar(
      response$manifest_json, type = "bytes"))) &&
    identical(response$manifest_sha256, digest::digest(
      response$manifest_json, algo = "sha256", serialize = FALSE)) &&
    identical(response$durable_memoization, TRUE) &&
    identical(response$deterministic_replay, TRUE) &&
    identical(response$data_access, FALSE) &&
    identical(response$operation_limit, FALSE) &&
    identical(response$request_limit, FALSE) &&
    identical(response$history_can_deny_operation, FALSE)
  if (!isTRUE(valid)) {
    stop("A peer returned an invalid authoritative capsule manifest",
         call. = FALSE)
  }
  .dsvert_dp_capsule_manifest_verify(response, "build", peer, context)
  manifest <- .dsvert_dp_capsule_source_manifest(
    response$manifest_json, context)
  artifact_index <- .dsvert_dp_capsule_artifact_commitment_index_client(
    manifest$value, reference, response$manifest_sha256)
  attestation <- manifest$value$workload$schema_attestation
  attested_signers <- tryCatch(
    .dsvert_dp_capsule_manifest_string_array(
      attestation$signers, "schema attestation signers"),
    error = function(error) character())
  if (!identical(manifest$capsule_id, response$capsule_id) ||
      !identical(
        attestation$manifest_sha256,
        schema$sha256) ||
      !identical(
        .dsvert_joint_dp_client_canonical(attestation$signatures),
        .dsvert_joint_dp_client_canonical(schema$signatures)) ||
      !identical(attested_signers,
                 sort(context$servers, method = "radix")) ||
      !identical(
        .dsvert_joint_dp_client_canonical(
          manifest$value$logical_snapshot),
        schema$logical_snapshot) ||
      !identical(
        .dsvert_joint_dp_client_canonical(response$artifact_commitments),
        artifact_index$value) ||
      !identical(as.numeric(response$artifact_commitment_count),
                 artifact_index$count) ||
      !identical(response$artifact_commitments_root,
                 artifact_index$root)) {
    stop("The authoritative manifest is not bound to its signed schema",
         call. = FALSE)
  }
  list(value = response, manifest = manifest,
       artifact_index = artifact_index)
}

#' Build one server-authoritative biomedical capsule manifest
#'
#' This internal orchestration has no analyst-controlled schema arguments.  It
#' validates signed policy drafts from the complete pinset, derives the only
#' accepted global schema and logical snapshot, collects every schema
#' signature, and requires every peer to return the same signed manifest bytes.
#'
#' @param datasources Complete named DataSHIELD connection set.
#' @param status Optional already validated reusable-capsule status handshake.
#' @param .aggregate Injectable DSI aggregate implementation for tests.
#' @return A list containing the canonical signed schema, canonical manifest,
#'   capsule identity and validated connection context.
#' @keywords internal
.dsvert_dp_capsule_manifest_build <- function(
    datasources, status = NULL, .aggregate = DSI::datashield.aggregate) {
  context <- .dsvert_joint_dp_client_context(
    datasources, status = status, .aggregate = .aggregate)
  peers <- context$servers
  draft_calls <- stats::setNames(lapply(peers, function(peer) {
    call(name = "dsvertDPCapsuleManifestDraftDS")
  }), peers)
  draft_json <- .dsvert_fanout_by_site(
    context$all_conns, draft_calls,
    operation = "custodian capsule-policy draft fan-out",
    .aggregate = .aggregate)
  drafts <- stats::setNames(lapply(peers, function(peer) {
    .dsvert_dp_capsule_manifest_draft(
      draft_json[[peer]], peer, context)
  }), peers)
  schema <- .dsvert_dp_capsule_manifest_schema(drafts, context)

  sign_calls <- stats::setNames(lapply(peers, function(peer) call(
    name = "dsvertDPCapsuleManifestSignDS",
    schema_json = schema$json,
    workload_contract_json = schema$workload_contract$json)), peers)
  sign_json <- .dsvert_fanout_by_site(
    context$all_conns, sign_calls,
    operation = "global capsule-schema signature fan-out",
    .aggregate = .aggregate)
  signatures <- stats::setNames(lapply(peers, function(peer) {
    .dsvert_dp_capsule_manifest_schema_signature(
      sign_json[[peer]], peer, context, schema)
  }), peers)
  schema$signatures <- signatures
  signed_schema <- c(schema$unsigned, list(signatures = signatures))
  signed_schema <- .dsvert_joint_dp_client_canonical(signed_schema)
  signed_schema_json <- .dsvert_joint_dp_client_json(signed_schema)

  build_calls <- stats::setNames(lapply(peers, function(peer) call(
    name = "dsvertDPCapsuleManifestBuildDS",
    schema_json = signed_schema_json,
    workload_contract_json = schema$workload_contract$json)), peers)
  build_json <- .dsvert_fanout_by_site(
    context$all_conns, build_calls,
    operation = "server-authoritative capsule manifest fan-out",
    .aggregate = .aggregate)
  builds <- stats::setNames(lapply(peers, function(peer) {
    .dsvert_dp_capsule_manifest_build_response(
      build_json[[peer]], peer, context, schema)
  }), peers)
  manifest_json <- builds[[1L]]$value$manifest_json
  if (!all(vapply(builds, function(value) {
        identical(value$value$manifest_json, manifest_json) &&
          identical(value$value$manifest_sha256,
                    builds[[1L]]$value$manifest_sha256) &&
          identical(value$value$capsule_id,
                    builds[[1L]]$value$capsule_id) &&
          identical(value$value$artifact_commitments_root,
                    builds[[1L]]$value$artifact_commitments_root) &&
          identical(as.numeric(value$value$artifact_commitment_count),
                    as.numeric(builds[[1L]]$value$
                      artifact_commitment_count)) &&
          identical(
            .dsvert_joint_dp_client_canonical(
              value$value$artifact_commitments),
            .dsvert_joint_dp_client_canonical(
              builds[[1L]]$value$artifact_commitments))
      }, logical(1L)))) {
    stop("Pinned peers built different biomedical capsule manifests",
         call. = FALSE)
  }
  list(
    schema_json = signed_schema_json,
    schema_sha256 = schema$sha256,
    workload_contract_json = schema$workload_contract$json,
    workload_contract_sha256 = schema$workload_contract$sha256,
    logical_snapshot = schema$logical_snapshot,
    manifest_json = manifest_json,
    manifest_sha256 = builds[[1L]]$value$manifest_sha256,
    capsule_id = builds[[1L]]$value$capsule_id,
    artifact_commitments = builds[[1L]]$artifact_index$value,
    artifact_commitment_count = builds[[1L]]$artifact_index$count,
    artifact_commitments_root = builds[[1L]]$artifact_index$root,
    manifest_build_receipts = stats::setNames(lapply(peers, function(peer) {
      builds[[peer]]$value[setdiff(
        names(builds[[peer]]$value),
        c("manifest_json", "artifact_commitments"))]
    }), peers),
    manifest_signatures = stats::setNames(lapply(peers, function(peer) {
      builds[[peer]]$value$signature
    }), peers),
    deterministic_replay = TRUE,
    operation_or_request_limit = FALSE,
    history_can_deny_operation = FALSE,
    context = context)
}
