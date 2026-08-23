# Pure client boundary for a synopsis vector that has already completed the
# bilateral server lifecycle.  The verified manifest bundle is deliberately
# retained as the trust anchor; a bare relay-supplied manifest is never enough.

.DSVERT_CLIENT_SYNOPSIS_COMPILE_VERSION <-
  "dsvert-stateless-catalog-synopsis-compile-receipt-v1"
.DSVERT_CLIENT_SYNOPSIS_COMPILE_SET_VERSION <-
  "dsvert-stateless-catalog-synopsis-compile-receipt-set-v1"
.DSVERT_CLIENT_SYNOPSIS_RELEASE_VERSION <-
  "dsvert-stateless-catalog-synopsis-release-v1"
.DSVERT_CLIENT_SYNOPSIS_REPLAY_VERSION <-
  "dsvert-stateless-catalog-synopsis-replay-v1"
.DSVERT_CLIENT_SYNOPSIS_PUBLIC_VERSION <-
  "dsvert-stateless-catalog-synopsis-public-chunk-v1"
.DSVERT_CLIENT_SYNOPSIS_PUBLIC_CHUNK_COORDINATES <- 8192L
.DSVERT_CLIENT_SYNOPSIS_EXACT_CHUNK_COORDINATES <- 64L
# Per-object availability bounds. They constrain one physical protocol object;
# they are not call, catalogue, rate or lifetime privacy limits.
.DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES <- 32L * 1024L^2
.DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES <- 2L * 1024L^2
.DSVERT_CLIENT_SYNOPSIS_PREPARE_MAX_OBJECT_BYTES <- 64L * 1024L
.DSVERT_CLIENT_SYNOPSIS_FRACTION_MAX_BYTES <- 512L

.dsvert_dp_synopsis_client_hash <- function(domain, value) {
  digest::digest(charToRaw(paste0(
    domain, .dsvert_joint_dp_client_json(
      .dsvert_joint_dp_client_canonical(value)))),
  algo = "sha256", serialize = FALSE)
}

.dsvert_dp_synopsis_client_json <- function(
    value, what,
    maximum_bytes = .DSVERT_CLIENT_SYNOPSIS_MAX_OBJECT_BYTES) {
  parsed <- .dsvert_joint_dp_client_decode(
    value, what, maximum_bytes)
  if (!identical(.dsvert_joint_dp_client_json(parsed), value)) {
    stop("The ", what, " is not canonical", call. = FALSE)
  }
  parsed
}

.dsvert_dp_synopsis_client_signature <- function(
    unsigned, signature, identity_pk, domain, what) {
  public <- .dsvert_joint_dp_client_b64url(
    identity_pk, 32L, paste(what, "identity key"))
  signature <- .dsvert_joint_dp_client_b64url(
    signature, 64L, paste(what, "signature"))
  message <- charToRaw(paste0(
    domain, .dsvert_joint_dp_client_json(unsigned)))
  valid <- tryCatch(openssl::ed25519_verify(
    message, signature, openssl::read_ed25519_pubkey(public)),
  error = function(error) FALSE)
  if (!isTRUE(valid)) stop("Invalid ", what, " signature", call. = FALSE)
  invisible(TRUE)
}

.dsvert_dp_synopsis_client_decimal_fraction <- function(value, what) {
  if (!.dsvert_vector_string(
        value, "^[0-9](?:\\.[0-9]+)?e[+-][0-9]{2,3}$", 64L)) {
    stop("Invalid synopsis ", what, call. = FALSE)
  }
  pieces <- strsplit(value, "e", fixed = TRUE)[[1L]]
  mantissa <- pieces[[1L]]
  exponent <- as.integer(pieces[[2L]])
  decimals <- if (grepl(".", mantissa, fixed = TRUE)) {
    nchar(sub("^[^.]*\\.", "", mantissa))
  } else 0L
  digits <- sub(".", "", mantissa, fixed = TRUE)
  numerator <- openssl::bignum(digits)
  power <- exponent - decimals
  if (power >= 0L) {
    numerator <- numerator * (openssl::bignum(10) ^ power)
    denominator <- openssl::bignum(1)
  } else {
    denominator <- openssl::bignum(10) ^ (-power)
  }
  list(numerator = numerator, denominator = denominator)
}

.dsvert_dp_synopsis_client_fraction_leq <- function(
    numerator, denominator, decimal, what) {
  if (!.dsvert_vector_integer_text(
        numerator, maximum_bytes = .DSVERT_CLIENT_SYNOPSIS_FRACTION_MAX_BYTES) ||
      !.dsvert_vector_integer_text(
        denominator, positive = TRUE,
        maximum_bytes = .DSVERT_CLIENT_SYNOPSIS_FRACTION_MAX_BYTES)) {
    stop("Invalid synopsis ", what, call. = FALSE)
  }
  right <- .dsvert_dp_synopsis_client_decimal_fraction(decimal, what)
  openssl::bignum(numerator) * right$denominator <=
    right$numerator * openssl::bignum(denominator)
}

.dsvert_dp_synopsis_client_lattice_scale <- function(bits) {
  if (!.dsvert_vector_whole(bits, 1, 62)) {
    stop("Invalid synopsis lattice scale", call. = FALSE)
  }
  as.character(openssl::bignum(2) ^ as.integer(bits))
}

.dsvert_dp_synopsis_client_bundle <- function(manifest_bundle, status) {
  fields <- c(
    "schema_json", "schema_sha256", "workload_contract_json",
    "workload_contract_sha256", "logical_snapshot", "manifest_json",
    "manifest_sha256", "capsule_id", "artifact_commitments",
    "artifact_commitment_count", "artifact_commitments_root",
    "manifest_build_receipts", "manifest_signatures",
    "deterministic_replay", "operation_or_request_limit",
    "history_can_deny_operation", "context")
  if (!.dsvert_dp_has_exact_names(manifest_bundle, fields) ||
      !is.list(manifest_bundle$context) ||
      !is.list(manifest_bundle$context$all_conns)) {
    stop("A trusted manifest bundle with a verified context is required",
         call. = FALSE)
  }
  if (inherits(status, "ds.vertDPSynopsisStatus")) {
    return(.dsvert_dp_synopsis_trusted_bundle_v1(manifest_bundle, status))
  }
  context <- manifest_bundle$context
  consensus <- .dsvert_joint_dp_capsule_status_consensus(
    status, context$all_conns)
  if (!identical(unclass(consensus), unclass(context$status)) ||
      !identical(context$servers,
                 sort(names(context$all_conns), method = "radix")) ||
      !identical(context$pinset,
                 consensus[[context$servers[[1L]]]]$policy$peer_pinset) ||
      !identical(context$designated,
                 consensus[[context$servers[[1L]]]]$policy$
                   designated_noise_peers)) {
    stop("The synopsis status differs from its trusted context",
         call. = FALSE)
  }
  if (!.dsvert_vector_hex(manifest_bundle$manifest_sha256) ||
      !identical(manifest_bundle$manifest_sha256, digest::digest(
        manifest_bundle$manifest_json, "sha256", serialize = FALSE)) ||
      !identical(manifest_bundle$deterministic_replay, TRUE) ||
      !identical(manifest_bundle$operation_or_request_limit, FALSE) ||
      !identical(manifest_bundle$history_can_deny_operation, FALSE)) {
    stop("The trusted synopsis manifest bundle is invalid", call. = FALSE)
  }

  signed_schema <- .dsvert_joint_dp_client_decode(
    manifest_bundle$schema_json, "signed synopsis capsule schema",
    .DSVERT_CLIENT_DP_CAPSULE_MANIFEST_MAX_SCHEMA_BYTES)
  signatures <- if (is.list(signed_schema)) signed_schema$signatures else NULL
  unsigned_schema <- signed_schema
  if (is.list(unsigned_schema)) unsigned_schema$signatures <- NULL
  if (!is.list(unsigned_schema) || !is.list(signatures) ||
      !identical(manifest_bundle$schema_sha256,
                 .dsvert_vector_hash(unsigned_schema)) ||
      !identical(manifest_bundle$workload_contract_sha256,
                 digest::digest(
                   manifest_bundle$workload_contract_json,
                   "sha256", serialize = FALSE)) ||
      is.null(names(signatures)) ||
      !setequal(names(signatures), context$servers)) {
    stop("The trusted signed synopsis schema is invalid", call. = FALSE)
  }
  for (peer in context$servers) {
    .dsvert_dp_capsule_schema_verify(
      unsigned_schema, signatures[[peer]], peer, context)
  }
  schema <- list(
    unsigned = unsigned_schema,
    sha256 = manifest_bundle$schema_sha256,
    logical_snapshot = manifest_bundle$logical_snapshot,
    workload_contract = list(
      json = manifest_bundle$workload_contract_json,
      sha256 = manifest_bundle$workload_contract_sha256),
    signatures = signatures)
  peers <- context$servers
  receipts <- manifest_bundle$manifest_build_receipts
  if (!is.list(receipts) || is.null(names(receipts)) ||
      !setequal(names(receipts), peers) || anyDuplicated(names(receipts))) {
    stop("The manifest build receipt coverage is invalid", call. = FALSE)
  }
  if (!is.list(manifest_bundle$manifest_signatures) ||
      is.null(names(manifest_bundle$manifest_signatures)) ||
      !setequal(names(manifest_bundle$manifest_signatures), peers) ||
      !all(vapply(peers, function(peer) identical(
        manifest_bundle$manifest_signatures[[peer]],
        receipts[[peer]]$signature), logical(1L)))) {
    stop("The manifest build signature set is inconsistent", call. = FALSE)
  }
  builds <- lapply(peers, function(peer) {
    receipt <- receipts[[peer]]
    receipt$manifest_json <- manifest_bundle$manifest_json
    receipt$artifact_commitments <- manifest_bundle$artifact_commitments
    .dsvert_dp_capsule_manifest_build_response(
      .dsvert_joint_dp_client_json(receipt), peer, context, schema)
  })
  reference <- builds[[1L]]
  if (!all(vapply(builds, function(build) {
        identical(build$value$manifest_sha256,
                  manifest_bundle$manifest_sha256) &&
          identical(build$value$capsule_id, manifest_bundle$capsule_id) &&
          identical(build$artifact_index$count,
                    manifest_bundle$artifact_commitment_count) &&
          identical(build$artifact_index$root,
                    manifest_bundle$artifact_commitments_root)
      }, logical(1L)))) {
    stop("The manifest build receipts do not agree", call. = FALSE)
  }
  list(
    context = context, status = consensus,
    manifest = reference$manifest$value,
    validated_manifest = reference$manifest)
}

.dsvert_dp_synopsis_client_catalog_ids <- function(value) {
  if (!is.list(value)) return(value)
  if (!is.null(names(value))) {
    changed <- names(value)
    changed[changed == "analysis_id"] <- "catalog_entry_id"
    if (anyDuplicated(changed)) {
      stop("The synopsis catalog ID namespace collides", call. = FALSE)
    }
    names(value) <- changed
  }
  lapply(value, .dsvert_dp_synopsis_client_catalog_ids)
}

.dsvert_dp_synopsis_client_projection <- function(policy, manifest, layout) {
  catalog <- list(
    version = "dsvert-stateless-catalog-synopsis-catalog-v1",
    domain = policy$domain, cohort_id = policy$cohort_id,
    peer_pinset_sha256 = policy$peer_pinset_sha256,
    logical_snapshot = manifest$logical_snapshot,
    capsule_schema = manifest$capsule_schema,
    schema_manifest_sha256 =
      manifest$workload$schema_attestation$manifest_sha256,
    admission = manifest$admission, bounds = manifest$bounds,
    families = manifest$workload$families,
    vertical_crosses = manifest$workload$vertical_crosses,
    primitive_scope = manifest$workload$primitive_scope,
    release_lattice = manifest$workload$release_lattice,
    sensitivity = manifest$workload$sensitivity,
    coordinate_count = manifest$workload$coordinate_count,
    coordinate_order_sha256 = layout$sha256,
    clipping_sha256 = manifest$workload$capsule_mechanism$clipping_hash)
  catalog <- .dsvert_dp_synopsis_client_catalog_ids(catalog)
  list(
    version = "dsvert-stateless-catalog-synopsis-projection-v1",
    sha256 = .dsvert_dp_synopsis_client_hash(
      "dsVert/stateless-catalog-synopsis/catalog/v1|", catalog),
    catalog = catalog)
}

.dsvert_dp_synopsis_client_lattice <- function(manifest, layout) {
  lattice <- manifest$workload$release_lattice
  bits <- as.integer(lattice$output_lattice_bits)
  if (!.dsvert_dp_is_integer(bits, 1L, 62L) ||
      !identical(as.numeric(lattice$output_lattice_scale), 2^bits)) {
    stop("The synopsis output lattice is invalid", call. = FALSE)
  }
  shifts <- rep(bits, layout$coordinate_count)
  upper <- numeric(layout$coordinate_count)
  for (block in layout$blocks) {
    indices <- seq.int(block$start, block$end)
    vector_bound <- FALSE
    if (identical(block$family, "numeric_moments")) {
      shifts[indices] <- c(bits, 0L, 0L)
      bound <- as.numeric(block$descriptor$statistic_maximum)
      vector_bound <- TRUE
    } else if (identical(block$family, "numeric_pair_moments")) {
      shifts[indices] <- c(bits, rep(0L, 5L))
      bound <- as.numeric(block$descriptor$statistic_maximum)
      vector_bound <- TRUE
    } else if (identical(block$family, "gaussian_models")) {
      already_scaled <- identical(
        block$descriptor$source_coordinate_scaling,
        "all_coordinates_already_on_common_numeric_lattice_v1")
      random_intercept <- identical(
        block$descriptor$version,
        "bounded-normalized-random-intercept-moments-v1")
      shifts[indices] <- if (random_intercept) {
        c(rep(bits, 3L), rep(0L, length(indices) - 3L))
      } else if (already_scaled) rep(0L, length(indices)) else
        c(bits, rep(0L, length(indices) - 1L))
      bound <- as.numeric(block$descriptor$statistic_maximum)
      vector_bound <- TRUE
    } else if (identical(block$family, "categorical_pairs") &&
               identical(block$descriptor$version,
                         "fixed-domain-categorical-cross-contingency-v1")) {
      shifts[indices] <- 0L
      bound <- as.numeric(block$descriptor$statistic_maximum)
      vector_bound <- TRUE
    } else {
      bound <- if (identical(block$family, "categorical_pairs")) {
        sets <- manifest$workload$families$categorical_pairs$sets
        matching <- names(sets)[vapply(sets, function(set) {
          columns <- vapply(set$columns, `[[`, character(1L), "column")
          identical(set$owner_peer, block$owner_peer) &&
            identical(set$dataset, block$dataset) &&
            all(c(block$descriptor$left$column,
                  block$descriptor$right$column) %in% columns)
        }, logical(1L))]
        if (length(matching) != 1L) {
          stop("The synopsis categorical bound is ambiguous", call. = FALSE)
        }
        sets[[matching]]$statistic_maximum
      } else block$descriptor$statistic_maximum
    }
    if (!is.numeric(bound) ||
        (isTRUE(vector_bound) && length(bound) != length(indices)) ||
        (!isTRUE(vector_bound) && length(bound) != 1L) ||
        anyNA(bound) || any(!is.finite(bound)) || any(bound < 0) ||
        any(bound != floor(bound))) {
      stop("The synopsis coordinate bounds are invalid", call. = FALSE)
    }
    upper[indices] <- bound
  }
  upper_text <- vapply(upper, function(value) format(
    value, scientific = FALSE, trim = TRUE, digits = 22L), character(1L))
  norm <- if (identical(
      manifest$workload$capsule_mechanism$mechanism,
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM)) "l2" else "l1"
  steps <- if (identical(norm, "l2")) {
    format(lattice$integer_l2_sensitivity_steps,
           scientific = TRUE, trim = TRUE, digits = 17L)
  } else {
    format(lattice$integer_l1_sensitivity_steps,
           scientific = FALSE, trim = TRUE, digits = 22L)
  }
  transform <- list(
    version = lattice$version, transform_rule = lattice$transform_rule,
    output_lattice_bits = bits, scale_shifts = as.list(shifts),
    raw_upper_bounds = as.list(upper_text),
    sensitivity_norm = norm, sensitivity_steps = steps)
  list(
    output_lattice_bits = bits,
    output_lattice_scale = as.numeric(lattice$output_lattice_scale),
    scale_shifts = shifts, raw_upper_bounds = upper_text,
    sensitivity_norm = norm, sensitivity_steps = steps,
    transform_sha256 = .dsvert_vector_hash(transform))
}

.dsvert_dp_synopsis_client_compile <- function(
    compilation, trusted, manifest_bundle) {
  fields <- c("version", "artifact", "receipts", "receipt_set_sha256")
  if (!.dsvert_dp_has_exact_names(compilation, fields) ||
      !identical(compilation$version,
                 .DSVERT_CLIENT_SYNOPSIS_COMPILE_SET_VERSION)) {
    stop("Invalid synopsis compile receipt set", call. = FALSE)
  }
  artifact <- compilation$artifact
  if (!.dsvert_dp_has_exact_names(
        artifact, c("semantic", "artifact_key", "physical_plan")) ||
      !.dsvert_dp_has_exact_names(artifact$semantic, c(
        "version", "catalog_projection", "source_claim_set_sha256",
        "noise_authority_roles", "privacy", "release", "public_shape")) ||
      !identical(artifact$semantic$version,
                 "dsvert-analysis-semantic-stateless-catalog-synopsis-v1") ||
      !.dsvert_vector_hex(artifact$semantic$source_claim_set_sha256) ||
      !identical(artifact$artifact_key,
        .dsvert_dp_synopsis_client_hash(
          "dsVert/analysis-artifact-key/v1|", artifact$semantic))) {
    stop("The synopsis artifact key is invalid", call. = FALSE)
  }
  physical <- artifact$physical_plan
  physical_fields <- c(
    "version", "request", "profile", "lattice", "backend_selection",
    "draw_law", "draw_law_sha256", "full_plan", "full_plan_sha256")
  if (!.dsvert_dp_has_exact_names(physical, physical_fields) ||
      !identical(physical$version,
                 "dsvert-stateless-catalog-synopsis-physical-plan-v1") ||
      !identical(physical$full_plan_sha256,
                 .dsvert_vector_hash(physical$full_plan)) ||
      !identical(physical$draw_law_sha256,
        .dsvert_dp_synopsis_client_hash(
          "dsVert/stateless-catalog-synopsis/draw-law/v1|",
          physical$draw_law)) ||
      !identical(
        .dsvert_joint_dp_client_json(artifact$semantic$release),
        .dsvert_joint_dp_client_json(physical[setdiff(
          physical_fields, c("full_plan", "full_plan_sha256"))]))) {
    stop("The synopsis physical plan is detached", call. = FALSE)
  }
  context <- trusted$context
  policy <- trusted$status[[context$servers[[1L]]]]$policy
  layout <- .dsvert_dp_capsule_vector_layout(trusted$manifest)
  projection <- .dsvert_dp_synopsis_client_projection(
    policy, trusted$manifest, layout)
  if (!identical(
        .dsvert_joint_dp_client_json(
          artifact$semantic$catalog_projection),
        .dsvert_joint_dp_client_json(projection))) {
    stop("The synopsis catalog projection changed", call. = FALSE)
  }
  expected_roles <- list(
    version = "dsvert-synopsis-noise-authority-roles-v1",
    role_order = list(
      "primary_noise_authority", "secondary_noise_authority"),
    authority_ids = as.list(unname(context$pinset[context$designated])))
  if (!identical(
      .dsvert_joint_dp_client_json(
        artifact$semantic$noise_authority_roles),
      .dsvert_joint_dp_client_json(expected_roles))) {
    stop("The synopsis noise-authority roles are invalid", call. = FALSE)
  }
  lattice <- .dsvert_dp_synopsis_client_lattice(trusted$manifest, layout)
  expected_lattice <- list(
    version = "dsvert-stateless-catalog-synopsis-lattice-v1",
    coordinate_count = as.integer(layout$coordinate_count),
    coordinate_order_sha256 = layout$sha256,
    clipping_sha256 =
      trusted$manifest$workload$capsule_mechanism$clipping_hash,
    transform_sha256 = lattice$transform_sha256,
    output_lattice_bits = lattice$output_lattice_bits,
    output_lattice_scale = lattice$output_lattice_scale,
    sensitivity_norm = lattice$sensitivity_norm,
    sensitivity_steps = lattice$sensitivity_steps,
    ring_bits = 128L, fractional_bits = 0L)
  if (!identical(
        .dsvert_joint_dp_client_json(physical$lattice),
        .dsvert_joint_dp_client_json(expected_lattice)) ||
      !identical(as.numeric(physical$request$total_coordinate_count),
                 as.numeric(layout$coordinate_count)) ||
      !identical(
        .dsvert_joint_dp_client_json(artifact$semantic$public_shape),
        .dsvert_joint_dp_client_json(list(
          version = "dsvert-stateless-catalog-synopsis-shape-v1",
          coordinates = as.integer(layout$coordinate_count))))) {
    stop("The synopsis lattice or public shape is invalid", call. = FALSE)
  }
  profile <- .dsvert_vector_profile(
    trusted$manifest$workload$capsule_mechanism,
    mechanism_selection =
      trusted$manifest$workload$mechanism_selection,
    backend = physical$profile$backend)
  expected_profile <- list(
    version = "dsvert-stateless-catalog-synopsis-backend-profile-v1",
    mechanism = profile$mechanism, backend = profile$backend,
    plan_version = profile$plan_version, sampler = profile$sampler,
    release_mechanism = profile$release_mechanism,
    draw_count = if (isTRUE(profile$exact_gc)) 1L else 2L,
    complete_epsilon_per_draw = TRUE,
    delta_aggregation = profile$delta_aggregation,
    adversary_model = if (isTRUE(profile$exact_gc)) {
      "analyst_plus_at_most_one_exact_gc_authority_v1"
    } else {
      paste0("analyst_plus_at_most_one_of_two_noncolluding_",
             "noise_authorities_v1")
    },
    output_transform = profile$postprocessing,
    commitment_purpose = if (isTRUE(profile$gaussian)) {
      "dyadic-discrete-gaussian-tv-bounded-v2"
    } else if (isTRUE(profile$exact_gc)) {
      "one-draw-exact-gc"
    } else "convolution")
  expected_selection <- if (isTRUE(profile$gaussian)) {
    list(
      version = "dsvert-stateless-catalog-synopsis-backend-selection-v1",
      rule = "manifest_certified_fixed_work_gaussian_v1",
      total_coordinate_count = as.integer(layout$coordinate_count),
      backend = profile$backend, selected_before_private_material = TRUE,
      retry_may_change_backend = FALSE)
  } else {
    promoted <- layout$coordinate_count <=
      .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_MAX_PROMOTED_COORDINATES
    list(
      version = "dsvert-stateless-catalog-synopsis-backend-selection-v1",
      rule = "public_coordinate_ceiling_v1",
      selected_before_private_material = TRUE,
      retry_may_change_backend = FALSE,
      policy_version =
        .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_COST_POLICY_VERSION,
      total_coordinate_count = as.integer(layout$coordinate_count),
      maximum_promoted_coordinates =
        .DSVERT_CLIENT_JOINT_DP_VECTOR_EXACT_GC_MAX_PROMOTED_COORDINATES,
      promoted = promoted,
      backend = if (promoted) .DSVERT_CLIENT_VECTOR_EXACT_BACKEND else
        .DSVERT_CLIENT_VECTOR_BACKEND,
      selection_reason = if (promoted) {
        "within_public_exact_gc_cost_ceiling"
      } else "above_public_exact_gc_cost_ceiling")
  }
  if (!identical(
        .dsvert_joint_dp_client_json(physical$profile),
        .dsvert_joint_dp_client_json(expected_profile)) ||
      !identical(
        .dsvert_joint_dp_client_json(physical$backend_selection),
        .dsvert_joint_dp_client_json(expected_selection))) {
    stop("The synopsis backend profile or selection is invalid",
         call. = FALSE)
  }
  draw_fields <- if (isTRUE(profile$gaussian)) {
    c(
      "version", "mechanism", "sampler", "reference",
      "total_coordinate_count", "request_binding_sha256",
      "epsilon_numerator", "epsilon_denominator",
      "allocated_delta_numerator", "allocated_delta_denominator",
      "core_delta_numerator", "core_delta_denominator",
      "tail_delta_numerator", "tail_delta_denominator",
      "l2_sensitivity_numerator", "l2_sensitivity_denominator",
      "rho_numerator", "rho_denominator", "zcdp_log_upper_integer",
      "zcdp_conversion_exponent_numerator",
      "zcdp_conversion_exponent_denominator",
      "sigma_squared_numerator", "sigma_squared_denominator",
      "proposal_scale", "maximum_noise_magnitude_per_peer",
      "maximum_noise_magnitude_two_peers",
      "tail_proof_exponent_numerator", "tail_proof_exponent_denominator",
      "tail_proof_target_numerator", "tail_proof_target_denominator",
      "vector_tail_tv_upper_numerator", "vector_tail_tv_upper_denominator",
      "vector_sampler_tv_upper_numerator",
      "vector_sampler_tv_upper_denominator",
      "vector_total_tv_upper_numerator", "vector_total_tv_upper_denominator",
      "per_peer_implementation_delta_numerator",
      "per_peer_implementation_delta_denominator",
      "sampler_candidate_count", "sampler_random_bits_per_coordinate",
      "sampler_table_precision_bits", "sampler_magnitude_count",
      "independent_noise_peer_count", "complete_epsilon_per_peer",
      "epsilon_divided_by_peer_count", "release_delta_aggregation",
      "nominal_variance_multiplier", "nominal_standard_deviation_factor",
      "at_least_one_honest_noise_peer", "maximum_colluding_noise_peers",
      "exact_rational_sampler", "finite_support_transfer_charged",
      "fixed_work_sampler", "sampler_branches_on_protected_values",
      "sampler_branches_on_private_randomness", "transcript_dp_claim",
      "logical_transcript_fixed_shape")
  } else {
    fields <- c(
      "version", "sampler", "stop_bits", "stop_numerator", "uniform_bits",
      "binary_geometric_bits", "bernoulli_thresholds", "sensitivity_steps",
      "total_coordinate_count", "epsilon_effective_upper_numerator",
      "epsilon_effective_upper_denominator", "one_geometric_tv_numerator",
      "one_geometric_tv_denominator", "tail_upper_numerator",
      "tail_upper_denominator", "rounding_upper_numerator",
      "rounding_upper_denominator", "implementation_delta_numerator",
      "implementation_delta_denominator", "maximum_noise_magnitude")
    if (!isTRUE(profile$exact_gc)) c(
      fields, "independent_noise_peer_count", "complete_epsilon_per_peer",
      "epsilon_divided_by_peer_count",
      "geometric_variables_per_peer_per_coordinate",
      "geometric_variables_total_per_coordinate",
      "per_peer_implementation_delta_numerator",
      "per_peer_implementation_delta_denominator",
      "release_implementation_delta_aggregation",
      "two_peer_ideal_transfer_delta_numerator",
      "two_peer_ideal_transfer_delta_denominator") else fields
  }
  expected_draw <- physical$full_plan[draw_fields]
  if ("bernoulli_thresholds" %in% draw_fields) {
    expected_draw$bernoulli_thresholds <- as.list(unname(unlist(
      expected_draw$bernoulli_thresholds, use.names = FALSE)))
  }
  if (!identical(
      .dsvert_joint_dp_client_json(physical$draw_law),
      .dsvert_joint_dp_client_json(expected_draw))) {
    stop("The synopsis draw law changed from its signed full plan",
         call. = FALSE)
  }
  privacy <- artifact$semantic$privacy
  mechanism <- if (is.list(privacy)) privacy$mechanism else NULL
  randomness <- if (is.list(mechanism)) mechanism$randomness else NULL
  lane <- if (is.list(randomness) && is.list(randomness$lanes)) {
    randomness$lanes$final_noise
  } else NULL
  epsilon_number <- if (is.list(privacy)) {
    .dsvert_vector_decimal(privacy$epsilon, 0, 8, open_minimum = TRUE)
  } else NULL
  delta_number <- if (is.list(privacy)) {
    .dsvert_vector_decimal(
      privacy$delta, 0, 1, open_minimum = isTRUE(profile$gaussian))
  } else NULL
  if (!isTRUE(profile$gaussian) && identical(delta_number, 0)) {
    stop(paste(
      "The finite synopsis Laplace v3 backend has positive implementation",
      "delta and cannot certify pure DP"), call. = FALSE)
  }
  family <- if (isTRUE(profile$gaussian)) "gaussian" else
    "discrete_laplace"
  expected_sensitivity <- list(
    version = "dsvert-sensitivity-v1", norm = lattice$sensitivity_norm,
    steps = if (isTRUE(profile$gaussian)) {
      physical$request$l2_sensitivity_steps
    } else physical$lattice$sensitivity_steps)
  expected_lane <- list(
    version = "dsvert-randomness-lane-v1",
    purpose = "privatize_final_vector", primitive = profile$sampler,
    coordinates = as.integer(layout$coordinate_count))
  privacy_valid <-
    .dsvert_dp_has_exact_names(privacy, c(
      "version", "adjacency", "privacy_unit", "epsilon", "delta",
      "mechanism")) &&
    identical(privacy$version, "dsvert-per-synopsis-dp-v1") &&
    identical(privacy$privacy_unit, "patient") &&
    identical(privacy$adjacency, projection$catalog$admission$adjacency) &&
    !is.null(epsilon_number) && !is.null(delta_number) &&
    .dsvert_dp_num_equal(epsilon_number, if (!is.null(
      policy$artifact_epsilon)) policy$artifact_epsilon else
        policy$capsule_epsilon, 128) &&
    .dsvert_dp_num_equal(delta_number, if (!is.null(
      policy$artifact_delta)) policy$artifact_delta else
        policy$capsule_delta, 128) &&
    identical(privacy$epsilon, formatC(
      epsilon_number, digits = 18L, format = "e", decimal.mark = ".")) &&
    identical(privacy$delta, formatC(
      delta_number, digits = 18L, format = "e", decimal.mark = ".")) &&
    .dsvert_dp_has_exact_names(mechanism, c(
      "version", "family", "sensitivity", "randomness")) &&
    identical(mechanism$version, profile$release_mechanism) &&
    identical(mechanism$family, family) &&
    identical(
      .dsvert_joint_dp_client_json(mechanism$sensitivity),
      .dsvert_joint_dp_client_json(expected_sensitivity)) &&
    .dsvert_dp_has_exact_names(randomness, c("version", "lanes")) &&
    identical(randomness$version, "dsvert-randomness-plan-v1") &&
    identical(names(randomness$lanes), "final_noise") &&
    identical(.dsvert_joint_dp_client_json(lane),
              .dsvert_joint_dp_client_json(expected_lane))
  expected_request <- if (isTRUE(profile$gaussian)) {
    trusted$manifest$workload$mechanism_selection$
      gaussian_calibration_request
  } else {
    list(
      epsilon = privacy$epsilon, delta = privacy$delta,
      sensitivity_steps = physical$lattice$sensitivity_steps,
      total_coordinate_count = as.integer(layout$coordinate_count))
  }
  if (!isTRUE(privacy_valid) ||
      !identical(.dsvert_joint_dp_client_json(physical$request),
                 .dsvert_joint_dp_client_json(expected_request))) {
    stop("The synopsis privacy contract is invalid", call. = FALSE)
  }
  .dsvert_vector_plan_validate(
    physical$full_plan, physical$full_plan_sha256, profile,
    layout$coordinate_count,
    if (isTRUE(profile$gaussian)) physical$request$l2_sensitivity_steps else
      lattice$sensitivity_steps)
  if (!isTRUE(profile$gaussian) &&
      (!.dsvert_dp_synopsis_client_fraction_leq(
        physical$draw_law$epsilon_effective_upper_numerator,
        physical$draw_law$epsilon_effective_upper_denominator,
        physical$request$epsilon, "effective epsilon") ||
       !.dsvert_dp_synopsis_client_fraction_leq(
        physical$draw_law$implementation_delta_numerator,
        physical$draw_law$implementation_delta_denominator,
        physical$request$delta, "implementation delta"))) {
    stop("The synopsis Laplace calibration exceeds its request",
         call. = FALSE)
  }

  peers <- context$servers
  receipts <- compilation$receipts
  if (!is.list(receipts) || length(receipts) != length(peers)) {
    stop("Synopsis compilation requires one receipt per peer", call. = FALSE)
  }
  verified <- lapply(receipts, function(receipt) {
    expected <- c(
      "version", "peer_name", "peer_identity_pk", "manifest_sha256",
      "artifact_key", "source_claim_set_sha256", "full_plan_sha256",
      "signature")
    if (!.dsvert_dp_has_exact_names(receipt, expected) ||
        !identical(receipt$version,
                   .DSVERT_CLIENT_SYNOPSIS_COMPILE_VERSION) ||
        !receipt$peer_name %in% peers ||
        !identical(receipt$peer_identity_pk,
                   unname(context$pinset[[receipt$peer_name]])) ||
        !identical(receipt$manifest_sha256,
                   manifest_bundle$manifest_sha256) ||
        !identical(receipt$artifact_key, artifact$artifact_key) ||
        !identical(receipt$source_claim_set_sha256,
                   artifact$semantic$source_claim_set_sha256) ||
        !identical(receipt$full_plan_sha256,
                   physical$full_plan_sha256)) {
      stop("Invalid synopsis compile receipt", call. = FALSE)
    }
    unsigned <- receipt[setdiff(names(receipt), "signature")]
    .dsvert_dp_synopsis_client_signature(
      unsigned, receipt$signature, receipt$peer_identity_pk,
      "dsVert/stateless-catalog-synopsis/compile-receipt/v1|",
      "synopsis compile receipt")
    receipt
  })
  names(verified) <- vapply(verified, `[[`, character(1L), "peer_name")
  if (anyDuplicated(names(verified)) || !setequal(names(verified), peers)) {
    stop("The synopsis compile receipt coverage is invalid", call. = FALSE)
  }
  unsigned <- unname(lapply(verified[peers], function(receipt) {
    receipt[setdiff(names(receipt), "signature")]
  }))
  expected_set <- .dsvert_dp_synopsis_client_hash(
    "dsVert/stateless-catalog-synopsis/compile-receipt-set/v1|", list(
      version = .DSVERT_CLIENT_SYNOPSIS_COMPILE_SET_VERSION,
      receipts = unsigned))
  if (!identical(compilation$receipt_set_sha256, expected_set)) {
    stop("The synopsis compile receipt-set hash is invalid", call. = FALSE)
  }
  list(
    artifact = artifact, physical = physical, profile = profile,
    layout = layout, lattice = lattice, receipts = verified)
}

.dsvert_dp_synopsis_client_execution <- function(compiled) {
  artifact <- compiled$artifact
  physical <- compiled$physical
  profile <- compiled$profile
  lattice <- compiled$lattice
  dimension <- as.integer(compiled$layout$coordinate_count)
  execution_id <- digest::digest(charToRaw(paste0(
    "dsVert/stateless-catalog-synopsis/execution-id/v1|",
    artifact$artifact_key)), "sha256", serialize = FALSE)
  support_field <- if (isTRUE(profile$gaussian)) {
    "maximum_noise_magnitude_two_peers"
  } else "maximum_noise_magnitude"
  support <- openssl::bignum(physical$full_plan[[support_field]])
  if (!isTRUE(profile$gaussian) && !isTRUE(profile$exact_gc)) {
    support <- support * openssl::bignum(2)
  }
  maximum <- openssl::bignum(0)
  for (index in seq_len(dimension)) {
    scaled <- openssl::bignum(compiled$lattice$raw_upper_bounds[[index]]) *
      (openssl::bignum(2) ^ compiled$lattice$scale_shifts[[index]])
    if (scaled > maximum) maximum <- scaled
  }
  positive_limit <- (openssl::bignum(2) ^ 127L) - 1L
  if (support > positive_limit || maximum + support > positive_limit) {
    stop("The synopsis Ring128 headroom is insufficient", call. = FALSE)
  }
  ring_unsigned <- list(
    version = "dsvert-stateless-catalog-synopsis-ring128-certificate-v1",
    ring_bits = 128L, fractional_bits = 0L,
    maximum_scaled_source_coordinate = as.character(maximum),
    maximum_release_noise_magnitude = as.character(support),
    positive_limit = as.character(positive_limit), no_wrap_certified = TRUE)
  ring <- c(ring_unsigned, list(sha256 = .dsvert_dp_synopsis_client_hash(
    "dsVert/stateless-catalog-synopsis/ring128-certificate/v1|",
    ring_unsigned)))
  public_coordinates <- min(
    .DSVERT_CLIENT_SYNOPSIS_PUBLIC_CHUNK_COORDINATES, dimension)
  list(
    execution_id = execution_id,
    geometry = list(
      coordinate_count = dimension,
      public_chunk_coordinates = public_coordinates,
      public_chunk_count = as.integer(ceiling(dimension / public_coordinates))),
    ring = ring)
}

.dsvert_dp_synopsis_client_release_set <- function(
    release_receipts, compiled, execution, trusted) {
  authorities <- trusted$context$designated
  if (!is.list(release_receipts) || length(release_receipts) != 2L) {
    stop("Exactly two synopsis RELEASE receipts are required", call. = FALSE)
  }
  fields <- c(
    "version", "phase", "execution_id", "artifact_key",
    "contract_sha256", "attempt_sha256", "source_contract_sha256",
    "result_set_sha256", "local_authority", "public_chunk_count",
    "final_chunk_hashes", "final_vector_root", "output_lattice_bits",
    "output_lattice_scale", "mechanism", "epsilon", "delta",
    "implementation_delta_numerator", "implementation_delta_denominator",
    "delta_aggregation", "postprocessing", "all_public_chunks_durable",
    "intermediate_payload_exposed", "durable_replay",
    "capability_available", "signature")
  plan <- compiled$physical$full_plan
  delta_fields <- if (isTRUE(compiled$profile$exact_gc)) {
    c("implementation_delta_numerator", "implementation_delta_denominator")
  } else c("per_peer_implementation_delta_numerator",
           "per_peer_implementation_delta_denominator")
  decoded <- lapply(release_receipts, function(json) {
    receipt <- .dsvert_dp_synopsis_client_json(
      json, "synopsis RELEASE",
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
    peer <- if (is.list(receipt$local_authority)) {
      receipt$local_authority$peer_name
    } else NULL
    role <- match(peer, authorities)
    hash_list <- receipt$final_chunk_hashes
    hashes <- if (is.list(hash_list) && is.null(names(hash_list)) &&
                 all(vapply(hash_list, function(value) {
                   is.character(value) && length(value) == 1L &&
                     !is.na(value)
                 }, logical(1L)))) {
      unlist(hash_list, use.names = FALSE)
    } else character()
    valid <- .dsvert_dp_has_exact_names(receipt, fields) &&
      identical(receipt$version,
                .DSVERT_CLIENT_SYNOPSIS_RELEASE_VERSION) &&
      identical(receipt$phase, "synopsis_released") &&
      length(role) == 1L && !is.na(role) &&
      identical(
        .dsvert_joint_dp_client_json(receipt$local_authority),
        .dsvert_joint_dp_client_json(list(
          peer_name = peer,
          identity_pk = unname(trusted$context$pinset[[peer]]),
          role = c("primary_noise_authority",
                   "secondary_noise_authority")[[role]]))) &&
      identical(receipt$execution_id, execution$execution_id) &&
      identical(receipt$artifact_key, compiled$artifact$artifact_key) &&
      .dsvert_vector_hex(receipt$contract_sha256) &&
      .dsvert_vector_hex(receipt$attempt_sha256) &&
      .dsvert_vector_hex(receipt$source_contract_sha256) &&
      .dsvert_vector_hex(receipt$result_set_sha256) &&
      identical(as.numeric(receipt$public_chunk_count),
                as.numeric(execution$geometry$public_chunk_count)) &&
      is.character(hashes) && length(hashes) == receipt$public_chunk_count &&
      all(vapply(hashes, .dsvert_vector_hex, logical(1L))) &&
      identical(receipt$final_vector_root,
                .dsvert_vector_merkle_root(hashes)) &&
      identical(as.numeric(receipt$output_lattice_bits),
                as.numeric(compiled$lattice$output_lattice_bits)) &&
      identical(receipt$output_lattice_scale,
                .dsvert_dp_synopsis_client_lattice_scale(
                  compiled$lattice$output_lattice_bits)) &&
      identical(receipt$mechanism, compiled$profile$release_mechanism) &&
      identical(receipt$epsilon, compiled$physical$request$epsilon) &&
      identical(receipt$delta, compiled$physical$request$delta) &&
      identical(receipt$implementation_delta_numerator,
                plan[[delta_fields[[1L]]]]) &&
      identical(receipt$implementation_delta_denominator,
                plan[[delta_fields[[2L]]]]) &&
      identical(receipt$delta_aggregation,
                compiled$profile$delta_aggregation) &&
      identical(receipt$postprocessing,
                compiled$profile$postprocessing) &&
      identical(receipt$all_public_chunks_durable, TRUE) &&
      identical(receipt$intermediate_payload_exposed, FALSE) &&
      identical(receipt$durable_replay, TRUE) &&
      identical(receipt$capability_available, TRUE)
    if (!isTRUE(valid)) stop("Invalid or misbound synopsis RELEASE",
                             call. = FALSE)
    unsigned <- receipt[setdiff(names(receipt), "signature")]
    .dsvert_dp_synopsis_client_signature(
      unsigned, receipt$signature, receipt$local_authority$identity_pk,
      "dsVert/stateless-catalog-synopsis/release/v1|",
      "synopsis RELEASE")
    receipt
  })
  names(decoded) <- vapply(decoded, function(value) {
    value$local_authority$peer_name
  }, character(1L))
  if (anyDuplicated(names(decoded)) || !setequal(names(decoded), authorities)) {
    stop("The synopsis RELEASE authority coverage is invalid", call. = FALSE)
  }
  decoded <- decoded[authorities]
  common <- setdiff(fields, c("local_authority", "signature"))
  if (!identical(
      .dsvert_joint_dp_client_json(decoded[[1L]][common]),
      .dsvert_joint_dp_client_json(decoded[[2L]][common]))) {
    stop("The two synopsis RELEASE receipts disagree", call. = FALSE)
  }
  list(receipts = decoded, reference = decoded[[1L]],
       json = stats::setNames(lapply(decoded, .dsvert_joint_dp_client_json),
                              authorities))
}

.dsvert_dp_synopsis_client_replay <- function(
    replay_responses, releases, compiled, execution, trusted) {
  reference <- releases$reference
  hashes <- unlist(reference$final_chunk_hashes, use.names = FALSE)
  chunk_count <- as.integer(reference$public_chunk_count)
  expected_names <- as.character(seq.int(0L, chunk_count - 1L))
  if (!is.list(replay_responses) || is.null(names(replay_responses)) ||
      anyDuplicated(names(replay_responses)) ||
      !setequal(names(replay_responses), expected_names)) {
    stop("The bilateral synopsis REPLAY chunk set is invalid", call. = FALSE)
  }
  authorities <- trusted$context$designated
  scaled <- character(compiled$layout$coordinate_count)
  replay_values <- vector("list", chunk_count)
  replay_fields <- c(
    "version", "phase", "execution_id", "artifact_key",
    "contract_sha256", "attempt_sha256", "source_contract_sha256",
    "result_set_sha256", "final_vector_root", "public_chunk_index",
    "public_chunk_count", "chunk_sha256", "chunk", "merkle_proof",
    "durable_replay", "source_store_read", "sampler_invoked",
    "finalizer_invoked", "transport_read")
  public_fields <- c(
    "version", "artifact_key", "execution_id", "contract_sha256",
    "attempt_sha256", "result_set_sha256", "public_chunk_index",
    "public_chunk_count", "coordinate_offset", "coordinate_count",
    "output_lattice_bits", "output_lattice_scale", "scaled_values",
    "value_encoding", "postprocessing", "source_values_exposed",
    "preclamp_values_exposed")
  client_limit <- (openssl::bignum(2) ^ 53L) - 1L
  for (index in seq.int(0L, chunk_count - 1L)) {
    bilateral <- replay_responses[[as.character(index)]]
    if (!is.list(bilateral) || is.null(names(bilateral)) ||
        anyDuplicated(names(bilateral)) ||
        !setequal(names(bilateral), authorities) ||
        !all(vapply(bilateral, .dsvert_dp_is_string, logical(1L))) ||
        length(unique(unname(unlist(bilateral, use.names = FALSE)))) != 1L) {
      stop("The bilateral synopsis REPLAY responses are different",
           call. = FALSE)
    }
    replay <- .dsvert_dp_synopsis_client_json(
      bilateral[[authorities[[1L]]]], "synopsis REPLAY",
      .DSVERT_CLIENT_SYNOPSIS_RECEIPT_MAX_OBJECT_BYTES)
    public <- replay$chunk
    offset <- index * execution$geometry$public_chunk_coordinates
    count <- min(
      execution$geometry$public_chunk_coordinates,
      compiled$layout$coordinate_count - offset)
    values <- if (is.list(public$scaled_values) &&
                  is.null(names(public$scaled_values)) &&
                  all(vapply(public$scaled_values, function(value) {
                    is.character(value) && length(value) == 1L &&
                      !is.na(value)
                  }, logical(1L)))) {
      unlist(public$scaled_values, use.names = FALSE)
    } else character()
    valid <- .dsvert_dp_has_exact_names(replay, replay_fields) &&
      .dsvert_dp_has_exact_names(public, public_fields) &&
      identical(replay$version, .DSVERT_CLIENT_SYNOPSIS_REPLAY_VERSION) &&
      identical(replay$phase, "synopsis_public_chunk_replayed") &&
      identical(replay$execution_id, execution$execution_id) &&
      identical(replay$artifact_key, compiled$artifact$artifact_key) &&
      identical(replay$contract_sha256, reference$contract_sha256) &&
      identical(replay$attempt_sha256, reference$attempt_sha256) &&
      identical(replay$source_contract_sha256,
                reference$source_contract_sha256) &&
      identical(replay$result_set_sha256, reference$result_set_sha256) &&
      identical(replay$final_vector_root, reference$final_vector_root) &&
      identical(as.numeric(replay$public_chunk_index), as.numeric(index)) &&
      identical(as.numeric(replay$public_chunk_count),
                as.numeric(chunk_count)) &&
      identical(replay$chunk_sha256, hashes[[index + 1L]]) &&
      identical(replay$chunk_sha256,
        .dsvert_dp_synopsis_client_hash(
          "dsVert/stateless-catalog-synopsis/public-chunk/v1|", public)) &&
      identical(.dsvert_joint_dp_client_canonical(replay$merkle_proof),
                .dsvert_joint_dp_client_canonical(
                  .dsvert_vector_merkle_proof(hashes, index))) &&
      identical(replay$durable_replay, TRUE) &&
      identical(replay$source_store_read, FALSE) &&
      identical(replay$sampler_invoked, FALSE) &&
      identical(replay$finalizer_invoked, FALSE) &&
      identical(replay$transport_read, FALSE) &&
      identical(public$version, .DSVERT_CLIENT_SYNOPSIS_PUBLIC_VERSION) &&
      identical(public$artifact_key, compiled$artifact$artifact_key) &&
      identical(public$execution_id, execution$execution_id) &&
      identical(public$contract_sha256, reference$contract_sha256) &&
      identical(public$attempt_sha256, reference$attempt_sha256) &&
      identical(public$result_set_sha256, reference$result_set_sha256) &&
      identical(as.numeric(public$public_chunk_index), as.numeric(index)) &&
      identical(as.numeric(public$public_chunk_count),
                as.numeric(chunk_count)) &&
      identical(as.numeric(public$coordinate_offset), as.numeric(offset)) &&
      identical(as.numeric(public$coordinate_count), as.numeric(count)) &&
      identical(as.numeric(public$output_lattice_bits),
                as.numeric(compiled$lattice$output_lattice_bits)) &&
      identical(as.numeric(public$output_lattice_scale),
                as.numeric(compiled$lattice$output_lattice_scale)) &&
      is.character(values) && length(values) == count && !anyNA(values) &&
      all(vapply(values, .dsvert_vector_integer_text, logical(1L))) &&
      identical(public$value_encoding,
                "nonnegative-decimal-integer-common-lattice-v1") &&
      identical(public$postprocessing, compiled$profile$postprocessing) &&
      identical(public$source_values_exposed, FALSE) &&
      identical(public$preclamp_values_exposed, FALSE)
    if (!isTRUE(valid)) stop("Invalid or misbound synopsis REPLAY chunk",
                             call. = FALSE)
    for (local in seq_along(values)) {
      coordinate <- offset + local
      integer <- tryCatch(openssl::bignum(values[[local]]),
                          error = function(error) NULL)
      upper <- openssl::bignum(
        compiled$lattice$raw_upper_bounds[[coordinate]]) *
        (openssl::bignum(2) ^
           compiled$lattice$scale_shifts[[coordinate]])
      if (is.null(integer) || integer > upper) {
        stop("A synopsis coordinate exceeds its signed bound", call. = FALSE)
      }
      if (integer > client_limit) {
        stop("A synopsis coordinate exceeds the exact client integer domain",
             call. = FALSE)
      }
    }
    scaled[seq.int(offset + 1L, length.out = count)] <- values
    replay_values[[index + 1L]] <- replay
  }
  list(scaled = scaled, replay = replay_values)
}

.dsvert_dp_synopsis_public_vector_v1 <- function(
    release_receipts, replay_responses, manifest_bundle, status, artifact) {
  trusted <- .dsvert_dp_synopsis_client_bundle(manifest_bundle, status)
  compiled <- .dsvert_dp_synopsis_client_compile(
    artifact, trusted, manifest_bundle)
  execution <- .dsvert_dp_synopsis_client_execution(compiled)
  releases <- .dsvert_dp_synopsis_client_release_set(
    release_receipts, compiled, execution, trusted)
  replay <- .dsvert_dp_synopsis_client_replay(
    replay_responses, releases, compiled, execution, trusted)
  reference <- releases$reference
  result <- list(
    version = "dsvert-stateless-synopsis-public-vector-client-v1",
    artifact_key = compiled$artifact$artifact_key,
    execution_id = execution$execution_id,
    manifest_sha256 = manifest_bundle$manifest_sha256,
    final_vector_root = reference$final_vector_root,
    result_set_sha256 = reference$result_set_sha256,
    contract_sha256 = reference$contract_sha256,
    attempt_sha256 = reference$attempt_sha256,
    source_contract_sha256 = reference$source_contract_sha256,
    coordinate_order_sha256 = compiled$layout$sha256,
    coordinate_count = as.integer(compiled$layout$coordinate_count),
    output_lattice_bits = as.integer(compiled$lattice$output_lattice_bits),
    output_lattice_scale = as.numeric(compiled$lattice$output_lattice_scale),
    scaled_values = replay$scaled,
    values = .dsvert_vector_scaled_to_double(
      replay$scaled, compiled$lattice$output_lattice_scale),
    backend = compiled$profile$backend,
    sampler = compiled$profile$sampler,
    mechanism_plan = compiled$physical$full_plan,
    plan_sha256 = compiled$physical$full_plan_sha256,
    postprocessing = compiled$profile$postprocessing,
    mechanism = reference$mechanism,
    epsilon = as.numeric(reference$epsilon),
    delta = as.numeric(reference$delta),
    implementation_delta = paste0(
      reference$implementation_delta_numerator, "/",
      reference$implementation_delta_denominator),
    delta_aggregation = reference$delta_aggregation,
    sticky_replay = TRUE, source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    manifest = trusted$manifest,
    signed_provenance = list(
      version = "dsvert-stateless-synopsis-public-provenance-v1",
      ordered_peer_pinset = as.list(trusted$context$pinset),
      peer_pinset_sha256 = trusted$status[[trusted$context$servers[[1L]]]]$
        policy$peer_pinset_sha256,
      designated_noise_peers = as.list(trusted$context$designated),
      artifact_key = compiled$artifact$artifact_key,
      execution_id = execution$execution_id,
      contract_sha256 = reference$contract_sha256,
      attempt_sha256 = reference$attempt_sha256,
      source_contract_sha256 = reference$source_contract_sha256,
      result_set_sha256 = reference$result_set_sha256,
      final_vector_root = reference$final_vector_root,
      full_plan_sha256 = compiled$physical$full_plan_sha256,
      execution_hash_binding = paste0(
        "both_pinned_designated_RELEASE_signatures;",
        "at_most_one_noise_authority_corrupt_v1"),
      execution_hashes_client_reconstructed = FALSE,
      compile_receipts = compiled$receipts,
      release_receipts = releases$receipts,
      replay_responses = replay$replay,
      protected_shares_included = FALSE,
      preclamp_values_included = FALSE,
      source_values_included = FALSE,
      intermediate_payload_exposed = FALSE,
      durable_replay = TRUE))
  class(result) <- c(
    "dsvert_synopsis_public_vector", "dsvert_joint_dp_vector", "list")
  result
}
