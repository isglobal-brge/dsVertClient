# Client orchestration for the single sticky biomedical capsule vector.  The
# relay validates only public signed receipts and final DP chunks. Exact source
# shares, sticky seeds, local noised shares and pre-clamp values never enter R.

.DSVERT_CLIENT_VECTOR_PREPARE_VERSION <-
  "dsvert-joint-dp-vector-prepare-v6"
.DSVERT_CLIENT_VECTOR_START_VERSION <-
  "dsvert-joint-dp-vector-chunk-start-v5"
.DSVERT_CLIENT_VECTOR_RESULT_VERSION <-
  "dsvert-joint-dp-vector-local-result-v5"
.DSVERT_CLIENT_VECTOR_RELEASE_VERSION <-
  "dsvert-joint-dp-vector-release-root-v5"
.DSVERT_CLIENT_VECTOR_CHUNK_VERSION <-
  "dsvert-joint-dp-vector-final-chunk-v4"
.DSVERT_CLIENT_VECTOR_ACK_VERSION <-
  "dsvert-joint-dp-vector-finalization-ack-v5"
.DSVERT_CLIENT_VECTOR_RELEASE_INSTANCE_VERSION <-
  "dsvert-joint-dp-vector-release-instance-v2"
.DSVERT_CLIENT_VECTOR_BACKEND <-
  "independent_full_global_draw_convolution_ring128_v3"
.DSVERT_CLIENT_VECTOR_SAMPLER <-
  "hkdf-sha256-chacha20-independent-full-draw-binary-geometric-tv-v3"
.DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM <-
  "two-independent-complete-vector-discrete-laplace-draws-v3"
.DSVERT_CLIENT_VECTOR_EXACT_BACKEND <-
  "exact_gc_one_joint_discrete_laplace_draw_ring128_v3"
.DSVERT_CLIENT_VECTOR_EXACT_SAMPLER <-
  "hkdf-sha256-chacha20-xor-binary-geometric-tv-v3"
.DSVERT_CLIENT_VECTOR_EXACT_RELEASE_MECHANISM <-
  "one-joint-complete-vector-discrete-laplace-draw-v3"
.DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM <-
  "dyadic_discrete_gaussian_truncated_tv_bounded"
.DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION <-
  "dsvert-joint-dp-vector-dyadic-discrete-gaussian-tv-bounded-plan-v2"
.DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND <-
  paste0("independent_full_global_dyadic_discrete_gaussian_",
         "tv_bounded_ring128_v2")
.DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER <-
  paste0("cks-target-outward-rational-dyadic-cdf-hkdf-sha256-",
         "chacha20-coordinate-domain-v2")
.DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM <-
  paste0("two-independent-complete-vector-dyadic-discrete-gaussian-",
         "tv-bounded-draws-v2")
.DSVERT_CLIENT_VECTOR_TYPED_CAPABILITY <-
  "blob.joint-dp.vector-final-share.v3"
.DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES <- 8192L
.DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES <- 2L * 1024L^2
.DSVERT_CLIENT_VECTOR_IMPLEMENTATION_DELTA_MAX_BYTES <- 1024L
.DSVERT_CLIENT_JOINT_DP_RECEIPT_DOMAIN <-
  "dsVert/joint-dp/signed-receipt/v1|"
.DSVERT_CLIENT_VECTOR_RETRY_CURRENT_INSTANCE_TOKEN <-
  "dsvert_retry_current_release_instance"

.dsvert_vector_profile <- function(capsule_mechanism,
                                   mechanism_selection = NULL,
                                   backend = NULL) {
  mechanism <- if (is.list(capsule_mechanism)) {
    capsule_mechanism$mechanism
  } else capsule_mechanism
  if (!.dsvert_vector_string(mechanism, maximum_bytes = 128L)) {
    stop("Invalid biomedical vector mechanism", call. = FALSE)
  }
  if (identical(mechanism, .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM)) {
    if (!is.null(backend) &&
        !identical(backend, .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND)) {
      stop("The signed Gaussian vector backend is invalid", call. = FALSE)
    }
    selection <- if (is.list(mechanism_selection)) {
      mechanism_selection
    } else if (is.list(capsule_mechanism)) {
      capsule_mechanism$certificate
    } else NULL
    return(list(
      gaussian = TRUE, exact_gc = FALSE, selection_bound = FALSE,
      mechanism = mechanism,
      plan_version = .DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION,
      backend = .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND,
      sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER,
      release_mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM,
      complete_epsilon_per_peer = TRUE,
      delta_aggregation = "max_per_peer_not_sum",
      postprocessing =
        "signed-Ring128-decode-then-fixed-public-coordinate-clamp-v1",
      manifest_selection = selection))
  }
  if (!identical(mechanism, "discrete-laplace")) {
    stop("Unsupported biomedical vector mechanism", call. = FALSE)
  }
  if (is.null(backend) ||
      identical(backend, .DSVERT_CLIENT_VECTOR_EXACT_BACKEND)) {
    list(
      gaussian = FALSE, exact_gc = TRUE, selection_bound = TRUE,
      mechanism = mechanism,
      plan_version = "dsvert-joint-dp-vector-laplace-plan-v3",
      backend = .DSVERT_CLIENT_VECTOR_EXACT_BACKEND,
      sampler = .DSVERT_CLIENT_VECTOR_EXACT_SAMPLER,
      release_mechanism = .DSVERT_CLIENT_VECTOR_EXACT_RELEASE_MECHANISM,
      complete_epsilon_per_peer = FALSE,
      delta_aggregation = "single_joint_draw",
      postprocessing = paste0(
        "fixed-public-coordinate-clamp-inside-exact-GC-before-",
        "selective-sharing-v1"),
      manifest_selection = if (is.list(capsule_mechanism)) {
        capsule_mechanism$certificate
      } else NULL)
  } else if (identical(backend, .DSVERT_CLIENT_VECTOR_BACKEND)) {
    list(
      gaussian = FALSE, exact_gc = FALSE, selection_bound = TRUE,
      mechanism = mechanism,
      plan_version = paste0(
        "dsvert-joint-dp-vector-independent-full-draw-",
        "convolution-plan-v3"),
      backend = .DSVERT_CLIENT_VECTOR_BACKEND,
      sampler = .DSVERT_CLIENT_VECTOR_SAMPLER,
      release_mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
      complete_epsilon_per_peer = TRUE,
      delta_aggregation = "max_per_peer_not_sum",
      postprocessing =
        "signed-Ring128-decode-then-fixed-public-coordinate-clamp-v1",
      manifest_selection = if (is.list(capsule_mechanism)) {
        capsule_mechanism$certificate
      } else NULL)
  } else {
    stop("The signed Laplace vector backend is invalid", call. = FALSE)
  }
}

.dsvert_vector_release_instance <- function(context, capsule_id) {
  peers <- context$designated
  if (!is.character(peers) || length(peers) != 2L || anyNA(peers) ||
      anyDuplicated(peers) ||
      !identical(peers, sort(peers, method = "radix")) ||
      !.dsvert_vector_hex(capsule_id)) {
    stop("Invalid designated-peer release-instance context", call. = FALSE)
  }
  roots <- stats::setNames(lapply(peers, function(peer) {
    root <- context$status[[peer]]$noise_root
    epoch <- if (is.list(root)) root$privacy_epoch else NULL
    key_id <- if (is.list(root)) root$key_id else NULL
    provider_id <- if (is.list(root)) root$provider_id else NULL
    release_domain <- context$status[[peer]]$release_domain
    release_domain_generation <- if (is.list(release_domain)) {
      release_domain$generation
    } else NULL
    release_domain_id <- if (is.list(release_domain)) {
      release_domain$domain_id
    } else NULL
    valid <- .dsvert_vector_whole(epoch, 1, 2^53 - 1) &&
      .dsvert_vector_string(key_id, maximum_bytes = 256L) &&
      .dsvert_vector_string(provider_id, maximum_bytes = 256L) &&
      .dsvert_vector_whole(
        release_domain_generation, 1, 2^53 - 1) &&
      .dsvert_vector_string(
        release_domain_id, "^rd_[0-9a-f]{64}$", 67L)
    if (!isTRUE(valid)) {
      stop("A designated peer has no valid signed noise-root status",
           call. = FALSE)
    }
    list(
      privacy_epoch = as.numeric(epoch), noise_key_id = key_id,
      provider_id = provider_id,
      release_domain_generation = as.numeric(release_domain_generation),
      release_domain_id = release_domain_id)
  }), peers)
  key_ids <- vapply(roots, `[[`, character(1L), "noise_key_id")
  if (anyDuplicated(key_ids)) {
    stop(structure(list(
      message = paste(
        "The two designated peers advertise the same noise-root key ID;",
        "ask the administrator of the duplicated peer to delete and",
        "regenerate only its privacy/noise_root, then retry."),
      call = NULL, code = "noise_root_not_independent"),
      class = c("dsvert_noise_root_not_independent", "error", "condition")))
  }
  value <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_VECTOR_RELEASE_INSTANCE_VERSION,
    capsule_id = capsule_id,
    peer_noise_roots = roots))
  list(
    value = value,
    json = .dsvert_joint_dp_client_json(value),
    id = .dsvert_vector_hash(value))
}

.dsvert_vector_exact_names <- function(value, expected) {
  is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && setequal(names(value), expected)
}

.dsvert_vector_string <- function(value, pattern = NULL,
                                  maximum_bytes = 1024L) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    nzchar(value) && nchar(value, type = "bytes") <= maximum_bytes &&
    (is.null(pattern) || grepl(pattern, value, perl = TRUE))
}

.dsvert_vector_hex <- function(value) {
  .dsvert_vector_string(value, "^[0-9a-f]{64}$", 64L)
}

.dsvert_vector_whole <- function(value, minimum = 0, maximum = 2^31 - 1) {
  is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value == floor(value) && value >= minimum &&
    value <= maximum
}

.dsvert_vector_integer_text <- function(value, positive = FALSE,
                                        maximum_bytes = 128L) {
  .dsvert_vector_string(
    value, if (isTRUE(positive)) "^[1-9][0-9]*$" else
      "^(0|[1-9][0-9]*)$", maximum_bytes)
}

.dsvert_vector_decimal <- function(value, minimum, maximum,
                                   open_minimum = FALSE) {
  if (!.dsvert_vector_string(value, maximum_bytes = 128L)) return(NULL)
  number <- suppressWarnings(as.numeric(value))
  if (length(number) != 1L || is.na(number) || !is.finite(number) ||
      number > maximum || if (isTRUE(open_minimum)) number <= minimum else
        number < minimum) return(NULL)
  number
}

.dsvert_vector_hash <- function(value) {
  digest::digest(.dsvert_joint_dp_client_json(value), algo = "sha256",
                 serialize = FALSE)
}

.dsvert_vector_gaussian_plan_fields <- function() {
  c(
    "version", "mechanism", "sampler", "reference",
    "total_coordinate_count", "maximum_chunk_coordinates",
    "request_binding_sha256",
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
    "vector_tail_tv_upper_numerator",
    "vector_tail_tv_upper_denominator",
    "per_peer_implementation_delta_numerator",
    "per_peer_implementation_delta_denominator",
    "vector_sampler_tv_upper_numerator",
    "vector_sampler_tv_upper_denominator",
    "vector_total_tv_upper_numerator",
    "vector_total_tv_upper_denominator",
    "sampler_candidate_count", "sampler_random_bits_per_coordinate",
    "sampler_random_bytes_per_coordinate", "sampler_table_precision_bits",
    "sampler_magnitude_count", "sampler_search_steps",
    "sampler_full_scan_steps", "sampler_cdf_table_bytes",
    "simultaneous_95_abs", "accuracy_accounting", "accounting",
    "privacy_theorem", "independent_noise_peer_count",
    "complete_epsilon_per_peer", "epsilon_divided_by_peer_count",
    "release_delta_aggregation", "nominal_variance_multiplier",
    "nominal_standard_deviation_factor", "at_least_one_honest_noise_peer",
    "maximum_colluding_noise_peers", "adversary_view",
    "adversary_view_privacy_argument", "source_share_hiding_precondition",
    "exact_rational_sampler", "finite_support_transfer_charged",
    "fixed_work_sampler", "sampler_branches_on_protected_values",
    "sampler_branches_on_private_randomness", "host_constant_time_claim",
    "transcript_dp_claim", "logical_transcript_fixed_shape",
    "physical_timing_dp_claim", "observable_worker_shape",
    "capability_available", "unavailable_reason")
}

.dsvert_vector_plan_validate <- function(
    plan, plan_sha256, profile, coordinate_count, sensitivity_steps) {
  if (isTRUE(profile$exact_gc)) {
    required <- c(
      "version", "sampler", "stop_bits", "stop_numerator",
      "uniform_bits", "binary_geometric_bits", "bernoulli_thresholds",
      "sensitivity_steps", "total_coordinate_count",
      "epsilon_effective_upper_numerator",
      "epsilon_effective_upper_denominator",
      "one_geometric_tv_numerator", "one_geometric_tv_denominator",
      "tail_upper_numerator", "tail_upper_denominator",
      "rounding_upper_numerator", "rounding_upper_denominator",
      "implementation_delta_numerator",
      "implementation_delta_denominator", "implementation_delta_bound",
      "maximum_noise_magnitude", "maximum_chunk_coordinates",
      "private_stream_bytes_per_coordinate", "accounting",
      "capability_available")
    integer_fields <- c(
      "stop_numerator", "sensitivity_steps",
      "epsilon_effective_upper_numerator",
      "epsilon_effective_upper_denominator",
      "one_geometric_tv_numerator", "one_geometric_tv_denominator",
      "tail_upper_numerator", "tail_upper_denominator",
      "rounding_upper_numerator", "rounding_upper_denominator",
      "implementation_delta_numerator",
      "implementation_delta_denominator", "maximum_noise_magnitude")
    thresholds <- if (is.list(plan)) {
      unlist(plan$bernoulli_thresholds, use.names = FALSE)
    } else NULL
    valid <- is.list(plan) && .dsvert_vector_exact_names(plan, required) &&
      .dsvert_vector_hex(plan_sha256) &&
      identical(.dsvert_vector_hash(plan), plan_sha256) &&
      identical(plan$version, profile$plan_version) &&
      identical(plan$sampler, profile$sampler) &&
      identical(as.numeric(plan$total_coordinate_count),
                as.numeric(coordinate_count)) &&
      identical(plan$sensitivity_steps, sensitivity_steps) &&
      identical(as.numeric(plan$stop_bits), 128) &&
      as.numeric(plan$uniform_bits) %in% c(128, 256) &&
      .dsvert_vector_whole(plan$binary_geometric_bits, 1, 63) &&
      is.character(thresholds) &&
      length(thresholds) == as.numeric(plan$binary_geometric_bits) &&
      !anyNA(thresholds) &&
      all(grepl("^(0|[1-9][0-9]*)$", thresholds)) &&
      all(vapply(plan[integer_fields], function(value) {
        .dsvert_vector_integer_text(
          value, maximum_bytes = .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES)
      }, logical(1L))) &&
      .dsvert_vector_whole(
        plan$maximum_chunk_coordinates, 1,
        .DSVERT_CLIENT_JOINT_DP_VECTOR_MAX_CHUNK) &&
      .dsvert_vector_whole(
        plan$private_stream_bytes_per_coordinate, 1, 2^31 - 1) &&
      .dsvert_vector_string(plan$implementation_delta_bound,
                            maximum_bytes =
                              .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES) &&
      .dsvert_vector_string(plan$accounting, maximum_bytes = 2048L) &&
      startsWith(plan$accounting, "global iid discrete Laplace") &&
      identical(plan$capability_available, TRUE)
    if (!isTRUE(valid)) {
      stop("The signed one-draw exact-GC vector plan is invalid",
           call. = FALSE)
    }
    return(invisible(.dsvert_joint_dp_client_canonical(plan)))
  }
  if (!is.list(plan) || is.null(names(plan)) || anyNA(names(plan)) ||
      anyDuplicated(names(plan)) || !.dsvert_vector_hex(plan_sha256) ||
      !identical(.dsvert_vector_hash(plan), plan_sha256) ||
      !is.list(profile) ||
      !identical(plan$version, profile$plan_version) ||
      !identical(plan$sampler, profile$sampler) ||
      !identical(as.numeric(plan$total_coordinate_count),
                 as.numeric(coordinate_count)) ||
      !identical(as.numeric(plan$maximum_chunk_coordinates),
                 as.numeric(min(.DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES,
                                coordinate_count))) ||
      !identical(as.numeric(plan$independent_noise_peer_count), 2) ||
      !identical(plan$complete_epsilon_per_peer, TRUE) ||
      !identical(plan$epsilon_divided_by_peer_count, FALSE) ||
      !identical(plan$capability_available, TRUE) ||
      !.dsvert_vector_integer_text(
        plan$per_peer_implementation_delta_numerator,
        maximum_bytes = .DSVERT_CLIENT_VECTOR_IMPLEMENTATION_DELTA_MAX_BYTES) ||
      !.dsvert_vector_integer_text(
        plan$per_peer_implementation_delta_denominator, TRUE,
        maximum_bytes = .DSVERT_CLIENT_VECTOR_IMPLEMENTATION_DELTA_MAX_BYTES)) {
    stop("The signed vector mechanism plan is invalid", call. = FALSE)
  }
  if (isTRUE(profile$gaussian)) {
    selection <- profile$manifest_selection
    request <- if (is.list(selection)) {
      selection$gaussian_calibration_request
    } else NULL
    expected_plan <- if (is.list(selection)) selection$gaussian_plan else NULL
    required_request <- c(
      "epsilon", "delta", "l2_sensitivity_steps",
      "total_coordinate_count")
    integer_scalar <- function(value, minimum = 0, maximum = Inf) {
      is.numeric(value) && length(value) == 1L && !is.na(value) &&
        is.finite(value) && value == floor(value) &&
        value >= minimum && value <= maximum
    }
    positive_fraction <- function(prefix) {
      .dsvert_vector_integer_text(
        plan[[paste0(prefix, "_numerator")]], positive = TRUE) &&
        .dsvert_vector_integer_text(
          plan[[paste0(prefix, "_denominator")]], positive = TRUE)
    }
    valid <- .dsvert_vector_exact_names(
        plan, .dsvert_vector_gaussian_plan_fields()) &&
      identical(plan$mechanism, profile$mechanism) &&
      identical(plan$release_delta_aggregation, "max_per_peer_not_sum") &&
      identical(plan$exact_rational_sampler, FALSE) &&
      identical(plan$finite_support_transfer_charged, TRUE) &&
      identical(plan$fixed_work_sampler, TRUE) &&
      identical(plan$sampler_branches_on_protected_values, FALSE) &&
      identical(plan$sampler_branches_on_private_randomness, FALSE) &&
      identical(plan$host_constant_time_claim, FALSE) &&
      identical(plan$transcript_dp_claim, TRUE) &&
      identical(plan$logical_transcript_fixed_shape, TRUE) &&
      identical(plan$physical_timing_dp_claim, FALSE) &&
      identical(as.numeric(plan$nominal_variance_multiplier), 2) &&
      identical(plan$nominal_standard_deviation_factor,
                "sqrt(2)_relative_to_one_full_draw") &&
      identical(plan$at_least_one_honest_noise_peer, TRUE) &&
      identical(as.numeric(plan$maximum_colluding_noise_peers), 1) &&
      identical(plan$adversary_view,
        paste0("analyst_plus_at_most_one_designated_noise_peer_including_",
               "its_seed_draw_source_share_and_protocol_transcript")) &&
      identical(plan$adversary_view_privacy_argument,
        paste0("conditioned_on_a_simulatable_own_share_and_fixed_corrupt_",
               "peer_view_the_other_independent_complete_epsilon_full_",
               "sensitivity_draw_is_an_epsilon_delta_DP_mechanism; own_",
               "draw_translation_second_draw_signed_decode_and_public_",
               "clamp_are_post_processing; release_delta_is_max_of_the_",
               "two_symmetric_conditional_guarantees")) &&
      identical(plan$source_share_hiding_precondition,
        paste0("each_single_pre_noise_aggregate_share_is_computationally_",
               "simulatable_without_the_protected_query_under_",
               "authenticated_semi_honest_fanin")) &&
      identical(as.numeric(plan$sampler_candidate_count), 1) &&
      integer_scalar(plan$sampler_random_bits_per_coordinate, 128, 16384) &&
      plan$sampler_random_bits_per_coordinate %% 8 == 0 &&
      identical(as.numeric(plan$sampler_random_bytes_per_coordinate),
                as.numeric(plan$sampler_random_bits_per_coordinate / 8 + 1)) &&
      integer_scalar(plan$sampler_table_precision_bits,
                     plan$sampler_random_bits_per_coordinate) &&
      integer_scalar(plan$sampler_magnitude_count, 1, 2^20 + 1) &&
      integer_scalar(plan$sampler_search_steps, 1, 21) &&
      integer_scalar(plan$sampler_full_scan_steps, 1, 2^20 + 1) &&
      identical(as.numeric(plan$sampler_full_scan_steps),
                as.numeric(plan$sampler_magnitude_count)) &&
      integer_scalar(plan$sampler_cdf_table_bytes, 1, 2^53 - 1) &&
      identical(as.numeric(plan$sampler_cdf_table_bytes),
                as.numeric(plan$sampler_magnitude_count) *
                  as.numeric(plan$sampler_random_bytes_per_coordinate)) &&
      positive_fraction("vector_tail_tv_upper") &&
      positive_fraction("vector_sampler_tv_upper") &&
      positive_fraction("vector_total_tv_upper") &&
      .dsvert_vector_string(plan$observable_worker_shape,
                            maximum_bytes = 1024L) &&
      .dsvert_vector_hex(plan$request_binding_sha256) &&
      .dsvert_vector_exact_names(request, required_request) &&
      identical(request$l2_sensitivity_steps, sensitivity_steps) &&
      identical(as.numeric(request$total_coordinate_count),
                as.numeric(coordinate_count)) &&
      is.list(expected_plan) &&
      identical(.dsvert_joint_dp_client_json(plan),
                .dsvert_joint_dp_client_json(expected_plan)) &&
      identical(selection$gaussian_plan_sha256, plan_sha256)
  } else {
    valid <- identical(plan$sensitivity_steps, sensitivity_steps) &&
      .dsvert_vector_integer_text(sensitivity_steps, TRUE) &&
      identical(plan$release_implementation_delta_aggregation,
                "max_per_peer_not_sum")
  }
  if (!isTRUE(valid)) {
    stop("The signed vector mechanism plan is misbound", call. = FALSE)
  }
  invisible(.dsvert_joint_dp_client_canonical(plan))
}

.dsvert_vector_implementation_delta <- function(contract) {
  exact_gc <- isTRUE(contract$profile$exact_gc)
  maximum_bytes <- if (exact_gc) {
    .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES
  } else .DSVERT_CLIENT_VECTOR_IMPLEMENTATION_DELTA_MAX_BYTES
  fields <- if (exact_gc) {
    c("implementation_delta_numerator", "implementation_delta_denominator")
  } else {
    c("per_peer_implementation_delta_numerator",
      "per_peer_implementation_delta_denominator")
  }
  values <- contract$mechanism_plan[fields]
  if (length(values) != 2L ||
      !.dsvert_vector_integer_text(
        values[[1L]], maximum_bytes = maximum_bytes) ||
      !.dsvert_vector_integer_text(
        values[[2L]], positive = TRUE, maximum_bytes = maximum_bytes)) {
    stop("The signed vector implementation-delta certificate is invalid",
         call. = FALSE)
  }
  paste0(values[[1L]], "/", values[[2L]])
}

.dsvert_vector_hex_raw <- function(value) {
  if (!.dsvert_vector_hex(value)) {
    stop("Invalid vector commitment", call. = FALSE)
  }
  as.raw(strtoi(substring(value, seq.int(1L, 63L, 2L),
                          seq.int(2L, 64L, 2L)), 16L))
}

.dsvert_vector_merkle_leaf <- function(value) {
  digest::digest(c(
    charToRaw("dsVert/joint-dp/vector-merkle-leaf/v3"), as.raw(0L),
    .dsvert_vector_hex_raw(value)), algo = "sha256", serialize = FALSE)
}

.dsvert_vector_merkle_parent <- function(left, right) {
  digest::digest(c(
    charToRaw("dsVert/joint-dp/vector-merkle-parent/v3"), as.raw(0L),
    .dsvert_vector_hex_raw(left), .dsvert_vector_hex_raw(right)),
    algo = "sha256", serialize = FALSE)
}

.dsvert_vector_merkle_root <- function(hashes) {
  if (!is.character(hashes) || !length(hashes) || anyNA(hashes) ||
      any(!vapply(hashes, .dsvert_vector_hex, logical(1L)))) {
    stop("Invalid vector Merkle leaf set", call. = FALSE)
  }
  nodes <- vapply(hashes, .dsvert_vector_merkle_leaf, character(1L))
  while (length(nodes) > 1L) {
    if (length(nodes) %% 2L) nodes <- c(nodes, tail(nodes, 1L))
    nodes <- vapply(seq.int(1L, length(nodes), by = 2L), function(index) {
      .dsvert_vector_merkle_parent(nodes[[index]], nodes[[index + 1L]])
    }, character(1L))
  }
  unname(nodes[[1L]])
}

.dsvert_vector_merkle_proof <- function(hashes, index) {
  if (!.dsvert_vector_whole(index, 0, length(hashes) - 1L)) {
    stop("Invalid vector Merkle index", call. = FALSE)
  }
  nodes <- vapply(hashes, .dsvert_vector_merkle_leaf, character(1L))
  position <- as.integer(index + 1L)
  proof <- list()
  while (length(nodes) > 1L) {
    if (length(nodes) %% 2L) nodes <- c(nodes, tail(nodes, 1L))
    sibling <- if (position %% 2L) position + 1L else position - 1L
    proof[[length(proof) + 1L]] <- list(
      side = if (position %% 2L) "right" else "left",
      sha256 = nodes[[sibling]])
    nodes <- vapply(seq.int(1L, length(nodes), by = 2L), function(node) {
      .dsvert_vector_merkle_parent(nodes[[node]], nodes[[node + 1L]])
    }, character(1L))
    position <- as.integer((position + 1L) %/% 2L)
  }
  proof
}

.dsvert_vector_receipt_fields <- function(version, profile) {
  common <- c(
    "version", "phase", "capsule_id", "release_instance_id",
    "release_instance", "release_contract_hash",
    "transcript_hash", "peer_name", "peer_identity_pk",
    "coordinate_count", "chunk_count", "backend", "sampler",
    "history_gate", "request_limit", "operation_limit", "signature")
  extras <- switch(version,
    "dsvert-joint-dp-vector-prepare-v6" = c(
      "source_contract_hash", "coordinate_order_sha256",
      "lattice_transform_sha256", "mechanism_plan", "plan_sha256", "epsilon",
      "allocated_delta", "sensitivity_steps", "commitment_context",
      "seed_commitment", "complete_epsilon_per_peer",
      "delta_aggregation", "capability_available",
      if (isTRUE(profile$selection_bound)) c(
        "backend_assessment", "backend_selection",
        "backend_selection_sha256",
        "one_joint_draw") else character()),
    "dsvert-joint-dp-vector-chunk-start-v5" = c(
      "chunk_index", "coordinate_offset", "coordinates_in_chunk",
      if (isTRUE(profile$exact_gc)) c(
        "backend_selection", "backend_selection_sha256", "binding",
        "binding_sha256", "operation_id", "purpose", "initialization",
        "source_share_exposed", "private_seed_exposed",
        "preclamp_values_exposed") else c(
        "noised_share_sha256", "sampler_contract_hash"),
      "implementation_delta_numerator",
      "implementation_delta_denominator", "intermediate_payload_exposed",
      "capability_available"),
    "dsvert-joint-dp-vector-local-result-v5" = c(
      "local_chunk_commitments", "local_chunk_set_root",
      "local_chunk_set_sha256", "all_chunks_durable",
      "intermediate_payload_exposed", "capability_available"),
    "dsvert-joint-dp-vector-release-root-v5" = c(
      "result_set_hash", "final_vector_root", "final_chunk_hashes",
      "output_lattice_bits", "output_lattice_scale", "mechanism",
      "epsilon", "delta", "implementation_delta_numerator",
      "implementation_delta_denominator", "delta_aggregation",
      "postprocessing", "intermediate_payload_exposed", "durable_replay",
      "capability_available"),
    "dsvert-joint-dp-vector-finalization-ack-v5" = c(
      "final_vector_root", "source_intermediates_compacted",
      "sampler_intermediates_compacted", "final_chunks_retained",
      "durable_replay_retained", "idempotent"),
    NULL)
  if (is.null(extras)) stop("Unknown vector receipt version", call. = FALSE)
  c(common, extras)
}

.dsvert_vector_verify_receipt <- function(
    receipt_json, version, phase, peer, context, contract = NULL,
    profile = NULL) {
  receipt <- .dsvert_joint_dp_client_decode(
    receipt_json, "vector signed receipt",
    .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES)
  expected_profile <- if (!is.null(contract)) contract$profile else profile
  expected <- .dsvert_vector_receipt_fields(version, expected_profile)
  valid <- is.list(expected_profile) &&
    .dsvert_vector_exact_names(receipt, expected) &&
    identical(receipt$version, version) &&
    identical(receipt$phase, phase) &&
    identical(receipt$peer_name, peer) &&
    identical(receipt$peer_identity_pk, unname(context$pinset[[peer]])) &&
    .dsvert_vector_hex(receipt$capsule_id) &&
    .dsvert_vector_hex(receipt$release_instance_id) &&
    is.list(receipt$release_instance) &&
    identical(.dsvert_vector_hash(receipt$release_instance),
              receipt$release_instance_id) &&
    .dsvert_vector_hex(receipt$release_contract_hash) &&
    .dsvert_vector_hex(receipt$transcript_hash) &&
    .dsvert_vector_whole(receipt$coordinate_count, 1) &&
    .dsvert_vector_whole(receipt$chunk_count, 1) &&
    identical(receipt$backend, expected_profile$backend) &&
    identical(receipt$sampler, expected_profile$sampler) &&
    identical(receipt$history_gate, TRUE) &&
    identical(receipt$request_limit, FALSE) &&
    identical(receipt$operation_limit, TRUE)
  if (!is.null(contract)) {
    valid <- valid && identical(receipt$capsule_id, contract$capsule_id) &&
      identical(receipt$release_instance_id,
                contract$release_instance_id) &&
      identical(.dsvert_joint_dp_client_canonical(
                  receipt$release_instance),
                .dsvert_joint_dp_client_canonical(
                  contract$release_instance)) &&
      identical(receipt$release_contract_hash,
                contract$release_contract_hash) &&
      identical(receipt$transcript_hash, contract$transcript_hash) &&
      identical(as.numeric(receipt$coordinate_count),
                as.numeric(contract$coordinate_count)) &&
      identical(as.numeric(receipt$chunk_count),
                as.numeric(contract$chunk_count))
  }
  if (!isTRUE(valid)) {
    stop("A pinned peer returned a misbound vector receipt", call. = FALSE)
  }
  unsigned <- receipt[setdiff(names(receipt), "signature")]
  message <- charToRaw(paste0(
    .DSVERT_CLIENT_JOINT_DP_RECEIPT_DOMAIN,
    .dsvert_joint_dp_client_json(unsigned)))
  public <- .dsvert_joint_dp_client_b64url(
    unname(context$pinset[[peer]]), 32L, "vector identity public key")
  signature <- .dsvert_joint_dp_client_b64url(
    receipt$signature, 64L, "vector receipt signature")
  verified <- tryCatch(openssl::ed25519_verify(
    message, signature, openssl::read_ed25519_pubkey(public)),
    error = function(error) FALSE)
  if (!isTRUE(verified)) {
    stop("A pinned peer returned an invalid vector signature", call. = FALSE)
  }
  receipt
}

.dsvert_vector_prepare_profile <- function(
    responses, peers, manifest, manifest_sha256) {
  declared <- .dsvert_vector_profile(
    manifest$workload$capsule_mechanism,
    manifest$workload$mechanism_selection)
  if (isTRUE(declared$gaussian)) return(declared)
  decoded <- lapply(peers, function(peer) {
    .dsvert_joint_dp_client_decode(
      responses[[peer]], "vector signed prepare receipt",
      .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES)
  })
  backends <- vapply(decoded, function(value) {
    if (!is.list(value) || !.dsvert_vector_string(
          value$backend, maximum_bytes = 256L)) {
      stop("A vector prepare has no signed backend", call. = FALSE)
    }
    value$backend
  }, character(1L))
  if (length(unique(backends)) != 1L) {
    stop("The designated peers disagree on the selected vector backend",
         call. = FALSE)
  }
  profile <- .dsvert_vector_profile(
    manifest$workload$capsule_mechanism,
    manifest$workload$mechanism_selection, backend = backends[[1L]])
  if (isTRUE(profile$selection_bound)) {
    invisible(lapply(decoded, function(value) {
      if (!is.list(value$backend_assessment) ||
          !is.list(value$backend_selection) ||
          !identical(value$backend, value$backend_selection$backend)) {
        stop("A Laplace vector prepare has no coherent backend selection",
             call. = FALSE)
      }
      .dsvert_joint_dp_vector_exact_gc_client_assessment(
        value$backend_assessment, manifest_sha256)
      .dsvert_joint_dp_vector_exact_gc_client_selection(
        value$backend_selection, manifest_sha256,
        require_exact = isTRUE(profile$exact_gc))
      if (!identical(value$backend_selection$assessment_sha256,
                     value$backend_assessment$assessment_sha256)) {
        stop("A vector backend selection is detached from its assessment",
             call. = FALSE)
      }
      TRUE
    }))
  }
  profile
}

.dsvert_vector_prepare_set <- function(
    responses, context, manifest, release_instance, manifest_sha256) {
  peers <- context$designated
  if (!is.list(responses) || !setequal(names(responses), peers)) {
    stop("Vector prepare did not cover both designated peers", call. = FALSE)
  }
  profile <- .dsvert_vector_prepare_profile(
    responses, peers, manifest, manifest_sha256)
  receipts <- stats::setNames(lapply(peers, function(peer) {
    .dsvert_vector_verify_receipt(
      responses[[peer]], .DSVERT_CLIENT_VECTOR_PREPARE_VERSION,
      "vector_prepared", peer, context, profile = profile)
  }), peers)
  stable <- c(
    "capsule_id", "release_instance_id", "release_instance",
    "release_contract_hash", "transcript_hash",
    "coordinate_count", "chunk_count", "backend", "sampler",
    "history_gate", "request_limit", "operation_limit",
    "source_contract_hash",
    "coordinate_order_sha256", "lattice_transform_sha256",
    "mechanism_plan", "plan_sha256",
    "epsilon", "allocated_delta", "sensitivity_steps",
    "complete_epsilon_per_peer", "delta_aggregation",
    "capability_available",
    if (isTRUE(profile$selection_bound)) c(
      "backend_assessment", "backend_selection",
      "backend_selection_sha256",
      "one_joint_draw") else character())
  if (length(unique(vapply(receipts, function(value) {
        .dsvert_joint_dp_client_json(value[stable])
      }, character(1L)))) != 1L) {
    stop("The designated peers disagree on the vector prepare contract",
         call. = FALSE)
  }
  reference <- receipts[[1L]]
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  if (!identical(reference$coordinate_order_sha256, layout$sha256)) {
    stop("The signed vector coordinate order does not match the manifest",
         call. = FALSE)
  }
  policy <- context$status[[peers[[1L]]]]$policy
  epsilon <- .dsvert_vector_decimal(reference$epsilon, 0, 8, TRUE)
  delta <- .dsvert_vector_decimal(reference$allocated_delta, 0, 1, TRUE)
  coordinate_count <- as.numeric(manifest$workload$coordinate_count)
  .dsvert_vector_plan_validate(
    reference$mechanism_plan, reference$plan_sha256, profile,
    coordinate_count, reference$sensitivity_steps)
  chunk_coordinates <- if (isTRUE(profile$exact_gc)) {
    as.integer(reference$mechanism_plan$maximum_chunk_coordinates)
  } else {
    .DSVERT_CLIENT_VECTOR_CHUNK_COORDINATES
  }
  chunk_count <- ceiling(
    coordinate_count / chunk_coordinates)
  valid <- identical(reference$capsule_id,
                     manifest$capsule_identity$capsule_id) &&
    is.list(release_instance) &&
    identical(reference$release_instance_id, release_instance$id) &&
    identical(.dsvert_joint_dp_client_canonical(reference$release_instance),
              release_instance$value) &&
    identical(as.numeric(reference$coordinate_count), coordinate_count) &&
    identical(as.numeric(reference$chunk_count), as.numeric(chunk_count)) &&
    .dsvert_vector_hex(reference$source_contract_hash) &&
    all(vapply(reference[c(
      "coordinate_order_sha256", "lattice_transform_sha256", "plan_sha256")],
      .dsvert_vector_hex, logical(1L))) &&
    (isTRUE(profile$gaussian) ||
       .dsvert_vector_integer_text(reference$sensitivity_steps, TRUE)) &&
    !is.null(epsilon) && !is.null(delta) &&
    .dsvert_dp_num_equal(epsilon, policy$capsule_epsilon) &&
    .dsvert_dp_num_equal(delta, policy$capsule_delta) &&
    identical(reference$complete_epsilon_per_peer,
              profile$complete_epsilon_per_peer) &&
    identical(reference$delta_aggregation,
              profile$delta_aggregation) &&
    identical(reference$capability_available, TRUE) &&
    all(vapply(receipts, function(value) {
      .dsvert_vector_hex(value$commitment_context) &&
        .dsvert_vector_hex(value$seed_commitment)
    }, logical(1L)))
  if (isTRUE(valid) && isTRUE(profile$selection_bound)) {
    valid <- .dsvert_vector_hex(manifest_sha256) &&
      identical(reference$one_joint_draw, isTRUE(profile$exact_gc)) &&
      identical(reference$backend_selection_sha256,
                reference$backend_selection$selection_sha256) &&
      identical(reference$backend_selection$assessment_sha256,
                reference$backend_assessment$assessment_sha256) &&
      identical(reference$backend_selection$exact_gc_plan_sha256,
                reference$backend_assessment$plan_sha256) &&
      identical(as.numeric(
        reference$backend_selection$total_coordinate_count),
        coordinate_count)
    if (isTRUE(valid)) {
      .dsvert_joint_dp_vector_exact_gc_client_assessment(
        reference$backend_assessment, manifest_sha256)
      .dsvert_joint_dp_vector_exact_gc_client_selection(
        reference$backend_selection, manifest_sha256,
        require_exact = isTRUE(profile$exact_gc))
    }
  }
  if (!isTRUE(valid)) {
    stop("The signed vector prepare contract is invalid", call. = FALSE)
  }
  list(
    receipts = receipts, json = responses[peers],
    contract = list(
      capsule_id = reference$capsule_id,
      release_instance_id = reference$release_instance_id,
      release_instance = reference$release_instance,
      release_contract_hash = reference$release_contract_hash,
      transcript_hash = reference$transcript_hash,
      manifest_sha256 = manifest_sha256,
      coordinate_count = as.integer(coordinate_count),
      chunk_coordinates = chunk_coordinates,
      chunk_count = as.integer(chunk_count),
      source_contract_hash = reference$source_contract_hash,
      coordinate_order_sha256 = layout$sha256,
      mechanism_plan = reference$mechanism_plan,
      plan_sha256 = reference$plan_sha256,
      backend_selection = if (isTRUE(profile$exact_gc)) {
        reference$backend_selection
      } else if (isTRUE(profile$selection_bound)) {
        reference$backend_selection
      } else NULL,
      backend_assessment = if (isTRUE(profile$selection_bound)) {
        reference$backend_assessment
      } else NULL,
      backend_selection_sha256 = if (isTRUE(profile$selection_bound)) {
        reference$backend_selection_sha256
      } else NULL,
      one_joint_draw = isTRUE(profile$exact_gc),
      capsule_mechanism = profile$mechanism,
      profile = profile,
      epsilon = reference$epsilon, delta = reference$allocated_delta,
      sensitivity_steps = reference$sensitivity_steps),
    layout = layout)
}

.dsvert_vector_chunk_geometry <- function(contract, index) {
  if (!.dsvert_vector_whole(index, 0, contract$chunk_count - 1L)) {
    stop("Invalid vector chunk index", call. = FALSE)
  }
  offset <- as.numeric(index) * contract$chunk_coordinates
  list(index = as.integer(index), offset = as.integer(offset),
       count = as.integer(min(contract$chunk_coordinates,
                              contract$coordinate_count - offset)))
}

.dsvert_vector_start_set <- function(
    responses, context, contract, chunk_index) {
  peers <- context$designated
  chunk <- .dsvert_vector_chunk_geometry(contract, chunk_index)
  delta_maximum_bytes <- if (isTRUE(contract$profile$exact_gc)) {
    .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES
  } else .DSVERT_CLIENT_VECTOR_IMPLEMENTATION_DELTA_MAX_BYTES
  receipts <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_vector_verify_receipt(
      responses[[peer]], .DSVERT_CLIENT_VECTOR_START_VERSION,
      "vector_chunk_noised", peer, context, contract)
    valid <- identical(as.numeric(value$chunk_index),
                       as.numeric(chunk$index)) &&
      identical(as.numeric(value$coordinate_offset),
                as.numeric(chunk$offset)) &&
      identical(as.numeric(value$coordinates_in_chunk),
                as.numeric(chunk$count)) &&
      .dsvert_vector_integer_text(
        value$implementation_delta_numerator,
        maximum_bytes = delta_maximum_bytes) &&
      .dsvert_vector_integer_text(
        value$implementation_delta_denominator, TRUE,
        maximum_bytes = delta_maximum_bytes) &&
      identical(value$intermediate_payload_exposed, FALSE) &&
      identical(value$capability_available, TRUE)
    if (isTRUE(valid) && isTRUE(contract$profile$exact_gc)) {
      valid <- identical(
          .dsvert_joint_dp_client_canonical(value$backend_selection),
          .dsvert_joint_dp_client_canonical(
            contract$backend_selection)) &&
        identical(value$backend_selection_sha256,
                  contract$backend_selection_sha256) &&
        identical(value$binding_sha256,
                  value$binding$binding_sha256) &&
        identical(value$operation_id, value$binding$operation_id) &&
        identical(value$purpose, value$binding$purpose) &&
        is.list(value$initialization) &&
        identical(value$source_share_exposed, FALSE) &&
        identical(value$private_seed_exposed, FALSE) &&
        identical(value$preclamp_values_exposed, FALSE)
      if (isTRUE(valid)) {
        .dsvert_joint_dp_vector_exact_gc_client_selection(
          value$backend_selection, contract$manifest_sha256)
        .dsvert_joint_dp_vector_exact_gc_client_binding(
          value$binding, contract$manifest_sha256,
          contract$release_contract_hash,
          contract$backend_selection_sha256,
          contract$transcript_hash, chunk$index, chunk$count)
      }
    } else if (isTRUE(valid)) {
      valid <- .dsvert_vector_hex(value$noised_share_sha256) &&
        .dsvert_vector_hex(value$sampler_contract_hash)
    }
    if (!isTRUE(valid)) stop("Invalid signed vector chunk receipt",
                             call. = FALSE)
    value
  }), peers)
  accounting <- vapply(receipts, function(value) paste0(
    value$implementation_delta_numerator, "/",
    value$implementation_delta_denominator), character(1L))
  expected_accounting <- .dsvert_vector_implementation_delta(contract)
  if (length(unique(accounting)) != 1L ||
      !identical(accounting[[1L]], expected_accounting)) {
    stop("The designated peers disagree on vector implementation delta",
         call. = FALSE)
  }
  if (isTRUE(contract$profile$exact_gc)) {
    bindings <- vapply(receipts, function(value) {
      .dsvert_joint_dp_client_json(value$binding)
    }, character(1L))
    if (length(unique(bindings)) != 1L) {
      stop("The designated peers disagree on the exact-GC binding",
           call. = FALSE)
    }
    binding <- receipts[[1L]]$binding
    initializations <-
      .dsvert_joint_dp_vector_exact_gc_initializations(
        receipts, peers, binding)
  } else {
    binding <- NULL
    initializations <- NULL
  }
  list(receipts = receipts, implementation_delta = accounting[[1L]],
       binding = binding, initializations = initializations)
}

.dsvert_vector_result_set <- function(responses, context, contract) {
  peers <- context$designated
  receipts <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_vector_verify_receipt(
      responses[[peer]], .DSVERT_CLIENT_VECTOR_RESULT_VERSION,
      "vector_local_result_committed", peer, context, contract)
    commitments <- unlist(value$local_chunk_commitments, use.names = FALSE)
    valid <- is.character(commitments) &&
      length(commitments) == contract$chunk_count && !anyNA(commitments) &&
      all(vapply(commitments, .dsvert_vector_hex, logical(1L))) &&
      identical(value$local_chunk_set_root,
                .dsvert_vector_merkle_root(commitments)) &&
      identical(value$local_chunk_set_sha256, .dsvert_vector_hash(list(
        protocol = "dsvert-joint-dp-vector-local-chunk-set-v3",
        peer_name = peer, commitments = as.list(commitments)))) &&
      identical(value$all_chunks_durable, TRUE) &&
      identical(value$intermediate_payload_exposed, FALSE) &&
      identical(value$capability_available, TRUE)
    if (!isTRUE(valid)) stop("Invalid signed vector result receipt",
                             call. = FALSE)
    value
  }), peers)
  list(
    receipts = receipts, json = responses[peers],
    result_set_hash = .dsvert_vector_hash(list(
      protocol = "dsvert-joint-dp-vector-result-set-v3",
      ordered_results = unname(receipts))))
}

.dsvert_vector_release_set <- function(
    responses, context, contract, result_set_hash,
    implementation_delta = NULL) {
  peers <- context$designated
  delta_maximum_bytes <- if (isTRUE(contract$profile$exact_gc)) {
    .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES
  } else .DSVERT_CLIENT_VECTOR_IMPLEMENTATION_DELTA_MAX_BYTES
  receipts <- stats::setNames(lapply(peers, function(peer) {
    value <- .dsvert_vector_verify_receipt(
      responses[[peer]], .DSVERT_CLIENT_VECTOR_RELEASE_VERSION,
      "vector_released", peer, context, contract)
    hashes <- unlist(value$final_chunk_hashes, use.names = FALSE)
    accounting <- paste0(value$implementation_delta_numerator, "/",
                         value$implementation_delta_denominator)
    valid <- identical(value$result_set_hash, result_set_hash) &&
      is.character(hashes) && length(hashes) == contract$chunk_count &&
      !anyNA(hashes) &&
      all(vapply(hashes, .dsvert_vector_hex, logical(1L))) &&
      identical(value$final_vector_root,
                .dsvert_vector_merkle_root(hashes)) &&
      .dsvert_vector_whole(value$output_lattice_bits, 1, 62) &&
      identical(as.numeric(value$output_lattice_scale),
                2^as.numeric(value$output_lattice_bits)) &&
      identical(value$mechanism, contract$profile$release_mechanism) &&
      identical(value$epsilon, contract$epsilon) &&
      identical(value$delta, contract$delta) &&
      .dsvert_vector_integer_text(
        value$implementation_delta_numerator,
        maximum_bytes = delta_maximum_bytes) &&
      .dsvert_vector_integer_text(
        value$implementation_delta_denominator, TRUE,
        maximum_bytes = delta_maximum_bytes) &&
      (is.null(implementation_delta) ||
         identical(accounting, implementation_delta)) &&
      identical(value$delta_aggregation,
                contract$profile$delta_aggregation) &&
      identical(value$postprocessing,
                contract$profile$postprocessing) &&
      identical(value$intermediate_payload_exposed, FALSE) &&
      identical(value$durable_replay, TRUE) &&
      identical(value$capability_available, TRUE)
    if (!isTRUE(valid)) stop("Invalid signed vector release receipt",
                             call. = FALSE)
    value
  }), peers)
  comparable <- vapply(receipts, function(value) {
    .dsvert_joint_dp_client_json(value[c(
      "result_set_hash", "final_vector_root", "final_chunk_hashes",
      "output_lattice_bits", "output_lattice_scale", "mechanism", "epsilon",
      "delta", "implementation_delta_numerator",
      "implementation_delta_denominator", "delta_aggregation",
      "postprocessing")])
  }, character(1L))
  if (length(unique(comparable)) != 1L) {
    stop("The two signed vector releases disagree", call. = FALSE)
  }
  list(receipts = receipts, json = responses[peers],
       reference = receipts[[1L]])
}

.dsvert_vector_final_share <- function(
    response, sender, recipient, context, contract, result_set_hash,
    chunk_index) {
  required <- c(
    "ciphertext", "transfer", "capsule_id", "release_contract_hash",
    "result_set_hash", "chunk_index", "intermediate_payload_exposed",
    "capability_available")
  valid <- sender %in% context$designated &&
    recipient %in% context$designated && !identical(sender, recipient) &&
    .dsvert_vector_exact_names(response, required) &&
    identical(response$capsule_id, contract$capsule_id) &&
    identical(response$release_contract_hash,
              contract$release_contract_hash) &&
    identical(response$result_set_hash, result_set_hash) &&
    identical(as.numeric(response$chunk_index), as.numeric(chunk_index)) &&
    identical(response$intermediate_payload_exposed, FALSE) &&
    identical(response$capability_available, TRUE)
  if (!isTRUE(valid)) {
    stop("A vector final-share response is misbound", call. = FALSE)
  }
  transfer <- .dsvert_validate_typed_blob_transfer(
    response$transfer, response$ciphertext, recipient)
  if (!identical(transfer$sender_name, sender) ||
      !identical(transfer$capability_id,
                 .DSVERT_CLIENT_VECTOR_TYPED_CAPABILITY)) {
    stop("A vector final-share transfer is misrouted", call. = FALSE)
  }
  list(ciphertext = response$ciphertext, transfer = transfer)
}

.dsvert_vector_replay_chunk <- function(
    responses, context, contract, release, index) {
  peers <- context$designated
  if (!is.list(responses) || !setequal(names(responses), peers) ||
      !identical(responses[[peers[[1L]]]], responses[[peers[[2L]]]])) {
    stop("The designated peers returned different vector replay chunks",
         call. = FALSE)
  }
  value <- .dsvert_joint_dp_client_decode(
    responses[[peers[[1L]]]], "vector replay chunk",
    .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES)
  expected <- c(
    "version", "capsule_id", "release_contract_hash", "result_set_hash",
    "final_vector_root", "chunk_hash", "chunk", "merkle_proof",
    "durable_replay", "source_store_read", "sampler_invoked",
    "history_gate", "request_limit", "operation_limit")
  chunk_fields <- c(
    "version", "capsule_id", "release_contract_hash", "chunk_index",
    "chunk_count", "coordinate_offset", "coordinates_in_chunk",
    "output_lattice_bits", "output_lattice_scale", "scaled_values",
    "value_encoding", "postprocessing", "source_values_exposed",
    "preclamp_values_exposed")
  chunk <- value$chunk
  geometry <- .dsvert_vector_chunk_geometry(contract, index)
  hashes <- unlist(release$final_chunk_hashes, use.names = FALSE)
  values <- if (is.list(chunk)) unlist(
    chunk$scaled_values, use.names = FALSE) else NULL
  valid <- .dsvert_vector_exact_names(value, expected) &&
    .dsvert_vector_exact_names(chunk, chunk_fields) &&
    identical(value$version, "dsvert-joint-dp-vector-replay-v4") &&
    identical(value$capsule_id, contract$capsule_id) &&
    identical(value$release_contract_hash,
              contract$release_contract_hash) &&
    identical(value$result_set_hash, release$result_set_hash) &&
    identical(value$final_vector_root, release$final_vector_root) &&
    identical(value$chunk_hash, hashes[[index + 1L]]) &&
    identical(value$chunk_hash, .dsvert_vector_hash(chunk)) &&
    identical(.dsvert_joint_dp_client_canonical(value$merkle_proof),
              .dsvert_joint_dp_client_canonical(
                .dsvert_vector_merkle_proof(hashes, index))) &&
    identical(value$durable_replay, TRUE) &&
    identical(value$source_store_read, FALSE) &&
    identical(value$sampler_invoked, FALSE) &&
    identical(value$history_gate, TRUE) &&
    identical(value$request_limit, FALSE) &&
    identical(value$operation_limit, TRUE) &&
    identical(chunk$version, .DSVERT_CLIENT_VECTOR_CHUNK_VERSION) &&
    identical(chunk$capsule_id, contract$capsule_id) &&
    identical(chunk$release_contract_hash,
              contract$release_contract_hash) &&
    identical(as.numeric(chunk$chunk_index), as.numeric(index)) &&
    identical(as.numeric(chunk$chunk_count),
              as.numeric(contract$chunk_count)) &&
    identical(as.numeric(chunk$coordinate_offset),
              as.numeric(geometry$offset)) &&
    identical(as.numeric(chunk$coordinates_in_chunk),
              as.numeric(geometry$count)) &&
    identical(as.numeric(chunk$output_lattice_bits),
              as.numeric(release$output_lattice_bits)) &&
    identical(as.numeric(chunk$output_lattice_scale),
              as.numeric(release$output_lattice_scale)) &&
    is.character(values) && length(values) == geometry$count &&
    !anyNA(values) &&
    all(vapply(values, .dsvert_vector_integer_text, logical(1L))) &&
    identical(chunk$value_encoding,
              "nonnegative-decimal-integer-common-lattice-v1") &&
    identical(chunk$postprocessing,
              contract$profile$postprocessing) &&
    identical(chunk$source_values_exposed, FALSE) &&
    identical(chunk$preclamp_values_exposed, FALSE)
  if (!isTRUE(valid)) stop("Invalid final DP vector replay chunk",
                           call. = FALSE)
  values
}

.dsvert_vector_scaled_to_double <- function(values, scale) {
  limit <- openssl::bignum("9007199254740992")
  result <- vapply(values, function(value) {
    integer <- tryCatch(openssl::bignum(value), error = function(error) NULL)
    if (is.null(integer) || integer > limit) {
      stop("A final DP coordinate exceeds the exact client integer domain",
           call. = FALSE)
    }
    number <- suppressWarnings(as.numeric(as.character(integer)))
    if (!is.finite(number) || number < 0 || number != floor(number)) {
      stop("A final DP coordinate is not exactly representable", call. = FALSE)
    }
    number / scale
  }, numeric(1L))
  unname(result)
}

.dsvert_vector_collect_replay <- function(contract, read_chunk) {
  if (!is.function(read_chunk)) {
    stop("Invalid final DP vector replay reader", call. = FALSE)
  }
  scaled <- character(contract$coordinate_count)
  for (index in seq.int(0L, contract$chunk_count - 1L)) {
    geometry <- .dsvert_vector_chunk_geometry(contract, index)
    values <- read_chunk(index)
    if (!is.character(values) || length(values) != geometry$count ||
        anyNA(values)) {
      stop("A final DP vector replay chunk has the wrong shape",
           call. = FALSE)
    }
    destination <- seq.int(
      geometry$offset + 1L, length.out = geometry$count)
    scaled[destination] <- values
  }
  scaled
}

.dsvert_vector_ack_set <- function(
    responses, context, contract, release) {
  servers <- context$servers
  acks <- stats::setNames(lapply(servers, function(peer) {
    value <- .dsvert_vector_verify_receipt(
      responses[[peer]], .DSVERT_CLIENT_VECTOR_ACK_VERSION,
      "vector_finalized_and_compacted", peer, context, contract)
    valid <- identical(value$final_vector_root,
                       release$final_vector_root) &&
      identical(value$source_intermediates_compacted, TRUE) &&
      identical(value$sampler_intermediates_compacted, TRUE) &&
      identical(value$final_chunks_retained, TRUE) &&
      identical(value$durable_replay_retained, TRUE) &&
      identical(value$idempotent, TRUE)
    if (!isTRUE(valid)) stop("Invalid vector finalization acknowledgement",
                             call. = FALSE)
    value
  }), servers)
  invisible(acks)
}

.dsvert_vector_try_phase <- function(code) {
  tryCatch(
    list(ok = TRUE, value = force(code)),
    dsvert_phase_not_ready = function(error) {
      list(ok = FALSE, value = NULL)
    })
}

# One complete vector attempt. The public internal wrapper below is solely
# responsible for the bounded pre-claim current-root refresh.
.dsvert_joint_dp_vector_capsule_once <- function(
    datasources, status = NULL, manifest_bundle = NULL,
    .aggregate = DSI::datashield.aggregate,
    .source_transport = .dsvert_dp_capsule_source_transport_context,
    .cross_orchestrate = .dsvert_dp_cross_orchestrate,
    .setup_transport = .dsvert_setup_peer_transport,
    .setup_exact = .dsvert_setup_exact_gc_transport,
    .run_exact = .dsvert_joint_dp_vector_exact_gc_run,
    .store_typed = .dsvert_store_typed_blob,
    .allocation = .dsvert_joint_dp_vector_allocation) {
  if (is.null(manifest_bundle)) {
    manifest_bundle <- .dsvert_dp_capsule_manifest_build(
      datasources, status = status, .aggregate = .aggregate)
  }
  context <- manifest_bundle$context
  if (!is.list(context) || !identical(names(context$all_conns),
                                      context$servers) ||
      !identical(names(context$conns), context$designated) ||
      !is.character(manifest_bundle$manifest_json)) {
    stop("Invalid biomedical capsule manifest context", call. = FALSE)
  }
  manifest <- .dsvert_dp_capsule_source_manifest(
    manifest_bundle$manifest_json, context)
  if (!identical(manifest$capsule_id, manifest_bundle$capsule_id)) {
    stop("The biomedical capsule identity changed after its signed manifest",
         call. = FALSE)
  }
  if (!is.function(.allocation)) {
    stop("Invalid biomedical capsule allocation implementation",
         call. = FALSE)
  }
  context <- .allocation(
    manifest_bundle$manifest_json, manifest_bundle$capsule_id,
    context, .aggregate)
  if (!is.list(context$allocation_openings) ||
      !identical(names(context$allocation_openings), context$designated) ||
      !is.list(context$allocation_certificate) ||
      !identical(context$allocation_certificate$capsule_id,
                 manifest_bundle$capsule_id) ||
      !identical(context$allocation_certificate$relay_is_authority, FALSE) ||
      !identical(context$allocation_certificate$capability_available, FALSE)) {
    stop("The biomedical capsule lacks a verified cross-signed allocation",
         call. = FALSE)
  }
  peers <- context$designated
  release_instance <- .dsvert_vector_release_instance(
    context, manifest_bundle$capsule_id)
  allocation_openings <- context$allocation_openings
  prepare_calls <- stats::setNames(lapply(peers, function(peer) call(
    name = "dsvertJointDPVectorPrepareDS",
    manifest_json = manifest_bundle$manifest_json,
    release_instance_json = release_instance$json,
    first_allocation_opening_json =
      allocation_openings[[peers[[1L]]]],
    second_allocation_opening_json =
      allocation_openings[[peers[[2L]]]])), peers)
  prepare_raw <- .dsvert_fanout_by_site(
    context$conns, prepare_calls, operation = "joint-DP vector prepare",
    .aggregate = .aggregate)
  prepared <- .dsvert_vector_prepare_set(
    prepare_raw, context, manifest$value, release_instance,
    manifest_bundle$manifest_sha256)
  contract <- prepared$contract

  session_id <- .dsvert_uuid4()
  exact_transport_ready <- FALSE
  compute_indices <- match(peers, context$servers)
  if (isTRUE(contract$profile$exact_gc) &&
      (length(compute_indices) != 2L || anyNA(compute_indices) ||
       anyDuplicated(compute_indices))) {
    stop("The exact-GC compute peers are absent from the DSI context",
         call. = FALSE)
  }
  result_calls <- stats::setNames(lapply(peers, function(peer) call(
    name = "dsvertJointDPVectorResultDS",
    manifest_json = manifest_bundle$manifest_json,
    first_prepare_json = prepared$json[[peers[[1L]]]],
    second_prepare_json = prepared$json[[peers[[2L]]]])), peers)
  existing_result <- .dsvert_vector_try_phase(.dsvert_fanout_by_site(
    context$conns, result_calls, operation = "joint-DP vector result lookup",
    .aggregate = .aggregate))
  implementation_delta <- .dsvert_vector_implementation_delta(contract)
  source_receipt <- NULL
  if (isTRUE(existing_result$ok)) {
    results <- .dsvert_vector_result_set(
      existing_result$value, context, contract)
  } else {
    if (!is.function(.source_transport)) {
      stop("Invalid biomedical source-transport implementation", call. = FALSE)
    }
    source_receipt <- .source_transport(
      manifest_bundle$manifest_json, context, .aggregate)
    cross_enabled <- length(
      .dsvert_dp_gaussian_cross_artifacts_client(manifest$value)) > 0L ||
      length(.dsvert_dp_categorical_cross_artifacts_client(
        manifest$value)) > 0L
    release_coordinate_count <- if (
        is.null(source_receipt$release_coordinate_count)) {
      source_receipt$coordinate_count
    } else {
      source_receipt$release_coordinate_count
    }
    valid_source <- is.list(source_receipt) &&
      identical(source_receipt$capsule_id, contract$capsule_id) &&
      identical(source_receipt$contract_hash,
                contract$source_contract_hash) &&
      identical(as.numeric(release_coordinate_count),
                as.numeric(contract$coordinate_count)) &&
      identical(source_receipt$payload_exposed, FALSE)
    if (!isTRUE(valid_source)) {
      stop("The biomedical capsule source is not bound to the vector prepare",
           call. = FALSE)
    }
    if (isTRUE(cross_enabled)) {
      if (!is.function(.cross_orchestrate)) {
        stop("Invalid cross-owner exact orchestration",
             call. = FALSE)
      }
      cross_receipt <- .cross_orchestrate(
        manifest_bundle$manifest_json, manifest$value, context,
        source_receipt, .aggregate)
      if (!is.list(cross_receipt) ||
          !identical(cross_receipt$enabled, TRUE) ||
          !identical(cross_receipt$sampler_handoff_ready, TRUE) ||
          !identical(cross_receipt$exact_intermediates_exposed, FALSE) ||
          !identical(cross_receipt$source_values_exposed, FALSE)) {
        stop("The exact cross-owner source is not ready for sampling",
             call. = FALSE)
      }
    } else if (
        !identical(as.numeric(source_receipt$coordinate_count),
                   as.numeric(contract$coordinate_count)) ||
        !identical(as.numeric(source_receipt$chunk_count),
                   as.numeric(contract$chunk_count)) ||
        !identical(source_receipt$sampler_handoff_ready, TRUE)) {
      stop("The biomedical capsule source is not ready for vector sampling",
           call. = FALSE)
    }
    if (isTRUE(contract$profile$exact_gc)) {
      if (!is.function(.setup_exact) || !is.function(.run_exact)) {
        stop("Invalid one-draw exact-GC transport implementation",
             call. = FALSE)
      }
      .setup_exact(
        context$all_conns, context$servers, compute_indices, session_id,
        .aggregate = .aggregate)
      exact_transport_ready <- TRUE
    }
    for (index in seq.int(0L, contract$chunk_count - 1L)) {
      start_calls <- stats::setNames(lapply(peers, function(peer) call(
        name = "dsvertJointDPVectorStartDS",
        manifest_json = manifest_bundle$manifest_json,
        first_prepare_json = prepared$json[[peers[[1L]]]],
        second_prepare_json = prepared$json[[peers[[2L]]]],
        chunk_index = as.integer(index),
        session_id = session_id)), peers)
      started <- .dsvert_vector_start_set(.dsvert_fanout_by_site(
        context$conns, start_calls,
        operation = "joint-DP vector sticky chunk sampling",
        .aggregate = .aggregate), context, contract, index)
      if (!identical(implementation_delta,
                     started$implementation_delta)) {
        stop("Vector sampler accounting changed between chunks",
             call. = FALSE)
      }
      if (isTRUE(contract$profile$exact_gc)) {
        .run_exact(
          context$all_conns, server_names = context$servers,
          servers = compute_indices, session_id = session_id,
          binding = started$binding,
          manifest_sha256 = contract$manifest_sha256,
          release_contract_hash = contract$release_contract_hash,
          selection_sha256 = contract$backend_selection_sha256,
          transcript_hash = contract$transcript_hash,
          chunk_index = as.integer(index),
          coordinate_count =
            .dsvert_vector_chunk_geometry(contract, index)$count,
          initialized = started$initializations,
          transport_ready = TRUE, .aggregate = .aggregate)
      }
    }
    commit_result_calls <- stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertJointDPVectorResultDS",
      manifest_json = manifest_bundle$manifest_json,
      first_prepare_json = prepared$json[[peers[[1L]]]],
      second_prepare_json = prepared$json[[peers[[2L]]]],
      session_id = session_id)), peers)
    results <- .dsvert_vector_result_set(.dsvert_fanout_by_site(
      context$conns, commit_result_calls,
      operation = "joint-DP vector result commit",
      .aggregate = .aggregate), context, contract)
  }

  release_calls <- stats::setNames(lapply(peers, function(peer) call(
    name = "dsvertJointDPVectorReleaseDS", session_id = session_id,
    manifest_json = manifest_bundle$manifest_json,
    first_result_json = results$json[[peers[[1L]]]],
    second_result_json = results$json[[peers[[2L]]]])), peers)
  existing_release <- .dsvert_vector_try_phase(.dsvert_fanout_by_site(
    context$conns, release_calls, operation = "joint-DP vector release lookup",
    .aggregate = .aggregate))
  if (isTRUE(existing_release$ok)) {
    releases <- .dsvert_vector_release_set(
      existing_release$value, context, contract, results$result_set_hash,
      implementation_delta)
  } else {
    if (!is.function(.store_typed) ||
        (!isTRUE(contract$profile$exact_gc) &&
         !is.function(.setup_transport)) ||
        (isTRUE(contract$profile$exact_gc) &&
         !is.function(.setup_exact))) {
      stop("Invalid vector peer-transport implementation", call. = FALSE)
    }
    if (isTRUE(contract$profile$exact_gc)) {
      if (!isTRUE(exact_transport_ready)) {
        .setup_exact(
          context$all_conns, context$servers, compute_indices, session_id,
          .aggregate = .aggregate)
        exact_transport_ready <- TRUE
      }
    } else {
      .setup_transport(
        context$all_conns, context$servers, context$servers, session_id,
        .aggregate = .aggregate)
    }
    for (index in seq.int(0L, contract$chunk_count - 1L)) {
      share_calls <- stats::setNames(lapply(peers, function(peer) call(
        name = "dsvertJointDPVectorFinalShareDS", session_id = session_id,
        manifest_json = manifest_bundle$manifest_json,
        first_result_json = results$json[[peers[[1L]]]],
        second_result_json = results$json[[peers[[2L]]]],
        chunk_index = as.integer(index))), peers)
      shares <- .dsvert_fanout_by_site(
        context$conns, share_calls,
        operation = "joint-DP vector encrypted final-share fan-out",
        .aggregate = .aggregate)
      for (sender in peers) {
        recipient <- setdiff(peers, sender)
        transfer <- .dsvert_vector_final_share(
          shares[[sender]], sender, recipient, context, contract,
          results$result_set_hash, index)
        .store_typed(
          transfer$ciphertext, transfer$transfer,
          context$all_conns[recipient], session_id,
          producer_conn = context$all_conns[sender],
          .aggregate = .aggregate)
      }
    }
    releases <- .dsvert_vector_release_set(.dsvert_fanout_by_site(
      context$conns, release_calls, operation = "joint-DP vector release",
      .aggregate = .aggregate), context, contract, results$result_set_hash,
      implementation_delta)
  }

  reference <- releases$reference
  scaled <- .dsvert_vector_collect_replay(contract, function(index) {
    replay_calls <- stats::setNames(lapply(peers, function(peer) call(
      name = "dsvertJointDPVectorReplayDS",
      manifest_json = manifest_bundle$manifest_json,
      first_release_json = releases$json[[peers[[1L]]]],
      second_release_json = releases$json[[peers[[2L]]]],
      chunk_index = as.integer(index))), peers)
    .dsvert_vector_replay_chunk(
      .dsvert_fanout_by_site(
        context$conns, replay_calls, operation = "joint-DP vector replay",
        .aggregate = .aggregate), context, contract, reference, index)
  })
  if (length(scaled) != contract$coordinate_count) {
    stop("The final DP vector has the wrong coordinate count", call. = FALSE)
  }

  ack_calls <- stats::setNames(lapply(context$servers, function(peer) call(
    name = "dsvertJointDPVectorFinalizeAckDS",
    manifest_json = manifest_bundle$manifest_json,
    first_release_json = releases$json[[peers[[1L]]]],
    second_release_json = releases$json[[peers[[2L]]]])), context$servers)
  finalization_receipts <- .dsvert_vector_ack_set(.dsvert_fanout_by_site(
    context$all_conns, ack_calls,
    operation = "joint-DP vector durable finalization",
    .aggregate = .aggregate), context, contract, reference)

  result <- list(
    version = "dsvert-joint-dp-biomedical-vector-client-v1",
    capsule_id = contract$capsule_id,
    release_instance_id = contract$release_instance_id,
    release_instance = contract$release_instance,
    manifest_sha256 = manifest_bundle$manifest_sha256,
    final_vector_root = reference$final_vector_root,
    coordinate_order_sha256 = contract$coordinate_order_sha256,
    coordinate_count = contract$coordinate_count,
    output_lattice_bits = as.integer(reference$output_lattice_bits),
    output_lattice_scale = as.numeric(reference$output_lattice_scale),
    values = .dsvert_vector_scaled_to_double(
      scaled, as.numeric(reference$output_lattice_scale)),
    backend = contract$profile$backend,
    sampler = contract$profile$sampler,
    one_joint_draw = isTRUE(contract$one_joint_draw),
    backend_selection = contract$backend_selection,
    backend_assessment = contract$backend_assessment,
    mechanism_plan = contract$mechanism_plan,
    plan_sha256 = contract$plan_sha256,
    postprocessing = contract$profile$postprocessing,
    mechanism = reference$mechanism,
    epsilon = as.numeric(reference$epsilon),
    delta = as.numeric(reference$delta),
    implementation_delta = paste0(
      reference$implementation_delta_numerator, "/",
      reference$implementation_delta_denominator),
    delta_aggregation = reference$delta_aggregation,
    sticky_replay = TRUE, source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    manifest = manifest$value,
    signed_provenance = list(
      version = "dsvert-joint-dp-vector-public-provenance-v1",
      ordered_peer_pinset = as.list(context$pinset),
      peer_pinset_sha256 =
        context$status[[context$servers[[1L]]]]$policy$peer_pinset_sha256,
      designated_noise_peers = as.list(context$designated),
      allocation_certificate = context$allocation_certificate,
      release_instance_id = contract$release_instance_id,
      release_instance = contract$release_instance,
      prepare_receipts = prepared$receipts,
      release_receipts = releases$receipts,
      finalization_receipts = finalization_receipts,
      protected_shares_included = FALSE,
      preclamp_values_included = FALSE,
      source_values_included = FALSE))
  class(result) <- c("dsvert_joint_dp_vector", "list")
  result
}

.dsvert_joint_dp_vector_retryable_instance_error <- function(error) {
  (inherits(error, "dsvert_release_instance_retry") ||
     inherits(error, "dsvert_dp_release_retry_current_instance")) && grepl(
    .DSVERT_CLIENT_VECTOR_RETRY_CURRENT_INSTANCE_TOKEN,
    conditionMessage(error), fixed = TRUE)
}

#' Materialize or replay the one server-authoritative DP capsule vector
#'
#' Internal orchestration. There is no analyst-selected epsilon, bound,
#' workload, peer, source, seed, allocation index, request counter or fallback
#' mechanism. An authenticated stale-root signal emitted before the first valid
#' START refreshes signed roots and retries the immutable handshake at most
#' once. After START has claimed an instance, another instance fails closed.
#'
#' @keywords internal
.dsvert_joint_dp_vector_capsule <- function(
    datasources, status = NULL, manifest_bundle = NULL,
    .aggregate = DSI::datashield.aggregate,
    .source_transport = .dsvert_dp_capsule_source_transport_context,
    .cross_orchestrate = .dsvert_dp_cross_orchestrate,
    .setup_transport = .dsvert_setup_peer_transport,
    .setup_exact = .dsvert_setup_exact_gc_transport,
    .run_exact = .dsvert_joint_dp_vector_exact_gc_run,
    .store_typed = .dsvert_store_typed_blob,
    .allocation = .dsvert_joint_dp_vector_allocation,
    .retry_refresh = NULL) {
  run <- function(current_status, current_manifest) {
    .dsvert_joint_dp_vector_capsule_once(
      datasources, status = current_status,
      manifest_bundle = current_manifest,
      .aggregate = .aggregate,
      .source_transport = .source_transport,
      .cross_orchestrate = .cross_orchestrate,
      .setup_transport = .setup_transport,
      .setup_exact = .setup_exact,
      .run_exact = .run_exact,
      .store_typed = .store_typed,
      .allocation = .allocation)
  }
  first <- tryCatch(run(status, manifest_bundle), error = identity)
  if (!inherits(first, "error")) return(first)
  if (!.dsvert_joint_dp_vector_retryable_instance_error(first)) {
    stop(first)
  }
  if (is.null(.retry_refresh)) {
    refreshed_status <- .dsvert_joint_dp_capsule_status_impl(
      datasources, .aggregate = .aggregate)
    refreshed <- list(
      status = refreshed_status, manifest_bundle = NULL)
  } else {
    if (!is.function(.retry_refresh)) {
      stop("Invalid vector release-instance refresh implementation",
           call. = FALSE)
    }
    refreshed <- .retry_refresh(datasources, .aggregate)
  }
  if (!is.list(refreshed) ||
      !setequal(names(refreshed), c("status", "manifest_bundle"))) {
    stop("The vector release-instance refresh returned an invalid context",
         call. = FALSE)
  }
  second <- tryCatch(
    run(refreshed$status, refreshed$manifest_bundle), error = identity)
  if (inherits(second, "error")) {
    if (.dsvert_joint_dp_vector_retryable_instance_error(second)) {
      stop("The vector release instance remained stale after one automatic root refresh: ",
           conditionMessage(second), call. = FALSE)
    }
    stop(second)
  }
  second
}
