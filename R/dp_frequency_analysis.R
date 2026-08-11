# Internal client preparation for one server-authoritative Frequency analysis.

.DSVERT_CLIENT_DP_FREQUENCY_CLAIM_VERSION <- "dsvert-dp-frequency-factor-claim-v1"
.DSVERT_CLIENT_DP_FREQUENCY_CONFIG_VERSION <- "dsvert-dp-frequency-config-v1"
.DSVERT_CLIENT_DP_FREQUENCY_RECEIPT_VERSION <- "dsvert-dp-frequency-receipt-v1"
.DSVERT_CLIENT_DP_FREQUENCY_PUBLIC_AUTH_VERSION <- "dsvert-dp-frequency-public-authorization-v1"
.DSVERT_CLIENT_DP_FREQUENCY_PREPARED_VERSION <- "dsvert-dp-frequency-prepared-v1"
.DSVERT_CLIENT_DP_FREQUENCY_FACTOR_VERSION <- "dsvert-psi-padded-factor-entry-v1"
.DSVERT_CLIENT_DP_FREQUENCY_CONFIG_MAX_BYTES <- 32L * 1024L^2
.DSVERT_CLIENT_DP_FREQUENCY_RECEIPTS_MAX_BYTES <- 32L * 1024L^2
.dsvert_dp_frequency_client_hash_v1 <- function(domain, value) {
    .dsvert_dp_analysis_frequency_hash_v1(domain, value)
}
.dsvert_dp_frequency_client_wire_json_v1 <- function(value) {
    as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", na = "null", digits = 17, pretty = FALSE))
}
.dsvert_dp_frequency_client_object_v1 <- function(value, fields) {
    is.list(value) && !is.null(names(value)) && !anyNA(names(value)) && !anyDuplicated(names(value)) && setequal(names(value), fields)
}
.dsvert_dp_frequency_client_hex_v1 <- function(value, what) {
    if (!.dsvert_dp_analysis_frequency_hex_v1(value)) {
        stop("Invalid Frequency ", what, ".", call. = FALSE)
    }
    value
}
.dsvert_dp_frequency_client_peer_v1 <- function(value) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$", value)) {
        stop("Invalid Frequency peer name.", call. = FALSE)
    }
    value
}
.dsvert_dp_frequency_client_pins_v1 <- function(value) {
    if (is.list(value))
        value <- unlist(value, use.names = TRUE)
    if (!is.character(value) || length(value) < 2L || length(value) > 4096L || is.null(names(value)) || anyNA(names(value)) || anyDuplicated(names(value))) {
        stop("Invalid Frequency peer pins.", call. = FALSE)
    }
    peers <- tryCatch(vapply(names(value), .dsvert_dp_frequency_client_peer_v1, character(1L)), error = function(error) character())
    pins <- tryCatch(vapply(value, .dsvert_dp_analysis_client_identity_pk, character(1L)), error = function(error) character())
    if (length(peers) != length(value) || length(pins) != length(value) || anyDuplicated(pins))
        stop("Invalid Frequency peer pins.", call. = FALSE)
    names(pins) <- peers
    pins[order(names(pins), method = "radix")]
}
.dsvert_dp_frequency_client_pinset_v1 <- function(pins) {
    pins <- .dsvert_dp_frequency_client_pins_v1(pins)
    encoded <- paste(vapply(names(pins), function(peer) paste0(nchar(peer, type = "bytes"), ":", peer, "=", nchar(pins[[peer]], type = "bytes"), ":", pins[[peer]]), character(1L)),
        collapse = "|")
    paste0("pinset_", digest::digest(encoded, "sha256", serialize = FALSE))
}
.dsvert_dp_frequency_client_text_v1 <- function(value, what) {
    bytes <- tryCatch(nchar(value, type = "bytes"), error = function(e) NA)
    if (!is.character(value) || length(value) != 1L || is.na(value) || is.na(bytes) || bytes < 1L || bytes > 1024L) {
        stop("Invalid Frequency ", what, ".", call. = FALSE)
    }
    value <- enc2utf8(value)
    if (!nzchar(value) || !isTRUE(validUTF8(value))) {
        stop("Invalid Frequency ", what, ".", call. = FALSE)
    }
    value
}
.dsvert_dp_frequency_client_factor_v1 <- function(value) {
    fields <- c("version", "variable_name", "variable_id", "levels", "dimension")
    if (!.dsvert_dp_frequency_client_object_v1(value, fields) || !identical(value$version, .DSVERT_CLIENT_DP_FREQUENCY_FACTOR_VERSION)) {
        stop("Invalid Frequency public factor entry.", call. = FALSE)
    }
    dimension <- .dsvert_dp_analysis_frequency_levels_dimension_v1(value$levels, value$dimension, value$variable_name)
    if (is.null(dimension))
        stop("Invalid Frequency public factor entry.", call. = FALSE)
    variable <- .dsvert_dp_frequency_client_text_v1(value$variable_name, "factor variable")
    levels <- tryCatch(vapply(value$levels, .dsvert_dp_frequency_client_text_v1, character(1L), what = "factor level"), error = function(error) character())
    expected_id <- paste0("var_", digest::digest(paste0("dsVert/psi-padded/factor-variable/v1|", .dsvert_joint_dp_client_json(list(variable_name = variable))), "sha256", serialize = FALSE))
    if (length(levels) != dimension || anyDuplicated(levels) || !identical(levels, sort(levels, method = "radix")) || !identical(value$variable_id, expected_id))
        stop("Invalid Frequency public factor entry.", call. = FALSE)
    list(version = .DSVERT_CLIENT_DP_FREQUENCY_FACTOR_VERSION, variable_name = variable, variable_id = expected_id, levels = as.list(unname(levels)), dimension = as.integer(dimension))
}
.dsvert_dp_frequency_client_factor_hash_v1 <- function(value) {
    entry <- .dsvert_dp_frequency_client_factor_v1(value)
    digest::digest(charToRaw(paste0("dsVert/psi-padded/factor-entry/v1|", .dsvert_dp_frequency_client_wire_json_v1(entry))), "sha256", serialize = FALSE)
}
.dsvert_dp_frequency_client_verify_v1 <- function(message, identity_pk, signature, what) {
    signature_raw <- .dsvert_joint_dp_client_b64url(signature, 64L, paste("Frequency", what, "signature"))
    public_raw <- .dsvert_joint_dp_client_b64url(identity_pk, 32L, paste("Frequency", what, "identity"))
    valid <- tryCatch(openssl::ed25519_verify(message, signature_raw, openssl::read_ed25519_pubkey(public_raw)), error = function(error) FALSE)
    if (!isTRUE(valid))
        stop("Frequency ", what, " signature verification failed.", call. = FALSE)
}
.dsvert_dp_frequency_client_claim_v1 <- function(value, source_owner, peer_pins = NULL) {
    fields <- c("version", "source_peer_name", "source_identity_pk", "psi_run_sha256", "attestation_id", "contract_hash", "source_binding_id", "alignment_hash", "alignment_purpose",
        "dataset_id", "dataset_version", "privacy_unit_column", "pinset_id", "capacity_bucket", "factor_entry", "factor_entry_sha256", "signature")
    if (!.dsvert_dp_frequency_client_object_v1(value, fields) || !identical(value$version, .DSVERT_CLIENT_DP_FREQUENCY_CLAIM_VERSION)) {
        stop("Invalid signed Frequency Claim.", call. = FALSE)
    }
    source_owner <- .dsvert_dp_frequency_client_peer_v1(source_owner)
    peer <- .dsvert_dp_frequency_client_peer_v1(value$source_peer_name)
    identity <- tryCatch(.dsvert_dp_analysis_client_identity_pk(value$source_identity_pk), error = function(error) NULL)
    factor <- tryCatch(.dsvert_dp_frequency_client_factor_v1(value$factor_entry), error = function(error) NULL)
    capacity <- suppressWarnings(as.numeric(value$capacity_bucket))
    ids <- c("alignment_purpose", "dataset_id", "dataset_version")
    source <- list(alignment_purpose = value$alignment_purpose, dataset_id = value$dataset_id, dataset_version = value$dataset_version, id_column = value$privacy_unit_column)
    source_binding <- paste0("source_", digest::digest(.dsvert_dp_frequency_client_wire_json_v1(source), "sha256", serialize = FALSE))
    valid <- identical(peer, source_owner) && !is.null(identity) && !is.null(factor) && all(vapply(value[ids], .dsvert_dp_analysis_client_scalar_id, logical(1L))) && is.character(value$privacy_unit_column) &&
        grepl("^[A-Za-z._][A-Za-z0-9._]{0,127}$", value$privacy_unit_column) && grepl("^attest_[0-9a-f]{64}$", value$attestation_id) && identical(value$source_binding_id, source_binding) &&
        grepl("^pinset_[0-9a-f]{64}$", value$pinset_id) && length(capacity) == 1L && is.finite(capacity) && capacity == floor(capacity) && capacity >= 64 && capacity <= 1048576 && bitwAnd(as.integer(capacity),
        as.integer(capacity) - 1L) == 0L
    for (field in c("psi_run_sha256", "contract_hash", "alignment_hash", "factor_entry_sha256")) {
        valid <- valid && .dsvert_dp_analysis_frequency_hex_v1(value[[field]])
    }
    if (!isTRUE(valid) || !identical(value$factor_entry_sha256, .dsvert_dp_frequency_client_factor_hash_v1(factor)))
        stop("Invalid signed Frequency Claim.", call. = FALSE)
    if (!is.null(peer_pins)) {
        pins <- .dsvert_dp_frequency_client_pins_v1(peer_pins)
        if (!peer %in% names(pins) || !identical(identity, unname(pins[[peer]])) || !identical(value$pinset_id, .dsvert_dp_frequency_client_pinset_v1(pins)))
            stop("The Frequency Claim source identity is not pinned.", call. = FALSE)
    }
    normalized <- value
    normalized$source_peer_name <- peer
    normalized$source_identity_pk <- identity
    normalized$capacity_bucket <- as.integer(capacity)
    normalized$factor_entry <- factor
    unsigned <- .dsvert_dp_analysis_client_canonical_value_v1(normalized[setdiff(names(normalized), "signature")])
    .dsvert_dp_frequency_client_verify_v1(charToRaw(paste0("dsVert/dp-frequency/factor-claim/v1|", .dsvert_joint_dp_client_json(unsigned))), identity, value$signature, "Claim")
    .dsvert_dp_analysis_client_canonical_value_v1(normalized)
}
.dsvert_dp_frequency_client_reduce_v1 <- function(numerator, denominator) {
    numerator <- openssl::bignum(numerator)
    denominator <- openssl::bignum(denominator)
    if (numerator >= denominator)
        return(list(numerator = "1", denominator = "1"))
    left <- numerator
    right <- denominator
    while (right != 0) {
        remainder <- left%%right
        left <- right
        right <- remainder
    }
    list(numerator = as.character(numerator%/%left), denominator = as.character(denominator%/%left))
}
.dsvert_dp_frequency_client_convolution_fields_v1 <- function() c("version", "sampler", "stop_bits", "stop_numerator", "uniform_bits", "binary_geometric_bits", "bernoulli_thresholds",
    "sensitivity_steps", "total_coordinate_count", "epsilon_effective_upper_numerator", "epsilon_effective_upper_denominator", "one_geometric_tv_numerator", "one_geometric_tv_denominator",
    "tail_upper_numerator", "tail_upper_denominator", "rounding_upper_numerator", "rounding_upper_denominator", "implementation_delta_numerator", "implementation_delta_denominator",
    "implementation_delta_bound", "maximum_noise_magnitude", "maximum_chunk_coordinates", "private_stream_bytes_per_coordinate", "accounting", "capability_available", "independent_noise_peer_count",
    "complete_epsilon_per_peer", "epsilon_divided_by_peer_count", "geometric_variables_per_peer_per_coordinate", "geometric_variables_total_per_coordinate", "per_peer_implementation_delta_numerator",
    "per_peer_implementation_delta_denominator", "per_peer_implementation_delta_bound", "release_implementation_delta_aggregation", "two_peer_ideal_transfer_delta_numerator", "two_peer_ideal_transfer_delta_denominator",
    "two_peer_ideal_transfer_delta_bound", "threat_model", "privacy_argument")
.dsvert_dp_frequency_client_accuracy_v1 <- function(certificate, plan, primitive, request, dimension) {
    fields <- c("primitive", "plan_sha256", "event", "method", "release_tv_upper_numerator", "release_tv_upper_denominator", "simultaneous_95_abs", "absolute_support")
    if (!.dsvert_dp_frequency_client_object_v1(certificate, fields) || !identical(certificate$primitive, primitive) || !identical(certificate$event, "max_j_abs_error_gt_radius") ||
        !identical(certificate$plan_sha256, .dsvert_dp_frequency_client_hash_v1("dsVert/frequency/full-plan/v1|", plan)))
        stop("Invalid Frequency selected accuracy certificate.", call. = FALSE)
    gaussian <- grepl("gaussian", primitive, fixed = TRUE)
    peer <- tryCatch(.dsvert_dp_analysis_frequency_uint_v1(plan[[if (gaussian) {
        "maximum_noise_magnitude_per_peer"
    }
    else "maximum_noise_magnitude"]], TRUE), error = function(error) NULL)
    radius <- tryCatch(.dsvert_dp_analysis_frequency_uint_v1(certificate$simultaneous_95_abs), error = function(error) NULL)
    support <- tryCatch(.dsvert_dp_analysis_frequency_uint_v1(certificate$absolute_support, TRUE), error = function(error) NULL)
    tv <- tryCatch(if (gaussian)
        .dsvert_dp_frequency_client_reduce_v1(2 * openssl::bignum(plan$vector_total_tv_upper_numerator), plan$vector_total_tv_upper_denominator)
    else .dsvert_dp_frequency_client_reduce_v1(4 * dimension * openssl::bignum(plan$one_geometric_tv_numerator), plan$one_geometric_tv_denominator), error = function(error) NULL)
    methods <- if (gaussian)
        "gaussian_plan_v2_subgaussian_mgf_tv_transfer"
    else c("exact_two_discrete_laplace_convolution_tail_v1", "dyadic_exponential_envelope_v1", "finite_support_v1")
    if (is.null(peer) || is.null(radius) || is.null(support) || is.null(tv) || !certificate$method %in% methods || !identical(support, 2 * peer) || radius > support || !identical(certificate$release_tv_upper_numerator,
        tv$numerator) || !identical(certificate$release_tv_upper_denominator, tv$denominator) || (gaussian && !identical(certificate$simultaneous_95_abs, plan$simultaneous_95_abs)) ||
        (!gaussian && identical(certificate$method, "finite_support_v1") && radius != support))
        stop("Invalid Frequency selected accuracy certificate.", call. = FALSE)
    if (!gaussian && !identical(certificate$method, "finite_support_v1")) {
        tv_upper <- .dsvert_dp_vector_fraction_upper(tv$numerator, tv$denominator)
        alpha <- .dsvert_dp_vector_next_down(.dsvert_dp_vector_next_down(0.05 - tv_upper))
        target <- if (alpha > 0)
            log(.dsvert_dp_vector_next_down(alpha/dimension))
        else -Inf
        tail <- if (identical(certificate$method, "dyadic_exponential_envelope_v1")) {
            context <- .dsvert_dp_vector_dyadic_tail_context(plan)
            function(value) .dsvert_dp_vector_plan_log_tail_upper(value, context, convolution = TRUE)
        }
        else function(value) .dsvert_dp_vector_convolution_log_tail(value, as.numeric(request$epsilon), as.numeric(request$sensitivity_steps))
        minimal <- is.finite(target) && tail(as.numeric(radius)) <= target && (radius == 0 || tail(as.numeric(radius - 1)) > target)
        if (!isTRUE(minimal) || radius >= support)
            stop("Invalid Frequency selected accuracy certificate.", call. = FALSE)
    }
    invisible(TRUE)
}
.dsvert_dp_frequency_client_plan_v1 <- function(config) {
    selection <- config$backend_selection
    fields <- c("summary", "selected_request", "selected_plan", "selected_accuracy_certificate", "selection_certificate")
    if (!.dsvert_dp_frequency_client_object_v1(selection, fields))
        stop("Invalid Frequency backend selection.", call. = FALSE)
    primitive <- tryCatch(selection$summary$selected_primitive, error = function(error) NULL)
    profile <- .dsvert_dp_analysis_frequency_profile_v1(primitive)
    kind <- if (is.list(profile) && isTRUE(profile$gaussian))
        "gaussian"
    else "convolution"
    requests <- .dsvert_dp_analysis_frequency_candidate_requests_v2(config$privacy, config$calibration, config$factor_domain$dimension)
    request <- requests[[kind]]
    plan <- selection$selected_plan
    plan_hash <- .dsvert_dp_frequency_client_hash_v1("dsVert/frequency/full-plan/v1|", plan)
    candidate <- tryCatch(selection$summary$candidates[[kind]], error = function(error) NULL)
    if (is.null(profile) || !identical(.dsvert_dp_analysis_client_canonical_value_v1(selection$selected_request), .dsvert_dp_analysis_client_canonical_value_v1(request)) || !is.list(candidate) || !identical(candidate$full_plan_sha256, plan_hash))
        stop("Invalid Frequency backend selection.", call. = FALSE)
    vector_profile <- list(exact_gc = FALSE, gaussian = isTRUE(profile$gaussian), plan_version = profile$plan, sampler = profile$sampler, mechanism = if (isTRUE(profile$gaussian)) "dyadic_discrete_gaussian_truncated_tv_bounded" else NULL,
        manifest_selection = if (isTRUE(profile$gaussian)) list(gaussian_calibration_request = request, gaussian_plan = plan, gaussian_plan_sha256 = .dsvert_vector_hash(plan)) else NULL)
    if (!isTRUE(profile$gaussian) && !.dsvert_vector_exact_names(plan, .dsvert_dp_frequency_client_convolution_fields_v1()))
        stop("Invalid Frequency selected physical plan.", call. = FALSE)
    tryCatch(.dsvert_vector_plan_validate(plan, .dsvert_vector_hash(plan), vector_profile, config$factor_domain$dimension, request[[if (isTRUE(profile$gaussian))
        "l2_sensitivity_steps"
    else "sensitivity_steps"]]), error = function(error) stop("Invalid Frequency selected physical plan.", call. = FALSE))
    .dsvert_dp_frequency_client_accuracy_v1(selection$selected_accuracy_certificate, plan, primitive, request, config$factor_domain$dimension)
    certificate_hash <- .dsvert_dp_frequency_client_hash_v1("dsVert/frequency/backend-selection/accuracy-certificate/v1|", selection$selected_accuracy_certificate)
    selection_certificate <- selection$selection_certificate
    expected_certificate <- list(version = "dsvert-joint-dp-frequency-backend-selection-certificate-v1", policy = "minimum_certified_simultaneous_95_abs_convolution_tie_v1", objective = "minimum_certified_simultaneous_95_abs",
        selected_primitive = primitive, selected_plan_sha256 = plan_hash, selected_simultaneous_95_abs = selection$selected_accuracy_certificate$simultaneous_95_abs, tie_break = "convolution_laplace_v3_on_equal_certified_radius",
        input_scope = paste0("public_adjacency_planner_requests_and_coordinate_upper_bound_only"), source_material_consulted = FALSE, private_randomness_consulted = FALSE, runtime_failure_consulted = FALSE,
        automatic_fallback = FALSE, utility_optimality_claimed = FALSE)
    if (!identical(.dsvert_dp_analysis_client_canonical_value_v1(selection_certificate), .dsvert_dp_analysis_client_canonical_value_v1(expected_certificate)) || !identical(candidate$accuracy_certificate_sha256, certificate_hash) || !identical(selection$summary$selection_certificate_sha256,
        .dsvert_dp_frequency_client_hash_v1("dsVert/frequency/backend-selection/certificate/v1|", selection_certificate)))
        stop("Invalid Frequency backend selection.", call. = FALSE)
    fraction <- function(numerator, denominator) list(numerator = as.character(numerator), denominator = as.character(denominator))
    allocated <- .dsvert_dp_analysis_frequency_decimal_fraction_v1(request$delta)
    implementation <- if (isTRUE(profile$gaussian))
        fraction(plan$per_peer_implementation_delta_numerator, plan$per_peer_implementation_delta_denominator)
    else fraction(plan$implementation_delta_numerator, plan$implementation_delta_denominator)
    core <- if (isTRUE(profile$gaussian))
        fraction(plan$core_delta_numerator, plan$core_delta_denominator)
    else fraction("0", "1")
    maximum_noise <- plan[[if (isTRUE(profile$gaussian))
        "maximum_noise_magnitude_per_peer"
    else "maximum_noise_magnitude"]]
    upper <- format(config$coordinate_upper_bound, scientific = FALSE, trim = TRUE)
    result <- list(version = "dsvert-frequency-plan-summary-v1", physical_plan_version = profile$plan, full_plan_sha256 = plan_hash, planner_request_sha256 = .dsvert_dp_frequency_client_hash_v1(profile$request_domain,
        request), coordinate_order_sha256 = .dsvert_dp_analysis_frequency_coordinate_order_sha256_v1(config$factor_domain$levels), d = config$factor_domain$dimension, chunk_coordinates = min(profile$max_chunk_coordinates,
        config$factor_domain$dimension), allocated_delta = fraction(allocated$numerator, allocated$denominator), core_delta = core, implementation_delta = implementation, maximum_noise_per_peer = maximum_noise,
        no_wrap_sha256 = .dsvert_dp_frequency_client_hash_v1("dsVert/frequency/ring128-no-wrap/v1|", list(version = "dsvert-frequency-ring128-no-wrap-v1", coordinate_upper_bound = upper,
            maximum_noise_per_peer = maximum_noise, maximum_noise_release = as.character(2 * openssl::bignum(maximum_noise)))), profile_sha256 = .dsvert_dp_frequency_client_hash_v1("dsVert/frequency/physical-profile/v1|",
            profile), backend_selection = selection$summary)
    sensitivity <- list(value = if (identical(config$privacy$adjacency, "replace_one_fixed_cohort")) {
        if (isTRUE(profile$gaussian)) sqrt(2) else 2
    } else 1)
    tryCatch(.dsvert_dp_analysis_frequency_plan_validate_v1(result, profile, config$privacy, sensitivity, config$factor_domain$dimension, config$coordinate_upper_bound, config$calibration),
        error = function(error) stop("Invalid Frequency backend selection.", call. = FALSE))
    result
}
.dsvert_dp_frequency_client_config_v1 <- function(value) {
    fields <- c("version", "domain", "cohort_id", "dataset_id", "dataset_version", "privacy_unit_column", "alignment_purpose", "source_owner", "source_binding_id", "factor_domain",
        "factor_entry_sha256", "coordinate_upper_bound", "max_records_per_unit", "repeated_record_policy", "overflow_policy", "missingness_policy", "privacy", "calibration", "peer_pins",
        "backend_build_sha256", "transport_chunk_coordinates", "backend_selection")
    if (!.dsvert_dp_frequency_client_object_v1(value, fields) || !identical(value$version, .DSVERT_CLIENT_DP_FREQUENCY_CONFIG_VERSION)) {
        stop("Invalid Frequency configuration.", call. = FALSE)
    }
    pins <- .dsvert_dp_frequency_client_pins_v1(value$peer_pins)
    owner <- value$source_owner
    factor <- .dsvert_dp_frequency_client_factor_v1(value$factor_domain)
    scalar <- function(x, lower, upper, integer = FALSE) is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= lower && x <= upper && (!integer || x == floor(x))
    ids <- c("domain", "cohort_id", "dataset_id", "dataset_version", "alignment_purpose")
    privacy <- value$privacy
    calibration <- value$calibration
    source <- list(alignment_purpose = value$alignment_purpose, dataset_id = value$dataset_id, dataset_version = value$dataset_version, id_column = value$privacy_unit_column)
    source_binding <- paste0("source_", digest::digest(.dsvert_dp_frequency_client_wire_json_v1(source), "sha256", serialize = FALSE))
    valid <- all(vapply(value[ids], .dsvert_dp_analysis_client_scalar_id, logical(1L))) && is.character(value$privacy_unit_column) && grepl("^[A-Za-z._][A-Za-z0-9._]{0,127}$", value$privacy_unit_column) &&
        identical(value$source_binding_id, source_binding) && .dsvert_dp_frequency_client_object_v1(owner, c("peer_name", "identity_pk")) && .dsvert_dp_frequency_client_object_v1(privacy,
        c("adjacency", "epsilon", "delta")) && .dsvert_dp_frequency_client_object_v1(calibration, "implementation_delta") && privacy$adjacency %in% c("add_remove_patient", "replace_one_fixed_cohort") &&
        scalar(privacy$epsilon, .Machine$double.xmin, 8) && scalar(privacy$delta, .Machine$double.xmin, 1 - .Machine$double.eps) && scalar(calibration$implementation_delta, .Machine$double.xmin,
        privacy$delta) && scalar(value$coordinate_upper_bound, 1, 1e+06, TRUE) && identical(as.numeric(value$max_records_per_unit), 1) && identical(value$repeated_record_policy, "psi_v4_first_eligible_source_record_per_privacy_unit_v1") &&
        identical(value$overflow_policy, "clip_to_psi_v4_first_eligible_source_record_v1") && identical(value$missingness_policy, "missing_or_out_of_domain_rows_are_ignored") && .dsvert_dp_analysis_frequency_hex_v1(value$backend_build_sha256) &&
        scalar(value$transport_chunk_coordinates, 1, .Machine$integer.max, TRUE) && identical(value$factor_entry_sha256, .dsvert_dp_frequency_client_factor_hash_v1(factor))
    owner_peer <- tryCatch(.dsvert_dp_frequency_client_peer_v1(owner$peer_name), error = function(error) NULL)
    owner_pk <- tryCatch(.dsvert_dp_analysis_client_identity_pk(owner$identity_pk), error = function(error) NULL)
    if (!isTRUE(valid) || is.null(owner_peer) || is.null(owner_pk) || !owner_peer %in% names(pins) || !identical(owner_pk, unname(pins[[owner_peer]])))
        stop("Invalid Frequency configuration.", call. = FALSE)
    result <- value
    result$source_owner <- list(peer_name = owner_peer, identity_pk = owner_pk)
    result$factor_domain <- factor
    result$coordinate_upper_bound <- as.numeric(value$coordinate_upper_bound)
    result$max_records_per_unit <- 1
    result$privacy <- list(adjacency = privacy$adjacency, epsilon = as.numeric(privacy$epsilon), delta = as.numeric(privacy$delta))
    result$calibration <- list(implementation_delta = as.numeric(calibration$implementation_delta))
    result$peer_pins <- pins
    result$transport_chunk_coordinates <- as.numeric(value$transport_chunk_coordinates)
    plan <- .dsvert_dp_frequency_client_plan_v1(result)
    if (!identical(result$transport_chunk_coordinates, as.numeric(plan$chunk_coordinates)))
        stop("Invalid Frequency transport geometry.", call. = FALSE)
    result
}
.dsvert_dp_frequency_client_config_hash_v1 <- function(value) {
    value <- .dsvert_dp_frequency_client_config_v1(value)
    value$peer_pins <- as.list(value$peer_pins)
    .dsvert_dp_frequency_client_hash_v1("dsVert/dp-frequency/config/v1|", value)
}
.dsvert_dp_frequency_client_claim_hash_v1 <- function(value, config) {
    value <- .dsvert_dp_frequency_client_claim_v1(value, config$source_owner$peer_name, config$peer_pins)
    .dsvert_dp_frequency_client_hash_v1("dsVert/dp-frequency/source-claim/v1|", value)
}
.dsvert_dp_frequency_client_receipt_v1 <- function(value, config, claim) {
    fields <- c("version", "peer_name", "peer_identity_pk", "config_sha256", "source_claim_sha256", "psi_run_sha256", "snapshot_commitment", "signature")
    if (!.dsvert_dp_frequency_client_object_v1(value, fields) || !identical(value$version, .DSVERT_CLIENT_DP_FREQUENCY_RECEIPT_VERSION))
        stop("Invalid signed Frequency receipt fields.", call. = FALSE)
    peer <- .dsvert_dp_frequency_client_peer_v1(value$peer_name)
    identity <- tryCatch(.dsvert_dp_analysis_client_identity_pk(value$peer_identity_pk), error = function(error) NULL)
    unsigned <- value[setdiff(names(value), "signature")]
    for (field in c("config_sha256", "source_claim_sha256", "psi_run_sha256", "snapshot_commitment")) {
        .dsvert_dp_frequency_client_hex_v1(unsigned[[field]], paste("receipt", field))
    }
    if (is.null(identity) || !peer %in% names(config$peer_pins) || !identical(identity, unname(config$peer_pins[[peer]])) || !identical(unsigned$config_sha256, .dsvert_dp_frequency_client_config_hash_v1(config)) ||
        !identical(unsigned$source_claim_sha256, .dsvert_dp_frequency_client_claim_hash_v1(claim, config)) || !identical(unsigned$psi_run_sha256, claim$psi_run_sha256))
        stop("Frequency receipt disagrees with its signed inputs.", call. = FALSE)
    unsigned$peer_name <- peer
    unsigned$peer_identity_pk <- identity
    unsigned <- .dsvert_dp_analysis_client_canonical_value_v1(unsigned)
    .dsvert_dp_frequency_client_verify_v1(charToRaw(paste0("dsVert/dp-frequency/receipt/v1|", .dsvert_joint_dp_client_json(unsigned))), identity, value$signature, "receipt")
    .dsvert_dp_analysis_client_canonical_value_v1(c(unsigned, list(signature = value$signature)))
}
.dsvert_dp_frequency_client_contract_v1 <- function(config, receipts) {
    plan <- .dsvert_dp_frequency_client_plan_v1(config)
    snapshots <- stats::setNames(lapply(receipts, function(receipt) list(version = .DSVERT_DP_ANALYSIS_SNAPSHOT_VERSION, dataset_id = config$dataset_id, dataset_version = config$dataset_version,
        snapshot_commitment = receipt$snapshot_commitment)), vapply(receipts, `[[`, character(1L), "peer_identity_pk"))
    snapshots <- snapshots[order(names(snapshots), method = "radix")]
    source <- config$source_owner$identity_pk
    secondary <- sort(setdiff(names(snapshots), source), method = "radix")[[1L]]
    roles <- list(version = "dsvert-frequency-noise-authority-roles-v1", role_order = list("source_owner", "secondary_noise_authority"), authority_ids = list(source, secondary))
    profile <- .dsvert_dp_analysis_frequency_profile_v1(config$backend_selection$summary$selected_primitive)
    sensitivity <- if (identical(config$privacy$adjacency, "replace_one_fixed_cohort")) {
        if (isTRUE(profile$gaussian))
            sqrt(2)
        else 2
    }
    else 1
    semantic <- list(version = .DSVERT_DP_ANALYSIS_FREQUENCY_SEMANTIC_VERSION, domain = config$domain, cohort_id = config$cohort_id, owner_snapshots = snapshots, noise_authority_roles = roles,
        analysis = list(primitive = config$backend_selection$summary$selected_primitive, formula = NULL, effective_arguments = list(version = "dsvert-fixed-domain-categorical-frequency-v2",
            statistic = "aligned_fixed_domain_categorical_frequency", source_owner = source, dataset_id = config$dataset_id, dataset_version = config$dataset_version, variable_id = config$factor_domain$variable_id,
            levels = config$factor_domain$levels, dimension = config$factor_domain$dimension, repeated_record_policy = config$repeated_record_policy, missingness_policy = config$missingness_policy,
            coordinate_bounds = list(lower = 0, upper = config$coordinate_upper_bound), sampler_plan = plan)), privacy = list(version = "dsvert-per-analysis-dp-v1", adjacency = config$privacy$adjacency,
            privacy_unit = "patient", contribution = list(version = "dsvert-contribution-policy-v1", max_records_per_unit = 1, overflow_policy = config$overflow_policy, constraints = list(version = "dsvert-contribution-constraints-v1",
                policy_sha256 = .dsvert_dp_analysis_frequency_contribution_sha256_v1())), mechanism = list(family = profile$mechanism_family, version = profile$mechanism, sensitivity = list(version = "dsvert-sensitivity-v1",
                norm = profile$sensitivity_norm, value = sensitivity), calibration = list(version = "dsvert-calibration-v1", sampler = profile$sampler, implementation_delta = config$calibration$implementation_delta),
                randomness = list(version = "dsvert-randomness-plan-v1", lanes = list(final_noise = list(version = "dsvert-randomness-lane-v1", purpose = "privatize_final_vector", primitive = profile$sampler,
                  coordinates = config$factor_domain$dimension)))), epsilon = config$privacy$epsilon, delta = config$privacy$delta), numeric = list(version = "dsvert-numeric-semantics-v1",
            value_bits = 128, fractional_bits = 0, rounding = "toward_zero", overflow = "reject", output_encoding = "twos_complement_integer_v1"), public_shape = list(counts = config$factor_domain$dimension))
    contract <- list(version = .DSVERT_DP_ANALYSIS_CONTRACT_VERSION, artifact_key = .dsvert_dp_analysis_artifact_key_v1(semantic), semantic = semantic, execution = list(version = .DSVERT_DP_ANALYSIS_EXECUTION_VERSION,
        peer_pins = as.list(config$peer_pins), backend = list(kernel = config$backend_selection$summary$selected_primitive, ring = "ring128", build_sha256 = config$backend_build_sha256),
        transport = list(chunk_coordinates = config$transport_chunk_coordinates)))
    .dsvert_dp_analysis_contract_validate_v1(contract)
}
.dsvert_dp_frequency_client_worker_v1 <- function(contract, config, binding) {
    plan <- contract$semantic$analysis$effective_arguments$sampler_plan
    if (!identical(plan, .dsvert_dp_analysis_client_canonical_value_v1(.dsvert_dp_frequency_client_plan_v1(config))))
        stop("Frequency worker contract and configuration disagree.", call. = FALSE)
    primitive <- contract$semantic$analysis$primitive
    profile <- .dsvert_dp_analysis_frequency_profile_v1(primitive)
    roles <- binding$binding$authority_roles
    tokens <- lapply(roles, .dsvert_exact_gc_identity_peer_id)
    purpose <- if (isTRUE(profile$gaussian))
        "dyadic-discrete-gaussian-tv-bounded-v2"
    else "convolution"
    d <- as.integer(plan$d)
    chunk <- as.integer(plan$chunk_coordinates)
    raw_bound <- list(lower = "0", upper = format(config$coordinate_upper_bound, scientific = FALSE, trim = TRUE), scale = 0L)
    release <- list(version = "dsvert-dp-frequency-worker-release-v1", artifact_key = contract$artifact_key, primitive = primitive, selected_plan_sha256 = plan$full_plan_sha256, coordinate_order_sha256 = plan$coordinate_order_sha256,
        d = d, chunk_coordinates = chunk, raw_bound = raw_bound, authority_roles = roles)
    release_hash <- .dsvert_dp_frequency_client_hash_v1("dsVert/dp-frequency/worker-release/v1|", release)
    transcript_hash <- .dsvert_dp_frequency_client_hash_v1("dsVert/dp-frequency/worker-transcript/v1|", list(version = "dsvert-dp-frequency-worker-transcript-v1", release_contract_hash = release_hash,
        analysis_binding_sha256 = binding$sha256, authority_tokens = tokens))
    raw_transcript <- as.raw(strtoi(substring(transcript_hash, seq.int(1L, 63L, 2L), seq.int(2L, 64L, 2L)), 16L))
    contexts <- lapply(tokens, function(token) digest::digest(c(charToRaw("dsvert-joint-dp-private-seed-commitment-v2"), as.raw(0L), raw_transcript, as.raw(0L), charToRaw(purpose),
        as.raw(0L), charToRaw(token)), "sha256", serialize = FALSE))
    list(version = "dsvert-dp-frequency-worker-static-v1", selected_primitive = primitive, selected_profile = profile, selected_request = .dsvert_dp_analysis_client_canonical_value_v1(config$backend_selection$selected_request),
        selected_plan = .dsvert_dp_analysis_client_canonical_value_v1(config$backend_selection$selected_plan), selected_plan_sha256 = plan$full_plan_sha256, ring_bits = 128L, frac_bits = 0L, output_lattice_bits = as.integer(profile$output_lattice_bits),
        d = d, chunk_coordinates = chunk, chunk_count = as.integer(ceiling(d/chunk)), raw_bound = raw_bound, authority_roles = roles, authority_tokens = tokens, release_contract_hash = release_hash,
        transcript_hash = transcript_hash, commitment_purpose = purpose, commitment_contexts = contexts, source_share_policy = list(source_owner = "private_frequency_vector_ring128_v1",
            secondary_noise_authority = "zero_vector_ring128_v1"))
}
.dsvert_dp_frequency_client_compile_v1 <- function(envelopes, server_names, claim) {
    if (!is.list(envelopes) || is.null(names(envelopes)) || anyNA(names(envelopes)) || anyDuplicated(names(envelopes)) || !setequal(names(envelopes), server_names))
        stop("Frequency compile results do not cover the complete federation.", call. = FALSE)
    envelopes <- envelopes[server_names]
    configs <- lapply(envelopes, function(envelope) {
        if (!.dsvert_dp_frequency_client_object_v1(envelope, c("config", "receipt")))
            stop("A peer returned an invalid Frequency compile envelope.", call. = FALSE)
        .dsvert_dp_frequency_client_config_v1(envelope$config)
    })
    config <- configs[[1L]]
    if (!all(vapply(configs, identical, logical(1L), config)) || !setequal(names(config$peer_pins), server_names))
        stop("Frequency peers disagree on one full-K configuration.", call. = FALSE)
    claim <- .dsvert_dp_frequency_client_claim_v1(claim, config$source_owner$peer_name, config$peer_pins)
    claim_binding <- list(dataset_id = claim$dataset_id, dataset_version = claim$dataset_version, privacy_unit_column = claim$privacy_unit_column, alignment_purpose = claim$alignment_purpose,
        source_binding_id = claim$source_binding_id, factor_entry_sha256 = claim$factor_entry_sha256)
    config_binding <- config[names(claim_binding)]
    if (!identical(config$source_owner, list(peer_name = claim$source_peer_name, identity_pk = claim$source_identity_pk)) || !identical(config_binding, claim_binding) || !identical(.dsvert_dp_analysis_client_canonical_value_v1(config$factor_domain),
        claim$factor_entry))
        stop("Frequency configuration disagrees with its signed Claim.", call. = FALSE)
    receipts <- Map(function(envelope, peer) {
        receipt <- .dsvert_dp_frequency_client_receipt_v1(envelope$receipt, config, claim)
        if (!identical(receipt$peer_name, peer))
            stop("Frequency receipt does not match its connected peer.", call. = FALSE)
        receipt
    }, envelopes, server_names)
    names(receipts) <- vapply(receipts, `[[`, character(1L), "peer_name")
    if (anyDuplicated(names(receipts)) || !setequal(names(receipts), names(config$peer_pins)))
        stop("Invalid Frequency receipt coverage.", call. = FALSE)
    receipts <- receipts[names(config$peer_pins)]
    common <- function(field) length(unique(vapply(receipts, `[[`, character(1L), field))) == 1L
    if (!all(vapply(c("config_sha256", "source_claim_sha256", "psi_run_sha256"), common, logical(1L))))
        stop("Frequency receipts lack exact configuration/Claim/PSI consensus.", call. = FALSE)
    contract <- .dsvert_dp_frequency_client_contract_v1(config, receipts)
    binding <- .dsvert_exact_gc_frequency_analysis_binding(contract)
    worker <- .dsvert_dp_frequency_client_worker_v1(contract, config, binding)
    role_ids <- unlist(binding$binding$authority_roles, use.names = TRUE)
    authorities <- names(config$peer_pins)[match(role_ids, unname(config$peer_pins))]
    if (length(authorities) != 2L || anyNA(authorities) || anyDuplicated(authorities))
        stop("Frequency authorities do not match the full-K pinset.", call. = FALSE)
    list(claim = claim, config = config, receipts = receipts, contract = contract, binding = binding, worker_static = worker, authorities = unname(authorities), config_sha256 = .dsvert_dp_frequency_client_config_hash_v1(config),
        source_claim_sha256 = receipts[[1L]]$source_claim_sha256, receipt_set_sha256 = .dsvert_dp_frequency_client_hash_v1("dsVert/dp-frequency/receipt-set/v1|", receipts), psi_run_sha256 = receipts[[1L]]$psi_run_sha256,
        contract_sha256 = .dsvert_dp_frequency_client_hash_v1("dsVert/dp-frequency/compiled-contract/v1|", contract), worker_static_sha256 = .dsvert_dp_frequency_client_hash_v1("dsVert/dp-frequency/worker-static/v1|",
            worker))
}
.dsvert_dp_frequency_client_public_auth_v1 <- function(value, peer, session_id, compiled) {
    fields <- c("version", "session_id", "artifact_key", "config_sha256", "source_claim_sha256", "receipt_set_sha256", "psi_run_sha256", "contract_sha256", "analysis_binding_sha256",
        "worker_static_sha256", "local_authority", "commitment_context", "seed_commitment", "authorization_sha256", "signature")
    if (!.dsvert_dp_frequency_client_object_v1(value, fields) || !identical(value$version, .DSVERT_CLIENT_DP_FREQUENCY_PUBLIC_AUTH_VERSION) || !identical(value$session_id, session_id))
        stop("A Frequency authority returned an invalid public authorization.", call. = FALSE)
    peer <- .dsvert_dp_frequency_client_peer_v1(peer)
    roles <- compiled$binding$binding$authority_roles
    role <- names(roles)[match(unname(compiled$config$peer_pins[[peer]]), unlist(roles, use.names = FALSE))]
    local <- .dsvert_dp_analysis_client_canonical_value_v1(list(peer_name = peer, identity_pk = unname(compiled$config$peer_pins[[peer]]), role = unname(role)))
    common <- c(artifact_key = compiled$contract$artifact_key, config_sha256 = compiled$config_sha256, source_claim_sha256 = compiled$source_claim_sha256, receipt_set_sha256 = compiled$receipt_set_sha256,
        psi_run_sha256 = compiled$psi_run_sha256, contract_sha256 = compiled$contract_sha256, analysis_binding_sha256 = compiled$binding$sha256, worker_static_sha256 = compiled$worker_static_sha256)
    if (!identical(unlist(value[names(common)], use.names = TRUE), common) || !identical(value$local_authority, local) || !identical(value$commitment_context, compiled$worker_static$commitment_contexts[[role]]))
        stop("A Frequency authority returned a misbound public authorization.", call. = FALSE)
    for (field in c(names(common), "commitment_context", "seed_commitment", "authorization_sha256")) {
        .dsvert_dp_frequency_client_hex_v1(value[[field]], paste("authorization", field))
    }
    private <- list(version = "dsvert-dp-frequency-session-authorization-v1", session_id = session_id, artifact_key = compiled$contract$artifact_key, config = compiled$config, config_sha256 = compiled$config_sha256,
        source_claim_sha256 = compiled$source_claim_sha256, receipt_peers = as.list(names(compiled$receipts)), receipt_set_sha256 = compiled$receipt_set_sha256, psi_run_sha256 = compiled$psi_run_sha256,
        contract = compiled$contract, contract_sha256 = compiled$contract_sha256, analysis_binding = compiled$binding$binding, analysis_binding_sha256 = compiled$binding$sha256, worker_static = compiled$worker_static,
        worker_static_sha256 = compiled$worker_static_sha256, local_authority = local)
    wire_private <- private
    wire_private$config$peer_pins <- as.list(wire_private$config$peer_pins)
    expected_hash <- .dsvert_dp_frequency_client_hash_v1("dsVert/dp-frequency/session-authorization/v1|", wire_private)
    if (!identical(value$authorization_sha256, expected_hash))
        stop("A Frequency authority returned a corrupt public authorization.", call. = FALSE)
    unsigned <- .dsvert_dp_analysis_client_canonical_value_v1(value[setdiff(names(value), "signature")])
    .dsvert_dp_frequency_client_verify_v1(charToRaw(paste0("dsVert/dp-frequency/public-authorization-signature/v1|", .dsvert_joint_dp_client_json(unsigned))), local$identity_pk, value$signature,
        "public authorization")
    .dsvert_dp_analysis_client_canonical_value_v1(value)
}
.dsvert_dp_frequency_client_authorizations_v1 <- function(values, session_id, compiled) {
    authorities <- compiled$authorities
    if (!is.list(values) || is.null(names(values)) || anyNA(names(values)) || anyDuplicated(names(values)) || !setequal(names(values), authorities))
        stop("Frequency public authorizations do not cover both authorities.", call. = FALSE)
    values <- values[authorities]
    verified <- Map(function(value, peer) {
        .dsvert_dp_frequency_client_public_auth_v1(value, peer, session_id, compiled)
    }, values, authorities)
    names(verified) <- authorities
    roles <- vapply(verified, function(value) value$local_authority$role, character(1L))
    commitments <- vapply(verified, `[[`, character(1L), "seed_commitment")
    if (!identical(unname(roles), c("source_owner", "secondary_noise_authority")) || anyDuplicated(commitments))
        stop("Frequency public authorizations have invalid role coverage.", call. = FALSE)
    verified
}
.dsvert_dp_frequency_prepare_v1 <- function(data_name, variable_name, source_owner, datasources, .aggregate = DSI::datashield.aggregate, .session_id = .dsvert_uuid4, .setup_frequency = .dsvert_setup_frequency_transport) {
    if (!is.character(data_name) || length(data_name) != 1L || is.na(data_name) || !grepl("^[A-Za-z._][A-Za-z0-9._]*$", data_name) || is.null(source_owner))
        stop("Frequency requires an explicit source owner before DSI.", call. = FALSE)
    variable_name <- .dsvert_dp_frequency_client_text_v1(variable_name, "factor variable")
    datasources <- .dsvert_dp_datasources(datasources)
    if (length(datasources) < 2L)
        stop("Frequency requires the complete named K>=2 datasource set.", call. = FALSE)
    source_owner <- .dsvert_dp_server(source_owner, datasources)
    server_names <- names(datasources)
    claim_values <- .dsvert_aggregate_strict(datasources[source_owner], call(name = "dsvertDPFrequencyClaimDS", data_name = data_name, variable_name = variable_name), operation = "Frequency source Claim",
        .aggregate = .aggregate)
    claim <- .dsvert_dp_frequency_client_claim_v1(claim_values[[source_owner]], source_owner)
    if (!identical(claim$factor_entry$variable_name, variable_name))
        stop("Frequency Claim names a different factor variable.", call. = FALSE)
    claim_json <- .dsvert_joint_dp_client_json(claim)
    if (nchar(claim_json, type = "bytes") > .DSVERT_CLIENT_DP_FREQUENCY_CONFIG_MAX_BYTES)
        stop("Frequency Claim exceeds its bounded wire contract.", call. = FALSE)
    compiled_values <- .dsvert_aggregate_strict(datasources, call(name = "dsvertDPFrequencyCompileDS", data_name = data_name, source_claim_json = .dsvert_dsi_text_encode(claim_json,
        "Frequency source Claim")), operation = "Frequency K-wide signed compilation", .aggregate = .aggregate)
    compiled <- .dsvert_dp_frequency_client_compile_v1(compiled_values, server_names, claim)
    session_id <- .session_id()
    if (!is.character(session_id) || length(session_id) != 1L || is.na(session_id) || !grepl(paste0("^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-", "[89ab][0-9a-f]{3}-[0-9a-f]{12}$"), session_id))
        stop("Invalid Frequency client session identifier.", call. = FALSE)
    wire_config <- compiled$config
    wire_config$peer_pins <- as.list(wire_config$peer_pins)
    config_json <- .dsvert_joint_dp_client_json(wire_config)
    receipts_json <- .dsvert_joint_dp_client_json(unname(compiled$receipts))
    if (nchar(config_json, type = "bytes") > .DSVERT_CLIENT_DP_FREQUENCY_CONFIG_MAX_BYTES || nchar(receipts_json, type = "bytes") > .DSVERT_CLIENT_DP_FREQUENCY_RECEIPTS_MAX_BYTES)
        stop("Frequency authorization input exceeds its bounded wire contract.", call. = FALSE)
    authority_conns <- datasources[compiled$authorities]
    completed <- FALSE
    on.exit(if (!completed) .dsvert_cleanup_best_effort(authority_conns,
        call(name = "dsvertDPFrequencyCleanupDS", session_id = session_id), .aggregate = .aggregate), add = TRUE)
    authorized <- .dsvert_aggregate_strict(authority_conns, call(name = "dsvertDPFrequencyAuthorizeDS", config_json = .dsvert_dsi_text_encode(config_json, "Frequency authorization configuration"),
        receipts_json = .dsvert_dsi_text_encode(receipts_json, "Frequency authorization receipts"), source_claim_json = .dsvert_dsi_text_encode(claim_json, "Frequency source Claim"),
        session_id = session_id), operation = "Frequency two-authority authorization", .aggregate = .aggregate)
    authorized <- .dsvert_dp_frequency_client_authorizations_v1(authorized, session_id, compiled)
    transport <- .setup_frequency(datasources = datasources, server_names = server_names, servers = match(compiled$authorities, server_names), session_id = session_id, frequency_analysis = compiled,
        .aggregate = .aggregate)
    result <- list(version = .DSVERT_CLIENT_DP_FREQUENCY_PREPARED_VERSION, session_id = session_id, source_owner = source_owner, claim = claim, config = compiled$config, receipts = compiled$receipts,
        contract = compiled$contract, worker_static = compiled$worker_static, authorities = compiled$authorities, authorizations = authorized, transport = transport)
    completed <- TRUE
    result
}
