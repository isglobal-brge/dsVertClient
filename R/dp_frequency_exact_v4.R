# Public, exact arithmetic checks for the unbounded v4 Laplace contract.
.DSVERT_JOINT_DP_PURE_CERTIFICATE_FIELDS <- c(
  "guarantee", "randomness", "noise_support", "noise_commitment",
  "representability_bound", "wrap_bound_certified", "wrap_bound",
  "wrap_bound_event", "admitted_ranges")

.dsvert_joint_dp_exact_exponential_bound <- function(numerator, denominator,
                                                    multiplicity) {
  # exp(3)>10; decimal exponents stay integers, never floating probabilities.
  power <- 0L
  factor <- 1
  while (factor < multiplicity) {
    power <- power + 1L
    factor <- factor * 10
  }
  exponent <- numerator %/% (3 * denominator)
  if (exponent <= power) return("1")
  paste0("1e-", as.character(exponent - power))
}

.dsvert_joint_dp_exact_decimal_fraction <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      nchar(value) > 512L ||
      !grepl("^[0-9]+(\\.[0-9]+)?([eE][+-]?[0-9]+)?$", value)) {
    stop("Invalid exact Laplace decimal.", call. = FALSE)
  }
  parts <- strsplit(tolower(value), "e", fixed = TRUE)[[1L]]
  exponent <- if (length(parts) == 2L) as.numeric(parts[[2L]]) else 0
  if (!is.finite(exponent) || abs(exponent) > 10000) {
    stop("Invalid exact Laplace decimal exponent.", call. = FALSE)
  }
  normalized <- paste0(parts[[1L]], "e", if (exponent >= 0) "+" else "-",
                       abs(exponent))
  .dsvert_dp_analysis_frequency_decimal_fraction_v1(normalized)
}

.dsvert_joint_dp_pure_request <- function(request) {
  fail <- function() stop("Invalid exact Laplace admitted ranges.", call. = FALSE)
  if (!is.list(request) || !setequal(names(request), c(
      "epsilon", "delta", "sensitivity_steps", "total_coordinate_count")) ||
      !is.character(request$sensitivity_steps) ||
      length(request$sensitivity_steps) != 1L ||
      is.na(request$sensitivity_steps) ||
      nchar(request$sensitivity_steps) > 512L ||
      !grepl("^[1-9][0-9]*$", request$sensitivity_steps) ||
      !is.numeric(request$total_coordinate_count) ||
      length(request$total_coordinate_count) != 1L ||
      is.na(request$total_coordinate_count) ||
      request$total_coordinate_count < 1 ||
      request$total_coordinate_count > 1000000 ||
      request$total_coordinate_count != floor(request$total_coordinate_count)) fail()
  epsilon <- tryCatch(.dsvert_joint_dp_exact_decimal_fraction(
    request$epsilon), error = function(e) NULL)
  delta <- tryCatch(.dsvert_joint_dp_exact_decimal_fraction(
    request$delta), error = function(e) NULL)
  if (is.null(epsilon) || is.null(delta) || epsilon$numerator == 0 ||
      epsilon$numerator > 10000 * epsilon$denominator ||
      delta$numerator >= delta$denominator) fail()
  sensitivity <- openssl::bignum(request$sensitivity_steps)
  denominator <- epsilon$denominator * sensitivity
  if (epsilon$numerator * (openssl::bignum(2)^100) < denominator) fail()
  list(numerator = epsilon$numerator, denominator = denominator,
       epsilon = epsilon, dimension = request$total_coordinate_count)
}

.dsvert_joint_dp_pure_certificate <- function(request) {
  rate <- .dsvert_joint_dp_pure_request(request)
  bound <- .dsvert_joint_dp_exact_exponential_bound(
    rate$numerator * (openssl::bignum(2)^127), rate$denominator,
    2 * rate$dimension)
  list(guarantee = "pure-dp-under-ideal-bits",
    randomness = "keyed-stream-computational",
    noise_support = "unbounded_integer", noise_commitment = "modulo_2^128",
    representability_bound = bound, wrap_bound_certified = TRUE,
    wrap_bound = bound,
    wrap_bound_event = "at_least_one_peer_draw_outside_signed_Ring128",
    admitted_ranges = list(epsilon_minimum_exclusive = "0",
      epsilon_maximum = "10000", sensitivity_steps = "positive_integer",
      epsilon_over_sensitivity_minimum = "2^-100",
      total_coordinate_count_minimum = 1L,
      total_coordinate_count_maximum = 1000000L, ring_bits = 128L))
}

.dsvert_joint_dp_pure_plan_validate <- function(plan, request) {
  fail <- function() stop("Invalid exact Laplace plan certificate.", call. = FALSE)
  expected <- .dsvert_joint_dp_pure_certificate(request)
  if (!is.list(plan) || anyDuplicated(names(plan)) ||
      !all(names(expected) %in% names(plan)) ||
      !identical(.dsvert_dp_analysis_client_canonical_value_v1(plan[names(expected)]),
                 .dsvert_dp_analysis_client_canonical_value_v1(expected)) ||
      !identical(plan$version, "dsvert-joint-dp-vector-independent-full-draw-convolution-plan-v4") ||
      !identical(plan$sampler, "hkdf-sha256-chacha20-independent-full-draw-exact-geometric-v4") ||
      !identical(plan$sensitivity_steps, request$sensitivity_steps) ||
      !identical(as.numeric(plan$total_coordinate_count),
                 as.numeric(request$total_coordinate_count)) ||
      !identical(plan$maximum_noise_magnitude, "unbounded") ||
      !identical(plan$complete_epsilon_per_peer, TRUE) ||
      !identical(plan$epsilon_divided_by_peer_count, FALSE) ||
      !identical(plan$release_implementation_delta_aggregation, "max_per_peer_not_sum") ||
      !identical(as.numeric(plan$stop_bits), 0) ||
      !identical(plan$stop_numerator, "0") ||
      !identical(as.numeric(plan$uniform_bits), 0) ||
      !identical(as.numeric(plan$binary_geometric_bits), 0) ||
      length(plan$bernoulli_thresholds) != 0L ||
      !identical(as.numeric(plan$independent_noise_peer_count), 2) ||
      !identical(as.numeric(plan$geometric_variables_per_peer_per_coordinate), 2) ||
      !identical(as.numeric(plan$geometric_variables_total_per_coordinate), 4)) fail()
  for (prefix in c("one_geometric_tv", "tail_upper", "rounding_upper",
      "implementation_delta", "per_peer_implementation_delta",
      "two_peer_ideal_transfer_delta")) {
    if (!identical(plan[[paste0(prefix, "_numerator")]], "0") ||
        !identical(plan[[paste0(prefix, "_denominator")]], "1")) fail()
  }
  epsilon <- .dsvert_joint_dp_pure_request(request)$epsilon
  if (!identical(plan$epsilon_effective_upper_numerator,
                 as.character(epsilon$numerator)) ||
      !identical(plan$epsilon_effective_upper_denominator,
                 as.character(epsilon$denominator))) fail()
  invisible(.dsvert_dp_analysis_client_canonical_value_v1(plan))
}

.dsvert_joint_dp_pure_sum_certificate <- function(plan, maximum) {
  boundary <- openssl::bignum(2)^127
  if (maximum >= boundary) stop("Invalid exact Laplace source bound.", call. = FALSE)
  threshold <- (boundary - maximum) %/% openssl::bignum(2)
  numerator <- openssl::bignum(plan$epsilon_effective_upper_numerator)
  denominator <- openssl::bignum(plan$epsilon_effective_upper_denominator) *
    openssl::bignum(plan$sensitivity_steps)
  list(sum_wrap_bound = .dsvert_joint_dp_exact_exponential_bound(
    numerator * threshold, denominator, 2 * plan$total_coordinate_count),
    sum_wrap_threshold = as.character(threshold))
}

.dsvert_joint_dp_pure_output_certificate <- function(plan, bounds, shifts) {
  maximum <- openssl::bignum(0)
  for (index in seq_along(bounds)) {
    value <- openssl::bignum(bounds[[index]]) * (openssl::bignum(2)^shifts[[index]])
    if (value > maximum) maximum <- value
  }
  c(plan[.DSVERT_JOINT_DP_PURE_CERTIFICATE_FIELDS],
    .dsvert_joint_dp_pure_sum_certificate(plan, maximum))
}

.dsvert_dp_analysis_frequency_reduce_v2 <- function(numerator, denominator) {
  left <- numerator
  right <- denominator
  while (!identical(as.character(right), "0")) {
    remainder <- left %% right
    left <- right
    right <- remainder
  }
  list(numerator = as.character(numerator %/% left),
       denominator = as.character(denominator %/% left))
}


# Frequency selector generation 2: exact modular Laplace and optional Gaussian.
# Generation 1 readers remain in dpAnalysisContract.R for historical artifacts.
.dsvert_dp_analysis_frequency_plan_fields_v1 <- function() list(
    convolution = c(
      "version", "sampler", "stop_bits", "stop_numerator", "uniform_bits",
      "binary_geometric_bits", "bernoulli_thresholds", "sensitivity_steps",
      "total_coordinate_count", "epsilon_effective_upper_numerator",
      "epsilon_effective_upper_denominator", "one_geometric_tv_numerator",
      "one_geometric_tv_denominator", "tail_upper_numerator",
      "tail_upper_denominator", "rounding_upper_numerator",
      "rounding_upper_denominator", "implementation_delta_numerator",
      "implementation_delta_denominator", "implementation_delta_bound",
      "maximum_noise_magnitude", "maximum_chunk_coordinates",
      "private_stream_bytes_per_coordinate", "accounting",
      "capability_available", "independent_noise_peer_count",
      "complete_epsilon_per_peer", "epsilon_divided_by_peer_count",
      "geometric_variables_per_peer_per_coordinate",
      "geometric_variables_total_per_coordinate",
      "per_peer_implementation_delta_numerator",
      "per_peer_implementation_delta_denominator",
      "per_peer_implementation_delta_bound",
      "release_implementation_delta_aggregation",
      "two_peer_ideal_transfer_delta_numerator",
      "two_peer_ideal_transfer_delta_denominator",
      "two_peer_ideal_transfer_delta_bound", "threat_model",
      "privacy_argument"),
    gaussian = c(
      "version", "mechanism", "sampler", "reference",
      "total_coordinate_count", "maximum_chunk_coordinates",
      "request_binding_sha256", "epsilon_numerator", "epsilon_denominator",
      "allocated_delta_numerator", "allocated_delta_denominator",
      "core_delta_numerator", "core_delta_denominator",
      "tail_delta_numerator", "tail_delta_denominator",
      "l2_sensitivity_numerator", "l2_sensitivity_denominator",
      "rho_numerator", "rho_denominator", "zcdp_log_upper_integer",
      "zcdp_conversion_exponent_numerator",
      "zcdp_conversion_exponent_denominator", "sigma_squared_numerator",
      "sigma_squared_denominator", "proposal_scale",
      "maximum_noise_magnitude_per_peer",
      "maximum_noise_magnitude_two_peers", "tail_proof_exponent_numerator",
      "tail_proof_exponent_denominator", "tail_proof_target_numerator",
      "tail_proof_target_denominator", "vector_tail_tv_upper_numerator",
      "vector_tail_tv_upper_denominator",
      "vector_sampler_tv_upper_numerator",
      "vector_sampler_tv_upper_denominator",
      "vector_total_tv_upper_numerator",
      "vector_total_tv_upper_denominator",
      "per_peer_implementation_delta_numerator",
      "per_peer_implementation_delta_denominator", "simultaneous_95_abs",
      "sampler_candidate_count", "sampler_random_bits_per_coordinate",
      "sampler_random_bytes_per_coordinate", "sampler_table_precision_bits",
      "sampler_magnitude_count", "sampler_search_steps",
      "sampler_full_scan_steps", "sampler_cdf_table_bytes",
      "accuracy_accounting", "accounting", "privacy_theorem",
      "independent_noise_peer_count", "complete_epsilon_per_peer",
      "epsilon_divided_by_peer_count", "release_delta_aggregation",
      "nominal_variance_multiplier", "nominal_standard_deviation_factor",
      "at_least_one_honest_noise_peer", "maximum_colluding_noise_peers",
      "adversary_view", "adversary_view_privacy_argument",
      "source_share_hiding_precondition", "exact_rational_sampler",
      "finite_support_transfer_charged", "fixed_work_sampler",
      "sampler_branches_on_protected_values",
      "sampler_branches_on_private_randomness", "host_constant_time_claim",
      "transcript_dp_claim", "logical_transcript_fixed_shape",
      "physical_timing_dp_claim", "observable_worker_shape",
      "capability_available", "unavailable_reason"))

.dsvert_dp_analysis_frequency_primitives_v3 <- function() list(
  convolution = "independent_full_global_draw_convolution_ring128_v4",
  gaussian = "independent_full_global_dyadic_discrete_gaussian_tv_bounded_ring128_v2")

.dsvert_dp_analysis_frequency_candidate_requests_v3 <- function(
    privacy, calibration, dimension) {
  values <- .dsvert_dp_analysis_frequency_candidate_requests_v2(
    privacy, calibration, dimension)
  .dsvert_joint_dp_pure_request(values$convolution)
  values
}

.dsvert_dp_analysis_frequency_policy_sha256_v3 <- function() {
  .dsvert_dp_analysis_frequency_hash_v1(
    "dsVert/frequency/backend-selection-policy/v3|", list(
      version = "dsvert-frequency-backend-selection-policy-v3",
      oracle_policy = "minimum_certified_simultaneous_95_abs_convolution_tie_v2",
      candidate_primitives = unname(.dsvert_dp_analysis_frequency_primitives_v3()),
      objective = "minimum_certified_simultaneous_95_abs",
      accuracy_event = "max_j_abs_error_gt_radius",
      tie_break = "convolution_laplace_v4_on_equal_certified_radius",
      input_scope = "public_adjacency_planner_requests_and_coordinate_upper_bound_only",
      source_material_consulted = FALSE, private_randomness_consulted = FALSE,
      runtime_failure_consulted = FALSE, automatic_fallback = FALSE,
      utility_optimality_claimed = FALSE))
}

.dsvert_dp_analysis_frequency_exact_radius_v4 <- function(request, bound) {
  rate <- .dsvert_joint_dp_pure_request(request)
  k <- 0L
  power <- 1
  while (power < 80 * rate$dimension) {
    k <- k + 1L
    power <- power * 10
  }
  numerator <- 3 * k * rate$denominator
  threshold <- (numerator + rate$numerator - 1) %/% rate$numerator
  upper <- openssl::bignum(format(bound, scientific = FALSE, trim = TRUE))
  radius <- 2 * threshold
  if (radius >= upper || upper + radius > (openssl::bignum(2)^127) - 1)
    radius <- upper
  as.character(radius)
}

.dsvert_dp_analysis_frequency_accuracy_validate_v3 <- function(
    certificate, plan, primitive, request, bound) {
  fail <- function() stop("Invalid Frequency accuracy certificate.", call. = FALSE)
  profile <- .dsvert_dp_analysis_frequency_profile_v1(primitive)
  fields <- c("primitive", "plan_sha256", "event", "method",
    "release_tv_upper_numerator", "release_tv_upper_denominator",
    "simultaneous_95_abs", "absolute_support")
  plan_fields <- .dsvert_dp_analysis_frequency_plan_fields_v1()[[
    if (isTRUE(profile$gaussian)) "gaussian" else "convolution"]]
  if (isTRUE(profile$exact)) plan_fields <- c(plan_fields, .DSVERT_JOINT_DP_PURE_CERTIFICATE_FIELDS)
  if (is.null(profile) || !.dsvert_dp_analysis_frequency_object_v1(plan, plan_fields) ||
      anyDuplicated(names(plan)) || !.dsvert_dp_analysis_frequency_object_v1(
      certificate, fields, list(primitive = primitive,
        event = "max_j_abs_error_gt_radius",
        plan_sha256 = .dsvert_dp_analysis_frequency_hash_v1(
          .DSVERT_DP_ANALYSIS_FREQUENCY_PLAN_DOMAIN_V1, plan))) ||
      !identical(plan$version, profile$plan) ||
      !identical(plan$sampler, profile$sampler) ||
      !identical(as.numeric(plan$total_coordinate_count),
                 as.numeric(request$total_coordinate_count)) ||
      !identical(as.numeric(plan$maximum_chunk_coordinates),
                 min(profile$max_chunk_coordinates, request$total_coordinate_count)) ||
      !identical(plan$capability_available, TRUE) ||
      !identical(as.numeric(plan$independent_noise_peer_count), 2) ||
      !identical(plan$complete_epsilon_per_peer, TRUE) ||
      !identical(plan$epsilon_divided_by_peer_count, FALSE)) fail()
  if (isTRUE(profile$exact)) {
    .dsvert_joint_dp_pure_plan_validate(plan, request)
    if (!identical(certificate$method,
        "exact_two_draw_exponential_union_modular_clamp_v1") ||
        !identical(certificate$release_tv_upper_numerator, "0") ||
        !identical(certificate$release_tv_upper_denominator, "1") ||
        !identical(certificate$absolute_support,
          format(bound, scientific = FALSE, trim = TRUE)) ||
        !identical(certificate$simultaneous_95_abs,
          .dsvert_dp_analysis_frequency_exact_radius_v4(request, bound))) fail()
  } else {
    peer <- .dsvert_dp_analysis_frequency_uint_v1(
      plan$maximum_noise_magnitude_per_peer, TRUE)
    support <- .dsvert_dp_analysis_frequency_uint_v1(
      plan$maximum_noise_magnitude_two_peers, TRUE)
    tv <- .dsvert_dp_analysis_frequency_fraction_v1(list(
      numerator = plan$vector_total_tv_upper_numerator,
      denominator = plan$vector_total_tv_upper_denominator), TRUE)
    released_tv <- .dsvert_dp_analysis_frequency_reduce_v2(
      2 * tv$numerator, tv$denominator)
    sensitivity <- .dsvert_dp_analysis_frequency_decimal_fraction_v1(
      request$l2_sensitivity_steps)
    if (!identical(plan$mechanism, "dyadic_discrete_gaussian_truncated_tv_bounded") ||
        !identical(plan$release_delta_aggregation, "max_per_peer_not_sum") ||
        !identical(plan$l2_sensitivity_numerator, as.character(sensitivity$numerator)) ||
        !identical(plan$l2_sensitivity_denominator, as.character(sensitivity$denominator)) ||
        !isTRUE(support == 2 * peer) || !isTRUE(40 * tv$numerator < tv$denominator) ||
        !identical(certificate$method, "gaussian_plan_v2_subgaussian_mgf_tv_transfer") ||
        !identical(certificate$release_tv_upper_numerator, released_tv$numerator) ||
        !identical(certificate$release_tv_upper_denominator, released_tv$denominator) ||
        !identical(certificate$absolute_support, as.character(support)) ||
        !identical(certificate$simultaneous_95_abs, plan$simultaneous_95_abs)) fail()
  }
  radius <- .dsvert_dp_analysis_frequency_uint_v1(certificate$simultaneous_95_abs)
  support <- .dsvert_dp_analysis_frequency_uint_v1(certificate$absolute_support, TRUE)
  if (radius > support || (!isTRUE(profile$exact) &&
      openssl::bignum(format(bound, scientific = FALSE, trim = TRUE)) + support >
        (openssl::bignum(2)^127) - 1)) fail()
  invisible(TRUE)
}

.dsvert_dp_analysis_frequency_selection_certificate_v3 <- function(
    primitive, plan_hash, radius) list(
  version = "dsvert-joint-dp-frequency-backend-selection-certificate-v2",
  policy = "minimum_certified_simultaneous_95_abs_convolution_tie_v2",
  objective = "minimum_certified_simultaneous_95_abs",
  selected_primitive = primitive, selected_plan_sha256 = plan_hash,
  selected_simultaneous_95_abs = radius,
  tie_break = "convolution_laplace_v4_on_equal_certified_radius",
  input_scope = "public_adjacency_planner_requests_and_coordinate_upper_bound_only",
  source_material_consulted = FALSE, private_randomness_consulted = FALSE,
  runtime_failure_consulted = FALSE, automatic_fallback = FALSE,
  utility_optimality_claimed = FALSE)

.dsvert_dp_analysis_frequency_selection_validate_v3 <- function(
    selection, selected_primitive, plan, privacy, calibration, dimension, bound) {
  fail <- function() stop("Invalid Frequency backend selection", call. = FALSE)
  fields <- c("version", "policy_sha256", "selection_certificate_sha256",
    "objective", "tie_break", "candidates", "selected_primitive", "selected_simultaneous_95_abs")
  if (!.dsvert_dp_analysis_frequency_object_v1(selection, fields, list(
      version = "dsvert-frequency-backend-selection-v3",
      policy_sha256 = .dsvert_dp_analysis_frequency_policy_sha256_v3(),
      objective = "minimum_certified_simultaneous_95_abs",
      tie_break = "convolution_laplace_v4_on_equal_certified_radius"))) fail()
  requests <- .dsvert_dp_analysis_frequency_candidate_requests_v3(privacy, calibration, dimension)
  primitives <- .dsvert_dp_analysis_frequency_primitives_v3()
  if (!.dsvert_dp_analysis_frequency_object_v1(selection$candidates, names(primitives))) fail()
  for (kind in names(primitives)) {
    candidate <- selection$candidates[[kind]]
    profile <- .dsvert_dp_analysis_frequency_profile_v1(primitives[[kind]])
    if (kind == "gaussian" && identical(candidate$available, FALSE)) {
      if (!.dsvert_dp_analysis_frequency_object_v1(candidate,
          c("available", "planner_request_sha256", "unavailable_reason")) ||
          !is.character(candidate$unavailable_reason) ||
          length(candidate$unavailable_reason) != 1L || is.na(candidate$unavailable_reason) ||
          !nzchar(candidate$unavailable_reason)) fail()
    } else {
      if (!.dsvert_dp_analysis_frequency_object_v1(candidate, c("available",
          "planner_request_sha256", "full_plan_sha256", "accuracy_certificate_sha256",
          "simultaneous_95_abs", "absolute_support"), list(available = TRUE)) ||
          !.dsvert_dp_analysis_frequency_hex_v1(candidate$full_plan_sha256) ||
          !.dsvert_dp_analysis_frequency_hex_v1(candidate$accuracy_certificate_sha256)) fail()
      radius <- .dsvert_dp_analysis_frequency_uint_v1(candidate$simultaneous_95_abs)
      support <- .dsvert_dp_analysis_frequency_uint_v1(candidate$absolute_support, TRUE)
      if (radius > support) fail()
      if (kind == "convolution" &&
          (!identical(candidate$absolute_support, format(bound, scientific = FALSE, trim = TRUE)) ||
           !identical(candidate$simultaneous_95_abs,
             .dsvert_dp_analysis_frequency_exact_radius_v4(requests$convolution, bound)))) fail()
      if (kind == "gaussian" && (calibration$implementation_delta <= 0 ||
          openssl::bignum(format(bound, scientific = FALSE, trim = TRUE)) + support >
            (openssl::bignum(2)^127) - 1)) fail()
    }
    if (!identical(candidate$planner_request_sha256,
        .dsvert_dp_analysis_frequency_hash_v1(profile$request_domain, requests[[kind]]))) fail()
  }
  candidates <- selection$candidates
  winner <- if (isTRUE(candidates$gaussian$available) &&
      openssl::bignum(candidates$gaussian$simultaneous_95_abs) <
      openssl::bignum(candidates$convolution$simultaneous_95_abs)) "gaussian" else "convolution"
  selected <- candidates[[winner]]
  expected <- .dsvert_dp_analysis_frequency_selection_certificate_v3(
    primitives[[winner]], selected$full_plan_sha256, selected$simultaneous_95_abs)
  if (!identical(selection$selected_primitive, primitives[[winner]]) ||
      !identical(selected_primitive, primitives[[winner]]) ||
      !identical(selection$selected_simultaneous_95_abs, selected$simultaneous_95_abs) ||
      !identical(plan$full_plan_sha256, selected$full_plan_sha256) ||
      !identical(plan$planner_request_sha256, selected$planner_request_sha256) ||
      !identical(selection$selection_certificate_sha256,
        .dsvert_dp_analysis_frequency_hash_v1(.DSVERT_DP_ANALYSIS_FREQUENCY_SELECTION_DOMAIN_V1, expected))) fail()
  invisible(TRUE)
}

.dsvert_dp_analysis_frequency_plan_summary_v2 <- function(config) {
  fail <- function() stop("Invalid Frequency backend selection", call. = FALSE)
  selection <- config$backend_selection
  if (!.dsvert_dp_analysis_frequency_object_v1(selection, c("summary", "selected_request",
      "selected_plan", "selected_accuracy_certificate", "selection_certificate"))) fail()
  primitive <- selection$summary$selected_primitive
  profile <- .dsvert_dp_analysis_frequency_profile_v1(primitive)
  if (is.null(profile)) fail()
  kind <- if (profile$gaussian) "gaussian" else "convolution"
  request <- .dsvert_dp_analysis_frequency_candidate_requests_v3(
    config$privacy, config$calibration, config$factor_domain$dimension)[[kind]]
  if (!identical(.dsvert_dp_analysis_frequency_hash_v1("", request),
                 .dsvert_dp_analysis_frequency_hash_v1("", selection$selected_request))) fail()
  full <- selection$selected_plan
  certificate <- selection$selected_accuracy_certificate
  .dsvert_dp_analysis_frequency_accuracy_validate_v3(
    certificate, full, primitive, request, config$coordinate_upper_bound)
  candidate <- selection$summary$candidates[[kind]]
  if (!identical(candidate$accuracy_certificate_sha256,
      .dsvert_dp_analysis_frequency_hash_v1(.DSVERT_DP_ANALYSIS_FREQUENCY_ACCURACY_DOMAIN_V1, certificate))) fail()
  expected_selection <- .dsvert_dp_analysis_frequency_selection_certificate_v3(
    primitive, certificate$plan_sha256, certificate$simultaneous_95_abs)
  if (!identical(.dsvert_dp_analysis_frequency_hash_v1("", selection$selection_certificate),
                 .dsvert_dp_analysis_frequency_hash_v1("", expected_selection))) fail()
  fraction <- function(n, d) list(numerator = as.character(n), denominator = as.character(d))
  allocated <- .dsvert_dp_analysis_frequency_decimal_fraction_v1(request$delta)
  implementation <- if (profile$gaussian) fraction(full$per_peer_implementation_delta_numerator,
    full$per_peer_implementation_delta_denominator) else fraction("0", "1")
  maximum <- if (profile$gaussian) full$maximum_noise_magnitude_per_peer else "unbounded"
  exact <- isTRUE(profile$exact)
  result <- list(version = if (exact) "dsvert-frequency-plan-summary-v2" else "dsvert-frequency-plan-summary-v1",
    physical_plan_version = profile$plan, full_plan_sha256 = certificate$plan_sha256,
    planner_request_sha256 = .dsvert_dp_analysis_frequency_hash_v1(profile$request_domain, request),
    coordinate_order_sha256 = .dsvert_dp_analysis_frequency_coordinate_order_sha256_v1(config$factor_domain$levels),
    d = config$factor_domain$dimension, chunk_coordinates = min(profile$max_chunk_coordinates, config$factor_domain$dimension),
    allocated_delta = fraction(allocated$numerator, allocated$denominator),
    core_delta = if (profile$gaussian) fraction(full$core_delta_numerator, full$core_delta_denominator) else fraction("0", "1"),
    implementation_delta = implementation, maximum_noise_per_peer = maximum,
    profile_sha256 = .dsvert_dp_analysis_frequency_hash_v1("dsVert/frequency/physical-profile/v1|", profile),
    backend_selection = selection$summary)
  upper <- format(config$coordinate_upper_bound, scientific = FALSE, trim = TRUE)
  if (exact) result$wrap_certificate <- .dsvert_joint_dp_pure_output_certificate(full, upper, 0L) else
    result$no_wrap_sha256 <- .dsvert_dp_analysis_frequency_hash_v1("dsVert/frequency/ring128-no-wrap/v1|",
      list(version = "dsvert-frequency-ring128-no-wrap-v1", coordinate_upper_bound = upper,
        maximum_noise_per_peer = maximum, maximum_noise_release = as.character(2 * openssl::bignum(maximum))))
  sensitivity <- list(value = if (config$privacy$adjacency == "replace_one_fixed_cohort")
    if (profile$gaussian) sqrt(2) else 2 else 1)
  .dsvert_dp_analysis_frequency_plan_validate_v1(result, profile, config$privacy,
    sensitivity, config$factor_domain$dimension, config$coordinate_upper_bound, config$calibration)
  if (!is.null(config$transport_chunk_coordinates) &&
      as.numeric(config$transport_chunk_coordinates) != result$chunk_coordinates) fail()
  result
}
