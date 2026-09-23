# Public exact arithmetic checks for the unbounded v4 Laplace contract.
.DSVERT_CLIENT_VECTOR_PURE_CERTIFICATE_FIELDS <- c(
  "guarantee", "randomness", "noise_support", "noise_commitment",
  "representability_bound", "wrap_bound_certified", "wrap_bound",
  "wrap_bound_event", "admitted_ranges")

.dsvert_vector_exact_exponential_bound <- function(numerator, denominator,
                                                   multiplicity) {
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

.dsvert_vector_exact_plan_valid <- function(plan, sensitivity_steps,
                                            coordinate_count) {
  tryCatch({
    if (!.dsvert_vector_integer_text(sensitivity_steps, TRUE, 512L) ||
        !identical(plan$sensitivity_steps, sensitivity_steps) ||
        !.dsvert_vector_whole(coordinate_count, 1, 1000000) ||
        !.dsvert_vector_integer_text(
          plan$epsilon_effective_upper_numerator, TRUE, 512L) ||
        !.dsvert_vector_integer_text(
          plan$epsilon_effective_upper_denominator, TRUE, 512L)) return(FALSE)
    numerator <- openssl::bignum(plan$epsilon_effective_upper_numerator)
    epsilon_denominator <- openssl::bignum(
      plan$epsilon_effective_upper_denominator)
    denominator <- epsilon_denominator * openssl::bignum(sensitivity_steps)
    if (numerator > 10000 * epsilon_denominator ||
        numerator * (openssl::bignum(2)^100) < denominator) return(FALSE)
    bound <- .dsvert_vector_exact_exponential_bound(
      numerator * (openssl::bignum(2)^127), denominator,
      2 * coordinate_count)
    expected <- list(guarantee = "pure-dp-under-ideal-bits",
      randomness = "keyed-stream-computational",
      noise_support = "unbounded_integer", noise_commitment = "modulo_2^128",
      representability_bound = bound, wrap_bound_certified = TRUE,
      wrap_bound = bound,
      wrap_bound_event = "at_least_one_peer_draw_outside_signed_Ring128",
      admitted_ranges = list(epsilon_minimum_exclusive = "0",
        epsilon_maximum = "10000", sensitivity_steps = "positive_integer",
        epsilon_over_sensitivity_minimum = "2^-100",
        total_coordinate_count_minimum = 1,
        total_coordinate_count_maximum = 1000000, ring_bits = 128))
    valid <- all(names(expected) %in% names(plan)) &&
      identical(.dsvert_joint_dp_client_json(plan[names(expected)]),
                .dsvert_joint_dp_client_json(expected)) &&
      identical(plan$maximum_noise_magnitude, "unbounded") &&
      identical(plan$complete_epsilon_per_peer, TRUE) &&
      identical(plan$epsilon_divided_by_peer_count, FALSE) &&
      identical(plan$release_implementation_delta_aggregation,
                "max_per_peer_not_sum") &&
      !any(c("no_wrap_headroom_certified", "no_wrap_certified") %in% names(plan)) &&
      identical(as.numeric(plan$stop_bits), 0) &&
      identical(plan$stop_numerator, "0") &&
      identical(as.numeric(plan$uniform_bits), 0) &&
      identical(as.numeric(plan$binary_geometric_bits), 0) &&
      length(plan$bernoulli_thresholds) == 0L &&
      identical(as.numeric(plan$independent_noise_peer_count), 2) &&
      identical(as.numeric(plan$geometric_variables_per_peer_per_coordinate), 2) &&
      identical(as.numeric(plan$geometric_variables_total_per_coordinate), 4)
    for (prefix in c("one_geometric_tv", "tail_upper", "rounding_upper",
        "implementation_delta", "per_peer_implementation_delta",
        "two_peer_ideal_transfer_delta")) {
      valid <- valid && identical(plan[[paste0(prefix, "_numerator")]], "0") &&
        identical(plan[[paste0(prefix, "_denominator")]], "1")
    }
    for (field in c("implementation_delta_bound",
        "per_peer_implementation_delta_bound", "two_peer_ideal_transfer_delta_bound")) {
      valid <- valid && identical(plan[[field]], "0")
    }
    isTRUE(valid)
  }, error = function(error) FALSE)
}

.dsvert_vector_exact_sum_certificate <- function(plan, maximum, dimension) {
  threshold <- ((openssl::bignum(2)^127) - maximum) %/% 2
  if (threshold < 0) stop("Invalid exact Laplace Ring128 source bound", call. = FALSE)
  numerator <- openssl::bignum(plan$epsilon_effective_upper_numerator)
  denominator <- openssl::bignum(plan$epsilon_effective_upper_denominator) *
    openssl::bignum(plan$sensitivity_steps)
  list(sum_wrap_threshold = as.character(threshold),
    sum_wrap_bound = .dsvert_vector_exact_exponential_bound(
      numerator * threshold, denominator, 2 * dimension))
}

.dsvert_vector_exact_probability_upper <- function(value) {
  if (!.dsvert_vector_string(value, "^(1|1e-[1-9][0-9]*)$", 512L)) {
    stop("Invalid exact Laplace utility bound", call. = FALSE)
  }
  number <- suppressWarnings(as.numeric(value))
  if (number == 0) .Machine$double.xmin else
    min(1, .dsvert_dp_vector_next_up(number))
}
