.DSVERT_DP_MINIMUM_EPSILON <- 2^-50
.DSVERT_DP_MAXIMUM_EPSILON <- 2^40
.DSVERT_DP_DEFAULT_CAPSULE_DELTA <- 2^-100
.DSVERT_DP_TRANSCRIPT_POLICY_VALUE <- paste0(
  "successful_release_semantic_messages_public_shape_or_dp_",
  "postprocessing_v1_physical_timing_availability_excluded")

#' Differential-privacy policy status
#'
#' Queries the reusable capsule control plane on every connected server. It
#' verifies the fixed per-capsule epsilon/delta, the authenticated lifetime
#' privacy bound, adjacency, exact
#' logical-name-to-Ed25519 pin map, privacy epoch, and the two designated noise
#' peers. Exact replay is unlimited and there is no request counter or accuracy
#' decay, while a new distinct capsule consumes one authenticated,
#' non-refundable distinct-capsule reservation unit at allocator commit. This
#' call releases no protected statistic and consumes no privacy allocation.
#' Printing the result distinguishes the remaining reservation units from
#' request quotas.
#'
#' @param datasources DataSHIELD connections. If `NULL`, use the active set.
#' @return A named list of server statuses with class `ds.vertDPStatus`.
#' @export
ds.vertDPStatus <- function(datasources = NULL) {
  .dsvert_joint_dp_capsule_status_impl(
    datasources, DSI::datashield.aggregate)
}

#' @export
print.ds.vertDPStatus <- function(x, ...) {
  reference <- x[[1L]]
  designated <- reference$policy$designated_noise_peers
  capsule <- x[[designated[[1L]]]]$composition_telemetry
  release <- x[[designated[[1L]]]]$release_instance_telemetry
  number <- function(value) {
    format(value, digits = 8L, scientific = TRUE, trim = TRUE)
  }

  cat("dsVert reusable joint-DP capsule status\n")
  cat("peers:", length(x), "| designated noise peers:",
      paste(designated, collapse = ", "), "\n")
  cat("per capsule: epsilon=", number(reference$policy$capsule_epsilon),
      ", delta=", number(reference$policy$capsule_delta),
      " | adjacency: ", reference$policy$adjacency, "\n", sep = "")
  cat("allocator-committed reservation units:",
      number(capsule$capsules_created),
      "| basic composition upper bound: epsilon=",
      number(capsule$cumulative_epsilon_upper_bound), ", delta=",
      number(capsule$cumulative_delta_upper_bound), "\n")
  cat("published release instances:", number(release$releases_published),
      "| basic composition upper bound: epsilon=",
      number(release$cumulative_epsilon_upper_bound), ", delta=",
      number(release$cumulative_delta_upper_bound), "\n")
  cat("allocator-committed reservation units remaining:",
      number(capsule$remaining_distinct_capsules), "of",
      number(capsule$lifetime_max_distinct_capsules), "\n")
  cat(paste(
    "new-capsule admission gate: enforced; request quota: none; exact replay",
    "of the same release is unlimited post-processing.\n"))
  invisible(x)
}

#' Preview DP utility candidates without accessing protected data
#'
#' Computes a public, deterministic utility preview for granular Laplace and
#' an analytic Gaussian approximation. Both the joint discrete-Laplace and
#' fixed-work dyadic discrete-Gaussian mechanisms are deployed. The exact
#' Gaussian planner uses conservative rational arithmetic in the server
#' worker, so this client-only preview never presents its approximation as the
#' executed choice. The signed capsule manifest is the sole authority for the
#' selected mechanism and its exact accuracy certificate.
#'
#' @param capsule_epsilon Fixed consortium epsilon assigned once to each
#'   immutable capsule being considered. Reusing a capsule does not change it.
#' @param peer_count Number of pinned vertical peers in the deployment being
#'   planned; at least two are required. This is metadata only: capsule epsilon
#'   is never divided by the peer count.
#' @param sensitivity Per-coordinate L1 sensitivity. For bounded moments,
#'   account for the three-way epsilon split by using effective sensitivities
#'   `3`, `3 * numeric_grid_scale`, and `3 * numeric_grid_scale` separately.
#' @param confidence Central probability used for the separate planning error
#'   radius. Candidate selection always reproduces the deployed 95% selector.
#' @param capsule_delta Fixed consortium delta assigned to each immutable
#'   capsule. Zero disables Gaussian.
#' @param coordinate_count Public vector dimension used by the simultaneous
#'   selector and its conservative
#'   `(1 + exp(epsilon)) * d * 2^-40` DP-transfer reserve.
#' @param gaussian_l2_sensitivity Public joint L2 sensitivity. Defaults to
#'   `sensitivity` for scalar queries.
#' @param objective Selector objective, or `"auto"` for marginal on scalars
#'   and simultaneous 95% radius on vectors.
#' @return A data frame with effective epsilon, planning approximations,
#'   upstream-compatible confidence radii, sampler-domain preview flags, and
#'   an explicit marker that the signed server manifest is authoritative.
#' @export
ds.vertDPCalibrate <- function(capsule_epsilon = c(1, 3),
                               peer_count = 2L,
                               sensitivity = 1,
                               confidence = 0.95,
                               capsule_delta = .DSVERT_DP_DEFAULT_CAPSULE_DELTA,
                               coordinate_count = 1L,
                               gaussian_l2_sensitivity = sensitivity,
                               objective = "auto") {
  if (!is.numeric(capsule_epsilon) || !length(capsule_epsilon) ||
      anyNA(capsule_epsilon) || any(!is.finite(capsule_epsilon)) ||
      any(capsule_epsilon <= 0)) {
    stop("capsule_epsilon must contain positive finite values",
         call. = FALSE)
  }
  if (!is.numeric(peer_count) || length(peer_count) != 1L ||
      is.na(peer_count) || !is.finite(peer_count) || peer_count < 2 ||
      peer_count > .Machine$integer.max || peer_count != floor(peer_count)) {
    stop("peer_count must be one positive integer of at least two",
         call. = FALSE)
  }
  peer_count <- as.integer(peer_count)
  if (!is.numeric(sensitivity) || length(sensitivity) != 1L ||
      is.na(sensitivity) || !is.finite(sensitivity) || sensitivity <= 0 ||
      sensitivity != floor(sensitivity) || sensitivity > 2^53 - 1) {
    stop("sensitivity must be one positive exactly representable integer ",
         "not exceeding 2^53 - 1", call. = FALSE)
  }
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      is.na(confidence) || !is.finite(confidence) ||
      confidence <= 0 || confidence >= 1) {
    stop("confidence must lie strictly between zero and one", call. = FALSE)
  }
  if (!is.numeric(capsule_delta) || length(capsule_delta) != 1L ||
      is.na(capsule_delta) || !is.finite(capsule_delta) ||
      capsule_delta < 0 || capsule_delta >= 1) {
    stop("capsule_delta must be one finite number in [0, 1)",
         call. = FALSE)
  }
  if (!is.numeric(coordinate_count) || length(coordinate_count) != 1L ||
      is.na(coordinate_count) || !is.finite(coordinate_count) ||
      coordinate_count < 1 || coordinate_count > .DSVERT_DP_MAX_COORDINATES ||
      coordinate_count != floor(coordinate_count)) {
    stop("coordinate_count must be one supported positive integer",
         call. = FALSE)
  }
  coordinate_count <- as.integer(coordinate_count)
  if (!is.numeric(gaussian_l2_sensitivity) ||
      length(gaussian_l2_sensitivity) != 1L ||
      is.na(gaussian_l2_sensitivity) ||
      !is.finite(gaussian_l2_sensitivity) ||
      gaussian_l2_sensitivity <= 0 ||
      gaussian_l2_sensitivity > 2^53 - 1) {
    stop("gaussian_l2_sensitivity must be positive, finite, and not exceed ",
         "2^53 - 1", call. = FALSE)
  }
  if (!is.character(objective) || length(objective) != 1L ||
      is.na(objective) ||
      !objective %in% c("auto", "marginal_95_abs",
                        "simultaneous_95_abs")) {
    stop("objective must be auto, marginal_95_abs, or simultaneous_95_abs",
         call. = FALSE)
  }
  if (identical(objective, "auto")) {
    objective <- if (coordinate_count == 1L) {
      "marginal_95_abs"
    } else {
      "simultaneous_95_abs"
    }
  }
  if (identical(objective, "marginal_95_abs") && coordinate_count != 1L) {
    stop("marginal_95_abs is valid only for coordinate_count = 1",
         call. = FALSE)
  }

  grid <- data.frame(
    capsule_epsilon = as.numeric(capsule_epsilon),
    stringsAsFactors = FALSE)
  grid$capsule_delta <- capsule_delta
  metrics <- .dsvert_dp_google_calibration_metrics(
    grid$capsule_epsilon / sensitivity, confidence,
    release_epsilon = grid$capsule_epsilon)
  grid$expected_absolute_error <- metrics$expected_absolute_error
  grid$rmse <- metrics$rmse
  grid$error_radius <- metrics$error_radius
  grid$sampler_supported <- metrics$sampler_supported
  # The deployed selector contract is specifically a conservative 95% rule.
  # `confidence` controls only the separate planning/reporting radius above;
  # it must never change which candidate this 95%-named selector predicts.
  selector_alpha <- 0.05 /
    if (identical(objective, "simultaneous_95_abs")) coordinate_count else 1
  grid$laplace_selection_radius <- vapply(
    seq_len(nrow(grid)), function(index) {
      .dsvert_dp_google_ci_radius(
        grid$capsule_epsilon[[index]], sensitivity,
        alpha = selector_alpha)
    }, numeric(1L))
  implementation_delta <- vapply(
    grid$capsule_epsilon,
    function(epsilon) {
      .dsvert_dp_gaussian_implementation_delta_bound(
        coordinate_count, epsilon)
    },
    numeric(1L))
  grid$gaussian_implementation_delta_bound <- implementation_delta
  grid$gaussian_analytic_delta <- pmax(
    0, grid$capsule_delta - implementation_delta)
  gaussian <- lapply(seq_len(nrow(grid)), function(index) {
    epsilon <- grid$capsule_epsilon[[index]]
    analytic_delta <- grid$gaussian_analytic_delta[[index]]
    if (!is.finite(epsilon) || epsilon < .DSVERT_DP_MINIMUM_EPSILON ||
        epsilon > .DSVERT_DP_MAXIMUM_EPSILON || analytic_delta <= 0) {
      return(c(sigma = NA_real_, radius = Inf, supported = 0))
    }
    sigma <- .dsvert_dp_gaussian_sigma_for_calibration(
      gaussian_l2_sensitivity, epsilon, analytic_delta)
    granularity <- if (is.finite(sigma) && sigma > 0) {
      2^ceiling(log2(sigma / 2^56))
    } else {
      Inf
    }
    radius <- .dsvert_dp_gaussian_accuracy_radius(sigma, selector_alpha)
    supported <- is.finite(sigma) && sigma > 0 &&
      is.finite(granularity) && granularity > 0 && granularity <= 1 &&
      is.finite(radius)
    c(sigma = sigma, radius = radius, supported = as.numeric(supported))
  })
  gaussian <- do.call(rbind, gaussian)
  grid$gaussian_sigma <- gaussian[, "sigma"]
  grid$gaussian_selection_radius <- gaussian[, "radius"]
  grid$gaussian_sampler_supported <- gaussian[, "supported"] == 1
  gaussian_wins <- grid$gaussian_sampler_supported &
    (!grid$sampler_supported |
       grid$gaussian_selection_radius < grid$laplace_selection_radius)
  grid$utility_candidate <- ifelse(
    gaussian_wins, "gaussian",
    ifelse(grid$sampler_supported, "laplace", "none"))
  grid$preview_candidate <- grid$utility_candidate
  # The server's fixed-work Gaussian planner uses exact rational arithmetic
  # and can reject public parameter combinations that this floating-point
  # preview cannot certify.  NA is intentional: reporting the preview as the
  # executed selection would create a false certificate.
  grid$selected_candidate <- NA_character_
  grid$selected_delta <- NA_real_
  grid$formal_backend <- "signed_server_capsule_manifest"
  grid$deployed_backends <- paste(
    "joint-discrete-laplace-v3",
    "dyadic-discrete-gaussian-tv-bounded-v2", sep = ";")
  grid$gaussian_backend_deployed <- TRUE
  grid$gaussian_preview_supported <- grid$gaussian_sampler_supported
  grid$selection_authority <- "signed_server_capsule_manifest"
  grid$deployment_decision <-
    "preview_only_server_plan_authoritative"
  grid$selector <- .DSVERT_DP_NOISE_SELECTOR
  grid$selector_objective <- objective
  grid$calibration_model <-
    "google_dp_v4.1.0_ci_with_continuous_laplace_moment_proxies"
  grid$metric_guarantee <- paste(
    "error_radius matches ComputeConfidenceIntervalInt64 when supported;",
    "expected_absolute_error and rmse are planning approximations")
  grid$confidence <- confidence
  grid$sensitivity <- sensitivity
  grid$gaussian_l2_sensitivity <- gaussian_l2_sensitivity
  grid$coordinate_count <- coordinate_count
  grid$peer_count <- peer_count
  grid$privacy_parameters_are_fixed <- TRUE
  grid$operation_accounting <- "none"
  grid$operation_limit <- FALSE
  grid$request_limit <- FALSE
  grid$history_can_deny_operation <- FALSE
  grid$accuracy_depends_on_request_history <- FALSE
  grid <- grid[order(grid$capsule_epsilon), , drop = FALSE]
  rownames(grid) <- NULL
  class(grid) <- c("ds.vertDPCalibration", class(grid))
  grid
}

.dsvert_dp_gaussian_sigma_for_calibration <- function(
    sensitivity, epsilon, delta) {
  lower <- 0
  upper <- sensitivity
  for (iteration in seq_len(2048L)) {
    achieved <- .dsvert_dp_gaussian_achieved_delta(
      upper, sensitivity, epsilon)
    if (!is.finite(achieved)) return(NA_real_)
    if (achieved <= delta) break
    lower <- upper
    upper <- 2 * upper
    if (!is.finite(upper)) return(NA_real_)
  }
  if (.dsvert_dp_gaussian_achieved_delta(
        upper, sensitivity, epsilon) > delta) return(NA_real_)
  for (iteration in seq_len(2048L)) {
    if (upper - lower <= 1e-3 * lower) break
    middle <- 0.5 * lower + 0.5 * upper
    achieved <- .dsvert_dp_gaussian_achieved_delta(
      middle, sensitivity, epsilon)
    if (!is.finite(achieved)) return(NA_real_)
    if (achieved > delta) lower <- middle else upper <- middle
  }
  upper
}

.dsvert_dp_gaussian_sigma_matches_calibration <- function(
    sigma, sensitivity, epsilon, delta) {
  expected <- .dsvert_dp_gaussian_sigma_for_calibration(
    sensitivity, epsilon, delta)
  is.finite(expected) && expected > 0 && is.finite(sigma) && sigma > 0 &&
    abs(sigma - expected) <= 0.002 * max(sigma, expected)
}

.dsvert_dp_google_calibration_metrics <- function(epsilon_over_sensitivity,
                                                  confidence = 0.95,
                                                  release_epsilon = NULL) {
  t <- as.numeric(epsilon_over_sensitivity)
  if (is.null(release_epsilon)) release_epsilon <- rep(Inf, length(t))
  release_epsilon <- rep_len(as.numeric(release_epsilon), length(t))
  expected_absolute_error <- rep(Inf, length(t))
  rmse <- rep(Inf, length(t))
  error_radius <- rep(Inf, length(t))
  sampler_supported <- rep(FALSE, length(t))

  positive <- !is.na(t) & t > 0
  if (any(positive)) {
    tp <- t[positive]
    lambda <- 1 / tp
    expected_absolute_error[positive] <- lambda
    rmse[positive] <- sqrt(2) * lambda
    supported <- is.finite(lambda) & lambda > 0 & lambda <= 2^40 &
      is.finite(release_epsilon[positive]) &
      release_epsilon[positive] >= .DSVERT_DP_MINIMUM_EPSILON &
      release_epsilon[positive] <= .DSVERT_DP_MAXIMUM_EPSILON
    sampler_supported[positive] <- supported
    radius <- rep(NA_real_, length(tp))
    radius[supported] <- floor(
      -lambda[supported] * log1p(-confidence) + 0.5)
    error_radius[positive] <- radius
  }
  list(expected_absolute_error = expected_absolute_error,
       rmse = rmse, error_radius = error_radius,
       sampler_supported = sampler_supported)
}

.dsvert_dp_datasources <- function(datasources) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (!is.list(datasources) || !length(datasources) ||
      is.null(names(datasources)) || anyNA(names(datasources)) ||
      any(!nzchar(names(datasources))) || anyDuplicated(names(datasources))) {
    stop("datasources must be a non-empty, uniquely named connection list",
         call. = FALSE)
  }
  datasources
}

.dsvert_dp_has_exact_names <- function(value, expected) {
  is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
    !anyDuplicated(names(value)) && length(value) == length(expected) &&
    setequal(names(value), expected)
}

.dsvert_dp_is_string <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    nzchar(value)
}

.dsvert_dp_is_number <- function(value, lower = -Inf, upper = Inf,
                                 lower_open = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value)) return(FALSE)
  lower_ok <- if (lower_open) value > lower else value >= lower
  isTRUE(lower_ok && value <= upper)
}

.dsvert_dp_is_integer <- function(value, lower = 0,
                                  upper = 2^53 - 1) {
  .dsvert_dp_is_number(value, lower, upper) && value == floor(value)
}

.DSVERT_DP_MAX_COORDINATES <- 1000000L

.dsvert_dp_num_equal <- function(left, right, multiplier = 128) {
  if (!.dsvert_dp_is_number(left) || !.dsvert_dp_is_number(right)) {
    return(FALSE)
  }
  scale <- max(abs(left), abs(right), .Machine$double.xmin)
  abs(left - right) <= multiplier * .Machine$double.eps * scale
}

.dsvert_dp_validate_named_domains <- function(value) {
  if (!is.list(value)) return(FALSE)
  if (!length(value)) return(is.null(names(value)))
  if (is.null(names(value)) || anyNA(names(value)) ||
      any(!nzchar(names(value))) || anyDuplicated(names(value)) ||
      !identical(names(value), sort(names(value), method = "radix"))) {
    return(FALSE)
  }
  all(vapply(value, function(domain) {
    is.character(domain) && length(domain) > 0L &&
      length(domain) <= .DSVERT_DP_MAX_COORDINATES && !anyNA(domain) &&
      all(nzchar(domain)) && !anyDuplicated(domain) &&
      identical(domain, sort(domain, method = "radix"))
  }, logical(1L)))
}

.dsvert_dp_validate_numeric_bounds <- function(value) {
  if (!is.list(value)) return(FALSE)
  if (!length(value)) return(is.null(names(value)))
  if (is.null(names(value)) || anyNA(names(value)) ||
      any(!nzchar(names(value))) || anyDuplicated(names(value)) ||
      !identical(names(value), sort(names(value), method = "radix"))) {
    return(FALSE)
  }
  all(vapply(value, function(bounds) {
    if (!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds) ||
        any(!is.finite(bounds)) || bounds[[1L]] >= bounds[[2L]]) {
      return(FALSE)
    }
    width <- bounds[[2L]] - bounds[[1L]]
    is.finite(width) && width > 0 && is.finite(width * width) &&
      width * width > 0
  }, logical(1L)))
}

.dsvert_dp_validate_datasets <- function(value) {
  if (!is.list(value) || !length(value) || is.null(names(value)) ||
      anyNA(names(value)) || any(!nzchar(names(value))) ||
      anyDuplicated(names(value)) ||
      !identical(names(value), sort(names(value), method = "radix"))) {
    return(FALSE)
  }
  identities <- character(length(value))
  valid <- vapply(seq_along(value), function(index) {
    descriptor <- value[[index]]
    if (!.dsvert_dp_has_exact_names(descriptor, c(
          "id", "version", "alignment_attested",
          "alignment_protocol_version")) ||
        !.dsvert_dp_is_string(descriptor$id) ||
        !.dsvert_dp_is_string(descriptor$version) ||
        !identical(descriptor$alignment_attested, TRUE) ||
        !.dsvert_dp_is_integer(
          descriptor$alignment_protocol_version, 1,
          .Machine$integer.max)) return(FALSE)
    pattern <- "^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$"
    if (!grepl(pattern, descriptor$id) ||
        !grepl(pattern, descriptor$version)) return(FALSE)
    identities[[index]] <<- paste(descriptor$id, descriptor$version,
                                  sep = "@")
    TRUE
  }, logical(1L))
  all(valid) && !anyDuplicated(identities)
}

.dsvert_dp_normalize_identity_pk <- function(value) {
  if (!.dsvert_dp_is_string(value) ||
      !grepl("^[A-Za-z0-9_-]{43}$", value)) return(NULL)
  standard <- chartr("-_", "+/", paste0(value, "="))
  decoded <- tryCatch(jsonlite::base64_dec(standard),
                      error = function(e) raw(0L))
  if (!is.raw(decoded) || length(decoded) != 32L) return(NULL)
  canonical <- chartr(
    "+/", "-_",
    sub("=+$", "", gsub("[\r\n]", "", jsonlite::base64_enc(decoded)),
        perl = TRUE))
  if (!identical(canonical, value)) NULL else canonical
}

.dsvert_dp_pinset_hash <- function(pinset) {
  canonical <- as.character(jsonlite::toJSON(
    as.list(pinset), auto_unbox = TRUE, null = "null", na = "null",
    digits = 17, pretty = FALSE))
  digest <- openssl::sha256(charToRaw(canonical))
  paste(format(digest), collapse = "")
}

.dsvert_dp_expected_anchor_id <- function(policy) {
  canonical <- as.character(jsonlite::toJSON(list(
    domain = policy$domain,
    cohort_id = policy$cohort_id,
    peer_name = policy$peer_name,
    pinset_sha256 = policy$peer_pinset_sha256
  ), auto_unbox = TRUE, null = "null", na = "null", digits = 17,
  pretty = FALSE))
  digest <- openssl::sha256(charToRaw(canonical))
  paste0("dpa1_", paste(format(digest), collapse = ""))
}

.dsvert_dp_validate_pinset <- function(value) {
  if (!is.character(value) || !length(value) || is.null(names(value)) ||
      anyNA(names(value)) || any(!nzchar(names(value))) ||
      any(!grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", names(value))) ||
      anyDuplicated(names(value)) || anyDuplicated(unname(value)) ||
      !identical(names(value), sort(names(value), method = "radix"))) {
    return(FALSE)
  }
  all(vapply(value, function(pk) {
    !is.null(.dsvert_dp_normalize_identity_pk(pk))
  }, logical(1L)))
}

.dsvert_dp_validate_rollback_policy <- function(value) {
  expected <- c(
    "schema_version", "mode", "anchor_id", "provider_contract",
    "external", "monotonic")
  common <- .dsvert_dp_has_exact_names(value, expected) &&
    .dsvert_dp_is_integer(value$schema_version, 2, 2) &&
    .dsvert_dp_is_string(value$anchor_id) &&
    grepl("^dpa1_[0-9a-f]{64}$", value$anchor_id) &&
    identical(value$monotonic, TRUE)
  if (!isTRUE(common)) return(FALSE)
  if (identical(value$mode, "external_cas")) {
    identical(value$provider_contract,
              "custodian_injected_durable_linearizable_receipt_cas_v2") &&
      identical(value$external, TRUE)
  } else if (identical(
      value$mode, "pinned_peer_global_allocator_pending")) {
    identical(value$provider_contract,
              "pinned_ed25519_cross_signed_global_allocator_v1") &&
      identical(value$external, FALSE)
  } else if (identical(
      value$mode,
      "pinned_peer_global_allocator_pending_plus_external_cas")) {
    identical(
      value$provider_contract,
      paste0(
        "pinned_ed25519_cross_signed_global_allocator_v1_",
        "plus_external_cas_v2")) &&
      identical(value$external, TRUE)
  } else {
    FALSE
  }
}

.dsvert_dp_validate_noise_selection_policy <- function(value) {
  expected <- c(
    "schema_version", "selector", "candidates", "scalar_objective",
    "vector_objective", "confidence", "tie_break", "delta_accounting")
  .dsvert_dp_has_exact_names(value, expected) &&
    .dsvert_dp_is_integer(value$schema_version, 1, 1) &&
    identical(value$selector, "minimum_conservative_95_radius_v3") &&
    identical(as.character(value$candidates), c(
      "deterministic_granular_laplace_int64_l1",
      paste0("deterministic_approximate_gaussian_int64_l2_",
             "dp_transfer_tv_accounted"))) &&
    identical(value$scalar_objective, "marginal_95_abs") &&
    identical(value$vector_objective, "simultaneous_95_abs") &&
    .dsvert_dp_num_equal(value$confidence, 0.95) &&
    identical(
      value$tie_break,
      "laplace_unless_gaussian_strictly_improves") &&
    identical(
      value$delta_accounting,
      "laplace_zero_gaussian_scheduled_allocation_only")
}

.dsvert_dp_validate_policy <- function(policy, server) {
  invalid <- function() {
    stop("Server '", server, "' returned an invalid DP v7 policy",
         call. = FALSE)
  }
  if (!is.list(policy)) invalid()
  global_policy <-
    .dsvert_dp_is_integer(policy$schema_version, 7, 7) &&
    identical(
      policy$policy_contract,
      "single_safe_policy_v4_global_allocator")
  legacy_test_policy <-
    .dsvert_dp_is_integer(policy$schema_version, 6, 6) &&
    identical(policy$policy_contract, "single_safe_policy_v3")
  if (!global_policy && !legacy_test_policy) invalid()
  accounting_fields <- if (global_policy) {
    c("global_total_epsilon", "global_total_delta", "peer_count",
      "designated_noise_peers")
  } else {
    c("global_total_epsilon", "global_total_delta",
      "composition_partitions", "local_total_epsilon",
      "local_total_delta")
  }
  expected <- c(
    "schema_version", "mechanism_version", "policy_contract", "domain",
    "cohort_id",
    "peer_name",
    "own_identity_pk", "logical_peers", "peer_pinset",
    "peer_pinset_sha256", accounting_fields, "decay",
    "adjacency", "numeric_grid_bits",
    "patient_column", "unit_capacity", "fixed_cohort_size",
    "max_records_per_unit", "overflow_policy", "noise_selection",
    "contingency_unit_aggregation_policy",
    "transcript_privacy",
    "snapshot_binding", "alignment_binding",
    "rollback_protection", "categorical_levels", "numeric_bounds",
    "datasets")
  if (!.dsvert_dp_has_exact_names(policy, expected) ||
      !identical(policy$mechanism_version,
                 "dsvert-dp-v7-contingency-unit-aggregation-1") ||
      !.dsvert_dp_is_string(policy$domain) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", policy$domain) ||
      !.dsvert_dp_is_string(policy$cohort_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", policy$cohort_id) ||
      !.dsvert_dp_is_string(policy$peer_name) ||
      !policy$peer_name %in% names(policy$peer_pinset) ||
      !.dsvert_dp_is_string(policy$own_identity_pk) ||
      !.dsvert_dp_validate_pinset(policy$peer_pinset) ||
      !is.character(policy$logical_peers) ||
      !identical(policy$logical_peers, names(policy$peer_pinset)) ||
      !identical(policy$own_identity_pk,
                 unname(policy$peer_pinset[[policy$peer_name]])) ||
      !.dsvert_dp_is_string(policy$peer_pinset_sha256) ||
      !grepl("^[0-9a-f]{64}$", policy$peer_pinset_sha256) ||
      !identical(.dsvert_dp_pinset_hash(policy$peer_pinset),
                 policy$peer_pinset_sha256) ||
      !.dsvert_dp_is_number(policy$global_total_epsilon, 0, Inf, TRUE) ||
      !.dsvert_dp_is_number(policy$global_total_delta, 0, 1) ||
      policy$global_total_delta >= 1 ||
      !.dsvert_dp_is_number(policy$decay, 0, 1, TRUE) ||
      policy$decay >= 1 ||
      !.dsvert_dp_is_integer(policy$numeric_grid_bits, 8, 18) ||
      !.dsvert_dp_is_string(policy$adjacency) ||
      !policy$adjacency %in% c(
        "add_remove_patient", "replace_one_fixed_cohort") ||
      !.dsvert_dp_is_string(policy$patient_column) ||
      !.dsvert_dp_is_integer(
        policy$unit_capacity, 1, .DSVERT_DP_MAX_COORDINATES) ||
      !.dsvert_dp_is_integer(
        policy$max_records_per_unit, 1, .DSVERT_DP_MAX_COORDINATES) ||
      !identical(policy$overflow_policy, "reject_snapshot") ||
      !identical(
        policy$contingency_unit_aggregation_policy,
        "consistent_cell_else_exclude_v1") ||
      !.dsvert_dp_validate_noise_selection_policy(
        policy$noise_selection) ||
      !identical(
        policy$transcript_privacy,
        .DSVERT_DP_TRANSCRIPT_POLICY_VALUE) ||
      !identical(policy$snapshot_binding,
                 "required_private_custodian_sha256") ||
      !identical(policy$alignment_binding,
                 "required_private_ordered_psi_manifest") ||
      !.dsvert_dp_validate_named_domains(policy$categorical_levels) ||
      !.dsvert_dp_validate_numeric_bounds(policy$numeric_bounds) ||
      !.dsvert_dp_validate_datasets(policy$datasets) ||
      !.dsvert_dp_validate_rollback_policy(
        policy$rollback_protection)) invalid()

  if (!identical(policy$rollback_protection$anchor_id,
                 .dsvert_dp_expected_anchor_id(policy))) invalid()

  if (identical(policy$adjacency, "replace_one_fixed_cohort")) {
    if (!.dsvert_dp_is_integer(
          policy$fixed_cohort_size, 1, .DSVERT_DP_MAX_COORDINATES) ||
        policy$fixed_cohort_size != policy$unit_capacity) invalid()
  } else if (!is.null(policy$fixed_cohort_size)) {
    invalid()
  }
  if (global_policy) {
    if (!.dsvert_dp_is_integer(
          policy$peer_count, 2, .Machine$integer.max) ||
        policy$peer_count != length(policy$peer_pinset) ||
        !is.character(policy$designated_noise_peers) ||
        length(policy$designated_noise_peers) != 2L ||
        anyNA(policy$designated_noise_peers) ||
        anyDuplicated(policy$designated_noise_peers) ||
        !identical(
          policy$designated_noise_peers,
          sort(policy$designated_noise_peers, method = "radix")) ||
        !all(policy$designated_noise_peers %in% names(policy$peer_pinset)) ||
        !policy$rollback_protection$mode %in% c(
          "pinned_peer_global_allocator_pending",
          "pinned_peer_global_allocator_pending_plus_external_cas") ||
        policy$global_total_epsilon * (1 - policy$decay) >
          .DSVERT_DP_MAXIMUM_EPSILON) invalid()
  } else if (!.dsvert_dp_is_integer(
        policy$composition_partitions, 1, .Machine$integer.max) ||
      policy$composition_partitions != length(policy$peer_pinset) ||
      !.dsvert_dp_is_number(policy$local_total_epsilon, 0, Inf, TRUE) ||
      !.dsvert_dp_is_number(policy$local_total_delta, 0, 1) ||
      !.dsvert_dp_num_equal(
        policy$local_total_epsilon,
        policy$global_total_epsilon / policy$composition_partitions) ||
      !.dsvert_dp_num_equal(
        policy$local_total_delta,
        policy$global_total_delta / policy$composition_partitions) ||
      policy$local_total_epsilon * (1 - policy$decay) >
        .DSVERT_DP_MAXIMUM_EPSILON ||
      !identical(policy$rollback_protection$mode, "external_cas")) {
    invalid()
  }
  invisible(TRUE)
}

.dsvert_dp_validate_noise_root <- function(value) {
  expected <- c(
    "protocol", "provider_id", "key_id", "privacy_epoch", "external",
    "storage", "automatic_generation", "automatic_recovery",
    "automatic_rotation", "rotation_count", "key_material_exposed")
  .dsvert_dp_has_exact_names(value, expected) &&
    identical(value$protocol, "dsvert-dp-noise-root-v1") &&
    .dsvert_dp_is_string(value$provider_id) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value$provider_id) &&
    .dsvert_dp_is_string(value$key_id) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value$key_id) &&
    .dsvert_dp_is_integer(value$privacy_epoch, 1) &&
    is.logical(value$external) && length(value$external) == 1L &&
    !is.na(value$external) &&
    is.logical(value$automatic_generation) &&
    length(value$automatic_generation) == 1L &&
    !is.na(value$automatic_generation) &&
    is.logical(value$automatic_recovery) &&
    length(value$automatic_recovery) == 1L &&
    !is.na(value$automatic_recovery) &&
    is.logical(value$automatic_rotation) &&
    length(value$automatic_rotation) == 1L &&
    !is.na(value$automatic_rotation) &&
    .dsvert_dp_is_integer(value$rotation_count, 0, 2^53 - 1) &&
    .dsvert_dp_is_string(value$storage) &&
    value$storage %in% c("hsm_kms_provider", "owner_only_file") &&
    identical(value$external, identical(value$storage, "hsm_kms_provider")) &&
    identical(
      value$automatic_generation,
      identical(value$storage, "owner_only_file")) &&
    (!isTRUE(value$automatic_recovery) ||
       identical(value$storage, "owner_only_file")) &&
    (identical(value$storage, "owner_only_file") ||
       (!isTRUE(value$automatic_rotation) && value$rotation_count == 0)) &&
    identical(value$key_material_exposed, FALSE)
}

.dsvert_dp_server <- function(server, datasources) {
  if (is.null(server)) return(sort(names(datasources), method = "radix")[[1L]])
  if (!is.character(server) || length(server) != 1L || is.na(server) ||
      !nzchar(server) || !server %in% names(datasources)) {
    stop("server must name one connected datasource", call. = FALSE)
  }
  server
}

.dsvert_dp_find_server <- function(data_name, variables, server, datasources,
                                   .aggregate) {
  if (!is.null(server)) return(.dsvert_dp_server(server, datasources))
  columns <- .dsvert_aggregate_strict(
    conns = datasources,
    expr = call(name = "dsvertColNamesDS", data_name = data_name),
    operation = "DP variable-location handshake", .aggregate = .aggregate)
  candidates <- names(datasources)[vapply(names(datasources), function(name) {
    response <- columns[[name]]
    is.list(response) && all(variables %in% response$columns)
  }, logical(1L))]
  if (!length(candidates)) {
    stop("No connected server holds every requested variable", call. = FALSE)
  }
  if (length(candidates) > 1L) {
    stop("The requested variables occur on multiple servers; specify server",
         call. = FALSE)
  }
  candidates[[1L]]
}

# Schema-v6 parsing remains only to regression-test retired certificates and to
# report a clear status for old deployments.  The release allocator gate above
# rejects schema v6, so these epsilon/K helpers cannot authorize a DSI call.
.dsvert_dp_release_epsilon_is_valid <- function(epsilon, policy) {
  if (!.dsvert_dp_is_number(epsilon, 0, Inf, TRUE) ||
      epsilon > policy$local_total_epsilon) return(FALSE)
  first <- policy$local_total_epsilon * (1 - policy$decay)
  if (!is.finite(first) || first <= 0 || epsilon > first &&
      !.dsvert_dp_num_equal(epsilon, first)) return(FALSE)
  index <- log(epsilon / first) / log(policy$decay)
  if (!is.finite(index) || index < -128 * .Machine$double.eps) return(FALSE)
  index <- max(0, round(index))
  .dsvert_dp_num_equal(
    epsilon, first * policy$decay^index, multiplier = 1024)
}

.dsvert_dp_validate_release_common <- function(result, status, server) {
  policy <- status$policy
  invalid <- function() {
    stop("Server '", server, "' returned an invalid DP v7 release",
         call. = FALSE)
  }
  if (!is.logical(result$released) || length(result$released) != 1L ||
      is.na(result$released) ||
      !identical(result$adjacency, policy$adjacency) ||
      !.dsvert_dp_is_number(result$epsilon, 0) ||
      !.dsvert_dp_is_number(result$delta, 0, 1)) invalid()
  public_fixed_count <- isTRUE(result$released) && identical(
    policy$adjacency, "replace_one_fixed_cohort") &&
    identical(result$mechanism, "public_fixed_cohort_size_v1") &&
    .dsvert_dp_num_equal(result$epsilon, 0) &&
    .dsvert_dp_num_equal(result$delta, 0)
  if (isTRUE(public_fixed_count)) {
    if (!.dsvert_dp_is_integer(result$peer_count, 2,
                               .Machine$integer.max) ||
        !identical(as.numeric(result$peer_count),
                   as.numeric(length(policy$peer_pinset)))) invalid()
    return(invisible(TRUE))
  }
  if (!.dsvert_dp_is_integer(result$composition_partitions, 1,
                             .Machine$integer.max) ||
      !.dsvert_dp_num_equal(result$composition_partitions,
                            policy$composition_partitions)) invalid()
  if (isTRUE(result$released)) {
    if (isTRUE(public_fixed_count)) return(invisible(TRUE))
    scheduled_delta <- result$epsilon *
      policy$local_total_delta / policy$local_total_epsilon
    valid_delta <- .dsvert_dp_num_equal(result$delta, 0) ||
      .dsvert_dp_num_equal(result$delta, scheduled_delta,
                           multiplier = 2048)
    if (!.dsvert_dp_release_epsilon_is_valid(result$epsilon, policy) ||
        !valid_delta ||
        !.dsvert_dp_is_integer(result$privacy_epoch, 1) ||
        !.dsvert_dp_num_equal(
          result$privacy_epoch, status$noise_root$privacy_epoch) ||
        !identical(result$noise_key_id, status$noise_root$key_id) ||
        !identical(result$sticky_noise, "dsvert-sticky-noise-v1")) invalid()
  } else if (!.dsvert_dp_num_equal(result$epsilon, 0) ||
             !.dsvert_dp_num_equal(result$delta, 0)) {
    invalid()
  }
  invisible(TRUE)
}

.dsvert_dp_noise_mechanism <-
  "dsvert_dp_v1_deterministic_granular_laplace_int64"
.dsvert_dp_noise_implementation <-
  paste0("dsVert adapted Google Differential Privacy v4.1.0 ",
         "granular Laplace integer mechanism")
.dsvert_dp_gaussian_mechanism <-
  "dsvert_dp_v3_deterministic_approximate_gaussian_int64"
.dsvert_dp_gaussian_implementation <-
  paste0("dsVert adapted Google Differential Privacy v4.1.0 ",
         "approximate Gaussian integer mechanism with a published ",
         "2^-40 scalar total-variation bound")
.DSVERT_DP_NOISE_SELECTOR <- "minimum_conservative_95_radius_v3"
.DSVERT_DP_NOISE_TIE_BREAK <-
  "laplace_unless_gaussian_strictly_improves"
.DSVERT_DP_GAUSSIAN_TV_BOUND_PER_COORDINATE <- 2^-40
.DSVERT_DP_GAUSSIAN_ACCOUNTING_RULE <-
  "analytic_gaussian_delta_plus_dp_transfer_from_total_variation_bound"
.DSVERT_DP_GAUSSIAN_ACCURACY_RULE <-
  "gaussian_tail_alpha_minus_total_variation_union_bound"

.dsvert_dp_sampler_metadata_is_valid <- function(result, mechanism,
                                                 coordinate_count) {
  identical(result$mechanism, mechanism) &&
    identical(result$implementation, .dsvert_dp_noise_implementation) &&
    identical(result$sampler, "deterministic_two_sided_geometric") &&
    identical(result$randomness, "HMAC-SHA256/ChaCha20") &&
    .dsvert_dp_is_integer(result$clipped_coordinates, 0,
                          coordinate_count)
}

.dsvert_dp_google_ci_radius <- function(epsilon, sensitivity,
                                       alpha = NULL, confidence = NULL) {
  if (!is.null(alpha) && !is.null(confidence)) return(NA_real_)
  if (is.null(alpha)) {
    alpha <- if (is.null(confidence)) 0.05 else 1 - confidence
  }
  lambda <- sensitivity / epsilon
  if (!.dsvert_dp_is_number(
        epsilon, .DSVERT_DP_MINIMUM_EPSILON,
        .DSVERT_DP_MAXIMUM_EPSILON) ||
      !.dsvert_dp_is_number(alpha, 0, 1, lower_open = TRUE) || alpha >= 1 ||
      !.dsvert_dp_is_number(lambda, 0, 2^40, lower_open = TRUE)) {
    return(NA_real_)
  }
  floor(-lambda * log(alpha) + 0.5)
}

.dsvert_dp_gaussian_achieved_delta <- function(sigma, sensitivity,
                                                epsilon) {
  if (!.dsvert_dp_is_number(sigma, 0, lower_open = TRUE) ||
      !.dsvert_dp_is_number(sensitivity, 0, lower_open = TRUE) ||
      !.dsvert_dp_is_number(epsilon, 0, lower_open = TRUE)) return(NA_real_)
  exponential <- exp(epsilon)
  if (is.infinite(exponential)) return(0)
  a <- sensitivity / (2 * sigma)
  b <- epsilon * sigma / sensitivity
  achieved <- stats::pnorm(a - b) - exponential * stats::pnorm(-a - b)
  if (!is.finite(achieved) || achieved < -1e-15) return(NA_real_)
  max(0, achieved)
}

.dsvert_dp_gaussian_implementation_delta_bound <- function(
    coordinate_count, epsilon) {
  if (!.dsvert_dp_is_integer(
        coordinate_count, 1, .DSVERT_DP_MAX_COORDINATES) ||
      !.dsvert_dp_is_number(epsilon, 0, Inf, lower_open = TRUE)) {
    return(Inf)
  }
  log_bound <- log(coordinate_count) +
    log(.DSVERT_DP_GAUSSIAN_TV_BOUND_PER_COORDINATE) +
    epsilon + log1p(exp(-epsilon))
  if (!is.finite(log_bound) || log_bound >= 0) return(Inf)
  # The Go runtime rounds this DP transfer bound outwards. Keep client-side
  # planning at least as conservative so it cannot select a Gaussian release
  # that the server rejects only because of floating-point roundoff.
  bound <- exp(log_bound) * (1 + 64 * .Machine$double.eps)
  if (!is.finite(bound) || bound <= 0 || bound >= 1) Inf else bound
}

.dsvert_dp_gaussian_accuracy_radius <- function(sigma, alpha) {
  adjusted <- alpha - .DSVERT_DP_GAUSSIAN_TV_BOUND_PER_COORDINATE
  if (!.dsvert_dp_is_number(sigma, 0, lower_open = TRUE) ||
      !.dsvert_dp_is_number(adjusted, 0, 1, lower_open = TRUE)) {
    return(NA_real_)
  }
  radius <- ceiling(abs(stats::qnorm(adjusted / 2)) * sigma)
  if (!is.finite(radius) || radius < 0 || radius > 2^53 - 1) NA_real_ else
    radius
}

.dsvert_dp_validate_noise_selection_certificate <- function(
    certificate, coordinate_count, objective, laplace_epsilon,
    laplace_l1_sensitivity, gaussian_epsilon, gaussian_delta,
    gaussian_l2_sensitivity) {
  top_names <- c(
    "schema_version", "selector", "objective", "coordinate_count",
    "laplace", "gaussian", "winner", "winner_mechanism",
    "winning_metric_abs", "winner_delta", "tie_break")
  candidate_names <- c(
    "available", "mechanism", "epsilon", "delta", "analytic_delta",
    "implementation_delta_bound", "accounting_rule",
    "accuracy_accounting", "sensitivity_norm", "sensitivity",
    "marginal_95_abs", "simultaneous_95_abs", "nominal_rmse", "sigma",
    "granularity", "analytic_accounting_verified", "unavailable_reason")
  if (!.dsvert_dp_has_exact_names(certificate, top_names) ||
      !.dsvert_dp_is_integer(certificate$schema_version, 2, 2) ||
      !identical(certificate$selector, .DSVERT_DP_NOISE_SELECTOR) ||
      !identical(certificate$objective, objective) ||
      !.dsvert_dp_is_integer(
        certificate$coordinate_count, coordinate_count, coordinate_count) ||
      !identical(certificate$tie_break, .DSVERT_DP_NOISE_TIE_BREAK) ||
      !certificate$winner %in% c("laplace", "gaussian", "none") ||
      !.dsvert_dp_is_integer(certificate$winning_metric_abs, 0) ||
      !.dsvert_dp_is_number(certificate$winner_delta, 0, 1)) return(NULL)

  candidate_common <- function(value, mechanism, epsilon, delta,
                               norm, sensitivity) {
    .dsvert_dp_has_exact_names(value, candidate_names) &&
      is.logical(value$available) && length(value$available) == 1L &&
      !is.na(value$available) && identical(value$mechanism, mechanism) &&
      .dsvert_dp_is_number(value$epsilon, 0, Inf, TRUE) &&
      .dsvert_dp_num_equal(value$epsilon, epsilon, multiplier = 2048) &&
      .dsvert_dp_is_number(value$delta, 0, 1) &&
      .dsvert_dp_num_equal(value$delta, delta, multiplier = 2048) &&
      .dsvert_dp_is_number(value$analytic_delta, 0, 1) &&
      .dsvert_dp_is_number(value$implementation_delta_bound, 0, 1) &&
      .dsvert_dp_is_string(value$accounting_rule) &&
      .dsvert_dp_is_string(value$accuracy_accounting) &&
      identical(value$sensitivity_norm, norm) &&
      .dsvert_dp_is_number(value$sensitivity, 0, Inf, TRUE) &&
      .dsvert_dp_num_equal(value$sensitivity, sensitivity,
                           multiplier = 2048) &&
      .dsvert_dp_is_integer(value$marginal_95_abs, 0) &&
      .dsvert_dp_is_integer(value$simultaneous_95_abs, 0) &&
      value$simultaneous_95_abs >= value$marginal_95_abs &&
      .dsvert_dp_is_number(value$nominal_rmse, 0) &&
      .dsvert_dp_is_number(value$sigma, 0) &&
      .dsvert_dp_is_number(value$granularity, 0) &&
      is.logical(value$analytic_accounting_verified) &&
      length(value$analytic_accounting_verified) == 1L &&
      !is.na(value$analytic_accounting_verified) &&
      is.character(value$unavailable_reason) &&
      length(value$unavailable_reason) == 1L &&
      !is.na(value$unavailable_reason)
  }
  laplace <- certificate$laplace
  gaussian <- certificate$gaussian
  if (!candidate_common(
        laplace, .dsvert_dp_noise_mechanism, laplace_epsilon, 0,
        "l1", laplace_l1_sensitivity) ||
      !candidate_common(
        gaussian, .dsvert_dp_gaussian_mechanism, gaussian_epsilon,
        gaussian_delta, "l2", gaussian_l2_sensitivity)) return(NULL)

  laplace_available <- isTRUE(laplace$available)
  laplace_marginal <- .dsvert_dp_google_ci_radius(
    laplace_epsilon, laplace_l1_sensitivity)
  laplace_simultaneous <- .dsvert_dp_google_ci_radius(
    laplace_epsilon, laplace_l1_sensitivity,
    alpha = 0.05 / coordinate_count)
  laplace_granularity <- 2^ceiling(log2(
    (laplace_l1_sensitivity / laplace_epsilon) / 2^40))
  laplace_valid <-
    .dsvert_dp_num_equal(laplace$analytic_delta, 0) &&
    .dsvert_dp_num_equal(laplace$implementation_delta_bound, 0) &&
    identical(laplace$accounting_rule,
              "pure_dp_no_implementation_slack") &&
    identical(laplace$accuracy_accounting,
              "exact_granular_laplace_confidence_interval") &&
    if (laplace_available) {
      identical(laplace$analytic_accounting_verified, TRUE) &&
        identical(laplace$unavailable_reason, "") &&
        is.finite(laplace_marginal) && is.finite(laplace_simultaneous) &&
        .dsvert_dp_num_equal(laplace$marginal_95_abs, laplace_marginal) &&
        .dsvert_dp_num_equal(
          laplace$simultaneous_95_abs, laplace_simultaneous) &&
        .dsvert_dp_num_equal(
          laplace$nominal_rmse,
          sqrt(2) * laplace_l1_sensitivity / laplace_epsilon,
          multiplier = 4096) &&
        .dsvert_dp_num_equal(laplace$sigma, 0) &&
        laplace_granularity > 0 && laplace_granularity <= 1 &&
        .dsvert_dp_num_equal(
          laplace$granularity, laplace_granularity, multiplier = 4096)
    } else {
      identical(laplace$analytic_accounting_verified, FALSE) &&
        identical(
          laplace$unavailable_reason,
          "laplace_calibration_not_representable") &&
        all(vapply(c(
          "marginal_95_abs", "simultaneous_95_abs", "nominal_rmse",
          "sigma", "granularity"), function(name) {
            .dsvert_dp_num_equal(laplace[[name]], 0)
          }, logical(1L)))
    }
  if (!isTRUE(laplace_valid)) return(NULL)

  implementation_delta <-
    .dsvert_dp_gaussian_implementation_delta_bound(
      coordinate_count, gaussian_epsilon)
  gaussian_available <- isTRUE(gaussian$available)
  gaussian_valid <-
    .dsvert_dp_num_equal(
      gaussian$implementation_delta_bound, implementation_delta,
      multiplier = 2048) &&
    identical(gaussian$accounting_rule,
              .DSVERT_DP_GAUSSIAN_ACCOUNTING_RULE) &&
    identical(gaussian$accuracy_accounting,
              .DSVERT_DP_GAUSSIAN_ACCURACY_RULE) &&
    if (gaussian_available) {
      analytic_delta <- gaussian_delta - implementation_delta
      achieved_delta <- .dsvert_dp_gaussian_achieved_delta(
        gaussian$sigma, gaussian_l2_sensitivity, gaussian_epsilon)
      calibrated_sigma <- .dsvert_dp_gaussian_sigma_matches_calibration(
        gaussian$sigma, gaussian_l2_sensitivity, gaussian_epsilon,
        analytic_delta)
      marginal <- .dsvert_dp_gaussian_accuracy_radius(
        gaussian$sigma, 0.05)
      simultaneous <- .dsvert_dp_gaussian_accuracy_radius(
        gaussian$sigma, 0.05 / coordinate_count)
      granularity <- 2^ceiling(log2(gaussian$sigma / 2^56))
      analytic_delta > 0 &&
        .dsvert_dp_num_equal(
          gaussian$analytic_delta, analytic_delta, multiplier = 2048) &&
        gaussian$analytic_delta + gaussian$implementation_delta_bound <=
          gaussian$delta &&
        is.finite(achieved_delta) &&
        achieved_delta <= gaussian$analytic_delta &&
        isTRUE(calibrated_sigma) &&
        identical(gaussian$analytic_accounting_verified, TRUE) &&
        identical(gaussian$unavailable_reason, "") &&
        .dsvert_dp_num_equal(gaussian$nominal_rmse, gaussian$sigma,
                             multiplier = 2048) &&
        granularity > 0 && granularity <= 1 &&
        .dsvert_dp_num_equal(
          gaussian$granularity, granularity, multiplier = 4096) &&
        .dsvert_dp_num_equal(gaussian$marginal_95_abs, marginal) &&
        .dsvert_dp_num_equal(
          gaussian$simultaneous_95_abs, simultaneous)
    } else {
      identical(gaussian$analytic_accounting_verified, FALSE) &&
        .dsvert_dp_num_equal(gaussian$analytic_delta, 0) &&
        gaussian$unavailable_reason %in% c(
          "gaussian_delta_is_zero",
          "gaussian_implementation_bound_not_representable",
          "gaussian_delta_does_not_cover_implementation_bound",
          "gaussian_calibration_not_representable",
          "gaussian_accuracy_not_representable",
          "gaussian_analytic_accounting_not_verified") &&
        all(vapply(c(
          "marginal_95_abs", "simultaneous_95_abs", "nominal_rmse",
          "sigma", "granularity"), function(name) {
            .dsvert_dp_num_equal(gaussian[[name]], 0)
          }, logical(1L)))
    }
  if (!isTRUE(gaussian_valid)) return(NULL)

  metric <- if (identical(objective, "marginal_95_abs")) {
    "marginal_95_abs"
  } else {
    "simultaneous_95_abs"
  }
  expected_winner <- if (gaussian_available &&
      (!laplace_available || gaussian[[metric]] < laplace[[metric]])) {
    "gaussian"
  } else if (laplace_available) {
    "laplace"
  } else {
    "none"
  }
  winner_valid <- if (identical(expected_winner, "none")) {
    identical(certificate$winner_mechanism, "") &&
      .dsvert_dp_num_equal(certificate$winning_metric_abs, 0) &&
      .dsvert_dp_num_equal(certificate$winner_delta, 0)
  } else {
    selected <- certificate[[expected_winner]]
    identical(certificate$winner_mechanism, selected$mechanism) &&
      .dsvert_dp_num_equal(
        certificate$winning_metric_abs, selected[[metric]]) &&
      .dsvert_dp_num_equal(certificate$winner_delta, selected$delta)
  }
  if (!identical(certificate$winner, expected_winner) || !winner_valid) {
    return(NULL)
  }
  list(winner = expected_winner,
       selected = if (identical(expected_winner, "none")) NULL else
         certificate[[expected_winner]],
       certificate = certificate)
}

.dsvert_dp_selected_sampler_metadata_is_valid <- function(
    result, selection, coordinate_count) {
  if (!is.list(selection) ||
      !selection$winner %in% c("laplace", "gaussian") ||
      !is.list(selection$selected) ||
      !identical(result$mechanism, selection$selected$mechanism) ||
      !identical(result$randomness, "HMAC-SHA256/ChaCha20") ||
      !.dsvert_dp_is_integer(
        result$clipped_coordinates, 0, coordinate_count)) return(FALSE)
  if (identical(selection$winner, "laplace")) {
    identical(result$implementation, .dsvert_dp_noise_implementation) &&
      identical(result$sampler, "deterministic_two_sided_geometric")
  } else {
    identical(result$implementation, .dsvert_dp_gaussian_implementation) &&
      identical(result$sampler, "deterministic_symmetric_binomial")
  }
}

.dsvert_dp_validate_count_release <- function(result, policy, invalid) {
  if (identical(policy$adjacency, "replace_one_fixed_cohort")) {
    expected <- c(
      "released", "value", "mechanism", "implementation", "sampler",
      "randomness", "sensitivity", "postprocessing",
      "clipped_coordinates", "accuracy_95_abs", "data_dependency",
      "epsilon", "delta", "adjacency", "peer_count")
    if (!.dsvert_dp_has_exact_names(result, expected) ||
        !.dsvert_dp_is_integer(result$value, policy$fixed_cohort_size,
                               policy$fixed_cohort_size) ||
        !identical(result$mechanism, "public_fixed_cohort_size_v1") ||
        !identical(result$implementation,
                   "custodian_owned_policy_constant") ||
        !identical(result$sampler, "none") ||
        !identical(result$randomness, "none") ||
        !.dsvert_dp_num_equal(result$sensitivity, 0) ||
        !identical(result$postprocessing, "none_public_policy_value") ||
        !.dsvert_dp_is_integer(result$clipped_coordinates, 0, 0) ||
        !.dsvert_dp_is_integer(result$accuracy_95_abs, 0, 0) ||
        !identical(result$data_dependency,
                   "none_public_fixed_cohort_policy") ||
        !.dsvert_dp_num_equal(result$epsilon, 0) ||
        !.dsvert_dp_num_equal(result$delta, 0) ||
        !.dsvert_dp_is_integer(result$peer_count, 2,
                               .Machine$integer.max) ||
        !identical(as.numeric(result$peer_count),
                   as.numeric(length(policy$peer_pinset)))) invalid()
    return(invisible(TRUE))
  }
  expected <- c(
    "released", "value", "mechanism", "implementation", "sampler",
    "randomness", "sensitivity", "l1_sensitivity", "l2_sensitivity",
    "noise_selection", "postprocessing", "clipped_coordinates",
    "accuracy_95_abs", "privacy_epoch", "noise_key_id", "sticky_noise",
    "epsilon", "delta", "adjacency",
    "composition_partitions")
  scheduled_delta <- result$epsilon *
    policy$local_total_delta / policy$local_total_epsilon
  selection <- .dsvert_dp_validate_noise_selection_certificate(
    result$noise_selection, 1L, "marginal_95_abs",
    result$epsilon, 1, result$epsilon, scheduled_delta, 1)
  if (!.dsvert_dp_has_exact_names(result, expected) ||
      !.dsvert_dp_is_integer(result$value) ||
      is.null(selection) ||
      !.dsvert_dp_selected_sampler_metadata_is_valid(
        result, selection, 1L) ||
      !.dsvert_dp_num_equal(result$sensitivity, 1) ||
      !.dsvert_dp_num_equal(result$l1_sensitivity, 1) ||
      !.dsvert_dp_num_equal(result$l2_sensitivity, 1) ||
      !.dsvert_dp_num_equal(result$delta, selection$selected$delta) ||
      !identical(result$postprocessing, "nonnegative_integer") ||
      !.dsvert_dp_is_integer(result$accuracy_95_abs) ||
      !.dsvert_dp_num_equal(
        result$accuracy_95_abs, selection$selected$marginal_95_abs) ||
      !identical(policy$adjacency, "add_remove_patient")) invalid()
}

.dsvert_dp_validate_table_release <- function(result, policy, row_var,
                                               col_var, invalid) {
  expected <- c(
    "released", "counts", "row_levels", "col_levels", "nrow", "ncol",
    "unit_aggregation_policy",
    "mechanism", "implementation", "sampler", "randomness", "sensitivity",
    "l1_sensitivity", "l2_sensitivity", "noise_selection",
    "postprocessing", "clipped_coordinates", "accuracy_95_abs_per_cell",
    "accuracy_simultaneous_95_abs", "accuracy_simultaneous_confidence",
    "accuracy_simultaneous_method",
    "privacy_epoch", "noise_key_id", "sticky_noise",
    "epsilon", "delta", "adjacency",
    "composition_partitions")
  row_levels <- policy$categorical_levels[[row_var]]
  col_levels <- policy$categorical_levels[[col_var]]
  if (!.dsvert_dp_has_exact_names(result, expected) ||
      is.null(row_levels) || is.null(col_levels)) invalid()
  n_row <- length(row_levels)
  n_col <- length(col_levels)
  if (n_row < 1L || n_col < 1L ||
      n_row > .DSVERT_DP_MAX_COORDINATES ||
      n_col > floor(.DSVERT_DP_MAX_COORDINATES / n_row)) invalid()
  coordinate_count <- as.integer(n_row * n_col)
  sensitivity <- if (identical(
    policy$adjacency, "add_remove_patient")) 1L else 2L
  l2_sensitivity <- sqrt(sensitivity)
  scheduled_delta <- result$epsilon *
    policy$local_total_delta / policy$local_total_epsilon
  selection <- .dsvert_dp_validate_noise_selection_certificate(
    result$noise_selection, coordinate_count, "simultaneous_95_abs",
    result$epsilon, sensitivity, result$epsilon, scheduled_delta,
    l2_sensitivity)
  if (!identical(result$row_levels, row_levels) ||
      !identical(result$col_levels, col_levels) ||
      !identical(
        result$unit_aggregation_policy,
        policy$contingency_unit_aggregation_policy) ||
      !.dsvert_dp_is_integer(result$nrow, 1, .Machine$integer.max) ||
      !.dsvert_dp_is_integer(result$ncol, 1, .Machine$integer.max) ||
      result$nrow != length(row_levels) ||
      result$ncol != length(col_levels) ||
      !is.numeric(result$counts) ||
      length(result$counts) != length(row_levels) * length(col_levels) ||
      anyNA(result$counts) || any(!is.finite(result$counts)) ||
      any(result$counts < 0) || any(result$counts != floor(result$counts)) ||
      any(result$counts > 2^53 - 1) ||
      is.null(selection) ||
      !.dsvert_dp_selected_sampler_metadata_is_valid(
        result, selection, coordinate_count) ||
      !.dsvert_dp_num_equal(result$sensitivity, sensitivity) ||
      !.dsvert_dp_num_equal(result$l1_sensitivity, sensitivity) ||
      !.dsvert_dp_num_equal(result$l2_sensitivity, l2_sensitivity,
                            multiplier = 2048) ||
      !.dsvert_dp_num_equal(result$delta, selection$selected$delta) ||
      !identical(result$postprocessing, "cellwise_nonnegative_integer") ||
      !is.numeric(result$accuracy_95_abs_per_cell) ||
      length(result$accuracy_95_abs_per_cell) != coordinate_count ||
      anyNA(result$accuracy_95_abs_per_cell) ||
      any(!is.finite(result$accuracy_95_abs_per_cell)) ||
      any(result$accuracy_95_abs_per_cell < 0) ||
      any(result$accuracy_95_abs_per_cell !=
            floor(result$accuracy_95_abs_per_cell)) ||
      any(result$accuracy_95_abs_per_cell > 2^53 - 1) ||
      any(result$accuracy_95_abs_per_cell !=
            selection$selected$marginal_95_abs) ||
      !.dsvert_dp_is_integer(result$accuracy_simultaneous_95_abs) ||
      !.dsvert_dp_num_equal(result$accuracy_simultaneous_95_abs,
                            selection$selected$simultaneous_95_abs) ||
      result$accuracy_simultaneous_95_abs <
        max(result$accuracy_95_abs_per_cell) ||
      !.dsvert_dp_num_equal(result$accuracy_simultaneous_confidence, 0.95) ||
      !identical(result$accuracy_simultaneous_method, "union_bound")) {
    invalid()
  }
}

.dsvert_dp_validate_meanvar_release <- function(result, policy, variable,
                                                invalid) {
  expected <- c(
    "released", "status", "reason", "n", "n_definition",
    "effective_count_95_lower_bound", "mean", "variance",
    "lower_bound", "upper_bound",
    "variance_definition",
    "mechanism", "implementation", "sampler", "randomness",
    "epsilon_allocation", "submechanism_count", "composition_rule",
    "coordinate_l1_sensitivity", "l2_sensitivity", "noise_selection",
    "numeric_grid_bits", "numeric_grid_scale",
    "max_abs_quantization_per_unit", "postprocessing",
    "clipped_coordinates",
    "accuracy_95_abs_count",
    "accuracy_95_abs_normalized_sum_noise_only",
    "accuracy_95_abs_normalized_sumsq_noise_only",
    "privacy_epoch", "noise_key_id", "sticky_noise", "epsilon", "delta",
    "adjacency", "composition_partitions")
  bounds <- policy$numeric_bounds[[variable]]
  if (!.dsvert_dp_has_exact_names(result, expected) || is.null(bounds)) {
    invalid()
  }
  grid_scale <- 2^policy$numeric_grid_bits
  epsilon_coordinate <- result$epsilon / 3
  expected_sensitivity <- c(1, grid_scale, grid_scale)
  l2_sensitivity <- sqrt(1 + 2 * grid_scale^2)
  scheduled_delta <- result$epsilon *
    policy$local_total_delta / policy$local_total_epsilon
  selection <- .dsvert_dp_validate_noise_selection_certificate(
    result$noise_selection, 3L, "simultaneous_95_abs",
    epsilon_coordinate, grid_scale, result$epsilon, scheduled_delta,
    l2_sensitivity)
  gaussian_selected <- !is.null(selection) &&
    identical(selection$winner, "gaussian")
  count_radius <- if (gaussian_selected) {
    selection$selected$marginal_95_abs
  } else {
    .dsvert_dp_google_ci_radius(epsilon_coordinate, 1)
  }
  moment_radius <- if (gaussian_selected) {
    selection$selected$marginal_95_abs / grid_scale
  } else {
    .dsvert_dp_google_ci_radius(
      epsilon_coordinate, grid_scale) / grid_scale
  }
  numeric_vector_equal <- function(value, expected) {
    is.numeric(value) && length(value) == length(expected) &&
      !anyNA(value) && all(is.finite(value)) &&
      all(vapply(seq_along(value), function(index) {
        .dsvert_dp_num_equal(value[[index]], expected[[index]],
                             multiplier = 1024)
      }, logical(1L)))
  }
  width <- bounds[[2L]] - bounds[[1L]]
  scale <- max(1, abs(bounds), width^2)
  tolerance <- 1024 * .Machine$double.eps * scale
  status_ok <- identical(result$status, "ok")
  status_point <- identical(
    result$status, "dp_point_available_count_not_certified_positive")
  status_nonestimable <- identical(
    result$status, "dp_effective_count_not_certified_positive")
  expected_count_lower <- if (.dsvert_dp_is_integer(result$n, 0) &&
      .dsvert_dp_is_integer(result$accuracy_95_abs_count, 0)) {
    max(0, result$n - result$accuracy_95_abs_count)
  } else {
    NA_real_
  }
  point_geometry <- FALSE
  if ((status_ok || status_point) &&
      ((status_ok && is.null(result$reason)) ||
       (status_point && identical(
         result$reason,
         "dp_noisy_effective_count_lower_bound_is_zero"))) &&
      .dsvert_dp_is_number(result$mean) &&
      .dsvert_dp_is_number(result$variance, 0)) {
    feasible_variance <-
      (result$mean - bounds[[1L]]) * (bounds[[2L]] - result$mean)
    point_geometry <-
      result$n > 0 &&
      ((status_ok && result$effective_count_95_lower_bound > 0) ||
       (status_point &&
          identical(as.numeric(result$effective_count_95_lower_bound), 0))) &&
      result$mean >= bounds[[1L]] - tolerance &&
      result$mean <= bounds[[2L]] + tolerance &&
      result$variance <= feasible_variance + tolerance
  }
  nonestimable_geometry <- status_nonestimable &&
    identical(result$reason,
              "dp_noisy_effective_count_lower_bound_is_zero") &&
    is.null(result$mean) && is.null(result$variance) &&
    identical(as.numeric(result$n), 0) &&
    identical(as.numeric(result$effective_count_95_lower_bound), 0)
  if (!(status_ok || status_point || status_nonestimable) ||
      !.dsvert_dp_is_integer(result$n, 0) ||
      !identical(
        result$n_definition,
        "nonnegative_postprocessed_dp_noisy_effective_unit_count") ||
      !.dsvert_dp_is_integer(result$effective_count_95_lower_bound, 0) ||
      !.dsvert_dp_num_equal(
        result$effective_count_95_lower_bound, expected_count_lower) ||
      !(point_geometry || nonestimable_geometry) ||
      !identical(
        result$variance_definition,
        "population_central_second_moment_denominator_n") ||
      !.dsvert_dp_num_equal(result$lower_bound, bounds[[1L]]) ||
      !.dsvert_dp_num_equal(result$upper_bound, bounds[[2L]]) ||
      is.null(selection) ||
      !.dsvert_dp_selected_sampler_metadata_is_valid(
        result, selection, 3L) ||
      !numeric_vector_equal(
        result$epsilon_allocation,
        if (gaussian_selected) result$epsilon else
          rep(epsilon_coordinate, 3L)) ||
      !.dsvert_dp_is_integer(
        result$submechanism_count,
        if (gaussian_selected) 1 else 3,
        if (gaussian_selected) 1 else 3) ||
      !identical(
        result$composition_rule,
        if (gaussian_selected) "single_joint_l2_gaussian_vector" else
          "three_sequential_laplace_coordinates") ||
      !numeric_vector_equal(result$coordinate_l1_sensitivity,
                            expected_sensitivity) ||
      !.dsvert_dp_num_equal(result$l2_sensitivity, l2_sensitivity,
                            multiplier = 2048) ||
      !.dsvert_dp_num_equal(result$delta, selection$selected$delta) ||
      !.dsvert_dp_is_integer(result$numeric_grid_bits, 8, 18) ||
      result$numeric_grid_bits != policy$numeric_grid_bits ||
      !.dsvert_dp_num_equal(result$numeric_grid_scale, grid_scale) ||
      !.dsvert_dp_num_equal(result$max_abs_quantization_per_unit,
                            0.5 / grid_scale) ||
      !identical(result$postprocessing,
                 "bounded_consistent_quantized_moments") ||
      !.dsvert_dp_is_integer(result$accuracy_95_abs_count) ||
      !.dsvert_dp_num_equal(result$accuracy_95_abs_count, count_radius) ||
      !.dsvert_dp_is_number(
        result$accuracy_95_abs_normalized_sum_noise_only, 0) ||
      !.dsvert_dp_num_equal(
        result$accuracy_95_abs_normalized_sum_noise_only, moment_radius,
        multiplier = 1024) ||
      !.dsvert_dp_is_number(
        result$accuracy_95_abs_normalized_sumsq_noise_only, 0) ||
      !.dsvert_dp_num_equal(
        result$accuracy_95_abs_normalized_sumsq_noise_only,
        moment_radius,
        multiplier = 1024)) invalid()
}

.dsvert_dp_validate_release <- function(result, server, status,
                                        kind = c("count", "table", "meanvar"),
                                        row_var = NULL, col_var = NULL,
                                        variable = NULL) {
  kind <- match.arg(kind)
  invalid <- function() {
    stop("Server '", server, "' returned an invalid DP v7 release",
         call. = FALSE)
  }
  if (!is.list(result)) invalid()
  if (length(result) == 1L && is.list(result[[1L]]) &&
      !"released" %in% names(result)) {
    if (is.null(names(result)) || !identical(names(result), server)) invalid()
    result <- result[[1L]]
  }
  if (!is.list(result) || !is.list(status) || !is.list(status$policy) ||
      !"released" %in% names(result)) invalid()

  if (identical(result$released, FALSE)) {
    expected <- c("released", "value", "reason", "epsilon", "delta",
                  "adjacency", "composition_partitions")
    if (identical(kind, "table")) {
      expected <- c(expected, "unit_aggregation_policy")
    }
    if (!.dsvert_dp_has_exact_names(result, expected) ||
        !is.null(result$value) ||
        !identical(result$reason,
                   "privacy_allocation_not_representable") ||
        (identical(kind, "table") && !identical(
          result$unit_aggregation_policy,
          status$policy$contingency_unit_aggregation_policy))) invalid()
    .dsvert_dp_validate_release_common(result, status, server)
  } else {
    if (!identical(result$released, TRUE)) invalid()
    .dsvert_dp_validate_release_common(result, status, server)
    switch(kind,
      count = .dsvert_dp_validate_count_release(
        result, status$policy, invalid),
      table = .dsvert_dp_validate_table_release(
        result, status$policy, row_var, col_var, invalid),
      meanvar = .dsvert_dp_validate_meanvar_release(
        result, status$policy, variable, invalid))
  }
  result$server <- server
  result
}

.dsvert_dp_vector_context <- function(run) {
  if (!is.list(run) || !inherits(run$release, "dsvert_joint_dp_vector") ||
      !is.list(run$layout) || !is.list(run$status) || !length(run$status) ||
      is.null(names(run$status)) || anyNA(names(run$status)) ||
      anyDuplicated(names(run$status)) ||
      !identical(as.numeric(run$layout$coordinate_count),
                 as.numeric(run$release$coordinate_count)) ||
      !is.list(run$release$manifest) ||
      !identical(run$release$history_gate, TRUE) ||
      !identical(run$release$request_limit, FALSE) ||
      !identical(run$release$operation_limit, TRUE)) {
    stop("The joint biomedical DP vector context is invalid", call. = FALSE)
  }
  adjacency <- tryCatch(vapply(run$status, function(value) {
    value$policy$adjacency
  }, character(1L)), error = function(error) character())
  if (length(adjacency) != length(run$status) ||
      !all(adjacency %in% c(
        "add_remove_patient", "replace_one_fixed_cohort")) ||
      length(unique(adjacency)) != 1L) {
    stop("The joint biomedical DP vector has inconsistent adjacency",
         call. = FALSE)
  }
  lattice <- run$release$manifest$workload$release_lattice
  if (!is.list(lattice) ||
      !.dsvert_dp_is_number(lattice$natural_l1_sensitivity, 0, Inf,
                            lower_open = TRUE) ||
      !.dsvert_dp_is_number(lattice$natural_l2_sensitivity, 0, Inf,
                            lower_open = TRUE) ||
      !.dsvert_dp_is_integer(lattice$output_lattice_bits, 1, 62) ||
      !identical(as.numeric(lattice$output_lattice_scale),
                 2^as.numeric(lattice$output_lattice_bits))) {
    stop("The signed biomedical DP release lattice is invalid",
         call. = FALSE)
  }
  list(release = run$release, layout = run$layout, status = run$status,
       manifest_bundle = run$manifest_bundle,
       manifest = run$release$manifest, lattice = lattice,
       adjacency = adjacency[[1L]])
}

.dsvert_dp_capsule_security_claim <- function() {
  list(
    version = "dsvert-capsule-security-claim-v3",
    privacy_definition = "bounded_lifetime_epsilon_delta_dp",
    operation_accounting =
      "one_per_distinct_capsule_allocator_commit",
    one_public_release_instance_per_capsule = TRUE,
    history_can_deny_new_capsule = TRUE,
    request_limit = FALSE,
    same_release_replay_is_postprocessing = TRUE,
    observable_scope =
      paste0(
        "successful_authenticated_semantic_messages_public_shape_or_",
        "dp_output_postprocessing_only"),
    timing_availability_and_traffic_flow_in_scope = FALSE,
    analyst_relay_trusted = FALSE,
    relay_tamper_behavior =
      "signature_or_AEAD_failure_aborts_without_an_accepted_result",
    peer_execution_model =
      "authenticated_protocol_compliant_semi_honest",
    designated_noise_peer_count = 2L,
    minimum_noncolluding_designated_noise_peers = 1L,
    authenticated_history_retention_assumption =
      paste0("at_least_one_noncolluding_designated_noise_peer_retains_and_",
             "uses_complete_authenticated_monotonic_history"),
    privacy_accountant_namespace_assumption = paste0(
      "one_stable_unique_namespace_across_domain_cohort_policy_pinset_",
      "and_ledger_reconfiguration_per_protected_privacy_universe"),
    simultaneous_designated_history_rollback_protection =
      "not_claimed_without_external_linearizable_cas",
    malicious_peer_security = FALSE,
    host_compromise_in_scope = FALSE,
    all_designated_noise_peer_collusion_in_scope = FALSE,
    unconditional_non_reconstruction_guarantee = FALSE)
}

.dsvert_dp_vector_public_metadata <- function(context) {
  release <- context$release
  lattice <- context$lattice
  profile <- .dsvert_dp_vector_release_profile(release, context$manifest)
  epochs <- vapply(context$status, function(value) {
    root <- value$noise_root
    epoch <- if (is.list(root) && !is.null(root$epoch)) {
      root$epoch
    } else if (is.list(root)) {
      root$privacy_epoch
    } else NULL
    if (!is.numeric(epoch) || length(epoch) != 1L || is.na(epoch) ||
        !is.finite(epoch) || epoch < 1 || epoch != floor(epoch)) NA_real_ else
      as.numeric(epoch)
  }, numeric(1L))
  key_ids <- vapply(context$status, function(value) {
    key_id <- tryCatch(value$noise_root$key_id, error = function(error) NULL)
    if (.dsvert_dp_is_string(key_id)) key_id else NA_character_
  }, character(1L))
  list(
    released = TRUE,
    mechanism = release$mechanism,
    implementation = profile$backend,
    backend = "exact_signed_Ring128_global_vector",
    sampler = profile$sampler,
    randomness = paste(
      "independent pinned-peer HKDF-SHA256/ChaCha20 streams;",
      "no analyst-controlled seed"),
    postprocessing = profile$postprocessing,
    one_joint_draw = isTRUE(profile$exact_gc),
    mechanism_plan = release$mechanism_plan,
    plan_sha256 = release$plan_sha256,
    backend_selection = release$backend_selection,
    backend_assessment = release$backend_assessment,
    epsilon = as.numeric(release$epsilon),
    delta = as.numeric(release$delta),
    implementation_delta = release$implementation_delta,
    delta_aggregation = release$delta_aggregation,
    adjacency = context$adjacency,
    sensitivity = as.numeric(if (isTRUE(profile$gaussian)) {
      lattice$natural_l2_sensitivity
    } else {
      lattice$natural_l1_sensitivity
    }),
    sensitivity_norm = if (isTRUE(profile$gaussian)) "l2" else "l1",
    l1_sensitivity = as.numeric(lattice$natural_l1_sensitivity),
    l2_sensitivity = as.numeric(lattice$natural_l2_sensitivity),
    sensitivity_scope = "complete_signed_biomedical_capsule_vector",
    capsule_mechanism = context$manifest$workload$capsule_mechanism,
    mechanism_selection = context$manifest$workload$mechanism_selection,
    output_lattice_bits = as.integer(lattice$output_lattice_bits),
    output_lattice_scale = as.numeric(lattice$output_lattice_scale),
    capsule_id = release$capsule_id,
    manifest_sha256 = release$manifest_sha256,
    final_vector_root = release$final_vector_root,
    coordinate_order_sha256 = release$coordinate_order_sha256,
    capsule_coordinate_count = as.integer(release$coordinate_count),
    sticky_noise = "immutable_capsule_durable_replay_v3",
    sticky_replay = isTRUE(release$sticky_replay),
    privacy_epochs = epochs,
    noise_key_ids = key_ids,
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    history_gate = TRUE,
    request_limit = FALSE,
    operation_limit = TRUE,
    composition_rule = "one_full_global_vector_release_no_coordinate_split",
    data_dependency = "immutable_signed_biomedical_capsule_snapshot",
    security_claim = .dsvert_dp_capsule_security_claim(),
    clipped_coordinates = NA_integer_,
    clamp_activation_disclosed = FALSE)
}

.dsvert_dp_vector_server_filter <- function(server, datasources) {
  if (is.null(server)) return(NULL)
  .dsvert_dp_server(server, datasources)
}

.dsvert_dp_vector_block_capacity <- function(block, numeric_moment = FALSE) {
  maximum <- tryCatch(block$descriptor$statistic_maximum,
                      error = function(error) NULL)
  maximum <- .dsvert_dp_capsule_manifest_numbers(
    maximum, "coordinate maximum")
  expected <- if (isTRUE(numeric_moment)) 3L else 1L
  if (length(maximum) != expected || any(maximum < 1) ||
      any(maximum != floor(maximum)) || any(maximum > 2^53 - 1)) {
    stop("Invalid signed capsule coordinate maximum", call. = FALSE)
  }
  if (isTRUE(numeric_moment)) {
    scale <- 2^as.numeric(block$descriptor$numeric_grid_bits)
    if (!identical(maximum[[2L]], maximum[[1L]] * scale) ||
        !identical(maximum[[3L]], maximum[[1L]] * scale)) {
      stop("Invalid signed numeric-moment coordinate maximum",
           call. = FALSE)
    }
  }
  as.numeric(maximum[[1L]])
}

#' Differentially private patient count
#'
#' Under add/remove adjacency this is local post-processing of the admitted-
#' count coordinate in the one signed, sticky biomedical capsule vector. Two
#' designated pinned peers contribute independent complete-vector geometric
#' draws; Ring128 decoding and the public coordinate clamp happen before the
#' client sees the vector. Under fixed-cohort adjacency, the only exception is
#' a public policy constant whose declared sensitivity, epsilon and delta are
#' all zero.
#'
#' @param data_name Name of the protected data frame.
#' @param server Optional returned-source label for the fixed-cohort
#'   public-policy path; every connected peer must return the same value and it
#'   must equal the consensus policy capacity. On the joint add/remove path it
#'   is only a compatibility assertion; the signed manifest selects the
#'   coordinate owner.
#' @param datasources DataSHIELD connections.
#' @return A DP count release with accuracy and composition metadata.
#' @export
ds.vertDPCount <- function(data_name, server = NULL, datasources = NULL) {
  .dsvert_dp_count_impl(
    data_name, server, datasources, DSI::datashield.aggregate)
}

.dsvert_dp_count_impl <- function(data_name, server = NULL,
                                  datasources = NULL, .aggregate) {
  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) || !nzchar(data_name)) {
    stop("data_name must be one non-empty string", call. = FALSE)
  }
  datasources <- .dsvert_dp_datasources(datasources)
  status <- .dsvert_joint_dp_capsule_status_impl(datasources, .aggregate)
  server <- .dsvert_dp_server(server, datasources)
  joint_add_remove <- all(vapply(status, function(value) {
    identical(value$policy$adjacency, "add_remove_patient")
  }, logical(1L)))
  if (isTRUE(joint_add_remove)) {
    # `server` is retained as a compatibility assertion only.  The signed
    # capsule descriptor, not the analyst, chooses the count owner.
    if (!is.null(server)) .dsvert_dp_server(server, datasources)
    run <- .dsvert_dp_capsule_vector_run(
      datasources, status = status, .aggregate = .aggregate)
    context <- .dsvert_dp_vector_context(run)
    block <- .dsvert_dp_capsule_single_block(
      context$layout, "admitted_count", dataset = data_name,
      description = paste0("signed admitted-count block for dataset '",
                           data_name, "'"))
    value <- .dsvert_dp_capsule_vector_values(context$release, block)
    capacity <- .dsvert_dp_vector_block_capacity(block)
    if (length(value) != 1L || value[[1L]] < 0 || value[[1L]] > capacity) {
      stop("The released admitted count violates its signed bounds",
           call. = FALSE)
    }
    accuracy <- .dsvert_dp_vector_accuracy_radius(
      context$release, context$manifest, coordinate_count = 1L,
      confidence = 0.95, maximum_error = capacity)
    result <- c(.dsvert_dp_vector_public_metadata(context), list(
      value = unname(value[[1L]]), server = block$owner_peer,
      coordinate_family = "admitted_count",
      coordinate_descriptor = block$descriptor,
      lower_bound = 0, upper_bound = capacity,
      artifact_l1_sensitivity = 1,
      accuracy_95_abs = accuracy$radius,
      accuracy_95_confidence = accuracy$confidence,
      accuracy_95_method = accuracy$method,
      accuracy_implementation_tv_upper_bound =
        accuracy$implementation_tv_upper_bound,
      accuracy_additional_privacy_cost =
        accuracy$additional_privacy_cost))
  } else {
    if (!all(vapply(status, function(value) {
      identical(value$policy$adjacency, "replace_one_fixed_cohort")
    }, logical(1L)))) {
      stop("Connected peers do not share a supported Count adjacency",
           call. = FALSE)
    }
    response <- .dsvert_aggregate_strict(
      conns = datasources,
      expr = call(
        name = "dsvertPublicFixedCohortCountDS", data_name = data_name),
      operation = "DP count release", .aggregate = .aggregate)
    fields <- c(
      "released", "value", "mechanism", "implementation", "sampler",
      "randomness", "sensitivity", "postprocessing",
      "clipped_coordinates", "accuracy_95_abs", "data_dependency",
      "epsilon", "delta", "adjacency", "peer_count")
    valid <- vapply(names(datasources), function(peer) {
      result <- response[[peer]]
      policy <- status[[peer]]$policy
      is.list(result) && !is.null(names(result)) &&
        !anyNA(names(result)) && !anyDuplicated(names(result)) &&
        setequal(names(result), fields) && identical(result$released, TRUE) &&
        .dsvert_dp_is_integer(
          result$value, policy$unit_capacity, policy$unit_capacity) &&
        identical(result$mechanism, "public_fixed_cohort_size_v1") &&
        identical(result$implementation,
                  "custodian_owned_policy_constant") &&
        identical(result$sampler, "none") &&
        identical(result$randomness, "none") &&
        identical(as.numeric(result$sensitivity), 0) &&
        identical(result$postprocessing, "none_public_policy_value") &&
        identical(as.numeric(result$clipped_coordinates), 0) &&
        identical(as.numeric(result$accuracy_95_abs), 0) &&
        identical(result$data_dependency,
                  "none_public_fixed_cohort_policy") &&
        identical(as.numeric(result$epsilon), 0) &&
        identical(as.numeric(result$delta), 0) &&
        identical(result$adjacency, "replace_one_fixed_cohort") &&
        .dsvert_dp_is_integer(
          result$peer_count, 2, .Machine$integer.max) &&
        identical(as.numeric(result$peer_count),
                  as.numeric(policy$peer_count))
    }, logical(1L))
    valid <- all(valid) &&
      length(unique(vapply(response, `[[`, numeric(1L), "value"))) == 1L
    if (!isTRUE(valid)) {
      stop("The peers returned an invalid public fixed-cohort Count",
           call. = FALSE)
    }
    result <- response[[server]]
    result$server <- server
    result$history_gate <- FALSE
    result$request_limit <- FALSE
    result$operation_limit <- FALSE
    result$composition_rule <- "public_policy_constant_no_dp_composition"
  }
  result$uncertainty_scope <- if (isTRUE(result$released) &&
      identical(result$mechanism, "public_fixed_cohort_size_v1")) {
    paste(
      "Public fixed-cohort policy value; no DP mechanism noise;",
      "sampling uncertainty excluded")
  } else {
    "DP mechanism noise only; sampling uncertainty excluded"
  }
  result$inferential_scope <- paste(
    "Finite-dataset DP release; no sampling confidence interval or",
    "population-supermodel inference is provided")
  class(result) <- c("ds.vertDPCount", "list")
  result
}

#' Differentially private fixed-domain contingency table
#'
#' The table domain, orientation and one-contribution-per-unit rule come from
#' the signed capsule descriptor. This function only reshapes the already-DP
#' block in the global sticky vector; it does not apply an
#' ordinary chi-square or Fisher p-value, because those reference
#' distributions are not calibrated for DP-noised counts. Under the current
#' repeated-record contract, a unit contributes once only when all its valid
#' fixed-domain rows occupy the same cell. Conflicting valid cells contribute
#' zero, while missing and out-of-domain rows are ignored for consistency.
#' This public, custodian-owned rule is returned with every release.
#'
#' @param data_name Name of the protected data frame.
#' @param row_var,col_var Fixed-domain categorical variables. They may be held
#'   by one custodian or by two different custodians named in a signed
#'   `vertical_cross_specs` `categorical_pair` descriptor.
#' @param server Optional owner assertion. For a cross-owner table it may name
#'   either declared source owner; it never triggers column discovery.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPContingency` object containing the noisy table and
#'   privacy/accuracy metadata.
#' @export
ds.vertDPContingency <- function(data_name, row_var, col_var, server = NULL,
                                 datasources = NULL) {
  .dsvert_dp_contingency_impl(
    data_name, row_var, col_var, server, datasources,
    DSI::datashield.aggregate)
}

.dsvert_dp_contingency_impl <- function(data_name, row_var, col_var,
                                        server = NULL, datasources = NULL,
                                        .aggregate) {
  for (value in list(data_name = data_name, row_var = row_var,
                     col_var = col_var)) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop("data_name, row_var and col_var must be non-empty strings",
           call. = FALSE)
    }
  }
  if (identical(row_var, col_var)) {
    stop("row_var and col_var must identify two distinct variables",
         call. = FALSE)
  }
  datasources <- .dsvert_dp_datasources(datasources)
  owner_assertion <- .dsvert_dp_vector_server_filter(server, datasources)
  run <- .dsvert_dp_capsule_vector_run(
    datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run)
  requested <- sort(c(row_var, col_var), method = "radix")
  predicate <- function(candidate) {
    descriptor <- candidate$descriptor
    left <- tryCatch(descriptor$left$column,
                     error = function(error) NULL)
    right <- tryCatch(descriptor$right$column,
                      error = function(error) NULL)
    columns_match <- is.character(left) && length(left) == 1L && !is.na(left) &&
      is.character(right) && length(right) == 1L && !is.na(right) &&
      identical(sort(c(left, right), method = "radix"), requested)
    if (!isTRUE(columns_match)) return(FALSE)
    cross <- identical(
      descriptor$version,
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION)
    if (isTRUE(cross)) {
      datasets <- c(descriptor$left$dataset, descriptor$right$dataset)
      owners <- c(descriptor$left$owner_peer, descriptor$right$owner_peer)
      identical(length(datasets), 2L) && data_name %in% datasets &&
        (is.null(owner_assertion) || owner_assertion %in% owners)
    } else {
      identical(candidate$dataset, data_name) &&
        (is.null(owner_assertion) ||
           identical(candidate$owner_peer, owner_assertion))
    }
  }
  block <- .dsvert_dp_capsule_single_block(
    context$layout, "categorical_pairs", predicate = predicate,
    description = paste0("signed categorical-pair block for '", row_var,
                         "' and '", col_var, "'"))
  descriptor <- block$descriptor
  cross_owner <- identical(
    descriptor$version,
    .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION)
  left_column <- descriptor$left$column
  right_column <- descriptor$right$column
  left_levels <- .dsvert_dp_capsule_manifest_strings(
    descriptor$left$levels, "left categorical levels", sorted = TRUE)
  right_levels <- .dsvert_dp_capsule_manifest_strings(
    descriptor$right$levels, "right categorical levels", sorted = TRUE)
  counts <- .dsvert_dp_capsule_vector_values(context$release, block)
  expected <- as.double(length(left_levels)) * length(right_levels)
  capacity <- if (isTRUE(cross_owner)) {
    value <- .dsvert_dp_capsule_manifest_numbers(
      context$manifest$admission$unit_capacity, "unit capacity")
    if (length(value) != 1L || value < 1 ||
        value != floor(value) || value > 2^53 - 1) {
      stop("The signed capsule unit capacity is invalid", call. = FALSE)
    }
    as.numeric(value)
  } else {
    count_block <- .dsvert_dp_capsule_single_block(
      context$layout, "admitted_count",
      description = "signed admitted-count capacity block")
    .dsvert_dp_vector_block_capacity(count_block)
  }
  if (length(counts) != expected || any(counts < 0) ||
      any(counts > capacity)) {
    stop("The released contingency block violates its signed domain",
         call. = FALSE)
  }
  canonical <- matrix(
    counts, nrow = length(left_levels), ncol = length(right_levels),
    dimnames = list(left_levels, right_levels))
  if (identical(left_column, row_var) && identical(right_column, col_var)) {
    table <- canonical
  } else if (identical(left_column, col_var) &&
             identical(right_column, row_var)) {
    table <- t(canonical)
  } else {
    stop("The signed contingency descriptor changed during mapping",
         call. = FALSE)
  }
  simultaneous <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest, coordinate_count = length(table),
    confidence = 0.95, maximum_error = capacity)
  marginal <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest, coordinate_count = 1L,
    confidence = 0.95, maximum_error = capacity)
  owners <- if (isTRUE(cross_owner)) {
    c(descriptor$left$owner_peer, descriptor$right$owner_peer)
  } else {
    block$owner_peer
  }
  datasets <- if (isTRUE(cross_owner)) {
    c(descriptor$left$dataset, descriptor$right$dataset)
  } else {
    block$dataset
  }
  owners <- sort(unique(owners), method = "radix")
  datasets <- sort(unique(datasets), method = "radix")
  result <- c(.dsvert_dp_vector_public_metadata(context), list(
    server = if (isTRUE(cross_owner)) {
      paste(owners, collapse = "::")
    } else {
      block$owner_peer
    },
    servers = unname(owners), datasets = unname(datasets),
    cross_owner = isTRUE(cross_owner),
    coordinate_family = "categorical_pairs",
    coordinate_descriptor = descriptor,
    row_var = row_var, col_var = col_var,
    row_levels = unname(rownames(table)),
    col_levels = unname(colnames(table)),
    nrow = as.integer(nrow(table)), ncol = as.integer(ncol(table)),
    counts = unname(as.numeric(table)), table = table,
    coordinate_maximum = capacity,
    unit_aggregation_policy = descriptor$repeated_record_policy,
    missingness_policy = descriptor$missingness_policy,
    artifact_l1_sensitivity = if (isTRUE(cross_owner)) {
      as.numeric(descriptor$selected_l1_sensitivity)
    } else if (identical(
      context$adjacency, "add_remove_patient")) 1 else 2,
    artifact_l2_sensitivity = if (isTRUE(cross_owner)) {
      as.numeric(descriptor$selected_l2_sensitivity)
    } else if (identical(
      context$adjacency, "add_remove_patient")) 1 else sqrt(2),
    accuracy_95_abs_per_cell = rep(
      marginal$radius, length(table)),
    accuracy_simultaneous_95_abs = simultaneous$radius,
    accuracy_simultaneous_confidence = simultaneous$confidence,
    accuracy_simultaneous_method = simultaneous$method,
    accuracy_implementation_tv_upper_bound =
      simultaneous$implementation_tv_upper_bound,
    accuracy_additional_privacy_cost =
      simultaneous$additional_privacy_cost))
  result$uncertainty_scope <-
    "DP mechanism noise only; sampling uncertainty excluded"
  result$inferential_scope <- paste(
    "DP-noised cell estimates; no sampling confidence interval or ordinary",
    "chi-square/Fisher p-value is provided")
  class(result) <- c("ds.vertDPContingency", "list")
  result
}

.dsvert_dp_meanvar_simultaneous_radii <- function(result) {
  signed <- result$accuracy_simultaneous_95_abs_raw_coordinates
  if (identical(result$backend, "exact_signed_Ring128_global_vector") &&
      is.numeric(signed) && length(signed) == 3L &&
      !is.null(names(signed)) &&
      setequal(names(signed), c("count", "normalized_sum",
                                "normalized_sumsq")) &&
      !anyNA(signed) && all(is.finite(signed)) && all(signed >= 0)) {
    signed <- signed[c("count", "normalized_sum", "normalized_sumsq")]
    names(signed) <- c("count", "sum", "sumsq")
    return(signed)
  }
  grid_scale <- result$numeric_grid_scale
  if (identical(result$noise_selection$winner, "gaussian")) {
    radius <- result$noise_selection$gaussian$simultaneous_95_abs
    return(c(count = radius, sum = radius, sumsq = radius))
  }
  epsilon_coordinate <- result$epsilon / 3
  c(
    count = .dsvert_dp_google_ci_radius(
      epsilon_coordinate, 1, alpha = 0.05 / 3),
    sum = .dsvert_dp_google_ci_radius(
      epsilon_coordinate, grid_scale, alpha = 0.05 / 3),
    sumsq = .dsvert_dp_google_ci_radius(
      epsilon_coordinate, grid_scale, alpha = 0.05 / 3))
}

.dsvert_dp_meanvar_projection_interval <- function(
    value, lower, upper, radius) {
  guard <- 256 * .Machine$double.eps *
    max(1, abs(value), abs(lower), abs(upper), radius)
  if (upper - lower <= guard) {
    return(list(interval = c(lower = -Inf, upper = Inf),
                state = "collapsed_physical_projection"))
  }
  if (value < lower - guard || value > upper + guard) {
    return(list(interval = c(lower = -Inf, upper = Inf),
                state = "unverifiable_projection_geometry"))
  }
  if (value <= lower + guard) {
    return(list(
      interval = c(lower = -Inf, upper = lower + radius + guard),
      state = "lower_censored_projection"))
  }
  if (value >= upper - guard) {
    return(list(
      interval = c(lower = upper - radius - guard, upper = Inf),
      state = "upper_censored_projection"))
  }
  if (abs(value - round(value)) > guard) {
    return(list(interval = c(lower = -Inf, upper = Inf),
                state = "unverifiable_projection_geometry"))
  }
  list(
    interval = c(
      lower = value - radius - guard,
      upper = value + radius + guard),
    state = "interior_projection")
}

.dsvert_dp_meanvar_mechanism_region <- function(result, raw_radius = NULL) {
  if (is.null(raw_radius)) {
    raw_radius <- .dsvert_dp_meanvar_simultaneous_radii(result)
  }
  if (!is.numeric(raw_radius) || length(raw_radius) != 3L ||
      is.null(names(raw_radius)) ||
      !setequal(names(raw_radius), c("count", "sum", "sumsq")) ||
      anyNA(raw_radius) || any(!is.finite(raw_radius)) ||
      any(raw_radius < 0)) {
    stop("The simultaneous DP moment radii are invalid", call. = FALSE)
  }
  raw_radius <- raw_radius[c("count", "sum", "sumsq")]
  maximum_integer <- 2^53 - 1
  n_lower <- max(0, result$n - raw_radius[["count"]])
  n_upper <- min(maximum_integer, result$n + raw_radius[["count"]])
  count_region <- c(lower = ceiling(n_lower), upper = floor(n_upper))
  width <- result$upper_bound - result$lower_bound
  full_region <- function(status, projection_state) {
    has_positive <- count_region[["upper"]] >= 1
    list(
      regions = list(
        effective_count = count_region,
        mean = if (has_positive) {
          c(lower = result$lower_bound, upper = result$upper_bound)
        } else c(lower = NA_real_, upper = NA_real_),
        variance = if (has_positive) {
          c(lower = 0, upper = width^2 / 4)
        } else c(lower = NA_real_, upper = NA_real_)),
      status = status,
      includes_non_estimable = count_region[["lower"]] == 0,
      projection_state = projection_state,
      tightness = "full_physical_region_after_projection_information_loss",
      raw_simultaneous_radius = raw_radius)
  }
  if (count_region[["upper"]] < 1) {
    return(full_region(
      "dp_region_has_no_positive_effective_count",
      c(sum = "not_available", sumsq = "not_available")))
  }
  if (result$n <= 0 || is.null(result$mean) || is.null(result$variance)) {
    return(full_region(
      "dp_region_includes_zero_effective_count",
      c(sum = "point_not_available", sumsq = "point_not_available")))
  }

  mean_normalized <-
    (result$mean - result$lower_bound) / width
  variance_normalized <- result$variance / width^2
  projected_sum <- mean_normalized * result$n * result$numeric_grid_scale
  sum_projection <- .dsvert_dp_meanvar_projection_interval(
    projected_sum, 0, result$n * result$numeric_grid_scale,
    raw_radius[["sum"]])
  projected_sumsq <-
    (variance_normalized + mean_normalized^2) *
    result$n * result$numeric_grid_scale
  sumsq_projection <- .dsvert_dp_meanvar_projection_interval(
    projected_sumsq,
    projected_sum^2 /
      (result$n * result$numeric_grid_scale),
    projected_sum,
    raw_radius[["sumsq"]])

  physical_coordinate_upper <- min(
    maximum_integer,
    count_region[["upper"]] * result$numeric_grid_scale)
  sum_interval <- c(
    lower = max(0, sum_projection$interval[["lower"]]),
    upper = min(
      physical_coordinate_upper,
      sum_projection$interval[["upper"]]))
  sumsq_interval <- c(
    lower = max(0, sumsq_projection$interval[["lower"]]),
    upper = min(
      physical_coordinate_upper,
      sumsq_projection$interval[["upper"]]))
  projection_state <- c(
    sum = sum_projection$state,
    sumsq = sumsq_projection$state)
  if (any(projection_state == "unverifiable_projection_geometry")) {
    return(full_region(
      "dp_projection_preimage_geometry_unverifiable", projection_state))
  }
  if (anyNA(c(sum_interval, sumsq_interval)) ||
      any(!is.finite(c(sum_interval, sumsq_interval))) ||
      sum_interval[["lower"]] > sum_interval[["upper"]] ||
      sumsq_interval[["lower"]] > sumsq_interval[["upper"]]) {
    return(full_region(
      "dp_projection_preimage_geometry_unverifiable", projection_state))
  }

  center <- function(interval) {
    0.5 * interval[["lower"]] + 0.5 * interval[["upper"]]
  }
  half_width <- function(interval) {
    0.5 * (interval[["upper"]] - interval[["lower"]])
  }
  count_center <- center(count_region)
  count_half_width <- half_width(count_region)
  propagated <- .dsvert_dp_describe_moment_region(
    n_dp = count_center,
    qsum_dp = center(sum_interval),
    qsumsq_dp = center(sumsq_interval),
    count_radius = count_half_width,
    sum_radius = half_width(sum_interval),
    sumsq_radius = half_width(sumsq_interval),
    grid_scale = result$numeric_grid_scale,
    lower_bound = result$lower_bound,
    upper_bound = result$upper_bound)
  propagated$effective_count <- count_region
  status <- if (count_region[["lower"]] == 0) {
    "dp_region_includes_zero_effective_count"
  } else {
    propagated$status
  }
  list(
    regions = list(
      effective_count = propagated$effective_count,
      mean = propagated$mean,
      variance = propagated$variance),
    status = status,
    includes_non_estimable = count_region[["lower"]] == 0,
    projection_state = projection_state,
    tightness = if (all(projection_state == "interior_projection")) {
      "analytic_independent_coordinate_box_with_physical_moment_constraints"
    } else {
      paste(
        "conservative_projection_preimage_box; cross-coordinate integer",
        "feasibility is relaxed")
    },
    raw_simultaneous_radius = raw_radius)
}

.dsvert_dp_meanvar_postprocess <- function(result) {
  if (!isTRUE(result$released)) {
    result$mechanism_regions <- NULL
    result$mechanism_region_status <- NULL
    return(result)
  }
  region <- .dsvert_dp_meanvar_mechanism_region(result)
  result$mechanism_regions <- region$regions
  result$mechanism_region_status <- region$status
  result$mechanism_region_includes_non_estimable <-
    region$includes_non_estimable
  result$mechanism_region_projection_state <- region$projection_state
  result$mechanism_region_tightness <- region$tightness
  result$accuracy_simultaneous_95_abs_raw_coordinates <-
    region$raw_simultaneous_radius
  result$mechanism_region_confidence <- 0.95
  result$mechanism_region_method <- paste(
    "simultaneous count/sum/sumsq noise box; inverse physical projection;",
    "bounded-moment and per-unit quantisation constraints")
  result$mechanism_region_scope <- paste(
    "DP mechanism noise plus deterministic quantisation enclosure;",
    "sampling uncertainty excluded")
  result$mechanism_region_additional_privacy_cost <-
    c(epsilon = 0, delta = 0)
  result$mechanism_region_additional_server_calls <- 0L
  result
}

#' Differentially private bounded mean and variance
#'
#' The signed capsule clips one contribution per protected unit to immutable
#' bounds and releases count, normalized quantized sum and normalized
#' quantized sum of squares as part of one global Ring128 vector. There is no
#' coordinate-wise epsilon split. The client applies only bounded-moment
#' post-processing and the exact affine conversion back to the original
#' scale. The reported variance is the population central second moment with
#' the DP-noisy denominator eqn{n}, not an eqn{n-1} sampling estimator.
#' Simultaneous regions cover mechanism noise and deterministic quantization;
#' they are not sampling confidence intervals.
#'
#' @param data_name Name of the protected data frame.
#' @param variable Numeric variable with custodian-configured bounds.
#' @param server Optional server name; auto-detected when unambiguous.
#' @param datasources DataSHIELD connections.
#' @return A DP moment release with bounds and accuracy metadata.
#' @export
ds.vertDPMeanVar <- function(data_name, variable, server = NULL,
                             datasources = NULL) {
  .dsvert_dp_meanvar_impl(
    data_name, variable, server, datasources, DSI::datashield.aggregate)
}

.dsvert_dp_meanvar_impl <- function(data_name, variable, server = NULL,
                                    datasources = NULL, .aggregate) {
  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) || !nzchar(data_name) ||
      !is.character(variable) || length(variable) != 1L ||
      is.na(variable) || !nzchar(variable)) {
    stop("data_name and variable must be non-empty strings", call. = FALSE)
  }
  datasources <- .dsvert_dp_datasources(datasources)
  owner <- .dsvert_dp_vector_server_filter(server, datasources)
  run <- .dsvert_dp_capsule_vector_run(
    datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run)
  block <- .dsvert_dp_capsule_single_block(
    context$layout, "numeric_moments", dataset = data_name,
    owner_peer = owner,
    predicate = function(candidate) {
      identical(candidate$descriptor$column, variable)
    },
    description = paste0("signed numeric-moment block for '", variable,
                         "'"))
  descriptor <- block$descriptor
  bounds <- .dsvert_dp_capsule_manifest_numbers(
    list(descriptor$lower, descriptor$upper), "numeric bounds")
  grid_bits <- descriptor$numeric_grid_bits
  if (length(bounds) != 2L || bounds[[1L]] >= bounds[[2L]] ||
      !.dsvert_dp_is_integer(grid_bits, 8, 18) ||
      !identical(as.numeric(grid_bits),
                 as.numeric(context$lattice$output_lattice_bits))) {
    stop("The signed numeric-moment descriptor is invalid", call. = FALSE)
  }
  grid_bits <- as.integer(grid_bits)
  grid_scale <- 2^grid_bits
  capacity <- .dsvert_dp_vector_block_capacity(
    block, numeric_moment = TRUE)
  coordinates <- .dsvert_dp_capsule_vector_values(context$release, block)
  if (length(coordinates) != 3L || any(coordinates < 0) ||
      any(coordinates > capacity)) {
    stop("The released numeric moments violate their signed bounds",
         call. = FALSE)
  }
  n_dp <- unname(coordinates[[1L]])
  normalized_sum_dp <- unname(coordinates[[2L]])
  normalized_sumsq_dp <- unname(coordinates[[3L]])
  mean_value <- variance_value <- sum_value <- sumsq_value <- NULL
  projected_sum <- projected_sumsq <- NULL
  if (n_dp > 0) {
    projected_sum <- min(n_dp, max(0, normalized_sum_dp))
    projected_sumsq <- min(
      projected_sum,
      max(projected_sum^2 / n_dp, normalized_sumsq_dp))
    normalized_mean <- projected_sum / n_dp
    normalized_variance <- max(
      0, projected_sumsq / n_dp - normalized_mean^2)
    width <- bounds[[2L]] - bounds[[1L]]
    mean_value <- bounds[[1L]] + width * normalized_mean
    variance_value <- width^2 * normalized_variance
    sum_value <- n_dp * bounds[[1L]] + width * projected_sum
    sumsq_value <- n_dp * bounds[[1L]]^2 +
      2 * bounds[[1L]] * width * projected_sum +
      width^2 * projected_sumsq
  }

  marginal <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest, coordinate_count = 1L,
    confidence = 0.95, maximum_error = capacity)
  simultaneous <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest, coordinate_count = 3L,
    confidence = 0.95, maximum_error = capacity)
  profile <- .dsvert_vector_profile(
    context$manifest$workload$capsule_mechanism,
    context$manifest$workload$mechanism_selection)
  region <- .dsvert_dp_describe_moment_region(
    n_dp = n_dp,
    qsum_dp = normalized_sum_dp * grid_scale,
    qsumsq_dp = normalized_sumsq_dp * grid_scale,
    count_radius = simultaneous$radius,
    sum_radius = simultaneous$radius * grid_scale,
    sumsq_radius = simultaneous$radius * grid_scale,
    grid_scale = grid_scale,
    lower_bound = bounds[[1L]], upper_bound = bounds[[2L]])
  count_certified <- region$effective_count[["lower"]] > 0
  point_available <- n_dp > 0
  release_status <- if (!point_available) {
    "dp_effective_count_not_certified_positive"
  } else if (!count_certified) {
    "dp_point_available_count_not_certified_positive"
  } else {
    "ok"
  }
  reason <- if (identical(release_status, "ok")) NULL else
    "dp_noisy_effective_count_lower_bound_is_zero"
  result <- c(.dsvert_dp_vector_public_metadata(context), list(
    server = block$owner_peer,
    coordinate_family = "numeric_moments",
    coordinate_descriptor = descriptor,
    variable = variable,
    status = release_status, reason = reason,
    n = n_dp,
    n_definition = paste(
      "nonnegative fixed-lattice DP noisy effective-unit count;",
      "not rounded by the client"),
    effective_count_95_lower_bound =
      unname(region$effective_count[["lower"]]),
    mean = mean_value, variance = variance_value,
    sum = sum_value, sumsq = sumsq_value,
    normalized_sum_dp = normalized_sum_dp,
    normalized_sumsq_dp = normalized_sumsq_dp,
    normalized_sum_projected = projected_sum,
    normalized_sumsq_projected = projected_sumsq,
    lower_bound = bounds[[1L]], upper_bound = bounds[[2L]],
    variance_definition =
      "population_central_second_moment_denominator_dp_noisy_n",
    numeric_grid_bits = grid_bits, numeric_grid_scale = grid_scale,
    max_abs_quantization_per_unit = 0.5 / grid_scale,
    moment_postprocessing = paste(
      "bounded consistent normalized moments, then exact affine conversion",
      "to the signed original scale"),
    artifact_l1_sensitivity = 3,
    artifact_l2_sensitivity = sqrt(3),
    coordinate_natural_l1_sensitivity = c(
      count = 1, normalized_sum = 1, normalized_sumsq = 1),
    epsilon_allocation = c(global_capsule_vector = context$release$epsilon),
    submechanism_count = 1L,
    noise_selection = list(
      winner = if (isTRUE(profile$gaussian)) {
        "two_peer_full_vector_dyadic_discrete_gaussian_tv_bounded"
      } else {
        "two_peer_full_vector_convolution"
      },
      mechanism = context$release$mechanism,
      objective = if (isTRUE(profile$gaussian)) {
        "signed_global_capsule_l2_sensitivity"
      } else {
        "signed_global_capsule_l1_sensitivity"
      },
      sensitivity_norm = if (isTRUE(profile$gaussian)) "l2" else "l1",
      coordinate_epsilon_split = FALSE),
    accuracy_95_abs_count = marginal$radius,
    accuracy_95_abs_normalized_sum_noise_only = marginal$radius,
    accuracy_95_abs_normalized_sumsq_noise_only = marginal$radius,
    accuracy_simultaneous_95_abs_raw_coordinates = c(
      count = simultaneous$radius,
      normalized_sum = simultaneous$radius,
      normalized_sumsq = simultaneous$radius),
    mechanism_regions = list(
      effective_count = region$effective_count,
      mean = region$mean, variance = region$variance),
    mechanism_region_status = region$status,
    mechanism_region_includes_non_estimable =
      region$effective_count[["lower"]] == 0,
    mechanism_region_confidence = simultaneous$confidence,
    mechanism_region_method = paste(
      simultaneous$method,
      "with bounded-moment and public quantisation propagation"),
    mechanism_region_scope = paste(
      "DP mechanism noise plus deterministic quantisation enclosure;",
      "sampling uncertainty excluded"),
    mechanism_region_additional_privacy_cost =
      simultaneous$additional_privacy_cost,
    mechanism_region_additional_server_calls = 0L,
    accuracy_implementation_tv_upper_bound =
      simultaneous$implementation_tv_upper_bound))
  result$uncertainty_scope <- paste(
    "Simultaneous regions cover DP mechanism noise only, with a deterministic",
    "quantisation enclosure; clipping changes the finite-snapshot estimand",
    "and sampling uncertainty is excluded")
  result$inferential_scope <- paste(
    "DP-noised bounded point estimates with conservative mechanism regions;",
    "no sampling confidence interval or hypothesis test is provided")
  class(result) <- c("ds.vertDPMeanVar", "list")
  result
}

#' @export
print.ds.vertDPCount <- function(x, ...) {
  if (!isTRUE(x$released)) {
    cat("dsVert DP count: suppressed (", x$reason, ")\n", sep = "")
  } else {
    cat("dsVert DP count:", x$value, "[", x$server, "]\n")
    cat("epsilon:", format(x$epsilon), " | 95% noise radius:",
        x$accuracy_95_abs, "\n")
    if (isTRUE(x$clipped_coordinates > 0)) {
      cat("WARNING: the public release clamp affected this coordinate.\n")
    }
    cat(x$uncertainty_scope, "\n")
  }
  invisible(x)
}

#' @export
print.ds.vertDPContingency <- function(x, ...) {
  if (!isTRUE(x$released)) {
    cat("dsVert DP contingency: suppressed (", x$reason, ")\n", sep = "")
  } else {
    print(x$table, ...)
    cat("epsilon:", format(x$epsilon),
        " | simultaneous 95% max-cell noise radius:",
        x$accuracy_simultaneous_95_abs, "(",
        x$accuracy_simultaneous_method, ")\n")
    if (isTRUE(x$clipped_coordinates > 0)) {
      cat("WARNING:", x$clipped_coordinates,
          "public release coordinates were clamped.\n")
    }
    cat(x$uncertainty_scope, "\n")
  }
  invisible(x)
}

#' @export
print.ds.vertDPMeanVar <- function(x, ...) {
  if (!isTRUE(x$released)) {
    cat("dsVert DP mean/variance: suppressed (", x$reason, ")\n", sep = "")
  } else if (identical(
      x$status, "dp_point_available_count_not_certified_positive")) {
    cat("dsVert DP mean:", format(x$mean), " | variance:",
        format(x$variance), "[", x$server, "]\n")
    cat("WARNING: point estimate available, but the 95% lower bound for ",
        "the DP noisy effective count is zero (n_dp=", x$n, ").\n",
        sep = "")
  } else if (!identical(x$status, "ok")) {
    cat("dsVert DP mean/variance:", x$status, "[", x$server, "]\n")
    cat("DP noisy effective n:", x$n, " | 95% lower bound:",
        x$effective_count_95_lower_bound, "\n")
  } else {
    cat("dsVert DP mean:", format(x$mean), " | variance:",
        format(x$variance), "[", x$server, "]\n")
    cat("bounds: [", x$lower_bound, ", ", x$upper_bound,
        "] | epsilon: ", format(x$epsilon), "\n", sep = "")
    if (isTRUE(x$clipped_coordinates > 0)) {
      cat("WARNING:", x$clipped_coordinates,
          "public release coordinates were clamped.\n")
    }
  }
  if (isTRUE(x$released)) {
    cat("simultaneous 95% mechanism-region status:",
        x$mechanism_region_status, "\n")
    if (!is.null(x$mechanism_regions)) {
      print(do.call(rbind, x$mechanism_regions), ...)
    }
    cat(x$uncertainty_scope, "\n")
  }
  invisible(x)
}
