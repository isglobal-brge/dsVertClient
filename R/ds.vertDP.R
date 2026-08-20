.DSVERT_DP_MINIMUM_EPSILON <- 2^-50
.DSVERT_DP_MAXIMUM_EPSILON <- 2^40
.DSVERT_DP_DEFAULT_CAPSULE_DELTA <- 2^-100
.DSVERT_DP_NOISE_SELECTOR <- "minimum_conservative_95_radius_v3"
.DSVERT_DP_GAUSSIAN_TV_BOUND_PER_COORDINATE <- 2^-40

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

.dsvert_dp_vector_context <- function(run, allow_synopsis = FALSE) {
  synopsis <- is.list(run) && is.list(run$release) &&
    inherits(run$release, "dsvert_synopsis_public_vector") &&
    identical(
      run$release$signed_provenance$version,
      "dsvert-stateless-synopsis-public-provenance-v1")
  if (isTRUE(synopsis) && !isTRUE(allow_synopsis)) {
    stop("This biomedical consumer does not yet accept synopsis provenance",
         call. = FALSE)
  }
  if (!is.list(run) || !inherits(run$release, "dsvert_joint_dp_vector") ||
      !is.list(run$layout) || !is.list(run$status) || !length(run$status) ||
      is.null(names(run$status)) || anyNA(names(run$status)) ||
      anyDuplicated(names(run$status)) ||
      !identical(as.numeric(run$layout$coordinate_count),
                 as.numeric(run$release$coordinate_count)) ||
      !is.list(run$release$manifest) ||
      (!isTRUE(synopsis) &&
       (!identical(run$release$history_gate, TRUE) ||
        !identical(run$release$request_limit, FALSE) ||
        !identical(run$release$operation_limit, TRUE)))) {
    stop("The joint biomedical DP vector context is invalid", call. = FALSE)
  }
  if (isTRUE(synopsis) &&
      any(c("capsule_id", "privacy_epoch", "noise_key_id",
            "history_gate", "request_limit", "operation_limit") %in%
          names(run$release))) {
    stop("The synopsis vector contains incompatible legacy provenance",
         call. = FALSE)
  }
  if (isTRUE(synopsis)) {
    provenance <- run$release$signed_provenance
    bindings <- c(
      "artifact_key", "execution_id", "contract_sha256",
      "attempt_sha256", "source_contract_sha256", "result_set_sha256",
      "final_vector_root")
    if (!all(vapply(bindings, function(field) {
      .dsvert_vector_hex(run$release[[field]]) &&
        identical(run$release[[field]], provenance[[field]])
    }, logical(1L)))) {
      stop("The synopsis vector provenance is detached", call. = FALSE)
    }
  }
  status_adjacency <- tryCatch(vapply(run$status, function(value) {
    value$policy$adjacency
  }, character(1L)), error = function(error) character())
  adjacency <- if (isTRUE(synopsis)) {
    run$release$manifest$admission$adjacency
  } else if (length(status_adjacency)) status_adjacency[[1L]] else NULL
  if (length(status_adjacency) != length(run$status) ||
      !all(status_adjacency %in% c(
        "add_remove_patient", "replace_one_fixed_cohort")) ||
      length(unique(status_adjacency)) != 1L ||
      !is.character(adjacency) || length(adjacency) != 1L ||
      is.na(adjacency) ||
      !adjacency %in% c("add_remove_patient", "replace_one_fixed_cohort") ||
      (isTRUE(synopsis) &&
       !identical(adjacency, status_adjacency[[1L]]))) {
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
       verification_compilation = run$verification_compilation,
       manifest = run$release$manifest, lattice = lattice,
       adjacency = adjacency, synopsis = synopsis)
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

.dsvert_dp_synopsis_vector_public_metadata <- function(context) {
  if (!is.list(context) || !isTRUE(context$synopsis)) {
    stop("The Synopsis biomedical DP vector context is invalid",
         call. = FALSE)
  }
  release <- context$release
  if (!identical(release$sticky_replay, TRUE) ||
      !identical(release$signed_provenance$durable_replay, TRUE)) {
    stop("The Synopsis biomedical DP vector has no durable replay binding",
         call. = FALSE)
  }
  lattice <- context$lattice
  profile <- .dsvert_dp_vector_release_profile(release, context$manifest)
  list(
    released = TRUE,
    mechanism = release$mechanism,
    implementation = profile$backend,
    backend = "exact_signed_Ring128_canonical_synopsis_vector",
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
    sensitivity_scope = "complete_signed_canonical_synopsis_vector",
    synopsis_mechanism = context$manifest$workload$capsule_mechanism,
    mechanism_selection = context$manifest$workload$mechanism_selection,
    output_lattice_bits = as.integer(lattice$output_lattice_bits),
    output_lattice_scale = as.numeric(lattice$output_lattice_scale),
    artifact_key = release$artifact_key,
    execution_id = release$execution_id,
    manifest_sha256 = release$manifest_sha256,
    contract_sha256 = release$contract_sha256,
    attempt_sha256 = release$attempt_sha256,
    source_contract_sha256 = release$source_contract_sha256,
    result_set_sha256 = release$result_set_sha256,
    final_vector_root = release$final_vector_root,
    coordinate_order_sha256 = release$coordinate_order_sha256,
    synopsis_coordinate_count = as.integer(release$coordinate_count),
    sticky_noise = "one_immutable_canonical_synopsis_artifact",
    sticky_replay = isTRUE(release$sticky_replay),
    unlimited_replay = TRUE,
    release_provenance = release$signed_provenance,
    privacy = list(
      version = "dsvert-per-synopsis-dp-v1",
      per_artifact_epsilon = release$epsilon,
      per_artifact_delta = release$delta,
      sticky_noise = TRUE,
      unlimited_replay = TRUE,
      replay_is_postprocessing = TRUE,
      public_openings = 1L,
      distinct_artifacts_compose = TRUE,
      finite_global_composition_claim = FALSE),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    composition_rule = "one_sticky_release_per_canonical_signed_artifact",
    data_dependency = "immutable_signed_canonical_synopsis_snapshot",
    security_claim = list(
      version = "dsvert-synopsis-security-claim-v1",
      privacy_definition = "per_synopsis_epsilon_delta_dp",
      two_pinned_noise_authorities = TRUE,
      maximum_colluding_noise_authorities = 1L,
      analyst_relay_trusted = FALSE,
      replay_is_postprocessing = TRUE,
      allocation_openings_used = FALSE,
      finite_global_composition_claim = FALSE),
    clipped_coordinates = NA_integer_,
    clamp_activation_disclosed = FALSE)
}

.dsvert_dp_vector_public_metadata <- function(context) {
  if (isTRUE(context$synopsis)) {
    return(.dsvert_dp_synopsis_vector_public_metadata(context))
  }
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

.dsvert_dp_count_synopsis_result_v1 <- function(
    data_name, server, datasources, .aggregate,
    .run = .dsvert_dp_synopsis_vector_run) {
  run <- .run(datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count",
    description = "signed admitted-count block")
  descriptor <- block$descriptor
  if (!identical(block$dataset, data_name) ||
      !identical(descriptor$dataset, block$dataset) ||
      !identical(descriptor$owner_peer, block$owner_peer) ||
      !identical(as.numeric(descriptor$statistic_minimum), 0) ||
      !block$owner_peer %in% names(datasources)) {
    stop("The signed admitted-count block does not match data_name",
         call. = FALSE)
  }
  capacity <- .dsvert_dp_vector_block_capacity(block)
  value <- .dsvert_dp_capsule_vector_values(context$release, block)
  if (length(value) != 1L || !is.finite(value) || value < 0 ||
      value > capacity) {
    stop("The released admitted-count block violates its signed domain",
         call. = FALSE)
  }
  metadata <- .dsvert_dp_vector_public_metadata(context)
  provenance <- metadata$release_provenance
  peers <- tryCatch(names(provenance$ordered_peer_pinset),
                    error = function(error) NULL)
  pinset <- tryCatch(.dsvert_dp_synopsis_client_pinset_v1(
    provenance$ordered_peer_pinset, peers), error = function(error) NULL)
  designated <- tryCatch(.dsvert_dp_synopsis_client_string_list_v1(
    provenance$designated_noise_peers, "Count noise authority set", 2L),
    error = function(error) NULL)
  if (is.null(pinset) || length(pinset) < 2L || is.null(designated) ||
      length(designated) != 2L || anyDuplicated(designated) ||
      any(!designated %in% names(pinset)) ||
      !identical(provenance$peer_pinset_sha256,
                 .dsvert_vector_hash(as.list(pinset)))) {
    stop("The signed Count noise authorities are invalid", call. = FALSE)
  }
  finalizer <- designated[[2L]]
  if (!is.null(server) && !identical(server, finalizer)) {
    stop("server must match the signed Count finalizer", call. = FALSE)
  }
  accuracy <- .dsvert_dp_vector_accuracy_radius(
    context$release, context$manifest, coordinate_count = 1L,
    confidence = 0.95, maximum_error = capacity)
  result <- c(metadata, list(
    value = unname(value), server = finalizer,
    source_owner = block$owner_peer, dataset = block$dataset,
    coordinate_family = "admitted_count", coordinate_descriptor = descriptor,
    coordinate_maximum = capacity,
    primitive = "signed_admitted_count_v1",
    artifact_l1_sensitivity = if (identical(
      context$adjacency, "add_remove_patient")) 1 else 2,
    artifact_l2_sensitivity = if (identical(
      context$adjacency, "add_remove_patient")) 1 else sqrt(2),
    accuracy_95_abs = accuracy$radius,
    accuracy_95_confidence = accuracy$confidence,
    accuracy_95_method = accuracy$method,
    accuracy_implementation_tv_upper_bound =
      accuracy$implementation_tv_upper_bound,
    accuracy_additional_privacy_cost = accuracy$additional_privacy_cost,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    operation_id = metadata$execution_id,
    release_sha256 = metadata$final_vector_root,
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    inferential_scope = paste(
      "Finite-snapshot DP count; no sampling confidence interval or",
      "population-supermodel inference is provided")))
  class(result) <- c("ds.vertDPCount", "list")
  result
}

#' Differentially private privacy-unit count
#'
#' The signed durable Synopsis contains one bounded \code{admitted_count}
#' coordinate for its canonical aligned snapshot.
#' Every connected peer validates that same artifact, while exactly two pinned
#' authorities supply the sticky mechanism randomness. This function reads and
#' validates that coordinate only: it neither draws additional noise nor opens
#' a second MPC output. Replaying the artifact is post-processing. Reported
#' accuracy covers mechanism noise, not population sampling uncertainty.
#'
#' @param data_name Name of the registered protected data frame.
#' @param server Optional connected-peer assertion. When supplied it must match
#'   the signed Synopsis finalizer.
#' @param datasources Named DataSHIELD connections.
#' @return A signed Count release with mechanism, accuracy, and per-artifact
#'   privacy metadata.
#' @export
ds.vertDPCount <- function(data_name, server = NULL, datasources = NULL) {
  resolved <- .dsvert_federation_argument(data_name, datasources)
  .dsvert_dp_count_impl(
    resolved$value, server, resolved$datasources, DSI::datashield.aggregate)
}

.dsvert_dp_count_legacy_impl <- function(data_name, server = NULL,
                                         datasources = NULL, .aggregate,
                                         .execute = .dsvert_dp_count_execute_v1) {
  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) || !nzchar(data_name)) {
    stop("data_name must be one non-empty string", call. = FALSE)
  }
  datasources <- .dsvert_dp_datasources(datasources)
  selected_server <- .dsvert_dp_server(server, datasources)
  execution <- .execute(
    data_name, datasources, .aggregate = .aggregate)
  if (!is.list(execution) || is.null(names(execution)) ||
      anyNA(names(execution)) || anyDuplicated(names(execution)) ||
      !setequal(names(execution), c("version", "mode", "payload")) ||
      !identical(execution$version,
                 .DSVERT_CLIENT_DP_COUNT_EXECUTION_RESULT_VERSION) ||
      !is.character(execution$mode) || length(execution$mode) != 1L ||
      is.na(execution$mode) ||
      !execution$mode %in% c("add_remove_dp", "fixed_cohort_public") ||
      !is.list(execution$payload)) {
    stop("Invalid closed Count execution result", call. = FALSE)
  }
  payload <- execution$payload
  if (identical(execution$mode, "add_remove_dp")) {
    payload_fields <- c(
      "release", "finalizer_peer", "accuracy_95_abs",
      "accuracy_95_confidence", "accuracy_95_method")
    release <- payload$release
    release_fields <- c(
      "version", "artifact_key", "contract_sha256",
      "analysis_binding_sha256", "worker_static_sha256", "circuit",
      "mechanism", "bounds", "value", "source_identity_pk",
      "finalizer_identity_pk", "backend", "postprocessing",
      "intermediate_values_exposed", "public_openings", "release_sha256",
      "signature")
    mechanism_fields <- c(
      "family", "version", "sampler", "epsilon", "delta",
      "implementation_delta", "sensitivity_l1")
    valid <- !is.null(names(payload)) && !anyNA(names(payload)) &&
      !anyDuplicated(names(payload)) &&
      setequal(names(payload), payload_fields) &&
      is.list(release) && !is.null(names(release)) &&
      !anyNA(names(release)) && !anyDuplicated(names(release)) &&
      setequal(names(release), release_fields) &&
      is.list(release$mechanism) &&
      setequal(names(release$mechanism), mechanism_fields) &&
      is.list(release$bounds) &&
      setequal(names(release$bounds), c("lower", "upper")) &&
      payload$finalizer_peer %in% names(datasources) &&
      .dsvert_dp_is_integer(payload$accuracy_95_abs, 0, 4096) &&
      identical(as.numeric(payload$accuracy_95_confidence), 0.95) &&
      identical(payload$accuracy_95_method,
                .DSVERT_CLIENT_DP_COUNT_ACCURACY_95_METHOD) &&
      identical(release$intermediate_values_exposed, FALSE) &&
      identical(as.numeric(release$public_openings), 1)
    lower <- suppressWarnings(as.numeric(release$bounds$lower))
    upper <- suppressWarnings(as.numeric(release$bounds$upper))
    value <- suppressWarnings(as.numeric(release$value))
    valid <- isTRUE(valid) && .dsvert_dp_is_integer(lower, 0, 0) &&
      .dsvert_dp_is_integer(upper, 1, 1000000) &&
      .dsvert_dp_is_integer(value, lower, upper)
    if (!isTRUE(valid)) {
      stop("Invalid add/remove Count execution payload", call. = FALSE)
    }
    result <- list(
      released = TRUE,
      value = value,
      server = payload$finalizer_peer,
      mechanism = release$mechanism$version,
      implementation = release$backend,
      sampler = release$mechanism$sampler,
      randomness = "two_persistent_identity_seeds_joint_exact_gc_v1",
      sensitivity = as.numeric(release$mechanism$sensitivity_l1),
      artifact_l1_sensitivity = 1,
      postprocessing = release$postprocessing,
      clipped_coordinates = NA_integer_,
      clamp_activation_disclosed = FALSE,
      accuracy_95_abs = min(as.numeric(payload$accuracy_95_abs), upper),
      accuracy_95_confidence = 0.95,
      accuracy_95_method = payload$accuracy_95_method,
      accuracy_additional_privacy_cost = 0,
      epsilon = as.numeric(release$mechanism$epsilon),
      delta = as.numeric(release$mechanism$delta),
      implementation_delta =
        as.numeric(release$mechanism$implementation_delta),
      adjacency = "add_remove_patient",
      lower_bound = lower,
      upper_bound = upper,
      artifact_key = release$artifact_key,
      release_sha256 = release$release_sha256,
      signed_release = release,
      sticky_replay = TRUE,
      intermediate_values_exposed = FALSE,
      privacy = list(
        per_artifact_epsilon = as.numeric(release$mechanism$epsilon),
        per_artifact_delta = as.numeric(release$mechanism$delta),
        sticky_noise = TRUE,
        finite_global_composition_claim = FALSE,
        distinct_artifacts_compose = TRUE,
        public_openings = 1L),
      composition_rule =
        "one_sticky_release_per_canonical_signed_artifact",
      data_dependency = "current_signed_aligned_snapshot")
  } else {
    if (is.null(names(payload)) || anyNA(names(payload)) ||
        anyDuplicated(names(payload)) || !setequal(names(payload), c(
          "declaration", "receipt_set_sha256", "peer_count"))) {
      stop("Invalid fixed-cohort Count execution payload", call. = FALSE)
    }
    declaration <- .dsvert_dp_count_client_fixed_declaration_v1(
      payload$declaration, names(datasources))
    receipt_set_sha256 <- .dsvert_dp_count_client_hex_v1(
      payload$receipt_set_sha256, "fixed receipt-set hash")
    if (!identical(as.numeric(payload$peer_count),
                   as.numeric(length(datasources)))) {
      stop("Invalid fixed-cohort Count execution payload", call. = FALSE)
    }
    result <- list(
      released = TRUE,
      value = as.numeric(declaration$fixed_cohort_size),
      server = selected_server,
      mechanism = "public_fixed_cohort_size_v1",
      implementation = "custodian_owned_signed_K_consensus",
      sampler = "none",
      randomness = "none",
      sensitivity = 0,
      postprocessing = "none_public_signed_declaration",
      clipped_coordinates = 0L,
      accuracy_95_abs = 0,
      accuracy_95_confidence = 0.95,
      accuracy_95_method = "exact_public_fixed_cohort_value",
      accuracy_additional_privacy_cost = 0,
      epsilon = 0,
      delta = 0,
      implementation_delta = 0,
      adjacency = "replace_one_fixed_cohort",
      peer_count = as.integer(payload$peer_count),
      declaration = declaration,
      receipt_set_sha256 = receipt_set_sha256,
      privacy = list(
        per_artifact_epsilon = 0,
        per_artifact_delta = 0,
        sticky_noise = FALSE,
        finite_global_composition_claim = FALSE,
        distinct_artifacts_compose = TRUE,
        public_openings = 1L),
      composition_rule = "public_signed_constant_no_dp_composition",
      data_dependency =
        "public_fixed_contract_validated_against_current_aligned_snapshot")
  }
  result$uncertainty_scope <- if (identical(
      execution$mode, "fixed_cohort_public")) {
    paste(
      "Public fixed-cohort signed value; no DP mechanism noise;",
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

.dsvert_dp_count_impl <- function(
    data_name, server = NULL, datasources = NULL, .aggregate,
    .execute = NULL, .run = .dsvert_dp_synopsis_vector_run) {
  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) || !nzchar(data_name)) {
    stop("data_name must be one non-empty string", call. = FALSE)
  }
  datasources <- .dsvert_dp_datasources(datasources)
  selected_server <- .dsvert_dp_vector_server_filter(server, datasources)
  if (is.null(.execute)) {
    return(.dsvert_dp_count_synopsis_result_v1(
      data_name, selected_server, datasources, .aggregate, .run = .run))
  }
  if (!is.function(.execute)) {
    stop("Invalid Count client dependency", call. = FALSE)
  }
  .dsvert_dp_count_legacy_impl(
    data_name, selected_server, datasources, .aggregate, .execute = .execute)
}

#' Differentially private fixed-domain contingency table
#'
#' The table domain, orientation and one-contribution-per-unit rule come from
#' signed catalog metadata. The canonical selector binds a reduced sticky
#' Synopsis artifact containing only the authorized same-owner or cross-owner
#' pair and its two marginals. Cross-owner products remain private exact-MPC
#' shares until the single noisy Synopsis release. This function only reshapes
#' that already-DP block; it does not apply an
#' ordinary chi-square or Fisher p-value, because those reference
#' distributions are not calibrated for DP-noised counts. Under the current
#' repeated-record contract, a unit contributes once only when all its valid
#' fixed-domain rows occupy the same cell. Conflicting valid cells contribute
#' zero, while missing and out-of-domain rows are ignored for consistency.
#' This public, custodian-owned rule is returned with every release.
#'
#' @param data_name Name of the protected data frame or a version-2
#'   `ds.vertFederation` carrying its public schema.
#' @param row_var,col_var Fixed-domain categorical variables authorized as one
#'   unambiguous signed same-owner or cross-owner pair. With a federation,
#'   ambiguous bare names require an exact `server$column` reference.
#' @param server Optional same-owner assertion; it never triggers column
#'   discovery.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPContingency` object containing the noisy table and
#'   privacy/accuracy metadata.
#' @export
ds.vertDPContingency <- function(data_name, row_var, col_var, server = NULL,
                                 datasources = NULL) {
  display_variables <- c(row_var, col_var)
  if (inherits(data_name, "ds.vertFederation") &&
      identical(data_name$version, 2L)) {
    references <- .dsvert_dp_synopsis_pair_federation_references_v1(
      data_name, display_variables)
    row_var <- references[[1L]]$reference
    col_var <- references[[2L]]$reference
  }
  resolved <- .dsvert_federation_argument(data_name, datasources)
  result <- .dsvert_dp_contingency_impl(
    resolved$value, row_var, col_var, server, resolved$datasources,
    DSI::datashield.aggregate)
  result$row_var <- display_variables[[1L]]
  result$col_var <- display_variables[[2L]]
  result
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
  requested_references <- lapply(
    list(row_var, col_var), .dsvert_dp_synopsis_pair_reference_v1,
    what = "A contingency variable")
  if (identical(
        requested_references[[1L]]$reference,
        requested_references[[2L]]$reference)) {
    stop("row_var and col_var must identify two distinct variables",
         call. = FALSE)
  }
  datasources <- .dsvert_dp_datasources(datasources)
  owner_assertion <- .dsvert_dp_vector_server_filter(server, datasources)
  selector_request <- list(
    version = .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION,
    family = "categorical_pair", dataset = data_name,
    columns = unname(sort(vapply(
      requested_references, `[[`, character(1L), "reference"),
      method = "radix")),
    owner_peer = owner_assertion)
  run <- .dsvert_dp_synopsis_vector_run(
    datasources, local_projection = selector_request,
    .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  accepts <- function(reference, side) {
    identical(reference$column, side$column) &&
      (is.null(reference$owner_peer) ||
         identical(reference$owner_peer, side$owner_peer))
  }
  predicate <- function(candidate) {
    descriptor <- candidate$descriptor
    sides <- tryCatch(list(
      descriptor$left[c("dataset", "column", "owner_peer")],
      descriptor$right[c("dataset", "column", "owner_peer")]),
      error = function(error) list())
    bindings <- if (length(sides) == 2L) c(
      accepts(requested_references[[1L]], sides[[1L]]) &&
        accepts(requested_references[[2L]], sides[[2L]]),
      accepts(requested_references[[1L]], sides[[2L]]) &&
        accepts(requested_references[[2L]], sides[[1L]])) else FALSE
    if (sum(bindings) != 1L) return(FALSE)
    cross <- identical(
      descriptor$version,
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION)
    if (isTRUE(cross)) {
      datasets <- c(descriptor$left$dataset, descriptor$right$dataset)
      is.null(owner_assertion) && data_name %in% datasets
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
  count_block <- .dsvert_dp_capsule_single_block(
    context$layout, "admitted_count",
    description = "signed admitted-count capacity block")
  capacity <- .dsvert_dp_vector_block_capacity(count_block)
  if (length(counts) != expected || any(counts < 0) ||
      any(counts > capacity)) {
    stop("The released contingency block violates its signed domain",
         call. = FALSE)
  }
  canonical <- matrix(
    counts, nrow = length(left_levels), ncol = length(right_levels),
    dimnames = list(left_levels, right_levels))
  left_side <- descriptor$left[c("dataset", "column", "owner_peer")]
  right_side <- descriptor$right[c("dataset", "column", "owner_peer")]
  if (accepts(requested_references[[1L]], left_side) &&
      accepts(requested_references[[2L]], right_side)) {
    table <- canonical
  } else if (accepts(requested_references[[1L]], right_side) &&
             accepts(requested_references[[2L]], left_side)) {
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
  } else block$owner_peer
  datasets <- if (isTRUE(cross_owner)) {
    c(descriptor$left$dataset, descriptor$right$dataset)
  } else block$dataset
  owners <- sort(unique(owners), method = "radix")
  datasets <- sort(unique(datasets), method = "radix")
  result <- c(.dsvert_dp_vector_public_metadata(context), list(
    server = if (isTRUE(cross_owner)) {
      paste(owners, collapse = "+")
    } else block$owner_peer,
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
    artifact_l1_sensitivity = if (identical(
      context$adjacency, "add_remove_patient")) 1 else 2,
    artifact_l2_sensitivity = if (identical(
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

#' Differentially private bounded mean and variance
#'
#' One canonical signed Synopsis artifact clips one contribution per protected
#' unit to immutable bounds and releases count, normalized quantized sum and
#' normalized quantized sum of squares in a common Ring128 vector. There is no
#' coordinate-wise epsilon split. Replay of the same durable artifact is
#' unlimited post-processing and is never denied by lifetime history. Distinct
#' artifacts compose and no finite global composition bound is claimed. The
#' client applies only bounded-moment post-processing and the exact affine
#' conversion back to the original scale. The reported variance is the
#' population central second moment with the DP-noisy denominator eqn{n}, not
#' an eqn{n-1} sampling estimator. Simultaneous regions cover mechanism noise
#' and deterministic quantization; they are not sampling confidence intervals.
#'
#' @param data_name Name of the protected data frame.
#' @param variable Numeric variable with custodian-configured bounds.
#' @param server Optional server name; auto-detected when unambiguous.
#' @param datasources DataSHIELD connections.
#' @return A DP moment release with bounds, accuracy, signed Synopsis
#'   provenance, and per-artifact privacy metadata.
#' @export
ds.vertDPMeanVar <- function(data_name, variable, server = NULL,
                             datasources = NULL) {
  resolved <- .dsvert_federation_argument(data_name, datasources)
  .dsvert_dp_meanvar_impl(
    resolved$value, variable, server, resolved$datasources,
    DSI::datashield.aggregate)
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
  run <- .dsvert_dp_synopsis_vector_run(
    datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  public_metadata <- .dsvert_dp_vector_public_metadata(context)
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
  result <- c(public_metadata, list(
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
    epsilon_allocation = c(
      canonical_synopsis_artifact = context$release$epsilon),
    submechanism_count = 1L,
    noise_selection = list(
      winner = if (isTRUE(profile$gaussian)) {
        "two_peer_full_vector_dyadic_discrete_gaussian_tv_bounded"
      } else {
        "two_peer_full_vector_convolution"
      },
      mechanism = context$release$mechanism,
      objective = if (isTRUE(profile$gaussian)) {
        "signed_canonical_synopsis_l2_sensitivity"
      } else {
        "signed_canonical_synopsis_l1_sensitivity"
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

#' @rdname ds.vertDPCount
#' @param x A `ds.vertDPCount` object.
#' @param ... Additional print arguments.
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
