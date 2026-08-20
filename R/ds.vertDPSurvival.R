.dsvert_dp_survival_normalize_hazards <- function(cause_hazard) {
  totals <- rowSums(cause_hazard)
  scale <- ifelse(totals > 1, 1 / totals, 1)
  cause_hazard <- cause_hazard * scale
  list(cause_hazard = cause_hazard,
       all_hazard = pmax(0, rowSums(cause_hazard)))
}

.dsvert_dp_survival_curves <- function(time_grid, entry_counts,
                                       exit_counts, causes) {
  time_count <- length(time_grid)
  if (!is.matrix(exit_counts) || nrow(exit_counts) != time_count ||
      ncol(exit_counts) != length(causes) + 1L ||
      !is.numeric(exit_counts) || anyNA(exit_counts) ||
      any(!is.finite(exit_counts)) || any(exit_counts < 0) ||
      (!is.null(entry_counts) &&
       (!is.numeric(entry_counts) || length(entry_counts) != time_count ||
        anyNA(entry_counts) || any(!is.finite(entry_counts)) ||
        any(entry_counts < 0)))) {
    stop("The DP survival histogram cannot be post-processed", call. = FALSE)
  }
  censor <- exit_counts[, 1L]
  events <- exit_counts[, -1L, drop = FALSE]
  exit_total <- rowSums(exit_counts)
  if (is.null(entry_counts)) {
    at_risk <- rev(cumsum(rev(exit_total)))
  } else {
    entered <- cumsum(entry_counts)
    exited_before <- c(0, head(cumsum(exit_total), -1L))
    at_risk <- pmax(0, exit_total, entered - exited_before)
  }

  cause_hazard <- matrix(
    0, nrow = time_count, ncol = length(causes),
    dimnames = list(NULL, causes))
  positive <- at_risk > 0
  if (any(positive)) {
    cause_hazard[positive, ] <-
      events[positive, , drop = FALSE] / at_risk[positive]
  }
  normalized_hazard <- .dsvert_dp_survival_normalize_hazards(cause_hazard)
  cause_hazard <- normalized_hazard$cause_hazard
  all_hazard <- normalized_hazard$all_hazard
  survival <- numeric(time_count)
  cumulative_hazard <- numeric(time_count)
  cumulative_incidence <- matrix(
    0, nrow = time_count, ncol = length(causes),
    dimnames = list(NULL, causes))
  survival_previous <- 1
  hazard_previous <- 0
  incidence_previous <- numeric(length(causes))
  for (index in seq_len(time_count)) {
    incidence_previous <- incidence_previous +
      survival_previous * cause_hazard[index, ]
    survival_previous <- survival_previous * (1 - all_hazard[[index]])
    hazard_previous <- hazard_previous + all_hazard[[index]]
    survival[[index]] <- max(0, min(1, survival_previous))
    cumulative_hazard[[index]] <- max(0, hazard_previous)
    cumulative_incidence[index, ] <-
      pmax(0, pmin(1, incidence_previous))
  }
  list(
    at_risk = at_risk,
    censor = censor,
    events = events,
    cause_hazard = cause_hazard,
    all_cause_hazard = all_hazard,
    survival = survival,
    cumulative_hazard = cumulative_hazard,
    cumulative_incidence = cumulative_incidence)
}

.dsvert_dp_survival_probability_bounds <- function(
    event_lower, event_upper, risk_lower, risk_upper) {
  if (risk_upper <= 0 || event_upper <= 0) {
    return(c(lower = 0, upper = 0))
  }
  complement_lower <- max(0, risk_lower - event_upper)
  complement_upper <- max(0, risk_upper - event_lower)
  lower_denominator <- event_lower + complement_upper
  upper_denominator <- event_upper + complement_lower
  c(
    lower = if (lower_denominator > 0) {
      event_lower / lower_denominator
    } else 0,
    upper = if (upper_denominator > 0) {
      event_upper / upper_denominator
    } else 1)
}

.dsvert_dp_survival_mechanism_bands <- function(
    time_grid, entry_lower, entry_upper,
    exit_lower, exit_upper, causes) {
  time_count <- length(time_grid)
  outcome_count <- length(causes) + 1L
  valid_matrix <- function(value) {
    is.matrix(value) && is.numeric(value) &&
      identical(dim(value), c(time_count, outcome_count)) &&
      !anyNA(value) && all(is.finite(value)) && all(value >= 0)
  }
  valid_entry <- function(value) {
    is.numeric(value) && length(value) == time_count &&
      !anyNA(value) && all(is.finite(value)) && all(value >= 0)
  }
  delayed_entry <- !is.null(entry_lower) || !is.null(entry_upper)
  if (!is.numeric(time_grid) || !length(time_grid) ||
      anyNA(time_grid) || any(!is.finite(time_grid)) ||
      any(diff(time_grid) <= 0) ||
      !is.character(causes) || !length(causes) || anyNA(causes) ||
      any(!nzchar(causes)) || anyDuplicated(causes) ||
      !valid_matrix(exit_lower) || !valid_matrix(exit_upper) ||
      any(exit_lower > exit_upper) ||
      (delayed_entry &&
       (!valid_entry(entry_lower) || !valid_entry(entry_upper) ||
        any(entry_lower > entry_upper)))) {
    stop("The DP survival histogram box is invalid", call. = FALSE)
  }

  exit_total_lower <- rowSums(exit_lower)
  exit_total_upper <- rowSums(exit_upper)
  if (!delayed_entry) {
    risk_lower <- rev(cumsum(rev(exit_total_lower)))
    risk_upper <- rev(cumsum(rev(exit_total_upper)))
  } else {
    prior_exit_upper <- c(0, head(cumsum(exit_total_upper), -1L))
    prior_exit_lower <- c(0, head(cumsum(exit_total_lower), -1L))
    flow_lower <- cumsum(entry_lower) - prior_exit_upper
    flow_upper <- cumsum(entry_upper) - prior_exit_lower
    risk_lower <- pmax(0, exit_total_lower, flow_lower)
    risk_upper <- pmax(0, exit_total_upper, flow_upper)
  }

  event_lower <- exit_lower[, -1L, drop = FALSE]
  event_upper <- exit_upper[, -1L, drop = FALSE]
  total_event_lower <- rowSums(event_lower)
  total_event_upper <- rowSums(event_upper)
  all_hazard <- t(vapply(seq_len(time_count), function(index) {
    .dsvert_dp_survival_probability_bounds(
      total_event_lower[[index]], total_event_upper[[index]],
      risk_lower[[index]], risk_upper[[index]])
  }, numeric(2L)))
  colnames(all_hazard) <- c("lower", "upper")
  cause_hazard_lower <- cause_hazard_upper <- matrix(
    0, nrow = time_count, ncol = length(causes),
    dimnames = list(NULL, causes))
  for (cause_index in seq_along(causes)) {
    bounds <- t(vapply(seq_len(time_count), function(index) {
      .dsvert_dp_survival_probability_bounds(
        event_lower[index, cause_index],
        event_upper[index, cause_index],
        risk_lower[[index]], risk_upper[[index]])
    }, numeric(2L)))
    cause_hazard_lower[, cause_index] <- bounds[, "lower"]
    cause_hazard_upper[, cause_index] <- bounds[, "upper"]
  }

  survival_lower <- survival_upper <- numeric(time_count)
  cumulative_hazard_lower <- cumulative_hazard_upper <-
    numeric(time_count)
  cumulative_incidence_lower <- cumulative_incidence_upper <- matrix(
    0, nrow = time_count, ncol = length(causes),
    dimnames = list(NULL, causes))
  lower_survival_previous <- upper_survival_previous <- 1
  lower_hazard_previous <- upper_hazard_previous <- 0
  lower_incidence_previous <- upper_incidence_previous <-
    numeric(length(causes))
  for (index in seq_len(time_count)) {
    lower_incidence_previous <- lower_incidence_previous +
      lower_survival_previous * cause_hazard_lower[index, ]
    upper_incidence_previous <- upper_incidence_previous +
      upper_survival_previous * cause_hazard_upper[index, ]
    lower_survival_previous <- lower_survival_previous *
      (1 - all_hazard[index, "upper"])
    upper_survival_previous <- upper_survival_previous *
      (1 - all_hazard[index, "lower"])
    lower_hazard_previous <-
      lower_hazard_previous + all_hazard[index, "lower"]
    upper_hazard_previous <-
      upper_hazard_previous + all_hazard[index, "upper"]
    survival_lower[[index]] <-
      max(0, min(1, lower_survival_previous))
    survival_upper[[index]] <-
      max(survival_lower[[index]], min(1, upper_survival_previous))
    cumulative_hazard_lower[[index]] <- max(0, lower_hazard_previous)
    cumulative_hazard_upper[[index]] <- max(
      cumulative_hazard_lower[[index]], upper_hazard_previous)
    total_upper <- 1 - survival_lower[[index]]
    upper_incidence_previous <- pmin(
      pmax(0, upper_incidence_previous), total_upper)
    lower_incidence_previous <- pmin(
      pmax(0, lower_incidence_previous), upper_incidence_previous)
    cumulative_incidence_lower[index, ] <- lower_incidence_previous
    cumulative_incidence_upper[index, ] <- upper_incidence_previous
  }
  total_incidence_lower <- 1 - survival_upper
  total_incidence_upper <- 1 - survival_lower
  list(
    at_risk = list(lower = risk_lower, upper = risk_upper),
    all_cause_hazard = list(
      lower = all_hazard[, "lower"], upper = all_hazard[, "upper"]),
    cause_specific_hazard = list(
      lower = cause_hazard_lower, upper = cause_hazard_upper),
    kaplan_meier = list(lower = survival_lower, upper = survival_upper),
    nelson_aalen = list(
      lower = cumulative_hazard_lower,
      upper = cumulative_hazard_upper),
    cumulative_incidence = list(
      lower = cumulative_incidence_lower,
      upper = cumulative_incidence_upper),
    total_cumulative_incidence = list(
      lower = total_incidence_lower,
      upper = total_incidence_upper),
    method = paste(
      "simultaneous histogram-coordinate box with interval propagation",
      "through risk sets, hazards, product limits, and Aalen-Johansen"),
    tightness = paste(
      "conservative rectangular relaxation; cross-time flow and cross-cause",
      "feasibility are not jointly optimized"))
}

.dsvert_dp_survival_release_bands <- function(result) {
  radius <- result$accuracy_simultaneous_95_abs
  histogram_lower <- pmax(0, result$histogram - radius)
  histogram_upper <- result$histogram + radius
  time_count <- length(result$time_grid)
  outcome_count <- length(result$causes) + 1L
  cursor <- 1L
  entry_lower <- entry_upper <- NULL
  if (isTRUE(result$delayed_entry)) {
    entry_lower <- histogram_lower[cursor:(cursor + time_count - 1L)]
    entry_upper <- histogram_upper[cursor:(cursor + time_count - 1L)]
    cursor <- cursor + time_count
  }
  exit_length <- time_count * outcome_count
  exit_lower <- matrix(
    histogram_lower[cursor:(cursor + exit_length - 1L)],
    nrow = time_count, ncol = outcome_count)
  exit_upper <- matrix(
    histogram_upper[cursor:(cursor + exit_length - 1L)],
    nrow = time_count, ncol = outcome_count)
  list(
    histogram_lower = histogram_lower,
    histogram_upper = histogram_upper,
    bands = .dsvert_dp_survival_mechanism_bands(
      result$time_grid, entry_lower, entry_upper,
      exit_lower, exit_upper, result$causes))
}

.dsvert_dp_survival_postprocess <- function(result) {
  if (!isTRUE(result$released)) {
    result$entry_histogram <- NULL
    result$exit_histogram <- NULL
    result$not_in_analysis <- NULL
    result$curve <- NULL
    result$cumulative_incidence <- NULL
    result$mechanism_bands <- NULL
    return(result)
  }
  time_count <- length(result$time_grid)
  outcome_levels <- c(result$censor_level, result$causes)
  cursor <- 1L
  entry_counts <- NULL
  if (isTRUE(result$delayed_entry)) {
    entry_counts <- result$histogram[cursor:(cursor + time_count - 1L)]
    cursor <- cursor + time_count
  }
  exit_length <- time_count * length(outcome_levels)
  exit_counts <- matrix(
    result$histogram[cursor:(cursor + exit_length - 1L)],
    nrow = time_count, ncol = length(outcome_levels),
    dimnames = list(result$time_grid, outcome_levels))
  cursor <- cursor + exit_length
  not_in_analysis <- result$histogram[[cursor]]
  derived <- .dsvert_dp_survival_curves(
    result$time_grid, entry_counts, exit_counts, result$causes)
  band_release <- .dsvert_dp_survival_release_bands(result)
  bands <- band_release$bands
  curve <- data.frame(
    time = result$time_grid,
    at_risk_dp = unname(derived$at_risk),
    event_dp = unname(rowSums(derived$events)),
    censor_dp = unname(derived$censor),
    kaplan_meier = unname(derived$survival),
    kaplan_meier_mechanism_lower_95 =
      unname(bands$kaplan_meier$lower),
    kaplan_meier_mechanism_upper_95 =
      unname(bands$kaplan_meier$upper),
    nelson_aalen = unname(derived$cumulative_hazard),
    nelson_aalen_mechanism_lower_95 =
      unname(bands$nelson_aalen$lower),
    nelson_aalen_mechanism_upper_95 =
      unname(bands$nelson_aalen$upper),
    total_cif_mechanism_lower_95 =
      unname(bands$total_cumulative_incidence$lower),
    total_cif_mechanism_upper_95 =
      unname(bands$total_cumulative_incidence$upper),
    stringsAsFactors = FALSE)
  status <- if (isTRUE(result$clipped_coordinates > 0)) {
    "dp_sampler_coordinates_clipped"
  } else if (!any(derived$at_risk > 0)) {
    "dp_curve_empty_after_postprocessing"
  } else if (is.na(result$clipped_coordinates)) {
    "fixed_public_clamp_applied_preclamp_state_not_released"
  } else {
    "ok"
  }
  result$status <- status
  result$entry_histogram <- entry_counts
  result$exit_histogram <- exit_counts
  result$not_in_analysis <- not_in_analysis
  result$curve <- curve
  result$cause_specific_hazard <- derived$cause_hazard
  result$cumulative_incidence <- derived$cumulative_incidence
  result$cumulative_incidence_mechanism_lower_95 <-
    bands$cumulative_incidence$lower
  result$cumulative_incidence_mechanism_upper_95 <-
    bands$cumulative_incidence$upper
  result$histogram_mechanism_lower_95 <-
    band_release$histogram_lower
  result$histogram_mechanism_upper_95 <-
    band_release$histogram_upper
  result$mechanism_bands <- bands
  result$mechanism_band_confidence <- 0.95
  result$mechanism_band_method <- bands$method
  result$mechanism_band_tightness <- bands$tightness
  result$mechanism_band_competing_risk_constraint <- paste(
    "At every grid time the true sum of cause-specific CIFs equals 1-S",
    "and lies in total_cif_mechanism_lower_95/upper_95; causewise interval",
    "corners are not asserted to be jointly attainable")
  result$mechanism_band_scope <- paste(
    "DP mechanism noise only; sampling uncertainty and public-grid",
    "discretisation error excluded")
  result$mechanism_band_additional_privacy_cost <-
    c(epsilon = 0, delta = 0)
  result$mechanism_band_additional_server_calls <- 0L
  result$curve_postprocessing <- paste(
    "risk=max(noisy exits, noisy entry-exit flow);",
    "KM, Nelson-Aalen, and Aalen-Johansen CIF")
  result$statistical_inference <- paste(
    "DP-noised fixed-grid point estimates; no sampling confidence",
    "intervals or p-values are provided")
  result$discretisation_error <- list(
    included_in_mechanism_bands = FALSE,
    status = "not_quantified",
    scope = paste(
      "The fixed public grid changes the time-resolved estimand; no",
      "between-grid interpolation or continuous-time error bound is claimed"))
  result$grid_error_scope <- paste(
    "public time discretisation error is separate from and not included",
    "in the simultaneous DP mechanism bands")
  result
}

.dsvert_dp_survival_vector_result <- function(
    capsule, data_name, analysis_id, server = NULL,
    allow_synopsis = FALSE) {
  capsule <- .dsvert_dp_vector_context(
    capsule, allow_synopsis = allow_synopsis)
  public_metadata <- if (isTRUE(capsule$synopsis)) {
    .dsvert_dp_vector_public_metadata(capsule)
  } else NULL
  release <- capsule$release
  manifest <- release$manifest
  families <- manifest$workload$families
  artifact <- families$survival_artifacts[[analysis_id]]
  scalar <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value)
  }
  time_grid <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$time_grid, "survival time grid"),
    error = function(error) numeric())
  time_bounds <- tryCatch(.dsvert_dp_capsule_manifest_numbers(
    artifact$time_bounds, "survival time bounds"),
    error = function(error) numeric())
  causes <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    artifact$causes, "survival causes"),
    error = function(error) character())
  valid_artifact <- is.list(artifact) &&
    identical(artifact$dataset, data_name) &&
    scalar(artifact$version) && scalar(artifact$owner_peer) &&
    scalar(artifact$time) && scalar(artifact$event) &&
    scalar(artifact$entry) && scalar(artifact$censor) &&
    length(time_grid) > 0L && length(time_bounds) == 2L &&
    time_bounds[[1L]] < time_grid[[1L]] &&
    identical(time_bounds[[2L]], time_grid[[length(time_grid)]]) &&
    all(diff(time_grid) > 0) && length(causes) > 0L &&
    !artifact$censor %in% causes &&
    artifact$owner_peer %in% names(capsule$status)
  if (!isTRUE(valid_artifact)) {
    stop("The signed biomedical release does not contain the requested survival ",
         "artifact", call. = FALSE)
  }
  owner <- artifact$owner_peer
  if (!is.null(server) &&
      (!is.character(server) || length(server) != 1L || is.na(server) ||
       !identical(server, owner))) {
    stop("server does not own the signed survival artifact", call. = FALSE)
  }

  block <- .dsvert_dp_capsule_single_block(
    capsule$layout, "survival_artifacts", dataset = data_name,
    owner_peer = owner,
    predicate = function(candidate) identical(candidate$key, analysis_id),
    description = paste0("survival block for '", analysis_id, "'"))
  delayed_entry <- !identical(artifact$entry, "none")
  time_count <- length(time_grid)
  coordinate_count <- time_count * (length(causes) + 1L) + 1L +
    if (delayed_entry) time_count else 0L
  adjacency <- capsule$status[[owner]]$policy$adjacency
  multiplier <- if (identical(adjacency, "add_remove_patient")) 1L else 2L
  block_l1 <- (if (delayed_entry) 2L else 1L) * multiplier
  block_l2 <- sqrt((if (delayed_entry) 2L else 1L) * multiplier)
  capacity <- as.numeric(capsule$status[[owner]]$policy$unit_capacity)
  artifact_valid <- identical(
    .dsvert_joint_dp_client_json(block$descriptor),
    .dsvert_joint_dp_client_json(artifact)) &&
    block$length == coordinate_count &&
    identical(as.numeric(artifact$coordinate_count),
              as.numeric(coordinate_count)) &&
    identical(artifact$coordinate_order,
              paste0("entry_bins_if_any_then_exit_time_within_outcome_",
                     "then_invalid_bin")) &&
    identical(artifact$repeated_record_policy,
              paste0("earliest_event_else_latest_censor_then_cause_",
                     "then_entry_deterministic_v2")) &&
    identical(artifact$missingness_policy,
              paste0("NA_NaN_Inf_or_out_of_domain_selected_unit_enters_",
                     "invalid_bin")) &&
    .dsvert_dp_num_equal(artifact$l1_sensitivity, block_l1) &&
    .dsvert_dp_num_equal(artifact$l2_sensitivity, block_l2, 2048) &&
    .dsvert_dp_num_equal(artifact$statistic_maximum, capacity)
  if (!isTRUE(artifact_valid)) {
    stop("The signed survival coordinate contract is inconsistent",
         call. = FALSE)
  }
  histogram <- .dsvert_dp_capsule_vector_values(release, block)
  if (any(histogram < 0) || any(histogram > capacity)) {
    stop("The final survival vector violates its signed public clamp",
         call. = FALSE)
  }
  marginal <- .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = 1L, maximum_error = capacity)
  simultaneous <- .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = coordinate_count,
    maximum_error = capacity)
  lattice <- manifest$workload$release_lattice
  if (isTRUE(capsule$synopsis)) {
    return(c(public_metadata, list(
      analysis_id = analysis_id,
      analysis_version = artifact$version,
      coordinate_family = "survival_artifacts",
      coordinate_descriptor = artifact,
      time_grid = time_grid,
      time_lower_bound = time_bounds[[1L]],
      time_upper_bound = time_bounds[[2L]],
      interval_semantics =
        "(previous_endpoint,current_endpoint] after public-bound clipping",
      unit_collapse = "first_event_else_latest_censor_public_tiebreak",
      censor_level = artifact$censor,
      causes = causes,
      delayed_entry = delayed_entry,
      histogram = histogram,
      coordinate_count = coordinate_count,
      histogram_layout = if (delayed_entry) {
        "entry[T],exit[T x (censor+causes) column-major],not_in_analysis"
      } else {
        "exit[T x (censor+causes) column-major],not_in_analysis"
      },
      not_in_analysis_definition = if (delayed_entry) {
        paste("DP-noisy bin for unknown outcome, non-finite time,",
              "non-finite entry, or entry after exit")
      } else {
        "DP-noisy bin for unknown outcome or non-finite time"
      },
      invalid_unit_rule = "invalid_patient_ids_rejected_by_admission",
      artifact_l1_sensitivity = block_l1,
      artifact_l2_sensitivity = block_l2,
      max_histogram_cells_per_unit = if (delayed_entry) 2L else 1L,
      contribution_unit = adjacency,
      clipping_observable = FALSE,
      accuracy_95_abs_per_coordinate = marginal$radius,
      accuracy_simultaneous_95_abs = simultaneous$radius,
      accuracy_simultaneous_confidence = simultaneous$confidence,
      accuracy_simultaneous_method = simultaneous$method,
      accuracy_implementation_tv_upper_bound =
        simultaneous$implementation_tv_upper_bound,
      accuracy_additional_privacy_cost =
        simultaneous$additional_privacy_cost,
      uncertainty_scope =
        "DP mechanism noise only; sampling uncertainty excluded",
      server = owner)))
  }
  noise_root <- capsule$status[[owner]]$noise_root
  profile <- .dsvert_vector_profile(
    manifest$workload$capsule_mechanism,
    manifest$workload$mechanism_selection)
  list(
    released = TRUE, analysis_id = analysis_id,
    analysis_version = artifact$version, time_grid = time_grid,
    time_lower_bound = time_bounds[[1L]],
    time_upper_bound = time_bounds[[2L]],
    interval_semantics =
      "(previous_endpoint,current_endpoint] after public-bound clipping",
    unit_collapse = "first_event_else_latest_censor_public_tiebreak",
    censor_level = artifact$censor, causes = causes,
    delayed_entry = delayed_entry, histogram = histogram,
    coordinate_count = coordinate_count,
    histogram_layout = if (delayed_entry) {
      "entry[T],exit[T x (censor+causes) column-major],not_in_analysis"
    } else {
      "exit[T x (censor+causes) column-major],not_in_analysis"
    },
    not_in_analysis_definition = if (delayed_entry) {
      paste("DP-noisy bin for unknown outcome, non-finite time,",
            "non-finite entry, or entry after exit")
    } else {
      "DP-noisy bin for unknown outcome or non-finite time"
    },
    invalid_unit_rule = "invalid_patient_ids_rejected_by_admission",
    mechanism = release$mechanism,
    implementation = "two pinned peers; Ring128 exact signed finalizer",
    sampler = profile$sampler,
    randomness = "independent HMAC/HKDF/ChaCha20 streams at both peers",
    l1_sensitivity = block_l1, l2_sensitivity = block_l2,
    global_l1_sensitivity =
      as.numeric(lattice$natural_l1_sensitivity),
    global_l2_sensitivity =
      as.numeric(lattice$natural_l2_sensitivity),
    sensitivity_scope = paste(
      "l1/l2 are this signed survival block; global_l1/global_l2",
      "calibrate the one complete capsule vector"),
    max_histogram_cells_per_unit = if (delayed_entry) 2L else 1L,
    contribution_unit = adjacency,
    postprocessing = "fixed public per-coordinate clamp",
    clipped_coordinates = NA_integer_, clipping_observable = FALSE,
    accuracy_95_abs_per_coordinate = marginal$radius,
    accuracy_simultaneous_95_abs = simultaneous$radius,
    accuracy_simultaneous_confidence = simultaneous$confidence,
    accuracy_simultaneous_method = simultaneous$method,
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    privacy_epoch = noise_root$privacy_epoch,
    noise_key_id = noise_root$key_id,
    sticky_noise = "one immutable capsule vector; unlimited replay",
    epsilon = release$epsilon, delta = release$delta,
    implementation_delta = release$implementation_delta,
    adjacency = adjacency,
    capsule_id = release$capsule_id,
    final_vector_root = release$final_vector_root,
    coordinate_order_sha256 = release$coordinate_order_sha256,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    security_claim = .dsvert_dp_capsule_security_claim(),
    server = owner)
}

#' Differentially private non-parametric survival curves
#'
#' Reads one canonical signed Synopsis release of a fixed, custodian-owned time
#' histogram and derives
#' Kaplan--Meier survival, Nelson--Aalen cumulative hazard, and competing-risk
#' cumulative incidence by post-processing that same release. Repeating the
#' same analysis id over the same snapshot replays the same sticky result
#' without lifetime admission.
#'
#' @param data_name Name of the registered protected data frame.
#' @param analysis_id Custodian-owned survival specification id.
#' @param server Optional datasource name. If omitted, the lexicographically
#'   first connected datasource is selected deterministically.
#' @param datasources DataSHIELD connections.
#' @return A `ds.vertDPSurvival` object with signed Synopsis provenance and
#'   per-artifact privacy metadata. Accuracy covers DP mechanism noise in
#'   histogram coordinates only, not sampling uncertainty or public-grid
#'   discretisation.
#' @export
ds.vertDPSurvival <- function(data_name, analysis_id, server = NULL,
                              datasources = NULL) {
  resolved <- .dsvert_federation_argument(data_name, datasources)
  .dsvert_dp_survival_impl(
    resolved$value, analysis_id, server, resolved$datasources,
    DSI::datashield.aggregate)
}

.dsvert_dp_survival_impl <- function(data_name, analysis_id, server = NULL,
                                     datasources = NULL, .aggregate) {
  for (value in list(data_name, analysis_id)) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop("data_name and analysis_id must be non-empty strings",
           call. = FALSE)
    }
  }
  datasources <- .dsvert_dp_datasources(datasources)
  capsule <- .dsvert_dp_synopsis_vector_run(
    datasources, .aggregate = .aggregate)
  result <- .dsvert_dp_survival_vector_result(
    capsule, data_name, analysis_id, server, allow_synopsis = TRUE)
  result <- .dsvert_dp_survival_postprocess(result)
  class(result) <- c("ds.vertDPSurvival", "list")
  result
}

.dsvert_dp_survival_is_synopsis <- function(x) {
  is.list(x) && is.list(x$release_provenance) &&
    identical(
      x$release_provenance$version,
      "dsvert-stateless-synopsis-public-provenance-v1")
}

.dsvert_dp_survival_synopsis_binding_valid <- function(x) {
  if (!.dsvert_dp_survival_is_synopsis(x)) return(FALSE)
  bindings <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  all(vapply(bindings, function(field) {
    .dsvert_vector_hex(x[[field]]) &&
      identical(x[[field]], x$release_provenance[[field]])
  }, logical(1L))) &&
    identical(x$release_provenance$durable_replay, TRUE) &&
    is.list(x$privacy) && identical(x$privacy$unlimited_replay, TRUE) &&
    !any(c(
      "capsule_id", "privacy_epoch", "noise_key_id", "history_gate",
      "request_limit", "operation_limit") %in% names(x))
}

.dsvert_dp_survival_source_provenance <- function(x) {
  fields <- if (.dsvert_dp_survival_is_synopsis(x)) {
    c(
      "analysis_id", "analysis_version", "server", "artifact_key",
      "execution_id", "manifest_sha256", "contract_sha256",
      "attempt_sha256", "source_contract_sha256", "result_set_sha256",
      "final_vector_root", "coordinate_order_sha256", "release_provenance",
      "privacy", "mechanism", "implementation", "sampler", "epsilon",
      "delta", "implementation_delta", "adjacency", "time_grid",
      "time_lower_bound", "time_upper_bound", "security_claim")
  } else {
    c(
      "analysis_id", "analysis_version", "server", "capsule_id",
      "final_vector_root", "coordinate_order_sha256", "privacy_epoch",
      "noise_key_id", "mechanism", "implementation", "sampler", "epsilon",
      "delta", "implementation_delta", "adjacency", "time_grid",
      "time_lower_bound", "time_upper_bound", "security_claim")
  }
  c(list(source_class = "ds.vertDPSurvival"), x[fields])
}

.dsvert_dp_survival_object <- function(x) {
  invalid <- function() {
    stop("x must be a validated released ds.vertDPSurvival object",
         call. = FALSE)
  }
  scalar_number <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value)
  }
  current_vector <- is.list(x) && x$mechanism %in% c(
    .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
    .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
  synopsis <- .dsvert_dp_survival_is_synopsis(x)
  if (!inherits(x, "ds.vertDPSurvival") || !is.list(x) ||
      !isTRUE(x$released) || !is.numeric(x$time_grid) ||
      !length(x$time_grid) || anyNA(x$time_grid) ||
      any(!is.finite(x$time_grid)) || any(diff(x$time_grid) <= 0) ||
      !scalar_number(x$time_lower_bound) ||
      !scalar_number(x$time_upper_bound) ||
      x$time_lower_bound >= x$time_grid[[1L]] ||
      !.dsvert_dp_num_equal(
        x$time_upper_bound, x$time_grid[[length(x$time_grid)]]) ||
      !is.character(x$causes) || !length(x$causes) || anyNA(x$causes) ||
      any(!nzchar(x$causes)) || anyDuplicated(x$causes) ||
      !is.character(x$censor_level) || length(x$censor_level) != 1L ||
      is.na(x$censor_level) || !nzchar(x$censor_level) ||
      x$censor_level %in% x$causes ||
      !is.logical(x$delayed_entry) || length(x$delayed_entry) != 1L ||
      is.na(x$delayed_entry) ||
      !is.numeric(x$histogram) || anyNA(x$histogram) ||
      any(!is.finite(x$histogram)) || any(x$histogram < 0) ||
      (isTRUE(current_vector) && !isTRUE(synopsis) &&
       (!identical(x$history_gate, TRUE) ||
        !identical(x$request_limit, FALSE) ||
        !identical(x$operation_limit, TRUE))) ||
      (isTRUE(synopsis) &&
       !.dsvert_dp_survival_synopsis_binding_valid(x)) ||
      !scalar_number(x$accuracy_simultaneous_95_abs) ||
      x$accuracy_simultaneous_95_abs < 0 ||
      !is.data.frame(x$curve) || !is.matrix(x$cumulative_incidence)) {
    invalid()
  }
  time_count <- length(x$time_grid)
  expected_coordinates <- time_count * (length(x$causes) + 1L) + 1L +
    if (isTRUE(x$delayed_entry)) time_count else 0L
  if (length(x$histogram) != expected_coordinates) invalid()

  # The source release is authenticated before this R object is created. This
  # second check protects client-only algebra against accidental or partial
  # in-memory edits by requiring every consumed derived field to remain the
  # deterministic post-processing of that release. It is not a replacement
  # for the signed release boundary.
  recomputed <- tryCatch(
    .dsvert_dp_survival_postprocess(x),
    error = function(error) NULL)
  derived_fields <- c(
    "status", "entry_histogram", "exit_histogram", "not_in_analysis",
    "curve", "cause_specific_hazard", "cumulative_incidence",
    "cumulative_incidence_mechanism_lower_95",
    "cumulative_incidence_mechanism_upper_95",
    "histogram_mechanism_lower_95", "histogram_mechanism_upper_95",
    "mechanism_bands", "mechanism_band_confidence",
    "mechanism_band_method", "mechanism_band_tightness",
    "mechanism_band_competing_risk_constraint", "mechanism_band_scope",
    "mechanism_band_additional_privacy_cost",
    "mechanism_band_additional_server_calls", "curve_postprocessing",
    "statistical_inference", "discretisation_error", "grid_error_scope")
  consistent <- !is.null(recomputed) && all(vapply(
    derived_fields,
    function(field) identical(x[[field]], recomputed[[field]]),
    logical(1L)))
  if (!isTRUE(consistent)) {
    invalid()
  }
  x
}

#' Kaplan--Meier curve from one DP survival release
#'
#' @param x A released `ds.vertDPSurvival` object.
#' @return A data frame containing public grid time, DP-derived survival, and
#'   conservative simultaneous mechanism-band limits.
#' @export
ds.vertDPKaplanMeier <- function(x) {
  x <- .dsvert_dp_survival_object(x)
  result <- x$curve[c(
    "time", "kaplan_meier",
    "kaplan_meier_mechanism_lower_95",
    "kaplan_meier_mechanism_upper_95")]
  attr(result, "uncertainty_scope") <- x$uncertainty_scope
  attr(result, "mechanism_band_scope") <- x$mechanism_band_scope
  attr(result, "mechanism_band_tightness") <- x$mechanism_band_tightness
  attr(result, "grid_error_scope") <- x$grid_error_scope
  attr(result, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  attr(result, "additional_server_calls") <- 0L
  result
}

#' Nelson--Aalen curve from one DP survival release
#'
#' @param x A released `ds.vertDPSurvival` object.
#' @return A data frame containing public grid time, cumulative hazard, and
#'   conservative simultaneous mechanism-band limits.
#' @export
ds.vertDPNelsonAalen <- function(x) {
  x <- .dsvert_dp_survival_object(x)
  result <- x$curve[c(
    "time", "nelson_aalen",
    "nelson_aalen_mechanism_lower_95",
    "nelson_aalen_mechanism_upper_95")]
  attr(result, "uncertainty_scope") <- x$uncertainty_scope
  attr(result, "mechanism_band_scope") <- x$mechanism_band_scope
  attr(result, "mechanism_band_tightness") <- x$mechanism_band_tightness
  attr(result, "grid_error_scope") <- x$grid_error_scope
  attr(result, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  attr(result, "additional_server_calls") <- 0L
  result
}

#' Competing-risks cumulative incidence from one DP survival release
#'
#' @param x A released `ds.vertDPSurvival` object.
#' @param cause Optional public cause label. `NULL` returns every cause.
#' @return A long data frame of public grid time, cause, cumulative incidence,
#'   and conservative simultaneous mechanism-band limits.
#' @export
ds.vertDPCumulativeIncidence <- function(x, cause = NULL) {
  x <- .dsvert_dp_survival_object(x)
  causes <- colnames(x$cumulative_incidence)
  if (!is.null(cause)) {
    if (!is.character(cause) || length(cause) != 1L || is.na(cause) ||
        !cause %in% causes) {
      stop("cause must name one released event cause", call. = FALSE)
    }
    causes <- cause
  }
  selected <- x$cumulative_incidence[, causes, drop = FALSE]
  result <- data.frame(
    time = rep(x$time_grid, times = length(causes)),
    cause = rep(causes, each = length(x$time_grid)),
    cumulative_incidence = as.vector(selected),
    cumulative_incidence_mechanism_lower_95 = as.vector(
      x$cumulative_incidence_mechanism_lower_95[, causes, drop = FALSE]),
    cumulative_incidence_mechanism_upper_95 = as.vector(
      x$cumulative_incidence_mechanism_upper_95[, causes, drop = FALSE]),
    total_cif_mechanism_lower_95 = rep(
      x$curve$total_cif_mechanism_lower_95, times = length(causes)),
    total_cif_mechanism_upper_95 = rep(
      x$curve$total_cif_mechanism_upper_95, times = length(causes)),
    stringsAsFactors = FALSE)
  attr(result, "uncertainty_scope") <- x$uncertainty_scope
  attr(result, "mechanism_band_scope") <- x$mechanism_band_scope
  attr(result, "mechanism_band_tightness") <- x$mechanism_band_tightness
  attr(result, "competing_risk_constraint") <-
    x$mechanism_band_competing_risk_constraint
  attr(result, "grid_error_scope") <- x$grid_error_scope
  attr(result, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  attr(result, "additional_server_calls") <- 0L
  result
}

.dsvert_dp_rmst_step_weights <- function(time_grid, time_lower_bound, tau) {
  starts <- c(time_lower_bound, head(time_grid, -1L))
  pmax(0, pmin(time_grid, tau) - starts)
}

#' Restricted mean survival time from one DP survival release
#'
#' Computes restricted mean survival time (RMST) by integrating the
#' left-continuous product-limit step curve on the release's fixed public time
#' grid. This is pure post-processing: it performs no server call and consumes
#' no additional privacy budget. Its limits propagate the simultaneous DP
#' mechanism band; they are not sampling confidence intervals and do not cover
#' error from replacing continuous event times by the public grid. The result
#' carries a copy of the source release's signed DP provenance.
#'
#' @param x A released `ds.vertDPSurvival` object.
#' @param tau One or more public finite restriction times greater than the
#'   release's lower time bound and no greater than its upper time bound.
#'   `NULL` uses the public upper bound.
#' @return A data frame containing RMST, conservative simultaneous
#'   DP-mechanism limits, and source-release provenance for every requested
#'   `tau`.
#' @export
ds.vertDPRMST <- function(x, tau = NULL) {
  x <- .dsvert_dp_survival_object(x)
  if (is.null(tau)) tau <- x$time_upper_bound
  if (!is.numeric(tau) || !length(tau) || anyNA(tau) ||
      any(!is.finite(tau)) || any(tau <= x$time_lower_bound) ||
      any(tau > x$time_upper_bound)) {
    stop(
      "tau must contain finite times above the public lower bound and at or below the public upper bound",
      call. = FALSE)
  }

  # With events assigned to public interval endpoints, survival is one before
  # the first endpoint and S_j on [t_j, t_{j+1}). The final grid-point value
  # has zero integration width at the public upper bound.
  point_levels <- c(1, head(x$curve$kaplan_meier, -1L))
  lower_levels <- c(
    1, head(x$curve$kaplan_meier_mechanism_lower_95, -1L))
  upper_levels <- c(
    1, head(x$curve$kaplan_meier_mechanism_upper_95, -1L))
  values <- lapply(as.numeric(tau), function(limit) {
    weights <- .dsvert_dp_rmst_step_weights(
      x$time_grid, x$time_lower_bound, limit)
    c(
      rmst = sum(weights * point_levels),
      lower = sum(weights * lower_levels),
      upper = sum(weights * upper_levels))
  })
  values <- do.call(rbind, values)
  result <- data.frame(
    time_lower_bound = rep(x$time_lower_bound, length(tau)),
    tau = as.numeric(tau),
    rmst = values[, "rmst"],
    rmst_mechanism_lower_95 = values[, "lower"],
    rmst_mechanism_upper_95 = values[, "upper"],
    stringsAsFactors = FALSE)
  attr(result, "uncertainty_scope") <- x$uncertainty_scope
  attr(result, "mechanism_band_scope") <- x$mechanism_band_scope
  attr(result, "mechanism_band_tightness") <-
    x$mechanism_band_tightness
  attr(result, "mechanism_band_confidence") <-
    x$mechanism_band_confidence
  attr(result, "mechanism_band_method") <- x$mechanism_band_method
  attr(result, "grid_error_scope") <- x$grid_error_scope
  attr(result, "statistical_inference") <- x$statistical_inference
  attr(result, "discretisation_error") <- x$discretisation_error
  attr(result, "integration_rule") <- paste(
    "left-continuous product-limit step curve; events occur at public",
    "interval endpoints")
  attr(result, "additional_privacy_cost") <- c(epsilon = 0, delta = 0)
  attr(result, "additional_server_calls") <- 0L
  attr(result, "source_release_provenance") <-
    .dsvert_dp_survival_source_provenance(x)
  class(result) <- c("ds.vertDPRMST", class(result))
  result
}

#' Restricted mean time lost from one DP survival release
#'
#' Computes restricted mean time lost (RMTL) as the exact complement of
#' `ds.vertDPRMST()` on the release's public interval:
#' `RMTL = (tau - time_lower_bound) - RMST`. The simultaneous mechanism limits
#' are transformed by the same identity with their order reversed. No curve is
#' re-estimated, no server is contacted, and no additional privacy is spent.
#' The inherited limits exclude sampling uncertainty and public-grid
#' discretisation error.
#'
#' @param x A validated released `ds.vertDPSurvival` object.
#' @param tau One or more public finite restriction times greater than the
#'   release's lower time bound and no greater than its upper time bound.
#'   `NULL` uses the public upper bound.
#' @return A `ds.vertDPRMTL` data frame retaining the complete RMST result and
#'   provenance, plus restriction width, RMTL, and its conservative
#'   simultaneous DP-mechanism limits.
#' @export
ds.vertDPRMTL <- function(x, tau = NULL) {
  result <- ds.vertDPRMST(x, tau)
  restriction_width <- result$tau - result$time_lower_bound
  result$restriction_width <- restriction_width
  result$rmtl <- restriction_width - result$rmst
  result$rmtl_mechanism_lower_95 <-
    restriction_width - result$rmst_mechanism_upper_95
  result$rmtl_mechanism_upper_95 <-
    restriction_width - result$rmst_mechanism_lower_95
  attr(result, "complement_identity") <- paste(
    "RMTL = restriction_width - RMST, where restriction_width =",
    "tau - time_lower_bound; mechanism limits use width-upper and",
    "width-lower without recalibration")
  class(result) <- c("ds.vertDPRMTL", class(result))
  result
}

#' @export
print.ds.vertDPSurvival <- function(x, ...) {
  if (!isTRUE(x$released)) {
    cat("dsVert DP survival: suppressed (", x$reason, ")\n", sep = "")
  } else {
    cat("dsVert DP survival:", x$analysis_id, "[", x$server, "]\n")
    cat("grid points:", length(x$time_grid), "| causes:",
        length(x$causes), "| status:", x$status, "\n")
    cat("epsilon:", format(x$epsilon),
        "| simultaneous 95% histogram-noise radius:",
        x$accuracy_simultaneous_95_abs, "\n")
    cat("curve bands:", x$mechanism_band_tightness, "\n")
    cat("Sampling confidence intervals and grid-discretisation error are ",
        "not included.\n", sep = "")
  }
  invisible(x)
}
