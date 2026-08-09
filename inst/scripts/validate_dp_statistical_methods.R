#!/usr/bin/env Rscript

# Reproducible finite-dataset validation of stable DP post-processors.
#
# This is a distributional simulation, not a DataSHIELD/DSI integration test.
# It uses the deployed signed-vector accuracy calculation and the ideal
# two-peer two-sided-geometric convolution to generate mechanism noise. The
# productive HMAC/HKDF/ChaCha20 finite sampler is not reimplemented here; its
# separately certified TV bound is retained by the deployed radius calculator.
# Gaussian regression coefficients and PCA loading directions are validated as
# point estimators only because the public API deliberately claims no regions
# for either nonlinear quantity.

.dv_cli_args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L && "--help" %in% .dv_cli_args) {
  cat(paste(
    "Usage: validate_dp_statistical_methods.R [--quick] [--replicates=N]",
    "[--output-dir=PATH] [--preflight]\n"
  ))
  quit(save = "no", status = 0L)
}
if (!requireNamespace("dsVertClient", quietly = TRUE)) {
  stop("Install dsVertClient before running this script", call. = FALSE)
}
if (sys.nframe() == 0L && "--preflight" %in% .dv_cli_args) {
  cat("dsVertClient preflight OK: ",
      as.character(utils::packageVersion("dsVertClient")), "\n", sep = "")
  quit(save = "no", status = 0L)
}

.dv_namespace <- asNamespace("dsVertClient")
.dv_get <- function(name) get(name, envir = .dv_namespace, inherits = FALSE)
.dv_laplace_mechanism <- .dv_get(
  ".DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM")
.dv_laplace_sampler <- .dv_get(".DSVERT_CLIENT_VECTOR_SAMPLER")
.dv_laplace_backend <- .dv_get(".DSVERT_CLIENT_VECTOR_BACKEND")
.dv_implementation_delta <-
  "1/1267650600228229401496703205376"
.dv_level <- 0.95
.dv_numeric_tolerance <- 1e-10

.dv_accuracy_contract <- function(
    family, epsilon, natural_l1_sensitivity, coordinate_count,
    scale = 256, maximum_error = Inf) {
  manifest <- list(workload = list(
    coordinate_count = as.integer(coordinate_count),
    capsule_mechanism = list(
      mechanism = "discrete-laplace", sensitivity_norm = "l1"),
    mechanism_selection = NULL,
    release_lattice = list(
      output_lattice_bits = as.integer(log2(scale)),
      output_lattice_scale = as.numeric(scale),
      natural_l1_sensitivity = as.numeric(natural_l1_sensitivity),
      integer_l1_sensitivity_steps =
        as.numeric(natural_l1_sensitivity * scale),
      natural_l2_sensitivity = sqrt(natural_l1_sensitivity),
      integer_l2_sensitivity_steps =
        sqrt(natural_l1_sensitivity) * scale)))
  release <- list(
    epsilon = as.numeric(epsilon),
    mechanism = .dv_laplace_mechanism,
    implementation_delta = .dv_implementation_delta)
  radius <- .dv_get(".dsvert_dp_vector_accuracy_radius")(
    release, manifest, coordinate_count = coordinate_count,
    confidence = .dv_level, maximum_error = maximum_error)
  steps <- radius$radius * scale
  if (abs(steps - round(steps)) >
      64 * .Machine$double.eps * max(1, abs(steps))) {
    stop("The deployed accuracy radius is not on the release lattice",
         call. = FALSE)
  }
  log_tail <- .dv_get(".dsvert_dp_vector_convolution_log_tail")(
    round(steps), epsilon, natural_l1_sensitivity * scale)
  failure_upper <- coordinate_count * exp(log_tail) +
    radius$implementation_tv_upper_bound
  list(
    family = family,
    epsilon = epsilon,
    natural_l1_sensitivity = natural_l1_sensitivity,
    coordinate_count = as.integer(coordinate_count),
    scale = scale,
    radius = radius$radius,
    accuracy_method = radius$method,
    implementation_tv_upper_bound =
      radius$implementation_tv_upper_bound,
    certified_failure_probability_upper = failure_upper,
    certified_coverage_lower = max(0, 1 - failure_upper))
}

.dv_draw_one_peer <- function(n, contract) {
  success_probability <- -expm1(
    -contract$epsilon /
      (contract$natural_l1_sensitivity * contract$scale))
  zero_probability <- success_probability / (2 - success_probability)
  nonzero <- stats::runif(n) >= zero_probability
  result <- numeric(n)
  count <- sum(nonzero)
  if (count) {
    magnitude <- stats::rgeom(count, prob = success_probability) + 1
    sign <- ifelse(stats::runif(count) < 0.5, -1, 1)
    result[nonzero] <- sign * magnitude
  }
  result
}

.dv_draw_noise <- function(n, contract) {
  (.dv_draw_one_peer(n, contract) +
     .dv_draw_one_peer(n, contract)) / contract$scale
}

.dv_clamp <- function(value, maximum) pmin(maximum, pmax(0, value))

.dv_inside <- function(truth, lower, upper) {
  tolerance <- .dv_numeric_tolerance *
    pmax(1, abs(truth), abs(lower), abs(upper))
  truth >= lower - tolerance & truth <= upper + tolerance
}

.dv_record <- function(
    family, method, estimand, replicate, truth, estimate, lower, upper,
    coordinate_event, mechanism_region_available = TRUE) {
  n <- max(
    length(estimand), length(truth), length(estimate), length(lower),
    length(upper))
  data.frame(
    family = rep_len(family, n),
    method = rep_len(method, n),
    estimand = rep_len(estimand, n),
    replicate = rep_len(as.integer(replicate), n),
    truth = rep_len(as.numeric(truth), n),
    estimate = rep_len(as.numeric(estimate), n),
    lower = rep_len(as.numeric(lower), n),
    upper = rep_len(as.numeric(upper), n),
    coordinate_event = rep_len(isTRUE(coordinate_event), n),
    mechanism_region_available = rep_len(
      isTRUE(mechanism_region_available), n),
    stringsAsFactors = FALSE)
}

.dv_describe_exact <- function(values, lower, upper, grid, scale) {
  bounded <- pmin(upper, pmax(lower, values))
  normalized <- (bounded - lower) / (upper - lower)
  histogram <- tabulate(
    as.integer(findInterval(bounded, grid, left.open = TRUE) + 1L),
    nbins = length(grid))
  list(
    bounded = bounded,
    natural_coordinates = c(
      length(bounded), sum(round(normalized * scale)) / scale,
      sum(round(normalized^2 * scale)) / scale,
      histogram, 0),
    histogram = histogram,
    mean = mean(bounded),
    variance = mean((bounded - mean(bounded))^2))
}

.dv_describe_release <- function(
    natural_coordinates, lower, upper, grid, capacity, contract,
    probs = c(0.25, 0.5, 0.75)) {
  scale <- contract$scale
  grid_length <- length(grid)
  marginal <- .dv_accuracy_contract(
    paste0(contract$family, "_marginal"), contract$epsilon,
    contract$natural_l1_sensitivity, 1L, scale, capacity)
  statistics <- c(
    round(natural_coordinates[[1L]]),
    round(natural_coordinates[[2L]] * scale),
    round(natural_coordinates[[3L]] * scale),
    natural_coordinates[4L:(grid_length + 4L)])
  result <- list(
    released = TRUE, analysis_id = "validation_describe",
    analysis_version = "v1", variables = "marker",
    variable_count = 1L, lower_bounds = lower, upper_bounds = upper,
    grid_lengths = as.integer(grid_length), grid_values = grid,
    histogram_semantics =
      "(previous_endpoint,current_endpoint] plus fixed invalid bin",
    unit_collapse =
      "mean_of_finite_rows_after_public_bound_clipping",
    count_definition =
      "DP-noisy effective units with at least one finite bounded value",
    invalid_unit_rule = "invalid_patient_ids_rejected_by_admission",
    statistics = statistics,
    coordinate_count = as.integer(length(statistics)),
    coordinate_layout = paste(
      "referenced capsule blocks per variable[count,qsum,qsumsq,",
      "histogram[grid_bins+invalid]]"),
    numeric_grid_bits = as.integer(log2(scale)),
    numeric_grid_scale = scale,
    quantization = "round(z*scale) and round(z^2*scale) independently",
    max_abs_normalized_quantization_per_unit = 0.5 / scale,
    allocation_names = c("count", "sum", "sumsq", "histogram"),
    allocation_weights = rep(0.25, 4L),
    mechanism = .dv_laplace_mechanism,
    implementation = "two pinned peers; Ring128 exact signed finalizer",
    sampler = .dv_laplace_sampler,
    randomness = "independent HMAC/HKDF/ChaCha20 streams at both peers",
    postprocessing = "fixed public per-coordinate clamp",
    clipped_coordinates = NA_integer_, clipping_observable = FALSE,
    accuracy_95_abs_by_variable_family = c(
      marginal$radius, marginal$radius * scale,
      marginal$radius * scale, marginal$radius),
    accuracy_simultaneous_95_abs_by_variable_family = c(
      contract$radius, contract$radius * scale,
      contract$radius * scale, contract$radius),
    accuracy_simultaneous_confidence = .dv_level,
    accuracy_simultaneous_method = contract$accuracy_method,
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    privacy_epoch = 1, noise_key_id = "validation-noise-key-a",
    sticky_noise = "one immutable capsule vector; unlimited replay",
    epsilon = contract$epsilon, delta = 2^-100,
    implementation_delta = .dv_implementation_delta,
    adjacency = "add_remove_patient",
    capsule_id = strrep("a", 64L), final_vector_root = strrep("b", 64L),
    coordinate_order_sha256 = strrep("c", 64L), server = "site_a",
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE)
  result <- .dv_get(".dsvert_dp_describe_postprocess")(result, probs)
  class(result) <- c("ds.vertDPDescribe", "list")
  result
}

.dv_run_describe <- function(replicates) {
  seed <- 20260802L
  set.seed(seed)
  lower <- 0
  upper <- 4
  grid <- seq(0.5, 4, by = 0.5)
  capacity <- 600
  # Keep the main quantile ranks away from exact CDF/bin boundaries. The
  # discontinuous exact-tie rule is exercised separately as a boundary case.
  bin_counts <- c(38L, 56L, 49L, 61L, 43L, 57L, 39L, 57L)
  values <- rep(grid - 0.25, times = bin_counts)
  exact <- .dv_describe_exact(values, lower, upper, grid, 256)
  contract <- .dv_accuracy_contract(
    "describe_quantile_median", epsilon = 6,
    natural_l1_sensitivity = 5, coordinate_count =
      length(exact$natural_coordinates), maximum_error = capacity)
  probs <- c(0.25, 0.5, 0.75)
  cumulative <- cumsum(exact$histogram)
  oracle_bins <- vapply(probs, function(probability) {
    which(cumulative >= probability * sum(exact$histogram))[[1L]]
  }, integer(1L))
  oracle_endpoints <- grid[oracle_bins]

  no_noise <- .dv_describe_release(
    exact$natural_coordinates, lower, upper, grid, capacity, contract, probs)
  no_noise_quantile <- .dv_get("ds.vertDPQuantile")(no_noise, probs)
  quantized_mean <- lower + (upper - lower) *
    exact$natural_coordinates[[2L]] / length(values)
  quantized_second <- exact$natural_coordinates[[3L]] / length(values)
  quantized_variance <- (upper - lower)^2 * max(
    0, quantized_second -
      (exact$natural_coordinates[[2L]] / length(values))^2)
  oracle_error <- max(
    abs(no_noise$descriptives$mean - quantized_mean),
    abs(no_noise$descriptives$variance - quantized_variance),
    abs(no_noise_quantile$bin_upper - oracle_endpoints))
  quantization_coverage <-
    .dv_inside(
      exact$mean,
      no_noise$descriptives$mean_mechanism_grid_lower_95,
      no_noise$descriptives$mean_mechanism_grid_upper_95) &&
    .dv_inside(
      exact$variance,
      no_noise$descriptives$variance_mechanism_grid_lower_95,
      no_noise$descriptives$variance_mechanism_grid_upper_95)

  records <- vector("list", replicates * 6L)
  identity_error <- 0
  position <- 0L
  started <- proc.time()[["elapsed"]]
  for (iteration in seq_len(replicates)) {
    noisy <- .dv_clamp(
      exact$natural_coordinates +
        .dv_draw_noise(length(exact$natural_coordinates), contract),
      capacity)
    coordinate_event <- all(
      abs(noisy - exact$natural_coordinates) <=
        contract$radius + .dv_numeric_tolerance)
    release <- .dv_describe_release(
      noisy, lower, upper, grid, capacity, contract, probs)
    quantiles <- .dv_get("ds.vertDPQuantile")(release, probs)
    median <- .dv_get("ds.vertDPMedian")(release)
    identity_error <- max(
      identity_error,
      abs(median$bin_index -
            quantiles$bin_index[quantiles$probability == 0.5]))

    position <- position + 1L
    records[[position]] <- .dv_record(
      "describe", "ds.vertDPDescribe", "bounded_mean", iteration,
      exact$mean, release$descriptives$mean,
      release$descriptives$mean_mechanism_grid_lower_95,
      release$descriptives$mean_mechanism_grid_upper_95,
      coordinate_event)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "describe", "ds.vertDPDescribe", "bounded_population_variance",
      iteration, exact$variance, release$descriptives$variance,
      release$descriptives$variance_mechanism_grid_lower_95,
      release$descriptives$variance_mechanism_grid_upper_95,
      coordinate_event)
    for (probability_index in seq_along(probs)) {
      position <- position + 1L
      row <- quantiles[probability_index, ]
      records[[position]] <- .dv_record(
        "describe", "ds.vertDPQuantile",
        paste0("q", 100 * probs[[probability_index]],
               "_public_bin_endpoint"), iteration,
        oracle_endpoints[[probability_index]], row$bin_upper,
        row$mechanism_grid_lower_95, row$mechanism_grid_upper_95,
        coordinate_event)
    }
    position <- position + 1L
    records[[position]] <- .dv_record(
      "describe", "ds.vertDPMedian", "median_public_bin_endpoint",
      iteration, oracle_endpoints[[2L]], median$bin_upper,
      median$mechanism_grid_lower_95, median$mechanism_grid_upper_95,
      coordinate_event)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    records = do.call(rbind, records[seq_len(position)]),
    contract = contract,
    elapsed = elapsed,
    seed = seed,
    oracle_error = oracle_error,
    quantization_coverage = quantization_coverage,
    identity_error = identity_error)
}

.dv_table_release <- function(table, contract, coordinate_maximum = 2000) {
  table <- as.matrix(table)
  result <- list(
    released = TRUE, mechanism = .dv_laplace_mechanism,
    implementation = .dv_laplace_backend,
    backend = "exact_signed_Ring128_global_vector",
    sampler = .dv_laplace_sampler,
    randomness = paste(
      "independent pinned-peer HKDF-SHA256/ChaCha20 streams;",
      "no analyst-controlled seed"),
    epsilon = contract$epsilon, delta = 2^-100,
    implementation_delta = .dv_implementation_delta,
    adjacency = "add_remove_patient",
    sensitivity = contract$natural_l1_sensitivity,
    sensitivity_norm = "l1",
    l1_sensitivity = contract$natural_l1_sensitivity,
    l2_sensitivity = sqrt(2),
    sensitivity_scope = "complete_signed_biomedical_capsule_vector",
    output_lattice_bits = as.integer(log2(contract$scale)),
    output_lattice_scale = contract$scale,
    sticky_noise = "immutable_capsule_durable_replay_v3",
    sticky_replay = TRUE, privacy_epochs = c(1, 1),
    noise_key_ids = c("validation-noise-a", "validation-noise-b"),
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    clipped_coordinates = NA_integer_, clamp_activation_disclosed = FALSE,
    table = table, counts = unname(as.numeric(table)),
    nrow = as.integer(nrow(table)), ncol = as.integer(ncol(table)),
    row_levels = unname(rownames(table)),
    col_levels = unname(colnames(table)),
    coordinate_maximum = coordinate_maximum,
    artifact_l1_sensitivity = 1, artifact_l2_sensitivity = 1,
    unit_aggregation_policy = "consistent_joint_cell_else_exclude_v1",
    server = "site_a", accuracy_simultaneous_confidence = .dv_level,
    accuracy_simultaneous_method = contract$accuracy_method,
    accuracy_additional_privacy_cost = c(epsilon = 0, delta = 0),
    accuracy_simultaneous_95_abs = contract$radius)
  class(result) <- c("ds.vertDPContingency", "list")
  result
}

.dv_epi_oracle <- function(table, weights, expected_rates) {
  exposed_risk <- table["exposed", "event"] / sum(table["exposed", ])
  unexposed_risk <-
    table["unexposed", "event"] / sum(table["unexposed", ])
  a <- table["exposed", "event"]
  b <- table["exposed", "nonevent"]
  c0 <- table["unexposed", "event"]
  d <- table["unexposed", "nonevent"]
  population_risk <- (a + c0) / sum(table)
  diagnostic <- c(
    sensitivity = exposed_risk,
    specificity = d / (c0 + d),
    ppv = a / (a + c0),
    npv = d / (b + d),
    prevalence = (a + b) / sum(table),
    accuracy = (a + d) / sum(table),
    balanced_accuracy = 0.5 * (
      exposed_risk + d / (c0 + d)),
    f1_score = 2 * a / (2 * a + b + c0),
    lr_positive = exposed_risk / (c0 / (c0 + d)),
    lr_negative = (b / (a + b)) / (d / (c0 + d)),
    diagnostic_odds_ratio = a * d / (b * c0))
  ordered_cases <- table["exposed", c("nonevent", "event")]
  ordered_controls <- table["unexposed", c("nonevent", "event")]
  controls_below <- c(0, head(cumsum(ordered_controls), -1L))
  roc_auc <- sum(ordered_cases *
    (controls_below + 0.5 * ordered_controls)) /
    (sum(ordered_cases) * sum(ordered_controls))
  list(
    epi = c(
      risk_exposed = exposed_risk,
      risk_unexposed = unexposed_risk,
      population_risk = population_risk,
      risk_difference = exposed_risk - unexposed_risk,
      risk_ratio = exposed_risk / unexposed_risk,
      odds_ratio = a * d / (b * c0),
      attributable_fraction_exposed =
        (exposed_risk - unexposed_risk) / exposed_risk,
      population_attributable_fraction =
        (population_risk - unexposed_risk) / population_risk),
    diagnostic = diagnostic, roc_auc = roc_auc,
    causal = c(
      risk_treated = exposed_risk,
      risk_control = unexposed_risk,
      risk_difference = exposed_risk - unexposed_risk,
      risk_ratio = exposed_risk / unexposed_risk,
      odds_ratio = a * d / (b * c0)),
    direct = sum(weights * c(unexposed_risk, exposed_risk)),
    indirect = sum(table[, "event"]) /
      sum(expected_rates * rowSums(table)))
}

.dv_extract_named_regions <- function(regions, names) {
  lower <- vapply(regions[names], `[[`, numeric(1L), "lower")
  upper <- vapply(regions[names], `[[`, numeric(1L), "upper")
  list(lower = unname(lower), upper = unname(upper))
}

.dv_run_epidemiology <- function(replicates) {
  seed <- 20260803L
  set.seed(seed)
  exact_table <- matrix(
    c(80, 320, 160, 240), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("event", "nonevent")))
  weights <- c(unexposed = 0.4, exposed = 0.6)
  expected_rates <- c(unexposed = 0.15, exposed = 0.30)
  contract <- .dv_accuracy_contract(
    "epidemiology_2x2", epsilon = 3,
    natural_l1_sensitivity = 2, coordinate_count = 4L,
    maximum_error = 2000)
  oracle <- .dv_epi_oracle(exact_table, weights, expected_rates)

  no_noise_release <- .dv_table_release(exact_table, contract)
  no_noise <- list(
    epi = .dv_get("ds.vertDPEpi2x2")(
      no_noise_release, exposed = "exposed", event = "event"),
    diagnostic = .dv_get("ds.vertDPDiagnostic2x2")(
      no_noise_release, disease_positive = "exposed",
      test_positive = "event"),
    roc = .dv_get("ds.vertDPROC")(
      no_noise_release, disease_positive = "exposed",
      score_order = c("nonevent", "event")),
    direct = .dv_get("ds.vertDPDirectStandardization")(
      no_noise_release, weights, event = "event"),
    causal = .dv_get("ds.vertDPCausalStandardization")(
      no_noise_release,
      strata = c("all", "all"),
      treatment = c("unexposed", "exposed"),
      treated = "exposed", standard_weights = c(all = 1),
      event = "event"),
    indirect = .dv_get("ds.vertDPIndirectStandardization")(
      no_noise_release, expected_rates, event = "event"))
  oracle_error <- max(
    abs(unlist(no_noise$epi$point_estimates) - oracle$epi),
    abs(no_noise$diagnostic$estimates - oracle$diagnostic),
    abs(no_noise$roc$auc - oracle$roc_auc),
    abs(unlist(no_noise$causal$point_estimates) - oracle$causal),
    abs(no_noise$direct$estimate - oracle$direct),
    abs(no_noise$indirect$estimate - oracle$indirect))

  records <- vector("list", replicates * 19L)
  position <- 0L
  identity_error <- 0
  started <- proc.time()[["elapsed"]]
  for (iteration in seq_len(replicates)) {
    noisy_table <- matrix(
      .dv_clamp(
        as.numeric(exact_table) + .dv_draw_noise(4L, contract), 2000),
      nrow = 2L, dimnames = dimnames(exact_table))
    coordinate_event <- all(
      abs(noisy_table - exact_table) <=
        contract$radius + .dv_numeric_tolerance)
    release <- .dv_table_release(noisy_table, contract)
    epi <- .dv_get("ds.vertDPEpi2x2")(
      release, exposed = "exposed", event = "event")
    diagnostic <- .dv_get("ds.vertDPDiagnostic2x2")(
      release, disease_positive = "exposed", test_positive = "event")
    roc <- .dv_get("ds.vertDPROC")(
      release, disease_positive = "exposed",
      score_order = c("nonevent", "event"))
    direct <- .dv_get("ds.vertDPDirectStandardization")(
      release, weights, event = "event")
    causal <- .dv_get("ds.vertDPCausalStandardization")(
      release,
      strata = c("all", "all"),
      treatment = c("unexposed", "exposed"),
      treated = "exposed", standard_weights = c(all = 1),
      event = "event")
    indirect <- .dv_get("ds.vertDPIndirectStandardization")(
      release, expected_rates, event = "event")

    identity_error <- max(
      identity_error,
      abs(epi$point_estimates$risk_exposed -
            diagnostic$estimates[["sensitivity"]]),
      abs(1 - epi$point_estimates$risk_unexposed -
            diagnostic$estimates[["specificity"]]),
      abs(epi$point_estimates$odds_ratio -
            diagnostic$estimates[["diagnostic_odds_ratio"]]),
      abs(epi$point_estimates$attributable_fraction_exposed -
            (1 - 1 / epi$point_estimates$risk_ratio)),
      abs(epi$point_estimates$population_attributable_fraction -
            (1 - epi$point_estimates$risk_unexposed /
              epi$point_estimates$population_risk)),
      abs(roc$auc - 0.5 * (
        diagnostic$estimates[["sensitivity"]] +
          diagnostic$estimates[["specificity"]])),
      abs(diagnostic$estimates[["balanced_accuracy"]] - 0.5 * (
        diagnostic$estimates[["sensitivity"]] +
          diagnostic$estimates[["specificity"]])),
      abs(diagnostic$estimates[["f1_score"]] -
        2 * diagnostic$estimates[["ppv"]] *
          diagnostic$estimates[["sensitivity"]] /
          (diagnostic$estimates[["ppv"]] +
             diagnostic$estimates[["sensitivity"]])),
      abs(direct$estimate -
            sum(direct$weights * direct$stratum_estimates)),
      abs(causal$point_estimates$risk_treated -
            epi$point_estimates$risk_exposed),
      abs(causal$point_estimates$risk_control -
            epi$point_estimates$risk_unexposed),
      abs(causal$point_estimates$risk_difference -
            epi$point_estimates$risk_difference),
      abs(causal$point_estimates$risk_ratio -
            epi$point_estimates$risk_ratio),
      abs(causal$point_estimates$odds_ratio -
            epi$point_estimates$odds_ratio))

    epi_names <- names(oracle$epi)
    epi_regions <- .dv_extract_named_regions(
      epi$mechanism_regions, epi_names)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "epidemiology", "ds.vertDPEpi2x2", epi_names, iteration,
      oracle$epi, unlist(epi$point_estimates[epi_names]),
      epi_regions$lower, epi_regions$upper, coordinate_event)

    diagnostic_names <- names(oracle$diagnostic)
    diagnostic_regions <- .dv_extract_named_regions(
      diagnostic$mechanism_regions, diagnostic_names)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "epidemiology", "ds.vertDPDiagnostic2x2", diagnostic_names,
      iteration, oracle$diagnostic,
      diagnostic$estimates[diagnostic_names],
      diagnostic_regions$lower, diagnostic_regions$upper,
      coordinate_event)

    position <- position + 1L
    records[[position]] <- .dv_record(
      "epidemiology", "ds.vertDPROC", "auc", iteration,
      oracle$roc_auc, roc$auc,
      roc$auc_mechanism_region[["lower"]],
      roc$auc_mechanism_region[["upper"]], coordinate_event)

    position <- position + 1L
    records[[position]] <- .dv_record(
      "epidemiology", "ds.vertDPDirectStandardization",
      "direct_standardized_risk", iteration, oracle$direct,
      direct$estimate, direct$mechanism_region[["lower"]],
      direct$mechanism_region[["upper"]], coordinate_event)
    causal_names <- names(oracle$causal)
    causal_regions <- .dv_extract_named_regions(
      causal$mechanism_regions, causal_names)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "epidemiology", "ds.vertDPCausalStandardization",
      causal_names, iteration, oracle$causal,
      unlist(causal$point_estimates[causal_names]),
      causal_regions$lower, causal_regions$upper, coordinate_event)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "epidemiology", "ds.vertDPIndirectStandardization",
      "observed_expected_ratio", iteration, oracle$indirect,
      indirect$estimate, indirect$mechanism_region[["lower"]],
      indirect$mechanism_region[["upper"]], coordinate_event)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    records = do.call(rbind, records[seq_len(position)]),
    contract = contract, elapsed = elapsed, seed = seed,
    oracle_error = oracle_error, identity_error = identity_error)
}

.dv_mh_reference <- function(table) {
  oriented <- array(
    0, dim = c(2L, 2L, nrow(table)),
    dimnames = list(
      exposure = c("exposed", "unexposed"),
      outcome = c("event", "no_event"),
      stratum = rownames(table)))
  oriented[1L, 1L, ] <- table[, "exposed_event"]
  oriented[1L, 2L, ] <- table[, "exposed_nonevent"]
  oriented[2L, 1L, ] <- table[, "unexposed_event"]
  oriented[2L, 2L, ] <- table[, "unexposed_nonevent"]
  unname(stats::mantelhaen.test(
    oriented, correct = FALSE, exact = FALSE)$estimate)
}

.dv_run_mantel_haenszel <- function(replicates) {
  seed <- 20260807L
  set.seed(seed)
  roles <- c(
    "exposed_event", "exposed_nonevent",
    "unexposed_event", "unexposed_nonevent")
  exact_table <- matrix(
    c(
      120, 80, 50, 150,
      200, 100, 100, 300,
      80, 120, 40, 160),
    nrow = 3L, byrow = TRUE,
    dimnames = list(c("young", "middle", "old"), roles))
  contract <- .dv_accuracy_contract(
    "mantel_haenszel", epsilon = 3,
    natural_l1_sensitivity = 2, coordinate_count = length(exact_table),
    maximum_error = 2000)
  oracle <- .dv_mh_reference(exact_table)
  no_noise <- .dv_get("ds.vertDPMantelHaenszel")(
    .dv_table_release(exact_table, contract))
  oracle_error <- abs(no_noise$estimate - oracle)
  zero_cost_contract <-
    identical(no_noise$additional_server_calls, 0L) &&
    identical(no_noise$additional_privacy_cost,
              c(epsilon = 0, delta = 0)) &&
    !any(c("statistic", "p_value") %in% names(no_noise))

  records <- vector("list", replicates)
  started <- proc.time()[["elapsed"]]
  for (iteration in seq_len(replicates)) {
    noisy_table <- matrix(
      .dv_clamp(
        as.numeric(exact_table) +
          .dv_draw_noise(length(exact_table), contract),
        2000),
      nrow = nrow(exact_table), dimnames = dimnames(exact_table))
    coordinate_event <- all(
      abs(noisy_table - exact_table) <=
        contract$radius + .dv_numeric_tolerance)
    fit <- .dv_get("ds.vertDPMantelHaenszel")(
      .dv_table_release(noisy_table, contract))
    records[[iteration]] <- .dv_record(
      "mantel_haenszel", "ds.vertDPMantelHaenszel",
      "common_odds_ratio", iteration, oracle, fit$estimate,
      fit$mechanism_region[["lower"]],
      fit$mechanism_region[["upper"]], coordinate_event)
    zero_cost_contract <- zero_cost_contract &&
      identical(fit$additional_server_calls, 0L) &&
      identical(fit$additional_privacy_cost,
                c(epsilon = 0, delta = 0)) &&
      !any(c("statistic", "p_value") %in% names(fit))
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    records = do.call(rbind, records), contract = contract,
    elapsed = elapsed, seed = seed, oracle_error = oracle_error,
    identity_error = oracle_error,
    zero_cost_contract = zero_cost_contract)
}

.dv_survival_oracle <- function(time_grid, exit_counts, causes) {
  events <- exit_counts[, -1L, drop = FALSE]
  at_risk <- rev(cumsum(rev(rowSums(exit_counts))))
  cause_hazard <- events / at_risk
  all_hazard <- rowSums(cause_hazard)
  survival <- cumprod(1 - all_hazard)
  nelson_aalen <- cumsum(all_hazard)
  cif <- matrix(0, nrow(exit_counts), length(causes),
                dimnames = list(NULL, causes))
  previous_survival <- 1
  previous_cif <- numeric(length(causes))
  for (index in seq_along(time_grid)) {
    previous_cif <- previous_cif +
      previous_survival * cause_hazard[index, ]
    previous_survival <- previous_survival * (1 - all_hazard[[index]])
    cif[index, ] <- previous_cif
  }
  list(
    at_risk = at_risk, survival = survival,
    nelson_aalen = nelson_aalen, cumulative_incidence = cif)
}

.dv_survival_release <- function(
    histogram, time_grid, causes, capacity, contract) {
  result <- list(
    released = TRUE, analysis_id = "validation_survival",
    analysis_version = "v1", time_grid = time_grid,
    time_lower_bound = 0, time_upper_bound = tail(time_grid, 1L),
    interval_semantics =
      "(previous_endpoint,current_endpoint] after public-bound clipping",
    unit_collapse = "first_event_else_latest_censor_public_tiebreak",
    censor_level = "0", causes = causes, delayed_entry = FALSE,
    histogram = histogram, coordinate_count = as.integer(length(histogram)),
    histogram_layout =
      "exit[T x (censor+causes) column-major],not_in_analysis",
    not_in_analysis_definition =
      "DP-noisy bin for unknown outcome or non-finite time",
    invalid_unit_rule = "invalid_patient_ids_rejected_by_admission",
    mechanism = .dv_laplace_mechanism,
    implementation = "two pinned peers; Ring128 exact signed finalizer",
    sampler = .dv_laplace_sampler,
    randomness = "independent HMAC/HKDF/ChaCha20 streams at both peers",
    l1_sensitivity = 1L, l2_sensitivity = 1,
    global_l1_sensitivity = contract$natural_l1_sensitivity,
    global_l2_sensitivity = sqrt(2),
    sensitivity_scope = paste(
      "l1/l2 are this signed survival block; global_l1/global_l2",
      "calibrate the one complete capsule vector"),
    max_histogram_cells_per_unit = 1L,
    contribution_unit = "add_remove_patient",
    postprocessing = "fixed public per-coordinate clamp",
    clipped_coordinates = NA_integer_, clipping_observable = FALSE,
    accuracy_95_abs_per_coordinate = .dv_accuracy_contract(
      paste0(contract$family, "_marginal"), contract$epsilon,
      contract$natural_l1_sensitivity, 1L, contract$scale,
      capacity)$radius,
    accuracy_simultaneous_95_abs = contract$radius,
    accuracy_simultaneous_confidence = .dv_level,
    accuracy_simultaneous_method = contract$accuracy_method,
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    privacy_epoch = 1, noise_key_id = "validation-noise-key-a",
    sticky_noise = "one immutable capsule vector; unlimited replay",
    epsilon = contract$epsilon, delta = 2^-100,
    implementation_delta = .dv_implementation_delta,
    adjacency = "add_remove_patient",
    capsule_id = strrep("d", 64L), final_vector_root = strrep("e", 64L),
    coordinate_order_sha256 = strrep("f", 64L), server = "site_a",
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE)
  result <- .dv_get(".dsvert_dp_survival_postprocess")(result)
  class(result) <- c("ds.vertDPSurvival", "list")
  result
}

.dv_fixed_grid_rmst <- function(time_grid, lower_bound, survival, tau) {
  starts <- c(lower_bound, head(time_grid, -1L))
  weights <- pmax(0, pmin(time_grid, tau) - starts)
  sum(weights * c(1, head(survival, -1L)))
}

.dv_run_survival <- function(replicates) {
  seed <- 20260804L
  set.seed(seed)
  time_grid <- as.numeric(1:6)
  causes <- c("A", "B")
  exit_counts <- cbind(
    `0` = c(5, 10, 15, 20, 30, 75),
    A = c(10, 8, 6, 4, 2, 0),
    B = c(5, 4, 3, 2, 1, 0))
  exact_histogram <- c(as.vector(exit_counts), 0)
  capacity <- 250
  contract <- .dv_accuracy_contract(
    "survival_rmst", epsilon = 3,
    natural_l1_sensitivity = 2,
    coordinate_count = length(exact_histogram),
    maximum_error = capacity)
  oracle <- .dv_survival_oracle(time_grid, exit_counts, causes)
  oracle_rmst <- .dv_fixed_grid_rmst(
    time_grid, 0, oracle$survival, tau = 6)
  quantile_probability <- 0.01
  oracle_quantile <- .dv_get(
    ".dsvert_dp_survival_quantile_endpoint")(
      time_grid, oracle$survival, quantile_probability)

  no_noise <- .dv_survival_release(
    exact_histogram, time_grid, causes, capacity, contract)
  no_noise_rmst <- .dv_get("ds.vertDPRMST")(no_noise, 6)
  no_noise_quantile <- .dv_get("ds.vertDPSurvivalQuantile")(
    no_noise, quantile_probability)
  oracle_error <- max(
    abs(no_noise$curve$kaplan_meier - oracle$survival),
    abs(no_noise$curve$nelson_aalen - oracle$nelson_aalen),
    abs(no_noise$cumulative_incidence - oracle$cumulative_incidence),
    abs(no_noise_rmst$rmst - oracle_rmst),
    abs(no_noise_quantile$quantile - oracle_quantile))

  records <- vector("list", replicates * 5L)
  position <- 0L
  identity_error <- 0
  started <- proc.time()[["elapsed"]]
  for (iteration in seq_len(replicates)) {
    noisy <- .dv_clamp(
      exact_histogram + .dv_draw_noise(length(exact_histogram), contract),
      capacity)
    coordinate_event <- all(
      abs(noisy - exact_histogram) <=
        contract$radius + .dv_numeric_tolerance)
    release <- .dv_survival_release(
      noisy, time_grid, causes, capacity, contract)
    rmst <- .dv_get("ds.vertDPRMST")(release, 6)
    survival_quantile <- .dv_get("ds.vertDPSurvivalQuantile")(
      release, quantile_probability)
    identity_error <- max(
      identity_error,
      abs(release$curve$kaplan_meier +
            rowSums(release$cumulative_incidence) - 1),
      pmax(0, diff(release$curve$kaplan_meier)),
      pmax(0, -diff(release$curve$nelson_aalen)),
      pmax(0, -apply(
        release$cumulative_incidence, 2L, diff)))

    position <- position + 1L
    records[[position]] <- .dv_record(
      "survival", "ds.vertDPSurvival", "kaplan_meier_fixed_grid",
      iteration, oracle$survival, release$curve$kaplan_meier,
      release$curve$kaplan_meier_mechanism_lower_95,
      release$curve$kaplan_meier_mechanism_upper_95,
      coordinate_event)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "survival", "ds.vertDPSurvival", "nelson_aalen_fixed_grid",
      iteration, oracle$nelson_aalen, release$curve$nelson_aalen,
      release$curve$nelson_aalen_mechanism_lower_95,
      release$curve$nelson_aalen_mechanism_upper_95,
      coordinate_event)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "survival", "ds.vertDPSurvival", "cumulative_incidence_fixed_grid",
      iteration, as.numeric(oracle$cumulative_incidence),
      as.numeric(release$cumulative_incidence),
      as.numeric(release$cumulative_incidence_mechanism_lower_95),
      as.numeric(release$cumulative_incidence_mechanism_upper_95),
      coordinate_event)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "survival", "ds.vertDPRMST", "rmst_tau_6_fixed_grid",
      iteration, oracle_rmst, rmst$rmst,
      rmst$rmst_mechanism_lower_95,
      rmst$rmst_mechanism_upper_95, coordinate_event)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "survival", "ds.vertDPSurvivalQuantile",
      "event_probability_0.01_first_grid_endpoint",
      iteration, oracle_quantile, survival_quantile$quantile,
      survival_quantile$quantile_mechanism_lower_95,
      survival_quantile$quantile_mechanism_upper_95,
      coordinate_event)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    records = do.call(rbind, records[seq_len(position)]),
    contract = contract, elapsed = elapsed, seed = seed,
    oracle_error = oracle_error, identity_error = identity_error)
}

.dv_gaussian_fixture <- function() {
  grid <- seq(0, 1, length.out = 30L)
  x1 <- rep(grid, each = length(grid))
  x2 <- rep(grid, times = length(grid))
  y <- 0.1 + 0.3 * x1 + 0.4 * x2
  design <- cbind(`(Intercept)` = 1, x1 = x1, x2 = x2)
  scale <- 256
  capacity <- 1000
  artifact <- list(
    design_terms = colnames(design), intercept = TRUE,
    predictor_order = c("x1", "x2"),
    outcome = list(column = "y", lower = 0, upper = 1),
    predictors = list(
      x1 = list(column = "x1", lower = 0, upper = 1),
      x2 = list(column = "x2", lower = 0, upper = 1)))
  coordinates <- length(y)
  for (right in seq_len(ncol(design))) {
    for (left in seq_len(right)) {
      coordinates <- c(coordinates, sum(round(
        design[, left] * design[, right] * scale)) / scale)
    }
  }
  coordinates <- c(
    coordinates,
    vapply(seq_len(ncol(design)), function(column) {
      sum(round(design[, column] * y * scale)) / scale
    }, numeric(1L)),
    sum(round(y^2 * scale)) / scale)
  contract <- .dv_accuracy_contract(
    "gaussian_correlation_pca", epsilon = 6,
    natural_l1_sensitivity = length(coordinates) + 1L,
    coordinate_count = length(coordinates) + 1L,
    scale = scale, maximum_error = capacity)
  central <- stats::lm.fit(design, y)$coefficients
  names(central) <- colnames(design)
  oracle_correlation <- stats::cor(cbind(x1 = x1, x2 = x2, y = y))
  oracle_eigen <- eigen(
    oracle_correlation, symmetric = TRUE, only.values = TRUE)$values
  list(
    artifact = artifact, coordinates = coordinates, capacity = capacity,
    scale = scale, contract = contract, coefficients = central,
    correlation = oracle_correlation, eigenvalues = oracle_eigen)
}

.dv_gaussian_fit <- function(coordinates, fixture) {
  moment <- .dv_get(".dsvert_dp_gaussian_unpack")(
    coordinates, fixture$artifact, fixture$capacity)
  fit <- .dv_get(".dsvert_dp_gaussian_solve")(
    moment, fixture$artifact, ridge = 0)
  coefficients <- .dv_get(".dsvert_dp_gaussian_original_coefficients")(
    fit$coefficients, fixture$artifact)
  list(coefficients = coefficients, moment = moment, fit = fit)
}

.dv_run_gaussian <- function(replicates) {
  seed <- 20260805L
  set.seed(seed)
  fixture <- .dv_gaussian_fixture()
  no_noise <- .dv_gaussian_fit(fixture$coordinates, fixture)
  oracle_error <- max(abs(
    no_noise$coefficients - fixture$coefficients))
  records <- vector("list", replicates)
  identity_error <- 0
  started <- proc.time()[["elapsed"]]
  for (iteration in seq_len(replicates)) {
    noisy <- .dv_clamp(
      fixture$coordinates + .dv_draw_noise(
        length(fixture$coordinates), fixture$contract),
      fixture$capacity)
    coordinate_event <- all(
      abs(noisy - fixture$coordinates) <=
        fixture$contract$radius + .dv_numeric_tolerance)
    fit <- .dv_gaussian_fit(noisy, fixture)
    identity_error <- max(
      identity_error,
      abs(fit$fit$residual_second_moment -
            max(0, fit$moment$outcome_square -
              2 * crossprod(fit$fit$coefficients, fit$moment$cross) +
              crossprod(
                fit$fit$coefficients,
                fit$moment$gram %*% fit$fit$coefficients))))
    records[[iteration]] <- .dv_record(
      "gaussian", "ds.vertDPGaussian", names(fixture$coefficients),
      iteration, fixture$coefficients, fit$coefficients,
      NA_real_, NA_real_, coordinate_event,
      mechanism_region_available = FALSE)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    records = do.call(rbind, records), contract = fixture$contract,
    elapsed = elapsed, seed = seed, oracle_error = oracle_error,
    identity_error = identity_error)
}

.dv_spectral_release <- function(coordinates, fixture, radius) {
  artifact <- fixture$artifact
  variables <- c(artifact$predictor_order, artifact$outcome$column)
  released <- .dv_get(".dsvert_dp_cor_gaussian_coordinates")(
    coordinates, artifact)
  moment <- .dv_get(".dsvert_dp_gaussian_unpack")(
    coordinates, artifact, fixture$capacity)
  augmented_names <- c(artifact$design_terms, artifact$outcome$column)
  augmented <- moment$augmented_projected
  dimnames(augmented) <- list(augmented_names, augmented_names)
  selected <- match(variables, augmented_names)
  intercept <- match("(Intercept)", augmented_names)
  mass <- augmented[intercept, intercept]
  sums <- augmented[intercept, selected]
  second <- augmented[selected, selected, drop = FALSE]
  centered <- (second - outer(sums, sums) / mass)
  centered <- (centered + t(centered)) / 2
  raw <- centered / sqrt(outer(diag(centered), diag(centered)))
  raw <- (raw + t(raw)) / 2
  diag(raw) <- 1
  dimnames(raw) <- list(variables, variables)
  psd <- .dv_get(".dsvert_dp_cor_psd")(raw)
  lower <- upper <- raw
  for (left_index in seq_len(length(variables) - 1L)) {
    for (right_index in seq.int(left_index + 1L, length(variables))) {
      left <- variables[[left_index]]
      right <- variables[[right_index]]
      interval <- .dv_get(".dsvert_dp_cor_interval")(c(
        released$mass, released$sums[[left]], released$sums[[right]],
        released$second[left, left], released$second[right, right],
        released$second[left, right]), radius, fixture$capacity,
        fixture$scale)
      lower[left_index, right_index] <- lower[right_index, left_index] <-
        interval$correlation[["lower"]]
      upper[left_index, right_index] <- upper[right_index, left_index] <-
        interval$correlation[["upper"]]
    }
  }
  projected_lower <- projected_upper <- psd$matrix
  for (row in seq_along(variables)) {
    for (column in seq_along(variables)) {
      enclosure_radius <- max(
        abs(psd$matrix[row, column] - lower[row, column]),
        abs(psd$matrix[row, column] - upper[row, column]))
      projected_lower[row, column] <- max(
        -1, psd$matrix[row, column] - enclosure_radius)
      projected_upper[row, column] <- min(
        1, psd$matrix[row, column] + enclosure_radius)
    }
  }
  complete_case_n <- matrix(
    mass, length(variables), length(variables),
    dimnames = list(variables, variables))
  correlation <- structure(list(
    released = TRUE, source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE, legacy_exact_route_called = FALSE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    psd_projection_applied = TRUE, source_artifact_family = "gaussian_models",
    estimand_missingness = "complete_case_joint", pca_eligible = TRUE,
    correlation_raw_pairwise = NULL,
    correlation_raw_complete_case = raw,
    complete_case_n = complete_case_n,
    capsule_id = strrep("a", 64L), analysis_id = "validation_spectral",
    correlation = psd$matrix, var_names = variables, n_obs = mass,
    psd_projection = psd,
    correlation_95_enclosure_raw_estimand_around_projected_release = list(
      lower = projected_lower, upper = projected_upper),
    epsilon = fixture$contract$epsilon, delta = 2^-100,
    mechanism = .dv_laplace_mechanism,
    cross_owner_state = "same_owner_materialized"),
    class = c("ds.vertDPCor", "ds.cor", "list"))
  list(correlation = correlation, lower = lower, upper = upper)
}

.dv_loading_sine <- function(estimate, truth) {
  estimate <- estimate / sqrt(sum(estimate^2))
  truth <- truth / sqrt(sum(truth^2))
  sqrt(max(0, 1 - min(1, abs(sum(estimate * truth)))^2))
}

.dv_run_spectral <- function(replicates) {
  seed <- 20260806L
  set.seed(seed)
  fixture <- .dv_gaussian_fixture()
  oracle_decomposition <- eigen(fixture$correlation, symmetric = TRUE)
  no_noise <- .dv_spectral_release(
    fixture$coordinates, fixture, radius = 0)
  no_noise_pca <- .dv_get("ds.vertPCA")(
    cor_result = no_noise$correlation, verbose = FALSE)
  oracle_error <- max(
    abs(no_noise$correlation$correlation_raw_complete_case -
          fixture$correlation),
    abs(no_noise_pca$eigenvalues - fixture$eigenvalues),
    .dv_loading_sine(
      no_noise_pca$loadings[, 1L], oracle_decomposition$vectors[, 1L]))
  records <- vector("list", replicates * 3L)
  identity_error <- 0
  position <- 0L
  pairs <- which(upper.tri(fixture$correlation), arr.ind = TRUE)
  pair_names <- paste0(
    rownames(fixture$correlation)[pairs[, 1L]], "_",
    colnames(fixture$correlation)[pairs[, 2L]])
  started <- proc.time()[["elapsed"]]
  for (iteration in seq_len(replicates)) {
    noisy <- .dv_clamp(
      fixture$coordinates + .dv_draw_noise(
        length(fixture$coordinates), fixture$contract),
      fixture$capacity)
    coordinate_event <- all(
      abs(noisy - fixture$coordinates) <=
        fixture$contract$radius + .dv_numeric_tolerance)
    released <- .dv_spectral_release(
      noisy, fixture, fixture$contract$radius)
    pca <- .dv_get("ds.vertPCA")(
      cor_result = released$correlation, verbose = FALSE)
    identity_error <- max(
      identity_error,
      abs(sum(pca$eigenvalues) - nrow(fixture$correlation)),
      pmax(0, -eigen(
        released$correlation$correlation,
        symmetric = TRUE, only.values = TRUE)$values))
    position <- position + 1L
    records[[position]] <- .dv_record(
      "correlation", "ds.vertCor", pair_names, iteration,
      fixture$correlation[pairs],
      released$correlation$correlation_raw_complete_case[pairs],
      released$lower[pairs], released$upper[pairs], coordinate_event)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "pca", "ds.vertPCA", paste0("eigenvalue_PC", 1:3), iteration,
      fixture$eigenvalues, pca$eigenvalues,
      pca$eigenvalue_95_mechanism_regions[, "lower"],
      pca$eigenvalue_95_mechanism_regions[, "upper"], coordinate_event)
    position <- position + 1L
    records[[position]] <- .dv_record(
      "pca", "ds.vertPCA", "PC1_sine_angle", iteration,
      0, .dv_loading_sine(
        pca$loadings[, 1L], oracle_decomposition$vectors[, 1L]),
      NA_real_, NA_real_, coordinate_event,
      mechanism_region_available = FALSE)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    records = do.call(rbind, records[seq_len(position)]),
    contract = fixture$contract, elapsed = elapsed, seed = seed,
    oracle_error = oracle_error, identity_error = identity_error)
}

.dv_edge_cases <- function(describe_contract, epi_contract,
                           survival_contract) {
  empty_describe <- .dv_describe_release(
    numeric(length = 12L), 0, 4, seq(0.5, 4, 0.5), 600,
    describe_contract, c(0, 0.5, 1))
  empty_quantile <- .dv_get("ds.vertDPQuantile")(
    empty_describe, c(0, 0.5, 1))
  empty_median <- .dv_get("ds.vertDPMedian")(empty_describe)
  tied_describe <- .dv_describe_release(
    c(400, 200, 120, rep(50, 8), 0), 0, 4,
    seq(0.5, 4, 0.5), 600, describe_contract,
    c(0, 0.25, 0.5, 0.75, 1))
  tied_quantile <- .dv_get("ds.vertDPQuantile")(
    tied_describe, c(0, 0.25, 0.5, 0.75, 1))

  empty_table <- matrix(
    0, 2L, 2L, dimnames = list(c("unexposed", "exposed"),
                               c("event", "nonevent")))
  empty_release <- .dv_table_release(empty_table, epi_contract)
  empty_epi <- .dv_get("ds.vertDPEpi2x2")(
    empty_release, "exposed", "event")
  empty_diagnostic <- .dv_get("ds.vertDPDiagnostic2x2")(
    empty_release, "exposed", "event")
  empty_direct <- .dv_get("ds.vertDPDirectStandardization")(
    empty_release, c(0.5, 0.5), "event")
  empty_causal <- .dv_get("ds.vertDPCausalStandardization")(
    empty_release,
    strata = c("all", "all"),
    treatment = c("unexposed", "exposed"),
    treated = "exposed", standard_weights = c(all = 1),
    event = "event")
  empty_indirect <- .dv_get("ds.vertDPIndirectStandardization")(
    empty_release, c(0.1, 0.2), "event")

  perfect <- matrix(
    c(50, 0, 0, 50), 2L, byrow = TRUE,
    dimnames = list(c("disease_positive", "disease_negative"),
                    c("test_positive", "test_negative")))
  perfect_diagnostic <- .dv_get("ds.vertDPDiagnostic2x2")(
    .dv_table_release(perfect, epi_contract),
    "disease_positive", "test_positive")

  zero_rate <- matrix(
    c(10, 0, 0, 0), 2L, byrow = TRUE,
    dimnames = list(c("zero_rate", "positive_rate"),
                    c("event", "nonevent")))
  infinite_indirect <- .dv_get("ds.vertDPIndirectStandardization")(
    .dv_table_release(zero_rate, epi_contract),
    c(zero_rate = 0, positive_rate = 0.2), "event")

  empty_survival <- .dv_survival_release(
    numeric(19L), as.numeric(1:6), c("A", "B"), 250,
    survival_contract)
  empty_rmst <- .dv_get("ds.vertDPRMST")(empty_survival, 6)
  empty_survival_quantile <- .dv_get("ds.vertDPSurvivalQuantile")(
    empty_survival, 0.5)
  empty_median_survival <- .dv_get("ds.vertDPMedianSurvival")(
    empty_survival)

  data.frame(
    case = c(
      "describe_empty_histogram", "quantile_probabilities_0_and_1",
      "quantile_nonempty_endpoint_probabilities",
      "quantile_exact_cdf_tie_first_crossing",
      "median_empty_histogram", "epi_empty_groups",
      "diagnostic_empty_table", "direct_empty_strata",
      "causal_empty_stratum_arms", "indirect_zero_expected_denominator",
      "diagnostic_perfect_boundary",
      "indirect_infinite_boundary", "survival_empty_curve",
      "rmst_empty_curve", "survival_quantile_beyond_grid",
      "median_survival_beyond_grid"),
    pass = c(
      identical(empty_describe$status,
                "fixed_public_clamp_applied_preclamp_state_not_released"),
      all(empty_quantile$status == "dp_projected_histogram_empty") &&
        identical(empty_quantile$probability, c(0, 0.5, 1)),
      identical(tied_quantile$bin_index[c(1L, 5L)], c(1L, 8L)),
      identical(tied_quantile$bin_index[c(2L, 3L, 4L)], c(2L, 4L, 6L)),
      identical(empty_median$status, "dp_projected_histogram_empty"),
      identical(empty_epi$status, "dp_point_non_estimable"),
      identical(empty_diagnostic$status, "non_estimable"),
      identical(empty_direct$status, "dp_point_non_estimable"),
      identical(empty_causal$status, "dp_point_non_estimable"),
      identical(empty_indirect$status, "dp_point_non_estimable"),
      perfect_diagnostic$status %in% c("boundary", "ok"),
      identical(infinite_indirect$status, "boundary_infinite"),
      identical(empty_survival$status,
                "dp_curve_empty_after_postprocessing"),
      identical(empty_rmst$rmst, 6),
      identical(
        empty_survival_quantile$status,
        "mechanism_region_extends_beyond_grid") &&
        identical(empty_survival_quantile$point_status,
                  "not_reached_by_grid_end") &&
        is.na(empty_survival_quantile$quantile) &&
        is.infinite(empty_survival_quantile$quantile_mechanism_upper_95),
      identical(
        empty_median_survival$status,
        "mechanism_region_extends_beyond_grid") &&
        identical(empty_median_survival$point_status,
                  "not_reached_by_grid_end") &&
        is.na(empty_median_survival$quantile)),
    detail = c(
      empty_describe$status,
      paste(empty_quantile$status, collapse = ","),
      paste0("bin_index=", paste(
        tied_quantile$bin_index[c(1L, 5L)], collapse = ",")),
      paste0("bin_index=", paste(
        tied_quantile$bin_index[c(2L, 3L, 4L)], collapse = ",")),
      empty_median$status,
      empty_epi$status,
      empty_diagnostic$status,
      empty_direct$status,
      empty_causal$status,
      empty_indirect$status,
      perfect_diagnostic$status,
      infinite_indirect$status,
      empty_survival$status,
      paste0("rmst=", format(empty_rmst$rmst)),
      empty_survival_quantile$status,
      empty_median_survival$status),
    stringsAsFactors = FALSE)
}

.dv_summarize <- function(records, elapsed_by_family, replicates) {
  records$covered <- NA
  region <- records$mechanism_region_available
  records$covered[region] <- .dv_inside(
    records$truth[region], records$lower[region], records$upper[region])
  groups <- split(
    seq_len(nrow(records)),
    interaction(records$family, records$method, records$estimand,
                drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(index) {
    value <- records[index, , drop = FALSE]
    available <- value$mechanism_region_available
    width <- value$upper[available] - value$lower[available]
    conditional <- value$covered[available & value$coordinate_event]
    data.frame(
      family = value$family[[1L]], method = value$method[[1L]],
      estimand = value$estimand[[1L]], observations = nrow(value),
      bias = mean(value$estimate - value$truth),
      rmse = sqrt(mean((value$estimate - value$truth)^2)),
      mechanism_region_available = any(available),
      mechanism_coverage = if (any(available)) {
        mean(value$covered[available])
      } else NA_real_,
      conditional_coordinate_coverage = if (length(conditional)) {
        mean(conditional)
      } else NA_real_,
      mean_region_width = if (length(width)) mean(width) else NA_real_,
      coordinate_event_rate = mean(value$coordinate_event),
      estimate_finite_rate = mean(is.finite(value$estimate)),
      finite_rate = if (any(available)) {
        mean(is.finite(value$estimate[available]) &
               is.finite(value$lower[available]) &
               is.finite(value$upper[available]))
      } else mean(is.finite(value$estimate)),
      elapsed_seconds_family = elapsed_by_family[[value$family[[1L]]]],
      replicates = as.integer(replicates),
      stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  result <- result[order(result$family, result$method, result$estimand), ]
  rownames(result) <- NULL
  result
}

.dv_contract_frame <- function(contracts) {
  do.call(rbind, lapply(contracts, function(value) data.frame(
    family = value$family, epsilon = value$epsilon,
    natural_l1_sensitivity = value$natural_l1_sensitivity,
    coordinate_count = value$coordinate_count, scale = value$scale,
    simultaneous_radius_95 = value$radius,
    implementation_tv_upper_bound =
      value$implementation_tv_upper_bound,
    certified_failure_probability_upper =
      value$certified_failure_probability_upper,
    certified_coverage_lower = value$certified_coverage_lower,
    accuracy_method = value$accuracy_method,
    stringsAsFactors = FALSE)))
}

.dv_postprocessors_have_no_dsi <- function() {
  methods <- c(
    "ds.vertDPQuantile", "ds.vertDPMedian", "ds.vertDPEpi2x2",
    "ds.vertDPEpi2x2Inference",
    "ds.vertDPMantelHaenszel",
    "ds.vertDPDiagnostic2x2", "ds.vertDPDiagnostic2x2Inference",
    "ds.vertDPROC",
    "ds.vertDPDirectStandardization",
    "ds.vertDPDirectStandardizationInference",
    "ds.vertDPCausalStandardization",
    "ds.vertDPCausalStandardizationInference",
    "ds.vertDPIndirectStandardization",
    "ds.vertDPIndirectStandardizationInference", "ds.vertDPRMST",
    "ds.vertDPRMTL", "ds.vertDPSurvivalQuantile",
    "ds.vertDPMedianSurvival")
  !any(vapply(methods, function(method) {
    text <- paste(deparse(body(.dv_get(method))), collapse = "\n")
    grepl("DSI::|datashield[.]", text, perl = TRUE)
  }, logical(1L)))
}

dsvert_run_dp_statistical_validation <- function(replicates = 1000L) {
  if (!is.numeric(replicates) || length(replicates) != 1L ||
      is.na(replicates) || !is.finite(replicates) || replicates < 1 ||
      replicates > 100000L || replicates != floor(replicates)) {
    stop("replicates must be one integer in [1, 100000]", call. = FALSE)
  }
  replicates <- as.integer(replicates)
  old_kind <- RNGkind()
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                            inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")

  describe <- .dv_run_describe(replicates)
  epidemiology <- .dv_run_epidemiology(replicates)
  mantel_haenszel <- .dv_run_mantel_haenszel(replicates)
  survival <- .dv_run_survival(replicates)
  gaussian <- .dv_run_gaussian(replicates)
  spectral <- .dv_run_spectral(replicates)
  suites <- list(
    describe, epidemiology, mantel_haenszel, survival, gaussian, spectral)
  records <- do.call(rbind, lapply(suites, `[[`, "records"))
  elapsed <- c(
    describe = describe$elapsed,
    epidemiology = epidemiology$elapsed,
    mantel_haenszel = mantel_haenszel$elapsed,
    survival = survival$elapsed,
    gaussian = gaussian$elapsed,
    correlation = spectral$elapsed,
    pca = spectral$elapsed)
  summary <- .dv_summarize(records, elapsed, replicates)
  contracts <- .dv_contract_frame(list(
    describe$contract, epidemiology$contract, mantel_haenszel$contract,
    survival$contract,
    gaussian$contract))
  edges <- .dv_edge_cases(
    describe$contract, epidemiology$contract, survival$contract)

  records$covered <- NA
  region <- records$mechanism_region_available
  records$covered[region] <- .dv_inside(
    records$truth[region], records$lower[region], records$upper[region])
  conditional <- records$covered[region & records$coordinate_event]
  region_width <- records$upper[region] - records$lower[region]
  gates <- data.frame(
    gate = c(
      "deployed_accuracy_contract_at_least_95_percent",
      "truth_covered_whenever_certified_coordinate_box_holds",
      "all_main_regions_ordered_and_finite",
      "all_point_estimates_finite",
      "central_oracle_zero_noise_identity",
      "gaussian_quantized_oracle_fidelity",
      "correlation_pca_quantized_oracle_fidelity",
      "describe_quantization_region_contains_bounded_truth",
      "quantile_median_identity",
      "epidemiology_diagnostic_and_causal_standardization_identities",
      "mantel_haenszel_zero_cost_no_classical_dp_inference",
      "survival_competing_risk_and_monotonicity_identities",
      "gaussian_solve_algebra_identity",
      "pca_trace_and_psd_identities",
      "all_edge_cases_typed",
      "release_only_postprocessors_have_no_dsi_calls"),
    pass = c(
      all(contracts$certified_coverage_lower >=
            .dv_level - 256 * .Machine$double.eps),
      length(conditional) > 0L && all(conditional),
      all(is.finite(records$truth[region])) &&
        all(is.finite(records$estimate[region])) &&
        all(is.finite(records$lower[region])) &&
        all(is.finite(records$upper[region])) &&
        all(records$lower[region] <= records$upper[region]),
      all(is.finite(records$truth)) && all(is.finite(records$estimate)),
      max(vapply(
        list(describe, epidemiology, mantel_haenszel, survival),
        `[[`, numeric(1L),
        "oracle_error")) <= 1e-12,
      gaussian$oracle_error <= 1 / 256,
      spectral$oracle_error <= 1 / 256,
      isTRUE(describe$quantization_coverage),
      describe$identity_error <= 1e-12,
      epidemiology$identity_error <= 1e-12,
      mantel_haenszel$identity_error <= 1e-12 &&
        isTRUE(mantel_haenszel$zero_cost_contract),
      survival$identity_error <= 1e-12,
      gaussian$identity_error <= 1e-10,
      spectral$identity_error <= 1e-10,
      all(edges$pass),
      .dv_postprocessors_have_no_dsi()),
    detail = c(
      paste0("minimum=", format(min(contracts$certified_coverage_lower),
                                digits = 12L)),
      paste0("checked=", length(conditional)),
      paste0(
        "rows=", nrow(records), "; width_range=[",
        format(min(region_width), digits = 8L), ",",
        format(max(region_width), digits = 8L), "]"),
      paste0("rows=", nrow(records), "; all finite"),
      paste0("max_abs_error=", format(max(vapply(
        list(describe, epidemiology, mantel_haenszel, survival),
        `[[`, numeric(1L),
        "oracle_error")), digits = 6L)),
      paste0("max_abs_error=", format(gaussian$oracle_error, digits = 6L)),
      paste0("max_abs_error=", format(spectral$oracle_error, digits = 6L)),
      "mean and variance bounded truth inside no-noise quantization region",
      paste0("max_abs_error=", format(describe$identity_error)),
      paste0("max_abs_error=", format(epidemiology$identity_error)),
      paste0(
        "oracle_error=", format(mantel_haenszel$identity_error),
        "; zero_cost_no_classical=",
        mantel_haenszel$zero_cost_contract),
      paste0("max_abs_error=", format(survival$identity_error)),
      paste0("max_abs_error=", format(gaussian$identity_error)),
      paste0("max_abs_error=", format(spectral$identity_error)),
      paste0("passed=", sum(edges$pass), "/", nrow(edges)),
      "static body audit of release-only public postprocessors"),
    stringsAsFactors = FALSE)
  if (!all(gates$pass)) {
    failed <- gates$gate[!gates$pass]
    stop("DP statistical validation failed: ",
         paste(failed, collapse = ", "), call. = FALSE)
  }
  list(
    summary = summary, contracts = contracts, gates = gates,
    edge_cases = edges,
    metadata = list(
      schema_version = "dsvert-dp-statistical-validation-v2",
      execution_scope = paste(
        "fixed finite-dataset distributional simulation using the ideal",
        "two-peer geometric convolution; not DSI E2E and not a productive",
        "cryptographic-sampler replay"),
      uncertainty_scope = paste(
        "DP mechanism noise and public grid/quantization only; no sampling",
        "or population confidence intervals"),
      replicates = replicates,
      seeds = c(
        describe = describe$seed, epidemiology = epidemiology$seed,
        mantel_haenszel = mantel_haenszel$seed,
        survival = survival$seed, gaussian = gaussian$seed,
        spectral = spectral$seed),
      r_version = R.version.string,
      package_version = as.character(
        utils::packageVersion("dsVertClient"))))
}

.dv_markdown_table <- function(value, digits = 6L) {
  printable <- value
  numeric_columns <- vapply(printable, is.numeric, logical(1L))
  printable[numeric_columns] <- lapply(
    printable[numeric_columns], function(column) {
      format(column, digits = digits, scientific = FALSE, trim = TRUE)
    })
  lines <- c(
    paste0("| ", paste(names(printable), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(printable)), collapse = " | "),
           " |"))
  body <- apply(printable, 1L, function(row) {
    paste0("| ", paste(row, collapse = " | "), " |")
  })
  c(lines, body)
}

.dv_write_validation <- function(result, output_dir) {
  if (!dir.exists(output_dir)) {
    stop("output_dir must already exist", call. = FALSE)
  }
  summary_path <- file.path(
    output_dir, "dp_statistical_validation_20260802.csv")
  contracts_path <- file.path(
    output_dir, "dp_statistical_validation_contracts_20260802.csv")
  edges_path <- file.path(
    output_dir, "dp_statistical_validation_edges_20260802.csv")
  report_path <- file.path(
    output_dir, "dp_statistical_validation_20260802.md")
  utils::write.csv(result$summary, summary_path, row.names = FALSE,
                   na = "NA")
  utils::write.csv(result$contracts, contracts_path, row.names = FALSE,
                   na = "NA")
  utils::write.csv(result$edge_cases, edges_path, row.names = FALSE,
                   na = "NA")
  report <- c(
    "# dsVert DP statistical validation — 2026-08-02",
    "",
    "## Scope",
    "",
    paste0("- ", result$metadata$execution_scope, "."),
    paste0("- ", result$metadata$uncertainty_scope, "."),
    paste0("- Fixed mechanism replicates per family: `",
           result$metadata$replicates, "`."),
    paste0("- Seeds: `", paste(
      names(result$metadata$seeds), result$metadata$seeds,
      sep = "=", collapse = "`, `"), "`."),
    "- Runtime is descriptive only and has no pass/fail threshold.",
    "- The central oracle is the bounded finite dataset: pre-quantization",
    "  moments for fidelity, and fixed public-bin endpoints for quantiles.",
    "  Zero-noise identities separately target the deployed quantized",
    "  estimand, and its reported grid region must contain the bounded truth.",
    "  No population or sampling confidence interval is invented.",
    "",
    "## Gates",
    "",
    .dv_markdown_table(result$gates),
    "",
    "## Deployed accuracy contracts",
    "",
    .dv_markdown_table(result$contracts),
    "",
    "The analytic gate combines the exact ideal convolution tail, the union",
    "bound over the released coordinates, and the productive sampler's",
    "published two-peer total-variation allowance. Region coverage is also",
    "required deterministically whenever every exact coordinate lies inside",
    "that certified simultaneous box.",
    "",
    "## Bias, RMSE, mechanism coverage, width, and runtime",
    "",
    .dv_markdown_table(result$summary),
    "",
    "`mechanism_coverage` is conditional on the fixed finite dataset and",
    "covers only DP mechanism noise plus the explicitly reported public-grid",
    "or quantization component. It is not a sampling-coverage estimate.",
    "Gaussian coefficients and the PC1 loading-angle diagnostic are marked",
    "as point-only because dsVert does not claim nonlinear coefficient or",
    "loading regions; their coverage and width fields are therefore `NA`.",
    paste0(
      "The minimum recorded fixed-seed mechanism coverage is `",
      format(min(result$summary$mechanism_coverage, na.rm = TRUE),
             digits = 8L),
      "`; this Monte Carlo diagnostic is reported without replacing the",
      " analytic and deterministic coverage gates above."),
    "",
    "## Boundary cases",
    "",
    .dv_markdown_table(result$edge_cases),
    "",
    "## Interpretation",
    "",
    "The run exercises the real dsVertClient postprocessors and the deployed",
    "signed-vector radius calculation. Noise draws come from a base-R sampler",
    "for the same ideal two-sided-geometric convolution distribution; they do",
    "not exercise HMAC/HKDF/ChaCha20, Ring128, peer pinning, DSI transport,",
    "sticky replay, server admission, or signed receipt verification. Those",
    "remain separate protocol/E2E validation obligations.")
  writeLines(report, report_path, useBytes = TRUE)
  invisible(c(summary_path, contracts_path, edges_path, report_path))
}

.dv_main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  quick <- "--quick" %in% args
  replicate_arg <- args[startsWith(args, "--replicates=")]
  replicates <- if (length(replicate_arg)) {
    as.integer(sub("^--replicates=", "", replicate_arg[[1L]]))
  } else if (quick) {
    25L
  } else {
    1000L
  }
  output_arg <- args[startsWith(args, "--output-dir=")]
  output_dir <- if (length(output_arg)) {
    normalizePath(
      sub("^--output-dir=", "", output_arg[[1L]]), mustWork = TRUE)
  } else {
    getwd()
  }
  result <- dsvert_run_dp_statistical_validation(replicates)
  paths <- .dv_write_validation(result, output_dir)
  cat("DP statistical validation passed\n")
  cat(paste(paths, collapse = "\n"), "\n")
}

if (sys.nframe() == 0L) .dv_main()
