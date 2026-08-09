.dp_describe_client_status <- function(
    adjacency = "add_remove_patient") {
  list(
    policy = list(
      adjacency = adjacency,
      composition_partitions = 2L,
      local_total_epsilon = 0.5,
      local_total_delta = 0,
      decay = 0.5,
      numeric_grid_bits = 8L,
      numeric_bounds = list(x = c(0, 4), y = c(0, 8))),
    noise_root = list(
      privacy_epoch = 1,
      key_id = "test-noise-key-v1"))
}

.dp_describe_client_release <- function(
    adjacency = "add_remove_patient") {
  epsilon <- 0.25
  weights <- c(0.2, 0.3, 0.3, 0.2)
  family_epsilon <- epsilon * weights / 2
  sensitivity <- c(
    1, 256, 256,
    if (identical(adjacency, "add_remove_patient")) 1 else 2)
  coordinates <- 16L
  marginal <- vapply(seq_len(4L), function(index) {
    .dsvert_dp_google_ci_radius(
      family_epsilon[[index]], sensitivity[[index]])
  }, numeric(1L))
  simultaneous <- vapply(seq_len(4L), function(index) {
    .dsvert_dp_google_ci_radius(
      family_epsilon[[index]], sensitivity[[index]],
      alpha = 0.05 / coordinates)
  }, numeric(1L))
  l2_sensitivity <- sqrt(2 * (2 * 256^2 + if (identical(
    adjacency, "add_remove_patient")) 2 else 3))
  worst <- which.max(sensitivity / family_epsilon)
  laplace <- list(
    available = TRUE,
    mechanism = "dsvert_dp_v1_deterministic_granular_laplace_int64",
    epsilon = family_epsilon[[worst]], delta = 0, analytic_delta = 0,
    implementation_delta_bound = 0,
    accounting_rule = "pure_dp_no_implementation_slack",
    accuracy_accounting = "exact_granular_laplace_confidence_interval",
    sensitivity_norm = "l1", sensitivity = sensitivity[[worst]],
    marginal_95_abs = marginal[[worst]],
    simultaneous_95_abs = simultaneous[[worst]],
    nominal_rmse = sqrt(2) * sensitivity[[worst]] /
      family_epsilon[[worst]], sigma = 0,
    granularity = 2^ceiling(log2(
      (sensitivity[[worst]] / family_epsilon[[worst]]) / 2^40)),
    analytic_accounting_verified = TRUE, unavailable_reason = "")
  gaussian <- list(
    available = FALSE,
    mechanism = "dsvert_dp_v3_deterministic_approximate_gaussian_int64",
    epsilon = epsilon, delta = 0, analytic_delta = 0,
    implementation_delta_bound =
      .dsvert_dp_gaussian_implementation_delta_bound(coordinates, epsilon),
    accounting_rule =
      "analytic_gaussian_delta_plus_dp_transfer_from_total_variation_bound",
    accuracy_accounting =
      "gaussian_tail_alpha_minus_total_variation_union_bound",
    sensitivity_norm = "l2", sensitivity = l2_sensitivity,
    marginal_95_abs = 0, simultaneous_95_abs = 0,
    nominal_rmse = 0, sigma = 0, granularity = 0,
    analytic_accounting_verified = FALSE,
    unavailable_reason = "gaussian_delta_is_zero")
  noise_selection <- list(
    schema_version = 2L,
    selector = "minimum_conservative_95_radius_v3",
    objective = "simultaneous_95_abs", coordinate_count = coordinates,
    laplace = laplace, gaussian = gaussian, winner = "laplace",
    winner_mechanism = laplace$mechanism,
    winning_metric_abs = laplace$simultaneous_95_abs,
    winner_delta = 0,
    tie_break = "laplace_unless_gaussian_strictly_improves")
  block <- c(4, 512, 336, 1, 1, 1, 1, 1)
  list(
    released = TRUE,
    analysis_id = "primary",
    analysis_version = "v1",
    variables = c("x", "y"),
    variable_count = 2L,
    lower_bounds = c(0, 0),
    upper_bounds = c(4, 8),
    grid_lengths = c(4L, 4L),
    grid_values = c(1, 2, 3, 4, 2, 4, 6, 8),
    histogram_semantics =
      "(previous_endpoint,current_endpoint] plus fixed invalid bin",
    unit_collapse =
      "mean_of_finite_rows_after_public_bound_clipping",
    count_definition =
      "DP-noisy effective units with at least one finite bounded value",
    invalid_unit_rule =
      "invalid_patient_ids_rejected_by_admission",
    statistics = c(block, block),
    coordinate_count = coordinates,
    coordinate_layout = paste(
      "per_variable[count,qsum,qsumsq,",
      "histogram[grid_bins+invalid]]"),
    numeric_grid_bits = 8L,
    numeric_grid_scale = 256,
    quantization =
      "round(z*scale) and round(z^2*scale) independently",
    max_abs_normalized_quantization_per_unit = 0.5 / 256,
    allocation_names = c("count", "sum", "sumsq", "histogram"),
    allocation_weights = weights,
    epsilon_per_variable_family = family_epsilon,
    epsilon_allocation_sum = epsilon,
    submechanism_count = 8L,
    composition_rule = paste(
      "sequential composition across variables and count/sum/sumsq/",
      "histogram families; each histogram is one L1 vector mechanism"),
    calibration_rule = paste(
      "each family is calibrated with its own epsilon and sensitivity;",
      "the summed sensitivity bound is descriptive, not a shared scale"),
    family_sensitivity = sensitivity,
    sum_family_l1_sensitivity_bound = 2 * sum(sensitivity),
    l2_sensitivity_bound = l2_sensitivity,
    noise_selection = noise_selection,
    mechanism = "dsvert_dp_v1_deterministic_granular_laplace_int64",
    implementation = paste0(
      "dsVert adapted Google Differential Privacy v4.1.0 ",
      "granular Laplace integer mechanism"),
    sampler = "deterministic_two_sided_geometric",
    randomness = "HMAC-SHA256/ChaCha20",
    postprocessing = "coordinatewise_nonnegative_integer",
    clipped_coordinates = 0L,
    accuracy_95_abs_by_variable_family = rep(marginal, times = 2L),
    accuracy_simultaneous_95_abs_by_variable_family =
      rep(simultaneous, times = 2L),
    accuracy_simultaneous_confidence = 0.95,
    accuracy_simultaneous_method = "union_bound",
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    privacy_epoch = 1,
    noise_key_id = "test-noise-key-v1",
    sticky_noise = "dsvert-sticky-noise-v1",
    epsilon = epsilon,
    delta = 0,
    adjacency = adjacency,
    composition_partitions = 2L)
}

.dp_describe_vector_capsule <- function(decoded = TRUE, gaussian = FALSE) {
  numeric <- list(artifacts = list(
    x = list(
      owner_peer = "site_a", dataset = "protected", column = "x",
      lower = 0, upper = 4, numeric_grid_bits = 8L,
      statistic_maximum = c(100, 25600, 25600)),
    y = list(
      owner_peer = "site_a", dataset = "protected", column = "y",
      lower = 0, upper = 8, numeric_grid_bits = 8L,
      statistic_maximum = c(100, 25600, 25600))))
  histograms <- list(artifacts = list(
    h_x = list(
      owner_peer = "site_a", dataset = "protected", column = "x",
      grid = c(1, 2, 3, 4), coordinate_count = 5L,
      statistic_maximum = 100),
    h_y = list(
      owner_peer = "site_a", dataset = "protected", column = "y",
      grid = c(2, 4, 6, 8), coordinate_count = 5L,
      statistic_maximum = 100)))
  describe <- list(primary = list(
    version = "v1", dataset = "protected", owner_peer = "site_a",
    variables = c("x", "y"),
    histogram_grids = list(x = c(1, 2, 3, 4), y = c(2, 4, 6, 8)),
    histogram_references = list(
      list(family = "fixed_numeric_histograms", primitive_id = "h_x"),
      list(family = "fixed_numeric_histograms", primitive_id = "h_y")),
    allocation_weights = c(0.2, 0.3, 0.3, 0.2)))
  manifest <- list(workload = list(
    coordinate_count = 17L,
    release_lattice = list(
      output_lattice_bits = 8L, output_lattice_scale = 256,
      natural_l1_sensitivity = 10,
      integer_l1_sensitivity_steps = 2560,
      natural_l2_sensitivity = sqrt(10),
      integer_l2_sensitivity_steps = sqrt(10) * 256),
    capsule_mechanism = list(
      mechanism = if (isTRUE(gaussian)) {
        .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM
      } else {
        "discrete-laplace"
      },
      sensitivity_norm = if (isTRUE(gaussian)) "l2" else "l1"),
    families = list(
      admitted_count = list(
        owner_peer = "site_a", dataset = "protected",
        statistic_maximum = 100),
      numeric_moments = numeric,
      numeric_pair_moments = list(artifacts = list()),
      gaussian_models = list(artifacts = list()),
      fixed_numeric_histograms = histograms,
      categorical_marginals = list(artifacts = list()),
      categorical_pairs = list(sets = list()), correlation_artifacts = list(),
      describe_artifacts = describe, survival_artifacts = list())))
  if (isTRUE(gaussian)) {
    sensitivity_steps <- format(
      sqrt(10) * 256, digits = 17L, scientific = TRUE, trim = TRUE)
    plan <- list(
      version = .DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION,
      mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
      sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER,
      request_binding_sha256 = digest::digest(
        "describe-gaussian-request", "sha256", serialize = FALSE),
      total_coordinate_count = 17,
      maximum_chunk_coordinates = 17,
      independent_noise_peer_count = 2,
      complete_epsilon_per_peer = TRUE,
      epsilon_divided_by_peer_count = FALSE,
      capability_available = TRUE,
      release_delta_aggregation = "max_per_peer_not_sum",
      exact_rational_sampler = FALSE,
      finite_support_transfer_charged = TRUE,
      fixed_work_sampler = TRUE,
      sampler_branches_on_protected_values = FALSE,
    sampler_branches_on_private_randomness = FALSE,
      host_constant_time_claim = FALSE,
    transcript_dp_claim = TRUE,
      sampler_candidate_count = 1,
      sampler_random_bits_per_coordinate = 128,
      sampler_random_bytes_per_coordinate = 17,
      sampler_table_precision_bits = 192,
      sampler_magnitude_count = 33,
      sampler_search_steps = 6,
      vector_tail_tv_upper_numerator = "1",
      vector_tail_tv_upper_denominator =
        "1000000000000000000000000000000",
      vector_sampler_tv_upper_numerator = "1",
      vector_sampler_tv_upper_denominator =
        "1000000000000000000000000000000",
      vector_total_tv_upper_numerator = "2",
      vector_total_tv_upper_denominator =
        "1000000000000000000000000000000",
      observable_worker_shape = "fixed dyadic CDF fixture",
      per_peer_implementation_delta_numerator = "1",
      per_peer_implementation_delta_denominator =
        "1000000000000000000000000000000",
      simultaneous_95_abs = "20")
    plan <- .client_complete_gaussian_plan_v2(plan)
    manifest$workload$mechanism_selection <- list(
      gaussian_calibration_request = list(
        epsilon = "1e+00", delta = format(
          2^-100, digits = 17L, scientific = TRUE, trim = TRUE),
        l2_sensitivity_steps = sensitivity_steps,
        total_coordinate_count = 17),
      gaussian_plan = plan,
      gaussian_plan_sha256 = .dsvert_vector_hash(plan))
  }
  if (isTRUE(decoded)) {
    manifest <- jsonlite::fromJSON(
      .dsvert_joint_dp_client_json(manifest), simplifyVector = FALSE)
  }
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  release <- list(
    capsule_id = strrep("a", 64L), final_vector_root = strrep("b", 64L),
    coordinate_order_sha256 = layout$sha256,
    coordinate_count = 17L,
    values = c(
      4,
      4, 2, 1.3125,
      4, 2, 1.3125,
      rep(1, 5L), rep(1, 5L)),
    epsilon = 1, delta = 2^-100,
    implementation_delta =
      "1/1000000000000000000000000000000",
    mechanism = if (isTRUE(gaussian)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM
    } else {
      .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM
    },
    manifest = manifest,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE)
  class(release) <- c("dsvert_joint_dp_vector", "list")
  status <- list(site_a = list(
    policy = list(adjacency = "add_remove_patient", unit_capacity = 100),
    noise_root = list(privacy_epoch = 1, key_id = "noise-key-a")))
  list(release = release, layout = layout, status = status)
}

test_that("describe maps only signed final vector blocks", {
  capsule <- .dp_describe_vector_capsule(decoded = TRUE)
  result <- .dsvert_dp_describe_vector_result(
    capsule, "protected", "primary")
  expect_identical(result$server, "site_a")
  expect_true(result$history_gate)
  expect_false(result$request_limit)
  expect_true(result$operation_limit)
  expect_identical(result$statistics,
                   c(4, 512, 336, rep(1, 5L),
                     4, 512, 336, rep(1, 5L)))
  expect_true(is.na(result$clipped_coordinates))
  expect_false(result$clipping_observable)
  expect_identical(result$security_claim,
                   .dsvert_dp_capsule_security_claim())
  expect_false(result$security_claim$analyst_relay_trusted)
  expect_false(
    result$security_claim$unconditional_non_reconstruction_guarantee)
  postprocessed <- .dsvert_dp_describe_postprocess(result, probs = 0.5)
  expect_equal(postprocessed$descriptives$mean, c(2, 4), tolerance = 0)
  expect_equal(postprocessed$descriptives$variance, c(1.25, 5),
               tolerance = 0)
  expect_identical(
    postprocessed$status,
    "fixed_public_clamp_applied_preclamp_state_not_released")
  expect_error(.dsvert_dp_describe_vector_result(
    capsule, "protected", "primary", server = "site_b"), "does not own")
  stale_gate <- capsule
  stale_gate$release$history_gate <- FALSE
  expect_error(.dsvert_dp_describe_vector_result(
    stale_gate, "protected", "primary"), "vector context is invalid")
})

test_that("describe propagates the signed Gaussian L2 plan", {
  capsule <- .dp_describe_vector_capsule(
    decoded = FALSE, gaussian = TRUE)
  result <- .dsvert_dp_describe_vector_result(
    capsule, "protected", "primary")
  expect_identical(result$mechanism,
                   .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
  expect_identical(result$sampler, .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER)
  expect_match(result$accuracy_simultaneous_method,
               "discrete-Gaussian plan", fixed = TRUE)
  expect_equal(
    result$accuracy_simultaneous_95_abs_by_variable_family,
    rep(c(20 / 256, 20, 20, 20 / 256), 2L), tolerance = 0)
  postprocessed <- .dsvert_dp_describe_postprocess(result, probs = 0.5)
  expect_identical(postprocessed$mechanism,
                   .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
})

test_that("client validates explicit sequential epsilon composition", {
  result <- .dsvert_dp_validate_describe_release(
    .dp_describe_client_release(), "primary", "site_a",
    .dp_describe_client_status())
  expect_equal(
    result$variable_count * sum(result$epsilon_per_variable_family),
    result$epsilon, tolerance = 1e-15)
  expect_equal(
    result$epsilon_per_variable_family,
    c(0.025, 0.0375, 0.0375, 0.025))
  expect_match(result$composition_rule, "each histogram is one L1 vector")
  expect_false(any(grepl(
    "observed_min|observed_max|exact|raw", names(result),
    ignore.case = TRUE)))
})

test_that("client enforces replacement histogram sensitivity", {
  status <- .dp_describe_client_status("replace_one_fixed_cohort")
  release <- .dp_describe_client_release("replace_one_fixed_cohort")
  result <- .dsvert_dp_validate_describe_release(
    release, "primary", "site_a", status)
  expect_identical(result$family_sensitivity, c(1, 256, 256, 2))
  expect_identical(result$sum_family_l1_sensitivity_bound, 1030)

  release$family_sensitivity[[4L]] <- 1
  expect_error(
    .dsvert_dp_validate_describe_release(
      release, "primary", "site_a", status),
    "invalid DP describe release")
})

test_that("no-noise descriptives and quantiles match the binned oracle", {
  result <- .dp_describe_client_release()
  result$accuracy_95_abs_by_variable_family[] <- 0
  result$accuracy_simultaneous_95_abs_by_variable_family[] <- 0
  result$server <- "site_a"
  result <- .dsvert_dp_describe_postprocess(
    result, probs = c(0.25, 0.5, 0.75))
  expect_equal(result$descriptives$n_dp, c(4, 4))
  expect_equal(result$descriptives$mean, c(2, 4), tolerance = 1e-14)
  expect_equal(result$descriptives$variance, c(1.25, 5),
               tolerance = 1e-14)
  expect_equal(result$descriptives$sd, sqrt(c(1.25, 5)),
               tolerance = 1e-14)
  expect_equal(
    result$quantiles$estimate[result$quantiles$variable == "x"],
    c(1, 2, 3))
  expect_equal(
    result$quantiles$estimate[result$quantiles$variable == "y"],
    c(2, 4, 6))
  expect_equal(result$histograms$x$counts, rep(1, 4L))
  expect_identical(result$histograms$x$invalid_dp, 1)
})

test_that("moment projection remains bounded and internally consistent", {
  result <- .dp_describe_client_release()
  result$statistics[1:3] <- c(2, 1e6, 0)
  result <- .dsvert_dp_describe_postprocess(result, probs = 0.5)
  x <- result$descriptives[result$descriptives$variable == "x", ]
  expect_gte(x$mean, x$lower_bound)
  expect_lte(x$mean, x$upper_bound)
  expect_gte(x$variance, 0)
  expect_lte(
    x$variance,
    (x$mean - x$lower_bound) * (x$upper_bound - x$mean) + 1e-12)
  expect_equal(x$sd^2, x$variance, tolerance = 1e-14)
})

test_that("moment release regions cover quantisation and the simultaneous box", {
  region <- .dsvert_dp_describe_moment_region(
    n_dp = 3, qsum_dp = 7, qsumsq_dp = 5,
    count_radius = 1, sum_radius = 2, sumsq_radius = 2,
    grid_scale = 4, lower_bound = -2, upper_bound = 6)

  candidates <- list()
  values <- seq(0, 1, by = 0.1)
  for (n in 1:4) {
    indices <- expand.grid(rep(list(seq_along(values)), n),
                           KEEP.OUT.ATTRS = FALSE)
    for (row in seq_len(nrow(indices))) {
      z <- values[as.integer(indices[row, ])]
      sufficient <- c(
        n, sum(round(z * 4)), sum(round(z^2 * 4)))
      if (all(abs(sufficient - c(3, 7, 5)) <= c(1, 2, 2))) {
        x <- -2 + 8 * z
        candidates[[length(candidates) + 1L]] <- c(
          mean = mean(x), variance = mean((x - mean(x))^2))
      }
    }
  }
  candidates <- do.call(rbind, candidates)
  expect_gt(nrow(candidates), 0L)
  expect_true(all(candidates[, "mean"] >=
                    region$mean[["lower"]] - 1e-12))
  expect_true(all(candidates[, "mean"] <=
                    region$mean[["upper"]] + 1e-12))
  expect_true(all(candidates[, "variance"] >=
                    region$variance[["lower"]] - 1e-12))
  expect_true(all(candidates[, "variance"] <=
                    region$variance[["upper"]] + 1e-12))
  expect_true(region$status %in% c(
    "ok", "dp_region_includes_zero_effective_count"))
})

test_that("moment regions cover bounded datasets after postprocessed noise", {
  set.seed(20260801)
  covered <- logical(250L)
  for (iteration in seq_along(covered)) {
    n <- sample.int(6L, 1L)
    z <- sample(seq(0, 1, by = 0.01), n, replace = TRUE)
    exact <- c(n, sum(round(z * 16)), sum(round(z^2 * 16)))
    radii <- c(2, 4, 4)
    noise <- vapply(radii, function(radius) {
      sample.int(2L * radius + 1L, 1L) - radius - 1L
    }, numeric(1L))
    observed <- pmax(0, exact + noise)
    region <- .dsvert_dp_describe_moment_region(
      observed[[1L]], observed[[2L]], observed[[3L]],
      radii[[1L]], radii[[2L]], radii[[3L]], 16, -10, 7)
    original <- -10 + 17 * z
    truth <- c(mean = mean(original),
               variance = mean((original - mean(original))^2))
    covered[[iteration]] <-
      truth[["mean"]] >= region$mean[["lower"]] - 1e-12 &&
      truth[["mean"]] <= region$mean[["upper"]] + 1e-12 &&
      truth[["variance"]] >= region$variance[["lower"]] - 1e-12 &&
      truth[["variance"]] <= region$variance[["upper"]] + 1e-12
  }
  expect_true(all(covered))
})

test_that("moment regions cover randomized scales, domains, and noise boxes", {
  withr::local_seed(20260802)
  covered <- logical(1000L)
  for (iteration in seq_along(covered)) {
    n <- sample.int(20L, 1L)
    z <- runif(n)
    grid_scale <- 2^sample(8:18, 1L)
    exact <- c(
      n,
      sum(round(z * grid_scale)),
      sum(round(z^2 * grid_scale)))
    radii <- c(
      sample(0:5, 1L),
      sample.int(2L * grid_scale + 1L, 1L) - 1L,
      sample.int(2L * grid_scale + 1L, 1L) - 1L)
    noise <- vapply(radii, function(radius) {
      if (!radius) return(0)
      sample.int(2L * radius + 1L, 1L) - radius - 1L
    }, numeric(1L))
    observed <- pmax(0, exact + noise)
    lower <- runif(1L, -1000, 1000)
    width <- 10^runif(1L, -4, 4)
    upper <- lower + width
    region <- .dsvert_dp_describe_moment_region(
      observed[[1L]], observed[[2L]], observed[[3L]],
      radii[[1L]], radii[[2L]], radii[[3L]],
      grid_scale, lower, upper)
    original <- lower + width * z
    truth <- c(
      mean = mean(original),
      variance = mean((original - mean(original))^2))
    tolerance <- 1e-8 * max(1, abs(truth), width^2)
    covered[[iteration]] <-
      truth[["mean"]] >= region$mean[["lower"]] - tolerance &&
      truth[["mean"]] <= region$mean[["upper"]] + tolerance &&
      truth[["variance"]] >= region$variance[["lower"]] - tolerance &&
      truth[["variance"]] <= region$variance[["upper"]] + tolerance
  }
  expect_true(all(covered))
})

test_that("describe exposes typed non-sampling moment regions", {
  result <- .dp_describe_client_release()
  result <- .dsvert_dp_describe_postprocess(result, probs = 0.5)
  expected <- c(
    "mean_mechanism_grid_lower_95", "mean_mechanism_grid_upper_95",
    "variance_mechanism_grid_lower_95",
    "variance_mechanism_grid_upper_95", "moment_region_status")
  expect_true(all(expected %in% names(result$descriptives)))
  expect_true(all(
    result$descriptives$mean_mechanism_grid_lower_95 <=
      result$descriptives$mean_mechanism_grid_upper_95))
  expect_true(all(
    result$descriptives$variance_mechanism_grid_lower_95 <=
      result$descriptives$variance_mechanism_grid_upper_95))
  expect_identical(result$moment_region_confidence, 0.95)
  expect_match(result$moment_region_scope, "quantization")
  expect_match(result$moment_region_scope, "sampling uncertainty excluded")
})

test_that("uncertified positive noisy counts retain an explicit warning", {
  result <- .dsvert_dp_validate_describe_release(
    .dp_describe_client_release(), "primary", "site_a",
    .dp_describe_client_status())
  result <- .dsvert_dp_describe_postprocess(result, probs = 0.5)
  expect_true(all(result$descriptives$status ==
    "dp_point_available_count_not_certified_positive"))
  expect_true(all(is.finite(result$descriptives$mean)))
  expect_true(all(result$descriptives$n_simultaneous_95_lower == 0))
})

test_that("quantile points and DP-grid bands are monotone", {
  result <- .dp_describe_client_release()
  result$server <- "site_a"
  result <- .dsvert_dp_describe_postprocess(
    result, probs = c(0.9, 0.1, 0.5, 0.5))
  expect_identical(unique(result$quantiles$probability), c(0.1, 0.5, 0.9))
  for (variable in result$variables) {
    q <- result$quantiles[result$quantiles$variable == variable, ]
    expect_true(all(diff(q$estimate) >= 0))
    expect_true(all(diff(q$dp_grid_lower) >= 0))
    expect_true(all(diff(q$dp_grid_upper) >= 0))
    expect_true(all(q$dp_grid_lower <= q$dp_grid_upper))
  }
  expect_match(result$quantile_band_scope, "sampling uncertainty excluded")
  expect_match(result$statistical_inference, "no sampling confidence")
})

test_that("quantile bands contain every table in the simultaneous box", {
  observed <- c(2, 1, 3)
  radius <- 1
  grid <- c(1, 2, 3)
  probs <- c(0.25, 0.5, 0.75)
  band <- .dsvert_dp_describe_quantiles(
    observed, grid, lower_bound = 0, probs, radius)
  lower <- pmax(0, observed - radius)
  upper <- observed + radius
  tables <- expand.grid(lapply(seq_along(observed), function(index) {
    seq.int(lower[[index]], upper[[index]])
  }), KEEP.OUT.ATTRS = FALSE)
  covered <- TRUE
  checked <- 0L
  for (row in seq_len(nrow(tables))) {
    counts <- as.numeric(tables[row, ])
    total <- sum(counts)
    if (!total) next
    cumulative <- cumsum(counts)
    for (index in seq_along(probs)) {
      quantile_bin <- which(cumulative >= probs[[index]] * total)[[1L]]
      quantile_endpoint <- grid[[quantile_bin]]
      covered <- covered &&
        quantile_endpoint >= band$dp_grid_lower[[index]] &&
        quantile_endpoint <= band$dp_grid_upper[[index]]
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 0L)
  expect_true(covered)
})

test_that("quantile bands exhaustively cover small simultaneous boxes", {
  grid <- c(1, 2, 3)
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  observed_tables <- expand.grid(rep(list(0:2), length(grid)),
                                 KEEP.OUT.ATTRS = FALSE)
  covered <- TRUE
  checked <- 0L
  failure <- ""
  for (observed_index in seq_len(nrow(observed_tables))) {
    observed <- as.numeric(observed_tables[observed_index, ])
    for (radius in 0:2) {
      band <- .dsvert_dp_describe_quantiles(
        observed, grid, lower_bound = 0, probs, radius)
      lower <- pmax(0, observed - radius)
      upper <- observed + radius
      candidate_tables <- expand.grid(
        Map(seq.int, lower, upper), KEEP.OUT.ATTRS = FALSE)
      for (candidate_index in seq_len(nrow(candidate_tables))) {
        counts <- as.numeric(candidate_tables[candidate_index, ])
        total <- sum(counts)
        if (!total) next
        cumulative <- cumsum(counts)
        endpoints <- vapply(probs, function(probability) {
          grid[[which(cumulative >= probability * total)[[1L]]]]
        }, numeric(1L))
        inside <- endpoints >= band$dp_grid_lower &
          endpoints <= band$dp_grid_upper
        if (!all(inside)) {
          covered <- FALSE
          failure <- paste(
            "observed", paste(observed, collapse = ","),
            "radius", radius,
            "candidate", paste(counts, collapse = ","))
        }
        checked <- checked + length(probs)
      }
    }
  }
  expect_identical(checked, 11155L)
  expect_true(covered, info = failure)
})

test_that("empty noisy histograms return typed full-domain quantile bands", {
  result <- .dp_describe_client_release()
  result$statistics[c(4:7, 12:15)] <- 0
  result$accuracy_simultaneous_95_abs_by_variable_family[] <- 0
  result <- .dsvert_dp_describe_postprocess(result, probs = c(0.25, 0.75))
  expect_true(all(is.na(result$quantiles$estimate)))
  expect_true(all(result$quantiles$status ==
                  "dp_histogram_empty_after_postprocessing"))
  x <- result$quantiles[result$quantiles$variable == "x", ]
  expect_equal(x$dp_grid_lower, c(0, 0))
  expect_equal(x$dp_grid_upper, c(4, 4))
})

test_that("client rejects malformed, fractional, or expanded releases", {
  status <- .dp_describe_client_status()
  base <- .dp_describe_client_release()
  malformed <- list(
    utils::modifyList(base, list(epsilon_allocation_sum = 0.3)),
    utils::modifyList(base, list(statistics = c(4.5, base$statistics[-1L]))),
    utils::modifyList(base, list(grid_lengths = c(5L, 4L))),
    utils::modifyList(base, list(
      composition_rule = "full epsilon independently for every statistic")),
    c(base, list(observed_minimum = 0)))
  for (value in malformed) {
    expect_error(
      .dsvert_dp_validate_describe_release(
        value, "primary", "site_a", status),
      "invalid DP describe release")
  }
})

test_that("describe validation rejects non-finite and overflowing coordinates", {
  status <- .dp_describe_client_status()
  for (bad in c(NA_real_, NaN, Inf, -Inf, 2^53)) {
    release <- .dp_describe_client_release()
    release$statistics[[1L]] <- bad
    expect_error(
      .dsvert_dp_validate_describe_release(
        release, "primary", "site_a", status),
      "invalid DP describe release")
  }
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    release <- .dp_describe_client_release()
    release$grid_values[[1L]] <- bad
    expect_error(
      .dsvert_dp_validate_describe_release(
        release, "primary", "site_a", status),
      "invalid DP describe release")
  }
})

test_that("probs never alter the canonical sticky server request", {
  conns <- list(site_a = structure(1, class = "fake"))
  calls <- 0L
  evaluate <- function(probs) {
    .dsvert_dp_describe_impl(
      "protected", "primary", probs, "site_a", conns,
      function(...) stop("unexpected direct aggregate"))
  }
  result <- testthat::with_mocked_bindings(
    list(evaluate(0.5), evaluate(c(0.1, 0.9))),
    .dsvert_dp_capsule_vector_run = function(
        datasources, status = NULL, .aggregate) {
      calls <<- calls + 1L
      .dp_describe_vector_capsule(decoded = TRUE)
    },
    .package = "dsVertClient")
  expect_identical(calls, 2L)
  expect_identical(result[[1L]]$final_vector_root,
                   result[[2L]]$final_vector_root)
  expect_identical(result[[1L]]$statistics, result[[2L]]$statistics)
  expect_length(result[[1L]]$quantiles$probability, 2L)
  expect_length(result[[2L]]$quantiles$probability, 4L)
})

test_that("describe release rejects DSI failures without raw text or retry", {
  conns <- list(site_a = structure(1, class = "fake"))
  for (kind in c("throw", "partial", "callback")) {
    attempts <- 0L
    aggregate <- function(conns, expr, error = NULL,
                          errors.print = TRUE, ...) {
      attempts <<- attempts + 1L
      if (kind == "throw") stop("SECRET_REMOTE_DETAIL", call. = FALSE)
      if (kind == "callback") {
        error("site_a", "SECRET_REMOTE_DETAIL")
      }
      setNames(list(NULL), names(conns))
    }
    evaluate <- function() {
      .dsvert_dp_describe_impl(
        "protected", "primary", 0.5, "site_a", conns, aggregate)
    }
    captured <- tryCatch(evaluate(), error = function(e) conditionMessage(e))
    expect_true(grepl(
      "DataSHIELD transport failed during 'reusable joint-DP capsule status'",
      captured, fixed = TRUE))
    expect_false(grepl("SECRET_REMOTE_DETAIL", captured, fixed = TRUE))
    expect_identical(attempts, 1L)
  }
})
