.dp_describe_client_release <- function() {
  .dsvert_dp_describe_vector_result(
    .dp_describe_vector_capsule(decoded = TRUE),
    "protected", "primary")
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
    backend = if (isTRUE(gaussian)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND
    } else {
      .DSVERT_CLIENT_VECTOR_BACKEND
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

test_that("describe admits synopsis provenance without legacy claims", {
  capsule <- .dp_describe_vector_capsule(decoded = TRUE)
  capsule$release$artifact_key <- strrep("c", 64L)
  capsule$release$execution_id <- strrep("d", 64L)
  capsule$release$contract_sha256 <- strrep("e", 64L)
  capsule$release$attempt_sha256 <- strrep("f", 64L)
  capsule$release$source_contract_sha256 <- strrep("1", 64L)
  capsule$release$result_set_sha256 <- strrep("2", 64L)
  capsule$release$signed_provenance <- c(list(
    version = "dsvert-stateless-synopsis-public-provenance-v1"),
    capsule$release[c(
      "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
      "source_contract_sha256", "result_set_sha256", "final_vector_root")])
  class(capsule$release) <- c(
    "dsvert_synopsis_public_vector", "dsvert_joint_dp_vector", "list")
  capsule$release[c(
    "capsule_id", "history_gate", "request_limit", "operation_limit")] <-
    NULL
  capsule$release$manifest$admission <- list(
    adjacency = "add_remove_patient", unit_capacity = 100)
  capsule$status[[1L]]$noise_root <- NULL

  expect_error(.dsvert_dp_vector_context(capsule), "does not yet accept")
  result <- .dsvert_dp_describe_vector_result(
    capsule, "protected", "primary")
  expect_identical(result$statistics,
                   c(4, 512, 336, rep(1, 5L),
                     4, 512, 336, rep(1, 5L)))
  expect_identical(result$artifact_key, strrep("c", 64L))
  expect_identical(result$execution_id, strrep("d", 64L))
  expect_identical(result$contract_sha256, strrep("e", 64L))
  expect_identical(result$attempt_sha256, strrep("f", 64L))
  expect_identical(result$source_contract_sha256, strrep("1", 64L))
  expect_identical(result$result_set_sha256, strrep("2", 64L))
  expect_identical(result$privacy$version, "dsvert-per-synopsis-dp-v1")
  expect_false(result$privacy$finite_global_composition_claim)
  expect_false(result$security_claim$allocation_openings_used)
  expect_false(any(c(
    "capsule_id", "privacy_epoch", "noise_key_id", "history_gate",
    "request_limit", "operation_limit") %in% names(result)))
  postprocessed <- .dsvert_dp_describe_postprocess(result, probs = 0.5)
  class(postprocessed) <- c("ds.vertDPDescribe", "list")
  expect_no_error(ds.vertDPQuantile(postprocessed, 0.5))
  expect_equal(postprocessed$descriptives$mean, c(2, 4), tolerance = 0)
  expect_equal(postprocessed$descriptives$variance, c(1.25, 5),
               tolerance = 0)
  capsule$release$manifest$admission$unit_capacity <- 99L
  expect_error(.dsvert_dp_describe_vector_result(
    capsule, "protected", "primary"), "capacity")
  capsule$release$manifest$admission$unit_capacity <- 100L
  capsule$release$artifact_key <- strrep("a", 64L)
  expect_error(.dsvert_dp_describe_vector_result(
    capsule, "protected", "primary"), "provenance")
})

test_that("describe vector contract rejects shape and lattice tampering", {
  base <- .dp_describe_vector_capsule(decoded = TRUE)
  reject <- function(change, pattern) {
    candidate <- change(base)
    expect_error(
      .dsvert_dp_describe_vector_result(
        candidate, "protected", "primary"),
      pattern)
  }
  reject(function(x) {
    x$release$manifest$workload$families$describe_artifacts$primary$
      histogram_references[[1L]]$primitive_id <- "missing"
    x
  }, "does not contain exactly one fixed histogram block")
  reject(function(x) {
    x$layout$blocks[["fixed_numeric_histograms::h_x"]]$descriptor$
      grid[[1L]] <- 0.5
    x
  }, "histogram reference is inconsistent")
  reject(function(x) {
    x$layout$blocks[["fixed_numeric_histograms::h_x"]]$length <- 4L
    x
  }, "histogram reference is inconsistent")
  reject(function(x) {
    x$release$manifest$workload$release_lattice$
      natural_l1_sensitivity <- 0
    x
  }, "release lattice is invalid")
  reject(function(x) {
    start <- x$layout$blocks[["numeric_moments::x"]]$start
    x$release$values[[start + 1L]] <-
      x$release$values[[start + 1L]] + 1 / 512
    x
  }, "not on its signed lattice")
})

test_that("describe preserves a fractional noisy count on the common lattice", {
  capsule <- .dp_describe_vector_capsule(decoded = TRUE)
  block <- capsule$layout$blocks[["numeric_moments::x"]]
  capsule$release$values[block$start:block$end] <-
    c(46.6953125, 24.89453125, 17.0625)

  result <- .dsvert_dp_describe_vector_result(
    capsule, "protected", "primary")
  expect_identical(result$statistics[1:3], c(46.6953125, 6373, 4368))

  postprocessed <- .dsvert_dp_describe_postprocess(result, probs = 0.5)
  x <- postprocessed$descriptives[
    postprocessed$descriptives$variable == "x", ]
  expect_identical(x$n_dp, 46.6953125)
  expect_true(all(is.finite(c(x$mean, x$variance, x$sd))))
  expect_true(x$mean >= x$lower_bound && x$mean <= x$upper_bound)
  expect_true(x$variance >= 0 &&
                x$variance <= (x$upper_bound - x$lower_bound)^2 / 4)
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
  result <- .dp_describe_client_release()
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

test_that("probs never alter the canonical sticky server request", {
  conns <- list(site_a = structure(1, class = "fake"))
  synopsis_calls <- 0L
  legacy_calls <- 0L
  evaluate <- function(probs) {
    .dsvert_dp_describe_impl(
      "protected", "primary", probs, "site_a", conns,
      function(...) stop("unexpected direct aggregate"))
  }
  result <- testthat::with_mocked_bindings(
    list(evaluate(0.5), evaluate(c(0.1, 0.9))),
    .dsvert_dp_synopsis_vector_run = function(
        datasources, status, .aggregate) {
      expect_null(status)
      synopsis_calls <<- synopsis_calls + 1L
      .dp_describe_vector_capsule(decoded = TRUE)
    },
    .dsvert_dp_describe_resume_token_v1 = function(...) structure(
      list(), class = c("dsvertDPDescribeResume", "list")),
    .dsvert_dp_capsule_vector_run = function(...) {
      legacy_calls <<- legacy_calls + 1L
      stop("legacy runner reached", call. = FALSE)
    },
    .package = "dsVertClient")
  expect_identical(synopsis_calls, 2L)
  expect_identical(legacy_calls, 0L)
  expect_identical(result[[1L]]$final_vector_root,
                   result[[2L]]$final_vector_root)
  expect_identical(result[[1L]]$statistics, result[[2L]]$statistics)
  expect_length(result[[1L]]$quantiles$probability, 2L)
  expect_length(result[[2L]]$quantiles$probability, 4L)
})

test_that("Describe exposes a portable bound resume token and cleanup state", {
  capsule <- .dp_describe_vector_capsule(decoded = TRUE)
  capsule$release$artifact_key <- strrep("c", 64L)
  capsule$release$execution_id <- strrep("d", 64L)
  capsule$release$contract_sha256 <- strrep("e", 64L)
  capsule$release$attempt_sha256 <- strrep("f", 64L)
  capsule$release$source_contract_sha256 <- strrep("1", 64L)
  capsule$release$result_set_sha256 <- strrep("2", 64L)
  capsule$release$signed_provenance <- c(list(
    version = "dsvert-stateless-synopsis-public-provenance-v1"),
    capsule$release[c(
      "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
      "source_contract_sha256", "result_set_sha256", "final_vector_root")])
  class(capsule$release) <- c(
    "dsvert_synopsis_public_vector", "dsvert_joint_dp_vector", "list")
  capsule$release[c(
    "capsule_id", "history_gate", "request_limit", "operation_limit")] <-
    NULL
  capsule$release$manifest$admission <- list(
    adjacency = "add_remove_patient", unit_capacity = 100)
  capsule$status[[1L]]$noise_root <- NULL
  class(capsule$status) <- c("ds.vertDPSynopsisStatus", "list")
  manifest_json <- .dsvert_joint_dp_client_json(capsule$release$manifest)
  conns <- list(site_a = structure(list(peer = "site_a"),
                                     class = "fake_connection"))
  context <- list(
    status = capsule$status, servers = "site_a", pinset = c(site_a = "pin"),
    designated = "site_a", policy = list(version = "test"),
    conns = conns, all_conns = conns)
  capsule$manifest_bundle <- list(
    manifest_json = manifest_json,
    manifest_sha256 = digest::digest(
      manifest_json, "sha256", serialize = FALSE),
    capsule_id = strrep("a", 64L), context = context)
  capsule$cleanup_pending <- TRUE

  statuses <- list()
  calls <- 0L
  evaluate <- function(resume = NULL) .dsvert_dp_describe_impl(
    "protected", "primary", 0.5, "site_a", conns,
    function(...) stop("unexpected aggregate"), resume = resume)
  run <- function(resume = NULL) testthat::with_mocked_bindings(
    evaluate(resume),
    .dsvert_dp_synopsis_vector_run = function(
      datasources, status, .aggregate) {
      calls <<- calls + 1L
      statuses[calls] <<- list(status)
      current <- capsule
      current$cleanup_pending <- calls == 1L
      current
    }, .package = "dsVertClient")

  first <- run()
  expect_s3_class(first$resume, "dsvertDPDescribeResume")
  expect_true(first$cleanup_pending)
  expect_null(statuses[[1L]])
  connection_reference <- function(value) {
    if (inherits(value, "fake_connection")) return(TRUE)
    is.list(value) && any(vapply(value, connection_reference, logical(1L)))
  }
  expect_false(connection_reference(first$resume))

  portable <- unserialize(serialize(first$resume, NULL))
  expect_false(connection_reference(portable))
  second <- run(portable)
  expect_false(second$cleanup_pending)
  expect_s3_class(statuses[[2L]], "dsvert_synopsis_bootstrap_v1")
  expect_identical(second$statistics, first$statistics)
  expect_no_error(run(first))

  tampered <- first$resume
  tampered$binding_sha256 <- strrep("0", 64L)
  expect_error(run(tampered), "invalid|misbound", ignore.case = TRUE)
  tampered <- first$resume
  tampered$method <- "ds.vertDPMeanVar"
  expect_error(run(tampered), "invalid|misbound", ignore.case = TRUE)
  tampered <- first$resume
  tampered$artifact_sha256 <- strrep("0", 64L)
  expect_error(run(tampered), "invalid|misbound", ignore.case = TRUE)
  tampered <- first$resume
  tampered$bootstrap$manifest_bundle$manifest_sha256 <- strrep("0", 64L)
  expect_error(run(tampered), "invalid|misbound", ignore.case = TRUE)
})

test_that("describe release rejects DSI failures without raw text or retry", {
  conns <- list(
    site_a = structure(1, class = "fake"),
    site_b = structure(2, class = "fake"))
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
      "DataSHIELD transport failed during 'synopsis connection identity fan-out'",
      captured, fixed = TRUE))
    expect_false(grepl("SECRET_REMOTE_DETAIL", captured, fixed = TRUE))
    expect_identical(attempts, 1L)
  }
})
