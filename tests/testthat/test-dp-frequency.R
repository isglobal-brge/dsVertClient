.dp_frequency_fixture <- function(gaussian = FALSE) {
  scale <- 2^8
  marginal <- list(
    dataset = "cohort", column = "status", owner_peer = "site_a",
    levels = as.list(c("none", "mild", "severe")),
    repeated_record_policy = "consistent_level_else_exclude_v1",
    missingness_policy = "missing_or_out_of_domain_rows_are_ignored",
    statistic_maximum = 100)
  mechanism_plan <- if (isTRUE(gaussian)) .client_complete_gaussian_plan_v2(list(
    version = .DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION,
    mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
    sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER,
    request_binding_sha256 = digest::digest(
      "frequency-gaussian-request", "sha256", serialize = FALSE),
    total_coordinate_count = 4,
    maximum_chunk_coordinates = 4,
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
      "1267650600228229401496703205376",
    vector_sampler_tv_upper_numerator = "1",
    vector_sampler_tv_upper_denominator =
      "1267650600228229401496703205376",
    vector_total_tv_upper_numerator = "2",
    vector_total_tv_upper_denominator =
      "1267650600228229401496703205376",
    observable_worker_shape = "fixed dyadic CDF fixture",
    per_peer_implementation_delta_numerator = "1",
    per_peer_implementation_delta_denominator =
      "1267650600228229401496703205376",
    simultaneous_95_abs = "25")) else NULL
  capsule_mechanism <- if (isTRUE(gaussian)) list(
    mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
    sensitivity_norm = "l2",
    certificate = list(
      gaussian_calibration_request = list(
        epsilon = "1e+00",
        delta = format(2^-100, digits = 17L,
                       scientific = TRUE, trim = TRUE),
        l2_sensitivity_steps = format(
          sqrt(2) * scale, digits = 17L,
          scientific = TRUE, trim = TRUE),
        total_coordinate_count = 4),
      gaussian_plan = mechanism_plan,
      gaussian_plan_sha256 = .dsvert_vector_hash(mechanism_plan))) else list(
        mechanism = "discrete-laplace", sensitivity_norm = "l1")
  manifest <- list(workload = list(
    coordinate_count = 4,
    release_lattice = list(
      version = "biomedical-capsule-common-lattice-v1",
      output_lattice_bits = 8L, output_lattice_scale = scale,
      natural_l1_sensitivity = 2,
      integer_l1_sensitivity_steps = 2 * scale,
      natural_l2_sensitivity = sqrt(2),
      integer_l2_sensitivity_steps = sqrt(2) * scale),
    capsule_mechanism = capsule_mechanism,
    families = list(
      admitted_count = list(
        version = "admitted-unit-count-v2", owner_peer = "site_a",
        dataset = "cohort", statistic_minimum = 0,
        statistic_maximum = 100, l1_sensitivity = 1,
        l2_sensitivity = 1),
      numeric_moments = list(artifacts = list()),
      numeric_pair_moments = list(artifacts = list()),
      gaussian_models = list(artifacts = list()),
      fixed_numeric_histograms = list(artifacts = list()),
      categorical_marginals = list(artifacts = list(status = marginal)),
      categorical_pairs = list(sets = list()),
      correlation_artifacts = list(), describe_artifacts = list(),
      survival_artifacts = list())))
  if (isTRUE(gaussian)) {
    manifest$workload$mechanism_selection <-
      capsule_mechanism$certificate
  }
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  release <- list(
    version = "dsvert-joint-dp-biomedical-vector-client-v1",
    capsule_id = strrep("1", 64L), manifest_sha256 = strrep("2", 64L),
    final_vector_root = strrep("3", 64L),
    coordinate_order_sha256 = layout$sha256, coordinate_count = 4L,
    output_lattice_bits = 8L, output_lattice_scale = scale,
    values = c(42, 11.25, 20.5, 7.75),
    mechanism = if (isTRUE(gaussian)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM
    } else {
      .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM
    },
    mechanism_plan = mechanism_plan,
    plan_sha256 = if (isTRUE(gaussian)) {
      .dsvert_vector_hash(mechanism_plan)
    } else NULL,
    epsilon = 1, delta = 2^-100,
    implementation_delta = if (isTRUE(gaussian)) {
      "1/1267650600228229401496703205376"
    } else "0/1",
    delta_aggregation = if (isTRUE(gaussian)) {
      "max_per_peer_not_sum"
    } else "zero",
    sticky_replay = TRUE, source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    manifest = manifest)
  class(release) <- c("dsvert_joint_dp_vector", "list")
  status <- stats::setNames(lapply(c("site_a", "site_b"), function(peer) {
    list(
      policy = list(adjacency = "add_remove_patient", peer_count = 2L),
      noise_root = list(epoch = 2, key_id = paste0("noise-", peer)))
  }), c("site_a", "site_b"))
  list(
    run = list(release = release, layout = layout, status = status,
               manifest_bundle = list()),
    conns = list(site_a = structure(1, class = "fake"),
                 site_b = structure(2, class = "fake")))
}

test_that("DP frequency reads one signed marginal with one capsule call", {
  fixture <- .dp_frequency_fixture()
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) {
      calls <<- calls + 1L
      fixture$run
    },
    .package = "dsVertClient")
  result <- .dsvert_dp_frequency_impl(
    "cohort", "status", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))

  expect_s3_class(result, "ds.vertDPFrequency")
  expect_identical(calls, 1L)
  expect_identical(result$server, "site_a")
  expect_identical(result$levels, c("none", "mild", "severe"))
  expect_equal(result$counts, c(none = 11.25, mild = 20.5, severe = 7.75))
  expect_equal(sum(result$proportions), 1)
  expect_gt(result$accuracy_simultaneous_95_abs, 0)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_false(result$request_limit)
})

test_that("DP frequency validates the finite-support Gaussian certificate", {
  fixture <- .dp_frequency_fixture(gaussian = TRUE)
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) fixture$run,
    .package = "dsVertClient")
  result <- .dsvert_dp_frequency_impl(
    "cohort", "status", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_s3_class(result, "ds.vertDPFrequency")
  expect_match(result$mechanism, "gaussian", ignore.case = TRUE)
  expect_match(result$accuracy_simultaneous_method,
               "fixed-work dyadic discrete-Gaussian")
  expect_gt(result$accuracy_implementation_tv_upper_bound, 0)
})

test_that("frequency mechanism regions contain every feasible count vector", {
  counts <- c(a = 3, b = 5, c = 2)
  regions <- .dsvert_dp_frequency_regions(counts, radius = 2, capacity = 10)
  grids <- expand.grid(a = 1:5, b = 3:7, c = 0:4)
  grids <- grids[rowSums(grids) <= 10 & rowSums(grids) > 0, , drop = FALSE]
  for (index in seq_len(nrow(grids))) {
    value <- unlist(grids[index, ], use.names = FALSE)
    proportion <- value / sum(value)
    expect_true(all(proportion >= regions$proportion[, "lower"] - 1e-15))
    expect_true(all(proportion <= regions$proportion[, "upper"] + 1e-15))
  }
})

test_that("frequency inference is zero-call and covers every CP box corner", {
  fixture <- .dp_frequency_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) fixture$run,
    .package = "dsVertClient")
  frequency <- .dsvert_dp_frequency_impl(
    "cohort", "status", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  inference <- ds.vertDPFrequencyInference(
    frequency, level = 0.95, dp_fraction = 0.4)

  expect_s3_class(inference, "ds.vertDPFrequencyInference")
  expect_identical(inference$additional_server_calls, 0L)
  expect_identical(inference$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_equal(inference$coverage_lower_bound, 0.95)
  expect_true(all(inference$intervals[, "lower"] >= 0))
  expect_true(all(inference$intervals[, "upper"] <= 1))
  expect_true(all(inference$intervals[, "lower"] <=
                  inference$intervals[, "upper"]))

  tampered <- frequency
  tampered$counts[[1L]] <- tampered$coordinate_maximum + 1
  expect_error(
    ds.vertDPFrequencyInference(tampered),
    "released, validated ds.vertDPFrequency")
})

test_that("frequency CP envelopes exhaust every small integer box", {
  alpha <- 0.01
  for (event_lower in 0:3) {
    for (event_upper in event_lower:4) {
      for (failure_lower in 0:3) {
        for (failure_upper in failure_lower:4) {
          envelope <- .dsvert_dp_cp_union_over_box(
            event_lower, event_upper, failure_lower, failure_upper, alpha)
          for (events in event_lower:event_upper) {
            for (failures in failure_lower:failure_upper) {
              interval <- .dsvert_dp_cp_interval(events, failures, alpha)
              expect_lte(envelope[["lower"]], interval[["lower"]])
              expect_gte(envelope[["upper"]], interval[["upper"]])
            }
          }
        }
      }
    }
  }
})

test_that("frequency rejects absent, ambiguous, or wrong-owner artifacts", {
  fixture <- .dp_frequency_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) fixture$run,
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_frequency_impl(
      "cohort", "unknown", NULL, fixture$conns, function(...) NULL),
    "exactly one|signed categorical")
  expect_error(
    .dsvert_dp_frequency_impl(
      "cohort", "status", "site_b", fixture$conns, function(...) NULL),
    "exactly one|signed categorical")
})
