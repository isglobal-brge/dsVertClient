.dp_cor_pair_descriptor <- function(analysis_id, left, right, capacity, scale) {
  list(
    version = "pairwise-complete-normalized-moments-v1",
    analysis_id = analysis_id, dataset = "cohort", owner_peer = "site_a",
    left = list(column = left, lower = 0, upper = 1),
    right = list(column = right, lower = 0, upper = 1),
    numeric_grid_bits = 8L, coordinate_count = 6L,
    coordinate_order = paste0(
      "count,quantized_sum_left,quantized_sum_right,",
      "quantized_sumsq_left,quantized_sumsq_right,",
      "quantized_cross_product"),
    repeated_record_policy = paste0(
      "clip_finite_rows_then_mean_each_variable_once_per_",
      "admitted_unit_v1"),
    missingness_policy =
      "pairwise_complete_units_with_both_collapsed_values_v1",
    statistic_maximum = c(capacity, rep(capacity * scale, 5L)),
    source_raw_l1_sensitivity = 1 + 5 * scale,
    source_raw_l2_sensitivity = sqrt(1 + 5 * scale^2) *
      (1 + 64 * .Machine$double.eps),
    natural_l1_sensitivity = 6,
    natural_l2_sensitivity = sqrt(6) *
      (1 + 64 * .Machine$double.eps),
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste0(
      "zero_missing_pair_vs_complete_unit_is_worst_case_for_",
      "add_remove_and_replace_one"))
}

.dp_cor_fixture <- function(gaussian = FALSE) {
  scale <- 256
  capacity <- 100
  analysis_id <- "cohort::site_a"
  descriptors <- list(
    pair_xy = .dp_cor_pair_descriptor(
      analysis_id, "x", "y", capacity, scale),
    pair_xz = .dp_cor_pair_descriptor(
      analysis_id, "x", "z", capacity, scale),
    pair_yz = .dp_cor_pair_descriptor(
      analysis_id, "y", "z", capacity, scale))
  correlation <- list(
    version = "same-owner-pairwise-correlation-artifact-v1",
    analysis_id = analysis_id, dataset = "cohort", owner_peer = "site_a",
    variables = c("x", "y", "z"),
    pair_references = names(descriptors), pair_count = 3L,
    coordinate_count = 18L,
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
  count <- list(
    owner_peer = "site_a", dataset = "cohort",
    statistic_maximum = capacity, l1_sensitivity = 1)
  capsule_mechanism <- list(
    mechanism = if (isTRUE(gaussian)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM
    } else {
      "discrete-laplace"
    },
    sensitivity_norm = if (isTRUE(gaussian)) "l2" else "l1")
  manifest <- list(workload = list(
    coordinate_count = 19,
    release_lattice = list(
      version = "biomedical-capsule-common-lattice-v1",
      output_lattice_bits = 8L, output_lattice_scale = scale,
      natural_l1_sensitivity = 19,
      integer_l1_sensitivity_steps = 19 * scale,
      natural_l2_sensitivity = sqrt(19),
      integer_l2_sensitivity_steps = sqrt(19) * scale),
    capsule_mechanism = capsule_mechanism,
    families = list(
      admitted_count = count,
      numeric_moments = list(artifacts = list()),
      numeric_pair_moments = list(
        artifacts = descriptors, natural_l1_sensitivity = 18),
      gaussian_models = list(
        artifacts = list(), natural_l1_sensitivity = 0),
      fixed_numeric_histograms = list(artifacts = list()),
      categorical_marginals = list(artifacts = list()),
      categorical_pairs = list(sets = list()),
      correlation_artifacts = stats::setNames(list(correlation), analysis_id),
      describe_artifacts = list(), survival_artifacts = list())))
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  # Each pair has n=80, sums=(40,40), squared sums=(30,30). Cross
  # sums 28,12,28 give correlations 0.8,-0.8,0.8. This valid pairwise matrix
  # is deliberately indefinite, as pairwise-complete matrices can be.
  values <- c(
    80,
    80, 40, 40, 30, 30, 28,
    80, 40, 40, 30, 30, 12,
    80, 40, 40, 30, 30, 28)
  if (isTRUE(gaussian)) {
    sensitivity_steps <- format(
      sqrt(19) * scale, digits = 17L, scientific = TRUE, trim = TRUE)
    plan <- list(
      version = .DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION,
      mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
      sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER,
      request_binding_sha256 = digest::digest(
        "correlation-gaussian-request", "sha256", serialize = FALSE),
      total_coordinate_count = 19,
      maximum_chunk_coordinates = 19,
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
      simultaneous_95_abs = "32")
    plan <- .client_complete_gaussian_plan_v2(plan)
    manifest$workload$mechanism_selection <- list(
      gaussian_calibration_request = list(
        epsilon = "1e+00", delta = format(
          2^-100, digits = 17L, scientific = TRUE, trim = TRUE),
        l2_sensitivity_steps = sensitivity_steps,
        total_coordinate_count = 19),
      gaussian_plan = plan,
      gaussian_plan_sha256 = .dsvert_vector_hash(plan))
  }
  release <- list(
    capsule_id = strrep("a", 64L), manifest_sha256 = strrep("b", 64L),
    final_vector_root = strrep("c", 64L),
    coordinate_order_sha256 = layout$sha256,
    coordinate_count = 19L, values = values,
    mechanism = if (isTRUE(gaussian)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM
    } else {
      .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM
    },
    epsilon = 1, delta = 2^-100,
    implementation_delta =
      "1/1267650600228229401496703205376",
    delta_aggregation = "max_per_peer_not_sum", sticky_replay = TRUE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    manifest = manifest)
  class(release) <- c("dsvert_joint_dp_vector", "list")
  status <- list(
    site_a = list(
      policy = list(adjacency = "add_remove_patient"),
      noise_root = list(epoch = 1, key_id = "root-a")),
    site_b = list(
      policy = list(adjacency = "add_remove_patient"),
      noise_root = list(epoch = 1, key_id = "root-b")))
  conns <- list(
    site_a = structure(list(), class = "fake"),
    site_b = structure(list(), class = "fake"))
  list(
    analysis_id = analysis_id, manifest = manifest, layout = layout,
    release = release, status = status, conns = conns,
    run = list(
      release = release, layout = layout, status = status,
      manifest_bundle = list()))
}

.dp_cor_gaussian_artifact <- function(capacity = 100, scale = 256,
                                      intercept = TRUE) {
  predictors <- list(
    x = list(column = "x", lower = 0, upper = 1),
    z = list(column = "z", lower = 0, upper = 1))
  predictor_order <- names(predictors)
  design_terms <- c(if (isTRUE(intercept)) "(Intercept)" else character(),
                    predictor_order)
  q <- length(design_terms)
  coordinate_count <- q * (q + 1) / 2 + q + 2
  list(
    version = "bounded-normalized-gaussian-sufficient-statistics-v1",
    spec_version = "v1", analysis_id = "gaussian_cor",
    dataset = "cohort", owner_peer = "site_a",
    outcome = list(column = "y", lower = 0, upper = 1),
    predictors = predictors, predictor_order = predictor_order,
    intercept = intercept, design_terms = design_terms,
    numeric_grid_bits = 8L, coordinate_count = coordinate_count,
    coordinate_order = paste0(
      "n_then_xtx_upper_column_major_then_xty_design_order_then_yty_",
      "v1"),
    repeated_record_policy = paste0(
      "clip_finite_rows_then_mean_each_variable_once_per_admitted_unit_",
      "v1"),
    missingness_policy =
      "complete_case_across_outcome_and_all_predictors_v1",
    contribution_domain =
      "one_vector_of_normalized_monomials_in_closed_unit_interval_v1",
    count_gram_intercept_policy = paste0(
      "n_is_complete_case_count_and_moment_upper_bound_gram11_governs_",
      "the_solve_no_averaging_v1"),
    statistic_maximum = c(
      capacity, rep(capacity * scale, coordinate_count - 1L)),
    source_raw_l1_sensitivity = 1 + (coordinate_count - 1L) * scale,
    source_raw_l2_sensitivity =
      sqrt(1 + (coordinate_count - 1L) * scale^2) *
      (1 + 64 * .Machine$double.eps),
    natural_l1_sensitivity = coordinate_count,
    natural_l2_sensitivity = sqrt(coordinate_count) *
      (1 + 64 * .Machine$double.eps),
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste0(
      "zero_missing_complete_case_vs_all_one_complete_unit_is_worst_",
      "case_for_add_remove_and_replace_one"),
    regularization_policy =
      "none_in_release_explicit_client_postprocessing_only_v1",
    implementation_state = "same_owner_materialized",
    cross_owner_state = "reserved_not_materialized")
}

.dp_cor_gaussian_cross_artifact <- function(k = 3L, capacity = 100,
                                            scale = 256) {
  stopifnot(k %in% 2:3)
  predictor_owner <- if (k == 2L) "site_a" else "site_c"
  computation <- c("site_a", "site_b")
  coordinate_count <- 7L
  product_error <- 2 + 1 / (4 * scale)
  list(
    version = .DSVERT_CLIENT_DP_GAUSSIAN_CROSS_ARTIFACT_VERSION,
    spec_version = "v2", analysis_id = "cross_cor",
    dataset = "outcome_data", owner_peer = "site_b",
    outcome = list(
      column = "y", dataset = "outcome_data", owner_peer = "site_b",
      lower = 0, upper = 1),
    predictors = list(x = list(
      column = "x", dataset = "predictor_data",
      owner_peer = predictor_owner, lower = 0, upper = 1)),
    predictor_order = list("x"), input_variable_order = list("x", "y"),
    participating_peers = as.list(sort(
      c("site_b", predictor_owner), method = "radix")),
    computation_peers = as.list(computation), intercept = TRUE,
    design_terms = list("(Intercept)", "x"), numeric_grid_bits = 8,
    coordinate_count = coordinate_count,
    coordinate_order = paste0(
      "n_then_xtx_upper_column_major_then_xty_design_order_then_yty_v2"),
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    private_input_layout = paste0(
      "capacity_padded_value_then_validity_per_signed_variable_",
      "manifest_order_v1"),
    repeated_record_policy = paste0(
      "clip_finite_rows_then_mean_each_variable_once_per_admitted_unit_v1"),
    missingness_policy = paste0(
      "complete_case_mask_remains_secret_shared_through_joint_noise_v1"),
    contribution_domain = paste0(
      "round_normalized_inputs_then_exact_floor_ring128_products_on_",
      "closed_unit_interval_v1"),
    count_gram_intercept_policy = paste0(
      "n_and_all_moments_share_one_secret_complete_case_mask_and_are_",
      "released_only_after_joint_dp_v1"),
    statistic_maximum = as.list(rep(
      capacity * scale, coordinate_count)),
    source_raw_l1_sensitivity = coordinate_count * scale,
    source_raw_l2_sensitivity = sqrt(coordinate_count) * scale,
    natural_l1_sensitivity = coordinate_count,
    natural_l2_sensitivity = sqrt(coordinate_count),
    adjacency = "add_remove_patient",
    adjacency_sensitivity_basis = paste0(
      "zero_missing_complete_case_vs_all_one_complete_unit_is_worst_",
      "case_for_add_remove_and_replace_one"),
    quantization_contract = list(
      input_rounding = "nearest_integer_ties_to_even_r_v1",
      product_rounding = "exact_signed_floor_after_division_by_scale_v1",
      per_product_max_abs_error_lattice_steps = product_error,
      per_product_max_abs_error_normalized = product_error / scale,
      per_sum_max_abs_error_lattice_steps = capacity * product_error,
      same_owner_v1_numerically_identical = FALSE),
    numeric_certificate = list(
      version = "dsvert-cross-gaussian-numeric-certificate-v1",
      ring_bits = 128, frac_bits = 8, required_signed_bits = 26,
      operand_maximum = scale, raw_product_maximum = scale^2,
      accumulated_coordinate_maximum = capacity * scale,
      truncation = "exact_signed_floor_gc_ot_or_direct_wide_v1",
      comparison = "not_used_after_custodian_bound_clipping",
      modular_wrap_proved_absent = TRUE,
      overflow_behavior = "typed_abort_before_commit"),
    transcript = list(
      version = "dsvert-cross-gaussian-fixed-transcript-v1",
      padded_units = capacity, variable_count = 2,
      validity_product_rounds = 1, masked_value_rounds = 1,
      moment_product_rounds = 1, data_dependent_branches = 0,
      exact_intermediate_release_count = 0),
    alignment_contract = list(
      version = "private-psi-ordered-manifest-consensus-v1",
      public_patient_dependent_hash = FALSE,
      mismatch_behavior = "typed_non_prealigned_cohort_failure"),
    regularization_policy =
      "none_in_release_explicit_client_postprocessing_only_v1",
    implementation_state = "cross_owner_exact_gc_materialized",
    cross_owner_state = "exact_gc_to_joint_dp_vector_v1")
}

.dp_cor_complete_case_fixture <- function(k = 1L, intercept = TRUE) {
  stopifnot(k %in% 1:3)
  capacity <- 100
  scale <- 256
  artifact <- if (k == 1L) {
    .dp_cor_gaussian_artifact(capacity, scale, intercept = intercept)
  } else {
    .dp_cor_gaussian_cross_artifact(k, capacity, scale)
  }
  peers <- c("site_a", "site_b", "site_c")[seq_len(max(2L, k))]
  count <- list(
    owner_peer = "site_a", dataset = "count_data",
    statistic_maximum = capacity, l1_sensitivity = 1)
  manifest <- list(workload = list(
    coordinate_count = 1L + artifact$coordinate_count,
    release_lattice = list(
      version = "biomedical-capsule-common-lattice-v1",
      output_lattice_bits = 8L, output_lattice_scale = scale,
      natural_l1_sensitivity = 1 + artifact$coordinate_count,
      integer_l1_sensitivity_steps =
        (1 + artifact$coordinate_count) * scale,
      natural_l2_sensitivity = sqrt(1 + artifact$coordinate_count),
      integer_l2_sensitivity_steps =
        sqrt(1 + artifact$coordinate_count) * scale),
    families = list(
      admitted_count = count,
      numeric_moments = list(artifacts = list()),
      numeric_pair_moments = list(
        artifacts = list(), natural_l1_sensitivity = 0),
      gaussian_models = list(
        artifacts = stats::setNames(list(artifact), artifact$analysis_id),
        natural_l1_sensitivity = artifact$coordinate_count),
      fixed_numeric_histograms = list(artifacts = list()),
      categorical_marginals = list(artifacts = list()),
      categorical_pairs = list(sets = list()),
      correlation_artifacts = list(), describe_artifacts = list(),
      survival_artifacts = list())))
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  gaussian_values <- if (k == 1L && isTRUE(intercept)) {
    # Four complete units: x=(0,.25,.75,1), z=(0,.5,.5,1),
    # y=(.1,.4,.6,.9). Coordinates follow the signed Gaussian order.
    c(4, 4, 2, 1.625, 2, 1.5, 1.5, 2, 1.45, 1.4, 1.34)
  } else if (k == 1L) {
    c(4, 1.625, 1.5, 1.5, 1.45, 1.4, 1.34)
  } else {
    c(4, 4, 2, 1.625, 2, 1.45, 1.34)
  }
  release <- list(
    capsule_id = strrep("a", 64L), values = c(4, gaussian_values),
    coordinate_count = 1L + artifact$coordinate_count,
    mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
    epsilon = 1, delta = 2^-100)
  class(release) <- c("dsvert_joint_dp_vector", "list")
  conns <- stats::setNames(lapply(peers, function(peer) {
    structure(list(peer = peer), class = "fake")
  }), peers)
  context <- list(
    release = release, manifest = manifest, layout = layout,
    lattice = manifest$workload$release_lattice,
    adjacency = "add_remove_patient")
  metadata <- list(
    released = TRUE, source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    capsule_id = release$capsule_id, epsilon = release$epsilon,
    delta = release$delta, mechanism = release$mechanism,
    implementation = .DSVERT_CLIENT_VECTOR_BACKEND,
    sampler = .DSVERT_CLIENT_VECTOR_SAMPLER,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    sensitivity_norm = "l1", mechanism_selection = list())
  list(
    artifact = artifact, manifest = manifest, layout = layout,
    release = release, conns = conns, context = context,
    metadata = metadata, run = list(fixture = TRUE))
}

.with_dp_cor_complete_case_mocks <- function(fixture, code) {
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) {
      calls <<- calls + 1L
      fixture$run
    },
    .dsvert_dp_vector_context = function(run) fixture$context,
    .dsvert_dp_vector_public_metadata = function(context) fixture$metadata,
    .dsvert_dp_vector_accuracy_radius = function(...) list(
      radius = 0.01, confidence = 0.95, method = "fixture_exact_radius",
      implementation_tv_upper_bound = 0,
      sampler_tv_upper_bound = 0,
      additional_privacy_cost = c(epsilon = 0, delta = 0)),
    .package = "dsVertClient")
  result <- force(code)
  attr(result, "capsule_calls") <- calls
  result
}

test_that("DP correlation matches the central pairwise-moment oracle", {
  fixture <- .dp_cor_fixture()
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) {
      calls <<- calls + 1L
      fixture$run
    },
    .dsvert_aggregate_strict = function(...) {
      stop("legacy aggregate route called", call. = FALSE)
    },
    .dsvert_fanout_by_site = function(...) {
      stop("legacy fanout route called", call. = FALSE)
    },
    .package = "dsVertClient")

  result <- .dsvert_dp_cor_impl(
    "cohort", fixture$analysis_id, c("x", "y", "z"), "site_a",
    fixture$conns, function(...) stop("raw DSI call", call. = FALSE))
  expect_identical(calls, 1L)
  expect_equal(result$correlation_raw_pairwise, matrix(
    c(1, .8, -.8, .8, 1, .8, -.8, .8, 1), 3L, 3L,
    dimnames = list(c("x", "y", "z"), c("x", "y", "z"))),
    tolerance = 1e-12)
  expect_true(isTRUE(all.equal(result$correlation, t(result$correlation),
                               tolerance = 1e-12)))
  expect_equal(unname(diag(result$correlation)), rep(1, 3L),
               tolerance = 1e-12)
  expect_gte(min(eigen(result$correlation, symmetric = TRUE,
                       only.values = TRUE)$values), -1e-12)
  expect_lt(result$psd_projection$input_min_eigenvalue, 0)
  expect_true(result$psd_projection_applied)
  expect_false(result$psd_projection$exact_nearest_correlation_matrix)
  expect_true(result$psd_projection_changes_pairwise_estimand)
  expect_identical(result$cross_owner_state, "reserved_not_materialized")
  expect_identical(result$legacy_exact_route_called, FALSE)
  expect_identical(result$source_artifact_family, "correlation_artifacts")
  expect_identical(result$estimand_missingness, "pairwise_complete")
  expect_identical(result$pca_eligible, FALSE)
  expect_identical(result$source_values_exposed, FALSE)
  expect_identical(result$intermediate_values_exposed, FALSE)
  expect_false(any(c("data", "exact", "raw_rows") %in% names(result)))
  expect_true(all(
    result$correlation_95_interval_raw_pairwise$lower <=
      result$correlation_raw_pairwise + 1e-12))
  expect_true(all(
    result$correlation_95_interval_raw_pairwise$upper >=
      result$correlation_raw_pairwise - 1e-12))
  expect_match(
    result$projected_enclosure_semantics,
    "not an interval for the projection", fixed = TRUE)
  expect_identical(result$accuracy_additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
})

test_that("same capsule permits flexible signed subsets without a new release", {
  fixture <- .dp_cor_fixture()
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) {
      calls <<- calls + 1L
      fixture$run
    },
    .package = "dsVertClient")
  subset <- .dsvert_dp_cor_impl(
    "cohort", fixture$analysis_id, c("z", "x"), NULL,
    fixture$conns, function(...) NULL)
  expect_identical(calls, 1L)
  expect_identical(subset$var_names, c("z", "x"))
  expect_equal(subset$correlation_raw_pairwise[1, 2], -.8,
               tolerance = 1e-12)
  expect_identical(subset$accuracy_additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
})

test_that("explicit pairwise correlation preserves the signed L2 plan", {
  fixture <- .dp_cor_fixture(gaussian = TRUE)
  calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) {
      calls <<- calls + 1L
      fixture$run
    },
    .dsvert_aggregate_strict = function(...) {
      stop("legacy aggregate route called", call. = FALSE)
    },
    .package = "dsVertClient")

  result <- .dsvert_dp_cor_impl(
    "cohort", fixture$analysis_id, c("x", "y", "z"), NULL,
    fixture$conns, function(...) stop("raw DSI call", call. = FALSE))
  expect_identical(calls, 1L)
  expect_identical(result$mechanism,
                   .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
  expect_identical(result$implementation,
                   .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND)
  expect_identical(result$sampler, .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER)
  expect_identical(result$sensitivity_norm, "l2")
  expect_equal(result$accuracy_simultaneous_95_abs_raw_coordinates,
               32 / 256, tolerance = 0)
  expect_match(result$accuracy_simultaneous_method,
               "discrete-Gaussian plan", fixed = TRUE)
  expect_true(all(
    result$correlation_95_interval_raw_pairwise$lower <=
      result$correlation_raw_pairwise + 1e-12))
  expect_true(all(
    result$correlation_95_interval_raw_pairwise$upper >=
      result$correlation_raw_pairwise - 1e-12))

  expect_error(
    ds.vertPCA(cor_result = result, n_components = 2L, verbose = FALSE),
    "complete-case ds.vertCor")
})

test_that("complete-case correlation matches the central joint oracle", {
  fixture <- .dp_cor_complete_case_fixture()
  result <- .with_dp_cor_complete_case_mocks(
    fixture,
    .dsvert_dp_cor_gaussian_impl(
      "cohort", "gaussian_cor", c("x", "z", "y"), NULL,
      fixture$conns, function(...) stop("raw DSI call", call. = FALSE)))
  oracle <- stats::cor(cbind(
    x = c(0, .25, .75, 1), z = c(0, .5, .5, 1),
    y = c(.1, .4, .6, .9)))

  expect_identical(attr(result, "capsule_calls"), 1L)
  expect_equal(result$correlation_raw_complete_case, oracle,
               tolerance = 1e-12)
  expect_equal(result$correlation, oracle, tolerance = 1e-12)
  expect_identical(result$source_artifact_family, "gaussian_models")
  expect_identical(result$estimand_missingness, "complete_case_joint")
  expect_identical(result$pca_eligible, TRUE)
  expect_identical(result$correlation_raw_pairwise, NULL)
  expect_equal(result$complete_case_n, matrix(
    4, 3L, 3L,
    dimnames = list(c("x", "z", "y"), c("x", "z", "y"))))
  expect_identical(result$pairwise_n, NULL)
  expect_identical(result$additional_server_calls_after_capsule, 0L)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_false(result$inference$sampling_inference_available)
  expect_identical(result$inference$p_values, NULL)
  expect_true(all(
    result$correlation_95_interval_complete_case$lower <=
      result$correlation_raw_complete_case + 1e-12))
  expect_true(all(
    result$correlation_95_interval_complete_case$upper >=
      result$correlation_raw_complete_case - 1e-12))
})

test_that("correlation interval arithmetic rounds outward", {
  x <- c(0.1, 0.2, 0.8, 0.9)
  y <- c(0.15, 0.4, 0.7, 0.95)
  coordinates <- c(
    length(x), sum(x), sum(y), sum(x^2), sum(y^2), sum(x * y))
  interval <- .dsvert_dp_cor_interval(
    coordinates, radius = 0, capacity = 100, scale = 256,
    quantization = 0)
  oracle <- stats::cor(x, y)
  expect_lte(interval$correlation[["lower"]], oracle)
  expect_gte(interval$correlation[["upper"]], oracle)
})

test_that("complete-case correlation requires signed intercept sums", {
  fixture <- .dp_cor_complete_case_fixture(intercept = FALSE)
  condition <- tryCatch(
    .with_dp_cor_complete_case_mocks(
      fixture,
      .dsvert_dp_cor_gaussian_impl(
        "cohort", "gaussian_cor", c("x", "z", "y"), NULL,
        fixture$conns, function(...) NULL)),
    non_identifiable = function(condition) condition)
  expect_s3_class(condition, "non_identifiable")
  expect_match(condition$message, "without an intercept")
  expect_identical(
    condition$reason, "pearson_requires_intercept_marginal_sums")
})

test_that("cross-owner complete-case Cor supports K2 and K3 role layouts", {
  for (k in 2:3) {
    fixture <- .dp_cor_complete_case_fixture(k)
    owners <- if (k == 2L) {
      list(site_a = "x", site_b = "y")
    } else {
      list(site_c = "x", site_b = "y")
    }
    result <- .with_dp_cor_complete_case_mocks(
      fixture,
      .dsvert_dp_cor_gaussian_impl(
        "outcome_data", "cross_cor", owners, NULL,
        fixture$conns, function(...) stop("raw DSI call", call. = FALSE)))
    expect_identical(attr(result, "capsule_calls"), 1L)
    expect_identical(result$participating_servers,
                     unlist(fixture$artifact$participating_peers,
                            use.names = FALSE))
    expect_identical(result$computation_servers, c("site_a", "site_b"))
    expect_identical(result$cross_owner_state,
                     "exact_gc_to_joint_dp_vector_v1")
    expect_identical(result$var_names, c("x", "y"))
    expect_identical(result$additional_server_calls_after_capsule, 0L)
    if (k == 3L) {
      expect_false("site_c" %in% result$computation_servers)
      expect_true("site_c" %in% result$participating_servers)
    }
  }
})

test_that("cross-owner requests fail before any capsule or legacy call", {
  fixture <- .dp_cor_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) {
      stop("capsule should not run", call. = FALSE)
    },
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_cor_impl(
      "cohort", fixture$analysis_id,
      list(site_a = "x", site_b = "z"), NULL,
      fixture$conns, function(...) NULL),
    "reserved_not_materialized")
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) fixture$run,
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_cor_impl(
      "cohort", fixture$analysis_id, c("x", "remote_z"), NULL,
      fixture$conns, function(...) NULL),
    "reserved_not_materialized")
  expect_error(
    ds.vertCor("cohort", c("x", "y"), analysis_id = NULL,
               verbose = FALSE, datasources = fixture$conns),
    "analysis_id")
})

test_that("ds.vertCor front door uses Gaussian and never pairwise fallback", {
  fixture <- .dp_cor_complete_case_fixture()
  result <- .with_dp_cor_complete_case_mocks(
    fixture,
    ds.vertCor(
      "cohort", list(site_a = c("x", "y")),
      analysis_id = "gaussian_cor", verbose = FALSE,
      datasources = fixture$conns))
  expect_identical(attr(result, "capsule_calls"), 1L)
  expect_identical(result$source_artifact_family, "gaussian_models")
  expect_identical(result$var_names, c("x", "y"))
})

test_that("PCA is zero-cost post-processing with honest eigengap diagnostics", {
  fixture <- .dp_cor_complete_case_fixture()
  correlation <- .with_dp_cor_complete_case_mocks(
    fixture,
    .dsvert_dp_cor_gaussian_impl(
      "cohort", "gaussian_cor", c("x", "z", "y"), NULL,
      fixture$conns, function(...) NULL))
  pca <- ds.vertPCA(
    cor_result = correlation, n_components = 2L, verbose = FALSE)
  expect_s3_class(pca, "ds.pca")
  expect_identical(dim(pca$loadings), c(3L, 2L))
  expect_identical(pca$scores, NULL)
  expect_false(pca$scores_available)
  expect_identical(pca$additional_server_calls, 0L)
  expect_identical(pca$additional_server_calls_after_capsule, 0L)
  expect_identical(pca$capsule_read_performed, FALSE)
  expect_identical(pca$capsule_workflow_count, 0L)
  expect_identical(pca$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_true(pca$history_gate)
  expect_false(pca$request_limit)
  expect_true(pca$operation_limit)
  expect_true(all(pca$eigenvalue_95_mechanism_regions[, "lower"] <=
                    pca$eigenvalues))
  expect_true(all(pca$eigenvalue_95_mechanism_regions[, "upper"] >=
                    pca$eigenvalues))
  enclosure <- correlation[[
    "correlation_95_enclosure_raw_estimand_around_projected_release"]]
  naive_spectral <- sqrt(sum(pmax(
    abs(correlation$correlation - enclosure$lower),
    abs(enclosure$upper - correlation$correlation))^2))
  expect_gte(pca$spectral_mechanism_radius_95, naive_spectral)
  expect_match(
    pca$uncertainty_scope, "raw complete-case estimand", fixed = TRUE)
  expect_identical(pca$source_artifact_family, "gaussian_models")
  expect_identical(pca$estimand_missingness, "complete_case_joint")
  expect_identical(pca$correlation_raw_pairwise, NULL)
  expect_true(pca$loading_identifiability$sign_indeterminate)
  expect_named(pca$loading_identifiability,
               c("sign_indeterminate", "adjacent_observed_gaps",
                 "adjacent_gap_95_lower", "adjacent_gap_95_upper",
                 "separated_within_mechanism_region",
                 "davis_kahan_angle_upper_radians",
                 "individual_directions_identifiable",
                 "tied_or_unresolved_eigenspaces_identifiable"))
  expect_identical(pca$legacy_exact_route_called, FALSE)

  for (field in c("history_gate", "request_limit", "operation_limit")) {
    tampered <- correlation
    tampered[[field]] <- !tampered[[field]]
    expect_error(
      ds.vertPCA(cor_result = tampered, verbose = FALSE),
      "authenticated complete-case ds.vertCor")
  }

  legacy <- structure(list(correlation = diag(2), var_names = c("a", "b")),
                      class = c("ds.cor", "list"))
  expect_error(
    ds.vertPCA(cor_result = legacy, verbose = FALSE),
    "authenticated complete-case ds.vertCor")

  pairwise <- .dp_cor_fixture()
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) pairwise$run,
    .package = "dsVertClient")
  pairwise_result <- .dsvert_dp_cor_impl(
    "cohort", pairwise$analysis_id, c("x", "y"), NULL,
    pairwise$conns, function(...) NULL)
  expect_error(
    ds.vertPCA(cor_result = pairwise_result, verbose = FALSE),
    "pairwise")
})

test_that("direct PCA front door performs exactly one capsule read", {
  fixture <- .dp_cor_complete_case_fixture()
  pca <- .with_dp_cor_complete_case_mocks(
    fixture,
    ds.vertPCA(
      data_name = "cohort", variables = c("x", "z", "y"),
      n_components = 2L, analysis_id = "gaussian_cor",
      verbose = FALSE, datasources = fixture$conns))
  expect_identical(attr(pca, "capsule_calls"), 1L)
  expect_identical(pca$additional_server_calls, 0L)
  expect_identical(pca$additional_server_calls_after_capsule, 0L)
  expect_identical(pca$capsule_read_performed, TRUE)
  expect_identical(pca$capsule_workflow_count, 1L)
  expect_identical(pca$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_identical(pca$scores, NULL)
})

test_that("pairwise denominator validation never invents a correlation", {
  fixture <- .dp_cor_fixture()
  fixture$release$values[5L] <- 20
  fixture$release$values[6L] <- 20
  fixture$run$release <- fixture$release
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) fixture$run,
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_cor_impl(
      "cohort", fixture$analysis_id, c("x", "y"), NULL,
      fixture$conns, function(...) NULL),
    "non_identifiable")
})

test_that("correlation rejects altered signed pair semantics", {
  fixture <- .dp_cor_fixture()
  fixture$manifest$workload$families$numeric_pair_moments$
    artifacts$pair_xy$missingness_policy <- "complete_case_all_variables"
  fixture$release$manifest <- fixture$manifest
  fixture$run$release <- fixture$release
  fixture$run$layout <- .dsvert_dp_capsule_vector_layout(fixture$manifest)
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) fixture$run,
    .package = "dsVertClient")
  expect_error(
    .dsvert_dp_cor_impl(
      "cohort", fixture$analysis_id, c("x", "y"), NULL,
      fixture$conns, function(...) NULL),
    "numeric-pair descriptor is invalid")
})

test_that("correlation and PCA front doors cannot reach legacy endpoints", {
  namespace <- asNamespace("dsVertClient")
  reachable <- character()
  queue <- c("ds.vertCor", "ds.vertDPCor", "ds.vertPCA")
  while (length(queue)) {
    name <- queue[[1L]]
    queue <- queue[-1L]
    if (name %in% reachable ||
        !exists(name, namespace, mode = "function", inherits = FALSE)) next
    value <- get(name, namespace, inherits = FALSE)
    reachable <- c(reachable, name)
    globals <- tryCatch(
      unique(unlist(codetools::findGlobals(value, merge = FALSE),
                    use.names = FALSE)),
      error = function(error) character())
    queue <- unique(c(queue, intersect(
      globals, ls(namespace, all.names = TRUE))))
  }
  forbidden <- c(
    "dsvertColNamesDS", "getObsCountDS", "glmStandardizeDS", "localCorDS",
    "k2ShareInputDS", "k2GradientR1DS", "k2GradientR2DS",
    "glmRing63CorSetColDS", "glmRing63CorSetZeroYDS")
  bodies <- paste(vapply(reachable, function(name) {
    paste(deparse(body(get(name, namespace, inherits = FALSE))),
          collapse = "\n")
  }, character(1L)), collapse = "\n")
  expect_contains(reachable, ".dsvert_dp_capsule_vector_run")
  expect_length(intersect(reachable, forbidden), 0L)
  for (endpoint in forbidden) {
    expect_false(grepl(endpoint, bodies, fixed = TRUE), info = endpoint)
  }

  source_files <- list.files(
    testthat::test_path("..", "..", "R"), pattern = "[.]R$",
    full.names = TRUE)
  source_text <- paste(vapply(source_files, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1L)), collapse = "\n")
  expect_false(grepl("glmRing63CorSetColDS", source_text, fixed = TRUE))
  expect_false(grepl("glmRing63CorSetZeroYDS", source_text, fixed = TRUE))
})
