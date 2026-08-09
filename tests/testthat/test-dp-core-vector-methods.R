.dp_core_vector_fixture <- function(gaussian = FALSE) {
  scale <- 2^8
  count <- list(
    version = "admitted-unit-count-v2", owner_peer = "site_b",
    dataset = "cohort", statistic_minimum = 0,
    statistic_maximum = 100, l1_sensitivity = 1,
    l2_sensitivity = 1)
  numeric <- list(artifacts = list(age = list(
    owner_peer = "site_b", dataset = "cohort", column = "age",
    lower = -10, upper = 30, numeric_grid_bits = 8L,
    statistic_maximum = as.list(c(100, 100 * scale, 100 * scale)),
    repeated_record_policy =
      "clip_finite_rows_then_mean_once_per_admitted_unit_v1",
    missingness_policy =
      "NA_NaN_Inf_have_no_numeric_value_finite_out_of_bounds_clip")))
  pairs <- list(sets = list(`cohort::site_b` = list(
    owner_peer = "site_b", dataset = "cohort",
    columns = list(
      list(column = "exposure", levels = as.list(c("no", "yes"))),
      list(column = "disease", levels = as.list(c("no", "yes")))),
    included_pairs = list(c("disease", "exposure")),
    repeated_record_policy = "consistent_joint_cell_else_exclude_v1",
    missingness_policy = "missing_or_out_of_domain_rows_are_ignored",
    coordinate_count = 4L, pair_count = 1L,
    statistic_maximum = 100)))
  manifest <- list(
    workload = list(
      coordinate_count = 8,
      release_lattice = list(
        version = "biomedical-capsule-common-lattice-v1",
        output_lattice_bits = 8L, output_lattice_scale = scale,
        natural_l1_sensitivity = 5,
        integer_l1_sensitivity_steps = 5 * scale,
        natural_l2_sensitivity = sqrt(5),
        integer_l2_sensitivity_steps = sqrt(5) * scale),
      capsule_mechanism = list(
        mechanism = if (isTRUE(gaussian)) {
          .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM
        } else {
          "discrete-laplace"
        },
        sensitivity_norm = if (isTRUE(gaussian)) "l2" else "l1"),
      families = list(
        admitted_count = count,
        numeric_moments = numeric,
        numeric_pair_moments = list(artifacts = list()),
        gaussian_models = list(artifacts = list()),
        fixed_numeric_histograms = list(artifacts = list()),
        categorical_marginals = list(artifacts = list()),
        categorical_pairs = pairs, correlation_artifacts = list(),
        describe_artifacts = list(), survival_artifacts = list())))
  if (isTRUE(gaussian)) {
    sensitivity_steps <- format(
      sqrt(5) * scale, digits = 17L, scientific = TRUE, trim = TRUE)
    plan <- list(
      version = .DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION,
      mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
      sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER,
      request_binding_sha256 = digest::digest(
        "core-vector-gaussian-request", "sha256", serialize = FALSE),
      total_coordinate_count = 8,
      maximum_chunk_coordinates = 8,
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
      simultaneous_95_abs = "25")
    plan <- .client_complete_gaussian_plan_v2(plan)
    manifest$workload$mechanism_selection <- list(
      gaussian_calibration_request = list(
        epsilon = "1e+00", delta = format(
          2^-100, digits = 17L, scientific = TRUE, trim = TRUE),
        l2_sensitivity_steps = sensitivity_steps,
        total_coordinate_count = 8),
      gaussian_plan = plan,
      gaussian_plan_sha256 = .dsvert_vector_hash(plan))
  }
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  release <- list(
    version = "dsvert-joint-dp-biomedical-vector-client-v1",
    capsule_id = paste(rep("1", 64L), collapse = ""),
    manifest_sha256 = paste(rep("2", 64L), collapse = ""),
    final_vector_root = paste(rep("3", 64L), collapse = ""),
    coordinate_order_sha256 = layout$sha256,
    coordinate_count = 8L,
    output_lattice_bits = 8L, output_lattice_scale = scale,
    # count; n, normalized sum, normalized sumsq; canonical pair cells.
    values = c(42.25, 4, 1.5, 0.75, 8.5, 2.25, 3.75, 9.5),
    mechanism = if (isTRUE(gaussian)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM
    } else {
      .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM
    },
    epsilon = 1, delta = 2^-100,
    implementation_delta =
      "1/1267650600228229401496703205376",
    delta_aggregation = "max_per_peer_not_sum",
    sticky_replay = TRUE, source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    manifest = manifest)
  class(release) <- c("dsvert_joint_dp_vector", "list")
  status <- stats::setNames(lapply(c("site_a", "site_b"), function(peer) {
    list(
      policy = list(adjacency = "add_remove_patient", peer_count = 2L),
      noise_root = list(epoch = 4, key_id = paste0("noise-", peer)))
  }), c("site_a", "site_b"))
  list(
    run = list(release = release, layout = layout, status = status,
               manifest_bundle = list()),
    status = status,
    conns = list(site_a = structure(1, class = "fake"),
                 site_b = structure(2, class = "fake")))
}

.dp_core_vector_mock_runner <- function(fixture, counter) {
  force(fixture)
  force(counter)
  function(datasources, status = NULL, .aggregate) {
    counter$calls <- counter$calls + 1L
    if (!is.null(status)) {
      expect_identical(status, fixture$status)
    }
    fixture$run
  }
}

test_that("DP count reads only the signed vector count coordinate", {
  fixture <- .dp_core_vector_fixture()
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_joint_dp_capsule_status_impl = function(...) fixture$status,
    .dsvert_dp_capsule_vector_run =
      .dp_core_vector_mock_runner(fixture, counter),
    .package = "dsVertClient")

  result <- .dsvert_dp_count_impl(
    "cohort", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_s3_class(result, "ds.vertDPCount")
  expect_identical(counter$calls, 1L)
  expect_identical(result$value, 42.25)
  expect_identical(result$server, "site_b")
  expect_identical(result$mechanism,
    "two-independent-complete-vector-discrete-laplace-draws-v3")
  expect_identical(result$implementation,
                   .DSVERT_CLIENT_VECTOR_BACKEND)
  expect_identical(result$backend, "exact_signed_Ring128_global_vector")
  expect_true(result$history_gate)
  expect_false(result$request_limit)
  expect_true(result$operation_limit)
  expect_identical(result$security_claim, list(
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
    unconditional_non_reconstruction_guarantee = FALSE))
  expect_gt(result$accuracy_95_abs, 0)
  expect_lte(result$accuracy_95_abs, 100)
  expect_match(result$accuracy_95_method, "convolution")
})

test_that("DP contingency respects signed column-major orientation", {
  fixture <- .dp_core_vector_fixture()
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run =
      .dp_core_vector_mock_runner(fixture, counter),
    .package = "dsVertClient")

  result <- .dsvert_dp_contingency_impl(
    "cohort", "exposure", "disease", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_s3_class(result, "ds.vertDPContingency")
  expect_identical(counter$calls, 1L)
  expect_equal(result$table, matrix(
    c(8.5, 3.75, 2.25, 9.5), 2L, 2L,
    dimnames = list(c("no", "yes"), c("no", "yes"))))
  expect_identical(result$counts, c(8.5, 3.75, 2.25, 9.5))
  expect_identical(result$row_var, "exposure")
  expect_identical(result$col_var, "disease")
  expect_identical(
    result$unit_aggregation_policy,
    "consistent_joint_cell_else_exclude_v1")
  expect_gt(result$accuracy_simultaneous_95_abs,
            result$accuracy_95_abs_per_cell[[1L]])
  expect_identical(result$accuracy_additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
})

test_that("vector statistics support fixed-cohort replacement adjacency", {
  fixture <- .dp_core_vector_fixture()
  for (peer in names(fixture$status)) {
    fixture$status[[peer]]$policy$adjacency <-
      "replace_one_fixed_cohort"
  }
  fixture$run$status <- fixture$status
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) fixture$run,
    .package = "dsVertClient")

  table <- .dsvert_dp_contingency_impl(
    "cohort", "exposure", "disease", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  moments <- .dsvert_dp_meanvar_impl(
    "cohort", "age", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_identical(table$adjacency, "replace_one_fixed_cohort")
  expect_identical(table$artifact_l1_sensitivity, 2)
  expect_equal(table$artifact_l2_sensitivity, sqrt(2), tolerance = 0)
  expect_identical(moments$adjacency, "replace_one_fixed_cohort")
  expect_identical(moments$artifact_l1_sensitivity, 3)
})

test_that("DP mean/variance converts normalized moments to natural scale", {
  fixture <- .dp_core_vector_fixture()
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run =
      .dp_core_vector_mock_runner(fixture, counter),
    .package = "dsVertClient")

  result <- .dsvert_dp_meanvar_impl(
    "cohort", "age", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_s3_class(result, "ds.vertDPMeanVar")
  expect_identical(counter$calls, 1L)
  expect_identical(result$n, 4)
  expect_equal(result$mean, 5, tolerance = 0)
  expect_equal(result$variance, 75, tolerance = 0)
  expect_equal(result$sum, 20, tolerance = 0)
  expect_equal(result$sumsq, 400, tolerance = 0)
  expect_identical(result$normalized_sum_dp, 1.5)
  expect_identical(result$normalized_sumsq_dp, 0.75)
  expect_identical(result$submechanism_count, 1L)
  expect_false(result$noise_selection$coordinate_epsilon_split)
  expect_identical(
    names(result$mechanism_regions),
    c("effective_count", "mean", "variance"))
  expect_match(result$mechanism_region_scope,
               "sampling uncertainty excluded")
  expect_identical(result$mechanism_region_additional_server_calls, 0L)
})

test_that("core methods preserve one signed Gaussian L2 release", {
  fixture <- .dp_core_vector_fixture(gaussian = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_joint_dp_capsule_status_impl = function(...) fixture$status,
    .dsvert_dp_capsule_vector_run =
      .dp_core_vector_mock_runner(fixture, counter),
    .package = "dsVertClient")

  count <- .dsvert_dp_count_impl(
    "cohort", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  table <- .dsvert_dp_contingency_impl(
    "cohort", "exposure", "disease", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  moments <- .dsvert_dp_meanvar_impl(
    "cohort", "age", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_identical(counter$calls, 3L)
  for (result in list(count, table, moments)) {
    expect_identical(result$mechanism,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
    expect_identical(result$implementation,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND)
    expect_identical(result$sampler,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER)
    expect_identical(result$sensitivity_norm, "l2")
    expect_identical(result$mechanism_selection,
                     fixture$run$release$manifest$workload$mechanism_selection)
  }
  expect_equal(count$accuracy_95_abs, 25 / 256, tolerance = 0)
  expect_equal(table$accuracy_simultaneous_95_abs,
               25 / 256, tolerance = 0)
  expect_equal(moments$accuracy_95_abs_count, 25 / 256, tolerance = 0)
  expect_identical(moments$noise_selection$sensitivity_norm, "l2")
  expect_match(moments$mechanism_region_method,
               "discrete-Gaussian plan", fixed = TRUE)
  expect_identical(
    moments$accuracy_simultaneous_95_abs_raw_coordinates,
    c(count = 25 / 256, normalized_sum = 25 / 256,
      normalized_sumsq = 25 / 256))
})

test_that("core vector methods reject missing and ambiguous signed blocks", {
  fixture <- .dp_core_vector_fixture()
  runner <- function(datasources, status = NULL, .aggregate) fixture$run
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = runner,
    .package = "dsVertClient")
  expect_error(.dsvert_dp_meanvar_impl(
    "cohort", "missing", NULL, fixture$conns, function(...) NULL),
    "exactly one signed numeric-moment block")
  expect_error(.dsvert_dp_contingency_impl(
    "cohort", "exposure", "missing", NULL, fixture$conns,
    function(...) NULL), "exactly one signed categorical-pair block")
  expect_error(.dsvert_dp_contingency_impl(
    "cohort", "exposure", "disease", "site_a", fixture$conns,
    function(...) NULL), "exactly one signed categorical-pair block")

  ambiguous <- fixture
  duplicate <- ambiguous$run$layout$blocks[[
    "numeric_moments::age"]]
  ambiguous$run$layout$blocks[["numeric_moments::age_alias"]] <- duplicate
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) ambiguous$run,
    .package = "dsVertClient")
  expect_error(.dsvert_dp_meanvar_impl(
    "cohort", "age", NULL, fixture$conns, function(...) NULL),
    "exactly one signed numeric-moment block")
})

test_that("fixed-cohort count remains only a zero-sensitivity public value", {
  fixture <- .dp_core_vector_fixture()
  fixed_status <- fixture$status
  for (peer in names(fixed_status)) {
    fixed_status[[peer]]$policy$adjacency <- "replace_one_fixed_cohort"
    fixed_status[[peer]]$policy$unit_capacity <- 100L
    fixed_status[[peer]]$policy$peer_count <- 2L
  }
  public <- list(
    released = TRUE, value = 100,
    mechanism = "public_fixed_cohort_size_v1",
    implementation = "custodian_owned_policy_constant",
    sampler = "none", randomness = "none", sensitivity = 0,
    postprocessing = "none_public_policy_value",
    clipped_coordinates = 0L, accuracy_95_abs = 0,
    data_dependency = "none_public_fixed_cohort_policy",
    epsilon = 0, delta = 0, adjacency = "replace_one_fixed_cohort",
    peer_count = 2L)
  aggregate <- function(conns, expr, error = NULL,
                        errors.print = TRUE, ...) {
    expect_false(errors.print)
    expect_identical(as.character(expr[[1L]]),
                     "dsvertPublicFixedCohortCountDS")
    stats::setNames(rep(list(public), length(conns)), names(conns))
  }
  testthat::local_mocked_bindings(
    .dsvert_joint_dp_capsule_status_impl = function(...) fixed_status,
    .dsvert_dp_capsule_vector_run = function(...) {
      stop("fixed public Count must not invoke the DP vector", call. = FALSE)
    },
    .package = "dsVertClient")

  result <- .dsvert_dp_count_impl(
    "cohort", "site_b", fixture$conns, aggregate)
  expect_identical(result$value, 100)
  expect_identical(result$sensitivity, 0)
  expect_identical(result$epsilon, 0)
  expect_identical(result$delta, 0)

  forged <- public
  forged$sensitivity <- 1
  expect_error(.dsvert_dp_count_impl(
    "cohort", "site_b", fixture$conns,
    function(conns, expr, ...) {
      stats::setNames(rep(list(forged), length(conns)), names(conns))
    }),
    "invalid public fixed-cohort Count")

  split <- stats::setNames(rep(list(public), 2L), names(fixture$conns))
  split$site_a$value <- 99
  expect_error(.dsvert_dp_count_impl(
    "cohort", "site_b", fixture$conns,
    function(...) split), "invalid public fixed-cohort Count")

  wrong_consensus <- public
  wrong_consensus$value <- 99
  expect_error(.dsvert_dp_count_impl(
    "cohort", "site_b", fixture$conns,
    function(conns, expr, ...) {
      stats::setNames(
        rep(list(wrong_consensus), length(conns)), names(conns))
    }), "invalid public fixed-cohort Count")
})

test_that("public vector metadata follows the signed backend selection", {
  plan <- list(version = "test-plan")
  plan_hash <- paste(rep("a", 64L), collapse = "")
  manifest_hash <- paste(rep("b", 64L), collapse = "")
  base_release <- list(
    backend = .DSVERT_CLIENT_VECTOR_BACKEND,
    mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
    mechanism_plan = plan, plan_sha256 = plan_hash,
    backend_selection = NULL, backend_assessment = NULL,
    epsilon = 1, delta = 1e-8, implementation_delta = "1/1000000000",
    delta_aggregation = "max_per_peer_not_sum",
    capsule_id = paste(rep("c", 64L), collapse = ""),
    manifest_sha256 = manifest_hash,
    final_vector_root = paste(rep("d", 64L), collapse = ""),
    coordinate_order_sha256 = paste(rep("e", 64L), collapse = ""),
    coordinate_count = 2L, sticky_replay = TRUE)
  context <- list(
    release = base_release,
    lattice = list(
      natural_l1_sensitivity = 2,
      natural_l2_sensitivity = sqrt(2),
      output_lattice_bits = 8L,
      output_lattice_scale = 256),
    manifest = list(workload = list(
      capsule_mechanism = list(mechanism = "discrete-laplace"),
      mechanism_selection = list())),
    status = list(
      peer_a = list(noise_root = list(
        privacy_epoch = 1, key_id = "root-a")),
      peer_b = list(noise_root = list(
        privacy_epoch = 1, key_id = "root-b"))),
    adjacency = "add_remove_patient")

  convolution <- .dsvert_dp_vector_public_metadata(context)
  expect_identical(convolution$implementation,
                   .DSVERT_CLIENT_VECTOR_BACKEND)
  expect_identical(convolution$mechanism,
                   .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM)
  expect_identical(convolution$sampler, .DSVERT_CLIENT_VECTOR_SAMPLER)
  expect_identical(
    convolution$postprocessing,
    "signed-Ring128-decode-then-fixed-public-coordinate-clamp-v1")
  expect_false(convolution$one_joint_draw)
  expect_identical(convolution$mechanism_plan, plan)
  expect_identical(convolution$plan_sha256, plan_hash)

  context$release$backend <- .DSVERT_CLIENT_VECTOR_EXACT_BACKEND
  context$release$mechanism <- .DSVERT_CLIENT_VECTOR_EXACT_RELEASE_MECHANISM
  context$release$delta_aggregation <- "single_joint_draw"
  exact <- .dsvert_dp_vector_public_metadata(context)
  expect_identical(exact$implementation,
                   .DSVERT_CLIENT_VECTOR_EXACT_BACKEND)
  expect_identical(exact$mechanism,
                   .DSVERT_CLIENT_VECTOR_EXACT_RELEASE_MECHANISM)
  expect_identical(exact$sampler, .DSVERT_CLIENT_VECTOR_EXACT_SAMPLER)
  expect_identical(
    exact$postprocessing,
    paste0("fixed-public-coordinate-clamp-inside-exact-GC-before-",
           "selective-sharing-v1"))
  expect_true(exact$one_joint_draw)
})
