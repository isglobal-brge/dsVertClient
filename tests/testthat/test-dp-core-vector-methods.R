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

.dp_meanvar_synopsis_fixture <- function(k = 2L, gaussian = FALSE) {
  fixture <- .dp_core_vector_fixture(gaussian = gaussian)
  peers <- paste0("site_", letters[seq_len(k)])
  fixture$conns <- stats::setNames(lapply(seq_along(peers), function(index) {
    structure(index, class = "fake")
  }), peers)
  fixture$status <- stats::setNames(lapply(peers, function(peer) {
    list(policy = list(
      adjacency = "add_remove_patient", peer_count = as.integer(k),
      unit_capacity = 100))
  }), peers)

  release <- fixture$run$release
  release$version <- "dsvert-stateless-synopsis-public-vector-client-v1"
  release$backend <- if (isTRUE(gaussian)) {
    .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND
  } else {
    .DSVERT_CLIENT_VECTOR_BACKEND
  }
  release$backend_selection <- NULL
  release$backend_assessment <- NULL
  if (isTRUE(gaussian)) {
    release$mechanism_plan <-
      release$manifest$workload$mechanism_selection$gaussian_plan
    release$plan_sha256 <-
      release$manifest$workload$mechanism_selection$gaussian_plan_sha256
  }
  release$manifest$admission <- list(
    adjacency = "add_remove_patient", unit_capacity = 100)
  release[c(
    "capsule_id", "history_gate", "request_limit", "operation_limit"
  )] <- NULL
  bindings <- c(
    artifact_key = "a", execution_id = "b", contract_sha256 = "c",
    attempt_sha256 = "d", source_contract_sha256 = "e",
    result_set_sha256 = "f", final_vector_root = "3")
  for (field in names(bindings)) {
    release[[field]] <- strrep(bindings[[field]], 64L)
  }
  release$signed_provenance <- c(list(
    version = "dsvert-stateless-synopsis-public-provenance-v1",
    ordered_peer_pinset = as.list(stats::setNames(
      paste0("pin-", peers), peers)),
    designated_noise_peers = as.list(peers[1:2])),
    release[names(bindings)], list(
      compile_receipts = stats::setNames(vector("list", k), peers),
      release_receipts = stats::setNames(vector("list", 2L), peers[1:2]),
      replay_responses = stats::setNames(vector("list", 2L), peers[1:2]),
      protected_shares_included = FALSE,
      source_values_included = FALSE,
      intermediate_payload_exposed = FALSE,
      durable_replay = TRUE))
  class(release) <- c(
    "dsvert_synopsis_public_vector", "dsvert_joint_dp_vector", "list")
  fixture$run <- list(
    release = release, layout = fixture$run$layout, status = fixture$status,
    manifest_bundle = list(manifest_sha256 = release$manifest_sha256))
  fixture
}

.dp_cross_contingency_synopsis_fixture <- function(
    k = 2L, duplicate_physical = FALSE) {
  fixture <- .dp_meanvar_synopsis_fixture(k = k)
  peers <- names(fixture$conns)
  authorities <- peers[1:2]
  scale <- 2^8
  artifact <- list(
    version = .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION,
    spec_version = "v2", analysis_id = "cross_table",
    implementation_state = "cross_owner_exact_gc_materialized",
    alignment_group = "cohort-v1",
    left = list(
      dataset = "leftdata", column = "disease",
      owner_peer = authorities[[1L]], levels = as.list(c("no", "yes"))),
    right = list(
      dataset = "rightdata", column = "exposure",
      owner_peer = authorities[[2L]], levels = as.list(c("high", "low"))),
    participating_peers = as.list(authorities),
    computation_peers = as.list(authorities), coordinate_count = 4L,
    coordinate_order = paste0(
      "canonical_left_level_rows_then_canonical_right_level_columns_",
      "column_major_v1"),
    source_coordinate_scaling =
      "all_coordinates_already_on_common_numeric_lattice_v1",
    private_input_layout = paste0(
      "capacity_padded_one_hot_by_public_level_then_side_",
      "manifest_order_v1"),
    repeated_record_policy =
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_UNIT_POLICY,
    missingness_policy = paste0(
      "missing_out_of_domain_or_conflicting_side_is_all_zero_and_",
      "contributes_no_joint_cell_v1"),
    statistic_maximum = as.list(rep(100 * scale, 4L)),
    selected_l1_sensitivity = 1, selected_l2_sensitivity = 1,
    numeric_certificate = list(frac_bits = 8L),
    transcript = list(padded_units = 100L))
  if (isTRUE(duplicate_physical)) {
    artifact$left$column <- "marker"
    artifact$right$column <- "marker"
  }
  manifest <- fixture$run$release$manifest
  manifest$workload$coordinate_count <- 9L
  manifest$workload$release_lattice$natural_l1_sensitivity <- 4
  manifest$workload$release_lattice$integer_l1_sensitivity_steps <-
    4 * scale
  manifest$workload$release_lattice$natural_l2_sensitivity <- 2
  manifest$workload$release_lattice$integer_l2_sensitivity_steps <-
    2 * scale
  manifest$workload$families$admitted_count <- list(
    version = "admitted-unit-count-v2", owner_peer = authorities[[1L]],
    dataset = "leftdata", statistic_minimum = 0,
    statistic_maximum = 100)
  for (family in c(
      "numeric_moments", "numeric_pair_moments", "gaussian_models",
      "fixed_numeric_histograms")) {
    manifest$workload$families[[family]] <- list(artifacts = list())
  }
  manifest$workload$families$categorical_marginals <- list(artifacts = list(
    `leftdata::disease` = list(
      owner_peer = authorities[[1L]], dataset = "leftdata",
      column = "disease", levels = as.list(c("no", "yes"))),
    `rightdata::exposure` = list(
      owner_peer = authorities[[2L]], dataset = "rightdata",
      column = "exposure", levels = as.list(c("high", "low")))))
  manifest$workload$families$categorical_pairs <- list(
    sets = list(), cross_artifacts = list(cross_table = artifact))
  manifest$workload$families$correlation_artifacts <- list()
  manifest$workload$families$describe_artifacts <- list()
  manifest$workload$families$survival_artifacts <- list()
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  release <- fixture$run$release
  release$manifest <- manifest
  release$coordinate_count <- layout$coordinate_count
  release$coordinate_order_sha256 <- layout$sha256
  # admitted; left/right marginals; canonical disease-by-exposure cells.
  release$values <- c(42, 12, 30, 20, 22, 8.5, 3.75, 2.25, 9.5)
  fixture$run$release <- release
  fixture$run$layout <- layout
  fixture
}

.dp_contingency_public_federation <- function(k = 2L) {
  sites <- paste0("site_", letters[seq_len(k)])
  id_columns <- stats::setNames(paste0("pid_", letters[seq_len(k)]), sites)
  rows <- lapply(seq_along(sites), function(index) {
    site <- sites[[index]]
    columns <- c(
      id_columns[[index]], "marker", paste0("unique_", letters[[index]]),
      paste0("numeric_", letters[[index]]))
    data.frame(
      server = rep(site, length(columns)), column = columns,
      kind = c("identifier", "categorical", "categorical", "numeric"),
      role = c("id", "data", "data", "data"),
      stringsAsFactors = FALSE)
  })
  structure(list(
    version = 2L, symbol = "DA", sites = sites,
    source_symbols = stats::setNames(paste0("source_", letters[seq_len(k)]),
                                     sites),
    id_columns = id_columns,
    public_schema = do.call(rbind, rows), attestation = list()),
    class = c("ds.vertFederation", "list"))
}

.dp_meanvar_synopsis_mock_runner <- function(fixture, counter) {
  force(fixture)
  force(counter)
  function(datasources, status = NULL, local_projection = NULL, .aggregate) {
    counter$synopsis <- counter$synopsis + 1L
    expect_identical(datasources, fixture$conns)
    expect_null(status)
    fixture$run
  }
}

.dp_count_public_add_execution <- function(
    value = "7", upper = "10", radius = 13L) {
  release <- list(
    version = "dsvert-dp-count-release-v1",
    artifact_key = strrep("1", 64L),
    contract_sha256 = strrep("2", 64L),
    analysis_binding_sha256 = strrep("3", 64L),
    worker_static_sha256 = strrep("4", 64L),
    circuit = paste0("joint-dp-laplace-v2/", strrep("5", 64L)),
    mechanism = list(
      family = "discrete_laplace",
      version = .DSVERT_DP_ANALYSIS_COUNT_TV_MECHANISM,
      sampler = .DSVERT_DP_ANALYSIS_COUNT_TV_SAMPLER,
      epsilon = 1,
      delta = 1e-6,
      implementation_delta = 1e-9,
      sensitivity_l1 = 1),
    bounds = list(lower = "0", upper = upper),
    value = value,
    source_identity_pk = strrep("A", 43L),
    finalizer_identity_pk = strrep("B", 43L),
    backend = "exact-gc-joint-dp-laplace-ring127-v2",
    postprocessing = "one-joint-noise-draw-and-one-clamp-inside-exact-gc",
    intermediate_values_exposed = FALSE,
    public_openings = 1,
    release_sha256 = strrep("6", 64L),
    signature = "signed-release-fixture")
  list(
    version = "dsvert-dp-count-execution-result-v1",
    mode = "add_remove_dp",
    payload = list(
      release = release,
      finalizer_peer = "site_b",
      accuracy_95_abs = radius,
      accuracy_95_confidence = 0.95,
      accuracy_95_method =
        "conservative_truncated_dyadic_two_geometric_tail_bound_v1"))
}

.dp_count_public_pin <- function(byte) {
  chartr("+/", "-_", sub(
    "=+$", "", gsub("[\r\n]", "", jsonlite::base64_enc(
      as.raw(rep(byte, 32L)))), perl = TRUE))
}

test_that("legacy Count validator is explicit-test-only", {
  fixture <- .dp_core_vector_fixture()
  expect_identical(names(formals(ds.vertDPCount)),
                   c("data_name", "server", "datasources"))
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_count_execute_v1 = function(
        data_name, datasources, .aggregate) {
      counter$calls <- counter$calls + 1L
      expect_identical(data_name, "cohort")
      expect_identical(datasources, fixture$conns)
      .dp_count_public_add_execution()
    },
    .package = "dsVertClient")

  result <- .dsvert_dp_count_impl(
    "cohort", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE),
    .execute = .dsvert_dp_count_execute_v1)
  expect_s3_class(result, "ds.vertDPCount")
  expect_true(all(c(
    "released", "server", "epsilon", "accuracy_95_abs",
    "uncertainty_scope") %in% names(result)))
  expect_output(print(result), "dsVert DP count: 7 [ site_b ]", fixed = TRUE)
  expect_identical(counter$calls, 1L)
  expect_identical(result$value, 7)
  expect_identical(result$server, "site_b")
  expect_identical(result$mechanism,
                   .DSVERT_DP_ANALYSIS_COUNT_TV_MECHANISM)
  expect_identical(result$implementation,
                   "exact-gc-joint-dp-laplace-ring127-v2")
  expect_identical(result$accuracy_95_abs, 10)
  expect_identical(result$epsilon, 1)
  expect_identical(result$randomness,
                   "two_persistent_identity_seeds_joint_exact_gc_v1")
  expect_identical(result$composition_rule,
                   "one_sticky_release_per_canonical_signed_artifact")
  expect_identical(result$artifact_l1_sensitivity, 1)
  expect_identical(result$privacy, list(
    per_artifact_epsilon = 1,
    per_artifact_delta = 1e-6,
    sticky_noise = TRUE,
    finite_global_composition_claim = FALSE,
    distinct_artifacts_compose = TRUE,
    public_openings = 1L))
  expect_false(any(c(
    "history_gate", "request_limit", "operation_limit", "capsule_id",
    "lifetime_budget") %in% names(result)))
  expect_match(result$uncertainty_scope, "mechanism noise")

  split_mode <- .dp_count_public_add_execution()
  split_mode$mode <- c("add_remove_dp", "fixed_cohort_public")
  testthat::local_mocked_bindings(
    .dsvert_dp_count_execute_v1 = function(...) split_mode,
    .package = "dsVertClient")
  expect_error(.dsvert_dp_count_impl(
    "cohort", NULL, fixture$conns, function(...) NULL,
    .execute = .dsvert_dp_count_execute_v1),
    "Invalid closed Count execution result")
})

test_that("DP Count reads a signed durable Synopsis at K=2, K=3, and K=5", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .dp_meanvar_synopsis_fixture(k = k)
    peers <- names(fixture$conns)
    pins <- stats::setNames(vapply(seq_along(peers), function(index) {
      chartr("+/", "-_", sub("=+$", "", jsonlite::base64_enc(
        as.raw(rep(index, 32L)))))
    }, character(1L)), peers)
    fixture$run$release$signed_provenance$ordered_peer_pinset <- as.list(pins)
    fixture$run$release$signed_provenance$peer_pinset_sha256 <-
      .dsvert_vector_hash(as.list(pins))
    calls <- 0L
    run <- function(datasources, status = NULL, local_projection = NULL,
                    .aggregate) {
      calls <<- calls + 1L
      expect_identical(datasources, fixture$conns)
      expect_null(status)
      expect_null(local_projection)
      fixture$run
    }
    result <- .dsvert_dp_count_impl(
      "cohort", "site_b", fixture$conns,
      function(...) stop("unexpected raw DSI call", call. = FALSE),
      .run = run)
    expect_s3_class(result, "ds.vertDPCount")
    expect_identical(calls, 1L)
    expect_identical(result$value, 42.25)
    expect_identical(result$source_owner, "site_b")
    expect_identical(result$server, "site_b")
    expect_identical(result$coordinate_family, "admitted_count")
    expect_identical(result$coordinate_maximum, 100)
    expect_true(is.finite(result$accuracy_95_abs))
    expect_gte(result$accuracy_95_abs, 0)
    expect_lte(result$accuracy_95_abs, 100)
    expect_identical(result$release_provenance$designated_noise_peers,
                     as.list(peers[1:2]))
    expect_length(result$release_provenance$ordered_peer_pinset, k)
    expect_identical(result$privacy$finite_global_composition_claim, FALSE)
    expect_identical(result$privacy$unlimited_replay, TRUE)
    expect_false(any(c(
      "history_gate", "request_limit", "operation_limit", "capsule_id",
      "lifetime_budget", "quota") %in% names(result)))

    tampered <- fixture
    tampered$run$release$values[[1L]] <- 101
    expect_error(.dsvert_dp_count_impl(
      "cohort", "site_b", tampered$conns, function(...) NULL,
      .run = function(...) tampered$run), "violates its signed domain")
    expect_error(.dsvert_dp_count_impl(
      "other", "site_b", fixture$conns, function(...) NULL,
      .run = run), "does not match data_name")
    expect_error(.dsvert_dp_count_impl(
      "cohort", "site_a", fixture$conns, function(...) NULL,
      .run = run), "signed Count finalizer")
    wrong_pinset <- fixture
    wrong_pinset$run$release$signed_provenance$peer_pinset_sha256 <-
      paste(rep("0", 64L), collapse = "")
    expect_error(.dsvert_dp_count_impl(
      "cohort", "site_b", wrong_pinset$conns, function(...) NULL,
      .run = function(...) wrong_pinset$run), "noise authorities")
  }
})

test_that("DP contingency respects signed column-major orientation", {
  fixture <- .dp_meanvar_synopsis_fixture()
  counter <- new.env(parent = emptyenv())
  counter$synopsis <- counter$legacy <- counter$sampling <-
    counter$claim <- 0L
  counter$selectors <- list()
  testthat::local_mocked_bindings(
    .dsvert_dp_synopsis_vector_run = function(
        datasources, status = NULL, local_projection = NULL, .aggregate) {
      counter$synopsis <- counter$synopsis + 1L
      counter$selectors[[counter$synopsis]] <- local_projection
      if (counter$synopsis == 1L) {
        counter$claim <- counter$claim + 1L
        counter$sampling <- counter$sampling + 1L
      }
      expect_identical(local_projection$family, "categorical_pair")
      expect_identical(local_projection$dataset, "cohort")
      expect_identical(local_projection$columns, c("disease", "exposure"))
      fixture$run
    },
    .dsvert_dp_capsule_vector_run = function(...) {
      counter$legacy <- counter$legacy + 1L
      stop("legacy capsule runner reached", call. = FALSE)
    },
    .package = "dsVertClient")

  result <- .dsvert_dp_contingency_impl(
    "cohort", "exposure", "disease", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_s3_class(result, "ds.vertDPContingency")
  expect_identical(counter$synopsis, 1L)
  expect_identical(counter$legacy, 0L)
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
  expect_identical(.dsvert_dp_table_contract(result), result)

  reversed <- .dsvert_dp_contingency_impl(
    "cohort", "disease", "exposure", NULL, fixture$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_identical(counter$selectors[[2L]], counter$selectors[[1L]])
  expect_identical(counter$claim, 1L)
  expect_identical(counter$sampling, 1L)
  expect_identical(counter$legacy, 0L)
  expect_identical(
    reversed$artifact_key, result$artifact_key)
  expect_identical(
    reversed$final_vector_root, result$final_vector_root)
  expect_equal(reversed$table, t(result$table), tolerance = 0)
  expect_identical(reversed$row_levels, result$col_levels)
  expect_identical(reversed$col_levels, result$row_levels)
})

test_that("DP contingency consumes one projected Synopsis for K=2,3,5", {
  forbidden <- c(
    "capsule_id", "privacy_epoch", "privacy_epochs", "noise_key_id",
    "noise_key_ids", "history_gate", "request_limit", "operation_limit",
    "lifetime_budget", "lifetime_composition", "privacy_accountant",
    "release_instance", "release_instance_id", "allocation_certificate",
    "capsule_coordinate_count")
  for (k in c(2L, 3L, 5L)) {
    fixture <- .dp_meanvar_synopsis_fixture(k = k)
    counter <- new.env(parent = emptyenv())
    counter$synopsis <- counter$legacy <- 0L
    runner <- function(
        datasources, status = NULL, local_projection = NULL, .aggregate) {
      counter$synopsis <- counter$synopsis + 1L
      expect_null(status)
      expect_identical(local_projection, list(
        version = .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION,
        family = "categorical_pair", dataset = "cohort",
        columns = c("disease", "exposure"), owner_peer = NULL))
      fixture$run
    }
    evaluate <- function() testthat::with_mocked_bindings(
      .dsvert_dp_contingency_impl(
        "cohort", "exposure", "disease", NULL, fixture$conns,
        function(...) stop("unexpected raw DSI call", call. = FALSE)),
      .dsvert_dp_synopsis_vector_run = runner,
      .dsvert_dp_capsule_vector_run = function(...) {
        counter$legacy <- counter$legacy + 1L
        stop("legacy capsule runner reached", call. = FALSE)
      }, .package = "dsVertClient")
    first <- evaluate()
    expect_identical(counter$synopsis, 1L, info = paste("K =", k))
    expect_identical(counter$legacy, 0L, info = paste("K =", k))
    expect_equal(first$table, matrix(
      c(8.5, 3.75, 2.25, 9.5), 2L, 2L,
      dimnames = list(c("no", "yes"), c("no", "yes"))))
    expect_false(first$cross_owner)
    expect_true(first$privacy$unlimited_replay)
    expect_true(first$privacy$replay_is_postprocessing)
    paths <- names(unlist(first, recursive = TRUE, use.names = TRUE))
    leaves <- sub("^.*[.]", "", paths)
    expect_length(intersect(forbidden, c(names(first), leaves)), 0L)

    chisq <- .dsvert_dp_chisq_from_release(
      first, simulations = 31L, mc_confidence = 0.9)
    fisher <- .dsvert_dp_fisher_from_release(
      first, simulations = 31L, mc_confidence = 0.9)
    epi <- ds.vertDPEpi2x2(first, exposed = "yes", event = "yes")
    direct <- ds.vertDPDirectStandardization(
      first, standard_weights = c(no = 0.4, yes = 0.6), event = "yes")
    indirect <- ds.vertDPIndirectStandardization(
      first, expected_rates = c(no = 0.2, yes = 0.3), event = "yes")
    expect_identical(counter$synopsis, 1L, info = paste("K =", k))
    expect_true(all(vapply(
      list(chisq, fisher, epi, direct, indirect),
      function(value) identical(value$additional_server_calls, 0L),
      logical(1L))), info = paste("K =", k))
    expect_identical(chisq$source_release$artifact_key, first$artifact_key)
    expect_false("capsule_id" %in% names(chisq$source_release))
    expect_identical(fisher$source_release$artifact_key, first$artifact_key)
    expect_false("capsule_id" %in% names(fisher$source_release))
    expect_equal(epi$point_estimates$risk_exposed,
                 9.5 / (9.5 + 3.75), tolerance = 0)
    expect_equal(direct$estimate,
      0.4 * 2.25 / (8.5 + 2.25) +
        0.6 * 9.5 / (3.75 + 9.5), tolerance = 1e-15)
    expect_true(is.finite(indirect$estimate))

    second <- evaluate()
    expect_identical(serialize(second, NULL, version = 3L),
                     serialize(first, NULL, version = 3L))
    expect_identical(counter$synopsis, 2L, info = paste("K =", k))
    expect_identical(counter$legacy, 0L, info = paste("K =", k))
  }
})

test_that(paste(
  "cross-owner contingency uses one canonical Synopsis draw for K=2,3,5",
  "and only transposes in the client"), {
  forbidden <- c(
    "capsule_id", "privacy_epoch", "privacy_epochs", "noise_key_id",
    "noise_key_ids", "history_gate", "request_limit", "operation_limit",
    "lifetime_budget", "lifetime_composition", "privacy_accountant",
    "release_instance", "release_instance_id", "allocation_certificate",
    "ledger", "reservation", "rate_limit", "catalog_limit", "quota")
  for (k in c(2L, 3L, 5L)) {
    fixture <- .dp_cross_contingency_synopsis_fixture(k)
    counter <- new.env(parent = emptyenv())
    counter$synopsis <- counter$claim <- counter$noise <- counter$legacy <- 0L
    counter$selectors <- list()
    runner <- function(
        datasources, status = NULL, local_projection = NULL, .aggregate) {
      counter$synopsis <- counter$synopsis + 1L
      counter$selectors[[counter$synopsis]] <- local_projection
      if (counter$synopsis == 1L) {
        counter$claim <- counter$claim + 1L
        counter$noise <- counter$noise + 1L
      }
      expect_identical(datasources, fixture$conns)
      expect_null(status)
      expect_identical(local_projection$version,
        .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION)
      expect_identical(local_projection$family, "categorical_pair")
      expect_identical(local_projection$columns, c("disease", "exposure"))
      fixture$run
    }
    evaluate <- function(dataset, row, column) {
      testthat::with_mocked_bindings(
        .dsvert_dp_contingency_impl(
          dataset, row, column, NULL, fixture$conns,
          function(...) stop("unexpected raw DSI call", call. = FALSE)),
        .dsvert_dp_synopsis_vector_run = runner,
        .dsvert_dp_capsule_vector_run = function(...) {
          counter$legacy <- counter$legacy + 1L
          stop("legacy capsule runner reached", call. = FALSE)
        }, .package = "dsVertClient")
    }

    first <- evaluate("leftdata", "disease", "exposure")
    expect_true(first$cross_owner, info = paste("K =", k))
    expect_identical(first$servers, fixture$run$status |> names() |> head(2L))
    expect_identical(first$datasets, c("leftdata", "rightdata"))
    expect_identical(
      first$unit_aggregation_policy,
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_UNIT_POLICY)
    expect_equal(first$table, matrix(
      c(8.5, 3.75, 2.25, 9.5), 2L, 2L,
      dimnames = list(c("no", "yes"), c("high", "low"))))
    expect_identical(.dsvert_dp_table_contract(first), first)
    expect_identical(counter$claim, 1L)
    expect_identical(counter$noise, 1L)
    expect_identical(counter$legacy, 0L)

    reverse <- evaluate("rightdata", "exposure", "disease")
    expect_identical(counter$selectors[[1L]]$columns,
                     counter$selectors[[2L]]$columns)
    expect_identical(counter$selectors[[1L]]$dataset, "leftdata")
    expect_identical(counter$selectors[[2L]]$dataset, "rightdata")
    expect_identical(reverse$artifact_key, first$artifact_key)
    expect_identical(reverse$final_vector_root, first$final_vector_root)
    expect_equal(reverse$table, t(first$table), tolerance = 0)
    expect_identical(counter$claim, 1L, info = paste("sticky Claim K =", k))
    expect_identical(counter$noise, 1L, info = paste("sticky draw K =", k))
    expect_identical(counter$legacy, 0L, info = paste("no capsule K =", k))

    chisq <- ds.vertChisqCross(
      reverse, fisher = TRUE, verbose = FALSE,
      simulations = 31L, mc_confidence = 0.9)
    epi <- ds.vertDPEpi2x2(reverse, exposed = "high", event = "yes")
    expect_identical(chisq$source_dp_release$artifact_key,
                     first$artifact_key)
    expect_identical(chisq$additional_server_calls, 0L)
    expect_identical(chisq$fisher$additional_server_calls, 0L)
    expect_identical(epi$additional_server_calls, 0L)
    expect_identical(counter$synopsis, 2L)
    paths <- names(unlist(first, recursive = TRUE, use.names = TRUE))
    leaves <- sub("^.*[.]", "", paths)
    expect_length(intersect(forbidden, c(names(first), leaves)), 0L)
  }
})

test_that(paste(
  "qualified duplicate physical columns map one cross Synopsis for K=2,3,5",
  "and bare duplicates fail before the runner"), {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .dp_cross_contingency_synopsis_fixture(
      k, duplicate_physical = TRUE)
    calls <- new.env(parent = emptyenv())
    calls$runner <- calls$legacy <- 0L
    runner <- function(
        datasources, status = NULL, local_projection = NULL, .aggregate) {
      calls$runner <- calls$runner + 1L
      expect_identical(
        local_projection$columns,
        c("site_a$marker", "site_b$marker"), info = paste("K =", k))
      fixture$run
    }
    evaluate <- function(row, column, dataset = "leftdata") {
      testthat::with_mocked_bindings(
        .dsvert_dp_contingency_impl(
          dataset, row, column, NULL, fixture$conns,
          function(...) stop("protected DSI reached", call. = FALSE)),
        .dsvert_dp_synopsis_vector_run = runner,
        .dsvert_dp_capsule_vector_run = function(...) {
          calls$legacy <- calls$legacy + 1L
          stop("legacy capsule runner reached", call. = FALSE)
        }, .package = "dsVertClient")
    }

    first <- evaluate("site_a$marker", "site_b$marker")
    reverse <- evaluate(
      "site_b$marker", "site_a$marker", dataset = "rightdata")
    expect_true(first$cross_owner, info = paste("K =", k))
    expect_equal(reverse$table, t(first$table), tolerance = 0,
                 info = paste("transpose K =", k))
    expect_identical(reverse$artifact_key, first$artifact_key)
    expect_identical(reverse$final_vector_root, first$final_vector_root)
    expect_identical(calls$runner, 2L, info = paste("runner K =", k))
    expect_identical(calls$legacy, 0L, info = paste("legacy K =", k))

    before <- calls$runner
    expect_error(evaluate("marker", "marker"), "distinct variables",
                 info = paste("bare duplicate K =", k))
    expect_identical(calls$runner, before,
                     info = paste("pre-runner K =", k))
  }
})

test_that(paste(
  "federation categorical references resolve locally before DSI for K=2,3,5"),
  {
    for (k in c(2L, 3L, 5L)) {
      federation <- .dp_contingency_public_federation(k)
      conns <- stats::setNames(rep(list(list()), k), federation$sites)
      calls <- new.env(parent = emptyenv())
      calls$aggregate <- calls$runner <- 0L
      received <- list()
      resolve <- function(value, datasources) {
        calls$aggregate <- calls$aggregate + 1L
        expect_identical(value, federation)
        expect_identical(datasources, conns)
        list(value = "leftdata", datasources = conns)
      }
      backend <- function(
          data_name, row_var, col_var, server, datasources, .aggregate) {
        calls$runner <- calls$runner + 1L
        received[[calls$runner]] <<- list(
          data_name = data_name, row_var = row_var, col_var = col_var,
          datasources = datasources)
        list(row_var = row_var, col_var = col_var)
      }
      evaluate <- function(row, column) testthat::with_mocked_bindings(
        ds.vertDPContingency(
          federation, row, column, datasources = conns),
        .dsvert_federation_argument = resolve,
        .dsvert_dp_contingency_impl = backend,
        .package = "dsVertClient")

      qualified <- evaluate("site_a$marker", "site_b$marker")
      expect_identical(
        received[[1L]][c("row_var", "col_var")],
        list(row_var = "site_a$marker", col_var = "site_b$marker"),
        info = paste("qualified K =", k))
      expect_identical(qualified$row_var, "site_a$marker")
      expect_identical(qualified$col_var, "site_b$marker")

      bare <- evaluate("unique_a", "unique_b")
      expect_identical(
        received[[2L]][c("row_var", "col_var")],
        list(row_var = "site_a$unique_a", col_var = "site_b$unique_b"),
        info = paste("bare-to-qualified K =", k))
      expect_identical(bare$row_var, "unique_a")
      expect_identical(bare$col_var, "unique_b")
      expect_identical(calls$aggregate, 2L)
      expect_identical(calls$runner, 2L)

      invalid <- list(
        ambiguous = list(c("marker", "marker"), "multiple owners"),
        wrong_owner = list(
          c("site_b$unique_a", "site_a$marker"), "not owned"),
        missing = list(
          c("site_a$missing", "site_b$marker"), "missing from public schema"),
        wrong_kind = list(
          c("site_a$numeric_a", "site_b$marker"), "expected categorical"),
        unknown_owner = list(
          c("site_z$marker", "site_b$marker"), "Unknown contingency server"),
        nested = list(
          c("site_a$dataset$marker", "site_b$marker"),
          "exact server\\$column"))
      for (name in names(invalid)) {
        before <- c(aggregate = calls$aggregate, runner = calls$runner)
        expect_error(
          evaluate(invalid[[name]][[1L]][[1L]],
                   invalid[[name]][[1L]][[2L]]),
          invalid[[name]][[2L]], info = paste(name, "K =", k))
        expect_identical(
          c(aggregate = calls$aggregate, runner = calls$runner), before,
          info = paste(name, "pre-DSI K =", k))
      }
    }
  })

test_that("cross chi-square preserves a federation for its single front door", {
  federation <- .dp_contingency_public_federation(2L)
  conns <- stats::setNames(rep(list(list()), 2L), federation$sites)
  delegated <- 0L
  expect_error(testthat::with_mocked_bindings(
    ds.vertChisqCross(
      federation, "site_a$marker", "site_b$marker", verbose = FALSE,
      datasources = conns, simulations = 31L),
    .dsvert_federation_argument = function(...) {
      stop("unexpected early federation unwrap", call. = FALSE)
    },
    ds.vertDPContingency = function(
        data_name, row_var, col_var, server, datasources) {
      delegated <<- delegated + 1L
      expect_identical(data_name, federation)
      expect_identical(row_var, "site_a$marker")
      expect_identical(col_var, "site_b$marker")
      expect_null(server)
      expect_identical(datasources, conns)
      stop("CONTINGENCY_FRONTDOOR", call. = FALSE)
    }, .package = "dsVertClient"), "CONTINGENCY_FRONTDOOR")
  expect_identical(delegated, 1L)
})

test_that("cross-owner contingency rejects all seven detached roots pre-map", {
  fixture <- .dp_cross_contingency_synopsis_fixture(k = 3L)
  bindings <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  for (field in bindings) {
    tampered <- fixture
    tampered$run$release[[field]] <- strrep("0", 64L)
    postprocess_calls <- 0L
    expect_error(testthat::with_mocked_bindings(
      .dsvert_dp_contingency_impl(
        "leftdata", "disease", "exposure", NULL, tampered$conns,
        function(...) stop("unexpected raw DSI call", call. = FALSE)),
      .dsvert_dp_synopsis_vector_run = function(...) tampered$run,
      .dsvert_dp_capsule_vector_run = function(...) {
        stop("legacy capsule runner reached", call. = FALSE)
      }, .dsvert_dp_capsule_single_block = function(...) {
        postprocess_calls <<- postprocess_calls + 1L
        stop("table mapping reached", call. = FALSE)
      }, .package = "dsVertClient"), "provenance is detached", info = field)
    expect_identical(postprocess_calls, 0L, info = field)
  }
})

test_that("cross-owner contingency fails closed on missing or ambiguous blocks", {
  fixture <- .dp_cross_contingency_synopsis_fixture(k = 3L)
  cross_name <- names(fixture$run$layout$blocks)[vapply(
    fixture$run$layout$blocks, function(block) {
      identical(block$family, "categorical_pairs")
    }, logical(1L))]
  expect_length(cross_name, 1L)
  malformed <- list(missing = fixture, ambiguous = fixture)
  malformed$missing$run$layout$blocks[[cross_name]] <- NULL
  malformed$ambiguous$run$layout$blocks[[paste0(cross_name, "::duplicate")]] <-
    malformed$ambiguous$run$layout$blocks[[cross_name]]
  for (name in names(malformed)) {
    candidate <- malformed[[name]]
    expect_error(testthat::with_mocked_bindings(
      .dsvert_dp_contingency_impl(
        "leftdata", "disease", "exposure", NULL, candidate$conns,
        function(...) stop("unexpected raw DSI call", call. = FALSE)),
      .dsvert_dp_synopsis_vector_run = function(...) candidate$run,
      .dsvert_dp_capsule_vector_run = function(...) {
        stop("legacy capsule runner reached", call. = FALSE)
      }, .package = "dsVertClient"),
      if (identical(name, "missing")) "does not contain" else "exactly one",
      info = name)
  }
})

test_that("DP contingency rejects detached Synopsis roots before table mapping", {
  fixture <- .dp_meanvar_synopsis_fixture(k = 3L)
  bindings <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  for (field in bindings) {
    tampered <- fixture
    tampered$run$release[[field]] <- strrep("0", 64L)
    postprocess_calls <- 0L
    expect_error(testthat::with_mocked_bindings(
      .dsvert_dp_contingency_impl(
        "cohort", "exposure", "disease", NULL, tampered$conns,
        function(...) NULL),
      .dsvert_dp_synopsis_vector_run = function(...) tampered$run,
      .dsvert_dp_capsule_vector_run = function(...) {
        stop("legacy capsule runner reached", call. = FALSE)
      },
      .dsvert_dp_capsule_single_block = function(...) {
        postprocess_calls <<- postprocess_calls + 1L
        stop("table mapping reached", call. = FALSE)
      }, .package = "dsVertClient"), "provenance is detached", info = field)
    expect_identical(postprocess_calls, 0L, info = field)
  }
})

test_that("vector statistics support fixed-cohort replacement adjacency", {
  synopsis <- .dp_meanvar_synopsis_fixture()
  for (peer in names(synopsis$status)) {
    synopsis$status[[peer]]$policy$adjacency <-
      "replace_one_fixed_cohort"
  }
  synopsis$run$status <- synopsis$status
  synopsis$run$release$manifest$admission$adjacency <-
    "replace_one_fixed_cohort"
  testthat::local_mocked_bindings(
    .dsvert_dp_synopsis_vector_run = function(...) synopsis$run,
    .dsvert_dp_capsule_vector_run = function(...) {
      stop("legacy capsule runner reached", call. = FALSE)
    },
    .package = "dsVertClient")

  table <- .dsvert_dp_contingency_impl(
    "cohort", "exposure", "disease", NULL, synopsis$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  moments <- .dsvert_dp_meanvar_impl(
    "cohort", "age", NULL, synopsis$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_identical(table$adjacency, "replace_one_fixed_cohort")
  expect_identical(table$artifact_l1_sensitivity, 2)
  expect_equal(table$artifact_l2_sensitivity, sqrt(2), tolerance = 0)
  expect_identical(moments$adjacency, "replace_one_fixed_cohort")
  expect_identical(moments$artifact_l1_sensitivity, 3)
})

test_that("DP mean/variance consumes one no-lifetime Synopsis for K=2,3,5", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .dp_meanvar_synopsis_fixture(k = k)
    counter <- new.env(parent = emptyenv())
    counter$synopsis <- 0L
    counter$legacy <- 0L
    result <- testthat::with_mocked_bindings(
      .dsvert_dp_meanvar_impl(
        "cohort", "age", NULL, fixture$conns,
        function(...) stop("unexpected raw DSI call", call. = FALSE)),
      .dsvert_dp_synopsis_vector_run =
        .dp_meanvar_synopsis_mock_runner(fixture, counter),
      .dsvert_dp_capsule_vector_run = function(...) {
        counter$legacy <- counter$legacy + 1L
        stop("legacy capsule runner reached", call. = FALSE)
      }, .package = "dsVertClient")

    expect_s3_class(result, "ds.vertDPMeanVar")
    expect_identical(counter$synopsis, 1L, info = paste("K =", k))
    expect_identical(counter$legacy, 0L, info = paste("K =", k))
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
    expect_identical(result$artifact_key,
                     fixture$run$release$artifact_key)
    expect_identical(result$execution_id,
                     fixture$run$release$execution_id)
    expect_identical(result$contract_sha256,
                     fixture$run$release$contract_sha256)
    expect_identical(result$attempt_sha256,
                     fixture$run$release$attempt_sha256)
    expect_identical(result$source_contract_sha256,
                     fixture$run$release$source_contract_sha256)
    expect_identical(result$result_set_sha256,
                     fixture$run$release$result_set_sha256)
    expect_identical(result$final_vector_root,
                     fixture$run$release$final_vector_root)
    expect_identical(result$release_provenance,
                     fixture$run$release$signed_provenance)
    expect_true(result$privacy$unlimited_replay)
    expect_true(result$privacy$replay_is_postprocessing)
    expect_identical(names(result$privacy), c(
      "version", "per_artifact_epsilon", "per_artifact_delta",
      "sticky_noise", "unlimited_replay", "replay_is_postprocessing",
      "public_openings", "distinct_artifacts_compose",
      "finite_global_composition_claim"))
    expect_identical(
      result$composition_rule,
      "one_sticky_release_per_canonical_signed_artifact")
    expect_false(result$privacy$finite_global_composition_claim)
    expect_true(result$privacy$distinct_artifacts_compose)
    legacy_fields <- c(
      "capsule_id", "privacy_epoch", "privacy_epochs", "noise_key_id",
      "noise_key_ids", "history_gate", "request_limit", "operation_limit",
      "lifetime_budget", "lifetime_composition", "privacy_accountant",
      "release_instance", "release_instance_id", "allocation_certificate",
      "capsule_coordinate_count")
    paths <- names(unlist(result, recursive = TRUE, use.names = TRUE))
    leaf_names <- sub("^.*[.]", "", paths)
    expect_length(intersect(legacy_fields, c(names(result), leaf_names)), 0L)
  }
})

test_that("DP mean/variance replay is byte-identical", {
  fixture <- .dp_meanvar_synopsis_fixture(k = 5L)
  counter <- new.env(parent = emptyenv())
  counter$synopsis <- counter$legacy <- 0L
  evaluate <- function() testthat::with_mocked_bindings(
    .dsvert_dp_meanvar_impl(
      "cohort", "age", NULL, fixture$conns,
      function(...) stop("unexpected raw DSI call", call. = FALSE)),
    .dsvert_dp_synopsis_vector_run =
      .dp_meanvar_synopsis_mock_runner(fixture, counter),
    .dsvert_dp_capsule_vector_run = function(...) {
      counter$legacy <- counter$legacy + 1L
      stop("legacy capsule runner reached", call. = FALSE)
    }, .package = "dsVertClient")
  first <- evaluate()
  second <- evaluate()
  expect_identical(serialize(second, NULL, version = 3L),
                   serialize(first, NULL, version = 3L))
  expect_identical(counter$synopsis, 2L)
  expect_identical(counter$legacy, 0L)
})

test_that("DP mean/variance rejects detached Synopsis hashes before postprocess", {
  fixture <- .dp_meanvar_synopsis_fixture(k = 3L)
  bindings <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  for (field in bindings) {
    tampered <- fixture
    tampered$run$release[[field]] <- strrep("0", 64L)
    postprocess_calls <- 0L
    expect_error(testthat::with_mocked_bindings(
      .dsvert_dp_meanvar_impl(
        "cohort", "age", NULL, tampered$conns, function(...) NULL),
      .dsvert_dp_synopsis_vector_run = function(...) tampered$run,
      .dsvert_dp_capsule_vector_run = function(...) {
        stop("legacy capsule runner reached", call. = FALSE)
      },
      .dsvert_dp_capsule_single_block = function(...) {
        postprocess_calls <<- postprocess_calls + 1L
        stop("postprocess reached", call. = FALSE)
      }, .package = "dsVertClient"), "provenance is detached",
      info = field)
    expect_identical(postprocess_calls, 0L, info = field)
  }

  tampered <- fixture
  tampered$run$release$signed_provenance$durable_replay <- FALSE
  postprocess_calls <- 0L
  expect_error(testthat::with_mocked_bindings(
    .dsvert_dp_meanvar_impl(
      "cohort", "age", NULL, tampered$conns, function(...) NULL),
    .dsvert_dp_synopsis_vector_run = function(...) tampered$run,
    .dsvert_dp_capsule_single_block = function(...) {
      postprocess_calls <<- postprocess_calls + 1L
      stop("postprocess reached", call. = FALSE)
    }, .package = "dsVertClient"), "durable replay binding")
  expect_identical(postprocess_calls, 0L)
})

test_that("remaining core vector methods preserve one signed Gaussian L2 release", {
  synopsis <- .dp_meanvar_synopsis_fixture(k = 2L, gaussian = TRUE)
  counter <- new.env(parent = emptyenv())
  counter$capsule <- counter$synopsis <- 0L
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) {
      counter$capsule <- counter$capsule + 1L
      stop("legacy capsule runner reached", call. = FALSE)
    },
    .dsvert_dp_synopsis_vector_run = function(...) {
      counter$synopsis <- counter$synopsis + 1L
      synopsis$run
    },
    .package = "dsVertClient")

  table <- .dsvert_dp_contingency_impl(
    "cohort", "exposure", "disease", NULL, synopsis$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  moments <- .dsvert_dp_meanvar_impl(
    "cohort", "age", NULL, synopsis$conns,
    function(...) stop("unexpected raw DSI call", call. = FALSE))
  expect_identical(counter$capsule, 0L)
  expect_identical(counter$synopsis, 2L)
  for (result in list(table, moments)) {
    expect_identical(result$mechanism,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
    expect_identical(result$implementation,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND)
    expect_identical(result$sampler,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER)
    expect_identical(result$sensitivity_norm, "l2")
    expect_identical(result$mechanism_selection,
                     synopsis$run$release$manifest$workload$mechanism_selection)
  }
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
  synopsis <- .dp_meanvar_synopsis_fixture()
  runner <- function(datasources, status = NULL, .aggregate) fixture$run
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = runner,
    .dsvert_dp_synopsis_vector_run = function(...) synopsis$run,
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

  ambiguous <- synopsis
  duplicate <- ambiguous$run$layout$blocks[[
    "numeric_moments::age"]]
  ambiguous$run$layout$blocks[["numeric_moments::age_alias"]] <- duplicate
  testthat::local_mocked_bindings(
    .dsvert_dp_synopsis_vector_run = function(...) ambiguous$run,
    .package = "dsVertClient")
  expect_error(.dsvert_dp_meanvar_impl(
    "cohort", "age", NULL, fixture$conns, function(...) NULL),
    "exactly one signed numeric-moment block")
})

test_that("fixed-cohort Count exposes only the signed K-consensus value", {
  fixture <- .dp_core_vector_fixture()
  pins <- list(
    site_a = .dp_count_public_pin(1L),
    site_b = .dp_count_public_pin(2L))
  declaration <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = "dsvert-fixed-cohort-count-declaration-v1",
    domain = "study-domain",
    cohort_id = "cohort-v1",
    dataset_id = "cohort-table",
    dataset_version = "v1",
    privacy_unit_column = "patient_id",
    alignment_purpose = "patient-record-alignment-v1",
    adjacency = "replace_one_fixed_cohort",
    fixed_cohort_size = 100L,
    peer_pins = pins))
  execution <- list(
    version = "dsvert-dp-count-execution-result-v1",
    mode = "fixed_cohort_public",
    payload = list(
      declaration = declaration,
      receipt_set_sha256 = strrep("f", 64L),
      peer_count = 2L))
  testthat::local_mocked_bindings(
    .dsvert_dp_count_execute_v1 = function(...) execution,
    .package = "dsVertClient")

  result <- .dsvert_dp_count_impl(
    "cohort", "site_b", fixture$conns,
    function(...) stop("fixed Count used an obsolete endpoint"),
    .execute = .dsvert_dp_count_execute_v1)
  expect_identical(result$value, 100)
  expect_identical(result$server, "site_b")
  expect_identical(result$mechanism, "public_fixed_cohort_size_v1")
  expect_identical(result$implementation,
                   "custodian_owned_signed_K_consensus")
  expect_identical(result$sensitivity, 0)
  expect_identical(result$epsilon, 0)
  expect_identical(result$delta, 0)
  expect_identical(result$accuracy_95_abs, 0)
  expect_identical(result$data_dependency,
    "public_fixed_contract_validated_against_current_aligned_snapshot")
  expect_identical(result$declaration, declaration)
  expect_identical(result$privacy$finite_global_composition_claim, FALSE)
  expect_identical(result$privacy$sticky_noise, FALSE)
  expect_false(grepl(
    "capsule|status|lifetime",
    paste(deparse(body(.dsvert_dp_count_legacy_impl)), collapse = " "),
    ignore.case = TRUE))
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
