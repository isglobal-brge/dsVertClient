.dp_survival_client_release <- function(delayed_entry = FALSE) {
  .dsvert_dp_survival_vector_result(
    .dp_survival_vector_capsule(delayed_entry, decoded = TRUE),
    "protected", "primary")
}
.dp_survival_vector_capsule <- function(delayed_entry = FALSE,
                                        decoded = TRUE, gaussian = FALSE) {
  grid <- as.numeric(1:4)
  exit <- cbind(
    `0` = c(0, 1, 0, 0),
    A = c(1, 0, 0, 1),
    B = c(0, 0, 1, 0))
  entry <- if (delayed_entry) c(2, 1, 1, 0) else NULL
  histogram <- c(entry, as.vector(exit), 0)
  block_l1 <- if (delayed_entry) 2L else 1L
  block_l2 <- sqrt(block_l1)
  natural_l1 <- 1L + block_l1
  artifact <- list(
    version = "v1", dataset = "protected", owner_peer = "site_a",
    time = "time", event = "event",
    entry = if (delayed_entry) "entry" else "none",
    censor = "0", causes = c("A", "B"), time_grid = grid,
    time_bounds = c(0, 4),
    coordinate_order = paste0(
      "entry_bins_if_any_then_exit_time_within_outcome_then_invalid_bin"),
    repeated_record_policy = paste0(
      "earliest_event_else_latest_censor_then_cause_then_entry_",
      "deterministic_v2"),
    missingness_policy = paste0(
      "NA_NaN_Inf_or_out_of_domain_selected_unit_enters_invalid_bin"),
    coordinate_count = length(histogram), l1_sensitivity = block_l1,
    l2_sensitivity = block_l2, statistic_maximum = 100)
  manifest <- list(workload = list(
    coordinate_count = 1L + length(histogram),
    release_lattice = list(
      output_lattice_bits = 8L, output_lattice_scale = 256,
      natural_l1_sensitivity = natural_l1,
      integer_l1_sensitivity_steps = natural_l1 * 256,
      natural_l2_sensitivity = sqrt(1 + block_l2^2),
      integer_l2_sensitivity_steps = sqrt(1 + block_l2^2) * 256),
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
      numeric_moments = list(artifacts = list()),
      numeric_pair_moments = list(artifacts = list()),
      gaussian_models = list(artifacts = list()),
      fixed_numeric_histograms = list(artifacts = list()),
      categorical_marginals = list(artifacts = list()),
      categorical_pairs = list(sets = list()), correlation_artifacts = list(),
      describe_artifacts = list(),
      survival_artifacts = list(primary = artifact))))
  if (isTRUE(gaussian)) {
    total_coordinates <- 1L + length(histogram)
    global_l2 <- sqrt(1 + block_l2^2)
    sensitivity_steps <- format(
      global_l2 * 256, digits = 17L, scientific = TRUE, trim = TRUE)
    plan <- list(
      version = .DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION,
      mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
      sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER,
      request_binding_sha256 = digest::digest(
        paste0("survival-gaussian-request-", delayed_entry),
        "sha256", serialize = FALSE),
      total_coordinate_count = as.numeric(total_coordinates),
      maximum_chunk_coordinates = as.numeric(total_coordinates),
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
      simultaneous_95_abs = "24")
    plan <- .client_complete_gaussian_plan_v2(plan)
    manifest$workload$mechanism_selection <- list(
      gaussian_calibration_request = list(
        epsilon = "1e+00", delta = format(
          2^-100, digits = 17L, scientific = TRUE, trim = TRUE),
        l2_sensitivity_steps = sensitivity_steps,
        total_coordinate_count = as.numeric(total_coordinates)),
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
    coordinate_count = 1L + length(histogram),
    values = c(4, histogram), epsilon = 1, delta = 2^-100,
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

.dp_survival_synopsis_run <- function(k = 2L, delayed_entry = FALSE,
                                      gaussian = FALSE) {
  run <- .dp_survival_vector_capsule(
    delayed_entry = delayed_entry, decoded = FALSE, gaussian = gaussian)
  peers <- paste0("site_", letters[seq_len(k)])
  run$status <- stats::setNames(lapply(peers, function(peer) {
    list(policy = list(
      adjacency = "add_remove_patient", peer_count = as.integer(k),
      unit_capacity = 100))
  }), peers)

  release <- run$release
  release$version <- "dsvert-stateless-synopsis-public-vector-client-v1"
  release$backend <- if (isTRUE(gaussian)) {
    .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND
  } else {
    .DSVERT_CLIENT_VECTOR_BACKEND
  }
  release$backend_selection <- NULL
  release$backend_assessment <- NULL
  release$manifest_sha256 <- strrep("b", 64L)
  release$sticky_replay <- TRUE
  release$delta_aggregation <- "max_per_peer_not_sum"
  release$manifest$admission <- list(
    adjacency = "add_remove_patient", unit_capacity = 100)
  if (isTRUE(gaussian)) {
    release$mechanism_plan <-
      release$manifest$workload$mechanism_selection$gaussian_plan
    release$plan_sha256 <-
      release$manifest$workload$mechanism_selection$gaussian_plan_sha256
  }
  release[c(
    "capsule_id", "history_gate", "request_limit", "operation_limit"
  )] <- NULL
  bindings <- c(
    artifact_key = "1", execution_id = "2", contract_sha256 = "3",
    attempt_sha256 = "4", source_contract_sha256 = "5",
    result_set_sha256 = "6", final_vector_root = "7")
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
  run$release <- release
  run$manifest_bundle <- list(manifest_sha256 = release$manifest_sha256)
  run$conns <- stats::setNames(lapply(seq_along(peers), function(index) {
    structure(list(index = index), class = "fake")
  }), peers)
  run
}

test_that("survival maps only its signed final-vector block", {
  for (delayed in c(FALSE, TRUE)) {
    capsule <- .dp_survival_vector_capsule(delayed, decoded = TRUE)
    result <- .dsvert_dp_survival_vector_result(
      capsule, "protected", "primary")
    expect_identical(result$server, "site_a")
    expect_identical(result$histogram,
                     capsule$release$values[-1L])
    expect_true(is.na(result$clipped_coordinates))
    expect_false(result$clipping_observable)
    expect_identical(result$security_claim,
                     .dsvert_dp_capsule_security_claim())
    expect_false(result$security_claim$malicious_peer_security)
    expect_false(
      result$security_claim$unconditional_non_reconstruction_guarantee)
    expect_true(result$history_gate)
    expect_false(result$request_limit)
    expect_true(result$operation_limit)
    expect_identical(result$global_l1_sensitivity,
                     if (delayed) 3 else 2)
    postprocessed <- .dsvert_dp_survival_postprocess(result)
    expect_equal(postprocessed$curve$kaplan_meier,
                 if (delayed) c(0.5, 0.5, 0.25, 0) else
                   c(0.75, 0.75, 0.375, 0))
    expect_identical(
      postprocessed$status,
      "fixed_public_clamp_applied_preclamp_state_not_released")
  }
  capsule <- .dp_survival_vector_capsule(FALSE, decoded = TRUE)
  stale_gate <- capsule
  stale_gate$release$history_gate <- FALSE
  expect_error(.dsvert_dp_survival_vector_result(
    stale_gate, "protected", "primary"), "vector context is invalid")
  expect_error(.dsvert_dp_survival_vector_result(
    capsule, "protected", "primary", server = "site_b"), "does not own")
  capsule$release$manifest$workload$families$survival_artifacts$
    primary$coordinate_count <- 12L
  expect_error(.dsvert_dp_survival_vector_result(
    capsule, "protected", "primary"), "coordinate contract")
})

test_that("survival vector contract rejects semantic and clamp tampering", {
  base <- .dp_survival_vector_capsule(FALSE, decoded = TRUE)
  forged <- list(
    coordinate_order = "forged",
    repeated_record_policy = "forged",
    missingness_policy = "forged",
    l1_sensitivity = 0.5,
    l2_sensitivity = 0.5,
    statistic_maximum = 99,
    coordinate_count = 12L)
  for (field in names(forged)) {
    candidate <- base
    candidate$release$manifest$workload$families$survival_artifacts$
      primary[[field]] <- forged[[field]]
    expect_error(
      .dsvert_dp_survival_vector_result(
        candidate, "protected", "primary"),
      "coordinate contract is inconsistent", info = field)
  }

  candidate <- base
  start <- candidate$layout$blocks[["survival_artifacts::primary"]]$start
  candidate$release$values[[start]] <- 101
  expect_error(
    .dsvert_dp_survival_vector_result(
      candidate, "protected", "primary"),
    "violates its signed public clamp")

  for (delayed in c(FALSE, TRUE)) {
    candidate <- .dp_survival_vector_capsule(delayed, decoded = TRUE)
    baseline <- .dsvert_dp_survival_postprocess(
      .dsvert_dp_survival_vector_result(
        candidate, "protected", "primary"))
    block <- candidate$layout$blocks[["survival_artifacts::primary"]]
    candidate$release$values[[block$end]] <- 7
    result <- .dsvert_dp_survival_postprocess(
      .dsvert_dp_survival_vector_result(
        candidate, "protected", "primary"))
    expect_identical(result$not_in_analysis, 7)
    expect_identical(result$curve, baseline$curve)
    expect_identical(result$cumulative_incidence,
                     baseline$cumulative_incidence)
  }
})

test_that("survival propagates the signed Gaussian L2 plan", {
  for (delayed in c(FALSE, TRUE)) {
    capsule <- .dp_survival_vector_capsule(
      delayed, decoded = FALSE, gaussian = TRUE)
    result <- .dsvert_dp_survival_vector_result(
      capsule, "protected", "primary")
    expect_identical(result$mechanism,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
    expect_identical(result$sampler,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER)
    expect_equal(result$accuracy_95_abs_per_coordinate,
                 24 / 256, tolerance = 0)
    expect_equal(result$accuracy_simultaneous_95_abs,
                 24 / 256, tolerance = 0)
    expect_match(result$accuracy_simultaneous_method,
                 "discrete-Gaussian plan", fixed = TRUE)
    postprocessed <- .dsvert_dp_survival_postprocess(result)
    expect_identical(postprocessed$mechanism,
                     .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
  }
})

test_that("RMTL accepts the formal vector survival object and binds provenance", {
  capsule <- .dp_survival_vector_capsule(
    delayed_entry = FALSE, decoded = FALSE, gaussian = TRUE)
  result <- .dsvert_dp_survival_vector_result(
    capsule, "protected", "primary")
  result <- .dsvert_dp_survival_postprocess(result)
  class(result) <- c("ds.vertDPSurvival", "list")

  rmtl <- ds.vertDPRMTL(result, tau = 2)
  provenance <- attr(rmtl, "source_release_provenance")
  expect_identical(provenance$capsule_id, capsule$release$capsule_id)
  expect_identical(
    provenance$final_vector_root, capsule$release$final_vector_root)
  expect_identical(
    provenance$coordinate_order_sha256,
    capsule$release$coordinate_order_sha256)
  expect_identical(provenance$privacy_epoch, 1)
  expect_identical(provenance$noise_key_id, "noise-key-a")
  expect_identical(attr(rmtl, "additional_server_calls"), 0L)
  stale_gate <- result
  stale_gate$history_gate <- FALSE
  expect_error(ds.vertDPRMTL(stale_gate, tau = 2), "validated released")
})

.dp_survival_client_object <- function(delayed_entry = FALSE) {
  release <- .dp_survival_client_release(delayed_entry)
  result <- .dsvert_dp_survival_postprocess(release)
  class(result) <- c("ds.vertDPSurvival", "list")
  result
}

.dp_survival_pure_views <- function(result) {
  list(
    kaplan_meier = ds.vertDPKaplanMeier(result),
    nelson_aalen = ds.vertDPNelsonAalen(result),
    cumulative_incidence = ds.vertDPCumulativeIncidence(result),
    rmst = ds.vertDPRMST(result, c(1, 2, 4)),
    rmtl = ds.vertDPRMTL(result, c(1, 2, 4)),
    quantile = ds.vertDPSurvivalQuantile(result, c(0.25, 0.5, 0.75)),
    median = ds.vertDPMedianSurvival(result),
    survival_contrast = ds.vertDPSurvivalContrast(
      result, result, "comparison", "reference"),
    rmst_contrast = ds.vertDPRMSTContrast(
      result, result, c(1, 2, 4), "comparison", "reference"))
}

test_that("survival uses one no-lifetime Synopsis for K=2,3,5", {
  baseline <- .dp_survival_client_object(FALSE)
  expected_views <- .dp_survival_pure_views(baseline)
  legacy_fields <- c(
    "capsule_id", "privacy_epoch", "privacy_epochs", "noise_key_id",
    "noise_key_ids", "history_gate", "request_limit", "operation_limit",
    "lifetime_budget", "lifetime_composition", "privacy_accountant",
    "release_instance", "release_instance_id", "allocation_certificate")
  legacy_names <- function(value) {
    paths <- names(unlist(value, recursive = TRUE, use.names = TRUE))
    intersect(legacy_fields, c(names(value), sub("^.*[.]", "", paths)))
  }

  for (k in c(2L, 3L, 5L)) {
    fixture <- .dp_survival_synopsis_run(k = k)
    calls <- new.env(parent = emptyenv())
    calls$synopsis <- calls$capsule <- 0L
    result <- testthat::with_mocked_bindings(
      .dsvert_dp_survival_impl(
        "protected", "primary", "site_a", fixture$conns,
        function(...) stop("raw DSI call", call. = FALSE)),
      .dsvert_dp_synopsis_vector_run = function(
          datasources, status = NULL, .aggregate) {
        calls$synopsis <- calls$synopsis + 1L
        expect_identical(datasources, fixture$conns)
        expect_null(status)
        fixture[c("release", "layout", "status", "manifest_bundle")]
      },
      .dsvert_dp_capsule_vector_run = function(...) {
        calls$capsule <- calls$capsule + 1L
        stop("legacy capsule runner reached", call. = FALSE)
      },
      .dsvert_aggregate_strict = function(...) {
        stop("legacy aggregate route reached", call. = FALSE)
      },
      .dsvert_fanout_by_site = function(...) {
        stop("legacy fanout route reached", call. = FALSE)
      }, .package = "dsVertClient")

    expect_identical(calls$synopsis, 1L, info = paste("K =", k))
    expect_identical(calls$capsule, 0L, info = paste("K =", k))
    expect_equal(result$curve, baseline$curve, tolerance = 0)
    expect_equal(result$cumulative_incidence, baseline$cumulative_incidence,
                 tolerance = 0)
    expect_identical(result$causes, c("A", "B"))
    expect_identical(result$artifact_key, fixture$release$artifact_key)
    expect_identical(result$execution_id, fixture$release$execution_id)
    expect_identical(result$contract_sha256, fixture$release$contract_sha256)
    expect_identical(result$attempt_sha256, fixture$release$attempt_sha256)
    expect_identical(result$source_contract_sha256,
                     fixture$release$source_contract_sha256)
    expect_identical(result$result_set_sha256,
                     fixture$release$result_set_sha256)
    expect_identical(result$final_vector_root,
                     fixture$release$final_vector_root)
    expect_identical(result$release_provenance,
                     fixture$release$signed_provenance)
    expect_true(result$privacy$unlimited_replay)
    expect_length(legacy_names(result), 0L)

    views <- testthat::with_mocked_bindings(
      .dp_survival_pure_views(result),
      .dsvert_dp_synopsis_vector_run = function(...) {
        stop("Synopsis rerun from pure postprocessing", call. = FALSE)
      },
      .dsvert_dp_capsule_vector_run = function(...) {
        stop("capsule rerun from pure postprocessing", call. = FALSE)
      },
      .dsvert_aggregate_strict = function(...) {
        stop("aggregate from pure postprocessing", call. = FALSE)
      },
      .dsvert_fanout_by_site = function(...) {
        stop("fanout from pure postprocessing", call. = FALSE)
      }, .package = "dsVertClient")
    for (name in names(views)) {
      expect_equal(
        as.matrix(views[[name]]), as.matrix(expected_views[[name]]),
        tolerance = 0, info = paste(name, "K =", k))
    }
    provenances <- list(
      attr(views$rmst, "source_release_provenance"),
      attr(views$rmtl, "source_release_provenance"),
      attr(views$quantile, "source_release_provenance"),
      attr(views$median, "source_release_provenance"),
      attr(views$survival_contrast, "source_release_provenance")$comparison,
      attr(views$rmst_contrast, "source_release_provenance")$comparison)
    for (provenance in provenances) {
      expect_identical(provenance$artifact_key, result$artifact_key)
      expect_identical(provenance$release_provenance,
                       result$release_provenance)
      expect_length(legacy_names(provenance), 0L)
    }
  }
})

test_that("survival Synopsis replay is byte-identical", {
  fixture <- .dp_survival_synopsis_run(k = 5L, delayed_entry = TRUE)
  calls <- new.env(parent = emptyenv())
  calls$synopsis <- calls$capsule <- 0L
  evaluate <- function() testthat::with_mocked_bindings(
    .dsvert_dp_survival_impl(
      "protected", "primary", "site_a", fixture$conns,
      function(...) stop("raw DSI call", call. = FALSE)),
    .dsvert_dp_synopsis_vector_run = function(...) {
      calls$synopsis <- calls$synopsis + 1L
      fixture[c("release", "layout", "status", "manifest_bundle")]
    },
    .dsvert_dp_capsule_vector_run = function(...) {
      calls$capsule <- calls$capsule + 1L
      stop("legacy capsule runner reached", call. = FALSE)
    }, .package = "dsVertClient")
  first <- evaluate()
  second <- evaluate()
  expect_identical(serialize(second, NULL, version = 3L),
                   serialize(first, NULL, version = 3L))
  expect_identical(calls$synopsis, 2L)
  expect_identical(calls$capsule, 0L)
})

test_that("survival rejects detached Synopsis bindings before block access", {
  fixture <- .dp_survival_synopsis_run(k = 3L)
  bindings <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  for (field in bindings) {
    tampered <- fixture
    tampered$release[[field]] <- strrep("0", 64L)
    block_calls <- postprocess_calls <- 0L
    expect_error(testthat::with_mocked_bindings(
      .dsvert_dp_survival_impl(
        "protected", "primary", "site_a", tampered$conns,
        function(...) NULL),
      .dsvert_dp_synopsis_vector_run = function(...) {
        tampered[c("release", "layout", "status", "manifest_bundle")]
      },
      .dsvert_dp_capsule_vector_run = function(...) {
        stop("legacy capsule runner reached", call. = FALSE)
      },
      .dsvert_dp_capsule_single_block = function(...) {
        block_calls <<- block_calls + 1L
        stop("block access reached", call. = FALSE)
      },
      .dsvert_dp_survival_postprocess = function(...) {
        postprocess_calls <<- postprocess_calls + 1L
        stop("postprocess reached", call. = FALSE)
      }, .package = "dsVertClient"), "provenance is detached", info = field)
    expect_identical(block_calls, 0L, info = field)
    expect_identical(postprocess_calls, 0L, info = field)
  }
})

test_that("survival Synopsis rejects missing ambiguous and detached blocks", {
  fixture <- .dp_survival_synopsis_run(k = 3L)
  missing <- fixture
  missing$layout$blocks[["survival_artifacts::primary"]] <- NULL
  ambiguous <- fixture
  ambiguous$layout$blocks[["survival_artifacts::duplicate"]] <-
    ambiguous$layout$blocks[["survival_artifacts::primary"]]
  detached <- fixture
  detached$layout$blocks[[
    "survival_artifacts::primary"]]$descriptor$causes <- c("B", "A")

  candidates <- list(
    missing = list(run = missing, message = "does not contain exactly one"),
    ambiguous = list(
      run = ambiguous, message = "does not contain exactly one"),
    detached = list(run = detached, message = "coordinate contract"))
  for (name in names(candidates)) {
    candidate <- candidates[[name]]
    postprocess_calls <- 0L
    expect_error(testthat::with_mocked_bindings(
      .dsvert_dp_survival_impl(
        "protected", "primary", "site_a", candidate$run$conns,
        function(...) NULL),
      .dsvert_dp_synopsis_vector_run = function(...) {
        candidate$run[c("release", "layout", "status", "manifest_bundle")]
      },
      .dsvert_dp_capsule_vector_run = function(...) {
        stop("legacy capsule runner reached", call. = FALSE)
      },
      .dsvert_dp_survival_postprocess = function(...) {
        postprocess_calls <<- postprocess_calls + 1L
        stop("postprocess reached", call. = FALSE)
      }, .package = "dsVertClient"), candidate$message, info = name)
    expect_identical(postprocess_calls, 0L, info = name)
  }
})

test_that("KM, Nelson-Aalen and CIF match the no-noise central oracle", {
  result <- .dp_survival_client_object(FALSE)
  expect_identical(
    result$status,
    "fixed_public_clamp_applied_preclamp_state_not_released")
  expect_equal(result$curve$at_risk_dp, c(4, 3, 2, 1))
  expect_equal(result$curve$event_dp, c(1, 0, 1, 1))
  expect_equal(result$curve$censor_dp, c(0, 1, 0, 0))
  expect_equal(result$curve$kaplan_meier, c(0.75, 0.75, 0.375, 0))
  expect_equal(result$curve$nelson_aalen, c(0.25, 0.25, 0.75, 1.75))
  expect_equal(result$cumulative_incidence[, "A"],
               c(0.25, 0.25, 0.25, 0.625))
  expect_equal(result$cumulative_incidence[, "B"],
               c(0, 0, 0.375, 0.375))
  expect_equal(
    result$curve$kaplan_meier + rowSums(result$cumulative_incidence),
    rep(1, 4L), tolerance = 1e-14)
  expect_true(all(
    result$curve$kaplan_meier >=
      result$curve$kaplan_meier_mechanism_lower_95 &
    result$curve$kaplan_meier <=
      result$curve$kaplan_meier_mechanism_upper_95))
  expect_true(all(
    result$curve$nelson_aalen >=
      result$curve$nelson_aalen_mechanism_lower_95 &
    result$curve$nelson_aalen <=
      result$curve$nelson_aalen_mechanism_upper_95))
  expect_equal(
    result$mechanism_band_additional_privacy_cost,
    c(epsilon = 0, delta = 0))
  expect_identical(result$mechanism_band_additional_server_calls, 0L)
  expect_false(
    result$discretisation_error$included_in_mechanism_bands)
  expect_identical(result$discretisation_error$status, "not_quantified")
})

test_that("delayed entry uses noisy flow without invalid hazards", {
  result <- .dp_survival_client_object(TRUE)
  expect_equal(result$curve$at_risk_dp, c(2, 2, 2, 1))
  expect_true(all(diff(result$curve$kaplan_meier) <= 0))
  expect_true(all(diff(result$curve$nelson_aalen) >= 0))
  expect_true(all(result$curve$kaplan_meier >= 0 &
                  result$curve$kaplan_meier <= 1))
  expect_true(all(apply(result$cumulative_incidence, 2L, function(x) {
    all(diff(x) >= 0) && all(x >= 0 & x <= 1)
  })))
  expect_equal(
    result$curve$kaplan_meier + rowSums(result$cumulative_incidence),
    rep(1, 4L), tolerance = 1e-14)
  expect_identical(result$l1_sensitivity, 2L)
  expect_identical(result$coordinate_count, 17L)
})

test_that("oversized noisy cause hazards are jointly normalized", {
  normalized <- .dsvert_dp_survival_normalize_hazards(matrix(
    c(0.8, 0.7, 0.2, 0.3), nrow = 2L, byrow = TRUE,
    dimnames = list(NULL, c("A", "B"))))
  expect_equal(rowSums(normalized$cause_hazard), c(1, 0.5))
  expect_equal(normalized$all_hazard, c(1, 0.5))
  expect_true(all(normalized$cause_hazard >= 0 &
                  normalized$cause_hazard <= 1))

  curves <- .dsvert_dp_survival_curves(
    time_grid = 1:2, entry_counts = c(1, 0),
    exit_counts = matrix(
      c(0, 80, 70, 0, 2, 3), nrow = 2L, byrow = TRUE,
      dimnames = list(NULL, c("0", "A", "B"))),
    causes = c("A", "B"))
  total_cif <- rowSums(curves$cumulative_incidence)
  expect_true(all(total_cif <= 1 + 1e-14))
  expect_equal(curves$survival + total_cif, rep(1, 2L),
               tolerance = 1e-14)
  expect_true(all(diff(curves$survival) <= 0))
  expect_true(all(apply(curves$cumulative_incidence, 2L, function(value) {
    all(diff(value) >= 0)
  })))
})

test_that("RMST is zero-cost fixed-grid post-processing with valid bands", {
  release <- .dp_survival_client_release(delayed_entry = FALSE)
  result <- .dsvert_dp_survival_postprocess(release)
  class(result) <- c("ds.vertDPSurvival", "list")

  rmst <- ds.vertDPRMST(result, tau = c(0.5, 1, 1.5, 4))
  expect_s3_class(rmst, "ds.vertDPRMST")
  expect_equal(rmst$rmst, c(0.5, 1, 1.375, 2.875), tolerance = 1e-14)
  expect_true(all(rmst$rmst_mechanism_lower_95 <= rmst$rmst))
  expect_true(all(rmst$rmst <= rmst$rmst_mechanism_upper_95))
  expect_identical(
    attr(rmst, "additional_privacy_cost"), c(epsilon = 0, delta = 0))
  expect_identical(attr(rmst, "additional_server_calls"), 0L)
  expect_match(attr(rmst, "integration_rule"), "left-continuous")
  expect_match(attr(rmst, "mechanism_band_scope"), "sampling uncertainty")
  expect_match(attr(rmst, "grid_error_scope"), "discretisation")

  default <- ds.vertDPRMST(result)
  expect_equal(default$tau, result$time_upper_bound)
  expect_equal(default$rmst, 2.875, tolerance = 1e-14)
})

test_that("RMST validates public restriction times and release provenance", {
  release <- .dp_survival_client_release(delayed_entry = TRUE)
  result <- .dsvert_dp_survival_postprocess(release)
  class(result) <- c("ds.vertDPSurvival", "list")

  for (invalid in list(numeric(), NA_real_, Inf, -1, 0, 4.1, "4")) {
    expect_error(ds.vertDPRMST(result, invalid), "tau must contain")
  }
  expect_error(ds.vertDPRMST(unclass(result), 2),
               "released ds.vertDPSurvival")
})

test_that("RMTL is the exact zero-cost complement of released RMST", {
  result <- .dp_survival_client_object(FALSE)
  tau <- c(0.5, 1, 1.5, 4)
  rmst <- ds.vertDPRMST(result, tau)
  rmtl <- ds.vertDPRMTL(result, tau)
  width <- tau - result$time_lower_bound

  expect_s3_class(rmtl, "ds.vertDPRMTL")
  expect_s3_class(rmtl, "ds.vertDPRMST")
  for (column in names(rmst)) {
    expect_identical(rmtl[[column]], rmst[[column]])
  }
  expect_identical(rmtl$restriction_width, width)
  expect_identical(rmtl$rmtl, width - rmst$rmst)
  expect_identical(
    rmtl$rmtl_mechanism_lower_95,
    width - rmst$rmst_mechanism_upper_95)
  expect_identical(
    rmtl$rmtl_mechanism_upper_95,
    width - rmst$rmst_mechanism_lower_95)
  expect_equal(rmtl$rmtl, c(0, 0, 0.125, 1.125), tolerance = 1e-14)
  expect_true(all(rmtl$rmtl_mechanism_lower_95 <= rmtl$rmtl))
  expect_true(all(rmtl$rmtl <= rmtl$rmtl_mechanism_upper_95))
  expect_true(all(rmtl$rmtl_mechanism_lower_95 >= -1e-14))
  expect_true(all(rmtl$rmtl_mechanism_upper_95 <= width + 1e-14))
  for (attribute in c(
      "uncertainty_scope", "mechanism_band_scope",
      "mechanism_band_tightness", "grid_error_scope", "integration_rule",
      "additional_privacy_cost", "additional_server_calls",
      "source_release_provenance", "statistical_inference",
      "discretisation_error")) {
    expect_identical(attr(rmtl, attribute), attr(rmst, attribute))
  }
  expect_match(attr(rmtl, "complement_identity"),
               "restriction_width - RMST", fixed = TRUE)
  expect_match(attr(rmtl, "statistical_inference"),
               "no sampling confidence", ignore.case = TRUE)
  expect_match(attr(rmtl, "grid_error_scope"), "not included")
})

test_that("RMTL uses the released lower bound rather than assuming time zero", {
  release <- .dp_survival_client_release(FALSE)
  release$time_lower_bound <- 0.25
  result <- .dsvert_dp_survival_postprocess(release)
  class(result) <- c("ds.vertDPSurvival", "list")
  tau <- c(0.5, 1.5, 4)

  rmst <- ds.vertDPRMST(result, tau)
  rmtl <- ds.vertDPRMTL(result, tau)
  width <- tau - 0.25
  expect_identical(rmtl$restriction_width, width)
  expect_identical(rmtl$rmtl, width - rmst$rmst)
  expect_false(isTRUE(all.equal(rmtl$rmtl, tau - rmst$rmst)))
})

test_that("survival algebra rejects inconsistent post-release tampering", {
  result <- .dp_survival_client_object(FALSE)
  result$curve$kaplan_meier[[1L]] <-
    result$curve$kaplan_meier[[1L]] / 2

  expect_error(ds.vertDPRMST(result, 2), "validated released")
  expect_error(ds.vertDPRMTL(result, 2), "validated released")
})

test_that("RMTL performs no DSI work after the survival release", {
  result <- .dp_survival_client_object(FALSE)
  fail <- function(...) stop("unexpected DSI call", call. = FALSE)

  rmtl <- testthat::with_mocked_bindings(
    ds.vertDPRMTL(result, c(1, 4)),
    .dsvert_dp_datasources = fail,
    .dsvert_dp_capsule_vector_run = fail,
    .dsvert_aggregate_strict = fail,
    .dsvert_fanout_by_site = fail,
    .package = "dsVertClient")

  expect_identical(attr(rmtl, "additional_server_calls"), 0L)
  expect_identical(
    attr(rmtl, "additional_privacy_cost"), c(epsilon = 0, delta = 0))
  expect_identical(
    attr(rmtl, "source_release_provenance")$analysis_id,
    result$analysis_id)
})

test_that("survival contrasts preserve only defensible joint coverage", {
  capsule <- .dp_survival_vector_capsule(FALSE)
  comparison <- .dsvert_dp_survival_vector_result(
    capsule, "protected", "primary")
  comparison <- .dsvert_dp_survival_postprocess(comparison)
  class(comparison) <- c("ds.vertDPSurvival", "list")

  same <- ds.vertDPSurvivalContrast(
    comparison, comparison, "treated", "control")
  expect_s3_class(same, "ds.vertDPSurvivalContrast")
  expect_equal(same$survival_difference, rep(0, nrow(same)))
  positive_reference <- same$reference_survival > 0
  expect_equal(same$survival_ratio[positive_reference],
               rep(1, sum(positive_reference)))
  expect_true(all(is.na(same$survival_ratio[!positive_reference])))
  expect_true(all(
    same$survival_difference_mechanism_lower <= 0 &
      same$survival_difference_mechanism_upper >= 0))
  expect_true(all(
    same$survival_ratio_mechanism_lower <= 1 &
      same$survival_ratio_mechanism_upper >= 1))
  expect_identical(attr(same, "joint_mechanism_confidence"), 0.95)
  expect_identical(attr(same, "joint_event"),
                   "same_signed_survival_artifact")

  distinct <- comparison
  distinct$final_vector_root <- strrep("c", 64L)
  contrast <- ds.vertDPSurvivalContrast(comparison, distinct)
  expect_equal(attr(contrast, "joint_mechanism_confidence"), 0.9,
               tolerance = 1e-15)
  expect_identical(attr(contrast, "joint_event"),
                   "bonferroni_across_distinct_releases")
  expect_match(attr(contrast, "statistical_inference"),
               "no sampling confidence", ignore.case = TRUE)
  expect_identical(attr(contrast, "additional_privacy_cost"),
                   c(epsilon = 0, delta = 0))
  expect_identical(attr(contrast, "additional_server_calls"), 0L)
})

test_that("survival contrasts type zero-denominator ratios", {
  comparison <- .dp_survival_client_object(FALSE)
  reference <- comparison

  contrast <- ds.vertDPSurvivalContrast(comparison, reference)
  zero <- reference$curve$kaplan_meier == 0
  expect_true(any(zero))
  expect_true(all(is.na(contrast$survival_ratio[zero])))
  expect_true(all(
    contrast$survival_ratio_status[zero] ==
      "reference_point_zero_not_estimable"))
  expect_true(all(is.infinite(
    contrast$survival_ratio_mechanism_upper[
      contrast$reference_mechanism_interval_includes_zero &
        contrast$comparison_survival_mechanism_upper > 0])))
})

test_that("RMST contrasts are exact zero-call post-processing", {
  comparison <- .dp_survival_client_object(FALSE)
  reference <- comparison
  result <- ds.vertDPRMSTContrast(
    comparison, reference, c(1, 2, 4), "exposed", "unexposed")

  expect_s3_class(result, "ds.vertDPRMSTContrast")
  expect_equal(result$rmst_difference, rep(0, 3L))
  expect_equal(result$rmst_ratio, rep(1, 3L))
  expect_true(all(result$rmst_difference_mechanism_lower <= 0))
  expect_true(all(result$rmst_difference_mechanism_upper >= 0))
  expect_true(all(result$rmst_ratio_mechanism_lower <= 1))
  expect_true(all(result$rmst_ratio_mechanism_upper >= 1))
  expect_match(attr(result, "grid_error_scope"), "not included")
  expect_identical(attr(result, "additional_privacy_cost"),
                   c(epsilon = 0, delta = 0))
  expect_identical(attr(result, "additional_server_calls"), 0L)
})

test_that("survival contrast rejects incompatible releases and bad labels", {
  comparison <- .dp_survival_client_object(FALSE)
  incompatible <- comparison
  incompatible$time_grid[[2L]] <- 2.5
  incompatible$time_upper_bound <- 4
  incompatible <- .dsvert_dp_survival_postprocess(incompatible)
  class(incompatible) <- c("ds.vertDPSurvival", "list")

  expect_error(
    ds.vertDPSurvivalContrast(comparison, incompatible),
    "same signed public time grid")
  expect_error(
    ds.vertDPRMSTContrast(comparison, incompatible, 2),
    "same signed public time grid")
  expect_error(
    ds.vertDPSurvivalContrast(comparison, comparison, "same", "same"),
    "distinct non-empty strings")
  expect_error(
    ds.vertDPRMSTContrast(comparison, comparison, 2, NA_character_, "b"),
      "distinct non-empty strings")
})

test_that("survival contrasts reject incompatible signed estimands", {
  comparison <- .dp_survival_client_object(FALSE)
  incompatible <- list()

  incompatible$censor <- comparison
  incompatible$censor$censor_level <- "censored"
  incompatible$censor <- .dsvert_dp_survival_postprocess(
    incompatible$censor)
  class(incompatible$censor) <- c("ds.vertDPSurvival", "list")

  incompatible$causes <- comparison
  incompatible$causes$causes <- c("C", "D")
  incompatible$causes <- .dsvert_dp_survival_postprocess(
    incompatible$causes)
  class(incompatible$causes) <- c("ds.vertDPSurvival", "list")

  incompatible$delayed_entry <- .dp_survival_client_object(TRUE)

  for (reference in incompatible) {
    expect_error(
      ds.vertDPSurvivalContrast(comparison, reference),
      "same signed public time grid")
    expect_error(
      ds.vertDPRMSTContrast(comparison, reference, 2),
      "same signed public time grid")
  }
})

test_that("copied identity metadata cannot bypass survival Bonferroni", {
  capsule <- .dp_survival_vector_capsule(FALSE)
  comparison <- .dsvert_dp_survival_vector_result(
    capsule, "protected", "primary")
  comparison <- .dsvert_dp_survival_postprocess(comparison)
  class(comparison) <- c("ds.vertDPSurvival", "list")

  reference <- comparison
  first_event <- length(reference$time_grid) + 1L
  reference$histogram[[first_event]] <-
    reference$histogram[[first_event]] + 1
  reference <- .dsvert_dp_survival_postprocess(reference)
  class(reference) <- c("ds.vertDPSurvival", "list")

  expect_false(identical(comparison$histogram, reference$histogram))
  identity_fields <- c(
    "analysis_id", "analysis_version", "server", "capsule_id",
    "final_vector_root", "coordinate_order_sha256", "privacy_epoch")
  expect_true(all(vapply(identity_fields, function(field) {
    identical(comparison[[field]], reference[[field]])
  }, logical(1L))))

  survival <- ds.vertDPSurvivalContrast(comparison, reference)
  rmst <- ds.vertDPRMSTContrast(comparison, reference, 2)
  for (result in list(survival, rmst)) {
    expect_identical(attr(result, "joint_event"),
                     "bonferroni_across_distinct_releases")
    expect_equal(attr(result, "joint_mechanism_confidence"), 0.9,
                 tolerance = 1e-15)
  }
})

test_that("survival contrast implementations contain no DSI route", {
  source <- paste(
    readLines(.dsvert_client_source_file(
                "ds.vertDPSurvivalContrast.R"),
              warn = FALSE),
    collapse = "\n")
  expect_false(grepl(
    "datashield\\.aggregate|\\.dsvert_fanout|\\.dsvert_dp_datasources|DSI::",
    source))
})

test_that("fixed-grid survival quantiles and median match the central curve", {
  result <- .dp_survival_client_object(FALSE)
  result$accuracy_simultaneous_95_abs <- 0
  result <- .dsvert_dp_survival_postprocess(result)
  class(result) <- c("ds.vertDPSurvival", "list")

  quantiles <- ds.vertDPSurvivalQuantile(result, c(0.25, 0.5, 0.75))
  median <- ds.vertDPMedianSurvival(result)

  expect_s3_class(quantiles, "ds.vertDPSurvivalQuantile")
  expect_s3_class(median, "ds.vertDPMedianSurvival")
  expect_identical(quantiles$probability, c(0.25, 0.5, 0.75))
  expect_identical(quantiles$survival_probability, c(0.75, 0.5, 0.25))
  expect_identical(quantiles$quantile, c(1, 3, 4))
  expect_identical(
    quantiles$quantile_mechanism_lower_95, quantiles$quantile)
  expect_identical(
    quantiles$quantile_mechanism_upper_95, quantiles$quantile)
  expect_identical(median$quantile, 3)
  expect_identical(median$probability, 0.5)
  expect_true(all(quantiles$status == "ok"))
  expect_false(any(quantiles$mechanism_region_includes_beyond_grid))
  expect_identical(attr(quantiles, "additional_server_calls"), 0L)
  expect_identical(
    attr(quantiles, "additional_privacy_cost"),
    c(epsilon = 0, delta = 0))
  expect_identical(
    attr(quantiles, "source_release_provenance")$analysis_id,
    result$analysis_id)
  expect_match(attr(quantiles, "grid_quantile_definition"),
               "first public endpoint")
})

test_that("survival quantiles type horizons that do not reach the target", {
  release <- .dp_survival_client_release(FALSE)
  release$histogram <- c(c(1, 1, 1, 1), rep(0, 8), 0)
  release$accuracy_simultaneous_95_abs <- 0
  result <- .dsvert_dp_survival_postprocess(release)
  class(result) <- c("ds.vertDPSurvival", "list")

  quantile <- ds.vertDPSurvivalQuantile(result, 0.5)
  expect_true(is.na(quantile$quantile))
  expect_identical(quantile$point_status, "not_reached_by_grid_end")
  expect_true(is.infinite(quantile$quantile_mechanism_lower_95))
  expect_true(is.infinite(quantile$quantile_mechanism_upper_95))
  expect_identical(
    quantile$status, "not_reached_by_any_curve_in_mechanism_band")
  expect_false(quantile$mechanism_region_has_finite_values)
  expect_true(quantile$mechanism_region_includes_beyond_grid)
})

test_that("survival quantile band inversion covers every monotone curve", {
  grid <- c(1, 2, 3)
  lower <- c(0.5, 0.25, 0)
  upper <- c(1, 0.75, 0.5)
  levels <- seq(0, 1, by = 0.25)
  candidates <- expand.grid(rep(list(levels), length(grid)))
  candidates <- candidates[apply(candidates, 1L, function(value) {
    all(diff(as.numeric(value)) <= 0) &&
      all(as.numeric(value) >= lower) && all(as.numeric(value) <= upper)
  }), , drop = FALSE]
  expect_gt(nrow(candidates), 0L)

  for (probability in c(0.25, 0.5, 0.75)) {
    region <- .dsvert_dp_survival_quantile_region(
      grid, lower, upper, probability)
    candidate_quantiles <- apply(candidates, 1L, function(value) {
      .dsvert_dp_survival_quantile_endpoint(
        grid, as.numeric(value), probability)
    })
    expect_true(all(candidate_quantiles >= region[["lower"]]))
    expect_true(all(candidate_quantiles <= region[["upper"]]))
  }
})

test_that("survival quantiles reject invalid inputs, tampering, and DSI", {
  result <- .dp_survival_client_object(FALSE)
  for (bad in list(numeric(), 0, 1, -0.1, 1.1, NA_real_, NaN, Inf,
                   "0.5")) {
    expect_error(ds.vertDPSurvivalQuantile(result, bad),
                 "probabilities")
  }

  tampered <- result
  tampered$curve$kaplan_meier[[1L]] <- 0
  expect_error(ds.vertDPSurvivalQuantile(tampered, 0.5),
               "validated released")

  fail <- function(...) stop("unexpected DSI call", call. = FALSE)
  quantile <- testthat::with_mocked_bindings(
    ds.vertDPSurvivalQuantile(result, c(0.25, 0.5)),
    .dsvert_dp_datasources = fail,
    .dsvert_dp_capsule_vector_run = fail,
    .dsvert_aggregate_strict = fail,
    .dsvert_fanout_by_site = fail,
    .package = "dsVertClient")
  expect_identical(attr(quantile, "additional_server_calls"), 0L)
})

test_that("random survival histograms always yield coherent finite curves", {
  withr::local_seed(20260803)
  valid <- logical(1000L)
  for (iteration in seq_along(valid)) {
    time_count <- sample.int(12L, 1L)
    cause_count <- sample.int(4L, 1L)
    causes <- paste0("cause", seq_len(cause_count))
    exit_counts <- matrix(
      sample(0:100, time_count * (cause_count + 1L), replace = TRUE),
      nrow = time_count)
    entry_counts <- if (sample(c(FALSE, TRUE), 1L)) {
      sample(0:100, time_count, replace = TRUE)
    } else NULL
    curves <- .dsvert_dp_survival_curves(
      seq_len(time_count), entry_counts, exit_counts, causes)
    radius <- sample(0:20, 1L)
    bands <- .dsvert_dp_survival_mechanism_bands(
      seq_len(time_count),
      if (is.null(entry_counts)) NULL else pmax(0, entry_counts - radius),
      if (is.null(entry_counts)) NULL else entry_counts + radius,
      pmax(exit_counts - radius, 0), exit_counts + radius, causes)
    total_incidence <- rowSums(curves$cumulative_incidence)
    valid[[iteration]] <-
      all(is.finite(unlist(curves, use.names = FALSE))) &&
      all(curves$at_risk >= 0) &&
      all(curves$cause_hazard >= 0) &&
      all(curves$all_cause_hazard >= 0 &
            curves$all_cause_hazard <= 1 + 1e-14) &&
      all(curves$survival >= 0 & curves$survival <= 1) &&
      all(diff(curves$survival) <= 1e-14) &&
      all(diff(curves$cumulative_hazard) >= -1e-14) &&
      all(apply(curves$cumulative_incidence, 2L, function(value) {
        all(diff(value) >= -1e-14)
      })) &&
      all(total_incidence >= 0 & total_incidence <= 1 + 1e-12) &&
      all(abs(curves$survival + total_incidence - 1) <= 1e-12) &&
      all(curves$survival >= bands$kaplan_meier$lower - 1e-12 &
            curves$survival <= bands$kaplan_meier$upper + 1e-12) &&
      all(curves$cumulative_hazard >=
            bands$nelson_aalen$lower - 1e-12 &
            curves$cumulative_hazard <=
            bands$nelson_aalen$upper + 1e-12) &&
      all(curves$cumulative_incidence >=
            bands$cumulative_incidence$lower - 1e-12 &
            curves$cumulative_incidence <=
            bands$cumulative_incidence$upper + 1e-12)
  }
  expect_true(all(valid))
})

test_that("survival mechanism bands exhaust every small histogram box", {
  configurations <- list(
    list(delayed = FALSE, causes = "A"),
    list(delayed = FALSE, causes = c("A", "B")),
    list(delayed = TRUE, causes = "A"))
  covered <- invariants <- TRUE
  failure <- ""
  checked <- 0L
  for (configuration_index in seq_along(configurations)) {
    configuration <- configurations[[configuration_index]]
    time_count <- 2L
    exit_length <- time_count * (length(configuration$causes) + 1L)
    entry_length <- if (configuration$delayed) time_count else 0L
    coordinate_count <- entry_length + exit_length
    observed_boxes <- expand.grid(
      rep(list(0:1), coordinate_count), KEEP.OUT.ATTRS = FALSE)
    for (box_index in seq_len(nrow(observed_boxes))) {
      observed <- as.numeric(observed_boxes[box_index, ])
      lower <- pmax(0, observed - 1)
      upper <- observed + 1
      cursor <- 1L
      entry_lower <- entry_upper <- NULL
      if (configuration$delayed) {
        entry_lower <- lower[cursor:(cursor + time_count - 1L)]
        entry_upper <- upper[cursor:(cursor + time_count - 1L)]
        cursor <- cursor + time_count
      }
      exit_lower <- matrix(
        lower[cursor:(cursor + exit_length - 1L)],
        nrow = time_count)
      exit_upper <- matrix(
        upper[cursor:(cursor + exit_length - 1L)],
        nrow = time_count)
      bands <- .dsvert_dp_survival_mechanism_bands(
        1:2, entry_lower, entry_upper,
        exit_lower, exit_upper, configuration$causes)

      invariants <- invariants &&
        all(bands$kaplan_meier$lower >= 0) &&
        all(bands$kaplan_meier$upper <= 1) &&
        all(bands$kaplan_meier$lower <= bands$kaplan_meier$upper) &&
        all(diff(bands$kaplan_meier$lower) <= 1e-14) &&
        all(diff(bands$kaplan_meier$upper) <= 1e-14) &&
        all(bands$nelson_aalen$lower >= 0) &&
        all(bands$nelson_aalen$lower <= bands$nelson_aalen$upper) &&
        all(diff(bands$nelson_aalen$lower) >= -1e-14) &&
        all(diff(bands$nelson_aalen$upper) >= -1e-14) &&
        all(bands$cumulative_incidence$lower >= 0) &&
        all(bands$cumulative_incidence$upper <= 1) &&
        all(bands$cumulative_incidence$lower <=
              bands$cumulative_incidence$upper) &&
        all(apply(
          bands$cumulative_incidence$lower, 2L,
          function(value) all(diff(value) >= -1e-14))) &&
        all(apply(
          bands$cumulative_incidence$upper, 2L,
          function(value) all(diff(value) >= -1e-14))) &&
        isTRUE(all.equal(
          bands$total_cumulative_incidence$lower,
          1 - bands$kaplan_meier$upper, tolerance = 1e-14)) &&
        isTRUE(all.equal(
          bands$total_cumulative_incidence$upper,
          1 - bands$kaplan_meier$lower, tolerance = 1e-14))

      candidates <- expand.grid(
        Map(seq.int, lower, upper), KEEP.OUT.ATTRS = FALSE)
      for (candidate_index in seq_len(nrow(candidates))) {
        candidate <- as.numeric(candidates[candidate_index, ])
        cursor <- 1L
        entry <- NULL
        if (configuration$delayed) {
          entry <- candidate[cursor:(cursor + time_count - 1L)]
          cursor <- cursor + time_count
        }
        exit <- matrix(
          candidate[cursor:(cursor + exit_length - 1L)],
          nrow = time_count)
        point <- .dsvert_dp_survival_curves(
          1:2, entry, exit, configuration$causes)
        total_incidence <- rowSums(point$cumulative_incidence)
        inside <-
          all(point$survival >= bands$kaplan_meier$lower - 1e-12 &
                point$survival <= bands$kaplan_meier$upper + 1e-12) &&
          all(point$cumulative_hazard >=
                bands$nelson_aalen$lower - 1e-12 &
                point$cumulative_hazard <=
                bands$nelson_aalen$upper + 1e-12) &&
          all(point$cumulative_incidence >=
                bands$cumulative_incidence$lower - 1e-12 &
                point$cumulative_incidence <=
                bands$cumulative_incidence$upper + 1e-12) &&
          all(total_incidence >=
                bands$total_cumulative_incidence$lower - 1e-12 &
                total_incidence <=
                bands$total_cumulative_incidence$upper + 1e-12)
        if (!inside) {
          covered <- FALSE
          failure <- paste(
            "configuration", configuration_index,
            "box", box_index, "candidate", candidate_index)
        }
        checked <- checked + 1L
      }
    }
  }
  expect_identical(checked, 31875L)
  expect_true(covered, info = failure)
  expect_true(invariants)
})

test_that("empty survival releases remain finite and explicitly typed", {
  for (delayed in c(FALSE, TRUE)) {
    release <- .dp_survival_client_release(delayed)
    release$histogram[] <- 0
    result <- .dsvert_dp_survival_postprocess(release)
    expect_identical(result$status, "dp_curve_empty_after_postprocessing")
    expect_true(all(result$curve$at_risk_dp == 0))
    expect_true(all(result$curve$kaplan_meier == 1))
    expect_true(all(result$curve$nelson_aalen == 0))
    expect_true(all(result$cumulative_incidence == 0))
    expect_false(anyNA(unlist(
      list(result$curve, result$cumulative_incidence),
      use.names = FALSE)))
  }
})

test_that("named methods are post-processing of the same release", {
  result <- .dp_survival_client_object(FALSE)
  km <- ds.vertDPKaplanMeier(result)
  na <- ds.vertDPNelsonAalen(result)
  cif <- ds.vertDPCumulativeIncidence(result)
  expect_equal(km$kaplan_meier, result$curve$kaplan_meier)
  expect_equal(na$nelson_aalen, result$curve$nelson_aalen)
  expect_equal(
    cif$cumulative_incidence,
    as.vector(result$cumulative_incidence))
  expect_equal(
    km$kaplan_meier_mechanism_lower_95,
    result$curve$kaplan_meier_mechanism_lower_95)
  expect_equal(
    na$nelson_aalen_mechanism_upper_95,
    result$curve$nelson_aalen_mechanism_upper_95)
  expect_equal(
    cif$cumulative_incidence_mechanism_lower_95,
    as.vector(result$cumulative_incidence_mechanism_lower_95))
  expect_match(attr(km, "uncertainty_scope"), "DP mechanism noise only")
  expect_match(attr(km, "mechanism_band_scope"), "grid")
  expect_match(attr(cif, "mechanism_band_tightness"), "rectangular")
  for (view in list(km, na, cif)) {
    expect_identical(
      attr(view, "additional_privacy_cost"),
      c(epsilon = 0, delta = 0))
    expect_identical(attr(view, "additional_server_calls"), 0L)
  }
  expect_match(result$statistical_inference, "no sampling confidence")
  expect_false(any(grepl("confidence|p_value|std_error",
                         names(result$curve), ignore.case = TRUE)))
  expect_error(
    ds.vertDPCumulativeIncidence(result, "unknown"),
    "released event cause")
})

test_that("repeated survival requests replay the one Synopsis vector", {
  fixture <- .dp_survival_synopsis_run(k = 2L)
  calls <- new.env(parent = emptyenv())
  calls$synopsis <- calls$capsule <- 0L
  evaluate <- function() {
    .dsvert_dp_survival_impl(
      "protected", "primary", "site_a", fixture$conns,
      function(...) stop("unexpected direct aggregate"))
  }
  result <- testthat::with_mocked_bindings(
    list(evaluate(), evaluate()),
    .dsvert_dp_synopsis_vector_run = function(
        datasources, status = NULL, .aggregate) {
      calls$synopsis <- calls$synopsis + 1L
      fixture[c("release", "layout", "status", "manifest_bundle")]
    },
    .dsvert_dp_capsule_vector_run = function(...) {
      calls$capsule <- calls$capsule + 1L
      stop("legacy capsule runner reached", call. = FALSE)
    },
    .package = "dsVertClient")
  expect_identical(result[[1L]], result[[2L]])
  expect_identical(calls$synopsis, 2L)
  expect_identical(calls$capsule, 0L)
  expect_identical(result[[1L]]$final_vector_root,
                   result[[2L]]$final_vector_root)
})

test_that("survival release rejects partial/raw DSI errors without retry", {
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
      .dsvert_dp_survival_impl(
        "protected", "primary", "site_a", conns, aggregate)
    }
    captured <- tryCatch(evaluate(), error = function(e) conditionMessage(e))
    expect_match(
      captured,
      "DataSHIELD transport failed during 'synopsis connection identity fan-out'",
      fixed = TRUE)
    expect_false(grepl("SECRET_REMOTE_DETAIL", captured, fixed = TRUE))
    expect_identical(attempts, 1L)
  }
})
