test_that("exported DP methods cannot reach retired DP endpoints", {
  namespace <- asNamespace("dsVertClient")
  roots <- sort(
    getNamespaceExports("dsVertClient")[
      startsWith(getNamespaceExports("dsVertClient"), "ds.vertDP")],
    method = "radix")
  retired <- c(
    "dsvertDPStatusDS", "dsvertDPCountDS", "dsvertDPContingencyDS",
    "dsvertDPMeanVarDS", "dsvertDPDescribeDS", "dsvertDPSurvivalDS",
    "dsvertJointDPPrepareDS", "dsvertJointDPCommitDS",
    "dsvertJointDPAuthorizeDS", "dsvertJointDPOpenDS",
    "dsvertJointDPResultReceiptDS", "dsvertJointDPDeliveryDS",
    "dsvertJointDPDeliveryContractDS", "dsvertJointDPCountReplayDS",
    "dsvertJointDPCountProposalDS", "dsvertJointDPCountSourceDS",
    "dsvertJointDPCountBackendPrepareDS",
    "dsvertJointDPCountBackendTokenDS", "dsvertJointDPCountStartDS",
    "dsvertJointDPCountResultDS", "dsvertJointDPCountFinalShareDS",
    "dsvertJointDPCountReleaseDS")

  offenders <- character()
  for (root in roots) {
    queue <- root
    seen <- character()
    while (length(queue)) {
      name <- queue[[1L]]
      queue <- queue[-1L]
      if (name %in% seen ||
          !exists(name, envir = namespace, inherits = FALSE)) {
        next
      }
      value <- get(name, envir = namespace, inherits = FALSE)
      if (!is.function(value)) next
      seen <- c(seen, name)
      source <- paste(
        deparse(body(value), width.cutoff = 500L), collapse = "\n")
      hits <- retired[vapply(
        retired, grepl, logical(1L), x = source, fixed = TRUE)]
      if (length(hits)) {
        offenders <- c(
          offenders, paste(root, name, hits, sep = " -> "))
      }
      globals <- tryCatch(
        unique(unlist(
          codetools::findGlobals(value, merge = FALSE),
          use.names = FALSE)),
        error = function(error) character())
      internal <- globals[vapply(globals, function(global) {
        exists(global, envir = namespace, inherits = FALSE) &&
          is.function(get(global, envir = namespace, inherits = FALSE))
      }, logical(1L))]
      queue <- unique(c(queue, setdiff(internal, seen)))
    }
  }

  expect_gt(length(roots), 0L)
  expect_false(
    length(offenders) > 0L, info = paste(offenders, collapse = "\n"))
})

test_that("retired local schema-6 validators are absent", {
  namespace <- asNamespace("dsVertClient")
  retired <- c(
    ".DSVERT_DP_TRANSCRIPT_POLICY_VALUE",
    ".dsvert_dp_validate_named_domains",
    ".dsvert_dp_validate_numeric_bounds",
    ".dsvert_dp_validate_datasets",
    ".dsvert_dp_expected_anchor_id",
    ".dsvert_dp_validate_rollback_policy",
    ".dsvert_dp_validate_noise_selection_policy",
    ".dsvert_dp_validate_policy",
    ".dsvert_dp_find_server",
    ".dsvert_dp_release_epsilon_is_valid",
    ".dsvert_dp_validate_release_common",
    ".dsvert_dp_sampler_metadata_is_valid",
    ".dsvert_dp_validate_noise_selection_certificate",
    ".dsvert_dp_selected_sampler_metadata_is_valid",
    ".dsvert_dp_validate_count_release",
    ".dsvert_dp_validate_table_release",
    ".dsvert_dp_validate_meanvar_release",
    ".dsvert_dp_validate_release",
    ".dsvert_dp_validate_describe_release",
    ".dsvert_dp_validate_survival_release",
    ".dsvert_dp_simultaneous_radius",
    ".dsvert_dp_gaussian_sigma_matches_calibration",
    ".dsvert_dp_meanvar_simultaneous_radii",
    ".dsvert_dp_meanvar_projection_interval",
    ".dsvert_dp_meanvar_mechanism_region",
    ".dsvert_dp_meanvar_postprocess")

  present <- retired[vapply(
    retired, exists, logical(1L), envir = namespace, inherits = FALSE)]
  expect_length(present, 0L)
})

test_that("DP calibration uses fixed per-capsule privacy parameters", {
  default_calibration <- ds.vertDPCalibrate(capsule_epsilon = 1)
  expect_identical(
    default_calibration$capsule_delta,
    .DSVERT_DP_DEFAULT_CAPSULE_DELTA)
  expect_identical(
    default_calibration$capsule_delta,
    7.888609052210118e-31)

  calibration <- ds.vertDPCalibrate(
    capsule_epsilon = c(1, 3), peer_count = 2L, sensitivity = 1,
    confidence = 0.95)
  expect_s3_class(calibration, "ds.vertDPCalibration")
  expect_equal(nrow(calibration), 2L)
  first_strict <- calibration[
    calibration$capsule_epsilon == 1, ]
  expect_equal(first_strict$capsule_epsilon, 1)
  expect_equal(first_strict$expected_absolute_error, 1)
  expect_equal(first_strict$rmse, sqrt(2))
  expect_equal(first_strict$error_radius,
               floor(-log1p(-0.95) + 0.5))
  expect_true(first_strict$sampler_supported)
  expect_match(first_strict$metric_guarantee, "planning approximations")
  strict <- calibration[calibration$capsule_epsilon == 1, ]
  balanced <- calibration[calibration$capsule_epsilon == 3, ]
  expect_true(all(balanced$expected_absolute_error <
                    strict$expected_absolute_error))
  expect_identical(unique(strict$expected_absolute_error), 1)
  expect_true(all(calibration$privacy_parameters_are_fixed))
  expect_identical(unique(calibration$operation_accounting), "none")
  expect_false(any(calibration$operation_limit))
  expect_false(any(calibration$request_limit))
  expect_false(any(calibration$history_can_deny_operation))
  expect_false(any(calibration$accuracy_depends_on_request_history))
  expect_false(any(grepl(
    "remaining|exhaust|budget|decay|queries_remaining",
    names(calibration), ignore.case = TRUE)))
  expect_identical(names(formals(ds.vertDPCalibrate)), c(
    "capsule_epsilon", "peer_count", "sensitivity", "confidence",
    "capsule_delta", "coordinate_count", "gaussian_l2_sensitivity",
    "objective"))
  for (retired in c(
      "decay", "release_indices", "composition_partitions",
      "total_epsilon", "total_delta")) {
    call <- as.call(c(
      list(quote(ds.vertDPCalibrate)), stats::setNames(list(1), retired)))
    expect_error(eval(call), "unused argument")
  }

  expect_error(ds.vertDPCalibrate(
    capsule_epsilon = 1, sensitivity = 1e308),
    "exactly representable integer")
  expect_error(ds.vertDPCalibrate(
    capsule_epsilon = 1, sensitivity = 1.5),
    "exactly representable integer")
  expect_error(ds.vertDPCalibrate(
    capsule_epsilon = 1,
    sensitivity = .Machine$double.xmin),
    "exactly representable integer")
  expect_error(ds.vertDPCalibrate(
    capsule_epsilon = 1, sensitivity = 2^53),
    "2\\^53 - 1")
  expect_silent(maximum <- ds.vertDPCalibrate(
    capsule_epsilon = 1, sensitivity = 2^53 - 1))
  expect_false(maximum$sampler_supported)

  expect_true(ds.vertDPCalibrate(
    capsule_epsilon = 2^40, peer_count = 2L,
    sensitivity = 1)$sampler_supported)
  expect_false(ds.vertDPCalibrate(
    capsule_epsilon = 2^40 + 1, peer_count = 2L,
    sensitivity = 1)$sampler_supported)

  expect_error(
    ds.vertDPCalibrate(capsule_epsilon = 1, peer_count = 1L),
    "at least two")
  k2 <- ds.vertDPCalibrate(
    capsule_epsilon = 1, peer_count = 2L, sensitivity = 1)
  k4 <- ds.vertDPCalibrate(
    capsule_epsilon = 1, peer_count = 4L, sensitivity = 1)
  expect_identical(k2$capsule_epsilon, k4$capsule_epsilon)
  expect_identical(k2$capsule_delta, k4$capsule_delta)
  expect_identical(k4$peer_count, 4L)
  expect_false(any(c(
    "composition_partitions", "local_lifetime_epsilon",
    "local_lifetime_delta") %in% names(k4)))
  selected <- ds.vertDPCalibrate(
    capsule_epsilon = 2, capsule_delta = 0.02,
    peer_count = 2L,
    sensitivity = 100, gaussian_l2_sensitivity = 1,
    coordinate_count = 64L)
  expect_identical(selected$selector,
                   "minimum_conservative_95_radius_v3")
  expect_identical(selected$selector_objective,
                   "simultaneous_95_abs")
  expect_identical(selected$utility_candidate, "gaussian")
  expect_identical(selected$preview_candidate, "gaussian")
  expect_true(is.na(selected$selected_candidate))
  expect_true(selected$gaussian_backend_deployed)
  expect_true(selected$gaussian_preview_supported)
  expect_identical(
    selected$selection_authority, "signed_server_capsule_manifest")
  expect_identical(
    selected$formal_backend, "signed_server_capsule_manifest")
  expect_match(selected$deployment_decision, "server_plan_authoritative")
  expect_equal(
    selected$gaussian_analytic_delta +
      selected$gaussian_implementation_delta_bound,
    selected$capsule_delta, tolerance = 0)
  expect_lt(selected$gaussian_selection_radius,
            selected$laplace_selection_radius)
  selected_90 <- ds.vertDPCalibrate(
    capsule_epsilon = 2, capsule_delta = 0.02,
    peer_count = 2L,
    sensitivity = 100, gaussian_l2_sensitivity = 1,
    coordinate_count = 64L, confidence = 0.90)
  selected_99 <- ds.vertDPCalibrate(
    capsule_epsilon = 2, capsule_delta = 0.02,
    peer_count = 2L,
    sensitivity = 100, gaussian_l2_sensitivity = 1,
    coordinate_count = 64L, confidence = 0.99)
  expect_identical(selected_90$selected_candidate,
                   selected$selected_candidate)
  expect_identical(selected_99$selected_candidate,
                   selected$selected_candidate)
  expect_identical(selected_90$laplace_selection_radius,
                   selected$laplace_selection_radius)
  expect_identical(selected_99$gaussian_selection_radius,
                   selected$gaussian_selection_radius)

  transfer_reference <- exp(
    log(64) + log(.DSVERT_DP_GAUSSIAN_TV_BOUND_PER_COORDINATE) +
      2 + log1p(exp(-2)))
  transfer_planner <-
    .dsvert_dp_gaussian_implementation_delta_bound(64L, 2)
  expect_gte(transfer_planner, transfer_reference)
  expect_lte(
    transfer_planner,
    transfer_reference * (1 + 128 * .Machine$double.eps))

  boundary <- ds.vertDPCalibrate(
    capsule_epsilon = 2,
    capsule_delta =
      .dsvert_dp_gaussian_implementation_delta_bound(64L, 2),
    peer_count = 2L,
    sensitivity = 100, gaussian_l2_sensitivity = 1,
    coordinate_count = 64L)
  expect_false(boundary$gaussian_sampler_supported)
  expect_true(is.na(boundary$selected_candidate))

})

test_that("calibration has no operation-history input or accuracy schedule", {
  calibration <- ds.vertDPCalibrate(
    capsule_epsilon = rep(1, 10000L), capsule_delta = 0)
  expect_equal(nrow(calibration), 10000L)
  expect_identical(unique(calibration$capsule_epsilon), 1)
  expect_identical(unique(calibration$capsule_delta), 0)
  expect_identical(unique(calibration$expected_absolute_error), 1)
  expect_identical(unique(calibration$operation_accounting), "none")
  expect_false(any(calibration$operation_limit))
  expect_false(any(calibration$request_limit))
  expect_false(any(calibration$history_can_deny_operation))
  expect_false(any(calibration$accuracy_depends_on_request_history))
})

test_that("random calibrations never present a preview as deployed selection", {
  withr::local_seed(20260807)
  valid <- logical(75L)
  for (iteration in seq_along(valid)) {
    coordinate_count <- sample(c(1L, 2L, 4L, 16L, 64L, 1000L), 1L)
    sensitivity <- sample.int(1000000L, 1L)
    capsule_delta <- sample(c(0, 10^runif(1L, -10, -1)), 1L)
    result <- ds.vertDPCalibrate(
      capsule_epsilon = 2^runif(1L, -45, 45),
      peer_count = sample(2:5, 1L),
      sensitivity = sensitivity,
      capsule_delta = capsule_delta,
      coordinate_count = coordinate_count,
      gaussian_l2_sensitivity = runif(1L, 1, sensitivity))
    utility_expected <- if (result$gaussian_sampler_supported &&
                    (!result$sampler_supported ||
                     result$gaussian_selection_radius <
                       result$laplace_selection_radius)) {
      "gaussian"
    } else if (result$sampler_supported) {
      "laplace"
    } else "none"
    valid[[iteration]] <-
      identical(result$utility_candidate, utility_expected) &&
      identical(result$preview_candidate, utility_expected) &&
      is.na(result$selected_candidate) &&
      is.na(result$selected_delta) &&
      identical(result$gaussian_backend_deployed, TRUE) &&
      identical(result$gaussian_preview_supported,
                result$gaussian_sampler_supported) &&
      identical(result$selection_authority,
                "signed_server_capsule_manifest") &&
      !is.nan(result$laplace_selection_radius) &&
      !is.nan(result$gaussian_selection_radius) &&
      identical(
        result$selector, "minimum_conservative_95_radius_v3") &&
      identical(
        result$selector_objective,
        if (coordinate_count == 1L) {
          "marginal_95_abs"
        } else "simultaneous_95_abs")
  }
  expect_true(all(valid))
})

test_that("postprocessed count boxes cover all bounded integer noise draws", {
  withr::local_seed(20260808)
  covered <- logical(2000L)
  for (iteration in seq_along(covered)) {
    truth <- sample(0:1000000, sample.int(32L, 1L), replace = TRUE)
    radius <- sample(0:10000, 1L)
    noise <- if (radius) {
      sample.int(2L * radius + 1L, length(truth), replace = TRUE) -
        radius - 1L
    } else rep(0, length(truth))
    released <- pmax(0, truth + noise)
    lower <- pmax(0, released - radius)
    upper <- released + radius
    covered[[iteration]] <- all(truth >= lower & truth <= upper)
  }
  expect_true(all(covered))
})
