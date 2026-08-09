.dp_prevalence_table_fixture <- function(table, epsilon = 3) {
  table <- as.matrix(table)
  result <- list(
    released = TRUE,
    mechanism =
      "two-independent-complete-vector-discrete-laplace-draws-v3",
    implementation = .DSVERT_CLIENT_VECTOR_BACKEND,
    backend = "exact_signed_Ring128_global_vector",
    sampler = .DSVERT_CLIENT_VECTOR_SAMPLER,
    randomness = paste(
      "independent pinned-peer HKDF-SHA256/ChaCha20 streams;",
      "no analyst-controlled seed"),
    epsilon = epsilon, delta = 2^-100,
    implementation_delta =
      "1/1267650600228229401496703205376",
    adjacency = "add_remove_patient",
    sensitivity = 2, sensitivity_norm = "l1",
    l1_sensitivity = 2, l2_sensitivity = sqrt(2),
    sensitivity_scope = "complete_signed_biomedical_capsule_vector",
    output_lattice_bits = 8L, output_lattice_scale = 256,
    sticky_noise = "immutable_capsule_durable_replay_v3",
    sticky_replay = TRUE, privacy_epochs = c(1, 1),
    noise_key_ids = c("noise-a", "noise-b"),
    source_values_exposed = FALSE, intermediate_values_exposed = FALSE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    clipped_coordinates = NA_integer_, clamp_activation_disclosed = FALSE,
    table = table, counts = unname(as.numeric(table)),
    nrow = as.integer(nrow(table)), ncol = as.integer(ncol(table)),
    row_levels = unname(rownames(table)),
    col_levels = unname(colnames(table)), coordinate_maximum = 1000,
    artifact_l1_sensitivity = 1, artifact_l2_sensitivity = 1,
    unit_aggregation_policy = "consistent_joint_cell_else_exclude_v1",
    server = "site-a",
    accuracy_simultaneous_confidence = 0.95,
    accuracy_simultaneous_method = paste(
      "exact ideal two-sided-geometric convolution tail with union bound;",
      "two-peer finite-sampler TV deducted; fixed-clamp range applied"),
    accuracy_additional_privacy_cost = c(epsilon = 0, delta = 0))
  result$accuracy_simultaneous_95_abs <-
    .dsvert_dp_vector_table_radius(result, 0.95)
  class(result) <- c("ds.vertDPContingency", "list")
  result
}

.dp_prevalence_fixture <- function(table = matrix(
    c(40, 20, 10, 30), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("prevalent", "not_prevalent")))) {
  .dp_prevalence_table_fixture(table)
}

test_that("prevalence view is a numeric identity of DP 2x2 epidemiology", {
  release <- .dp_prevalence_fixture()
  base <- ds.vertDPEpi2x2(
    release, exposed = "exposed", event = "prevalent", level = 0.95)
  view <- ds.vertDPPrevalenceRatio(
    release, exposed = "exposed", prevalent = "prevalent", level = 0.95)

  expect_s3_class(view, "ds.vertDPPrevalenceRatio")
  expect_s3_class(view, "ds.vertDPEpi2x2")
  expect_identical(view$point_estimates, base$point_estimates)
  expect_identical(view$point_status, base$point_status)
  expect_identical(view$mechanism_regions, base$mechanism_regions)
  expect_identical(view$mechanism_region_types,
                   base$mechanism_region_types)
  expect_identical(view$number_needed, base$number_needed)
  expect_identical(view$simultaneous_radius, base$simultaneous_radius)
  expect_identical(view$coverage_method, base$coverage_method)
  expect_identical(view$epsilon, base$epsilon)
  expect_identical(view$delta, base$delta)
  expect_identical(view$server, base$server)
  expect_identical(view$additional_server_calls, 0L)
  expect_identical(view$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))

  map <- view$prevalence_measure_mapping
  expect_identical(
    view$prevalence_point_estimates,
    stats::setNames(base$point_estimates[unname(map)], names(map)))
  expect_identical(
    view$prevalence_mechanism_regions,
    stats::setNames(base$mechanism_regions[unname(map)], names(map)))
  expect_identical(view$number_needed_from_prevalence_difference,
                   base$number_needed)
})

test_that("prevalence inference view preserves every coverage object", {
  release <- .dp_prevalence_fixture()
  base <- ds.vertDPEpi2x2Inference(
    release, exposed = "exposed", event = "prevalent", level = 0.9,
    mechanism_alpha_share = 0.3)
  view <- ds.vertDPPrevalenceRatioInference(
    release, exposed = "exposed", prevalent = "prevalent", level = 0.9,
    mechanism_alpha_share = 0.3)

  expect_s3_class(view, "ds.vertDPPrevalenceRatioInference")
  expect_s3_class(view, "ds.vertDPEpi2x2Inference")
  for (field in c(
      "point_estimates", "point_status", "combined_regions",
      "combined_region_types", "mechanism_regions",
      "confidential_count_integer_box", "number_needed", "level",
      "coverage_lower_bound", "mechanism_level",
      "sampling_familywise_level", "base_sampling_interval_level",
      "alpha_allocation", "coverage_method", "sampling_model",
      "uncertainty_scope", "inferential_scope", "epsilon", "delta",
      "server")) {
    expect_identical(view[[field]], base[[field]], info = field)
  }
  expect_identical(view$additional_server_calls, 0L)
  expect_identical(view$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))

  map <- view$prevalence_measure_mapping
  expect_identical(
    view$prevalence_combined_regions,
    stats::setNames(base$combined_regions[unname(map)], names(map)))
  expect_identical(
    view$prevalence_mechanism_regions,
    stats::setNames(base$mechanism_regions[unname(map)], names(map)))
  expect_identical(view$number_needed_from_prevalence_difference,
                   base$number_needed)
})

test_that("cross-sectional design and both orientations are explicit", {
  release <- .dp_prevalence_fixture()
  expect_error(ds.vertDPPrevalenceRatio(release),
               "exposed and prevalent must be supplied explicitly")
  expect_error(
    ds.vertDPPrevalenceRatio(release, exposed = "exposed"),
    "exposed and prevalent must be supplied explicitly")
  expect_error(ds.vertDPPrevalenceRatioInference(release),
               "exposed and prevalent must be supplied explicitly")
  expect_error(
    ds.vertDPPrevalenceRatioInference(
      release, exposed = "exposed", prevalent = "missing"),
    "Unknown prevalent level")

  view <- ds.vertDPPrevalenceRatio(
    release, exposed = "exposed", prevalent = "prevalent")
  expect_identical(view$study_design, "cross_sectional")
  expect_identical(view$study_design_inferred_from_table, FALSE)
  expect_match(view$study_design_source, "caller_declared")
  expect_identical(view$orientation$exposed, "exposed")
  expect_identical(view$orientation$prevalent, "prevalent")
})

test_that("prevalence wrappers make no DSI or capsule retrieval call", {
  release <- .dp_prevalence_fixture()
  reached_remote <- function(...) {
    stop("test observed a remote operation", call. = FALSE)
  }
  testthat::local_mocked_bindings(
    .dsvert_dp_datasources = reached_remote,
    .dsvert_dp_capsule_vector_run = reached_remote,
    .dsvert_aggregate_strict = reached_remote,
    .dsvert_fanout_by_site = reached_remote,
    .package = "dsVertClient")

  point <- ds.vertDPPrevalenceRatio(
    release, exposed = "exposed", prevalent = "prevalent")
  inference <- ds.vertDPPrevalenceRatioInference(
    release, exposed = "exposed", prevalent = "prevalent")
  expect_identical(point$additional_server_calls, 0L)
  expect_identical(inference$additional_server_calls, 0L)
})

test_that("tamper and boundary behavior are inherited without narrowing", {
  release <- .dp_prevalence_fixture()
  tampered <- release
  tampered$request_limit <- TRUE
  expect_error(
    ds.vertDPPrevalenceRatio(
      tampered, exposed = "exposed", prevalent = "prevalent"),
    "validated")
  expect_error(
    ds.vertDPPrevalenceRatioInference(
      tampered, exposed = "exposed", prevalent = "prevalent"),
    "validated")

  boundary_table <- matrix(
    c(0, 20, 10, 10), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("prevalent", "not_prevalent")))
  boundary <- .dp_prevalence_fixture(boundary_table)
  base <- ds.vertDPEpi2x2(
    boundary, exposed = "exposed", event = "prevalent")
  view <- ds.vertDPPrevalenceRatio(
    boundary, exposed = "exposed", prevalent = "prevalent")
  expect_identical(view$point_status$risk_ratio,
                   base$point_status$risk_ratio)
  expect_identical(view$prevalence_mechanism_regions$prevalence_ratio,
                   base$mechanism_regions$risk_ratio)
  expect_identical(view$number_needed_from_prevalence_difference,
                   base$number_needed)

  empty_table <- matrix(
    c(0, 0, 10, 10), nrow = 2L, byrow = TRUE,
    dimnames = dimnames(boundary_table))
  empty <- .dp_prevalence_fixture(empty_table)
  empty_base <- ds.vertDPEpi2x2(
    empty, exposed = "exposed", event = "prevalent")
  empty_view <- ds.vertDPPrevalenceRatio(
    empty, exposed = "exposed", prevalent = "prevalent")
  expect_identical(empty_view$point_estimates, empty_base$point_estimates)
  expect_identical(empty_view$mechanism_regions,
                   empty_base$mechanism_regions)
})

test_that("prevalence print methods disclose the semantic non-claim", {
  release <- .dp_prevalence_fixture()
  point <- ds.vertDPPrevalenceRatio(
    release, exposed = "exposed", prevalent = "prevalent")
  inference <- ds.vertDPPrevalenceRatioInference(
    release, exposed = "exposed", prevalent = "prevalent")

  point_output <- capture.output(print(point))
  inference_output <- capture.output(print(inference))
  expect_true(any(grepl("cross-sectional", point_output, fixed = TRUE)))
  expect_true(any(grepl("not inferred from the table", point_output,
                        fixed = TRUE)))
  expect_true(any(grepl("prevalence ratio", point_output,
                        ignore.case = TRUE)))
  expect_true(any(grepl("cross-sectional", inference_output,
                        fixed = TRUE)))
  expect_true(any(grepl("not inferred from the table", inference_output,
                        fixed = TRUE)))
})
