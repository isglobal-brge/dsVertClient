.mh_cells <- function() {
  matrix(
    c(
      12, 8, 5, 15,
      20, 10, 10, 30,
      8, 12, 4, 16),
    nrow = 3L, byrow = TRUE,
    dimnames = list(
      c("young", "middle", "old"),
      c(
        "exposed_event", "exposed_nonevent",
        "unexposed_event", "unexposed_nonevent")))
}

.mh_array <- function(cells = .mh_cells()) {
  out <- array(
    0, dim = c(2L, 2L, nrow(cells)),
    dimnames = list(
      exposure = c("unexposed", "exposed"),
      outcome = c("no_event", "event"),
      stratum = rownames(cells)))
  out["exposed", "event", ] <- cells[, "exposed_event"]
  out["exposed", "no_event", ] <- cells[, "exposed_nonevent"]
  out["unexposed", "event", ] <- cells[, "unexposed_event"]
  out["unexposed", "no_event", ] <- cells[, "unexposed_nonevent"]
  out
}

.mh_dp_fixture <- function(table, epsilon = 3) {
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
    col_levels = unname(colnames(table)), coordinate_maximum = 100,
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

test_that("exact Mantel-Haenszel output matches stats for authorised arrays", {
  input <- .mh_array()
  fit <- ds.vertMantelHaenszel(
    input, exposed = "exposed", event = "event", correct = TRUE)
  reference <- stats::mantelhaen.test(
    fit$oriented_tables, correct = TRUE, exact = FALSE,
    conf.level = 0.95)

  expect_s3_class(fit, "ds.vertMantelHaenszel")
  expect_identical(fit$status, "ok")
  expect_identical(fit$common_odds_ratio$estimate_type, "finite")
  expect_equal(
    fit$common_odds_ratio$estimate, unname(reference$estimate),
    tolerance = 1e-15)
  expect_equal(
    fit$mantel_haenszel_test$statistic, unname(reference$statistic),
    tolerance = 1e-15)
  expect_equal(
    fit$mantel_haenszel_test$p_value, unname(reference$p.value),
    tolerance = 1e-15)
  expect_equal(
    as.numeric(fit$mantel_haenszel_test$confidence_interval),
    as.numeric(reference$conf.int), tolerance = 1e-15)
  expect_identical(fit$additional_server_calls, 0L)
  expect_match(fit$inferential_scope, "already-authorized", fixed = TRUE)
  expect_output(print(fit), "classical p-value")

  uncorrected <- ds.vertMantelHaenszel(
    input, exposed = "exposed", event = "event", correct = FALSE)
  reference_uncorrected <- stats::mantelhaen.test(
    uncorrected$oriented_tables, correct = FALSE, exact = FALSE)
  expect_equal(
    uncorrected$mantel_haenszel_test$statistic,
    unname(reference_uncorrected$statistic), tolerance = 1e-15)
})

test_that("matrix cell mapping is explicit, order-invariant, and guarded", {
  canonical <- .mh_cells()
  shuffled <- canonical[, c(3L, 1L, 4L, 2L), drop = FALSE]
  colnames(shuffled) <- c("ue", "ee", "un", "en")
  map <- c(
    exposed_event = "ee", exposed_nonevent = "en",
    unexposed_event = "ue", unexposed_nonevent = "un")

  expected <- ds.vertMantelHaenszel(canonical)
  actual <- ds.vertMantelHaenszel(shuffled, cell_map = map)
  indexed <- ds.vertMantelHaenszel(
    unname(shuffled),
    cell_map = c(
      exposed_event = 2, exposed_nonevent = 4,
      unexposed_event = 1, unexposed_nonevent = 3))
  expect_equal(actual$common_odds_ratio, expected$common_odds_ratio)
  expect_equal(indexed$common_odds_ratio$estimate,
               expected$common_odds_ratio$estimate)
  expect_identical(actual$orientation$strata, rownames(canonical))

  checked <- list(
    observed = canonical,
    disclosure_guard = list(passed = TRUE, threshold = 5L))
  checked_fit <- ds.vertMantelHaenszel(checked)
  expect_identical(
    checked_fit$input_provenance,
    "disclosure_checked_authorized_aggregate")
  expect_identical(checked_fit$disclosure_guard, checked$disclosure_guard)

  checked$disclosure_guard$passed <- FALSE
  expect_error(ds.vertMantelHaenszel(checked), "disclosure guard")
  expect_error(ds.vertMantelHaenszel(unname(canonical)), "cell_map")
  expect_error(
    ds.vertMantelHaenszel(
      shuffled, cell_map = c(
        exposed_event = "ee", exposed_nonevent = "ee",
        unexposed_event = "ue", unexposed_nonevent = "un")),
    "distinct")
})

test_that("exact Mantel-Haenszel validation and boundary types are explicit", {
  expect_error(
    ds.vertMantelHaenszel(matrix(1:9, nrow = 3L)),
    "exactly four columns")
  invalid <- .mh_cells()
  invalid[1L, 1L] <- -1
  expect_error(ds.vertMantelHaenszel(invalid), "non-negative")
  invalid <- .mh_cells()
  invalid[1L, 1L] <- 1.5
  expect_error(ds.vertMantelHaenszel(invalid), "whole-number")
  duplicate_strata <- .mh_cells()
  rownames(duplicate_strata) <- c("same", "same", "other")
  expect_error(ds.vertMantelHaenszel(duplicate_strata), "unique")
  expect_error(
    ds.vertMantelHaenszel(.mh_array(), exposed = "missing",
                          event = "event"),
    "Unknown exposed")
  dp_chisq <- structure(
    list(observed = matrix(1:4, 2L), calibration = "dp-bootstrap-v1"),
    class = c("ds.vertChisq", "list"))
  expect_error(ds.vertMantelHaenszel(dp_chisq),
               "ds.vertDPMantelHaenszel", fixed = TRUE)

  zero <- matrix(
    c(0, 10, 5, 0), nrow = 1L,
    dimnames = list("s", .dsvert_mh_cell_roles()))
  infinite <- matrix(
    c(5, 0, 0, 10), nrow = 1L,
    dimnames = list("s", .dsvert_mh_cell_roles()))
  undefined <- matrix(
    c(5, 0, 0, 0), nrow = 1L,
    dimnames = list("s", .dsvert_mh_cell_roles()))
  expect_identical(
    ds.vertMantelHaenszel(zero)$common_odds_ratio$estimate_type, "zero")
  expect_identical(
    ds.vertMantelHaenszel(infinite)$common_odds_ratio$estimate_type,
    "infinite")
  undefined_fit <- ds.vertMantelHaenszel(undefined)
  expect_identical(
    undefined_fit$common_odds_ratio$estimate_type, "non_estimable")
  expect_true(is.na(undefined_fit$common_odds_ratio$estimate))
  expect_identical(
    undefined_fit$mantel_haenszel_test$status,
    "not_available_single_stratum")
  unnamed_single <- ds.vertMantelHaenszel(matrix(c(8, 2, 4, 6), 2L))
  expect_identical(
    unnamed_single$mantel_haenszel_test$status,
    "not_available_single_stratum")
  expect_identical(unnamed_single$orientation$strata, "stratum_1")
})

test_that("DP Mantel-Haenszel is zero-cost post-processing of one table", {
  table <- .mh_cells()
  release <- .mh_dp_fixture(table)
  fit <- ds.vertDPMantelHaenszel(release)
  exact <- .dsvert_mh_point(table)

  expect_s3_class(fit, "ds.vertDPMantelHaenszel")
  expect_identical(fit$status, "ok")
  expect_equal(fit$estimate, exact$estimate, tolerance = 1e-15)
  expect_identical(fit$estimate_type, "finite")
  expect_lte(fit$mechanism_region[["lower"]], fit$estimate)
  expect_gte(fit$mechanism_region[["upper"]], fit$estimate)
  expect_identical(fit$additional_server_calls, 0L)
  expect_identical(
    fit$additional_privacy_cost, c(epsilon = 0, delta = 0))
  expect_match(fit$inferential_scope, "no classical CMH p-value",
               fixed = TRUE)
  expect_false(any(c("statistic", "p_value") %in% names(fit)))
  expect_match(fit$interval_construction, "interval division")
  expect_match(fit$contribution_contract, "exactly one global")
  expect_output(print(fit), "common odds ratio")

  shuffled <- table[, c(3L, 1L, 4L, 2L), drop = FALSE]
  colnames(shuffled) <- c("ue", "ee", "un", "en")
  mapped <- ds.vertDPMantelHaenszel(
    .mh_dp_fixture(shuffled),
    cell_map = c(
      exposed_event = "ee", exposed_nonevent = "en",
      unexposed_event = "ue", unexposed_nonevent = "un"))
  expect_equal(mapped$estimate, fit$estimate, tolerance = 0)
  expect_identical(mapped$orientation$strata, rownames(table))

  replace_release <- release
  replace_release$adjacency <- "replace_one_fixed_cohort"
  replace_release$artifact_l1_sensitivity <- 2
  replace_release$artifact_l2_sensitivity <- sqrt(2)
  expect_error(
    ds.vertDPMantelHaenszel(replace_release),
    "block L1 sensitivity 1", fixed = TRUE)
})

test_that("DP Mantel-Haenszel point boundaries are typed without tests", {
  make <- function(values) {
    table <- matrix(
      values, nrow = 1L,
      dimnames = list("s", .dsvert_mh_cell_roles()))
    ds.vertDPMantelHaenszel(.mh_dp_fixture(table))
  }
  finite <- make(c(5, 2, 3, 7))
  zero <- make(c(0, 2, 3, 7))
  infinite <- make(c(5, 0, 0, 7))
  undefined <- make(c(5, 0, 0, 0))

  expect_identical(finite$estimate_type, "finite")
  expect_identical(zero$estimate_type, "zero")
  expect_equal(zero$estimate, 0)
  expect_identical(infinite$estimate_type, "infinite")
  expect_true(is.infinite(infinite$estimate))
  expect_identical(undefined$estimate_type, "non_estimable")
  expect_true(is.na(undefined$estimate))
  for (fit in list(finite, zero, infinite, undefined)) {
    expect_identical(fit$additional_server_calls, 0L)
    expect_false(any(c("statistic", "p_value") %in% names(fit)))
  }
})

.mh_validate_all_small_boxes <- function(strata, maximum) {
  roles <- .dsvert_mh_cell_roles()
  coordinates <- 4L * strata
  point_grid <- as.matrix(expand.grid(
    rep(list(0:maximum), coordinates), KEEP.OUT.ATTRS = FALSE))
  point_values <- vapply(seq_len(nrow(point_grid)), function(index) {
    cells <- matrix(
      point_grid[index, ], nrow = strata, byrow = TRUE,
      dimnames = list(paste0("s", seq_len(strata)), roles))
    .dsvert_mh_point(cells)$estimate
  }, numeric(1L))
  pairs <- as.matrix(expand.grid(
    lower = 0:maximum, upper = 0:maximum,
    KEEP.OUT.ATTRS = FALSE))
  pairs <- pairs[pairs[, "lower"] <= pairs[, "upper"], , drop = FALSE]
  box_grid <- as.matrix(expand.grid(
    rep(list(seq_len(nrow(pairs))), coordinates),
    KEEP.OUT.ATTRS = FALSE))

  for (box_index in seq_len(nrow(box_grid))) {
    selected <- box_grid[box_index, ]
    lower_vector <- pairs[selected, "lower"]
    upper_vector <- pairs[selected, "upper"]
    lower <- matrix(
      lower_vector, nrow = strata, byrow = TRUE,
      dimnames = list(paste0("s", seq_len(strata)), roles))
    upper <- matrix(
      upper_vector, nrow = strata, byrow = TRUE,
      dimnames = list(paste0("s", seq_len(strata)), roles))
    region <- .dsvert_dp_mh_region(lower, upper)
    inside <- rowSums(sweep(point_grid, 2L, lower_vector, `>=`) &
                        sweep(point_grid, 2L, upper_vector, `<=`)) ==
      coordinates
    values <- point_values[inside]
    finite <- values[is.finite(values)]
    if (length(finite) &&
        (any(finite < region$interval[["lower"]] - 1e-14) ||
         any(finite > region$interval[["upper"]] + 1e-14))) {
      stop("finite point escaped MH interval in box ", box_index,
           call. = FALSE)
    }
    if (any(values == 0, na.rm = TRUE) && !region$includes_zero) {
      stop("zero boundary was not flagged in box ", box_index,
           call. = FALSE)
    }
    if (any(is.infinite(values)) &&
        (!region$includes_infinite ||
         !is.infinite(region$interval[["upper"]]))) {
      stop("infinite boundary was not covered in box ", box_index,
           call. = FALSE)
    }
    if (anyNA(values) && !region$includes_non_estimable) {
      stop("non-estimable boundary was not flagged in box ", box_index,
           call. = FALSE)
    }
    if (any(!is.na(values)) && !region$has_estimable_values) {
      stop("estimable point was omitted in box ", box_index,
           call. = FALSE)
    }
  }
  TRUE
}

test_that("DP MH interval arithmetic exhaustively encloses small boxes", {
  expect_true(.mh_validate_all_small_boxes(strata = 1L, maximum = 2L))
  expect_true(.mh_validate_all_small_boxes(strata = 2L, maximum = 1L))
})
