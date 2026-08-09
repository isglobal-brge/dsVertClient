.dp_table_fixture <- function(table, epsilon = 1) {
  table <- as.matrix(table)
  simultaneous_radius <- .dsvert_dp_google_ci_radius(
    epsilon, 1, confidence = 1 - 0.05 / length(table))
  marginal_radius <- .dsvert_dp_google_ci_radius(epsilon, 1)
  laplace <- list(
    available = TRUE,
    mechanism = "dsvert_dp_v1_deterministic_granular_laplace_int64",
    epsilon = epsilon, delta = 0, analytic_delta = 0,
    implementation_delta_bound = 0,
    accounting_rule = "pure_dp_no_implementation_slack",
    accuracy_accounting = "exact_granular_laplace_confidence_interval",
    sensitivity_norm = "l1", sensitivity = 1,
    marginal_95_abs = marginal_radius,
    simultaneous_95_abs = simultaneous_radius,
    nominal_rmse = sqrt(2) / epsilon, sigma = 0,
    granularity = 2^ceiling(log2((1 / epsilon) / 2^40)),
    analytic_accounting_verified = TRUE, unavailable_reason = "")
  gaussian <- list(
    available = FALSE,
    mechanism = "dsvert_dp_v3_deterministic_approximate_gaussian_int64",
    epsilon = epsilon, delta = 0, analytic_delta = 0,
    implementation_delta_bound =
      .dsvert_dp_gaussian_implementation_delta_bound(
        length(table), epsilon),
    accounting_rule =
      "analytic_gaussian_delta_plus_dp_transfer_from_total_variation_bound",
    accuracy_accounting =
      "gaussian_tail_alpha_minus_total_variation_union_bound",
    sensitivity_norm = "l2", sensitivity = 1,
    marginal_95_abs = 0, simultaneous_95_abs = 0,
    nominal_rmse = 0, sigma = 0, granularity = 0,
    analytic_accounting_verified = FALSE,
    unavailable_reason = "gaussian_delta_is_zero")
  noise_selection <- list(
    schema_version = 2L,
    selector = "minimum_conservative_95_radius_v3",
    objective = "simultaneous_95_abs", coordinate_count = length(table),
    laplace = laplace, gaussian = gaussian, winner = "laplace",
    winner_mechanism = laplace$mechanism,
    winning_metric_abs = simultaneous_radius, winner_delta = 0,
    tie_break = "laplace_unless_gaussian_strictly_improves")
  structure(list(
    released = TRUE, table = table, counts = as.numeric(table),
    nrow = nrow(table), ncol = ncol(table),
    row_levels = rownames(table), col_levels = colnames(table),
    mechanism = "dsvert_dp_v1_deterministic_granular_laplace_int64",
    implementation = paste0(
      "dsVert adapted Google Differential Privacy v4.1.0 ",
      "granular Laplace integer mechanism"),
    sampler = "deterministic_two_sided_geometric",
    randomness = "HMAC-SHA256/ChaCha20",
    privacy_epoch = 1, noise_key_id = "test-noise-key-v1",
    sticky_noise = "dsvert-sticky-noise-v1",
    sensitivity = 1, l1_sensitivity = 1, l2_sensitivity = 1,
    noise_selection = noise_selection, clipped_coordinates = 0L,
    accuracy_simultaneous_95_abs = simultaneous_radius,
    accuracy_simultaneous_confidence = 0.95,
    accuracy_simultaneous_method = "union_bound",
    unit_aggregation_policy = "consistent_cell_else_exclude_v1",
    epsilon = epsilon, delta = 0, server = "site-a"),
    class = c("ds.vertDPContingency", "list"))
}

.dp_vector_table_fixture <- function(table, epsilon = 3) {
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

.dp_vector_table_gaussian_fixture <- function(table, epsilon = 3) {
  result <- .dp_vector_table_fixture(table, epsilon)
  coordinate_count <- length(result$table)
  sensitivity_steps <- format(
    result$l2_sensitivity * result$output_lattice_scale,
    digits = 17L, scientific = TRUE, trim = TRUE)
  plan <- list(
    version = .DSVERT_CLIENT_VECTOR_GAUSSIAN_PLAN_VERSION,
    mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
    sampler = .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER,
    request_binding_sha256 = digest::digest(
      "epidemiology-gaussian-request", "sha256", serialize = FALSE),
    total_coordinate_count = as.numeric(coordinate_count),
    maximum_chunk_coordinates = as.numeric(coordinate_count),
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
    sigma_squared_numerator = "25",
    sigma_squared_denominator = "2",
    maximum_noise_magnitude_two_peers = "64",
    observable_worker_shape = "fixed dyadic CDF fixture",
    per_peer_implementation_delta_numerator = "1",
    per_peer_implementation_delta_denominator =
      "1267650600228229401496703205376",
    simultaneous_95_abs = "17")
  plan <- .client_complete_gaussian_plan_v2(plan)
  request <- list(
    epsilon = format(epsilon, digits = 17L, scientific = TRUE, trim = TRUE),
    delta = format(result$delta, digits = 17L,
                   scientific = TRUE, trim = TRUE),
    l2_sensitivity_steps = sensitivity_steps,
    total_coordinate_count = as.numeric(coordinate_count))
  result$capsule_mechanism <- list(
    mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
    sensitivity_norm = "l2")
  result$mechanism_selection <- list(
    gaussian_calibration_request = request,
    gaussian_plan = plan,
    gaussian_plan_sha256 = .dsvert_vector_hash(plan))
  result$mechanism <- .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM
  result$implementation <- .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND
  result$sampler <- .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER
  result$sensitivity <- result$l2_sensitivity
  result$sensitivity_norm <- "l2"
  result$capsule_coordinate_count <- as.integer(coordinate_count)
  accuracy <- .dsvert_dp_vector_table_radius(result, 0.95)
  result$accuracy_simultaneous_95_abs <- accuracy
  result$accuracy_simultaneous_method <- paste(
    "signed fixed-work dyadic discrete-Gaussian plan v2 simultaneous 95% bound;",
    "tail and CDF TV transfers already charged; fixed-clamp range applied")
  result
}

test_that("epidemiology accepts the signed global-vector contingency", {
  table <- matrix(
    c(40, 20, 10, 30), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("event", "nonevent")))
  released <- .dp_vector_table_fixture(table)
  expect_null(released$noise_selection)
  expect_identical(.dsvert_dp_table_contract(released), released)

  epi <- ds.vertDPEpi2x2(
    released, exposed = "exposed", event = "event")
  diagnostic <- ds.vertDPDiagnostic2x2(
    released, disease_positive = "exposed", test_positive = "event")
  direct <- ds.vertDPDirectStandardization(
    released, c(unexposed = 0.4, exposed = 0.6), event = "event")
  indirect <- ds.vertDPIndirectStandardization(
    released, c(unexposed = 0.2, exposed = 0.3), event = "event")
  for (value in list(epi, diagnostic, direct, indirect)) {
    expect_match(value$coverage_method, "convolution tail", fixed = TRUE)
    expect_true(is.finite(value$simultaneous_radius))
  }
  expect_gte(.dsvert_dp_table_simultaneous_radius(released, 0.99),
             .dsvert_dp_table_simultaneous_radius(released, 0.95))

  forged <- released
  forged$accuracy_simultaneous_95_abs <-
    forged$accuracy_simultaneous_95_abs - 1
  expect_error(.dsvert_dp_table_contract(forged), "validated")
})

test_that("epidemiology consumes the signed Gaussian plan without rerouting", {
  table <- matrix(
    c(40, 20, 10, 30), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("event", "nonevent")))
  released <- .dp_vector_table_gaussian_fixture(table)
  expect_identical(.dsvert_dp_table_contract(released), released)
  expect_equal(
    .dsvert_dp_vector_table_radius(released, 0.95), 17 / 256,
    tolerance = 0)

  results <- list(
    ds.vertDPEpi2x2(released, exposed = "exposed", event = "event"),
    ds.vertDPDiagnostic2x2(
      released, disease_positive = "exposed", test_positive = "event"),
    ds.vertDPDirectStandardization(
      released, c(unexposed = 0.4, exposed = 0.6), event = "event"),
    ds.vertDPIndirectStandardization(
      released, c(unexposed = 0.2, exposed = 0.3), event = "event"))
  for (value in results) {
    expect_match(value$coverage_method, "discrete-Gaussian plan", fixed = TRUE)
    expect_equal(value$simultaneous_radius, 17 / 256, tolerance = 0)
  }
  expect_gte(
    .dsvert_dp_table_simultaneous_radius(released, 0.99),
    .dsvert_dp_table_simultaneous_radius(released, 0.95))
  expect_equal(
    .dsvert_dp_table_simultaneous_radius(released, 0.99), 19 / 256,
    tolerance = 1e-15)
  confidence_grid <- c(0.5, 0.9, 0.95, 0.975, 0.99, 0.999)
  radii <- vapply(
    confidence_grid,
    function(level) .dsvert_dp_table_simultaneous_radius(released, level),
    numeric(1L))
  expect_true(all(diff(radii) >= 0))

  support_plan <- released$mechanism_selection$gaussian_plan
  support_plan$vector_total_tv_upper_numerator <- "1"
  support_plan$vector_total_tv_upper_denominator <- "100"
  support_bound <- .dsvert_dp_vector_gaussian_accuracy_steps(
    support_plan, coordinate_count = 4L, confidence = 0.99)
  expect_identical(support_bound$finite_support, TRUE)
  expect_gte(support_bound$steps, 64)

  incomplete_plan <- released$mechanism_selection$gaussian_plan
  incomplete_plan$sigma_squared_numerator <- NULL
  expect_error(
    .dsvert_dp_vector_gaussian_accuracy_steps(
      incomplete_plan, coordinate_count = 4L, confidence = 0.99),
    "accuracy certificate is invalid")

  forged <- released
  forged$mechanism_selection$gaussian_plan$simultaneous_95_abs <- "18"
  expect_error(.dsvert_dp_table_contract(forged), "validated")
})

test_that("DP 2x2 regions contain effects at the released table", {
  tab <- matrix(c(40, 20, 10, 30), nrow = 2L, byrow = TRUE,
                dimnames = list(c("no", "yes"), c("event", "nonevent")))
  result <- ds.vertDPEpi2x2(
    .dp_table_fixture(tab, epsilon = 2),
    exposed = "yes", event = "event")
  expect_s3_class(result, "ds.vertDPEpi2x2")
  expect_identical(result$status, "ok")
  expect_gte(result$simultaneous_radius, 0)
  expect_lte(result$mechanism_regions$risk_difference[["lower"]],
             result$point_estimates[["risk_difference"]])
  expect_gte(result$mechanism_regions$risk_difference[["upper"]],
             result$point_estimates[["risk_difference"]])
  expect_lte(result$mechanism_regions$risk_ratio[["lower"]],
             result$point_estimates[["risk_ratio"]])
  expect_gte(result$mechanism_regions$risk_ratio[["upper"]],
             result$point_estimates[["risk_ratio"]])
  expect_match(result$uncertainty_scope, "sampling uncertainty excluded")
})

test_that("DP diagnostic ROC is zero-cost post-processing of one ordered table", {
  tab <- matrix(
    c(40, 30, 20, 10,
       5, 10, 20, 40),
    nrow = 2L, byrow = TRUE,
    dimnames = list(c("control", "case"), paste0("q", 1:4)))
  released <- .dp_vector_table_fixture(tab, epsilon = 3)

  roc <- ds.vertDPROC(
    released, disease_positive = "case", direction = "higher")
  expect_s3_class(roc, "ds.vertDPROC")
  expect_identical(roc$status, "ok")
  expect_equal(roc$curve$false_positive_rate[c(1L, nrow(roc$curve))],
               c(0, 1), tolerance = 0)
  expect_equal(roc$curve$sensitivity[c(1L, nrow(roc$curve))],
               c(0, 1), tolerance = 0)
  expect_true(all(diff(roc$curve$false_positive_rate) >= 0))
  expect_true(all(diff(roc$curve$sensitivity) >= 0))

  cases <- tab["case", ]
  controls <- tab["control", ]
  expected_auc <- sum(vapply(seq_along(cases), function(index) {
    cases[[index]] *
      (sum(controls[seq_len(index - 1L)]) + 0.5 * controls[[index]])
  }, numeric(1L))) / (sum(cases) * sum(controls))
  expect_equal(roc$auc, expected_auc, tolerance = 1e-15)
  expect_lte(roc$auc_mechanism_region[["lower"]], roc$auc)
  expect_gte(roc$auc_mechanism_region[["upper"]], roc$auc)
  expect_identical(roc$additional_server_calls, 0L)
  expect_identical(roc$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_match(roc$uncertainty_scope, "sampling uncertainty excluded",
               fixed = TRUE)

  reverse <- ds.vertDPROC(
    released, disease_positive = "case", direction = "lower")
  expect_equal(roc$auc + reverse$auc, 1, tolerance = 1e-15)
  expect_identical(reverse$score_order, rev(roc$score_order))
})

test_that("DP ROC mechanism box encloses every admissible table", {
  tab <- matrix(
    c(40, 30, 20, 10,
       5, 10, 20, 40),
    nrow = 2L, byrow = TRUE,
    dimnames = list(c("control", "case"), paste0("q", 1:4)))
  roc <- ds.vertDPROC(
    .dp_vector_table_fixture(tab, epsilon = 3),
    disease_positive = "case", direction = "higher")

  set.seed(90210)
  for (iteration in seq_len(250L)) {
    candidate <- matrix(
      runif(length(tab), as.numeric(roc$count_lower),
            as.numeric(roc$count_upper)),
      nrow = 2L, byrow = FALSE,
      dimnames = dimnames(tab))
    candidate_auc <- .dsvert_dp_auc_from_ordered_bins(
      candidate["case", ], candidate["control", ])
    if (is.finite(candidate_auc)) {
      expect_gte(candidate_auc,
                 roc$auc_mechanism_region[["lower"]] - 1e-12)
      expect_lte(candidate_auc,
                 roc$auc_mechanism_region[["upper"]] + 1e-12)
    }
  }
})

test_that("DP ROC validates order and handles non-estimable classes explicitly", {
  tab <- matrix(
    c(10, 8, 6, 4, 0, 0, 0, 0), nrow = 2L, byrow = TRUE,
    dimnames = list(c("control", "case"), paste0("q", 1:4)))
  released <- .dp_vector_table_fixture(tab, epsilon = 2)
  result <- ds.vertDPROC(released, disease_positive = "case")
  expect_identical(result$status, "non_estimable_zero_disease_total")
  expect_true(is.na(result$auc))
  expect_true(result$auc_mechanism_region_includes_non_estimable)
  expect_equal(result$auc_mechanism_region, c(lower = 0, upper = 1))

  expect_error(
    ds.vertDPROC(released, disease_positive = "case",
                 score_order = c("q1", "q1", "q3", "q4")),
    "permutation")
  expect_error(
    ds.vertDPROC(released, disease_positive = "case", direction = "sideways"),
    "arg")
  three_rows <- rbind(tab, unknown = c(1, 1, 1, 1))
  expect_error(
    ds.vertDPROC(.dp_table_fixture(three_rows), disease_positive = "case"),
    "exactly two disease-status rows")
})

test_that("DP 2x2 includes attributable fractions and typed number needed", {
  tab <- matrix(
    c(40, 20, 10, 30), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("event", "nonevent")))
  result <- ds.vertDPEpi2x2(
    .dp_vector_table_fixture(tab, epsilon = 3),
    exposed = "exposed", event = "event")
  risk_exposed <- 10 / 40
  risk_unexposed <- 40 / 60
  risk_difference <- risk_exposed - risk_unexposed
  population_risk <- 50 / 100

  expect_equal(
    result$point_estimates$attributable_fraction_exposed,
    risk_difference / risk_exposed, tolerance = 1e-15)
  expect_lte(
    result$mechanism_regions$attributable_fraction_exposed[["lower"]],
    result$point_estimates$attributable_fraction_exposed)
  expect_gte(
    result$mechanism_regions$attributable_fraction_exposed[["upper"]],
    result$point_estimates$attributable_fraction_exposed)
  expect_equal(
    result$point_estimates$population_attributable_fraction,
    (population_risk - risk_unexposed) / population_risk,
    tolerance = 1e-15)
  expect_lte(
    result$mechanism_regions$population_attributable_fraction[["lower"]],
    result$point_estimates$population_attributable_fraction)
  expect_gte(
    result$mechanism_regions$population_attributable_fraction[["upper"]],
    result$point_estimates$population_attributable_fraction)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_identical(result$additional_server_calls, 0L)
  expect_identical(result$number_needed$point_direction, "benefit")
  expect_equal(result$number_needed$point_estimate,
               1 / abs(risk_difference), tolerance = 1e-15)
  expect_lte(result$number_needed$mechanism_regions$benefit[["lower"]],
             result$number_needed$point_estimate)
  expect_gte(result$number_needed$mechanism_regions$benefit[["upper"]],
             result$number_needed$point_estimate)

  equal_risks <- matrix(
    c(20, 20, 10, 10), nrow = 2L, byrow = TRUE,
    dimnames = dimnames(tab))
  boundary <- ds.vertDPEpi2x2(
    .dp_vector_table_fixture(equal_risks, epsilon = 3),
    exposed = "exposed", event = "event")
  expect_identical(boundary$number_needed$point_direction, "none")
  expect_true(is.infinite(boundary$number_needed$point_estimate))
  expect_true(boundary$number_needed$mechanism_region_includes_infinite)
})

test_that("random 2x2 orientations agree with direct epidemiology and diagnostic oracles", {
  withr::local_seed(20260805)
  valid <- logical(500L)
  for (iteration in seq_along(valid)) {
    table <- matrix(
      sample.int(100L, 4L, replace = TRUE), nrow = 2L,
      dimnames = list(c("row_a", "row_b"), c("col_a", "col_b")))
    positive_row <- sample.int(2L, 1L)
    positive_col <- sample.int(2L, 1L)
    negative_row <- setdiff(1:2, positive_row)
    negative_col <- setdiff(1:2, positive_col)

    diagnostic <- ds.vertDPDiagnostic2x2(
      .dp_table_fixture(table, epsilon = 3),
      disease_positive = positive_row, test_positive = positive_col)
    cells <- c(
      true_positive = table[positive_row, positive_col],
      false_negative = table[positive_row, negative_col],
      false_positive = table[negative_row, positive_col],
      true_negative = table[negative_row, negative_col])
    diagnostic_oracle <- .dsvert_dp_diagnostic_values(cells)

    epi <- ds.vertDPEpi2x2(
      .dp_table_fixture(table, epsilon = 3),
      exposed = positive_row, event = positive_col)
    exposed_risk <- table[positive_row, positive_col] /
      sum(table[positive_row, ])
    unexposed_risk <- table[negative_row, positive_col] /
      sum(table[negative_row, ])
    epi_oracle <- c(
      risk_exposed = exposed_risk,
      risk_unexposed = unexposed_risk,
      risk_difference = exposed_risk - unexposed_risk,
      risk_ratio = exposed_risk / unexposed_risk,
      odds_ratio = table[positive_row, positive_col] *
        table[negative_row, negative_col] /
        (table[positive_row, negative_col] *
           table[negative_row, positive_col]),
      attributable_fraction_exposed =
        (exposed_risk - unexposed_risk) / exposed_risk,
      population_attributable_fraction =
        (sum(table[, positive_col]) / sum(table) - unexposed_risk) /
        (sum(table[, positive_col]) / sum(table)))
    epi_point <- unlist(epi$point_estimates, use.names = TRUE)

    valid[[iteration]] <-
      isTRUE(all.equal(
        diagnostic$estimates, diagnostic_oracle,
        tolerance = 1e-14, check.attributes = TRUE)) &&
      isTRUE(all.equal(
        epi_point[names(epi_oracle)], epi_oracle,
        tolerance = 1e-14, check.attributes = TRUE))
  }
  expect_true(all(valid))
})

test_that("DP 2x2 fails closed on unproven or malformed tables", {
  raw <- matrix(c(10, 20, 30, 40), 2L)
  expect_error(ds.vertDPEpi2x2(raw), "validated")
  x <- .dp_table_fixture(raw)
  x$mechanism <- "unknown"
  expect_error(ds.vertDPEpi2x2(x), "validated")
  expect_error(
    ds.vertDPEpi2x2(.dp_table_fixture(matrix(1:6, 2L))), "2-by-2")
})

test_that("DP direct standardisation is bounded post-processing", {
  tab <- matrix(c(10, 90, 20, 80, 30, 70), nrow = 3L, byrow = TRUE,
                dimnames = list(c("young", "middle", "old"),
                                c("event", "nonevent")))
  result <- ds.vertDPDirectStandardization(
    .dp_table_fixture(tab, epsilon = 3),
    c(old = 0.5, young = 0.2, middle = 0.3), event = "event")
  expect_s3_class(result, "ds.vertDPStandardization")
  expect_equal(result$estimate, 0.23)
  expect_lte(result$mechanism_region[["lower"]], result$estimate)
  expect_gte(result$mechanism_region[["upper"]], result$estimate)
  expect_true(all(result$mechanism_region >= 0 &
                  result$mechanism_region <= 1))
  expect_equal(sum(result$weights), 1)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_identical(result$additional_server_calls, 0L)
})

test_that("DP direct standardisation normalises extreme public weights", {
  tab <- matrix(c(10, 90, 20, 80), nrow = 2L, byrow = TRUE,
                dimnames = list(c("young", "old"),
                                c("event", "nonevent")))
  result <- ds.vertDPDirectStandardization(
    .dp_table_fixture(tab, epsilon = 3), c(1e308, 1e308), event = "event")
  expect_identical(result$status, "ok")
  expect_equal(result$weights, c(0.5, 0.5))
  expect_equal(result$estimate, 0.15)
})

test_that("DP indirect standardisation is zero-cost O/E post-processing", {
  tab <- matrix(c(10, 90, 20, 80, 30, 70), nrow = 3L, byrow = TRUE,
                dimnames = list(c("young", "middle", "old"),
                                c("event", "nonevent")))
  result <- ds.vertDPIndirectStandardization(
    .dp_table_fixture(tab, epsilon = 3),
    c(old = 0.20, young = 0.05, middle = 0.10), event = "event")
  expect_s3_class(result, "ds.vertDPIndirectStandardization")
  expect_equal(result$observed_events_dp, 60)
  expect_equal(result$expected_events_dp, 35)
  expect_equal(result$estimate, 60 / 35)
  expect_lte(result$mechanism_region[["lower"]], result$estimate)
  expect_gte(result$mechanism_region[["upper"]], result$estimate)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_identical(result$additional_server_calls, 0L)
  expect_match(result$inferential_scope, "no Garwood")
})

test_that("DP indirect-standardisation region covers every small table", {
  tab <- matrix(c(2, 3, 1, 4), nrow = 2L, byrow = TRUE,
                dimnames = list(c("s1", "s2"), c("event", "nonevent")))
  rates <- c(0.1, 0.4)
  result <- ds.vertDPIndirectStandardization(
    .dp_table_fixture(tab, epsilon = 3), rates, event = "event")
  lower <- as.integer(result$count_lower)
  upper <- as.integer(result$count_upper)
  candidates <- expand.grid(lapply(seq_along(lower), function(index) {
    seq.int(lower[[index]], upper[[index]])
  }), KEEP.OUT.ATTRS = FALSE)
  ratios <- apply(candidates, 1L, function(row) {
    candidate <- matrix(as.numeric(row), nrow = 2L)
    numerator <- sum(candidate[, 1L])
    denominator <- sum(rates * rowSums(candidate))
    if (denominator > 0) numerator / denominator else NA_real_
  })
  ratios <- ratios[is.finite(ratios)]
  expect_true(length(ratios) > 0)
  expect_gte(min(ratios), result$mechanism_region[["lower"]] - 1e-12)
  expect_lte(max(ratios), result$mechanism_region[["upper"]] + 1e-12)
})

test_that("DP indirect standardisation types denominator boundaries", {
  zero <- matrix(c(0, 0, 0, 0), nrow = 2L, byrow = TRUE,
                 dimnames = list(c("s1", "s2"), c("event", "nonevent")))
  result <- ds.vertDPIndirectStandardization(
    .dp_table_fixture(zero, epsilon = 3), c(0.1, 0.2), event = "event")
  expect_identical(result$status, "dp_point_non_estimable")
  expect_true(result$mechanism_region_includes_non_estimable)

  expect_error(
    ds.vertDPIndirectStandardization(
      .dp_table_fixture(zero), c(0, 0), event = "event"),
    "expected_rates")
  expect_error(
    ds.vertDPIndirectStandardization(
      .dp_table_fixture(zero), c(s1 = 0.1, wrong = 0.2), event = "event"),
    "match")
})

test_that("random indirect-standardisation boxes conservatively enclose ratios", {
  withr::local_seed(20260807)
  valid <- logical(500L)
  for (iteration in seq_along(valid)) {
    strata <- sample.int(2L, 1L)
    event_lower <- sample(0:2, strata, replace = TRUE)
    event_upper <- event_lower + sample(0:1, strata, replace = TRUE)
    nonevent_lower <- sample(0:2, strata, replace = TRUE)
    nonevent_upper <- nonevent_lower + sample(0:1, strata, replace = TRUE)
    rates <- sample(c(0, 0.05, 0.2, 0.5, 1), strata, replace = TRUE)
    if (!any(rates > 0)) rates[[1L]] <- 0.2
    region <- .dsvert_dp_indirect_box_region(
      event_lower, event_upper, nonevent_lower, nonevent_upper, rates)
    candidates <- expand.grid(c(
      Map(seq.int, event_lower, event_upper),
      Map(seq.int, nonevent_lower, nonevent_upper)),
      KEEP.OUT.ATTRS = FALSE)
    numerator <- rowSums(candidates[, seq_len(strata), drop = FALSE])
    denominator <- vapply(seq_len(nrow(candidates)), function(index) {
      event <- as.numeric(candidates[index, seq_len(strata)])
      nonevent <- as.numeric(
        candidates[index, strata + seq_len(strata)])
      sum(rates * (event + nonevent))
    }, numeric(1L))
    estimable <- denominator > 0
    ratios <- numerator[estimable] / denominator[estimable]
    valid[[iteration]] <-
      identical(region$estimable_values, any(estimable)) &&
      (!length(ratios) ||
       (min(ratios) >= region$interval[["lower"]] - 1e-12 &&
        max(ratios) <= region$interval[["upper"]] + 1e-12)) &&
      (!any(!estimable & numerator == 0) ||
       isTRUE(region$includes_non_estimable)) &&
      (!any(!estimable & numerator > 0) ||
       isTRUE(region$includes_infinite))
  }
  expect_true(all(valid))
})

test_that("DP direct standardisation ignores empty zero-weight strata only", {
  tab <- matrix(c(0, 0, 20, 80), nrow = 2L, byrow = TRUE,
                dimnames = list(c("empty", "observed"),
                                c("event", "nonevent")))
  released <- .dp_table_fixture(tab, epsilon = 3)
  ignored <- ds.vertDPDirectStandardization(
    released, c(empty = 0, observed = 1), event = "event")
  expect_identical(ignored$status, "ok")
  expect_equal(ignored$estimate, 0.2)

  required <- ds.vertDPDirectStandardization(
    released, c(empty = 1, observed = 1), event = "event")
  expect_identical(required$status, "dp_point_non_estimable")
  expect_true(is.na(required$estimate))
})

test_that("DP 2x2 reports ratio boundaries instead of an ok non-finite point", {
  tab <- matrix(c(0, 10, 5, 5), nrow = 2L, byrow = TRUE,
                dimnames = list(c("unexposed", "exposed"),
                                c("event", "nonevent")))
  result <- ds.vertDPEpi2x2(
    .dp_table_fixture(tab, epsilon = 3),
    exposed = "exposed", event = "event")
  expect_false(identical(result$status, "ok"))
  expect_null(result$point_estimates$risk_ratio)
  expect_null(result$point_estimates$odds_ratio)
  expect_true(result$point_status$risk_ratio %in%
                c("boundary_infinite", "non_estimable_undefined_ratio"))
  expect_match(result$mechanism_region_types[["risk_ratio"]], "unbounded")
})

test_that("DP 2x2 retains each estimable point when one group is empty", {
  tab <- matrix(c(0, 0, 4, 6), nrow = 2L, byrow = TRUE,
                dimnames = list(c("unexposed", "exposed"),
                                c("event", "nonevent")))
  result <- ds.vertDPEpi2x2(
    .dp_table_fixture(tab, epsilon = 3),
    exposed = "exposed", event = "event")

  expect_identical(result$status, "dp_point_non_estimable")
  expect_equal(result$point_estimates$risk_exposed, 0.4)
  expect_null(result$point_estimates$risk_unexposed)
  expect_null(result$point_estimates$risk_difference)
  expect_null(result$point_estimates$risk_ratio)
  expect_null(result$point_estimates$odds_ratio)
  expect_identical(result$point_status$risk_exposed, "ok")
  expect_identical(
    result$point_status$risk_unexposed,
    "non_estimable_zero_group_total")
  expect_match(result$inferential_scope, "no sampling")
})

test_that("DP 2x2 gives conservative regions for undefined zero ratios", {
  tab <- matrix(c(0, 10, 0, 10), nrow = 2L, byrow = TRUE,
                dimnames = list(c("unexposed", "exposed"),
                                c("event", "nonevent")))
  result <- ds.vertDPEpi2x2(
    .dp_table_fixture(tab, epsilon = 3),
    exposed = "exposed", event = "event")
  expect_identical(result$status, "dp_point_non_estimable")
  expect_null(result$point_estimates$risk_ratio)
  expect_null(result$point_estimates$odds_ratio)
  expect_equal(result$mechanism_regions$risk_ratio[["lower"]], 0)
  expect_true(is.infinite(result$mechanism_regions$risk_ratio[["upper"]]))
})

test_that("DP uncertainty radius responds monotonically", {
  expect_gte(.dsvert_dp_simultaneous_radius(1, 10, 0.95),
             .dsvert_dp_simultaneous_radius(1, 4, 0.95))
  expect_gte(.dsvert_dp_simultaneous_radius(0.5, 4, 0.95),
             .dsvert_dp_simultaneous_radius(2, 4, 0.95))
  expect_error(.dsvert_dp_simultaneous_radius(0, 4, 0.95), "Invalid")
  expect_true(is.finite(.dsvert_dp_simultaneous_radius(
    1, 1000000, 1 - 1e-12)))
})

test_that("DP 2x2 mechanism regions enclose every integer table in the box", {
  tab <- matrix(c(3, 4, 2, 5), nrow = 2L, byrow = TRUE,
                dimnames = list(c("unexposed", "exposed"),
                                c("nonevent", "event")))
  x <- .dp_table_fixture(tab, epsilon = 3)
  result <- ds.vertDPEpi2x2(x, exposed = "exposed", event = "event")
  lower <- as.integer(result$count_lower)
  upper <- as.integer(result$count_upper)
  candidates <- expand.grid(lapply(seq_along(lower), function(index) {
    seq.int(lower[[index]], upper[[index]])
  }))
  values <- t(apply(candidates, 1L, function(row) {
    candidate <- matrix(as.numeric(row), nrow = 2L)
    exposed_risk <- if (sum(candidate[2L, ]) > 0) {
      candidate[2L, 2L] / sum(candidate[2L, ])
    } else NA_real_
    unexposed_risk <- if (sum(candidate[1L, ]) > 0) {
      candidate[1L, 2L] / sum(candidate[1L, ])
    } else NA_real_
    c(
      risk_difference = exposed_risk - unexposed_risk,
      risk_ratio = exposed_risk / unexposed_risk,
      odds_ratio = candidate[2L, 2L] * candidate[1L, 1L] /
        (candidate[2L, 1L] * candidate[1L, 2L]))
  }))
  for (effect in colnames(values)) {
    finite <- values[, effect][is.finite(values[, effect])]
    expect_gte(min(finite),
               result$mechanism_regions[[effect]][["lower"]] - 1e-12)
    expect_lte(max(finite),
               result$mechanism_regions[[effect]][["upper"]] + 1e-12)
  }
})

test_that("DP standardisation region encloses every integer table in the box", {
  tab <- matrix(c(2, 4, 3, 5), nrow = 2L, byrow = TRUE,
                dimnames = list(c("young", "old"), c("nonevent", "event")))
  x <- .dp_table_fixture(tab, epsilon = 3)
  result <- ds.vertDPDirectStandardization(
    x, c(young = 0.25, old = 0.75), event = "event")
  lower <- as.integer(result$count_lower %||% (x$table - result$simultaneous_radius))
  upper <- as.integer(result$count_upper %||% (x$table + result$simultaneous_radius))
  lower[lower < 0L] <- 0L
  candidates <- expand.grid(lapply(seq_along(lower), function(index) {
    seq.int(lower[[index]], upper[[index]])
  }))
  estimates <- apply(candidates, 1L, function(row) {
    candidate <- matrix(as.numeric(row), nrow = 2L)
    denominators <- rowSums(candidate)
    if (any(denominators == 0)) return(NA_real_)
    sum(c(0.25, 0.75) * candidate[, 2L] / denominators)
  })
  estimates <- estimates[is.finite(estimates)]
  expect_gte(min(estimates), result$mechanism_region[["lower"]] - 1e-12)
  expect_lte(max(estimates), result$mechanism_region[["upper"]] + 1e-12)
  expect_equal(
    as.numeric(result$count_lower),
    as.numeric(pmax(0, x$table - result$simultaneous_radius)))
  expect_equal(
    as.numeric(result$count_upper),
    as.numeric(x$table + result$simultaneous_radius))
  expect_match(result$inferential_scope, "no sampling")
})

test_that("random direct-standardisation boxes cover every stratum risk", {
  withr::local_seed(20260806)
  valid <- logical(250L)
  for (iteration in seq_along(valid)) {
    strata <- sample.int(5L, 1L)
    table <- matrix(
      sample(0:5, 2L * strata, replace = TRUE), nrow = strata,
      dimnames = list(paste0("s", seq_len(strata)),
                      c("outcome_a", "outcome_b")))
    event <- sample.int(2L, 1L)
    nonevent <- setdiff(1:2, event)
    event_selector <- if (iteration %% 2L) {
      colnames(table)[[event]]
    } else event
    weights <- runif(strata)
    result <- ds.vertDPDirectStandardization(
      .dp_table_fixture(table, epsilon = 3),
      weights, event = event_selector)
    finite_min <- finite_max <- numeric(strata)
    stratum_ok <- TRUE
    for (index in seq_len(strata)) {
      candidates <- expand.grid(
        event = seq.int(result$count_lower[index, event],
                        result$count_upper[index, event]),
        nonevent = seq.int(result$count_lower[index, nonevent],
                           result$count_upper[index, nonevent]))
      denominator <- rowSums(candidates)
      risks <- candidates$event[denominator > 0] /
        denominator[denominator > 0]
      finite_min[[index]] <- min(risks)
      finite_max[[index]] <- max(risks)
      stratum_ok <- stratum_ok &&
        finite_min[[index]] >=
          result$stratum_regions[index, "lower"] - 1e-12 &&
        finite_max[[index]] <=
          result$stratum_regions[index, "upper"] + 1e-12
    }
    global_min <- sum(result$weights * finite_min)
    global_max <- sum(result$weights * finite_max)
    valid[[iteration]] <- stratum_ok &&
      global_min >= result$mechanism_region[["lower"]] - 1e-12 &&
      global_max <= result$mechanism_region[["upper"]] + 1e-12
  }
  expect_true(all(valid))
})

test_that("epidemiology regions and states cover every small oriented table", {
  tables <- expand.grid(rep(list(0:2), 4L), KEEP.OUT.ATTRS = FALSE)
  covered <- states_valid <- TRUE
  coverage_failure <- state_failure <- ""
  checked <- 0L
  for (table_index in seq_len(nrow(tables))) {
    table <- matrix(
      as.numeric(tables[table_index, ]), nrow = 2L,
      dimnames = list(c("r0", "r1"), c("c0", "c1")))
    for (exposed in 1:2) {
      for (event in 1:2) {
        unexposed <- setdiff(1:2, exposed)
        nonevent <- setdiff(1:2, event)
        result <- ds.vertDPEpi2x2(
          .dp_table_fixture(table, epsilon = 3), exposed, event)
        candidates <- expand.grid(lapply(seq_along(result$count_lower),
                                         function(index) {
          seq.int(result$count_lower[[index]], result$count_upper[[index]])
        }), KEEP.OUT.ATTRS = FALSE)
        values <- t(apply(candidates, 1L, function(candidate) {
          value <- matrix(as.numeric(candidate), nrow = 2L)
          totals <- rowSums(value)
          exposed_risk <- if (totals[[exposed]] > 0) {
            value[exposed, event] / totals[[exposed]]
          } else NA_real_
          unexposed_risk <- if (totals[[unexposed]] > 0) {
            value[unexposed, event] / totals[[unexposed]]
          } else NA_real_
          population_risk <- if (sum(value) > 0) {
            sum(value[, event]) / sum(value)
          } else NA_real_
          c(
            risk_exposed = exposed_risk,
            risk_unexposed = unexposed_risk,
            risk_difference = exposed_risk - unexposed_risk,
            risk_ratio = exposed_risk / unexposed_risk,
            odds_ratio = value[exposed, event] *
              value[unexposed, nonevent] /
              (value[exposed, nonevent] * value[unexposed, event]),
            attributable_fraction_exposed =
              (exposed_risk - unexposed_risk) / exposed_risk,
            population_attributable_fraction =
              (population_risk - unexposed_risk) / population_risk)
        }))
        for (metric in colnames(values)) {
          observed <- values[, metric]
          interval <- result$mechanism_regions[[metric]]
          finite <- observed[is.finite(observed)]
          if ((length(finite) &&
               (min(finite) < interval[["lower"]] - 1e-12 ||
                max(finite) > interval[["upper"]] + 1e-12)) ||
              (any(is.infinite(observed) & observed > 0) &&
               !is.infinite(interval[["upper"]])) ||
              (any(is.infinite(observed) & observed < 0) &&
               !is.infinite(interval[["lower"]]))) {
            covered <- FALSE
            coverage_failure <- paste(
              "table", table_index, "row", exposed,
              "column", event, "metric", metric)
          }
        }

        exposed_total <- sum(table[exposed, ])
        unexposed_total <- sum(table[unexposed, ])
        groups_estimable <- exposed_total > 0 && unexposed_total > 0
        exposed_risk <- if (exposed_total > 0) {
          table[exposed, event] / exposed_total
        } else NA_real_
        unexposed_risk <- if (unexposed_total > 0) {
          table[unexposed, event] / unexposed_total
        } else NA_real_
        expected_ratio_status <- if (!groups_estimable) {
          "non_estimable_zero_group_total"
        } else if (unexposed_risk == 0 && exposed_risk == 0) {
          "non_estimable_undefined_ratio"
        } else if (unexposed_risk == 0) {
          "boundary_infinite"
        } else "ok"
        odds_numerator <- table[exposed, event] *
          table[unexposed, nonevent]
        odds_denominator <- table[exposed, nonevent] *
          table[unexposed, event]
        expected_odds_status <- if (!groups_estimable) {
          "non_estimable_zero_group_total"
        } else if (odds_denominator == 0 && odds_numerator == 0) {
          "non_estimable_undefined_ratio"
        } else if (odds_denominator == 0) {
          "boundary_infinite"
        } else "ok"
        if (!identical(
              result$point_status$risk_ratio,
              expected_ratio_status) ||
            !identical(
              result$point_status$odds_ratio,
              expected_odds_status) ||
            any(vapply(
              result$mechanism_regions,
              function(interval) anyNA(interval) ||
                interval[["lower"]] > interval[["upper"]],
              logical(1L)))) {
          states_valid <- FALSE
          state_failure <- paste(
            "table", table_index, "row", exposed, "column", event)
        }
        checked <- checked + 1L
      }
    }
  }
  expect_identical(checked, 324L)
  expect_true(covered, info = coverage_failure)
  expect_true(states_valid, info = state_failure)
})

test_that("DP diagnostic 2x2 uses explicit disease-by-test orientation", {
  tab <- matrix(
    c(2, 8, 9, 1), nrow = 2L, byrow = TRUE,
    dimnames = list(
      disease = c("present", "absent"),
      test = c("negative", "positive")))
  result <- ds.vertDPDiagnostic2x2(
    .dp_table_fixture(tab, epsilon = 3),
    disease_positive = "present", test_positive = "positive")

  expect_s3_class(result, "ds.vertDPDiagnostic2x2")
  expect_identical(result$orientation$row_role, "disease_status")
  expect_identical(result$orientation$column_role, "test_result")
  expect_identical(
    result$orientation$disease_positive,
    list(index = 1L, level = "present"))
  expect_identical(
    result$orientation$test_positive,
    list(index = 2L, level = "positive"))
  expect_equal(result$estimates[c("sensitivity", "specificity")],
               c(sensitivity = 0.8, specificity = 0.9))
  expect_equal(result$estimates[["ppv"]], 8 / 9)
  expect_equal(result$estimates[["npv"]], 9 / 11)
  expect_equal(result$estimates[["prevalence"]], 0.5)
  expect_equal(result$estimates[["accuracy"]], 0.85)
  expect_equal(result$estimates[["balanced_accuracy"]], 0.85)
  expect_equal(result$estimates[["f1_score"]], 16 / 19)
  expect_equal(result$estimates[["lr_positive"]], 8)
  expect_equal(result$estimates[["lr_negative"]], 2 / 9)
  expect_equal(result$estimates[["diagnostic_odds_ratio"]], 36)
  expect_null(result$p_value)
  expect_identical(
    result$oriented_table,
    matrix(c(8, 2, 1, 9), nrow = 2L, byrow = TRUE,
           dimnames = list(
             disease = c("positive", "negative"),
             test = c("positive", "negative"))))
  expect_match(result$uncertainty_scope, "sampling uncertainty excluded")
  expect_match(result$inferential_scope, "no hypothesis test")
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_identical(result$additional_server_calls, 0L)
  expect_output(print(result), "rows=disease status; columns=test result")
})

test_that("DP diagnostic 2x2 types zero, infinity, and non-estimability", {
  perfect <- matrix(
    c(10, 0, 0, 10), nrow = 2L, byrow = TRUE,
    dimnames = list(c("disease", "healthy"), c("positive", "negative")))
  result <- ds.vertDPDiagnostic2x2(
    .dp_table_fixture(perfect, epsilon = 3),
    disease_positive = 1L, test_positive = 1L)
  expect_true(is.infinite(result$estimates[["lr_positive"]]))
  expect_true(is.infinite(result$estimates[["diagnostic_odds_ratio"]]))
  expect_identical(result$continuity_correction, "none")
  expect_identical(
    result$point_status[["lr_positive"]], "boundary_infinite")
  expect_identical(
    result$point_status[["lr_negative"]], "boundary_zero")
  expect_identical(
    result$point_status[["diagnostic_odds_ratio"]], "boundary_infinite")

  empty <- perfect * 0
  undefined <- ds.vertDPDiagnostic2x2(
    .dp_table_fixture(empty, epsilon = 3),
    disease_positive = "disease", test_positive = "positive")
  expect_true(all(is.na(undefined$estimates)))
  expect_true(all(grepl("^non_estimable", undefined$point_status)))

  reversed <- matrix(
    c(0, 10, 10, 0), nrow = 2L, byrow = TRUE,
    dimnames = dimnames(perfect))
  reversed_result <- ds.vertDPDiagnostic2x2(
    .dp_table_fixture(reversed, epsilon = 3),
    disease_positive = 1L, test_positive = 1L)
  expect_identical(
    reversed_result$point_status[["lr_positive"]], "boundary_zero")
  expect_identical(
    reversed_result$point_status[["lr_negative"]], "boundary_infinite")
  expect_identical(
    reversed_result$point_status[["diagnostic_odds_ratio"]],
    "boundary_zero")
})

test_that("DP diagnostic 2x2 rejects implicit or malformed orientation", {
  tab <- matrix(
    c(8, 2, 1, 9), nrow = 2L, byrow = TRUE,
    dimnames = list(c("disease", "healthy"), c("positive", "negative")))
  x <- .dp_table_fixture(tab)
  expect_error(
    ds.vertDPDiagnostic2x2(x),
    "disease_positive")
  expect_error(
    ds.vertDPDiagnostic2x2(
      x, disease_positive = "missing", test_positive = "positive"),
    "Unknown disease_positive")
  expect_error(
    ds.vertDPDiagnostic2x2(
      x, disease_positive = "disease", test_positive = 3L),
    "test_positive")
  x$unit_aggregation_policy <- NULL
  expect_error(
    ds.vertDPDiagnostic2x2(
      x, disease_positive = 1L, test_positive = 1L),
    "validated")
  expect_false(grepl(
    "datashield|DSI::", paste(deparse(body(ds.vertDPDiagnostic2x2)),
                               collapse = "")))
})

test_that("diagnostic regions exhaustively cover every small integer box", {
  intervals <- do.call(rbind, lapply(0:2, function(lower) {
    cbind(lower = lower, upper = lower:2)
  }))
  boxes <- expand.grid(rep(list(seq_len(nrow(intervals))), 4L))
  metric_names <- names(.dsvert_dp_diagnostic_values(c(
    true_positive = 1, false_negative = 1,
    false_positive = 1, true_negative = 1)))
  coverage_ok <- flags_ok <- TRUE
  coverage_failure <- flags_failure <- ""

  for (box_index in seq_len(nrow(boxes))) {
    interval_index <- as.integer(boxes[box_index, ])
    lower <- setNames(intervals[interval_index, "lower"], c(
      "true_positive", "false_negative",
      "false_positive", "true_negative"))
    upper <- setNames(intervals[interval_index, "upper"], names(lower))
    built <- .dsvert_dp_diagnostic_regions(lower, upper)
    candidates <- expand.grid(Map(seq.int, lower, upper))
    values <- vapply(seq_len(nrow(candidates)), function(index) {
      .dsvert_dp_diagnostic_values(
        setNames(as.numeric(candidates[index, ]), names(lower)))
    }, numeric(length(metric_names)))
    rownames(values) <- metric_names

    for (metric in metric_names) {
      observed <- values[metric, ]
      finite <- observed[is.finite(observed)]
      interval <- built$regions[[metric]]
      if (length(finite) &&
          (min(finite) < interval[["lower"]] - 1e-12 ||
           max(finite) > interval[["upper"]] + 1e-12)) {
        coverage_ok <- FALSE
        coverage_failure <- paste("box", box_index, "metric", metric)
      }
      expected_flags <- c(
        includes_zero = any(is.finite(observed) & observed == 0),
        includes_infinite = any(is.infinite(observed) & observed > 0),
        includes_non_estimable = any(is.na(observed)))
      actual_flags <- unlist(
        built$flags[metric, names(expected_flags), drop = TRUE],
        use.names = TRUE)
      if (!identical(actual_flags, expected_flags)) {
        flags_ok <- FALSE
        flags_failure <- paste("box", box_index, "metric", metric)
      }
    }
  }
  expect_true(coverage_ok, info = coverage_failure)
  expect_true(flags_ok, info = flags_failure)
})

test_that("combined DP epidemiology uses exact binomial limits", {
  alpha <- 0.037
  expect_equal(
    .dsvert_dp_cp_interval(0, 0, alpha),
    c(lower = 0, upper = 1), tolerance = 0)
  for (events in 0:8) {
    for (nonevents in 0:8) {
      if (events + nonevents == 0L) next
      expected <- stats::binom.test(
        events, events + nonevents, conf.level = 1 - alpha)$conf.int
      expect_equal(
        as.numeric(.dsvert_dp_cp_interval(events, nonevents, alpha)),
        as.numeric(expected), tolerance = 1e-14)
    }
  }
})

test_that("exact-binomial union encloses every integer table in a count box", {
  alpha <- 0.041
  checked <- 0L
  for (event_lower in 0:3) {
    for (event_upper in event_lower:3) {
      for (nonevent_lower in 0:3) {
        for (nonevent_upper in nonevent_lower:3) {
          union <- .dsvert_dp_cp_union_over_box(
            event_lower, event_upper, nonevent_lower, nonevent_upper, alpha)
          for (events in event_lower:event_upper) {
            for (nonevents in nonevent_lower:nonevent_upper) {
              candidate <- .dsvert_dp_cp_interval(events, nonevents, alpha)
              expect_lte(union[["lower"]], candidate[["lower"]])
              expect_gte(union[["upper"]], candidate[["upper"]])
              checked <- checked + 1L
            }
          }
        }
      }
    }
  }
  expect_identical(checked, 400L)
})

test_that("exact-binomial base intervals meet their finite-sample coverage", {
  alpha <- 0.05
  minimum_coverage <- 1
  for (size in 1:25) {
    intervals <- vapply(0:size, function(events) {
      .dsvert_dp_cp_interval(events, size - events, alpha)
    }, numeric(2L))
    for (probability in seq(0, 1, by = 0.025)) {
      covered <- intervals["lower", ] <= probability &
        intervals["upper", ] >= probability
      coverage <- sum(stats::dbinom(0:size, size, probability)[covered])
      minimum_coverage <- min(minimum_coverage, coverage)
      expect_gte(coverage, 1 - alpha - 1e-14)
    }
  }
  expect_gte(minimum_coverage, 1 - alpha - 1e-14)
})

test_that("combined DP epidemiology reports its full coverage contract", {
  tab <- matrix(
    c(40, 20, 10, 30), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("event", "nonevent")))
  released <- .dp_table_fixture(tab, epsilon = 2)
  result <- ds.vertDPEpi2x2Inference(
    released, exposed = "exposed", event = "event")

  expect_s3_class(result, "ds.vertDPEpi2x2Inference")
  expect_identical(result$status, "ok")
  expect_identical(result$combined_region_status, "ok")
  expect_equal(result$coverage_lower_bound, 0.95, tolerance = 0)
  expect_equal(result$mechanism_level, 0.975, tolerance = 1e-15)
  expect_equal(sum(result$alpha_allocation[c(
    "mechanism", "sampling_familywise")]),
    result$alpha_allocation[["total"]], tolerance = 1e-15)
  expect_equal(
    3 * result$alpha_allocation[["each_of_three_sampling_intervals"]],
    result$alpha_allocation[["sampling_familywise"]], tolerance = 1e-15)
  expect_equal(
    result$base_sampling_interval_level,
    1 - result$alpha_allocation[["each_of_three_sampling_intervals"]],
    tolerance = 1e-15)
  expect_false(result$confidential_count_integer_box$empty)
  expect_true(all(
    result$confidential_count_integer_box$lower <= tab &
      result$confidential_count_integer_box$upper >= tab))
  expect_true(all(result$combined_regions$risk_exposed >= 0))
  expect_lte(result$combined_regions$risk_exposed[["upper"]], 1)
  expect_gte(result$combined_regions$risk_difference[["lower"]], -1)
  expect_lte(result$combined_regions$risk_difference[["upper"]], 1)
  expect_gte(result$combined_regions$risk_ratio[["lower"]], 0)
  expect_identical(result$additional_server_calls, 0L)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_match(result$uncertainty_scope, "sampling uncertainty")
  expect_match(result$inferential_scope, "no hypothesis test")
  expect_false(grepl(
    "datashield|DSI::|aggregate\\(",
    paste(deparse(body(ds.vertDPEpi2x2Inference)), collapse = "")))
})

test_that("combined epidemiology transforms simultaneous base-risk regions", {
  exposed <- c(lower = 0.2, upper = 0.7)
  unexposed <- c(lower = 0.1, upper = 0.4)
  risk_difference <- c(
    lower = exposed[["lower"]] - unexposed[["upper"]],
    upper = exposed[["upper"]] - unexposed[["lower"]])
  risk_ratio <- .dsvert_dp_probability_ratio_region(exposed, unexposed)
  odds_ratio <- .dsvert_dp_nonnegative_ratio_region(
    .dsvert_dp_odds_region(exposed), .dsvert_dp_odds_region(unexposed))

  for (risk_exposed in seq(exposed[["lower"]], exposed[["upper"]],
                           length.out = 21L)) {
    for (risk_unexposed in seq(unexposed[["lower"]],
                               unexposed[["upper"]], length.out = 21L)) {
      difference <- risk_exposed - risk_unexposed
      ratio <- risk_exposed / risk_unexposed
      odds <- (risk_exposed / (1 - risk_exposed)) /
        (risk_unexposed / (1 - risk_unexposed))
      expect_gte(difference, risk_difference[["lower"]] - 1e-15)
      expect_lte(difference, risk_difference[["upper"]] + 1e-15)
      expect_gte(ratio, risk_ratio[["lower"]] - 1e-15)
      expect_lte(ratio, risk_ratio[["upper"]] + 1e-15)
      expect_gte(odds, odds_ratio[["lower"]] - 1e-15)
      expect_lte(odds, odds_ratio[["upper"]] + 1e-15)
    }
  }
})

test_that("combined inference derives arbitrary Gaussian confidence safely", {
  tab <- matrix(
    c(40, 20, 10, 30), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("event", "nonevent")))
  released <- .dp_vector_table_gaussian_fixture(tab)
  default <- ds.vertDPEpi2x2Inference(
    released, exposed = "exposed", event = "event")
  expect_equal(default$mechanism_level, 0.975, tolerance = 1e-15)
  expect_identical(default$combined_region_status, "ok")

  result <- ds.vertDPEpi2x2Inference(
    released, exposed = "exposed", event = "event", level = 0.90)
  expect_equal(result$mechanism_level, 0.95, tolerance = 1e-15)
  expect_identical(result$combined_region_status, "ok")

  fractional <- released
  fractional$table[1L, 1L] <- fractional$counts[[1L]] <- 40.5
  vacuous <- ds.vertDPEpi2x2Inference(
    fractional, exposed = "exposed", event = "event", level = 0.90)
  expect_identical(
    vacuous$combined_region_status,
    "vacuous_empty_integer_mechanism_box")
  expect_equal(vacuous$combined_regions$risk_exposed,
               c(lower = 0, upper = 1), tolerance = 0)
  expect_equal(vacuous$combined_regions$risk_unexposed,
               c(lower = 0, upper = 1), tolerance = 0)
  expect_equal(vacuous$combined_regions$risk_difference,
               c(lower = -1, upper = 1), tolerance = 0)
  expect_equal(vacuous$combined_regions$risk_ratio,
               c(lower = 0, upper = Inf), tolerance = 0)
})

test_that("combined inference rejects malformed releases and allocations", {
  tab <- matrix(
    c(4, 2, 1, 3), nrow = 2L, byrow = TRUE,
    dimnames = list(c("unexposed", "exposed"),
                    c("event", "nonevent")))
  released <- .dp_table_fixture(tab)
  forged <- released
  forged$randomness <- "analyst-controlled"
  expect_error(ds.vertDPEpi2x2Inference(forged), "validated")
  expect_error(
    ds.vertDPEpi2x2Inference(released, level = 1), "in \\(0, 1\\)")
  expect_error(
    ds.vertDPEpi2x2Inference(released, mechanism_alpha_share = 0),
    "in \\(0, 1\\)")
  expect_error(
    .dsvert_dp_integer_count_range(0.2, 0.3),
    "no representable integer")
})

test_that("combined diagnostic inference covers all reported point metrics", {
  tab <- matrix(
    c(400, 100, 80, 420), nrow = 2L, byrow = TRUE,
    dimnames = list(c("disease", "healthy"),
                    c("positive", "negative")))
  released <- .dp_table_fixture(tab, epsilon = 3)
  result <- ds.vertDPDiagnostic2x2Inference(
    released, disease_positive = "disease", test_positive = "positive")

  expect_s3_class(result, "ds.vertDPDiagnostic2x2Inference")
  expect_identical(result$status, "ok")
  expect_identical(result$combined_region_status, "ok")
  expect_equal(result$coverage_lower_bound, 0.95, tolerance = 0)
  expect_equal(result$mechanism_level, 0.975, tolerance = 1e-15)
  expect_equal(
    6 * result$alpha_allocation[["each_of_six_sampling_intervals"]],
    result$alpha_allocation[["sampling_familywise"]], tolerance = 1e-15)
  expect_identical(names(result$combined_regions),
                   names(result$point_estimates))
  for (metric in names(result$point_estimates)) {
    point <- result$point_estimates[[metric]]
    region <- result$combined_regions[[metric]]
    expect_false(is.na(point), info = metric)
    expect_gte(point, region[["lower"]] - 1e-14)
    expect_lte(point, region[["upper"]] + 1e-14)
  }
  expect_false(any(result$combined_region_includes_non_estimable))
  expect_identical(result$additional_server_calls, 0L)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_match(result$uncertainty_scope, "sampling uncertainty")
  expect_match(result$inferential_scope, "no hypothesis test")
  expect_false(grepl(
    "datashield|DSI::|aggregate\\(",
    paste(deparse(body(ds.vertDPDiagnostic2x2Inference)), collapse = "")))
})

test_that("combined diagnostic transformations enclose their base regions", {
  sensitivity <- c(lower = 0.2, upper = 0.8)
  specificity <- c(lower = 0.3, upper = 0.9)
  ppv <- c(lower = 0.1, upper = 0.7)
  balanced <- .dsvert_dp_probability_mean_region(
    sensitivity, specificity)
  f1 <- .dsvert_dp_probability_harmonic_region(sensitivity, ppv)
  fpr <- .dsvert_dp_probability_complement_region(specificity)
  fnr <- .dsvert_dp_probability_complement_region(sensitivity)
  lr_positive <- .dsvert_dp_probability_ratio_region(sensitivity, fpr)
  lr_negative <- .dsvert_dp_probability_ratio_region(fnr, specificity)
  diagnostic_odds_ratio <- .dsvert_dp_nonnegative_ratio_region(
    .dsvert_dp_odds_region(sensitivity), .dsvert_dp_odds_region(fpr))

  for (sens in seq(0.2, 0.8, length.out = 13L)) {
    for (spec in seq(0.3, 0.9, length.out = 13L)) {
      expect_true((sens + spec) / 2 >= balanced[["lower"]] - 1e-15)
      expect_true((sens + spec) / 2 <= balanced[["upper"]] + 1e-15)
      expect_true(sens / (1 - spec) >= lr_positive[["lower"]] - 1e-15)
      expect_true(sens / (1 - spec) <= lr_positive[["upper"]] + 1e-15)
      expect_true((1 - sens) / spec >= lr_negative[["lower"]] - 1e-15)
      expect_true((1 - sens) / spec <= lr_negative[["upper"]] + 1e-15)
      value <- (sens * spec) / ((1 - sens) * (1 - spec))
      expect_true(value >= diagnostic_odds_ratio[["lower"]] - 1e-15)
      expect_true(value <= diagnostic_odds_ratio[["upper"]] + 1e-15)
    }
    for (precision in seq(0.1, 0.7, length.out = 13L)) {
      value <- 2 * sens * precision / (sens + precision)
      expect_true(value >= f1[["lower"]] - 1e-15)
      expect_true(value <= f1[["upper"]] + 1e-15)
    }
  }
})

test_that("combined diagnostic inference types zero denominators", {
  tab <- matrix(
    c(0, 0, 0, 20), nrow = 2L, byrow = TRUE,
    dimnames = list(c("disease", "healthy"),
                    c("positive", "negative")))
  result <- ds.vertDPDiagnostic2x2Inference(
    .dp_table_fixture(tab, epsilon = 3),
    disease_positive = "disease", test_positive = "positive")
  expect_true(result$combined_region_includes_non_estimable[["sensitivity"]])
  expect_true(result$combined_region_includes_non_estimable[["ppv"]])
  expect_true(result$combined_region_includes_non_estimable[["f1_score"]])
  expect_match(result$combined_region_status, "non_estimable")
  expect_equal(result$combined_regions$sensitivity,
               c(lower = 0, upper = 1), tolerance = 0)
  expect_gte(result$combined_regions$lr_positive[["lower"]], 0)
})

test_that("combined diagnostic inference types every likelihood-ratio 0/0", {
  cases <- list(
    no_test_positives = list(
      table = matrix(
        c(0, 5, 0, 5), nrow = 2L, byrow = TRUE,
        dimnames = list(c("disease", "healthy"),
                        c("positive", "negative"))),
      undefined = c("lr_positive", "diagnostic_odds_ratio")),
    no_test_negatives = list(
      table = matrix(
        c(5, 0, 5, 0), nrow = 2L, byrow = TRUE,
        dimnames = list(c("disease", "healthy"),
                        c("positive", "negative"))),
      undefined = c("lr_negative", "diagnostic_odds_ratio")))

  for (case in cases) {
    result <- ds.vertDPDiagnostic2x2Inference(
      .dp_table_fixture(case$table, epsilon = 3),
      disease_positive = "disease", test_positive = "positive")
    expect_true(all(
      result$combined_region_includes_non_estimable[case$undefined]),
      info = paste(case$undefined, collapse = ", "))
    expect_match(result$combined_region_status, "non_estimable")
  }
})

test_that("diagnostic ratio non-estimability flags match an exact oracle", {
  metrics <- c("lr_positive", "lr_negative", "diagnostic_odds_ratio")
  candidates <- expand.grid(rep(list(0:2), 4L))
  for (index in seq_len(nrow(candidates))) {
    cells <- as.numeric(candidates[index, ])
    names(cells) <- c(
      "true_positive", "false_negative", "false_positive", "true_negative")
    expected <- is.na(.dsvert_dp_diagnostic_values(cells)[metrics])
    observed <- .dsvert_dp_diagnostic_nonestimable_flags(cells)[metrics]
    expect_identical(
      unname(observed), unname(expected),
      info = paste(cells, collapse = ","))
  }
})

test_that("combined diagnostic inference derives Gaussian confidence", {
  tab <- matrix(
    c(40, 10, 8, 42), nrow = 2L, byrow = TRUE,
    dimnames = list(c("disease", "healthy"),
                    c("positive", "negative")))
  released <- .dp_vector_table_gaussian_fixture(tab)
  default <- ds.vertDPDiagnostic2x2Inference(
    released, disease_positive = "disease", test_positive = "positive")
  expect_equal(default$mechanism_level, 0.975, tolerance = 1e-15)
  result <- ds.vertDPDiagnostic2x2Inference(
    released, disease_positive = "disease", test_positive = "positive",
    level = 0.90)
  expect_equal(result$mechanism_level, 0.95, tolerance = 1e-15)

  fractional <- released
  fractional$table[1L, 1L] <- fractional$counts[[1L]] <- 40.5
  vacuous <- ds.vertDPDiagnostic2x2Inference(
    fractional, disease_positive = "disease", test_positive = "positive",
    level = 0.90)
  expect_identical(
    vacuous$combined_region_status,
    "vacuous_empty_integer_mechanism_box")
  expect_true(all(vapply(
    vacuous$base_sampling_regions,
    identical, logical(1L), c(lower = 0, upper = 1))))
  expect_equal(vacuous$combined_regions$lr_positive,
               c(lower = 0, upper = Inf), tolerance = 0)
})

test_that("combined diagnostic inference validates release and orientation", {
  tab <- matrix(
    c(8, 2, 1, 9), nrow = 2L, byrow = TRUE,
    dimnames = list(c("disease", "healthy"),
                    c("positive", "negative")))
  released <- .dp_table_fixture(tab)
  forged <- released
  forged$randomness <- "analyst-controlled"
  expect_error(
    ds.vertDPDiagnostic2x2Inference(
      forged, disease_positive = "disease", test_positive = "positive"),
    "validated")
  expect_error(
    ds.vertDPDiagnostic2x2Inference(
      released, disease_positive = "unknown", test_positive = "positive"),
    "Unknown disease_positive")
  expect_error(
    ds.vertDPDiagnostic2x2Inference(
      released, disease_positive = "disease", test_positive = "positive",
      mechanism_alpha_share = 1),
    "in \\(0, 1\\)")
})

test_that("combined direct standardization propagates exact stratum intervals", {
  tab <- matrix(
    c(80, 120, 50, 150, 30, 170), nrow = 3L, byrow = TRUE,
    dimnames = list(c("young", "middle", "old"),
                    c("event", "nonevent")))
  released <- .dp_table_fixture(tab, epsilon = 3)
  weights <- c(young = 0.2, middle = 0.3, old = 0.5)
  result <- ds.vertDPDirectStandardizationInference(
    released, standard_weights = weights, event = "event")

  expect_s3_class(result, "ds.vertDPDirectStandardizationInference")
  expect_identical(result$status, "ok")
  expect_identical(result$combined_region_status, "ok")
  expect_equal(result$coverage_lower_bound, 0.95, tolerance = 0)
  expect_equal(result$mechanism_level, 0.975, tolerance = 1e-15)
  expect_identical(result$positive_weight_stratum_count, 3L)
  expect_equal(
    3 * result$alpha_allocation[["each_positive_weight_stratum"]],
    result$alpha_allocation[["sampling_familywise"]], tolerance = 1e-15)
  expect_equal(
    result$combined_region,
    c(
      lower = sum(result$weights *
                    result$stratum_combined_regions[, "lower"]),
      upper = sum(result$weights *
                    result$stratum_combined_regions[, "upper"])),
    tolerance = 0)
  expect_gte(result$estimate, result$combined_region[["lower"]] - 1e-14)
  expect_lte(result$estimate, result$combined_region[["upper"]] + 1e-14)
  expect_false(any(result$positive_weight_zero_sample_possible))
  expect_identical(result$additional_server_calls, 0L)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_match(result$uncertainty_scope, "sampling uncertainty")
  expect_match(result$inferential_scope, "no hypothesis test")
  expect_false(grepl(
    "datashield|DSI::|aggregate\\(",
    paste(deparse(body(ds.vertDPDirectStandardizationInference)),
          collapse = "")))
})

test_that("zero-weight standardization strata consume no sampling alpha", {
  tab <- matrix(
    c(20, 30, 10, 40, 0, 0), nrow = 3L, byrow = TRUE,
    dimnames = list(c("a", "b", "excluded"),
                    c("event", "nonevent")))
  result <- ds.vertDPDirectStandardizationInference(
    .dp_table_fixture(tab, epsilon = 3),
    standard_weights = c(a = 1, b = 1, excluded = 0), event = "event")
  expect_identical(result$positive_weight_stratum_count, 2L)
  expect_equal(result$weights, c(0.5, 0.5, 0), tolerance = 0)
  expect_equal(result$stratum_combined_regions["excluded", ],
               c(lower = 0, upper = 1), tolerance = 0)
  expect_equal(
    2 * result$alpha_allocation[["each_positive_weight_stratum"]],
    result$alpha_allocation[["sampling_familywise"]], tolerance = 1e-15)

  empty_positive <- ds.vertDPDirectStandardizationInference(
    .dp_table_fixture(tab, epsilon = 3),
    standard_weights = c(a = 0, b = 0, excluded = 1), event = "event")
  expect_true(empty_positive$positive_weight_zero_sample_possible[[1L]])
  expect_match(empty_positive$combined_region_status, "zero_stratum")
  expect_equal(empty_positive$combined_region,
               c(lower = 0, upper = 1), tolerance = 0)
})

test_that("combined direct standardization derives Gaussian confidence", {
  tab <- matrix(
    c(40, 20, 10, 30), nrow = 2L, byrow = TRUE,
    dimnames = list(c("a", "b"), c("event", "nonevent")))
  released <- .dp_vector_table_gaussian_fixture(tab)
  default <- ds.vertDPDirectStandardizationInference(
    released, c(a = 0.4, b = 0.6), event = "event")
  expect_equal(default$mechanism_level, 0.975, tolerance = 1e-15)
  result <- ds.vertDPDirectStandardizationInference(
    released, c(a = 0.4, b = 0.6), event = "event", level = 0.90)
  expect_equal(result$mechanism_level, 0.95, tolerance = 1e-15)

  fractional <- released
  fractional$table[1L, 1L] <- fractional$counts[[1L]] <- 40.5
  vacuous <- ds.vertDPDirectStandardizationInference(
    fractional, c(a = 0.4, b = 0.6),
    event = "event", level = 0.90)
  expect_identical(
    vacuous$combined_region_status,
    "vacuous_empty_integer_mechanism_box")
  expect_equal(vacuous$combined_region,
               c(lower = 0, upper = 1), tolerance = 0)
})

test_that("combined direct standardization validates release and inputs", {
  tab <- matrix(
    c(8, 2, 1, 9), nrow = 2L, byrow = TRUE,
    dimnames = list(c("a", "b"), c("event", "nonevent")))
  released <- .dp_table_fixture(tab)
  forged <- released
  forged$randomness <- "analyst-controlled"
  expect_error(
    ds.vertDPDirectStandardizationInference(
      forged, c(a = 0.5, b = 0.5), event = "event"),
    "validated")
  expect_error(
    ds.vertDPDirectStandardizationInference(
      released, c(a = 1), event = "event"),
    "match the strata")
  expect_error(
    ds.vertDPDirectStandardizationInference(
      released, c(a = 0.5, b = 0.5), event = "event",
      mechanism_alpha_share = 0),
      "in \\(0, 1\\)")
})

test_that("combined indirect standardization matches the Garwood envelope", {
  tab <- matrix(
    c(10, 90, 20, 80, 30, 70), nrow = 3L, byrow = TRUE,
    dimnames = list(c("young", "middle", "old"),
                    c("event", "nonevent")))
  released <- .dp_table_fixture(tab, epsilon = 3)
  released$capsule_id <- strrep("a", 64L)
  released$manifest_sha256 <- strrep("b", 64L)
  released$final_vector_root <- strrep("c", 64L)
  result <- ds.vertDPIndirectStandardizationInference(
    released, c(old = 0.20, young = 0.05, middle = 0.10),
    event = "event")

  expect_s3_class(result, "ds.vertDPIndirectStandardizationInference")
  expect_identical(result$status, "ok")
  expect_identical(result$combined_region_status, "ok")
  expect_equal(result$expected_rates, c(0.05, 0.10, 0.20), tolerance = 0)
  expect_equal(result$coverage_lower_bound, 0.95, tolerance = 0)
  expect_equal(result$mechanism_level, 0.975, tolerance = 1e-15)
  expect_equal(result$poisson_sampling_level, 0.975, tolerance = 1e-15)
  expect_equal(
    sum(result$alpha_allocation[c("mechanism", "sampling_familywise")]),
    result$alpha_allocation[["total"]], tolerance = 1e-15)
  expect_equal(
    result$alpha_allocation[["one_poisson_garwood_interval"]],
    result$alpha_allocation[["sampling_familywise"]], tolerance = 0)

  observed <- result$confidential_observed_events_integer_range
  expected <- result$confidential_expected_events_range
  sampling_alpha <- result$alpha_allocation[["sampling_familywise"]]
  oracle <- c(
    lower = 0.5 * stats::qchisq(
      sampling_alpha / 2, 2 * observed[["lower"]]) /
      expected[["upper"]],
    upper = 0.5 * stats::qchisq(
      1 - sampling_alpha / 2, 2 * (observed[["upper"]] + 1)) /
      expected[["lower"]])
  expect_equal(result$combined_region, oracle, tolerance = 1e-15)
  expect_gte(result$estimate, result$combined_region[["lower"]])
  expect_lte(result$estimate, result$combined_region[["upper"]])
  expect_identical(result$p_values, NULL)
  expect_identical(result$hypothesis_tests, NULL)
  expect_identical(result$additional_server_calls, 0L)
  expect_identical(result$additional_random_draws, 0L)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_false(result$new_release)
  expect_match(result$coverage_method, "Garwood")
  expect_match(result$inferential_scope, "no p-value", fixed = TRUE)
  expect_match(result$inferential_scope, "no causal", fixed = TRUE)
  expect_identical(
    result$source_release_provenance$capsule_id, released$capsule_id)
  expect_identical(
    result$source_release_provenance$manifest_sha256,
    released$manifest_sha256)
  expect_identical(
    result$source_release_provenance$final_vector_root,
    released$final_vector_root)
})

test_that("combined indirect Garwood envelope covers every small integer table", {
  tab <- matrix(
    c(5, 7, 4, 8), nrow = 2L, byrow = TRUE,
    dimnames = list(c("s1", "s2"), c("event", "nonevent")))
  rates <- c(s1 = 0.2, s2 = 0.4)
  result <- ds.vertDPIndirectStandardizationInference(
    .dp_table_fixture(tab, epsilon = 5), rates, event = "event")
  box <- result$confidential_count_integer_box
  candidates <- expand.grid(lapply(seq_along(box$lower), function(index) {
    seq.int(box$lower[[index]], box$upper[[index]])
  }), KEEP.OUT.ATTRS = FALSE)
  expect_lte(nrow(candidates), 10000L)
  alpha <- result$alpha_allocation[["sampling_familywise"]]

  for (index in seq_len(nrow(candidates))) {
    candidate <- matrix(as.numeric(candidates[index, ]), nrow = 2L)
    observed <- sum(candidate[, 1L])
    expected <- sum(rates * rowSums(candidate))
    oracle <- c(
      lower = if (observed == 0) 0 else
        0.5 * stats::qchisq(alpha / 2, 2 * observed) / expected,
      upper = 0.5 * stats::qchisq(
        1 - alpha / 2, 2 * (observed + 1)) / expected)
    expect_true(
      oracle[["lower"]] >= result$combined_region[["lower"]] - 1e-14,
      info = paste(as.numeric(candidate), collapse = ","))
    expect_true(
      oracle[["upper"]] <= result$combined_region[["upper"]] + 1e-14,
      info = paste(as.numeric(candidate), collapse = ","))
  }
})

test_that("combined indirect inference types zero and vacuous boundaries", {
  zero_observed <- matrix(
    c(0, 20, 0, 30), nrow = 2L, byrow = TRUE,
    dimnames = list(c("s1", "s2"), c("event", "nonevent")))
  zero <- ds.vertDPIndirectStandardizationInference(
    .dp_table_fixture(zero_observed, epsilon = 20),
    c(s1 = 0.1, s2 = 0.2), event = "event")
  expect_identical(zero$status, "boundary_zero")
  expect_identical(zero$combined_region_status,
                   "ok_includes_zero_observed_events")
  expect_true(zero$observed_count_box_includes_zero)
  expect_false(zero$expected_denominator_box_includes_zero)
  expect_equal(zero$combined_region[["lower"]], 0, tolerance = 0)
  expect_true(is.finite(zero$combined_region[["upper"]]))

  zero_expected <- matrix(
    c(0, 0, 5, 5), nrow = 2L, byrow = TRUE,
    dimnames = list(c("s1", "s2"), c("event", "nonevent")))
  denominator <- ds.vertDPIndirectStandardizationInference(
    .dp_table_fixture(zero_expected, epsilon = 5),
    c(s1 = 1, s2 = 0), event = "event")
  expect_true(denominator$expected_denominator_box_includes_zero)
  expect_identical(
    denominator$combined_region_status,
    "vacuous_expected_denominator_includes_zero")
  expect_equal(denominator$combined_region,
               c(lower = 0, upper = Inf), tolerance = 0)

  released <- .dp_vector_table_gaussian_fixture(zero_observed)
  released$table[1L, 1L] <- released$counts[[1L]] <- 0.5
  empty <- ds.vertDPIndirectStandardizationInference(
    released, c(s1 = 0.1, s2 = 0.2), event = "event", level = 0.90)
  expect_identical(
    empty$combined_region_status,
    "vacuous_empty_integer_mechanism_box")
  expect_equal(empty$combined_region,
               c(lower = 0, upper = Inf), tolerance = 0)
})

test_that("combined indirect inference validates inputs and has no remote route", {
  tab <- matrix(
    c(8, 2, 1, 9), nrow = 2L, byrow = TRUE,
    dimnames = list(c("s1", "s2"), c("event", "nonevent")))
  released <- .dp_table_fixture(tab)
  forged <- released
  forged$randomness <- "analyst-controlled"
  expect_error(
    ds.vertDPIndirectStandardizationInference(
      forged, c(s1 = 0.1, s2 = 0.2), event = "event"),
    "validated")
  expect_error(
    ds.vertDPIndirectStandardizationInference(
      released, c(s1 = 0.1, wrong = 0.2), event = "event"),
    "match")
  expect_error(
    ds.vertDPIndirectStandardizationInference(
      released, c(s1 = 0.1, s2 = 0.2), event = "event",
      mechanism_alpha_share = 1),
    "in \\(0, 1\\)")

  source <- paste(readLines(.dsvert_client_source_file(
    "ds.vertDPStandardizationInference.R"),
    warn = FALSE), collapse = "\n")
  expect_false(grepl(
    "datashield\\.aggregate|\\.dsvert_fanout|\\.dsvert_dp_datasources|DSI::",
    source))
  expect_false(grepl("runif\\(|rnorm\\(|sample\\(", source))
})

test_that("stratified causal standardization matches the saturated g-formula", {
  tab <- matrix(
    c(10, 90, 20, 80, 30, 70, 50, 50), nrow = 4L, byrow = TRUE,
    dimnames = list(
      c("s1_control", "s1_treated", "s2_control", "s2_treated"),
      c("event", "nonevent")))
  strata <- c("s1", "s1", "s2", "s2")
  treatment <- c("control", "treated", "control", "treated")
  weights <- c(s1 = 0.25, s2 = 0.75)
  result <- ds.vertDPCausalStandardization(
    .dp_vector_table_fixture(tab), strata, treatment, "treated", weights,
    event = "event")

  expected <- c(
    risk_treated = 0.25 * 0.2 + 0.75 * 0.5,
    risk_control = 0.25 * 0.1 + 0.75 * 0.3)
  expected <- c(
    expected,
    risk_difference = expected[["risk_treated"]] -
      expected[["risk_control"]],
    risk_ratio = expected[["risk_treated"]] /
      expected[["risk_control"]],
    odds_ratio =
      (expected[["risk_treated"]] / (1 - expected[["risk_treated"]])) /
      (expected[["risk_control"]] / (1 - expected[["risk_control"]])))
  expect_equal(unlist(result$point_estimates), expected, tolerance = 1e-15)
  expect_identical(result$status, "ok")
  expect_identical(
    result$identification_assumptions,
    c("consistency", "conditional_exchangeability_within_public_strata",
      "positivity", "no_interference", "correct_public_row_mapping",
      "fixed_valid_target_population_weights"))
  expect_identical(result$additional_server_calls, 0L)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  for (name in names(expected)) {
    expect_lte(result$mechanism_regions[[name]][["lower"]], expected[[name]])
    expect_gte(result$mechanism_regions[[name]][["upper"]], expected[[name]])
  }
})

test_that("causal mechanism regions cover every compatible integer table", {
  tab <- matrix(
    c(2, 3, 3, 2, 1, 4, 4, 1), nrow = 4L, byrow = TRUE,
    dimnames = list(paste0("r", 1:4), c("event", "nonevent")))
  strata <- c("s1", "s1", "s2", "s2")
  treatment <- c("c", "t", "c", "t")
  result <- ds.vertDPCausalStandardization(
    .dp_table_fixture(tab, epsilon = 20), strata, treatment, "t",
    c(s1 = 0.4, s2 = 0.6), event = "event")
  domains <- Map(
    seq, as.integer(ceiling(result$count_lower)),
    as.integer(floor(result$count_upper)))
  candidates <- expand.grid(domains, KEEP.OUT.ATTRS = FALSE)
  expect_lte(nrow(candidates), 10000L)
  covered <- TRUE
  failure <- ""
  for (index in seq_len(nrow(candidates))) {
    candidate <- matrix(
      as.numeric(candidates[index, ]), nrow = 4L,
      dimnames = dimnames(tab))
    row_total <- rowSums(candidate)
    if (any(row_total == 0)) next
    row_risk <- candidate[, "event"] / row_total
    treated <- 0.4 * row_risk[[2L]] + 0.6 * row_risk[[4L]]
    control <- 0.4 * row_risk[[1L]] + 0.6 * row_risk[[3L]]
    truth <- c(
      risk_treated = treated, risk_control = control,
      risk_difference = treated - control,
      risk_ratio = if (control > 0) treated / control else Inf,
      odds_ratio = if (control > 0 && control < 1) {
        (treated / (1 - treated)) / (control / (1 - control))
      } else if (treated == control) {
        NA_real_
      } else {
        Inf
      })
    for (name in names(truth)[!is.na(truth)]) {
      interval <- result$mechanism_regions[[name]]
      if (truth[[name]] < interval[["lower"]] - 1e-14 ||
          truth[[name]] > interval[["upper"]] + 1e-14) {
        covered <- FALSE
        failure <- paste("candidate", index, "estimand", name)
        break
      }
    }
    if (!covered) break
  }
  expect_true(covered, info = failure)
})

test_that("causal combined inference reports the full honest contract", {
  tab <- matrix(
    c(10, 90, 20, 80, 30, 70, 50, 50), nrow = 4L, byrow = TRUE,
    dimnames = list(
      c("s1_control", "s1_treated", "s2_control", "s2_treated"),
      c("event", "nonevent")))
  strata <- c("s1", "s1", "s2", "s2")
  treatment <- c("control", "treated", "control", "treated")
  weights <- c(s1 = 0.25, s2 = 0.75)
  result <- ds.vertDPCausalStandardizationInference(
    .dp_vector_table_fixture(tab), strata, treatment, "treated", weights,
    event = "event")
  expect_equal(result$coverage_lower_bound, 0.95, tolerance = 0)
  expect_equal(result$mechanism_level, 0.975, tolerance = 1e-15)
  expect_equal(sum(result$alpha_allocation[c(
    "mechanism", "sampling_familywise")]), 0.05, tolerance = 1e-15)
  expect_equal(
    result$base_sampling_interval_level,
    1 - result$alpha_allocation[[
      "each_positive_stratum_treatment_interval"]], tolerance = 1e-15)
  expect_identical(result$additional_server_calls, 0L)
  expect_identical(result$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_match(result$inferential_scope, "no p-value", fixed = TRUE)
  expect_identical(result$combined_region_status, "ok")

  gaussian <- ds.vertDPCausalStandardizationInference(
    .dp_vector_table_gaussian_fixture(tab), strata, treatment, "treated",
    weights, event = "event")
  expect_equal(gaussian$mechanism_level, 0.975, tolerance = 1e-15)
  expect_identical(gaussian$combined_region_status, "ok")
})

test_that("causal inference is vacuous on an integer-empty relevant DP box", {
  tab <- matrix(
    c(10, 90, 20, 80, 30, 70, 50, 50), nrow = 4L, byrow = TRUE,
    dimnames = list(paste0("r", 1:4), c("event", "nonevent")))
  released <- .dp_vector_table_gaussian_fixture(tab)
  released$table[1L, 1L] <- released$counts[[1L]] <- 10.5
  result <- ds.vertDPCausalStandardizationInference(
    released, c("s1", "s1", "s2", "s2"),
    c("c", "t", "c", "t"), "t", c(s1 = 0.5, s2 = 0.5),
    event = "event", level = 0.90)
  expect_identical(
    result$combined_region_status,
    "vacuous_empty_integer_mechanism_box")
  expect_equal(result$combined_regions$risk_treated,
               c(lower = 0, upper = 1), tolerance = 0)
  expect_equal(result$combined_regions$risk_control,
               c(lower = 0, upper = 1), tolerance = 0)
  expect_equal(result$combined_regions$risk_difference,
               c(lower = -1, upper = 1), tolerance = 0)
})

test_that("causal standardization rejects ambiguous public designs", {
  tab <- matrix(
    c(1, 9, 2, 8, 3, 7, 4, 6), nrow = 4L, byrow = TRUE,
    dimnames = list(paste0("r", 1:4), c("event", "nonevent")))
  released <- .dp_table_fixture(tab)
  valid_strata <- c("s1", "s1", "s2", "s2")
  valid_treatment <- c("c", "t", "c", "t")
  valid_weights <- c(s1 = 0.5, s2 = 0.5)
  expect_error(ds.vertDPCausalStandardization(
    released, valid_strata, rep("t", 4L), "t", valid_weights),
    "exactly two")
  expect_error(ds.vertDPCausalStandardization(
    released, valid_strata, c("c", "c", "t", "t"), "t",
    valid_weights), "combination")
  expect_error(ds.vertDPCausalStandardization(
    released, valid_strata, valid_treatment, "t", c(0.5, 0.5)),
    "named public weights")
  expect_error(ds.vertDPCausalStandardization(
    released, valid_strata, valid_treatment, "t",
    c(s1 = 0.5, wrong = 0.5)), "named public weights")
  forged <- released
  forged$noise_key_id <- ""
  expect_error(ds.vertDPCausalStandardization(
    forged, valid_strata, valid_treatment, "t", valid_weights),
    "released, validated")
})

test_that("causal standardization types zero and infinite odds boundaries", {
  all_event <- matrix(
    c(10, 0, 10, 0), nrow = 2L, byrow = TRUE,
    dimnames = list(c("control", "treated"),
                    c("event", "nonevent")))
  undefined <- ds.vertDPCausalStandardization(
    .dp_table_fixture(all_event, epsilon = 5), c("s", "s"),
    c("c", "t"), "t", c(s = 1), event = "event")
  expect_null(undefined$point_estimates$odds_ratio)
  expect_identical(
    undefined$point_status[["odds_ratio"]],
    "non_estimable_undefined_ratio")

  separated <- matrix(
    c(0, 10, 10, 0), nrow = 2L, byrow = TRUE,
    dimnames = list(c("control", "treated"),
                    c("event", "nonevent")))
  infinite <- ds.vertDPCausalStandardization(
    .dp_table_fixture(separated, epsilon = 5), c("s", "s"),
    c("c", "t"), "t", c(s = 1), event = "event")
  expect_identical(infinite$point_estimates$odds_ratio, Inf)
  expect_identical(infinite$point_status[["odds_ratio"]],
                   "boundary_infinite")
  expect_identical(infinite$status, "boundary")
})
