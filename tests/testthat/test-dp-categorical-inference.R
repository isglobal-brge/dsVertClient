.dp_categorical_exact_plan <- function(coordinate_count = 4L,
                                       sensitivity_steps = "256") {
  denominator <- "1267650600228229401496703205376"
  list(
    version = "dsvert-joint-dp-vector-laplace-plan-v3",
    sampler = .DSVERT_CLIENT_VECTOR_EXACT_SAMPLER,
    stop_bits = 128,
    stop_numerator = "1326635224458652993228324758901174720",
    uniform_bits = 128, binary_geometric_bits = 1,
    bernoulli_thresholds = list("1"),
    sensitivity_steps = sensitivity_steps,
    total_coordinate_count = as.numeric(coordinate_count),
    epsilon_effective_upper_numerator = "1",
    epsilon_effective_upper_denominator = "1",
    one_geometric_tv_numerator = "1",
    one_geometric_tv_denominator = denominator,
    tail_upper_numerator = "1", tail_upper_denominator = denominator,
    rounding_upper_numerator = "1",
    rounding_upper_denominator = denominator,
    implementation_delta_numerator = "8",
    implementation_delta_denominator = denominator,
    implementation_delta_bound = "6.310887241768095e-30",
    maximum_noise_magnitude = "1048576",
    maximum_chunk_coordinates = as.numeric(coordinate_count),
    private_stream_bytes_per_coordinate = 64,
    accounting = paste(
      "global iid discrete Laplace calibrated once to the workload joint",
      "L1 sensitivity; exact binary-geometric coupling"),
    capability_available = TRUE)
}

.dp_categorical_release <- function(
    table = matrix(c(28, 22, 22, 28), 2L, 2L,
                   dimnames = list(c("no", "yes"), c("no", "yes"))),
    capacity = 1000, root_digit = "3", exact = FALSE) {
  plan <- if (isTRUE(exact)) {
    .dp_categorical_exact_plan(length(table))
  } else NULL
  value <- list(
    released = TRUE,
    mechanism = if (isTRUE(exact)) {
      .DSVERT_CLIENT_VECTOR_EXACT_RELEASE_MECHANISM
    } else .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
    implementation = if (isTRUE(exact)) {
      .DSVERT_CLIENT_VECTOR_EXACT_BACKEND
    } else .DSVERT_CLIENT_VECTOR_BACKEND,
    backend = "exact_signed_Ring128_global_vector",
    sampler = if (isTRUE(exact)) {
      .DSVERT_CLIENT_VECTOR_EXACT_SAMPLER
    } else .DSVERT_CLIENT_VECTOR_SAMPLER,
    randomness = paste(
      "independent pinned-peer HKDF-SHA256/ChaCha20 streams;",
      "no analyst-controlled seed"),
    epsilon = 1, delta = 2^-100,
    implementation_delta = if (isTRUE(exact)) paste0(
      plan$implementation_delta_numerator, "/",
      plan$implementation_delta_denominator) else
        "1/1267650600228229401496703205376",
    one_joint_draw = isTRUE(exact), mechanism_plan = plan,
    plan_sha256 = if (isTRUE(exact)) .dsvert_vector_hash(plan) else NULL,
    capsule_mechanism = if (isTRUE(exact)) {
      list(mechanism = "discrete-laplace")
    } else NULL,
    adjacency = "add_remove_patient",
    sensitivity = 1, sensitivity_norm = "l1",
    l1_sensitivity = 1, l2_sensitivity = 1,
    sensitivity_scope = "complete_signed_biomedical_capsule_vector",
    output_lattice_bits = 8L, output_lattice_scale = 256,
    capsule_id = paste(rep("1", 64L), collapse = ""),
    manifest_sha256 = paste(rep("2", 64L), collapse = ""),
    final_vector_root = paste(rep(root_digit, 64L), collapse = ""),
    coordinate_order_sha256 = paste(rep("4", 64L), collapse = ""),
    capsule_coordinate_count = 4L,
    sticky_noise = "immutable_capsule_durable_replay_v3",
    sticky_replay = TRUE,
    privacy_epochs = c(site_a = 1, site_b = 1),
    noise_key_ids = c(site_a = "noise-a", site_b = "noise-b"),
    source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE,
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    clipped_coordinates = NA_integer_, clamp_activation_disclosed = FALSE,
    server = "site_b", row_var = "exposure", col_var = "outcome",
    row_levels = unname(rownames(table)),
    col_levels = unname(colnames(table)),
    nrow = as.integer(nrow(table)), ncol = as.integer(ncol(table)),
    counts = unname(as.numeric(table)), table = table,
    coordinate_maximum = capacity,
    unit_aggregation_policy = "consistent_joint_cell_else_exclude_v1",
    missingness_policy = "missing_or_out_of_domain_rows_are_ignored",
    artifact_l1_sensitivity = 1, artifact_l2_sensitivity = 1,
    accuracy_simultaneous_confidence = 0.95,
    accuracy_simultaneous_method = if (isTRUE(exact)) paste(
      "exact ideal one-draw two-sided-geometric tail with union bound;",
      "signed vector sampler TV deducted once; clamp inside exact GC applied")
    else paste(
      "exact ideal two-sided-geometric convolution tail with union bound;",
      "two-peer finite-sampler TV deducted; fixed-clamp range applied"),
    accuracy_additional_privacy_cost = c(epsilon = 0, delta = 0))
  class(value) <- c("ds.vertDPContingency", "list")
  value$accuracy_simultaneous_95_abs <-
    .dsvert_dp_vector_table_radius(value, 0.95)
  value
}

test_that("two-draw sampler agrees with a small brute-force oracle", {
  noise <- list(
    stop_probability = 0.5, continuation_probability = 0.5, scale = 1)
  support <- -40:40
  normalizer <- (1 - noise$continuation_probability) /
    (1 + noise$continuation_probability)
  one <- normalizer * noise$continuation_probability^abs(support)
  oracle <- vapply(-3:3, function(value) {
    sum(one * normalizer *
          noise$continuation_probability^abs(value - support))
  }, numeric(1L))
  sampled <- .dsvert_dp_chisq_with_seed(
    76123L, .dsvert_dp_chisq_sample_noise(200000L, noise))
  empirical <- vapply(-3:3, function(value) mean(sampled == value),
                      numeric(1L))
  expect_equal(empirical, oracle, tolerance = 0.002)
})

test_that("one-draw exact-GC sampler agrees with its signed ideal law", {
  noise <- list(
    stop_probability = 0.5, continuation_probability = 0.5,
    scale = 1, draw_count = 1L)
  support <- -40:40
  normalizer <- (1 - noise$continuation_probability) /
    (1 + noise$continuation_probability)
  oracle <- normalizer * noise$continuation_probability^abs(-3:3)
  sampled <- .dsvert_dp_chisq_with_seed(
    76123L, .dsvert_dp_chisq_sample_noise(200000L, noise))
  empirical <- vapply(-3:3, function(value) mean(sampled == value),
                      numeric(1L))
  expect_equal(empirical, oracle, tolerance = 0.002)

  legacy <- noise[names(noise) != "draw_count"]
  explicit_two <- c(legacy, list(draw_count = 2L))
  expect_identical(
    .dsvert_dp_chisq_with_seed(
      991L, .dsvert_dp_chisq_sample_noise(1000L, legacy)),
    .dsvert_dp_chisq_with_seed(
      991L, .dsvert_dp_chisq_sample_noise(1000L, explicit_two)))
})

test_that("exact-GC categorical inference uses one signed draw", {
  exact <- .dp_categorical_release(exact = TRUE)
  convolution <- .dp_categorical_release()
  expect_silent(.dsvert_dp_table_contract(exact))
  expect_identical(
    .dsvert_dp_chisq_mechanism_dispatch(exact),
    "signed_one_draw_discrete_laplace")
  expect_identical(
    .dsvert_dp_fisher_mechanism_dispatch(exact),
    "signed_one_draw_discrete_laplace")

  exact_noise <- .dsvert_dp_chisq_noise_contract(exact)
  convolution_noise <- .dsvert_dp_chisq_noise_contract(convolution)
  expect_identical(exact_noise$draw_count, 1L)
  expect_identical(convolution_noise$draw_count, 2L)
  expect_identical(exact_noise$stop_probability_source,
                   "signed_exact_gc_dyadic_plan")
  expect_identical(convolution_noise$stop_probability_source,
                   "declared_epsilon_exponential_law")
  signed_q <- .dsvert_dp_vector_dyadic_tail_context(
    exact$mechanism_plan)$q
  expect_identical(exact_noise$stop_probability, signed_q)
  expect_gte(exact_noise$stop_probability,
             exact_noise$stop_probability_interval[["lower"]])
  expect_lte(exact_noise$stop_probability,
             exact_noise$stop_probability_interval[["upper"]])
  expect_lt(exact_noise$stop_probability,
            convolution_noise$stop_probability)
  exact_formula <- 2 * (1 - exact_noise$variance_stop_probability) /
    exact_noise$variance_stop_probability^2 / exact$output_lattice_scale^2
  convolution_formula <- 4 * convolution_noise$continuation_probability /
    convolution_noise$stop_probability^2 /
    convolution$output_lattice_scale^2
  expect_equal(exact_noise$variance, exact_formula, tolerance = 1e-12)
  expect_gte(exact_noise$variance, exact_formula)
  expect_equal(convolution_noise$variance, convolution_formula,
               tolerance = 0)
  expect_true(exact_noise$variance_is_upper_bound)
  expect_false(convolution_noise$variance_is_upper_bound)
  expect_identical(
    exact_noise$vector_tv_upper_bound,
    .dsvert_dp_vector_sampler_tv_upper(exact$mechanism_plan, TRUE))
  expect_gt(exact_noise$numeric_reference_tv_upper_bound, 0)
  expect_gte(
    exact_noise$one_transfer_tv_upper_bound,
    exact_noise$vector_tv_upper_bound +
      exact_noise$numeric_reference_tv_upper_bound)
  expect_gte(
    exact_noise$calibration_tv_upper_bound,
    2 * exact_noise$one_transfer_tv_upper_bound)

  exactly_representable <- exact$mechanism_plan
  exactly_representable$stop_numerator <-
    "170141183460469231731687303715884105728"
  representable <- .dsvert_dp_chisq_exact_stop_probability(
    exactly_representable)
  expect_equal(representable$value, 0.5, tolerance = 0)
  expect_lt(representable$numeric_vector_tv_upper_bound, 1e-12)

  chisq <- ds.vertChisq(exact, correct = FALSE, simulations = 199L)
  fisher <- ds.vertFisher(exact, simulations = 199L)
  expect_identical(chisq$calibration,
                   .DSVERT_DP_CHISQ_EXACT_BOOTSTRAP_VERSION)
  expect_identical(fisher$calibration,
                   .DSVERT_DP_FISHER_EXACT_BOOTSTRAP_VERSION)
  expect_match(chisq$reference_mechanism, "one joint complete")
  expect_match(fisher$reference_mechanism, "one joint complete")
  expect_match(chisq$finite_sampler_handling,
               "R-parameter conversion TV propagated twice")
  expect_match(fisher$finite_sampler_handling,
               "R-parameter conversion TV propagated twice")
  expect_identical(
    chisq$mechanism_vector_tv_upper_bound,
    exact_noise$vector_tv_upper_bound)
  expect_identical(
    chisq$mechanism_numeric_reference_tv_upper_bound,
    exact_noise$numeric_reference_tv_upper_bound)
  expect_identical(
    chisq$mechanism_one_transfer_tv_upper_bound,
    exact_noise$one_transfer_tv_upper_bound)
  expect_identical(
    chisq$mechanism_reference_tv_upper_bound,
    exact_noise$calibration_tv_upper_bound)
  expect_identical(
    fisher$mechanism_numeric_reference_tv_upper_bound,
    exact_noise$numeric_reference_tv_upper_bound)
  expect_identical(
    fisher$mechanism_reference_tv_upper_bound,
    exact_noise$calibration_tv_upper_bound)
  expect_false(identical(
    chisq$bootstrap_seed_id,
    ds.vertChisq(convolution, correct = FALSE,
                 simulations = 199L)$bootstrap_seed_id))

  epi <- ds.vertDPEpi2x2(exact)
  expect_match(epi$coverage_method, "one-draw")

  wrong_exact_method <- exact
  wrong_exact_method$accuracy_simultaneous_method <-
    convolution$accuracy_simultaneous_method
  expect_error(.dsvert_dp_table_contract(wrong_exact_method),
               "released, validated ds.vertDPContingency", fixed = TRUE)
  wrong_convolution_method <- convolution
  wrong_convolution_method$accuracy_simultaneous_method <-
    exact$accuracy_simultaneous_method
  expect_error(.dsvert_dp_table_contract(wrong_convolution_method),
               "released, validated ds.vertDPContingency", fixed = TRUE)

  unsigned <- exact
  unsigned$mechanism_plan <- NULL
  expect_error(.dsvert_dp_chisq_noise_contract(unsigned),
               "lacks its signed finite-sampler plan", fixed = TRUE)
  detached <- exact
  detached$mechanism_plan$stop_bits <- 127
  expect_error(.dsvert_dp_chisq_noise_contract(detached),
               "stop probability is not representable", fixed = TRUE)
})

test_that("table-vector accuracy method text is canonical for every profile", {
  exact <- .dp_categorical_release(exact = TRUE)
  convolution <- .dp_categorical_release()
  gaussian <- paste(
    "signed fixed-work dyadic discrete-Gaussian plan v2 simultaneous",
    "95% bound; tail and CDF TV transfers already charged;",
    "fixed-clamp range applied")
  methods <- list(
    gaussian = gaussian,
    exact = exact$accuracy_simultaneous_method,
    convolution = convolution$accuracy_simultaneous_method)
  profiles <- list(
    gaussian = list(gaussian = TRUE, exact_gc = FALSE),
    exact = list(gaussian = FALSE, exact_gc = TRUE),
    convolution = list(gaussian = FALSE, exact_gc = FALSE))
  for (name in names(methods)) {
    value <- list(accuracy_simultaneous_method = methods[[name]])
    expect_true(.dsvert_dp_table_vector_accuracy_method_is_valid(
      value, profiles[[name]]), info = name)
    for (forged in c(
        paste("forged prefix", methods[[name]]),
        paste(methods[[name]], "forged suffix"))) {
      value$accuracy_simultaneous_method <- forged
      expect_false(.dsvert_dp_table_vector_accuracy_method_is_valid(
        value, profiles[[name]]), info = name)
    }
  }
  both <- paste(methods$exact, methods$convolution)
  for (profile in profiles) {
    expect_false(.dsvert_dp_table_vector_accuracy_method_is_valid(
      list(accuracy_simultaneous_method = both), profile))
  }

  forged <- exact
  forged$accuracy_simultaneous_method <- both
  consumers <- c(
    "ds.vertDPCausalStandardization",
    "ds.vertDPCausalStandardizationInference",
    "ds.vertDPROC",
    "ds.vertDPDiagnostic2x2Inference",
    "ds.vertDPDiagnostic2x2",
    "ds.vertDPEpi2x2",
    "ds.vertDPDirectStandardization",
    "ds.vertDPIndirectStandardization",
    "ds.vertDPEpi2x2Inference",
    "ds.vertDPDirectStandardizationInference",
    "ds.vertDPIndirectStandardizationInference",
    "ds.vertDPMantelHaenszel")
  expect_length(consumers, 12L)
  for (consumer in consumers) {
    expect_error(
      do.call(get(consumer, envir = asNamespace("dsVertClient")),
              list(x = forged)),
      "released, validated ds.vertDPContingency", fixed = TRUE,
      info = consumer)
  }
})

test_that("the public per-cell clamp is reproduced exactly", {
  expect_identical(
    .dsvert_dp_chisq_apply_release_clamp(
      c(0, 5, 10), c(-1.25, 1.5, 2), 10),
    c(0, 6.5, 10))
})

test_that("DP-aware inference is deterministic and preserves caller RNG", {
  release <- .dp_categorical_release()
  set.seed(90210)
  before <- .Random.seed
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(...) {
      stop("an existing release must not cause a DSI round trip")
    },
    .package = "dsVertClient")
  first <- ds.vertChisq(release, correct = FALSE, simulations = 399L)
  expect_identical(.Random.seed, before)
  second <- ds.vertChisq(release, correct = FALSE, simulations = 399L)
  expect_identical(.Random.seed, before)
  expect_identical(first$p_value, second$p_value)
  expect_identical(first$exceedances, second$exceedances)
  expect_identical(first$bootstrap_seed_id, second$bootstrap_seed_id)
  expect_identical(first$additional_server_calls, 0L)
  expect_identical(first$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_identical(first$source_dp_release, release)
  expect_identical(first$source_release$final_vector_root,
                   release$final_vector_root)
})

test_that("the character front door obtains one capsule artifact only", {
  release <- .dp_categorical_release()
  calls <- 0L
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(...) {
      calls <<- calls + 1L
      release
    },
    .package = "dsVertClient")
  result <- ds.vertChisq(
    "cohort", "exposure", "outcome", simulations = 199L,
    correct = FALSE)
  expect_identical(calls, 1L)
  expect_identical(result$additional_server_calls, 0L)
  frontdoor <- paste(deparse(body(ds.vertChisq)), collapse = "\n")
  expect_false(grepl("dsvertContingencyDS", frontdoor, fixed = TRUE))
  expect_false(grepl("pchisq", frontdoor, fixed = TRUE))
})

test_that("the DP-aware result reports MC and finite-sampler uncertainty", {
  result <- ds.vertChisq(
    .dp_categorical_release(), correct = FALSE, simulations = 499L)
  expect_s3_class(result, "ds.vertChisq")
  expect_identical(result$status, "ok")
  expect_true(result$decision_available)
  expect_gte(result$p_value, 1 / 500)
  expect_lte(result$p_value, 1)
  expect_lte(result$p_value_mc_interval[["lower"]], result$p_value)
  expect_gte(result$p_value_mc_interval[["upper"]], result$p_value)
  expect_lte(result$p_value_calibration_interval[["lower"]],
             result$p_value_mc_interval[["lower"]])
  expect_gte(result$p_value_calibration_interval[["upper"]],
             result$p_value_mc_interval[["upper"]])
  expect_match(result$inferential_contract, "not a finite-sample exact")
  expect_match(result$reference_mechanism, "two independent complete")
  expect_false(result$finite_sampler_bit_exact_reference)
})

test_that("association evidence grows with effective sample size", {
  small <- .dp_categorical_release(
    matrix(c(30, 20, 20, 30), 2L, 2L,
           dimnames = list(c("no", "yes"), c("no", "yes"))),
    root_digit = "5")
  large <- .dp_categorical_release(
    matrix(c(300, 200, 200, 300), 2L, 2L,
           dimnames = list(c("no", "yes"), c("no", "yes"))),
    root_digit = "6")
  small_result <- ds.vertChisq(
    small, correct = FALSE, simulations = 999L)
  large_result <- ds.vertChisq(
    large, correct = FALSE, simulations = 999L)
  expect_gt(large_result$statistic, small_result$statistic)
  expect_lte(large_result$p_value, small_result$p_value)
})

test_that("weakly estimable tables return no misleading p-value", {
  release <- .dp_categorical_release(
    matrix(c(2, 0, 0, 2), 2L, 2L,
           dimnames = list(c("no", "yes"), c("no", "yes"))))
  result <- ds.vertChisq(release, simulations = 99L)
  expect_false(result$decision_available)
  expect_true(is.na(result$p_value))
  expect_match(result$status, "not_tested_")
  expect_identical(result$source_dp_release, release)
  expect_output(print(result), "no p-value")
})

test_that("tampered DP artifacts and analyst seed-like inputs are rejected", {
  release <- .dp_categorical_release()
  release$counts[[1L]] <- release$counts[[1L]] + 1
  expect_error(ds.vertChisq(release, simulations = 19L),
               "released, validated ds.vertDPContingency")
  expect_error(ds.vertChisq(
    .dp_categorical_release(), server = "site_b", simulations = 19L),
    "must be omitted")
  expect_error(ds.vertChisq(
    .dp_categorical_release(), simulations = 0L),
    "positive integer")
})

test_that("same-owner cross routing cannot append ordinary Fisher inference", {
  cross <- paste(deparse(body(ds.vertChisqCross)), collapse = "\n")
  expect_false(grepl("stats::fisher.test(out$observed)", cross,
                     fixed = TRUE))
})

test_that("lowercase and epidemiology front doors reuse one validated release", {
  release <- .dp_categorical_release(
    matrix(c(40, 20, 10, 30), 2L, byrow = TRUE,
           dimnames = list(c("unexposed", "exposed"),
                           c("event", "nonevent"))))
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(...) {
      stop("an existing release must not cause a DSI round trip")
    },
    .package = "dsVertClient")
  chisq <- ds.vert.chisq(release, simulations = 19L)
  expect_s3_class(chisq, "ds.vertChisq")
  expect_identical(chisq$frontdoor, "ds.vert.chisq")
  epi <- ds.vertEpi2x2(
    chisq, exposed = "exposed", event = "event")
  expect_s3_class(epi, "ds.vertDPEpi2x2")
  expect_identical(epi$uncertainty_scope,
                   "DP mechanism noise only; sampling uncertainty excluded")
  expect_error(
    ds.vertEpi2x2(chisq, zero_correction = "always"),
    "not available for a DP contingency release")
  tampered <- chisq
  tampered$calibration <- "unknown"
  expect_error(ds.vertEpi2x2(tampered),
               "invalid DP contingency provenance")
})

test_that("DP-aware Fisher reuses one sticky release and preserves RNG", {
  release <- .dp_categorical_release(
    matrix(c(60, 20, 20, 60), 2L, 2L,
           dimnames = list(c("unexposed", "exposed"),
                           c("event", "nonevent"))))
  set.seed(20403)
  before <- .Random.seed
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(...) {
      stop("an existing release must not cause a DSI round trip")
    },
    .package = "dsVertClient")
  first <- ds.vertFisher(release, simulations = 399L)
  expect_identical(.Random.seed, before)
  second <- ds.vertFisher(release, simulations = 399L)
  expect_identical(.Random.seed, before)
  expect_identical(first$p_value, second$p_value)
  expect_identical(first$exceedances, second$exceedances)
  expect_identical(first$bootstrap_seed_id, second$bootstrap_seed_id)
  expect_identical(first$additional_server_calls, 0L)
  expect_identical(first$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))
  expect_identical(first$source_dp_release, release)
  expect_false(first$finite_sample_exact)
  expect_false(first$raw_table_fisher_exact)
  expect_null(first$conf_int)
  expect_false(first$confidence_interval_available)
  expect_match(first$inferential_contract, "not Fisher-exact")
  expect_lte(first$p_value_mc_interval[["lower"]], first$p_value)
  expect_gte(first$p_value_mc_interval[["upper"]], first$p_value)
  expect_lte(first$p_value_calibration_interval[["lower"]],
             first$p_value_mc_interval[["lower"]])
  expect_gte(first$p_value_calibration_interval[["upper"]],
             first$p_value_mc_interval[["upper"]])
})

test_that("Fisher character and lowercase front doors obtain one capsule", {
  release <- .dp_categorical_release()
  calls <- 0L
  testthat::local_mocked_bindings(
    ds.vertDPContingency = function(...) {
      calls <<- calls + 1L
      release
    },
    .package = "dsVertClient")
  direct <- ds.vertFisher(
    "cohort", "exposure", "outcome", simulations = 49L)
  expect_identical(calls, 1L)
  expect_identical(direct$additional_server_calls, 0L)

  testthat::local_mocked_bindings(
    .dsvert_datasources = function(...) {
      stop("an existing release must not resolve DSI connections")
    },
    ds.vertDPContingency = function(...) {
      stop("an existing release must not cause a DSI round trip")
    },
    .package = "dsVertClient")
  alias <- ds.vert.fisher(release, simulations = 49L)
  expect_s3_class(alias, "ds.vertFisher")
  expect_identical(alias$frontdoor, "ds.vert.fisher")
})

test_that("conditional latent sampler has the declared exact margins", {
  table <- matrix(c(45, 25, 15, 35), 2L, 2L,
                  dimnames = list(c("r0", "r1"), c("c0", "c1")))
  fit <- .dsvert_dp_fisher_fit(table, 1000)
  expect_true(fit$ok)
  draws <- .dsvert_dp_chisq_with_seed(
    94117L, .dsvert_dp_fisher_latent_tables(50000L, fit))
  expect_true(all(draws[1L, ] + draws[3L, ] ==
                    fit$integer_row_margins[[1L]]))
  expect_true(all(draws[2L, ] + draws[4L, ] ==
                    fit$integer_row_margins[[2L]]))
  expect_true(all(draws[1L, ] + draws[2L, ] ==
                    fit$integer_col_margins[[1L]]))
  expect_true(all(draws[3L, ] + draws[4L, ] ==
                    fit$integer_col_margins[[2L]]))
  oracle_mean <- fit$integer_row_margins[[1L]] *
    fit$integer_col_margins[[1L]] / fit$n
  expect_equal(mean(draws[1L, ]), oracle_mean, tolerance = 0.03)
  expect_equal(stats::var(draws[1L, ]), fit$hypergeometric_variance,
               tolerance = 0.08)
})

test_that("DP-aware conditional test separates a strong association", {
  null <- .dp_categorical_release(
    matrix(c(50, 50, 50, 50), 2L, 2L,
           dimnames = list(c("no", "yes"), c("no", "yes"))),
    root_digit = "7")
  associated <- .dp_categorical_release(
    matrix(c(80, 20, 20, 80), 2L, 2L,
           dimnames = list(c("no", "yes"), c("no", "yes"))),
    root_digit = "8")
  null_result <- ds.vertFisher(null, simulations = 999L)
  associated_result <- ds.vertFisher(associated, simulations = 999L)
  expect_gt(null_result$p_value, 0.05)
  expect_lt(associated_result$p_value, null_result$p_value)
  expect_lte(associated_result$p_value, 0.01)
  expect_gt(associated_result$signed_root_pearson, 0)
  less_result <- ds.vertFisher(
    associated, alternative = "less", simulations = 499L)
  expect_gt(less_result$p_value, associated_result$p_value)
})

test_that("conditional calibration has reasonable null behavior", {
  release <- .dp_categorical_release(
    matrix(c(50, 50, 50, 50), 2L, 2L,
           dimnames = list(c("no", "yes"), c("no", "yes"))))
  fit <- .dsvert_dp_fisher_fit(release$table, release$coordinate_maximum)
  noise <- .dsvert_dp_chisq_noise_contract(release)
  generated <- .dsvert_dp_chisq_with_seed(81277L, {
    latent <- .dsvert_dp_fisher_latent_tables(30L, fit)
    mechanism_noise <- matrix(
      .dsvert_dp_chisq_sample_noise(4L * 30L, noise), nrow = 4L)
    .dsvert_dp_chisq_apply_release_clamp(
      latent, mechanism_noise, release$coordinate_maximum)
  })
  p_values <- vapply(seq_len(ncol(generated)), function(index) {
    candidate <- release
    candidate$table <- matrix(
      generated[, index], 2L, 2L, dimnames = dimnames(release$table))
    candidate$counts <- unname(as.numeric(candidate$table))
    candidate$final_vector_root <- digest::digest(
      paste("dp-fisher-null-oracle", index, sep = "|"),
      algo = "sha256", serialize = FALSE)
    ds.vertFisher(candidate, simulations = 199L)$p_value
  }, numeric(1L))
  # This is a regression smoke test, not a proof of finite-sample level. The
  # generous binomial envelope catches grossly anti-conservative calibration
  # without making the suite depend on a fragile simulated percentage.
  expect_lte(sum(p_values <= 0.10), 8L)
  expect_gt(stats::median(p_values), 0.20)
})

test_that("Fisher fails honestly for unsupported or degenerate contracts", {
  degenerate <- .dp_categorical_release(
    matrix(c(100, 0, 0, 0), 2L, 2L,
           dimnames = list(c("no", "yes"), c("no", "yes"))))
  result <- ds.vertFisher(degenerate, simulations = 19L)
  expect_false(result$decision_available)
  expect_true(is.na(result$p_value))
  expect_match(result$status, "not_tested_")
  expect_output(print(result), "no p-value")

  non_2x2 <- .dp_categorical_release(
    matrix(c(10, 10, 10, 10, 10, 10), 3L, 2L,
           dimnames = list(c("a", "b", "c"), c("no", "yes"))))
  expect_error(
    ds.vertFisher(non_2x2, simulations = 19L),
    class = "dsvert_dp_fisher_unsupported_dimension")
  expect_error(
    .dsvert_dp_fisher_mechanism_dispatch(list(
      noise_selection = list(winner = "gaussian"))),
    class = "dsvert_dp_fisher_mechanism_not_certified")

  gaussian <- list(
    capsule_mechanism = list(
      mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_MECHANISM,
      sensitivity_norm = "l2"),
    mechanism_selection = list(),
    mechanism = .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
  expect_error(
    .dsvert_dp_chisq_mechanism_dispatch(gaussian),
    class = "dsvert_dp_chisq_mechanism_not_certified")
  expect_error(
    .dsvert_dp_fisher_mechanism_dispatch(gaussian),
    class = "dsvert_dp_fisher_mechanism_not_certified")
})

test_that("Fisher rejects tampering and has no legacy exact release route", {
  release <- .dp_categorical_release()
  release$table[1L, 1L] <- release$table[1L, 1L] + 1
  expect_error(ds.vertFisher(release, simulations = 19L),
               "released, validated ds.vertDPContingency")
  expect_error(ds.vertFisher(
    .dp_categorical_release(), server = "site_b", simulations = 19L),
    "must be omitted")
  fisher_surface <- paste(
    deparse(body(ds.vertFisher)),
    deparse(body(.dsvert_dp_fisher_from_release)), collapse = "\n")
  expect_false(grepl("dsvertContingencyDS", fisher_surface, fixed = TRUE))
  expect_false(grepl("dsvertColNamesDS", fisher_surface, fixed = TRUE))
  expect_false(grepl("fisher.test", fisher_surface, fixed = TRUE))
})
