.vector_layout_manifest <- function() {
  count <- list(owner_peer = "site_a", dataset = "cohort")
  numeric <- list(artifacts = list(
    z_key = list(owner_peer = "site_b", dataset = "remote", column = "z",
                 lower = -1, upper = 1, numeric_grid_bits = 8L),
    a_key = list(owner_peer = "site_a", dataset = "cohort", column = "age",
                 lower = 0, upper = 100, numeric_grid_bits = 8L)))
  histograms <- list(artifacts = list(
    age_hist = list(
      owner_peer = "site_a", dataset = "cohort", column = "age",
      grid = c(50, 100), coordinate_count = 3L)))
  marginals <- list(artifacts = list(
    sex = list(owner_peer = "site_a", dataset = "cohort", column = "sex",
               levels = c("female", "male")),
    status = list(
      owner_peer = "site_a", dataset = "cohort", column = "status",
      levels = c("censor", "event_a", "event_b"))))
  pairs <- list(sets = list(primary = list(
    owner_peer = "site_a", dataset = "cohort",
    columns = list(
      list(column = "status", levels = c("censor", "event_a", "event_b")),
      list(column = "sex", levels = c("female", "male"))),
    included_pairs = list(c("sex", "status")),
    repeated_record_policy = "consistent_cell_else_exclude_v1",
    missingness_policy = "missing_or_out_of_domain_excluded_v1",
    coordinate_count = 6L, pair_count = 1L)))
  survival <- list(primary = list(
    owner_peer = "site_a", dataset = "cohort", coordinate_count = 5L,
    version = "v1", time = "time", event = "status", entry = "none",
    censor = "censor", causes = c("event_a", "event_b"),
    time_grid = c(5, 10), time_bounds = c(0, 10)))
  list(workload = list(
    coordinate_count = 26,
    families = list(
      admitted_count = count, numeric_moments = numeric,
      numeric_pair_moments = list(artifacts = list()),
      gaussian_models = list(artifacts = list()),
      fixed_numeric_histograms = histograms,
      categorical_marginals = marginals, categorical_pairs = pairs,
      correlation_artifacts = list(), describe_artifacts = list(),
      survival_artifacts = survival)))
}

test_that("client reconstructs the canonical server coordinate order", {
  manifest <- .vector_layout_manifest()
  layout <- .dsvert_dp_capsule_vector_layout(manifest)
  expect_identical(layout$version,
                   .DSVERT_CLIENT_DP_CAPSULE_LAYOUT_VERSION)
  expect_identical(layout$coordinate_count, 26L)
  expect_match(layout$sha256, "^[0-9a-f]{64}$")
  expect_identical(unname(vapply(
    layout$blocks, `[[`, character(1L), "family")), c(
      "admitted_count", "numeric_moments", "numeric_moments",
      "fixed_numeric_histograms", "categorical_marginals",
      "categorical_marginals", "categorical_pairs", "survival_artifacts"))
  expect_identical(unname(vapply(
    layout$blocks, `[[`, integer(1L), "start")),
    c(1L, 2L, 5L, 8L, 11L, 13L, 16L, 22L))
  expect_identical(
    layout$blocks[["numeric_moments::a_key"]]$descriptor$column, "age")
  pair <- layout$blocks[["categorical_pairs::primary::sex::status"]]
  expect_identical(pair$descriptor$left$column, "sex")
  expect_identical(pair$descriptor$right$column, "status")
})

test_that("client rejects ambiguous or mismatched coordinate manifests", {
  bad_count <- .vector_layout_manifest()
  bad_count$workload$coordinate_count <- 25
  expect_error(
    .dsvert_dp_capsule_vector_layout(bad_count), "does not match")

  duplicate <- .vector_layout_manifest()
  duplicate$workload$families$categorical_pairs$sets$primary$columns[[2L]] <-
    duplicate$workload$families$categorical_pairs$sets$primary$columns[[1L]]
  expect_error(
    .dsvert_dp_capsule_vector_layout(duplicate), "pair layout is invalid")

  unknown <- .vector_layout_manifest()
  unknown$workload$families$categorical_pairs$sets$primary$
    included_pairs <- list(c("sex", "missing"))
  expect_error(
    .dsvert_dp_capsule_vector_layout(unknown), "pair layout is invalid")

  wrong_pair_count <- .vector_layout_manifest()
  wrong_pair_count$workload$families$categorical_pairs$sets$primary$
    pair_count <- 2L
  expect_error(
    .dsvert_dp_capsule_vector_layout(wrong_pair_count),
    "pair layout is invalid")

  wrong_coordinate_count <- .vector_layout_manifest()
  wrong_coordinate_count$workload$families$categorical_pairs$sets$primary$
    coordinate_count <- 5L
  expect_error(
    .dsvert_dp_capsule_vector_layout(wrong_coordinate_count),
    "pair layout is invalid")
})

test_that("client never expands an undeclared categorical pair", {
  manifest <- .vector_layout_manifest()
  set <- manifest$workload$families$categorical_pairs$sets$primary
  set$columns[[3L]] <- list(
    column = "site", levels = c("north", "south"))
  manifest$workload$families$categorical_pairs$sets$primary <- set
  layout <- .dsvert_dp_capsule_vector_layout(manifest)

  pair_blocks <- .dsvert_dp_capsule_vector_blocks(
    layout, "categorical_pairs")
  expect_length(pair_blocks, 1L)
  expect_true("categorical_pairs::primary::sex::status" %in%
                names(pair_blocks))
  expect_false(any(grepl("site", names(pair_blocks), fixed = TRUE)))
})

test_that("only final DP coordinates can be extracted from a layout block", {
  layout <- .dsvert_dp_capsule_vector_layout(.vector_layout_manifest())
  release <- list(values = as.numeric(seq_len(26L)), coordinate_count = 26L)
  class(release) <- c("dsvert_joint_dp_vector", "list")
  block <- layout$blocks[["fixed_numeric_histograms::age_hist"]]
  expect_identical(
    .dsvert_dp_capsule_vector_values(release, block), c(8, 9, 10))
  class(release) <- "list"
  expect_error(
    .dsvert_dp_capsule_vector_values(release, block), "cannot be mapped")
})

test_that("two-peer vector accuracy uses the exact convolution tail", {
  epsilon <- 1
  sensitivity <- 2
  radius <- 3L
  support <- -100L:100L
  probability <- exp(-epsilon * abs(support) / sensitivity)
  probability <- probability / sum(probability)
  brute <- sum(outer(
    support, support,
    function(left, right) abs(left + right) > radius) *
      outer(probability, probability))
  certified <- exp(.dsvert_dp_vector_convolution_log_tail(
    radius, epsilon, sensitivity))
  expect_equal(certified, brute, tolerance = 1e-14)

  manifest <- list(workload = list(release_lattice = list(
    output_lattice_bits = 4L, output_lattice_scale = 16,
    natural_l1_sensitivity = 2,
    integer_l1_sensitivity_steps = 32,
    natural_l2_sensitivity = sqrt(2),
    integer_l2_sensitivity_steps = sqrt(2) * 16),
    capsule_mechanism = list(
      mechanism = "discrete-laplace", sensitivity_norm = "l1")))
  release <- list(epsilon = 1, implementation_delta = paste0(
    "1/", "1000000000000000000000000000000"),
    mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM)
  accuracy <- .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = 5L)
  target <- log((0.05 - accuracy$implementation_tv_upper_bound) / 5)
  expect_lte(.dsvert_dp_vector_convolution_log_tail(
    floor(accuracy$radius * 16), 1, 32), target)
  expect_identical(accuracy$coordinate_count, 5L)
  expect_identical(accuracy$additional_privacy_cost,
                   c(epsilon = 0, delta = 0))

  clamped <- .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = 5L, maximum_error = 7)
  expect_identical(clamped$radius, 7)
})

test_that("L2 lattice validation tolerates only binary scale round-off", {
  integer_steps <- 6596.6585480835238
  natural_sensitivity <- 25.76819745345108
  expect_false(identical(integer_steps, natural_sensitivity * 256))
  expect_true(.dsvert_dp_vector_l2_lattice_consistent(
    integer_steps, natural_sensitivity, 256))
  expect_false(.dsvert_dp_vector_l2_lattice_consistent(
    integer_steps + 1e-6, natural_sensitivity, 256))
})

test_that("one-draw accuracy follows the signed dyadic law and separate TV", {
  denominator <- "1267650600228229401496703205376"
  plan <- list(
    version = "dsvert-joint-dp-vector-laplace-plan-v3",
    sampler = .DSVERT_CLIENT_VECTOR_EXACT_SAMPLER,
    stop_bits = 128,
    stop_numerator = "3385866123089935016502817816940606708",
    uniform_bits = 128, binary_geometric_bits = 1,
    bernoulli_thresholds = list("1"), sensitivity_steps = "100",
    total_coordinate_count = 1,
    epsilon_effective_upper_numerator = "1",
    epsilon_effective_upper_denominator = "1",
    one_geometric_tv_numerator = "1",
    one_geometric_tv_denominator = denominator,
    tail_upper_numerator = "1", tail_upper_denominator = denominator,
    rounding_upper_numerator = "1",
    rounding_upper_denominator = denominator,
    implementation_delta_numerator = "1024",
    implementation_delta_denominator = denominator,
    implementation_delta_bound = "8.077935669463161e-28",
    maximum_noise_magnitude = "1048576",
    maximum_chunk_coordinates = 128,
    private_stream_bytes_per_coordinate = 64,
    accounting = paste(
      "global iid discrete Laplace calibrated once to the workload joint",
      "L1 sensitivity; exact binary-geometric coupling"),
    capability_available = TRUE)
  manifest <- list(workload = list(
    coordinate_count = 1,
    release_lattice = list(
      output_lattice_bits = 4L, output_lattice_scale = 16,
      natural_l1_sensitivity = 6.25,
      integer_l1_sensitivity_steps = 100,
      natural_l2_sensitivity = 6.25,
      integer_l2_sensitivity_steps = 100),
    capsule_mechanism = list(
      mechanism = "discrete-laplace", sensitivity_norm = "l1")))
  release <- list(
    epsilon = 1,
    implementation_delta = paste0("1024/", denominator),
    backend = .DSVERT_CLIENT_VECTOR_EXACT_BACKEND,
    mechanism = .DSVERT_CLIENT_VECTOR_EXACT_RELEASE_MECHANISM,
    mechanism_plan = plan, plan_sha256 = .dsvert_vector_hash(plan))
  accuracy <- .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = 1L, confidence = 0.95)

  q <- .dsvert_dp_vector_dyadic_tail_context(plan)$q
  p <- 1 - q
  tv <- 2 / as.numeric(denominator)
  target <- 0.05 - tv
  oracle <- 0L
  while (2 * p^(oracle + 1) / (1 + p) > target) {
    oracle <- oracle + 1L
  }
  certified_steps <- accuracy$radius * 16
  expect_gte(certified_steps, oracle)
  expect_lte(exp(.dsvert_dp_vector_plan_log_tail_upper(
    certified_steps, .dsvert_dp_vector_dyadic_tail_context(plan))),
    target)
  expect_equal(
    accuracy$sampler_tv_upper_bound, tv,
    tolerance = 64 * .Machine$double.eps * tv)
  expect_gt(accuracy$implementation_delta_bound,
            accuracy$sampler_tv_upper_bound)
  expect_true(accuracy$accuracy_plan_certified)

  large <- plan
  large$stop_numerator <- "309485009821345068724781056"
  large$sensitivity_steps <- "999999999999"
  large_release <- release
  large_release$mechanism_plan <- large
  large_release$plan_sha256 <- .dsvert_vector_hash(large)
  large_manifest <- manifest
  large_manifest$workload$release_lattice$output_lattice_bits <- 1L
  large_manifest$workload$release_lattice$output_lattice_scale <- 2
  large_manifest$workload$release_lattice$natural_l1_sensitivity <-
    999999999999 / 2
  large_manifest$workload$release_lattice$integer_l1_sensitivity_steps <-
    999999999999
  large_manifest$workload$release_lattice$natural_l2_sensitivity <-
    999999999999 / 2
  large_manifest$workload$release_lattice$integer_l2_sensitivity_steps <-
    999999999999
  large_accuracy <- .dsvert_dp_vector_accuracy_radius(
    large_release, large_manifest, coordinate_count = 1L,
    confidence = 0.95)
  large_q <- 2^-40
  large_p <- 1 - large_q
  large_oracle <- ceiling(
    log(target * (1 + large_p) / 2) / log1p(-large_q)) - 1
  expect_gte(large_accuracy$radius * 2, large_oracle)
})

.vector_double_compare_rational <- function(value, numerator, denominator) {
  encoded <- sprintf("%a", value)
  matched <- regmatches(encoded, regexec(
    "^0x([0-9a-f]+)(?:\\.([0-9a-f]+))?p([+-][0-9]+)$",
    encoded, perl = TRUE))[[1L]]
  if (length(matched) != 4L) return(NA_integer_)
  fractional <- matched[[3L]]
  mantissa_hex <- paste0(matched[[2L]], fractional)
  if (nchar(mantissa_hex) %% 2L == 1L) {
    mantissa_hex <- paste0("0", mantissa_hex)
  }
  mantissa <- openssl::bignum(
    mantissa_hex, hex = TRUE)
  exponent <- as.integer(matched[[4L]]) - 4L * nchar(fractional)
  if (exponent >= 0L) {
    value_numerator <- mantissa * openssl::bignum(2)^exponent
    value_denominator <- openssl::bignum(1)
  } else {
    value_numerator <- mantissa
    value_denominator <- openssl::bignum(2)^(-exponent)
  }
  left <- value_numerator * openssl::bignum(denominator)
  right <- openssl::bignum(numerator) * value_denominator
  if (left < right) -1L else if (left > right) 1L else 0L
}

.vector_double_at_least_rational <- function(
    value, numerator, denominator) {
  .vector_double_compare_rational(value, numerator, denominator) >= 0L
}

.vector_double_at_most_rational <- function(
    value, numerator, denominator) {
  .vector_double_compare_rational(value, numerator, denominator) <= 0L
}

test_that("signed dyadic conversion encloses the exact rational", {
  denominator <- as.character(openssl::bignum("2")^128L)
  numerators <- c(
    "1326635224458652993228324758901174720",
    "170141183460469231731687303715884105728")
  withr::local_seed(20260809L)
  random_hex <- replicate(1000L, paste0(sprintf(
    "%02x", sample.int(256L, 16L, replace = TRUE) - 1L),
    collapse = ""))
  random <- vapply(random_hex, function(value) {
    as.character(openssl::bignum(value, hex = TRUE))
  }, character(1L))
  random[random == "0"] <- "1"
  numerators <- c(numerators, random)
  for (numerator in numerators) {
    interval <- .dsvert_dp_vector_dyadic_fraction_interval(numerator, 128L)
    expect_true(.vector_double_at_most_rational(
      interval$q_lower, numerator, denominator), info = numerator)
    expect_true(.vector_double_at_least_rational(
      interval$q_upper, numerator, denominator), info = numerator)
    expect_gte(interval$q, interval$q_lower)
    expect_lte(interval$q, interval$q_upper)
  }

  half <- .dsvert_dp_vector_dyadic_fraction_interval(
    "170141183460469231731687303715884105728", 128L)
  expect_identical(half$q, 0.5)
  subnormal <- .dsvert_dp_vector_dyadic_fraction_interval("1", 1023L)
  subnormal_denominator <- as.character(openssl::bignum("2")^1023L)
  expect_lt(subnormal$q, .Machine$double.xmin)
  expect_true(.vector_double_at_most_rational(
    subnormal$q_lower, "1", subnormal_denominator))
  expect_true(.vector_double_at_least_rational(
    subnormal$q_upper, "1", subnormal_denominator))
})

test_that("exact-GC variance is an outward bound on the signed rational law", {
  numerator <- openssl::bignum(
    "1326635224458652993228324758901174720")
  denominator <- openssl::bignum("2")^128L
  plan <- list(
    stop_bits = 128,
    stop_numerator = as.character(numerator),
    total_coordinate_count = 4,
    one_geometric_tv_numerator = "1",
    one_geometric_tv_denominator =
      "1267650600228229401496703205376")
  release <- list(
    capsule_mechanism = list(mechanism = "discrete-laplace"),
    mechanism = .DSVERT_CLIENT_VECTOR_EXACT_RELEASE_MECHANISM,
    implementation = .DSVERT_CLIENT_VECTOR_EXACT_BACKEND,
    output_lattice_scale = 256,
    l1_sensitivity = 1,
    epsilon = 1,
    mechanism_plan = plan)
  noise <- .dsvert_dp_chisq_noise_contract(release)
  variance_numerator <- openssl::bignum("2") *
    (denominator - numerator) * denominator
  variance_denominator <- numerator * numerator *
    openssl::bignum(as.character(256^2))
  expect_true(.vector_double_at_least_rational(
    noise$variance,
    as.character(variance_numerator),
    as.character(variance_denominator)))
  expect_true(noise$variance_is_upper_bound)
})

test_that("implementation-delta conversion is bounded and outward", {
  numerator <- paste0(
    "211558757151573125914622082483080874355185947817867791267922291",
    "1252012249457309860615360823402031939982373274520985549")
  denominator <- paste0(
    "258964424476258387776423832887791397529165256639620593669822365",
    "235784270237541625590084668100605512305241594195201437163147851",
    "6170478321664000000000000000000")
  certificate <- paste0(numerator, "/", denominator)
  expect_gt(nchar(denominator, type = "bytes"), 128L)

  bound <- .dsvert_dp_vector_fraction(certificate)
  expect_true(.vector_double_at_least_rational(
    bound, numerator, denominator))
  expect_identical(.dsvert_dp_vector_fraction("0/1"), 0)

  near_numerator <- strrep("9", 400L)
  near_denominator <- paste0("1", strrep("0", 400L))
  near_bound <- .dsvert_dp_vector_fraction(paste0(
    near_numerator, "/", near_denominator))
  expect_identical(near_bound, 1)
  expect_true(.vector_double_at_least_rational(
    near_bound, near_numerator, near_denominator))

  subnormal_denominator <- paste0("1", strrep("0", 400L))
  subnormal_bound <- .dsvert_dp_vector_fraction(paste0(
    "1/", subnormal_denominator))
  expect_gt(subnormal_bound, 0)
  expect_true(.vector_double_at_least_rational(
    subnormal_bound, "1", subnormal_denominator))

  malformed <- c("1/0", "1/1", "2/1", "01/2", "1//2", "1/2\n")
  for (value in malformed) {
    expect_error(.dsvert_dp_vector_fraction(value),
                 "Invalid vector implementation-delta certificate",
                 fixed = TRUE)
  }
  oversized <- paste0(
    strrep("1", .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES), "/2")
  expect_error(.dsvert_dp_vector_fraction(oversized),
               "Invalid vector implementation-delta certificate",
               fixed = TRUE)

  manifest <- list(workload = list(release_lattice = list(
    output_lattice_bits = 4L, output_lattice_scale = 16,
    natural_l1_sensitivity = 2,
    integer_l1_sensitivity_steps = 32,
    natural_l2_sensitivity = sqrt(2),
    integer_l2_sensitivity_steps = sqrt(2) * 16),
    capsule_mechanism = list(
      mechanism = "discrete-laplace", sensitivity_norm = "l1")))
  release <- list(
    epsilon = 1, implementation_delta = certificate,
    mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM)
  accuracy <- .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = 1L)
  expect_true(.vector_double_at_least_rational(
    accuracy$implementation_delta_bound, numerator, denominator))

  noise <- .dsvert_dp_chisq_noise_contract(list(
    output_lattice_scale = 16, l1_sensitivity = 2, epsilon = 1,
    implementation_delta = certificate))
  twice_numerator <- as.character(openssl::bignum(numerator) * 2)
  expect_true(.vector_double_at_least_rational(
    noise$vector_tv_upper_bound, twice_numerator, denominator))
})
