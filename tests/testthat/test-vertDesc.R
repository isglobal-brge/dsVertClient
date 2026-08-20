# Tests for .dsvert_interp_quantile (pure client-side helper for ds.vertDesc).

library(testthat)

.desc_dp_fixture <- function() {
  variables <- c("age", "bmi")
  descriptives <- data.frame(
    variable = variables,
    status = c("ok", "dp_point_available_count_not_certified_positive"),
    n_dp = c(101, 99),
    n_simultaneous_95_lower = c(90, 0),
    mean = c(52.5, 24.75),
    variance = c(100, 16),
    sd = c(10, 4),
    lower_bound = c(0, 10),
    upper_bound = c(100, 50),
    invalid_dp = c(3, 5),
    stringsAsFactors = FALSE)
  quantiles <- do.call(rbind, lapply(seq_along(variables), function(index) {
    data.frame(
      variable = variables[[index]],
      probability = c(0.25, 0.5, 0.75),
      estimate = c(40, 52, 65) + index,
      dp_grid_lower = c(35, 48, 60) + index,
      dp_grid_upper = c(45, 57, 70) + index,
      status = "ok",
      stringsAsFactors = FALSE)
  }))
  result <- list(
    released = TRUE,
    analysis_id = "baseline_v1",
    analysis_version = "1",
    variables = variables,
    grid_lengths = c(4L, 5L),
    descriptives = descriptives,
    quantiles = quantiles,
    server = "site_a",
    capsule_id = paste(rep("d", 64L), collapse = ""),
    final_vector_root = paste(rep("a", 64L), collapse = ""),
    coordinate_order_sha256 = paste(rep("b", 64L), collapse = ""),
    privacy_epoch = 1L,
    noise_key_id = paste(rep("c", 64L), collapse = ""),
    history_gate = TRUE,
    request_limit = FALSE,
    operation_limit = TRUE,
    epsilon = 1,
    delta = 1e-6,
    implementation_delta = 0,
    adjacency = "add_remove_patient",
    mechanism = "complete_discrete_laplace_two_peer",
    sampler = "hkdf_chacha20",
    sticky_noise = "one immutable capsule vector; unlimited replay",
    uncertainty_scope = "DP mechanism noise only; sampling uncertainty excluded",
    histogram_semantics = "fixed signed bins plus invalid bin",
    unit_collapse = "one bounded value per admitted unit",
    count_definition = "DP-noisy effective units",
    invalid_unit_rule = "invalid units enter the fixed invalid bin",
    quantization = "signed fixed lattice",
    postprocessing = "fixed public clamp",
    quantile_band_confidence = 0.95,
    quantile_band_scope = "mechanism noise plus grid; no sampling uncertainty",
    moment_region_confidence = 0.95,
    moment_region_method = "simultaneous coordinate-box propagation",
    moment_region_scope = "mechanism and quantisation; no sampling uncertainty",
    statistical_inference = "DP points; no sampling intervals or p-values")
  class(result) <- c("ds.vertDPDescribe", "list")
  result
}

.desc_synopsis_fixture <- function() {
  result <- .desc_dp_fixture()
  anchors <- list(
    artifact_key = strrep("e", 64L),
    execution_id = strrep("f", 64L),
    contract_sha256 = strrep("1", 64L),
    attempt_sha256 = strrep("2", 64L),
    source_contract_sha256 = strrep("3", 64L),
    result_set_sha256 = strrep("4", 64L),
    final_vector_root = result$final_vector_root)
  result[c(
    "capsule_id", "privacy_epoch", "noise_key_id", "history_gate",
    "request_limit", "operation_limit")] <- NULL
  result[names(anchors)[-length(anchors)]] <- anchors[-length(anchors)]
  result$release_provenance <- c(list(
    version = "dsvert-stateless-synopsis-public-provenance-v1"), anchors)
  result$privacy <- list(
    version = "dsvert-per-synopsis-dp-v1",
    per_artifact_epsilon = result$epsilon,
    per_artifact_delta = result$delta,
    sticky_noise = TRUE,
    public_openings = 1L,
    distinct_artifacts_compose = TRUE,
    finite_global_composition_claim = FALSE)
  result$composition_rule <-
    "one_sticky_release_per_canonical_signed_artifact"
  result$security_claim <- list(
    version = "dsvert-synopsis-security-claim-v1",
    finite_global_composition_claim = FALSE)
  result
}

test_that("ds.vertDesc is a one-release compatibility adapter", {
  fixture <- .desc_dp_fixture()
  calls <- new.env(parent = emptyenv())
  calls$count <- 0L
  calls$args <- NULL
  testthat::local_mocked_bindings(
    ds.vertDPDescribe = function(data_name, analysis_id, probs,
                                 datasources = NULL) {
      calls$count <- calls$count + 1L
      calls$args <- list(data_name = data_name, analysis_id = analysis_id,
                         probs = probs, datasources = datasources)
      fixture$quantiles <- fixture$quantiles[
        fixture$quantiles$probability %in% probs, , drop = FALSE]
      fixture
    },
    .package = "dsVertClient")

  conns <- list(site_a = structure(list(), class = "mock"))
  result <- ds.vertDesc(
    "cohort", variables = "bmi", probs = c(0.5, 0.25),
    verbose = FALSE, datasources = conns, analysis_id = "baseline_v1")

  expect_s3_class(result, "ds.vertDesc")
  expect_identical(calls$count, 1L)
  expect_identical(calls$args$analysis_id, "baseline_v1")
  expect_identical(calls$args$probs, c(0.25, 0.5))
  expect_identical(result$server, "site_a")
  expect_identical(result$variable, "bmi")
  expect_identical(result$n, 99)
  expect_identical(result$n_na, 5)
  expect_true(is.na(result$min) && is.na(result$max))
  expect_identical(result$range_low, 10)
  expect_identical(result$range_high, 50)
  expect_identical(result$range_method, "custodian_signed_fixed_grid")
  expect_identical(result$histogram_buckets, 5L)
  expect_identical(result$quantile_status, "ran_dp_fixed_grid")
  expect_identical(result$q25, 42)
  expect_identical(result$q50, 54)
  expect_identical(attr(result, "signed_histogram_buckets"), c(bmi = 5L))
  expect_true(attr(result, "dp_release")$formal_dp)
  expect_identical(
    attr(result, "dp_release")$capsule_id, strrep("d", 64L))
  expect_match(attr(result, "dp_release")$uncertainty_scope,
               "sampling uncertainty excluded")
  expect_identical(attr(result, "dp_release")$selected_variables, "bmi")
  expect_match(attr(result, "dp_release")$statistical_inference,
               "no sampling intervals")
  expect_match(attr(result, "compatibility_semantics")$sd,
               "not the usual sample standard deviation")
  expect_equal(nrow(attr(result, "dp_quantile_bands")), 2L)
})

test_that("ds.vertDesc propagates synopsis release provenance by branch", {
  fixture <- .desc_synopsis_fixture()
  testthat::local_mocked_bindings(
    ds.vertDPDescribe = function(...) fixture,
    .package = "dsVertClient")

  result <- ds.vertDesc(
    "cohort", variables = "bmi", probs = 0.5,
    verbose = FALSE, analysis_id = "baseline_v1")
  release <- attr(result, "dp_release")
  expect_identical(release$artifact_key, fixture$artifact_key)
  expect_identical(release$execution_id, fixture$execution_id)
  expect_identical(release$contract_sha256, fixture$contract_sha256)
  expect_identical(release$attempt_sha256, fixture$attempt_sha256)
  expect_identical(
    release$source_contract_sha256, fixture$source_contract_sha256)
  expect_identical(release$result_set_sha256, fixture$result_set_sha256)
  expect_identical(release$release_provenance, fixture$release_provenance)
  expect_identical(release$privacy, fixture$privacy)
  expect_identical(release$security_claim, fixture$security_claim)
  expect_false(any(c(
    "capsule_id", "privacy_epoch", "noise_key_id", "history_gate",
    "request_limit", "operation_limit") %in% names(release)))
})

test_that("ds.vertDesc rejects mixed describe provenance", {
  synopsis <- .desc_synopsis_fixture()
  synopsis$capsule_id <- strrep("5", 64L)
  testthat::local_mocked_bindings(
    ds.vertDPDescribe = function(...) synopsis,
    .package = "dsVertClient")
  expect_error(
    ds.vertDesc("cohort", analysis_id = "baseline_v1", verbose = FALSE),
    "intact released")
})

test_that("legacy knobs cannot alter the signed describe workload", {
  fixture <- .desc_dp_fixture()
  calls <- new.env(parent = emptyenv())
  calls$count <- 0L
  testthat::local_mocked_bindings(
    ds.vertDPDescribe = function(...) {
      calls$count <- calls$count + 1L
      fixture
    },
    .package = "dsVertClient")

  expect_error(
    ds.vertDesc("cohort", analysis_id = "baseline_v1",
                exact_extrema = TRUE, verbose = FALSE),
    "not disclosure-safe")
  expect_error(
    ds.vertDesc("cohort", analysis_id = "baseline_v1",
                range_sd = 4, verbose = FALSE),
    "cannot redefine")
  expect_error(
    ds.vertDesc("cohort", analysis_id = "baseline_v1",
                open_ended = TRUE, verbose = FALSE),
    "cannot redefine")
  expect_identical(calls$count, 0L)

  expect_error(
    ds.vertDesc("cohort", variables = "bmi", n_buckets = 4L,
                analysis_id = "baseline_v1", verbose = FALSE),
    "does not match the signed fixed grid")
  expect_identical(calls$count, 1L)
  matched <- ds.vertDesc(
    "cohort", variables = "bmi", n_buckets = 5L,
    analysis_id = "baseline_v1", verbose = FALSE)
  expect_identical(matched$histogram_buckets, 5L)
  expect_identical(calls$count, 2L)
})

test_that("analysis and variable selection fail honestly", {
  fixture <- .desc_dp_fixture()
  calls <- new.env(parent = emptyenv())
  calls$count <- 0L
  testthat::local_mocked_bindings(
    ds.vertDPDescribe = function(...) {
      calls$count <- calls$count + 1L
      fixture
    },
    .package = "dsVertClient")

  expect_error(ds.vertDesc("cohort", verbose = FALSE),
               "analysis_id is required")
  expect_identical(calls$count, 0L)
  expect_error(
    ds.vertDesc("cohort", variables = "glucose",
                analysis_id = "baseline_v1", verbose = FALSE),
    "not included in signed describe artifact")
  expect_error(
    ds.vertDesc("cohort", variables = list(site_b = "age"),
                analysis_id = "baseline_v1", verbose = FALSE),
    "owned by server 'site_a'")
  expect_identical(calls$count, 2L)
})

test_that("ds.vertDesc call graph cannot reach exact adaptive endpoints", {
  namespace <- asNamespace("dsVertClient")
  reachable <- character()
  queue <- "ds.vertDesc"
  while (length(queue)) {
    name <- queue[[1L]]
    queue <- queue[-1L]
    if (name %in% reachable ||
        !exists(name, namespace, mode = "function", inherits = FALSE)) {
      next
    }
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
    "dsvertColNamesDS", "dsvertLocalMomentsDS", "dsvertHistogramDS",
    "dsvertDPDescribeDS")
  bodies <- paste(vapply(reachable, function(name) {
    paste(deparse(body(get(name, namespace, inherits = FALSE))),
          collapse = "\n")
  }, character(1L)), collapse = "\n")
  expect_contains(reachable, "ds.vertDPDescribe")
  expect_length(intersect(reachable, forbidden), 0L)
  for (endpoint in forbidden) {
    expect_false(grepl(endpoint, bodies, fixed = TRUE), info = endpoint)
  }
})

test_that("interp_quantile matches known closed-form for uniform buckets", {
  # 100 observations perfectly uniform across [0, 1] in 10 buckets
  edges <- seq(0, 1, length.out = 11L)
  counts <- rep(10L, 10L)
  q <- dsVertClient:::.dsvert_interp_quantile(
    edges, counts, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  expect_equal(q, c(0.1, 0.25, 0.5, 0.75, 0.9), tolerance = 1e-12)
})

test_that("interp_quantile is within bucket width of exact quantile on normal", {
  set.seed(7)
  x <- rnorm(5000)
  edges <- seq(min(x) - 1e-9, max(x) + 1e-9, length.out = 201L)
  bucket <- findInterval(x, edges, rightmost.closed = TRUE)
  counts <- tabulate(bucket, nbins = 200L)

  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  bucket_width <- (edges[length(edges)] - edges[1]) / 200

  q_approx <- dsVertClient:::.dsvert_interp_quantile(edges, counts, probs)
  q_exact <- stats::quantile(x, probs = probs, names = FALSE)

  # With 200 buckets on a 5000-point normal sample, interpolation error
  # should be on the order of one bucket width or less.
  for (k in seq_along(probs)) {
    expect_true(abs(q_approx[k] - q_exact[k]) < 1.5 * bucket_width,
      info = sprintf("p=%.2f: approx=%.4f exact=%.4f diff=%.4f bucket_w=%.4f",
                     probs[k], q_approx[k], q_exact[k],
                     q_approx[k] - q_exact[k], bucket_width))
  }
})

test_that("interp_quantile handles under/overflow cleanly", {
  edges <- c(0, 0.5, 1.0)
  counts <- c(10L, 10L)
  # 20 observations below edges[1] and 10 above edges[3]
  # total = sum(counts) + below + above = 20 + 20 + 10 = 50
  q <- dsVertClient:::.dsvert_interp_quantile(
    edges, counts, probs = c(0.3, 0.5, 0.75, 0.9),
    below = 20L, above = 10L)
  # p=0.3 -> target 15 (<= below=20)                    -> edges[1] = 0
  # p=0.5 -> target 25, in bucket 1 [cum 20..30)         -> 0 + 0.5*0.5 = 0.25
  # p=0.75 -> target 37.5, in bucket 2 [cum 30..40)       -> 0.5 + 0.75*0.5 = 0.875
  # p=0.9 -> target 45 (>= below + sum(counts) = 40)      -> edges[K+1] = 1.0
  expect_equal(q[1], 0)
  expect_equal(q[2], 0.25, tolerance = 1e-12)
  expect_equal(q[3], 0.875, tolerance = 1e-12)
  expect_equal(q[4], 1.0, tolerance = 1e-12)
})

test_that("interp_quantile returns NA on empty input", {
  q <- dsVertClient:::.dsvert_interp_quantile(
    edges = c(0, 1), counts = c(0L), probs = c(0.5),
    below = 0, above = 0)
  expect_true(all(is.na(q)))
})

test_that("interp_quantile handles zero-count buckets gracefully", {
  # Buckets: [0,1) empty, [1,2) 10 obs, [2,3) empty
  edges <- c(0, 1, 2, 3)
  counts <- c(0L, 10L, 0L)
  q <- dsVertClient:::.dsvert_interp_quantile(
    edges, counts, probs = c(0.5))
  # 50% of 10 = 5 → midpoint of bucket 2
  expect_equal(q, 1.5, tolerance = 1e-12)
})
