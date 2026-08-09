.dp_quantile_release <- function(
    grid = c(1, 2, 3, 4), counts = c(1, 2, 3, 4), lower = 0,
    radius = 1) {
  stopifnot(length(grid) == length(counts), length(grid) >= 1L)
  scale <- 256
  statistics <- c(
    sum(counts), round(0.5 * sum(counts) * scale),
    round(0.35 * sum(counts) * scale), counts, 0)
  result <- list(
    released = TRUE,
    analysis_id = "primary",
    analysis_version = "v1",
    variables = "x",
    variable_count = 1L,
    lower_bounds = lower,
    upper_bounds = tail(grid, 1L),
    grid_lengths = as.integer(length(grid)),
    grid_values = grid,
    histogram_semantics =
      "(previous_endpoint,current_endpoint] plus fixed invalid bin",
    unit_collapse =
      "mean_of_finite_rows_after_public_bound_clipping",
    count_definition =
      "DP-noisy effective units with at least one finite bounded value",
    invalid_unit_rule = "invalid_patient_ids_rejected_by_admission",
    statistics = statistics,
    coordinate_count = as.integer(length(statistics)),
    coordinate_layout = paste(
      "referenced capsule blocks per variable[count,qsum,qsumsq,",
      "histogram[grid_bins+invalid]]"),
    numeric_grid_bits = 8L,
    numeric_grid_scale = scale,
    quantization = "round(z*scale) and round(z^2*scale) independently",
    max_abs_normalized_quantization_per_unit = 0.5 / scale,
    allocation_names = c("count", "sum", "sumsq", "histogram"),
    allocation_weights = rep(0.25, 4L),
    mechanism = .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM,
    implementation = "two pinned peers; Ring128 exact signed finalizer",
    sampler = .DSVERT_CLIENT_VECTOR_SAMPLER,
    randomness =
      "independent HMAC/HKDF/ChaCha20 streams at both peers",
    postprocessing = "fixed public per-coordinate clamp",
    clipped_coordinates = NA_integer_,
    clipping_observable = FALSE,
    accuracy_95_abs_by_variable_family =
      c(radius / 2, radius, radius, radius / 2),
    accuracy_simultaneous_95_abs_by_variable_family =
      c(radius, 2 * radius, 2 * radius, radius),
    accuracy_simultaneous_confidence = 0.95,
    accuracy_simultaneous_method =
      "exact convolution tail with union bound",
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    privacy_epoch = 1,
    noise_key_id = "noise-key-a",
    sticky_noise = "one immutable capsule vector; unlimited replay",
    epsilon = 1,
    delta = 2^-100,
    implementation_delta =
      "1/1267650600228229401496703205376",
    adjacency = "add_remove_patient",
    composition_partitions = 1L,
    capsule_id = strrep("a", 64L),
    final_vector_root = strrep("b", 64L),
    coordinate_order_sha256 = strrep("c", 64L),
    history_gate = TRUE, request_limit = FALSE, operation_limit = TRUE,
    server = "site_a")
  result <- .dsvert_dp_describe_postprocess(result, probs = 0.5)
  class(result) <- c("ds.vertDPDescribe", "list")
  result
}

.dp_quantile_candidate_index <- function(counts, probability) {
  total <- sum(counts)
  if (!total) return(NA_integer_)
  if (probability == 0) return(which(counts > 0)[[1L]])
  which(cumsum(counts) >= probability * total)[[1L]]
}

test_that("binned quantiles match a synthetic projected-histogram oracle", {
  release <- .dp_quantile_release()
  result <- ds.vertDPQuantile(
    release, probs = c(1, 0.5, 0.3, 0.1, 0, 0.5))

  expect_s3_class(result, "ds.vertDPQuantile")
  expect_identical(result$probability, c(0, 0.1, 0.3, 0.5, 1))
  expect_identical(result$bin_index, c(1L, 1L, 2L, 3L, 4L))
  expect_identical(
    result$bin_id,
    paste0("x::bin_", c(1L, 1L, 2L, 3L, 4L)))
  expect_identical(result$bin_lower, c(0, 0, 1, 2, 3))
  expect_identical(result$bin_upper, c(1, 1, 2, 3, 4))
  expect_true(result$bin_left_closed[[1L]])
  expect_false(any(result$bin_left_closed[-c(1L, 2L)]))
  expect_true(all(result$bin_right_closed))
  expect_false("estimate" %in% names(result))
  expect_match(attr(result, "estimate_scope"), "no exact sample")
  expect_match(attr(result, "estimate_scope"), "no within-bin")
})

test_that("all released variables share one deterministic probability order", {
  one <- .dp_quantile_release()
  two <- unclass(one)
  block <- one$statistics
  two$variables <- c("x", "y")
  two$variable_count <- 2L
  two$lower_bounds <- c(0, -4)
  two$upper_bounds <- c(4, 0)
  two$grid_lengths <- c(4L, 4L)
  two$grid_values <- c(1, 2, 3, 4, -3, -2, -1, 0)
  two$statistics <- c(block, block)
  two$coordinate_count <- 16L
  two$accuracy_95_abs_by_variable_family <-
    rep(one$accuracy_95_abs_by_variable_family, 2L)
  two$accuracy_simultaneous_95_abs_by_variable_family <-
    rep(one$accuracy_simultaneous_95_abs_by_variable_family, 2L)
  two <- .dsvert_dp_describe_postprocess(two, probs = 0.5)
  class(two) <- c("ds.vertDPDescribe", "list")

  result <- ds.vertDPQuantile(two, c(0.75, 0.25))
  expect_identical(result$variable, rep(c("x", "y"), each = 2L))
  expect_identical(result$probability, rep(c(0.25, 0.75), 2L))
  expect_identical(result$bin_index, rep(c(2L, 4L), 2L))
  expect_identical(result$bin_lower, c(1, 3, -3, -1))
  expect_identical(result$bin_upper, c(2, 4, -2, 0))
})

test_that("simultaneous quantile regions cover every small count box", {
  release <- .dp_quantile_release(
    grid = c(1, 2, 3), counts = c(2, 1, 3), radius = 1)
  probs <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  result <- ds.vertDPQuantile(release, probs)
  lower <- pmax(0, release$histograms$x$counts - 1)
  upper <- release$histograms$x$counts + 1
  candidates <- expand.grid(
    Map(seq.int, lower, upper), KEEP.OUT.ATTRS = FALSE)
  checked <- 0L
  for (row in seq_len(nrow(candidates))) {
    counts <- as.numeric(candidates[row, ])
    for (index in seq_along(probs)) {
      bin <- .dp_quantile_candidate_index(counts, probs[[index]])
      if (is.na(bin)) next
      endpoint <- release$histograms$x$grid[[bin]]
      expect_gte(
        endpoint, result$mechanism_grid_lower_95[[index]] - 1e-12)
      expect_lte(
        endpoint, result$mechanism_grid_upper_95[[index]] + 1e-12)
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 0L)
  expect_true(all(diff(result$bin_index) >= 0))
  expect_true(all(diff(result$mechanism_grid_lower_95) >= 0))
  expect_true(all(diff(result$mechanism_grid_upper_95) >= 0))
})

test_that("band-index oracle is exhaustive over endpoint and interior probs", {
  probs <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  observed <- expand.grid(rep(list(0:3), 3L), KEEP.OUT.ATTRS = FALSE)
  covered <- TRUE
  checked <- 0L
  for (observed_index in seq_len(nrow(observed))) {
    counts_dp <- as.numeric(observed[observed_index, ])
    for (radius in 0:3) {
      bands <- lapply(probs, function(probability) {
        .dsvert_dp_quantile_band_indices(
          counts_dp, probability, radius)
      })
      candidates <- expand.grid(
        Map(seq.int, pmax(0, counts_dp - radius), counts_dp + radius),
        KEEP.OUT.ATTRS = FALSE)
      for (candidate_index in seq_len(nrow(candidates))) {
        counts <- as.numeric(candidates[candidate_index, ])
        if (!sum(counts)) next
        for (probability_index in seq_along(probs)) {
          bin <- .dp_quantile_candidate_index(
            counts, probs[[probability_index]])
          band <- bands[[probability_index]]
          covered <- covered &&
            bin >= band[["lower"]] && bin <= band[["upper"]]
          checked <- checked + 1L
        }
      }
    }
  }
  expect_identical(checked, 117992L)
  expect_true(covered)
})

test_that("median is a zero-cost release-only quantile wrapper", {
  release <- .dp_quantile_release()
  median <- ds.vertDPMedian(release)
  quantile <- ds.vertDPQuantile(release, 0.5)

  expect_s3_class(median, "ds.vertDPMedian")
  expect_s3_class(median, "ds.vertDPQuantile")
  expect_identical(median$probability, 0.5)
  expect_identical(as.data.frame(median), as.data.frame(quantile))
  expect_identical(
    attr(median, "additional_privacy_cost"),
    c(epsilon = 0, delta = 0))
  expect_identical(attr(median, "additional_server_calls"), 0L)
  expect_true(attr(median, "postprocessing_only"))
  provenance <- attr(median, "source_provenance")
  expect_identical(provenance$analysis_id, release$analysis_id)
  expect_identical(provenance$capsule_id, release$capsule_id)
  expect_identical(provenance$final_vector_root, release$final_vector_root)
  expect_identical(provenance$mechanism, release$mechanism)
  expect_identical(provenance$epsilon, release$epsilon)
  expect_identical(provenance$delta, release$delta)
  expect_identical(provenance$adjacency, release$adjacency)
})

test_that("validated Laplace and Gaussian capsule profiles are explicit", {
  laplace <- .dp_quantile_release()
  expect_no_error(ds.vertDPQuantile(laplace, 0.5))

  gaussian <- laplace
  gaussian$mechanism <- .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM
  gaussian$sampler <- .DSVERT_CLIENT_VECTOR_GAUSSIAN_SAMPLER
  result <- ds.vertDPQuantile(gaussian, 0.5)
  expect_identical(
    attr(result, "source_provenance")$mechanism,
    .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)
})

test_that("empty and one-bin projected histograms are explicit", {
  empty <- .dp_quantile_release(counts = rep(0, 4L))
  result <- ds.vertDPQuantile(empty, c(0, 0.5, 1))
  expect_true(all(result$status == "dp_projected_histogram_empty"))
  expect_true(all(is.na(result$bin_id)))
  expect_true(all(is.na(result$bin_index)))
  expect_true(all(is.na(result$bin_lower)))
  expect_true(all(is.na(result$bin_upper)))
  expect_identical(result$projected_histogram_mass_dp, rep(0, 3L))
  expect_identical(result$mechanism_grid_lower_95, rep(0, 3L))
  expect_identical(result$mechanism_grid_upper_95, rep(4, 3L))

  one_bin <- .dp_quantile_release(grid = 4, counts = 7)
  degenerate <- ds.vertDPQuantile(one_bin, c(0, 0.5, 1))
  expect_identical(degenerate$bin_index, rep(1L, 3L))
  expect_identical(degenerate$bin_lower, rep(0, 3L))
  expect_identical(degenerate$bin_upper, rep(4, 3L))
  expect_true(all(degenerate$bin_left_closed))
})

test_that("common-lattice fractional DP mass is not silently renormalized", {
  release <- .dp_quantile_release(
    counts = c(0.25, 0.5, 1.25, 0), radius = 0.5)
  release$statistics[[1L]] <- 7
  release <- .dsvert_dp_describe_postprocess(release, probs = 0.5)
  class(release) <- c("ds.vertDPDescribe", "list")
  result <- ds.vertDPQuantile(release, 0.5)
  expect_identical(result$projected_histogram_mass_dp, 2)
  expect_identical(release$descriptives$n_dp, 7)
  expect_identical(result$bin_index, 3L)
  expect_identical(result$bin_lower, 2)
  expect_identical(result$bin_upper, 3)
  expect_match(attr(result, "histogram_mass_scope"), "invalid bin")
  expect_match(attr(result, "histogram_mass_scope"), "not silently forced")
})

test_that("quantile post-processing rejects malformed provenance and fields", {
  release <- .dp_quantile_release()
  malformed <- list(
    unclass(release),
    within(release, status <- "unchecked"),
    within(release, released <- FALSE),
    within(release, adjacency <- NULL),
    within(release, final_vector_root <- "not-a-root"),
    within(release, mechanism <- "unvalidated-noise"),
    within(release, sampler <- "base-r-rng"),
    within(release, history_gate <- FALSE),
    within(release, upper_bounds <- 5),
    within(release, histograms$x$counts[[1L]] <- 99),
    within(release, {
      statistics[[4L]] <- 0.3
      histograms$x$counts[[1L]] <- 0.3
    }),
    within(release, histograms$x$grid[[1L]] <- 0.5),
    within(release, descriptives$n_dp[[1L]] <- 999),
    within(release,
           histograms$x$cell_noise_radius_simultaneous_95 <- 999))
  for (value in malformed) {
    expect_error(
      ds.vertDPQuantile(value, 0.5),
      "intact released ds.vertDPDescribe capsule object")
  }
  for (probs in list(numeric(), NA_real_, NaN, Inf, -0.1, 1.1, "0.5")) {
    expect_error(ds.vertDPQuantile(release, probs), "probabilities in")
  }
})

test_that("quantile post-processing leaves RNG and DSI untouched", {
  release <- .dp_quantile_release()
  set.seed(20260802)
  before <- .Random.seed
  testthat::local_mocked_bindings(
    .dsvert_dp_capsule_vector_run = function(...) {
      stop("unexpected DSI capsule call", call. = FALSE)
    },
    .dsvert_dp_datasources = function(...) {
      stop("unexpected DSI connection lookup", call. = FALSE)
    },
    .package = "dsVertClient")

  result <- ds.vertDPQuantile(release, c(0, 0.5, 1))
  expect_identical(.Random.seed, before)
  expect_identical(attr(result, "additional_server_calls"), 0L)
  expect_identical(
    attr(result, "additional_privacy_cost"),
    c(epsilon = 0, delta = 0))
  expect_match(attr(result, "histogram_mass_scope"), "not silently forced")
  expect_match(attr(result, "uncertainty_scope"),
               "sampling uncertainty excluded")
})

test_that("quantile and median are registered as validated post-processors", {
  exports <- getNamespaceExports("dsVertClient")
  expect_contains(exports, c("ds.vertDPQuantile", "ds.vertDPMedian"))

  status <- ds.vertMethodStatus(c("ds.vertDPQuantile", "ds.vertDPMedian"))
  expect_true(all(status$status == "provisional"))
  expect_true(all(status$release_contract == "postprocessing_inherits_input"))
  expect_true(all(grepl("fixed public histogram bin",
                        status$principal_limitation, fixed = TRUE)))

  inventory <- .dsvert_capsule_method_inventory()
  rows <- inventory[inventory$method %in%
                      c("ds.vertDPQuantile", "ds.vertDPMedian"), ]
  expect_equal(nrow(rows), 2L)
  expect_true(all(rows$current_route_status ==
                    "client_only_validated_capsule_postprocess"))
  expect_false(any(rows$same_capsule_replay_history_can_deny))
  expect_false(any(rows$new_capsule_reservation_history_can_deny))
  expect_true(all(vapply(
    rows$legacy_remote_call_evidence, nrow, integer(1L)) == 0L))
})
