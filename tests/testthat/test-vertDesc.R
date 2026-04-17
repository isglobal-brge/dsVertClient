# Tests for .dsvert_interp_quantile (pure client-side helper for ds.vertDesc).

library(testthat)

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
