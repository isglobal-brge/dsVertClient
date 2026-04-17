# Tests for .dsvert_chisq_compute.

library(testthat)

test_that("chisq matches base R chisq.test on 2x2 without correction", {
  # Canonical Fisher / chi² 2x2 example
  observed <- matrix(c(12, 5, 2, 15), 2, 2, byrow = TRUE)
  row_m <- as.integer(rowSums(observed))
  col_m <- as.integer(colSums(observed))
  n <- sum(observed)

  res <- dsVertClient:::.dsvert_chisq_compute(
    observed, row_m, col_m, n, correct = FALSE)
  ref <- chisq.test(observed, correct = FALSE)

  expect_equal(res$statistic, unname(ref$statistic), tolerance = 1e-10)
  expect_equal(res$df, unname(ref$parameter))
  expect_equal(res$p_value, ref$p.value, tolerance = 1e-10)
  expect_equal(res$expected, ref$expected, tolerance = 1e-10)
})

test_that("chisq with Yates correction matches base R on 2x2", {
  observed <- matrix(c(12, 5, 2, 15), 2, 2, byrow = TRUE)
  row_m <- as.integer(rowSums(observed))
  col_m <- as.integer(colSums(observed))
  n <- sum(observed)

  res <- dsVertClient:::.dsvert_chisq_compute(
    observed, row_m, col_m, n, correct = TRUE)
  ref <- chisq.test(observed, correct = TRUE)

  expect_equal(res$statistic, unname(ref$statistic), tolerance = 1e-10)
  expect_equal(res$p_value, ref$p.value, tolerance = 1e-10)
})

test_that("chisq on larger table matches base R", {
  observed <- matrix(c(
    25, 30, 12,
    18, 22, 20,
     9, 15, 10), 3, 3, byrow = TRUE)
  row_m <- as.integer(rowSums(observed))
  col_m <- as.integer(colSums(observed))
  n <- sum(observed)

  res <- dsVertClient:::.dsvert_chisq_compute(
    observed, row_m, col_m, n, correct = TRUE)
  ref <- chisq.test(observed, correct = TRUE)  # correct ignored for 3x3

  expect_equal(res$statistic, unname(ref$statistic), tolerance = 1e-10)
  expect_equal(res$df, unname(ref$parameter))
  expect_equal(res$p_value, ref$p.value, tolerance = 1e-10)
  expect_false(res$correct)  # Yates only for 2x2
})

test_that("chisq residuals match (observed - expected) / sqrt(expected)", {
  observed <- matrix(c(12, 5, 2, 15), 2, 2, byrow = TRUE)
  row_m <- as.integer(rowSums(observed))
  col_m <- as.integer(colSums(observed))
  n <- sum(observed)

  res <- dsVertClient:::.dsvert_chisq_compute(
    observed, row_m, col_m, n, correct = FALSE)
  manual <- (observed - res$expected) / sqrt(res$expected)
  expect_equal(res$residuals, manual, tolerance = 1e-10)
})

test_that("chisq errors on empty table", {
  observed <- matrix(0L, 2, 2)
  expect_error(dsVertClient:::.dsvert_chisq_compute(
    observed, c(0L, 0L), c(0L, 0L), n = 0L),
    "no observations")
})
