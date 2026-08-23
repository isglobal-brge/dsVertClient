.lmm_test_summary <- function(y, cluster, lower = -3, upper = 5) {
  normalized <- (y - lower) / (upper - lower)
  sizes <- as.numeric(table(cluster))
  sums <- tapply(normalized, cluster, sum)
  list(
    n = length(normalized),
    clusters = length(sizes),
    cluster_size_sq = sum(sizes^2),
    sum_y_normalized = sum(normalized),
    sum_y_sq_normalized = sum(normalized^2),
    sum_cluster_mean_sq_normalized = sum(sums^2 / sizes))
}

test_that("bounded random-intercept moments recover the exact ANOVA estimate", {
  y <- c(-1.8, -1.1, -0.9, 0.4, 0.9, 1.1, 2.0, 1.7, 2.4)
  cluster <- c("a", "a", "a", "b", "b", "c", "c", "c", "c")
  summary <- .lmm_test_summary(y, cluster)
  result <- dsVertClient:::.dsvert_lmm_random_intercept_moments(
    summary, outcome_lower = -3, outcome_upper = 5,
    observation_capacity = 32, cluster_capacity = 8)

  sizes <- as.numeric(table(cluster))
  means <- as.numeric(tapply(y, cluster, mean))
  grand <- mean(y)
  within <- sum((y - ave(y, cluster, FUN = mean))^2)
  between <- sum(sizes * (means - grand)^2)
  msw <- within / (length(y) - length(sizes))
  msb <- between / (length(sizes) - 1)
  n0 <- (length(y) - sum(sizes^2) / length(y)) / (length(sizes) - 1)

  expect_identical(result$status, "ok")
  expect_false(result$projection_applied)
  expect_equal(unname(result$coefficient), grand, tolerance = 1e-12)
  expect_equal(result$sigma2, msw, tolerance = 1e-12)
  expect_equal(result$sigma_b2, max((msb - msw) / n0, 0), tolerance = 1e-12)
  expect_equal(result$icc,
               result$sigma_b2 / (result$sigma2 + result$sigma_b2),
               tolerance = 1e-12)
})

test_that("LMM moment projection produces bounded plausible components", {
  result <- dsVertClient:::.dsvert_lmm_random_intercept_moments(
    list(
      n = 1000, clusters = 3, cluster_size_sq = -1,
      sum_y_normalized = -2, sum_y_sq_normalized = 100,
      sum_cluster_mean_sq_normalized = -10),
    outcome_lower = 0, outcome_upper = 2,
    observation_capacity = 20, cluster_capacity = 4)

  expect_identical(result$status, "ok")
  expect_true(result$projection_applied)
  expect_true(all(is.finite(c(result$coefficient, result$sigma2,
                              result$sigma_b2, result$icc))))
  expect_gte(result$sigma2, 0)
  expect_gte(result$sigma_b2, 0)
  expect_true(result$icc >= 0 && result$icc <= 1)
  expect_lte(result$projected_summary[["n"]], 20)
  expect_lte(result$projected_summary[["cluster_size_sq"]],
             4 * result$projected_summary[["n"]])
})

test_that("LMM moment core fails closed without cluster degrees of freedom", {
  result <- dsVertClient:::.dsvert_lmm_random_intercept_moments(
    list(
      n = 1, clusters = 1, cluster_size_sq = 1,
      sum_y_normalized = 0.5, sum_y_sq_normalized = 0.25,
      sum_cluster_mean_sq_normalized = 0.25),
    outcome_lower = 0, outcome_upper = 1,
    observation_capacity = 8, cluster_capacity = 4)
  expect_identical(result$status, "non_identifiable")
  expect_null(result$coefficient)
  expect_error(
    dsVertClient:::.dsvert_lmm_random_intercept_moments(
      list(n = 1), 0, 1, 8, 4),
    "invalid schema")
})
