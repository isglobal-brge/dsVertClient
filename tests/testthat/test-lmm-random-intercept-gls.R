make_lmm_gls_summary <- function(y, x, cluster) {
  design <- cbind("(Intercept)" = 1, x = x)
  sizes <- tabulate(cluster)
  dimension <- ncol(design)
  by_size <- lapply(seq_len(max(sizes)), function(size) {
    members <- which(sizes == size)
    cluster_xtx <- matrix(0, dimension, dimension,
                          dimnames = list(colnames(design), colnames(design)))
    cluster_xty <- stats::setNames(rep(0, dimension), colnames(design))
    cluster_yty <- 0
    for (member in members) {
      rows <- which(cluster == member)
      sx <- colSums(design[rows, , drop = FALSE])
      sy <- sum(y[rows])
      cluster_xtx <- cluster_xtx + tcrossprod(sx)
      cluster_xty <- cluster_xty + sx * sy
      cluster_yty <- cluster_yty + sy^2
    }
    list(count = length(members), xtx = cluster_xtx,
         xty = cluster_xty, yty = cluster_yty)
  })
  list(
    global = list(n = length(y), xtx = crossprod(design),
                  xty = drop(crossprod(design, y)), yty = sum(y^2)),
    by_size = by_size)
}

explicit_lmm_gls <- function(y, x, cluster, ratio) {
  design <- cbind("(Intercept)" = 1, x = x)
  covariance <- diag(length(y))
  for (level in unique(cluster)) {
    rows <- which(cluster == level)
    covariance[rows, rows] <- covariance[rows, rows] + ratio
  }
  inverse <- solve(covariance)
  information <- crossprod(design, inverse %*% design)
  score <- crossprod(design, inverse %*% y)
  coefficients <- drop(solve(information, score))
  residual <- drop(crossprod(y - design %*% coefficients,
                             inverse %*% (y - design %*% coefficients)))
  list(coefficients = coefficients, sigma2 = residual / length(y),
       objective = determinant(covariance, logarithm = TRUE)$modulus +
         length(y) * log(residual / length(y)))
}

explicit_lmm_reml <- function(y, x, cluster, ratio) {
  design <- cbind("(Intercept)" = 1, x = x)
  covariance <- diag(length(y))
  for (level in unique(cluster)) {
    rows <- which(cluster == level)
    covariance[rows, rows] <- covariance[rows, rows] + ratio
  }
  inverse <- solve(covariance)
  information <- crossprod(design, inverse %*% design)
  score <- crossprod(design, inverse %*% y)
  coefficients <- drop(solve(information, score))
  residual <- drop(crossprod(y - design %*% coefficients,
                             inverse %*% (y - design %*% coefficients)))
  degrees <- length(y) - ncol(design)
  sigma2 <- residual / degrees
  list(coefficients = coefficients, sigma2 = sigma2,
       objective = determinant(covariance, logarithm = TRUE)$modulus +
         determinant(information, logarithm = TRUE)$modulus +
         degrees * log(sigma2))
}

test_that("random-intercept GLS uses only cluster sufficient statistics", {
  y <- c(1.0, 1.4, 2.7, 3.1, 2.8, 4.2, 4.7)
  x <- c(0.1, 0.5, 0.2, 0.6, 0.9, 0.3, 0.8)
  cluster <- c(1L, 1L, 2L, 2L, 2L, 3L, 3L)
  aggregate <- make_lmm_gls_summary(y, x, cluster)
  grid <- c(0, 0.25, 1, 4)

  fit <- dsVertClient:::.dsvert_lmm_random_intercept_gls(
    aggregate$global, aggregate$by_size, grid)
  expect_identical(fit$status, "ok")
  candidates <- lapply(grid, explicit_lmm_gls, y = y, x = x, cluster = cluster)
  objectives <- vapply(candidates, `[[`, numeric(1L), "objective")
  selected <- candidates[[which.min(objectives)]]
  expect_equal(fit$variance_ratio, grid[[which.min(objectives)]])
  expect_equal(unname(fit$coefficients), unname(selected$coefficients),
               tolerance = 1e-10)
  expect_equal(fit$sigma2, selected$sigma2, tolerance = 1e-10)
})

test_that("random-intercept GLS supports signed finite-grid REML", {
  y <- c(1.0, 1.4, 2.7, 3.1, 2.8, 4.2, 4.7)
  x <- c(0.1, 0.5, 0.2, 0.6, 0.9, 0.3, 0.8)
  cluster <- c(1L, 1L, 2L, 2L, 2L, 3L, 3L)
  aggregate <- make_lmm_gls_summary(y, x, cluster)
  grid <- c(0, 0.25, 1, 4)

  fit <- dsVertClient:::.dsvert_lmm_random_intercept_gls(
    aggregate$global, aggregate$by_size, grid, estimation_profile = "reml")
  expect_identical(fit$status, "ok")
  candidates <- lapply(grid, explicit_lmm_reml, y = y, x = x, cluster = cluster)
  objectives <- vapply(candidates, `[[`, numeric(1L), "objective")
  selected <- candidates[[which.min(objectives)]]
  expect_identical(fit$estimation_profile, "reml")
  expect_equal(fit$variance_ratio, grid[[which.min(objectives)]])
  expect_equal(unname(fit$coefficients), unname(selected$coefficients),
               tolerance = 1e-10)
  expect_equal(fit$sigma2, selected$sigma2, tolerance = 1e-10)
})

test_that("random-intercept GLS fails closed on malformed or singular summaries", {
  y <- c(1, 2, 3, 4)
  x <- c(0, 1, 0, 1)
  aggregate <- make_lmm_gls_summary(y, x, c(1L, 1L, 2L, 2L))
  aggregate$global$xtx[1, 2] <- aggregate$global$xtx[1, 2] + 1
  expect_error(
    dsVertClient:::.dsvert_lmm_random_intercept_gls(
      aggregate$global, aggregate$by_size, c(0, 1)),
    "global information")

  singular <- make_lmm_gls_summary(y, rep(1, length(y)), c(1L, 1L, 2L, 2L))
  fit <- dsVertClient:::.dsvert_lmm_random_intercept_gls(
    singular$global, singular$by_size, c(0, 1))
  expect_identical(fit$status, "non_identifiable")

  incomplete <- make_lmm_gls_summary(y, x, c(1L, 1L, 2L, 2L))
  incomplete$by_size[[2L]]$count <- 0L
  fit <- dsVertClient:::.dsvert_lmm_random_intercept_gls(
    incomplete$global, incomplete$by_size, c(0, 1))
  expect_identical(fit$status, "non_identifiable")
  expect_identical(fit$reason, "inconsistent_cluster_counts")
})
