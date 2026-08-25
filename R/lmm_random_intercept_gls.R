# Deterministic post-processing for a signed random-intercept LMM synopsis.
#
# The input is intentionally limited to aggregate cross-products.  It cannot
# accept rows, cluster labels, or a caller-selected optimiser trace.  A future
# sticky artifact supplies these aggregates; this helper only evaluates its
# signed, finite variance-ratio grid.

.dsvert_lmm_gls_non_identifiable <- function(reason) {
  structure(list(
    status = "non_identifiable",
    reason = as.character(reason),
    coefficients = NULL,
    sigma2 = NULL,
    sigma_b2 = NULL,
    icc = NULL,
    variance_ratio = NULL,
    objective = NULL),
    class = c("dsvert_lmm_gls_summary", "list"))
}

.dsvert_lmm_gls_square <- function(value, dimension, field) {
  valid <- is.numeric(value) && is.matrix(value) &&
    identical(dim(value), c(dimension, dimension)) &&
    all(is.finite(value)) &&
    isTRUE(all.equal(value, t(value), tolerance = 1e-10))
  if (!isTRUE(valid)) {
    stop(paste0("The signed LMM ", field, " is invalid."), call. = FALSE)
  }
  (value + t(value)) / 2
}

.dsvert_lmm_gls_vector <- function(value, dimension, field) {
  if (!is.numeric(value) || length(value) != dimension ||
      any(!is.finite(value))) {
    stop(paste0("The signed LMM ", field, " is invalid."), call. = FALSE)
  }
  as.numeric(value)
}

#' Evaluate a finite-grid random-intercept GLS synopsis
#'
#' @keywords internal
.dsvert_lmm_random_intercept_gls <- function(
    global, by_cluster_size, variance_ratio_grid) {
  required_global <- c("n", "xtx", "xty", "yty")
  required_size <- c("count", "xtx", "xty", "yty")
  if (!is.list(global) || !identical(names(global), required_global) ||
      !is.list(by_cluster_size) || !length(by_cluster_size) ||
      !is.numeric(variance_ratio_grid) || !length(variance_ratio_grid) ||
      anyNA(variance_ratio_grid) || any(!is.finite(variance_ratio_grid)) ||
      any(variance_ratio_grid < 0) || any(diff(variance_ratio_grid) <= 0)) {
    stop("The signed random-intercept GLS synopsis is invalid.", call. = FALSE)
  }
  n <- suppressWarnings(as.integer(global$n))
  if (!is.numeric(global$n) || length(global$n) != 1L || is.na(n) ||
      n < 2L || global$n != n || !is.matrix(global$xtx) ||
      nrow(global$xtx) != ncol(global$xtx) || nrow(global$xtx) < 1L) {
    stop("The signed random-intercept GLS cohort is invalid.", call. = FALSE)
  }
  dimension <- nrow(global$xtx)
  xtx <- .dsvert_lmm_gls_square(global$xtx, dimension, "global information")
  xty <- .dsvert_lmm_gls_vector(global$xty, dimension, "global score")
  yty <- suppressWarnings(as.numeric(global$yty))
  if (length(yty) != 1L || !is.finite(yty)) {
    stop("The signed LMM global sum of squares is invalid.", call. = FALSE)
  }

  sizes <- seq_along(by_cluster_size)
  counts <- numeric(length(sizes))
  size_xtx <- vector("list", length(sizes))
  size_xty <- vector("list", length(sizes))
  size_yty <- numeric(length(sizes))
  for (index in seq_along(by_cluster_size)) {
    component <- by_cluster_size[[index]]
    if (!is.list(component) || !identical(names(component), required_size)) {
      stop("The signed LMM cluster summary is invalid.", call. = FALSE)
    }
    count <- suppressWarnings(as.integer(component$count))
    if (!is.numeric(component$count) || length(component$count) != 1L ||
        is.na(count) || count < 0L || component$count != count) {
      stop("The signed LMM cluster count is invalid.", call. = FALSE)
    }
    counts[[index]] <- count
    size_xtx[[index]] <- .dsvert_lmm_gls_square(
      component$xtx, dimension, "cluster information")
    size_xty[[index]] <- .dsvert_lmm_gls_vector(
      component$xty, dimension, "cluster score")
    size_yty[[index]] <- suppressWarnings(as.numeric(component$yty))
    if (length(size_yty[[index]]) != 1L || !is.finite(size_yty[[index]])) {
      stop("The signed LMM cluster sum of squares is invalid.", call. = FALSE)
    }
  }
  if (sum(counts * sizes) != n || !any(counts > 0L)) {
    return(.dsvert_lmm_gls_non_identifiable("inconsistent_cluster_counts"))
  }

  candidates <- vector("list", length(variance_ratio_grid))
  for (index in seq_along(variance_ratio_grid)) {
    ratio <- variance_ratio_grid[[index]]
    weight <- ratio / (1 + ratio * sizes)
    information <- xtx
    score <- xty
    sumsq <- yty
    for (size_index in seq_along(sizes)) {
      information <- information - weight[[size_index]] * size_xtx[[size_index]]
      score <- score - weight[[size_index]] * size_xty[[size_index]]
      sumsq <- sumsq - weight[[size_index]] * size_yty[[size_index]]
    }
    inverse <- tryCatch(solve(information), error = function(error) NULL)
    if (is.null(inverse) || any(!is.finite(inverse))) next
    coefficients <- drop(inverse %*% score)
    residual <- sumsq - 2 * sum(coefficients * score) +
      drop(crossprod(coefficients, information %*% coefficients))
    if (!is.finite(residual) || residual <= 0) next
    sigma2 <- residual / n
    if (!is.finite(sigma2) || sigma2 <= 0) next
    objective <- sum(counts * log1p(ratio * sizes)) + n * log(sigma2)
    if (!is.finite(objective)) next
    candidates[[index]] <- list(
      variance_ratio = ratio, objective = objective,
      coefficients = stats::setNames(coefficients, colnames(xtx) %||%
        paste0("x", seq_len(dimension))),
      information = information, sigma2 = sigma2)
  }
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) {
    return(.dsvert_lmm_gls_non_identifiable("singular_or_nonpositive_profile"))
  }
  objectives <- vapply(candidates, `[[`, numeric(1L), "objective")
  selected <- candidates[[which.min(objectives)]]
  ratio <- selected$variance_ratio
  structure(list(
    status = "ok",
    reason = NULL,
    coefficients = selected$coefficients,
    information = selected$information,
    sigma2 = selected$sigma2,
    sigma_b2 = ratio * selected$sigma2,
    icc = ratio / (1 + ratio),
    variance_ratio = ratio,
    objective = selected$objective),
    class = c("dsvert_lmm_gls_summary", "list"))
}
