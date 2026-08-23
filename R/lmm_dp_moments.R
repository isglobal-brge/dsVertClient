# Internal projection and method-of-moments core for the future sticky-DP
# random-intercept LMM artifact.  It deliberately accepts only a fixed,
# already released six-coordinate summary; it has no DataSHIELD or source-data
# access path.

.dsvert_lmm_dp_scalar <- function(value, what) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value)) {
    stop(what, " must be one finite number", call. = FALSE)
  }
  as.numeric(value)
}

.dsvert_lmm_random_intercept_moments <- function(
    summary, outcome_lower, outcome_upper, observation_capacity,
    cluster_capacity) {
  fields <- c(
    "n", "clusters", "cluster_size_sq", "sum_y_normalized",
    "sum_y_sq_normalized", "sum_cluster_mean_sq_normalized")
  if (!is.list(summary) || is.null(names(summary)) || anyNA(names(summary)) ||
      anyDuplicated(names(summary)) || !setequal(names(summary), fields)) {
    stop("The LMM moment summary has an invalid schema", call. = FALSE)
  }
  values <- vapply(fields, function(field) {
    .dsvert_lmm_dp_scalar(summary[[field]], paste("LMM", field))
  }, numeric(1L))
  outcome_lower <- .dsvert_lmm_dp_scalar(outcome_lower, "outcome_lower")
  outcome_upper <- .dsvert_lmm_dp_scalar(outcome_upper, "outcome_upper")
  observation_capacity <- .dsvert_lmm_dp_scalar(
    observation_capacity, "observation_capacity")
  cluster_capacity <- .dsvert_lmm_dp_scalar(
    cluster_capacity, "cluster_capacity")
  if (outcome_lower >= outcome_upper || observation_capacity < 2 ||
      observation_capacity != floor(observation_capacity) ||
      cluster_capacity < 2 || cluster_capacity != floor(cluster_capacity) ||
      cluster_capacity > observation_capacity) {
    stop("The LMM moment bounds are invalid", call. = FALSE)
  }

  n <- min(max(values[["n"]], 0), observation_capacity)
  clusters <- min(max(values[["clusters"]], 0), n)
  size_sq <- min(max(values[["cluster_size_sq"]], n),
                 min(cluster_capacity * n, n^2))
  sum_y <- min(max(values[["sum_y_normalized"]], 0), n)
  min_second <- if (n > 0) sum_y^2 / n else 0
  sum_y_sq <- min(max(values[["sum_y_sq_normalized"]], min_second), sum_y)
  cluster_mean_sq <- min(max(values[["sum_cluster_mean_sq_normalized"]],
                             min_second), sum_y_sq)
  projected <- c(
    n = n, clusters = clusters, cluster_size_sq = size_sq,
    sum_y_normalized = sum_y, sum_y_sq_normalized = sum_y_sq,
    sum_cluster_mean_sq_normalized = cluster_mean_sq)
  projection_applied <- !isTRUE(all.equal(values, projected,
                                           tolerance = 1e-12))

  within_df <- n - clusters
  between_df <- clusters - 1
  effective_cluster_size <- if (between_df > 0) {
    (n - size_sq / n) / between_df
  } else 0
  if (n <= 0 || clusters < 2 || within_df <= 0 || between_df <= 0 ||
      !is.finite(effective_cluster_size) || effective_cluster_size <= 0) {
    return(structure(list(
      status = "non_identifiable", projected_summary = projected,
      projection_applied = projection_applied,
      reason = "insufficient_cluster_information"),
      class = c("dsvert_lmm_dp_moments", "list")))
  }

  within_ss <- max(sum_y_sq - cluster_mean_sq, 0)
  between_ss <- max(cluster_mean_sq - sum_y^2 / n, 0)
  sigma_e2_normalized <- within_ss / within_df
  sigma_b2_normalized <- max(
    (between_ss / between_df - sigma_e2_normalized) /
      effective_cluster_size,
    0)
  span <- outcome_upper - outcome_lower
  sigma_e2 <- span^2 * sigma_e2_normalized
  sigma_b2 <- span^2 * sigma_b2_normalized
  total_variance <- sigma_e2 + sigma_b2
  if (!is.finite(sigma_e2) || !is.finite(sigma_b2) || sigma_e2 < 0 ||
      sigma_b2 < 0 || !is.finite(total_variance)) {
    stop("The projected LMM moments are numerically invalid", call. = FALSE)
  }
  structure(list(
    status = "ok",
    coefficient = stats::setNames(
      outcome_lower + span * sum_y / n, "(Intercept)"),
    sigma2 = sigma_e2,
    sigma_b2 = sigma_b2,
    icc = if (total_variance > 0) sigma_b2 / total_variance else 0,
    effective_cluster_size = effective_cluster_size,
    projected_summary = projected,
    projection_applied = projection_applied,
    estimand = "bounded_random_intercept_method_of_moments"),
    class = c("dsvert_lmm_dp_moments", "list"))
}
