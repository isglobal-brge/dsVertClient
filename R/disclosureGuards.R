#' Disclosure guard for low-dimensional aggregate matrix releases
#'
#' @param n_obs Number of observations contributing to the aggregate.
#' @param p Number of released variables.
#' @param what Human-readable aggregate name for error messages.
#' @return A list with the guard settings used.
#' @keywords internal
#' @noRd
.dsvert_guard_matrix_release <- function(n_obs, p, what = "aggregate matrix") {
  n_obs <- as.integer(n_obs)
  p <- as.integer(p)
  if (!is.finite(n_obs) || n_obs < 1L || !is.finite(p) || p < 1L) {
    stop("invalid dimensions for ", what, " disclosure guard",
         call. = FALSE)
  }
  min_n <- as.integer(getOption(
    "dsvert.min_aggregate_n",
    getOption("datashield.privacyLevel", 5L)))
  min_n_per_variable <- as.numeric(getOption(
    "dsvert.min_n_per_released_variable", 5))
  max_p_over_n <- as.numeric(getOption("dsvert.max_p_over_n", 0.2))
  allow_high_dim <- isTRUE(getOption("dsvert.allow_high_dim_aggregates",
                                     FALSE))

  required_n <- max(min_n, ceiling(min_n_per_variable * p))
  p_over_n <- p / n_obs
  if (!allow_high_dim &&
      (n_obs < required_n || p_over_n > max_p_over_n)) {
    stop(sprintf(
      paste0(
        "%s release blocked by disclosure guard: n=%d, p=%d, p/n=%.3f. ",
        "Require n >= max(%d, %.1f*p) and p/n <= %.3f. ",
        "Set options(dsvert.allow_high_dim_aggregates=TRUE) only for ",
        "controlled diagnostics."
      ),
      what, n_obs, p, p_over_n, min_n, min_n_per_variable, max_p_over_n),
      call. = FALSE)
  }
  list(
    n_obs = n_obs,
    p = p,
    p_over_n = p_over_n,
    min_n = min_n,
    min_n_per_variable = min_n_per_variable,
    max_p_over_n = max_p_over_n,
    allow_high_dim = allow_high_dim)
}

#' Disclosure guard for aggregate contingency-table releases
#'
#' @param counts Observed count matrix.
#' @param row_margins Optional row sums. If \code{NULL}, computed from
#'   \code{counts}.
#' @param col_margins Optional column sums. If \code{NULL}, computed from
#'   \code{counts}.
#' @param n Optional total count. If \code{NULL}, computed from counts.
#' @param what Human-readable aggregate name for error messages.
#' @return A list with guard metadata.
#' @keywords internal
#' @noRd
.dsvert_guard_table_release <- function(counts,
                                        row_margins = NULL,
                                        col_margins = NULL,
                                        n = NULL,
                                        what = "contingency table") {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("invalid counts for ", what, " disclosure guard", call. = FALSE)
  }
  row_m <- if (is.null(row_margins)) rowSums(counts) else as.numeric(row_margins)
  col_m <- if (is.null(col_margins)) colSums(counts) else as.numeric(col_margins)
  total <- if (is.null(n)) sum(counts) else as.numeric(n)
  if (any(!is.finite(row_m)) || any(row_m < 0) ||
      any(!is.finite(col_m)) || any(col_m < 0) ||
      !is.finite(total) || total < 0) {
    stop("invalid margins for ", what, " disclosure guard", call. = FALSE)
  }

  privacy_min <- suppressWarnings(as.numeric(getOption(
    "dsvert.min_cell_count",
    getOption("datashield.privacyLevel", 5L)))[1L])
  if (!is.finite(privacy_min)) privacy_min <- 0
  allow_small_cells <- isTRUE(getOption("dsvert.allow_small_cell_tables",
                                        FALSE))
  min_positive <- function(x) {
    x <- x[x > 0]
    if (length(x) == 0L) Inf else min(x)
  }
  guard <- list(
    min_cell_count = as.integer(privacy_min),
    allow_small_cell_tables = allow_small_cells,
    min_positive_cell = min_positive(counts),
    min_positive_row_margin = min_positive(row_m),
    min_positive_col_margin = min_positive(col_m),
    n = total)

  if (!allow_small_cells && is.finite(privacy_min) && privacy_min > 0L) {
    small_cell <- any(counts > 0 & counts < privacy_min)
    small_margin <- any(row_m > 0 & row_m < privacy_min) ||
      any(col_m > 0 & col_m < privacy_min) ||
      (total > 0 && total < privacy_min)
    if (small_cell || small_margin) {
      stop(sprintf(
        paste0(
          "%s release blocked by disclosure guard: positive cells and ",
          "margins must be either zero or >= %d. Set ",
          "options(dsvert.allow_small_cell_tables=TRUE) only for controlled ",
          "diagnostics."
        ),
        what, as.integer(privacy_min)),
        call. = FALSE)
    }
  }
  guard
}
