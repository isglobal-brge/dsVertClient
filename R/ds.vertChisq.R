#' @title Federated chi-square test on a 2-way contingency table
#' @description Compute Pearson chi-square (and optional Yates and LR
#'   variants) on a contingency table of two categorical variables held
#'   at the same server in a DataSHIELD federation. The server returns
#'   only aggregate cell counts and margins (with
#'   datashield.privacyLevel-aware small-cell suppression); chi-square
#'   statistics are then computed entirely client-side.
#'
#' @param data_name Character. Name of the aligned data frame on each server.
#' @param var1 Character. Row variable.
#' @param var2 Character. Column variable.
#' @param server Character. Optional name of the server that holds both
#'   variables. If \code{NULL} the server is auto-detected via
#'   \code{dsvertColNamesDS}.
#' @param correct Logical. If TRUE (default for 2x2 tables), apply Yates'
#'   continuity correction.
#' @param datasources DataSHIELD connections.
#'
#' @return An object of class \code{ds.vertChisq} with elements
#'   \code{statistic}, \code{df}, \code{p_value}, \code{observed},
#'   \code{expected}, \code{residuals}, \code{n}, \code{correct}.
#'
#' @details For the cross-server case (variables vertically partitioned
#'   across two or more servers) a separate entry point
#'   \code{ds.vertChisqCross} will reuse the Beaver dot-product infrastructure
#'   to compute cell counts on secret shares; that is tracked as a follow-on
#'   in the dsVert v2.0 implementation plan.
#'
#' @importFrom stats pchisq
#' @export
ds.vertChisq <- function(data_name, var1, var2, server = NULL,
                         correct = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (!is.character(var1) || length(var1) != 1L ||
      !is.character(var2) || length(var2) != 1L) {
    stop("var1 and var2 must be single character strings", call. = FALSE)
  }

  if (is.null(server)) {
    col_results <- DSI::datashield.aggregate(datasources,
      call(name = "dsvertColNamesDS", data_name = data_name))
    candidates <- character(0)
    for (srv in names(datasources)) {
      cols <- col_results[[srv]]$columns
      if (all(c(var1, var2) %in% cols)) {
        candidates <- c(candidates, srv)
      }
    }
    if (length(candidates) == 0L) {
      stop("No server holds both '", var1, "' and '", var2,
           "'. Cross-server chi-square is a separate routine (planned).",
           call. = FALSE)
    }
    if (length(candidates) > 1L) {
      stop("Variables found on multiple servers (",
           paste(candidates, collapse = ", "),
           "); specify `server = '...'` to disambiguate.",
           call. = FALSE)
    }
    server <- candidates
  }

  if (!server %in% names(datasources)) {
    stop("Server '", server, "' not in datasources", call. = FALSE)
  }

  ci <- which(names(datasources) == server)
  tab <- DSI::datashield.aggregate(datasources[ci],
    call(name = "dsvertContingencyDS", data_name = data_name,
         var1 = var1, var2 = var2,
         suppress_small_cells = TRUE))[[1]]
  disclosure_guard <- .dsvert_guard_table_release(
    tab$counts,
    row_margins = tab$row_margins,
    col_margins = tab$col_margins,
    n = tab$n,
    what = "same-server contingency table")

  stats <- .dsvert_chisq_compute(
    observed = tab$counts,
    row_margins = tab$row_margins,
    col_margins = tab$col_margins,
    n = tab$n,
    correct = correct)

  rownames(stats$observed) <- tab$row_levels
  colnames(stats$observed) <- tab$col_levels
  rownames(stats$expected) <- tab$row_levels
  colnames(stats$expected) <- tab$col_levels
  rownames(stats$residuals) <- tab$row_levels
  colnames(stats$residuals) <- tab$col_levels

  out <- c(stats, list(
    n_na = as.integer(tab$n_na),
    server = server,
    var1 = var1,
    var2 = var2,
    disclosure_guard = disclosure_guard))
  class(out) <- c("ds.vertChisq", "list")
  out
}

#' Pure helper: compute chi-square statistics from observed counts + margins.
#' Exposed as an internal function so unit tests can exercise the math path
#' without a DataSHIELD round trip.
#'
#' @param observed Integer matrix of cell counts.
#' @param row_margins Integer vector of row sums.
#' @param col_margins Integer vector of column sums.
#' @param n Total observation count.
#' @param correct Apply Yates continuity correction for 2x2 tables.
#' @return list with statistic, df, p_value, observed, expected,
#'   residuals, n, correct.
#' @keywords internal
.dsvert_chisq_compute <- function(observed, row_margins, col_margins, n,
                                  correct = TRUE) {
  observed <- as.matrix(observed)
  n <- as.numeric(n)
  if (n < 1 || sum(observed) < 1) {
    stop("Contingency table has no observations", call. = FALSE)
  }
  row_m <- as.numeric(row_margins)
  col_m <- as.numeric(col_margins)
  expected <- outer(row_m, col_m) / n

  dims <- dim(observed)
  is_2x2 <- identical(dims, c(2L, 2L))
  correct <- isTRUE(correct) && is_2x2

  eps <- .Machine$double.eps
  if (correct) {
    diff <- abs(observed - expected) - 0.5
    diff[diff < 0] <- 0
    stat <- sum(diff^2 / pmax(expected, eps))
  } else {
    stat <- sum((observed - expected)^2 / pmax(expected, eps))
  }
  df <- (dims[1] - 1L) * (dims[2] - 1L)
  p <- pchisq(stat, df = df, lower.tail = FALSE)
  residuals <- (observed - expected) / sqrt(pmax(expected, eps))

  list(
    statistic = stat,
    df = as.integer(df),
    p_value = p,
    observed = observed,
    expected = expected,
    residuals = residuals,
    n = n,
    correct = correct)
}

#' @title Federated Fisher exact test (same-server case)
#' @description Exact conditional p-value for a 2x2 (or small 2xK) contingency
#'   table of two variables held at the same server. Uses the
#'   \code{stats::fisher.test} implementation on the released count matrix
#'   (itself already subject to datashield.privacyLevel small-cell
#'   suppression), so this is a safe exact-inference complement to
#'   ds.vertChisq when any expected cell count is below the chi-square
#'   validity floor.
#'
#' @param data_name,var1,var2,server,datasources Same semantics as in
#'   \code{ds.vertChisq}.
#' @param alternative One of "two.sided" (default), "greater", "less"; only
#'   consulted for 2x2 tables.
#' @param conf.int Logical. Return an odds-ratio confidence interval for 2x2.
#' @param conf.level Confidence level for 2x2 OR interval.
#'
#' @return An object of class \code{ds.vertFisher} (a
#'   \code{htest}-compatible list plus table metadata).
#'
#' @importFrom stats fisher.test
#' @export
ds.vertFisher <- function(data_name, var1, var2, server = NULL,
                          alternative = c("two.sided", "greater", "less"),
                          conf.int = TRUE, conf.level = 0.95,
                          datasources = NULL) {
  alternative <- match.arg(alternative)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  if (is.null(server)) {
    col_results <- DSI::datashield.aggregate(datasources,
      call(name = "dsvertColNamesDS", data_name = data_name))
    candidates <- character(0)
    for (srv in names(datasources)) {
      cols <- col_results[[srv]]$columns
      if (all(c(var1, var2) %in% cols)) candidates <- c(candidates, srv)
    }
    if (length(candidates) == 0L) {
      stop("No server holds both '", var1, "' and '", var2, "'",
           call. = FALSE)
    }
    if (length(candidates) > 1L) {
      stop("Variables found on multiple servers; specify `server = '...'`",
           call. = FALSE)
    }
    server <- candidates
  }

  ci <- which(names(datasources) == server)
  tab <- DSI::datashield.aggregate(datasources[ci],
    call(name = "dsvertContingencyDS", data_name = data_name,
         var1 = var1, var2 = var2,
         suppress_small_cells = TRUE))[[1]]
  disclosure_guard <- .dsvert_guard_table_release(
    tab$counts,
    row_margins = tab$row_margins,
    col_margins = tab$col_margins,
    n = tab$n,
    what = "same-server Fisher table")

  observed <- tab$counts
  rownames(observed) <- tab$row_levels
  colnames(observed) <- tab$col_levels
  if (sum(observed) < 1) {
    stop("Contingency table has no observations after disclosure filtering",
         call. = FALSE)
  }

  is_2x2 <- identical(dim(observed), c(2L, 2L))
  fish <- if (is_2x2) {
    fisher.test(observed, alternative = alternative,
                conf.int = conf.int, conf.level = conf.level)
  } else {
    fisher.test(observed)
  }

  out <- list(
    p_value = fish$p.value,
    odds_ratio = if (!is.null(fish$estimate)) unname(fish$estimate) else NA_real_,
    conf_int = if (!is.null(fish$conf.int)) as.numeric(fish$conf.int) else NULL,
    alternative = alternative,
    observed = observed,
    n = as.integer(tab$n),
    n_na = as.integer(tab$n_na),
    server = server,
    var1 = var1,
    var2 = var2,
    disclosure_guard = disclosure_guard)
  class(out) <- c("ds.vertFisher", "list")
  out
}

#' @export
print.ds.vertFisher <- function(x, ...) {
  cat(sprintf("dsVert Fisher exact test on %s x %s (server: %s)\n",
              x$var1, x$var2, x$server))
  cat(sprintf("  n = %d (n_na = %d)\n", x$n, x$n_na))
  cat(sprintf("  p-value = %s  (alternative: %s)\n",
              format.pval(x$p_value, digits = 4L), x$alternative))
  if (!is.na(x$odds_ratio)) {
    cat(sprintf("  odds ratio = %.4f", x$odds_ratio))
    if (!is.null(x$conf_int)) {
      cat(sprintf("  CI: [%.4f, %.4f]", x$conf_int[1], x$conf_int[2]))
    }
    cat("\n")
  }
  cat("\nObserved counts:\n")
  print(x$observed)
  invisible(x)
}

#' @export
print.ds.vertChisq <- function(x, ...) {
  cat(sprintf("dsVert chi-square test on %s x %s (server: %s)\n",
              x$var1, x$var2, x$server))
  if (x$correct) cat("  (Yates continuity correction applied)\n")
  cat(sprintf("  n = %d (n_na = %d)\n", as.integer(x$n), x$n_na))
  cat(sprintf("  chi-sq = %.4f  df = %d  p-value = %s\n",
              x$statistic, x$df, format.pval(x$p_value, digits = 4L)))
  cat("\nObserved counts:\n")
  print(x$observed)
  invisible(x)
}
