#' @title Federated descriptive statistics with approximate quantiles
#' @description Compute a federated \code{summary()}-style table of numeric
#'   variables held across vertically partitioned servers. Exact means and
#'   standard deviations are obtained from each server's local moments;
#'   exact extrema are suppressed by default because they can disclose
#'   outliers. Quantiles are interpolated from disclosure-checked histogram
#'   bucket counts. No observation-level quantity is reconstructed at the
#'   client.
#'
#' @param data_name Character. Name of the aligned data frame held on each
#'   server.
#' @param variables Optional: a character vector of variables of interest,
#'   a named list mapping server -> variables, or \code{NULL} (default) to
#'   auto-detect.
#' @param probs Numeric vector of quantile probabilities to report.
#'   Defaults to the usual tertiles (0.25, 0.5, 0.75).
#' @param n_buckets Integer. Number of uniform histogram buckets used for
#'   quantile interpolation. Higher values give tighter quantile resolution
#'   at the cost of one extra aggregate call per variable; when
#'   \code{datashield.privacyLevel > 0}, the effective bucket count is capped
#'   by the cohort size and privacy threshold.
#' @param range_sd Numeric. When exact extrema are not released, histogram
#'   edges are built over \code{mean +/- range_sd * sd}. Values outside that
#'   range are retained as under/overflow aggregate counts.
#' @param exact_extrema Logical. Return exact min/max and use them for the
#'   histogram range. Defaults to
#'   \code{getOption("dsvert.allow_exact_extrema", FALSE)}.
#' @param open_ended Logical. When exact extrema are suppressed, use
#'   open-ended first/last histogram buckets so outliers are absorbed into
#'   coarse tail aggregates rather than disclosed as separate under/overflow
#'   counts.
#' @param verbose Logical. Print per-variable progress when TRUE.
#' @param datasources DataSHIELD connections. If \code{NULL}, auto-detected.
#'
#' @return A data frame with one row per variable containing columns:
#'   \itemize{
#'     \item \code{server}: which server holds the variable
#'     \item \code{variable}: column name
#'     \item \code{n}: number of non-missing observations
#'     \item \code{n_na}: number of missing observations
#'     \item \code{mean}, \code{sd}: local moments
#'     \item \code{min}, \code{max}: exact extrema only when
#'       \code{exact_extrema = TRUE}; otherwise \code{NA}
#'     \item \code{range_low}, \code{range_high}: histogram range used for
#'       approximate quantiles
#'     \item \code{quantile_status}: \code{"ran"} or the suppression reason
#'     \item one column per requested quantile (named \code{q25}, \code{q50},
#'       \code{q75} etc. by default)
#'   }
#'
#' @details
#' Per-variable flow:
#'   1. \code{dsvertLocalMomentsDS} returns scalar moments in a single
#'      aggregate call; exact extrema are omitted unless explicitly enabled.
#'   2. Edges for \code{dsvertHistogramDS} are derived from either exact
#'      extrema (opt-in) or a moment-bounded range
#'      \code{mean +/- range_sd * sd}.
#'   3. Client-side linear interpolation within the target bucket converts
#'      the cumulative count to the requested quantile.
#'
#' The quantile estimate inherits the histogram's resolution; for
#' \code{n_buckets = 100} on a roughly normal cohort, the median is typically
#' accurate within 0.5\% of the exact value, tightening with \code{n_buckets}
#' at linear extra cost.
#'
#' @examples
#' \dontrun{
#' conns <- DSI::datashield.login(logindata)
#' ds.vertDesc("DA", variables = c("age", "bmi", "glu", "bp"),
#'             probs = c(0.25, 0.5, 0.75, 0.9))
#' }
#' @export
ds.vertDesc <- function(data_name,
                        variables = NULL,
                        probs = c(0.25, 0.5, 0.75),
                        n_buckets = 100L,
                        range_sd = getOption("dsvert.desc_range_sd", 4),
                        exact_extrema = getOption("dsvert.allow_exact_extrema",
                                                  FALSE),
                        open_ended = getOption("dsvert.desc_open_ended", TRUE),
                        verbose = TRUE,
                        datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (!is.numeric(probs) || any(probs <= 0) || any(probs >= 1)) {
    stop("probs must be a numeric vector strictly between 0 and 1",
         call. = FALSE)
  }
  n_buckets <- as.integer(n_buckets)
  if (n_buckets < 2L) {
    stop("n_buckets must be >= 2", call. = FALSE)
  }
  range_sd <- as.numeric(range_sd)
  if (!is.finite(range_sd) || range_sd <= 0) {
    stop("range_sd must be a positive finite number", call. = FALSE)
  }
  exact_extrema <- isTRUE(exact_extrema)
  open_ended <- isTRUE(open_ended)

  # Build variable -> server map ---------------------------------------
  if (is.null(variables) || is.character(variables)) {
    user_vars <- variables
    col_results <- DSI::datashield.aggregate(
      datasources, call(name = "dsvertColNamesDS", data_name = data_name))
    var_map <- list()
    for (srv in names(datasources)) {
      feats <- setdiff(col_results[[srv]]$columns, c("id", "patient_id"))
      if (!is.null(user_vars)) feats <- intersect(feats, user_vars)
      if (length(feats) > 0) var_map[[srv]] <- feats
    }
  } else if (is.list(variables) && !is.null(names(variables))) {
    var_map <- variables
  } else {
    stop("variables must be NULL, a character vector, or a named list",
         call. = FALSE)
  }

  server_names <- names(datasources)
  rows <- list()

  q_names <- sprintf("q%02d", as.integer(round(100 * probs)))

  for (srv in names(var_map)) {
    ci <- which(server_names == srv)
    if (length(ci) == 0L) {
      stop("Server '", srv, "' not in datasources", call. = FALSE)
    }
    for (v in var_map[[srv]]) {
      if (verbose) message("[ds.vertDesc] ", srv, ":", v)

      moments <- DSI::datashield.aggregate(datasources[ci],
        call(name = "dsvertLocalMomentsDS", data_name = data_name,
             variable = v, return_extrema = exact_extrema))[[1]]

      row <- list(
        server = srv, variable = v,
        n = moments$n_total, n_na = moments$n_na,
        mean = moments$mean, sd = moments$sd,
        min = moments$min, max = moments$max,
        range_low = NA_real_, range_high = NA_real_,
        range_method = if (exact_extrema) "exact_extrema" else "moment_sd",
        histogram_buckets = NA_integer_,
        quantile_status = "not_run")

      # Quantiles via histogram interpolation --------------------------
      if (is.na(moments$mean) || moments$n_total < 2L ||
          !is.finite(moments$sd) || moments$sd <= 0) {
        for (qn in q_names) row[[qn]] <- NA_real_
        row$quantile_status <- "insufficient_variation"
      } else {
        if (exact_extrema &&
            is.finite(moments$min) && is.finite(moments$max) &&
            moments$min < moments$max) {
          lo <- moments$min
          hi <- moments$max
        } else {
          lo <- moments$mean - range_sd * moments$sd
          hi <- moments$mean + range_sd * moments$sd
        }
        if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
          for (qn in q_names) row[[qn]] <- NA_real_
          row$quantile_status <- "invalid_range"
        } else {
          privacy_min <- suppressWarnings(as.numeric(getOption(
            "dsvert.min_cell_count",
            getOption("datashield.privacyLevel", 5L)))[1L])
          effective_buckets <- n_buckets
          if (is.numeric(privacy_min) && is.finite(privacy_min) &&
              privacy_min > 0) {
            target_n_per_bucket <- as.numeric(getOption(
              "dsvert.desc_target_n_per_bucket", privacy_min + 1))
            effective_buckets <- min(
              effective_buckets,
              max(2L, floor(as.numeric(moments$n_total) /
                              target_n_per_bucket)))
          }
          row$range_low <- lo
          row$range_high <- hi
          make_edges <- function(bucket_count) {
            if (!exact_extrema && open_ended) {
              finite_grid <- seq(lo, hi, length.out = bucket_count + 1L)
              inner <- finite_grid[-c(1L, length(finite_grid))]
              c(-Inf, inner, Inf)
            } else {
              pad <- (hi - lo) * 1e-9
              seq(lo, hi + pad, length.out = bucket_count + 1L)
            }
          }
          candidate_buckets <- if (isTRUE(getOption(
            "dsvert.desc_adaptive_buckets", TRUE))) {
            seq(effective_buckets, 2L, by = -1L)
          } else {
            effective_buckets
          }
          hist_res <- NULL
          hist_err <- NULL
          edges <- NULL
          for (bucket_count in candidate_buckets) {
            trial_edges <- make_edges(bucket_count)
            trial_res <- tryCatch(
              DSI::datashield.aggregate(datasources[ci],
                call(name = "dsvertHistogramDS", data_name = data_name,
                     variable = v, edges = trial_edges,
                     suppress_small_cells = TRUE,
                     fail_on_small_cells = TRUE))[[1]],
              error = function(e) e)
            if (!inherits(trial_res, "error")) {
              hist_res <- trial_res
              edges <- trial_edges
              row$histogram_buckets <- as.integer(bucket_count)
              break
            }
            hist_err <- conditionMessage(trial_res)
          }
          if (!exact_extrema && open_ended) {
            row$range_method <- "moment_sd_open_ended"
          }
          if (inherits(hist_res, "error")) {
            for (qn in q_names) row[[qn]] <- NA_real_
            row$quantile_status <- conditionMessage(hist_res)
          } else if (is.null(hist_res)) {
            for (qn in q_names) row[[qn]] <- NA_real_
            row$quantile_status <- hist_err %||% "histogram_suppressed"
          } else {
            q_vals <- .dsvert_interp_quantile(
              edges = edges,
              counts = hist_res$counts,
              probs = probs,
              below = as.numeric(hist_res$below),
              above = as.numeric(hist_res$above))
            for (k in seq_along(probs)) row[[q_names[k]]] <- q_vals[k]
            row$quantile_status <- "ran"
          }
        }
      }
      rows[[length(rows) + 1L]] <- row
    }
  }

  # Assemble data frame
  if (length(rows) == 0L) {
    return(data.frame(server = character(0), variable = character(0),
                      n = integer(0), n_na = integer(0),
                      mean = numeric(0), sd = numeric(0),
                      min = numeric(0), max = numeric(0),
                      range_low = numeric(0), range_high = numeric(0),
                      range_method = character(0),
                      histogram_buckets = integer(0),
                      quantile_status = character(0)))
  }
  cols <- names(rows[[1L]])
  out <- lapply(cols, function(cname) {
    vals <- sapply(rows, function(r) r[[cname]])
    unname(vals)
  })
  names(out) <- cols
  out_df <- as.data.frame(out, stringsAsFactors = FALSE)
  class(out_df) <- c("ds.vertDesc", "data.frame")
  attr(out_df, "probs") <- probs
  attr(out_df, "n_buckets") <- n_buckets
  attr(out_df, "range_sd") <- range_sd
  attr(out_df, "exact_extrema") <- exact_extrema
  attr(out_df, "open_ended") <- open_ended
  out_df
}

#' Interpolate a quantile from histogram bucket counts.
#'
#' Pure helper function extracted for unit testing. Given bucket edges, bucket
#' counts (optionally with under/overflow), and target probabilities, returns
#' the linearly-interpolated quantile within the target bucket.
#'
#' @param edges Numeric vector of length K+1 defining the bucket edges.
#' @param counts Integer vector of length K with per-bucket counts.
#' @param probs Numeric vector of target probabilities (0, 1).
#' @param below Count of observations strictly below \code{edges[1]} (default 0).
#' @param above Count of observations strictly above \code{edges[K+1]} (default 0).
#' @return Numeric vector of the same length as probs.
#' @keywords internal
.dsvert_interp_quantile <- function(edges, counts, probs,
                                    below = 0L, above = 0L) {
  counts <- as.numeric(counts)
  K <- length(counts)
  total <- sum(counts) + below + above
  if (total < 1) return(rep(NA_real_, length(probs)))
  # Cumulative starts with `below` observations BEFORE the first bucket.
  cum_before <- c(below, below + cumsum(counts))
  # The total above edges[K+1] contributes `above` extra mass.
  out <- numeric(length(probs))
  for (k in seq_along(probs)) {
    target <- probs[k] * total
    if (target <= below) {
      out[k] <- edges[1L]
      next
    }
    if (target >= below + sum(counts)) {
      out[k] <- edges[K + 1L]
      next
    }
    bucket <- findInterval(target, cum_before, rightmost.closed = TRUE)
    if (bucket < 1L) bucket <- 1L
    if (bucket > K) bucket <- K
    bucket_total <- counts[bucket]
    if (bucket_total == 0) {
      out[k] <- edges[bucket]
    } else {
      left <- edges[bucket]
      right <- edges[bucket + 1L]
      if (!is.finite(left) && is.finite(right)) {
        out[k] <- right
        next
      }
      if (is.finite(left) && !is.finite(right)) {
        out[k] <- left
        next
      }
      if (!is.finite(left) || !is.finite(right)) {
        out[k] <- NA_real_
        next
      }
      frac <- (target - cum_before[bucket]) / bucket_total
      out[k] <- left + frac * (right - left)
    }
  }
  out
}

#' @export
print.ds.vertDesc <- function(x, ...) {
  cat("dsVert descriptive summary (", nrow(x), " variables, ",
      attr(x, "n_buckets"), " requested histogram buckets for quantiles)\n",
      sep = "")
  y <- x
  class(y) <- "data.frame"
  print(y, row.names = FALSE, digits = 5L)
  invisible(x)
}
