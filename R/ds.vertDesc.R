#' @title Disclosure-safe descriptive statistics compatibility adapter
#' @description Return the historical \code{ds.vertDesc} data-frame shape from
#'   one custodian-owned, sticky \code{ds.vertDPDescribe} capsule artifact.
#'   Counts, moments and quantiles are differentially private. Observed extrema
#'   and data-adaptive histogram queries are never requested.
#'
#' @param data_name Character. Name of the aligned data frame held on each
#'   server.
#' @param variables Optional character vector, or a named list mapping the
#'   artifact-owning server to variables, used only to filter the variables
#'   already present in the signed artifact. \code{NULL} returns every signed
#'   variable. A variable absent from the artifact is an error.
#' @param probs Numeric vector of quantile probabilities to report.
#'   Defaults to the usual quartiles (0.25, 0.5, 0.75).
#' @param n_buckets Legacy compatibility argument. If explicitly supplied, it
#'   must equal the fixed signed grid length of every selected variable. It
#'   never changes the server workload.
#' @param range_sd Legacy data-adaptive range argument. Its omitted default is
#'   ignored. Explicit values are rejected because ranges are fixed by the
#'   custodian-owned artifact.
#' @param exact_extrema Logical legacy argument. \code{TRUE} is always rejected:
#'   exact observed minima and maxima are not disclosure-safe releases.
#' @param open_ended Legacy adaptive-tail argument. Its omitted default is
#'   ignored. Explicit values are rejected because the signed public bounds,
#'   grids and invalid bin define the histogram semantics.
#' @param verbose Logical. Print per-variable progress when TRUE.
#' @param datasources DataSHIELD connections. If \code{NULL}, auto-detected.
#' @param analysis_id Custodian-owned describe specification id. It must name
#'   an existing artifact in the signed capsule. The client never discovers or
#'   chooses an id remotely.
#'
#' @return A data frame with one row per variable containing columns:
#'   \itemize{
#'     \item \code{server}: which server holds the variable
#'     \item \code{variable}: column name
#'     \item \code{n}: DP-noisy effective count
#'     \item \code{n_na}: DP-noisy invalid/missing-unit count
#'     \item \code{mean}, \code{sd}: bounded DP mean and square root of the
#'       DP population central second moment (not the usual sample SD)
#'     \item \code{min}, \code{max}: always \code{NA}
#'     \item \code{range_low}, \code{range_high}: signed public bounds
#'     \item \code{quantile_status}: fixed-grid DP post-processing status
#'     \item one column per requested quantile (named \code{q25}, \code{q50},
#'       \code{q75} etc. by default)
#'   }
#'
#' @details
#' The function makes one call to \code{ds.vertDPDescribe}. Quantile
#' probabilities and variable selection are then pure client-side
#' post-processing of that single immutable release. The returned object has
#' class \code{ds.vertDesc} for compatibility and carries \code{dp_release},
#' \code{dp_descriptives}, and \code{dp_quantile_bands} attributes with the
#' formal release and mechanism-uncertainty metadata. These regions exclude
#' sampling uncertainty. In particular, \code{n_na} is the noised fixed invalid
#' bin rather than an exact row-level NA count, and quantiles are fixed-grid
#' upper endpoints rather than interpolation from analyst-selected bins.
#'
#' @examples
#' \dontrun{
#' conns <- DSI::datashield.login(logindata)
#' ds.vertDesc("DA", variables = c("age", "bmi", "glu", "bp"),
#'             probs = c(0.25, 0.5, 0.75, 0.9),
#'             analysis_id = "baseline_describe_v1")
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
                        datasources = NULL,
                        analysis_id = NULL) {
  buckets_supplied <- !missing(n_buckets)
  range_supplied <- !missing(range_sd)
  tails_supplied <- !missing(open_ended)

  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) || !nzchar(data_name)) {
    stop("data_name must be a non-empty string", call. = FALSE)
  }
  if (!is.character(analysis_id) || length(analysis_id) != 1L ||
      is.na(analysis_id) || !nzchar(analysis_id)) {
    stop(paste(
      "analysis_id is required for ds.vertDesc.",
      "Ask the server administrator to provision a custodian-owned describe",
      "specification in the signed DP capsule, then pass that id explicitly."),
      call. = FALSE)
  }
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("probs must contain finite probabilities strictly inside (0,1)",
         call. = FALSE)
  }
  probs <- sort(unique(as.numeric(probs)), method = "radix")
  q_names <- sprintf("q%02d", as.integer(round(100 * probs)))
  if (anyDuplicated(q_names)) {
    stop("probs produce duplicate legacy quantile column names", call. = FALSE)
  }
  if (!is.logical(exact_extrema) || length(exact_extrema) != 1L ||
      is.na(exact_extrema)) {
    stop("exact_extrema must be TRUE or FALSE", call. = FALSE)
  }
  if (isTRUE(exact_extrema)) {
    stop(paste(
      "exact_extrema = TRUE is unavailable because observed minima and",
      "maxima are not disclosure-safe under repeated queries.",
      "Use the signed public bounds and fixed-grid DP quantiles instead."),
      call. = FALSE)
  }
  if (range_supplied) {
    if (!is.numeric(range_sd) || length(range_sd) != 1L ||
        is.na(range_sd) || !is.finite(range_sd) || range_sd <= 0) {
      stop("range_sd must be a positive finite number", call. = FALSE)
    }
    stop(paste(
      "range_sd cannot redefine a disclosure-safe signed workload.",
      "Ask the server administrator for a custodian-owned analysis_id whose",
      "fixed public bounds and histogram grid match the intended analysis."),
      call. = FALSE)
  }
  if (tails_supplied) {
    if (!is.logical(open_ended) || length(open_ended) != 1L ||
        is.na(open_ended)) {
      stop("open_ended must be TRUE or FALSE", call. = FALSE)
    }
    stop(paste(
      "open_ended cannot redefine a disclosure-safe signed workload.",
      "The custodian-owned analysis_id fixes public bounds, bins and the",
      "invalid-bin policy."), call. = FALSE)
  }
  if (buckets_supplied &&
      (!is.numeric(n_buckets) || length(n_buckets) != 1L ||
       is.na(n_buckets) || !is.finite(n_buckets) || n_buckets < 2 ||
       n_buckets != floor(n_buckets))) {
    stop("n_buckets must be one integer >= 2", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(variables)) {
    if (is.character(variables)) {
      if (!length(variables) || anyNA(variables) || any(!nzchar(variables)) ||
          anyDuplicated(variables)) {
        stop("variables must contain unique non-empty names", call. = FALSE)
      }
    } else if (is.list(variables) && !is.null(names(variables)) &&
               length(variables) && !anyNA(names(variables)) &&
               all(nzchar(names(variables))) && !anyDuplicated(names(variables)) &&
               all(vapply(variables, function(value) {
                 is.character(value) && length(value) && !anyNA(value) &&
                   all(nzchar(value)) && !anyDuplicated(value)
               }, logical(1L)))) {
      # Validated after the signed artifact identifies its single owner.
    } else {
      stop(paste(
        "variables must be NULL, unique variable names, or a named list of",
        "unique variable names for one artifact-owning server"), call. = FALSE)
    }
  }

  if (isTRUE(verbose)) {
    message("[ds.vertDesc] reading signed DP describe artifact '",
            analysis_id, "'")
  }
  release <- ds.vertDPDescribe(
    data_name = data_name, analysis_id = analysis_id, probs = probs,
    datasources = datasources)
  if (!is.list(release) || !isTRUE(release$released)) {
    stop(paste(
      "The signed DP describe artifact did not produce a release.",
      "Ask the server administrator to verify the custodian-owned",
      "analysis_id and its representable fixed privacy allocation."),
      call. = FALSE)
  }

  available <- release$variables
  selected <- available
  if (is.character(variables)) {
    selected <- variables
  } else if (is.list(variables)) {
    other_servers <- setdiff(names(variables), release$server)
    if (length(other_servers)) {
      stop("analysis_id is owned by server '", release$server,
           "'; variables cannot select another server", call. = FALSE)
    }
    selected <- variables[[release$server]]
  }
  absent <- setdiff(selected, available)
  if (length(absent)) {
    stop(paste0(
      "Variable(s) not included in signed describe artifact '", analysis_id,
      "': ", paste(absent, collapse = ", "), ". Ask the server administrator ",
      "for a custodian-owned analysis_id containing the required fixed ",
      "variables and grids."), call. = FALSE)
  }

  signed_buckets <- stats::setNames(
    as.integer(release$grid_lengths), available)
  selected_buckets <- signed_buckets[selected]
  if (buckets_supplied && any(selected_buckets != as.integer(n_buckets))) {
    details <- paste(paste0(names(selected_buckets), "=", selected_buckets),
                     collapse = ", ")
    stop(paste0(
      "n_buckets does not match the signed fixed grid (", details,
      "). It cannot change or reroll the capsule workload; ask the server ",
      "administrator for a custodian-owned analysis_id with the intended grid."),
      call. = FALSE)
  }

  descriptions <- release$descriptives[
    match(selected, release$descriptives$variable), , drop = FALSE]
  quantiles <- release$quantiles[
    release$quantiles$variable %in% selected, , drop = FALSE]
  rows <- lapply(seq_along(selected), function(index) {
    variable <- selected[[index]]
    description <- descriptions[index, , drop = FALSE]
    variable_quantiles <- quantiles[quantiles$variable == variable, ,
                                    drop = FALSE]
    positions <- match(probs, variable_quantiles$probability)
    if (anyNA(positions)) {
      stop("The DP describe artifact omitted a requested quantile",
           call. = FALSE)
    }
    variable_quantiles <- variable_quantiles[positions, , drop = FALSE]
    row <- list(
      server = release$server,
      variable = variable,
      n = description$n_dp,
      n_na = description$invalid_dp,
      mean = description$mean,
      sd = description$sd,
      min = NA_real_,
      max = NA_real_,
      range_low = description$lower_bound,
      range_high = description$upper_bound,
      range_method = "custodian_signed_fixed_grid",
      histogram_buckets = unname(selected_buckets[[variable]]),
      quantile_status = if (all(variable_quantiles$status == "ok")) {
        "ran_dp_fixed_grid"
      } else {
        paste(unique(variable_quantiles$status), collapse = ";")
      })
    for (quantile_index in seq_along(q_names)) {
      row[[q_names[[quantile_index]]]] <-
        variable_quantiles$estimate[[quantile_index]]
    }
    as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  })
  out_df <- do.call(rbind, rows)
  rownames(out_df) <- NULL
  class(out_df) <- c("ds.vertDesc", "data.frame")
  attr(out_df, "probs") <- probs
  attr(out_df, "n_buckets") <- if (length(unique(selected_buckets)) == 1L) {
    unname(selected_buckets[[1L]])
  } else selected_buckets
  attr(out_df, "signed_histogram_buckets") <- selected_buckets
  attr(out_df, "range_sd") <- NA_real_
  attr(out_df, "exact_extrema") <- FALSE
  attr(out_df, "open_ended") <- NA
  attr(out_df, "analysis_id") <- release$analysis_id
  attr(out_df, "compatibility_semantics") <- list(
    n = "DP-noisy effective privacy-unit count",
    n_na = "DP-noisy fixed invalid-bin count, not an exact row NA count",
    sd = paste(
      "square root of the bounded DP population central second moment;",
      "not the usual sample standard deviation"),
    min_max = "observed extrema are not released",
    quantiles = "upper endpoints of custodian-signed fixed histogram bins")
  attr(out_df, "dp_release") <- list(
    formal_dp = TRUE,
    source = "ds.vertDPDescribe",
    analysis_id = release$analysis_id,
    analysis_version = release$analysis_version,
    server = release$server,
    capsule_id = release$capsule_id,
    final_vector_root = release$final_vector_root,
    coordinate_order_sha256 = release$coordinate_order_sha256,
    privacy_epoch = release$privacy_epoch,
    noise_key_id = release$noise_key_id,
    epsilon = release$epsilon,
    delta = release$delta,
    implementation_delta = release$implementation_delta,
    adjacency = release$adjacency,
    mechanism = release$mechanism,
    sampler = release$sampler,
    sticky_noise = release$sticky_noise,
    uncertainty_scope = release$uncertainty_scope,
    histogram_semantics = release$histogram_semantics,
    unit_collapse = release$unit_collapse,
    count_definition = release$count_definition,
    invalid_unit_rule = release$invalid_unit_rule,
    quantization = release$quantization,
    postprocessing = release$postprocessing,
    artifact_variables = release$variables,
    selected_variables = selected,
    quantile_band_confidence = release$quantile_band_confidence,
    quantile_band_scope = release$quantile_band_scope,
    moment_region_confidence = release$moment_region_confidence,
    moment_region_method = release$moment_region_method,
    moment_region_scope = release$moment_region_scope,
    statistical_inference = release$statistical_inference)
  attr(out_df, "dp_descriptives") <- descriptions
  quantile_keys <- unlist(lapply(selected, function(variable) {
    paste(variable, probs, sep = "\r")
  }), use.names = FALSE)
  attr(out_df, "dp_quantile_bands") <- quantiles[
    match(quantile_keys,
          paste(quantiles$variable, quantiles$probability, sep = "\r")),
    , drop = FALSE]
  out_df
}

#' @title Interpolate a quantile from histogram bucket counts
#' @description Pure helper function extracted for unit testing. Given bucket
#'   edges, bucket counts (optionally with under/overflow), and target
#'   probabilities, returns the linearly-interpolated quantile within the target
#'   bucket.
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
  release <- attr(x, "dp_release")
  cat("dsVert DP descriptive summary (", nrow(x),
      " variables; custodian-signed fixed grids)\n", sep = "")
  cat("analysis_id: ", attr(x, "analysis_id"), " | epsilon: ",
      format(release$epsilon), " | delta: ", format(release$delta), "\n",
      sep = "")
  y <- x
  class(y) <- "data.frame"
  print(y, row.names = FALSE, digits = 5L)
  invisible(x)
}
