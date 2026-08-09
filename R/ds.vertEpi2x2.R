#' Epidemiological effect measures from a protected 2-by-2 table
#'
#' Computes risks, risk difference, risk ratio, odds ratio, attributable
#' fractions among the exposed and in the population, and number needed to
#' benefit/harm from an existing aggregate 2-by-2 table. This function is
#' entirely client-side: it performs
#' no additional DataSHIELD query and therefore does not enlarge the disclosure
#' surface of the table-producing method.
#'
#' @param x A 2-by-2 count matrix, or an object such as
#'   \code{ds.vertChisq} containing an \code{observed} matrix.
#' @param exposed Row identifying the exposed group, as a name or index.
#'   Defaults to the second row.
#' @param event Column identifying the event, as a name or index. Defaults to
#'   the second column.
#' @param level Confidence level for Wald intervals on an ordinary aggregate,
#'   or simultaneous DP-mechanism coverage when `x` is a DP-aware
#'   `ds.vertChisq` result.
#' @param zero_correction One of \code{"none"} (default), \code{"if_zero"},
#'   or \code{"always"}. The latter two explicitly request a 0.5
#'   Haldane--Anscombe correction for ratio estimates and report that choice
#'   in the result.
#'
#' @return For an ordinary aggregate, a `ds.vertEpi2x2` object with standard
#'   epidemiological estimates and Wald intervals. For a DP-aware
#'   `ds.vertChisq` result, a `ds.vertDPEpi2x2` object with simultaneous
#'   mechanism-noise regions and no ordinary sampling interval.
#'
#' @details A DP-aware `ds.vertChisq` object is routed through its validated
#'   source capsule. Ordinary Wald formulae and continuity corrections are not
#'   applied to DP-noised cells, and no additional server call is made.
#' @export
ds.vertEpi2x2 <- function(x, exposed = 2L, event = 2L, level = 0.95,
                          zero_correction = c("none", "if_zero", "always")) {
  zero_correction <- match.arg(zero_correction)
  carries_dp_provenance <- inherits(x, "ds.vertChisq") &&
    (!is.null(x$calibration) || !is.null(x$source_dp_release))
  is_dp_chisq <- carries_dp_provenance &&
    identical(x$calibration, .DSVERT_DP_CHISQ_BOOTSTRAP_VERSION) &&
    inherits(x$source_dp_release, "ds.vertDPContingency")
  if (carries_dp_provenance && !is_dp_chisq) {
    stop("x carries incomplete or invalid DP contingency provenance",
         call. = FALSE)
  }
  if (is_dp_chisq) {
    if (!identical(zero_correction, "none")) {
      stop("zero_correction is not available for a DP contingency release; ",
           "use its certified mechanism regions", call. = FALSE)
    }
    return(ds.vertDPEpi2x2(
      x$source_dp_release, exposed = exposed, event = event, level = level))
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }

  is_checked_object <- is.list(x) && !is.null(x$observed)
  disclosure_guard <- if (is_checked_object) x$disclosure_guard else NULL
  if (!is.null(disclosure_guard$passed) &&
      identical(disclosure_guard$passed, FALSE)) {
    stop("The input table did not pass its disclosure guard", call. = FALSE)
  }
  observed <- if (is_checked_object) x$observed else x
  observed <- as.matrix(observed)
  if (!identical(dim(observed), c(2L, 2L))) {
    stop("x must contain exactly 2 x 2 aggregate counts", call. = FALSE)
  }
  storage.mode(observed) <- "double"
  if (any(!is.finite(observed)) || any(observed < 0)) {
    stop("counts must be finite and non-negative", call. = FALSE)
  }
  if (any(abs(observed - round(observed)) > sqrt(.Machine$double.eps))) {
    stop("counts must be whole-number aggregate counts", call. = FALSE)
  }

  select_dimension <- function(value, labels, what) {
    if (is.character(value) && length(value) == 1L) {
      if (is.null(labels) || !value %in% labels) {
        stop("Unknown ", what, " level: '", value, "'", call. = FALSE)
      }
      return(match(value, labels))
    }
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != as.integer(value) || !value %in% 1:2) {
      stop(what, " must identify one row/column by name or index",
           call. = FALSE)
    }
    as.integer(value)
  }

  exposed_idx <- select_dimension(exposed, rownames(observed), "exposed")
  event_idx <- select_dimension(event, colnames(observed), "event")
  unexposed_idx <- setdiff(1:2, exposed_idx)
  nonevent_idx <- setdiff(1:2, event_idx)

  # Conventional layout after orientation:
  #                event  no event
  # exposed           a       b
  # unexposed         c       d
  a <- observed[exposed_idx, event_idx]
  b <- observed[exposed_idx, nonevent_idx]
  c_ <- observed[unexposed_idx, event_idx]
  d <- observed[unexposed_idx, nonevent_idx]
  oriented <- matrix(
    c(a, c_, b, d), nrow = 2L,
    dimnames = list(exposure = c("exposed", "unexposed"),
                    outcome = c("event", "no_event")))

  n_exposed <- a + b
  n_unexposed <- c_ + d
  if (n_exposed <= 0 || n_unexposed <= 0) {
    stop("Both exposure groups must contain observations", call. = FALSE)
  }
  risk_exposed <- a / n_exposed
  risk_unexposed <- c_ / n_unexposed
  population_risk <- (a + c_) / sum(observed)
  rd <- risk_exposed - risk_unexposed
  z <- stats::qnorm((1 + level) / 2)
  se_rd <- sqrt(risk_exposed * (1 - risk_exposed) / n_exposed +
                  risk_unexposed * (1 - risk_unexposed) / n_unexposed)

  apply_correction <- identical(zero_correction, "always") ||
    (identical(zero_correction, "if_zero") && any(oriented == 0))
  correction <- if (apply_correction) 0.5 else 0
  ac <- a + correction
  bc <- b + correction
  cc <- c_ + correction
  dc <- d + correction
  ne_c <- ac + bc
  nu_c <- cc + dc

  rr <- (ac / ne_c) / (cc / nu_c)
  or <- (ac * dc) / (bc * cc)
  se_log_rr <- sqrt(1 / ac - 1 / ne_c + 1 / cc - 1 / nu_c)
  se_log_or <- sqrt(1 / ac + 1 / bc + 1 / cc + 1 / dc)
  ratio_interval <- function(estimate, se) {
    if (!is.finite(estimate) || estimate <= 0 || !is.finite(se)) {
      return(c(lower = NA_real_, upper = NA_real_))
    }
    exp(log(estimate) + c(-1, 1) * z * se)
  }
  rr_ci <- ratio_interval(rr, se_log_rr)
  or_ci <- ratio_interval(or, se_log_or)
  rd_ci <- pmax(-1, pmin(1, rd + c(-1, 1) * z * se_rd))

  af <- if (is.finite(rr) && rr != 0) (rr - 1) / rr else NA_real_
  paf <- if (is.finite(population_risk) && population_risk > 0) {
    (population_risk - risk_unexposed) / population_risk
  } else {
    NA_real_
  }
  nnt <- if (rd == 0) Inf else 1 / rd
  number_needed <- if (rd == 0) Inf else 1 / abs(rd)
  number_needed_direction <- if (rd < 0) {
    "benefit"
  } else if (rd > 0) {
    "harm"
  } else {
    "none"
  }
  notes <- character(0)
  if (apply_correction) {
    notes <- paste0(
      "A 0.5 Haldane-Anscombe correction was applied to all four cells ",
      "for ratio estimates because zero_correction='", zero_correction, "'.")
  }

  out <- list(
    observed = observed,
    oriented = oriented,
    n = sum(observed),
    risks = c(exposed = risk_exposed, unexposed = risk_unexposed,
              population = population_risk),
    risk_difference = list(
      estimate = rd, std_error = se_rd,
      lower = unname(rd_ci[1L]), upper = unname(rd_ci[2L])),
    risk_ratio = list(
      estimate = rr, std_error_log = se_log_rr,
      lower = unname(rr_ci[1L]), upper = unname(rr_ci[2L])),
    odds_ratio = list(
      estimate = or, std_error_log = se_log_or,
      lower = unname(or_ci[1L]), upper = unname(or_ci[2L])),
    attributable_fraction_exposed = list(estimate = af),
    population_attributable_fraction = list(estimate = paf),
    nnt = list(estimate = nnt,
               interpretation = if (rd > 0) "harm" else if (rd < 0) "benefit" else "none"),
    number_needed = list(
      estimate = number_needed,
      direction = number_needed_direction,
      definition = "absolute reciprocal of the finite-snapshot risk difference"),
    level = level,
    correction = correction,
    notes = notes,
    disclosure_guard = disclosure_guard,
    input_provenance = if (is_checked_object) {
      "disclosure_checked_table"
    } else {
      "caller_supplied_aggregate"
    },
    method_status = "promoted",
    quality = .dsvert_quality(
      status = if (any(oriented == 0) && !apply_correction) {
        "degraded"
      } else {
        "ok"
      },
      warnings = if (any(oriented == 0) && !apply_correction) {
        paste0("At least one ratio estimate is on a zero-cell boundary; ",
               "no continuity correction was applied.")
      } else {
        character(0)
      },
      metrics = list(
        zero_cell_boundary = any(oriented == 0) && !apply_correction)))
  class(out) <- c("ds.vertEpi2x2", "list")
  out
}

#' @export
print.ds.vertEpi2x2 <- function(x, ...) {
  pct <- function(v) sprintf("%.2f%%", 100 * v)
  cat("dsVert epidemiological 2x2 measures\n")
  cat(sprintf("  Risks: exposed %s; unexposed %s\n",
              pct(x$risks[["exposed"]]), pct(x$risks[["unexposed"]])))
  cat(sprintf("  RD = %.6g  (%.1f%% CI %.6g to %.6g)\n",
              x$risk_difference$estimate, 100 * x$level,
              x$risk_difference$lower, x$risk_difference$upper))
  cat(sprintf("  RR = %.6g  (%.1f%% CI %.6g to %.6g)\n",
              x$risk_ratio$estimate, 100 * x$level,
              x$risk_ratio$lower, x$risk_ratio$upper))
  cat(sprintf("  OR = %.6g  (%.1f%% CI %.6g to %.6g)\n",
              x$odds_ratio$estimate, 100 * x$level,
              x$odds_ratio$lower, x$odds_ratio$upper))
  if (length(x$notes)) cat("  Note: ", x$notes, "\n", sep = "")
  invisible(x)
}
