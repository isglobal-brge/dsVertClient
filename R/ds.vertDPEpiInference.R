.dsvert_dp_cp_interval <- function(events, nonevents, alpha) {
  values <- c(events = events, nonevents = nonevents, alpha = alpha)
  if (any(!is.finite(values)) || events < 0 || nonevents < 0 ||
      events != floor(events) || nonevents != floor(nonevents) ||
      events + nonevents > 2^53 - 1 || alpha <= 0 || alpha >= 1) {
    stop("Invalid exact-binomial confidence-interval inputs", call. = FALSE)
  }
  lower <- if (events == 0) {
    0
  } else {
    stats::qbeta(alpha / 2, events, nonevents + 1)
  }
  upper <- if (nonevents == 0) {
    1
  } else {
    stats::qbeta(1 - alpha / 2, events + 1, nonevents)
  }
  if (!is.finite(lower) || !is.finite(upper) ||
      lower < 0 || upper > 1 || lower > upper) {
    stop("Exact-binomial confidence limits are not representable",
         call. = FALSE)
  }
  c(lower = lower, upper = upper)
}

.dsvert_dp_integer_count_range <- function(lower, upper, allow_empty = FALSE) {
  if (!is.numeric(lower) || !is.numeric(upper) ||
      length(lower) != 1L || length(upper) != 1L ||
      is.na(lower) || is.na(upper) || !is.finite(lower) ||
      !is.finite(upper) || lower < 0 || lower > upper ||
      upper > 2^53 - 1) {
    stop("The DP count box cannot be converted to exact integer counts",
         call. = FALSE)
  }
  lower_margin <- 128 * .Machine$double.eps * max(1, abs(lower))
  upper_margin <- 128 * .Machine$double.eps * max(1, abs(upper))
  result <- c(
    lower = ceiling(max(0, lower - lower_margin)),
    upper = min(2^53 - 1, floor(upper + upper_margin)))
  if (result[["lower"]] > result[["upper"]]) {
    if (isTRUE(allow_empty)) return(c(lower = NA_real_, upper = NA_real_))
    stop("The DP count box contains no representable integer count",
         call. = FALSE)
  }
  result
}

.dsvert_dp_integer_count_box <- function(lower, upper) {
  if (!is.matrix(lower) || !is.matrix(upper) ||
      !identical(dim(lower), dim(upper)) ||
      !identical(dimnames(lower), dimnames(upper))) {
    stop("The DP count bounds must be conformable named matrices",
         call. = FALSE)
  }
  ranges <- vapply(seq_along(lower), function(index) {
    .dsvert_dp_integer_count_range(
      lower[[index]], upper[[index]], allow_empty = TRUE)
  }, numeric(2L))
  integer_lower <- matrix(
    ranges["lower", ], nrow = nrow(lower), ncol = ncol(lower),
    dimnames = dimnames(lower))
  integer_upper <- matrix(
    ranges["upper", ], nrow = nrow(upper), ncol = ncol(upper),
    dimnames = dimnames(upper))
  list(
    lower = integer_lower,
    upper = integer_upper,
    empty_coordinates = which(
      is.na(integer_lower) | is.na(integer_upper), arr.ind = TRUE),
    empty = anyNA(integer_lower) || anyNA(integer_upper))
}

.dsvert_dp_cp_union_over_box <- function(
    event_lower, event_upper, nonevent_lower, nonevent_upper, alpha) {
  event <- .dsvert_dp_integer_count_range(event_lower, event_upper)
  nonevent <- .dsvert_dp_integer_count_range(
    nonevent_lower, nonevent_upper)
  lower <- .dsvert_dp_cp_interval(
    event[["lower"]], nonevent[["upper"]], alpha)[["lower"]]
  upper <- .dsvert_dp_cp_interval(
    event[["upper"]], nonevent[["lower"]], alpha)[["upper"]]
  c(lower = lower, upper = upper)
}

.dsvert_dp_probability_ratio_region <- function(numerator, denominator) {
  valid <- function(value) {
    is.numeric(value) && length(value) == 2L &&
      identical(names(value), c("lower", "upper")) &&
      !anyNA(value) && all(is.finite(value)) &&
      value[["lower"]] >= 0 && value[["upper"]] <= 1 &&
      value[["lower"]] <= value[["upper"]]
  }
  if (!valid(numerator) || !valid(denominator)) {
    stop("Invalid probability regions", call. = FALSE)
  }
  c(
    lower = if (denominator[["upper"]] > 0) {
      numerator[["lower"]] / denominator[["upper"]]
    } else if (numerator[["lower"]] > 0) {
      Inf
    } else {
      0
    },
    upper = if (denominator[["lower"]] > 0) {
      numerator[["upper"]] / denominator[["lower"]]
    } else if (numerator[["upper"]] > 0) {
      Inf
    } else {
      0
    })
}

.dsvert_dp_odds_region <- function(probability) {
  if (!is.numeric(probability) || length(probability) != 2L ||
      !identical(names(probability), c("lower", "upper")) ||
      anyNA(probability) || any(!is.finite(probability)) ||
      probability[["lower"]] < 0 || probability[["upper"]] > 1 ||
      probability[["lower"]] > probability[["upper"]]) {
    stop("Invalid probability region", call. = FALSE)
  }
  c(
    lower = if (probability[["lower"]] >= 1) {
      Inf
    } else {
      probability[["lower"]] / (1 - probability[["lower"]])
    },
    upper = if (probability[["upper"]] >= 1) {
      Inf
    } else {
      probability[["upper"]] / (1 - probability[["upper"]])
    })
}

.dsvert_dp_nonnegative_ratio_region <- function(numerator, denominator) {
  valid <- function(value) {
    is.numeric(value) && length(value) == 2L &&
      identical(names(value), c("lower", "upper")) &&
      !anyNA(value) && all(value >= 0) &&
      value[["lower"]] <= value[["upper"]]
  }
  if (!valid(numerator) || !valid(denominator)) {
    stop("Invalid non-negative ratio regions", call. = FALSE)
  }
  c(
    lower = if (is.infinite(denominator[["upper"]])) {
      0
    } else if (denominator[["upper"]] > 0) {
      numerator[["lower"]] / denominator[["upper"]]
    } else if (numerator[["lower"]] > 0) {
      Inf
    } else {
      0
    },
    upper = if (denominator[["lower"]] > 0) {
      numerator[["upper"]] / denominator[["lower"]]
    } else if (numerator[["upper"]] > 0) {
      Inf
    } else {
      0
    })
}

.dsvert_dp_inference_region_type <- function(interval) {
  if (!is.numeric(interval) || length(interval) != 2L ||
      !identical(names(interval), c("lower", "upper")) ||
      anyNA(interval) || interval[["lower"]] > interval[["upper"]]) {
    stop("Combined uncertainty interval construction failed", call. = FALSE)
  }
  if (all(is.finite(interval))) return("finite")
  if (is.finite(interval[["lower"]]) &&
      is.infinite(interval[["upper"]]) && interval[["upper"]] > 0) {
    return(if (interval[["lower"]] == 0) {
      "unbounded_nonnegative"
    } else {
      "unbounded_above"
    })
  }
  "unbounded"
}

#' Epidemiological 2-by-2 inference with DP and sampling uncertainty
#'
#' Builds conservative simultaneous confidence regions for two binomial risks
#' and their epidemiological contrasts from one validated DP 2-by-2 release.
#' It first forms a simultaneous DP count box, then unions exact
#' Clopper--Pearson intervals over every integer table compatible with that
#' box. Bonferroni allocation across the DP event and three sampling intervals
#' gives the reported joint lower coverage bound. No noisy cell is treated as
#' an exact confidential count.
#'
#' The sampling statement assumes independent, identically sampled privacy
#' units from a joint exposure/outcome population. Group-risk intervals are
#' conditionally binomial given the confidential group sizes; the population
#' risk interval is marginally binomial. Effects are deterministic
#' post-processing of their simultaneous base-risk regions. The construction
#' is intentionally conservative and provides no p-value.
#'
#' @param x A released `ds.vertDPContingency` with a 2-by-2 table.
#' @param exposed,event Exposed row and event column by name or index.
#' @param level Desired lower bound for joint sampling-plus-mechanism coverage.
#' @param mechanism_alpha_share Fraction of total non-coverage allocated to
#'   the DP mechanism event. The remainder is divided equally across exposed,
#'   unexposed, and population exact-binomial intervals.
#' @return A `ds.vertDPEpi2x2Inference` object. It makes no server call and
#'   consumes no additional privacy.
#' @references Clopper, C. J. and Pearson, E. S. (1934). The use of
#'   confidence or fiducial limits illustrated in the case of the binomial.
#'   Biometrika 26(4), 404--413. \doi{10.1093/biomet/26.4.404}.
#' @export
ds.vertDPEpi2x2Inference <- function(
    x, exposed = 2L, event = 2L, level = 0.95,
    mechanism_alpha_share = 0.5) {
  x <- .dsvert_dp_table_contract(x)
  if (!identical(dim(x$table), c(2L, 2L))) {
    stop("x must contain exactly a 2-by-2 DP table", call. = FALSE)
  }
  values <- c(level = level, mechanism_alpha_share = mechanism_alpha_share)
  if (any(!is.finite(values)) || level <= 0 || level >= 1 ||
      mechanism_alpha_share <= 0 || mechanism_alpha_share >= 1) {
    stop("level and mechanism_alpha_share must be finite numbers in (0, 1)",
         call. = FALSE)
  }
  exposed <- .dsvert_dp_dimension_index(
    exposed, rownames(x$table), 2L, "exposed")
  event <- .dsvert_dp_dimension_index(
    event, colnames(x$table), 2L, "event")
  unexposed <- setdiff(1:2, exposed)
  nonevent <- setdiff(1:2, event)

  alpha <- 1 - level
  mechanism_alpha <- alpha * mechanism_alpha_share
  sampling_alpha <- alpha - mechanism_alpha
  mechanism_level <- 1 - mechanism_alpha
  base_sampling_alpha <- sampling_alpha / 3
  mechanism <- ds.vertDPEpi2x2(
    x, exposed = exposed, event = event, level = mechanism_level)
  lower <- mechanism$count_lower
  upper <- mechanism$count_upper
  integer_box <- .dsvert_dp_integer_count_box(lower, upper)
  if (integer_box$empty) {
    risk_exposed <- risk_unexposed <- population_risk <-
      c(lower = 0, upper = 1)
  } else {
    risk_exposed <- .dsvert_dp_cp_union_over_box(
      integer_box$lower[exposed, event],
      integer_box$upper[exposed, event],
      integer_box$lower[exposed, nonevent],
      integer_box$upper[exposed, nonevent], base_sampling_alpha)
    risk_unexposed <- .dsvert_dp_cp_union_over_box(
      integer_box$lower[unexposed, event],
      integer_box$upper[unexposed, event],
      integer_box$lower[unexposed, nonevent],
      integer_box$upper[unexposed, nonevent], base_sampling_alpha)
    population_risk <- .dsvert_dp_cp_union_over_box(
      sum(integer_box$lower[, event]),
      sum(integer_box$upper[, event]),
      sum(integer_box$lower[, nonevent]),
      sum(integer_box$upper[, nonevent]), base_sampling_alpha)
  }

  risk_difference <- c(
    lower = risk_exposed[["lower"]] - risk_unexposed[["upper"]],
    upper = risk_exposed[["upper"]] - risk_unexposed[["lower"]])
  risk_ratio <- .dsvert_dp_probability_ratio_region(
    risk_exposed, risk_unexposed)
  odds_ratio <- .dsvert_dp_nonnegative_ratio_region(
    .dsvert_dp_odds_region(risk_exposed),
    .dsvert_dp_odds_region(risk_unexposed))
  attributable_fraction <- .dsvert_dp_attributable_fraction_exposed_region(
    risk_exposed, risk_unexposed)
  population_attributable_fraction <-
    .dsvert_dp_attributable_fraction_exposed_region(
      population_risk, risk_unexposed)
  regions <- list(
    risk_exposed = risk_exposed,
    risk_unexposed = risk_unexposed,
    population_risk = population_risk,
    risk_difference = risk_difference,
    risk_ratio = risk_ratio,
    odds_ratio = odds_ratio,
    attributable_fraction_exposed = attributable_fraction$interval,
    population_attributable_fraction =
      population_attributable_fraction$interval)
  region_types <- vapply(
    regions, .dsvert_dp_inference_region_type, character(1L))
  number_needed <- .dsvert_dp_number_needed_region(risk_difference)

  result <- list(
    status = mechanism$status,
    combined_region_status = if (integer_box$empty) {
      "vacuous_empty_integer_mechanism_box"
    } else {
      "ok"
    },
    point_estimates = mechanism$point_estimates,
    point_status = mechanism$point_status,
    combined_regions = regions,
    combined_region_types = region_types,
    mechanism_regions = mechanism$mechanism_regions,
    confidential_count_integer_box = integer_box,
    number_needed = list(
      combined_regions = number_needed$regions,
      combined_region_includes_infinite = number_needed$includes_infinite,
      combined_region_possible_directions = number_needed$possible_directions),
    level = level,
    coverage_lower_bound = level,
    mechanism_level = mechanism_level,
    sampling_familywise_level = 1 - sampling_alpha,
    base_sampling_interval_level = 1 - base_sampling_alpha,
    alpha_allocation = c(
      total = alpha, mechanism = mechanism_alpha,
      sampling_familywise = sampling_alpha,
      each_of_three_sampling_intervals = base_sampling_alpha),
    coverage_method = paste(
      "union bound over the signed DP simultaneous count-box event and",
      "three Bonferroni exact Clopper-Pearson intervals; effect regions",
      "are deterministic interval post-processing; an integer-empty",
      "mechanism box returns the full parameter space"),
    sampling_model = paste(
      "iid privacy units from a joint exposure/outcome population;",
      "group risks conditionally binomial given confidential group sizes",
      "and population risk marginally binomial"),
    uncertainty_scope =
      "joint DP-mechanism and superpopulation sampling uncertainty",
    inferential_scope = paste(
      "Conservative confidence regions for two binomial risks and derived",
      "epidemiological effects; no hypothesis test or p-value"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPEpi2x2Inference", "list")
  result
}

#' @export
print.ds.vertDPEpi2x2Inference <- function(x, ...) {
  cat("dsVert DP 2x2 sampling inference:", x$status, "\n")
  cat("Conservative joint coverage >= ", format(100 * x$level), "%",
      " (", x$combined_region_status, ")\n", sep = "")
  print(do.call(rbind, x$combined_regions), ...)
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}
