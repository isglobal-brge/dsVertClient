.dsvert_dp_causal_design <- function(
    x, strata, treatment, treated, standard_weights, event) {
  rows <- nrow(x$table)
  strata <- as.character(strata)
  if (length(strata) != rows || anyNA(strata) ||
      any(!nzchar(strata)) || length(unique(strata)) < 1L) {
    stop("strata must give one non-empty public stratum per table row",
         call. = FALSE)
  }
  treatment <- as.character(treatment)
  treated <- as.character(treated)
  if (length(treatment) != rows || anyNA(treatment) ||
      any(!nzchar(treatment)) || length(unique(treatment)) != 2L ||
      length(treated) != 1L || is.na(treated) || !nzchar(treated) ||
      !treated %in% treatment) {
    stop("treatment must define exactly two public levels and identify treated",
         call. = FALSE)
  }
  control <- setdiff(unique(treatment), treated)
  stratum_levels <- unique(strata)
  row_index <- matrix(
    NA_integer_, nrow = length(stratum_levels), ncol = 2L,
    dimnames = list(stratum_levels, c("control", "treated")))
  for (index in seq_len(rows)) {
    column <- if (identical(treatment[[index]], treated)) {
      "treated"
    } else {
      "control"
    }
    row <- match(strata[[index]], stratum_levels)
    if (!is.na(row_index[row, column])) {
      stop("Each public stratum-treatment combination must occur once",
           call. = FALSE)
    }
    row_index[row, column] <- index
  }
  if (anyNA(row_index)) {
    stop("Every public stratum must contain treated and control rows",
         call. = FALSE)
  }
  if (!is.numeric(standard_weights) ||
      length(standard_weights) != length(stratum_levels) ||
      is.null(names(standard_weights)) || anyNA(standard_weights) ||
      any(!is.finite(standard_weights)) || any(standard_weights < 0) ||
      !any(standard_weights > 0) || anyDuplicated(names(standard_weights)) ||
      !setequal(names(standard_weights), stratum_levels)) {
    stop(paste(
      "standard_weights must be finite, non-negative, named public weights",
      "for every stratum"), call. = FALSE)
  }
  standard_weights <- standard_weights[stratum_levels]
  scale <- max(standard_weights)
  scaled <- standard_weights / scale
  total <- sum(scaled)
  if (!is.finite(total) || total <= 0) {
    stop("standard_weights cannot be normalised safely", call. = FALSE)
  }
  event <- .dsvert_dp_dimension_index(
    event, colnames(x$table), 2L, "event")
  list(
    row_index = row_index,
    strata = stratum_levels,
    treated = treated,
    control = control,
    weights = unname(scaled / total),
    event = event,
    nonevent = setdiff(1:2, event))
}

.dsvert_dp_causal_effect_regions <- function(treated, control) {
  difference <- c(
    lower = treated[["lower"]] - control[["upper"]],
    upper = treated[["upper"]] - control[["lower"]])
  list(
    risk_treated = treated,
    risk_control = control,
    risk_difference = difference,
    risk_ratio = .dsvert_dp_probability_ratio_region(treated, control),
    odds_ratio = .dsvert_dp_nonnegative_ratio_region(
      .dsvert_dp_odds_region(treated), .dsvert_dp_odds_region(control)))
}

.dsvert_dp_causal_point <- function(treated, control) {
  estimable <- length(treated) == 1L && length(control) == 1L &&
    is.finite(treated) && is.finite(control)
  if (!estimable) {
    return(list(
      values = list(
        risk_treated = NULL, risk_control = NULL,
        risk_difference = NULL, risk_ratio = NULL, odds_ratio = NULL),
      status = c(
        risk_treated = "non_estimable_zero_stratum_arm",
        risk_control = "non_estimable_zero_stratum_arm",
        risk_difference = "non_estimable_zero_stratum_arm",
        risk_ratio = "non_estimable_zero_stratum_arm",
        odds_ratio = "non_estimable_zero_stratum_arm")))
  }
  difference <- treated - control
  ratio <- if (control > 0) treated / control else if (treated > 0) Inf else NULL
  treated_odds <- if (treated < 1) treated / (1 - treated) else Inf
  control_odds <- if (control < 1) control / (1 - control) else Inf
  odds_ratio <- if (is.finite(control_odds) && control_odds > 0) {
    treated_odds / control_odds
  } else if (is.infinite(control_odds)) {
    if (is.infinite(treated_odds)) NULL else 0
  } else if (treated_odds > 0) {
    Inf
  } else {
    NULL
  }
  values <- list(
    risk_treated = treated, risk_control = control,
    risk_difference = difference, risk_ratio = ratio,
    odds_ratio = odds_ratio)
  status <- vapply(values, function(value) {
    if (is.null(value)) "non_estimable_undefined_ratio" else
      if (is.infinite(value)) "boundary_infinite" else "ok"
  }, character(1L))
  list(values = values, status = status)
}

#' DP stratified causal standardisation
#'
#' Computes a saturated, stratum-standardised g-formula from one already
#' released DP table whose rows are the public Cartesian product of strata and
#' a binary treatment and whose columns are a binary outcome. Public fixed
#' target-population weights are used; no propensity score is estimated.
#'
#' A causal interpretation additionally requires consistency, conditional
#' exchangeability within every declared stratum, positivity, no interference,
#' correct public row mapping, and scientifically valid target weights. DP
#' protects the release; it does not establish those identification assumptions.
#'
#' @param x A released `ds.vertDPContingency` with two outcome columns.
#' @param strata Public stratum label for every table row.
#' @param treatment Public binary treatment label for every table row.
#' @param treated The treatment level interpreted as treated.
#' @param standard_weights Named, non-negative public weights for all strata.
#' @param event Event column by name or index.
#' @param level Simultaneous DP-mechanism coverage.
#' @return A `ds.vertDPCausalStandardization` object. It makes no server call
#'   and consumes no additional privacy.
#' @export
ds.vertDPCausalStandardization <- function(
    x, strata, treatment, treated, standard_weights,
    event = 2L, level = 0.95) {
  x <- .dsvert_dp_table_contract(x)
  if (ncol(x$table) != 2L || nrow(x$table) < 2L ||
      !.dsvert_dp_is_number(level, 0, 1, lower_open = TRUE) || level >= 1) {
    stop("x must be a non-empty row-by-binary-outcome DP table and level in (0,1)",
         call. = FALSE)
  }
  design <- .dsvert_dp_causal_design(
    x, strata, treatment, treated, standard_weights, event)
  radius <- .dsvert_dp_table_simultaneous_radius(x, level)
  if (identical(level, 0.95) &&
      !.dsvert_dp_table_published_accuracy_matches(x, radius)) {
    stop("x does not carry a valid simultaneous DP accuracy certificate",
         call. = FALSE)
  }
  box <- .dsvert_dp_count_box(x$table, radius)
  row_regions <- t(vapply(seq_len(nrow(x$table)), function(index) {
    .dsvert_dp_risk_bounds(
      box$lower[index, design$event], box$upper[index, design$event],
      box$lower[index, design$nonevent], box$upper[index, design$nonevent])
  }, numeric(2L)))
  row_total <- rowSums(x$table)
  row_risk <- ifelse(
    row_total > 0, x$table[, design$event] / row_total, NA_real_)
  positive <- design$weights > 0
  standardized <- function(column, source) {
    rows <- design$row_index[, column]
    if (any(!is.finite(source[rows][positive]))) return(NA_real_)
    sum(design$weights[positive] * source[rows][positive])
  }
  point_treated <- standardized("treated", row_risk)
  point_control <- standardized("control", row_risk)
  point <- .dsvert_dp_causal_point(point_treated, point_control)
  region_for <- function(column) {
    rows <- design$row_index[, column]
    c(
      lower = sum(design$weights * row_regions[rows, "lower"]),
      upper = sum(design$weights * row_regions[rows, "upper"]))
  }
  regions <- .dsvert_dp_causal_effect_regions(
    region_for("treated"), region_for("control"))
  number_needed <- .dsvert_dp_number_needed_region(
    regions$risk_difference)

  result <- list(
    status = if (all(point$status == "ok")) "ok" else
      if (any(point$status == "boundary_infinite")) "boundary" else
        "dp_point_non_estimable",
    point_estimates = point$values,
    point_status = point$status,
    mechanism_regions = regions,
    mechanism_region_types = vapply(
      regions, .dsvert_dp_inference_region_type, character(1L)),
    number_needed = list(
      mechanism_regions = number_needed$regions,
      mechanism_region_includes_infinite = number_needed$includes_infinite,
      mechanism_region_possible_directions =
        number_needed$possible_directions),
    stratum_treatment_row_index = design$row_index,
    stratum_row_estimates = row_risk,
    stratum_row_regions = row_regions,
    standard_weights = stats::setNames(design$weights, design$strata),
    treated = design$treated,
    control = design$control,
    count_lower = box$lower,
    count_upper = box$upper,
    level = level,
    simultaneous_radius = radius,
    coverage_method = .dsvert_dp_table_coverage_method(x),
    uncertainty_scope = "DP mechanism noise only; sampling uncertainty excluded",
    identification_assumptions = c(
      "consistency", "conditional_exchangeability_within_public_strata",
      "positivity", "no_interference", "correct_public_row_mapping",
      "fixed_valid_target_population_weights"),
    inferential_scope = paste(
      "Saturated stratum-standardised finite-release contrast; causal",
      "interpretation is conditional on every stated identification assumption"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPCausalStandardization", "list")
  result
}

#' Joint DP and sampling inference for causal standardisation
#'
#' Unions the signed DP count-box event with exact Clopper--Pearson intervals
#' for every positive-weight stratum-treatment risk, then standardises those
#' intervals with fixed public weights. It adds sampling uncertainty to
#' `ds.vertDPCausalStandardization()` without another release.
#'
#' @inheritParams ds.vertDPCausalStandardization
#' @param mechanism_alpha_share Fraction of total non-coverage assigned to the
#'   simultaneous DP mechanism event.
#' @return A `ds.vertDPCausalStandardizationInference` object.
#' @export
ds.vertDPCausalStandardizationInference <- function(
    x, strata, treatment, treated, standard_weights,
    event = 2L, level = 0.95, mechanism_alpha_share = 0.5) {
  x <- .dsvert_dp_table_contract(x)
  values <- c(level = level, mechanism_alpha_share = mechanism_alpha_share)
  if (ncol(x$table) != 2L || nrow(x$table) < 2L ||
      any(!is.finite(values)) || level <= 0 || level >= 1 ||
      mechanism_alpha_share <= 0 || mechanism_alpha_share >= 1) {
    stop("x, level, and mechanism_alpha_share are invalid", call. = FALSE)
  }
  design <- .dsvert_dp_causal_design(
    x, strata, treatment, treated, standard_weights, event)
  alpha <- 1 - level
  mechanism_alpha <- alpha * mechanism_alpha_share
  sampling_alpha <- alpha - mechanism_alpha
  mechanism_level <- 1 - mechanism_alpha
  mechanism <- ds.vertDPCausalStandardization(
    x, strata, treatment, treated, standard_weights,
    event = event, level = mechanism_level)
  integer_box <- .dsvert_dp_integer_count_box(
    mechanism$count_lower, mechanism$count_upper)
  positive_rows <- unname(design$row_index[design$weights > 0, , drop = FALSE])
  empty <- anyNA(integer_box$lower[positive_rows, , drop = FALSE]) ||
    anyNA(integer_box$upper[positive_rows, , drop = FALSE])
  base_alpha <- sampling_alpha / length(positive_rows)
  row_regions <- matrix(
    c(0, 1), nrow = nrow(x$table), ncol = 2L, byrow = TRUE,
    dimnames = list(rownames(x$table), c("lower", "upper")))
  if (!empty) {
    for (row in positive_rows) {
      row_regions[row, ] <- .dsvert_dp_cp_union_over_box(
        integer_box$lower[row, design$event],
        integer_box$upper[row, design$event],
        integer_box$lower[row, design$nonevent],
        integer_box$upper[row, design$nonevent], base_alpha)
    }
  }
  region_for <- function(column) {
    rows <- design$row_index[, column]
    c(
      lower = sum(design$weights * row_regions[rows, "lower"]),
      upper = sum(design$weights * row_regions[rows, "upper"]))
  }
  regions <- .dsvert_dp_causal_effect_regions(
    region_for("treated"), region_for("control"))
  number_needed <- .dsvert_dp_number_needed_region(regions$risk_difference)
  result <- list(
    status = mechanism$status,
    combined_region_status = if (empty) {
      "vacuous_empty_integer_mechanism_box"
    } else {
      "ok"
    },
    point_estimates = mechanism$point_estimates,
    point_status = mechanism$point_status,
    combined_regions = regions,
    combined_region_types = vapply(
      regions, .dsvert_dp_inference_region_type, character(1L)),
    mechanism_regions = mechanism$mechanism_regions,
    confidential_count_integer_box = integer_box,
    stratum_treatment_sampling_regions = row_regions,
    standard_weights = mechanism$standard_weights,
    number_needed = list(
      combined_regions = number_needed$regions,
      combined_region_includes_infinite = number_needed$includes_infinite,
      combined_region_possible_directions =
        number_needed$possible_directions),
    level = level,
    coverage_lower_bound = level,
    mechanism_level = mechanism_level,
    sampling_familywise_level = 1 - sampling_alpha,
    base_sampling_interval_level = 1 - base_alpha,
    alpha_allocation = c(
      total = alpha, mechanism = mechanism_alpha,
      sampling_familywise = sampling_alpha,
      each_positive_stratum_treatment_interval = base_alpha),
    coverage_method = paste(
      "union bound over the signed DP count-box event and Bonferroni exact",
      "Clopper-Pearson intervals for every positive-weight public",
      "stratum-treatment risk; effects are interval post-processing"),
    sampling_model = paste(
      "conditionally binomial outcomes within each fixed public",
      "stratum-treatment arm"),
    uncertainty_scope =
      "joint DP-mechanism and superpopulation sampling uncertainty",
    identification_assumptions = mechanism$identification_assumptions,
    inferential_scope = paste(
      "Conservative standardized risk/effect regions; causal interpretation",
      "requires every stated identification assumption; no p-value"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPCausalStandardizationInference", "list")
  result
}

#' @export
print.ds.vertDPCausalStandardization <- function(x, ...) {
  cat("dsVert DP stratified causal standardisation:", x$status, "\n")
  print(unlist(x$point_estimates), ...)
  cat(x$inferential_scope, "\n")
  invisible(x)
}

#' @export
print.ds.vertDPCausalStandardizationInference <- function(x, ...) {
  cat("dsVert DP causal standardisation inference:", x$status, "\n")
  cat("Conservative joint coverage >= ", format(100 * x$level), "%\n",
      sep = "")
  print(do.call(rbind, x$combined_regions), ...)
  cat(x$inferential_scope, "\n")
  invisible(x)
}
