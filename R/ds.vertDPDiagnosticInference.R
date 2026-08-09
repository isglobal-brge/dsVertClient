.dsvert_dp_cp_union_named <- function(
    lower, upper, events, nonevents, alpha) {
  expected <- c(
    "true_positive", "false_negative",
    "false_positive", "true_negative")
  if (!is.numeric(lower) || !is.numeric(upper) ||
      !identical(names(lower), expected) ||
      !identical(names(upper), expected) ||
      !length(events) || !length(nonevents) ||
      anyDuplicated(c(events, nonevents)) ||
      any(!c(events, nonevents) %in% expected)) {
    stop("Invalid diagnostic integer count box", call. = FALSE)
  }
  .dsvert_dp_cp_union_over_box(
    sum(lower[events]), sum(upper[events]),
    sum(lower[nonevents]), sum(upper[nonevents]), alpha)
}

.dsvert_dp_probability_complement_region <- function(probability) {
  if (!is.numeric(probability) || length(probability) != 2L ||
      !identical(names(probability), c("lower", "upper")) ||
      anyNA(probability) || any(!is.finite(probability)) ||
      probability[["lower"]] < 0 || probability[["upper"]] > 1 ||
      probability[["lower"]] > probability[["upper"]]) {
    stop("Invalid probability region", call. = FALSE)
  }
  c(
    lower = 1 - probability[["upper"]],
    upper = 1 - probability[["lower"]])
}

.dsvert_dp_probability_mean_region <- function(left, right) {
  valid <- function(value) {
    is.numeric(value) && length(value) == 2L &&
      identical(names(value), c("lower", "upper")) &&
      !anyNA(value) && all(is.finite(value)) &&
      value[["lower"]] >= 0 && value[["upper"]] <= 1 &&
      value[["lower"]] <= value[["upper"]]
  }
  if (!valid(left) || !valid(right)) {
    stop("Invalid probability regions", call. = FALSE)
  }
  c(
    lower = (left[["lower"]] + right[["lower"]]) / 2,
    upper = (left[["upper"]] + right[["upper"]]) / 2)
}

.dsvert_dp_probability_harmonic_region <- function(left, right) {
  valid <- function(value) {
    is.numeric(value) && length(value) == 2L &&
      identical(names(value), c("lower", "upper")) &&
      !anyNA(value) && all(is.finite(value)) &&
      value[["lower"]] >= 0 && value[["upper"]] <= 1 &&
      value[["lower"]] <= value[["upper"]]
  }
  if (!valid(left) || !valid(right)) {
    stop("Invalid probability regions", call. = FALSE)
  }
  harmonic <- function(x, y) {
    if (x + y > 0) 2 * x * y / (x + y) else 0
  }
  c(
    lower = harmonic(left[["lower"]], right[["lower"]]),
    upper = harmonic(left[["upper"]], right[["upper"]]))
}

.dsvert_dp_diagnostic_integer_vectors <- function(integer_box) {
  if (!is.list(integer_box) || !identical(integer_box$empty, FALSE) ||
      !is.matrix(integer_box$lower) || !is.matrix(integer_box$upper) ||
      !identical(dim(integer_box$lower), c(2L, 2L)) ||
      !identical(dim(integer_box$upper), c(2L, 2L))) {
    stop("Invalid diagnostic integer count box", call. = FALSE)
  }
  extract <- function(value) {
    c(
      true_positive = value[1L, 1L],
      false_negative = value[1L, 2L],
      false_positive = value[2L, 1L],
      true_negative = value[2L, 2L])
  }
  list(lower = extract(integer_box$lower),
       upper = extract(integer_box$upper))
}

.dsvert_dp_diagnostic_nonestimable_flags <- function(lower) {
  if (!is.numeric(lower) || anyNA(lower) || any(lower < 0) ||
      !setequal(names(lower), c(
        "true_positive", "false_negative",
        "false_positive", "true_negative"))) {
    stop("Invalid diagnostic integer lower bounds", call. = FALSE)
  }
  tp <- lower[["true_positive"]]
  fn <- lower[["false_negative"]]
  fp <- lower[["false_positive"]]
  tn <- lower[["true_negative"]]
  disease_zero <- tp + fn == 0
  nondisease_zero <- fp + tn == 0
  test_positive_zero <- tp + fp == 0
  test_negative_zero <- tn + fn == 0
  total_zero <- tp + fn + fp + tn == 0
  c(
    sensitivity = disease_zero,
    specificity = nondisease_zero,
    ppv = test_positive_zero,
    npv = test_negative_zero,
    prevalence = total_zero,
    accuracy = total_zero,
    balanced_accuracy = disease_zero || nondisease_zero,
    f1_score = 2 * tp + fp + fn == 0,
    lr_positive = disease_zero || nondisease_zero || test_positive_zero,
    lr_negative = disease_zero || nondisease_zero || test_negative_zero,
    diagnostic_odds_ratio = disease_zero || nondisease_zero ||
      test_positive_zero || test_negative_zero)
}

#' Diagnostic-accuracy inference with DP and sampling uncertainty
#'
#' Builds conservative simultaneous confidence regions from one validated
#' disease-by-test DP 2-by-2 release. It combines the signed simultaneous DP
#' count-box event with six Bonferroni exact Clopper--Pearson intervals for
#' sensitivity, specificity, PPV, NPV, prevalence, and accuracy. Balanced
#' accuracy, F1, likelihood ratios, and the diagnostic odds ratio are
#' deterministic interval post-processing of those base regions.
#'
#' The sampling statement assumes independent, identically sampled privacy
#' units from one joint disease/test population. Conditional binomial models
#' are used for sensitivity, specificity, PPV and NPV; prevalence and accuracy
#' use marginal binomial models. No independence among the six intervals is
#' assumed. Possible zero-denominator states are reported explicitly.
#'
#' @param x A released `ds.vertDPContingency` with a 2-by-2 table.
#' @param disease_positive Positive disease-status row by name or index.
#' @param test_positive Positive diagnostic-test column by name or index.
#' @param level Desired lower bound for joint sampling-plus-mechanism coverage.
#' @param mechanism_alpha_share Fraction of total non-coverage allocated to
#'   the DP mechanism event. The remainder is divided equally among the six
#'   exact-binomial base intervals.
#' @return A `ds.vertDPDiagnostic2x2Inference` object. It makes no server call
#'   and consumes no additional privacy.
#' @references Clopper, C. J. and Pearson, E. S. (1934). The use of
#'   confidence or fiducial limits illustrated in the case of the binomial.
#'   Biometrika 26(4), 404--413. \doi{10.1093/biomet/26.4.404}.
#' @export
ds.vertDPDiagnostic2x2Inference <- function(
    x, disease_positive, test_positive, level = 0.95,
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
  alpha <- 1 - level
  mechanism_alpha <- alpha * mechanism_alpha_share
  sampling_alpha <- alpha - mechanism_alpha
  mechanism_level <- 1 - mechanism_alpha
  base_sampling_alpha <- sampling_alpha / 6
  mechanism <- ds.vertDPDiagnostic2x2(
    x, disease_positive = disease_positive,
    test_positive = test_positive, level = mechanism_level)
  integer_box <- .dsvert_dp_integer_count_box(
    mechanism$count_lower, mechanism$count_upper)

  if (integer_box$empty) {
    base <- rep(list(c(lower = 0, upper = 1)), 6L)
    names(base) <- c(
      "sensitivity", "specificity", "ppv", "npv",
      "prevalence", "accuracy")
    includes_non_estimable <- setNames(
      rep(TRUE, length(mechanism$estimates)), names(mechanism$estimates))
  } else {
    cells <- .dsvert_dp_diagnostic_integer_vectors(integer_box)
    cp <- function(events, nonevents) {
      .dsvert_dp_cp_union_named(
        cells$lower, cells$upper, events, nonevents,
        base_sampling_alpha)
    }
    base <- list(
      sensitivity = cp("true_positive", "false_negative"),
      specificity = cp("true_negative", "false_positive"),
      ppv = cp("true_positive", "false_positive"),
      npv = cp("true_negative", "false_negative"),
      prevalence = cp(
        c("true_positive", "false_negative"),
        c("false_positive", "true_negative")),
      accuracy = cp(
        c("true_positive", "true_negative"),
        c("false_positive", "false_negative")))
    includes_non_estimable <-
      .dsvert_dp_diagnostic_nonestimable_flags(cells$lower)
  }

  false_positive_rate <-
    .dsvert_dp_probability_complement_region(base$specificity)
  false_negative_rate <-
    .dsvert_dp_probability_complement_region(base$sensitivity)
  regions <- c(base, list(
    balanced_accuracy = .dsvert_dp_probability_mean_region(
      base$sensitivity, base$specificity),
    f1_score = .dsvert_dp_probability_harmonic_region(
      base$sensitivity, base$ppv),
    lr_positive = .dsvert_dp_probability_ratio_region(
      base$sensitivity, false_positive_rate),
    lr_negative = .dsvert_dp_probability_ratio_region(
      false_negative_rate, base$specificity),
    diagnostic_odds_ratio = .dsvert_dp_nonnegative_ratio_region(
      .dsvert_dp_odds_region(base$sensitivity),
      .dsvert_dp_odds_region(false_positive_rate))))
  regions <- regions[names(mechanism$estimates)]
  region_types <- vapply(
    regions, .dsvert_dp_inference_region_type, character(1L))

  result <- list(
    status = mechanism$status,
    combined_region_status = if (integer_box$empty) {
      "vacuous_empty_integer_mechanism_box"
    } else if (any(includes_non_estimable)) {
      "ok_includes_non_estimable_states"
    } else {
      "ok"
    },
    orientation = mechanism$orientation,
    point_estimates = mechanism$estimates,
    point_status = mechanism$point_status,
    base_sampling_regions = base,
    combined_regions = regions,
    combined_region_types = region_types,
    combined_region_includes_non_estimable = includes_non_estimable,
    mechanism_regions = mechanism$mechanism_regions,
    confidential_count_integer_box = integer_box,
    level = level,
    coverage_lower_bound = level,
    mechanism_level = mechanism_level,
    sampling_familywise_level = 1 - sampling_alpha,
    base_sampling_interval_level = 1 - base_sampling_alpha,
    alpha_allocation = c(
      total = alpha, mechanism = mechanism_alpha,
      sampling_familywise = sampling_alpha,
      each_of_six_sampling_intervals = base_sampling_alpha),
    coverage_method = paste(
      "union bound over the signed DP simultaneous count-box event and",
      "six Bonferroni exact Clopper-Pearson intervals; derived diagnostic",
      "regions are deterministic interval post-processing; an integer-empty",
      "mechanism box returns the full parameter space"),
    sampling_model = paste(
      "iid privacy units from one joint disease/test population;",
      "conditional binomial sensitivity, specificity, PPV and NPV, with",
      "marginal binomial prevalence and accuracy"),
    uncertainty_scope =
      "joint DP-mechanism and superpopulation sampling uncertainty",
    inferential_scope = paste(
      "Conservative confidence regions for diagnostic accuracy;",
      "no hypothesis test or p-value"),
    continuity_correction = "none",
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPDiagnostic2x2Inference", "list")
  result
}

#' @export
print.ds.vertDPDiagnostic2x2Inference <- function(x, ...) {
  cat("dsVert DP diagnostic sampling inference:", x$status, "\n")
  cat("Conservative joint coverage >= ", format(100 * x$level), "%",
      " (", x$combined_region_status, ")\n", sep = "")
  print(do.call(rbind, x$combined_regions), ...)
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}
