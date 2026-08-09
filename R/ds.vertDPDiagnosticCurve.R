.dsvert_dp_auc_from_ordered_bins <- function(cases, controls) {
  valid <- function(value) {
    is.numeric(value) && length(value) > 1L && !anyNA(value) &&
      all(is.finite(value)) && all(value >= 0)
  }
  if (!valid(cases) || !valid(controls) || length(cases) != length(controls)) {
    stop("Ordered diagnostic-bin counts are invalid", call. = FALSE)
  }
  case_total <- sum(cases)
  control_total <- sum(controls)
  if (case_total <= 0 || control_total <= 0) return(NA_real_)
  controls_below <- c(0, head(cumsum(controls), -1L))
  numerator <- sum(cases * (controls_below + 0.5 * controls))
  max(0, min(1, numerator / (case_total * control_total)))
}

.dsvert_dp_auc_ordered_bin_region <- function(
    case_lower, case_upper, control_lower, control_upper) {
  values <- list(case_lower, case_upper, control_lower, control_upper)
  valid <- vapply(values, function(value) {
    is.numeric(value) && length(value) > 1L && !anyNA(value) &&
      all(is.finite(value)) && all(value >= 0)
  }, logical(1L))
  lengths <- vapply(values, length, integer(1L))
  if (!all(valid) || length(unique(lengths)) != 1L ||
      any(case_lower > case_upper) ||
      any(control_lower > control_upper)) {
    stop("The ordered diagnostic DP count box is invalid", call. = FALSE)
  }
  numerator_bound <- function(cases, controls) {
    controls_below <- c(0, head(cumsum(controls), -1L))
    sum(cases * (controls_below + 0.5 * controls))
  }
  numerator_lower <- numerator_bound(case_lower, control_lower)
  numerator_upper <- numerator_bound(case_upper, control_upper)
  denominator_lower <- sum(case_lower) * sum(control_lower)
  denominator_upper <- sum(case_upper) * sum(control_upper)
  estimable_values <- denominator_upper > 0
  includes_non_estimable <- denominator_lower <= 0
  interval <- if (!estimable_values) {
    c(lower = 0, upper = 1)
  } else {
    c(
      lower = max(0, numerator_lower / denominator_upper),
      upper = if (denominator_lower > 0) {
        min(1, numerator_upper / denominator_lower)
      } else {
        1
      })
  }
  list(
    interval = interval,
    includes_non_estimable = includes_non_estimable,
    estimable_values = estimable_values,
    tightness = paste(
      "conservative rectangular count-box propagation; numerator and",
      "denominator dependence is not jointly optimized"))
}

.dsvert_dp_diagnostic_score_order <- function(score_order, labels, size) {
  if (is.null(labels) || length(labels) != size || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("The diagnostic score bins must have unique non-empty labels",
         call. = FALSE)
  }
  if (is.null(score_order)) return(seq_len(size))
  if (is.character(score_order) && length(score_order) == size &&
      !anyNA(score_order) && !anyDuplicated(score_order) &&
      setequal(score_order, labels)) {
    return(match(score_order, labels))
  }
  if (is.numeric(score_order) && length(score_order) == size &&
      !anyNA(score_order) && all(is.finite(score_order)) &&
      all(score_order == floor(score_order)) &&
      identical(sort(as.integer(score_order)), seq_len(size))) {
    return(as.integer(score_order))
  }
  stop("score_order must be a permutation of all diagnostic score bins",
       call. = FALSE)
}

.dsvert_dp_diagnostic_curve_row <- function(
    cases, controls, case_lower, case_upper,
    control_lower, control_upper, positive) {
  negative <- setdiff(seq_along(cases), positive)
  cells <- c(
    true_positive = sum(cases[positive]),
    false_negative = sum(cases[negative]),
    false_positive = sum(controls[positive]),
    true_negative = sum(controls[negative]))
  lower <- c(
    true_positive = sum(case_lower[positive]),
    false_negative = sum(case_lower[negative]),
    false_positive = sum(control_lower[positive]),
    true_negative = sum(control_lower[negative]))
  upper <- c(
    true_positive = sum(case_upper[positive]),
    false_negative = sum(case_upper[negative]),
    false_positive = sum(control_upper[positive]),
    true_negative = sum(control_upper[negative]))
  point <- .dsvert_dp_diagnostic_values(cells)
  region <- .dsvert_dp_diagnostic_regions(lower, upper)$regions
  c(
    sensitivity = point[["sensitivity"]],
    false_positive_rate = if (is.finite(point[["specificity"]])) {
      1 - point[["specificity"]]
    } else {
      NA_real_
    },
    specificity = point[["specificity"]],
    sensitivity_lower = region$sensitivity[["lower"]],
    sensitivity_upper = region$sensitivity[["upper"]],
    false_positive_rate_lower =
      1 - region$specificity[["upper"]],
    false_positive_rate_upper =
      1 - region$specificity[["lower"]])
}

#' Diagnostic ROC curve and AUC from one ordered DP table
#'
#' Purely post-processes one validated disease-status by ordered-test-bin DP
#' contingency table. It makes no server call and consumes no additional
#' privacy. The reported regions enclose DP mechanism noise under the table's
#' simultaneous certificate; they are not sampling confidence intervals.
#'
#' @param x A released `ds.vertDPContingency` whose two rows encode disease
#'   status and whose columns are ordered diagnostic-score bins.
#' @param disease_positive Positive disease-status row by name or index.
#' @param score_order Optional permutation of every column, from low to high
#'   score. By default the released column order is used.
#' @param direction Whether higher or lower scores indicate disease.
#' @param level Simultaneous coverage for DP mechanism noise.
#' @return A `ds.vertDPROC` object containing the finite-snapshot ROC curve,
#'   tie-adjusted AUC and conservative simultaneous mechanism regions.
#' @export
ds.vertDPROC <- function(x, disease_positive, score_order = NULL,
                         direction = c("higher", "lower"), level = 0.95) {
  x <- .dsvert_dp_table_contract(x)
  if (nrow(x$table) != 2L) {
    stop("x must contain exactly two disease-status rows", call. = FALSE)
  }
  if (ncol(x$table) < 2L) {
    stop("x must contain at least two ordered diagnostic score bins",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }
  direction <- match.arg(direction)
  disease_positive <- .dsvert_dp_dimension_index(
    disease_positive, rownames(x$table), 2L, "disease_positive")
  disease_negative <- setdiff(1:2, disease_positive)
  order_index <- .dsvert_dp_diagnostic_score_order(
    score_order, colnames(x$table), ncol(x$table))
  if (identical(direction, "lower")) order_index <- rev(order_index)
  labels <- colnames(x$table)[order_index]

  radius <- .dsvert_dp_table_simultaneous_radius(x, level)
  if (identical(level, 0.95) &&
      !.dsvert_dp_table_published_accuracy_matches(x, radius)) {
    stop("x does not carry a valid simultaneous DP accuracy certificate",
         call. = FALSE)
  }
  box <- .dsvert_dp_count_box(x$table, radius)
  cases <- unname(x$table[disease_positive, order_index])
  controls <- unname(x$table[disease_negative, order_index])
  case_lower <- unname(box$lower[disease_positive, order_index])
  case_upper <- unname(box$upper[disease_positive, order_index])
  control_lower <- unname(box$lower[disease_negative, order_index])
  control_upper <- unname(box$upper[disease_negative, order_index])

  bin_count <- length(cases)
  curve_values <- t(vapply(0:bin_count, function(count) {
    positive <- if (count == 0L) integer() else
      seq.int(bin_count - count + 1L, bin_count)
    .dsvert_dp_diagnostic_curve_row(
      cases, controls, case_lower, case_upper,
      control_lower, control_upper, positive)
  }, numeric(7L)))
  cutoff_index <- c(NA_integer_, rev(seq_len(bin_count)))
  cutoff <- c(NA_character_, labels[cutoff_index[-1L]])
  rule <- ifelse(
    is.na(cutoff), "no score bin classified positive",
    paste0(if (identical(direction, "higher")) ">= " else "<= ",
           cutoff))
  curve <- data.frame(
    threshold_bin = cutoff, positive_rule = rule,
    curve_values, row.names = NULL, check.names = FALSE,
    stringsAsFactors = FALSE)

  auc <- .dsvert_dp_auc_from_ordered_bins(cases, controls)
  auc_region <- .dsvert_dp_auc_ordered_bin_region(
    case_lower, case_upper, control_lower, control_upper)
  status <- if (sum(cases) <= 0) {
    "non_estimable_zero_disease_total"
  } else if (sum(controls) <= 0) {
    "non_estimable_zero_nondisease_total"
  } else {
    "ok"
  }
  row_label <- function(index) {
    labels <- rownames(x$table)
    if (is.null(labels)) NA_character_ else labels[[index]]
  }
  result <- list(
    status = status,
    orientation = list(
      row_role = "disease_status",
      disease_positive = list(
        index = disease_positive, level = row_label(disease_positive)),
      disease_negative = list(
        index = disease_negative, level = row_label(disease_negative)),
      direction = direction),
    score_order = labels,
    curve = curve,
    auc = auc,
    auc_definition = paste(
      "finite-snapshot empirical rank AUC with one-half credit for ties",
      "between public score bins"),
    auc_mechanism_region = auc_region$interval,
    auc_mechanism_region_includes_non_estimable =
      auc_region$includes_non_estimable,
    auc_mechanism_region_has_estimable_values =
      auc_region$estimable_values,
    auc_mechanism_region_tightness = auc_region$tightness,
    noisy_table = x$table,
    count_lower = box$lower,
    count_upper = box$upper,
    level = level,
    simultaneous_radius = radius,
    coverage_method = .dsvert_dp_table_coverage_method(x),
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    inferential_scope = paste(
      "Finite-snapshot discrimination from one DP-noised ordered table;",
      "no population ROC confidence band, hypothesis test, p-value or",
      "calibration claim is provided"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon, delta = x$delta, server = x$server)
  class(result) <- c("ds.vertDPROC", "list")
  result
}

#' @export
print.ds.vertDPROC <- function(x, ...) {
  cat("dsVert DP diagnostic ROC:", x$status, "\n")
  cat("AUC:", format(x$auc), " | DP-noise region: [",
      format(x$auc_mechanism_region[["lower"]]), ", ",
      format(x$auc_mechanism_region[["upper"]]), "]\n", sep = "")
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}
