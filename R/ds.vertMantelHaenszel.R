.dsvert_mh_cell_roles <- function() {
  c(
    "exposed_event", "exposed_nonevent",
    "unexposed_event", "unexposed_nonevent")
}

.dsvert_mh_validate_counts <- function(x, whole_number, what) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
    stop(what, " must contain finite non-negative counts", call. = FALSE)
  }
  if (isTRUE(whole_number) &&
      any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop(what, " must contain whole-number aggregate counts", call. = FALSE)
  }
  invisible(x)
}

.dsvert_mh_dimension_index <- function(value, labels, what) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    if (is.null(labels) || !value %in% labels) {
      stop("Unknown ", what, " level: '", value, "'", call. = FALSE)
    }
    return(match(value, labels))
  }
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || !value %in% 1:2) {
    stop(what, " must identify one row/column by name or index",
         call. = FALSE)
  }
  as.integer(value)
}

.dsvert_mh_cell_map <- function(table, cell_map = NULL) {
  roles <- .dsvert_mh_cell_roles()
  if (!is.matrix(table) || ncol(table) != 4L) {
    stop("A strata-by-cells input must have exactly four columns",
         call. = FALSE)
  }
  labels <- colnames(table)
  if (is.null(cell_map)) {
    if (is.null(labels) || anyNA(labels) || any(!nzchar(labels)) ||
        anyDuplicated(labels) || !setequal(labels, roles)) {
      stop(
        "cell_map is required unless columns use the four canonical cell names",
        call. = FALSE)
    }
    indices <- match(roles, labels)
  } else {
    if (!is.atomic(cell_map) || length(cell_map) != 4L ||
        is.null(names(cell_map)) || anyNA(names(cell_map)) ||
        anyDuplicated(names(cell_map)) || !setequal(names(cell_map), roles)) {
      stop("cell_map must name each canonical 2-by-2 cell exactly once",
           call. = FALSE)
    }
    cell_map <- cell_map[roles]
    if (is.character(cell_map)) {
      if (is.null(labels) || anyNA(cell_map) || any(!nzchar(cell_map)) ||
          anyDuplicated(cell_map) || any(!cell_map %in% labels)) {
        stop("character cell_map values must be distinct table column names",
             call. = FALSE)
      }
      indices <- match(cell_map, labels)
    } else if (is.numeric(cell_map)) {
      if (anyNA(cell_map) || any(!is.finite(cell_map)) ||
          any(cell_map != floor(cell_map)) ||
          any(!cell_map %in% seq_len(ncol(table))) ||
          anyDuplicated(cell_map)) {
        stop("numeric cell_map values must be distinct column indices",
             call. = FALSE)
      }
      indices <- as.integer(cell_map)
    } else {
      stop("cell_map values must be all column names or all column indices",
           call. = FALSE)
    }
  }
  names(indices) <- roles
  indices
}

.dsvert_mh_strata_labels <- function(labels, strata) {
  if (is.null(labels)) return(paste0("stratum_", seq_len(strata)))
  if (!is.character(labels) || length(labels) != strata || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("stratum labels must be unique non-empty public names",
         call. = FALSE)
  }
  labels
}

.dsvert_mh_matrix_cells <- function(table, cell_map = NULL,
                                    whole_number = FALSE) {
  table <- as.matrix(table)
  .dsvert_mh_validate_counts(
    table, whole_number, "The strata-by-cells table")
  indices <- .dsvert_mh_cell_map(table, cell_map)
  strata <- .dsvert_mh_strata_labels(rownames(table), nrow(table))
  cells <- table[, indices, drop = FALSE]
  storage.mode(cells) <- "double"
  colnames(cells) <- .dsvert_mh_cell_roles()
  rownames(cells) <- strata
  list(
    cells = cells,
    strata = strata,
    cell_map = stats::setNames(
      if (is.null(colnames(table))) indices else colnames(table)[indices],
      .dsvert_mh_cell_roles()),
    cell_indices = indices,
    layout = "rows=strata; columns=publicly mapped 2-by-2 cells")
}

.dsvert_mh_array_cells <- function(table, exposed, event,
                                   whole_number = TRUE) {
  if (!is.array(table) || length(dim(table)) != 3L ||
      !identical(as.integer(dim(table)[1:2]), c(2L, 2L)) ||
      dim(table)[3L] < 1L) {
    stop("x must be a 2-by-2-by-K array or a K-by-4 matrix",
         call. = FALSE)
  }
  .dsvert_mh_validate_counts(table, whole_number, "x")
  labels <- dimnames(table)
  row_labels <- if (is.null(labels)) NULL else labels[[1L]]
  col_labels <- if (is.null(labels)) NULL else labels[[2L]]
  stratum_labels <- if (is.null(labels)) NULL else labels[[3L]]
  exposed_index <- .dsvert_mh_dimension_index(
    exposed, row_labels, "exposed")
  event_index <- .dsvert_mh_dimension_index(event, col_labels, "event")
  unexposed_index <- setdiff(1:2, exposed_index)
  nonevent_index <- setdiff(1:2, event_index)
  strata <- .dsvert_mh_strata_labels(stratum_labels, dim(table)[3L])
  cells <- cbind(
    exposed_event = as.numeric(table[exposed_index, event_index, ]),
    exposed_nonevent = as.numeric(
      table[exposed_index, nonevent_index, ]),
    unexposed_event = as.numeric(
      table[unexposed_index, event_index, ]),
    unexposed_nonevent = as.numeric(
      table[unexposed_index, nonevent_index, ]))
  rownames(cells) <- strata
  list(
    cells = cells,
    strata = strata,
    cell_map = stats::setNames(.dsvert_mh_cell_roles(),
                               .dsvert_mh_cell_roles()),
    cell_indices = NULL,
    layout = "dimensions=exposure,outcome,stratum",
    exposed = list(
      index = exposed_index,
      level = if (is.null(row_labels)) exposed_index else
        row_labels[[exposed_index]]),
    event = list(
      index = event_index,
      level = if (is.null(col_labels)) event_index else
        col_labels[[event_index]]))
}

.dsvert_mh_point <- function(cells) {
  roles <- .dsvert_mh_cell_roles()
  if (!is.matrix(cells) || nrow(cells) < 1L || ncol(cells) != 4L ||
      is.null(colnames(cells)) || !identical(colnames(cells), roles)) {
    stop("The oriented Mantel-Haenszel cells are invalid", call. = FALSE)
  }
  .dsvert_mh_validate_counts(cells, FALSE, "The oriented cells")
  a <- cells[, "exposed_event"]
  b <- cells[, "exposed_nonevent"]
  c0 <- cells[, "unexposed_event"]
  d <- cells[, "unexposed_nonevent"]
  totals <- a + b + c0 + d
  if (any(!is.finite(totals))) {
    stop("The stratum totals are not numerically representable",
         call. = FALSE)
  }
  diagonal_terms <- ifelse(totals > 0, (a / totals) * d, 0)
  off_diagonal_terms <- ifelse(totals > 0, (b / totals) * c0, 0)
  diagonal <- sum(diagonal_terms)
  off_diagonal <- sum(off_diagonal_terms)
  if (!is.finite(diagonal) || !is.finite(off_diagonal)) {
    stop("The Mantel-Haenszel sums are not numerically representable",
         call. = FALSE)
  }
  if (off_diagonal > 0) {
    estimate <- diagonal / off_diagonal
    estimate_type <- if (estimate == 0) "zero" else "finite"
  } else if (diagonal > 0) {
    estimate <- Inf
    estimate_type <- "infinite"
  } else {
    estimate <- NA_real_
    estimate_type <- "non_estimable"
  }
  list(
    estimate = estimate,
    estimate_type = estimate_type,
    diagonal_sum = diagonal,
    off_diagonal_sum = off_diagonal,
    diagonal_terms = diagonal_terms,
    off_diagonal_terms = off_diagonal_terms,
    stratum_totals = totals)
}

.dsvert_mh_oriented_array <- function(cells) {
  strata <- rownames(cells)
  out <- array(
    0, dim = c(2L, 2L, nrow(cells)),
    dimnames = list(
      exposure = c("exposed", "unexposed"),
      outcome = c("event", "no_event"),
      stratum = strata))
  out[1L, 1L, ] <- cells[, "exposed_event"]
  out[1L, 2L, ] <- cells[, "exposed_nonevent"]
  out[2L, 1L, ] <- cells[, "unexposed_event"]
  out[2L, 2L, ] <- cells[, "unexposed_nonevent"]
  out
}

.dsvert_mh_classical_inference <- function(cells, level, correct) {
  unavailable <- function(status, reason) {
    list(
      status = status,
      reason = reason,
      statistic = NA_real_, degrees_of_freedom = NA_real_,
      p_value = NA_real_, confidence_interval = c(
        lower = NA_real_, upper = NA_real_),
      level = level, continuity_correction = correct,
      reference = "stats::mantelhaen.test(exact = FALSE)")
  }
  if (nrow(cells) < 2L) {
    return(unavailable(
      "not_available_single_stratum",
      "Classical Mantel-Haenszel inference requires at least two strata."))
  }
  totals <- rowSums(cells)
  if (any(totals <= 1)) {
    return(unavailable(
      "not_available_small_stratum",
      "Every stratum must contain at least two observations."))
  }
  warnings <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      stats::mantelhaen.test(
        .dsvert_mh_oriented_array(cells), correct = correct,
        exact = FALSE, conf.level = level),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }),
    error = function(condition) condition)
  if (inherits(fit, "error")) {
    return(unavailable(
      "not_available_degenerate_table", conditionMessage(fit)))
  }
  interval <- unname(fit$conf.int[1:2])
  names(interval) <- c("lower", "upper")
  list(
    status = if (length(warnings)) "available_with_warning" else "available",
    reason = if (length(warnings)) paste(unique(warnings), collapse = "; ")
      else "",
    statistic = unname(fit$statistic),
    degrees_of_freedom = unname(fit$parameter),
    p_value = unname(fit$p.value),
    confidence_interval = interval,
    level = level,
    continuity_correction = correct,
    reference = "stats::mantelhaen.test(exact = FALSE)")
}

.dsvert_mh_product_over_total <- function(x, y, other_1, other_2) {
  total <- x + y + other_1 + other_2
  if (total > 0) (x / total) * y else 0
}

.dsvert_dp_mh_region <- function(lower, upper) {
  roles <- .dsvert_mh_cell_roles()
  valid <- function(x) {
    is.matrix(x) && nrow(x) >= 1L && ncol(x) == 4L &&
      identical(colnames(x), roles) && is.numeric(x) && !anyNA(x) &&
      all(is.finite(x)) && all(x >= 0)
  }
  if (!valid(lower) || !valid(upper) || !identical(dim(lower), dim(upper)) ||
      any(lower > upper)) {
    stop("The Mantel-Haenszel DP count box is invalid", call. = FALSE)
  }
  term_bounds <- function(diagonal) {
    if (diagonal) {
      x <- "exposed_event"
      y <- "unexposed_nonevent"
      other <- c("exposed_nonevent", "unexposed_event")
    } else {
      x <- "exposed_nonevent"
      y <- "unexposed_event"
      other <- c("exposed_event", "unexposed_nonevent")
    }
    term_lower <- vapply(seq_len(nrow(lower)), function(index) {
      .dsvert_mh_product_over_total(
        lower[index, x], lower[index, y],
        upper[index, other[[1L]]], upper[index, other[[2L]]])
    }, numeric(1L))
    term_upper <- vapply(seq_len(nrow(lower)), function(index) {
      .dsvert_mh_product_over_total(
        upper[index, x], upper[index, y],
        lower[index, other[[1L]]], lower[index, other[[2L]]])
    }, numeric(1L))
    c(lower = sum(term_lower), upper = sum(term_upper))
  }
  diagonal <- term_bounds(TRUE)
  off_diagonal <- term_bounds(FALSE)
  if (any(!is.finite(c(diagonal, off_diagonal)))) {
    stop("The Mantel-Haenszel DP interval is not numerically representable",
         call. = FALSE)
  }

  has_estimable_values <- off_diagonal[["upper"]] > 0 ||
    diagonal[["upper"]] > 0
  includes_zero <- diagonal[["lower"]] == 0 &&
    off_diagonal[["upper"]] > 0
  includes_infinite <- diagonal[["upper"]] > 0 &&
    off_diagonal[["lower"]] == 0
  includes_non_estimable <- diagonal[["lower"]] == 0 &&
    off_diagonal[["lower"]] == 0

  if (!has_estimable_values) {
    interval <- c(lower = NA_real_, upper = NA_real_)
    region_type <- "non_estimable"
  } else {
    lower_bound <- if (off_diagonal[["upper"]] > 0) {
      diagonal[["lower"]] / off_diagonal[["upper"]]
    } else {
      Inf
    }
    upper_bound <- if (off_diagonal[["lower"]] > 0) {
      diagonal[["upper"]] / off_diagonal[["lower"]]
    } else if (diagonal[["upper"]] > 0) {
      Inf
    } else {
      0
    }
    interval <- c(lower = lower_bound, upper = upper_bound)
    region_type <- if (is.infinite(lower_bound)) {
      "infinite_boundary"
    } else if (is.infinite(upper_bound)) {
      "unbounded_above"
    } else if (lower_bound == 0) {
      "finite_including_zero"
    } else {
      "finite_positive"
    }
  }
  list(
    interval = interval,
    region_type = region_type,
    includes_zero = includes_zero,
    includes_infinite = includes_infinite,
    includes_non_estimable = includes_non_estimable,
    has_estimable_values = has_estimable_values,
    diagonal_sum_bounds = diagonal,
    off_diagonal_sum_bounds = off_diagonal,
    construction = paste(
      "monotone per-stratum product-over-total bounds followed by",
      "non-negative interval division over the simultaneous cell box"))
}

#' Common Mantel-Haenszel odds ratio from authorised stratified tables
#'
#' Purely post-processes an already disclosure-authorised set of stratified
#' 2-by-2 aggregate tables. A three-dimensional input uses dimensions
#' exposure by outcome by stratum. A matrix uses rows as strata and four
#' columns mapped to the canonical cell roles. Bare arrays and matrices cannot
#' carry independently verifiable DataSHIELD provenance; supplying them is an
#' explicit caller attestation that the aggregates were authorised upstream.
#'
#' @param x A 2-by-2-by-K count array, a K-by-4 count matrix, or a checked
#'   object containing one of these in `observed`. A single 2-by-2 matrix is
#'   accepted as one descriptive stratum.
#' @param exposed,event For array inputs, the exposed row and event column by
#'   name or index. Ignored for a K-by-4 input.
#' @param cell_map For a K-by-4 matrix, a named character or numeric vector
#'   mapping `exposed_event`, `exposed_nonevent`, `unexposed_event`, and
#'   `unexposed_nonevent` to distinct columns. It may be omitted only when the
#'   four columns already use those exact canonical names.
#' @param level Confidence level for the asymptotic classical interval.
#' @param correct Logical; use the continuity correction in the classical
#'   Mantel-Haenszel test.
#'
#' @return A `ds.vertMantelHaenszel` object containing the descriptive common
#'   odds ratio and its explicit finite/zero/infinite/non-estimable type. When
#'   the authorised aggregate satisfies the classical test requirements, the
#'   result also contains the `stats::mantelhaen.test(exact = FALSE)` statistic,
#'   p-value and asymptotic confidence interval. No server call is made.
#'
#' @details Classical inference is conditional on the caller-provided table,
#'   independent observational units, a scientifically valid common-odds-ratio
#'   model and valid public strata. It is not a privacy guarantee and cannot
#'   repair selection, sparse-stratum, heterogeneity or provenance problems.
#' @export
ds.vertMantelHaenszel <- function(
    x, exposed = 2L, event = 2L, cell_map = NULL, level = 0.95,
    correct = TRUE) {
  carries_dp_provenance <- inherits(x, "ds.vertDPContingency") ||
    (inherits(x, "ds.vertChisq") && is.list(x) &&
       (!is.null(x$calibration) || !is.null(x$source_dp_release)))
  if (carries_dp_provenance) {
    stop("Use ds.vertDPMantelHaenszel() for a DP contingency release",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }
  if (!is.logical(correct) || length(correct) != 1L || is.na(correct)) {
    stop("correct must be TRUE or FALSE", call. = FALSE)
  }
  is_checked_object <- is.list(x) && !is.null(x$observed)
  disclosure_guard <- if (is_checked_object) x$disclosure_guard else NULL
  if (!is.null(disclosure_guard$passed) &&
      identical(disclosure_guard$passed, FALSE)) {
    stop("The input tables did not pass their disclosure guard",
         call. = FALSE)
  }
  table <- if (is_checked_object) x$observed else x
  if (is.matrix(table) && identical(dim(table), c(2L, 2L))) {
    if (!is.null(cell_map)) {
      stop("cell_map is not used for a single 2-by-2 matrix",
           call. = FALSE)
    }
    table_names <- dimnames(table)
    if (is.null(table_names)) table_names <- list(NULL, NULL)
    single <- array(
      table, dim = c(2L, 2L, 1L),
      dimnames = list(
        table_names[[1L]], table_names[[2L]], "stratum_1"))
    prepared <- .dsvert_mh_array_cells(single, exposed, event, TRUE)
  } else if (is.array(table) && length(dim(table)) == 3L) {
    if (!is.null(cell_map)) {
      stop("cell_map is used only for a K-by-4 matrix", call. = FALSE)
    }
    prepared <- .dsvert_mh_array_cells(table, exposed, event, TRUE)
  } else if (is.matrix(table)) {
    prepared <- .dsvert_mh_matrix_cells(table, cell_map, TRUE)
  } else {
    stop("x must be a 2-by-2-by-K array or a K-by-4 matrix",
         call. = FALSE)
  }
  point <- .dsvert_mh_point(prepared$cells)
  classical <- .dsvert_mh_classical_inference(
    prepared$cells, level, correct)
  estimate_ci <- classical$confidence_interval
  result <- list(
    status = if (point$estimate_type == "non_estimable") {
      "non_estimable"
    } else if (point$estimate_type %in% c("zero", "infinite")) {
      "boundary"
    } else {
      "ok"
    },
    common_odds_ratio = list(
      estimate = point$estimate,
      estimate_type = point$estimate_type,
      lower = estimate_ci[["lower"]],
      upper = estimate_ci[["upper"]]),
    mantel_haenszel_test = classical,
    oriented_tables = .dsvert_mh_oriented_array(prepared$cells),
    orientation = prepared[intersect(
      c("layout", "cell_map", "strata", "exposed", "event"),
      names(prepared))],
    stratum_components = data.frame(
      stratum = prepared$strata,
      total = point$stratum_totals,
      diagonal_over_total = point$diagonal_terms,
      off_diagonal_over_total = point$off_diagonal_terms,
      stringsAsFactors = FALSE),
    input_provenance = if (is_checked_object) {
      "disclosure_checked_authorized_aggregate"
    } else {
      "caller_attested_authorized_aggregate"
    },
    disclosure_guard = disclosure_guard,
    inferential_scope = paste(
      "Classical sampling inference is conditional on the already-authorized",
      "aggregate and the stated model assumptions; it is not DP inference"),
    additional_server_calls = 0L,
    method_status = "promoted")
  class(result) <- c("ds.vertMantelHaenszel", "list")
  result
}

#' DP common Mantel-Haenszel odds ratio from one fixed capsule table
#'
#' Purely post-processes one validated `ds.vertDPContingency` whose rows are
#' public strata and whose four columns are publicly mapped 2-by-2 cells. It
#' computes the finite-snapshot common Mantel-Haenszel odds ratio and a
#' conservative simultaneous mechanism-noise region by interval arithmetic.
#' It never computes a classical CMH p-value or sampling confidence interval.
#' The capsule must use add/remove adjacency with the sole
#' `consistent_joint_cell_else_exclude_v1` rule: each admitted unit contributes
#' to exactly one global stratum-by-cell coordinate, giving block L1
#' sensitivity one.
#'
#' @param x One released `ds.vertDPContingency` with K rows and four columns.
#' @param cell_map A public canonical cell mapping as described for
#'   `ds.vertMantelHaenszel()`; optional only for canonical column names.
#' @param level Simultaneous coverage for the DP mechanism-noise cell box.
#'
#' @return A `ds.vertDPMantelHaenszel` object. The estimate is explicitly typed
#'   as finite, zero, infinite, or non-estimable. The function performs zero
#'   server calls and consumes zero additional privacy budget.
#' @export
ds.vertDPMantelHaenszel <- function(x, cell_map = NULL, level = 0.95) {
  x <- .dsvert_dp_table_contract(x)
  if (!identical(x$adjacency, "add_remove_patient") ||
      !.dsvert_dp_num_equal(x$artifact_l1_sensitivity, 1) ||
      !identical(x$unit_aggregation_policy,
                 "consistent_joint_cell_else_exclude_v1")) {
    stop(paste(
      "x must use the one-cell add/remove contribution contract with",
      "block L1 sensitivity 1"), call. = FALSE)
  }
  if (!is.matrix(x$table) || nrow(x$table) < 1L || ncol(x$table) != 4L) {
    stop("x must contain a strata-by-four-cells DP table", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be one finite number in (0, 1)", call. = FALSE)
  }
  prepared <- .dsvert_mh_matrix_cells(x$table, cell_map, FALSE)
  radius <- .dsvert_dp_table_simultaneous_radius(x, level)
  if (identical(level, 0.95) &&
      !.dsvert_dp_table_published_accuracy_matches(x, radius)) {
    stop("x does not carry a valid simultaneous DP accuracy certificate",
         call. = FALSE)
  }
  box <- .dsvert_dp_count_box(x$table, radius)
  lower <- box$lower[, prepared$cell_indices, drop = FALSE]
  upper <- box$upper[, prepared$cell_indices, drop = FALSE]
  colnames(lower) <- colnames(upper) <- .dsvert_mh_cell_roles()
  rownames(lower) <- rownames(upper) <- prepared$strata
  point <- .dsvert_mh_point(prepared$cells)
  region <- .dsvert_dp_mh_region(lower, upper)
  result <- list(
    status = if (point$estimate_type == "non_estimable") {
      "non_estimable"
    } else if (point$estimate_type %in% c("zero", "infinite")) {
      "boundary"
    } else {
      "ok"
    },
    estimate = point$estimate,
    estimate_type = point$estimate_type,
    mechanism_region = region$interval,
    mechanism_region_type = region$region_type,
    mechanism_region_includes_zero = region$includes_zero,
    mechanism_region_includes_infinite = region$includes_infinite,
    mechanism_region_includes_non_estimable =
      region$includes_non_estimable,
    mechanism_region_has_estimable_values = region$has_estimable_values,
    interval_construction = region$construction,
    component_bounds = list(
      diagonal_sum = region$diagonal_sum_bounds,
      off_diagonal_sum = region$off_diagonal_sum_bounds),
    noisy_strata_table = prepared$cells,
    count_lower = lower,
    count_upper = upper,
    orientation = prepared[c("layout", "cell_map", "strata")],
    contribution_contract = paste(
      "each admitted privacy unit contributes to exactly one global",
      "stratum-by-cell coordinate or is excluded; add/remove block L1",
      "sensitivity = 1"),
    level = level,
    simultaneous_radius = radius,
    coverage_method = .dsvert_dp_table_coverage_method(x),
    uncertainty_scope =
      "DP mechanism noise only; sampling uncertainty excluded",
    inferential_scope = paste(
      "Finite-dataset common Mantel-Haenszel odds ratio from one DP-noised",
      "table; no classical CMH p-value or sampling confidence interval is",
      "provided"),
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    additional_server_calls = 0L,
    epsilon = x$epsilon,
    delta = x$delta,
    server = x$server)
  class(result) <- c("ds.vertDPMantelHaenszel", "list")
  result
}

#' @export
print.ds.vertMantelHaenszel <- function(x, ...) {
  cat("dsVert Mantel-Haenszel post-processing:", x$status, "\n")
  cat("common odds ratio:",
      format(x$common_odds_ratio$estimate),
      "(", x$common_odds_ratio$estimate_type, ")\n")
  if (grepl("^available", x$mantel_haenszel_test$status)) {
    cat("classical p-value:",
        format(x$mantel_haenszel_test$p_value), "\n")
  } else {
    cat("classical inference:", x$mantel_haenszel_test$status, "\n")
  }
  invisible(x)
}

#' @export
print.ds.vertDPMantelHaenszel <- function(x, ...) {
  cat("dsVert DP Mantel-Haenszel post-processing:", x$status, "\n")
  cat("common odds ratio:", format(x$estimate),
      "(", x$estimate_type, ")\n")
  cat("simultaneous DP-noise region: [",
      format(x$mechanism_region[["lower"]]), ", ",
      format(x$mechanism_region[["upper"]]), "]\n", sep = "")
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}
