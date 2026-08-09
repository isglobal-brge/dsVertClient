.DSVERT_DP_PREVALENCE_MEASURE_MAP <- c(
  exposed_prevalence = "risk_exposed",
  unexposed_prevalence = "risk_unexposed",
  population_prevalence = "population_risk",
  prevalence_difference = "risk_difference",
  prevalence_ratio = "risk_ratio",
  prevalence_odds_ratio = "odds_ratio",
  attributable_prevalence_fraction_exposed =
    "attributable_fraction_exposed",
  population_attributable_prevalence_fraction =
    "population_attributable_fraction")

.dsvert_dp_prevalence_orientation <- function(x, exposed, prevalent) {
  x <- .dsvert_dp_table_contract(x)
  if (!identical(dim(x$table), c(2L, 2L))) {
    stop("x must contain exactly a 2-by-2 DP table", call. = FALSE)
  }
  list(
    x = x,
    exposed = .dsvert_dp_dimension_index(
      exposed, rownames(x$table), 2L, "exposed"),
    prevalent = .dsvert_dp_dimension_index(
      prevalent, colnames(x$table), 2L, "prevalent"))
}

.dsvert_dp_prevalence_aliases <- function(values) {
  mapped <- values[unname(.DSVERT_DP_PREVALENCE_MEASURE_MAP)]
  names(mapped) <- names(.DSVERT_DP_PREVALENCE_MEASURE_MAP)
  mapped
}

.dsvert_dp_prevalence_view <- function(
    result, exposed, prevalent, inference = FALSE) {
  result$study_design <- "cross_sectional"
  result$study_design_source <-
    "caller_declared_by_selecting_the_cross_sectional_prevalence_view"
  result$study_design_inferred_from_table <- FALSE
  result$study_design_nonclaim <- paste(
    "The released table cannot establish study design, temporality,",
    "representativeness, exchangeability, or a causal interpretation.")
  result$orientation <- list(
    exposed = exposed,
    prevalent = prevalent,
    source = "caller_supplied_explicit_orientation")
  result$prevalence_measure_mapping <-
    .DSVERT_DP_PREVALENCE_MEASURE_MAP
  result$prevalence_point_estimates <-
    .dsvert_dp_prevalence_aliases(result$point_estimates)
  result$prevalence_point_status <-
    .dsvert_dp_prevalence_aliases(result$point_status)
  result$prevalence_mechanism_regions <-
    .dsvert_dp_prevalence_aliases(result$mechanism_regions)
  if (isTRUE(inference)) {
    result$prevalence_combined_regions <-
      .dsvert_dp_prevalence_aliases(result$combined_regions)
    result$prevalence_combined_region_types <-
      .dsvert_dp_prevalence_aliases(result$combined_region_types)
  } else {
    result$prevalence_mechanism_region_types <-
      .dsvert_dp_prevalence_aliases(result$mechanism_region_types)
  }
  result$number_needed_from_prevalence_difference <- result$number_needed
  result$number_needed_scope <- paste(
    "Descriptive reciprocal of the finite-snapshot prevalence difference;",
    "it is not a therapeutic NNT without additional causal assumptions.")
  result
}

#' Cross-sectional prevalence effects from one DP 2-by-2 release
#'
#' Provides an explicitly cross-sectional naming view of
#' [ds.vertDPEpi2x2()]. It performs no new statistical calculation: prevalence
#' difference, prevalence ratio, prevalence odds ratio, attributable
#' prevalence fractions, and the reciprocal prevalence-difference summary are
#' aliases of the corresponding risk-named values and regions in the same
#' validated result.
#'
#' @param x A released `ds.vertDPContingency` with a 2-by-2 table.
#' @param exposed Required row name or index for the exposed group. There is no
#'   default because orientation is a scientific declaration.
#' @param prevalent Required column name or index for prevalent disease or
#'   status. There is no default because orientation is a scientific
#'   declaration.
#' @param level Simultaneous coverage for DP mechanism noise, passed unchanged
#'   to [ds.vertDPEpi2x2()].
#'
#' @return A `ds.vertDPPrevalenceRatio` object that retains every field of the
#'   underlying `ds.vertDPEpi2x2` result and adds prevalence-named aliases.
#'   The aliases are numerically identical to the source fields. No server call
#'   or additional privacy release is made.
#'
#' @details Calling this function declares a cross-sectional interpretation;
#' the design is not and cannot be inferred from the released table. The
#' method adds no causal, temporal, or population-transportability claim.
#' Its mechanism regions exclude population-sampling uncertainty; use
#' [ds.vertDPPrevalenceRatioInference()] when its iid cross-sectional sampling
#' model is scientifically justified.
#' @export
ds.vertDPPrevalenceRatio <- function(x, exposed, prevalent, level = 0.95) {
  if (missing(exposed) || missing(prevalent)) {
    stop("exposed and prevalent must be supplied explicitly", call. = FALSE)
  }
  orientation <- .dsvert_dp_prevalence_orientation(
    x, exposed, prevalent)
  result <- ds.vertDPEpi2x2(
    orientation$x, exposed = orientation$exposed,
    event = orientation$prevalent, level = level)
  result <- .dsvert_dp_prevalence_view(
    result, exposed = exposed, prevalent = prevalent)
  class(result) <- c("ds.vertDPPrevalenceRatio", class(result))
  result
}

#' @rdname ds.vertDPPrevalenceRatio
#' @param mechanism_alpha_share Fraction of total non-coverage allocated to
#'   the DP mechanism event, passed unchanged to
#'   [ds.vertDPEpi2x2Inference()].
#' @return `ds.vertDPPrevalenceRatioInference()` returns a
#'   `ds.vertDPPrevalenceRatioInference` object retaining every field and
#'   coverage statement of the underlying `ds.vertDPEpi2x2Inference` result,
#'   plus prevalence-named aliases. It makes no server call and consumes no
#'   additional privacy.
#' @export
ds.vertDPPrevalenceRatioInference <- function(
    x, exposed, prevalent, level = 0.95, mechanism_alpha_share = 0.5) {
  if (missing(exposed) || missing(prevalent)) {
    stop("exposed and prevalent must be supplied explicitly", call. = FALSE)
  }
  orientation <- .dsvert_dp_prevalence_orientation(
    x, exposed, prevalent)
  result <- ds.vertDPEpi2x2Inference(
    orientation$x, exposed = orientation$exposed,
    event = orientation$prevalent, level = level,
    mechanism_alpha_share = mechanism_alpha_share)
  result <- .dsvert_dp_prevalence_view(
    result, exposed = exposed, prevalent = prevalent, inference = TRUE)
  class(result) <- c(
    "ds.vertDPPrevalenceRatioInference", class(result))
  result
}

#' @export
print.ds.vertDPPrevalenceRatio <- function(x, ...) {
  cat("dsVert DP cross-sectional prevalence ratio:", x$status, "\n")
  cat("Design: cross-sectional, caller-declared; not inferred from the table.\n")
  cat("Prevalence point estimates:\n")
  print(x$prevalence_point_estimates, ...)
  cat("Simultaneous DP-mechanism regions:\n")
  print(do.call(rbind, x$prevalence_mechanism_regions), ...)
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}

#' @export
print.ds.vertDPPrevalenceRatioInference <- function(x, ...) {
  cat("dsVert DP cross-sectional prevalence-ratio inference:",
      x$status, "\n")
  cat("Design: cross-sectional, caller-declared; not inferred from the table.\n")
  cat("Conservative joint coverage >= ", format(100 * x$level), "% (",
      x$combined_region_status, ")\n", sep = "")
  print(do.call(rbind, x$prevalence_combined_regions), ...)
  cat(x$uncertainty_scope, "\n")
  invisible(x)
}
