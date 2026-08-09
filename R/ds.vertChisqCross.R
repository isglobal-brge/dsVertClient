#' DP-aware inference for a signed cross-owner categorical table
#'
#' This compatibility entry point obtains exactly one fixed-domain joint-DP
#' contingency release and performs only local post-processing. It never
#' discovers columns, constructs analyst-addressable one-hot objects, opens an
#' exact table, or invokes the retired cross-count endpoints. The row/column
#' ownership and domains must already be present in the custodian-signed
#' capsule `vertical_cross_specs` contract.
#'
#' @param data One of the two signed dataset names, or an existing
#'   `ds.vertDPContingency` release.
#' @param var1,var2 Row and column variables. When `data` is an existing
#'   release these are optional orientation assertions.
#' @param correct Apply the DP-aware Yates-style correction for a 2-by-2 table.
#' @param fisher Also compute the DP-aware conditional 2-by-2 calibration from
#'   the same release. No second DSI request is made.
#' @param datasources DataSHIELD connections. Omit for an existing release.
#' @param verbose Print one non-sensitive progress message.
#' @param simulations Monte Carlo replicates for each requested calibration.
#' @param mc_confidence Confidence level for its Monte Carlo interval.
#' @return A `ds.vertChisq` result. If `fisher=TRUE`, the same object also
#'   contains `fisher`, `fisher_p`, and `source_dp_release`.
#' @export
ds.vertChisqCross <- function(
    data, var1 = NULL, var2 = NULL, correct = TRUE, fisher = FALSE,
    datasources = NULL, verbose = TRUE, simulations = 9999L,
    mc_confidence = 0.95) {
  if (!is.logical(correct) || length(correct) != 1L || is.na(correct) ||
      !is.logical(fisher) || length(fisher) != 1L || is.na(fisher) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("correct, fisher and verbose must be non-missing logical scalars",
         call. = FALSE)
  }
  .dsvert_dp_chisq_validate_simulation_inputs(
    simulations, mc_confidence)

  if (inherits(data, "ds.vertDPContingency")) {
    if (!is.null(datasources)) {
      stop("datasources must be omitted for an existing DP release",
           call. = FALSE)
    }
    release <- data
    assertions <- list(row = var1, column = var2)
    expected <- list(row = release$row_var, column = release$col_var)
    for (name in names(assertions)) {
      value <- assertions[[name]]
      if (!is.null(value) &&
          (!is.character(value) || length(value) != 1L || is.na(value) ||
           !identical(value, expected[[name]]))) {
        stop(name, " variable assertion does not match the DP release",
             call. = FALSE)
      }
    }
  } else {
    for (value in list(data = data, var1 = var1, var2 = var2)) {
      if (!is.character(value) || length(value) != 1L || is.na(value) ||
          !nzchar(value)) {
        stop("data, var1 and var2 must be non-empty strings",
             call. = FALSE)
      }
    }
    release <- ds.vertDPContingency(
      data_name = data, row_var = var1, col_var = var2,
      datasources = datasources)
  }
  if (isTRUE(verbose)) {
    message("[ds.vertChisqCross] DP-aware post-processing of one capsule release")
  }
  chisq <- .dsvert_dp_chisq_from_release(
    release, correct = correct, simulations = simulations,
    mc_confidence = mc_confidence)
  chisq$source_dp_release <- release
  chisq$cross_owner <- isTRUE(release$cross_owner)
  chisq$servers <- release$servers
  if (isTRUE(fisher)) {
    fisher_result <- .dsvert_dp_fisher_from_release(
      release, simulations = simulations, mc_confidence = mc_confidence)
    chisq$fisher <- fisher_result
    chisq$fisher_p <- fisher_result$p_value
  } else {
    chisq$fisher_p <- NULL
  }
  chisq
}
