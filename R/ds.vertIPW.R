#' @title Federated inverse-probability-weighted GLM (two-stage)
#' @description Convenience wrapper implementing the two-stage IPW
#'   (inverse-probability weighting) workflow on vertically partitioned
#'   DataSHIELD data.
#'
#' Stage 1: fit the propensity model via \code{ds.vertGLM} with a binomial
#'   family on the supplied \code{propensity_formula}.
#' Stage 2: fit the outcome model via \code{ds.vertGLM} with user-supplied
#'   \code{weights_column} (already holding IPW weights on the server
#'   that owns the treatment / outcome variables).
#'
#' @param outcome_formula   Formula for the outcome model.
#' @param propensity_formula Formula for the propensity (binary treatment) model.
#' @param data              Name of the aligned data frame on each server.
#' @param weights_column    Name of the column on the outcome server that
#'   holds 1/P(T=1|X) for treated units and 1/(1-P(T=1|X)) for controls,
#'   or the stabilised analogue. Pre-computed server-side (the dsVert
#'   pipeline for ON-SHARE propensity score computation is a planned
#'   follow-on; see V2_PROGRESS.md).
#' @param outcome_family    Family for the outcome model. Default "gaussian".
#' @param verbose Logical. Print stage-by-stage progress (default TRUE).
#' @param datasources DataSHIELD connections; if NULL, uses
#'   \code{DSI::datashield.connections_find()}.
#' @param ...               Passed through to both underlying
#'   \code{ds.vertGLM} calls.
#' @return list of class \code{ds.vertIPW} with \code{propensity} and
#'   \code{outcome} ds.glm objects.
#' @export
ds.vertIPW <- function(outcome_formula, propensity_formula, data = NULL,
                      weights_column = "ipw",
                      outcome_family = "gaussian",
                      verbose = TRUE, datasources = NULL, ...) {
  if (verbose) message("[ds.vertIPW] Stage 1: propensity model")
  prop_fit <- ds.vertGLM(propensity_formula, data = data,
                         family = "binomial",
                         verbose = verbose,
                         datasources = datasources, ...)

  if (verbose) message("[ds.vertIPW] Stage 2: outcome model (weighted)")
  outcome_fit <- ds.vertGLM(outcome_formula, data = data,
                            family = outcome_family,
                            weights = weights_column,
                            verbose = verbose,
                            datasources = datasources, ...)

  out <- list(
    propensity = prop_fit,
    outcome = outcome_fit,
    weights_column = weights_column,
    call = match.call())
  class(out) <- c("ds.vertIPW", "list")
  out
}

#' @export
print.ds.vertIPW <- function(x, ...) {
  cat("dsVert IPW estimator\n")
  cat("\nStage 1 -- Propensity (binomial):\n")
  print(x$propensity)
  cat("\nStage 2 -- Outcome (", x$outcome$family, ", weighted):\n", sep = "")
  print(x$outcome)
  invisible(x)
}
