#' @title Quarantined IPW compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and computes or returns no propensity model, weights, effect
#'   estimate, or diagnostic. Retained two-stage code after the gate is
#'   unreachable through this public frontdoor and carries no disclosure, DP,
#'   causal-identification, accuracy, or availability claim.
#' @details Promotion requires signed treatment/outcome binding, a complete
#'   propensity and outcome workflow, bounded weights and contributions,
#'   explicit estimand and identification assumptions, and validated
#'   uncertainty.
#' @param outcome_formula,propensity_formula,data,weights_column,outcome_family,verbose,datasources
#'   Retained compatibility arguments. They are not evaluated because the
#'   public frontdoor fails locally.
#' @param ... Retained compatibility arguments; not evaluated.
#' @return No fitted object. The function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertIPW <- function(outcome_formula, propensity_formula, data = NULL,
                      weights_column = "ipw",
                      outcome_family = "gaussian",
                      verbose = TRUE, datasources = NULL, ...) {
  .dsvert_block_retired_remote_route("ipw")
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
