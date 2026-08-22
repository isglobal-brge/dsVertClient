#' @title Quarantined IPW compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and computes or returns no propensity model, weights, effect
#'   estimate, or diagnostic. The former two-stage implementation was removed:
#'   it could not establish a signed treatment/outcome binding, bounded weights,
#'   or a non-rerollable causal release.
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
