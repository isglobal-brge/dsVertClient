#' @title Signed finite-grid Cox frontdoor
#' @description With \code{analysis_id}, selects one point from a signed
#'   same-owner Breslow Cox partial-likelihood beta grid. With
#'   \code{formal_analysis_id}, reads one already completed,
#'   two-authority-signed sticky formal Cox certificate. It returns only the
#'   certified coefficient and hazard-ratio point/range values and never
#'   exposes source records, shares or intermediate state. The finite-grid
#'   route is the only fresh Cox computation admitted here;
#'   \code{fresh_formal_analysis_id} remains sealed.
#'
#' @details The finite-grid route requires one bounded row per patient, a
#'   binary categorical event, no delayed entry, a fixed public time grid and
#'   bounded numeric predictors. It selects only from its signed beta
#'   candidates; it has no arbitrary optimizer. The read-only formal route
#'   needs an explicit formula, data name and custodian-owned analysis id. Both
#'   routes have no analyst-owned privacy parameters, ring, retry, source or
#'   MPC settings. Covariance, standard errors, p-values, baseline hazards,
#'   residuals, individual predictions and new analysis runs remain unavailable
#'   until separately attested artifacts exist.
#' @param formula,data Explicit formula and aligned data name selecting the
#'   custodian-configured completed release.
#' @param analysis_id Custodian-configured same-owner Cox partial-likelihood
#'   finite-grid selector. It requires \code{Surv(time,event) ~ predictors}
#'   and rejects legacy optimizer controls.
#' @param formal_analysis_id Custodian-configured formal Cox public certificate
#'   selector. It is read-only and cannot create a new release.
#' @param fresh_formal_analysis_id Retained compatibility alias for
#'   \code{analysis_id}. It is currently unavailable before DSI.
#' @param max_iter,tol,max_event_times,newton,ridge_eps,debug_trace,verbose,datasources
#'   Legacy compatibility arguments. They are rejected when a formal id is
#'   supplied; without it the route fails before DSI.
#' @return With \code{analysis_id}, a coefficient-only
#'   \code{dsvert_dp_cox_partial_grid} object. With
#'   \code{formal_analysis_id}, a coefficient-only
#'   \code{dsvert_formal_dp_cox} object. Otherwise a typed unavailable error.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertCox <- function(formula, data = NULL,
                       max_iter = 30L, tol = 1e-4,
                       max_event_times = NULL,
                       newton = TRUE,
                       ridge_eps = 1e-6,
                       debug_trace = FALSE, verbose = TRUE,
                       datasources = NULL, analysis_id = NULL,
                       formal_analysis_id = NULL,
                       fresh_formal_analysis_id = NULL) {
  explicit <- names(match.call(expand.dots = FALSE))[-1L]
  if (sum(!vapply(list(analysis_id, formal_analysis_id,
                       fresh_formal_analysis_id), is.null, logical(1L))) > 1L) {
    stop(paste(
      "analysis_id, formal_analysis_id and fresh_formal_analysis_id are",
      "mutually exclusive"),
         call. = FALSE)
  }
  if (!is.null(analysis_id)) {
    return(.dsvert_dp_cox_partial_grid_adapter(
      explicit, formula, data, verbose, datasources, analysis_id))
  }
  if (!is.null(formal_analysis_id)) {
    return(.dsvert_formal_cox_frontdoor_adapter(
      explicit, formula, data, verbose, datasources, formal_analysis_id))
  }
  if (!is.null(fresh_formal_analysis_id)) {
    .dsvert_block_retired_remote_route("cox", .allow_test = FALSE)
  }
  .dsvert_block_retired_remote_route("cox")
}

#' @export
print.ds.vertCox <- function(x, ...) {
  if (inherits(x, "dsvert_dp_cox_partial_grid")) {
    cat("dsVert signed Breslow Cox partial-likelihood finite grid\n")
    print(round(data.frame(coef = x$coefficients,
      `exp(coef)` = x$hazard_ratio, check.names = FALSE), 5L))
    cat("No covariance, p-values, baseline hazard or source values are released.\n")
    return(invisible(x))
  }
  if (inherits(x, "dsvert_formal_dp_cox")) {
    cat("dsVert formal Cox public release\n")
    print(round(data.frame(
      coef = x$coefficients,
      `exp(coef)` = x$hazard_ratio,
      `HR lower` = x$hazard_ratio_lower,
      `HR upper` = x$hazard_ratio_upper,
      check.names = FALSE), 5L))
    cat("No covariance, p-values, baseline hazard, predictions or source values are released.\n")
    return(invisible(x))
  }
  cat("dsVert Cox proportional hazards\n")
  cat(sprintf("  N = %d, events = %d\n", x$n_obs, x$n_events))
  cat(sprintf("  converged: %s (iterations = %d)\n",
              x$converged, x$iterations))
  df <- data.frame(
    coef       = x$coefficients,
    `exp(coef)` = exp(x$coefficients),
    check.names = FALSE)
  if (!all(is.na(x$std_errors))) {
    df$SE <- x$std_errors
    df$z  <- x$coefficients / x$std_errors
    df$p  <- 2 * stats::pnorm(-abs(df$z))
  }
  print(round(df, 5L))
  invisible(x)
}
