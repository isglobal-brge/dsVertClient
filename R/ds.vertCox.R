#' @title Quarantined Cox compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and computes or returns no model, statistic, diagnostic, or
#'   patient-derived metadata. Retained implementation code after the local
#'   migration gate is unreachable through this frontdoor and is not evidence
#'   of disclosure safety, differential privacy, numerical accuracy, or
#'   availability.
#'
#' @details Promotion requires a purpose-bound contribution-limited Cox
#'   capsule, private risk-set evaluation, certified ties/strata/delayed-entry
#'   semantics, covariance and identifiability certificates, and independent
#'   multi-host validation.
#' @param formula,data,max_iter,tol,max_event_times,newton,ridge_eps,debug_trace,verbose,datasources
#'   Retained compatibility arguments. They are not evaluated because the
#'   public frontdoor fails locally.
#' @return No fitted object. The function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertCox <- function(formula, data = NULL,
                       max_iter = 30L, tol = 1e-4,
                       max_event_times = NULL,
                       newton = TRUE,
                       ridge_eps = 1e-6,
                       debug_trace = FALSE,
                       verbose = TRUE, datasources = NULL) {
  .dsvert_block_retired_remote_route("cox")
  fit <- ds.vertCoxProfileNonDisclosive(
    formula = formula,
    data = data,
    max_event_times = max_event_times,
    max_iter = max_iter,
    tol = tol,
    newton = newton,
    ridge_eps = ridge_eps,
    debug_trace = debug_trace,
    verbose = verbose,
    datasources = datasources)
  fit$method <- "profile_nd"
  fit$iterations <- fit$n_iter %||% fit$iterations %||% NA_integer_
  fit$std_errors <- fit$std_errors %||%
    stats::setNames(rep(NA_real_, length(fit$coefficients)),
                    names(fit$coefficients))
  fit$covariance <- fit$covariance %||% NULL
  fit$loglik <- fit$loglik %||% NA_real_
  fit$call <- match.call()
  class(fit) <- unique(c("ds.vertCox", class(fit), "list"))
  fit
}

#' @export
print.ds.vertCox <- function(x, ...) {
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
