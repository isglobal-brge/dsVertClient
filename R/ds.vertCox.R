#' @title Federated Cox proportional-hazards regression
#' @description User-facing Cox PH wrapper. Dispatches to
#'   \code{\link{ds.vertCoxProfileNonDisclosive}}, the non-disclosive
#'   Breslow profile route validated for K=2 and K>=3. The historical
#'   rank/permutation and person-time Poisson Cox routes were removed from the
#'   product package because they exposed event-rank metadata or required
#'   person-time-expanded inputs.
#'
#' @param formula Formula of the form \code{Surv(time, event) ~ x1 + ...}.
#' @param data Aligned data-frame name on each server.
#' @param max_iter Maximum Newton iterations for the profile route.
#' @param tol Convergence tolerance on max |delta beta|.
#' @param max_event_times Integer runtime guard passed to
#'   \code{\link{ds.vertCoxProfileNonDisclosive}}.
#' @param newton,ridge_eps,debug_trace Parameters passed through to
#'   \code{\link{ds.vertCoxProfileNonDisclosive}}.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertCox} object: \code{coefficients},
#'   \code{std_errors}, \code{covariance}, \code{loglik}
#'   (partial log-likelihood), \code{n_obs}, \code{n_events},
#'   \code{iterations}, \code{converged}.
#' @export
ds.vertCox <- function(formula, data = NULL,
                       max_iter = 30L, tol = 1e-4,
                       max_event_times = NULL,
                       newton = TRUE,
                       ridge_eps = 1e-6,
                       debug_trace = FALSE,
                       verbose = TRUE, datasources = NULL) {
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
