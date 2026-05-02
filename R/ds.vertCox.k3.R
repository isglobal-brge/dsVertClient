#' @title Federated Cox proportional-hazards via Allison-1982 Poisson trick (K=3)
#' @description Discrete-time Allison-1982 / Whitehead-1980 / Prentice-Gloeckler
#'   1978 Sec.2 equivalence: the partial-likelihood Cox estimator is the
#'   maximum-likelihood Poisson regression on person-time-expanded data
#'   with one row per (subject, at-risk-interval) pair, log-link, and a
#'   piecewise-constant baseline hazard \eqn{\alpha_t} as a factor
#'   covariate. The slopes \eqn{\beta} of the Poisson fit are the
#'   Cox PH slopes (asymptotically).
#'
#'   Implementation reduces to a single
#'   \code{\link{ds.vertGLM}} call on the K=3 secure-aggregation path
#'   with \code{family="poisson"}, \code{offset=log(at_risk * width)},
#'   and the baseline factor \code{factor(t_bin)} included as design
#'   columns. No new MPC primitives -- full inheritance of the K=3 GLM
#'   threat model (D-INV-1..5 preserved; honest-majority secure-agg
#'   envelope; OT-Beaver dishonest-majority at the DCF-pair level).
#'
#'   The expanded design is expected to live on the cluster under a
#'   user-supplied symbol (default "DA_pt"), constructed by the caller
#'   prior to invocation. The validation probe in
#'   \code{scripts/run_opal_demo_probes.R::probe_cox_k3} pre-expands
#'   the lung K=3 split client-side and uploads
#'   \code{lung_pt_s1/s2/s3} to the test cluster; production
#'   deployments would do the expansion server-side via a non-
#'   disclosive helper (see \code{coxDiscreteShareDS.R} for the K=2
#'   share-mask precedent -- the K=3 generalisation is straightforward
#'   3-party additive sharing).
#'
#' @param formula One-sided formula listing the slope covariates only --
#'   \emph{not} a Surv(...) formula. Example: \code{~ age + sex_num + bmi}.
#'   The outcome column is taken to be the user-supplied \code{event_col}.
#' @param data Aligned data-frame name on each server (e.g. "DA_pt").
#' @param event_col Outcome column on the outcome server holding the
#'   per-row 0/1 event indicator (1 = event occurred in this person-time
#'   row).
#' @param offset_col Offset column on the outcome server holding
#'   \eqn{\log(\Delta t)} where \eqn{\Delta t} is the at-risk width of the
#'   row. For not-at-risk rows the row is excluded prior to upload (or a
#'   weights column is set to 0).
#' @param baseline_col Column name on the outcome server holding the
#'   discrete-time bin index (e.g. \code{"t_bin"}). Will be included via
#'   \code{factor(t_bin)} in the design.
#' @param verbose Print progress.
#' @param datasources DataSHIELD K=3 connections.
#' @return list of class \code{ds.vertCox.k3} with \code{coefficients}
#'   (slopes only -- baseline factor levels dropped), \code{n_pp}
#'   (number of person-time rows), and the underlying
#'   \code{ds.glm} object as \code{$fit}.
#' @references
#' Allison, P. D. (1982). Discrete-time methods for the analysis of event
#'   histories. \emph{Sociological Methodology}, 13, 61-98.
#' Prentice, R. L. & Gloeckler, L. A. (1978). Regression analysis of
#'   grouped survival data with application to breast cancer data.
#'   \emph{Biometrics}, 34(1), 57-67.
#' Whitehead, J. (1980). Fitting Cox's regression model to survival
#'   data using GLIM. \emph{Applied Statistics}, 29(3), 268-275.
#' Andreux, M. \emph{et al.} (2020). Federated survival analysis with
#'   discrete-time Cox models. arXiv:2006.08997.
#' @seealso \code{\link{ds.vertGLM}}, \code{\link{ds.vertCox}} (K=2
#'   masked-Newton path).
#' @export
ds.vertCox.k3 <- function(formula, data, event_col,
                           offset_col, baseline_col,
                           verbose = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (length(datasources) < 3L)
    stop("ds.vertCox.k3 requires K=3 connections (got ",
         length(datasources), ")", call. = FALSE)
  if (!inherits(formula, "formula"))
    stop("formula must be a one-sided ~ x1 + x2 + ... formula", call. = FALSE)
  rhs <- attr(terms(formula), "term.labels")
  if (length(rhs) == 0L)
    stop("formula must list at least one slope covariate", call. = FALSE)

  ## Build the Poisson formula: event ~ x1 + x2 + ... + factor(t_bin)
  ## with offset on the outcome server. ds.vertGLM auto-detects which
  ## server holds each variable.
  poisson_formula <- as.formula(
    sprintf("%s ~ %s + factor(%s)",
            event_col,
            paste(rhs, collapse = " + "),
            baseline_col))

  if (verbose) message("[ds.vertCox.k3] Allison-1982 Poisson trick -- ",
                        "K=3 GLM family=poisson on person-time expansion.")

  fit <- ds.vertGLM(formula = poisson_formula,
                    data = data,
                    family = "poisson",
                    offset = offset_col,
                    verbose = verbose,
                    datasources = datasources)

  ## Drop the baseline factor levels from the reported slopes --
  ## they are nuisance parameters representing the piecewise-constant
  ## log-baseline-hazard and are not part of the Cox beta.
  bn <- names(fit$coefficients)
  is_baseline <- grepl(sprintf("^factor\\(%s\\)", baseline_col), bn)
  is_intercept <- bn %in% c("(Intercept)")
  slopes <- fit$coefficients[!is_baseline & !is_intercept]

  out <- list(
    coefficients = slopes,
    fit = fit,
    n_pp = fit$n_obs,
    family = "Cox PH (Allison-1982 Poisson trick, K=3)",
    call = match.call())
  class(out) <- c("ds.vertCox.k3", "list")
  out
}

#' @export
print.ds.vertCox.k3 <- function(x, ...) {
  cat("dsVert Cox PH (K=3, Allison-1982 Poisson trick)\n")
  cat(sprintf("  Person-time rows: %d\n", x$n_pp))
  cat("\nSlope coefficients (beta):\n")
  print(round(x$coefficients, 5))
  cat("\nReference: survival::coxph(Surv(time, status) ~ ...). ",
      "Baseline strata factor levels are dropped from the slopes ",
      "above; they live in the underlying $fit object.\n", sep = "")
  invisible(x)
}
