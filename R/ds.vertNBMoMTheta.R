#' @title NB regression with Method-of-Moments theta-estimator (K=2-safe)
#' @description Negative-binomial GLM with Anscombe 1950 / Saha-Paul
#'   2005 Method-of-Moments dispersion estimator. Replaces the digamma-
#'   based MLE theta-Newton (\code{\link{ds.vertNB}}) with the closed-form
#'   sample-moment estimator
#'     \deqn{\hat\theta_{\mathrm{MoM}} = \bar y^2 / (s^2 - \bar y)}
#'   where \eqn{s^2 = (\sum y_i^2 - n\bar y^2)/(n-1)} is the bias-corrected
#'   sample variance. Under the iid-mu approximation (mu == ybar) this is the
#'   reduction of the regression Saha-Paul 2005 Sec.3 moment equation, and is
#'   psi-free -- so it has ZERO new MPC primitive cost beyond the four
#'   y-aggregates already in the existing iid-mu disclosure budget.
#'
#' @details Disclosure-equivalent to \code{\link{ds.vertNB}}: only y
#'   aggregates (Sumy, Sumy^2, n, ybar) are revealed at the coordinator, which
#'   is identical to the \code{dsvertNBProfileSumsDS} cost. No
#'   per-patient mu disclosure (which would require eta-reveal at OS,
#'   currently outside the K=2-safe budget). Trade-off vs MLE: ~5-10%
#'   asymptotic efficiency loss for moderate overdispersion (Lloyd-Smith
#'   2007 PLoS ONE 2(2):e180); but consistent under correct model
#'   specification (Saha-Paul 2005 Sec.3 Theorem 1) and computable in
#'   closed form (no Newton iteration on theta).
#'
#'   Honesty note: under the iid-mu approximation, this estimator
#'   propagates the same heteroscedasticity bias as iid-mu MLE-theta; the
#'   structural full-regression MoM (Saha-Paul Method 2 with per-patient
#'   mu) requires eta at OS plaintext, which is the same disclosure
#'   pattern that was disabled in ord_joint K=2-safe. Future work:
#'   share-space mu_i Beaver vecmul to recover the full-regression form
#'   without eta-reveal.
#'
#'   References:
#'   - Anscombe 1950 *Biometrika* 37(3-4):358-382 (NB moment estimators).
#'   - Saha & Paul 2005 *Biometrics* 61(1):179-185 Sec.3 (bias-corrected
#'     regression MoM).
#'   - Lloyd-Smith 2007 *PLoS ONE* 2(2):e180 (MLE-vs-MoM efficiency
#'     comparison for overdispersed counts).
#'
#' @inheritParams ds.vertNB
#' @return Object of class \code{c("ds.vertNBMoMTheta", "ds.vertNB", "ds.glm")}
#'   compatible with \code{\link{ds.vertNB}} consumers. \code{$theta}
#'   carries the MoM estimate; \code{$theta_method = "mom"};
#'   \code{$theta_mom_underdispersed} flags pathological s^2 <= ybar cases.
#' @export
#' @seealso \code{\link{ds.vertNB}}, \code{\link{ds.vertNBFullRegTheta}}
ds.vertNBMoMTheta <- function(formula, data = NULL,
                               verbose = TRUE,
                               datasources = NULL, ...) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (verbose) message("[ds.vertNBMoMTheta] Stage 1: Poisson GLM point estimator")
  fit <- ds.vertGLM(formula, data = data, family = "poisson",
                    verbose = FALSE, datasources = datasources, ...)

  y_var <- .ds_gee_extract_lhs(formula)
  server_names <- names(datasources)
  y_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                        data, y_var)
  if (is.null(y_srv)) {
    stop("could not locate outcome server for y='", y_var, "'",
         call. = FALSE)
  }
  conn_idx <- which(server_names == y_srv)

  # Aggregate y-moment sums (Sumy, Sumy^2, n, ybar, s^2) via the new
  # dsvertNBMomentSumsDS DS function. Disclosure: 4 floats per call
  # (same magnitude as dsvertNBProfileSumsDS for MLE-theta).
  sums <- tryCatch({
    r <- DSI::datashield.aggregate(
      datasources[conn_idx],
      call("dsvertNBMomentSumsDS", data_name = data, variable = y_var))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    r
  }, error = function(e) {
    stop("dsvertNBMomentSumsDS failed: ", conditionMessage(e),
         call. = FALSE)
  })
  if (!is.finite(sums$y_mean) || !is.finite(sums$y_var)) {
    stop("dsvertNBMomentSumsDS returned non-finite y_mean / y_var (privacy floor or empty cohort)",
         call. = FALSE)
  }
  ybar <- as.numeric(sums$y_mean)
  yvar <- as.numeric(sums$y_var)
  n_total <- as.integer(sums$n_total)

  # Anscombe 1950 closed-form theta_MoM = ybar^2 / (s^2 - ybar). For under-
  # dispersed samples (s^2 <= ybar) the moment equation has no positive
  # solution; flag with NA + diagnostic and fall back to large theta
  # (approxPoisson) to keep downstream printing/coercion stable.
  underdispersed <- (yvar - ybar) <= 0
  theta_mom <- if (!underdispersed) ybar * ybar / (yvar - ybar) else NA_real_
  if (verbose) {
    if (underdispersed)
      message(sprintf("[ds.vertNBMoMTheta] underdispersed sample: s^2=%.4g <= ybar=%.4g; theta_MoM undefined, using Poisson limit",
                       yvar, ybar))
    else
      message(sprintf("[ds.vertNBMoMTheta] theta_MoM = ybar^2/(s^2-ybar) = %.4g^2/(%.4g-%.4g) = %.4g (n=%d)",
                       ybar, yvar, ybar, theta_mom, n_total))
  }

  out <- fit
  out$theta <- if (is.finite(theta_mom) && theta_mom > 0) theta_mom else Inf
  out$theta_method <- "mom"
  out$theta_mom <- theta_mom
  out$theta_mom_underdispersed <- underdispersed
  out$theta_mom_n <- n_total
  out$theta_mom_y_mean <- ybar
  out$theta_mom_y_var <- yvar
  out$theta_mom_sum_y <- as.numeric(sums$sum_y)
  out$theta_mom_sum_y_sq <- as.numeric(sums$sum_y_sq)
  out$family <- "negative_binomial_mom"
  class(out) <- unique(c("ds.vertNBMoMTheta", "ds.vertNB", class(out)))
  out
}

#' @export
print.ds.vertNBMoMTheta <- function(x, ...) {
  cat("dsVert negative-binomial regression (MoM-theta; Anscombe 1950 / Saha-Paul 2005 Sec.3)\n")
  cat(sprintf("  theta_MoM = %.4g  (sample n=%d, ybar=%.4g, s^2=%.4g)\n",
              x$theta_mom %||% NA_real_, x$theta_mom_n %||% NA_integer_,
              x$theta_mom_y_mean %||% NA_real_,
              x$theta_mom_y_var %||% NA_real_))
  if (isTRUE(x$theta_mom_underdispersed))
    cat("  WARNING: sample under-dispersed (s^2 <= ybar); theta_MoM undefined.\n")
  invisible(x)
}
