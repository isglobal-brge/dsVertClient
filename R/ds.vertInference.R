#' @title Likelihood-ratio test on two nested ds.vertGLM fits
#' @description Compute the standard LR chi-square statistic between a reduced
#'   and a full ds.vertGLM fit. Both fits must be on the same cohort and
#'   the reduced model must be nested within the full model.
#'
#' @param reduced ds.glm object with fewer coefficients.
#' @param full ds.glm object with more coefficients.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{statistic}: 2 * (reduced$deviance - full$deviance)
#'     \item \code{df}: full$n_vars - reduced$n_vars
#'     \item \code{p_value}: upper-tail chi-square p-value
#'     \item \code{deviance_reduced}, \code{deviance_full}
#'   }
#'
#' @details This is a pure client-side computation over the scalar
#'   deviances already returned by ds.vertGLM; no additional MPC round
#'   is performed and no observation-level information is exposed.
#' @export
ds.vertLR <- function(reduced, full) {
  if (!inherits(reduced, "ds.glm")) {
    stop("`reduced` must be a ds.glm object (from ds.vertGLM)", call. = FALSE)
  }
  if (!inherits(full, "ds.glm")) {
    stop("`full` must be a ds.glm object (from ds.vertGLM)", call. = FALSE)
  }
  if (!identical(reduced$family, full$family)) {
    stop("LR test requires both fits to share the same family", call. = FALSE)
  }
  if (!isTRUE(reduced$n_obs == full$n_obs)) {
    stop("LR test requires both fits on the same cohort (n_obs differs)",
         call. = FALSE)
  }
  if (!all(names(reduced$coefficients) %in% names(full$coefficients))) {
    stop("Coefficients of `reduced` must be a subset of `full` (nested)",
         call. = FALSE)
  }

  df <- full$n_vars - reduced$n_vars
  if (df <= 0L) {
    stop("Full model must have strictly more coefficients than reduced",
         call. = FALSE)
  }

  stat <- as.numeric(reduced$deviance - full$deviance)
  if (!is.finite(stat) || stat < 0) {
    # Numerical noise can push a well-nested fit slightly negative; floor at 0
    stat <- max(stat, 0)
  }
  p <- pchisq(stat, df = df, lower.tail = FALSE)

  out <- list(
    statistic = stat,
    df = as.integer(df),
    p_value = p,
    deviance_reduced = as.numeric(reduced$deviance),
    deviance_full = as.numeric(full$deviance),
    family = full$family,
    n_obs = full$n_obs)
  class(out) <- c("ds.vertLR", "list")
  out
}

#' @export
print.ds.vertLR <- function(x, ...) {
  cat("dsVert likelihood-ratio test\n")
  cat(sprintf("  Family : %s\n", x$family))
  cat(sprintf("  N      : %d\n", x$n_obs))
  cat(sprintf("  Deviance: reduced=%.4f  full=%.4f\n",
              x$deviance_reduced, x$deviance_full))
  cat(sprintf("  LR chi-sq = %.4f on %d df,  p-value = %s\n",
              x$statistic, x$df, format.pval(x$p_value, digits = 4L)))
  invisible(x)
}


#' @title Wald confidence intervals for ds.vertGLM coefficients
#' @description Return Wald-type confidence intervals using the
#'   finite-difference standard errors already stored in a ds.glm object.
#'   Observation-level quantities are never touched; this is a scalar
#'   client-side transformation.
#'
#' @param fit A ds.glm object.
#' @param parm Optional character vector of coefficient names to report;
#'   default all.
#' @param level Confidence level, default 0.95.
#'
#' @return A data frame with columns \code{estimate}, \code{std_error},
#'   \code{lower}, \code{upper}; row names are coefficient names.
#' @export
ds.vertConfint <- function(fit, parm = NULL, level = 0.95) {
  if (!inherits(fit, "ds.glm")) {
    stop("`fit` must be a ds.glm object", call. = FALSE)
  }
  if (!is.numeric(level) || level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1)", call. = FALSE)
  }
  coefs <- fit$coefficients
  ses <- fit$std_errors
  if (is.null(names(coefs)) || is.null(names(ses))) {
    stop("fit must expose named coefficients and std_errors", call. = FALSE)
  }
  if (!is.null(parm)) {
    missing_parm <- setdiff(parm, names(coefs))
    if (length(missing_parm) > 0L) {
      stop("Unknown coefficient(s): ", paste(missing_parm, collapse = ", "),
           call. = FALSE)
    }
    idx <- match(parm, names(coefs))
  } else {
    idx <- seq_along(coefs)
  }
  z <- qnorm((1 + level) / 2)
  est <- as.numeric(coefs[idx])
  se <- as.numeric(ses[idx])
  out <- data.frame(
    estimate = est,
    std_error = se,
    lower = est - z * se,
    upper = est + z * se,
    stringsAsFactors = FALSE)
  rownames(out) <- names(coefs)[idx]
  attr(out, "level") <- level
  out
}


#' @title Univariate Wald test for a single ds.vertGLM coefficient
#' @description Test H0: beta_j = null against a two-sided alternative using
#'   the diagonal Wald statistic (estimate - null) / SE. A convenience
#'   wrapper over the z_values / p_values already stored in a ds.glm
#'   object, extended to non-zero nulls.
#'
#' @param fit A ds.glm object.
#' @param parm Single coefficient name.
#' @param null Null value (default 0).
#'
#' @return List with estimate, SE, z, p_value, null.
#' @export
ds.vertWald <- function(fit, parm, null = 0) {
  if (!inherits(fit, "ds.glm")) {
    stop("`fit` must be a ds.glm object", call. = FALSE)
  }
  if (!is.character(parm) || length(parm) != 1L) {
    stop("`parm` must be a single coefficient name", call. = FALSE)
  }
  if (!parm %in% names(fit$coefficients)) {
    stop("Coefficient '", parm, "' not in fit", call. = FALSE)
  }
  est <- as.numeric(fit$coefficients[parm])
  se <- as.numeric(fit$std_errors[parm])
  z <- (est - null) / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  out <- list(
    parm = parm,
    estimate = est,
    std_error = se,
    null = null,
    z = z,
    p_value = p)
  class(out) <- c("ds.vertWald", "list")
  out
}

#' @export
print.ds.vertWald <- function(x, ...) {
  cat(sprintf("dsVert Wald test: H0: %s = %g  vs  two-sided\n",
              x$parm, x$null))
  cat(sprintf("  estimate = %.6f  SE = %.6f  z = %.4f\n",
              x$estimate, x$std_error, x$z))
  cat(sprintf("  p-value  = %s\n", format.pval(x$p_value, digits = 4L)))
  invisible(x)
}
