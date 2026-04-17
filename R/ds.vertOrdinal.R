#' @title Federated ordinal logistic regression (naive cumulative binomials)
#' @description Fit a K-category ordinal logistic regression by training
#'   K-1 cumulative binomial logistics of the form
#'   \eqn{P(Y \leq k | X) = \sigma(\theta_k - X^\top \beta)}, one for
#'   each threshold level. This is the NAIVE approach: each cumulative
#'   binomial is fit independently, giving per-level beta estimates that
#'   may differ across thresholds. Proper proportional-odds fitting with
#'   a single shared beta and K-1 threshold parameters jointly optimised
#'   requires a dedicated joint L-BFGS objective (Month 2 follow-on).
#'
#'   The naive version is useful for exploratory analyses and for
#'   testing the proportional-odds assumption by comparing the beta
#'   estimates across threshold levels.
#'
#' @param formula R formula with the ORDERED outcome on the LHS (passed
#'   through as a factor level name in the per-threshold formulas).
#' @param data Name of the aligned data frame on each server.
#' @param levels_ordered Character vector of the ordered levels,
#'   smallest-to-largest.
#' @param cumulative_template String format to build cumulative
#'   indicator column names, e.g. \code{"%s_leq"} means the column
#'   \code{<level>_leq} is 1 when the patient's outcome is <= that
#'   level. Columns must already exist server-side. Default "\%s_leq".
#' @param ...  Passed to each underlying \code{ds.vertGLM} call.
#' @return \code{ds.vertOrdinal} object: per-threshold fits + consolidated
#'   beta matrix + threshold-parameter vector.
#' @export
ds.vertOrdinal <- function(formula, data = NULL, levels_ordered,
                           cumulative_template = "%s_leq",
                           verbose = TRUE, datasources = NULL, ...) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be an R formula", call. = FALSE)
  }
  if (!is.character(levels_ordered) || length(levels_ordered) < 2L) {
    stop("levels_ordered must be a character vector with >= 2 entries",
         call. = FALSE)
  }
  rhs <- attr(terms(formula), "term.labels")

  # Drop the topmost level (cumulative P(Y <= top) = 1 always)
  thresholds <- head(levels_ordered, -1L)
  if (length(thresholds) < 2L) {
    stop("Need at least 2 non-trivial thresholds (3+ ordered levels)",
         call. = FALSE)
  }

  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  fits <- list()
  for (k in thresholds) {
    ind_col <- sprintf(cumulative_template, k)
    if (verbose) {
      message(sprintf("[ds.vertOrdinal] Threshold Y <= '%s' (indicator '%s')",
                       k, ind_col))
    }
    fm <- as.formula(paste(ind_col, "~", paste(rhs, collapse = " + ")))
    fits[[k]] <- ds.vertGLM(fm, data = data, family = "binomial",
                             verbose = verbose,
                             datasources = datasources, ...)
  }

  # Extract intercepts as threshold parameters and non-intercept betas
  coef_names <- names(fits[[1]]$coefficients)
  int_idx <- which(coef_names == "(Intercept)")
  theta_hat <- sapply(fits, function(f) f$coefficients[int_idx])
  names(theta_hat) <- thresholds

  beta_mat <- sapply(fits, function(f) f$coefficients[-int_idx])
  if (is.null(dim(beta_mat))) beta_mat <- matrix(beta_mat, ncol = length(thresholds),
                                                  dimnames = list(coef_names[-int_idx], thresholds))

  # Proportional odds diagnostic: if the assumption holds, rows of
  # beta_mat should be roughly constant across columns.
  po_diag <- apply(beta_mat, 1, function(row) max(row) - min(row))

  out <- list(
    fits = fits,
    thresholds = theta_hat,
    beta = beta_mat,
    levels = levels_ordered,
    proportional_odds_range = po_diag,
    n_obs = fits[[1]]$n_obs,
    family = "ordinal (naive cumulative)",
    call = match.call())
  class(out) <- c("ds.vertOrdinal", "list")
  out
}

#' @export
print.ds.vertOrdinal <- function(x, ...) {
  cat(sprintf("dsVert ordinal logistic regression (naive, %d levels)\n",
              length(x$levels)))
  cat(sprintf("  N = %d\n\n", x$n_obs))
  cat("Threshold parameters (per-level intercepts):\n")
  print(round(x$thresholds, 4L))
  cat("\nBeta estimates (rows = predictors, columns = threshold levels):\n")
  print(round(x$beta, 4L))
  cat("\nProportional-odds diagnostic (per-predictor range across thresholds):\n")
  print(round(x$proportional_odds_range, 4L))
  cat("\nSmall ranges suggest the proportional-odds assumption holds;\n")
  cat("large ranges flag violations that would benefit from a joint fit.\n")
  invisible(x)
}
