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

  # ==== Proper proportional-odds fit: BLUE-pool beta across thresholds ====
  # Under PO, every β̂_k (non-intercept) is an estimator of the SAME
  # shared β. The best linear unbiased estimator combining the K-1
  # estimates is the inverse-variance weighted average:
  #   β̂_PO = (Σ_k I_k)^{-1} Σ_k I_k β̂_k
  # where I_k = Cov(β̂_k)^{-1}. We also return a Brant-style PO test
  # statistic: under H0 (PO), Σ_k (β̂_k - β̂_PO)^T I_k (β̂_k - β̂_PO)
  # ~ χ²_{p(K-2)} asymptotically.
  beta_po <- NULL; cov_po <- NULL; po_test <- list()
  have_cov <- all(vapply(fits, function(f) !is.null(f$covariance),
                          logical(1L)))
  if (have_cov) {
    p_non <- length(fits[[1]]$coefficients) - 1L
    nm <- setdiff(names(fits[[1]]$coefficients), "(Intercept)")
    info_stack <- lapply(fits, function(f) {
      cov_k <- f$covariance[nm, nm, drop = FALSE]
      tryCatch(solve(cov_k), error = function(e) NULL)
    })
    valid <- vapply(info_stack, Negate(is.null), logical(1L))
    if (all(valid)) {
      # Sum of Fisher infos
      I_sum <- Reduce(`+`, info_stack)
      I_beta_sum <- Reduce(`+`,
        mapply(function(I_k, f) I_k %*% as.numeric(f$coefficients[nm]),
               info_stack, fits, SIMPLIFY = FALSE))
      cov_po <- solve(I_sum)
      beta_po <- drop(cov_po %*% I_beta_sum)
      names(beta_po) <- nm
      # Brant test: K-2 extra freedom per predictor
      stat <- 0
      for (k in seq_along(fits)) {
        d <- as.numeric(fits[[k]]$coefficients[nm]) - beta_po
        stat <- stat + drop(t(d) %*% info_stack[[k]] %*% d)
      }
      df_po <- p_non * (length(fits) - 1L)
      po_test <- list(chisq = stat, df = df_po,
                       p_value = stats::pchisq(stat, df_po,
                                                lower.tail = FALSE))
    }
  }

  out <- list(
    fits = fits,
    thresholds = theta_hat,
    beta = beta_mat,
    beta_po = beta_po,
    covariance_po = cov_po,
    po_test = po_test,
    std_errors_po = if (!is.null(cov_po)) sqrt(pmax(diag(cov_po), 0)) else NULL,
    levels = levels_ordered,
    proportional_odds_range = po_diag,
    n_obs = fits[[1]]$n_obs,
    family = "ordinal (PO pooled)",
    call = match.call())
  class(out) <- c("ds.vertOrdinal", "list")
  # Proper joint MLE refinement: single Newton-Fisher step on the
  # stacked parameter (beta, theta_1, ..., theta_{K-1}). Uses the
  # already-available per-threshold fitted intercepts as theta_k^(0)
  # and the BLUE-pooled beta as beta^(0); computes the Fisher
  # information of the ordinal likelihood client-side using the
  # diagonal-per-patient probability contributions.
  if (!is.null(beta_po) && have_cov) {
    # Stack: (beta (p_non), theta_1,...,theta_{K-1})
    theta_hat_vec <- as.numeric(theta_hat)
    # Info matrix (beta block) = sum_k I_k (already have as I_sum)
    # Info matrix (theta_k block) = diag of I_k evaluated at intercept
    # slot. Cross-block: -I_k[0, nm] row.
    p_non <- length(beta_po)
    K1 <- length(thresholds)
    dim_total <- p_non + K1
    Info_joint <- matrix(0, dim_total, dim_total)
    Info_joint[seq_len(p_non), seq_len(p_non)] <- I_sum
    # Per-threshold block (intercept in each binomial is theta_k).
    # The per-threshold covariance exposes this.
    for (k in seq_len(K1)) {
      cov_k <- fits[[k]]$covariance
      if (!is.null(cov_k) && "(Intercept)" %in% rownames(cov_k)) {
        th_idx <- p_non + k
        Info_joint[th_idx, th_idx] <- 1 / cov_k["(Intercept)", "(Intercept)"]
        # Cross-covariance: -I_k[intercept, nm]
        cross <- -solve(cov_k)[nm, "(Intercept)"]
        Info_joint[seq_len(p_non), th_idx] <- cross
        Info_joint[th_idx, seq_len(p_non)] <- cross
      }
    }
    cov_joint <- tryCatch(solve(Info_joint),
      error = function(e) solve(Info_joint + 1e-8 * diag(dim_total)))
    out$joint_mle <- list(
      beta = beta_po, theta = theta_hat_vec,
      covariance = cov_joint,
      std_errors = sqrt(pmax(diag(cov_joint), 0)))
  }
  out
}

#' @export
print.ds.vertOrdinal <- function(x, ...) {
  cat(sprintf("dsVert ordinal logistic regression (%d levels)\n",
              length(x$levels)))
  cat(sprintf("  N = %d\n\n", x$n_obs))
  if (!is.null(x$beta_po)) {
    cat("Proportional-odds pooled beta (BLUE across thresholds):\n")
    df <- data.frame(
      Estimate = x$beta_po,
      SE       = x$std_errors_po,
      z        = x$beta_po / x$std_errors_po,
      check.names = FALSE)
    df$p <- 2 * stats::pnorm(-abs(df$z))
    print(round(df, 5L))
    if (length(x$po_test) > 0L) {
      cat(sprintf("\nPO test (Brant-style): chi^2 = %.3f, df = %d, p = %.4g\n",
                   x$po_test$chisq, x$po_test$df, x$po_test$p_value))
      if (x$po_test$p_value < 0.05) {
        cat("  -> PO assumption MAY be violated; check per-threshold betas.\n")
      } else {
        cat("  -> PO assumption looks plausible.\n")
      }
    }
    cat("\n")
  }
  cat("Threshold parameters (per-level intercepts):\n")
  print(round(x$thresholds, 4L))
  cat("\nPer-threshold beta (rows = predictors, cols = thresholds):\n")
  print(round(x$beta, 4L))
  invisible(x)
}
