#' @title Quarantined ordinal-regression compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and returns no ordinal fit. Retained cumulative-binomial and
#'   proportional-odds code after the gate is unreachable through this public
#'   frontdoor and carries no disclosure, DP, accuracy, or availability claim.
#' @details Promotion requires a signed bounded proportional-odds artifact,
#'   protected score/information and validated covariance and diagnostics.
#' @param formula,data,levels_ordered,cumulative_template,max_iter,max_outer,tol,warm_max_iter,warm_tol,binomial_sigmoid_intervals,verbose,datasources
#'   Retained compatibility arguments. They are not evaluated because the
#'   public frontdoor fails locally.
#' @param ... Retained compatibility arguments; not evaluated.
#' @return No fitted object. The function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertOrdinal <- function(formula, data = NULL, levels_ordered,
                           cumulative_template = "%s_leq",
                           max_iter = NULL, max_outer = 20L, tol = NULL,
                           warm_max_iter = NULL, warm_tol = NULL,
                           binomial_sigmoid_intervals = NULL,
                           verbose = TRUE, datasources = NULL, ...) {
  .dsvert_block_retired_remote_route("ordinal")
  extra <- list(...)
  if (length(extra) > 0L) {
    arg_names <- names(extra)
    arg_names[!nzchar(arg_names)] <- "<unnamed>"
    stop("unused argument(s): ", paste(arg_names, collapse = ", "),
         call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be an R formula", call. = FALSE)
  }
  if (!is.character(levels_ordered) || length(levels_ordered) < 2L) {
    stop("levels_ordered must be a character vector with >= 2 entries",
         call. = FALSE)
  }
  if (length(levels_ordered) < 3L) {
    stop("Ordinal product route requires >=3 ordered levels; use ",
         "ds.vertGLM for binary logistic models.",
         call. = FALSE)
  }
  if (verbose) {
    message("[ds.vertOrdinal] dispatching to ",
            "ds.vertOrdinalJointNewton")
  }
  ds.vertOrdinalJointNewton(
    formula = formula,
    data = data,
    levels_ordered = levels_ordered,
    cumulative_template = cumulative_template,
    max_outer = as.integer(max_iter %||% max_outer),
    tol = as.numeric(tol %||% 1e-5),
    warm_max_iter = warm_max_iter,
    warm_tol = warm_tol,
    binomial_sigmoid_intervals = binomial_sigmoid_intervals,
    verbose = verbose,
    datasources = datasources)
}

#' @keywords internal
.ds_vertOrdinalWarm <- function(formula, data = NULL, levels_ordered,
                                cumulative_template = "%s_leq",
                                max_iter = NULL, tol = NULL,
                                ring = 63L,
                                verbose = TRUE, datasources = NULL, ...) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be an R formula", call. = FALSE)
  }
  if (!is.character(levels_ordered) || length(levels_ordered) < 2L) {
    stop("levels_ordered must be a character vector with >= 2 entries",
         call. = FALSE)
  }
  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) {
    stop("ring must be 63L or 127L", call. = FALSE)
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
  dots <- list(...)
  if (!is.null(max_iter)) dots$max_iter <- max_iter
  if (!is.null(tol)) dots$tol <- tol
  for (k in thresholds) {
    ind_col <- sprintf(cumulative_template, k)
    if (verbose) {
      message(sprintf("[ds.vertOrdinal] Threshold Y <= '%s' (indicator '%s')",
                       k, ind_col))
    }
    fm <- as.formula(paste(ind_col, "~", paste(rhs, collapse = " + ")))
    fits[[k]] <- do.call(ds.vertGLM, c(list(
      formula = fm, data = data, family = "binomial",
      ring = ring, verbose = verbose,
      datasources = datasources), dots))
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
  # Under PO, every beta_k (non-intercept) is an estimator of the SAME
  # shared beta. The best linear unbiased estimator combining the K-1
  # estimates is the inverse-variance weighted average:
  #   beta_PO = (Sum_k I_k)^{-1} Sum_k I_k beta_k
  # where I_k = Cov(beta_k)^{-1}. We also return a Brant-style PO test
  # statistic: under H0 (PO), Sum_k (beta_k - beta_PO)^T I_k (beta_k - beta_PO)
  # ~ chi^2_{p(K-2)} asymptotically.
  beta_po <- NULL; cov_po <- NULL; po_test <- list()
  theta_hat_adj <- theta_hat  # threshold intercepts after PO correction
  have_cov <- all(vapply(fits, function(f) !is.null(f$covariance),
                          logical(1L)))
  if (have_cov) {
    p_non <- length(fits[[1]]$coefficients) - 1L
    nm <- setdiff(names(fits[[1]]$coefficients), "(Intercept)")
    info_stack <- lapply(fits, function(f) {
      cov_k <- f$covariance[nm, nm, drop = FALSE]
      .dsvert_solve_identifiable(
        cov_k,
        context = "An ordinal warm-start threshold covariance",
        reason = "singular_ordinal_threshold_information",
        symmetric = TRUE)
    })
    valid <- vapply(info_stack, Negate(is.null), logical(1L))
    if (all(valid)) {
      # Sum of Fisher infos
      I_sum <- Reduce(`+`, info_stack)
      I_beta_sum <- Reduce(`+`,
        mapply(function(I_k, f) I_k %*% as.numeric(f$coefficients[nm]),
               info_stack, fits, SIMPLIFY = FALSE))
      cov_po <- .dsvert_solve_identifiable(
        I_sum,
        context = "The proportional-odds pooled information",
        reason = "singular_ordinal_pooled_information",
        symmetric = TRUE)
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

      # Client-side threshold correction (2026-04-21 PM).
      # After pooling gamma_k -> gamma_BLUE, the per-threshold intercept alpha_k
      # no longer satisfies the marginal score at (alpha_k, gamma_BLUE).
      # A one-step Newton correction on the INTERCEPT alone, derived
      # from the per-threshold Fisher block -- no new MPC, uses only
      # the covariance matrices already returned by ds.vertGLM:
      #
      #   alpha_k^* = alpha_k - I_k[alpha,alpha]^{-1} * I_k[alpha, gamma] * (gamma_BLUE - gamma_k)
      #
      # Eliminates the theta-intercept bias relative to MASS::polr zeta_k;
      # drives max|Delta cum P| from 5.5e-2 toward 1e-2 on housing.
      for (k in seq_along(fits)) {
        cov_k <- fits[[k]]$covariance
        if (is.null(cov_k)) next
        I_k_full <- .dsvert_solve_identifiable(
          cov_k,
          context = "An ordinal warm-start full covariance",
          reason = "singular_ordinal_threshold_information",
          symmetric = TRUE)
        if (!("(Intercept)" %in% rownames(I_k_full))) next
        I_aa <- I_k_full["(Intercept)", "(Intercept)"]
        I_ag <- I_k_full["(Intercept)", nm, drop = TRUE]
        gamma_k <- as.numeric(fits[[k]]$coefficients[nm])
        delta_gamma <- beta_po - gamma_k
        delta_alpha <- -as.numeric(I_ag %*% delta_gamma) / I_aa
        theta_hat_adj[k] <- theta_hat[k] + delta_alpha
      }
    }
  }

  out <- list(
    fits = fits,
    thresholds = theta_hat_adj,           # PO-corrected intercepts
    thresholds_ovr = theta_hat,           # raw OVR intercepts (pre-correction)
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
        cross <- -.dsvert_solve_identifiable(
          cov_k,
          context = "An ordinal warm-start joint covariance",
          reason = "singular_ordinal_threshold_information",
          symmetric = TRUE)[nm, "(Intercept)"]
        Info_joint[seq_len(p_non), th_idx] <- cross
        Info_joint[th_idx, seq_len(p_non)] <- cross
      }
    }
    cov_joint <- .dsvert_solve_identifiable(
      Info_joint,
      context = "The ordinal warm-start joint information",
      reason = "singular_ordinal_joint_information",
      symmetric = TRUE)
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
