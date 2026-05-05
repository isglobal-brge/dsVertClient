#' @title Federated multiple imputation with Rubin pooling
#' @description Fit a GLM on vertically partitioned DataSHIELD data with
#'   multiple imputation of missing values. Imputations stay on the
#'   server holding the missing variable (via
#'   \code{dsvertImputeColumnDS}); the client only ever sees the
#'   \eqn{M} pooled coefficient vectors and covariance matrices and
#'   applies Rubin's rules client-side.
#'
#'   Protocol for each of \eqn{m = 1..M}:
#'     1. On every server holding a column with missingness, call
#'        \code{dsvertImputeColumnDS} with seed \eqn{s_m}. The server
#'        draws a local imputation using a Bayesian-ridge model
#'        conditional on the other complete-case columns available on
#'        that server. The imputed column is written back into the
#'        aligned data frame under a per-round name (e.g.
#'        \code{__dsvert_imp_<var>_<m>}).
#'     2. Run \code{\link{ds.vertGLM}} on the imputed data, collect
#'        \code{beta_m} and \code{Cov(beta_m)}.
#'     3. Client accumulates \eqn{(\beta_m, \mathrm{Cov}_m)}.
#'
#'   Rubin's pooling rules are applied client-side:
#'     \deqn{\bar\beta = \frac{1}{M}\sum_m \beta_m}
#'     \deqn{W = \frac{1}{M}\sum_m \mathrm{Cov}_m}
#'     \deqn{B = \frac{1}{M-1} \sum_m (\beta_m-\bar\beta)(\beta_m-\bar\beta)^T}
#'     \deqn{T = W + (1 + 1/M) B}
#'
#' @param formula Model formula.
#' @param data Aligned data-frame name.
#' @param impute_columns Character vector of column names with
#'   missingness that should be imputed (on whichever server holds
#'   them). Per-server column presence is auto-detected.
#' @param m Number of imputations (default 20).
#' @param max_iter Inner \code{ds.vertGLM} \code{max_iter}.
#' @param tol Convergence tolerance for inner fits.
#' @param family GLM family.
#' @param lambda L2 regularisation for inner fits.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connection object.
#' @param seed RNG seed (default 1L). Per-round seed = \code{seed + m}.
#' @return A \code{ds.vertMI} object with fields \code{coefficients},
#'   \code{covariance} (Rubin total variance T), \code{std_errors},
#'   \code{within}, \code{between}, \code{fmi} (fraction of missing
#'   information), \code{m}, \code{family}, \code{fits} (list of the M
#'   inner \code{ds.glm} fits).
#' @export
ds.vertMI <- function(formula, data = NULL, impute_columns = NULL,
                      m = 20L, family = "gaussian",
                      max_iter = 50L, tol = 1e-4, lambda = 1e-4,
                      verbose = TRUE, datasources = NULL, seed = 1L) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)
  if (is.null(impute_columns) || length(impute_columns) == 0L) {
    stop("impute_columns required: name at least one column with missingness",
         call. = FALSE)
  }
  m <- as.integer(m)
  if (m < 2L) stop("m must be >= 2", call. = FALSE)

  # Auto-detect which server holds each impute_column.
  col_locs <- list()
  for (srv in server_names) {
    ci <- which(server_names == srv)
    cols <- tryCatch(
      DSI::datashield.aggregate(datasources[ci],
        call(name = "dsvertColNamesDS", data_name = data))[[1]]$columns,
      error = function(e) character(0))
    for (v in intersect(impute_columns, cols)) col_locs[[v]] <- srv
  }
  missing_vars <- setdiff(impute_columns, names(col_locs))
  if (length(missing_vars) > 0L) {
    stop("impute_columns not found on any server: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  if (verbose) {
    message(sprintf("[ds.vertMI] M=%d imputations, variables: %s",
                     m, paste(impute_columns, collapse = ", ")))
  }

  fits <- vector("list", m)
  imputation_log <- list()
  for (mi in seq_len(m)) {
    if (verbose) message(sprintf("[ds.vertMI] Imputation round %d/%d", mi, m))
    round_tag <- sprintf("__mi_%d", mi)
    for (v in names(col_locs)) {
      srv <- col_locs[[v]]
      ci <- which(server_names == srv)
      imp_res <- tryCatch(
        DSI::datashield.aggregate(datasources[ci],
          call(name = "dsvertImputeColumnDS",
               data_name = data,
               impute_column = v,
               output_column = paste0(v, round_tag),
               seed = as.integer(seed + mi))),
        error = function(e) {
          stop("dsvertImputeColumnDS failed on server '", srv, "': ",
               conditionMessage(e), "\nEnsure dsVert >= 1.1.0 is ",
               "deployed (provides the server-side imputation helper).",
               call. = FALSE)
        })
      if (is.list(imp_res) && length(imp_res) == 1L) imp_res <- imp_res[[1L]]
      imputation_log[[length(imputation_log) + 1L]] <- data.frame(
        round = mi,
        variable = v,
        server = srv,
        n_imputed = as.integer(imp_res$n_imputed %||% NA_integer_),
        n_observed = as.integer(imp_res$n_observed %||% NA_integer_),
        method = as.character(imp_res$method %||% NA_character_),
        n_predictors = as.integer(imp_res$n_predictors %||% NA_integer_),
        intercept_only = as.logical(imp_res$intercept_only %||% NA),
        stringsAsFactors = FALSE)
    }
    # Swap each impute_column for its imputed twin in the formula.
    f_txt <- deparse(formula)
    for (v in impute_columns) {
      f_txt <- gsub(paste0("\\b", v, "\\b"),
                    paste0(v, round_tag), f_txt)
    }
    f_mi <- stats::as.formula(paste(f_txt, collapse = " "))
    fits[[mi]] <- ds.vertGLM(formula = f_mi, data = data, family = family,
                              max_iter = max_iter, tol = tol,
                              lambda = lambda, verbose = FALSE,
                              datasources = datasources)
  }

  # Rubin pooling.
  betas <- do.call(cbind, lapply(fits, function(f) as.numeric(f$coefficients)))
  covs <- lapply(fits, function(f) {
    if (is.null(f$covariance)) {
      stop("Inner ds.vertGLM did not expose covariance; refit with ",
           "dsVert >= 8bb7902.", call. = FALSE)
    }
    as.matrix(f$covariance)
  })
  beta_bar <- rowMeans(betas)
  W <- Reduce(`+`, covs) / m
  dev <- sweep(betas, 1L, beta_bar)
  B <- (dev %*% t(dev)) / max(m - 1L, 1L)
  Tmat <- W + (1 + 1 / m) * B
  se <- sqrt(pmax(diag(Tmat), 0))
  nm <- names(fits[[1]]$coefficients)
  restore_impute_names <- function(x) {
    for (v in impute_columns) {
      x <- sub(paste0("^", v, "__mi_[0-9]+$"), v, x)
    }
    x
  }
  nm <- restore_impute_names(nm)
  names(beta_bar) <- names(se) <- nm
  dimnames(W) <- dimnames(B) <- dimnames(Tmat) <- list(nm, nm)
  target_order <- c("(Intercept)", attr(stats::terms(formula), "term.labels"))
  target_order <- target_order[target_order %in% nm]
  ord <- match(c(target_order, setdiff(nm, target_order)), nm)
  if (length(ord) == length(nm) && all(!is.na(ord))) {
    beta_bar <- beta_bar[ord]
    se <- se[ord]
    W <- W[ord, ord, drop = FALSE]
    B <- B[ord, ord, drop = FALSE]
    Tmat <- Tmat[ord, ord, drop = FALSE]
    nm <- nm[ord]
  }
  # Fraction of missing information (per-coefficient).
  lambda_hat <- (1 + 1 / m) * diag(B) / pmax(diag(Tmat), 1e-30)
  r <- (1 + 1 / m) * diag(B) / pmax(diag(W), 1e-30)
  df <- (m - 1) * (1 + 1 / r)^2
  fmi <- (r + 2 / (df + 3)) / (r + 1)
  names(lambda_hat) <- names(fmi) <- nm
  imputation_log <- if (length(imputation_log)) {
    do.call(rbind, imputation_log)
  } else {
    data.frame()
  }
  quality_warnings <- character(0)
  quality_status <- "ok"
  if (nrow(imputation_log) > 0L) {
    if (all((imputation_log$n_imputed %||% 0L) == 0L, na.rm = TRUE)) {
      quality_warnings <- c(quality_warnings,
        "No missing cells were imputed. If missingness existed before alignment, run ds.psiAlign(..., na.action = 'none') before ds.vertMI().")
      quality_status <- "degraded"
    }
    aggregate_fallback <- imputation_log$intercept_only %in% TRUE &
      imputation_log$method %in% c("mean_intercept", "mode_intercept")
    unstable_intercept <- imputation_log$intercept_only %in% TRUE &
      !aggregate_fallback
    if (any(unstable_intercept)) {
      quality_warnings <- c(quality_warnings,
        "At least one imputation used an intercept-only local model because no complete numeric predictor was available on the imputation server.")
      quality_status <- "degraded"
    }
    if (any(aggregate_fallback)) {
      quality_warnings <- c(quality_warnings,
        "At least one imputation used an aggregate-only mean/mode fallback because no complete numeric predictor was available on the imputation server.")
      if (identical(quality_status, "ok")) quality_status <- "approximate"
    }
  }

  out <- list(
    coefficients = beta_bar,
    covariance   = Tmat,
    within       = W,
    between      = B,
    std_errors   = se,
    lambda_hat   = lambda_hat,
    fmi          = fmi,
    m            = m,
    family       = family,
    n_obs        = fits[[1]]$n_obs,
    imputation_log = imputation_log,
    quality      = list(status = quality_status,
                        warnings = quality_warnings),
    fits         = fits,
    call         = match.call())
  class(out) <- c("ds.vertMI", "list")
  out
}

#' @export
print.ds.vertMI <- function(x, ...) {
  cat("dsVert multiple-imputation GLM (Rubin-pooled)\n")
  cat(sprintf("  M = %d imputations   family = %s   N = %d\n",
              x$m, x$family, x$n_obs))
  if (!is.null(x$quality$status)) {
    cat(sprintf("  Quality: %s\n", x$quality$status))
    if (length(x$quality$warnings)) {
      for (w in x$quality$warnings) cat("  - ", w, "\n", sep = "")
    }
  }
  z <- x$coefficients / x$std_errors
  p <- 2 * stats::pnorm(-abs(z))
  df <- data.frame(
    Estimate = x$coefficients,
    SE       = x$std_errors,
    z        = z,
    p        = p,
    FMI      = x$fmi,
    check.names = FALSE)
  print(round(df, 5L))
  invisible(x)
}
