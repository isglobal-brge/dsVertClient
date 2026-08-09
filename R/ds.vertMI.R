#' @title Quarantined multiple-imputation compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and mutates no server data, draws no imputation, and returns no
#'   coefficients, counts, covariance, or diagnostic. Retained implementation
#'   code after the gate is unreachable through this public frontdoor and
#'   carries no disclosure, DP, accuracy, or availability claim.
#' @details Promotion requires a signed bounded imputation contract,
#'   non-rerollable cryptographic randomness, no exact per-round count release,
#'   and validated Rubin-rule inference.
#' @param formula,data,impute_columns,m,family,max_iter,tol,lambda,intercept_only,verbose,datasources,seed
#'   Retained compatibility arguments. They are not evaluated because the
#'   public frontdoor fails locally.
#' @return No fitted object. The function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertMI <- function(formula, data = NULL, impute_columns = NULL,
                      m = 20L, family = "gaussian",
                      max_iter = 50L, tol = 1e-4, lambda = 0,
                      intercept_only = c("error", "aggregate"),
                      verbose = TRUE, datasources = NULL, seed = 1L) {
  .dsvert_block_retired_remote_route("mi")
  intercept_only <- match.arg(intercept_only)
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
  column_results <- .dsvert_aggregate_strict(
    datasources,
    call(name = "dsvertColNamesDS", data_name = data),
    operation = "multiple-imputation column discovery")
  for (srv in server_names) {
    cols <- column_results[[srv]]$columns
    for (v in intersect(impute_columns, cols)) {
      if (!is.null(col_locs[[v]])) {
        stop("impute column '", v,
             "' is present on more than one server; choose an unambiguous vertical partition",
             call. = FALSE)
      }
      col_locs[[v]] <- srv
    }
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
        .dsvert_aggregate_strict(datasources[ci],
          call(name = "dsvertImputeColumnDS",
               data_name = data,
               impute_column = v,
               output_column = paste0(v, round_tag),
               seed = as.integer(seed + mi),
               allow_intercept_only = identical(intercept_only,
                                                  "aggregate")),
          operation = "protected multiple-imputation update")[[1L]],
        error = function(e) {
          stop("The protected multiple-imputation update failed; ",
               "no partial round was accepted.",
               call. = FALSE)
        })
      imputation_log[[length(imputation_log) + 1L]] <- data.frame(
        round = mi,
        variable = v,
        server = srv,
        n_imputed = as.integer(imp_res$n_imputed %||% NA_integer_),
        n_observed = as.integer(imp_res$n_observed %||% NA_integer_),
        method = as.character(imp_res$method %||% NA_character_),
        model_regularization = as.character(
          imp_res$model_regularization %||% NA_character_),
        regularization_alpha = as.numeric(
          imp_res$regularization_alpha %||% NA_real_),
        numerical_stabilization = as.character(
          imp_res$numerical_stabilization %||% NA_character_),
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
  if (is.numeric(lambda) && length(lambda) == 1L && is.finite(lambda) &&
      lambda > 0) {
    quality_warnings <- c(
      quality_warnings,
      "A positive lambda was explicitly requested; pooled coefficients target a ridge-penalized imputation-analysis model, not the unpenalized GLM estimand.")
    quality_status <- "approximate"
  }
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
    lambda       = as.numeric(lambda),
    estimand     = if (isTRUE(as.numeric(lambda) > 0)) {
      "explicit_ridge_penalized_mi_glm"
    } else {
      "unpenalized_mi_glm"
    },
    intercept_only_policy = intercept_only,
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
  if (!is.null(x$estimand)) {
    cat(sprintf("  Estimand: %s (lambda = %.4g)\n",
                x$estimand, x$lambda %||% 0))
  }
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
