#' @title Federated generalised estimating equations
#' @description Fit a GLM/GEE for vertically partitioned DataSHIELD data and
#'   return sandwich (robust) standard errors alongside the usual
#'   model-based ones. For \code{corstr = "independence"}, the point estimate
#'   \eqn{\hat{\beta}} is obtained by a single call to
#'   \code{\link{ds.vertGLM}}. For Gaussian \code{corstr = "exchangeable"},
#'   dsVert promotes the point estimate to a protected cluster-level
#'   exchangeable GLS/GEE update. For binomial and Poisson
#'   \code{corstr = "exchangeable"}, dsVert uses a Ring127 protected
#'   Pearson-score exchangeable GEE update. For Gaussian, binomial, and
#'   Poisson \code{corstr = "ar1"}, dsVert uses guarded within-cluster order
#'   metadata on the DCF parties and returns only low-dimensional AR1
#'   aggregates.
#'   For Gaussian, binomial, and Poisson independence models the sandwich meat
#'   is computed in the share domain. If
#'   \code{id_col} is supplied, dsVert computes the clustered meat
#'   \eqn{\sum_c S_c S_c^\top}, where
#'   \eqn{S_c=\sum_{i\in c} X_i r_i}; otherwise it computes the row-level
#'   HC0 meat \eqn{X^\top \mathrm{diag}(r^2) X}. Cluster membership is
#'   transport-encrypted between DCF parties and clusters below
#'   \code{datashield.privacyLevel} fail closed.
#'
#'   The client receives only low-dimensional aggregate matrices/vectors.
#'   It never materialises the \eqn{n}-length residual, squared-residual, or
#'   weighted-column vectors.
#'
#'   Formula:
#'     \deqn{V_{sand} = \mathrm{Cov}(\hat\beta) \, A \, \mathrm{Cov}(\hat\beta)}
#'     \deqn{A = \sum_c S_c S_c^\top}
#'   for clustered Gaussian/binomial/Poisson fits, or
#'     \deqn{A = X^T \, \mathrm{diag}(r^2) \, X}
#'   for row-level HC0 when no \code{id_col} is supplied.
#'
#' @param formula A model formula passed through to \code{ds.vertGLM}.
#' @param data Character. Aligned data-frame name on each server.
#' @param family One of \code{"gaussian"}, \code{"binomial"},
#'   \code{"poisson"}.
#' @param id_col Optional character. For Gaussian/binomial/Poisson models this
#'   enables the cluster-robust sandwich meat; the cluster column must live with
#'   the outcome and all clusters must pass \code{datashield.privacyLevel}.
#' @param order_col Optional character. Required for \code{corstr = "ar1"}.
#'   The order column must live with the outcome and is used only
#'   server-to-server to derive guarded adjacent records.
#' @param corstr Working correlation. \code{"independence"} is available for
#'   Gaussian/binomial/Poisson. \code{"exchangeable"} currently fits true
#'   exchangeable Gaussian, binomial, and Poisson GEE coefficients from guarded
#'   cluster-level sufficient statistics. \code{"ar1"} fits Gaussian,
#'   binomial, and Poisson GEE from guarded adjacent-pair sufficient
#'   statistics.
#' @param max_iter,tol,lambda,verbose Passed to \code{ds.vertGLM}.
#' @param working_max_iter Optional integer. Maximum iterations for
#'   exchangeable/AR1 working-correlation updates. Defaults to \code{max_iter}.
#' @param ring Integer 63 or 127. Binomial/Poisson exchangeable and all AR1
#'   routes are automatically run in Ring127 because protected nonlinear link
#'   operations and adjacent-product statistics need high precision.
#' @param binomial_sigmoid_intervals Optional integer. Number of DCF spline
#'   intervals for protected binomial sigmoid evaluations used by the
#'   underlying GLM fit and GEE sandwich/working-correlation updates. When
#'   \code{NULL}, \code{dsvert.gee_binomial_sigmoid_intervals},
#'   \code{dsvert.glm_num_intervals_binomial}, or the 100-interval default is
#'   used.
#' @param datasources DataSHIELD connection object.
#' @return An object of class \code{ds.vertGEE} with components
#'   \code{coefficients}, \code{model_se} (sqrt of
#'   \code{diag(Cov(beta))}), \code{robust_se} (sandwich SEs),
#'   \code{robust_covariance}, \code{corstr}, and \code{fit}
#'   (the underlying \code{ds.glm} object).
#' @export
ds.vertGEE <- function(formula, data = NULL,
                       family = c("gaussian", "binomial", "poisson"),
                       id_col = NULL,
                       order_col = NULL,
                       corstr = c("independence", "exchangeable", "ar1"),
                       max_iter = 100L, tol = 1e-4, lambda = 1e-4,
                       working_max_iter = NULL,
                       ring = 63L,
                       binomial_sigmoid_intervals = NULL,
                       verbose = TRUE, datasources = NULL) {
  family <- match.arg(family)
  corstr <- match.arg(corstr)
  if (!is.null(binomial_sigmoid_intervals)) {
    binomial_sigmoid_intervals <- as.integer(binomial_sigmoid_intervals)
    if (length(binomial_sigmoid_intervals) != 1L ||
        !is.finite(binomial_sigmoid_intervals) ||
        binomial_sigmoid_intervals < 10L) {
      stop("binomial_sigmoid_intervals must be NULL or an integer >= 10",
           call. = FALSE)
    }
  }
  effective_binomial_sigmoid_intervals <- if (identical(family, "binomial")) {
    if (!is.null(binomial_sigmoid_intervals)) {
      binomial_sigmoid_intervals
    } else {
      suppressWarnings(as.integer(getOption(
        "dsvert.gee_binomial_sigmoid_intervals",
        getOption("dsvert.glm_num_intervals_binomial",
                  getOption("dsvert.glm_num_intervals", 100L)))[[1L]]))
    }
  } else {
    NA_integer_
  }
  if (identical(family, "binomial") &&
      (!is.finite(effective_binomial_sigmoid_intervals) ||
       effective_binomial_sigmoid_intervals < 10L)) {
    effective_binomial_sigmoid_intervals <- 100L
  }
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)
  if (identical(corstr, "exchangeable") &&
      !(family %in% c("gaussian", "binomial", "poisson"))) {
    stop("ds.vertGEE corstr='exchangeable' is currently implemented only ",
         "for family='gaussian', family='binomial', and family='poisson'. Use ",
         "corstr='independence' for ",
         family, " until the non-Gaussian estimating equations are added.",
         call. = FALSE)
  }
  # Cluster-ID broadcast is used only inside the share-domain
  # sandwich path below. The analyst client relays opaque encrypted blobs
  # and does not receive row-level cluster membership.
  if (corstr %in% c("exchangeable", "ar1")) {
    if (is.null(id_col)) {
      stop("corstr='", corstr, "' requires id_col", call. = FALSE)
    }
    if (identical(corstr, "ar1") &&
        (is.null(order_col) || !is.character(order_col) ||
         length(order_col) != 1L || !nzchar(order_col))) {
      stop("corstr='ar1' requires order_col", call. = FALSE)
    }
    if (verbose) {
      message("[ds.vertGEE] corstr='", corstr,
              "' - using protected working-correlation path")
    }
  }
  if (verbose && !is.null(id_col)) {
    if (family %in% c("gaussian", "binomial", "poisson")) {
      message("[ds.vertGEE] id_col='", id_col,
              "' enables share-domain cluster sandwich")
    } else {
      message("[ds.vertGEE] id_col='", id_col,
              "' recorded; cluster sandwich currently implemented for ",
              "Gaussian/binomial/Poisson models only")
    }
  }

  ring_use <- as.integer(ring %||% 63L)
  if (!(ring_use %in% c(63L, 127L))) {
    stop("ring must be 63 or 127", call. = FALSE)
  }
  working_max_iter_use <- as.integer(working_max_iter %||% max_iter)
  if (!is.finite(working_max_iter_use) || working_max_iter_use < 1L) {
    stop("working_max_iter must be a positive integer", call. = FALSE)
  }
  if (identical(corstr, "exchangeable") &&
      family %in% c("binomial", "poisson") &&
      ring_use != 127L) {
    if (verbose) {
      message("[ds.vertGEE] ", family, " exchangeable promoted to Ring127")
    }
    ring_use <- 127L
  }
  if (identical(corstr, "ar1") && ring_use != 127L) {
    if (verbose) {
      message("[ds.vertGEE] AR1 promoted to Ring127")
    }
    ring_use <- 127L
  }

  # Stage 1: ordinary fit -> betahat, Cov_model (bread).
  if (verbose) message("[ds.vertGEE] Stage 1: ds.vertGLM point estimate")
  fit <- ds.vertGLM(formula = formula, data = data,
                    family = family, max_iter = max_iter, tol = tol,
                    lambda = lambda, ring = ring_use, verbose = verbose,
                    datasources = datasources,
                    binomial_sigmoid_intervals =
                      if (identical(family, "binomial")) {
                        effective_binomial_sigmoid_intervals
                      } else {
                        NULL
                      },
                    keep_session = TRUE)
  if (!is.null(fit$session_id)) {
    on.exit({
      for (.srv in fit$server_list %||% names(datasources)) {
        .ci <- which(names(datasources) == .srv)
        if (length(.ci) == 1L) {
          tryCatch(DSI::datashield.aggregate(datasources[.ci],
            call(name = "mpcCleanupDS", session_id = fit$session_id)),
            error = function(e) NULL)
          tryCatch(DSI::datashield.aggregate(datasources[.ci],
            call(name = "mpcGcDS")), error = function(e) NULL)
        }
      }
    }, add = TRUE)
  }
  if (is.null(fit$covariance)) {
    stop("ds.vertGEE: the underlying ds.vertGLM fit does not expose ",
         "Cov(beta). Refit with dsVert >= 8bb7902 / dsVertClient >= ",
         "8bb7902.", call. = FALSE)
  }
  Cov_model <- as.matrix(fit$covariance)
  Cov_model <- (Cov_model + t(Cov_model)) / 2
  Cov_model_info <- as.matrix(fit$covariance_information %||%
                                fit$covariance_unscaled %||%
                                fit$covariance)
  Cov_model_info <- (Cov_model_info + t(Cov_model_info)) / 2
  coef_out <- fit$coefficients
  model_se_out <- fit$std_errors
  working_correlation <- list(corstr = corstr, alpha = NA_real_,
                              estimator = "independence")

  # Stage 2: share-domain sandwich meat. Gaussian/binomial/Poisson models use
  # the retained GLM session so residuals and X columns remain additive shares.
  # A legacy weighted-r2 fallback is kept for older servers only.
  {
    if (verbose) message("[ds.vertGEE] Stage 2: weighted fit for sandwich meat (family=",
                          family, ")")
    # Request residuals-squared column on the outcome server. This uses
    # the `dsvertPearsonR2ColDS` server helper (Gaussian-only in the
    # first pass): given data_name, y, x_vars, and a plaintext betahat,
    # it writes a new column `__dsvert_r2` into the aligned data frame.
    # The client then re-fits with weights="__dsvert_r2".
    # Discover y_server
    y_var <- .ds_gee_extract_lhs(formula)
    y_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                          data, y_var)
    robust_method <- "model_based"
    if (is.null(y_srv)) {
      warning("ds.vertGEE: could not locate outcome server for y='",
              y_var, "'. Returning model-based SE.", call. = FALSE)
      robust_se <- fit$std_errors
      robust_cov <- Cov_model
    } else {
      exchangeable_fit <- NULL
      exchangeable_error <- NULL
      if (identical(corstr, "exchangeable")) {
        exchangeable_fit <- tryCatch(
          if (identical(family, "gaussian")) {
            .ds_gee_secure_gaussian_exchangeable(
              fit = fit, datasources = datasources,
              server_names = server_names, lambda = lambda,
              data = data, id_col = id_col,
              max_iter = working_max_iter_use,
              tol = tol, verbose = verbose)
          } else if (family %in% c("binomial", "poisson")) {
            .ds_gee_secure_poisson_exchangeable(
              fit = fit, datasources = datasources,
              server_names = server_names, lambda = lambda,
              data = data, id_col = id_col,
              max_iter = working_max_iter_use,
              tol = tol, verbose = verbose,
              binomial_sigmoid_intervals =
                effective_binomial_sigmoid_intervals)
          } else {
            stop("unsupported exchangeable GEE family: ", family,
                 call. = FALSE)
          },
          error = function(e) {
            exchangeable_error <<- e
            if (verbose) {
              message("[ds.vertGEE] exchangeable path unavailable: ",
                      conditionMessage(e))
            }
            NULL
          })
      }
      if (!is.null(exchangeable_fit)) {
        working_fit <- exchangeable_fit
      } else {
        working_fit <- NULL
      }
      ar1_error <- NULL
      if (identical(corstr, "ar1")) {
        working_fit <- tryCatch(
          if (identical(family, "gaussian")) {
            .ds_gee_secure_gaussian_ar1(
              fit = fit, datasources = datasources,
              server_names = server_names, lambda = lambda,
              data = data, id_col = id_col, order_col = order_col,
              max_iter = working_max_iter_use, tol = tol, verbose = verbose)
          } else if (family %in% c("binomial", "poisson")) {
            .ds_gee_secure_poisson_exchangeable(
              fit = fit, datasources = datasources,
              server_names = server_names, lambda = lambda,
              data = data, id_col = id_col, order_col = order_col,
              corstr = "ar1", max_iter = working_max_iter_use,
              tol = tol, verbose = verbose,
              binomial_sigmoid_intervals =
                effective_binomial_sigmoid_intervals)
          } else {
            stop("unsupported AR1 GEE family: ", family, call. = FALSE)
          },
          error = function(e) {
            ar1_error <<- e
            if (verbose) {
              message("[ds.vertGEE] AR1 path unavailable: ",
                      conditionMessage(e))
            }
            NULL
          })
      }
      if (!is.null(working_fit)) {
        coef_out <- working_fit$coefficients
        model_se_out <- working_fit$model_se
        robust_cov <- working_fit$robust_covariance
        robust_se <- working_fit$robust_se
        robust_method <- working_fit$method
        working_correlation <- working_fit$working_correlation
      } else {
        if (identical(corstr, "exchangeable")) {
          msg <- if (is.null(exchangeable_error)) {
            "unknown error"
          } else {
            conditionMessage(exchangeable_error)
          }
          stop("exchangeable GEE unavailable: ", msg,
               call. = FALSE)
        }
        if (identical(corstr, "ar1")) {
          msg <- if (is.null(ar1_error)) {
            "unknown error"
          } else {
            conditionMessage(ar1_error)
          }
          stop("AR1 GEE unavailable: ", msg, call. = FALSE)
        }
        secure_hc0 <- NULL
        secure_error <- NULL
        if (family %in% c("gaussian", "binomial", "poisson")) {
          secure_hc0 <- tryCatch(
            .ds_gee_secure_hc0(
              fit = fit, datasources = datasources,
              server_names = server_names, lambda = lambda,
              data = data, id_col = id_col,
              verbose = verbose,
              binomial_sigmoid_intervals =
                effective_binomial_sigmoid_intervals),
            error = function(e) {
              secure_error <<- e
              if (verbose) {
                message("[ds.vertGEE] secure sandwich path unavailable: ",
                        conditionMessage(e))
              }
              NULL
            })
        }
        if (!is.null(secure_hc0)) {
          robust_cov <- secure_hc0$robust_covariance
          robust_se <- secure_hc0$robust_se
          robust_method <- secure_hc0$method
        } else {
        if (family %in% c("gaussian", "binomial", "poisson") &&
            !is.null(id_col) &&
            is.character(id_col) && length(id_col) == 1L && nzchar(id_col)) {
          msg <- if (is.null(secure_error)) {
            "unknown error"
          } else {
            conditionMessage(secure_error)
          }
          stop("clustered GEE sandwich unavailable: ", msg, call. = FALSE)
        }
        x_all <- names(fit$coefficients)
        x_all <- setdiff(x_all, "(Intercept)")
        r2_ok <- tryCatch({
          DSI::datashield.aggregate(
            datasources[which(server_names == y_srv)],
            call(name = "dsvertPearsonR2ColDS",
                 data_name = data, y_var = y_var,
                 x_names = x_all,
                 betahat = as.numeric(fit$coefficients[x_all]),
                 intercept = as.numeric(fit$coefficients["(Intercept)"]),
                 family = family,
                 r2_column = "__dsvert_r2"))
          TRUE
        }, error = function(e) {
          if (verbose) message("[ds.vertGEE] dsvertPearsonR2ColDS not ",
                               "available yet: ", conditionMessage(e),
                               "\n  Returning model-based SE.")
          FALSE
        })
        if (!isTRUE(r2_ok)) {
          robust_se <- fit$std_errors
          robust_cov <- Cov_model
        } else {
        fit_w <- ds.vertGLM(formula = formula, data = data,
                             family = family, max_iter = max_iter,
                             tol = tol, lambda = lambda,
                             weights = "__dsvert_r2",
                             verbose = verbose,
                             datasources = datasources)
        if (is.null(fit_w$covariance)) {
          warning("weighted fit did not expose covariance; returning ",
                  "model-based SE.", call. = FALSE)
          robust_se <- fit$std_errors
          robust_cov <- Cov_model
        } else {
          # A ~ inverse-Fisher of weighted fit (meat). Sandwich:
          #   V_sand = Bread * Fisher_weighted * Bread
          #          ~= Cov_info_model * solve(Cov_info_weighted) *
          #              Cov_info_model
          # where Cov_info_* is inverse Fisher, not a Gaussian sigma^2-scaled
          # model covariance.
          Cov_w_info <- as.matrix(fit_w$covariance_information %||%
                                     fit_w$covariance_unscaled %||%
                                     fit_w$covariance)
          Cov_w_info <- (Cov_w_info + t(Cov_w_info)) / 2
          A <- tryCatch(solve(Cov_w_info), error = function(e) NULL)
          if (is.null(A)) {
            warning("meat matrix singular; returning model-based SE.",
                    call. = FALSE)
            robust_se <- fit$std_errors
            robust_cov <- Cov_model
          } else {
            robust_cov <- Cov_model_info %*% A %*% Cov_model_info
            dimnames(robust_cov) <- list(names(fit$coefficients),
                                          names(fit$coefficients))
            robust_se <- sqrt(pmax(diag(robust_cov), 0))
            names(robust_se) <- names(fit$coefficients)
            robust_method <- "weighted_r2_fallback"
          }
        }
        }
        }
      }
    }
  }

  out <- list(
    coefficients       = coef_out,
    model_se           = model_se_out,
    robust_se          = robust_se,
    robust_covariance  = robust_cov,
    robust_method      = robust_method,
    corstr             = corstr,
    working_correlation = working_correlation,
    family             = family,
    binomial_sigmoid_intervals =
      if (identical(family, "binomial")) {
        effective_binomial_sigmoid_intervals
      } else {
        NA_integer_
      },
    id_col             = id_col,
    order_col          = if (identical(corstr, "ar1")) order_col else NULL,
    n_obs              = fit$n_obs,
    fit                = fit,
    call               = match.call())
  class(out) <- c("ds.vertGEE", "list")
  out
}

#' @keywords internal
.ds_gee_extract_lhs <- function(formula) {
  if (inherits(formula, "formula")) {
    as.character(attr(terms(formula), "variables")[[2]])
  } else if (is.character(formula) && grepl("~", formula)) {
    as.character(attr(terms(as.formula(formula)), "variables")[[2]])
  } else stop("cannot extract y from formula", call. = FALSE)
}

#' @keywords internal
.ds_gee_find_server_holding <- function(datasources, server_names,
                                        data_name, var) {
  for (srv in server_names) {
    ci <- which(server_names == srv)
    cols <- tryCatch(
      DSI::datashield.aggregate(datasources[ci],
        call(name = "dsvertColNamesDS", data_name = data_name))[[1]]$columns,
      error = function(e) NULL)
    if (!is.null(cols) && var %in% cols) return(srv)
  }
  NULL
}

#' @keywords internal
.ds_gee_standardized_parameters <- function(fit, features) {
  beta_orig <- as.numeric(fit$coefficients[features])
  names(beta_orig) <- features
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])
  names(x_sds) <- names(x_means) <- features

  if (fit$family == "gaussian" && isTRUE(fit$standardize_y) &&
      !is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
    beta_std <- beta_orig * x_sds / fit$y_sd
    intercept_std <- (as.numeric(fit$coefficients["(Intercept)"]) +
                        sum(beta_orig * x_means) - fit$y_mean) / fit$y_sd
  } else {
    beta_std <- beta_orig * x_sds
    intercept_std <- as.numeric(fit$coefficients["(Intercept)"]) +
                       sum(beta_orig * x_means)
  }
  names(beta_std) <- features
  list(intercept = as.numeric(intercept_std), beta = beta_std)
}

#' @keywords internal
.ds_gee_standardization_jacobian <- function(fit, features) {
  p <- length(features)
  out_names <- c("(Intercept)", features)
  J <- diag(p + 1L)
  dimnames(J) <- list(out_names, out_names)
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])

  if (fit$family == "gaussian" && isTRUE(fit$standardize_y) &&
      !is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
    for (jj in seq_len(p)) {
      J[jj + 1L, jj + 1L] <- fit$y_sd / x_sds[jj]
      J[1L, jj + 1L] <- -fit$y_sd * x_means[jj] / x_sds[jj]
    }
    J[1L, 1L] <- fit$y_sd
  } else {
    for (jj in seq_len(p)) {
      J[jj + 1L, jj + 1L] <- 1 / x_sds[jj]
      J[1L, jj + 1L] <- -x_means[jj] / x_sds[jj]
    }
  }
  J
}

#' @keywords internal
.ds_gee_unstandardized_parameters <- function(fit, beta_std, features) {
  beta_std <- as.numeric(beta_std[c("(Intercept)", features)])
  names(beta_std) <- c("(Intercept)", features)
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])
  names(x_sds) <- names(x_means) <- features

  if (fit$family == "gaussian" && isTRUE(fit$standardize_y) &&
      !is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
    beta_orig <- beta_std[features] * fit$y_sd / x_sds
    intercept <- beta_std["(Intercept)"] * fit$y_sd + fit$y_mean -
      sum(beta_orig * x_means)
  } else {
    beta_orig <- beta_std[features] / x_sds
    intercept <- beta_std["(Intercept)"] -
      sum(beta_orig * x_means)
  }
  out <- c("(Intercept)" = as.numeric(intercept), beta_orig)
  names(out) <- c("(Intercept)", features)
  out
}

#' @keywords internal
.ds_gee_solve_exchangeable_from_cluster_stats <- function(
    sx, sy, sxy, xx, syy, beta_start, max_iter = 100L, tol = 1e-6) {
  q <- ncol(sx)
  n_clusters <- nrow(sx)
  cluster_sizes <- as.numeric(sx[, 1L])
  if (q < 1L || n_clusters < 1L) {
    stop("empty cluster statistics", call. = FALSE)
  }
  if (any(!is.finite(cluster_sizes)) || any(cluster_sizes < 2)) {
    stop("invalid exchangeable cluster sizes", call. = FALSE)
  }
  n_obs <- sum(cluster_sizes)
  beta <- as.numeric(beta_start)
  if (length(beta) != q || any(!is.finite(beta))) {
    beta <- rep(0, q)
  }
  alpha <- 0
  lower_alpha <- -1 / (max(cluster_sizes) - 1) + 1e-6
  upper_alpha <- 0.95

  safe_solve <- function(A, b = NULL) {
    A <- (A + t(A)) / 2
    out <- tryCatch(if (is.null(b)) solve(A) else solve(A, b),
                    error = function(e) NULL)
    if (!is.null(out) && all(is.finite(out))) return(out)
    ridge <- 1e-8 * max(1, mean(abs(diag(A))))
    tryCatch(if (is.null(b)) solve(A + diag(ridge, nrow(A)))
             else solve(A + diag(ridge, nrow(A)), b),
             error = function(e) {
               stop("exchangeable GLS system is singular", call. = FALSE)
             })
  }

  residual_summaries <- function(beta) {
    rss <- numeric(n_clusters)
    rsum <- numeric(n_clusters)
    for (cc in seq_len(n_clusters)) {
      xxy <- xx[cc, , , drop = FALSE][1L, , ]
      rss[cc] <- syy[cc] - 2 * sum(beta * sxy[cc, ]) +
        drop(t(beta) %*% xxy %*% beta)
      rsum[cc] <- sy[cc] - sum(sx[cc, ] * beta)
    }
    rss <- pmax(rss, 0)
    list(rss = rss, rsum = rsum)
  }

  update_alpha <- function(beta) {
    rs <- residual_summaries(beta)
    phi <- sum(rs$rss) / max(n_obs - q, 1)
    if (!is.finite(phi) || phi <= 0) return(0)
    pair_num <- sum((rs$rsum^2 - rs$rss) / 2)
    pair_den <- sum(cluster_sizes * (cluster_sizes - 1) / 2) * phi
    if (!is.finite(pair_den) || pair_den <= 0) return(0)
    alpha_new <- pair_num / pair_den
    min(max(alpha_new, lower_alpha), upper_alpha)
  }

  fixed_alpha_fit <- function(alpha) {
    A <- matrix(0, q, q)
    rhs <- numeric(q)
    for (cc in seq_len(n_clusters)) {
      m <- cluster_sizes[cc]
      a <- 1 / (1 - alpha)
      b <- alpha / ((1 - alpha) * (1 - alpha + alpha * m))
      A <- A + a * xx[cc, , ] - b * tcrossprod(sx[cc, ])
      rhs <- rhs + a * sxy[cc, ] - b * sx[cc, ] * sy[cc]
    }
    list(beta = as.numeric(safe_solve(A, rhs)), A = A)
  }

  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(as.integer(max_iter))) {
    beta_old <- beta
    alpha_old <- alpha
    alpha <- update_alpha(beta)
    fa <- fixed_alpha_fit(alpha)
    beta <- fa$beta
    if (max(abs(beta - beta_old), abs(alpha - alpha_old)) <= tol) {
      converged <- TRUE
      break
    }
  }
  fa <- fixed_alpha_fit(alpha)
  beta <- fa$beta
  A <- (fa$A + t(fa$A)) / 2
  bread <- safe_solve(A)

  rs <- residual_summaries(beta)
  phi_num <- 0
  scores <- matrix(0, nrow = n_clusters, ncol = q)
  for (cc in seq_len(n_clusters)) {
    m <- cluster_sizes[cc]
    a <- 1 / (1 - alpha)
    b <- alpha / ((1 - alpha) * (1 - alpha + alpha * m))
    rsum <- rs$rsum[cc]
    phi_num <- phi_num + a * rs$rss[cc] - b * rsum^2
    scores[cc, ] <- a * (sxy[cc, ] - xx[cc, , ] %*% beta) -
      b * sx[cc, ] * rsum
  }
  phi <- phi_num / max(n_obs - q, 1)
  if (!is.finite(phi) || phi <= 0) {
    phi <- sum(rs$rss) / max(n_obs - q, 1)
  }
  model_cov <- as.numeric(phi) * bread
  robust_cov <- bread %*% crossprod(scores) %*% bread
  model_cov <- (model_cov + t(model_cov)) / 2
  robust_cov <- (robust_cov + t(robust_cov)) / 2

  list(beta = beta, alpha = alpha, phi = as.numeric(phi),
       model_cov = model_cov, robust_cov = robust_cov,
       iterations = as.integer(iter), converged = converged,
       scores = scores)
}

#' @keywords internal
.ds_gee_residual_quad_from_stats <- function(stats, beta) {
  stats$yy - 2 * sum(beta * stats$xy) +
    drop(t(beta) %*% stats$xx %*% beta)
}

#' @keywords internal
.ds_gee_adjacent_residual_cross <- function(adj, beta) {
  adj$yy - sum(beta * adj$xy_forward) - sum(beta * adj$xy_backward) +
    drop(t(beta) %*% adj$xx %*% beta)
}

#' @keywords internal
.ds_gee_solve_ar1_from_stats <- function(
    total, interior, nonlast, adj, beta_start, max_iter = 100L,
    tol = 1e-6, lag_stats = NULL) {
  q <- ncol(total$xx)
  n_obs <- as.numeric(total$xx[1L, 1L])
  if (q < 1L || !is.finite(n_obs) || n_obs < q + 1L) {
    stop("invalid AR1 sufficient statistics", call. = FALSE)
  }
  beta <- as.numeric(beta_start)
  if (length(beta) != q || any(!is.finite(beta))) beta <- rep(0, q)
  rho <- 0
  lower_rho <- -0.95
  upper_rho <- 0.95

  safe_solve <- function(A, b = NULL) {
    A <- (A + t(A)) / 2
    out <- tryCatch(if (is.null(b)) solve(A) else solve(A, b),
                    error = function(e) NULL)
    if (!is.null(out) && all(is.finite(out))) return(out)
    ridge <- 1e-8 * max(1, mean(abs(diag(A))))
    tryCatch(if (is.null(b)) solve(A + diag(ridge, nrow(A)))
             else solve(A + diag(ridge, nrow(A)), b),
             error = function(e) {
               stop("AR1 GLS system is singular", call. = FALSE)
             })
  }

  if (is.null(lag_stats) || length(lag_stats) == 0L) {
    lag_stats <- list(adj)
  }

  update_rho <- function(beta_current) {
    if (length(lag_stats) == 1L) {
      num <- .ds_gee_adjacent_residual_cross(adj, beta_current)
      den <- .ds_gee_residual_quad_from_stats(nonlast, beta_current)
      if (!is.finite(num) || !is.finite(den) || den <= 0) return(0)
      return(min(max(num / den, lower_rho), upper_rho))
    }
    phi <- .ds_gee_residual_quad_from_stats(total, beta_current) /
      max(n_obs - q, 1)
    if (!is.finite(phi) || phi <= 0) return(0)
    counts <- vapply(lag_stats, function(x) as.numeric(x$count %||%
                                                        x$xx[1L, 1L]),
                     numeric(1L))
    rho_targets <- vapply(lag_stats, function(x) {
      .ds_gee_adjacent_residual_cross(x, beta_current) /
        max(as.numeric(x$count %||% x$xx[1L, 1L]) * phi, .Machine$double.eps)
    }, numeric(1L))
    keep <- is.finite(counts) & counts > 0 & is.finite(rho_targets)
    if (!any(keep)) return(0)
    counts <- counts[keep]
    rho_targets <- rho_targets[keep]
    lags <- seq_along(lag_stats)[keep]
    objective <- function(rho_current) {
      sum(counts * (rho_targets - rho_current^lags)^2)
    }
    out <- tryCatch(stats::optimize(objective,
                                    interval = c(lower_rho, upper_rho))$minimum,
                    error = function(e) 0)
    if (!is.finite(out)) out <- 0
    min(max(out, lower_rho), upper_rho)
  }

  fixed_rho_fit <- function(rho_current) {
    cfac <- 1 / (1 - rho_current^2)
    A <- cfac * (total$xx + rho_current^2 * interior$xx -
                   rho_current * (adj$xx + t(adj$xx)))
    rhs <- cfac * (total$xy + rho_current^2 * interior$xy -
                     rho_current * (adj$xy_forward + adj$xy_backward))
    list(beta = as.numeric(safe_solve(A, rhs)), A = A)
  }

  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(as.integer(max_iter))) {
    beta_old <- beta
    rho_old <- rho
    rho <- update_rho(beta)
    fr <- fixed_rho_fit(rho)
    beta <- fr$beta
    if (max(abs(beta - beta_old), abs(rho - rho_old)) <= tol) {
      converged <- TRUE
      break
    }
  }
  fr <- fixed_rho_fit(rho)
  beta <- fr$beta
  A <- (fr$A + t(fr$A)) / 2
  bread <- safe_solve(A)

  rss_total <- .ds_gee_residual_quad_from_stats(total, beta)
  phi <- rss_total / max(n_obs - q, 1)
  if (!is.finite(phi) || phi <= 0) phi <- 1
  model_cov <- as.numeric(phi) * bread
  model_cov <- (model_cov + t(model_cov)) / 2
  list(beta = beta, rho = rho, phi = as.numeric(phi), A = A,
       bread = bread, model_cov = model_cov,
       iterations = as.integer(iter), converged = converged)
}

#' @keywords internal
.ds_gee_secure_gaussian_ar1 <- function(
    fit, datasources, server_names, lambda = 0, data = NULL, id_col = NULL,
    order_col = NULL, max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  if (!identical(fit$family, "gaussian")) return(NULL)
  if (is.null(id_col) || length(id_col) != 1L ||
      !is.character(id_col) || !nzchar(id_col)) {
    stop("Gaussian AR1 GEE requires id_col", call. = FALSE)
  }
  if (is.null(order_col) || length(order_col) != 1L ||
      !is.character(order_col) || !nzchar(order_col)) {
    stop("Gaussian AR1 GEE requires order_col", call. = FALSE)
  }
  required <- c("session_id", "transport_pks", "server_list", "x_vars",
                "y_server", "x_sds", "x_means")
  missing_req <- required[vapply(required, function(nm) is.null(fit[[nm]]),
                                 logical(1L))]
  if (length(missing_req) > 0L) {
    stop("GLM session metadata missing: ", paste(missing_req, collapse = ", "),
         call. = FALSE)
  }

  session_id <- fit$session_id
  server_list <- fit$server_list
  x_vars <- fit$x_vars
  coordinator <- fit$y_server
  n_obs <- as.integer(fit$n_obs)
  ring <- as.integer(fit$ring %||% 127L)
  if (ring != 127L) {
    stop("Gaussian AR1 GEE requires Ring127", call. = FALSE)
  }
  frac_bits <- 50L
  ring_tag <- "ring127"
  transport_pks <- fit$transport_pks

  .to_b64url <- function(x) gsub("+", "-", gsub("/", "_",
    gsub("=+$", "", x, perl = TRUE), fixed = TRUE), fixed = TRUE)
  .fp_const <- function(x) {
    .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
      values = array(as.numeric(x), dim = 1L),
      frac_bits = frac_bits, ring = ring_tag))$fp_data)
  }
  .dsAgg <- function(conns, expr, ...) {
    tryCatch(
      DSI::datashield.aggregate(conns = conns, expr = expr, ...),
      error = function(e) {
        fn_name <- if (is.call(expr)) as.character(expr[[1]]) else "?"
        srv_name <- tryCatch(names(conns)[1], error = function(x) "?")
        stop(sprintf("%s on %s failed: %s", fn_name, srv_name,
                     conditionMessage(e)), call. = FALSE)
      })
  }
  .sendBlob <- function(blob, key, conn_idx) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        .dsAgg(datasources[conn_idx],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               session_id = session_id))
      } else {
        .dsAgg(datasources[conn_idx],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               chunk_index = chunk_idx, n_chunks = n_chunks,
               session_id = session_id))
      }
    })
  }

  if (fit$eta_privacy == "k2_beaver") {
    nl <- setdiff(server_list, coordinator)
    if (length(nl) != 1L) {
      stop("K=2 Gaussian AR1 GEE requires exactly one non-outcome server",
           call. = FALSE)
    }
    nl <- nl[[1L]]
    dcf_parties <- c(coordinator, nl)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    dealer_conn <- dcf_conns[[2L]]
    canonical_features <- c(x_vars[[coordinator]], x_vars[[nl]])
  } else if (fit$eta_privacy == "secure_agg") {
    fusion <- .k3_select_fusion_server(server_list, coordinator, x_vars)
    dcf_parties <- c(fusion, coordinator)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    non_dcf <- setdiff(server_list, dcf_parties)
    dealer <- if (length(non_dcf) > 0L) non_dcf[[1L]] else fusion
    dealer_conn <- which(server_names == dealer)
    canonical_features <- c(x_vars[[coordinator]], x_vars[[fusion]])
    for (srv in non_dcf) canonical_features <- c(canonical_features, x_vars[[srv]])
  } else {
    stop("unsupported eta_privacy for Gaussian AR1 GEE: ",
         fit$eta_privacy, call. = FALSE)
  }

  cluster_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                             data, id_col)
  if (is.null(cluster_srv)) {
    stop("id_col '", id_col, "' not found on any server", call. = FALSE)
  }
  order_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                           data, order_col)
  if (is.null(order_srv)) {
    stop("order_col '", order_col, "' not found on any server",
         call. = FALSE)
  }
  if (!identical(cluster_srv, coordinator) ||
      !identical(order_srv, coordinator)) {
    stop("Gaussian AR1 GEE requires id_col and order_col to live on the ",
         "outcome server", call. = FALSE)
  }

  target_features <- unlist(x_vars[server_list], use.names = FALSE)
  target_order <- c("(Intercept)", target_features)
  canonical_order <- c("(Intercept)", canonical_features)
  p_total <- length(canonical_features)
  q <- p_total + 1L
  std <- .ds_gee_standardized_parameters(fit, target_features)

  .restore_eta_shape <- function(beta_std) {
    if (fit$eta_privacy == "k2_beaver") {
      b_coord <- as.numeric(beta_std[x_vars[[coordinator]]])
      b_nl <- as.numeric(beta_std[x_vars[[dcf_parties[[2L]]]]])
      for (i in seq_along(dcf_parties)) {
        srv <- dcf_parties[[i]]
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "dsvertGEERestoreFeatureShapeDS",
               p_own = as.integer(length(x_vars[[srv]])),
               p_peer = as.integer(length(x_vars[[setdiff(dcf_parties, srv)]])),
               session_id = session_id))
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "k2ComputeEtaShareDS",
               beta_coord = b_coord, beta_nl = b_nl,
               intercept = if (srv == coordinator) {
                 as.numeric(beta_std["(Intercept)"])
               } else {
                 0
               },
               is_coordinator = (srv == coordinator),
               output_key = "gee_ar1_eta",
               session_id = session_id))
      }
    } else {
      fusion <- dcf_parties[[1L]]
      non_dcf <- setdiff(server_list, dcf_parties)
      p_coord <- length(x_vars[[coordinator]])
      p_fusion <- length(x_vars[[fusion]])
      p_extras <- sum(vapply(non_dcf, function(s) length(x_vars[[s]]),
                             integer(1L)))
      for (i in seq_along(dcf_parties)) {
        srv <- dcf_parties[[i]]
        is_coord <- srv == coordinator
        p_peer <- if (is_coord) p_fusion + p_extras else p_coord + p_extras
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "dsvertGEERestoreFeatureShapeDS",
               p_own = as.integer(length(x_vars[[srv]])),
               p_peer = as.integer(p_peer),
               session_id = session_id))
        b_coord <- as.numeric(beta_std[x_vars[[coordinator]]])
        if (is_coord) {
          b_nl <- c(as.numeric(beta_std[x_vars[[fusion]]]))
          for (ns in non_dcf) b_nl <- c(b_nl, as.numeric(beta_std[x_vars[[ns]]]))
        } else {
          b_nl <- numeric(0)
          for (ns in non_dcf) b_nl <- c(b_nl, as.numeric(beta_std[x_vars[[ns]]]))
          b_nl <- c(b_nl, as.numeric(beta_std[x_vars[[fusion]]]))
        }
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "k2ComputeEtaShareDS",
               beta_coord = b_coord, beta_nl = b_nl,
               intercept = if (is_coord) {
                 as.numeric(beta_std["(Intercept)"])
               } else {
                 0
               },
               is_coordinator = is_coord,
               output_key = "gee_ar1_eta",
               session_id = session_id))
        if (!is_coord && p_extras > 0L) {
          .dsAgg(datasources[dcf_conns[[i]]],
            call(name = "glmRing63ReorderXFullDS",
                 p_coord = as.integer(p_coord),
                 p_fusion = as.integer(p_fusion),
                 p_extras = as.integer(p_extras),
                 session_id = session_id))
        }
      }
    }
    invisible(TRUE)
  }

  .ring_affine <- function(a_key = NULL, b_key = NULL,
                           sign_a = 0L, sign_b = 0L,
                           output_key, public_const_fp = NULL) {
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2Ring127AffineCombineDS",
             a_key = a_key, b_key = b_key,
             sign_a = as.numeric(sign_a), sign_b = as.numeric(sign_b),
             public_const_fp = public_const_fp,
             is_party0 = identical(srv, coordinator),
             output_key = output_key, n = as.numeric(n_obs),
             session_id = session_id))
    }
    invisible(output_key)
  }
  .ring_scale <- function(in_key, scalar, output_key) {
    scalar_fp <- .fp_const(scalar)
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2Ring127LocalScaleDS",
             in_key = in_key, scalar_fp = scalar_fp,
             output_key = output_key, n = as.numeric(n_obs),
             session_id = session_id,
             is_party0 = identical(srv, coordinator)))
    }
    invisible(output_key)
  }
  .vecmul <- function(x_key, y_key, output_key) {
    tri <- .dsAgg(datasources[dealer_conn],
      call(name = "k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           n = as.numeric(n_obs), session_id = session_id,
           frac_bits = frac_bits, ring = ring))
    if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple", dcf_conns[[1L]])
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple", dcf_conns[[2L]])
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2BeaverVecmulConsumeTripleDS",
             session_id = session_id))
    }
    r1 <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      peer <- dcf_parties[[3L - i]]
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR1DS",
             peer_pk = transport_pks[[peer]],
             x_key = x_key, y_key = y_key,
             n = as.numeric(n_obs), session_id = session_id,
             frac_bits = frac_bits, ring = ring))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[i]] <- r
    }
    .sendBlob(r1[[1L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              dcf_conns[[2L]])
    .sendBlob(r1[[2L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              dcf_conns[[1L]])
    for (i in seq_along(dcf_parties)) {
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR2DS",
             is_party0 = (i == 1L),
             x_key = x_key, y_key = y_key,
             output_key = output_key,
             n = as.numeric(n_obs), session_id = session_id,
             frac_bits = frac_bits, ring = ring))
    }
    invisible(output_key)
  }
  .transform <- function(source_key, output_key, transform) {
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "dsvertGEEAR1TransformShareDS",
             source_key = source_key, output_key = output_key,
             transform = transform, session_id = session_id,
             frac_bits = frac_bits, ring = ring))
    }
    invisible(output_key)
  }
  .sum_scalar <- function(key) {
    parts <- vector("list", length(dcf_parties))
    for (i in seq_along(dcf_parties)) {
      part <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverSumShareDS",
             source_key = key, session_id = session_id,
             frac_bits = frac_bits, ring = ring_tag))
      if (is.list(part) && length(part) == 1L) part <- part[[1L]]
      parts[[i]] <- part
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = parts[[1L]]$sum_share_fp,
      share_b = parts[[2L]]$sum_share_fp,
      frac_bits = frac_bits, ring = ring_tag))
    as.numeric(agg$values[1L])
  }
  .prod_sum <- function(x_key, y_key, output_key) {
    .vecmul(x_key, y_key, output_key)
    .sum_scalar(output_key)
  }
  .cluster_sum <- function(key) {
    parts <- vector("list", length(dcf_parties))
    for (i in seq_along(dcf_parties)) {
      part <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "dsvertPerClusterSumShareDS",
             share_key = key, session_id = session_id,
             frac_bits = frac_bits, ring = ring))
      if (is.list(part) && length(part) == 1L) part <- part[[1L]]
      parts[[i]] <- part
    }
    n_clusters <- length(parts[[1L]]$per_cluster_fp)
    out <- numeric(n_clusters)
    for (ck in seq_len(n_clusters)) {
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = parts[[1L]]$per_cluster_fp[[ck]],
        share_b = parts[[2L]]$per_cluster_fp[[ck]],
        frac_bits = frac_bits, ring = ring_tag))
      out[ck] <- as.numeric(agg$values[1L])
    }
    coord_i <- match(coordinator, dcf_parties)
    list(values = out,
         cluster_sizes = as.integer(parts[[coord_i]]$cluster_sizes))
  }

  beta_start <- c(std$intercept, as.numeric(std$beta[target_features]))
  names(beta_start) <- target_order
  .restore_eta_shape(beta_start)

  coord_i <- match(coordinator, dcf_parties)
  peer_i <- setdiff(seq_along(dcf_parties), coord_i)
  ob <- .dsAgg(datasources[dcf_conns[[coord_i]]],
    call(name = "dsvertGEEAR1OrderBroadcastDS",
         data_name = data, cluster_col = id_col, order_col = order_col,
         peer_pk = transport_pks[[dcf_parties[[peer_i]]]],
         session_id = session_id))
  if (is.list(ob) && length(ob) == 1L) ob <- ob[[1L]]
  .sendBlob(ob$peer_blob, "dsvert_gee_ar1_order_blob", dcf_conns[[peer_i]])
  .dsAgg(datasources[dcf_conns[[peer_i]]],
    call(name = "dsvertGEEAR1OrderReceiveDS", session_id = session_id))

  cb <- .dsAgg(datasources[dcf_conns[[coord_i]]],
    call(name = "dsvertClusterIDsBroadcastDS",
         data_name = data, cluster_col = id_col,
         peer_pk = transport_pks[[dcf_parties[[peer_i]]]],
         session_id = session_id))
  if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
  .sendBlob(cb$peer_blob, "dsvert_cluster_ids_blob", dcf_conns[[peer_i]])
  .dsAgg(datasources[dcf_conns[[peer_i]]],
    call(name = "dsvertClusterIDsReceiveDS", session_id = session_id))

  for (i in seq_along(dcf_parties)) {
    .dsAgg(datasources[dcf_conns[[i]]],
      call(name = "dsvertGEEInterceptShareDS",
           output_key = "gee_ar1_x_col_0",
           n = as.numeric(n_obs),
           is_party0 = (i == 1L),
           session_id = session_id,
           frac_bits = frac_bits, ring = ring))
  }
  if (p_total > 0L) {
    for (j in seq_len(p_total)) {
      for (ci in dcf_conns) {
        .dsAgg(datasources[ci],
          call(name = "k2BeaverExtractColumnDS",
               source_key = "k2_x_full_fp",
               n = as.numeric(n_obs), K = as.numeric(p_total),
               col_index = as.numeric(j),
               output_key = paste0("gee_ar1_x_col_", j),
               session_id = session_id,
               frac_bits = frac_bits, ring = ring_tag))
      }
    }
  }

  x_keys <- paste0("gee_ar1_x_col_", 0:p_total)
  y_key <- "k2_y_share_fp_original"
  nonlast_x_keys <- paste0("gee_ar1_nonlast_x_col_", 0:p_total)
  interior_x_keys <- paste0("gee_ar1_interior_x_col_", 0:p_total)
  max_lag <- max(1L, as.integer(ob$max_lag %||% 1L))
  max_lag_cap <- getOption("dsvert.gee_ar1_max_lag", max_lag)
  max_lag <- min(max_lag, max(1L, as.integer(max_lag_cap)))
  lead_x_by_lag <- vector("list", max_lag)
  lead_y_by_lag <- character(max_lag)
  lead_x_by_lag[[1L]] <- paste0("gee_ar1_lead1_x_col_", 0:p_total)
  lead_y_by_lag[[1L]] <- "gee_ar1_lead1_y"
  for (j in seq_len(q)) {
    .transform(x_keys[[j]], lead_x_by_lag[[1L]][[j]], "lead")
    .transform(x_keys[[j]], nonlast_x_keys[[j]], "nonlast")
    .transform(x_keys[[j]], interior_x_keys[[j]], "interior")
  }
  .transform(y_key, lead_y_by_lag[[1L]], "lead")
  .transform(y_key, "gee_ar1_nonlast_y", "nonlast")
  .transform(y_key, "gee_ar1_interior_y", "interior")
  if (max_lag > 1L) {
    for (lag in 2:max_lag) {
      lead_x_by_lag[[lag]] <- paste0("gee_ar1_lead", lag, "_x_col_",
                                     0:p_total)
      lead_y_by_lag[[lag]] <- paste0("gee_ar1_lead", lag, "_y")
      for (j in seq_len(q)) {
        .transform(lead_x_by_lag[[lag - 1L]][[j]],
                   lead_x_by_lag[[lag]][[j]], "lead")
      }
      .transform(lead_y_by_lag[[lag - 1L]], lead_y_by_lag[[lag]], "lead")
    }
  }

  build_mask_stats <- function(prefix, x_masked_keys, y_masked_key) {
    xy <- numeric(q)
    xy[1L] <- .sum_scalar(y_masked_key)
    if (q > 1L) {
      for (j in 2:q) {
        xy[j] <- .prod_sum(x_masked_keys[[j]], y_key,
                           paste0(prefix, "_xy_", j))
      }
    }
    xx <- matrix(0, q, q)
    for (j in seq_len(q)) {
      for (k in j:q) {
        val <- if (j == 1L && k == 1L) {
          .sum_scalar(x_masked_keys[[1L]])
        } else if (j == 1L) {
          .sum_scalar(x_masked_keys[[k]])
        } else {
          .prod_sum(x_masked_keys[[j]], x_keys[[k]],
                    paste0(prefix, "_xx_", j, "_", k))
        }
        xx[j, k] <- val
        xx[k, j] <- val
      }
    }
    yy <- .prod_sum(y_masked_key, y_key, paste0(prefix, "_yy"))
    dimnames(xx) <- list(canonical_order, canonical_order)
    names(xy) <- canonical_order
    list(xx = xx, xy = xy, yy = yy)
  }

  if (verbose) {
    message("[ds.vertGEE] Gaussian AR1: protected adjacent-pair stats")
  }
  total <- build_mask_stats("gee_ar1_total", x_keys, y_key)
  nonlast <- build_mask_stats("gee_ar1_nonlast", nonlast_x_keys,
                              "gee_ar1_nonlast_y")
  interior <- build_mask_stats("gee_ar1_interior", interior_x_keys,
                               "gee_ar1_interior_y")

  build_lag_adj <- function(lag) {
    lead_x_keys <- lead_x_by_lag[[lag]]
    lead_y_key <- lead_y_by_lag[[lag]]
    prefix <- paste0("gee_ar1_adj_lag", lag)
    adj_xx <- matrix(0, q, q)
    for (j in seq_len(q)) {
      for (k in seq_len(q)) {
        adj_xx[j, k] <- if (j == 1L && k == 1L) {
          .sum_scalar(lead_x_keys[[1L]])
        } else if (j == 1L) {
          .sum_scalar(lead_x_keys[[k]])
        } else if (k == 1L) {
          .prod_sum(x_keys[[j]], lead_x_keys[[1L]],
                    paste0(prefix, "_xx_", j, "_", k))
        } else {
          .prod_sum(x_keys[[j]], lead_x_keys[[k]],
                    paste0(prefix, "_xx_", j, "_", k))
        }
      }
    }
    adj_xy_forward <- numeric(q)
    adj_xy_backward <- numeric(q)
    adj_xy_forward[1L] <- .sum_scalar(lead_y_key)
    adj_xy_backward[1L] <- .prod_sum(lead_x_keys[[1L]], y_key,
                                     paste0(prefix, "_xy_b_1"))
    if (q > 1L) {
      for (j in 2:q) {
        adj_xy_forward[j] <- .prod_sum(x_keys[[j]], lead_y_key,
                                       paste0(prefix, "_xy_f_", j))
        adj_xy_backward[j] <- .prod_sum(lead_x_keys[[j]], y_key,
                                        paste0(prefix, "_xy_b_", j))
      }
    }
    adj_yy <- .prod_sum(y_key, lead_y_key, paste0(prefix, "_yy"))
    dimnames(adj_xx) <- list(canonical_order, canonical_order)
    names(adj_xy_forward) <- names(adj_xy_backward) <- canonical_order
    list(xx = adj_xx, xy_forward = adj_xy_forward,
         xy_backward = adj_xy_backward, yy = adj_yy,
         count = as.numeric(adj_xx[1L, 1L]), lag = as.integer(lag))
  }
  lag_adjs <- lapply(seq_len(max_lag), build_lag_adj)
  adj <- lag_adjs[[1L]]

  perm <- match(target_order, canonical_order)
  if (anyNA(perm)) {
    stop("could not align Gaussian AR1 GEE order", call. = FALSE)
  }
  perm_stats <- function(x) {
    out <- list(xx = x$xx[perm, perm, drop = FALSE],
                xy = x$xy[perm], yy = x$yy)
    dimnames(out$xx) <- list(target_order, target_order)
    names(out$xy) <- target_order
    out
  }
  total_t <- perm_stats(total)
  nonlast_t <- perm_stats(nonlast)
  interior_t <- perm_stats(interior)
  adj_t <- list(
    xx = adj$xx[perm, perm, drop = FALSE],
    xy_forward = adj$xy_forward[perm],
    xy_backward = adj$xy_backward[perm],
    yy = adj$yy,
    count = adj$count,
    lag = adj$lag)
  dimnames(adj_t$xx) <- list(target_order, target_order)
  names(adj_t$xy_forward) <- names(adj_t$xy_backward) <- target_order
  perm_adj <- function(x) {
    out <- list(
      xx = x$xx[perm, perm, drop = FALSE],
      xy_forward = x$xy_forward[perm],
      xy_backward = x$xy_backward[perm],
      yy = x$yy,
      count = x$count,
      lag = x$lag)
    dimnames(out$xx) <- list(target_order, target_order)
    names(out$xy_forward) <- names(out$xy_backward) <- target_order
    out
  }
  lag_adjs_t <- lapply(lag_adjs, perm_adj)

  fit_ar1 <- .ds_gee_solve_ar1_from_stats(
    total = total_t, interior = interior_t, nonlast = nonlast_t,
    adj = adj_t, beta_start = beta_start,
    max_iter = max_iter, tol = tol, lag_stats = lag_adjs_t)
  names(fit_ar1$beta) <- target_order
  dimnames(fit_ar1$model_cov) <- list(target_order, target_order)
  dimnames(fit_ar1$bread) <- list(target_order, target_order)

  .restore_eta_shape(fit_ar1$beta)
  .ring_affine(a_key = "gee_ar1_eta", sign_a = 1L,
               output_key = "secure_mu_share")
  for (ci in dcf_conns) {
    .dsAgg(datasources[ci],
      call(name = "k2PrepareWeightedResidualShareDS",
           session_id = session_id))
  }
  .ring_affine(a_key = "k2_weight_residual_share_fp", sign_a = -1L,
               output_key = "gee_ar1_resid")
  .transform("gee_ar1_resid", "gee_ar1_resid_lead", "lead")
  .transform("gee_ar1_resid", "gee_ar1_resid_lag", "lag")
  .transform("gee_ar1_resid", "gee_ar1_resid_interior", "interior")
  rho <- as.numeric(fit_ar1$rho)
  cfac <- 1 / (1 - rho^2)
  .ring_scale("gee_ar1_resid", cfac, "gee_ar1_z_base")
  .ring_scale("gee_ar1_resid_interior", cfac * rho^2,
              "gee_ar1_z_interior")
  .ring_scale("gee_ar1_resid_lead", -cfac * rho, "gee_ar1_z_lead")
  .ring_scale("gee_ar1_resid_lag", -cfac * rho, "gee_ar1_z_lag")
  .ring_affine(a_key = "gee_ar1_z_base", b_key = "gee_ar1_z_interior",
               sign_a = 1L, sign_b = 1L, output_key = "gee_ar1_z_tmp1")
  .ring_affine(a_key = "gee_ar1_z_tmp1", b_key = "gee_ar1_z_lead",
               sign_a = 1L, sign_b = 1L, output_key = "gee_ar1_z_tmp2")
  .ring_affine(a_key = "gee_ar1_z_tmp2", b_key = "gee_ar1_z_lag",
               sign_a = 1L, sign_b = 1L, output_key = "gee_ar1_z")

  S_canonical <- NULL
  cluster_sizes <- NULL
  for (j in seq_len(q)) {
    score_key <- paste0("gee_ar1_score_x_", j)
    .vecmul(x_keys[[j]], "gee_ar1_z", score_key)
    cs <- .cluster_sum(score_key)
    if (is.null(S_canonical)) {
      S_canonical <- matrix(0, nrow = length(cs$values), ncol = q)
      cluster_sizes <- cs$cluster_sizes
    }
    S_canonical[, j] <- cs$values
    if (verbose) {
      message(sprintf("[ds.vertGEE] Gaussian AR1 score column %d/%d", j, q))
    }
  }
  S_std <- S_canonical[, perm, drop = FALSE]
  dimnames(S_std) <- list(NULL, target_order)
  robust_std <- fit_ar1$bread %*% crossprod(S_std) %*% fit_ar1$bread
  robust_std <- (robust_std + t(robust_std)) / 2
  dimnames(robust_std) <- list(target_order, target_order)

  J <- .ds_gee_standardization_jacobian(fit, target_features)
  model_cov <- J %*% fit_ar1$model_cov %*% t(J)
  robust_cov <- J %*% robust_std %*% t(J)
  model_cov <- (model_cov + t(model_cov)) / 2
  robust_cov <- (robust_cov + t(robust_cov)) / 2
  dimnames(model_cov) <- list(target_order, target_order)
  dimnames(robust_cov) <- list(target_order, target_order)

  coefficients <- .ds_gee_unstandardized_parameters(
    fit, fit_ar1$beta, target_features)
  model_se <- sqrt(pmax(diag(model_cov), 0))
  robust_se <- sqrt(pmax(diag(robust_cov), 0))
  names(model_se) <- names(robust_se) <- target_order

  list(coefficients = coefficients,
       model_se = model_se,
       robust_se = robust_se,
       model_covariance = model_cov,
       robust_covariance = robust_cov,
       method = "share_domain_ar1_gaussian",
       working_correlation = list(
         corstr = "ar1",
         alpha = as.numeric(fit_ar1$rho),
         phi = as.numeric(fit_ar1$phi),
         iterations = as.integer(fit_ar1$iterations),
         converged = isTRUE(fit_ar1$converged),
         estimator = "gaussian_ar1_adjacent_pair_stats",
         disclosure = paste0(
           "guarded adjacent-pair sufficient statistics and cluster scores ",
           "only; no row-level order, residuals, fitted values, or visit ",
           "labels returned")),
       cluster_sizes = cluster_sizes)
}

#' @keywords internal
.ds_gee_secure_gaussian_exchangeable <- function(
    fit, datasources, server_names, lambda = 0, data = NULL, id_col = NULL,
    max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  if (!identical(fit$family, "gaussian")) return(NULL)
  if (is.null(id_col) || length(id_col) != 1L ||
      !is.character(id_col) || !nzchar(id_col)) {
    stop("Gaussian exchangeable GEE requires id_col", call. = FALSE)
  }
  required <- c("session_id", "transport_pks", "server_list", "x_vars",
                "y_server", "x_sds", "x_means")
  missing_req <- required[vapply(required, function(nm) is.null(fit[[nm]]),
                                 logical(1L))]
  if (length(missing_req) > 0L) {
    stop("GLM session metadata missing: ", paste(missing_req, collapse = ", "),
         call. = FALSE)
  }

  session_id <- fit$session_id
  server_list <- fit$server_list
  x_vars <- fit$x_vars
  coordinator <- fit$y_server
  n_obs <- as.integer(fit$n_obs)
  ring <- as.integer(fit$ring %||% 63L)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  frac_bits <- if (ring == 127L) 50L else 20L
  ring_tag <- if (ring == 127L) "ring127" else "ring63"
  transport_pks <- fit$transport_pks

  .to_b64url <- function(x) gsub("+", "-", gsub("/", "_",
    gsub("=+$", "", x, perl = TRUE), fixed = TRUE), fixed = TRUE)
  .dsAgg <- function(conns, expr, ...) {
    tryCatch(
      DSI::datashield.aggregate(conns = conns, expr = expr, ...),
      error = function(e) {
        fn_name <- if (is.call(expr)) as.character(expr[[1]]) else "?"
        srv_name <- tryCatch(names(conns)[1], error = function(x) "?")
        stop(sprintf("%s on %s failed: %s", fn_name, srv_name,
                     conditionMessage(e)), call. = FALSE)
      })
  }
  .sendBlob <- function(blob, key, conn_idx) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        .dsAgg(datasources[conn_idx],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               session_id = session_id))
      } else {
        .dsAgg(datasources[conn_idx],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               chunk_index = chunk_idx, n_chunks = n_chunks,
               session_id = session_id))
      }
    })
  }

  target_features <- unlist(x_vars[server_list], use.names = FALSE)
  target_order <- c("(Intercept)", target_features)
  std <- .ds_gee_standardized_parameters(fit, target_features)

  if (fit$eta_privacy == "k2_beaver") {
    nl <- setdiff(server_list, coordinator)
    if (length(nl) != 1L) {
      stop("K=2 Gaussian exchangeable GEE requires exactly one non-outcome ",
           "server", call. = FALSE)
    }
    nl <- nl[[1L]]
    dcf_parties <- c(coordinator, nl)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    dealer_conn <- dcf_conns[[2L]]
    canonical_features <- c(x_vars[[coordinator]], x_vars[[nl]])
    b_coord <- as.numeric(std$beta[x_vars[[coordinator]]])
    b_nl <- as.numeric(std$beta[x_vars[[nl]]])
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "dsvertGEERestoreFeatureShapeDS",
             p_own = as.integer(length(x_vars[[srv]])),
             p_peer = as.integer(length(x_vars[[setdiff(dcf_parties, srv)]])),
             session_id = session_id))
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2ComputeEtaShareDS",
             beta_coord = b_coord, beta_nl = b_nl,
             intercept = if (srv == coordinator) std$intercept else 0,
             is_coordinator = (srv == coordinator),
             session_id = session_id))
    }
  } else if (fit$eta_privacy == "secure_agg") {
    fusion <- .k3_select_fusion_server(server_list, coordinator, x_vars)
    dcf_parties <- c(fusion, coordinator)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    non_dcf <- setdiff(server_list, dcf_parties)
    dealer <- if (length(non_dcf) > 0L) non_dcf[[1L]] else fusion
    dealer_conn <- which(server_names == dealer)
    p_coord <- length(x_vars[[coordinator]])
    p_fusion <- length(x_vars[[fusion]])
    p_extras <- sum(vapply(non_dcf, function(s) length(x_vars[[s]]),
                           integer(1L)))
    canonical_features <- c(x_vars[[coordinator]], x_vars[[fusion]])
    for (srv in non_dcf) {
      canonical_features <- c(canonical_features, x_vars[[srv]])
    }
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      is_coord <- srv == coordinator
      p_peer <- if (is_coord) p_fusion + p_extras else p_coord + p_extras
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "dsvertGEERestoreFeatureShapeDS",
             p_own = as.integer(length(x_vars[[srv]])),
             p_peer = as.integer(p_peer),
             session_id = session_id))
      b_coord <- as.numeric(std$beta[x_vars[[coordinator]]])
      if (is_coord) {
        b_nl <- c(as.numeric(std$beta[x_vars[[fusion]]]))
        for (ns in non_dcf) b_nl <- c(b_nl, as.numeric(std$beta[x_vars[[ns]]]))
      } else {
        b_nl <- numeric(0)
        for (ns in non_dcf) b_nl <- c(b_nl, as.numeric(std$beta[x_vars[[ns]]]))
        b_nl <- c(b_nl, as.numeric(std$beta[x_vars[[fusion]]]))
      }
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2ComputeEtaShareDS",
             beta_coord = b_coord, beta_nl = b_nl,
             intercept = if (is_coord) std$intercept else 0,
             is_coordinator = is_coord,
             session_id = session_id))
      if (!is_coord && p_extras > 0L) {
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "glmRing63ReorderXFullDS",
               p_coord = as.integer(p_coord),
               p_fusion = as.integer(p_fusion),
               p_extras = as.integer(p_extras),
               session_id = session_id))
      }
    }
  } else {
    stop("unsupported eta_privacy for Gaussian exchangeable GEE: ",
         fit$eta_privacy, call. = FALSE)
  }

  cluster_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                             data, id_col)
  if (is.null(cluster_srv)) {
    stop("id_col '", id_col, "' not found on any server", call. = FALSE)
  }
  if (!identical(cluster_srv, coordinator)) {
    stop("Gaussian exchangeable GEE requires id_col to live on the ",
         "outcome server; found on '", cluster_srv, "'", call. = FALSE)
  }
  coord_i <- match(coordinator, dcf_parties)
  peer_i <- setdiff(seq_along(dcf_parties), coord_i)
  cb <- .dsAgg(datasources[dcf_conns[[coord_i]]],
    call(name = "dsvertClusterIDsBroadcastDS",
         data_name = data, cluster_col = id_col,
         peer_pk = transport_pks[[dcf_parties[[peer_i]]]],
         session_id = session_id))
  if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
  .sendBlob(cb$peer_blob, "dsvert_cluster_ids_blob", dcf_conns[[peer_i]])
  .dsAgg(datasources[dcf_conns[[peer_i]]],
    call(name = "dsvertClusterIDsReceiveDS", session_id = session_id))

  p_total <- length(canonical_features)
  canonical_order <- c("(Intercept)", canonical_features)
  q <- p_total + 1L

  .vecmul <- function(x_key, y_key, output_key) {
    tri <- .dsAgg(datasources[dealer_conn],
      call(name = "k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           n = as.numeric(n_obs), session_id = session_id,
           frac_bits = frac_bits, ring = ring))
    if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple", dcf_conns[[1L]])
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple", dcf_conns[[2L]])
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2BeaverVecmulConsumeTripleDS",
             session_id = session_id))
    }
    r1 <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      peer <- dcf_parties[[3L - i]]
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR1DS",
             peer_pk = transport_pks[[peer]],
             x_key = x_key, y_key = y_key,
             n = as.numeric(n_obs), session_id = session_id,
             frac_bits = frac_bits, ring = ring))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[i]] <- r
    }
    .sendBlob(r1[[1L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              dcf_conns[[2L]])
    .sendBlob(r1[[2L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              dcf_conns[[1L]])
    for (i in seq_along(dcf_parties)) {
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR2DS",
             is_party0 = (i == 1L),
             x_key = x_key, y_key = y_key,
             output_key = output_key,
             n = as.numeric(n_obs), session_id = session_id,
             frac_bits = frac_bits, ring = ring))
    }
    invisible(output_key)
  }

  .cluster_sum <- function(key) {
    parts <- vector("list", length(dcf_parties))
    for (i in seq_along(dcf_parties)) {
      part <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "dsvertPerClusterSumShareDS",
             share_key = key, session_id = session_id,
             frac_bits = frac_bits, ring = ring))
      if (is.list(part) && length(part) == 1L) part <- part[[1L]]
      parts[[i]] <- part
    }
    n_clusters <- length(parts[[1L]]$per_cluster_fp)
    out <- numeric(n_clusters)
    for (ck in seq_len(n_clusters)) {
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = parts[[1L]]$per_cluster_fp[[ck]],
        share_b = parts[[2L]]$per_cluster_fp[[ck]],
        frac_bits = frac_bits, ring = ring_tag))
      out[ck] <- as.numeric(agg$values[1L])
    }
    list(values = out,
         cluster_sizes = as.integer(parts[[coord_i]]$cluster_sizes))
  }

  for (i in seq_along(dcf_parties)) {
    .dsAgg(datasources[dcf_conns[[i]]],
      call(name = "dsvertGEEInterceptShareDS",
           output_key = "gee_ex_x_col_0",
           n = as.numeric(n_obs),
           is_party0 = (i == 1L),
           session_id = session_id,
           frac_bits = frac_bits, ring = ring))
  }
  for (j in seq_len(p_total)) {
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2BeaverExtractColumnDS",
             source_key = "k2_x_full_fp",
             n = as.numeric(n_obs), K = as.numeric(p_total),
             col_index = as.numeric(j),
             output_key = paste0("gee_ex_x_col_", j),
             session_id = session_id,
             frac_bits = frac_bits, ring = ring_tag))
    }
  }

  if (verbose) {
    message("[ds.vertGEE] Gaussian exchangeable: protected cluster stats")
  }
  x_keys <- paste0("gee_ex_x_col_", 0:p_total)
  sx <- NULL
  cluster_sizes <- NULL
  for (j in seq_len(q)) {
    cs <- .cluster_sum(x_keys[[j]])
    if (is.null(sx)) {
      sx <- matrix(0, nrow = length(cs$values), ncol = q)
      cluster_sizes <- cs$cluster_sizes
    }
    sx[, j] <- cs$values
  }
  sy <- .cluster_sum("k2_y_share_fp_original")$values

  sxy <- matrix(0, nrow = nrow(sx), ncol = q)
  sxy[, 1L] <- sy
  if (q > 1L) {
    for (j in 2:q) {
      key <- paste0("gee_ex_xy_", j - 1L)
      .vecmul(x_keys[[j]], "k2_y_share_fp_original", key)
      sxy[, j] <- .cluster_sum(key)$values
    }
  }
  .vecmul("k2_y_share_fp_original", "k2_y_share_fp_original", "gee_ex_y2")
  syy <- .cluster_sum("gee_ex_y2")$values

  xx <- array(0, dim = c(nrow(sx), q, q))
  for (j in seq_len(q)) {
    for (k in j:q) {
      if (j == 1L && k == 1L) {
        vals <- as.numeric(cluster_sizes)
      } else if (j == 1L) {
        vals <- sx[, k]
      } else if (k == 1L) {
        vals <- sx[, j]
      } else {
        key <- paste0("gee_ex_xx_", j - 1L, "_", k - 1L)
        .vecmul(x_keys[[j]], x_keys[[k]], key)
        vals <- .cluster_sum(key)$values
      }
      xx[, j, k] <- vals
      xx[, k, j] <- vals
    }
    if (verbose) {
      message(sprintf("[ds.vertGEE] Gaussian exchangeable stats column %d/%d",
                      j, q))
    }
  }

  perm <- match(target_order, canonical_order)
  if (anyNA(perm)) {
    stop("could not align Gaussian exchangeable GEE order", call. = FALSE)
  }
  sx <- sx[, perm, drop = FALSE]
  sxy <- sxy[, perm, drop = FALSE]
  xx <- xx[, perm, perm, drop = FALSE]
  dimnames(sx) <- list(NULL, target_order)
  dimnames(sxy) <- list(NULL, target_order)
  dimnames(xx) <- list(NULL, target_order, target_order)

  beta_start <- c(std$intercept, as.numeric(std$beta[target_features]))
  names(beta_start) <- target_order
  fit_ex <- .ds_gee_solve_exchangeable_from_cluster_stats(
    sx = sx, sy = sy, sxy = sxy, xx = xx, syy = syy,
    beta_start = beta_start, max_iter = max_iter, tol = tol)
  names(fit_ex$beta) <- target_order
  dimnames(fit_ex$model_cov) <- list(target_order, target_order)
  dimnames(fit_ex$robust_cov) <- list(target_order, target_order)

  J <- .ds_gee_standardization_jacobian(fit, target_features)
  model_cov <- J %*% fit_ex$model_cov %*% t(J)
  robust_cov <- J %*% fit_ex$robust_cov %*% t(J)
  model_cov <- (model_cov + t(model_cov)) / 2
  robust_cov <- (robust_cov + t(robust_cov)) / 2
  dimnames(model_cov) <- list(target_order, target_order)
  dimnames(robust_cov) <- list(target_order, target_order)

  coefficients <- .ds_gee_unstandardized_parameters(
    fit, fit_ex$beta, target_features)
  model_se <- sqrt(pmax(diag(model_cov), 0))
  robust_se <- sqrt(pmax(diag(robust_cov), 0))
  names(model_se) <- names(robust_se) <- target_order

  list(coefficients = coefficients,
       model_se = model_se,
       robust_se = robust_se,
       model_covariance = model_cov,
       robust_covariance = robust_cov,
       method = "share_domain_exchangeable_gaussian",
       working_correlation = list(
         corstr = "exchangeable",
         alpha = as.numeric(fit_ex$alpha),
         phi = as.numeric(fit_ex$phi),
         iterations = as.integer(fit_ex$iterations),
         converged = isTRUE(fit_ex$converged),
         estimator = "gaussian_exchangeable_cluster_stats",
         disclosure = paste0(
           "guarded cluster-level sufficient statistics only; no row-level ",
           "residuals, fitted values, scores, or cluster labels returned")),
       cluster_sizes = cluster_sizes)
}

#' @keywords internal
.ds_gee_secure_poisson_exchangeable <- function(
    fit, datasources, server_names, lambda = 0, data = NULL, id_col = NULL,
    order_col = NULL, corstr = c("exchangeable", "ar1"),
    max_iter = 100L, tol = 1e-6, verbose = FALSE,
    binomial_sigmoid_intervals = NULL) {
  corstr <- match.arg(corstr)
  if (!(fit$family %in% c("binomial", "poisson"))) return(NULL)
  gee_family <- fit$family
  cor_label <- if (identical(corstr, "ar1")) "AR1" else "exchangeable"
  family_title <- paste0(toupper(substr(gee_family, 1L, 1L)),
                         substr(gee_family, 2L, nchar(gee_family)))
  if (is.null(id_col) || length(id_col) != 1L ||
      !is.character(id_col) || !nzchar(id_col)) {
    stop(family_title, " ", cor_label, " GEE requires id_col",
         call. = FALSE)
  }
  if (identical(corstr, "ar1") &&
      (is.null(order_col) || length(order_col) != 1L ||
       !is.character(order_col) || !nzchar(order_col))) {
    stop(family_title, " AR1 GEE requires order_col", call. = FALSE)
  }
  required <- c("session_id", "transport_pks", "server_list", "x_vars",
                "y_server", "x_sds", "x_means")
  missing_req <- required[vapply(required, function(nm) is.null(fit[[nm]]),
                                 logical(1L))]
  if (length(missing_req) > 0L) {
    stop("GLM session metadata missing: ", paste(missing_req, collapse = ", "),
         call. = FALSE)
  }

  session_id <- fit$session_id
  server_list <- fit$server_list
  x_vars <- fit$x_vars
  coordinator <- fit$y_server
  n_obs <- as.integer(fit$n_obs)
  ring <- as.integer(fit$ring %||% 63L)
  if (ring != 127L) {
    stop(family_title, " ", cor_label, " GEE requires Ring127",
         call. = FALSE)
  }
  frac_bits <- 50L
  ring_tag <- "ring127"
  transport_pks <- fit$transport_pks
  sigmoid_intervals <- if (identical(gee_family, "binomial")) {
    suppressWarnings(as.integer(
      (binomial_sigmoid_intervals %||%
         getOption("dsvert.gee_binomial_sigmoid_intervals",
                   getOption("dsvert.glm_num_intervals_binomial",
                             getOption("dsvert.glm_num_intervals", 100L))))[[1L]]
    ))
  } else {
    NA_integer_
  }
  if (identical(gee_family, "binomial") &&
      (!is.finite(sigmoid_intervals) || sigmoid_intervals < 10L)) {
    sigmoid_intervals <- 100L
  }

  .to_b64url <- function(x) gsub("+", "-", gsub("/", "_",
    gsub("=+$", "", x, perl = TRUE), fixed = TRUE), fixed = TRUE)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-", "+", gsub("_", "/", x, fixed = TRUE), fixed = TRUE)
    pad <- nchar(x) %% 4
    if (pad == 2) x <- paste0(x, "==")
    if (pad == 3) x <- paste0(x, "=")
    x
  }
  .fp_const <- function(x) {
    .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
      values = array(as.numeric(x), dim = 1L),
      frac_bits = frac_bits, ring = ring_tag))$fp_data)
  }
  .dsAgg <- function(conns, expr, ...) {
    tryCatch(
      DSI::datashield.aggregate(conns = conns, expr = expr, ...),
      error = function(e) {
        fn_name <- if (is.call(expr)) as.character(expr[[1]]) else "?"
        srv_name <- tryCatch(names(conns)[1], error = function(x) "?")
        stop(sprintf("%s on %s failed: %s", fn_name, srv_name,
                     conditionMessage(e)), call. = FALSE)
      })
  }
  .sendBlob <- function(blob, key, conn_idx) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        .dsAgg(datasources[conn_idx],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               session_id = session_id))
      } else {
        .dsAgg(datasources[conn_idx],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               chunk_index = chunk_idx, n_chunks = n_chunks,
               session_id = session_id))
      }
    })
  }

  if (fit$eta_privacy == "k2_beaver") {
    nl <- setdiff(server_list, coordinator)
    if (length(nl) != 1L) {
      stop("K=2 ", gee_family, " exchangeable GEE requires exactly one non-outcome ",
           "server", call. = FALSE)
    }
    nl <- nl[[1L]]
    dcf_parties <- c(coordinator, nl)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    dealer_conn <- dcf_conns[[2L]]
    canonical_features <- c(x_vars[[coordinator]], x_vars[[nl]])
  } else if (fit$eta_privacy == "secure_agg") {
    fusion <- .k3_select_fusion_server(server_list, coordinator, x_vars)
    dcf_parties <- c(fusion, coordinator)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    non_dcf <- setdiff(server_list, dcf_parties)
    dealer <- if (length(non_dcf) > 0L) non_dcf[[1L]] else fusion
    dealer_conn <- which(server_names == dealer)
    canonical_features <- c(x_vars[[coordinator]], x_vars[[fusion]])
    for (srv in non_dcf) {
      canonical_features <- c(canonical_features, x_vars[[srv]])
    }
  } else {
    stop("unsupported eta_privacy for ", gee_family, " ", cor_label,
         " GEE: ",
         fit$eta_privacy, call. = FALSE)
  }

  cluster_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                             data, id_col)
  if (is.null(cluster_srv)) {
    stop("id_col '", id_col, "' not found on any server", call. = FALSE)
  }
  if (!identical(cluster_srv, coordinator)) {
    stop(family_title, " ", cor_label,
         " GEE requires id_col to live on the ",
         "outcome server; found on '", cluster_srv, "'", call. = FALSE)
  }
  if (identical(corstr, "ar1")) {
    order_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                             data, order_col)
    if (is.null(order_srv)) {
      stop("order_col '", order_col, "' not found on any server",
           call. = FALSE)
    }
    if (!identical(order_srv, coordinator)) {
      stop(family_title, " AR1 GEE requires order_col to live on the ",
           "outcome server; found on '", order_srv, "'", call. = FALSE)
    }
  }

  target_features <- unlist(x_vars[server_list], use.names = FALSE)
  target_order <- c("(Intercept)", target_features)
  canonical_order <- c("(Intercept)", canonical_features)
  p_total <- length(canonical_features)
  q <- p_total + 1L
  std <- .ds_gee_standardized_parameters(fit, target_features)

  .restore_eta_shape <- function(beta_std) {
    if (fit$eta_privacy == "k2_beaver") {
      b_coord <- as.numeric(beta_std[x_vars[[coordinator]]])
      b_nl <- as.numeric(beta_std[x_vars[[dcf_parties[[2L]]]]])
      for (i in seq_along(dcf_parties)) {
        srv <- dcf_parties[[i]]
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "dsvertGEERestoreFeatureShapeDS",
               p_own = as.integer(length(x_vars[[srv]])),
               p_peer = as.integer(length(x_vars[[setdiff(dcf_parties, srv)]])),
               session_id = session_id))
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "k2ComputeEtaShareDS",
               beta_coord = b_coord, beta_nl = b_nl,
               intercept = if (srv == coordinator) {
                 as.numeric(beta_std["(Intercept)"])
               } else {
                 0
               },
               is_coordinator = (srv == coordinator),
               output_key = "gee_poisson_eta",
               session_id = session_id))
      }
    } else {
      fusion <- dcf_parties[[1L]]
      non_dcf <- setdiff(server_list, dcf_parties)
      p_coord <- length(x_vars[[coordinator]])
      p_fusion <- length(x_vars[[fusion]])
      p_extras <- sum(vapply(non_dcf, function(s) length(x_vars[[s]]),
                             integer(1L)))
      for (i in seq_along(dcf_parties)) {
        srv <- dcf_parties[[i]]
        is_coord <- srv == coordinator
        p_peer <- if (is_coord) p_fusion + p_extras else p_coord + p_extras
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "dsvertGEERestoreFeatureShapeDS",
               p_own = as.integer(length(x_vars[[srv]])),
               p_peer = as.integer(p_peer),
               session_id = session_id))
        b_coord <- as.numeric(beta_std[x_vars[[coordinator]]])
        if (is_coord) {
          b_nl <- c(as.numeric(beta_std[x_vars[[fusion]]]))
          for (ns in non_dcf) b_nl <- c(b_nl, as.numeric(beta_std[x_vars[[ns]]]))
        } else {
          b_nl <- numeric(0)
          for (ns in non_dcf) b_nl <- c(b_nl, as.numeric(beta_std[x_vars[[ns]]]))
          b_nl <- c(b_nl, as.numeric(beta_std[x_vars[[fusion]]]))
        }
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "k2ComputeEtaShareDS",
               beta_coord = b_coord, beta_nl = b_nl,
               intercept = if (is_coord) {
                 as.numeric(beta_std["(Intercept)"])
               } else {
                 0
               },
               is_coordinator = is_coord,
               output_key = "gee_poisson_eta",
               session_id = session_id))
        if (!is_coord && p_extras > 0L) {
          .dsAgg(datasources[dcf_conns[[i]]],
            call(name = "glmRing63ReorderXFullDS",
                 p_coord = as.integer(p_coord),
                 p_fusion = as.integer(p_fusion),
                 p_extras = as.integer(p_extras),
                 session_id = session_id))
        }
      }
    }
    invisible(TRUE)
  }

  .copy_key <- function(source_key, target_key) {
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2Ring127AffineCombineDS",
             a_key = source_key, b_key = NULL,
             sign_a = 1L, sign_b = 0L,
             public_const_fp = NULL,
             is_party0 = identical(srv, coordinator),
             output_key = target_key, n = as.numeric(n_obs),
             session_id = session_id))
    }
    invisible(target_key)
  }
  .ring_affine <- function(a_key = NULL, b_key = NULL,
                           sign_a = 0L, sign_b = 0L,
                           output_key, public_const_fp = NULL) {
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2Ring127AffineCombineDS",
             a_key = a_key, b_key = b_key,
             sign_a = as.numeric(sign_a), sign_b = as.numeric(sign_b),
             public_const_fp = public_const_fp,
             is_party0 = identical(srv, coordinator),
             output_key = output_key, n = as.numeric(n_obs),
             session_id = session_id))
    }
    invisible(output_key)
  }
  .ring_scale <- function(in_key, scalar, output_key) {
    scalar_fp <- .fp_const(scalar)
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2Ring127LocalScaleDS",
             in_key = in_key, scalar_fp = scalar_fp,
             output_key = output_key, n = as.numeric(n_obs),
             session_id = session_id,
             is_party0 = identical(srv, coordinator)))
    }
    invisible(output_key)
  }
  .vecmul <- function(x_key, y_key, output_key) {
    tri <- .dsAgg(datasources[dealer_conn],
      call(name = "k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           n = as.numeric(n_obs), session_id = session_id,
           frac_bits = frac_bits, ring = ring))
    if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple", dcf_conns[[1L]])
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple", dcf_conns[[2L]])
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2BeaverVecmulConsumeTripleDS",
             session_id = session_id))
    }
    r1 <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      peer <- dcf_parties[[3L - i]]
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR1DS",
             peer_pk = transport_pks[[peer]],
             x_key = x_key, y_key = y_key,
             n = as.numeric(n_obs), session_id = session_id,
             frac_bits = frac_bits, ring = ring))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[i]] <- r
    }
    .sendBlob(r1[[1L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              dcf_conns[[2L]])
    .sendBlob(r1[[2L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              dcf_conns[[1L]])
    for (i in seq_along(dcf_parties)) {
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR2DS",
             is_party0 = (i == 1L),
             x_key = x_key, y_key = y_key,
             output_key = output_key,
             n = as.numeric(n_obs), session_id = session_id,
             frac_bits = frac_bits, ring = ring))
    }
    invisible(output_key)
  }
  .transform <- function(source_key, output_key, transform) {
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "dsvertGEEAR1TransformShareDS",
             source_key = source_key, output_key = output_key,
             transform = transform, session_id = session_id,
             frac_bits = frac_bits, ring = ring))
    }
    invisible(output_key)
  }
  .sum_scalar <- function(key) {
    parts <- vector("list", length(dcf_parties))
    for (i in seq_along(dcf_parties)) {
      part <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverSumShareDS",
             source_key = key, session_id = session_id,
             frac_bits = frac_bits, ring = ring_tag))
      if (is.list(part) && length(part) == 1L) part <- part[[1L]]
      parts[[i]] <- part
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = parts[[1L]]$sum_share_fp,
      share_b = parts[[2L]]$sum_share_fp,
      frac_bits = frac_bits, ring = ring_tag))
    as.numeric(agg$values[1L])
  }
  .prod_sum <- function(x_key, y_key, output_key) {
    .vecmul(x_key, y_key, output_key)
    .sum_scalar(output_key)
  }
  .cluster_sum <- function(key) {
    parts <- vector("list", length(dcf_parties))
    for (i in seq_along(dcf_parties)) {
      part <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "dsvertPerClusterSumShareDS",
             share_key = key, session_id = session_id,
             frac_bits = frac_bits, ring = ring))
      if (is.list(part) && length(part) == 1L) part <- part[[1L]]
      parts[[i]] <- part
    }
    n_clusters <- length(parts[[1L]]$per_cluster_fp)
    out <- numeric(n_clusters)
    for (ck in seq_len(n_clusters)) {
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = parts[[1L]]$per_cluster_fp[[ck]],
        share_b = parts[[2L]]$per_cluster_fp[[ck]],
        frac_bits = frac_bits, ring = ring_tag))
      out[ck] <- as.numeric(agg$values[1L])
    }
    coord_i <- match(coordinator, dcf_parties)
    list(values = out,
         cluster_sizes = as.integer(parts[[coord_i]]$cluster_sizes))
  }
  dcf_cache_key <- NULL
  dcf_blob_cache <- new.env(parent = emptyenv())
  .ensure_dcf_keys <- function(spline_family, num_intervals) {
    key <- paste(spline_family, as.integer(num_intervals), sep = ":")
    if (identical(dcf_cache_key, key)) return(invisible(TRUE))
    dcf <- dcf_blob_cache[[key]]
    if (is.null(dcf)) {
      dcf <- .dsAgg(datasources[dealer_conn],
        call(name = "glmRing63GenDcfKeysDS",
             dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
             dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
             family = spline_family, n = as.integer(n_obs),
             frac_bits = frac_bits,
             num_intervals = as.integer(num_intervals),
             ring = ring, session_id = session_id))
      if (is.list(dcf) && length(dcf) == 1L) dcf <- dcf[[1L]]
      dcf <- list(blob0 = dcf$dcf_blob_0, blob1 = dcf$dcf_blob_1)
      dcf_blob_cache[[key]] <- dcf
    }
    .sendBlob(dcf$blob0, "k2_dcf_keys_persistent", dcf_conns[[1L]])
    .sendBlob(dcf$blob1, "k2_dcf_keys_persistent", dcf_conns[[2L]])
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2StoreDcfKeysPersistentDS", session_id = session_id))
    }
    dcf_cache_key <<- key
    invisible(TRUE)
  }
  .wide_spline <- function(input_key, output_key, spline_family,
                           num_intervals = 100L) {
    .ensure_dcf_keys(spline_family, num_intervals)
    .copy_key(input_key, "k2_eta_share_fp")
    spline_t <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenSplineTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           n = as.integer(n_obs), frac_bits = frac_bits,
           ring = ring, session_id = session_id))
    if (is.list(spline_t) && length(spline_t) == 1L) spline_t <- spline_t[[1L]]
    .sendBlob(spline_t$spline_blob_0, "k2_spline_triples", dcf_conns[[1L]])
    .sendBlob(spline_t$spline_blob_1, "k2_spline_triples", dcf_conns[[2L]])
    for (ph in 1:4) {
      pr <- vector("list", 2L)
      for (i in seq_along(dcf_parties)) {
        r <- .dsAgg(datasources[dcf_conns[[i]]],
          call(paste0("k2WideSplinePhase", ph, "DS"),
               party_id = as.integer(i - 1L),
               family = spline_family,
               num_intervals = as.integer(num_intervals),
               frac_bits = frac_bits,
               ring = ring,
               session_id = session_id))
        if (is.list(r) && length(r) == 1L) r <- r[[1L]]
        pr[[i]] <- r
      }
      if (ph == 1L) {
        .sendBlob(pr[[1L]]$dcf_masked, "k2_peer_dcf_masked",
                  dcf_conns[[2L]])
        .sendBlob(pr[[2L]]$dcf_masked, "k2_peer_dcf_masked",
                  dcf_conns[[1L]])
      } else if (ph == 2L) {
        for (i in 1:2) {
          peer_i <- 3L - i
          pk_b64 <- .b64url_to_b64(transport_pks[[dcf_parties[[peer_i]]]])
          payload <- jsonlite::toJSON(list(
            and_xma = pr[[i]]$and_xma, and_ymb = pr[[i]]$and_ymb,
            had1_xma = pr[[i]]$had1_xma, had1_ymb = pr[[i]]$had1_ymb),
            auto_unbox = TRUE)
          sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
            data = jsonlite::base64_enc(charToRaw(payload)),
            recipient_pk = pk_b64))
          .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1",
                    dcf_conns[[peer_i]])
        }
      } else if (ph == 3L) {
        for (i in 1:2) {
          peer_i <- 3L - i
          pk_b64 <- .b64url_to_b64(transport_pks[[dcf_parties[[peer_i]]]])
          payload <- jsonlite::toJSON(list(
            had2_xma = pr[[i]]$had2_xma,
            had2_ymb = pr[[i]]$had2_ymb),
            auto_unbox = TRUE)
          sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
            data = jsonlite::base64_enc(charToRaw(payload)),
            recipient_pk = pk_b64))
          .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1",
                    dcf_conns[[peer_i]])
        }
      }
    }
    .copy_key("secure_mu_share", output_key)
    invisible(output_key)
  }

  coord_i <- match(coordinator, dcf_parties)
  peer_i <- setdiff(seq_along(dcf_parties), coord_i)
  order_meta <- NULL
  if (identical(corstr, "ar1")) {
    order_meta <- .dsAgg(datasources[dcf_conns[[coord_i]]],
      call(name = "dsvertGEEAR1OrderBroadcastDS",
           data_name = data, cluster_col = id_col, order_col = order_col,
           peer_pk = transport_pks[[dcf_parties[[peer_i]]]],
           session_id = session_id))
    if (is.list(order_meta) && length(order_meta) == 1L) {
      order_meta <- order_meta[[1L]]
    }
    .sendBlob(order_meta$peer_blob, "dsvert_gee_ar1_order_blob",
              dcf_conns[[peer_i]])
    .dsAgg(datasources[dcf_conns[[peer_i]]],
      call(name = "dsvertGEEAR1OrderReceiveDS", session_id = session_id))
  }
  cb <- .dsAgg(datasources[dcf_conns[[coord_i]]],
    call(name = "dsvertClusterIDsBroadcastDS",
         data_name = data, cluster_col = id_col,
         peer_pk = transport_pks[[dcf_parties[[peer_i]]]],
         session_id = session_id))
  if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
  .sendBlob(cb$peer_blob, "dsvert_cluster_ids_blob", dcf_conns[[peer_i]])
  .dsAgg(datasources[dcf_conns[[peer_i]]],
    call(name = "dsvertClusterIDsReceiveDS", session_id = session_id))

  beta <- c(std$intercept, as.numeric(std$beta[target_features]))
  names(beta) <- target_order
  perm <- match(target_order, canonical_order)
  if (anyNA(perm)) {
    stop("could not align ", gee_family, " ", cor_label, " GEE order",
         call. = FALSE)
  }
  beta <- beta[canonical_order]
  beta[is.na(beta)] <- 0
  .restore_eta_shape(beta)

  for (i in seq_along(dcf_parties)) {
    .dsAgg(datasources[dcf_conns[[i]]],
      call(name = "dsvertGEEInterceptShareDS",
           output_key = "gee_px_col_0",
           n = as.numeric(n_obs),
           is_party0 = (i == 1L),
           session_id = session_id,
           frac_bits = frac_bits, ring = ring))
  }
  for (j in seq_len(p_total)) {
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2BeaverExtractColumnDS",
             source_key = "k2_x_full_fp",
             n = as.numeric(n_obs), K = as.numeric(p_total),
             col_index = as.numeric(j),
             output_key = paste0("gee_px_col_", j),
             session_id = session_id,
             frac_bits = frac_bits, ring = ring_tag))
    }
  }
  x_keys <- paste0("gee_px_col_", 0:p_total)

  xx_base_keys <- matrix("", nrow = q, ncol = q)
  for (j in seq_len(q)) {
    for (k in j:q) {
      if (j == 1L && k == 1L) {
        xx_base_keys[j, k] <- ""
      } else if (j == 1L) {
        xx_base_keys[j, k] <- x_keys[[k]]
      } else {
        key <- paste0("gee_px_xxbase_", j - 1L, "_", k - 1L)
        .vecmul(x_keys[[j]], x_keys[[k]], key)
        xx_base_keys[j, k] <- key
      }
      xx_base_keys[k, j] <- xx_base_keys[j, k]
    }
  }

  safe_solve <- function(A, b = NULL) {
    A <- (A + t(A)) / 2
    out <- tryCatch(if (is.null(b)) solve(A) else solve(A, b),
                    error = function(e) NULL)
    if (!is.null(out) && all(is.finite(out))) return(out)
    ridge <- 1e-8 * max(1, mean(abs(diag(A))))
    tryCatch(if (is.null(b)) solve(A + diag(ridge, nrow(A)))
             else solve(A + diag(ridge, nrow(A)), b),
             error = function(e) {
               stop(family_title, " ", cor_label, " GEE system is singular",
                    call. = FALSE)
             })
  }

  if (identical(corstr, "ar1")) {
    rho <- 0
    converged <- FALSE
    cluster_sizes <- NULL
    last_stats <- NULL
    lower_rho <- -0.95
    upper_rho <- 0.95

    prepare_family_shares <- function(beta_current) {
      beta_named <- beta_current
      names(beta_named) <- canonical_order
      .restore_eta_shape(beta_named)

      if (identical(gee_family, "poisson")) {
        .wide_spline("gee_poisson_eta", "gee_px_mu",
                     spline_family = "poisson", num_intervals = 100L)
        .ring_scale("gee_poisson_eta", 0.5, "gee_px_eta_half")
        .wide_spline("gee_px_eta_half", "gee_px_sqrt_var",
                     spline_family = "poisson", num_intervals = 100L)
        .ring_scale("gee_poisson_eta", -0.5, "gee_px_eta_neghalf")
        .wide_spline("gee_px_eta_neghalf", "gee_px_inv_sqrt_var",
                     spline_family = "poisson", num_intervals = 100L)
        .copy_key("gee_px_mu", "secure_mu_share")
        for (ci in dcf_conns) {
          .dsAgg(datasources[ci],
            call(name = "k2PrepareWeightedResidualShareDS",
                 session_id = session_id))
        }
        .ring_affine(a_key = "k2_weight_residual_share_fp",
                     sign_a = -1L, output_key = "gee_px_y_minus_mu")
        .vecmul("k2_y_share_fp_original", "gee_px_inv_sqrt_var",
                "gee_px_y_inv_sqrt_var")
        .ring_affine(a_key = "gee_px_y_inv_sqrt_var",
                     b_key = "gee_px_sqrt_var",
                     sign_a = 1L, sign_b = -1L,
                     output_key = "gee_px_pearson")
        sqrt_var_key <- "gee_px_sqrt_var"
      } else {
        .wide_spline("gee_poisson_eta", "gee_px_mu",
                     spline_family = "sigmoid",
                     num_intervals = sigmoid_intervals)
        .ring_scale("gee_poisson_eta", 0.5, "gee_px_eta_half")
        .wide_spline("gee_px_eta_half", "gee_px_h",
                     spline_family = "poisson", num_intervals = 100L)
        .ring_scale("gee_poisson_eta", -0.5, "gee_px_eta_neghalf")
        .wide_spline("gee_px_eta_neghalf", "gee_px_inv_h",
                     spline_family = "poisson", num_intervals = 100L)
        .copy_key("gee_px_mu", "secure_mu_share")
        for (ci in dcf_conns) {
          .dsAgg(datasources[ci],
            call(name = "k2PrepareWeightedResidualShareDS",
                 session_id = session_id))
        }
        .ring_affine(a_key = "k2_weight_residual_share_fp",
                     sign_a = -1L, output_key = "gee_px_y_minus_mu")
        for (i in seq_along(dcf_parties)) {
          .dsAgg(datasources[dcf_conns[[i]]],
            call(name = "dsvertGLMMOneMinusMuDS",
                 output_key = "gee_px_one_minus_mu",
                 is_party0 = (i == 1L),
                 session_id = session_id,
                 frac_bits = frac_bits, ring = ring))
        }
        .vecmul("gee_px_mu", "gee_px_one_minus_mu", "gee_px_var")
        .vecmul("gee_px_mu", "gee_px_inv_h", "gee_px_sqrt_var")
        .ring_affine(a_key = "gee_px_h", b_key = "gee_px_inv_h",
                     sign_a = 1L, sign_b = 1L,
                     output_key = "gee_px_inv_sqrt_var")
        .vecmul("k2_y_share_fp_original", "gee_px_inv_sqrt_var",
                "gee_px_y_inv_sqrt_var")
        .ring_affine(a_key = "gee_px_y_inv_sqrt_var",
                     b_key = "gee_px_h",
                     sign_a = 1L, sign_b = -1L,
                     output_key = "gee_px_pearson")
        sqrt_var_key <- "gee_px_sqrt_var"
      }
      sqrt_var_key
    }

    compute_ar1_stats <- function(beta_current, rho_current) {
      sqrt_var_key <- prepare_family_shares(beta_current)
      .transform("gee_px_pearson", "gee_ar1_px_pearson_lead", "lead")
      .transform("gee_px_pearson", "gee_ar1_px_pearson_lag", "lag")
      .transform("gee_px_pearson", "gee_ar1_px_pearson_interior",
                 "interior")
      .transform("gee_px_pearson", "gee_ar1_px_pearson_nonlast",
                 "nonlast")

      e2 <- .prod_sum("gee_px_pearson", "gee_px_pearson",
                      "gee_ar1_px_e2")
      enonlast2 <- .prod_sum("gee_ar1_px_pearson_nonlast",
                             "gee_px_pearson",
                             "gee_ar1_px_enonlast2")
      eadj <- .prod_sum("gee_px_pearson", "gee_ar1_px_pearson_lead",
                        "gee_ar1_px_eadj")

      wx_keys <- paste0("gee_ar1_px_wx_", seq_len(q) - 1L)
      wx_lead_keys <- paste0("gee_ar1_px_wx_lead_", seq_len(q) - 1L)
      wx_interior_keys <- paste0("gee_ar1_px_wx_interior_", seq_len(q) - 1L)
      .copy_key(sqrt_var_key, wx_keys[[1L]])
      if (q > 1L) {
        for (j in 2:q) {
          .vecmul(x_keys[[j]], sqrt_var_key, wx_keys[[j]])
        }
      }
      for (j in seq_len(q)) {
        .transform(wx_keys[[j]], wx_lead_keys[[j]], "lead")
        .transform(wx_keys[[j]], wx_interior_keys[[j]], "interior")
      }

      total <- interior <- adj <- matrix(0, q, q)
      for (j in seq_len(q)) {
        for (k in j:q) {
          total[j, k] <- .prod_sum(wx_keys[[j]], wx_keys[[k]],
                                   paste0("gee_ar1_px_total_", j, "_", k))
          interior[j, k] <- .prod_sum(
            wx_interior_keys[[j]], wx_keys[[k]],
            paste0("gee_ar1_px_interior_", j, "_", k))
          total[k, j] <- total[j, k]
          interior[k, j] <- interior[j, k]
        }
        for (k in seq_len(q)) {
          adj[j, k] <- .prod_sum(wx_keys[[j]], wx_lead_keys[[k]],
                                 paste0("gee_ar1_px_adj_", j, "_", k))
        }
        if (verbose) {
          message(sprintf("[ds.vertGEE] %s AR1 bread column %d/%d",
                          family_title, j, q))
        }
      }
      cfac <- 1 / (1 - rho_current^2)
      A <- cfac * (total + rho_current^2 * interior -
                     rho_current * (adj + t(adj)))
      A <- (A + t(A)) / 2

      .ring_scale("gee_px_pearson", cfac, "gee_ar1_px_z_base")
      .ring_scale("gee_ar1_px_pearson_interior", cfac * rho_current^2,
                  "gee_ar1_px_z_interior")
      .ring_scale("gee_ar1_px_pearson_lead", -cfac * rho_current,
                  "gee_ar1_px_z_lead")
      .ring_scale("gee_ar1_px_pearson_lag", -cfac * rho_current,
                  "gee_ar1_px_z_lag")
      .ring_affine(a_key = "gee_ar1_px_z_base",
                   b_key = "gee_ar1_px_z_interior",
                   sign_a = 1L, sign_b = 1L,
                   output_key = "gee_ar1_px_z_tmp1")
      .ring_affine(a_key = "gee_ar1_px_z_tmp1",
                   b_key = "gee_ar1_px_z_lead",
                   sign_a = 1L, sign_b = 1L,
                   output_key = "gee_ar1_px_z_tmp2")
      .ring_affine(a_key = "gee_ar1_px_z_tmp2",
                   b_key = "gee_ar1_px_z_lag",
                   sign_a = 1L, sign_b = 1L,
                   output_key = "gee_ar1_px_rinv_pearson")
      .vecmul(sqrt_var_key, "gee_ar1_px_rinv_pearson",
              "gee_ar1_px_score_z")

      scores <- NULL
      for (j in seq_len(q)) {
        score_key <- if (j == 1L) {
          "gee_ar1_px_score_z"
        } else {
          key <- paste0("gee_ar1_px_score_x_", j)
          .vecmul(x_keys[[j]], "gee_ar1_px_score_z", key)
          key
        }
        cs <- .cluster_sum(score_key)
        if (is.null(scores)) {
          scores <- matrix(0, nrow = length(cs$values), ncol = q)
          cluster_sizes <<- cs$cluster_sizes
        }
        scores[, j] <- cs$values
        if (verbose) {
          message(sprintf("[ds.vertGEE] %s AR1 score column %d/%d",
                          family_title, j, q))
        }
      }
      list(A = A, U = colSums(scores), scores = scores,
           e2 = e2, enonlast2 = enonlast2, eadj = eadj,
           scale = e2 / max(n_obs, 1),
           cluster_sizes = cluster_sizes)
    }

    max_iter <- as.integer(max_iter)
    iter <- 0L
    for (iter in seq_len(max_iter)) {
      beta_old <- beta
      rho_old <- rho
      last_stats <- compute_ar1_stats(beta, rho)
      rho <- if (is.finite(last_stats$enonlast2) &&
                 last_stats$enonlast2 > 0) {
        last_stats$eadj / last_stats$enonlast2
      } else {
        0
      }
      rho <- min(max(rho, lower_rho), upper_rho)
      step <- as.numeric(safe_solve(last_stats$A, last_stats$U))
      if (any(!is.finite(step))) {
        stop("non-finite ", gee_family, " AR1 Newton step",
             call. = FALSE)
      }
      step_norm <- max(abs(step))
      if (step_norm > 1) step <- step / step_norm
      beta <- beta + step
      if (verbose) {
        message(sprintf(
          "[ds.vertGEE] %s AR1 iter %d rho=%.5g step=%.3e",
          family_title, iter, rho, max(abs(step))))
      }
      if (max(abs(beta - beta_old), abs(rho - rho_old)) <= tol) {
        converged <- TRUE
        break
      }
    }
    if (is.null(last_stats)) last_stats <- compute_ar1_stats(beta, rho)
    final_stats <- compute_ar1_stats(beta, rho)
    bread <- safe_solve(final_stats$A)
    scale <- final_stats$scale
    if (!is.finite(scale) || scale <= 0) scale <- 1
    model_cov_std <- as.numeric(scale) * bread
    robust_cov_std <- bread %*% crossprod(final_stats$scores) %*% bread
    model_cov_std <- (model_cov_std + t(model_cov_std)) / 2
    robust_cov_std <- (robust_cov_std + t(robust_cov_std)) / 2

    beta_target <- beta[perm]
    names(beta_target) <- target_order
    dimnames(model_cov_std) <- list(canonical_order, canonical_order)
    dimnames(robust_cov_std) <- list(canonical_order, canonical_order)
    model_cov_std <- model_cov_std[perm, perm, drop = FALSE]
    robust_cov_std <- robust_cov_std[perm, perm, drop = FALSE]

    J <- .ds_gee_standardization_jacobian(fit, target_features)
    model_cov <- J %*% model_cov_std %*% t(J)
    robust_cov <- J %*% robust_cov_std %*% t(J)
    model_cov <- (model_cov + t(model_cov)) / 2
    robust_cov <- (robust_cov + t(robust_cov)) / 2
    dimnames(model_cov) <- list(target_order, target_order)
    dimnames(robust_cov) <- list(target_order, target_order)

    coefficients <- .ds_gee_unstandardized_parameters(
      fit, beta_target, target_features)
    model_se <- sqrt(pmax(diag(model_cov), 0))
    robust_se <- sqrt(pmax(diag(robust_cov), 0))
    names(model_se) <- names(robust_se) <- target_order

    return(list(
      coefficients = coefficients,
      model_se = model_se,
      robust_se = robust_se,
      model_covariance = model_cov,
      robust_covariance = robust_cov,
      method = paste0("share_domain_ar1_", gee_family),
      working_correlation = list(
        corstr = "ar1",
        alpha = as.numeric(rho),
        phi = as.numeric(scale),
        iterations = as.integer(iter),
        converged = isTRUE(converged),
        estimator = paste0(gee_family,
                           "_ar1_guarded_adjacent_pearson_stats"),
        disclosure = paste0(
          "guarded adjacent-pair Pearson score sufficient statistics only; ",
          "no row-level mu, residuals, scores, visit labels, order vectors, ",
          "or cluster labels returned")),
      cluster_sizes = cluster_sizes))
  }

  alpha <- 0
  converged <- FALSE
  cluster_sizes <- NULL
  last_stats <- NULL
  lower_alpha <- -Inf
  upper_alpha <- 0.95

  compute_stats <- function(beta_current) {
    beta_named <- beta_current
    names(beta_named) <- canonical_order
    .restore_eta_shape(beta_named)

    if (identical(gee_family, "poisson")) {
      .wide_spline("gee_poisson_eta", "gee_px_mu",
                   spline_family = "poisson", num_intervals = 100L)
      .ring_scale("gee_poisson_eta", 0.5, "gee_px_eta_half")
      .wide_spline("gee_px_eta_half", "gee_px_sqrt_var",
                   spline_family = "poisson", num_intervals = 100L)
      .ring_scale("gee_poisson_eta", -0.5, "gee_px_eta_neghalf")
      .wide_spline("gee_px_eta_neghalf", "gee_px_inv_sqrt_var",
                   spline_family = "poisson", num_intervals = 100L)
      .copy_key("gee_px_mu", "secure_mu_share")
      for (ci in dcf_conns) {
        .dsAgg(datasources[ci],
          call(name = "k2PrepareWeightedResidualShareDS",
               session_id = session_id))
      }
      .ring_affine(a_key = "k2_weight_residual_share_fp",
                   sign_a = -1L, output_key = "gee_px_y_minus_mu")
      .vecmul("k2_y_share_fp_original", "gee_px_inv_sqrt_var",
              "gee_px_y_inv_sqrt_var")
      .ring_affine(a_key = "gee_px_y_inv_sqrt_var",
                   b_key = "gee_px_sqrt_var",
                   sign_a = 1L, sign_b = -1L,
                   output_key = "gee_px_pearson")
      var_key <- "gee_px_mu"
    } else {
      .wide_spline("gee_poisson_eta", "gee_px_mu",
                   spline_family = "sigmoid",
                   num_intervals = sigmoid_intervals)
      .ring_scale("gee_poisson_eta", 0.5, "gee_px_eta_half")
      .wide_spline("gee_px_eta_half", "gee_px_h",
                   spline_family = "poisson", num_intervals = 100L)
      .ring_scale("gee_poisson_eta", -0.5, "gee_px_eta_neghalf")
      .wide_spline("gee_px_eta_neghalf", "gee_px_inv_h",
                   spline_family = "poisson", num_intervals = 100L)
      .copy_key("gee_px_mu", "secure_mu_share")
      for (ci in dcf_conns) {
        .dsAgg(datasources[ci],
          call(name = "k2PrepareWeightedResidualShareDS",
               session_id = session_id))
      }
      .ring_affine(a_key = "k2_weight_residual_share_fp",
                   sign_a = -1L, output_key = "gee_px_y_minus_mu")
      for (i in seq_along(dcf_parties)) {
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "dsvertGLMMOneMinusMuDS",
               output_key = "gee_px_one_minus_mu",
               is_party0 = (i == 1L),
               session_id = session_id,
               frac_bits = frac_bits, ring = ring))
      }
      .vecmul("gee_px_mu", "gee_px_one_minus_mu", "gee_px_var")
      .vecmul("gee_px_mu", "gee_px_inv_h", "gee_px_sqrt_var")
      .ring_affine(a_key = "gee_px_h", b_key = "gee_px_inv_h",
                   sign_a = 1L, sign_b = 1L,
                   output_key = "gee_px_inv_sqrt_var")
      .vecmul("k2_y_share_fp_original", "gee_px_inv_sqrt_var",
              "gee_px_y_inv_sqrt_var")
      .ring_affine(a_key = "gee_px_y_inv_sqrt_var",
                   b_key = "gee_px_h",
                   sign_a = 1L, sign_b = -1L,
                   output_key = "gee_px_pearson")
      var_key <- "gee_px_var"
    }
    .vecmul("gee_px_pearson", "gee_px_pearson", "gee_px_pearson2")

    s1 <- sxv <- matrix(0, nrow = 1L, ncol = q)
    for (j in seq_len(q)) {
      key_s1 <- paste0("gee_px_s1_", j)
      .vecmul(x_keys[[j]], "gee_px_y_minus_mu", key_s1)
      cs1 <- .cluster_sum(key_s1)
      if (j == 1L) {
        s1 <- matrix(0, nrow = length(cs1$values), ncol = q)
        sxv <- matrix(0, nrow = length(cs1$values), ncol = q)
        cluster_sizes <<- cs1$cluster_sizes
      }
      s1[, j] <- cs1$values
      key_sxv <- paste0("gee_px_sxv_", j)
      .vecmul(x_keys[[j]], "gee_px_sqrt_var", key_sxv)
      sxv[, j] <- .cluster_sum(key_sxv)$values
    }
    se <- .cluster_sum("gee_px_pearson")$values
    e2 <- .cluster_sum("gee_px_pearson2")$values
    var_sum <- .cluster_sum(var_key)$values

    vxx <- array(0, dim = c(nrow(s1), q, q))
    for (j in seq_len(q)) {
      for (k in j:q) {
        if (j == 1L && k == 1L) {
          vals <- var_sum
        } else {
          key <- paste0("gee_px_vxx_", j, "_", k)
          .vecmul(xx_base_keys[j, k], var_key, key)
          vals <- .cluster_sum(key)$values
        }
        vxx[, j, k] <- vals
        vxx[, k, j] <- vals
      }
    }
    list(s1 = s1, sxv = sxv, se = se, e2 = e2, vxx = vxx,
         cluster_sizes = cluster_sizes)
  }

  build_system <- function(stats, alpha_current) {
    sizes <- as.numeric(stats$cluster_sizes)
    U <- numeric(q)
    A <- matrix(0, q, q)
    scores <- matrix(0, nrow = nrow(stats$s1), ncol = q)
    for (cc in seq_len(nrow(stats$s1))) {
      m <- sizes[cc]
      a <- 1 / (1 - alpha_current)
      b <- alpha_current / ((1 - alpha_current) *
                              (1 + alpha_current * (m - 1)))
      sc <- a * stats$s1[cc, ] - b * stats$sxv[cc, ] * stats$se[cc]
      Ac <- a * stats$vxx[cc, , ] - b * tcrossprod(stats$sxv[cc, ])
      U <- U + sc
      A <- A + Ac
      scores[cc, ] <- sc
    }
    lam <- as.numeric(lambda %||% 0)
    if (is.finite(lam) && lam > 0) {
      U <- U - lam * beta
      A <- A + lam * diag(q)
    }
    list(U = U, A = (A + t(A)) / 2, scores = scores)
  }

  max_iter <- as.integer(max_iter)
  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    alpha_old <- alpha
    last_stats <- compute_stats(beta)
    sizes <- as.numeric(last_stats$cluster_sizes)
    lower_alpha <- -1 / (max(sizes) - 1) + 1e-6
    scale <- sum(last_stats$e2) / max(n_obs, 1)
    pair_num <- sum((last_stats$se^2 - last_stats$e2) / 2)
    pair_den <- sum(sizes * (sizes - 1) / 2) * scale
    alpha <- if (is.finite(pair_den) && pair_den > 0) {
      pair_num / pair_den
    } else {
      0
    }
    alpha <- min(max(alpha, lower_alpha), upper_alpha)
    sys <- build_system(last_stats, alpha)
    step <- as.numeric(safe_solve(sys$A, sys$U))
    if (any(!is.finite(step))) {
      stop("non-finite ", gee_family, " exchangeable Newton step",
           call. = FALSE)
    }
    step_norm <- max(abs(step))
    if (step_norm > 1) step <- step / step_norm
    beta <- beta + step
    if (verbose) {
      message(sprintf(
        "[ds.vertGEE] %s exchangeable iter %d alpha=%.5g step=%.3e",
        family_title, iter, alpha, max(abs(step))))
    }
    if (max(abs(beta - beta_old), abs(alpha - alpha_old)) <= tol) {
      converged <- TRUE
      break
    }
  }
  if (is.null(last_stats)) last_stats <- compute_stats(beta)
  sys <- build_system(last_stats, alpha)
  bread <- safe_solve(sys$A)
  scale <- sum(last_stats$e2) / max(n_obs, 1)
  model_cov_std <- as.numeric(scale) * bread
  robust_cov_std <- bread %*% crossprod(sys$scores) %*% bread
  model_cov_std <- (model_cov_std + t(model_cov_std)) / 2
  robust_cov_std <- (robust_cov_std + t(robust_cov_std)) / 2

  beta_target <- beta[perm]
  names(beta_target) <- target_order
  dimnames(model_cov_std) <- list(canonical_order, canonical_order)
  dimnames(robust_cov_std) <- list(canonical_order, canonical_order)
  model_cov_std <- model_cov_std[perm, perm, drop = FALSE]
  robust_cov_std <- robust_cov_std[perm, perm, drop = FALSE]

  J <- .ds_gee_standardization_jacobian(fit, target_features)
  model_cov <- J %*% model_cov_std %*% t(J)
  robust_cov <- J %*% robust_cov_std %*% t(J)
  model_cov <- (model_cov + t(model_cov)) / 2
  robust_cov <- (robust_cov + t(robust_cov)) / 2
  dimnames(model_cov) <- list(target_order, target_order)
  dimnames(robust_cov) <- list(target_order, target_order)

  coefficients <- .ds_gee_unstandardized_parameters(
    fit, beta_target, target_features)
  model_se <- sqrt(pmax(diag(model_cov), 0))
  robust_se <- sqrt(pmax(diag(robust_cov), 0))
  names(model_se) <- names(robust_se) <- target_order

  list(coefficients = coefficients,
       model_se = model_se,
       robust_se = robust_se,
       model_covariance = model_cov,
       robust_covariance = robust_cov,
       method = paste0("share_domain_exchangeable_", gee_family),
       working_correlation = list(
         corstr = "exchangeable",
         alpha = as.numeric(alpha),
         phi = as.numeric(scale),
         iterations = as.integer(iter),
         converged = isTRUE(converged),
         estimator = paste0(gee_family,
                            "_exchangeable_pearson_cluster_stats"),
         disclosure = paste0(
           "guarded cluster-level Pearson score sufficient statistics only; ",
           "no row-level mu, residuals, scores, or cluster labels returned")),
       cluster_sizes = cluster_sizes)
}

#' @keywords internal
.ds_gee_secure_hc0 <- function(fit, datasources, server_names,
                               lambda = 0, data = NULL, id_col = NULL,
                               verbose = FALSE,
                               binomial_sigmoid_intervals = NULL) {
  if (!fit$family %in% c("gaussian", "binomial", "poisson")) return(NULL)
  required <- c("session_id", "transport_pks", "server_list", "x_vars",
                "y_server", "hessian_std", "x_sds", "x_means")
  missing_req <- required[vapply(required, function(nm) is.null(fit[[nm]]),
                                 logical(1L))]
  if (length(missing_req) > 0L) {
    stop("GLM session metadata missing: ", paste(missing_req, collapse = ", "),
         call. = FALSE)
  }

  session_id <- fit$session_id
  server_list <- fit$server_list
  x_vars <- fit$x_vars
  coordinator <- fit$y_server
  n_obs <- as.integer(fit$n_obs)
  ring <- as.integer(fit$ring %||% 63L)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  frac_bits <- if (ring == 127L) 50L else 20L
  ring_tag <- if (ring == 127L) "ring127" else "ring63"
  transport_pks <- fit$transport_pks
  sigmoid_intervals <- if (identical(fit$family, "binomial")) {
    suppressWarnings(as.integer(
      (binomial_sigmoid_intervals %||%
         getOption("dsvert.gee_binomial_sigmoid_intervals",
                   getOption("dsvert.glm_num_intervals_binomial",
                             getOption("dsvert.glm_num_intervals", 100L))))[[1L]]
    ))
  } else {
    NA_integer_
  }
  if (identical(fit$family, "binomial") &&
      (!is.finite(sigmoid_intervals) || sigmoid_intervals < 10L)) {
    sigmoid_intervals <- 100L
  }

  .to_b64url <- function(x) gsub("+", "-", gsub("/", "_",
    gsub("=+$", "", x, perl = TRUE), fixed = TRUE), fixed = TRUE)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-", "+", gsub("_", "/", x, fixed = TRUE), fixed = TRUE)
    pad <- nchar(x) %% 4
    if (pad == 2) x <- paste0(x, "==")
    if (pad == 3) x <- paste0(x, "=")
    x
  }
  .dsAgg <- function(conns, expr, ...) {
    tryCatch(
      DSI::datashield.aggregate(conns = conns, expr = expr, ...),
      error = function(e) {
        fn_name <- if (is.call(expr)) as.character(expr[[1]]) else "?"
        srv_name <- tryCatch(names(conns)[1], error = function(x) "?")
        stop(sprintf("%s on %s failed: %s", fn_name, srv_name,
                     conditionMessage(e)), call. = FALSE)
      })
  }
  .sendBlob <- function(blob, key, conn_idx) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        .dsAgg(datasources[conn_idx],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               session_id = session_id))
      } else {
        .dsAgg(datasources[conn_idx],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               chunk_index = chunk_idx, n_chunks = n_chunks,
               session_id = session_id))
      }
    })
  }

  target_features <- unlist(x_vars[server_list], use.names = FALSE)
  target_order <- c("(Intercept)", target_features)
  std <- .ds_gee_standardized_parameters(fit, target_features)
  cluster_enabled <- !is.null(id_col) && length(id_col) == 1L &&
    is.character(id_col) && nzchar(id_col)
  cluster_sizes <- NULL

  if (fit$eta_privacy == "k2_beaver") {
    nl <- setdiff(server_list, coordinator)
    if (length(nl) != 1L) {
      stop("K=2 GEE HC0 requires exactly one non-outcome server",
           call. = FALSE)
    }
    nl <- nl[[1L]]
    dcf_parties <- c(coordinator, nl)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    dealer_conn <- dcf_conns[[2L]]
    canonical_features <- c(x_vars[[coordinator]], x_vars[[nl]])
    b_coord <- as.numeric(std$beta[x_vars[[coordinator]]])
    b_nl <- as.numeric(std$beta[x_vars[[nl]]])
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "dsvertGEERestoreFeatureShapeDS",
             p_own = as.integer(length(x_vars[[srv]])),
             p_peer = as.integer(length(x_vars[[setdiff(dcf_parties, srv)]])),
             session_id = session_id))
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2ComputeEtaShareDS",
             beta_coord = b_coord, beta_nl = b_nl,
             intercept = if (srv == coordinator) std$intercept else 0,
             is_coordinator = (srv == coordinator),
             session_id = session_id))
    }
  } else if (fit$eta_privacy == "secure_agg") {
    fusion <- .k3_select_fusion_server(server_list, coordinator, x_vars)
    dcf_parties <- c(fusion, coordinator)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    non_dcf <- setdiff(server_list, dcf_parties)
    dealer <- if (length(non_dcf) > 0L) non_dcf[[1L]] else fusion
    dealer_conn <- which(server_names == dealer)
    p_coord <- length(x_vars[[coordinator]])
    p_fusion <- length(x_vars[[fusion]])
    p_extras <- sum(vapply(non_dcf, function(s) length(x_vars[[s]]),
                           integer(1L)))
    canonical_features <- c(x_vars[[coordinator]], x_vars[[fusion]])
    for (srv in non_dcf) canonical_features <- c(canonical_features, x_vars[[srv]])
    for (i in seq_along(dcf_parties)) {
      srv <- dcf_parties[[i]]
      is_coord <- srv == coordinator
      p_peer <- if (is_coord) p_fusion + p_extras else p_coord + p_extras
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "dsvertGEERestoreFeatureShapeDS",
             p_own = as.integer(length(x_vars[[srv]])),
             p_peer = as.integer(p_peer),
             session_id = session_id))
      b_coord <- as.numeric(std$beta[x_vars[[coordinator]]])
      if (is_coord) {
        b_nl <- c(as.numeric(std$beta[x_vars[[fusion]]]))
        for (ns in non_dcf) b_nl <- c(b_nl, as.numeric(std$beta[x_vars[[ns]]]))
      } else {
        b_nl <- numeric(0)
        for (ns in non_dcf) b_nl <- c(b_nl, as.numeric(std$beta[x_vars[[ns]]]))
        b_nl <- c(b_nl, as.numeric(std$beta[x_vars[[fusion]]]))
      }
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2ComputeEtaShareDS",
             beta_coord = b_coord, beta_nl = b_nl,
             intercept = if (is_coord) std$intercept else 0,
             is_coordinator = is_coord,
             session_id = session_id))
      if (!is_coord && p_extras > 0L) {
        .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "glmRing63ReorderXFullDS",
               p_coord = as.integer(p_coord),
               p_fusion = as.integer(p_fusion),
               p_extras = as.integer(p_extras),
               session_id = session_id))
      }
    }
  } else {
    stop("unsupported eta_privacy for secure GEE HC0: ", fit$eta_privacy,
         call. = FALSE)
  }

  if (identical(fit$family, "gaussian")) {
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2IdentityLinkDS", session_id = session_id))
    }
  } else if (fit$family %in% c("binomial", "poisson")) {
    spline_family <- fit$family
    dcf_family <- if (identical(fit$family, "poisson")) "poisson" else "sigmoid"
    num_intervals <- if (identical(fit$family, "poisson")) 100L else sigmoid_intervals
    dcf <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenDcfKeysDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           family = dcf_family, n = as.integer(n_obs),
           frac_bits = frac_bits, num_intervals = num_intervals,
           ring = ring, session_id = session_id))
    if (is.list(dcf) && length(dcf) == 1L) dcf <- dcf[[1L]]
    .sendBlob(dcf$dcf_blob_0, "k2_dcf_keys_persistent", dcf_conns[[1L]])
    .sendBlob(dcf$dcf_blob_1, "k2_dcf_keys_persistent", dcf_conns[[2L]])
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2StoreDcfKeysPersistentDS", session_id = session_id))
    }

    spline_t <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenSplineTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           n = as.integer(n_obs), frac_bits = frac_bits,
           ring = ring, session_id = session_id))
    if (is.list(spline_t) && length(spline_t) == 1L) {
      spline_t <- spline_t[[1L]]
    }
    .sendBlob(spline_t$spline_blob_0, "k2_spline_triples", dcf_conns[[1L]])
    .sendBlob(spline_t$spline_blob_1, "k2_spline_triples", dcf_conns[[2L]])
    for (ph in 1:4) {
      pr <- vector("list", 2L)
      for (i in seq_along(dcf_parties)) {
        r <- .dsAgg(datasources[dcf_conns[[i]]],
          call(paste0("k2WideSplinePhase", ph, "DS"),
               party_id = as.integer(i - 1L),
               family = spline_family,
               num_intervals = num_intervals,
               frac_bits = frac_bits,
               ring = ring,
               session_id = session_id))
        if (is.list(r) && length(r) == 1L) r <- r[[1L]]
        pr[[i]] <- r
      }
      if (ph == 1L) {
        .sendBlob(pr[[1L]]$dcf_masked, "k2_peer_dcf_masked",
                  dcf_conns[[2L]])
        .sendBlob(pr[[2L]]$dcf_masked, "k2_peer_dcf_masked",
                  dcf_conns[[1L]])
      } else if (ph == 2L) {
        for (i in 1:2) {
          peer_i <- 3L - i
          pk_b64 <- .b64url_to_b64(transport_pks[[dcf_parties[[peer_i]]]])
          payload <- jsonlite::toJSON(list(
            and_xma = pr[[i]]$and_xma, and_ymb = pr[[i]]$and_ymb,
            had1_xma = pr[[i]]$had1_xma, had1_ymb = pr[[i]]$had1_ymb),
            auto_unbox = TRUE)
          sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
            data = jsonlite::base64_enc(charToRaw(payload)),
            recipient_pk = pk_b64))
          .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1",
                    dcf_conns[[peer_i]])
        }
      } else if (ph == 3L) {
        for (i in 1:2) {
          peer_i <- 3L - i
          pk_b64 <- .b64url_to_b64(transport_pks[[dcf_parties[[peer_i]]]])
          payload <- jsonlite::toJSON(list(
            had2_xma = pr[[i]]$had2_xma,
            had2_ymb = pr[[i]]$had2_ymb),
            auto_unbox = TRUE)
          sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
            data = jsonlite::base64_enc(charToRaw(payload)),
            recipient_pk = pk_b64))
          .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1",
                    dcf_conns[[peer_i]])
        }
      }
    }
  }

  if (cluster_enabled) {
    if (is.null(data) || !is.character(data) || length(data) != 1L ||
        !nzchar(data)) {
      stop("data name is required for clustered GEE sandwich",
           call. = FALSE)
    }
    cluster_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                               data, id_col)
    if (is.null(cluster_srv)) {
      stop("id_col '", id_col, "' not found on any server", call. = FALSE)
    }
    if (!identical(cluster_srv, coordinator)) {
      stop("clustered GEE sandwich requires id_col to live on the ",
           "outcome server; found on '", cluster_srv, "'", call. = FALSE)
    }
    coord_i <- match(coordinator, dcf_parties)
    peer_i <- setdiff(seq_along(dcf_parties), coord_i)
    cb <- .dsAgg(datasources[dcf_conns[[coord_i]]],
      call(name = "dsvertClusterIDsBroadcastDS",
           data_name = data, cluster_col = id_col,
           peer_pk = transport_pks[[dcf_parties[[peer_i]]]],
           session_id = session_id))
    if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
    .sendBlob(cb$peer_blob, "dsvert_cluster_ids_blob", dcf_conns[[peer_i]])
    .dsAgg(datasources[dcf_conns[[peer_i]]],
      call(name = "dsvertClusterIDsReceiveDS", session_id = session_id))
  }

  p_total <- length(canonical_features)
  if (p_total < 1L) stop("no predictors available for GEE HC0", call. = FALSE)
  canonical_order <- c("(Intercept)", canonical_features)

  .vecmul <- function(x_key, y_key, output_key) {
    tri <- .dsAgg(datasources[dealer_conn],
      call(name = "k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           n = as.numeric(n_obs), session_id = session_id,
           frac_bits = frac_bits, ring = ring))
    if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple", dcf_conns[[1L]])
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple", dcf_conns[[2L]])
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2BeaverVecmulConsumeTripleDS",
             session_id = session_id))
    }
    r1 <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      peer <- dcf_parties[[3L - i]]
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR1DS",
             peer_pk = transport_pks[[peer]],
             x_key = x_key, y_key = y_key,
             n = as.numeric(n_obs), session_id = session_id,
             frac_bits = frac_bits, ring = ring))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[i]] <- r
    }
    .sendBlob(r1[[1L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              dcf_conns[[2L]])
    .sendBlob(r1[[2L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              dcf_conns[[1L]])
    for (i in seq_along(dcf_parties)) {
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR2DS",
             is_party0 = (i == 1L),
             x_key = x_key, y_key = y_key,
             output_key = output_key,
             n = as.numeric(n_obs), session_id = session_id,
             frac_bits = frac_bits, ring = ring))
    }
    invisible(output_key)
  }

  .gradient <- function() {
    gt <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenGradTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           n = as.integer(n_obs), p = as.integer(p_total),
           ring = ring, session_id = session_id))
    if (is.list(gt) && length(gt) == 1L) gt <- gt[[1L]]
    .sendBlob(gt$grad_blob_0, "k2_grad_triple_fp", dcf_conns[[1L]])
    .sendBlob(gt$grad_blob_1, "k2_grad_triple_fp", dcf_conns[[2L]])
    r1 <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      peer <- dcf_parties[[3L - i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2GradientR1DS",
             peer_pk = transport_pks[[peer]],
             session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[i]] <- r
    }
    .sendBlob(r1[[1L]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[[2L]])
    .sendBlob(r1[[2L]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[[1L]])
    r2 <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2GradientR2DS",
             party_id = as.integer(i - 1L), session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r2[[i]] <- r
    }
    g <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r2[[1L]]$gradient_fp,
      share_b = r2[[2L]]$gradient_fp,
      frac_bits = frac_bits, ring = ring_tag))
    s <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r1[[1L]]$sum_residual_fp,
      share_b = r1[[2L]]$sum_residual_fp,
      frac_bits = frac_bits, ring = ring_tag))
    c(s$values[1L], as.numeric(g$values))
  }

  for (ci in dcf_conns) {
    .dsAgg(datasources[ci],
      call(name = "k2PrepareWeightedResidualShareDS",
           session_id = session_id))
  }
  if (!cluster_enabled) {
    .vecmul("k2_weight_residual_share_fp",
            "k2_weight_residual_share_fp",
            "gee_r2_share")
  }

  for (i in seq_along(dcf_parties)) {
    .dsAgg(datasources[dcf_conns[[i]]],
      call(name = "dsvertGEEInterceptShareDS",
           output_key = "gee_x_col_0",
           n = as.numeric(n_obs),
           is_party0 = (i == 1L),
           session_id = session_id,
           frac_bits = frac_bits, ring = ring))
  }
  for (j in seq_len(p_total)) {
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2BeaverExtractColumnDS",
             source_key = "k2_x_full_fp",
             n = as.numeric(n_obs), K = as.numeric(p_total),
             col_index = as.numeric(j),
             output_key = paste0("gee_x_col_", j),
             session_id = session_id,
             frac_bits = frac_bits, ring = ring_tag))
    }
  }

  meat_method <- "share_domain_hc0"
  if (cluster_enabled) {
    S_canonical <- NULL
    for (j in 0:p_total) {
      x_key <- paste0("gee_x_col_", j)
      rx_key <- paste0("gee_r_x_col_", j)
      .vecmul("k2_weight_residual_share_fp", x_key, rx_key)
      parts <- vector("list", length(dcf_parties))
      for (i in seq_along(dcf_parties)) {
        part <- .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "dsvertPerClusterSumShareDS",
               share_key = rx_key, session_id = session_id,
               frac_bits = frac_bits, ring = ring))
        if (is.list(part) && length(part) == 1L) part <- part[[1L]]
        parts[[i]] <- part
      }
      n_clusters <- length(parts[[1L]]$per_cluster_fp)
      if (is.null(S_canonical)) {
        S_canonical <- matrix(0, nrow = n_clusters, ncol = p_total + 1L)
      }
      for (ck in seq_len(n_clusters)) {
        agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
          share_a = parts[[1L]]$per_cluster_fp[[ck]],
          share_b = parts[[2L]]$per_cluster_fp[[ck]],
          frac_bits = frac_bits, ring = ring_tag))
        S_canonical[ck, j + 1L] <- as.numeric(agg$values[1L])
      }
      coord_i <- match(coordinator, dcf_parties)
      cluster_sizes <- as.integer(parts[[coord_i]]$cluster_sizes)
      if (verbose) {
        message(sprintf("[ds.vertGEE] secure cluster meat column %d/%d",
                        j + 1L, p_total + 1L))
      }
    }
    A_canonical <- crossprod(S_canonical)
    meat_method <- "share_domain_cluster"
  } else {
    A_canonical <- matrix(0, nrow = p_total + 1L, ncol = p_total + 1L)
    for (j in 0:p_total) {
      x_key <- paste0("gee_x_col_", j)
      wx_key <- paste0("gee_r2_x_col_", j)
      .vecmul("gee_r2_share", x_key, wx_key)
      for (ci in dcf_conns) {
        .dsAgg(datasources[ci],
          call(name = "k2FinalizeWeightedResidualShareDS",
               input_key = wx_key, session_id = session_id))
      }
      A_canonical[, j + 1L] <- .gradient()
      if (verbose) {
        message(sprintf("[ds.vertGEE] secure HC0 meat column %d/%d",
                        j + 1L, p_total + 1L))
      }
    }
    A_canonical <- (A_canonical + t(A_canonical)) / 2
  }
  dimnames(A_canonical) <- list(canonical_order, canonical_order)

  perm <- match(target_order, canonical_order)
  if (anyNA(perm)) {
    stop("could not align GEE HC0 meat order", call. = FALSE)
  }
  A_std <- A_canonical[perm, perm, drop = FALSE]
  dimnames(A_std) <- list(target_order, target_order)

  H <- as.matrix(fit$hessian_std)
  if (!is.null(rownames(H))) {
    H <- H[target_order, target_order, drop = FALSE]
  }
  H_adj <- H - as.numeric(lambda %||% 0) * diag(nrow(H))
  fisher_std <- n_obs * H_adj
  cov_std <- tryCatch(solve(fisher_std), error = function(e) NULL)
  if (is.null(cov_std)) stop("standardized Fisher matrix is singular",
                             call. = FALSE)
  dimnames(cov_std) <- list(target_order, target_order)

  robust_std <- cov_std %*% A_std %*% cov_std
  J <- .ds_gee_standardization_jacobian(fit, target_features)
  robust_cov <- J %*% robust_std %*% t(J)
  robust_cov <- (robust_cov + t(robust_cov)) / 2
  dimnames(robust_cov) <- list(target_order, target_order)
  robust_se <- sqrt(pmax(diag(robust_cov), 0))
  names(robust_se) <- target_order
  list(robust_covariance = robust_cov,
       robust_se = robust_se,
       meat_std = A_std,
       cluster_sizes = cluster_sizes,
       method = meat_method)
}

#' @export
print.ds.vertGEE <- function(x, ...) {
  cat("dsVert GEE (", x$corstr, " working correlation)\n", sep = "")
  cat(sprintf("  Family: %s   N = %d\n", x$family, x$n_obs))
  cat("\nCoefficients:\n")
  df <- data.frame(
    Estimate = x$coefficients,
    `SE (model)` = x$model_se,
    `SE (robust)` = x$robust_se,
    z = x$coefficients / x$robust_se,
    check.names = FALSE)
  df$p <- 2 * stats::pnorm(-abs(df$z))
  print(round(df, 5L))
  invisible(x)
}
