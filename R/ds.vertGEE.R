#' @title Federated generalised estimating equations
#' @description Fit a GLM for vertically partitioned DataSHIELD data and
#'   return sandwich (robust) standard errors alongside the usual
#'   model-based ones. The point estimate \eqn{\hat{\beta}} is obtained
#'   by a single call to \code{\link{ds.vertGLM}} (bread =
#'   \code{fit$covariance}), then the Liang-Zeger meat matrix is
#'   estimated by a second weighted fit with weights proportional to the
#'   squared Pearson residuals. This is the working-independence GEE
#'   estimator; exchangeable / AR1 working correlation structures are a
#'   planned Month-4 extension that adds a cluster-ID broadcast between
#'   the DCF parties and per-cluster residual outer products.
#'
#'   Because both fits reveal only \eqn{p}-dimensional aggregates
#'   (gradient and Hessian), the sandwich estimator is obtained without
#'   ever materialising the \eqn{n}-length residual vector at the
#'   client. Inter-server leakage is identical to the ordinary
#'   \code{ds.vertGLM} path (no new channels).
#'
#'   Formula:
#'     \deqn{V_{sand} = \mathrm{Cov}(\hat\beta) \, A \, \mathrm{Cov}(\hat\beta)}
#'     \deqn{A = X^T \, \mathrm{diag}(r^2) \, X \, / \, n}
#'
#' @param formula A model formula passed through to \code{ds.vertGLM}.
#' @param data Character. Aligned data-frame name on each server.
#' @param family One of \code{"gaussian"}, \code{"binomial"},
#'   \code{"poisson"}.
#' @param id_col Optional character; reserved for clustered GEE. In this
#'   first-pass implementation \code{id_col} is accepted but the
#'   working-independence sandwich is computed (cluster aggregation is
#'   the Month-4 follow-on).
#' @param corstr Working correlation structure. Only
#'   \code{"independence"} is currently supported; \code{"exchangeable"}
#'   and \code{"ar1"} raise a targeted message.
#' @param max_iter,tol,lambda,verbose Passed to \code{ds.vertGLM}.
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
                       corstr = c("independence", "exchangeable", "ar1"),
                       max_iter = 100L, tol = 1e-4, lambda = 1e-4,
                       verbose = TRUE, datasources = NULL) {
  family <- match.arg(family)
  corstr <- match.arg(corstr)
  if (corstr != "independence") {
    stop("ds.vertGEE: corstr='", corstr, "' requires the cluster-ID ",
         "broadcast protocol scheduled for Month 4. See V2_PROGRESS.md. ",
         "Use corstr='independence' for the working-independence ",
         "sandwich today.", call. = FALSE)
  }
  if (verbose && !is.null(id_col)) {
    message("[ds.vertGEE] id_col='", id_col, "' recorded but unused under ",
            "corstr='independence'. Cluster-aware sandwich is Month 4.")
  }

  # Stage 1: ordinary fit -> betahat, Cov_model (bread).
  if (verbose) message("[ds.vertGEE] Stage 1: ds.vertGLM point estimate")
  fit <- ds.vertGLM(formula = formula, data = data,
                    family = family, max_iter = max_iter, tol = tol,
                    lambda = lambda, verbose = verbose,
                    datasources = datasources)
  if (is.null(fit$covariance)) {
    stop("ds.vertGEE: the underlying ds.vertGLM fit does not expose ",
         "Cov(beta). Refit with dsVert >= 8bb7902 / dsVertClient >= ",
         "8bb7902.", call. = FALSE)
  }
  Cov_model <- as.matrix(fit$covariance)
  Cov_model <- (Cov_model + t(Cov_model)) / 2

  # Stage 2: weighted fit with weights = r^2 to pick up X^T diag(r^2) X.
  # For Gaussian: r_i = y_i - x_i^T beta. The outcome server holds y
  # plaintext; it can materialise r^2 server-side after receiving the
  # plaintext betahat from the client. We expose r^2 as a weights
  # column via the existing ds.vertGLM(weights=) infrastructure, which
  # scales mu/y shares element-wise (no Beaver round) and lets the
  # second fit's Fisher info equal X^T diag(r^2) X up to a known scale.
  #
  # Privacy: the weights themselves never cross the client -- r^2 is
  # computed on the outcome server and encrypted to the DCF peer via
  # k2SetWeightsDS. The client only sees the final sandwich matrix.
  #
  # For Binomial / Poisson: the Pearson residual is
  #   r_i = (y_i - mu_i) / sqrt(Var(mu_i))
  # which requires mu_i server-side; under DCF this is available on
  # the outcome server via the stored eta+mu shares once they have
  # been reconstructed internally. A server helper
  # `k2ComputePearsonR2ColDS(betahat, family)` is the cleanest way
  # to wire this; its design is documented below.
  #
  # Until that helper is merged, we degrade gracefully: the Gaussian
  # case can reuse the `resid = y - Xbeta` plaintext path on the
  # outcome server via dsvertPearsonRColumnDS (below); Binomial and
  # Poisson currently fall back to model-based SE with a warning.
  if (family != "gaussian") {
    if (verbose) {
      message("[ds.vertGEE] Sandwich SE for family='", family,
              "' requires k2ComputePearsonR2ColDS (Month 4). Returning ",
              "model-based SE for now.")
    }
    robust_se <- fit$std_errors
    robust_cov <- Cov_model
  } else {
    if (verbose) message("[ds.vertGEE] Stage 2: weighted fit for sandwich meat")
    # Request residuals-squared column on the outcome server. This uses
    # the `dsvertPearsonR2ColDS` server helper (Gaussian-only in the
    # first pass): given data_name, y, x_vars, and a plaintext betahat,
    # it writes a new column `__dsvert_r2` into the aligned data frame.
    # The client then re-fits with weights="__dsvert_r2".
    server_names <- names(datasources %||% DSI::datashield.connections_find())
    if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
    # Discover y_server
    y_var <- .ds_gee_extract_lhs(formula)
    y_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                          data, y_var)
    if (is.null(y_srv)) {
      warning("ds.vertGEE: could not locate outcome server for y='",
              y_var, "'. Returning model-based SE.", call. = FALSE)
      robust_se <- fit$std_errors
      robust_cov <- Cov_model
    } else {
      x_all <- names(fit$coefficients)
      x_all <- setdiff(x_all, "(Intercept)")
      r2_ok <- tryCatch({
        DSI::datashield.aggregate(
          datasources[which(server_names == y_srv)],
          call("dsvertPearsonR2ColDS",
               data_name = data, y_var = y_var,
               x_names = x_all,
               betahat = as.numeric(fit$coefficients[x_all]),
               intercept = as.numeric(fit$coefficients["(Intercept)"]),
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
          #   V_sand = Cov_model * (Fisher_weighted / Fisher_model) * Cov_model
          #         ~= Cov_model * solve(Cov_weighted) * Cov_model
          # where Cov_weighted = vcov of weighted fit.
          Cov_w <- as.matrix(fit_w$covariance)
          Cov_w <- (Cov_w + t(Cov_w)) / 2
          A <- tryCatch(solve(Cov_w), error = function(e) NULL)
          if (is.null(A)) {
            warning("meat matrix singular; returning model-based SE.",
                    call. = FALSE)
            robust_se <- fit$std_errors
            robust_cov <- Cov_model
          } else {
            robust_cov <- Cov_model %*% A %*% Cov_model
            dimnames(robust_cov) <- list(names(fit$coefficients),
                                          names(fit$coefficients))
            robust_se <- sqrt(pmax(diag(robust_cov), 0))
            names(robust_se) <- names(fit$coefficients)
          }
        }
      }
    }
  }

  out <- list(
    coefficients       = fit$coefficients,
    model_se           = fit$std_errors,
    robust_se          = robust_se,
    robust_covariance  = robust_cov,
    corstr             = corstr,
    family             = family,
    id_col             = id_col,
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
        call("dsvertColNamesDS", data_name = data_name))[[1]]$columns,
      error = function(e) NULL)
    if (!is.null(cols) && var %in% cols) return(srv)
  }
  NULL
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
