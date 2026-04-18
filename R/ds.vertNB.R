#' @title Federated negative binomial regression with dispersion estimate
#' @description Fit a negative binomial GLM on vertically partitioned
#'   DataSHIELD data. Two-stage:
#'   \enumerate{
#'     \item Fit \eqn{\hat\beta} via the dsVert Poisson GLM (identical
#'           score under canonical log link so the point estimate is
#'           the same; only the covariance differs).
#'     \item Estimate \eqn{\hat\theta} server-side via method-of-moments
#'           on the outcome variable: \eqn{\hat\theta = \bar y^2 /
#'           (s_y^2 - \bar y)} when \eqn{s_y^2 > \bar y}, else
#'           \eqn{+\infty} (Poisson limit).
#'     \item Rescale the Poisson SE by \eqn{\sqrt{1 + \bar y / \hat\theta}}
#'           so the reported z-stats reflect NB variance inflation.
#'   }
#'   \eqn{\bar y} and \eqn{s_y^2} are already-aggregated moments (via
#'   the shipped \code{dsvertLocalMomentsDS}); no per-patient disclosure.
#' @export
ds.vertNB <- function(formula, data = NULL, theta = NULL,
                      verbose = TRUE, datasources = NULL, ...) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (verbose) message("[ds.vertNB] Stage 1: Poisson GLM point estimator")
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

  moments <- tryCatch({
    r <- DSI::datashield.aggregate(
      datasources[which(server_names == y_srv)],
      call("dsvertLocalMomentsDS", data_name = data, x_vars = y_var))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    r
  }, error = function(e) {
    stop("dsvertLocalMomentsDS failed: ", conditionMessage(e),
         call. = FALSE)
  })

  y_mean <- as.numeric(moments$means[[y_var]])
  y_sd <- as.numeric(moments$sds[[y_var]])
  y_var_val <- y_sd^2

  theta_est <- theta
  if (is.null(theta_est)) {
    if (is.finite(y_mean) && is.finite(y_var_val) &&
        y_var_val > y_mean + 1e-10) {
      theta_est <- y_mean^2 / (y_var_val - y_mean)
    } else {
      theta_est <- Inf
    }
  }

  var_inflation <- if (is.finite(theta_est) && theta_est > 0) {
    sqrt(1 + y_mean / theta_est)
  } else 1
  nb_se <- fit$std_errors * var_inflation
  nb_z <- fit$coefficients / nb_se
  nb_p <- 2 * stats::pnorm(-abs(nb_z))
  nb_cov <- if (!is.null(fit$covariance)) fit$covariance * var_inflation^2 else NULL

  out <- list(
    coefficients = fit$coefficients,
    std_errors   = nb_se,
    z_values     = nb_z,
    p_values     = nb_p,
    covariance   = nb_cov,
    theta        = theta_est,
    y_mean       = y_mean,
    y_var        = y_var_val,
    var_inflation = var_inflation,
    family       = "negbin",
    n_obs        = fit$n_obs,
    deviance     = fit$deviance,
    iterations   = fit$iterations,
    converged    = fit$converged,
    poisson_fit  = fit,
    call         = match.call())
  class(out) <- c("ds.vertNB", "ds.glm", "list")
  out
}

#' @export
print.ds.vertNB <- function(x, ...) {
  cat("dsVert negative-binomial regression\n")
  cat(sprintf("  N = %d   theta = %.4g (var inflation = %.3f)\n",
              x$n_obs, x$theta, x$var_inflation))
  df <- data.frame(
    Estimate = x$coefficients,
    SE       = x$std_errors,
    z        = x$z_values,
    p        = x$p_values,
    check.names = FALSE)
  print(round(df, 5L))
  invisible(x)
}
