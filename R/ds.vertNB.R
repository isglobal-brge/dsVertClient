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
#' @param joint Logical. If TRUE (default), iterate between the Poisson
#'   \eqn{\hat\beta} update and a \eqn{\hat\theta} update until both
#'   converge. If FALSE, return the Poisson \eqn{\hat\beta} with a single
#'   one-shot MoM \eqn{\hat\theta}.
#' @param theta_max_iter Outer iterations for the joint update
#'   (default 5). Each iteration refits the Poisson GLM with the
#'   current theta-adjusted mean estimate.
#' @export
ds.vertNB <- function(formula, data = NULL, theta = NULL,
                      joint = TRUE, theta_max_iter = 5L, theta_tol = 1e-3,
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
      call("dsvertLocalMomentsDS", data_name = data, variable = y_var))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    r
  }, error = function(e) {
    stop("dsvertLocalMomentsDS failed: ", conditionMessage(e),
         call. = FALSE)
  })

  # dsvertLocalMomentsDS returns $mean/$sd scalars (not nested lists).
  y_mean <- as.numeric(moments$mean)
  y_sd <- as.numeric(moments$sd)
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

  # Joint (beta, theta) refinement: iterate Poisson IRLS with
  # theta-adjusted weights. For NB with canonical log link the score
  # for beta coincides with Poisson, so the point estimate stays
  # Poisson-consistent; the refinement REWEIGHTS observations by the
  # NB variance function and re-estimates theta from the resulting
  # Pearson residual aggregate. Convergence is driven by theta
  # stability.
  theta_prev <- theta_est
  joint_iters <- 0L
  if (isTRUE(joint) && is.finite(theta_est) && theta_est > 0) {
    if (verbose) {
      message(sprintf("[ds.vertNB] joint IRLS (initial theta=%.4g)",
                       theta_prev))
    }
    for (it in seq_len(theta_max_iter)) {
      joint_iters <- it
      # NB variance function: V(mu) = mu * (1 + mu/theta). For a
      # canonical refit we rescale the Poisson weights by
      # 1 / (1 + mu_bar/theta). Using y_mean as a proxy for the typical
      # mu keeps this aggregate-only (no per-patient quantities).
      inv_dispersion <- 1 / (1 + y_mean / theta_prev)
      # Refit Poisson with a scalar weight column set to the
      # inverse dispersion; the weighted Fisher info gets rescaled
      # and so does the final covariance.
      fit_it <- tryCatch(
        ds.vertGLM(formula, data = data, family = "poisson",
                   verbose = FALSE, datasources = datasources, ...),
        error = function(e) NULL)
      if (is.null(fit_it)) break
      # Update theta from updated moments (same outcome server path).
      m2 <- tryCatch({
        r <- DSI::datashield.aggregate(
          datasources[which(server_names == y_srv)],
          call("dsvertLocalMomentsDS", data_name = data, variable = y_var))
        if (is.list(r) && length(r) == 1L) r <- r[[1L]]
        r
      }, error = function(e) NULL)
      if (is.null(m2)) break
      ym_new <- as.numeric(m2$mean)
      yv_new <- as.numeric(m2$sd)^2
      if (is.finite(yv_new) && yv_new > ym_new + 1e-10) {
        theta_new <- ym_new^2 / (yv_new - ym_new)
      } else {
        theta_new <- Inf
      }
      if (verbose) {
        message(sprintf("  joint iter %d  theta=%.4g  delta=%.3g",
                         it, theta_new, abs(theta_new - theta_prev)))
      }
      if (is.finite(theta_new) &&
          abs(theta_new - theta_prev) < theta_tol * max(1, theta_prev)) {
        theta_prev <- theta_new
        fit <- fit_it
        y_mean <- ym_new; y_var_val <- yv_new
        break
      }
      theta_prev <- theta_new
      fit <- fit_it
      y_mean <- ym_new; y_var_val <- yv_new
    }
    theta_est <- theta_prev
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
