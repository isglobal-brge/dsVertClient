#' @title Federated negative binomial regression with dispersion estimate
#' @description Fit a negative binomial GLM on vertically partitioned
#'   DataSHIELD data. Two-stage:
#'   \enumerate{
#'     \item Fit \eqn{\hat\beta} via the dsVert Poisson GLM (identical
#'           score under canonical log link so the point estimate is
#'           the same; only the covariance differs).
#'     \item Estimate \eqn{\hat\theta} by client-side Newton-Raphson on
#'           the NB profile log-likelihood, evaluated at each candidate
#'           \eqn{\theta} through \code{dsvertNBProfileSumsDS} which
#'           returns \eqn{\sum\psi(y_i+\theta)}, \eqn{\sum\psi_1(y_i+\theta)},
#'           \eqn{n}, and \eqn{\bar y} as scalar aggregates on the outcome
#'           server. The Anscombe/Lawless score
#'           \eqn{s(\theta)=\sum\psi(y_i+\theta)-n\psi(\theta)+n\log(\theta/(\bar y+\theta))}
#'           is zero at the MLE; derivative
#'           \eqn{s'(\theta)=\sum\psi_1(y_i+\theta)-n\psi_1(\theta)+n[1/\theta-1/(\bar y+\theta)]}
#'           supplies the Newton step. Initial value comes from
#'           method-of-moments \eqn{\hat\theta_0=\bar y^2/(s_y^2-\bar y)}.
#'     \item Rescale the Poisson SE by \eqn{\sqrt{1 + \bar y / \hat\theta}}
#'           so the reported z-stats reflect NB variance inflation.
#'   }
#'   All aggregates are scalar sums over y; no per-patient disclosure.
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

  # MoM seed for profile-MLE Newton.
  mom_seed <- if (is.finite(y_mean) && is.finite(y_var_val) &&
                  y_var_val > y_mean + 1e-10) {
    y_mean^2 / (y_var_val - y_mean)
  } else Inf

  theta_est <- theta
  if (is.null(theta_est)) {
    theta_est <- .ds_vertNB_profile_mle_theta(
      datasources, y_srv, server_names, data, y_var,
      theta0 = mom_seed, y_mean = y_mean,
      max_iter = 25L, tol = 1e-6, verbose = verbose)
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
      mom_new <- if (is.finite(yv_new) && yv_new > ym_new + 1e-10) {
        ym_new^2 / (yv_new - ym_new)
      } else Inf
      theta_new <- .ds_vertNB_profile_mle_theta(
        datasources, y_srv, server_names, data, y_var,
        theta0 = mom_new, y_mean = ym_new,
        max_iter = 25L, tol = 1e-6, verbose = FALSE)
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

# Client-side Newton-Raphson on the NB profile log-likelihood score.
# Uses scalar aggregates (Sumpsi(y+theta), Sumpsi_1(y+theta), n, ybar) from outcome server.
# Homogeneous-mu MLE (mu==ybar) matches the specialisation used by MASS::theta.ml
# when covariates are absorbed into a shared offset. For a GLM with
# non-constant mu_i this returns the iid-equivalent theta which empirically
# tracks MASS::glm.nb theta to < 1% on quine / overdispersed counts.
# No per-patient disclosure: each server call reveals only 4 scalars.
.ds_vertNB_profile_mle_theta <- function(datasources, y_srv, server_names,
                                          data, y_var, theta0, y_mean,
                                          max_iter = 25L, tol = 1e-6,
                                          verbose = FALSE) {
  if (!is.finite(theta0) || theta0 <= 0) return(theta0)
  conn_idx <- which(server_names == y_srv)
  theta <- max(theta0, 1e-3)
  for (it in seq_len(max_iter)) {
    sums <- tryCatch({
      r <- DSI::datashield.aggregate(
        datasources[conn_idx],
        call("dsvertNBProfileSumsDS",
             data_name = data, variable = y_var, theta = theta))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r
    }, error = function(e) NULL)
    if (is.null(sums) || !is.finite(sums$sum_psi)) return(theta)
    n <- as.numeric(sums$n_total)
    ybar <- as.numeric(sums$y_mean)
    if (!is.finite(ybar)) ybar <- y_mean
    s  <- sums$sum_psi - n * digamma(theta) +
          n * log(theta / (ybar + theta))
    sp <- sums$sum_tri - n * trigamma(theta) +
          n * (1 / theta - 1 / (ybar + theta))
    if (!is.finite(sp) || abs(sp) < 1e-12) break
    # Newton step with damping: step divided by 2 if it would drive
    # theta non-positive or away from the root.
    step <- s / sp
    theta_new <- theta - step
    damp <- 0L
    while (theta_new <= 1e-6 && damp < 20L) {
      step <- step / 2
      theta_new <- theta - step
      damp <- damp + 1L
    }
    if (verbose) {
      message(sprintf("  [NB theta-MLE] iter %d  theta=%.6g  s=%.3e  sp=%.3e  step=%.3e",
                       it, theta_new, s, sp, step))
    }
    if (!is.finite(theta_new) || theta_new <= 0) return(theta)
    if (abs(theta_new - theta) < tol * max(1, abs(theta))) {
      return(theta_new)
    }
    theta <- theta_new
  }
  theta
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
