#' @title Internal binomial GLMM Laplace development route
#' @description Fit a binomial random-intercept GLMM on vertically
#'   partitioned DataSHIELD data using a Laplace marginal-likelihood target.
#'
#'   This is an internal development route for an explicit
#'   \code{glmer}-like likelihood target. The public \code{ds.vert.glmm()}
#'   frontdoor uses the aggregate PQL route because it is the validated GLMM
#'   surface. The Laplace target is the random-intercept binomial
#'   approximation used by
#'   \code{lme4::glmer(..., family = binomial(), nAGQ = 1)}. The default
#'   estimator uses one protected Newton mode step per objective evaluation;
#'   setting \code{mode_max_iter > 1} trades runtime for more exact cluster
#'   mode solving.
#'
#'   \deqn{
#'     \ell(\beta,\sigma_b^2) =
#'       \sum_i \left[
#'         y_i^\top \eta_i - \sum_j \log(1 + e^{\eta_{ij}})
#'         - b_i^2/(2\sigma_b^2) - \frac12\log(\sigma_b^2)
#'         - \frac12\log\{ \sum_j p_{ij}(1-p_{ij}) + 1/\sigma_b^2 \}
#'       \right],
#'   }
#'
#'   where \eqn{\eta_{ij}=x_{ij}^\top\beta+b_i} and each \eqn{b_i} is the
#'   protected Newton mode for cluster \eqn{i}.
#'
#'   Privacy: per-patient linear predictors, probabilities, residuals,
#'   working responses, scores, and BLUPs never leave the share-domain
#'   servers. The client only orchestrates protected sigmoid, softplus,
#'   Beaver products, scalar likelihood pieces, and guarded per-cluster
#'   score/curvature sums. Per-cluster working modes are internal optimiser
#'   state and are deliberately not returned.
#'
#' @param formula Fixed-effects formula (0/1 binomial outcome on LHS).
#' @param data Aligned data-frame name.
#' @param cluster_col Cluster id column on the outcome server.
#' @param start Optional named vector of fixed-effect starting values.
#' @param max_outer Total protected objective-evaluation budget across
#'   random-effect standard-deviation starts. \code{NULL} chooses an adaptive
#'   budget from the number of fixed effects and starts.
#' @param mode_max_iter Newton updates for each cluster mode. The default
#'   \code{1} is the adaptive one-step Laplace route: it uses guarded
#'   score/curvature sums once per objective evaluation and is the practical
#'   explicit Laplace setting. Larger values move toward exact Laplace cluster
#'   modes at substantially higher protected-computation cost.
#' @param tol Relative optimiser tolerance.
#' @param mode_tol Cluster-mode Newton step tolerance.
#' @param lambda L2 penalty for the initial fixed-effect GLM prime.
#' @param ring Ring id. The Laplace route requires Ring127.
#' @param sigma_sd_starts Candidate random-intercept standard-deviation
#'   starts. \code{NULL} estimates one protected moment start from guarded
#'   cluster score/curvature sums. Passing a numeric vector enables explicit
#'   multi-start optimisation when extra runtime is acceptable.
#' @param sigma_sd_bounds Lower/upper bounds for the random-intercept
#'   standard deviation.
#' @param prime_iter PIRLS iterations for the initial ds.vertGLM prime.
#' @param compute_se Reserved for future protected Hessian support.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return \code{ds.vertGLMMLaplace} object with fixed effects, scalar variance
#'   component, protected objective diagnostics, and disclosure metadata.
#' @keywords internal
#' @noRd
ds.vertGLMMLaplace <- function(formula, data = NULL, cluster_col,
                               start = NULL,
                               max_outer = NULL, mode_max_iter = 1L,
                               tol = 1e-4, mode_tol = 1e-5,
                               lambda = 0, ring = NULL,
                               sigma_sd_starts = NULL,
                               sigma_sd_bounds = c(1e-4, 20),
                               prime_iter = 50L,
                               compute_se = FALSE,
                               verbose = TRUE,
                               datasources = NULL) {
  if (is.null(ring)) ring <- 127L
  ring <- as.integer(ring)
  if (ring != 127L) {
    stop("ds.vertGLMMLaplace Laplace currently requires ring=127", call. = FALSE)
  }
  if (isTRUE(compute_se)) {
    stop("Protected Hessian standard errors for ds.vertGLMMLaplace are not ",
         "implemented yet; use compute_se=FALSE.", call. = FALSE)
  }
  max_outer <- if (is.null(max_outer)) NA_integer_ else as.integer(max_outer)
  mode_max_iter <- as.integer(mode_max_iter)
  prime_iter <- as.integer(prime_iter)
  if (!is.na(max_outer) && max_outer < 1L) {
    stop("max_outer must be >= 1", call. = FALSE)
  }
  if (mode_max_iter < 1L) {
    stop("mode_max_iter must be >= 1", call. = FALSE)
  }
  if (!is.numeric(sigma_sd_bounds) || length(sigma_sd_bounds) != 2L ||
      any(!is.finite(sigma_sd_bounds)) || any(sigma_sd_bounds <= 0) ||
      sigma_sd_bounds[[1L]] >= sigma_sd_bounds[[2L]]) {
    stop("sigma_sd_bounds must be increasing positive finite values",
         call. = FALSE)
  }
  auto_sigma_starts <- is.null(sigma_sd_starts) ||
    identical(sigma_sd_starts, "auto")
  if (!auto_sigma_starts) {
    sigma_sd_starts <- as.numeric(sigma_sd_starts)
    sigma_sd_starts <- sigma_sd_starts[is.finite(sigma_sd_starts) &
                                         sigma_sd_starts > 0]
    sigma_sd_starts <- pmin(pmax(sigma_sd_starts, sigma_sd_bounds[[1L]]),
                            sigma_sd_bounds[[2L]])
    sigma_sd_starts <- unique(sigma_sd_starts)
    if (!length(sigma_sd_starts)) {
      sigma_sd_starts <- exp(mean(log(sigma_sd_bounds)))
    }
  }
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)
  y_var <- .ds_gee_extract_lhs(formula)
  y_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                        data, y_var)
  clust_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                            data, cluster_col)
  if (is.null(y_srv) || clust_srv != y_srv) {
    stop("cluster_col must live on the outcome server", call. = FALSE)
  }

  if (verbose) {
    message("[ds.vertGLMMLaplace] prime: binomial ds.vertGLM")
  }
  fit <- ds.vertGLM(formula = formula, data = data, family = "binomial",
                    max_iter = prime_iter, tol = tol, lambda = lambda,
                    ring = ring, compute_se = FALSE,
                    compute_deviance = FALSE,
                    verbose = isTRUE(verbose), datasources = datasources,
                    keep_session = TRUE)
  on.exit(.ds_glmm_cleanup_fit_session(fit, datasources), add = TRUE)

  feature_order <- .ds_glmm_feature_order(fit)
  coef_names <- c("(Intercept)", feature_order)
  beta0 <- fit$coefficients[coef_names]
  if (!is.null(start)) {
    if (is.null(names(start))) {
      stop("start must be a named fixed-effect vector", call. = FALSE)
    }
    missing_start <- setdiff(coef_names, names(start))
    if (length(missing_start)) {
      stop("start missing coefficient(s): ",
           paste(missing_start, collapse = ", "), call. = FALSE)
    }
    beta0 <- as.numeric(start[coef_names])
    names(beta0) <- coef_names
  }

  clust_info <- DSI::datashield.aggregate(
    datasources[which(server_names == y_srv)],
    call(name = "dsvertClusterSizesDS", data_name = data,
         cluster_col = cluster_col))
  if (is.list(clust_info) && length(clust_info) == 1L) {
    clust_info <- clust_info[[1L]]
  }
  n_clusters <- as.integer(clust_info$n_clusters %||%
                             length(clust_info$sizes))
  n_obs <- as.integer(clust_info$n_total %||% sum(clust_info$sizes))
  if (n_clusters < 2L) {
    stop("ds.vertGLMMLaplace needs at least two protected clusters", call. = FALSE)
  }
  if (auto_sigma_starts) {
    sigma_sd_starts <- .ds_glmm_laplace_auto_sd_starts(
      beta_orig = beta0, fit = fit, data = data, cluster_col = cluster_col,
      datasources = datasources, server_names = server_names,
      y_srv = y_srv, n_clusters = n_clusters,
      sigma_sd_bounds = sigma_sd_bounds, verbose = verbose)
  }
  initial_moments <- attr(sigma_sd_starts, "initial_moments", exact = TRUE)

  state <- new.env(parent = emptyenv())
  state$evals <- 0L
  state$history <- numeric(0)
  zero_b <- rep(0, n_clusters)
  state$cache <- new.env(parent = emptyenv())
  state$best_eval <- NULL

  objective <- function(par) {
    p <- length(coef_names)
    beta <- par[seq_len(p)]
    names(beta) <- coef_names
    log_sd <- par[[p + 1L]]
    key <- paste(format(signif(par, 11), scientific = TRUE), collapse = "|")
    cached <- get0(key, envir = state$cache, inherits = FALSE)
    if (!is.null(cached)) return(cached$value)
    initial_par <- c(beta0, log(sigma_sd_starts[[1L]]))
    use_initial_moments <- !is.null(initial_moments) &&
      length(sigma_sd_starts) == 1L &&
      length(initial_par) == length(par) &&
      max(abs(as.numeric(par) - as.numeric(initial_par))) < 1e-10

    val <- tryCatch(
      .ds_glmm_laplace_eval(
        beta_orig = beta,
        log_sd = log_sd,
        fit = fit, data = data, cluster_col = cluster_col,
        datasources = datasources, server_names = server_names,
        y_srv = y_srv, b_start = zero_b,
        mode_max_iter = mode_max_iter, mode_tol = mode_tol,
        initial_moments = if (use_initial_moments) initial_moments else NULL,
        verbose = FALSE),
      error = function(e) {
        list(value = 1e100, error = conditionMessage(e), b = zero_b,
             mode_max_step = Inf, mode_converged = FALSE)
      })
    state$evals <- state$evals + 1L
    state$history <- c(state$history, val$value)
    if (is.finite(val$value) && val$value < 1e90) {
      if (is.null(state$best_eval) || val$value < state$best_eval$value) {
        state$best_eval <- val
      }
    }
    assign(key, val, envir = state$cache)
    if (verbose) {
      msg <- sprintf(
        "[ds.vertGLMMLaplace] eval %d: nll=%.6g sigma_b=%.4g mode_step=%.3g",
        state$evals, val$value, exp(log_sd), val$mode_max_step)
      message(msg)
    }
    val$value
  }

  p <- length(coef_names)
  if (is.na(max_outer)) {
    max_outer <- max(6L, p + 1L)
  }
  eval_budget <- max(2L, ceiling(max_outer / length(sigma_sd_starts)))
  lower <- c(rep(-Inf, p), log(sigma_sd_bounds[[1L]]))
  upper <- c(rep(Inf, p), log(sigma_sd_bounds[[2L]]))
  fits <- vector("list", length(sigma_sd_starts))
  for (ii in seq_along(sigma_sd_starts)) {
    if (verbose) {
      message(sprintf("[ds.vertGLMMLaplace] Laplace start %d/%d: sd=%.4g",
                      ii, length(sigma_sd_starts), sigma_sd_starts[[ii]]))
    }
    par0 <- c(beta0, log(sigma_sd_starts[[ii]]))
    if (max_outer <= 1L) {
      fits[[ii]] <- list(
        par = par0,
        objective = objective(par0),
        convergence = 0L,
        message = "single-pass protected objective evaluation",
        evaluations = stats::setNames(c(1L, 0L), c("function", "gradient")))
    } else {
      fits[[ii]] <- stats::nlminb(
        start = par0, objective = objective,
        lower = lower, upper = upper,
        control = list(iter.max = eval_budget, eval.max = eval_budget,
                       rel.tol = tol, x.tol = tol))
    }
  }
  obj <- vapply(fits, function(x) as.numeric(x$objective), numeric(1L))
  best <- fits[[which.min(obj)]]
  beta_hat <- best$par[seq_len(p)]
  names(beta_hat) <- coef_names
  sigma_b2 <- exp(2 * best$par[[p + 1L]])

  best_key <- paste(format(signif(best$par, 11), scientific = TRUE),
                    collapse = "|")
  final_eval <- get0(best_key, envir = state$cache, inherits = FALSE)
  if (is.null(final_eval) || !is.finite(final_eval$value)) {
    final_eval <- .ds_glmm_laplace_eval(
      beta_orig = beta_hat,
      log_sd = best$par[[p + 1L]],
      fit = fit, data = data, cluster_col = cluster_col,
      datasources = datasources, server_names = server_names,
      y_srv = y_srv, b_start = zero_b,
      mode_max_iter = mode_max_iter, mode_tol = mode_tol,
      verbose = FALSE)
  }

  finite_hist <- state$history[is.finite(state$history) &
                                 state$history < 1e90]
  tail_hist <- utils::tail(finite_hist, min(6L, length(finite_hist)))
  tail_rel_range <- NA_real_
  if (length(tail_hist) >= 4L) {
    tail_rel_range <- (max(tail_hist) - min(tail_hist)) /
      max(1, abs(min(tail_hist)))
  }
  objective_floor_hit <- is.finite(tail_rel_range) &&
    tail_rel_range <= max(tol, 1e-6)
  optimiser_ok <- best$convergence == 0L || objective_floor_hit
  eval_ok <- is.null(final_eval$error) &&
    is.finite(final_eval$value) && final_eval$value < 1e90
  mode_ok <- mode_max_iter == 1L && is.finite(final_eval$mode_max_step) ||
    isTRUE(final_eval$mode_converged) ||
    is.finite(final_eval$mode_max_step) &&
      final_eval$mode_max_step <= max(10 * mode_tol, 1e-4)
  convergence_metric <- if (mode_max_iter == 1L) 0 else final_eval$mode_max_step
  if (is.finite(tail_rel_range)) {
    convergence_metric <- if (mode_max_iter == 1L) {
      tail_rel_range
    } else {
      max(tail_rel_range, final_eval$mode_max_step, na.rm = TRUE)
    }
  }
  quality <- .dsvert_quality_from_convergence(
    converged = eval_ok && optimiser_ok && mode_ok,
    metric = convergence_metric,
    tolerance = max(tol, mode_tol),
    label = "Laplace GLMM")
  quality$metrics <- c(quality$metrics, list(
    objective = final_eval$value,
    optimizer_convergence = best$convergence,
    optimizer_message = best$message,
    function_evals = as.integer(state$evals),
    tail_objective_relative_range = tail_rel_range,
    mode_strategy = if (mode_max_iter == 1L) "one_step" else "iterated",
    mode_max_step = final_eval$mode_max_step,
    mode_iterations = final_eval$mode_iterations,
    sigma_b2 = sigma_b2))
  if (best$convergence != 0L && objective_floor_hit) {
    quality$warnings <- setdiff(
      quality$warnings,
      "Laplace GLMM did not meet its convergence criterion.")
  }
  if (!eval_ok) {
    quality$status <- "failed"
    quality$warnings <- unique(c(
      quality$warnings,
      "Laplace GLMM protected objective evaluation failed closed."))
  }

  icc <- sigma_b2 / (sigma_b2 + pi^2 / 3)
  out <- list(
    coefficients = beta_hat,
    std_errors = stats::setNames(rep(NA_real_, length(beta_hat)),
                                 names(beta_hat)),
    sigma_b2 = sigma_b2,
    sigma_b = sqrt(sigma_b2),
    icc = icc,
    n_clusters = n_clusters,
    n_obs = n_obs,
    converged = eval_ok && optimiser_ok && mode_ok,
    optimizer_convergence = best$convergence,
    optimizer_message = best$message,
    objective = final_eval$value,
    evaluations = as.integer(state$evals),
    quality = quality,
    method = if (mode_max_iter == 1L) {
      "laplace_one_step_random_intercept"
    } else {
      "laplace_random_intercept"
    },
    family = "binomial (Laplace random-intercept GLMM)",
    disclosure = list(
      patient_level_returned = FALSE,
      random_effects_returned = FALSE,
      cluster_labels_returned = FALSE,
      cluster_vectors_returned = FALSE,
      opened = c("scalar log-likelihood pieces",
                 "guarded per-cluster score sums",
                 "guarded per-cluster curvature sums")),
    call = match.call())
  class(out) <- c("ds.vertGLMMLaplace", "list")
  out
}

#' @keywords internal
print.ds.vertGLMMLaplace <- function(x, ...) {
  cat("dsVert binomial GLMM-Laplace (random intercept)\n")
  cat(sprintf("  Clusters = %d    N = %d\n", x$n_clusters, x$n_obs))
  cat(sprintf("  sigma_b^2 = %.4g    ICC (latent) = %.3f\n",
              x$sigma_b2, x$icc))
  cat(sprintf("  Converged: %s    evaluations = %d\n",
              x$converged, x$evaluations))
  if (!is.null(x$quality$status)) {
    cat(sprintf("  Quality: %s\n", x$quality$status))
    if (length(x$quality$warnings)) {
      for (w in x$quality$warnings) cat("  - ", w, "\n", sep = "")
    }
  }
  cat("\nFixed effects (log-odds):\n")
  z <- x$coefficients / x$std_errors
  print(round(data.frame(
    Estimate = x$coefficients, SE = x$std_errors,
    z = z, p = 2 * stats::pnorm(-abs(z)),
    check.names = FALSE), 5L))
  invisible(x)
}

#' @keywords internal
.ds_glmm_laplace_eval <- function(beta_orig, log_sd, fit, data, cluster_col,
                                   datasources, server_names, y_srv, b_start,
                                   mode_max_iter, mode_tol,
                                   initial_moments = NULL,
                                   verbose = FALSE) {
  if (!is.finite(log_sd)) {
    return(list(value = 1e100, b = b_start, h = rep(Inf, length(b_start)),
                mode_converged = FALSE, mode_iterations = 0L,
                mode_max_step = Inf))
  }
  sigma2 <- exp(2 * log_sd)
  if (!is.finite(sigma2) || sigma2 <= 0 || sigma2 > 1e12) {
    return(list(value = 1e100, b = b_start, h = rep(Inf, length(b_start)),
                mode_converged = FALSE, mode_iterations = 0L,
                mode_max_step = Inf))
  }
  b <- as.numeric(b_start)
  if (!length(b) || any(!is.finite(b))) b <- rep(0, length(b_start))
  offset_col <- "__dsvert_glmm_laplace_b"
  fit$coefficients <- beta_orig
  h <- rep(Inf, length(b))
  mode_converged <- FALSE
  mode_max_step <- Inf
  last_comp <- NULL

  for (iter in seq_len(mode_max_iter)) {
    if (iter == 1L && !is.null(initial_moments) &&
        length(initial_moments$rsum) == length(b) &&
        length(initial_moments$vsum) == length(b) &&
        max(abs(b)) < 1e-12) {
      comp <- NULL
      rsum <- as.numeric(initial_moments$rsum)
      s <- pmax(as.numeric(initial_moments$vsum), 1e-12)
    } else {
      .ds_glmm_laplace_set_cluster_offset(
        data = data, cluster_col = cluster_col, weights = b,
        output_column = offset_col, datasources = datasources,
        server_names = server_names, y_srv = y_srv,
        session_id = fit$session_id)
      comp <- .ds_glmm_share_domain_moments(
        fit = fit, data = data, cluster_col = cluster_col,
        datasources = datasources, server_names = server_names,
        verbose = FALSE, laplace_components = FALSE)
      rsum <- as.numeric(comp$rsum_per_cluster)
      s <- pmax(as.numeric(comp$vsum_per_cluster), 1e-12)
    }
    if (length(rsum) != length(b) || length(s) != length(b)) {
      stop("protected GLMM-Laplace cluster aggregate length changed",
           call. = FALSE)
    }
    h <- s + 1 / sigma2
    score <- rsum - b / sigma2
    step <- pmax(pmin(score / h, 5), -5)
    step[!is.finite(step)] <- 0
    b <- b + step
    mode_max_step <- max(abs(step))
    last_comp <- comp
    if (verbose) {
      message(sprintf("[ds.vertGLMMLaplace] mode %d max_step=%.3g",
                      iter, mode_max_step))
    }
    if (is.finite(mode_max_step) && mode_max_step <= mode_tol) {
      mode_converged <- TRUE
      break
    }
  }

  .ds_glmm_laplace_set_cluster_offset(
    data = data, cluster_col = cluster_col, weights = b,
    output_column = offset_col, datasources = datasources,
    server_names = server_names, y_srv = y_srv,
    session_id = fit$session_id)
  final_comp <- .ds_glmm_share_domain_moments(
    fit = fit, data = data, cluster_col = cluster_col,
    datasources = datasources, server_names = server_names,
    verbose = FALSE, laplace_components = TRUE)
  lap <- final_comp$laplace_components
  s <- pmax(as.numeric(lap$s_by_cluster), 1e-12)
  h <- s + 1 / sigma2
  value <- -(as.numeric(lap$y_dot_eta) - as.numeric(lap$sum_softplus) -
      sum(b * b) / (2 * sigma2) -
      0.5 * length(b) * log(sigma2) -
      0.5 * sum(log(pmax(h, 1e-12))))
  if (!is.finite(value)) value <- 1e100
  list(value = value, b = b, h = h,
       mode_converged = mode_converged,
       mode_iterations = iter,
       mode_max_step = mode_max_step,
       components = NULL,
       last_components = is.null(last_comp))
}

#' @keywords internal
.ds_glmm_laplace_auto_sd_starts <- function(beta_orig, fit, data, cluster_col,
                                           datasources, server_names, y_srv,
                                           n_clusters, sigma_sd_bounds,
                                           verbose = FALSE) {
  offset_col <- "__dsvert_glmm_laplace_b"
  .ds_glmm_laplace_set_cluster_offset(
    data = data, cluster_col = cluster_col,
    weights = rep(0, n_clusters), output_column = offset_col,
    datasources = datasources, server_names = server_names, y_srv = y_srv,
    session_id = fit$session_id)
  fit$coefficients <- beta_orig
  comp <- .ds_glmm_share_domain_moments(
    fit = fit, data = data, cluster_col = cluster_col,
    datasources = datasources, server_names = server_names,
    verbose = FALSE, laplace_components = FALSE)
  rsum <- as.numeric(comp$rsum_per_cluster)
  s <- pmax(as.numeric(comp$vsum_per_cluster), 1e-8)
  start_prior_sd <- as.numeric(getOption(
    "dsvert.glmm_laplace_start_prior_sd", 0.5))
  if (!is.finite(start_prior_sd) || start_prior_sd <= 0) {
    start_prior_sd <- 0.5
  }
  start_precision <- 1 / (start_prior_sd * start_prior_sd)
  b_proxy <- pmax(pmin(rsum / (s + start_precision), 5), -5)
  sd0 <- stats::sd(b_proxy)
  if (!is.finite(sd0) || sd0 <= 0) {
    sd0 <- exp(mean(log(sigma_sd_bounds)))
  }
  sd0 <- pmin(pmax(sd0, sigma_sd_bounds[[1L]]), sigma_sd_bounds[[2L]])
  starts <- sd0
  starts <- pmin(pmax(starts, sigma_sd_bounds[[1L]]), sigma_sd_bounds[[2L]])
  starts <- unique(signif(starts[is.finite(starts) & starts > 0], 6L))
  if (!length(starts)) starts <- exp(mean(log(sigma_sd_bounds)))
  attr(starts, "initial_moments") <- list(rsum = rsum, vsum = s)
  if (verbose) {
    message("[ds.vertGLMMLaplace] adaptive sd starts: ",
            paste(format(starts, digits = 4), collapse = ", "))
  }
  starts
}

#' @keywords internal
.ds_glmm_laplace_set_cluster_offset <- function(data, cluster_col, weights,
                                                output_column, datasources,
                                                server_names, y_srv,
                                                session_id) {
  y_ci <- which(server_names == y_srv)
  tryCatch(
    DSI::datashield.aggregate(
      datasources[y_ci],
      call(name = "dsvertExpandClusterWeightsDS",
           data_name = data, cluster_col = cluster_col,
           weights_per_cluster = as.numeric(weights),
           output_column = output_column)),
    error = function(e) {
      stop("[GLMM-Laplace] offset expand failed: ",
           conditionMessage(e), call. = FALSE)
    })
  tryCatch(
    DSI::datashield.aggregate(
      datasources[y_ci],
      call(name = "k2SetOffsetDS", data_name = data,
           offset_column = output_column, session_id = session_id)),
    error = function(e) {
      stop("[GLMM-Laplace] offset registration failed: ",
           conditionMessage(e), call. = FALSE)
    })
  invisible(TRUE)
}
