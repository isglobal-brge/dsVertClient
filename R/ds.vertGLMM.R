#' @title Federated binomial GLMM via Laplace approximation
#' @description Fit a binomial generalised linear mixed model with a
#'   single random intercept on vertically partitioned DataSHIELD data
#'   via the Laplace-approximated marginal likelihood:
#'
#'   \deqn{\ell_L(\beta, \sigma_b^2) = \sum_i \left[\ell_i(\beta, b_i^*)
#'         - \tfrac{1}{2} \log |H_i(b_i^*)| \right]}
#'
#'   where \eqn{b_i^*} is the mode of the per-cluster penalised
#'   log-likelihood \eqn{\ell_i(\beta, b) = \sum_j \log f(y_{ij}|
#'   x_{ij}^T \beta + b) - b^2/(2\sigma_b^2)} and \eqn{H_i} is its
#'   Hessian at the mode.
#'
#'   Architecture:
#'   \enumerate{
#'     \item For each cluster (cluster IDs on outcome server), binomial
#'           score and information terms \eqn{\sum(y-p)} and
#'           \eqn{\sum p(1-p)} are computed from DCF shares. The outcome
#'           server broadcasts cluster membership only to the DCF peer; the
#'           client receives aggregate cluster sums only.
#'     \item The outer optimiser on \eqn{(\beta, \sigma_b^2)} is
#'           client-side, re-calling ds.vertGLM with the shrinkage-
#'           weight column derived from the current \eqn{b_i^*}.
#'     \item Variance-component update uses the moment-matching
#'           estimate \eqn{\hat\sigma_b^2 = \mathrm{var}(\hat b_i)}
#'           across clusters (Breslow-Clayton approximation).
#'   }
#'
#'   The scaffolding here reuses every primitive already shipped in
#'   dsVert 1.1.0+: the keep_session flag on ds.vertGLM, cluster
#'   membership broadcast, per-cluster share sums, Beaver vecmul, and the
#'   DCF sigmoid wide-spline.
#'
#'   Privacy: the returned object contains fixed effects, variance-component
#'   estimates, scalar diagnostics, and the final inner \code{ds.vertGLM} fit.
#'   Per-patient quantities never leave the DCF parties. Guarded per-cluster
#'   sufficient statistics are used internally for the random-intercept update,
#'   but per-cluster BLUPs and cluster-size vectors are not returned. Original
#'   cluster labels are not returned and clusters below
#'   \code{datashield.privacyLevel} fail closed.
#'
#'   Inter-server leakage: cluster membership (same tier as ds.vertLMM).
#'
#' @param formula Fixed-effects formula (binomial outcome on LHS).
#' @param data Aligned data-frame name.
#' @param cluster_col Cluster id column on the outcome server.
#' @param max_outer Outer (beta, sigma_b^2) iterations.
#' @param inner_iter Inner PIRLS iterations per cluster per outer step.
#' @param tol Outer convergence tolerance.
#' @param lambda L2 penalty passed to the inner binomial GLM fits.
#' @param method Character. \code{"pql_aggregate"} is the default accuracy
#'   route and runs the aggregate PQL weighted-LMM update using Ring127
#'   share-domain working statistics. \code{"em"} keeps the legacy
#'   offset-GLM plus BLUP/EM approximation.
#' @param use_pearson_cap Logical. If TRUE, cap the EM variance-component
#'   update by the first marginal Pearson cluster-moment estimate. This
#'   prevents BLUP posterior-variance inflation when the random-effect signal
#'   is weak or singular, using only the same guarded cluster aggregates as
#'   the default GLMM loop.
#' @param compute_se Logical. Compute GLM finite-difference standard errors
#'   for the inner fits. Set FALSE for coefficient/variance validation runs.
#' @param ring Integer 63 or 127. The PQL aggregate route requires Ring127.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return \code{ds.vertGLMM} object: fixed-effect coefficients,
#'   random-effect variance \eqn{\hat\sigma_b^2}, scalar fit diagnostics, and
#'   the converged binomial \code{fit}. Per-cluster BLUP vectors are internal
#'   working quantities and are not returned.
#' @export
ds.vertGLMM <- function(formula, data = NULL, cluster_col,
                        max_outer = 10L, inner_iter = 10L,
                        tol = 1e-3, lambda = 0,
                        method = c("pql_aggregate", "em"),
                        use_pearson_cap = getOption(
                          "dsvert.glmm_use_pearson_cap", TRUE),
                        compute_se = TRUE,
                        ring = NULL,
                        verbose = TRUE,
                        datasources = NULL) {
  method <- match.arg(method)
  if (is.null(ring)) ring <- if (method == "pql_aggregate") 127L else 63L
  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  if (method == "pql_aggregate" && ring != 127L) {
    stop("method='pql_aggregate' currently requires ring=127", call. = FALSE)
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

  # Prime with a straight binomial fit ignoring random effects.
  if (verbose) message("[ds.vertGLMM] prime: binomial ds.vertGLM")
  fit <- ds.vertGLM(formula = formula, data = data, family = "binomial",
                    max_iter = inner_iter, tol = tol, lambda = lambda,
                    ring = ring,
                    compute_se = isTRUE(compute_se),
                    verbose = FALSE,
                    datasources = datasources,
                    keep_session = TRUE)

  # Cluster sizes + per-cluster residuals (all aggregates).
  clust_info <- DSI::datashield.aggregate(
    datasources[which(server_names == y_srv)],
    call(name = "dsvertClusterSizesDS", data_name = data,
         cluster_col = cluster_col))
  if (is.list(clust_info) && length(clust_info) == 1L)
    clust_info <- clust_info[[1L]]
  n_i <- as.integer(clust_info$sizes)
  n_clusters <- length(n_i)

  if (method == "pql_aggregate") {
    return(.ds_glmm_pql_aggregate_loop(
      fit = fit, data = data, cluster_col = cluster_col,
      n_i = n_i, n_clusters = n_clusters, max_outer = max_outer,
      tol = tol, compute_se = compute_se,
      datasources = datasources, server_names = server_names,
      y_srv = y_srv, verbose = verbose, call = match.call()))
  }

  # Outer optimiser: alternate (i) compute cluster BLUPs \hat b_i via
  # Laplace-approximated inner PIRLS on y vs. binomial likelihood, and
  # (ii) refit the fixed effects with the outer-loop intercept offsets
  # b_i. Moment-match \sigma_b^2 = var(\hat b_i).
  sigma_b2 <- 1.0
  sigma_b2_anchor <- NA_real_
  sigma_b2_pearson_cap <- NA_real_
  b_hat <- rep(0, n_clusters)
  converged <- FALSE
  offset_col <- NULL
  for (outer in seq_len(max_outer)) {
    # Per-cluster binomial score aggregates at current fixed effects and
    # current random-intercept offset. The vertical path keeps p, y-p and
    # p(1-p) in the DCF share domain and releases only per-cluster sums.
    cl <- .ds_glmm_share_domain_moments(
      fit = fit, data = data, cluster_col = cluster_col,
      datasources = datasources, server_names = server_names,
      verbose = verbose)
    rsum <- as.numeric(cl$rsum_per_cluster)
    vsum <- as.numeric(cl$vsum_per_cluster)
    if (length(cl$n_per_cluster) == length(n_i)) {
      n_i <- as.integer(cl$n_per_cluster)
    }
    active <- n_i > 0L & is.finite(rsum) & is.finite(vsum)
    if (isTRUE(use_pearson_cap) && !is.finite(sigma_b2_pearson_cap) &&
        sum(active) >= 2L) {
      inv_v <- 1 / pmax(vsum[active], 1e-8)
      pearson_cluster <- rsum[active] * inv_v
      sigma_b2_pearson_cap <- max(
        stats::var(pearson_cluster) - mean(inv_v),
        1e-6)
    }
    score_mom_num <- sum(rsum[active]^2 - vsum[active])
    score_mom_den <- sum(vsum[active]^2)
    sigma_b2_score <- if (is.finite(score_mom_den) && score_mom_den > 0) {
      max(score_mom_num / score_mom_den, 1e-6)
    } else {
      1e-6
    }
    info <- rep(Inf, n_clusters)
    info[active] <- pmax(vsum[active], 1e-8) + 1 / sigma_b2
    score <- rsum - b_hat / sigma_b2
    step <- score / info
    step[!is.finite(step)] <- 0
    step <- pmax(pmin(step, 3), -3)
    b_hat_new <- b_hat + step
    b_hat_new[!active] <- 0
    if (any(active)) {
      # Random intercepts are identified with mean zero; otherwise their
      # common shift is confounded with the fixed intercept and the offset
      # refit can run away.
      b_hat_new[active] <- b_hat_new[active] - mean(b_hat_new[active])
    }
    # EM update for sigma_b^2 (task #99 fix 2026-04-21):
    # `var(b_hat)` is biased DOWN because b_hat is a shrunk BLUP:
    # E[b_hat_i^2] = shrinkage_i * sigma_b^2 + shrinkage_i^2 * posterior_var,
    # so var(b_hat) systematically underestimates sigma_b^2 by the
    # shrinkage factor -- on a 15x20 synth with sigma_b = 0.7 the
    # previous rule collapsed sigma_b^2 -> 0.001 instead of 0.49.
    #
    # The canonical Laird-Ware / Lindstrom-Bates EM update uses
    #   sigma_b^2_new = mean(b_hat_i^2 + posterior_var_i)
    # which is exactly the conditional-expectation of b_i^2 | y
    # under the current sigma_b^2. It is positive, non-decreasing from
    # var(b_hat), and is the fixed point of the EM iteration that
    # converges to the ML estimator.
    post_var <- 1 / info
    sigma_b2_em <- max(mean(b_hat_new[active]^2 + post_var[active]), 1e-6)
    use_score_anchor <- isTRUE(
      getOption("dsvert.glmm_use_score_anchor", FALSE)
    )
    if (use_score_anchor && !is.finite(sigma_b2_anchor)) {
      # The first marginal score moment prevents the EM/PQL variance update
      # from collapsing after BLUP offsets absorb cluster signal. Cap it
      # conservatively because fixed-point score moments are noisy in very
      # small clustered fixtures. It is opt-in until validated on a broader
      # grid because low-variance fixtures can otherwise be over-inflated.
      anchor_cap <- getOption("dsvert.glmm_sigma_anchor_cap",
                              2 * sigma_b2_em)
      anchor_cap <- suppressWarnings(as.numeric(anchor_cap)[1L])
      if (!is.finite(anchor_cap) || anchor_cap <= 0) anchor_cap <- Inf
      sigma_b2_anchor <- min(sigma_b2_score, anchor_cap)
    }
    sigma_b2_floor <- if (is.finite(sigma_b2_anchor)) sigma_b2_anchor else 1e-6
    sigma_b2_new <- max(sigma_b2_em, sigma_b2_floor)
    if (isTRUE(use_pearson_cap) && is.finite(sigma_b2_pearson_cap)) {
      sigma_b2_new <- max(sigma_b2_floor,
                          min(sigma_b2_new, sigma_b2_pearson_cap))
    }
    if (verbose) {
      message(sprintf(
        paste0("[GLMM] outer %d  sigma_b^2=%.4g  var(b_hat)=%.4g  ",
               "mean(post_var)=%.4g  score_anchor=%.4g  pearson_cap=%.4g"),
        outer, sigma_b2_new, stats::var(b_hat_new[active]),
        mean(post_var[active]), sigma_b2_anchor, sigma_b2_pearson_cap))
    }
    b_delta <- max(abs(b_hat_new - b_hat))
    sigma_delta <- abs(sigma_b2_new - sigma_b2)
    old_coef <- fit$coefficients
    b_hat <- b_hat_new
    sigma_b2 <- sigma_b2_new

    # Expand b_i into a per-patient offset column so the next
    # ds.vertGLM fit adds it through the existing offset= plumbing.
    tryCatch(
      DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call(name = "dsvertExpandClusterWeightsDS",
             data_name = data, cluster_col = cluster_col,
             weights_per_cluster = as.numeric(b_hat),
             output_column = "__dsvert_glmm_b")),
      error = function(e) {
        message("[GLMM] offset expand failed: ", conditionMessage(e))
      })
    offset_col <- "__dsvert_glmm_b"
    .ds_glmm_cleanup_fit_session(fit, datasources)
    new_fit <- tryCatch(
      ds.vertGLM(formula = formula, data = data, family = "binomial",
                 offset = offset_col,
                 max_iter = inner_iter, tol = tol, lambda = lambda,
                 ring = ring,
                 compute_se = isTRUE(compute_se),
                 verbose = FALSE,
                 datasources = datasources,
                 keep_session = TRUE),
      error = function(e) {
        message("[GLMM] inner binomial refit failed: ",
                conditionMessage(e)); NULL })
    if (is.null(new_fit)) break
    fit <- new_fit
    common_coef <- intersect(names(old_coef), names(fit$coefficients))
    beta_delta <- if (length(common_coef)) {
      max(abs(old_coef[common_coef] - fit$coefficients[common_coef]))
    } else Inf
    if (max(sigma_delta, b_delta, beta_delta) < tol * max(1, sigma_b2)) {
      converged <- TRUE
      break
    }
  }
  .ds_glmm_cleanup_fit_session(fit, datasources)

  icc <- sigma_b2 / (sigma_b2 + pi^2 / 3)  # logistic latent-variance
  out <- list(
    coefficients = fit$coefficients,
    std_errors   = fit$std_errors,
    sigma_b2     = sigma_b2,
    sigma_b2_anchor = sigma_b2_anchor,
    sigma_b2_pearson_cap = sigma_b2_pearson_cap,
    icc          = icc,
    n_clusters   = n_clusters,
    n_obs        = sum(n_i),
    converged    = converged,
    iterations   = outer,
    fit          = fit,
    family       = "binomial (Laplace-approximated)",
    call         = match.call())
  class(out) <- c("ds.vertGLMM", "list")
  out
}

#' @export
print.ds.vertGLMM <- function(x, ...) {
  cat("dsVert binomial GLMM (Laplace approximation)\n")
  n_obs <- if (!is.null(x$n_obs)) x$n_obs else sum(x$cluster_sizes %||% 0L)
  cat(sprintf("  Clusters = %d    N = %d\n",
              x$n_clusters, n_obs))
  cat(sprintf("  sigma_b^2 = %.4g    ICC (latent) = %.3f\n",
              x$sigma_b2, x$icc))
  cat(sprintf("  Converged: %s (%d outer iters)\n",
              x$converged, x$iterations))
  cat("\nFixed effects (log-odds):\n")
  z <- x$coefficients / x$std_errors
  print(round(data.frame(
    Estimate = x$coefficients, SE = x$std_errors,
    z = z, p = 2 * stats::pnorm(-abs(z)),
    check.names = FALSE), 5L))
  invisible(x)
}

#' @keywords internal
.ds_glmm_cleanup_fit_session <- function(fit, datasources) {
  if (is.null(fit$session_id) || is.null(fit$server_list)) return(invisible(NULL))
  for (srv in fit$server_list) {
    ci <- which(names(datasources) == srv)
    if (length(ci) == 1L) {
      tryCatch(DSI::datashield.aggregate(datasources[ci],
        call(name = "mpcCleanupDS", session_id = fit$session_id)),
        error = function(e) NULL)
      tryCatch(DSI::datashield.aggregate(datasources[ci],
        call(name = "mpcGcDS")), error = function(e) NULL)
    }
  }
  invisible(NULL)
}

#' @keywords internal
.ds_glmm_feature_order <- function(fit) {
  x_vars <- fit$x_vars
  server_list <- fit$server_list
  coordinator <- fit$y_server
  if (identical(fit$eta_privacy, "secure_agg")) {
    fusion <- .k3_select_fusion_server(server_list, coordinator, x_vars)
    non_dcf <- setdiff(server_list, c(coordinator, fusion))
    c(x_vars[[coordinator]], x_vars[[fusion]],
      unlist(x_vars[non_dcf], use.names = FALSE))
  } else {
    nl <- setdiff(server_list, coordinator)
    c(x_vars[[coordinator]], unlist(x_vars[nl], use.names = FALSE))
  }
}

#' @keywords internal
.ds_glmm_unstandardize_beta <- function(beta_std, fit, features) {
  beta_std <- as.numeric(beta_std)
  names(beta_std) <- c("(Intercept)", features)
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])
  names(x_sds) <- names(x_means) <- features
  beta_feat <- beta_std[features] / x_sds
  intercept <- beta_std["(Intercept)"] - sum(beta_feat * x_means)
  out <- c("(Intercept)" = as.numeric(intercept), beta_feat)
  names(out) <- c("(Intercept)", features)
  out
}

#' @keywords internal
.ds_glmm_safe_solve <- function(a, b) {
  out <- tryCatch(
    solve(a, b),
    error = function(e) qr.solve(a + diag(1e-10, nrow(a)), b))
  as.numeric(out)
}

#' @keywords internal
.ds_glmm_pql_solve_components <- function(comp,
                                          tau_bounds = c(1e-12, 1e6)) {
  cols <- comp$columns
  xwx <- as.matrix(comp$xwx)
  xwz <- as.numeric(comp$xwz)
  q <- as.matrix(comp$q_by_cluster)
  zsum <- as.numeric(comp$zsum_by_cluster)
  wz2 <- as.numeric(comp$wz2_by_cluster)
  s <- as.numeric(comp$s_by_cluster)
  sizes <- as.integer(comp$cluster_sizes)
  active <- sizes > 0L & is.finite(s) & s > 1e-10 &
    is.finite(zsum) & is.finite(wz2)
  if (sum(active) < 2L) {
    stop("PQL aggregate solve needs at least two active clusters",
         call. = FALSE)
  }
  q <- q[active, , drop = FALSE]
  zsum <- zsum[active]
  wz2 <- wz2[active]
  s <- s[active]
  n <- sum(sizes[active])
  p <- length(cols)

  stats_for_tau <- function(tau) {
    denom <- 1 / tau + s
    a <- xwx
    rhs <- xwz
    for (ii in seq_along(denom)) {
      qi <- q[ii, ]
      a <- a - tcrossprod(qi) / denom[ii]
      rhs <- rhs - qi * zsum[ii] / denom[ii]
    }
    beta <- .ds_glmm_safe_solve(a, rhs)
    quad <- as.numeric(crossprod(beta, xwx %*% beta))
    r_v0_r <- sum(wz2) - 2 * sum(beta * xwz) + quad
    adj <- numeric(length(denom))
    for (ii in seq_along(denom)) {
      zi <- zsum[ii] - sum(q[ii, ] * beta)
      adj[ii] <- zi^2 / denom[ii]
    }
    r_v0_r <- max(r_v0_r - sum(adj), 1e-12)
    phi <- max(r_v0_r / n, 1e-12)
    nll <- n * log(phi) + sum(log1p(tau * s)) + r_v0_r / phi
    list(beta = beta, phi = phi, r_v0_r = r_v0_r, nll = nll)
  }

  objective <- function(log_tau) stats_for_tau(exp(log_tau))$nll
  opt <- optimize(objective, interval = log(tau_bounds))
  tau <- exp(opt$minimum)
  fit <- stats_for_tau(tau)
  beta <- fit$beta
  names(beta) <- cols
  b_hat <- numeric(length(sizes))
  denom <- 1 / tau + s
  b_active <- numeric(length(s))
  for (ii in seq_along(s)) {
    b_active[ii] <- (zsum[ii] - sum(q[ii, ] * beta)) / denom[ii]
  }
  b_hat[active] <- b_active
  list(beta_std = beta,
       b_hat = b_hat,
       sigma_b2 = fit$phi * tau,
       phi = fit$phi,
       tau = tau,
       nll = fit$nll,
       active = active,
       optimizer = opt)
}

#' @keywords internal
.ds_glmm_pql_aggregate_loop <- function(fit, data, cluster_col,
                                        n_i, n_clusters, max_outer,
                                        tol, compute_se,
                                        datasources, server_names,
                                        y_srv, verbose, call) {
  if (as.integer(fit$ring %||% 63L) != 127L) {
    stop("aggregate PQL GLMM requires a Ring127 GLM session", call. = FALSE)
  }
  feature_order <- .ds_glmm_feature_order(fit)
  beta_orig <- fit$coefficients[c("(Intercept)", feature_order)]
  b_hat <- rep(0, n_clusters)
  sigma_b2 <- 1
  phi <- 1
  tau <- sigma_b2 / phi
  converged <- FALSE
  trace <- data.frame()
  offset_col <- "__dsvert_glmm_pql_b"
  on.exit(.ds_glmm_cleanup_fit_session(fit, datasources), add = TRUE)

  for (outer in seq_len(max_outer)) {
    tryCatch(
      DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call(name = "dsvertExpandClusterWeightsDS",
             data_name = data, cluster_col = cluster_col,
             weights_per_cluster = as.numeric(b_hat),
             output_column = offset_col)),
      error = function(e) {
        stop("[GLMM PQL] offset expand failed: ",
             conditionMessage(e), call. = FALSE)
      })
    tryCatch(
      DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call(name = "k2SetOffsetDS", data_name = data,
             offset_column = offset_col, session_id = fit$session_id)),
      error = function(e) {
        stop("[GLMM PQL] offset registration failed: ",
             conditionMessage(e), call. = FALSE)
      })

    fit$coefficients <- beta_orig
    cl <- .ds_glmm_share_domain_moments(
      fit = fit, data = data, cluster_col = cluster_col,
      datasources = datasources, server_names = server_names,
      verbose = verbose, pql_components = TRUE,
      pql_feature_order = feature_order)
    if (length(cl$n_per_cluster) == length(n_i)) {
      n_i <- as.integer(cl$n_per_cluster)
    }
    step <- .ds_glmm_pql_solve_components(cl$pql_components)
    beta_new <- .ds_glmm_unstandardize_beta(step$beta_std, fit,
                                            feature_order)
    common <- intersect(names(beta_orig), names(beta_new))
    beta_delta <- max(abs(beta_new[common] - beta_orig[common]))
    b_delta <- max(abs(step$b_hat - b_hat))
    sigma_delta <- abs(step$sigma_b2 - sigma_b2)

    beta_orig <- beta_new
    b_hat <- step$b_hat
    sigma_b2 <- step$sigma_b2
    phi <- step$phi
    tau <- step$tau
    trace <- rbind(trace, data.frame(
      outer = outer, beta_delta = beta_delta, b_delta = b_delta,
      sigma_delta = sigma_delta, sigma_b2 = sigma_b2, phi = phi,
      tau = tau, nll = step$nll))
    if (verbose) {
      message(sprintf(
        "[GLMM PQL] outer %d  sigma_b^2=%.4g  phi=%.4g  tau=%.4g  beta_delta=%.3e",
        outer, sigma_b2, phi, tau, beta_delta))
    }
    if (max(beta_delta, b_delta, sigma_delta) < tol) {
      converged <- TRUE
      break
    }
  }

  se <- if (isTRUE(compute_se)) {
    stats::setNames(rep(NA_real_, length(beta_orig)), names(beta_orig))
  } else {
    stats::setNames(rep(NA_real_, length(beta_orig)), names(beta_orig))
  }
  fit$coefficients <- beta_orig
  fit$std_errors <- se
  icc <- sigma_b2 / (sigma_b2 + pi^2 / 3)
  out <- list(
    coefficients = beta_orig,
    std_errors = se,
    sigma_b2 = sigma_b2,
    sigma_b2_anchor = NA_real_,
    sigma_b2_pearson_cap = NA_real_,
    phi = phi,
    tau = tau,
    icc = icc,
    n_clusters = n_clusters,
    n_obs = sum(n_i),
    converged = converged,
    iterations = outer,
    trace = trace,
    fit = fit,
    method = "pql_aggregate",
    family = "binomial (aggregate PQL weighted-LMM)",
    call = call)
  class(out) <- c("ds.vertGLMM", "list")
  out
}

#' @keywords internal
.ds_glmm_share_domain_moments <- function(fit, data, cluster_col,
                                          datasources, server_names,
                                          verbose = FALSE,
                                          pql_components = FALSE,
                                          pql_feature_order = NULL) {
  if (!identical(fit$family, "binomial")) {
    stop("GLMM share-domain moments require a binomial ds.vertGLM fit",
         call. = FALSE)
  }
  required <- c("session_id", "transport_pks", "server_list", "x_vars",
                "y_server", "x_sds", "x_means")
  missing_req <- required[vapply(required, function(nm) is.null(fit[[nm]]),
                                 logical(1L))]
  if (length(missing_req) > 0L) {
    stop("GLMM fit session metadata missing: ",
         paste(missing_req, collapse = ", "), call. = FALSE)
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
  target_features <- unlist(x_vars[server_list], use.names = FALSE)
  std <- .ds_gee_standardized_parameters(fit, target_features)

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

  if (fit$eta_privacy == "k2_beaver") {
    nl <- setdiff(server_list, coordinator)
    if (length(nl) != 1L) {
      stop("K=2 GLMM moments require exactly one non-outcome server",
           call. = FALSE)
    }
    nl <- nl[[1L]]
    dcf_parties <- c(coordinator, nl)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    dealer_conn <- dcf_conns[[2L]]
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
    stop("unsupported eta_privacy for GLMM moments: ", fit$eta_privacy,
         call. = FALSE)
  }

  sigmoid_intervals <- as.integer(getOption(
    "dsvert.glm_num_intervals_binomial", 100L))
  if (!is.finite(sigmoid_intervals) || sigmoid_intervals < 10L) {
    sigmoid_intervals <- 100L
  }

  dcf <- .dsAgg(datasources[dealer_conn],
    call(name = "glmRing63GenDcfKeysDS",
         dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
         dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
         family = "sigmoid", n = as.integer(n_obs),
         frac_bits = frac_bits, num_intervals = sigmoid_intervals,
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
  if (is.list(spline_t) && length(spline_t) == 1L) spline_t <- spline_t[[1L]]
  .sendBlob(spline_t$spline_blob_0, "k2_spline_triples", dcf_conns[[1L]])
  .sendBlob(spline_t$spline_blob_1, "k2_spline_triples", dcf_conns[[2L]])
  for (ph in 1:4) {
    pr <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(paste0("k2WideSplinePhase", ph, "DS"),
             party_id = as.integer(i - 1L),
             family = "binomial",
             num_intervals = sigmoid_intervals,
             frac_bits = frac_bits,
             ring = ring,
             session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      pr[[i]] <- r
    }
    if (ph == 1L) {
      .sendBlob(pr[[1L]]$dcf_masked, "k2_peer_dcf_masked", dcf_conns[[2L]])
      .sendBlob(pr[[2L]]$dcf_masked, "k2_peer_dcf_masked", dcf_conns[[1L]])
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

  peer_i <- if (dcf_parties[[1L]] == coordinator) 2L else 1L
  coord_i <- if (dcf_parties[[1L]] == coordinator) 1L else 2L
  cb <- .dsAgg(datasources[dcf_conns[[coord_i]]],
    call(name = "dsvertClusterIDsBroadcastDS",
         data_name = data, cluster_col = cluster_col,
         peer_pk = transport_pks[[dcf_parties[[peer_i]]]],
         session_id = session_id))
  if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
  .sendBlob(cb$peer_blob, "dsvert_cluster_ids_blob", dcf_conns[[peer_i]])
  .dsAgg(datasources[dcf_conns[[peer_i]]],
    call(name = "dsvertClusterIDsReceiveDS", session_id = session_id))

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

  .cluster_sum <- function(share_key) {
    parts <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "dsvertPerClusterSumShareDS",
             share_key = share_key,
             session_id = session_id,
             frac_bits = frac_bits, ring = ring))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      parts[[i]] <- r
    }
    K <- length(parts[[1L]]$per_cluster_fp)
    vals <- numeric(K)
    for (kk in seq_len(K)) {
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = parts[[1L]]$per_cluster_fp[[kk]],
        share_b = parts[[2L]]$per_cluster_fp[[kk]],
        frac_bits = frac_bits, ring = ring_tag))
      vals[[kk]] <- as.numeric(agg$values[1L])
    }
    list(values = vals,
         sizes = as.integer(parts[[coord_i]]$cluster_sizes))
  }

  for (ci in dcf_conns) {
    .dsAgg(datasources[ci],
      call(name = "k2PrepareWeightedResidualShareDS",
           session_id = session_id))
  }
  py <- .cluster_sum("k2_weight_residual_share_fp")
  rsum <- -py$values

  for (i in seq_along(dcf_parties)) {
    .dsAgg(datasources[dcf_conns[[i]]],
      call(name = "dsvertGLMMOneMinusMuDS",
           output_key = "glmm_one_minus_mu_share",
           is_party0 = (i == 1L),
           session_id = session_id,
           frac_bits = frac_bits, ring = ring))
  }
  .vecmul("secure_mu_share", "glmm_one_minus_mu_share", "glmm_v_share")
  vs <- .cluster_sum("glmm_v_share")
  if (verbose) {
    message(sprintf("[GLMM] share-domain cluster moments: K=%d",
                    length(rsum)))
  }

  pql <- NULL
  if (isTRUE(pql_components)) {
    if (ring != 127L) {
      stop("GLMM aggregate PQL components require ring=127", call. = FALSE)
    }
    feature_order <- pql_feature_order
    if (is.null(feature_order)) feature_order <- .ds_glmm_feature_order(fit)
    feature_order <- as.character(feature_order)
    p_feat <- length(feature_order)
    p_full <- p_feat + 1L
    col_names <- c("(Intercept)", feature_order)
    n_int <- as.integer(n_obs)
    ring_parties <- c(coordinator, setdiff(dcf_parties, coordinator)[1L])
    ring_nl <- ring_parties[[2L]]

    .ring_vecmul <- function(x_key, y_key, output_key) {
      .ring127_vecmul(
        x_key = x_key, y_key = y_key, output_key = output_key,
        n = n_int, datasources = datasources, dealer_ci = dealer_conn,
        server_list = ring_parties, server_names = server_names,
        y_server = coordinator, nl = ring_nl,
        transport_pks = transport_pks, session_id = session_id,
        .dsAgg = .dsAgg, .sendBlob = .sendBlob)
      invisible(output_key)
    }
    .ring_affine <- function(a_key = NULL, b_key = NULL,
                             sign_a = 0L, sign_b = 0L,
                             output_key, public_const_fp = NULL) {
      for (srv in ring_parties) {
        .dsAgg(datasources[which(server_names == srv)],
          call(name = "k2Ring127AffineCombineDS",
               a_key = a_key, b_key = b_key,
               sign_a = as.numeric(sign_a), sign_b = as.numeric(sign_b),
               public_const_fp = public_const_fp,
               is_party0 = identical(srv, coordinator),
               output_key = output_key, n = as.numeric(n_int),
               session_id = session_id))
      }
      invisible(output_key)
    }
    .ring_scale <- function(in_key, scalar, output_key) {
      scalar_fp <- dsVert:::.callMpcTool("k2-float-to-fp", list(
        values = array(as.numeric(scalar), dim = 1L),
        frac_bits = 50L, ring = "ring127"))$fp_data
      scalar_fp <- .to_b64url(scalar_fp)
      for (srv in ring_parties) {
        .dsAgg(datasources[which(server_names == srv)],
          call(name = "k2Ring127LocalScaleDS",
               in_key = in_key, scalar_fp = scalar_fp,
               output_key = output_key, n = as.numeric(n_int),
               session_id = session_id,
               is_party0 = identical(srv, coordinator)))
      }
      invisible(output_key)
    }
    .sum_share <- function(share_key) {
      parts <- vector("list", 2L)
      for (i in seq_along(dcf_parties)) {
        r <- .dsAgg(datasources[dcf_conns[[i]]],
          call(name = "k2BeaverSumShareDS",
               source_key = share_key,
               session_id = session_id,
               frac_bits = as.numeric(frac_bits), ring = ring_tag))
        if (is.list(r) && length(r) == 1L) r <- r[[1L]]
        parts[[i]] <- r
      }
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = parts[[1L]]$sum_share_fp,
        share_b = parts[[2L]]$sum_share_fp,
        frac_bits = frac_bits, ring = ring_tag))
      as.numeric(agg$values[1L])
    }

    for (i in seq_along(ring_parties)) {
      .dsAgg(datasources[which(server_names == ring_parties[[i]])],
        call(name = "dsvertGEEInterceptShareDS",
             output_key = "glmm_pql_x_col_0",
             n = as.numeric(n_int),
             is_party0 = identical(ring_parties[[i]], coordinator),
             session_id = session_id,
             frac_bits = as.numeric(frac_bits), ring = as.numeric(ring)))
    }
    x_keys <- stats::setNames(
      c("glmm_pql_x_col_0", paste0("glmm_pql_x_col_", seq_len(p_feat))),
      col_names)
    if (p_feat > 0L) {
      for (jj in seq_len(p_feat)) {
        for (srv in ring_parties) {
          .dsAgg(datasources[which(server_names == srv)],
            call(name = "k2BeaverExtractColumnDS",
                 source_key = "k2_x_full_fp",
                 n = as.numeric(n_int), K = as.numeric(p_feat),
                 col_index = as.numeric(jj),
                 output_key = x_keys[[jj + 1L]],
                 session_id = session_id,
                 frac_bits = as.numeric(frac_bits), ring = ring_tag))
        }
      }
    }

    .ring_vecmul("glmm_v_share", "k2_eta_share_fp",
                 "glmm_pql_weta_share")
    .ring_affine(a_key = "glmm_pql_weta_share",
                 b_key = "k2_weight_residual_share_fp",
                 sign_a = 1L, sign_b = -1L,
                 output_key = "glmm_pql_h_share")
    .ring_vecmul("glmm_pql_h_share", "glmm_pql_h_share",
                 "glmm_pql_h2_share")
    # w lies in (0, 0.25]. Scaling by 4096 maps ordinary GLMM working
    # weights into the Ring127 reciprocal core [1, 3000] without revealing
    # any patient-level range information. The final local scale restores
    # 1 / w = 4096 * (1 / (4096 * w)).
    .ring_scale("glmm_v_share", 4096, "glmm_pql_w_scaled_share")
    .ring127_recip_round_keyed(
      in_key = "glmm_pql_w_scaled_share",
      out_key = "glmm_pql_inv_w_scaled_share",
      n = n_int, datasources = datasources, dealer_ci = dealer_conn,
      server_list = ring_parties, server_names = server_names,
      y_server = coordinator, nl = ring_nl,
      transport_pks = transport_pks, session_id = session_id,
      .dsAgg = .dsAgg, .sendBlob = .sendBlob)
    .ring_scale("glmm_pql_inv_w_scaled_share", 4096,
                "glmm_pql_inv_w_share")
    .ring_vecmul("glmm_pql_h2_share", "glmm_pql_inv_w_share",
                 "glmm_pql_wz2_share")

    zsum <- .cluster_sum("glmm_pql_h_share")$values
    wz2 <- .cluster_sum("glmm_pql_wz2_share")$values
    q_by_cluster <- matrix(0, nrow = length(vs$values), ncol = p_full,
                           dimnames = list(NULL, col_names))
    q_by_cluster[, 1L] <- vs$values
    xw_keys <- stats::setNames(rep(NA_character_, p_full), col_names)
    xw_keys[1L] <- "glmm_v_share"
    xwz <- stats::setNames(numeric(p_full), col_names)
    xwz[1L] <- sum(zsum)
    if (p_feat > 0L) {
      for (jj in seq_len(p_feat)) {
        col_j <- col_names[[jj + 1L]]
        xw_key <- paste0("glmm_pql_xw_", jj)
        xw_keys[[jj + 1L]] <- xw_key
        .ring_vecmul(x_keys[[col_j]], "glmm_v_share", xw_key)
        q_by_cluster[, jj + 1L] <- .cluster_sum(xw_key)$values
        xh_key <- paste0("glmm_pql_xh_", jj)
        .ring_vecmul(x_keys[[col_j]], "glmm_pql_h_share", xh_key)
        xwz[[jj + 1L]] <- .sum_share(xh_key)
      }
    }

    xwx <- matrix(0, nrow = p_full, ncol = p_full,
                  dimnames = list(col_names, col_names))
    xwx[1L, 1L] <- sum(vs$values)
    if (p_feat > 0L) {
      for (jj in seq_len(p_feat)) {
        xwx[1L, jj + 1L] <- sum(q_by_cluster[, jj + 1L])
        xwx[jj + 1L, 1L] <- xwx[1L, jj + 1L]
      }
      for (jj in seq_len(p_feat)) {
        for (kk in jj:p_feat) {
          key <- paste0("glmm_pql_xwx_", jj, "_", kk)
          .ring_vecmul(xw_keys[[jj + 1L]], x_keys[[kk + 1L]], key)
          val <- .sum_share(key)
          xwx[jj + 1L, kk + 1L] <- val
          xwx[kk + 1L, jj + 1L] <- val
        }
      }
    }
    pql <- list(
      columns = col_names,
      feature_names = feature_order,
      xwx = xwx,
      xwz = xwz,
      q_by_cluster = q_by_cluster,
      zsum_by_cluster = zsum,
      wz2_by_cluster = wz2,
      s_by_cluster = vs$values,
      cluster_sizes = py$sizes,
      n = sum(py$sizes))
  }

  list(rsum_per_cluster = rsum,
       vsum_per_cluster = vs$values,
       n_per_cluster = py$sizes,
       pql_components = pql)
}
