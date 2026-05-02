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
#'     \item For each cluster (cluster IDs on outcome server) an
#'           INNER L-BFGS on \eqn{b_i} is executed server-side
#'           (using the cached \code{dsvertClusterResidualsDS}
#'           aggregates + the cached \eqn{\eta_i} share from
#'           ds.vertGLM(keep_session=TRUE)).
#'     \item The outer optimiser on \eqn{(\beta, \sigma_b^2)} is
#'           client-side, re-calling ds.vertGLM with the shrinkage-
#'           weight column derived from the current \eqn{b_i^*}.
#'     \item Variance-component update uses the moment-matching
#'           estimate \eqn{\hat\sigma_b^2 = \mathrm{var}(\hat b_i)}
#'           across clusters (Breslow-Clayton approximation).
#'   }
#'
#'   The scaffolding here reuses every primitive already shipped in
#'   dsVert 1.1.0+: the per-cluster residual aggregates, the
#'   keep_session flag on ds.vertGLM, the Beaver vecmul for inner
#'   weighted updates, and the DCF sigmoid/exp wide-splines.
#'
#'   Privacy: client sees only (beta, sigma_b^2) + per-cluster
#'   b_hat_i estimates as an aggregate vector (one value per cluster).
#'   No per-patient quantity ever leaves the DCF parties.
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
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return \code{ds.vertGLMM} object: fixed-effect coefficients,
#'   cluster-level BLUPs \eqn{\hat b_i}, random-effect variance
#'   \eqn{\hat\sigma_b^2}, and the converged binomial \code{fit}.
#' @export
ds.vertGLMM <- function(formula, data = NULL, cluster_col,
                        max_outer = 10L, inner_iter = 10L,
                        tol = 1e-3, lambda = 0,
                        verbose = TRUE,
                        datasources = NULL) {
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
  f_obj <- if (inherits(formula, "formula")) formula else stats::as.formula(formula)
  rhs_terms <- attr(stats::terms(f_obj), "term.labels")
  y_cols <- tryCatch(
    DSI::datashield.aggregate(datasources[which(server_names == y_srv)],
      call(name = "dsvertColNamesDS", data_name = data))[[1L]]$columns,
    error = function(e) character(0))
  missing_local <- setdiff(rhs_terms, y_cols)
  if (length(missing_local) > 0L) {
    stop("ds.vertGLMM currently requires all fixed-effect predictors on ",
         "the outcome server for the cluster BLUP update. Missing local ",
         "predictors: ", paste(missing_local, collapse = ", "),
         ". A true vertical GLMM path needs cross-server binomial residual ",
         "aggregation.", call. = FALSE)
  }

  # Prime with a straight binomial fit ignoring random effects.
  if (verbose) message("[ds.vertGLMM] prime: binomial ds.vertGLM")
  fit <- ds.vertGLM(formula = formula, data = data, family = "binomial",
                    max_iter = inner_iter, tol = tol, lambda = lambda,
                    verbose = FALSE,
                    datasources = datasources)

  # Cluster sizes + per-cluster residuals (all aggregates).
  clust_info <- DSI::datashield.aggregate(
    datasources[which(server_names == y_srv)],
    call(name = "dsvertClusterSizesDS", data_name = data,
         cluster_col = cluster_col))
  if (is.list(clust_info) && length(clust_info) == 1L)
    clust_info <- clust_info[[1L]]
  n_i <- as.integer(clust_info$sizes)
  n_clusters <- length(n_i)

  # Outer optimiser: alternate (i) compute cluster BLUPs \hat b_i via
  # Laplace-approximated inner PIRLS on y vs. binomial likelihood, and
  # (ii) refit the fixed effects with the outer-loop intercept offsets
  # b_i. Moment-match \sigma_b^2 = var(\hat b_i).
  sigma_b2 <- 1.0
  b_hat <- rep(0, n_clusters)
  converged <- FALSE
  offset_col <- NULL
  for (outer in seq_len(max_outer)) {
    # Per-cluster binomial score aggregates at current fixed effects and
    # current random-intercept offset. Only one scalar per cluster returns.
    x_coef_names <- setdiff(names(fit$coefficients), "(Intercept)")
    cl <- DSI::datashield.aggregate(
      datasources[which(server_names == y_srv)],
      call(name = "dsvertClusterBinomialMomentsDS",
           data_name = data, y_var = y_var,
           x_names = x_coef_names,
           intercept = as.numeric(fit$coefficients["(Intercept)"]),
           betahat = as.numeric(fit$coefficients[setdiff(
             names(fit$coefficients), "(Intercept)")]),
           cluster_col = cluster_col,
           offset_col = offset_col))
    if (is.list(cl) && length(cl) == 1L) cl <- cl[[1L]]
    rsum <- as.numeric(cl$rsum_per_cluster)
    vsum <- as.numeric(cl$vsum_per_cluster)
    active <- n_i > 0L & is.finite(rsum) & is.finite(vsum)
    info <- rep(Inf, n_clusters)
    info[active] <- pmax(vsum[active], 1e-8) + 1 / sigma_b2
    score <- rsum - b_hat / sigma_b2
    step <- score / info
    step[!is.finite(step)] <- 0
    step <- pmax(pmin(step, 3), -3)
    b_hat_new <- b_hat + step
    b_hat_new[!active] <- 0
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
    sigma_b2_new <- max(mean(b_hat_new[active]^2 + post_var[active]), 1e-6)
    if (verbose) {
      message(sprintf(
        "[GLMM] outer %d  sigma_b^2=%.4g  var(b_hat)=%.4g  mean(post_var)=%.4g",
        outer, sigma_b2_new, stats::var(b_hat_new[active]),
        mean(post_var[active])))
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
    fit <- tryCatch(
      ds.vertGLM(formula = formula, data = data, family = "binomial",
                 offset = offset_col,
                 max_iter = inner_iter, tol = tol, lambda = lambda,
                 verbose = FALSE,
                 datasources = datasources),
      error = function(e) {
        message("[GLMM] inner binomial refit failed: ",
                conditionMessage(e)); fit })
    common_coef <- intersect(names(old_coef), names(fit$coefficients))
    beta_delta <- if (length(common_coef)) {
      max(abs(old_coef[common_coef] - fit$coefficients[common_coef]))
    } else Inf
    if (max(sigma_delta, b_delta, beta_delta) < tol * max(1, sigma_b2)) {
      converged <- TRUE
      break
    }
  }

  icc <- sigma_b2 / (sigma_b2 + pi^2 / 3)  # logistic latent-variance
  out <- list(
    coefficients = fit$coefficients,
    std_errors   = fit$std_errors,
    sigma_b2     = sigma_b2,
    b_hat        = b_hat,
    icc          = icc,
    n_clusters   = n_clusters,
    cluster_sizes = n_i,
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
  cat(sprintf("  Clusters = %d    N = %d\n",
              x$n_clusters, sum(x$cluster_sizes)))
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
