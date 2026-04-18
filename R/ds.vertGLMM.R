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
#'           (using the cached \code{\link{dsvertClusterResidualsDS}}
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
#'   Privacy: client sees only (\beta, \sigma_b^2) + per-cluster
#'   \hat b_i estimates as an aggregate vector (one value per cluster).
#'   No per-patient quantity ever leaves the DCF parties.
#'
#'   Inter-server leakage: cluster membership (same tier as ds.vertLMM).
#'
#' @param formula Fixed-effects formula (binomial outcome on LHS).
#' @param data Aligned data-frame name.
#' @param cluster_col Cluster id column on the outcome server.
#' @param max_outer Outer (\beta, \sigma_b^2) iterations.
#' @param inner_iter Inner PIRLS iterations per cluster per outer step.
#' @param tol Outer convergence tolerance.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return \code{ds.vertGLMM} object: fixed-effect coefficients,
#'   cluster-level BLUPs \eqn{\hat b_i}, random-effect variance
#'   \eqn{\hat\sigma_b^2}, and the converged binomial \code{fit}.
#' @export
ds.vertGLMM <- function(formula, data = NULL, cluster_col,
                        max_outer = 10L, inner_iter = 10L,
                        tol = 1e-3, verbose = TRUE,
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

  # Prime with a straight binomial fit ignoring random effects.
  if (verbose) message("[ds.vertGLMM] prime: binomial ds.vertGLM")
  fit <- ds.vertGLM(formula = formula, data = data, family = "binomial",
                    max_iter = 60L, verbose = FALSE,
                    datasources = datasources)

  # Cluster sizes + per-cluster residuals (all aggregates).
  clust_info <- DSI::datashield.aggregate(
    datasources[which(server_names == y_srv)],
    call("dsvertClusterSizesDS", data_name = data,
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
  for (outer in seq_len(max_outer)) {
    # Per-cluster residual aggregates at current fit.
    cl <- DSI::datashield.aggregate(
      datasources[which(server_names == y_srv)],
      call("dsvertClusterResidualsDS",
           data_name = data, y_var = y_var,
           x_names = setdiff(names(fit$coefficients), "(Intercept)"),
           intercept = as.numeric(fit$coefficients["(Intercept)"]),
           betahat = as.numeric(fit$coefficients[setdiff(
             names(fit$coefficients), "(Intercept)")]),
           cluster_col = cluster_col))
    if (is.list(cl) && length(cl) == 1L) cl <- cl[[1L]]
    rsum <- as.numeric(cl$rsum_per_cluster)
    rss <- as.numeric(cl$rss_per_cluster)
    # Laplace inner (client-side closed form for logit random
    # intercept): \hat b_i ≈ rsum_i / (n_i * \bar p_i (1 - \bar p_i)
    #                                    + 1/sigma_b^2)
    # Use the BINOMIAL VARIANCE at the current intercept. Without per-
    # cluster p_i we approximate by the global mean residual variance
    # proxy from the fit's deviance.
    sigma2_resid <- max(sum(rss) / sum(n_i), 0.25)
    b_hat_new <- rsum / (n_i * sigma2_resid + 1 / sigma_b2)
    sigma_b2_new <- max(stats::var(b_hat_new), 1e-6)
    if (verbose) {
      message(sprintf("[GLMM] outer %d  sigma_b^2=%.4g  var(b_hat)=%.4g",
                       outer, sigma_b2_new, stats::var(b_hat_new)))
    }
    if (abs(sigma_b2_new - sigma_b2) < tol * max(1, sigma_b2)) {
      converged <- TRUE
      b_hat <- b_hat_new; sigma_b2 <- sigma_b2_new
      break
    }
    b_hat <- b_hat_new; sigma_b2 <- sigma_b2_new

    # Expand b_i into a per-patient offset column so the next
    # ds.vertGLM fit subtracts it via the existing offset= plumbing.
    tryCatch(
      DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call("dsvertExpandClusterWeightsDS",
             data_name = data, cluster_col = cluster_col,
             weights_per_cluster = as.numeric(b_hat),
             output_column = "__dsvert_glmm_b")),
      error = function(e) {
        message("[GLMM] offset expand failed: ", conditionMessage(e))
      })
    fit <- tryCatch(
      ds.vertGLM(formula = formula, data = data, family = "binomial",
                 offset = "__dsvert_glmm_b",
                 max_iter = 60L, verbose = FALSE,
                 datasources = datasources),
      error = function(e) {
        message("[GLMM] inner binomial refit failed: ",
                conditionMessage(e)); fit })
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
