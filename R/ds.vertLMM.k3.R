#' @title Federated linear mixed model via REML 1-D profile (K=3)
#' @description Random-intercept LMM
#'   \deqn{y_{ij} = X_{ij} \beta + b_i + \varepsilon_{ij},
#'         \quad b_i \sim \mathcal{N}(0, \sigma_b^2),\;
#'         \varepsilon_{ij} \sim \mathcal{N}(0, \sigma^2)}
#'   estimated by REML 1-D profile over the intra-cluster
#'   correlation \eqn{\rho = \sigma_b^2 / (\sigma_b^2 + \sigma^2/n_i)}.
#'   For balanced designs the GLS-transformed score equations are
#'   solved exactly by a weighted Gaussian \code{ds.vertGLM} call; the
#'   profile likelihood at each candidate \eqn{\rho} is
#'   evaluated from the inner fit's \code{$deviance} plus the
#'   per-cluster log-determinant closed form (Christensen 2019 §A.3,
#'   Pinheiro & Bates 2000 §2.4).
#'
#'   K=3 implementation reuses the secure-aggregation \code{ds.vertGLM}
#'   pipeline at each profile evaluation — no new MPC primitives. The
#'   per-cluster GLS weights are a public function of \eqn{\rho} and
#'   the cluster sizes (already disclosed under
#'   \code{datashield.privacyLevel}), so D-INV-1..5 are preserved
#'   identically to the K=2 \code{ds.vertLMM} path.
#'
#'   Outer optimisation is golden-section search over
#'   \eqn{\rho \in (\rho_{\min}, \rho_{\max})}; default
#'   \eqn{(0.001, 0.999)} so the search never enters the boundary
#'   singularities at \eqn{\rho = 0} (no random effect — degenerate
#'   GLS) or \eqn{\rho = 1} (infinite \eqn{\sigma_b^2}).
#'
#' @param formula Fixed-effects formula \code{y ~ X}.
#' @param data Aligned data-frame name on each server.
#' @param cluster_col Cluster id column on the outcome server.
#' @param rho_lo,rho_hi Profile search interval (default 0.001, 0.999).
#' @param tol Profile convergence tolerance on \eqn{\rho} (default 1e-3).
#' @param max_outer Maximum golden-section iterations.
#' @param verbose Print progress.
#' @param datasources DataSHIELD K=3 connections.
#' @return list of class \code{ds.vertLMM.k3} with
#'   \code{coefficients}, \code{sigma_b2}, \code{sigma2}, \code{rho_hat},
#'   \code{n_clusters}, \code{cluster_sizes}, and the inner ds.glm fit.
#' @references
#' Christensen, R. H. B. (2019). \emph{Linear Models}. \code{ordinal::clm.fit}
#'   §A.3 — REML profile likelihood for variance components.
#' Pinheiro, J. C. & Bates, D. M. (2000). \emph{Mixed-Effects Models in
#'   S/S-PLUS}, §2.4.
#' Laird, N. M. & Ware, J. H. (1982). Random-effects models for
#'   longitudinal data. \emph{Biometrics}, 38(4), 963-974.
#' Lindstrom, M. J. & Bates, D. M. (1990). Nonlinear mixed effects
#'   models for repeated measures data. \emph{Biometrics}, 46(3),
#'   673-687.
#' @seealso \code{\link{ds.vertGLM}}, \code{\link{ds.vertLMM}} (K=2
#'   exact closed-form path).
#' @export
ds.vertLMM.k3 <- function(formula, data, cluster_col,
                           rho_lo = 0.001, rho_hi = 0.999,
                           tol = 1e-3, max_outer = 30L,
                           verbose = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (length(datasources) < 3L)
    stop("ds.vertLMM.k3 requires K=3 connections (got ",
         length(datasources), ")", call. = FALSE)
  if (!inherits(formula, "formula"))
    stop("formula must be a y ~ X formula", call. = FALSE)

  ## Locate cluster_col on one of the servers (must be the outcome
  ## server — same constraint as ds.vertGLMM).
  server_names <- names(datasources)
  cluster_srv <- NULL
  for (.srv in server_names) {
    .ci <- which(server_names == .srv)
    cols <- tryCatch(
      DSI::datashield.aggregate(datasources[.ci],
        call("dsvertColNamesDS", data_name = data))[[1]]$columns,
      error = function(e) NULL)
    if (!is.null(cols) && cluster_col %in% cols) {
      cluster_srv <- .srv; break
    }
  }
  if (is.null(cluster_srv))
    stop("cluster_col '", cluster_col, "' not found on any server",
         call. = FALSE)

  ## Per-cluster sizes — public under privacy-level threshold.
  if (verbose) message("[ds.vertLMM.k3] Querying cluster sizes …")
  ci <- which(server_names == cluster_srv)
  ci_info <- DSI::datashield.aggregate(datasources[ci],
    call("dsvertClusterSizesDS", data_name = data,
         cluster_col = cluster_col))
  if (is.list(ci_info) && length(ci_info) == 1L) ci_info <- ci_info[[1]]
  n_i <- as.integer(ci_info$sizes)
  n_clusters <- length(n_i)
  n_total <- sum(n_i)

  ## REML profile log-likelihood. For balanced or near-balanced
  ## random-intercept design the per-row GLS weight is
  ##   w_ij(ρ) = 1 - ρ/(1 + ρ(n_i - 1))
  ## which gives ds.vertGLM the same fixed-effect estimator as the
  ## REML solution at that ρ. The profile-loglik-up-to-constants is
  ##   l_p(ρ) = -½ Σ_i (n_i - 1) log(1 - ρ) - ½ Σ_i log(1 + ρ(n_i - 1))
  ##           - (n - p)/2 · log(RSS_p(ρ) / (n - p))
  ## with RSS_p(ρ) recovered from the weighted GLM fit's deviance.
  ## (Pinheiro & Bates 2000 §2.4.1 Eq. 2.20.)

  .compute_weights_per_cluster <- function(rho) {
    ## w_i = 1 - ρ/(1 + ρ(n_i - 1)).
    1 - rho / (1 + rho * (n_i - 1))
  }

  .neg_profile_loglik <- function(rho) {
    if (verbose) message(sprintf("[ds.vertLMM.k3] inner fit at ρ=%.4f", rho))
    w_per_cluster <- .compute_weights_per_cluster(rho)
    DSI::datashield.aggregate(datasources[ci],
      call("dsvertExpandClusterWeightsDS",
           data_name = data, cluster_col = cluster_col,
           weights_per_cluster = as.numeric(w_per_cluster),
           output_column = "__lmm_k3_w"))
    fit <- ds.vertGLM(formula = formula, data = data,
                      family = "gaussian",
                      weights = "__lmm_k3_w",
                      max_iter = 30L, tol = 1e-6,
                      verbose = FALSE,
                      datasources = datasources)
    rss <- as.numeric(fit$deviance)
    p <- length(fit$coefficients)
    df_resid <- max(n_total - p, 1L)
    sigma2 <- rss / df_resid
    log_det <- sum((n_i - 1) * log(pmax(1 - rho, 1e-12))) +
               sum(log(1 + rho * (n_i - 1)))
    -(- 0.5 * log_det - 0.5 * df_resid * log(sigma2 / df_resid))
  }

  ## Golden-section search over ρ.
  if (verbose) message("[ds.vertLMM.k3] golden-section over ρ ∈ [",
                        rho_lo, ", ", rho_hi, "] …")
  phi <- (sqrt(5) - 1) / 2
  a <- rho_lo; b_ <- rho_hi
  c_ <- b_ - phi * (b_ - a); d_ <- a + phi * (b_ - a)
  fc <- .neg_profile_loglik(c_); fd <- .neg_profile_loglik(d_)
  for (it in seq_len(max_outer)) {
    if (abs(b_ - a) < tol) break
    if (fc < fd) {
      b_ <- d_; d_ <- c_; fd <- fc
      c_ <- b_ - phi * (b_ - a); fc <- .neg_profile_loglik(c_)
    } else {
      a <- c_; c_ <- d_; fc <- fd
      d_ <- a + phi * (b_ - a); fd <- .neg_profile_loglik(d_)
    }
    if (verbose) message(sprintf("  it=%d  ρ ∈ [%.4f, %.4f]", it, a, b_))
  }
  rho_hat <- (a + b_) / 2

  ## Final fit at ρ̂.
  w_final <- .compute_weights_per_cluster(rho_hat)
  DSI::datashield.aggregate(datasources[ci],
    call("dsvertExpandClusterWeightsDS",
         data_name = data, cluster_col = cluster_col,
         weights_per_cluster = as.numeric(w_final),
         output_column = "__lmm_k3_w"))
  fit_final <- ds.vertGLM(formula = formula, data = data,
                           family = "gaussian",
                           weights = "__lmm_k3_w",
                           max_iter = 60L, tol = 1e-7,
                           verbose = FALSE,
                           datasources = datasources)

  ## Variance components from the converged fit.
  rss <- as.numeric(fit_final$deviance)
  p <- length(fit_final$coefficients)
  df_resid <- max(n_total - p, 1L)
  sigma2_hat <- rss / df_resid
  ## σ_b² = ρ · σ² · (1 / (1 - ρ)) at the balanced-design fix point.
  sigma_b2_hat <- rho_hat * sigma2_hat / max(1 - rho_hat, 1e-12)

  out <- list(
    coefficients  = fit_final$coefficients,
    rho_hat       = rho_hat,
    sigma_b2      = sigma_b2_hat,
    sigma2        = sigma2_hat,
    n_clusters    = n_clusters,
    cluster_sizes = n_i,
    fit           = fit_final,
    family        = "Gaussian (REML 1-D profile, K=3)",
    call          = match.call())
  class(out) <- c("ds.vertLMM.k3", "list")
  out
}

#' @export
print.ds.vertLMM.k3 <- function(x, ...) {
  cat("dsVert LMM (K=3, REML 1-D profile)\n")
  cat(sprintf("  Clusters = %d   Total N = %d\n",
              x$n_clusters, sum(x$cluster_sizes)))
  cat(sprintf("  ρ̂ = %.4f   σ̂_b² = %.4g   σ̂² = %.4g\n",
              x$rho_hat, x$sigma_b2, x$sigma2))
  cat("\nFixed effects (β):\n")
  print(round(x$coefficients, 5))
  invisible(x)
}
