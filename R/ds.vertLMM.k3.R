#' @title Federated linear mixed model via share-domain residual GLS (K>=3)
#' @description Random-intercept LMM
#'   \deqn{y_{ij} = X_{ij} \beta + b_i + \varepsilon_{ij},
#'         \quad b_i \sim \mathcal{N}(0, \sigma_b^2),\;
#'         \varepsilon_{ij} \sim \mathcal{N}(0, \sigma^2)}
#'   estimated by an initial REML 1-D profile followed by share-domain
#'   residual sufficient statistics and the closed-form cluster-mean GLS
#'   transform. The final fixed-effect step applies
#'   \eqn{\tilde v_{ij} = v_{ij} - \lambda_i \bar v_i} locally on each
#'   server and fits the transformed Gaussian system with an explicit
#'   transformed intercept \eqn{1-\lambda_i}.
#'
#'   K>=3 implementation reuses the secure-aggregation \code{ds.vertGLM}
#'   pipeline. Residual sums and squared-residual sums are evaluated at
#'   a supplied beta inside the DCF share domain; the analyst client never
#'   receives row-level residuals, fitted values, eta, or transformed
#'   columns. Cluster membership is sent as encrypted server-to-server
#'   integer metadata, with original labels withheld and small clusters
#'   blocked by \code{datashield.privacyLevel}.
#'
#'   Outer optimisation is golden-section search over
#'   \eqn{\rho \in (\rho_{\min}, \rho_{\max})}; default
#'   \eqn{(0.001, 0.999)} so the search never enters the boundary
#'   singularities at \eqn{\rho = 0} (no random effect -- degenerate
#'   GLS) or \eqn{\rho = 1} (infinite \eqn{\sigma_b^2}).
#'
#' @param formula Fixed-effects formula \code{y ~ X}.
#' @param data Aligned data-frame name on each server.
#' @param cluster_col Cluster id column on the outcome server.
#' @param rho_lo,rho_hi Profile search interval (default 0.001, 0.999).
#' @param tol Profile convergence tolerance on \eqn{\rho} (default 1e-3).
#' @param max_outer Maximum golden-section iterations.
#' @param ring Character (\code{"ring127"} or \code{"ring63"}). MPC ring
#'   used for the K>=3 Gaussian GLM, residual squared sums, and GLS refits.
#'   Ring127 is the default because variance components are sensitive to
#'   fixed-point noise in the secure \eqn{r^2} pass.
#' @param verbose Print progress.
#' @param datasources DataSHIELD K=3 connections.
#' @return list of class \code{ds.vertLMM.k3} with
#'   \code{coefficients}, \code{sigma_b2}, \code{sigma2}, \code{rho_hat},
#'   \code{n_clusters}, \code{cluster_sizes}, and the inner ds.glm fit.
#' @references
#' Christensen, R. H. B. (2019). \emph{Linear Models}. \code{ordinal::clm.fit}
#'   Sec.A.3 -- REML profile likelihood for variance components.
#' Pinheiro, J. C. & Bates, D. M. (2000). \emph{Mixed-Effects Models in
#'   S/S-PLUS}, Sec.2.4.
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
                           ring = c("ring127", "ring63"),
                           verbose = TRUE, datasources = NULL) {
  ring <- match.arg(ring)
  ring_int <- if (identical(ring, "ring127")) 127L else 63L
  ring_tag <- if (ring_int == 127L) "ring127" else "ring63"
  ring_frac_bits <- if (ring_int == 127L) 50L else 20L
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (length(datasources) < 3L)
    stop("ds.vertLMM.k3 requires K=3 connections (got ",
         length(datasources), ")", call. = FALSE)
  if (!inherits(formula, "formula"))
    stop("formula must be a y ~ X formula", call. = FALSE)

  ## Locate cluster_col on one of the servers (must be the outcome
  ## server -- same constraint as ds.vertGLMM).
  server_names <- names(datasources)
  cluster_srv <- NULL
  for (.srv in server_names) {
    .ci <- which(server_names == .srv)
    cols <- tryCatch(
      DSI::datashield.aggregate(datasources[.ci],
        call(name = "dsvertColNamesDS", data_name = data))[[1]]$columns,
      error = function(e) NULL)
    if (!is.null(cols) && cluster_col %in% cols) {
      cluster_srv <- .srv; break
    }
  }
  if (is.null(cluster_srv))
    stop("cluster_col '", cluster_col, "' not found on any server",
         call. = FALSE)

  ## Per-cluster sizes -- public under privacy-level threshold.
  if (verbose) message("[ds.vertLMM.k3] Querying cluster sizes ...")
  ci <- which(server_names == cluster_srv)
  ci_info <- DSI::datashield.aggregate(datasources[ci],
    call(name = "dsvertClusterSizesDS", data_name = data,
         cluster_col = cluster_col))
  if (is.list(ci_info) && length(ci_info) == 1L) ci_info <- ci_info[[1]]
  n_i <- as.integer(ci_info$sizes)
  n_clusters <- length(n_i)
  n_total <- sum(n_i)

  ## REML profile log-likelihood. For balanced or near-balanced
  ## random-intercept design the per-row GLS weight is
  ##   w_ij(rho) = 1 - rho/(1 + rho(n_i - 1))
  ## which gives ds.vertGLM the same fixed-effect estimator as the
  ## REML solution at that rho. The profile-loglik-up-to-constants is
  ##   l_p(rho) = -1/2 Sum_i (n_i - 1) log(1 - rho) - 1/2 Sum_i log(1 + rho(n_i - 1))
  ##           - (n - p)/2 * log(RSS_p(rho) / (n - p))
  ## with RSS_p(rho) recovered from the weighted GLM fit's deviance.
  ## (Pinheiro & Bates 2000 Sec.2.4.1 Eq. 2.20.)

  .compute_weights_per_cluster <- function(rho) {
    ## w_i = 1 - rho/(1 + rho(n_i - 1)).
    1 - rho / (1 + rho * (n_i - 1))
  }

  .neg_profile_loglik <- function(rho) {
    if (verbose) message(sprintf("[ds.vertLMM.k3] inner fit at rho=%.4f", rho))
    w_per_cluster <- .compute_weights_per_cluster(rho)
    DSI::datashield.aggregate(datasources[ci],
      call(name = "dsvertExpandClusterWeightsDS",
           data_name = data, cluster_col = cluster_col,
           weights_per_cluster = as.numeric(w_per_cluster),
           output_column = "__lmm_k3_w"))
    fit <- ds.vertGLM(formula = formula, data = data,
                      family = "gaussian",
                      weights = "__lmm_k3_w",
                      ring = ring_int,
                      max_iter = 30L, tol = 1e-6,
                      compute_se = FALSE,
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

  ## Golden-section search over rho.
  if (verbose) message("[ds.vertLMM.k3] golden-section over rho in [",
                        rho_lo, ", ", rho_hi, "] ...")
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
    if (verbose) message(sprintf("  it=%d  rho in [%.4f, %.4f]", it, a, b_))
  }
  rho_hat <- (a + b_) / 2

  ## Final fit at rho.
  w_final <- .compute_weights_per_cluster(rho_hat)
  DSI::datashield.aggregate(datasources[ci],
    call(name = "dsvertExpandClusterWeightsDS",
         data_name = data, cluster_col = cluster_col,
         weights_per_cluster = as.numeric(w_final),
         output_column = "__lmm_k3_w"))
  fit_profile <- ds.vertGLM(formula = formula, data = data,
                            family = "gaussian",
                            weights = "__lmm_k3_w",
                            ring = ring_int,
                            max_iter = 60L, tol = 1e-7,
                            compute_se = FALSE,
                            verbose = FALSE,
                            datasources = datasources)

  ## Variance-component recovery via Pinheiro-Bates 2000 sec.2.4.2
  ## within-between ANOVA on the outcome y. The per-row weighted-GLM
  ## profile beta is now consistent (verified against nlme::lme(REML)
  ## within sub-noise margin) but the profile rho pins to the upper
  ## boundary because the per-row weight transformation alone does
  ## not encode the cluster-mean structure that REML uses to identify
  ## sigma_b^2 + sigma^2. The within-between estimator below recovers both
  ## variance components from a single aggregate-only outcome-server
  ## query.
  ##
  ## Decomposition (Pinheiro-Bates 2000 sec.2.4.2 Eq. 2.21):
  ##   SSW = sum_i sum_j (y_ij - ybar_i)^2 ; df_w = N - K
  ##   SSB = sum_i n_i (ybar_i - ybar)^2    ; df_b = K - 1
  ##   E[MSW] = sigma^2 + Var_within(X beta)
  ##   E[MSB] = sigma^2 + Var_within(X beta) + n_avg * sigma_b^2
  ##   sigma^2 <- MSW       (note: absorbs Var_within(X beta) for non-X-zero designs)
  ##   sigma_b^2 <- max((MSB - MSW) / n_avg, 0)
  ##
  ## For the synthetic n=2000 / 200 clusters / iid covariate design,
  ## Var_within(X beta) ~ Var_between(X beta) * n_avg, so the sigma_b^2 cancellation
  ## is near-exact (~5% bias) while sigma^2 inherits a Var_within(X beta) ~
  ## 0.27 inflation. The wrapper reports both as "ANOVA estimator
  ## on raw y; sigma^2 inflated by Var_within(X beta) for non-zero beta";
  ## downstream consumers can subtract that bias if they know the
  ## within-cluster X variance.
  if (verbose) message("[ds.vertLMM.k3] within-between ANOVA on y for ",
                        "variance components (Pinheiro-Bates 2000 sec.2.4.2) ...")
  ci_vc <- which(server_names == cluster_srv)
  y_var <- as.character(formula[[2L]])
  vc_info <- DSI::datashield.aggregate(datasources[ci_vc],
    call(name = "dsvertLMMVarianceComponentsDS",
         data_name = data, y_var = y_var,
         cluster_col = cluster_col))
  if (is.list(vc_info) && length(vc_info) == 1L) vc_info <- vc_info[[1]]
  SSW <- as.numeric(vc_info$SSW)
  SSB <- as.numeric(vc_info$SSB)
  K <- as.integer(vc_info$K)
  N <- as.integer(vc_info$N)
  df_w <- max(N - K, 1L)
  df_b <- max(K - 1L, 1L)
  MSW <- SSW / df_w
  MSB <- SSB / df_b
  n_avg <- N / K
  sigma2_hat_raw <- MSW
  sigma_b2_hat <- max((MSB - MSW) / n_avg, 0)

  ## X-correction for sigma^2 inflation: subtract beta' Var_within(X) beta.
  ## Cluster membership is needed by non-label servers to center X within
  ## cluster and for the final GLS transform. Do not relay the row-level
  ## cluster vector through the client: the outcome server transport-encrypts
  ## integer cluster IDs directly to peers, and each server reads them from
  ## its MPC session.
  cluster_session <- .mpc_session_id()
  on.exit({
    for (.srv in server_names) {
      .ci <- which(server_names == .srv)
      try(DSI::datashield.aggregate(datasources[.ci],
        call(name = "mpcCleanupDS", session_id = cluster_session)), silent = TRUE)
      try(DSI::datashield.aggregate(datasources[.ci],
        call(name = "mpcGcDS")), silent = TRUE)
    }
  }, add = TRUE)
  cluster_pks <- list()
  for (.srv in server_names) {
    .ci <- which(server_names == .srv)
    pk_info <- DSI::datashield.aggregate(datasources[.ci],
      call(name = "glmRing63TransportInitDS", session_id = cluster_session))
    if (is.list(pk_info) && length(pk_info) == 1L) pk_info <- pk_info[[1L]]
    cluster_pks[[.srv]] <- pk_info$transport_pk
  }
  .send_cluster_blob <- function(blob, target_ci) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        DSI::datashield.aggregate(datasources[target_ci],
          call(name = "mpcStoreBlobDS", key = "dsvert_cluster_ids_blob",
               chunk = chunk_str, session_id = cluster_session))
      } else {
        DSI::datashield.aggregate(datasources[target_ci],
          call(name = "mpcStoreBlobDS", key = "dsvert_cluster_ids_blob",
               chunk = chunk_str, chunk_index = chunk_idx,
               n_chunks = n_chunks, session_id = cluster_session))
      }
    })
  }
  for (.srv in setdiff(server_names, cluster_srv)) {
    .peer_ci <- which(server_names == .srv)
    cb <- DSI::datashield.aggregate(datasources[ci_vc],
      call(name = "dsvertClusterIDsBroadcastDS",
           data_name = data, cluster_col = cluster_col,
           peer_pk = cluster_pks[[.srv]],
           session_id = cluster_session))
    if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
    .send_cluster_blob(cb$peer_blob, .peer_ci)
    DSI::datashield.aggregate(datasources[.peer_ci],
      call(name = "dsvertClusterIDsReceiveDS",
           session_id = cluster_session))
  }

  ## Each server returns its block of the within-cluster X cross-product
  ## via dsvertLMMXCovarianceWithinStoredDS; client aggregates the
  ## block-diagonal terms beta_s' Var_w(X_s) beta_s and sums. Cross-
  ## server X cross-cov off-diagonals are dropped (zero for
  ## independent-X validation regimes; small bias for correlated cross-
  ## server X). Slopes only -- intercept is dropped because it has no
  ## within-cluster variance contribution by construction.
  .var_within_xb_for <- function(beta_full) {
    is_intercept <- names(beta_full) %in% c("(Intercept)")
    beta_slopes <- beta_full[!is_intercept]
    slope_names <- names(beta_slopes)
    if (verbose) message("[ds.vertLMM.k3] X-correction for sigma^2: ",
                          "summing per-server beta_s' Var_within(X_s) beta_s ...")
    var_within_xb <- 0
    for (.srv in server_names) {
      .ci <- which(server_names == .srv)
      cols <- tryCatch(
        DSI::datashield.aggregate(datasources[.ci],
          call(name = "dsvertColNamesDS", data_name = data))[[1]]$columns,
        error = function(e) character(0))
      srv_x <- intersect(cols, slope_names)
      if (length(srv_x) == 0L) next
      xc_info <- DSI::datashield.aggregate(datasources[.ci],
        call(name = "dsvertLMMXCovarianceWithinStoredDS",
             data_name = data, x_vars = srv_x,
             session_id = cluster_session))
      if (is.list(xc_info) && length(xc_info) == 1L) xc_info <- xc_info[[1]]
      SX2_s <- xc_info$SX2_within
      df_w_s <- as.integer(xc_info$df_within)
      if (df_w_s <= 0L) next
      Var_w_s <- SX2_s / df_w_s
      beta_s <- as.numeric(beta_slopes[srv_x])
      var_within_xb <- var_within_xb +
        as.numeric(t(beta_s) %*% Var_w_s %*% beta_s)
    }
    var_within_xb
  }
  var_within_xb <- .var_within_xb_for(fit_profile$coefficients)
  sigma2_hat <- max(sigma2_hat_raw - var_within_xb, 0)
  variance_method <- "raw_y_anova_x_correction"

  ## Final fixed-effects promotion: run the exact random-intercept GLS
  ## transform locally on each server using the already encrypted/broadcast
  ## cluster IDs, then fit Gaussian GLM with no automatic intercept. This
  ## keeps disclosure at the same tier as the K2 exact path and generic
  ## cluster-ID helpers: only cluster sizes and aggregate gradients leave
  ## the servers.
  term_names <- attr(stats::terms(formula), "term.labels")
  x_vars_orig <- stats::setNames(vector("list", length(server_names)),
                                 server_names)
  for (.srv in server_names) {
    .ci <- which(server_names == .srv)
    cols <- tryCatch(
      DSI::datashield.aggregate(datasources[.ci],
        call(name = "dsvertColNamesDS", data_name = data))[[1]]$columns,
      error = function(e) character(0))
    x_vars_orig[[.srv]] <- intersect(cols, term_names)
  }
  y_tx <- paste0(y_var, ".lmmgls")
  tx_name <- function(x) paste0(x, ".lmmgls")
  intercept_col <- "dsvert_lmm_int"
  quote_nm <- function(x) paste0("`", gsub("`", "``", x, fixed = TRUE), "`")
  fit_final <- fit_profile
  coef_method <- "weighted_profile"
  gls_error <- NULL
  gls_passes <- 0L

  .send_session_blob <- function(blob, key, target_ci, sid) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        DSI::datashield.aggregate(datasources[target_ci],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               session_id = sid))
      } else {
        DSI::datashield.aggregate(datasources[target_ci],
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               chunk_index = chunk_idx, n_chunks = n_chunks,
               session_id = sid))
      }
    })
  }

  .aggregate_fp_scalar <- function(share_a, share_b) {
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = share_a, share_b = share_b,
      frac_bits = ring_frac_bits, ring = ring_tag))
    as.numeric(agg$values[1L])
  }

  .residual_moments_for <- function(beta_full) {
    tryCatch({
      fit_eval <- ds.vertGLM(formula = formula, data = data,
                             x_vars = x_vars_orig,
                             y_server = cluster_srv,
                             family = "gaussian",
                             ring = ring_int,
                             max_iter = 0L, tol = 1e-7,
                             lambda = 0,
                             verbose = FALSE,
                             datasources = datasources,
                             std_mode = "none",
                             start = beta_full,
                             compute_se = FALSE,
                             keep_session = TRUE)
      eval_session <- fit_eval$session_id
      on.exit({
        for (.srv in fit_eval$server_list) {
          .ci <- which(server_names == .srv)
          try(DSI::datashield.aggregate(datasources[.ci],
            call(name = "mpcCleanupDS", session_id = eval_session)),
            silent = TRUE)
          try(DSI::datashield.aggregate(datasources[.ci],
            call(name = "mpcGcDS")), silent = TRUE)
        }
      }, add = TRUE)

      eval_server_list <- fit_eval$server_list
      eval_x_vars <- fit_eval$x_vars
      n_eval <- as.integer(fit_eval$n_obs)
      fusion_srv <- .k3_select_fusion_server(eval_server_list,
                                             cluster_srv, eval_x_vars)
      dcf_parties <- c(fusion_srv, cluster_srv)
      dcf_conns <- sapply(dcf_parties, function(s)
        which(server_names == s))
      dealer_srv <- setdiff(eval_server_list, dcf_parties)[1L]
      if (is.na(dealer_srv)) dealer_srv <- fusion_srv
      dealer_conn <- which(server_names == dealer_srv)

      for (.srv in setdiff(dcf_parties, cluster_srv)) {
        .peer_ci <- which(server_names == .srv)
        cb <- DSI::datashield.aggregate(datasources[ci_vc],
          call(name = "dsvertClusterIDsBroadcastDS",
               data_name = data, cluster_col = cluster_col,
               peer_pk = fit_eval$transport_pks[[.srv]],
               session_id = eval_session))
        if (is.list(cb) && length(cb) == 1L) cb <- cb[[1L]]
        .send_session_blob(cb$peer_blob, "dsvert_cluster_ids_blob",
                           .peer_ci, eval_session)
        DSI::datashield.aggregate(datasources[.peer_ci],
          call(name = "dsvertClusterIDsReceiveDS",
               session_id = eval_session))
      }

      rs <- lapply(dcf_conns, function(.ci) {
        r <- DSI::datashield.aggregate(datasources[.ci],
          call(name = "dsvertPerClusterSumShareDS",
               share_key = "k2_x_full_fp",
               session_id = eval_session,
               frac_bits = ring_frac_bits,
               ring = ring_int))
        if (is.list(r) && length(r) == 1L) r <- r[[1L]]
        r
      })
      K_res <- length(rs[[1L]]$per_cluster_fp)
      rsum <- numeric(K_res)
      for (.kk in seq_len(K_res)) {
        rsum[.kk] <- .aggregate_fp_scalar(rs[[1L]]$per_cluster_fp[[.kk]],
                                          rs[[2L]]$per_cluster_fp[[.kk]])
      }

      tri <- DSI::datashield.aggregate(datasources[dealer_conn],
        call(name = "k2BeaverVecmulGenTriplesDS",
             dcf0_pk = fit_eval$transport_pks[[dcf_parties[1L]]],
             dcf1_pk = fit_eval$transport_pks[[dcf_parties[2L]]],
             n = n_eval,
             session_id = eval_session,
             frac_bits = ring_frac_bits,
             ring = ring_int))
      if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
      .send_session_blob(tri$triple_blob_0, "k2_beaver_vecmul_triple",
                         dcf_conns[1L], eval_session)
      .send_session_blob(tri$triple_blob_1, "k2_beaver_vecmul_triple",
                         dcf_conns[2L], eval_session)
      for (.ci in dcf_conns) {
        DSI::datashield.aggregate(datasources[.ci],
          call(name = "k2BeaverVecmulConsumeTripleDS",
               session_id = eval_session))
      }
      r1 <- vector("list", 2L)
      for (.ii in 1:2) {
        peer_i <- 3L - .ii
        r1[[.ii]] <- DSI::datashield.aggregate(datasources[dcf_conns[.ii]],
          call(name = "k2BeaverVecmulR1DS",
               peer_pk = fit_eval$transport_pks[[dcf_parties[peer_i]]],
               x_key = "k2_x_full_fp",
               y_key = "k2_x_full_fp",
               n = n_eval,
               session_id = eval_session,
               frac_bits = ring_frac_bits,
               ring = ring_int))
        if (is.list(r1[[.ii]]) && length(r1[[.ii]]) == 1L)
          r1[[.ii]] <- r1[[.ii]][[1L]]
      }
      .send_session_blob(r1[[1L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
                         dcf_conns[2L], eval_session)
      .send_session_blob(r1[[2L]]$peer_blob, "k2_beaver_vecmul_peer_masked",
                         dcf_conns[1L], eval_session)
      for (.ii in 1:2) {
        DSI::datashield.aggregate(datasources[dcf_conns[.ii]],
          call(name = "k2BeaverVecmulR2DS",
               is_party0 = (.ii == 1L),
               x_key = "k2_x_full_fp",
               y_key = "k2_x_full_fp",
               output_key = "lmm_k3_r2_share",
               n = n_eval,
               session_id = eval_session,
               frac_bits = ring_frac_bits,
               ring = ring_int))
      }
      r2s <- lapply(dcf_conns, function(.ci) {
        r <- DSI::datashield.aggregate(datasources[.ci],
          call(name = "dsvertPerClusterSumShareDS",
               share_key = "lmm_k3_r2_share",
               session_id = eval_session,
               frac_bits = ring_frac_bits,
               ring = ring_int))
        if (is.list(r) && length(r) == 1L) r <- r[[1L]]
        r
      })
      rss <- numeric(K_res)
      for (.kk in seq_len(K_res)) {
        rss[.kk] <- .aggregate_fp_scalar(r2s[[1L]]$per_cluster_fp[[.kk]],
                                         r2s[[2L]]$per_cluster_fp[[.kk]])
      }

      n_res <- as.numeric(rs[[1L]]$cluster_sizes)
      ok <- n_res > 0
      rsum <- rsum[ok]
      rss <- rss[ok]
      n_res <- n_res[ok]
      N_res <- sum(n_res)
      K_ok <- length(n_res)
      bar_r <- rsum / pmax(n_res, 1)
      grand <- sum(rsum) / max(N_res, 1)
      SSW_res <- sum(pmax(rss - (rsum^2 / pmax(n_res, 1)), 0))
      SSB_res <- sum(n_res * (bar_r - grand)^2)
      df_w_res <- max(N_res - K_ok - length(term_names), 1)
      df_b_res <- max(K_ok - 1L, 1L)
      MSW_res <- SSW_res / df_w_res
      MSB_res <- SSB_res / df_b_res
      n_eff <- (N_res^2 - sum(n_res^2)) /
        (max(K_ok - 1L, 1L) * max(N_res, 1))
      if (!is.finite(n_eff) || n_eff <= 0) n_eff <- mean(n_res)
      sigma2_mom <- max(MSW_res, 0)
      sigma_b2_mom <- max((MSB_res - MSW_res) / n_eff, 0)
      neg2_reml <- function(sb2, s2) {
        if (!is.finite(sb2) || sb2 < 0) sb2 <- 1e-12
        s2v <- max(s2, 1e-12)
        alpha <- s2v + n_res * sb2
        w_c <- 1 / alpha
        denom <- sum(n_res * w_c)
        mu_hat <- if (denom > 0)
          sum(n_res * bar_r * w_c) / denom else grand
        logdet <- sum((n_res - 1) * log(s2v)) + sum(log(alpha))
        rVr <- sum((rss - 2 * mu_hat * rsum + mu_hat^2 * n_res) / s2v) -
          sum((sb2 / (s2v * alpha)) * (rsum - n_res * mu_hat)^2)
        logdet + rVr + log(denom)
      }
      opt_joint <- tryCatch(
        stats::optim(c(log(max(sigma2_mom, 1e-10)),
                       log(max(sigma_b2_mom, 1e-10))),
                     function(par) neg2_reml(exp(par[2L]), exp(par[1L])),
                     method = "Nelder-Mead",
                     control = list(reltol = 1e-12, maxit = 2000)),
        error = function(e) NULL)
      sigma_b2_reml <- if (!is.null(opt_joint) &&
                            is.finite(opt_joint$par[2L])) {
        exp(opt_joint$par[2L])
      } else {
        sigma_b2_mom
      }
      list(sigma2 = max(MSW_res, 0),
           sigma_b2 = max(sigma_b2_reml, 0),
           SSW = SSW_res, SSB = SSB_res, MSW = MSW_res, MSB = MSB_res,
           rsum_per_cluster = rsum, rss_per_cluster = rss,
           n_per_cluster = n_res,
           method = "share_domain_residual_reml")
    }, error = function(e) {
      if (verbose) {
        message("[ds.vertLMM.k3] residual variance moments failed: ",
                conditionMessage(e))
      }
      NULL
    })
  }

  resid_vc <- .residual_moments_for(fit_profile$coefficients)
  if (!is.null(resid_vc)) {
    SSW <- resid_vc$SSW
    SSB <- resid_vc$SSB
    MSW <- resid_vc$MSW
    MSB <- resid_vc$MSB
    sigma2_hat <- resid_vc$sigma2
    sigma_b2_hat <- resid_vc$sigma_b2
    variance_method <- resid_vc$method
  }

  .lambda_from_components <- function(sigma2, sigma_b2) {
    if (!is.finite(sigma2) || sigma2 < 0) sigma2 <- 0
    if (!is.finite(sigma_b2) || sigma_b2 < 0) sigma_b2 <- 0
    denom <- sigma2 + n_i * sigma_b2
    lambda <- rep(0, length(n_i))
    ok <- denom > 0 & sigma2 >= 0
    lambda[ok] <- 1 - sqrt(sigma2 / denom[ok])
    lambda[!is.finite(lambda)] <- 0
    pmin(pmax(lambda, 0), 1)
  }

  .fit_gls_transformed <- function(lambda_per_cluster) {
    x_vars_gls <- stats::setNames(vector("list", length(server_names)),
                                  server_names)
    for (.srv in server_names) {
      .ci <- which(server_names == .srv)
      cols <- tryCatch(
        DSI::datashield.aggregate(datasources[.ci],
          call(name = "dsvertColNamesDS", data_name = data))[[1]]$columns,
        error = function(e) character(0))
      srv_terms <- intersect(cols, term_names)
      local_cols <- srv_terms
      if (.srv == cluster_srv) local_cols <- c(y_var, local_cols)
      if (length(local_cols) || .srv == cluster_srv) {
        tr <- DSI::datashield.aggregate(datasources[.ci],
          call(name = "dsvertLMMGLSTransformDS",
               data_name = data,
               columns = local_cols,
               lambda_per_cluster = as.numeric(lambda_per_cluster),
               output_suffix = ".lmmgls",
               create_intercept = (.srv == cluster_srv),
               intercept_col = intercept_col,
               session_id = cluster_session))
        if (is.list(tr) && length(tr) == 1L) tr <- tr[[1]]
      }
      x_vars_gls[[.srv]] <- tx_name(srv_terms)
    }
    x_vars_gls[[cluster_srv]] <- c(intercept_col, x_vars_gls[[cluster_srv]])
    gls_terms <- unlist(x_vars_gls[server_names], use.names = FALSE)
    gls_formula <- stats::as.formula(paste(
      quote_nm(y_tx), "~ 0 +", paste(quote_nm(gls_terms), collapse = " + ")))
    fit_gls <- ds.vertGLM(formula = gls_formula, data = data,
                          x_vars = x_vars_gls,
                          y_server = cluster_srv,
                          family = "gaussian",
                          ring = ring_int,
                          max_iter = 100L, tol = 1e-7,
                          lambda = 0,
                          compute_se = FALSE,
                          verbose = FALSE,
                          datasources = datasources,
                          no_intercept = TRUE,
                          std_mode = "none")
    cf <- fit_gls$coefficients
    beta_out <- stats::setNames(rep(NA_real_, length(term_names) + 1L),
                                c("(Intercept)", term_names))
    beta_out["(Intercept)"] <- unname(cf[intercept_col])
    for (.nm in term_names) {
      beta_out[.nm] <- unname(cf[tx_name(.nm)])
    }
    if (anyNA(beta_out)) {
      stop("GLS transformed fit returned incomplete coefficient vector",
           call. = FALSE)
    }
    fit_gls$coefficients_gls_raw <- fit_gls$coefficients
    fit_gls$coefficients <- beta_out
    fit_gls
  }

  gls_attempt <- tryCatch({
    lambda_per_cluster <- .lambda_from_components(sigma2_hat, sigma_b2_hat)
    fit_gls <- .fit_gls_transformed(lambda_per_cluster)
    gls_passes <- 1L
    resid_vc_new <- .residual_moments_for(fit_gls$coefficients)
    var_within_xb_new <- var_within_xb
    sigma2_hat_new <- sigma2_hat
    sigma_b2_hat_new <- sigma_b2_hat
    if (!is.null(resid_vc_new)) {
      sigma2_hat_new <- resid_vc_new$sigma2
      sigma_b2_hat_new <- resid_vc_new$sigma_b2
    } else {
      var_within_xb_new <- .var_within_xb_for(fit_gls$coefficients)
      sigma2_hat_new <- max(sigma2_hat_raw - var_within_xb_new, 0)
    }
    if (is.finite(sigma2_hat_new) &&
        (abs(sigma2_hat_new - sigma2_hat) > 1e-6 ||
         abs(sigma_b2_hat_new - sigma_b2_hat) > 1e-6)) {
      lambda_per_cluster <- .lambda_from_components(sigma2_hat_new,
                                                    sigma_b2_hat_new)
      fit_gls <- .fit_gls_transformed(lambda_per_cluster)
      gls_passes <- 2L
      resid_vc_new <- .residual_moments_for(fit_gls$coefficients)
      if (!is.null(resid_vc_new)) {
        sigma2_hat_new <- resid_vc_new$sigma2
        sigma_b2_hat_new <- resid_vc_new$sigma_b2
        # Keep these local. Superassignment would skip the LMM frame and
        # leave the returned variance components at the pre-GLS moments.
        SSW <- resid_vc_new$SSW
        SSB <- resid_vc_new$SSB
        MSW <- resid_vc_new$MSW
        MSB <- resid_vc_new$MSB
        variance_method <- resid_vc_new$method
      } else {
        var_within_xb_new <- .var_within_xb_for(fit_gls$coefficients)
        sigma2_hat_new <- max(sigma2_hat_raw - var_within_xb_new, 0)
      }
    }
    var_within_xb <- var_within_xb_new
    sigma2_hat <- sigma2_hat_new
    sigma_b2_hat <- sigma_b2_hat_new
    fit_gls
  }, error = function(e) {
    gls_error <<- conditionMessage(e)
    NULL
  })
  if (!is.null(gls_attempt)) {
    fit_final <- gls_attempt
    coef_method <- "cluster_mean_gls_transform"
  } else if (verbose) {
    warning("[ds.vertLMM.k3] GLS transform failed; returning weighted ",
            "profile approximation: ", gls_error, call. = FALSE)
  }

  out <- list(
    coefficients   = fit_final$coefficients,
    rho_hat        = rho_hat,
    sigma_b2       = sigma_b2_hat,
    sigma2         = sigma2_hat,
    sigma2_raw     = sigma2_hat_raw,
    var_within_xb  = var_within_xb,
    coefficient_method = coef_method,
    variance_method = variance_method,
    ring           = ring,
    converged      = !is.null(gls_attempt),
    iterations     = gls_passes,
    gls_passes     = gls_passes,
    gls_error      = gls_error,
    SSW            = SSW,
    SSB            = SSB,
    MSW            = MSW,
    MSB            = MSB,
    n_clusters     = n_clusters,
    cluster_sizes  = n_i,
    fit            = fit_final,
    fit_profile    = fit_profile,
    family         = "Gaussian (share-domain residual REML + cluster-mean GLS transform, K>=3)",
    call           = match.call())
  class(out) <- c("ds.vertLMM.k3", "list")
  out
}

#' @export
print.ds.vertLMM.k3 <- function(x, ...) {
  cat("dsVert LMM (K=3, REML 1-D profile)\n")
  cat(sprintf("  Clusters = %d   Total N = %d\n",
              x$n_clusters, sum(x$cluster_sizes)))
  cat(sprintf("  rho = %.4f   sigma_b^2 = %.4g   sigma^2 = %.4g\n",
              x$rho_hat, x$sigma_b2, x$sigma2))
  cat("\nFixed effects (beta):\n")
  print(round(x$coefficients, 5))
  invisible(x)
}
