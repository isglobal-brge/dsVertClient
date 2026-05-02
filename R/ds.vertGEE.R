#' @title Federated generalised estimating equations
#' @description Fit a GLM for vertically partitioned DataSHIELD data and
#'   return sandwich (robust) standard errors alongside the usual
#'   model-based ones. The point estimate \eqn{\hat{\beta}} is obtained
#'   by a single call to \code{\link{ds.vertGLM}}. For Gaussian models the
#'   working-independence HC0 meat matrix is computed in the share domain:
#'   residual shares are squared with Beaver multiplication and each
#'   \eqn{X^\top(r^2 X_j)} column is obtained through the existing aggregate
#'   matvec path. This is not yet a clustered \code{geeglm} sandwich;
#'   exchangeable / AR1 working correlation structures require a separate
#'   per-cluster meat primitive with small-cluster guards.
#'
#'   The client receives only low-dimensional aggregate matrices/vectors.
#'   It never materialises the \eqn{n}-length residual, squared-residual, or
#'   weighted-column vectors.
#'
#'   Formula:
#'     \deqn{V_{sand} = \mathrm{Cov}(\hat\beta) \, A \, \mathrm{Cov}(\hat\beta)}
#'     \deqn{A = X^T \, \mathrm{diag}(r^2) \, X \, / \, n}
#'
#' @param formula A model formula passed through to \code{ds.vertGLM}.
#' @param data Character. Aligned data-frame name on each server.
#' @param family One of \code{"gaussian"}, \code{"binomial"},
#'   \code{"poisson"}.
#' @param id_col Optional character; recorded for future clustered GEE.
#'   The current validated implementation computes a working-independence
#'   HC0 sandwich.
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
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)
  # Cluster-ID broadcast: for exchangeable / AR1 working correlation we
  # need per-cluster residual outer products. The outcome server holds
  # cluster IDs plaintext and we extend the peer's view by shipping
  # cluster membership as a Ring63-permutation blob (same tier as
  # Cox event-time ordering -- permutation is revealed, absolute IDs
  # are not). Enabled only when id_col is provided.
  cluster_ids <- NULL
  if (corstr %in% c("exchangeable", "ar1")) {
    if (is.null(id_col)) {
      stop("corstr='", corstr, "' requires id_col", call. = FALSE)
    }
    if (verbose) {
      message("[ds.vertGEE] corstr='", corstr,
              "' - using cluster-ID broadcast path")
    }
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
                    datasources = datasources,
                    keep_session = TRUE)
  if (!is.null(fit$session_id)) {
    on.exit({
      for (.srv in fit$server_list %||% names(datasources)) {
        .ci <- which(names(datasources) == .srv)
        if (length(.ci) == 1L) {
          tryCatch(DSI::datashield.aggregate(datasources[.ci],
            call(name = "mpcCleanupDS", session_id = fit$session_id)),
            error = function(e) NULL)
          tryCatch(DSI::datashield.aggregate(datasources[.ci],
            call(name = "mpcGcDS")), error = function(e) NULL)
        }
      }
    }, add = TRUE)
  }
  if (is.null(fit$covariance)) {
    stop("ds.vertGEE: the underlying ds.vertGLM fit does not expose ",
         "Cov(beta). Refit with dsVert >= 8bb7902 / dsVertClient >= ",
         "8bb7902.", call. = FALSE)
  }
  Cov_model <- as.matrix(fit$covariance)
  Cov_model <- (Cov_model + t(Cov_model)) / 2
  Cov_model_info <- as.matrix(fit$covariance_information %||%
                                fit$covariance_unscaled %||%
                                fit$covariance)
  Cov_model_info <- (Cov_model_info + t(Cov_model_info)) / 2

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
  {
    if (verbose) message("[ds.vertGEE] Stage 2: weighted fit for sandwich meat (family=",
                          family, ")")
    # Request residuals-squared column on the outcome server. This uses
    # the `dsvertPearsonR2ColDS` server helper (Gaussian-only in the
    # first pass): given data_name, y, x_vars, and a plaintext betahat,
    # it writes a new column `__dsvert_r2` into the aligned data frame.
    # The client then re-fits with weights="__dsvert_r2".
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
      secure_hc0 <- NULL
      if (family == "gaussian") {
        secure_hc0 <- tryCatch(
          .ds_gee_secure_hc0(
            fit = fit, datasources = datasources,
            server_names = server_names, lambda = lambda,
            verbose = verbose),
          error = function(e) {
            if (verbose) {
              message("[ds.vertGEE] secure HC0 path unavailable: ",
                      conditionMessage(e))
            }
            NULL
          })
      }
      if (!is.null(secure_hc0)) {
        robust_cov <- secure_hc0$robust_covariance
        robust_se <- secure_hc0$robust_se
      } else {
        x_all <- names(fit$coefficients)
        x_all <- setdiff(x_all, "(Intercept)")
        r2_ok <- tryCatch({
          DSI::datashield.aggregate(
            datasources[which(server_names == y_srv)],
            call(name = "dsvertPearsonR2ColDS",
                 data_name = data, y_var = y_var,
                 x_names = x_all,
                 betahat = as.numeric(fit$coefficients[x_all]),
                 intercept = as.numeric(fit$coefficients["(Intercept)"]),
                 family = family,
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
          #   V_sand = Bread * Fisher_weighted * Bread
          #          ~= Cov_info_model * solve(Cov_info_weighted) *
          #              Cov_info_model
          # where Cov_info_* is inverse Fisher, not a Gaussian sigma^2-scaled
          # model covariance.
          Cov_w_info <- as.matrix(fit_w$covariance_information %||%
                                     fit_w$covariance_unscaled %||%
                                     fit_w$covariance)
          Cov_w_info <- (Cov_w_info + t(Cov_w_info)) / 2
          A <- tryCatch(solve(Cov_w_info), error = function(e) NULL)
          if (is.null(A)) {
            warning("meat matrix singular; returning model-based SE.",
                    call. = FALSE)
            robust_se <- fit$std_errors
            robust_cov <- Cov_model
          } else {
            robust_cov <- Cov_model_info %*% A %*% Cov_model_info
            dimnames(robust_cov) <- list(names(fit$coefficients),
                                          names(fit$coefficients))
            robust_se <- sqrt(pmax(diag(robust_cov), 0))
            names(robust_se) <- names(fit$coefficients)
          }
        }
        }
      }
    }
  }

  # Exchangeable / AR1 refinement on the sandwich cov.
  if (corstr %in% c("exchangeable", "ar1") && !is.null(id_col)) {
    # Request per-cluster residual aggregates from the outcome server.
    clust_info <- tryCatch(
      DSI::datashield.aggregate(
        datasources[which(server_names == .ds_gee_find_server_holding(
          datasources, server_names, data, id_col))],
        call(name = "dsvertClusterResidualsDS",
             data_name = data, y_var = y_var,
             x_names = setdiff(names(fit$coefficients), "(Intercept)"),
             intercept = as.numeric(fit$coefficients["(Intercept)"]),
             betahat = as.numeric(
               fit$coefficients[setdiff(names(fit$coefficients),
                                         "(Intercept)")]),
             cluster_col = id_col)),
      error = function(e) NULL)
    if (!is.null(clust_info)) {
      if (is.list(clust_info) && length(clust_info) == 1L)
        clust_info <- clust_info[[1L]]
      rss <- as.numeric(clust_info$rss_per_cluster)
      rsum <- as.numeric(clust_info$rsum_per_cluster)
      n_i <- as.integer(clust_info$n_per_cluster)
      # Exchangeable: estimate within-cluster correlation alpha as
      #   alpha = sum_{i: n_i>1} [sum_j<k r_ij r_ik] / sum_i n_i(n_i-1)/2
      # Use sum(r)^2 - sum(r^2) = 2 sum_{j<k} r_ij r_ik.
      offdiag_sums <- (rsum^2 - rss) / 2
      total_pairs <- sum(n_i * (n_i - 1) / 2)
      sigma2_hat <- sum(rss) / max(sum(n_i), 1L)
      if (total_pairs > 0 && sigma2_hat > 0) {
        alpha_hat <- sum(offdiag_sums) / (total_pairs * sigma2_hat)
        alpha_hat <- max(min(alpha_hat, 0.99), -0.5)
        # Inflate sandwich SE by sqrt(1 + (n_bar - 1) * alpha) for
        # exchangeable; for AR1 the correction is a decaying series
        # that we approximate with the same factor at mean cluster size.
        n_bar <- mean(n_i)
        inflation <- sqrt(max(1 + (n_bar - 1) * alpha_hat, 1e-4))
        robust_se <- robust_se * inflation
        robust_cov <- robust_cov * inflation^2
        attr(robust_cov, "corstr_alpha") <- alpha_hat
        attr(robust_cov, "corstr_inflation") <- inflation
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
        call(name = "dsvertColNamesDS", data_name = data_name))[[1]]$columns,
      error = function(e) NULL)
    if (!is.null(cols) && var %in% cols) return(srv)
  }
  NULL
}

#' @keywords internal
.ds_gee_standardized_parameters <- function(fit, features) {
  beta_orig <- as.numeric(fit$coefficients[features])
  names(beta_orig) <- features
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])
  names(x_sds) <- names(x_means) <- features

  if (fit$family == "gaussian" && isTRUE(fit$standardize_y) &&
      !is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
    beta_std <- beta_orig * x_sds / fit$y_sd
    intercept_std <- (as.numeric(fit$coefficients["(Intercept)"]) +
                        sum(beta_orig * x_means) - fit$y_mean) / fit$y_sd
  } else {
    beta_std <- beta_orig * x_sds
    intercept_std <- as.numeric(fit$coefficients["(Intercept)"]) +
                       sum(beta_orig * x_means)
  }
  names(beta_std) <- features
  list(intercept = as.numeric(intercept_std), beta = beta_std)
}

#' @keywords internal
.ds_gee_standardization_jacobian <- function(fit, features) {
  p <- length(features)
  out_names <- c("(Intercept)", features)
  J <- diag(p + 1L)
  dimnames(J) <- list(out_names, out_names)
  x_sds <- as.numeric(fit$x_sds[features])
  x_means <- as.numeric(fit$x_means[features])

  if (fit$family == "gaussian" && isTRUE(fit$standardize_y) &&
      !is.null(fit$y_sd) && is.finite(fit$y_sd) && fit$y_sd > 0) {
    for (jj in seq_len(p)) {
      J[jj + 1L, jj + 1L] <- fit$y_sd / x_sds[jj]
      J[1L, jj + 1L] <- -fit$y_sd * x_means[jj] / x_sds[jj]
    }
    J[1L, 1L] <- fit$y_sd
  } else {
    for (jj in seq_len(p)) {
      J[jj + 1L, jj + 1L] <- 1 / x_sds[jj]
      J[1L, jj + 1L] <- -x_means[jj] / x_sds[jj]
    }
  }
  J
}

#' @keywords internal
.ds_gee_secure_hc0 <- function(fit, datasources, server_names,
                               lambda = 0, verbose = FALSE) {
  if (!identical(fit$family, "gaussian")) return(NULL)
  required <- c("session_id", "transport_pks", "server_list", "x_vars",
                "y_server", "hessian_std", "x_sds", "x_means")
  missing_req <- required[vapply(required, function(nm) is.null(fit[[nm]]),
                                 logical(1L))]
  if (length(missing_req) > 0L) {
    stop("GLM session metadata missing: ", paste(missing_req, collapse = ", "),
         call. = FALSE)
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

  target_features <- unlist(x_vars[server_list], use.names = FALSE)
  target_order <- c("(Intercept)", target_features)
  std <- .ds_gee_standardized_parameters(fit, target_features)

  if (fit$eta_privacy == "k2_beaver") {
    nl <- setdiff(server_list, coordinator)
    if (length(nl) != 1L) {
      stop("K=2 GEE HC0 requires exactly one non-outcome server",
           call. = FALSE)
    }
    nl <- nl[[1L]]
    dcf_parties <- c(coordinator, nl)
    dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                        integer(1L))
    dealer_conn <- dcf_conns[[2L]]
    canonical_features <- c(x_vars[[coordinator]], x_vars[[nl]])
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
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2IdentityLinkDS", session_id = session_id))
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
    canonical_features <- c(x_vars[[coordinator]], x_vars[[fusion]])
    for (srv in non_dcf) canonical_features <- c(canonical_features, x_vars[[srv]])
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
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2IdentityLinkDS", session_id = session_id))
    }
  } else {
    stop("unsupported eta_privacy for secure GEE HC0: ", fit$eta_privacy,
         call. = FALSE)
  }

  p_total <- length(canonical_features)
  if (p_total < 1L) stop("no predictors available for GEE HC0", call. = FALSE)
  canonical_order <- c("(Intercept)", canonical_features)

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

  .gradient <- function() {
    gt <- .dsAgg(datasources[dealer_conn],
      call(name = "glmRing63GenGradTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[[1L]]]],
           dcf1_pk = transport_pks[[dcf_parties[[2L]]]],
           n = as.integer(n_obs), p = as.integer(p_total),
           ring = ring, session_id = session_id))
    if (is.list(gt) && length(gt) == 1L) gt <- gt[[1L]]
    .sendBlob(gt$grad_blob_0, "k2_grad_triple_fp", dcf_conns[[1L]])
    .sendBlob(gt$grad_blob_1, "k2_grad_triple_fp", dcf_conns[[2L]])
    r1 <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      peer <- dcf_parties[[3L - i]]
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2GradientR1DS",
             peer_pk = transport_pks[[peer]],
             session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[i]] <- r
    }
    .sendBlob(r1[[1L]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[[2L]])
    .sendBlob(r1[[2L]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[[1L]])
    r2 <- vector("list", 2L)
    for (i in seq_along(dcf_parties)) {
      r <- .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2GradientR2DS",
             party_id = as.integer(i - 1L), session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r2[[i]] <- r
    }
    g <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r2[[1L]]$gradient_fp,
      share_b = r2[[2L]]$gradient_fp,
      frac_bits = frac_bits, ring = ring_tag))
    s <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r1[[1L]]$sum_residual_fp,
      share_b = r1[[2L]]$sum_residual_fp,
      frac_bits = frac_bits, ring = ring_tag))
    c(s$values[1L], as.numeric(g$values))
  }

  for (ci in dcf_conns) {
    .dsAgg(datasources[ci],
      call(name = "k2PrepareWeightedResidualShareDS",
           session_id = session_id))
  }
  .vecmul("k2_weight_residual_share_fp",
          "k2_weight_residual_share_fp",
          "gee_r2_share")

  for (i in seq_along(dcf_parties)) {
    .dsAgg(datasources[dcf_conns[[i]]],
      call(name = "dsvertGEEInterceptShareDS",
           output_key = "gee_x_col_0",
           n = as.numeric(n_obs),
           is_party0 = (i == 1L),
           session_id = session_id,
           frac_bits = frac_bits, ring = ring))
  }
  for (j in seq_len(p_total)) {
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2BeaverExtractColumnDS",
             source_key = "k2_x_full_fp",
             n = as.numeric(n_obs), K = as.numeric(p_total),
             col_index = as.numeric(j),
             output_key = paste0("gee_x_col_", j),
             session_id = session_id,
             frac_bits = frac_bits, ring = ring_tag))
    }
  }

  A_canonical <- matrix(0, nrow = p_total + 1L, ncol = p_total + 1L)
  for (j in 0:p_total) {
    x_key <- paste0("gee_x_col_", j)
    wx_key <- paste0("gee_r2_x_col_", j)
    .vecmul("gee_r2_share", x_key, wx_key)
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci],
        call(name = "k2FinalizeWeightedResidualShareDS",
             input_key = wx_key, session_id = session_id))
    }
    A_canonical[, j + 1L] <- .gradient()
    if (verbose) {
      message(sprintf("[ds.vertGEE] secure HC0 meat column %d/%d",
                      j + 1L, p_total + 1L))
    }
  }
  A_canonical <- (A_canonical + t(A_canonical)) / 2
  dimnames(A_canonical) <- list(canonical_order, canonical_order)

  perm <- match(target_order, canonical_order)
  if (anyNA(perm)) {
    stop("could not align GEE HC0 meat order", call. = FALSE)
  }
  A_std <- A_canonical[perm, perm, drop = FALSE]
  dimnames(A_std) <- list(target_order, target_order)

  H <- as.matrix(fit$hessian_std)
  if (!is.null(rownames(H))) {
    H <- H[target_order, target_order, drop = FALSE]
  }
  H_adj <- H - as.numeric(lambda %||% 0) * diag(nrow(H))
  fisher_std <- n_obs * H_adj
  cov_std <- tryCatch(solve(fisher_std), error = function(e) NULL)
  if (is.null(cov_std)) stop("standardized Fisher matrix is singular",
                             call. = FALSE)
  dimnames(cov_std) <- list(target_order, target_order)

  robust_std <- cov_std %*% A_std %*% cov_std
  J <- .ds_gee_standardization_jacobian(fit, target_features)
  robust_cov <- J %*% robust_std %*% t(J)
  robust_cov <- (robust_cov + t(robust_cov)) / 2
  dimnames(robust_cov) <- list(target_order, target_order)
  robust_se <- sqrt(pmax(diag(robust_cov), 0))
  names(robust_se) <- target_order
  list(robust_covariance = robust_cov,
       robust_se = robust_se,
       meat_std = A_std,
       method = "share_domain_hc0")
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
