#' @title Signed binary random-effect GLMM
#' @description Fits a binary random-effect model from one signed, sticky
#'   Synopsis.  The intercept-only call retains its historical moment
#'   projection. Additive covariate calls select the minimum from the signed
#'   finite marginal-likelihood grid; they never call PQL or expose
#'   cluster-level statistics.
#' @details The signed outcome bounds must be exactly \code{[0, 1]}. For an
#'   intercept-only formula, \code{sigma_b2} is the conventional logistic
#'   latent-scale approximation to the released observed ICC. For additive
#'   covariates, it is the selected value from a custodian-signed finite random
#'   intercept variance grid or a two-effect covariance grid. Neither route
#'   supplies standard errors, p-values or sampling inference; a named random slope is available only when it
#'   exactly matches a signed finite two-effect grid. Interactions and
#'   unconstrained likelihood optimisation remain unavailable.
#' @param formula An intercept-only formula or additive bare column names.
#' @param data Signed protected dataset name or federation.
#' @param cluster_col Cluster column required to match the signed artifact.
#' @param analysis_id Custodian-configured signed random-intercept artifact id.
#' @param random_slopes Optional bare predictor name for a signed finite-grid
#'   random-slope artifact; it must match the artifact exactly.
#' @param max_outer,inner_iter,tol,ring,verbose Retained compatibility controls;
#'   they do not alter the signed estimand.
#' @param lambda Must be zero.
#' @param compute_se Must be \code{FALSE}.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertGLMM} object containing the certified public DP
#'   moment or finite-grid projection and no cluster-level statistics.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertGLMM <- function(formula, data = NULL, cluster_col,
                        analysis_id = NULL,
                        random_slopes = NULL,
                        max_outer = 30L, inner_iter = 50L,
                        tol = 1e-4, lambda = 0,
                        compute_se = FALSE,
                        ring = NULL,
                        verbose = TRUE,
                        datasources = NULL) {
  terms <- if (inherits(formula, "formula")) stats::terms(formula) else NULL
  predictors <- if (is.null(terms)) character() else
    attr(terms, "term.labels")
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]]) ||
      !identical(attr(terms, "intercept"), 1L) ||
      any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", predictors))) {
    stop("ds.vertGLMM requires an intercept and additive bare column names",
         call. = FALSE)
  }
  if (!is.character(cluster_col) || length(cluster_col) != 1L ||
      !nzchar(cluster_col) || !is.character(analysis_id) ||
      length(analysis_id) != 1L || !nzchar(analysis_id)) {
    stop("ds.vertGLMM requires cluster_col and signed analysis_id strings",
         call. = FALSE)
  }
  if (!is.null(random_slopes) && (!is.character(random_slopes) ||
      length(random_slopes) != 1L || is.na(random_slopes) ||
      !grepl("^[A-Za-z.][A-Za-z0-9._]*$", random_slopes))) {
    stop("random_slopes must be one bare signed predictor name or NULL",
         call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) ||
      !is.finite(lambda) || lambda != 0 || !identical(compute_se, FALSE)) {
    stop(paste(
      "ds.vertGLMM supports only signed binary random-intercept routes:",
      "lambda=0 and compute_se=FALSE"), call. = FALSE)
  }
  if (length(predictors)) {
    resolved <- .dsvert_federation_argument(data, datasources)
    result <- .dsvert_dp_glmm_grid_impl(
      formula, resolved$value, cluster_col, analysis_id,
      datasources = resolved$datasources, .aggregate = DSI::datashield.aggregate)
    signed_slopes <- result$signed_artifact$random_effect_order[-1L] %||% character()
    if (!identical(signed_slopes, random_slopes %||% character())) {
      stop("random_slopes must exactly match the signed GLMM artifact", call. = FALSE)
    }
    return(result)
  }
  if (!is.null(random_slopes)) {
    stop("random_slopes requires a signed finite-grid covariate formula", call. = FALSE)
  }
  signed <- ds.vertDPLMM(
    data_name = data, analysis_id = analysis_id, datasources = datasources)
  artifact <- signed$signed_artifact
  outcome <- as.character(formula[[2L]])
  if (!is.list(artifact) || !is.list(artifact$outcome) ||
      !is.list(artifact$cluster) ||
      !identical(artifact$outcome$column, outcome) ||
      !identical(artifact$cluster$column, cluster_col) ||
      !isTRUE(all.equal(as.numeric(artifact$outcome$lower), 0)) ||
      !isTRUE(all.equal(as.numeric(artifact$outcome$upper), 1))) {
    stop(paste(
      "formula, cluster_col and binary [0, 1] outcome bounds must match",
      "the signed random-intercept artifact"), call. = FALSE)
  }
  if (!identical(signed$status, "ok")) {
    signed$family <- "binomial"
    signed$estimand <- "binary_random_intercept_moment_non_identifiable"
    signed$coefficients <- signed$coefficient <- NULL
    signed$sigma2 <- signed$sigma_b2 <- signed$icc <- NULL
    signed$standard_errors <- NULL
    signed$legacy_fallback_called <- FALSE
    class(signed) <- unique(c("ds.vertGLMM", class(signed)))
    return(signed)
  }
  probability <- as.numeric(signed$coefficients[["(Intercept)"]])
  n_obs <- as.numeric(signed$n_obs)
  rho <- as.numeric(signed$icc)
  if (!is.finite(probability) || !is.finite(n_obs) || n_obs < 2 ||
      !is.finite(rho)) {
    stop("The signed binary random-intercept projection is invalid",
         call. = FALSE)
  }
  floor_probability <- 1 / (2 * n_obs)
  marginal_probability <- min(
    max(probability, floor_probability), 1 - floor_probability)
  observed_icc <- min(max(rho, 0), 1 - 1 / n_obs)
  latent_sigma_b2 <- if (observed_icc == 0) 0 else {
    pi^2 * observed_icc / (3 * (1 - observed_icc))
  }
  signed$family <- "binomial"
  signed$estimand <-
    "binary_random_intercept_population_average_moment_approximation"
  signed$coefficients <- signed$coefficient <- stats::setNames(
    stats::qlogis(marginal_probability), "(Intercept)")
  signed$marginal_probability <- marginal_probability
  signed$probability_projection_applied <-
    !isTRUE(all.equal(probability, marginal_probability, tolerance = 1e-12))
  signed$icc <- signed$icc_observed <- observed_icc
  signed$icc_scale <- "observed_pair_correlation"
  signed$sigma2 <- NULL
  signed$sigma_b2 <- signed$latent_sigma_b2_approx <- latent_sigma_b2
  signed$random_effect_scale <- "latent_logit_approximation"
  signed$standard_errors <- signed$p_values <- NULL
  signed$cluster_sizes <- NULL
  signed$iterations <- 0L
  signed$converged <- TRUE
  signed$legacy_fallback_called <- FALSE
  class(signed) <- unique(c("ds.vertGLMM", class(signed)))
  signed
}

#' @export
print.ds.vertGLMM <- function(x, ...) {
  if (!identical(x$status, "ok")) {
    cat("dsVert signed binary random-intercept GLMM\n")
    cat("  Status: ", x$status %||% "non_identifiable", "\n", sep = "")
    return(invisible(x))
  }
  if (!is.null(x$selected_candidate)) {
    cat("dsVert signed binary finite-grid GLMM\n")
    cat(sprintf("  Selected signed candidate = %d\n", x$selected_candidate))
    cat(sprintf("  Random-intercept variance = %.5g\n", x$sigma_b2))
    cat(sprintf("  DP negative log likelihood = %.5f\n",
                x$selected_dp_negative_log_likelihood))
    cat("  Coefficients:\n")
    print(round(x$coefficients, 5L))
    if (!is.null(x$random_effect_covariance)) {
      cat("  Random-effect covariance:\n")
      print(round(x$random_effect_covariance, 5L))
    }
    cat("  Standard errors and sampling inference are unavailable.\n")
    return(invisible(x))
  }
  cat("dsVert signed binary random-intercept moment approximation\n")
  cat(sprintf("  Population-average event probability = %.5f\n",
              x$marginal_probability))
  cat(sprintf("  Observed ICC = %.4f    latent-scale sigma_b^2 = %.4g\n",
              x$icc_observed, x$latent_sigma_b2_approx))
  cat("  Intercept (population-average log-odds):\n")
  print(round(x$coefficients, 5L))
  cat("  Standard errors and sampling inference are unavailable.\n")
  invisible(x)
}

#' @keywords internal
.ds_glmm_cleanup_fit_session <- function(fit, datasources) {
  if (is.null(fit$session_id) || is.null(fit$server_list)) return(invisible(NULL))
  cleanup_sites <- intersect(fit$server_list, names(datasources))
  cleanup_conns <- datasources[cleanup_sites]
  .dsvert_cleanup_best_effort(
    cleanup_conns,
    call(name = "mpcCleanupDS", session_id = fit$session_id))
  .dsvert_cleanup_best_effort(cleanup_conns, call(name = "mpcGcDS"))
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
.ds_glmm_share_domain_moments <- function(fit, data, cluster_col,
                                          datasources, server_names,
                                          verbose = FALSE,
                                          laplace_components = FALSE) {
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
    .dsvert_aggregate_strict(
      conns, expr, operation = "protected GLMM MPC phase")
  }
  .sendBlob <- function(blob, contract, conn_idx) {
    .dsvert_store_transfer_or_legacy(
      blob, contract, datasources[conn_idx], session_id,
      producer_conns = datasources,
      .aggregate = .dsAgg)
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

  # REVEAL-FREE link (F1 fix): k2_eta_share_fp -> secure_mu_share on ring127 shares
  # (.ring127_sigmoid_round_keyed). dcf_masked relayed 0x; dealer-free.
  do.call(.ring127_sigmoid_round_keyed, list(
    in_key = "k2_eta_share_fp", out_key = "secure_mu_share",
    n = n_obs, datasources = datasources, dealer_ci = dcf_conns[[2L]],
    server_list = dcf_parties, server_names = server_names,
    y_server = dcf_parties[[1L]], nl = dcf_parties[[2L]],
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob))

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
    if (identical(as.integer(ring), 127L)) {
      .ring127_vecmul(
        x_key, y_key, output_key, n_obs,
        datasources, dcf_conns[[2L]], dcf_parties, server_names,
        dcf_parties[[1L]], dcf_parties[[2L]], transport_pks, session_id,
        .dsAgg, .sendBlob)
      return(invisible(output_key))
    }
    .ot_beaver_prepare_vecmul(
      datasources = datasources,
      party_conns = dcf_conns,
      party_names = dcf_parties,
      transport_pks = transport_pks,
      session_id = session_id,
      n = n_obs,
      ring = ring,
      .dsAgg = .dsAgg,
      .sendBlob = .sendBlob)
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
    .sendBlob(r1[[1L]]$peer_blob, r1[[1L]]$peer_transfer,
              dcf_conns[[2L]])
    .sendBlob(r1[[2L]]$peer_blob, r1[[2L]]$peer_transfer,
              dcf_conns[[1L]])
    for (i in seq_along(dcf_parties)) {
      .dsAgg(datasources[dcf_conns[[i]]],
        call(name = "k2BeaverVecmulR2DS",
             is_party0 = (i == 1L),
             peer_name = dcf_parties[[3L - i]],
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

  laplace <- NULL
  if (isTRUE(laplace_components)) {
    if (ring != 127L) {
      stop("GLMM Laplace aggregate components require ring=127",
           call. = FALSE)
    }
    softplus_intervals <- as.integer(getOption(
      "dsvert.glmm_softplus_intervals", 80L))
    if (!is.finite(softplus_intervals) || softplus_intervals < 10L) {
      softplus_intervals <- 80L
    }

    # REVEAL-FREE link (F1 fix): k2_eta_share_fp -> softplus_share_fp on ring127 shares
    # (.ring127_softplus_round_keyed). dcf_masked relayed 0x; dealer-free.
    do.call(.ring127_softplus_round_keyed, list(
      in_key = "k2_eta_share_fp", out_key = "softplus_share_fp",
      n = n_obs, datasources = datasources, dealer_ci = dcf_conns[[2L]],
      server_list = dcf_parties, server_names = server_names,
      y_server = dcf_parties[[1L]], nl = dcf_parties[[2L]],
      transport_pks = transport_pks, session_id = session_id,
      .dsAgg = .dsAgg, .sendBlob = .sendBlob))

    .vecmul("k2_y_share_fp_original", "k2_eta_share_fp",
            "glmm_laplace_yeta_share")
    laplace <- list(
      sum_softplus = .sum_share("softplus_share_fp"),
      y_dot_eta = .sum_share("glmm_laplace_yeta_share"),
      s_by_cluster = vs$values,
      rsum_per_cluster = rsum,
      cluster_sizes = py$sizes,
      n = sum(py$sizes))
  }

  list(rsum_per_cluster = rsum,
       vsum_per_cluster = vs$values,
       n_per_cluster = py$sizes,
       laplace_components = laplace)
}
