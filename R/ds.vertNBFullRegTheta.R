#' @title Quarantined negative-binomial compatibility frontdoor
#' @description This exported name is retained for API compatibility. It
#'   raises a typed \code{dsvert_route_unavailable} condition before any DSI
#'   call and computes or returns no fit or dispersion estimate. Retained
#'   implementation code after the gate is unreachable through this public
#'   frontdoor and carries no disclosure, DP, accuracy, or availability claim.
#' @details Promotion requires a bounded NB2 score/information capsule, joint
#'   beta/theta inference, certified deviance and covariance, and independent
#'   multi-host validation.
#' @param formula,data,theta,joint,theta_max_iter,theta_tol,variant,beta_max_iter,beta_tol,compute_covariance,verbose,datasources
#'   Retained compatibility arguments. They are not evaluated because the
#'   public frontdoor fails locally.
#' @param ... Retained compatibility arguments; not evaluated.
#' @return No fitted object. The function raises
#'   \code{dsvert_route_unavailable} before DSI.
#' @seealso \code{\link{ds.vertMethodStatus}}
#' @export
ds.vertNBFullRegTheta <- function(formula, data = NULL, theta = NULL,
                                  joint = TRUE, theta_max_iter = 5L,
                                  theta_tol = 1e-3, variant = "full_reg_nd",
                                  beta_max_iter = 2L, beta_tol = 1e-4,
                                  compute_covariance = TRUE,
                                  verbose = TRUE, datasources = NULL, ...) {
  .dsvert_block_retired_remote_route("negative_binomial")
  if (!identical(variant, "full_reg_nd")) {
    stop("variant must be 'full_reg_nd'", call. = FALSE)
  }

  base_fit <- ds.vertNB(formula = formula, data = data, theta = theta,
                        joint = joint, theta_max_iter = theta_max_iter,
                        theta_tol = theta_tol, verbose = verbose,
                        datasources = datasources, ...)

  theta_iid <- base_fit$theta
  y_mean <- base_fit$y_mean
  y_var <- base_fit$y_var
  n <- base_fit$n_obs

  # ============================================================
  # Non-disclosive full-regression theta MLE -- share-domain pipeline
  # ============================================================
  # Closes D-INV-4 (per-patient eta^nl reveal at label, present in the
  # legacy "full_reg" variant). eta^nl stays in Ring127 additive secret
  # shares end-to-end through mu = exp(eta)_share, log(mu+theta)_share,
  # 1/(theta+mu)_share and (y+theta)*1/(theta+mu)_share via Beaver vecmul +
  # AffineCombine + Chebyshev-Clenshaw primitives. Only the four
  # scalar aggregates Sumlog(mu+theta), Sum1/(theta+mu), Sum(y+theta)/(theta+mu),
  # Sum(y+theta)/(theta+mu)^2 (plus label-only Sumpsi(y+theta), Sumpsi_1(y+theta)) are revealed
  # per Newton-theta iter.
  #
  # Refs: Lawless 1987 *Can. J. Statist.* 15(3):209-225 (NB profile-MLE
  # theta score); Venables-Ripley 2002 *MASS* Sec.7.4 (\code{glm.nb} Newton);
  # Catrina-Saxena 2010 *Financial Cryptography* Sec.3.3 (multiplicative-
  # depth ULP); Beaver 1991 *CRYPTO* Sec.3 (precomputed multiplication
  # triples); Demmler-Schneider-Zohner ABY 2015 Sec.III.B (K=2 OT-Beaver
  # dishonest-majority); Trefethen ATAP Sec.8 (Bernstein-ellipse Cheb).
  if (identical(variant, "full_reg_nd")) {
    if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
    server_names <- names(datasources)
    y_var_char <- .ds_gee_extract_lhs(formula)
    y_srv <- .ds_gee_find_server_holding(datasources, server_names, data, y_var_char)
    if (is.null(y_srv)) stop("outcome server not found", call. = FALSE)
    non_y_srv <- setdiff(server_names, y_srv)
    if (length(non_y_srv) < 1L)
      stop("full_reg_nd variant requires at least two servers", call. = FALSE)
    y_ci  <- which(server_names == y_srv)

    # Identify which features live on each server.
    rhs <- attr(stats::terms(formula), "term.labels")
    column_results <- .dsvert_aggregate_strict(
      datasources,
      call(name = "dsvertColNamesDS", data_name = data),
      operation = "negative-binomial column discovery")
    cols_by_server <- lapply(column_results, function(value) {
      if (is.list(value)) value$columns else value
    })
    x_vars_full <- lapply(server_names, function(srv)
      intersect(rhs, cols_by_server[[srv]]))
    names(x_vars_full) <- server_names
    x_label <- x_vars_full[[y_srv]]

    # beta-slice per server from the Poisson fit's revealed beta.
    beta_all <- base_fit$poisson_fit$coefficients
    int_val  <- beta_all[["(Intercept)"]]
    beta_label <- beta_all[x_label]

    session_id <- .dsvert_uuid4()

    # Closures to abstract the DataSHIELD aggregate / blob-relay calls
    # so the orchestrator can be invoked under both real-Opal and
    # local-harness test fixtures (mirrors ord_joint / mnl_joint .dsAgg
    # / .sendBlob plumbing).
    .dsAgg <- function(conns, expr) {
      .dsvert_aggregate_strict(
        conns, expr,
        operation = "negative-binomial protected MPC phase")
    }
    .sendBlob <- function(blob, contract, target_ci) {
      .dsvert_store_transfer_or_legacy(
        blob, contract, datasources[target_ci], session_id,
        producer_conns = datasources)
    }

    # Use the generic Ring127 input-sharing setup for K=2 and K>=3. The
    # older K=2-only eta-share setup was non-disclosive for theta, but it did
    # not keep X/y shares needed for the NB beta score.
    setup <- .nb_fullreg_nd_session_setup_k3(
      data = data,
      x_vars = x_vars_full,
      beta_all = beta_all,
      y_var_char = y_var_char,
      y_srv = y_srv,
      datasources = datasources,
      server_names = server_names,
      server_list = server_names,
      session_id = session_id,
      .dsAgg = .dsAgg,
      .sendBlob = .sendBlob,
      verbose = verbose)
    score_server_list <- setup$dcf_parties
    score_nl_srv <- setup$fusion_server
    score_nl_ci <- which(server_names == score_nl_srv)
    dealer_ci <- which(server_names == setup$dealer_servers[[1L]])
    yt_from_y_share <- TRUE
    n_obs        <- setup$n
    transport_pks <- setup$transport_pks
    p_total <- setup$p_total
    coef_order <- c("(Intercept)", setup$feature_order)
    beta_current <- base_fit$poisson_fit$coefficients
    beta_current <- beta_current[names(base_fit$poisson_fit$coefficients)]

    # Keep the existing helper index for DCF/key orchestration. Beaver
    # preprocessing is negotiated separately by the shared helper.

    # Newton-theta loop on the share-domain score. MoM warm-init from the
    # Poisson fit's residual moments (Venables-Ripley Sec.7.4 default seed).
    theta_iid <- base_fit$theta
    y_mean <- base_fit$y_mean
    y_var <- base_fit$y_var
    theta_mom <- if (is.finite(y_var) && y_var > y_mean + 1e-10)
      max(y_mean^2 / max(y_var - y_mean, 1e-6), 0.1) else NA_real_
    # The y-only iid profile systematically under-seeds regression-rich NB
    # fixtures because it collapses mu_i variation into ybar. A doubled iid
    # seed stays aggregate-only and cuts the Ring127 full-score evaluations
    # from ~5 to ~3 on the current validation fixtures.
    theta_regression_seed <- if (is.finite(theta_iid) && theta_iid > 0)
      2 * theta_iid else NA_real_
    theta_seed <- max(c(theta_mom, theta_iid, theta_regression_seed, 1e-3),
                      na.rm = TRUE)
    if (!is.finite(theta_seed) || theta_seed <= 0) theta_seed <- 1e-3
    theta_cur <- theta_seed

    score_eval <- function(th) {
      tryCatch(.nb_fullreg_nd_score(
        theta = th, n_obs = n_obs,
        datasources = datasources, dealer_ci = dealer_ci,
        server_list = score_server_list, server_names = server_names,
        y_server = y_srv, nl = score_nl_srv,
        ci_os = y_ci, ci_nl = score_nl_ci,
        transport_pks = transport_pks, session_id = session_id,
        .dsAgg = .dsAgg, .sendBlob = .sendBlob, verbose = verbose,
        yt_from_y_share = yt_from_y_share),
        error = function(e) {
          message(sprintf("[NBFullRegND] score eval ERR at theta=%.4f: %s",
                           th, conditionMessage(e)))
          list(score = NA_real_, deriv = NA_real_, n = n_obs)
        })
    }

    theta_max_iter <- max(1L, as.integer(theta_max_iter))
    theta_tol <- max(as.numeric(theta_tol), .Machine$double.eps)
    beta_max_iter <- max(0L, as.integer(beta_max_iter))
    beta_tol <- max(as.numeric(beta_tol), .Machine$double.eps)
    theta_trace <- data.frame(
      beta_iter = integer(0L),
      iter = integer(0L),
      theta = numeric(0L),
      score = numeric(0L),
      deriv = numeric(0L),
      theta_next = numeric(0L),
      stringsAsFactors = FALSE)
    beta_trace <- data.frame(
      iter = integer(0L),
      theta = numeric(0L),
      max_abs_score = numeric(0L),
      max_abs_step = numeric(0L),
      stringsAsFactors = FALSE)
    theta_converged <- FALSE
    beta_converged <- beta_max_iter == 0L
    beta_stats <- NULL

    for (b_it in seq.int(0L, beta_max_iter)) {
      if (b_it > 0L) {
        .nb_fullreg_nd_recompute_eta(
          beta_all = beta_current, setup = setup, y_srv = y_srv,
          datasources = datasources, server_names = server_names,
          session_id = session_id, .dsAgg = .dsAgg)
      }

      theta_converged_iter <- FALSE
      for (it in seq_len(theta_max_iter)) {
        s <- score_eval(theta_cur)
        if (anyNA(unlist(s[c("score","deriv")]))) break
        if (!is.finite(s$deriv) || abs(s$deriv) < 1e-12) break
        # Update on log(theta), not theta. This keeps positivity by
        # construction and is much better scaled when the fixed-mu theta is
        # several-fold above the y-only MoM/iid starting values.
        step <- s$score / (s$deriv * theta_cur)
        step <- max(min(step, log(4)), -log(4))
        theta_new <- theta_cur * exp(-step)
        damp <- 0L
        while (theta_new <= 1e-6 && damp < 20L) {
          step <- step / 2
          theta_new <- theta_cur * exp(-step)
          damp <- damp + 1L
        }
        if (!is.finite(theta_new) || theta_new <= 0) break
        theta_trace <- rbind(theta_trace, data.frame(
          beta_iter = b_it,
          iter = it,
          theta = theta_cur,
          score = as.numeric(s$score),
          deriv = as.numeric(s$deriv),
          theta_next = theta_new,
          stringsAsFactors = FALSE))
        if (abs(theta_new - theta_cur) < theta_tol * max(1, abs(theta_cur))) {
          theta_cur <- theta_new
          theta_converged <- TRUE
          theta_converged_iter <- TRUE
          break
        }
        theta_cur <- theta_new
      }
      if (!theta_converged_iter && nrow(theta_trace) == 0L) break

      if (b_it >= beta_max_iter && !isTRUE(compute_covariance)) break

      # Refresh mu and reciprocal at the final theta for this beta. The theta
      # loop evaluates the score before proposing theta_next, so the final
      # accepted value needs one explicit pass before beta score/covariance.
      s_final <- score_eval(theta_cur)
      if (anyNA(unlist(s_final[c("score","deriv")]))) break

      beta_stats <- .nb_fullreg_nd_beta_score_fisher(
        theta = theta_cur, n_obs = n_obs, p_total = p_total,
        datasources = datasources, dealer_ci = dealer_ci,
        server_list = score_server_list, server_names = server_names,
        y_server = y_srv, nl = score_nl_srv,
        transport_pks = transport_pks, session_id = session_id,
        .dsAgg = .dsAgg, .sendBlob = .sendBlob, verbose = verbose)
      names(beta_stats$score) <- coef_order
      dimnames(beta_stats$fisher) <- list(coef_order, coef_order)

      if (b_it >= beta_max_iter) break

      fisher <- beta_stats$fisher
      step <- .dsvert_solve_identifiable(
        fisher, beta_stats$score,
        context = "The protected negative-binomial beta system",
        reason = "singular_negative_binomial_information",
        symmetric = TRUE)
      if (any(!is.finite(step))) break
      max_step <- max(abs(step))
      if (max_step > 1) {
        step <- step / max_step
        max_step <- 1
      }
      beta_vec <- beta_current[coef_order] + step
      beta_current[coef_order] <- beta_vec
      beta_trace <- rbind(beta_trace, data.frame(
        iter = b_it + 1L,
        theta = theta_cur,
        max_abs_score = max(abs(beta_stats$score)),
        max_abs_step = max_step,
        stringsAsFactors = FALSE))
      if (max_step < beta_tol) {
        beta_converged <- TRUE
        break
      }
    }

    var_inflation <- if (is.finite(theta_cur) && theta_cur > 0)
      sqrt(1 + base_fit$y_mean / theta_cur) else 1
    pf <- base_fit$poisson_fit
    coef_names <- names(pf$coefficients)
    cov_beta <- NULL
    nb_se <- pf$std_errors * var_inflation
    if (isTRUE(compute_covariance) &&
        !is.null(beta_stats) && !is.null(beta_stats$fisher)) {
      cov_beta <- .dsvert_solve_identifiable(
        beta_stats$fisher,
        context = "The protected negative-binomial covariance system",
        reason = "singular_negative_binomial_information",
        symmetric = TRUE)
      cov_beta <- cov_beta[coef_names, coef_names, drop = FALSE]
      nb_se <- sqrt(pmax(diag(cov_beta), 0))
    }
    out <- base_fit
    out$coefficients <- beta_current[coef_names]
    out$theta <- theta_cur
    out$theta_iid <- theta_iid
    out$variance_correction <- NA_real_
    out$variant <- "full_reg_nd"
    out$theta_trace <- theta_trace
    out$theta_iter <- nrow(theta_trace)
    out$theta_converged <- theta_converged
    out$beta_trace <- beta_trace
	    out$beta_iter <- nrow(beta_trace)
	    out$beta_converged <- beta_converged
	    out$std_errors <- nb_se
	    out$z_values <- out$coefficients / out$std_errors
	    out$p_values <- 2 * stats::pnorm(-abs(out$z_values))
	    out$covariance <- cov_beta
	    out$var_inflation <- var_inflation
	    out$estimand <- "negative_binomial_full_regression_unpenalized"
	    out$optimizer_regularization <- "none"
	    out$quality <- .dsvert_quality_from_convergence(
	      theta_converged, metric = if (nrow(theta_trace)) {
	        abs(theta_trace$theta_next[nrow(theta_trace)] -
	              theta_trace$theta[nrow(theta_trace)])
	      } else Inf,
	      tolerance = theta_tol * max(1, abs(theta_cur)),
	      label = "NB full-regression theta")
	    if (!isTRUE(beta_converged) && beta_max_iter > 0L) {
	      out$quality$status <- "degraded"
	      out$quality$warnings <- unique(c(
	        out$quality$warnings,
	        "NB full-regression beta refinement reached beta_max_iter before convergence."))
	    }
	    out$quality$metrics$beta_converged <- beta_converged
	    out$quality$metrics$beta_iter <- nrow(beta_trace)
	    out$quality$metrics$theta_iter <- nrow(theta_trace)
	    class(out) <- c("ds.vertNBFullRegTheta", class(out))
	    return(out)
	  }

  if (identical(variant, "iid_mu")) {
    out <- base_fit
	    out$theta_iid <- theta_iid
	    out$variance_correction <- 0
	    out$variant <- "iid_mu"
	    out$quality <- .dsvert_quality(
	      status = "approximate",
	      warnings = "NB iid_mu variant keeps the outcome-only theta approximation; use method = 'accurate' for the non-disclosive full-regression route.",
	      metrics = list(theta = theta_iid, variance_correction = 0))
	    class(out) <- c("ds.vertNBFullRegTheta", class(out))
	    return(out)
	  }

  # Aggregate Var(mu) estimate from law of total variance:
  #   Var(y) = E[mu] + E[mu^2]/theta + Var(mu)  (NB conditional variance)
  # With iid-mu assumption E[mu^2] approx ybar^2:
  #   Var(mu) approx s_y^2 - ybar - ybar^2/theta_iid
  var_mu_hat <- max(0, y_var - y_mean - y_mean^2 / max(theta_iid, 1e-6))

  if (!is.finite(var_mu_hat) || var_mu_hat <= 0 || !is.finite(theta_iid) ||
      theta_iid <= 0) {
    if (isTRUE(verbose)) {
      message(sprintf("[ds.vertNBFullRegTheta] no correction (Var(mu)=%.3g, theta_iid=%.3g) -- returning iid-mu result",
                      var_mu_hat, theta_iid))
    }
    out <- base_fit
	    out$theta_iid <- theta_iid
	    out$variance_correction <- var_mu_hat
	    out$variant <- "iid_mu (fallback)"
	    out$quality <- .dsvert_quality(
	      status = "approximate",
	      warnings = "NB corrected route could not estimate a positive aggregate Var(mu); returned the iid-mu approximation.",
	      metrics = list(theta = theta_iid, variance_correction = var_mu_hat))
	    class(out) <- c("ds.vertNBFullRegTheta", class(out))
	    return(out)
	  }

  # Variance-corrected profile score: use outcome server's scalar
  # aggregates (Sumpsi(y+theta), Sumpsi_1(y+theta), n, ybar) via dsvertNBProfileSumsDS,
  # then add the Taylor correction to Sum log(mu+theta) around ybar.
  server_names <- names(datasources %||% DSI::datashield.connections_find())
  conns <- datasources %||% DSI::datashield.connections_find()
  y_var_name <- .ds_gee_extract_lhs(formula)
  y_srv <- .ds_gee_find_server_holding(conns, server_names, data, y_var_name)
  conn_idx <- which(server_names == y_srv)

  score_corrected <- function(th) {
    if (!is.finite(th) || th <= 0) return(NA_real_)
    sums <- tryCatch({
      r <- .dsvert_aggregate_strict(
        conns[conn_idx],
        call(name = "dsvertNBProfileSumsDS",
             data_name = data, variable = y_var_name, theta = th),
        operation = "negative-binomial protected profile release")[[1L]]
      r
    }, error = function(e) NULL)
    if (is.null(sums) || !is.finite(sums$sum_psi)) return(NA_real_)
    s_iid <- sums$sum_psi - n * digamma(th) +
             n * log(th / (y_mean + th))
    # Taylor correction: d/dtheta [-1/2 n V_mu / (ybar+theta)^2] = n V_mu / (ybar+theta)^3
    s_iid + n * var_mu_hat / (y_mean + th)^3
  }

  # Brent root-find on the corrected score. Bracket around theta_iid.
  lo <- max(1e-4, theta_iid * 0.25)
  hi <- theta_iid * 4
  s_lo <- score_corrected(lo)
  s_hi <- score_corrected(hi)
  # Expand bracket if needed
  while (is.finite(s_lo) && is.finite(s_hi) && sign(s_lo) == sign(s_hi) &&
         hi < 1e6) {
    hi <- hi * 2
    s_hi <- score_corrected(hi)
  }
  if (!is.finite(s_lo) || !is.finite(s_hi) || sign(s_lo) == sign(s_hi)) {
    if (isTRUE(verbose)) {
      message("[ds.vertNBFullRegTheta] could not bracket corrected root -- returning iid-mu")
    }
    out <- base_fit
	    out$theta_iid <- theta_iid
	    out$variance_correction <- var_mu_hat
	    out$variant <- "iid_mu (no bracket)"
	    out$quality <- .dsvert_quality(
	      status = "degraded",
	      warnings = "NB corrected route could not bracket the corrected theta root; returned the iid-mu approximation.",
	      metrics = list(theta = theta_iid, variance_correction = var_mu_hat,
	                     lower_score = s_lo, upper_score = s_hi))
	    class(out) <- c("ds.vertNBFullRegTheta", class(out))
	    return(out)
	  }
  theta_corr <- tryCatch(
    stats::uniroot(score_corrected, lower = lo, upper = hi,
                    tol = 1e-6, maxiter = 40L)$root,
    error = function(e) theta_iid)

  # Rescale SE using corrected theta.
  var_inflation <- if (is.finite(theta_corr) && theta_corr > 0) {
    sqrt(1 + y_mean / theta_corr)
  } else 1
  poisson_fit <- base_fit$poisson_fit
  nb_se <- poisson_fit$std_errors * var_inflation
  nb_z <- poisson_fit$coefficients / nb_se
  nb_p <- 2 * stats::pnorm(-abs(nb_z))
  nb_cov <- if (!is.null(poisson_fit$covariance))
    poisson_fit$covariance * var_inflation^2 else NULL

  out <- base_fit
  out$theta <- theta_corr
  out$theta_iid <- theta_iid
  out$variance_correction <- var_mu_hat
  out$variant <- "corrected"
  out$std_errors <- nb_se
  out$z_values <- nb_z
	  out$p_values <- nb_p
	  out$covariance <- nb_cov
	  out$var_inflation <- var_inflation
	  out$quality <- .dsvert_quality(
	    status = "approximate",
	    warnings = "NB corrected route uses an aggregate Taylor correction to the full-regression theta score.",
	    metrics = list(theta = theta_corr, theta_iid = theta_iid,
	                   variance_correction = var_mu_hat))
	  class(out) <- c("ds.vertNBFullRegTheta", class(out))
	  out
	}

#' @export
print.ds.vertNBFullRegTheta <- function(x, ...) {
  cat("dsVert NB regression (full-reg theta: variance-corrected profile MLE)\n")
  cat(sprintf("  N = %d   theta = %.4g   theta_iid = %.4g   variant = %s\n",
              x$n_obs, x$theta, x$theta_iid, x$variant))
  cat(sprintf("  Var(mu) estimate = %.4g   var-inflation = %.3f\n",
              x$variance_correction, x$var_inflation))
  if (!is.null(x$quality$status)) {
    cat(sprintf("  Quality: %s\n", x$quality$status))
    if (length(x$quality$warnings)) {
      for (w in x$quality$warnings) cat("  - ", w, "\n", sep = "")
    }
  }
  df <- data.frame(
    Estimate = x$coefficients,
    SE       = x$std_errors,
    z        = x$z_values,
    p        = x$p_values,
    check.names = FALSE)
  print(round(df, 5L))
  invisible(x)
}

# Null-coalescing helper -- internal, may already exist elsewhere but
# redefined here for standalone loading safety.
`%||%` <- function(a, b) if (is.null(a)) b else a
