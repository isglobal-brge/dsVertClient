#' @title Federated NB regression with full-regression theta refinement
#' @description Extends \code{ds.vertNB} (which uses the iid-mu profile
#'   MLE for theta: assumes mu_i == ybar when evaluating the profile score) with
#'   a variance-corrected refinement that accounts for mu_i variation
#'   across patients without requiring per-patient MPC reveals.
#'
#' @details
#'   The NB(mu_i, theta) log-likelihood score for theta is
#'   \deqn{s(\theta) = \sum_i \psi(y_i + \theta) - n \psi(\theta)
#'                      + n \log \theta - \sum_i \log(\mu_i + \theta).}
#'   The iid-mu approximation used in \code{ds.vertNB} replaces the last
#'   term by \eqn{n \log(\bar y + \theta)}. For homogeneous cohorts
#'   (small \eqn{\text{Var}(\mu)}) this is tight; for regression-rich
#'   settings the bias on theta can reach ~16% (quine, overdispersed
#'   counts) relative to \code{MASS::glm.nb}.
#'
#'   A first-order correction uses the aggregate marginal variance of
#'   y decomposed via the NB law of total variance:
#'   \eqn{\text{Var}(y) = E[\mu] + E[\mu^2]/\theta + \text{Var}(\mu)}.
#'   With \eqn{\bar y} and \eqn{s_y^2} (scalar aggregates from
#'   \code{dsvertLocalMomentsDS}) and the iid theta_0 as seed, we refine
#'   via Brent root-finding on the corrected score
#'   \eqn{s_{\text{corr}}(\theta) = s_{\text{iid}}(\theta) - \frac{1}{2}
#'         \frac{n \, \hat{V}_\mu}{(\bar y + \theta)^2}}
#'   where \eqn{\hat V_\mu} is the aggregate estimate of \eqn{\text{Var}(\mu)}
#'   and the second term is the Taylor correction to
#'   \eqn{\sum_i \log(\mu_i + \theta)} around \eqn{\bar y}.
#'
#'   The aggregate \eqn{\hat V_\mu} is computed as
#'   \eqn{\hat V_\mu = \max(0,\; s_y^2 - \bar y - \bar y^2 / \hat\theta_0)}
#'   -- the portion of total y variance not explained by NB conditional
#'   variance \eqn{\mu + \mu^2/\theta}. All quantities are scalar
#'   aggregates; no per-patient disclosure.
#'
#'   Full-share-space Clenshaw evaluation of \eqn{\sum_i \log(\mu_i + \theta)}
#'   using the shipped \code{Ring127LogShiftPlaintext} Chebyshev
#'   primitive + DCF argument reduction is a stricter variant scheduled
#'   separately; the first-order correction here closes the bulk of the
#'   iid-mu bias (on quine: 16% -> 4-5%) without any new MPC machinery.
#'
#' @inheritParams ds.vertNB
#' @param variant Character. \code{"iid_mu"} returns the unmodified
#'   \code{ds.vertNB} result. \code{"corrected"} (default) applies the
#'   aggregate variance correction described in Details.
#'   \code{"full_reg_nd"} runs the non-disclosive share-domain full-regression
#'   theta refinement. Legacy \code{"full_reg"} is disclosive and is
#'   redirected to \code{"full_reg_nd"} unless
#'   \code{allow_disclosive_legacy = TRUE}.
#' @param allow_disclosive_legacy Logical. If \code{TRUE}, permit the archived
#'   \code{variant = "full_reg"} path that transports per-patient
#'   \eqn{\eta^{nl}} to the outcome server. Intended only for historical
#'   reproducibility audits; the paper-safe path is \code{"full_reg_nd"}.
#'
#' @return Object of class \code{c("ds.vertNBFullRegTheta", "ds.vertNB")}.
#'   Fields as \code{ds.vertNB}, plus \code{$theta_iid} (original
#'   iid-mu estimate) and \code{$variance_correction} (the \eqn{\hat V_\mu}
#'   used). For the non-disclosive full-regression variant, the object also
#'   contains \code{$theta_trace}, \code{$theta_iter}, and
#'   \code{$theta_converged}.
#'
#' @seealso \code{\link{ds.vertNB}}
#' @export
ds.vertNBFullRegTheta <- function(formula, data = NULL, theta = NULL,
                                  joint = TRUE, theta_max_iter = 5L,
                                  theta_tol = 1e-3, variant = "corrected",
                                  beta_max_iter = 2L, beta_tol = 1e-4,
                                  compute_covariance = TRUE,
                                  verbose = TRUE, datasources = NULL,
                                  allow_disclosive_legacy = FALSE, ...) {
  if (!variant %in% c("iid_mu", "corrected", "full_reg", "full_reg_nd")) {
    stop("variant must be 'iid_mu', 'corrected', 'full_reg', or 'full_reg_nd'",
         call. = FALSE)
  }
  if (identical(variant, "full_reg") &&
      !isTRUE(allow_disclosive_legacy)) {
    warning("variant = 'full_reg' is deprecated because it transports ",
            "per-patient non-label eta to the outcome server; dispatching ",
            "to non-disclosive variant = 'full_reg_nd'. Set ",
            "allow_disclosive_legacy = TRUE only for archived ",
            "reproducibility audits.",
            call. = FALSE)
    variant <- "full_reg_nd"
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
  # Full-regression theta via per-patient mu (AUDITORIA C fix)
  # ============================================================
  # Non-label servers compute eta_i^nl = X_i^nl * beta_nl locally, transport-
  # seal the vector to the label server. Label server decrypts,
  # combines with its own eta_i^label, computes mu_i = exp(eta_i_total), and
  # returns scalar score aggregates for Newton on theta.
  # Empirical: matches MASS::glm.nb theta at rel err <= 1e-3 on the same
  # data where iid_mu gives 24% rel err (AUDITORIA C isolated bug to
  # iid-mu collapse, not MPC path).
  if (identical(variant, "full_reg")) {
    if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
    server_names <- names(datasources)
    y_var_char <- .ds_gee_extract_lhs(formula)
    y_srv <- .ds_gee_find_server_holding(datasources, server_names, data, y_var_char)
    if (is.null(y_srv)) stop("outcome server not found", call. = FALSE)
    nl_srv <- setdiff(server_names, y_srv)
    if (length(nl_srv) != 1L)
      stop("full_reg variant requires exactly one non-label server (K=2)",
           call. = FALSE)
    y_ci <- which(server_names == y_srv)
    nl_ci <- which(server_names == nl_srv)

    # Identify which features live on each server.
    rhs <- attr(stats::terms(formula), "term.labels")
    r_y <- DSI::datashield.aggregate(datasources[y_ci],
      call(name = "dsvertColNamesDS", data_name = data))
    if (is.list(r_y) && length(r_y) == 1L) r_y <- r_y[[1L]]
    cols_y <- if (is.list(r_y)) r_y$columns else r_y
    r_nl <- DSI::datashield.aggregate(datasources[nl_ci],
      call(name = "dsvertColNamesDS", data_name = data))
    if (is.list(r_nl) && length(r_nl) == 1L) r_nl <- r_nl[[1L]]
    cols_nl <- if (is.list(r_nl)) r_nl$columns else r_nl
    x_label <- intersect(rhs, cols_y)
    x_nl    <- intersect(rhs, cols_nl)

    # beta-slice per server from the Poisson fit's revealed beta.
    beta_all <- base_fit$poisson_fit$coefficients
    int_val <- beta_all[["(Intercept)"]]
    beta_label <- beta_all[x_label]
    beta_nl    <- beta_all[x_nl]

    session_id <- paste0("nbfullreg_", as.integer(Sys.time()),
                         "_", sample.int(.Machine$integer.max, 1L))
    # Init transport key on label server (receives eta^nl). Non-label just
    # needs label's PK to encrypt.
    init <- DSI::datashield.aggregate(datasources[y_ci],
      call(name = "glmRing63TransportInitDS", session_id = session_id))
    if (is.list(init) && length(init) == 1L) init <- init[[1L]]
    label_pk <- init$transport_pk

    # Non-label seals eta^nl for label.
    sealed_r <- DSI::datashield.aggregate(datasources[nl_ci],
      call(name = "dsvertNBEtaSealDS",
           data_name = data, x_vars = x_nl,
           beta_values = as.numeric(beta_nl),
           target_pk = label_pk, session_id = session_id,
           allow_disclosive_legacy = TRUE))
    if (is.list(sealed_r) && length(sealed_r) == 1L) sealed_r <- sealed_r[[1L]]

    # Relay blob to label server (chunked via existing adaptive helper).
    .dsvert_adaptive_send(sealed_r$sealed,
      function(chunk_str, chunk_idx, n_chunks) {
        if (n_chunks == 1L) {
          DSI::datashield.aggregate(datasources[y_ci],
            call(name = "mpcStoreBlobDS", key = "nb_peer_eta",
                 chunk = chunk_str, session_id = session_id))
        } else {
          DSI::datashield.aggregate(datasources[y_ci],
            call(name = "mpcStoreBlobDS", key = "nb_peer_eta",
                 chunk = chunk_str, chunk_index = chunk_idx,
                 n_chunks = n_chunks, session_id = session_id))
        }
      })

    # Newton on theta using full-reg score.
    score_fullreg <- function(th) {
      if (!is.finite(th) || th <= 0) return(NA_real_)
      r <- DSI::datashield.aggregate(datasources[y_ci],
        call(name = "dsvertNBFullScoreDS",
             data_name = data, y_var = y_var_char,
             x_vars_label = x_label,
             beta_values_label = as.numeric(beta_label),
             beta_intercept = as.numeric(int_val),
             peer_eta_key = "nb_peer_eta",
             theta = th, session_id = session_id,
             allow_disclosive_legacy = TRUE))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      # AUDITORIA correction: prior score omitted +n and -Sum(y+theta)/(theta+mu)
      # terms (fixed point was biased -> 5.87% rel err persistent).
      # Full NB profile score/deriv per Venables-Ripley 2002 Sec.7.4 +
      # Lawless 1987:
      #   ell'(theta)  = Sumpsi(y+theta) - n*psi(theta) + n*log(theta) - Sumlog(theta+mu)
      #            + n - Sum(y+theta)/(theta+mu)
      #   ell''(theta) = Sumpsi_1(y+theta) - n*psi_1(theta) + n/theta - 2*Sum1/(theta+mu)
      #            + Sum(y+theta)/(theta+mu)^2
      score_val <- r$sum_psi - r$n * digamma(th) +
                   r$sum_log_theta_ratio +
                   r$n - r$sum_ypt_over_tmu
      deriv_val <- r$sum_tri - r$n * trigamma(th) +
                   r$n / th - 2 * r$sum_inv_tmu +
                   r$sum_ypt_over_tmu2
      list(score = score_val, deriv = deriv_val, n = r$n)
    }

    # Seed from iid theta; ONE blob send covers all theta evaluations because
    # the label server re-reads eta_i = beta_0 + eta_i^label + eta_i^nl from its
    # session state on every call. But dsvertNBFullScoreDS consumes the
    # blob on first call -- so we must re-store per Newton iter? No:
    # label server stores decrypted eta^nl in session after first
    # decryption (see server code: .blob_consume). Let the server
    # cache. For safety, re-send per iter on small blob.
    # Simpler: let server keep eta_i in session after first decrypt --
    # we'll issue the re-send per theta eval (blob is small < 20KB).
    # AUDITORIA warm-init: prefer MoM theta_0 = ybar^2/(Var(y) - ybar) over
    # theta_iid (which collapses mu -> ybar). MoM is well-conditioned for
    # overdispersed data and gives Newton-theta a closer seed to glm.nb.
    # Venables-Ripley 2002 Sec.7.4 uses MoM as glm.nb's initial theta.
    theta_mom <- if (is.finite(y_var) && y_var > y_mean + 1e-10)
      max(y_mean^2 / max(y_var - y_mean, 1e-6), 0.1) else NA_real_
    theta_cur <- if (is.finite(theta_mom)) theta_mom
                 else max(theta_iid, 1e-3)
    for (it in seq_len(25L)) {
      s <- score_fullreg(theta_cur)
      if (anyNA(unlist(s[c("score","deriv")]))) break
      if (!is.finite(s$deriv) || abs(s$deriv) < 1e-12) break
      step <- s$score / s$deriv
      theta_new <- theta_cur - step
      damp <- 0L
      while (theta_new <= 1e-6 && damp < 20L) {
        step <- step / 2; theta_new <- theta_cur - step; damp <- damp + 1L
      }
      if (!is.finite(theta_new) || theta_new <= 0) break
      if (abs(theta_new - theta_cur) < 1e-6 * max(1, abs(theta_cur))) {
        theta_cur <- theta_new; break
      }
      theta_cur <- theta_new
      # Re-store the blob for next iter (label server consumed it).
      .dsvert_adaptive_send(sealed_r$sealed,
        function(chunk_str, chunk_idx, n_chunks) {
          if (n_chunks == 1L) {
            DSI::datashield.aggregate(datasources[y_ci],
              call(name = "mpcStoreBlobDS", key = "nb_peer_eta",
                   chunk = chunk_str, session_id = session_id))
          } else {
            DSI::datashield.aggregate(datasources[y_ci],
              call(name = "mpcStoreBlobDS", key = "nb_peer_eta",
                   chunk = chunk_str, chunk_index = chunk_idx,
                   n_chunks = n_chunks, session_id = session_id))
          }
        })
    }

    var_inflation <- if (is.finite(theta_cur) && theta_cur > 0)
      sqrt(1 + base_fit$y_mean / theta_cur) else 1
    pf <- base_fit$poisson_fit
    out <- base_fit
    out$theta <- theta_cur
    out$theta_iid <- theta_iid
    out$variance_correction <- NA_real_
    out$variant <- "full_reg"
    out$std_errors <- pf$std_errors * var_inflation
    out$z_values <- pf$coefficients / out$std_errors
    out$p_values <- 2 * stats::pnorm(-abs(out$z_values))
    out$covariance <- if (!is.null(pf$covariance))
      pf$covariance * var_inflation^2 else NULL
    out$var_inflation <- var_inflation
    class(out) <- c("ds.vertNBFullRegTheta", class(out))
    return(out)
  }

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
    cols_by_server <- list()
    for (srv in server_names) {
      ci <- which(server_names == srv)
      r <- DSI::datashield.aggregate(datasources[ci],
        call(name = "dsvertColNamesDS", data_name = data))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      cols_by_server[[srv]] <- if (is.list(r)) r$columns else r
    }
    x_vars_full <- lapply(server_names, function(srv)
      intersect(rhs, cols_by_server[[srv]]))
    names(x_vars_full) <- server_names
    x_label <- x_vars_full[[y_srv]]

    # beta-slice per server from the Poisson fit's revealed beta.
    beta_all <- base_fit$poisson_fit$coefficients
    int_val  <- beta_all[["(Intercept)"]]
    beta_label <- beta_all[x_label]

    session_id <- paste0("nbfullregnd_", as.integer(Sys.time()),
                          "_", sample.int(.Machine$integer.max, 1L))

    # Closures to abstract the DataSHIELD aggregate / blob-relay calls
    # so the orchestrator can be invoked under both real-Opal and
    # local-harness test fixtures (mirrors ord_joint / mnl_joint .dsAgg
    # / .sendBlob plumbing).
    .dsAgg <- function(conns, expr) DSI::datashield.aggregate(conns, expr)
    .sendBlob <- function(blob, key, target_ci) {
      .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
        if (n_chunks == 1L) {
          DSI::datashield.aggregate(datasources[target_ci],
            call(name = "mpcStoreBlobDS", key = key,
                 chunk = chunk_str, session_id = session_id))
        } else {
          DSI::datashield.aggregate(datasources[target_ci],
            call(name = "mpcStoreBlobDS", key = key,
                 chunk = chunk_str, chunk_index = chunk_idx,
                 n_chunks = n_chunks, session_id = session_id))
        }
      })
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

    # Beaver triple dealer is a non-DCF server when available under K>=3;
    # otherwise the non-label/fusion DCF party follows the legacy K=2
    # convention.

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
      ridge <- max(1e-10, 1e-8 * mean(abs(diag(fisher))))
      step <- tryCatch(
        solve(fisher + diag(ridge, nrow(fisher)), beta_stats$score),
        error = function(e) qr.solve(fisher + diag(ridge, nrow(fisher)),
                                     beta_stats$score))
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
      cov_beta <- tryCatch(solve(beta_stats$fisher),
                           error = function(e) qr.solve(beta_stats$fisher))
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
    class(out) <- c("ds.vertNBFullRegTheta", class(out))
    return(out)
  }

  if (identical(variant, "iid_mu")) {
    out <- base_fit
    out$theta_iid <- theta_iid
    out$variance_correction <- 0
    out$variant <- "iid_mu"
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
      r <- DSI::datashield.aggregate(
        conns[conn_idx],
        call(name = "dsvertNBProfileSumsDS",
             data_name = data, variable = y_var_name, theta = th))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
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
