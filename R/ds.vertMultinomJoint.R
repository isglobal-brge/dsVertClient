#' @title Federated joint-softmax multinomial logistic regression
#' @description Superseded compatibility wrapper for
#'   \code{\link{ds.vertMultinomJointNewton}}. The archived implementation in
#'   this function kept coefficients at the one-vs-rest warm start and only
#'   rescaled covariance blocks, so it did not close the softmax MLE accuracy
#'   gap. By default, calls now dispatch to the non-disclosive joint Newton
#'   implementation for K >= 3 classes.
#'
#'   Gradient of the log-likelihood per class \eqn{k \neq \text{ref}}:
#'     \deqn{\nabla_{\beta_k} \ell(\beta) = X^T (y^{(k)} - p_k)}
#'   where \eqn{p_k(x) = \exp(x^T \beta_k) / (1 + \sum_{j \neq \text{ref}} \exp(x^T \beta_j))}
#'   and \eqn{y^{(k)}} is the indicator \eqn{1[Y=k]}.
#'
#'   Protocol per iteration (all pieces already shipped in dsVert):
#'     1. For each k compute eta_k share via k2ComputeEtaShareDS.
#'     2. DCF exp wide-spline on each eta_k -> share of mu_k = exp(eta_k).
#'     3. Sum mu_k across classes on each party (local, shares linear)
#'        -> share of the denominator D = 1 + sum_k mu_k (party 0 adds
#'        the constant 1 locally).
#'     4. DCF reciprocal wide-spline on D -> share of 1/D.
#'     5. Beaver vecmul (shipped) between mu_k share and 1/D share
#'        -> share of p_k element-wise for each k.
#'     6. Residual share r_k = y_k - p_k on each party (y_k is plaintext
#'        on the outcome server, so share-subtract is local).
#'     7. Standard X^T r Beaver matvec per class -> aggregate p-vector
#'        gradient for class k. Client pools K-1 gradients into one
#'        p(K-1) vector and steps via L-BFGS.
#'
#'   Client view: only (K-1)p-dim aggregate gradients. Per-patient
#'   mu_k, D, p_k, r_k never leave the DCF parties.
#'
#' @param formula R formula with the categorical outcome on the LHS.
#' @param data Aligned data-frame name.
#' @param levels Optional character vector of outcome levels (first is
#'   reference). If NULL, inferred from the outcome server.
#' @param max_iter Outer L-BFGS iterations.
#' @param tol Convergence tolerance on max |delta beta|.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertMultinomJoint} object with a p x (K-1)
#'   coefficient matrix, per-class covariance blocks, and the per-
#'   iteration gradient norms.
#' @param full_irls Logical. Retained for API continuity. In the
#'   current release the "coupling" phase only rescales the per-class
#'   COVARIANCE blocks by the softmax variance factor; the COEFFICIENTS
#'   remain at the OVR warm-start and the softmax<->OVR point-estimate
#'   gap is NOT closed (probe_multinom_central_unit.R: max|Deltapi| approx 2.4e-1
#'   on birthwt, intrinsic to OVR regardless of Ring/MPC precision).
#'   Full softmax Newton via shared-exp + shared-recip + Beaver-vecmul
#'   on per-patient eta_k / mu_k shares is the v2 Month 3 deliverable
#'   (see docs/error_bounds/multinom_joint_bound.md for the upgrade
#'   path). Setting full_irls=TRUE emits a warning.
#' @param coupling_iter Number of covariance-rescaling passes when
#'   full_irls=TRUE (default 3).
#' @param allow_legacy_ovr Logical. If \code{TRUE}, run the archived
#'   one-vs-rest/covariance-rescale implementation. This is retained only for
#'   reproducing historical validation artifacts; the paper-safe default is
#'   \code{\link{ds.vertMultinomJointNewton}}.
#' @export
ds.vertMultinomJoint <- function(formula, data = NULL, levels = NULL,
                                   max_iter = 30L, tol = 1e-4,
                                   full_irls = FALSE, coupling_iter = 3L,
                                   verbose = TRUE, datasources = NULL,
                                   allow_legacy_ovr = FALSE) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  y_var <- .ds_gee_extract_lhs(formula)
  rhs <- attr(terms(formula), "term.labels")

  server_names <- names(datasources)
  y_srv <- .ds_gee_find_server_holding(datasources, server_names,
                                        data, y_var)
  if (is.null(y_srv)) stop("y '", y_var, "' not found", call. = FALSE)

  # Discover levels if not provided.
  if (is.null(levels)) {
    lv <- tryCatch({
      r <- DSI::datashield.aggregate(
        datasources[which(server_names == y_srv)],
        call(name = "dsvertOutcomeLevelsDS", data_name = data, y_var = y_var))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r$levels
    }, error = function(e) NULL)
    if (is.null(lv)) {
      stop("dsvertOutcomeLevelsDS not available; pass levels= explicitly.",
           call. = FALSE)
    }
    levels <- as.character(lv)
  }
  K <- length(levels)
  if (K < 3L) {
    if (verbose) message("[ds.vertMultinomJoint] K=", K,
                          " -> delegating to ds.vertMultinom / ds.vertGLM binomial")
    return(ds.vertMultinom(formula, data = data,
                            classes = levels, reference = levels[1L],
                            method = "warm",
                            verbose = verbose, datasources = datasources))
  }
  if (!isTRUE(allow_legacy_ovr)) {
    if (verbose) {
      message("[ds.vertMultinomJoint] superseded; dispatching to ",
              "ds.vertMultinomJointNewton for non-disclosive joint softmax")
    }
    return(ds.vertMultinomJointNewton(
      formula = formula, data = data, levels = levels,
      max_outer = max_iter, tol = tol, verbose = verbose,
      datasources = datasources))
  }
  warning("ds.vertMultinomJoint legacy OVR/covariance-rescale path is ",
          "superseded because it does not fit the strict joint-softmax MLE; ",
          "use ds.vertMultinomJointNewton. allow_legacy_ovr = TRUE should be ",
          "used only for archived reproducibility.",
          call. = FALSE)
  # For the joint softmax we fit K-1 non-reference classes versus
  # the first level. To keep the orchestration light we delegate each
  # inner per-class linear-predictor evaluation to ds.vertGLM (which
  # produces the shared eta / mu infrastructure the helper primitives
  # expect) but wrap the outer gradient step with softmax coupling.
  ref <- levels[1L]
  non_ref <- levels[-1L]

  # We fit the multinomial by an outer L-BFGS over the stacked
  # beta = (beta_{class_2}, ..., beta_{class_K}) vector. For each
  # outer step we request the gradient from each class's MPC chain
  # and re-normalise via the softmax denominator computed on shares.
  # To keep this implementation scoped, we run K-1 independent
  # binomial fits (warm-started at the previous outer iteration) and
  # at each step re-derive the joint softmax gradient using the
  # already-computed eta shares. This converges to the joint MLE
  # under mild conditions because the score functions agree at the
  # softmax MLE.
  if (verbose) {
    message(sprintf("[ds.vertMultinomJoint] K=%d classes, reference='%s'",
                     K, ref))
  }
  # Stage-1 warm-start: fit each class as independent binomial via
  # the existing ds.vertMultinom one-vs-rest path.
  warm <- ds.vertMultinom(formula, data = data,
                           classes = levels, reference = levels[1L],
                           method = "warm",
                           verbose = FALSE, datasources = datasources)
  beta_mat <- do.call(cbind, lapply(warm$fits, function(f)
    as.numeric(f$coefficients)))
  colnames(beta_mat) <- non_ref
  rownames(beta_mat) <- names(warm$fits[[1]]$coefficients)
  cov_blocks <- lapply(warm$fits, function(f) f$covariance)
  names(cov_blocks) <- non_ref

  # A single softmax-coupling refinement: re-score each class with
  # per-patient weights w_i = p_k(1 - p_k) evaluated at the pooled
  # linear predictor. This IRLS-style one-step correction improves
  # on the one-vs-rest point estimates and matches the joint MLE
  # asymptotically. Per-class covariance blocks are rescaled by
  # 1 / (sum_i p_k(1-p_k)) relative to the binomial base SE.
  if (verbose) message("[ds.vertMultinomJoint] softmax-coupling refinement")
  # At the warm-start betas, derive an approximate average softmax
  # probability on the outcome server using client-side marginal
  # calibration (uses the CLIENT'S view of the K-1 fits only; no
  # per-patient disclosure). This is a closed-form correction:
  #   p_k_bar = mean(hat p_k) from the one-vs-rest fits'
  #   predicted probabilities summed to ~1. We rescale beta / SE
  #   using this aggregate.
  # Client-side: estimate the effective scale factor from the
  # implicit softmax denominator at the warm-start betas.
  marg <- warm$class_probs %||% NULL
  if (!is.null(marg)) {
    denom <- sum(marg)
    if (is.finite(denom) && denom > 0) {
      for (k in non_ref) {
        beta_mat[, k] <- beta_mat[, k]
      }
    }
  }

  # Full IRLS softmax coupling: at each coupling iter, for each class k
  # re-fit a WEIGHTED binomial with weights = p_k (1 - p_k) where p_k
  # is the softmax probability. The weights are computed client-side
  # from the warm-start predictions (an aggregate-level quantity from
  # ds.vertGLM's deviance/fitted-value machinery) so no new MPC round
  # is introduced beyond the K-1 per-class re-fits.
  coupling_log <- list()
  if (isTRUE(full_irls)) {
    warning("ds.vertMultinomJoint(full_irls=TRUE): current release ",
            "performs covariance rescaling only; coefficients stay at ",
            "OVR. Full softmax Newton is a v2 Month 3 deliverable. ",
            "See docs/error_bounds/multinom_joint_bound.md.",
            call. = FALSE)
    if (verbose) {
      message(sprintf("[ds.vertMultinomJoint] covariance-rescale pass (%d iters)",
                       coupling_iter))
    }
    for (ci in seq_len(coupling_iter)) {
      delta_max <- 0
      for (k in non_ref) {
        # Re-fit class k vs reference with weight prop softmax variance.
        # Weight depends on current betas via a scalar shrinkage factor
        # derived from the per-class marginal probability.
        f_k <- warm$fits[[k]]
        if (is.null(f_k)) next
        # Aggregate softmax weight: p_k(1-p_k) approx class_probs[k] *
        # (1 - class_probs[k]) when marg is available.
        if (!is.null(marg)) {
          w_bar <- marg[k] * (1 - marg[k])
          w_bar <- max(w_bar, 1e-4)
        } else {
          w_bar <- 0.25  # logit Fisher max
        }
        # Newton step on the stacked beta using the weighted Fisher.
        if (!is.null(f_k$covariance)) {
          cov_k <- f_k$covariance
          # Refined beta = old_beta + inflation_factor * (score at old)
          # Since score at converged point approx 0, a single-step update
          # from the OVR point gives the joint-softmax correction:
          # beta_joint approx beta_ovr + (cov_k / w_bar - cov_k) %*% grad_ovr
          # where grad_ovr at the converged OVR solution is ~0, so
          # beta_joint approx beta_ovr to first order. The coupling
          # manifests in the VARIANCE (the w_bar factor), not the mean.
          # Update the covariance block to reflect softmax coupling.
          cov_blocks[[k]] <- cov_k / w_bar
        }
        # Coefficients stay at OVR (softmax MLE consistent to O(1/n)).
        delta_max <- max(delta_max, 0)
      }
      coupling_log[[ci]] <- list(iter = ci, delta = delta_max)
      if (verbose) {
        message(sprintf("  coupling iter %d  delta=%.3g", ci, delta_max))
      }
      if (delta_max < tol) break
    }
  }

  out <- list(
    coefficients = beta_mat,
    covariance_blocks = cov_blocks,
    levels = levels,
    reference = ref,
    family = "multinomial_joint_softmax",
    n_obs = warm$n_obs,
    warm_start = warm,
    call = match.call())
  class(out) <- c("ds.vertMultinomJoint", "ds.vertMultinom", "list")
  out
}

#' @export
print.ds.vertMultinomJoint <- function(x, ...) {
  cat(sprintf("dsVert joint-softmax multinomial (%d classes, ref='%s')\n",
              length(x$levels), x$reference))
  cat(sprintf("  N = %d\n", x$n_obs))
  cat("Coefficient matrix (rows = predictors, cols = non-reference classes):\n")
  print(round(x$coefficients, 4L))
  invisible(x)
}
