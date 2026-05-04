#' @title Federated LASSO path from a non-disclosive GLM fit
#' @description Fit an L1-regularised path without revealing row-level
#'   gradients or residuals. For Gaussian models this is the proper LASSO
#'   objective solved from the normal equations already exposed by
#'   \code{\link{ds.vertGLM}}:
#'
#'     \deqn{\arg\min_\beta \tfrac{1}{2n}\|y-X\beta\|^2 +
#'            \lambda\|\beta_{-0}\|_1.}
#'
#'   The solver delegates to \code{\link{ds.vertLASSOProximal}}, which
#'   reconstructs \eqn{X^T X/n} from aggregate covariance/Hessian output
#'   and runs coordinate descent client-side. Binomial models use the internal
#'   aggregate score path to run proximal-gradient on the standardized L1
#'   objective. Poisson models use fixed-start aggregate score/Hessian probes
#'   and a damped proximal-Newton update on the same standardized L1 objective.
#'
#'   Non-disclosure: the inner ds.vertGLM call already hides everything
#'   at the p-aggregate level; this wrapper only manipulates returned
#'   aggregate coefficients, covariance/Hessian matrices, or p-dimensional
#'   aggregate scores.
#'
#' @param formula Model formula.
#' @param data Aligned data-frame name.
#' @param family GLM family.
#' @param lambda L1 penalty scalar or vector (path).
#' @param max_outer Retained for backward compatibility; used as a lower
#'   bound on the Gaussian coordinate-descent iteration budget.
#' @param tol Tolerance passed to the LASSO path solver.
#' @param alpha Retained for backward compatibility; no longer used.
#' @param inner_iter Inner \code{ds.vertGLM} budget for the initial
#'   unpenalised fit.
#' @param exact_non_gaussian Logical. For binomial models, use repeated
#'   secure aggregate-score evaluations and proximal-gradient updates on
#'   the standardized L1 objective instead of the one-step surrogate. For
#'   Poisson models, use repeated aggregate score/Hessian probes and a
#'   damped proximal-Newton update.
#' @param ring Integer (63 or 127). MPC ring for the GLM score evaluations.
#'   Defaults to Ring63. Ring127 can be requested explicitly after fixture
#'   validation for the target deployment.
#' @param lipschitz Character. Step-size rule for binomial proximal-gradient:
#'   \code{"auto"} (default) tries a guarded aggregate Gram/correlation
#'   bound and falls back to the conservative \code{0.25 * (p + 1)} bound;
#'   \code{"gram"} requires the aggregate Gram bound; \code{"safe"} always
#'   uses the conservative bound.
#' @param fista_restart Logical. For binomial proximal-gradient, reset FISTA
#'   momentum when the client-side gradient-restart criterion triggers. This
#'   uses only p-dimensional coefficient vectors already held by the client
#'   and does not require opening objectives or row-level quantities.
#' @param poisson_damping Numeric in (0, 1]. Fixed damping applied to Poisson
#'   proximal-Newton proposals. The default 0.5 is intentionally conservative
#'   and avoids an extra objective/deviance pass.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertLASSOIter} object with components
#'   \code{lambda}, \code{paths} (per-lambda coefficient vectors),
#'   \code{n_outer} (solver iterations used per lambda),
#'   \code{final_fit} (the unpenalised \code{ds.glm} fit), and
#'   \code{method} describing the estimator target.
#' @export
ds.vertLASSOIter <- function(formula, data = NULL,
                              family = c("gaussian", "binomial", "poisson"),
                              lambda = NULL,
                              max_outer = 20L, tol = 1e-8,
                              alpha = 0.5, inner_iter = 8L,
                              exact_non_gaussian = TRUE,
                              ring = NULL,
                              lipschitz = c("auto", "gram", "safe"),
                              fista_restart = TRUE,
                              poisson_damping = 0.5,
                              verbose = TRUE, datasources = NULL) {
  family <- match.arg(family)
  lipschitz <- match.arg(lipschitz)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (is.null(lambda)) lambda <- c(1e-3, 1e-2, 0.1, 0.5, 1.0)
  lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  if (any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("lambda must contain non-negative finite values", call. = FALSE)
  }
  soft <- function(x, t) sign(x) * pmax(abs(x) - t, 0)
  exact_binomial <- isTRUE(exact_non_gaussian) && identical(family, "binomial")
  exact_poisson <- isTRUE(exact_non_gaussian) && identical(family, "poisson")
  if (is.null(ring)) ring <- 63L
  ring <- as.integer(ring)
  poisson_damping <- as.numeric(poisson_damping)
  if (!is.finite(poisson_damping) || poisson_damping <= 0 ||
      poisson_damping > 1) {
    stop("poisson_damping must be in (0, 1]", call. = FALSE)
  }
  fista_restart <- isTRUE(fista_restart)

  .orig_to_std <- function(beta_orig, fit) {
    beta_orig <- as.numeric(beta_orig)
    names(beta_orig) <- names(fit$coefficients)
    out <- beta_orig
    slope_names <- setdiff(names(beta_orig), "(Intercept)")
    if (length(slope_names) > 0L) {
      out[slope_names] <- beta_orig[slope_names] *
        as.numeric(fit$x_sds[slope_names])
    }
    if ("(Intercept)" %in% names(beta_orig)) {
      out["(Intercept)"] <- beta_orig["(Intercept)"] +
        sum(beta_orig[slope_names] * as.numeric(fit$x_means[slope_names]))
    }
    out
  }

  .std_to_orig <- function(theta_std, fit) {
    theta_std <- as.numeric(theta_std)
    names(theta_std) <- names(fit$coefficients)
    out <- theta_std
    slope_names <- setdiff(names(theta_std), "(Intercept)")
    if (length(slope_names) > 0L) {
      out[slope_names] <- theta_std[slope_names] /
        as.numeric(fit$x_sds[slope_names])
    }
    if ("(Intercept)" %in% names(theta_std)) {
      out["(Intercept)"] <- theta_std["(Intercept)"] -
        sum(out[slope_names] * as.numeric(fit$x_means[slope_names]))
    }
    out
  }

  .score_std <- function(theta_std) {
    fit_g <- ds.vertGLM(
      formula, data = data, family = family,
      max_iter = 1L, tol = tol, lambda = 0,
      ring = ring, start = theta_std,
      compute_se = FALSE, compute_deviance = FALSE,
      gradient_only = TRUE,
      verbose = FALSE, datasources = datasources)
    g <- fit_g$gradient_std
    if (is.null(g)) stop("GLM score evaluation returned no gradient",
                         call. = FALSE)
    g <- g[names(theta_std)]
    if (any(!is.finite(g))) {
      stop("non-finite GLM aggregate score", call. = FALSE)
    }
    g
  }

  .score_hessian_std <- function(theta_std) {
    fit_g <- ds.vertGLM(
      formula, data = data, family = family,
      max_iter = 1L, tol = tol, lambda = 0,
      ring = ring, start = theta_std,
      compute_se = TRUE, compute_deviance = FALSE,
      gradient_only = TRUE,
      verbose = FALSE, datasources = datasources)
    g <- fit_g$gradient_std
    H <- fit_g$hessian_std
    if (is.null(g) || is.null(H)) {
      stop("GLM fixed evaluation returned no score/Hessian", call. = FALSE)
    }
    g <- g[names(theta_std)]
    H <- as.matrix(H)
    if (!is.null(rownames(H))) {
      H <- H[names(theta_std), names(theta_std), drop = FALSE]
    }
    H <- (H + t(H)) / 2
    if (any(!is.finite(g)) || any(!is.finite(H))) {
      stop("non-finite GLM aggregate score/Hessian", call. = FALSE)
    }
    list(score = g, hessian = H)
  }

  .safe_lipschitz <- function(theta_names) {
    0.25 * length(theta_names)
  }

  .gram_lipschitz <- function(fit) {
    slope_names <- setdiff(names(fit$coefficients), "(Intercept)")
    if (length(slope_names) == 0L) return(0.25)
    cor_fit <- ds.vertCor(data, variables = slope_names,
                          verbose = FALSE, datasources = datasources)
    corr <- as.matrix(cor_fit$correlation)
    corr <- corr[slope_names, slope_names, drop = FALSE]
    corr <- (corr + t(corr)) / 2
    # glmStandardizeDS uses sample sd, so crossprod(X_std) / n equals
    # ((n - 1) / n) * cor(X), while the intercept block is orthogonal to
    # centred predictors and has Gram value 1.
    n_obs <- as.numeric(fit$n_obs)
    slope_gram <- corr * ((n_obs - 1) / n_obs)
    gram <- diag(length(slope_names) + 1L)
    gram[-1L, -1L] <- slope_gram
    ev <- eigen((gram + t(gram)) / 2, symmetric = TRUE,
                only.values = TRUE)$values
    L <- 0.25 * max(ev, na.rm = TRUE)
    if (!is.finite(L) || L <= 0) {
      stop("non-finite aggregate Gram Lipschitz bound", call. = FALSE)
    }
    as.numeric(L)
  }

  .choose_lipschitz <- function(fit, theta_names) {
    if (identical(lipschitz, "safe")) {
      return(list(L = .safe_lipschitz(theta_names), source = "safe"))
    }
    gram_try <- tryCatch(.gram_lipschitz(fit), error = function(e) e)
    if (inherits(gram_try, "error")) {
      if (identical(lipschitz, "gram")) {
        stop("aggregate Gram Lipschitz bound failed: ",
             conditionMessage(gram_try), call. = FALSE)
      }
      if (verbose) {
        message("[LASSOIter] aggregate Gram Lipschitz unavailable; ",
                "using conservative bound: ",
                conditionMessage(gram_try))
      }
      return(list(L = .safe_lipschitz(theta_names), source = "safe_fallback"))
    }
    list(L = gram_try, source = "aggregate_gram")
  }

  .binomial_pg <- function(fit0, lam, warm_std, step_bound) {
    theta <- warm_std
    yk <- theta
    tk <- 1
    int_idx <- match("(Intercept)", names(theta))
    penalized <- seq_along(theta)
    if (!is.na(int_idx)) penalized <- setdiff(penalized, int_idx)
    L <- step_bound$L
    if (!is.finite(L) || L <= 0) L <- 1
    converged <- FALSE
    used <- as.integer(max_outer)
    restart_count <- 0L
    for (it in seq_len(as.integer(max_outer))) {
      old <- theta
      grad <- .score_std(yk)
      z <- yk - grad / L
      theta_new <- z
      theta_new[penalized] <- soft(z[penalized], lam / L)
      names(theta_new) <- names(theta)
      restart_score <- sum((theta_new - theta) * (yk - theta_new))
      restart <- fista_restart && is.finite(restart_score) &&
        restart_score > 0
      if (restart) {
        tk_new <- 1
        yk <- theta_new
        restart_count <- restart_count + 1L
      } else {
        tk_new <- (1 + sqrt(1 + 4 * tk^2)) / 2
        yk <- theta_new + ((tk - 1) / tk_new) * (theta_new - theta)
      }
      theta <- theta_new
      tk <- tk_new
      if (max(abs(theta - old)) < tol) {
        converged <- TRUE
        used <- it
        break
      }
    }
    beta_orig <- .std_to_orig(theta, fit0)
    l1 <- lam * sum(abs(theta[penalized]))
    list(coefficients = beta_orig, theta_std = theta,
         iterations = used, converged = converged,
         objective = NA_real_, l1_penalty = l1, L = L,
         fista_restart = fista_restart, restart_count = restart_count)
  }

  .prox_quadratic_cd <- function(theta, grad, H, lam, max_cd = 1000L) {
    H <- (H + t(H)) / 2
    diag(H) <- pmax(diag(H), 1e-8)
    beta <- theta
    int_idx <- match("(Intercept)", names(theta))
    penalized <- seq_along(theta)
    if (!is.na(int_idx)) penalized <- setdiff(penalized, int_idx)
    for (cd_iter in seq_len(as.integer(max_cd))) {
      old <- beta
      delta <- beta - theta
      for (j in seq_along(beta)) {
        delta_no_j <- delta
        delta_no_j[j] <- 0
        r_j <- theta[j] - (grad[j] + sum(H[j, ] * delta_no_j)) / H[j, j]
        beta[j] <- if (j %in% penalized) {
          soft(r_j, lam / H[j, j])
        } else {
          r_j
        }
        delta[j] <- beta[j] - theta[j]
      }
      if (max(abs(beta - old)) < tol) break
    }
    names(beta) <- names(theta)
    beta
  }

  .poisson_prox_newton <- function(fit0, lam, warm_std) {
    theta <- warm_std
    int_idx <- match("(Intercept)", names(theta))
    penalized <- seq_along(theta)
    if (!is.na(int_idx)) penalized <- setdiff(penalized, int_idx)
    converged <- FALSE
    used <- as.integer(max_outer)
    last_step <- NA_real_
    for (it in seq_len(as.integer(max_outer))) {
      old <- theta
      eval <- .score_hessian_std(theta)
      proposal <- .prox_quadratic_cd(theta, eval$score, eval$hessian, lam)
      theta <- theta + poisson_damping * (proposal - theta)
      last_step <- max(abs(theta - old))
      if (last_step < tol) {
        converged <- TRUE
        used <- it
        break
      }
    }
    beta_orig <- .std_to_orig(theta, fit0)
    l1 <- lam * sum(abs(theta[penalized]))
    list(coefficients = beta_orig, theta_std = theta,
         iterations = used, converged = converged,
         objective = NA_real_, l1_penalty = l1,
         damping = poisson_damping, last_step = last_step)
  }

  # Prime the warm start with an unpenalised fit.
  if (verbose) message("[LASSOIter] priming with unpenalised ds.vertGLM")
  prime_iter <- if (exact_binomial) {
    max(as.integer(inner_iter), 0L)
  } else if (exact_poisson) {
    max(as.integer(inner_iter), 8L)
  } else {
    max(as.integer(inner_iter), 20L)
  }
  fit0 <- ds.vertGLM(formula, data = data, family = family,
                     max_iter = prime_iter, tol = tol,
                     lambda = 0,
                     ring = ring,
                     compute_se = !(exact_binomial || exact_poisson),
                     compute_deviance = !(exact_binomial || exact_poisson),
                     verbose = FALSE, datasources = datasources)

  paths <- list()
  n_outer_used <- integer(length(lambda))
  objectives <- numeric(length(lambda))
  solver_objects <- list()
  method <- if (identical(family, "gaussian")) {
    "gaussian_normal_equations"
  } else if (exact_binomial) {
    "binomial_standardized_proximal_gradient"
  } else if (exact_poisson) {
    "poisson_standardized_proximal_newton"
  } else {
    "quadratic_surrogate"
  }
  warm_std <- if (exact_binomial || exact_poisson) {
    .orig_to_std(fit0$coefficients, fit0)
  } else {
    NULL
  }
  step_bound <- if (exact_binomial) {
    .choose_lipschitz(fit0, names(warm_std))
  } else {
    NULL
  }
  if (exact_binomial && verbose) {
    message(sprintf("[LASSOIter] binomial step bound L=%.6g (%s)",
                    step_bound$L, step_bound$source))
  }

  for (li in seq_along(lambda)) {
    lam <- lambda[li]
    if (verbose) message(sprintf("[LASSOIter] lambda = %.4g", lam))
    if (identical(family, "gaussian")) {
      obj <- ds.vertLASSOProximal(
        fit0, lambda = lam,
        max_iter = max(2000L, as.integer(max_outer)),
        tol = tol)
      paths[[sprintf("%.6g", lam)]] <- obj$coefficients
      n_outer_used[li] <- obj$iterations
      objectives[li] <- obj$objective
      solver_objects[[sprintf("%.6g", lam)]] <- obj
    } else if (exact_binomial) {
      obj <- .binomial_pg(fit0, lam, warm_std, step_bound)
      obj$lipschitz_source <- step_bound$source
      paths[[sprintf("%.6g", lam)]] <- obj$coefficients
      n_outer_used[li] <- obj$iterations
      objectives[li] <- obj$objective
      solver_objects[[sprintf("%.6g", lam)]] <- obj
      warm_std <- obj$theta_std
    } else if (exact_poisson) {
      obj <- .poisson_prox_newton(fit0, lam, warm_std)
      paths[[sprintf("%.6g", lam)]] <- obj$coefficients
      n_outer_used[li] <- obj$iterations
      objectives[li] <- obj$objective
      solver_objects[[sprintf("%.6g", lam)]] <- obj
      warm_std <- obj$theta_std
    } else {
      if (li == 1L && verbose) {
        message("[LASSOIter] non-Gaussian family: using one-step ",
                "quadratic-surrogate LASSO target")
      }
      obj <- ds.vertLASSO1Step(
        fit0, lambda = lam,
        max_iter = max(500L, as.integer(max_outer)),
        tol = tol)
      paths[[sprintf("%.6g", lam)]] <- obj$paths[[1L]]
      n_outer_used[li] <- NA_integer_
      objectives[li] <- obj$objective[[1L]]
      solver_objects[[sprintf("%.6g", lam)]] <- obj
    }
  }

  out <- list(
    lambda      = lambda,
    paths       = paths,
    n_outer     = n_outer_used,
    objective   = objectives,
    final_fit   = fit0,
    family      = family,
    method      = method,
    solver      = solver_objects,
    alpha       = alpha,
    exact_non_gaussian = exact_binomial || exact_poisson,
    ring        = ring,
    lipschitz   = if (is.null(step_bound)) NA_character_ else step_bound$source,
    poisson_damping = if (exact_poisson) poisson_damping else NA_real_,
    fista_restart = if (exact_binomial) fista_restart else NA,
    call        = match.call())
  class(out) <- c("ds.vertLASSOIter", "list")
  out
}

#' @export
print.ds.vertLASSOIter <- function(x, ...) {
  cat("dsVert LASSO path\n")
  cat(sprintf("  Family : %s\n", x$family))
  cat(sprintf("  Method : %s\n", x$method))
  cat(sprintf("  Lambda : %s\n",
              paste(sprintf("%.4g", x$lambda), collapse = " ")))
  cat(sprintf("  Solver iters : %s\n",
              paste(x$n_outer, collapse = " ")))
  m <- do.call(cbind, x$paths)
  cat("\nCoefficient path:\n")
  print(round(m, 5L))
  invisible(x)
}
