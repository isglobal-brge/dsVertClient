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
#'   and runs coordinate descent client-side. For binomial and Poisson GLMs,
#'   a true proximal solver would need repeated secure gradient evaluations
#'   at sparse warm starts, which the current GLM API does not expose. Those
#'   families therefore return the documented \code{\link{ds.vertLASSO1Step}}
#'   quadratic-surrogate path instead of pretending to be exact.
#'
#'   Non-disclosure: the inner ds.vertGLM call already hides everything
#'   at the p-aggregate level; this wrapper only manipulates returned
#'   aggregate coefficients and covariance/Hessian matrices.
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
                              verbose = TRUE, datasources = NULL) {
  family <- match.arg(family)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (is.null(lambda)) lambda <- c(1e-3, 1e-2, 0.1, 0.5, 1.0)
  lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  if (any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("lambda must contain non-negative finite values", call. = FALSE)
  }

  # Prime the warm start with an unpenalised fit.
  if (verbose) message("[LASSOIter] priming with unpenalised ds.vertGLM")
  fit0 <- ds.vertGLM(formula, data = data, family = family,
                     max_iter = max(inner_iter, 20L), tol = tol,
                     lambda = 0,
                     verbose = FALSE, datasources = datasources)

  paths <- list()
  n_outer_used <- integer(length(lambda))
  objectives <- numeric(length(lambda))
  solver_objects <- list()
  method <- if (identical(family, "gaussian")) {
    "gaussian_normal_equations"
  } else {
    "quadratic_surrogate"
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
