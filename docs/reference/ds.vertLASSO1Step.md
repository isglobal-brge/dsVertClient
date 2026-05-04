# One-step LASSO via quadratic-surrogate proximal gradient

Fit a proper L1-regularised coefficient vector by expanding the
log-likelihood at the converged GLM point as a local quadratic and then
running proximal coordinate descent with soft- thresholding to minimise

    0.5 (beta - betahat)^T H (beta - betahat) + lambda ||beta||_1

where H = solve(Cov(betahat)) is the Fisher information matrix already
exposed by `ds.vertGLM` (`fit$covariance`). This yields a consistent,
efficient L1-regularised estimator in large samples without a second MPC
round: it uses only the already- returned coefficient vector and full
covariance matrix.

Useful for: - Regularised post-selection inference. - Model compression
while preserving the statistical geometry. - A baseline for the planned
Month 3 proper proximal-gradient LASSO that fits from scratch with
per-iteration MPC gradients.

## Usage

``` r
ds.vertLASSO1Step(
  fit,
  lambda,
  keep_intercept = TRUE,
  max_iter = 500L,
  tol = 1e-08
)
```

## Arguments

- fit:

  A `ds.glm` object with `fit$covariance` populated (commits \>=
  8bb7902).

- lambda:

  Numeric vector of L1 penalty values (a regularisation path).

- keep_intercept:

  Never penalise the intercept.

- max_iter:

  Coordinate-descent iterations per lambda.

- tol:

  Convergence tolerance on max \|Delta beta\|.

## Value

A ds.vertLASSO1Step object: per-lambda coefficient vectors, the penalty
path, the quadratic-surrogate objective at each lambda, and the fit.
