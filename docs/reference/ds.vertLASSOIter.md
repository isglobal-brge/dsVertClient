# Iterative proximal-gradient LASSO over the MPC GLM gradient

Proper proximal-gradient L1-regularised GLM fitting that invokes the
full MPC gradient pipeline at each outer step. At iteration \\t\\:
\$\$\nabla_t = \mathrm{ds.vertGLM}(\beta_t)\$\$ \$\$\beta\_{t+1} =
\mathrm{soft}(\beta_t - \alpha \nabla_t, \alpha \lambda)\$\$ where
\\\alpha\\ is a backtracking step size chosen to guarantee descent on
the smooth part of the objective, and the soft-threshold operator
enforces the L1 sparsity pattern. The intercept is not penalised.

Unlike
[`ds.vertLASSO1Step`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO1Step.md)
(which post-hoc soft-thresholds the converged GLM solution via a local
quadratic surrogate around \\\hat\beta\\), this routine re-evaluates the
gradient at the current sparse iterate \\\beta_t\\, so the final
estimate is the true proximal-gradient L1 solution rather than a
one-step surrogate. Costs M + 1 MPC GLM gradient invocations for M iter.

Non-disclosure: the inner ds.vertGLM call already hides everything at
the p-aggregate level; this outer loop only manipulates the returned
\\p\\-vector of coefficients.

## Usage

``` r
ds.vertLASSOIter(
  formula,
  data = NULL,
  family = c("gaussian", "binomial", "poisson"),
  lambda = NULL,
  max_outer = 20L,
  tol = 0.001,
  alpha = 0.5,
  inner_iter = 8L,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  Model formula.

- data:

  Aligned data-frame name.

- family:

  GLM family.

- lambda:

  L1 penalty scalar or vector (path).

- max_outer:

  Outer proximal-gradient iterations (default 20).

- tol:

  Outer tolerance on \\\\\beta_t - \beta\_{t-1}\\\_\infty\\.

- alpha:

  Initial step size (default 0.5); a simple halving line search is
  applied.

- inner_iter:

  Inner `ds.vertGLM` budget per outer step. A small value (5-10) is
  usually sufficient since we're only reading the gradient at the
  warm-start point.

- verbose:

  Print progress.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertLASSOIter` object with components `lambda`, `paths`
(per-lambda coefficient vectors), `n_outer` (outer iterations used per
lambda), `final_fit` (the last inner `ds.glm` fit).
