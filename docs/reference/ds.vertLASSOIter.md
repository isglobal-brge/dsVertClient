# Federated LASSO path from a non-disclosive GLM fit

Fit an L1-regularised path without revealing row-level gradients or
residuals. For Gaussian models this is the proper LASSO objective solved
from the normal equations already exposed by
[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md):

    \deqn{\arg\min_\beta \tfrac{1}{2n}\|y-X\beta\|^2 +
               \lambda\|\beta_{-0}\|_1.}

The solver delegates to
[`ds.vertLASSOProximal`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOProximal.md),
which reconstructs \\X^T X/n\\ from aggregate covariance/Hessian output
and runs coordinate descent client-side. Binomial models use the
internal aggregate score path to run proximal-gradient on the
standardized L1 objective. Poisson models use fixed-start aggregate
score/Hessian probes and a damped proximal-Newton update on the same
standardized L1 objective.

Non-disclosure: the inner ds.vertGLM call already hides everything at
the p-aggregate level; this wrapper only manipulates returned aggregate
coefficients, covariance/Hessian matrices, or p-dimensional aggregate
scores.

## Usage

``` r
ds.vertLASSOIter(
  formula,
  data = NULL,
  family = c("gaussian", "binomial", "poisson"),
  lambda = NULL,
  max_outer = 20L,
  tol = 1e-08,
  alpha = 0.5,
  inner_iter = 8L,
  exact_non_gaussian = TRUE,
  ring = NULL,
  lipschitz = c("auto", "gram", "safe"),
  fista_restart = TRUE,
  binomial_sigmoid_intervals = getOption("dsvert.lasso_binomial_sigmoid_intervals", 200L),
  poisson_damping = 0.5,
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

  Retained for backward compatibility; used as a lower bound on the
  Gaussian coordinate-descent iteration budget.

- tol:

  Tolerance passed to the LASSO path solver.

- alpha:

  Retained for backward compatibility; no longer used.

- inner_iter:

  Inner `ds.vertGLM` budget for the initial unpenalised fit.

- exact_non_gaussian:

  Logical. For binomial models, use repeated secure aggregate-score
  evaluations and proximal-gradient updates on the standardized L1
  objective instead of the one-step surrogate. For Poisson models, use
  repeated aggregate score/Hessian probes and a damped proximal-Newton
  update.

- ring:

  Integer (63 or 127). MPC ring for the GLM score evaluations. Defaults
  to Ring63. Ring127 can be requested explicitly after fixture
  validation for the target deployment.

- lipschitz:

  Character. Step-size rule for binomial proximal-gradient: `"auto"`
  (default) tries a guarded aggregate Gram/correlation bound and falls
  back to the conservative `0.25 * (p + 1)` bound; `"gram"` requires the
  aggregate Gram bound; `"safe"` always uses the conservative bound.

- fista_restart:

  Logical. For binomial proximal-gradient, reset FISTA momentum when the
  client-side gradient-restart criterion triggers. This uses only
  p-dimensional coefficient vectors already held by the client and does
  not require opening objectives or row-level quantities.

- binomial_sigmoid_intervals:

  Integer or NULL. Number of secure sigmoid spline intervals used by
  binomial aggregate-score evaluations. The default,
  `getOption("dsvert.lasso_binomial_sigmoid_intervals", 200L)`, is
  intentionally finer than the generic GLM default because L1
  optimisation repeatedly reuses the score oracle and amplifies the
  spline floor. Set NULL to inherit the ambient
  `dsvert.glm_num_intervals_binomial` option.

- poisson_damping:

  Numeric in (0, 1\]. Fixed damping applied to Poisson proximal-Newton
  proposals. The default 0.5 is intentionally conservative and avoids an
  extra objective/deviance pass.

- verbose:

  Print progress.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertLASSOIter` object with components `lambda`, `paths`
(per-lambda coefficient vectors), `n_outer` (solver iterations used per
lambda), `final_fit` (the unpenalised `ds.glm` fit), and `method`
describing the estimator target.
