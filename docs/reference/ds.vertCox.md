# Federated Cox proportional-hazards regression

Fit a Cox PH model on vertically partitioned DataSHIELD data using the
reverse-cumsum reformulation of the partial- likelihood score. The
client runs an L-BFGS outer loop; each step obtains the aggregate
gradient \$\$\nabla \ell(\beta) = \sum_j x_j (\delta_j - e^{\eta_j}
G_j)\$\$ \$\$G_j = \sum\_{i: \delta_i=1, t_i \le t_j} 1 / S(t_i)\$\$
\$\$S(t_i) = \sum\_{k: t_k \ge t_i} e^{\eta_k}\$\$ via the
already-deployed server helpers (k2SetCoxTimesDS /
k2ApplyCoxPermutationDS / k2CoxReverseCumsumSDS / k2StoreCoxRecipDS /
k2CoxForwardCumsumGDS) and the existing 4-phase DCF protocol with
family="exp" + family="reciprocal". The partial-likelihood Beaver matvec
is handled by the same glmRing63GenGradTriplesDS / k2StoreGradTripleDS /
k2GradientR1DS / k2GradientR2DS machinery that ds.vertGLM uses, so no
new cryptographic round is introduced.

Client view per iteration: the p-dimensional aggregate gradient and a
scalar partial log-likelihood (via the standard Beaver-sum path). The
client never sees \\\eta_j\\, \\S(t_j)\\, \\G_j\\, or any per-patient
quantity.

Inter-server disclosure: the DCF peer learns the ascending-time sort
permutation (ranking of event times) and the binary event indicator.
Absolute event times are NOT disclosed.

## Usage

``` r
ds.vertCox(
  formula,
  data = NULL,
  time_col = NULL,
  event_col = NULL,
  tstart_col = NULL,
  strata_col = NULL,
  max_iter = 30L,
  tol = 1e-04,
  lambda = 1e-04,
  compute_loglik = FALSE,
  compute_se = FALSE,
  num_intervals_exp = 75L,
  num_intervals_recip = 75L,
  one_step_newton = TRUE,
  newton_refine_iters = 5L,
  newton_refine_tol = 1e-05,
  ring = 127L,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  Formula of the form `Surv(time, event) ~ x1 + ...`. If the LHS is not
  a `Surv(...)` expression, supply `time_col` / `event_col` explicitly.

- data:

  Aligned data-frame name on each server.

- time_col, event_col:

  Column names on the outcome server.

- tstart_col:

  Optional left-truncation / start-time column for counting-process /
  time-varying covariates. When supplied, the data are interpreted as
  (tstart, time, event) triplets.

- strata_col:

  Optional stratifying column name (categorical). Each level fits a
  separate baseline hazard with shared beta.

- max_iter:

  Outer L-BFGS iterations (default 30).

- tol:

  Convergence tolerance on max \|delta beta\|.

- lambda:

  L2 regularisation.

- compute_loglik:

  Logical (default FALSE). If TRUE, run an extra DCF pipeline pass at
  beta_hat to evaluate the partial log-likelihood for AIC/BIC / LR
  tests. Adds ~60s per fit (Opal) or ~6 min (local).

- compute_se:

  Logical (default FALSE). If TRUE, evaluate the observed information
  matrix at beta_hat and return covariance / standard errors. Adds one
  full DCF pipeline pass.

- num_intervals_exp:

  Integer. DCF spline grid size for exp(eta); 75 is the sweet spot on
  Ring63 (error \< 1e-3 while ~25 percent faster than 100). Drop to 50
  for small-n speed.

- num_intervals_recip:

  Integer. DCF spline grid size for the reciprocal 1/S(t); same default
  rationale as `num_intervals_exp`.

- one_step_newton:

  Logical. If TRUE (default), use the bias-free one-step Newton at
  beta=0 followed by fixed-Fisher refinement; FALSE falls back to the
  legacy iterative gradient-descent loop.

- newton_refine_iters:

  Integer. Number of damped fixed-Fisher refinement iterations after the
  one-step Newton (HARD CAP 5 by P3 disclosure budget).

- newton_refine_tol:

  Numeric. Convergence tolerance on max \|delta beta\| during
  refinement.

- ring:

  Integer (63 or 127). Selects the MPC ring / fracBits pipeline. Default
  127 (Catrina-Saxena fracBits=50, ~1e-15 per-op). Pass 63 to force the
  legacy Ring63 pipeline.

- verbose:

  Print progress.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertCox` object: `coefficients`, `std_errors`, `covariance`,
`loglik` (partial log-likelihood), `n_obs`, `n_events`, `iterations`,
`converged`.
