# Federated binomial GLMM-PQL

Fit a binomial generalised linear mixed model with a single random
intercept on vertically partitioned DataSHIELD data. The paper-supported
route is aggregate PQL: working responses, weights, probabilities,
residuals, and row scores remain in Ring127 share-domain statistics and
the client receives fixed effects, scalar variance components, and
guarded diagnostics only.

Architecture:

1.  For each cluster (cluster IDs on outcome server), binomial score and
    information terms \\\sum(y-p)\\ and \\\sum p(1-p)\\ are computed
    from DCF shares. The outcome server broadcasts cluster membership
    only to the DCF peer; the client receives aggregate cluster sums
    only.

2.  The outer optimiser on \\(\beta, \sigma_b^2)\\ is client-side,
    re-calling ds.vertGLM with the shrinkage- weight column derived from
    the current \\b_i^\*\\.

3.  Variance-component update uses the moment-matching estimate
    \\\hat\sigma_b^2 = \mathrm{var}(\hat b_i)\\ across clusters
    (Breslow-Clayton approximation).

The scaffolding here reuses every primitive already shipped in dsVert
1.1.0+: the keep_session flag on ds.vertGLM, cluster membership
broadcast, per-cluster share sums, Beaver vecmul, and the DCF sigmoid
wide-spline.

Privacy: the returned object contains fixed effects, variance-component
estimates, scalar diagnostics, and the final inner `ds.vertGLM` fit.
Per-patient quantities never leave the DCF parties. Guarded per-cluster
sufficient statistics are used internally for the random-intercept
update, but per-cluster BLUPs and cluster-size vectors are not returned.
Original cluster labels are not returned and clusters below
`datashield.privacyLevel` fail closed.

Inter-server leakage: cluster membership (same tier as ds.vertLMM).

## Usage

``` r
ds.vertGLMM(
  formula,
  data = NULL,
  cluster_col,
  max_outer = 10L,
  inner_iter = 10L,
  tol = 0.001,
  lambda = 0,
  compute_se = TRUE,
  ring = NULL,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  Fixed-effects formula (binomial outcome on LHS).

- data:

  Aligned data-frame name.

- cluster_col:

  Cluster id column on the outcome server.

- max_outer:

  Outer (beta, sigma_b^2) iterations.

- inner_iter:

  Inner PIRLS iterations per cluster per outer step.

- tol:

  Outer convergence tolerance.

- lambda:

  L2 penalty passed to the inner binomial GLM fits.

- compute_se:

  Logical. Compute GLM finite-difference standard errors for the inner
  fits. Set FALSE for coefficient/variance validation runs.

- ring:

  Integer 63 or 127. The PQL aggregate route requires Ring127.

- verbose:

  Print progress.

- datasources:

  DataSHIELD connections.

## Value

`ds.vertGLMM` object: fixed-effect coefficients, random-effect variance
\\\hat\sigma_b^2\\, scalar fit diagnostics, and the converged binomial
`fit`. Per-cluster BLUP vectors are internal working quantities and are
not returned.
