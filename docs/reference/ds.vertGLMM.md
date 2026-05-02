# Federated binomial GLMM via Laplace approximation

Fit a binomial generalised linear mixed model with a single random
intercept on vertically partitioned DataSHIELD data via the
Laplace-approximated marginal likelihood:

\$\$\ell_L(\beta, \sigma_b^2) = \sum_i \left\[\ell_i(\beta, b_i^\*) -
\tfrac{1}{2} \log \|H_i(b_i^\*)\| \right\]\$\$

where \\b_i^\*\\ is the mode of the per-cluster penalised log-likelihood
\\\ell_i(\beta, b) = \sum_j \log f(y\_{ij}\| x\_{ij}^T \beta + b) -
b^2/(2\sigma_b^2)\\ and \\H_i\\ is its Hessian at the mode.

Architecture:

1.  For each cluster (cluster IDs on outcome server) an INNER L-BFGS on
    \\b_i\\ is executed server-side (using the cached
    `dsvertClusterResidualsDS` aggregates + the cached \\\eta_i\\ share
    from ds.vertGLM(keep_session=TRUE)).

2.  The outer optimiser on \\(\beta, \sigma_b^2)\\ is client-side,
    re-calling ds.vertGLM with the shrinkage- weight column derived from
    the current \\b_i^\*\\.

3.  Variance-component update uses the moment-matching estimate
    \\\hat\sigma_b^2 = \mathrm{var}(\hat b_i)\\ across clusters
    (Breslow-Clayton approximation).

The scaffolding here reuses every primitive already shipped in dsVert
1.1.0+: the per-cluster residual aggregates, the keep_session flag on
ds.vertGLM, the Beaver vecmul for inner weighted updates, and the DCF
sigmoid/exp wide-splines.

Privacy: client sees only (beta, sigma_b^2) + per-cluster b_hat_i
estimates as an aggregate vector (one value per cluster). No per-patient
quantity ever leaves the DCF parties.

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

- verbose:

  Print progress.

- datasources:

  DataSHIELD connections.

## Value

`ds.vertGLMM` object: fixed-effect coefficients, cluster-level BLUPs
\\\hat b_i\\, random-effect variance \\\hat\sigma_b^2\\, and the
converged binomial `fit`.
