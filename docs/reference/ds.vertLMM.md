# Federated linear mixed model with a single random intercept

Fit a random-intercept linear mixed model \\y\_{ij} = x\_{ij}^T \beta +
b_i + \varepsilon\_{ij}\\ on vertically partitioned DataSHIELD data,
where the cluster indicator `id` lives on the outcome server. The REML
log-likelihood profile is expressed in terms of a single variance ratio
\\\rho = \sigma_b^2 / (\sigma^2 + n_i\sigma_b^2)\\ so the outer
optimiser is one-dimensional. Each outer step calls
[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md)
with per-patient weights derived from the current \\\rho\\, reusing the
already-deployed `ds.vertGLM(weights=)` infrastructure.

Derivation (Laird-Ware compact form): \$\$V_i = \sigma^2 I + \sigma_b^2
\mathbf{1} \mathbf{1}^T\$\$ \$\$V_i^{-1} = \frac{1}{\sigma^2}\left(I -
\rho_i \mathbf{1}\mathbf{1}^T\right)\$\$ \$\$\log\|V_i\| = (n_i - 1)
\log \sigma^2 + \log(\sigma^2 + n_i\sigma_b^2)\$\$

Both summands are one-dimensional functions of \\\rho\\ that the client
evaluates on centralised aggregates (sum of \\n_i \rho_i\\ and sum of
\\\log(\sigma^2 + n_i \sigma_b^2)\\); per-cluster residual sums
\\\sum\_{ij} r\_{ij}\\ are returned by the outcome server as a single
aggregate vector (one scalar per cluster) under the already-documented
cluster-ID inter-server leakage tier.

Inter-server disclosure: the DCF peer learns integer cluster membership
(same class as Cox event-time ordering). Original cluster labels are not
returned, clusters below datashield.privacyLevel fail closed, and
individual observations are not revealed to the client.

## Usage

``` r
ds.vertLMM(
  formula,
  data = NULL,
  cluster_col,
  random_slopes = NULL,
  reml = TRUE,
  max_iter = 30L,
  inner_iter = 50L,
  tol = 1e-04,
  exact_cross_server = TRUE,
  sigma_b2_override = NULL,
  ring = c("ring63", "ring127"),
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  Fixed-effects formula.

- data:

  Aligned data-frame name.

- cluster_col:

  Column holding the cluster id (must be on the outcome server).

- random_slopes:

  Optional character vector of column names for random slopes (in
  addition to the random intercept). NULL fits an intercept-only
  random-effects structure.

- reml:

  Use REML (default TRUE) vs ML.

- max_iter:

  Outer variance-component iterations (default 30).

- inner_iter:

  Inner `ds.vertGLM` iteration budget.

- tol:

  Outer tolerance on \\\rho\\ change.

- exact_cross_server:

  Logical (default TRUE). If TRUE, use the exact cross-server Beaver
  vecmul gram pass; FALSE falls back to a diagonal-only approximation.

- sigma_b2_override:

  Numeric. Optional fixed value of the between-cluster variance,
  bypassing the iterative REML/ML variance-component update.

- ring:

  Character (`"ring63"` or `"ring127"`). Selects the Beaver vecmul
  pipeline ring; Ring127 (fracBits=50) is needed for STRICT closure on
  dense Gram matrices.

- verbose:

  Print progress.

- datasources:

  DataSHIELD connection object.

## Value

A `ds.vertLMM` object with components `coefficients`, `covariance`,
`std_errors`, `sigma2` (residual variance), `sigma_b2` (random-effect
variance), `icc`, `n_clusters`, `converged`, `iterations`, `fit` (final
inner `ds.glm`).
