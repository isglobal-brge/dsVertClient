# Federated generalised estimating equations

Fit a GLM/GEE for vertically partitioned DataSHIELD data and return
sandwich (robust) standard errors alongside the usual model-based ones.
For `corstr = "independence"`, the point estimate \\\hat{\beta}\\ is
obtained by a single call to
[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md).
For Gaussian `corstr = "exchangeable"`, dsVert promotes the point
estimate to a protected cluster-level exchangeable GLS/GEE update. For
binomial and Poisson `corstr = "exchangeable"`, dsVert uses a Ring127
protected Pearson-score exchangeable GEE update. For Gaussian, binomial,
and Poisson `corstr = "ar1"`, dsVert uses guarded within-cluster order
metadata on the DCF parties and returns only low-dimensional AR1
aggregates. For Gaussian, binomial, and Poisson independence models the
sandwich meat is computed in the share domain. If `id_col` is supplied,
dsVert computes the clustered meat \\\sum_c S_c S_c^\top\\, where
\\S_c=\sum\_{i\in c} X_i r_i\\; otherwise it computes the row-level HC0
meat \\X^\top \mathrm{diag}(r^2) X\\. Cluster membership is
transport-encrypted between DCF parties and clusters below
`datashield.privacyLevel` fail closed.

The client receives only low-dimensional aggregate matrices/vectors. It
never materialises the \\n\\-length residual, squared-residual, or
weighted-column vectors.

Formula: \$\$V\_{sand} = \mathrm{Cov}(\hat\beta) \\ A \\
\mathrm{Cov}(\hat\beta)\$\$ \$\$A = \sum_c S_c S_c^\top\$\$ for
clustered Gaussian/binomial/Poisson fits, or \$\$A = X^T \\
\mathrm{diag}(r^2) \\ X\$\$ for row-level HC0 when no `id_col` is
supplied.

## Usage

``` r
ds.vertGEE(
  formula,
  data = NULL,
  family = c("gaussian", "binomial", "poisson"),
  id_col = NULL,
  order_col = NULL,
  corstr = c("independence", "exchangeable", "ar1"),
  max_iter = 100L,
  tol = 1e-04,
  lambda = 1e-04,
  working_max_iter = NULL,
  ring = 63L,
  binomial_sigmoid_intervals = NULL,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  A model formula passed through to `ds.vertGLM`.

- data:

  Character. Aligned data-frame name on each server.

- family:

  One of `"gaussian"`, `"binomial"`, `"poisson"`.

- id_col:

  Optional character. For Gaussian/binomial/Poisson models this enables
  the cluster-robust sandwich meat; the cluster column must live with
  the outcome and all clusters must pass `datashield.privacyLevel`.

- order_col:

  Optional character. Required for `corstr = "ar1"`. The order column
  must live with the outcome and is used only server-to-server to derive
  guarded adjacent records.

- corstr:

  Working correlation. `"independence"` is available for
  Gaussian/binomial/Poisson. `"exchangeable"` currently fits true
  exchangeable Gaussian, binomial, and Poisson GEE coefficients from
  guarded cluster-level sufficient statistics. `"ar1"` fits Gaussian,
  binomial, and Poisson GEE from guarded adjacent-pair sufficient
  statistics.

- max_iter, tol, lambda, verbose:

  Passed to `ds.vertGLM`.

- working_max_iter:

  Optional integer. Maximum iterations for exchangeable/AR1
  working-correlation updates. Defaults to `max_iter`.

- ring:

  Integer 63 or 127. Binomial/Poisson exchangeable and all AR1 routes
  are automatically run in Ring127 because protected nonlinear link
  operations and adjacent-product statistics need high precision.

- binomial_sigmoid_intervals:

  Optional integer. Number of DCF spline intervals for protected
  binomial sigmoid evaluations used by the underlying GLM fit and GEE
  sandwich/working-correlation updates. When `NULL`,
  `dsvert.gee_binomial_sigmoid_intervals`,
  `dsvert.glm_num_intervals_binomial`, or the 100-interval default is
  used.

- datasources:

  DataSHIELD connection object.

## Value

An object of class `ds.vertGEE` with components `coefficients`,
`model_se` (sqrt of `diag(Cov(beta))`), `robust_se` (sandwich SEs),
`robust_covariance`, `corstr`, and `fit` (the underlying `ds.glm`
object).
