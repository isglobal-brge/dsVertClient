# Federated generalised estimating equations

Fit a GLM for vertically partitioned DataSHIELD data and return sandwich
(robust) standard errors alongside the usual model-based ones. The point
estimate \\\hat{\beta}\\ is obtained by a single call to
[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md)
(bread = `fit$covariance`), then the Liang-Zeger meat matrix is
estimated by a second weighted fit with weights proportional to the
squared Pearson residuals. This is the working-independence GEE
estimator; exchangeable / AR1 working correlation structures are a
planned Month-4 extension that adds a cluster-ID broadcast between the
DCF parties and per-cluster residual outer products.

Because both fits reveal only \\p\\-dimensional aggregates (gradient and
Hessian), the sandwich estimator is obtained without ever materialising
the \\n\\-length residual vector at the client. Inter-server leakage is
identical to the ordinary `ds.vertGLM` path (no new channels).

Formula: \$\$V\_{sand} = \mathrm{Cov}(\hat\beta) \\ A \\
\mathrm{Cov}(\hat\beta)\$\$ \$\$A = X^T \\ \mathrm{diag}(r^2) \\ X \\ /
\\ n\$\$

## Usage

``` r
ds.vertGEE(
  formula,
  data = NULL,
  family = c("gaussian", "binomial", "poisson"),
  id_col = NULL,
  corstr = c("independence", "exchangeable", "ar1"),
  max_iter = 100L,
  tol = 1e-04,
  lambda = 1e-04,
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

  Optional character; reserved for clustered GEE. In this first-pass
  implementation `id_col` is accepted but the working-independence
  sandwich is computed (cluster aggregation is the Month-4 follow-on).

- corstr:

  Working correlation structure. Only `"independence"` is currently
  supported; `"exchangeable"` and `"ar1"` raise a targeted message.

- max_iter, tol, lambda, verbose:

  Passed to `ds.vertGLM`.

- datasources:

  DataSHIELD connection object.

## Value

An object of class `ds.vertGEE` with components `coefficients`,
`model_se` (sqrt of `diag(Cov(beta))`), `robust_se` (sandwich SEs),
`robust_covariance`, `corstr`, and `fit` (the underlying `ds.glm`
object).
