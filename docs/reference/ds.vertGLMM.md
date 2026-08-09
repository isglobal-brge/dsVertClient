# Quarantined generalized-linear mixed-model frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no GLMM coefficient, random effect, variance component, cluster
statistic, or diagnostic. Retained implementation code after the gate is
unreachable through this public frontdoor and carries no disclosure, DP,
accuracy, or availability claim.

## Usage

``` r
ds.vertGLMM(
  formula,
  data = NULL,
  cluster_col,
  max_outer = 30L,
  inner_iter = 50L,
  tol = 1e-04,
  lambda = 0,
  compute_se = FALSE,
  ring = NULL,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula, data, cluster_col, max_outer, inner_iter, tol, lambda,
  compute_se, ring, verbose, datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires bounded cluster contributions, a precisely specified
target estimator and validated covariance and diagnostics.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
