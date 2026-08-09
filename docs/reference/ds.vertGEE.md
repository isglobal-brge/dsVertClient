# Quarantined GEE compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no GEE coefficient, working correlation, sandwich covariance,
cluster statistic, or diagnostic. Retained implementation code after the
gate is unreachable through this public frontdoor and carries no
disclosure, DP, accuracy, or availability claim.

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
  lambda = 0,
  working_max_iter = NULL,
  ring = 63L,
  binomial_sigmoid_intervals = NULL,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula, data, family, id_col, order_col, corstr, max_iter, tol,
  lambda, working_max_iter, ring, binomial_sigmoid_intervals, verbose,
  datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires contribution-bounded cluster score/meat artifacts, a
protected working-correlation contract and validated robust covariance.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
