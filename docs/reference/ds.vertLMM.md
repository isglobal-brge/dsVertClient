# Quarantined linear-mixed-model compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no mixed model, variance component, cluster statistic, or
diagnostic. Retained LMM implementation code after the gate is
unreachable through this public frontdoor and carries no disclosure, DP,
accuracy, or availability claim.

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

- formula, data, cluster_col, random_slopes, reml, max_iter, inner_iter,
  tol, exact_cross_server, sigma_b2_override, ring, verbose,
  datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires contribution-bounded cluster statistics, private
cluster handling, certified ML/REML and random-effects semantics,
covariance and identifiability certificates, and independent multi-host
validation.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
