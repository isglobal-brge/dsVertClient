# Quarantined deprecated K\>=3 LMM frontdoor

This exported deprecated name is retained for API compatibility. It
raises a typed `dsvert_route_unavailable` condition before any DSI call
and returns no mixed model. Retained K\>=3 implementation code after the
gate is unreachable through this public frontdoor and carries no
disclosure, DP, accuracy, or availability claim.

## Usage

``` r
ds.vertLMM.k3(
  formula,
  data,
  cluster_col,
  rho_lo = 0.001,
  rho_hi = 0.999,
  tol = 1e-04,
  max_outer = 30L,
  ring = c("ring127", "ring63"),
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula, data, cluster_col, rho_lo, rho_hi, tol, max_outer, ring,
  verbose, datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
