# Quarantined joint-multinomial compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and returns no
joint-softmax fit. Retained code after the gate is unreachable through
this public frontdoor and carries no disclosure, DP, accuracy, or
availability claim.

## Usage

``` r
ds.vertMultinomJoint(
  formula,
  data = NULL,
  levels = NULL,
  max_iter = 30L,
  tol = 1e-04,
  verbose = TRUE,
  datasources = NULL,
  design_analysis_id = NULL
)
```

## Arguments

- formula, data, levels, max_iter, tol, verbose, datasources,
  design_analysis_id:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
