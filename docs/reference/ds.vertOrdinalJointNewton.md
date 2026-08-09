# Quarantined proportional-odds Newton compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and returns no
ordinal fit. Retained Newton/MPC code after the gate is unreachable
through this public frontdoor and carries no disclosure, DP, accuracy,
or availability claim.

## Usage

``` r
ds.vertOrdinalJointNewton(
  formula,
  data = NULL,
  levels_ordered,
  cumulative_template = "%s_leq",
  max_outer = 20L,
  tol = 1e-05,
  warm_max_iter = NULL,
  warm_tol = NULL,
  binomial_sigmoid_intervals = NULL,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula, data, levels_ordered, cumulative_template, max_outer, tol,
  warm_max_iter, warm_tol, binomial_sigmoid_intervals, verbose,
  datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
