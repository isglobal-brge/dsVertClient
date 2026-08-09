# Quarantined joint-softmax Newton compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and returns no
multinomial fit. Retained Newton/MPC code after the gate is unreachable
through this public frontdoor and carries no disclosure, DP, accuracy,
or availability claim.

## Usage

``` r
ds.vertMultinomJointNewton(
  formula,
  data = NULL,
  levels,
  indicator_template = "%s_ind",
  max_outer = 20L,
  tol = 1e-05,
  warm_max_iter = NULL,
  warm_tol = NULL,
  binomial_sigmoid_intervals = NULL,
  verbose = TRUE,
  datasources = NULL,
  design_analysis_id = NULL
)
```

## Arguments

- formula, data, levels, indicator_template, max_outer, tol,
  warm_max_iter, warm_tol, binomial_sigmoid_intervals, verbose,
  datasources, design_analysis_id:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires a purpose-bound signed `multinomial_design_grams`
artifact over the exact bounded score design and a validated
joint-softmax inference contract.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
