# Quarantined multinomial compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no multinomial fit. Retained one-vs-rest and joint-softmax code
after the gate is unreachable through this public frontdoor and carries
no disclosure, DP, accuracy, or availability claim.

## Usage

``` r
ds.vertMultinom(
  formula,
  data = NULL,
  classes = NULL,
  reference = NULL,
  indicator_template = "%s_ind",
  max_iter = NULL,
  max_outer = 20L,
  tol = NULL,
  warm_max_iter = NULL,
  warm_tol = NULL,
  binomial_sigmoid_intervals = NULL,
  verbose = TRUE,
  datasources = NULL,
  ...,
  design_analysis_id = NULL
)
```

## Arguments

- formula, data, classes, reference, indicator_template, max_iter,
  max_outer, tol, warm_max_iter, warm_tol, binomial_sigmoid_intervals,
  verbose, datasources, design_analysis_id:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

- ...:

  Retained compatibility arguments; not evaluated.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires a signed bounded joint-softmax artifact over the
exact score design plus validated joint covariance and inference.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
