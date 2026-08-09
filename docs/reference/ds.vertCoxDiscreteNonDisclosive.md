# Quarantined discrete-hazard compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no model, statistic, diagnostic, or patient-derived metadata.
The retained pooled-logistic research code is unreachable through this
public frontdoor and is not a Cox fallback or a public disclosure, DP,
accuracy, or availability claim.

## Usage

``` r
ds.vertCoxDiscreteNonDisclosive(
  formula,
  data = NULL,
  J = 5L,
  bin_breaks = NULL,
  target = c("discrete_logit", "cox_profile"),
  max_event_times = NULL,
  max_iter = 20L,
  tol = 1e-06,
  newton = TRUE,
  ridge_eps = 1e-06,
  debug_trace = FALSE,
  verbose = FALSE,
  datasources = NULL
)
```

## Arguments

- formula, data, J, bin_breaks, target, max_event_times, max_iter, tol,
  newton, ridge_eps, debug_trace, verbose, datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

A future discrete-time route must use its own signed bounded hazard
estimand and must remain explicitly distinct from Cox partial
likelihood.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
