# Quarantined Cox compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no model, statistic, diagnostic, or patient-derived metadata.
Retained implementation code after the local migration gate is
unreachable through this frontdoor and is not evidence of disclosure
safety, differential privacy, numerical accuracy, or availability.

## Usage

``` r
ds.vertCox(
  formula,
  data = NULL,
  max_iter = 30L,
  tol = 1e-04,
  max_event_times = NULL,
  newton = TRUE,
  ridge_eps = 1e-06,
  debug_trace = FALSE,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula, data, max_iter, tol, max_event_times, newton, ridge_eps,
  debug_trace, verbose, datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires a purpose-bound contribution-limited Cox capsule,
private risk-set evaluation, certified ties/strata/delayed-entry
semantics, covariance and identifiability certificates, and independent
multi-host validation.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
