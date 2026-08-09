# Quarantined Cox-profile compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and returns no
Cox fit. The word “NonDisclosive” is a historical name, not a current
security or differential-privacy claim.

## Usage

``` r
ds.vertCoxProfileNonDisclosive(
  formula,
  data = NULL,
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

- formula, data, max_event_times, max_iter, tol, newton, ridge_eps,
  debug_trace, verbose, datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
