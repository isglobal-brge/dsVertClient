# Quarantined multiple-imputation compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and mutates no
server data, draws no imputation, and returns no coefficients, counts,
covariance, or diagnostic. Retained implementation code after the gate
is unreachable through this public frontdoor and carries no disclosure,
DP, accuracy, or availability claim.

## Usage

``` r
ds.vertMI(
  formula,
  data = NULL,
  impute_columns = NULL,
  m = 20L,
  family = "gaussian",
  max_iter = 50L,
  tol = 1e-04,
  lambda = 0,
  intercept_only = c("error", "aggregate"),
  verbose = TRUE,
  datasources = NULL,
  seed = 1L
)
```

## Arguments

- formula, data, impute_columns, m, family, max_iter, tol, lambda,
  intercept_only, verbose, datasources, seed:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires a signed bounded imputation contract, non-rerollable
cryptographic randomness, no exact per-round count release, and
validated Rubin-rule inference.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
