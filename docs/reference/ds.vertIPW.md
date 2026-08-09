# Quarantined IPW compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no propensity model, weights, effect estimate, or diagnostic.
Retained two-stage code after the gate is unreachable through this
public frontdoor and carries no disclosure, DP, causal-identification,
accuracy, or availability claim.

## Usage

``` r
ds.vertIPW(
  outcome_formula,
  propensity_formula,
  data = NULL,
  weights_column = "ipw",
  outcome_family = "gaussian",
  verbose = TRUE,
  datasources = NULL,
  ...
)
```

## Arguments

- outcome_formula, propensity_formula, data, weights_column,
  outcome_family, verbose, datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

- ...:

  Retained compatibility arguments; not evaluated.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires signed treatment/outcome binding, a complete
propensity and outcome workflow, bounded weights and contributions,
explicit estimand and identification assumptions, and validated
uncertainty.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
