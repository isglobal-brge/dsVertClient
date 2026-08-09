# Quarantined iterative-LASSO compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no coefficient path, score, Hessian, selection result, or
diagnostic. Retained iterative code after the gate is unreachable
through this public frontdoor and carries no disclosure, DP, accuracy,
model-selection, or availability claim.

## Usage

``` r
ds.vertLASSOIter(
  formula,
  data = NULL,
  family = c("gaussian", "binomial", "poisson"),
  lambda = NULL,
  max_outer = 20L,
  tol = 1e-08,
  alpha = 0.5,
  inner_iter = 8L,
  exact_non_gaussian = TRUE,
  ring = NULL,
  lipschitz = c("auto", "gram", "safe"),
  fista_restart = TRUE,
  binomial_sigmoid_intervals = getOption("dsvert.lasso_binomial_sigmoid_intervals", 200L),
  poisson_damping = 0.5,
  verbose = TRUE,
  datasources = NULL,
  cor_analysis_id = NULL
)
```

## Arguments

- formula, data, family, lambda, max_outer, tol, alpha, inner_iter,
  exact_non_gaussian, ring, lipschitz, fista_restart,
  binomial_sigmoid_intervals, poisson_damping, verbose, datasources,
  cor_analysis_id:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires a signed bounded artifact over the exact score
design, a whole-path privacy contract, KKT validation and DP-aware
selection and inference.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
