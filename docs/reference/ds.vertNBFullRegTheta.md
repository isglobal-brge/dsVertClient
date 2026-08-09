# Quarantined negative-binomial compatibility frontdoor

This exported name is retained for API compatibility. It raises a typed
`dsvert_route_unavailable` condition before any DSI call and computes or
returns no fit or dispersion estimate. Retained implementation code
after the gate is unreachable through this public frontdoor and carries
no disclosure, DP, accuracy, or availability claim.

## Usage

``` r
ds.vertNBFullRegTheta(
  formula,
  data = NULL,
  theta = NULL,
  joint = TRUE,
  theta_max_iter = 5L,
  theta_tol = 0.001,
  variant = "full_reg_nd",
  beta_max_iter = 2L,
  beta_tol = 1e-04,
  compute_covariance = TRUE,
  verbose = TRUE,
  datasources = NULL,
  ...
)
```

## Arguments

- formula, data, theta, joint, theta_max_iter, theta_tol, variant,
  beta_max_iter, beta_tol, compute_covariance, verbose, datasources:

  Retained compatibility arguments. They are not evaluated because the
  public frontdoor fails locally.

- ...:

  Retained compatibility arguments; not evaluated.

## Value

No fitted object. The function raises `dsvert_route_unavailable` before
DSI.

## Details

Promotion requires a bounded NB2 score/information capsule, joint
beta/theta inference, certified deviance and covariance, and independent
multi-host validation.

## See also

[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
