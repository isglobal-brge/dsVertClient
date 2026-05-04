# Federated Cox proportional-hazards regression

User-facing Cox PH wrapper. Dispatches to
[`ds.vertCoxProfileNonDisclosive`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCoxProfileNonDisclosive.md),
the non-disclosive Breslow profile route validated for K=2 and K\>=3.
The historical rank/permutation and person-time Poisson Cox routes were
removed from the product package because they exposed event-rank
metadata or required person-time-expanded inputs.

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

- formula:

  Formula of the form `Surv(time, event) ~ x1 + ...`.

- data:

  Aligned data-frame name on each server.

- max_iter:

  Maximum Newton iterations for the profile route.

- tol:

  Convergence tolerance on max \|delta beta\|.

- max_event_times:

  Integer runtime guard passed to
  [`ds.vertCoxProfileNonDisclosive`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCoxProfileNonDisclosive.md).

- newton, ridge_eps, debug_trace:

  Parameters passed through to
  [`ds.vertCoxProfileNonDisclosive`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCoxProfileNonDisclosive.md).

- verbose:

  Print progress.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertCox` object: `coefficients`, `std_errors`, `covariance`,
`loglik` (partial log-likelihood), `n_obs`, `n_events`, `iterations`,
`converged`.
