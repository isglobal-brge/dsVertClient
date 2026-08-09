# Preflight a vertically partitioned GLM numeric backend

Negotiates only server-published bounds and runtime capabilities. The
analyst may request a preferred backend/ring, but cannot supply or relax
a bound. A stronger backend is selected automatically when required.

## Usage

``` r
ds.vertNumericPreflight(
  n_obs,
  n_predictors,
  family = "gaussian",
  max_iter = 100L,
  compute_se = TRUE,
  compute_deviance = TRUE,
  weights_active = FALSE,
  offset_active = FALSE,
  backend = "auto",
  ring = 63L,
  datasources = NULL
)
```

## Arguments

- n_obs:

  Public observation count.

- n_predictors:

  Number of predictor columns.

- family:

  GLM family.

- max_iter:

  Maximum optimizer iterations.

- compute_se:

  Whether the workload includes Hessian evaluations.

- compute_deviance:

  Whether the workload includes deviance evaluation.

- weights_active:

  Whether protected weights are active.

- offset_active:

  Whether a protected offset is active.

- backend:

  Preferred backend: auto, ring63, ring127, exact_gc, or multiprecision.

- ring:

  Preferred fast ring, 63 or 127.

- datasources:

  DataSHIELD connections.

## Value

A numeric certificate.
