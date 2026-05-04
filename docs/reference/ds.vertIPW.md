# Federated inverse-probability-weighted GLM (two-stage)

Convenience wrapper implementing the two-stage IPW (inverse-probability
weighting) workflow on vertically partitioned DataSHIELD data.

Stage 1: fit the propensity model via `ds.vertGLM` with a binomial
family on the supplied `propensity_formula`. Stage 2: fit the outcome
model via `ds.vertGLM` with user-supplied `weights_column` (already
holding IPW weights on the server that owns the treatment / outcome
variables).

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

- outcome_formula:

  Formula for the outcome model.

- propensity_formula:

  Formula for the propensity (binary treatment) model.

- data:

  Name of the aligned data frame on each server.

- weights_column:

  Name of the column on the outcome server that holds 1/P(T=1\|X) for
  treated units and 1/(1-P(T=1\|X)) for controls, or the stabilised
  analogue. Pre-computed server-side (the dsVert pipeline for ON-SHARE
  propensity score computation is a planned follow-on; see
  V2_PROGRESS.md).

- outcome_family:

  Family for the outcome model. Default "gaussian".

- verbose:

  Logical. Print stage-by-stage progress (default TRUE).

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

- ...:

  Passed through to both underlying `ds.vertGLM` calls.

## Value

list of class `ds.vertIPW` with `propensity` and `outcome` ds.glm
objects.
