# Disclosure-safe compatibility adapter for correlation

This existing name consumes only a signed `gaussian_models` sufficient-
statistic artifact. A signed `analysis_id` is mandatory. Same-owner and
cross-owner variables share one secret complete-case mask; there is no
silent fallback to pairwise moments.

## Usage

``` r
ds.vertCor(
  data_name,
  variables = NULL,
  analysis_id = NULL,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- data_name:

  Signed protected dataset name.

- variables:

  Optional character subset, or a named list identifying the signed
  owner of each requested variable. At least two variables are required.

- analysis_id:

  Mandatory signed Gaussian artifact id.

- verbose:

  Logical progress flag retained for compatibility.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertDPCor`/`ds.cor` object.
