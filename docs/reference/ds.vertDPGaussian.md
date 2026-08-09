# Bounded Gaussian regression from the sticky DP capsule

Fits a descriptive Gaussian linear model from one signed sufficient-
statistic artifact. Outcome and predictors are clipped to
custodian-owned public bounds, collapsed once per admitted patient,
normalized to `[0,1]`, and restricted to a fixed complete-case estimand.
The function performs one sticky capsule retrieval and then only
deterministic client post-processing.

## Usage

``` r
ds.vertDPGaussian(
  data_name,
  analysis_id,
  ridge = 0,
  server = NULL,
  datasources = NULL
)
```

## Arguments

- data_name:

  Signed protected dataset name.

- analysis_id:

  Custodian-configured signed Gaussian artifact id.

- ridge:

  Explicit non-negative ridge penalty on normalized predictors. The
  default zero preserves the unpenalized bounded-statistic estimand.

- server:

  Optional expected signed outcome-owner server name.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertDPGaussian` object. It contains no classical standard errors,
p-values, individual fitted values, residuals, or scores.
