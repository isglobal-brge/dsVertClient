# Differentially private same-owner Pearson correlation

Reads pairwise-complete bounded moments from exactly one signed, sticky
biomedical capsule vector. The requested variables must all be
co-located in one dataset and owner. Cross-owner products remain
deliberately unavailable. The returned raw matrix uses pairwise-complete
DP moments and can therefore be indefinite; `correlation` is an
explicitly labelled eigenvalue-clipping and diagonal-rescaling
projection, not an exact correlation reconstruction or an exact
nearest-correlation solution.

## Usage

``` r
ds.vertDPCor(
  data_name,
  analysis_id,
  variables = NULL,
  server = NULL,
  datasources = NULL
)
```

## Arguments

- data_name:

  Signed protected dataset name.

- analysis_id:

  Signed correlation artifact id, normally
  `paste(data_name, owner, sep = "::")`.

- variables:

  Optional character subset, or a named list containing exactly one
  owner entry. At least two variables are required.

- server:

  Optional assertion of the single owner.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertDPCor`/`ds.cor` object containing raw and projected
correlations, pair counts, simultaneous mechanism/quantization regions,
and capsule provenance.
