# Principal components from a signed DP correlation artifact

Performs client-only eigen decomposition of the explicitly PSD-projected
complete-case matrix returned by `ds.vertCor`. It never accepts the
pairwise-complete correlation artifact, invokes the former Ring63/exact
correlation protocol, or returns individual component scores.

## Usage

``` r
ds.vertPCA(
  data_name = NULL,
  variables = NULL,
  n_components = NULL,
  analysis_id = NULL,
  cor_result = NULL,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- data_name:

  Signed protected dataset name. Ignored when `cor_result` is supplied.

- variables:

  Optional signed variable subset. Named owner lists are supported when
  they agree with the Gaussian artifact.

- n_components:

  Number of components, or `NULL` for all.

- analysis_id:

  Mandatory signed Gaussian artifact id when `cor_result` is not
  supplied.

- cor_result:

  An existing complete-case `ds.vertCor` result from the same sticky
  capsule. Pairwise, arbitrary, and legacy `ds.cor` objects are
  rejected.

- verbose:

  Logical progress flag.

- datasources:

  DataSHIELD connections.

## Value

A `ds.pca` object with loadings, eigenvalues, explicit spectral
mechanism diagnostics and inherited DP provenance. Scores are
unavailable.
