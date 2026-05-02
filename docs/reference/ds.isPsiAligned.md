# Check whether a data symbol is already PSI-aligned on every server

Ask each server whether the named symbol has the same row count. Used by
`ds.vertGLM`/LMM/GEE/Cox to SKIP a fresh PSI run when the caller already
aligned the cohort. Security- preserving: only queries shape.

## Usage

``` r
ds.isPsiAligned(newobj = "DA", datasources = NULL)
```

## Arguments

- newobj:

  Character. Symbol name to test on each server (default `"DA"`).

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

## Value

List with `aligned` (logical) and `n_common` (integer or NA_integer\_).
