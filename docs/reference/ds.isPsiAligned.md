# Check whether a data symbol is pinned-padded-PSI aligned

Each server validates a persistent, secret-token-authenticated manifest
against the current identifier values and row order. The client accepts
the result only when every named server returns the same fixed-schema
public attestation. No row count, identifier commitment or intersection
size is released.

## Usage

``` r
ds.isPsiAligned(newobj = "DA", datasources = NULL)
```

## Arguments

- newobj:

  Character. Symbol to validate on every server.

- datasources:

  Named DataSHIELD connections. When `NULL`, the current connections
  returned by
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html)
  are used.

## Value

A list with `aligned` and `n_common`. `n_common` is always
`NA_integer_`; the alignment-status route never releases cardinality.
