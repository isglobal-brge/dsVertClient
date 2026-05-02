# Federated chi-square test on a 2-way contingency table

Compute Pearson chi-square (and optional Yates and LR variants) on a
contingency table of two categorical variables held at the same server
in a DataSHIELD federation. The server returns only aggregate cell
counts and margins (with datashield.privacyLevel-aware small-cell
suppression); chi-square statistics are then computed entirely
client-side.

## Usage

``` r
ds.vertChisq(
  data_name,
  var1,
  var2,
  server = NULL,
  correct = TRUE,
  datasources = NULL
)
```

## Arguments

- data_name:

  Character. Name of the aligned data frame on each server.

- var1:

  Character. Row variable.

- var2:

  Character. Column variable.

- server:

  Character. Optional name of the server that holds both variables. If
  `NULL` the server is auto-detected via `dsvertColNamesDS`.

- correct:

  Logical. If TRUE (default for 2x2 tables), apply Yates' continuity
  correction.

- datasources:

  DataSHIELD connections.

## Value

An object of class `ds.vertChisq` with elements `statistic`, `df`,
`p_value`, `observed`, `expected`, `residuals`, `n`, `correct`.

## Details

For the cross-server case (variables vertically partitioned across two
or more servers) a separate entry point `ds.vertChisqCross` will reuse
the Beaver dot-product infrastructure to compute cell counts on secret
shares; that is tracked as a follow-on in the dsVert v2.0 implementation
plan.
