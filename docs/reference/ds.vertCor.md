# Ring63 Privacy-Preserving Correlation for Vertically Partitioned Data

Computes Pearson correlation matrix using Ring63 Beaver MPC. Only
aggregate scalars disclosed.

## Usage

``` r
ds.vertCor(data_name, variables = NULL, verbose = TRUE, datasources = NULL)
```

## Arguments

- data_name:

  Character. Aligned data frame name.

- variables:

  Named list: server -\> variable names.

- verbose:

  Logical. If TRUE (default), print progress messages.

- datasources:

  DataSHIELD connections.

## Value

List with correlation matrix, variable names, n_obs.
