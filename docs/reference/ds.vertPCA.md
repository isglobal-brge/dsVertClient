# Principal Component Analysis for Vertically Partitioned Data

Client-side function that performs PCA on vertically partitioned data
using distributed Block SVD.

## Usage

``` r
ds.vertPCA(data_name, variables, n_components = NULL, datasources = NULL)
```

## Arguments

- data_name:

  Character string. Name of the (aligned) data frame on each server.

- variables:

  A named list where each name corresponds to a server name and each
  element is a character vector of variable names from that server.

- n_components:

  Integer. Number of principal components to return. Default is NULL
  (returns all).

- datasources:

  DataSHIELD connection object or list of connections. If NULL, uses all
  available connections.

## Value

A list with class "ds.pca" containing:

- `scores`: Matrix of principal component scores (n_obs x n_components)

- `loadings`: Matrix of variable loadings (n_vars x n_components)

- `variance`: Variance explained by each component

- `variance_pct`: Percentage of variance explained

- `cumulative_pct`: Cumulative percentage explained

- `var_names`: Variable names

## Details

This function extends the Block SVD approach for correlation to perform
full PCA:

1.  Each server computes U\*D from SVD of its standardized variables

2.  Client combines and performs final SVD

3.  Principal components = U \* D (combined)

4.  Loadings = V (right singular vectors)

5.  Variance = D^2 / (n-1)

## See also

[`ds.vertCor`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md)
for correlation analysis

## Examples

``` r
if (FALSE) { # \dontrun{
vars <- list(
  server1 = c("age", "weight"),
  server2 = c("height", "bmi")
)

pca_result <- ds.vertPCA("D_aligned", vars, n_components = 3)

# Plot first two components
plot(pca_result$scores[, 1], pca_result$scores[, 2],
     xlab = paste0("PC1 (", round(pca_result$variance_pct[1], 1), "%)"),
     ylab = paste0("PC2 (", round(pca_result$variance_pct[2], 1), "%)"))
} # }
```
