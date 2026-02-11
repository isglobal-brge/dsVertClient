# Principal Component Analysis for Vertically Partitioned Data

Performs PCA on vertically partitioned data using the privacy-preserving
correlation matrix computed via Homomorphic Encryption.

## Usage

``` r
ds.vertPCA(
  data_name = NULL,
  variables = NULL,
  n_components = NULL,
  cor_result = NULL,
  log_n = 12,
  log_scale = 40,
  datasources = NULL
)
```

## Arguments

- data_name:

  Character string. Name of the (aligned) data frame on each server.
  Ignored if `cor_result` is provided.

- variables:

  A named list where each name corresponds to a server name and each
  element is a character vector of variable names from that server.
  Ignored if `cor_result` is provided.

- n_components:

  Integer. Number of principal components to return. Default is NULL
  (returns all).

- cor_result:

  An existing `ds.cor` object from
  [`ds.vertCor`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md).
  If provided, the MHE protocol is NOT re-run; the PCA is computed
  directly from this correlation matrix. This avoids running the
  expensive MHE protocol twice when you already have the correlation.

- log_n:

  Integer. CKKS ring dimension parameter for MHE. Default is 12. Ignored
  if `cor_result` is provided.

- log_scale:

  Integer. CKKS precision parameter for MHE. Default is 40. Ignored if
  `cor_result` is provided.

- datasources:

  DataSHIELD connection object or list of connections. If NULL, uses all
  available connections. Ignored if `cor_result` is provided.

## Value

A list with class "ds.pca" containing:

- `loadings`: Matrix of variable loadings (n_vars x n_components)

- `eigenvalues`: Eigenvalues for each component

- `variance_pct`: Percentage of variance explained

- `cumulative_pct`: Cumulative percentage explained

- `var_names`: Variable names

- `n_obs`: Number of observations

- `correlation`: The correlation matrix used for PCA

## Details

This function performs PCA using the correlation matrix obtained via
Multiparty Homomorphic Encryption (MHE). The approach is:

1.  Compute the privacy-preserving correlation matrix using
    [`ds.vertCor`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md)
    (or reuse an existing one via the `cor_result` parameter)

2.  Perform eigen decomposition on the correlation matrix

3.  Extract loadings (eigenvectors) and eigenvalues

Since PCA on standardized data is equivalent to eigen decomposition of
the correlation matrix, this gives correct loadings and variance
explained.

**Note on scores:** This function does NOT return principal component
scores because computing scores would require access to the raw data. If
you need scores, you would need to compute them on each server using the
loadings and aggregate the results (which is a separate operation).

### Interpreting Loadings

Each column of the loadings matrix represents a principal component. The
values show how much each variable contributes to that component:

- Values close to 1 or -1 indicate strong contribution

- Values close to 0 indicate weak contribution

- Sign indicates direction of relationship

## Security

This function inherits all security properties from
[`ds.vertCor`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md):

- Individual observations are never exposed

- The client cannot decrypt without all servers cooperating

- Only aggregate statistics (correlation matrix) are revealed

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

# Or reuse an existing correlation result (avoids running MHE again):
cor_result <- ds.vertCor("D_aligned", vars)
pca_result <- ds.vertPCA(cor_result = cor_result, n_components = 3)

# View variance explained
print(pca_result)

# Biplot of loadings for first two PCs
plot(pca_result$loadings[, 1], pca_result$loadings[, 2],
     xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = paste0("PC1 (", round(pca_result$variance_pct[1], 1), "%)"),
     ylab = paste0("PC2 (", round(pca_result$variance_pct[2], 1), "%)"))
text(pca_result$loadings[, 1], pca_result$loadings[, 2],
     labels = pca_result$var_names, pos = 3)
abline(h = 0, v = 0, lty = 2)
} # }
```
