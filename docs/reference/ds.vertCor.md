# Correlation Matrix for Vertically Partitioned Data

Client-side function that computes a correlation matrix across
vertically partitioned data using distributed Block SVD.

## Usage

``` r
ds.vertCor(data_name, variables, datasources = NULL)
```

## Arguments

- data_name:

  Character string. Name of the (aligned) data frame on each server.

- variables:

  A named list where each name corresponds to a server name and each
  element is a character vector of variable names from that server.

- datasources:

  DataSHIELD connection object or list of connections. If NULL, uses all
  available connections.

## Value

A correlation matrix with dimensions equal to the total number of
variables across all servers. Row and column names are the variable
names.

## Details

This function implements distributed correlation analysis using Block
Singular Value Decomposition:

1.  Each server computes U\*D from SVD of its standardized variables

2.  Client combines into a single matrix

3.  Client computes final SVD to get V and D

4.  Correlation matrix = V \* D^2 \* V' (normalized)

This method is privacy-preserving because:

- Individual U\*D matrices cannot reconstruct original data

- Only the combined correlation structure is revealed

## References

Iwen, M. & Ong, B.W. (2016). A distributed and incremental SVD algorithm
for agglomerative data analysis on large networks. SIAM Journal on
Matrix Analysis and Applications.

## See also

[`ds.vertPCA`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertPCA.md)
for principal component analysis

## Examples

``` r
if (FALSE) { # \dontrun{
# Define which variables are on which server
vars <- list(
  server1 = c("age", "weight"),
  server2 = c("height", "bmi"),
  server3 = c("glucose", "cholesterol")
)

# Compute correlation matrix
cor_matrix <- ds.vertCor("D_aligned", vars, datasources = conns)
print(cor_matrix)
} # }
```
