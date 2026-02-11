# dsVertClient - DataSHIELD Client Functions for Vertically Partitioned Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**dsVertClient** is a client-side DataSHIELD package that enables privacy-preserving statistical analysis on **vertically partitioned federated data**. In vertical partitioning, different data sources hold different variables (columns) for the same set of observations (rows).

This package provides user-friendly functions for:
- **Record Matching**: Align records across servers using secure hashing
- **Correlation Analysis**: Compute correlation matrices using distributed Block SVD
- **Principal Component Analysis**: Perform PCA across distributed data
- **Generalized Linear Models**: Fit GLMs using Block Coordinate Descent

## Installation

```r
# Install from GitHub (install dsVert first on servers)
devtools::install_github("isglobal-brge/dsVertClient")
```

## Quick Start

```r
library(dsVertClient)
library(DSI)

# Connect to DataSHIELD servers
conns <- datashield.login(logindata)

# 1. Align records across servers
ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
ds.alignRecords("D", "patient_id", ref_hashes$hashes, "D_aligned", datasources = conns)

# 2. Define which variables are on which server
variables <- list(
  server1 = c("age", "weight"),
  server2 = c("height", "bmi"),
  server3 = c("glucose", "cholesterol")
)

# 3. Compute correlation matrix
cor_matrix <- ds.vertCor("D_aligned", variables, datasources = conns)

# 4. Perform PCA
pca <- ds.vertPCA("D_aligned", variables, n_components = 3, datasources = conns)

# 5. Fit a GLM
model <- ds.vertGLM("D_aligned", "outcome", variables,
                    family = "gaussian", datasources = conns)

# Disconnect
datashield.logout(conns)
```

## Client-Side Functions

| Function | Description |
|----------|-------------|
| `ds.hashId` | Get hashed identifiers from a server |
| `ds.alignRecords` | Align records across servers to match hashes |
| `ds.vertCor` | Compute correlation matrix across servers |
| `ds.vertPCA` | Perform Principal Component Analysis |
| `ds.vertGLM` | Fit Generalized Linear Models (gaussian, binomial, poisson) |

## Record Matching Workflow

```r
# Step 1: Get reference hashes from one server
ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])

# Step 2: Align all servers to match these hashes
ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                newobj = "D_aligned", datasources = conns)

# Now all servers have aligned data in "D_aligned"
```

## Correlation and PCA

The package uses **Block Singular Value Decomposition** for distributed correlation and PCA:

```r
# Correlation
cor_matrix <- ds.vertCor("D_aligned", variables, datasources = conns)

# PCA with variance explained
pca <- ds.vertPCA("D_aligned", variables, n_components = 4, datasources = conns)
print(pca)  # Shows variance explained by each component

# Plot first two components
plot(pca$scores[,1], pca$scores[,2],
     xlab = paste0("PC1 (", round(pca$variance_pct[1], 1), "%)"),
     ylab = paste0("PC2 (", round(pca$variance_pct[2], 1), "%)"))
```

## Generalized Linear Models

Fit GLMs using **Block Coordinate Descent**:

```r
# Linear regression (Gaussian)
model_linear <- ds.vertGLM("D_aligned", "continuous_outcome", x_vars,
                           family = "gaussian", datasources = conns)

# Logistic regression (Binomial)
model_logistic <- ds.vertGLM("D_aligned", "binary_outcome", x_vars,
                             family = "binomial", datasources = conns)

# Poisson regression
model_poisson <- ds.vertGLM("D_aligned", "count_outcome", x_vars,
                            family = "poisson", datasources = conns)

# View results
print(model_linear)
coef(model_linear)
```

## Requirements

- R >= 4.0.0
- DSI package
- dsVert package installed on DataSHIELD servers

## Authors

- David Sarrat González (david.sarrat@isglobal.org)
- Miron Banjac (miron.banjac@isglobal.org)
- Juan R González (juanr.gonzalez@isglobal.org)

## References

- van Kesteren, E.J. et al. (2019). Privacy-preserving generalized linear models using distributed block coordinate descent. arXiv:1911.05935.
- Iwen, M. & Ong, B.W. (2016). A distributed and incremental SVD algorithm for agglomerative data analysis on large networks.
