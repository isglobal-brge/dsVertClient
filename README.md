# dsVertClient — DataSHIELD Client for Vertically Partitioned Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**dsVertClient** provides user-friendly R functions for privacy-preserving analysis on vertically partitioned data across DataSHIELD servers. The analyst calls simple functions; all cryptographic protocols run transparently.

## Quick Start

```r
library(DSI); library(DSOpal); library(dsVertClient)

# Connect to servers
conns <- DSI::datashield.login(logins = builder$build())
DSI::datashield.assign.table(conns, "D", "project.data")

# 1. Align records (find common patients across servers)
ds.psiAlign("D", "patient_id", "D_aligned", datasources = conns)

# 2. Fit a logistic regression
result <- ds.vertGLM(
  data_name = "D_aligned",
  y_var     = "diabetes",
  x_vars    = list(
    hospital_A = c("age", "bmi"),
    hospital_B = c("glu", "bp")
  ),
  y_server  = "hospital_B",
  family    = "binomial",
  datasources = conns
)

# 3. View results (coefficients + SE + p-values)
print(result)

# 4. Correlation matrix
cor <- ds.vertCor("D_aligned", variables = list(...), datasources = conns)

# 5. PCA
pca <- ds.vertPCA("D_aligned", variables = list(...), datasources = conns)

DSI::datashield.logout(conns)
```

## Functions

| Function | Description | Output |
|----------|-------------|--------|
| `ds.psiAlign()` | Private Set Intersection | Aligned data frame |
| `ds.vertGLM()` | GLM (Gaussian, Binomial, Poisson) | Coefficients, SE, p-values, deviance |
| `ds.vertCor()` | Pearson correlation matrix | p x p correlation matrix |
| `ds.vertPCA()` | Principal Component Analysis | Loadings, eigenvalues, variance % |

## GLM Output

```
Coefficients:
                Estimate Std.Error z value  Pr(>|z|)
(Intercept)     -9.2784    2.0981  -4.425  0.000010 ***
age              0.0733    0.0223   3.301  0.000993 ***
bmi              0.0798    0.0522   1.534  0.126024
glu              0.0290    0.0101   2.878  0.004046 **

Converged: TRUE (14 iterations)
Deviance: 22.09
```

## Supported Configurations

| | K=2 | K>=3 |
|--|-----|------|
| Gaussian | L-BFGS + identity link | L-BFGS + identity link |
| Binomial | L-BFGS + DCF sigmoid (50 intervals) | L-BFGS + DCF sigmoid |
| Poisson | L-BFGS + DCF exp (100 intervals) | L-BFGS + DCF exp |
| SE + p-values | Finite-difference Hessian | Central-difference Hessian |
| Deviance | Beaver dot-product (1 scalar) | Beaver dot-product |
| Correlation | — | Ring63 Beaver (p matvec) |
| PCA | — | Eigen of correlation |

## Security

- **Zero observation-level disclosure**: client sees only p-dimensional aggregates
- **Server-generated Beaver triples**: client never sees cryptographic material
- **Dealer rotation**: different server generates triples each iteration (K>=4)
- **Transport encryption**: X25519 + AES-256-GCM between servers
- **No CKKS/Lattigo**: pure Ring63 fixed-point MPC

## Performance (Pima diabetes, p=6)

| Analysis | K=2 | K=3 |
|----------|-----|-----|
| Key setup | 0.3s | 0.5s |
| Binomial GLM | ~330s (12 iters + SE) | ~360s (14 iters + SE) |
| Gaussian GLM | ~60s (11 iters + SE) | ~60s (11 iters + SE) |
| Poisson GLM | ~580s (11 iters + SE) | ~580s (12 iters + SE) |
| Correlation | — | 18s |
| PCA | — | 18s |

## Installation

```r
# From GitHub
devtools::install_github("isglobal-brge/dsVertClient")

# Or from source
R CMD build --no-build-vignettes .
R CMD INSTALL dsVertClient_2.0.0.tar.gz
```

Requires the server-side package [dsVert](https://github.com/isglobal-brge/dsVert) installed on all DataSHIELD servers.
