# dsVertClient: DataSHIELD Client Functions for Vertically Partitioned Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**dsVertClient** is a client-side DataSHIELD package that enables privacy-preserving statistical analysis on **vertically partitioned federated data**. In vertical partitioning, different data sources hold different variables (columns) for the same set of observations (rows).

This package provides user-friendly functions for:

- **ECDH-PSI Record Alignment**: Privacy-preserving record matching using P-256 elliptic curves (no dictionary attacks possible)
- **Correlation Analysis**: Compute correlation matrices using **Multiparty Homomorphic Encryption (MHE)** with threshold decryption
- **Principal Component Analysis**: Perform PCA from the MHE-based correlation matrix
- **Generalized Linear Models**: Fit GLMs using Block Coordinate Descent with encrypted labels (5 families)

## Security Guarantees

### Record Alignment (ECDH-PSI)

| Property | What it prevents |
|----------|-----------------|
| **Dictionary attack resistance** | Unlike SHA-256 hashing, the client CANNOT reverse-engineer patient IDs from the masked curve points (requires server's secret scalar) |
| **Scalar confidentiality** | Each server's P-256 random scalar stays on-server |
| **Unlinkability (DDH)** | The client CANNOT link single-masked points across servers |

### Statistical Analysis (MHE-CKKS)

| Property | What it prevents |
|----------|-----------------|
| **Client privacy** | The researcher CANNOT decrypt any individual data. Only aggregate statistics (correlations, coefficients) are revealed after all servers cooperate. |
| **Server privacy** | Each server's raw data never leaves the server. Other servers only see encrypted ciphertexts (opaque). |
| **Collusion resistance** | Even K-1 colluding servers cannot decrypt without the K-th server's key share. Full decryption requires ALL K servers. |

## Installation

```r
# Install from GitHub (install dsVert first on servers)
devtools::install_github("isglobal-brge/dsVertClient")
```

## Quick Start

```r
library(dsVertClient)
library(DSI)

# Connect to Opal/DataSHIELD servers
conns <- datashield.login(logindata)

# 1. Align records across servers using ECDH-PSI
ds.psiAlign("D", "patient_id", "D_aligned", datasources = conns)

# 2. Define which variables are on which server
variables <- list(
  hospital_A = c("age", "bmi"),
  hospital_B = c("glucose", "systolic_bp"),
  hospital_C = c("cholesterol", "hdl")
)

# 3. Compute correlation matrix (MHE with threshold decryption)
cor_result <- ds.vertCor("D_aligned", variables, datasources = conns)
print(cor_result)

# 4. Perform PCA
pca <- ds.vertPCA("D_aligned", variables, n_components = 3, datasources = conns)
print(pca)

# 5. Fit a GLM
model <- ds.vertGLM("D_aligned", "outcome", variables,
                    y_server = "hospital_B", family = "gaussian",
                    datasources = conns)
summary(model)

# Disconnect
datashield.logout(conns)
```

## Client-Side Functions

| Function | Description |
|----------|-------------|
| `ds.psiAlign` | Align records across servers using ECDH-PSI (recommended) |
| `ds.vertCor` | Compute cross-server correlation matrix via MHE threshold decryption |
| `ds.vertPCA` | Perform PCA from the MHE-based correlation matrix |
| `ds.vertGLM` | Fit Generalized Linear Models via encrypted-label BCD |
| `ds.validateIdFormat` | Validate identifier format consistency across servers |
| `ds.hashId` | Get hashed identifiers (deprecated, use `ds.psiAlign`) |
| `ds.alignRecords` | Hash-based alignment (deprecated, use `ds.psiAlign`) |

## Workflow: Align → Analyze

```
 ┌──────────────────┐   ┌───────────────────────────┐
 │  ds.psiAlign     │──▶│        Analyze            │
 │  (ECDH-PSI)      │   │  ds.vertCor (MHE)         │
 │                   │   │  ds.vertPCA               │
 │                   │   │  ds.vertGLM (BCD + MHE)   │
 └──────────────────┘   └───────────────────────────┘
```

## MHE Correlation: How It Works

The 6-phase threshold MHE protocol:

1. **Key Generation**: Each server generates its own secret key share and public key share. Party 0 creates the Common Reference Polynomial (CRP).
2. **Key Combination**: Public key shares are combined into a Collective Public Key (CPK). Encryption under the CPK requires ALL servers for decryption.
3. **Encryption**: Each server standardizes its data (Z-scores) and encrypts columns under the CPK.
4. **Local Correlation**: Within-server correlations are computed in plaintext.
5. **Cross-Server Correlation**: For each server pair (A, B): server A computes Z_A * Enc(Z_B) homomorphically. Each server provides a partial decryption share. The client fuses all shares to recover the inner product.
6. **Assembly**: The full p x p correlation matrix is assembled from local and cross-server blocks.

## Supported GLM Families

| Family | Link | Use Case |
|--------|------|----------|
| `gaussian` | Identity | Continuous outcomes (linear regression) |
| `binomial` | Logit | Binary outcomes (logistic regression) |
| `poisson` | Log | Count data |
| `Gamma` | Log | Positive continuous data (costs, times) |
| `inverse.gaussian` | Log | Positive continuous with high variance |

## Requirements

- R >= 4.0.0
- DSI package
- jsonlite package
- dsVert package installed on DataSHIELD servers (Opal/Rock)

## Documentation

- [Getting Started](https://isglobal-brge.github.io/dsVertClient/articles/a-getting-started.html): Setup, PSI record alignment, first analysis
- [Statistical Analysis](https://isglobal-brge.github.io/dsVertClient/articles/b-statistical-analysis.html): Correlation, PCA, and GLMs
- [Methodology](https://isglobal-brge.github.io/dsVertClient/articles/c-methodology.html): MHE protocol, BCD-IRLS, ECDH-PSI details
- [Validation](https://isglobal-brge.github.io/dsVertClient/articles/d-validation.html): Numerical comparison against local R
- [Security Model](https://isglobal-brge.github.io/dsVertClient/articles/e-security.html): Threat model, MHE guarantees, PSI security

## Authors

- David Sarrat Gonzalez
- Miron Banjac
- Juan R Gonzalez

## References

- De Cristofaro, E. & Tsudik, G. (2010). "Practical Private Set Intersection Protocols with Linear Complexity". *FC 2010*.
- Mouchet, C. et al. (2021). "Multiparty Homomorphic Encryption from Ring-Learning-With-Errors". *Proceedings on Privacy Enhancing Technologies (PETS)*.
- Cheon, J.H. et al. (2017). "Homomorphic Encryption for Arithmetic of Approximate Numbers". *ASIACRYPT 2017*.
- van Kesteren, E.J. et al. (2019). "Privacy-preserving generalized linear models using distributed block coordinate descent". arXiv:1911.03183.
- Lattigo v6: https://github.com/tuneinsight/lattigo
