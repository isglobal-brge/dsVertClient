# dsVertClient: DataSHIELD Client Functions for Vertically Partitioned Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**dsVertClient** is a client-side DataSHIELD package that enables privacy-preserving statistical analysis on **vertically partitioned federated data**. In vertical partitioning, different data sources hold different variables (columns) for the same set of observations (rows).

This package provides user-friendly functions for:

- **ECDH-PSI Record Alignment**: Privacy-preserving record matching using P-256 elliptic curves (no dictionary attacks possible)
- **Correlation Analysis**: Compute correlation matrices using **Multiparty Homomorphic Encryption (MHE)** with share-wrapped threshold decryption
- **Principal Component Analysis**: Perform PCA from the MHE-based correlation matrix
- **Generalized Linear Models**: Fit GLMs using IRLS with encrypted-label BCD and GLM secure routing (3 families)

## Security Guarantees

### Record Alignment (ECDH-PSI with Blind Relay)

| Property | What it prevents |
|----------|-----------------|
| **Dictionary attack resistance** | Unlike SHA-256 hashing, the client CANNOT reverse-engineer patient IDs from the masked curve points (requires server's secret scalar) |
| **Scalar confidentiality** | Each server's P-256 random scalar stays on-server |
| **Unlinkability (DDH)** | The client CANNOT link single-masked points across servers |
| **Blind relay** | Client sees ONLY opaque encrypted blobs (X25519 + AES-256-GCM). All EC point exchanges are transport-encrypted server-to-server. |
| **PSI Firewall** | Server-side FSM enforces phase ordering and one-shot semantics per target, preventing OPRF oracle attacks |
| **MITM-resistant mode** | Optional pre-shared key pinning validates transport PKs against administrator-configured peers, detecting key substitution |

### Statistical Analysis (MHE-CKKS)

| Property | What it prevents |
|----------|-----------------|
| **Client privacy** | The researcher CANNOT decrypt any individual data. Only aggregate statistics (correlations, coefficients) are revealed after all servers cooperate. |
| **Server privacy** | Each server's raw data never leaves the server. Other servers only see encrypted ciphertexts (opaque). |
| **Collusion resistance** | Even K-1 colluding servers cannot decrypt without the K-th server's key share. Full decryption requires ALL K servers. |

### Transport & Routing Security

| Property | What it prevents |
|----------|-----------------|
| **Share-wrapping** | Client CANNOT see or manipulate partial decryption shares. Each share is transport-encrypted under the fusion server's X25519 key before relay. |
| **GLM secure routing** | Client CANNOT see individual-level vectors (eta, mu, w, v). Only sees p_k-length coefficient vectors and opaque encrypted blobs. |
| **Protocol Firewall** | Only ciphertexts produced by authorized operations can be decrypted. Prevents decryption oracle attacks. |

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

# 3. Compute correlation matrix (MHE with share-wrapped threshold decryption)
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
| `ds.psiAlign` | Align records across servers using ECDH-PSI |
| `ds.vertCor` | Compute cross-server correlation matrix via MHE share-wrapped threshold decryption |
| `ds.vertPCA` | Perform PCA from the MHE-based correlation matrix |
| `ds.vertGLM` | Fit Generalized Linear Models via encrypted-label BCD with GLM secure routing or HE-Link |

## Workflow: Align -> Analyze

```
 ┌──────────────────┐   ┌───────────────────────────┐
 │  ds.psiAlign     │──▶│        Analyze            │
 │  (ECDH-PSI)      │   │  ds.vertCor (MHE)         │
 │                  │   │  ds.vertPCA               │
 │                  │   │  ds.vertGLM (BCD + MHE)   │
 └──────────────────┘   └───────────────────────────┘
```

## MHE Correlation: How It Works

The 6-phase threshold MHE protocol:

1. **Key Generation**: Each server generates its own secret key share and public key share. Party 0 creates the Common Reference Polynomial (CRP).
2. **Key Combination**: Public key shares are combined into a Collective Public Key (CPK). Encryption under the CPK requires ALL servers for decryption.
3. **Encryption**: Each server standardizes its data (Z-scores) and encrypts columns under the CPK.
4. **Local Correlation**: Within-server correlations are computed in plaintext.
5. **Cross-Server Correlation (share-wrapped threshold decryption)**: For each server pair (A, B): server A computes Z_A * Enc(Z_B) homomorphically. Each server produces a partial decryption share, **wrapped** (transport-encrypted) under the fusion server's X25519 public key. The client relays these opaque wrapped shares to the fusion server (party 0). The fusion server unwraps all shares, computes its own share, and returns only the final aggregate scalar. The client **never** sees raw partial shares.
6. **Assembly**: The full p x p correlation matrix is assembled from local and cross-server blocks.

## GLM Secure Routing (Coordinator Model)

The GLM protocol uses a **coordinator model** so that individual-level vectors never pass through the client:

| Role | Responsibilities |
|------|-----------------|
| **Label server (coordinator)** | Runs the IRLS loop. Encrypts individual-level vectors (mu, w, v) under each non-label server's public key and sends the ciphertexts through the client as opaque blobs. |
| **Non-label servers** | Decrypt locally, compute their gradient contribution (X_k' W v), encrypt eta_k for the coordinator. |
| **Client (blind relay)** | Routes opaque encrypted blobs between servers. Sees only p_k-length coefficient vectors (public model output) and encrypted payloads it cannot decrypt. |

This ensures the client acts as a transport layer with no access to observation-level information.

## HE-Link: Homomorphic Sigmoid for K=2

With K=2 servers, standard secure routing leaks the non-label server's `eta_nonlabel` to the label server (because `eta_total = eta_label + eta_nonlabel`). The `eta_privacy = "he_link"` mode fixes this by computing `mu = sigmoid(eta_total)` entirely in the encrypted domain using CKKS polynomial approximation:

- Each server encrypts `eta_k` under the CPK
- Coordinator evaluates degree-7 sigmoid polynomial homomorphically → `ct_mu`
- Gradient computed from `ct_y - ct_mu` (both encrypted)
- Only p_k-length gradient scalars are threshold-decrypted

## Secure Aggregation: Pairwise PRG Masks for K>=3

With K>=3 servers (>=2 non-label), the `eta_privacy = "secure_agg"` mode uses **pairwise PRG masks** so the coordinator only sees the aggregate `sum(eta_k)`, never individual per-server contributions. Masks are generated via ChaCha20 PRG from shared seeds (X25519 ECDH + HKDF) and cancel exactly thanks to fixed-point integer arithmetic. A server-side **FSM firewall** prevents protocol abuse (out-of-order calls, iteration replay, partial aggregation).

### Ring Topology

For K>=4 servers, the `topology = "ring"` option reduces seed derivations from O(K-1) to O(2) per server. Instead of all-pairs seed derivation, servers are arranged in a sorted circular order and each server only derives seeds with its two ring neighbors. For K=3, ring and pairwise topologies are identical.

### Peer Manifest Consensus

When `use_secure_agg = TRUE`, the client automatically runs a **peer manifest consensus** protocol (Phase 0.5) before the BCD loop. Each server computes `SHA-256(canonical_manifest_json)` and transport-encrypts its hash to every peer. Peers decrypt and compare. Any mismatch (e.g., phantom peer injection by a malicious client) triggers an immediate stop. This is enabled by default in Phase 0.5 when secure aggregation is active; server-side enforcement via `glmSecureAggInitDS` is controlled by the `dsvert.manifest_consensus` option.

### K=2 Gaussian Warning

When `eta_privacy = "auto"` and K=2 with Gaussian family, the auto-selection picks `"transport"` mode where the coordinator sees the non-label server's eta vector directly. A warning is emitted to alert users to this privacy limitation. Set `eta_privacy = "transport"` explicitly to suppress the warning.

## eta_privacy Modes

The `eta_privacy` parameter controls which privacy protection is used for the linear predictor exchange:

- `"auto"` (default): selects the strongest available mode:
  - K>=3 servers → `"secure_agg"` (pairwise PRG masks)
  - K=2 servers + binomial/poisson → `"he_link"` (homomorphic sigmoid)
  - K=2 servers + gaussian → `"transport"` (standard secure routing)
- `"secure_agg"`: forces pairwise PRG mask aggregation (requires K>=3)
- `"he_link"`: forces HE-Link mode (requires K=2, binomial, `log_n >= 14`)
- `"transport"`: always uses standard secure routing (eta visible to coordinator)

## Supported GLM Families

| Family | Link | Use Case |
|--------|------|----------|
| `gaussian` | Identity | Continuous outcomes (linear regression) |
| `binomial` | Logit | Binary outcomes (logistic regression) |
| `poisson` | Log | Count data |

## Session Isolation

Each call to `ds.vertGLM`, `ds.psiAlign`, or `ds.vertCor` automatically
generates a UUIDv4 session identifier and passes it to every server-side
function call. This ensures that multiple concurrent analyses use completely
independent server-side storage --- keys, ciphertexts, and protocol state
never collide between parallel jobs. Session cleanup happens automatically
at the end of each analysis (or via 24-hour TTL expiry for failed jobs).

## Chunked Transfer Protocol

DataSHIELD's R expression parser limits the size of arguments that can be
passed inline in `call()` expressions. Cryptographic objects (CKKS
ciphertexts, EC points, key shares, transport-encrypted blobs) routinely
exceed this limit. dsVertClient uses **adaptive chunking** (default 200 KB,
configurable via `options(dsvert.chunk_size = N)`) that automatically detects
each server's expression size limit and reduces chunk size by 25% on failure.
Payloads are sent via `mheStoreBlobDS`, which the server reassembles
transparently. The working chunk size is cached for the session and converges
to the minimum across all connected servers. This chunking is applied
uniformly across all protocols (PSI, MHE correlation, GLM) for any data that
scales with the number of observations, variables, or servers --- ensuring
that no DataSHIELD call can overflow regardless of dataset size or number of
parties.

## Performance Optimizations

dsVertClient implements several optimizations to reduce communication overhead:

- **Parallel DSI dispatch**: Key generation, CPK distribution, transport key exchange, PSI init, index collection, and cleanup are dispatched to all servers simultaneously using DSI's async aggregate calls.
- **Batched threshold decryption**: All ciphertexts for a cross-server correlation block are sent as blobs and processed in a single server-side call, reducing client-server round-trips from O(p_A × p_B × K) to O(K).
- **MHE context caching**: When `reuse_mhe = TRUE`, consecutive analyses with the same peer set and CKKS parameters reuse cached encryption keys, skipping the expensive key-generation and key-combination phases. Forward secrecy is preserved via fresh X25519 transport keys per job.
- **Binary wire formats**: PSI points use compressed P-256 binary packing (~66% smaller). GLM transport vectors use little-endian IEEE 754 binary format (~2× smaller than JSON).
- **Adaptive chunking**: Automatic chunk size probing and fallback for large payloads (default 200 KB, cached per session).

## Requirements

- R >= 4.0.0
- DSI package
- jsonlite package
- dsVert package installed on DataSHIELD servers (Opal/Rock)

## Documentation

- [Validation Study](https://isglobal-brge.github.io/dsVertClient/articles/validation.html): Full pipeline validation against centralized R

## Authors

- David Sarrat González (david.sarrat@isglobal.org)
- Miron Banjac (miron.banjac@isglobal.org)
- Juan R González (juanr.gonzalez@isglobal.org)
