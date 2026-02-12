# Security Model

## Overview

dsVertClient uses Multiparty Homomorphic Encryption (MHE) based on the
CKKS scheme to protect data during distributed computation. This
vignette describes the threat model, what each participant learns, and
the security guarantees.

------------------------------------------------------------------------

## Threat Model

The protocol assumes a **semi-honest** (honest-but-curious) adversary
model:

- All parties follow the protocol correctly
- Parties may attempt to infer information from messages they receive
- No party intentionally sends corrupted data or deviates from the
  protocol

This is the standard model for federated analytics in healthcare and
research settings, where institutions are contractually obligated to
follow protocols but may be curious about others’ data.

------------------------------------------------------------------------

## What MHE Protects

### Correlation Analysis

During distributed correlation, cross-server products are computed on
encrypted data. Each server encrypts its standardized columns under a
collective public key (CPK). No single server can decrypt the
ciphertexts.

| Step | Data Flow | Protection |
|----|----|----|
| Encryption | Server k encrypts Z_k under CPK | Only collective decryption possible |
| Cross-product | Server A computes z_A x Enc(z_B) | Server A sees Enc(z_B), not z_B |
| Partial decrypt | Each server provides share_k | Individual shares are meaningless |
| Fuse | Client combines all K shares | Only aggregate correlation revealed |

### GLM with Encrypted Labels

The response variable y never leaves the label server in plaintext.
Non-label servers compute gradient contributions using encrypted y.

| Step | Data Flow | Protection |
|----|----|----|
| Encrypt y | Label server encrypts y under CPK | ct_y is opaque to all others |
| Gradient | Non-label server computes X_k^T (ct_y - mu) | Result is encrypted |
| Threshold decrypt | All K servers provide partial shares | Only p_k-length gradient revealed |
| BCD update | Server k solves local system | Raw features X_k stay local |

------------------------------------------------------------------------

## What the Client Sees

The client (researcher coordinating the analysis) receives:

1.  **Collective Public Key (CPK)** – cannot decrypt anything alone
2.  **Encrypted ciphertexts** – opaque blobs, undecryptable without all
    K secret key shares
3.  **Partial decryption shares** – each share individually reveals
    nothing
4.  **Final aggregated results** – correlation coefficients, GLM
    coefficients, deviance

The client CANNOT see:

- Individual-level data from any server
- Raw feature values
- Response variable values
- Residuals or per-observation predictions

------------------------------------------------------------------------

## What Each Server Sees

| Server Role | Sees | Does NOT See |
|----|----|----|
| Any server | Its own raw data | Other servers’ raw data |
| Any server | Encrypted columns from others (opaque) | Plaintext of encrypted columns |
| Any server | mu and w vectors (derived from aggregate eta) | Which specific observations drive the aggregate |
| Label server | Response variable y in plaintext | Other servers’ feature matrices |
| Non-label server | Encrypted y (ct_y) | Response variable in plaintext |

------------------------------------------------------------------------

## Threshold Property

Decryption requires cooperation from ALL K parties. The collective
secret key is additively shared:

``` math
sk = sk_1 + sk_2 + \ldots + sk_K
```

Each party k computes a partial decryption share using its sk_k. Only
when ALL K shares are combined can the ciphertext be decrypted. This is
a **K-of-K threshold scheme**: every party must participate.

------------------------------------------------------------------------

## Collusion Resistance

Even if K-1 servers collude (share their secret key shares), they cannot
decrypt any ciphertext without the remaining server’s share. The RLWE
assumption guarantees that the partial information from K-1 shares is
computationally indistinguishable from random noise.

Example with 3 servers:

- Servers 1 and 2 collude: they hold sk_1 + sk_2, but need sk_3
- Servers 1 and 3 collude: they hold sk_1 + sk_3, but need sk_2
- Only if all 3 cooperate: sk_1 + sk_2 + sk_3 = sk, and decryption
  succeeds

------------------------------------------------------------------------

## What IS Revealed

The following aggregate statistics ARE revealed to the client:

| Analysis    | Revealed                                            |
|-------------|-----------------------------------------------------|
| Correlation | Pairwise correlation coefficients r_ij              |
| PCA         | Eigenvalues and loadings (derived from correlation) |
| GLM         | Regression coefficients, deviance, fit statistics   |

These are summary statistics over n observations. They do not directly
reveal individual data points. However, as with any statistical release,
care should be taken with very small samples or highly identifying
variable combinations.

------------------------------------------------------------------------

## PSI Alignment Security

Record alignment uses an ECDH-based Private Set Intersection (PSI)
protocol. Each server holds a secret P-256 scalar that never leaves the
server. The protocol has the following security properties:

### Dictionary Attack Resistance

Unlike SHA-256 hashing, where an attacker who obtains the hash list can
pre-compute hashes for a dictionary of plausible identifiers and match
them, ECDH-PSI masks each identifier with a server-held secret scalar.
Without knowledge of the scalar, the masked curve points are
computationally indistinguishable from random group elements. An
attacker who intercepts the masked points cannot mount a dictionary
attack because they cannot reproduce the server’s scalar multiplication.

### Server Scalar Confidentiality

Each server $`k`$ draws a secret scalar $`\alpha_k`$ uniformly at random
from the P-256 curve order. This scalar is generated and stored
exclusively on the server. It is never transmitted to the client or to
other servers. The masked points $`\alpha_k \cdot H(\text{id})`$ are
one-way under the elliptic-curve discrete logarithm assumption: given a
masked point, recovering $`\alpha_k`$ is computationally infeasible.

### Unlinkability Across Servers (DDH)

The client receives doubly-masked points from each server pair. Under
the Decisional Diffie–Hellman (DDH) assumption on P-256, the client
cannot link a singly-masked point from server $`A`$ to the corresponding
singly-masked point from server $`B`$. The client can only determine
whether two doubly-masked points (one from each direction) correspond to
the same identifier, which is exactly the information needed for
alignment and nothing more.

### Non-Matching Identifier Privacy

No party learns which identifiers are held by other servers but are
absent from its own dataset. The intersection protocol reveals only the
set of common identifiers (in masked form). Identifiers present on
server $`A`$ but absent from server $`B`$ produce masked points that
have no match in the doubly-masked output, and the client simply
discards them without learning anything about the underlying identifier
values.

------------------------------------------------------------------------

## Assumptions and Limitations

1.  **Semi-honest model**: Parties are assumed to follow the protocol
    correctly. A malicious party could send corrupted key shares or
    ciphertexts. Verifiable computation is not implemented.

2.  **RLWE hardness**: Security relies on the Ring Learning With Errors
    problem being computationally hard. RLWE is a leading candidate for
    post-quantum cryptography.

3.  **Network security**: The protocol assumes secure (TLS)
    communication channels between the client and servers.
    DataSHIELD/Opal provides this via HTTPS.

4.  **Aggregate leakage**: The revealed aggregate statistics
    (correlations, coefficients) could theoretically allow inference
    about individuals if the sample is very small or variables are
    quasi-identifiers. This is inherent to any statistical release, not
    specific to MHE.

5.  **No differential privacy**: The protocol provides cryptographic
    protection of intermediate computations but does not add noise to
    final outputs. Combining MHE with differential privacy is possible
    but not currently implemented.
