# ECDH-PSI Record Alignment (Blind Relay)

Privacy-preserving record alignment using Elliptic Curve Diffie-Hellman
Private Set Intersection (ECDH-PSI) with blind-relay transport
encryption. Aligns data frames across vertically partitioned DataSHIELD
servers so that rows correspond to the same individuals. The client
never sees raw EC points — only opaque encrypted blobs.

## Usage

``` r
ds.psiAlign(
  data_name,
  id_col,
  newobj = "D_aligned",
  ref_server = NULL,
  datasources = NULL
)
```

## Arguments

- data_name:

  Character string. Name of the data frame on each server.

- id_col:

  Character string. Name of the identifier column.

- newobj:

  Character string. Name for the aligned data frame on servers. Default
  is `"D_aligned"`.

- ref_server:

  Character string or NULL. Name of the reference server. If NULL
  (default), the first connection is used.

- datasources:

  DataSHIELD connection object or list of connections. If NULL, uses all
  available connections.

## Value

Invisibly returns a list with alignment statistics for each server:

- `n_matched`: Number of records matched

- `n_total`: Number of records on that server

## Details

This function performs privacy-preserving record alignment in a single
call, using ECDH-PSI with blind-relay transport encryption.

### Protocol overview

ECDH-PSI exploits the commutativity of elliptic curve scalar
multiplication: \\\alpha \cdot (\beta \cdot H(id)) = \beta \cdot (\alpha
\cdot H(id))\\.

All EC point exchanges are encrypted server-to-server (X25519 +
AES-256-GCM ECIES). The client acts as a blind relay, seeing only opaque
blobs.

1.  **Phase 0**: Each server generates an X25519 transport keypair.
    Public keys are exchanged via the client.

2.  **Phase 1**: The reference server masks IDs with scalar \\\alpha\\.
    Points are stored server-side (not returned to client).

3.  For each target server:

    - The reference encrypts masked points under the target's PK.

    - The target decrypts, generates scalar \\\beta\\, double-masks ref
      points (stores locally), masks own IDs, encrypts them under the
      ref's PK.

    - The reference decrypts, double-masks with \\\alpha\\, encrypts
      result under target's PK.

    - The target decrypts, matches double-masked sets, aligns data.

4.  A multi-server intersection ensures only records present on ALL
    servers are retained.

### Security (DDH assumption on P-256, malicious-client model)

- The client sees only opaque encrypted blobs — not EC points.

- Each server's scalar never leaves the server.

- PSI firewall: phase ordering + one-shot semantics prevent OPRF oracle
  attacks.

## See also

[`ds.vertCor`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md),
[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md)
for analysis functions that operate on aligned data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Align records across all servers using PSI
ds.psiAlign("D", "patient_id", "D_aligned", datasources = connections)

# Now "D_aligned" on all servers has matching, ordered observations
} # }
```
