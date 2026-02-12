# ECDH-PSI Record Alignment

Privacy-preserving record alignment using Elliptic Curve Diffie-Hellman
Private Set Intersection (ECDH-PSI). Aligns data frames across
vertically partitioned DataSHIELD servers so that rows correspond to the
same individuals, without exposing identifiers.

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
call, using ECDH-PSI instead of SHA-256 hashing for stronger privacy
guarantees.

### Protocol overview

ECDH-PSI exploits the commutativity of elliptic curve scalar
multiplication: \\\alpha \cdot (\beta \cdot H(id)) = \beta \cdot (\alpha
\cdot H(id))\\.

1.  The reference server hashes IDs to P-256 curve points and multiplies
    by a random scalar \\\alpha\\. Returns masked points to client.

2.  For each target server:

    - The target generates scalar \\\beta\\, double-masks ref points
      with \\\beta\\ (stores locally), and masks own IDs with \\\beta\\.

    - The reference double-masks target points with \\\alpha\\.

    - The target matches double-masked sets to find the intersection,
      then reorders its data to match the reference order.

3.  A multi-server intersection ensures only records present on ALL
    servers are retained.

### Security (DDH assumption on P-256)

- The client sees only opaque elliptic curve points — not reversible to
  identifiers, not vulnerable to dictionary attacks.

- Each server's scalar never leaves the server.

- The DDH assumption prevents the client from linking single-masked
  points across servers.

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
