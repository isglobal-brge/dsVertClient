# Hash Identifier Column (Deprecated)

**Deprecated**: Use
[`ds.psiAlign`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md)
instead.

Client-side function that retrieves hashed identifiers from a DataSHIELD
server. Used as the first step in record matching for vertically
partitioned data.

## Usage

``` r
ds.hashId(data_name, id_col, algo = "sha256", datasource = NULL)
```

## Arguments

- data_name:

  Character string. Name of the data frame on the server.

- id_col:

  Character string. Name of the identifier column to hash.

- algo:

  Character string. Hash algorithm. Default is "sha256".

- datasource:

  A DataSHIELD connection object (DSConnection). If NULL, uses the first
  available connection.

## Value

A list containing:

- `hashes`: Character vector of hashed identifiers

- `n`: Number of observations

## Details

**Deprecated**: This function uses SHA-256 hashing, which is vulnerable
to dictionary attacks when identifiers are predictable (e.g. sequential
patient IDs). Use
[`ds.psiAlign`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md)
instead, which provides ECDH-PSI alignment with stronger privacy
guarantees: the client only sees opaque elliptic curve points that
cannot be reversed to identifiers.

This function is typically called on the "reference" server in a
vertical partitioning scenario. The returned hashes are then passed to
[`ds.alignRecords`](https://isglobal-brge.github.io/dsVertClient/reference/ds.alignRecords.md)
on other servers to align the data.

## See also

[`ds.psiAlign`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md)
for the recommended replacement,
[`ds.alignRecords`](https://isglobal-brge.github.io/dsVertClient/reference/ds.alignRecords.md)
for legacy hash-based alignment

## Examples

``` r
if (FALSE) { # \dontrun{
# Get hashes from reference server
ref_hashes <- ds.hashId("D", "patient_id", datasource = conns$server1)

# Use hashes to align other servers
ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                datasources = conns[c("server2", "server3")])
} # }
```
