# Align Records Across Servers

Client-side function that aligns records on DataSHIELD servers to match
a reference set of hashed identifiers. This ensures observations are
properly matched across vertically partitioned data.

## Usage

``` r
ds.alignRecords(
  data_name,
  id_col,
  reference_hashes,
  newobj = "D_aligned",
  algo = "sha256",
  datasources = NULL
)
```

## Arguments

- data_name:

  Character string. Name of the data frame on each server.

- id_col:

  Character string. Name of the identifier column.

- reference_hashes:

  Character vector. Hashes from the reference server (obtained via
  [`ds.hashId`](https://isglobal-brge.github.io/dsVertClient/reference/ds.hashId.md)).

- newobj:

  Character string. Name for the aligned data frame on servers. Default
  is "D_aligned".

- algo:

  Character string. Hash algorithm (must match reference). Default is
  "sha256".

- datasources:

  DataSHIELD connection object or list of connections. If NULL, uses all
  available connections.

## Value

Invisibly returns a list with alignment statistics for each server:

- `n_matched`: Number of records matched

- `n_reference`: Number of reference hashes

## Details

This function performs record alignment for vertically partitioned data:

1.  Takes reference hashes (from `ds.hashId` on one server)

2.  For each target server:

    - Hashes local identifiers

    - Matches against reference hashes

    - Reorders data to match reference order

    - Creates new aligned data frame

After alignment, all servers will have:

- Same number of observations

- Observations in the same order (by identifier)

- Only observations present in all partitions

## See also

[`ds.hashId`](https://isglobal-brge.github.io/dsVertClient/reference/ds.hashId.md)
for obtaining reference hashes

## Examples

``` r
if (FALSE) { # \dontrun{
# Step 1: Get reference hashes from first server
ref <- ds.hashId("D", "patient_id", datasource = conns$server1)

# Step 2: Align all servers (including reference) to these hashes
ds.alignRecords("D", "patient_id", ref$hashes,
                newobj = "D_aligned", datasources = conns)

# Now "D_aligned" on all servers has matching, ordered observations
} # }
```
