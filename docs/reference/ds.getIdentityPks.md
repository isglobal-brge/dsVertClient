# Query Server Identity Public Keys

Queries all connected DataSHIELD servers for their Ed25519 identity
public keys. Used by administrators to discover PKs for configuring the
`dsvert.trusted_peers` option on each server.

## Usage

``` r
ds.getIdentityPks(datasources = NULL)
```

## Arguments

- datasources:

  DataSHIELD connections. If NULL, uses all available.

## Value

Named list: server_name -\> identity_pk (base64url string).

## Examples

``` r
if (FALSE) { # \dontrun{
pks <- ds.getIdentityPks(datasources = conns)
# Set trusted peers on each server via opalr:
# dsadmin.set_option(opal, "dsvert.trusted_peers", paste(pks, collapse=","))
} # }
```
