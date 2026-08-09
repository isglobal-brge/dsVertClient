# Audited maturity status of dsVert client methods

Returns the conservative method-maturity registry used to keep
experimental or compatibility routes out of the promoted path. The
status reflects both statistical validity and disclosure behavior under
a pinned, honest-but-curious peer model with an untrusted analyst/relay.

## Usage

``` r
ds.vertMethodStatus(method = NULL, status = NULL)
```

## Arguments

- method:

  Optional character vector of exported method names.

- status:

  Optional subset of `"promoted"`, `"provisional"`, `"compatibility"`,
  or `"quarantine"`.

## Value

A data frame with one row per public method, its canonical route,
maturity status, release contract, defensible scope, and principal
limitation. This registry is an auditable release contract, not a
cryptographic certification.
