# Verify a bounded Gaussian capsule certificate without DSI

Revalidates the signed schema, manifest build commitments, pinned-peer
vector receipts, Gaussian artifact membership and the intersecting
public DP chunks required by the current Merkle layout. A self-contained
certificate proves internal integrity. Strong peer authenticity
additionally requires a caller-supplied trusted pinset or the ephemeral
pinset cache populated by an online fit in this R session.

## Usage

``` r
ds.validateDPGaussianCertificate(x, trusted_pinset = NULL)
```

## Arguments

- x:

  A `ds.vertDPGaussian` result or its `provenance_certificate`.

- trusted_pinset:

  Optional trusted named peer-to-Ed25519-public-key map.

## Value

A verified public provenance view with separate `integrity_valid` and
`authenticity` fields.
