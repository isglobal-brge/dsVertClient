# Relay an immutable biomedical capsule source to its pinned peers

Internal orchestration only. It performs a complete reusable-capsule
status handshake, relays paired authenticated Ring128 ciphertexts in
canonical owner/chunk order, and returns public completion metadata. It
never returns a coordinate, share, mask, seed, plaintext, or
patient-derived digest.

## Usage

``` r
.dsvert_dp_capsule_source_transport(
  manifest_json,
  datasources,
  status = NULL,
  allocation_openings = NULL,
  .aggregate = DSI::datashield.aggregate
)
```

## Arguments

- manifest_json:

  Canonical immutable biomedical capsule manifest.

- datasources:

  Complete named DataSHIELD connection set.

- status:

  Optional already validated full capsule-status handshake.

- allocation_openings:

  Canonical cross-signed allocation openings, named exactly by the two
  designated peers.

- .aggregate:

  Injectable DSI aggregate implementation for tests.

## Value

A redacted internal sampler-handoff receipt.
