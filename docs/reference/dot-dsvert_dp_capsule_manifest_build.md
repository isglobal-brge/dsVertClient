# Build one server-authoritative biomedical capsule manifest

This internal orchestration has no analyst-controlled schema arguments.
It validates signed policy drafts from the complete pinset, derives the
only accepted global schema and logical snapshot, collects every schema
signature, and requires every peer to return the same signed manifest
bytes.

## Usage

``` r
.dsvert_dp_capsule_manifest_build(
  datasources,
  status = NULL,
  .aggregate = DSI::datashield.aggregate
)
```

## Arguments

- datasources:

  Complete named DataSHIELD connection set.

- status:

  Optional already validated reusable-capsule status handshake.

- .aggregate:

  Injectable DSI aggregate implementation for tests.

## Value

A list containing the canonical signed schema, canonical manifest,
capsule identity and validated connection context.
