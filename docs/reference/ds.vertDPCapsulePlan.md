# Dry-run one server-authoritative DP capsule

Performs only the reusable-capsule status handshake and the three signed
manifest phases (draft, global-schema signature, and byte-identical
build). The selection is owned by custodian policy and signed workload
specifications: this function has no analyst-controlled expansion input.
It validates the current manifest, primitive-scope, coordinate-layout,
and mechanism contracts, but never resolves a protected snapshot,
materialises coordinates, invokes a producer, or creates a DP release.

## Usage

``` r
ds.vertDPCapsulePlan(datasources = NULL, status = NULL)
```

## Arguments

- datasources:

  Complete named DataSHIELD connection set. If `NULL`, use the active
  connections.

- status:

  Optional result from
  [`ds.vertDPStatus()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPStatus.md).
  It is validated again against `datasources`; if omitted, the handshake
  is performed.

## Value

A compact `ds.vertDPCapsulePlan` object containing the signed capsule
identity, immutable primitive selection, projected coordinate and
sensitivity cost, family inventory, mechanism calibration, pinset, and
explicit zero-access/zero-release guarantees for this dry-run.
