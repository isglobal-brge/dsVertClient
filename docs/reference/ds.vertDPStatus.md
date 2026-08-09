# Differential-privacy policy status

Queries the reusable capsule control plane on every connected server. It
verifies the fixed per-capsule epsilon/delta, the authenticated lifetime
privacy bound, adjacency, exact logical-name-to-Ed25519 pin map, privacy
epoch, and the two designated noise peers. It also requires one
consensus `jdpc1_` privacy-accountant namespace protected by an
identity-bound immutable receipt. Exact replay is unlimited and there is
no request counter or accuracy decay, while a new distinct capsule
consumes one authenticated, non-refundable distinct-capsule reservation
unit at allocator commit. This call releases no protected statistic and
consumes no privacy allocation. Printing the result distinguishes the
remaining reservation units from request quotas.

## Usage

``` r
ds.vertDPStatus(datasources = NULL)
```

## Arguments

- datasources:

  DataSHIELD connections. If `NULL`, use the active set.

## Value

A named list of server statuses with class `ds.vertDPStatus`.
