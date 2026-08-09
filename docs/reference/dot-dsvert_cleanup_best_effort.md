# Final best-effort cleanup without exposing connector details (internal)

Cleanup is the sole phase allowed to ignore missing replies: it runs
only from an on-exit handler and cannot advance a statistical protocol.

## Usage

``` r
.dsvert_cleanup_best_effort(
  conns,
  expressions,
  .aggregate = DSI::datashield.aggregate
)
```
