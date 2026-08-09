# One idempotent DSI fan-out cycle with an explicit availability state

This helper does not retry and never advances protocol state. It is for
offset-addressed, idempotent exchanges whose caller can replay the
identical request. DSI 1.8 represents a missing site result as a named
NULL; that and a top-level transport exception are the only non-fatal
unavailable states. A present but misassociated outer result is a
protocol failure.

## Usage

``` r
.dsvert_fanout_cycle(
  conns,
  expressions,
  operation = "DataSHIELD exchange",
  .aggregate = DSI::datashield.aggregate
)
```
