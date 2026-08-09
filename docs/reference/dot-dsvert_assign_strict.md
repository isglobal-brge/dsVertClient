# Fail-closed DataSHIELD assignment fan-out (internal)

DSI assignments do not return a typed acknowledgement. This wrapper uses
the connector error callback as the completion contract, suppresses
remote error text, and rejects any partial fan-out. A named expression
list is mapped to sites in one asynchronous DSI phase.

## Usage

``` r
.dsvert_assign_strict(
  conns,
  symbol,
  values,
  operation = "DataSHIELD assignment",
  .assign = DSI::datashield.assign.expr
)
```
