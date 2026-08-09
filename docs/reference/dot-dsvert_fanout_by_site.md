# One-request named per-site DSI fan-out (internal)

DSI 1.8 maps a named expression list to connections by logical site name
and submits async-capable connectors before fetching their results. This
helper validates the map before invoking that path and then applies the
strict result contract above.

## Usage

``` r
.dsvert_fanout_by_site(
  conns,
  expressions,
  operation = "DataSHIELD fan-out",
  allow_null = character(),
  result_contract = c("non_null", "logical_true"),
  .aggregate = DSI::datashield.aggregate
)
```
