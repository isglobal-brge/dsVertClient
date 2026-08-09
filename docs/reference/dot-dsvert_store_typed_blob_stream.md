# Pump one producer-spooled typed payload without materialising it (internal)

Source and destination use immutable absolute offsets. The client
verifies every frame immediately. The recipient verifies the complete
ticket-bound SHA-256 before atomic commit and signs a receipt binding
that hash and length; the producer verifies that receipt before
releasing its source spool. Thus the untrusted relay needs no second
O(payload) integrity spool and its peak memory and disk use are both
O(frame).

## Usage

``` r
.dsvert_store_typed_blob_stream(
  transfer,
  conn,
  session_id,
  producer_conn,
  .aggregate = DSI::datashield.aggregate
)
```
