# Send one producer-ticketed opaque payload to its pinned recipient (internal)

The call intentionally has no key, filename, recipient or purpose
override. Ambiguous DSI responses replay the exact same request until a
monotonic availability deadline. The server commits by signed transfer
ID, absolute offset and payload hash; no request-count or rate quota is
consulted.

## Usage

``` r
.dsvert_store_typed_blob(
  blob,
  transfer,
  conn,
  session_id,
  producer_conn,
  .aggregate = DSI::datashield.aggregate
)
```
