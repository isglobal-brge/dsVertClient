# Store one acknowledged opaque blob on one DataSHIELD target (internal)

This is the only legacy client helper that enables automatic chunk
replay. Both supported server methods reject conflicting duplicates and
return an exact logical `TRUE` after either the initial commit or an
identical replay.

## Usage

``` r
.dsvert_store_blob(
  blob,
  key,
  conn,
  session_id,
  store_function = "mpcStoreBlobDS",
  .aggregate = DSI::datashield.aggregate
)
```
