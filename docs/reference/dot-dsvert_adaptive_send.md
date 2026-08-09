# Send chunked data via DataSHIELD with immutable geometry (internal)

The complete transfer geometry is fixed before the first request. A
callback exception or missing acknowledgement may be replayed when
`idempotent = TRUE`; every replay receives byte-for-byte identical chunk
data and identical indices until an availability deadline. There is no
request-attempt quota. Explicit negative or malformed replies are
protocol failures and are never retried.

## Usage

``` r
.dsvert_adaptive_send(data, send_one_chunk, target, idempotent = FALSE)
```

## Arguments

- data:

  Character. The full payload string to send.

- send_one_chunk:

  Function(chunk_str, chunk_index, n_chunks). Callback that sends one
  chunk and returns the unmodified result from
  [`DSI::datashield.aggregate`](https://datashield.github.io/DSI/reference/datashield.aggregate.html).

- target:

  Expected logical DataSHIELD connection name.

- idempotent:

  Whether exact replay is permitted after an ambiguous response. This
  must only be enabled for a server operation whose duplicate request is
  explicitly idempotent.

## Value

Integer. Number of chunks sent (invisible).
