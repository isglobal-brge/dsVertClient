# Send chunked data via DataSHIELD with adaptive chunk sizing (internal)

Wraps the standard chunk-and-send pattern with automatic fallback. On
the first large send in a session, probes the server with the initial
chunk size. If the probe fails (expression too large), reduces the chunk
size by 25\\ minimum is reached). The successful chunk size is cached
for all subsequent sends.

## Usage

``` r
.dsvert_adaptive_send(data, send_one_chunk, min_chunk_size = 10000L)
```

## Arguments

- data:

  Character. The full payload string to send.

- send_one_chunk:

  Function(chunk_str, chunk_index, n_chunks). Callback that sends a
  single chunk via
  [`DSI::datashield.aggregate`](https://datashield.github.io/DSI/reference/datashield.aggregate.html).
  Must throw an error on failure.

- min_chunk_size:

  Integer. Minimum chunk size before giving up. Default 10000 (10KB).

## Value

Integer. Number of chunks sent (invisible).

## Details

Small payloads (less than half the chunk size) are sent directly without
counting as a conclusive probe, since they don't test the actual limit.
