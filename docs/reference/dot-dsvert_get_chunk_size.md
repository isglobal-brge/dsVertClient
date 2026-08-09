# Get the current effective chunk size (internal)

Returns the cached effective chunk size if available (from a previous
successful transfer), otherwise returns the configured chunk size from
`getOption("dsvert.chunk_size")` or the portable default (640 KiB of
Base64 text, representing at most 480 KiB of raw bytes).

## Usage

``` r
.dsvert_get_chunk_size()
```

## Value

Integer. Chunk size in characters.
