# Get the current effective chunk size (internal)

Returns the cached effective chunk size if available (from a previous
successful probe), otherwise returns the initial chunk size from
`getOption("dsvert.chunk_size")` or the default (200KB).

## Usage

``` r
.dsvert_get_chunk_size()
```

## Value

Integer. Chunk size in characters.
