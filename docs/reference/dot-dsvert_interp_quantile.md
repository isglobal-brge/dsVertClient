# Interpolate a quantile from histogram bucket counts.

Pure helper function extracted for unit testing. Given bucket edges,
bucket counts (optionally with under/overflow), and target
probabilities, returns the linearly-interpolated quantile within the
target bucket.

## Usage

``` r
.dsvert_interp_quantile(edges, counts, probs, below = 0L, above = 0L)
```

## Arguments

- edges:

  Numeric vector of length K+1 defining the bucket edges.

- counts:

  Integer vector of length K with per-bucket counts.

- probs:

  Numeric vector of target probabilities (0, 1).

- below:

  Count of observations strictly below `edges[1]` (default 0).

- above:

  Count of observations strictly above `edges[K+1]` (default 0).

## Value

Numeric vector of the same length as probs.
