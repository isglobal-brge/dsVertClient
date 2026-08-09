# Binned DP quantiles from one validated describe release

Computes fixed-grid quantiles by deterministic client-side
post-processing of one already released
[`ds.vertDPDescribe()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPDescribe.md)
object. No server call is made and the additional privacy cost is
`(0, 0)`. The result identifies a public histogram bin and its interval;
it is not an exact sample quantile and does not interpolate within bins.

## Usage

``` r
ds.vertDPQuantile(x, probs = c(0.25, 0.5, 0.75))
```

## Arguments

- x:

  An intact released `ds.vertDPDescribe` capsule object.

- probs:

  Finite public probabilities in `[0, 1]`. Duplicates are removed and
  the result is sorted by probability within each variable.

## Value

A `ds.vertDPQuantile` data frame with the selected public bin, interval,
projected DP histogram mass, and a simultaneous 95-percent
mechanism/grid region. Sampling uncertainty is excluded.
