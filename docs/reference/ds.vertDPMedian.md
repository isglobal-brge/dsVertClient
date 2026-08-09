# Binned DP median from one validated describe release

A release-only convenience wrapper around
[`ds.vertDPQuantile()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPQuantile.md)
at probability `0.5`. It performs no DSI operation and has additional
privacy cost `(0, 0)`.

## Usage

``` r
ds.vertDPMedian(x)
```

## Arguments

- x:

  An intact released `ds.vertDPDescribe` capsule object.

## Value

A `ds.vertDPMedian` data frame with one binned median per released
variable and the same mechanism/grid metadata as
[`ds.vertDPQuantile()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPQuantile.md).
