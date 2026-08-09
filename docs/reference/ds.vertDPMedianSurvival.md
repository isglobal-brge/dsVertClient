# Median survival from one DP survival release

Exact convenience view of `ds.vertDPSurvivalQuantile(x, 0.5)`. It makes
no server call, draws no new noise, and retains the same fixed-grid and
uncertainty qualifications.

## Usage

``` r
ds.vertDPMedianSurvival(x)
```

## Arguments

- x:

  A validated released `ds.vertDPSurvival` object.

## Value

A one-row `ds.vertDPMedianSurvival` data frame retaining the full
survival-quantile provenance and mechanism-band contract.
