# Nelson–Aalen Curve from One DP Survival Release

Pure post-processing; this function makes no server request and consumes
no additional privacy allocation.

## Usage

``` r
ds.vertDPNelsonAalen(x)
```

## Arguments

- x:

  A released `ds.vertDPSurvival` object.

## Value

A data frame with public grid time and the DP-histogram-derived
cumulative hazard. No sampling confidence interval is claimed.
