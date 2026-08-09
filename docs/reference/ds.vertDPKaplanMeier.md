# Kaplan–Meier Curve from One DP Survival Release

Pure post-processing; this function makes no server request and consumes
no additional privacy allocation.

## Usage

``` r
ds.vertDPKaplanMeier(x)
```

## Arguments

- x:

  A released `ds.vertDPSurvival` object.

## Value

A data frame with public grid time and the DP-histogram-derived survival
estimate. Its `uncertainty_scope` attribute states that only DP
mechanism noise is covered.
