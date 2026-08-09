# Competing-Risks Cumulative Incidence from One DP Release

Pure post-processing of the same histogram used for Kaplan–Meier and
Nelson–Aalen. It makes no server request and consumes no additional
privacy allocation.

## Usage

``` r
ds.vertDPCumulativeIncidence(x, cause = NULL)
```

## Arguments

- x:

  A released `ds.vertDPSurvival` object.

- cause:

  Optional public event-cause label. `NULL` returns every released
  cause.

## Value

A long data frame with public grid time, cause, and the Aalen–Johansen
cumulative-incidence estimate derived from the DP histogram.
