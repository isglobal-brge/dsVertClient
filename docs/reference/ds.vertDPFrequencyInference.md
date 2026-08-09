# Conservative sampling regions for a DP frequency distribution

Combines one simultaneous DP count-box event with exact Clopper–Pearson
intervals for every public level. A Bonferroni union bound supplies the
advertised joint coverage. This is pure post-processing and makes no DSI
call.

## Usage

``` r
ds.vertDPFrequencyInference(x, level = 0.95, dp_fraction = 0.5)
```

## Arguments

- x:

  A validated `ds.vertDPFrequency` object.

- level:

  Requested joint coverage.

- dp_fraction:

  Fraction of the total error probability assigned to the simultaneous
  DP mechanism event; the remainder is split over categories.

## Value

A `ds.vertDPFrequencyInference` object.
