# Fixed-grid survival quantiles from one DP survival release

Inverts the Kaplan–Meier curve and its simultaneous DP-mechanism band
from one validated `ds.vertDPSurvival` object. The finite-snapshot
quantile is the first public interval endpoint at which the
product-limit survival is no greater than `1 - probability`. If the
target is not reached by the public horizon, the point is `NA` and the
mechanism-region endpoint is represented by `Inf` with an explicit
status flag. This is zero-call post-processing and consumes no
additional privacy. The limits exclude sampling uncertainty and
public-grid discretisation error.

## Usage

``` r
ds.vertDPSurvivalQuantile(x, probabilities = 0.5)
```

## Arguments

- x:

  A validated released `ds.vertDPSurvival` object.

- probabilities:

  One or more finite event-distribution probabilities strictly between
  zero and one.

## Value

A `ds.vertDPSurvivalQuantile` data frame containing fixed-grid
quantiles, inverted simultaneous DP-mechanism limits, estimability
flags, and complete source-release provenance.
