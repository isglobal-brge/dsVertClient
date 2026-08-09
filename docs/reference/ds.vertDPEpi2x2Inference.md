# Epidemiological 2-by-2 inference with DP and sampling uncertainty

Builds conservative simultaneous confidence regions for two binomial
risks and their epidemiological contrasts from one validated DP 2-by-2
release. It first forms a simultaneous DP count box, then unions exact
Clopper–Pearson intervals over every integer table compatible with that
box. Bonferroni allocation across the DP event and three sampling
intervals gives the reported joint lower coverage bound. No noisy cell
is treated as an exact confidential count.

## Usage

``` r
ds.vertDPEpi2x2Inference(
  x,
  exposed = 2L,
  event = 2L,
  level = 0.95,
  mechanism_alpha_share = 0.5
)
```

## Arguments

- x:

  A released `ds.vertDPContingency` with a 2-by-2 table.

- exposed, event:

  Exposed row and event column by name or index.

- level:

  Desired lower bound for joint sampling-plus-mechanism coverage.

- mechanism_alpha_share:

  Fraction of total non-coverage allocated to the DP mechanism event.
  The remainder is divided equally across exposed, unexposed, and
  population exact-binomial intervals.

## Value

A `ds.vertDPEpi2x2Inference` object. It makes no server call and
consumes no additional privacy.

## Details

The sampling statement assumes independent, identically sampled privacy
units from a joint exposure/outcome population. Group-risk intervals are
conditionally binomial given the confidential group sizes; the
population risk interval is marginally binomial. Effects are
deterministic post-processing of their simultaneous base-risk regions.
The construction is intentionally conservative and provides no p-value.

## References

Clopper, C. J. and Pearson, E. S. (1934). The use of confidence or
fiducial limits illustrated in the case of the binomial. Biometrika
26(4), 404–413.
[doi:10.1093/biomet/26.4.404](https://doi.org/10.1093/biomet/26.4.404) .
