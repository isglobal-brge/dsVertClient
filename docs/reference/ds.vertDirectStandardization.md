# Direct standardization of aggregate stratum-specific rates

Computes a directly standardized rate using caller-supplied aggregate
event counts, person-time, and an external standard population. The
function is entirely client-side: it makes no DSI or DataSHIELD calls
and exposes no new server-side values.

## Usage

``` r
ds.vertDirectStandardization(
  cases,
  person_time,
  standard_population,
  scale = 1e+05,
  level = 0.95,
  ci_method = c("poisson_gamma", "normal")
)
```

## Arguments

- cases:

  Numeric vector of non-negative, whole-number aggregate event counts by
  stratum.

- person_time:

  Numeric vector of strictly positive aggregate person-time by stratum.

- standard_population:

  Numeric vector of non-negative external standard population weights.
  Values are normalized internally to sum to one.

- scale:

  Positive multiplier for reported rates, commonly `1e5`.

- level:

  Confidence level.

- ci_method:

  `"poisson_gamma"` (default) for the Fay–Feuer interval or `"normal"`
  for a bounded Wald interval.

## Value

An object of class `ds.vertDirectStandardization`. Its estimate,
standard error, confidence interval, and stratum rates are on the scale
selected by `scale`. The object records aggregate provenance,
assumptions, and any interval correction.

## Details

The default confidence interval is the Poisson-gamma interval of Fay and
Feuer, including its one-event upper-tail adjustment. A normal interval
using the independent-Poisson variance is also available; its lower
bound is truncated at zero and the truncation is reported in
`correction`.

## References

Fay, M. P. and Feuer, E. J. (1997). Confidence intervals for directly
standardized rates: a method based on the gamma distribution.
*Statistics in Medicine*, 16, 791–801.
