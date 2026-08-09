# Indirect standardization from aggregate observed and expected events

Computes a standardized mortality ratio (SMR) or standardized incidence
ratio (SIR) from caller-supplied aggregate observed and expected events.
The two labels share the same estimator; `measure` records the intended
scientific interpretation. The confidence interval is the exact central
Garwood interval for a Poisson count with fixed expected events.

## Usage

``` r
ds.vertIndirectStandardization(
  observed,
  expected,
  measure = c("SMR", "SIR"),
  scale = 100,
  level = 0.95
)
```

## Arguments

- observed:

  Numeric vector of non-negative, whole-number aggregate observed event
  counts by stratum.

- expected:

  Numeric vector of non-negative expected events by stratum.

- measure:

  Either `"SMR"` or `"SIR"`.

- scale:

  Positive multiplier for the reported ratio. The conventional value
  `100` reports percent of expected events; use `1` for an unscaled
  ratio.

- level:

  Confidence level.

## Value

An object of class `ds.vertIndirectStandardization` containing the
observed-to-expected ratio, its scaled estimate, an exact Poisson
confidence interval, provenance, and assumptions.

## Details

This function is entirely client-side. It makes no DSI or DataSHIELD
calls and exposes no new server-side values.

## References

Garwood, F. (1936). Fiducial limits for the Poisson distribution.
*Biometrika*, 28, 437–442.
