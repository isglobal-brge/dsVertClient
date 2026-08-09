# Differentially private indirect standardisation

Computes an observed-to-expected ratio (for example, an SMR or SIR) from
one validated DP strata-by-binary-outcome table and public
stratum-specific expected event probabilities. It makes no server call
and consumes no additional privacy budget. The returned region
propagates simultaneous DP mechanism noise; it is not a Garwood,
sampling, or population confidence interval.

## Usage

``` r
ds.vertDPIndirectStandardization(x, expected_rates, event = 2L, level = 0.95)
```

## Arguments

- x:

  A released `ds.vertDPContingency` whose rows are strata and whose two
  columns are event/non-event outcomes.

- expected_rates:

  Public expected event probabilities, one per stratum. Named values are
  reordered to the released row levels.

- event:

  Event column by name or index.

- level:

  Simultaneous coverage for DP mechanism noise.

## Value

A `ds.vertDPIndirectStandardization` object. No server call is made.
