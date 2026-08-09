# Joint DP and sampling inference for causal standardisation

Unions the signed DP count-box event with exact Clopper–Pearson
intervals for every positive-weight stratum-treatment risk, then
standardises those intervals with fixed public weights. It adds sampling
uncertainty to
[`ds.vertDPCausalStandardization()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCausalStandardization.md)
without another release.

## Usage

``` r
ds.vertDPCausalStandardizationInference(
  x,
  strata,
  treatment,
  treated,
  standard_weights,
  event = 2L,
  level = 0.95,
  mechanism_alpha_share = 0.5
)
```

## Arguments

- x:

  A released `ds.vertDPContingency` with two outcome columns.

- strata:

  Public stratum label for every table row.

- treatment:

  Public binary treatment label for every table row.

- treated:

  The treatment level interpreted as treated.

- standard_weights:

  Named, non-negative public weights for all strata.

- event:

  Event column by name or index.

- level:

  Simultaneous DP-mechanism coverage.

- mechanism_alpha_share:

  Fraction of total non-coverage assigned to the simultaneous DP
  mechanism event.

## Value

A `ds.vertDPCausalStandardizationInference` object.
