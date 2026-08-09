# DP stratified causal standardisation

Computes a saturated, stratum-standardised g-formula from one already
released DP table whose rows are the public Cartesian product of strata
and a binary treatment and whose columns are a binary outcome. Public
fixed target-population weights are used; no propensity score is
estimated.

## Usage

``` r
ds.vertDPCausalStandardization(
  x,
  strata,
  treatment,
  treated,
  standard_weights,
  event = 2L,
  level = 0.95
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

## Value

A `ds.vertDPCausalStandardization` object. It makes no server call and
consumes no additional privacy.

## Details

A causal interpretation additionally requires consistency, conditional
exchangeability within every declared stratum, positivity, no
interference, correct public row mapping, and scientifically valid
target weights. DP protects the release; it does not establish those
identification assumptions.
