# Direct standardisation with simultaneous DP-noise bounds

Treats rows of an already released DP contingency table as strata and
one column as the event. Reference weights are public client-side
inputs.

## Usage

``` r
ds.vertDPDirectStandardization(x, standard_weights, event = 2L, level = 0.95)
```

## Arguments

- x:

  A released `ds.vertDPContingency` (strata by outcome).

- standard_weights:

  Non-negative public weights, optionally named by table row levels.

- event:

  Event column by name or index.

- level:

  Simultaneous coverage for all DP table cells.

## Value

A `ds.vertDPStandardization` object. No server call is made.
