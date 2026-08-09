# Contrast two differentially private survival releases

Computes fixed-grid survival differences and ratios from two already
released, validated `ds.vertDPSurvival` objects. It performs no server
call and consumes no additional privacy. Difference and ratio limits use
conservative interval arithmetic on the inputs' simultaneous
DP-mechanism bands. If the inputs are not the exact same signed survival
artifact, joint coverage is bounded by Bonferroni without assuming
independent DP noise. These are mechanism-only regions, not sampling
confidence intervals or hypothesis tests.

## Usage

``` r
ds.vertDPSurvivalContrast(
  comparison,
  reference,
  comparison_label = "comparison",
  reference_label = "reference"
)
```

## Arguments

- comparison:

  A validated released `ds.vertDPSurvival` object for the
  numerator/comparison group.

- reference:

  A validated released `ds.vertDPSurvival` object for the
  denominator/reference group on the identical public time grid.

- comparison_label:

  Non-empty label for the comparison group.

- reference_label:

  Distinct non-empty label for the reference group.

## Value

A `ds.vertDPSurvivalContrast` data frame with point contrasts,
conservative joint DP-mechanism limits, denominator status and complete
source-release provenance.
