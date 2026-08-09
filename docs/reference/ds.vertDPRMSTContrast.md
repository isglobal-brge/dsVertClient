# Contrast restricted mean survival time from two DP releases

Computes comparison-minus-reference RMST and comparison/reference RMST
from two compatible, already released DP survival objects. It reuses
[`ds.vertDPRMST()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPRMST.md)
and conservative interval arithmetic, makes no server call, draws no
noise and spends no additional privacy. The joint mechanism confidence
is retained only for the exact same signed artifact; otherwise it uses a
Bonferroni lower bound. No sampling inference is implied.

## Usage

``` r
ds.vertDPRMSTContrast(
  comparison,
  reference,
  tau = NULL,
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

- tau:

  One or more common finite restriction times within the signed public
  time interval. `NULL` uses the common upper bound.

- comparison_label:

  Non-empty label for the comparison group.

- reference_label:

  Distinct non-empty label for the reference group.

## Value

A `ds.vertDPRMSTContrast` data frame with RMST differences, ratios,
conservative joint mechanism limits and source-release provenance.
