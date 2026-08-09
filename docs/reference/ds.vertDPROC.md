# Diagnostic ROC curve and AUC from one ordered DP table

Purely post-processes one validated disease-status by ordered-test-bin
DP contingency table. It makes no server call and consumes no additional
privacy. The reported regions enclose DP mechanism noise under the
table's simultaneous certificate; they are not sampling confidence
intervals.

## Usage

``` r
ds.vertDPROC(
  x,
  disease_positive,
  score_order = NULL,
  direction = c("higher", "lower"),
  level = 0.95
)
```

## Arguments

- x:

  A released `ds.vertDPContingency` whose two rows encode disease status
  and whose columns are ordered diagnostic-score bins.

- disease_positive:

  Positive disease-status row by name or index.

- score_order:

  Optional permutation of every column, from low to high score. By
  default the released column order is used.

- direction:

  Whether higher or lower scores indicate disease.

- level:

  Simultaneous coverage for DP mechanism noise.

## Value

A `ds.vertDPROC` object containing the finite-snapshot ROC curve,
tie-adjusted AUC and conservative simultaneous mechanism regions.
