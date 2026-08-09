# Diagnostic-accuracy measures with simultaneous DP-noise regions

Purely post-processes one validated DP 2-by-2 histogram. Rows must
encode disease status and columns must encode the diagnostic-test
result. The positive disease row and positive test column are mandatory
explicit choices. No continuity correction, server call, p-value, or
classical sampling confidence interval is used. Measures include
sensitivity, specificity, predictive values, prevalence, accuracy,
balanced accuracy, F1 score, likelihood ratios, and the diagnostic odds
ratio.

## Usage

``` r
ds.vertDPDiagnostic2x2(x, disease_positive, test_positive, level = 0.95)
```

## Arguments

- x:

  A released `ds.vertDPContingency` with a 2-by-2 table.

- disease_positive:

  Positive disease-status row by name or index.

- test_positive:

  Positive diagnostic-test column by name or index.

- level:

  Simultaneous coverage for DP mechanism noise.

## Value

A `ds.vertDPDiagnostic2x2` object containing typed estimates and
simultaneous DP-mechanism regions. No server call is made.
