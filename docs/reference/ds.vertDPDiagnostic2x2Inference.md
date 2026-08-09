# Diagnostic-accuracy inference with DP and sampling uncertainty

Builds conservative simultaneous confidence regions from one validated
disease-by-test DP 2-by-2 release. It combines the signed simultaneous
DP count-box event with six Bonferroni exact Clopper–Pearson intervals
for sensitivity, specificity, PPV, NPV, prevalence, and accuracy.
Balanced accuracy, F1, likelihood ratios, and the diagnostic odds ratio
are deterministic interval post-processing of those base regions.

## Usage

``` r
ds.vertDPDiagnostic2x2Inference(
  x,
  disease_positive,
  test_positive,
  level = 0.95,
  mechanism_alpha_share = 0.5
)
```

## Arguments

- x:

  A released `ds.vertDPContingency` with a 2-by-2 table.

- disease_positive:

  Positive disease-status row by name or index.

- test_positive:

  Positive diagnostic-test column by name or index.

- level:

  Desired lower bound for joint sampling-plus-mechanism coverage.

- mechanism_alpha_share:

  Fraction of total non-coverage allocated to the DP mechanism event.
  The remainder is divided equally among the six exact-binomial base
  intervals.

## Value

A `ds.vertDPDiagnostic2x2Inference` object. It makes no server call and
consumes no additional privacy.

## Details

The sampling statement assumes independent, identically sampled privacy
units from one joint disease/test population. Conditional binomial
models are used for sensitivity, specificity, PPV and NPV; prevalence
and accuracy use marginal binomial models. No independence among the six
intervals is assumed. Possible zero-denominator states are reported
explicitly.

## References

Clopper, C. J. and Pearson, E. S. (1934). The use of confidence or
fiducial limits illustrated in the case of the binomial. Biometrika
26(4), 404–413.
[doi:10.1093/biomet/26.4.404](https://doi.org/10.1093/biomet/26.4.404) .
