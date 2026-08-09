# Direct-standardisation inference with DP and sampling uncertainty

Builds a conservative joint confidence region for a directly
standardised risk from one validated DP strata-by-binary-outcome release
and fixed public standard-population weights. It combines the signed
simultaneous DP count box with Bonferroni exact Clopper–Pearson
intervals for every positive-weight stratum. The final region is the
public weighted sum of those stratum intervals.

## Usage

``` r
ds.vertDPDirectStandardizationInference(
  x,
  standard_weights,
  event = 2L,
  level = 0.95,
  mechanism_alpha_share = 0.5
)
```

## Arguments

- x:

  A released `ds.vertDPContingency` whose rows are strata and whose two
  columns are event/non-event outcomes.

- standard_weights:

  Fixed public non-negative weights, optionally named by the released
  stratum levels.

- event:

  Event column by name or index.

- level:

  Desired lower bound for joint sampling-plus-mechanism coverage.

- mechanism_alpha_share:

  Fraction of total non-coverage allocated to the DP mechanism event.
  The remainder is divided equally across all positive-weight stratum
  intervals.

## Value

A `ds.vertDPDirectStandardizationInference` object. It makes no server
call and consumes no additional privacy.

## Details

The sampling statement assumes independent privacy units and
conditionally binomial outcomes within each public stratum. It
conditions on the confidential stratum sample sizes and treats the
supplied standard weights as fixed. No independence among stratum
intervals is needed for the union bound. A possible zero observed sample
in a positive-weight stratum returns an uninformative `[0, 1]` component
rather than a fabricated estimate.

## References

Clopper, C. J. and Pearson, E. S. (1934). The use of confidence or
fiducial limits illustrated in the case of the binomial. Biometrika
26(4), 404–413.
[doi:10.1093/biomet/26.4.404](https://doi.org/10.1093/biomet/26.4.404) .
