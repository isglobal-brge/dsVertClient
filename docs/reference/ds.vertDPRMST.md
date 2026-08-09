# Restricted mean survival time from one DP survival release

Computes restricted mean survival time (RMST) by integrating the
left-continuous product-limit step curve on the release's fixed public
time grid. This is pure post-processing: it performs no server call and
consumes no additional privacy budget. Its limits propagate the
simultaneous DP mechanism band; they are not sampling confidence
intervals and do not cover error from replacing continuous event times
by the public grid. The result carries a copy of the source release's DP
and capsule provenance.

## Usage

``` r
ds.vertDPRMST(x, tau = NULL)
```

## Arguments

- x:

  A released `ds.vertDPSurvival` object.

- tau:

  One or more public finite restriction times greater than the release's
  lower time bound and no greater than its upper time bound. `NULL` uses
  the public upper bound.

## Value

A data frame containing RMST, conservative simultaneous DP-mechanism
limits, and source-release provenance for every requested `tau`.
