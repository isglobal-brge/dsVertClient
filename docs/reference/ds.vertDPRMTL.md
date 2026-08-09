# Restricted mean time lost from one DP survival release

Computes restricted mean time lost (RMTL) as the exact complement of
[`ds.vertDPRMST()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPRMST.md)
on the release's public interval:
`RMTL = (tau - time_lower_bound) - RMST`. The simultaneous mechanism
limits are transformed by the same identity with their order reversed.
No curve is re-estimated, no server is contacted, and no additional
privacy is spent. The inherited limits exclude sampling uncertainty and
public-grid discretisation error.

## Usage

``` r
ds.vertDPRMTL(x, tau = NULL)
```

## Arguments

- x:

  A validated released `ds.vertDPSurvival` object.

- tau:

  One or more public finite restriction times greater than the release's
  lower time bound and no greater than its upper time bound. `NULL` uses
  the public upper bound.

## Value

A `ds.vertDPRMTL` data frame retaining the complete RMST result and
provenance, plus restriction width, RMTL, and its conservative
simultaneous DP-mechanism limits.
