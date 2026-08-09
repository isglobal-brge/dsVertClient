# DP-aware inference for a signed cross-owner categorical table

This compatibility entry point obtains exactly one fixed-domain joint-DP
contingency release and performs only local post-processing. It never
discovers columns, constructs analyst-addressable one-hot objects, opens
an exact table, or invokes the retired cross-count endpoints. The
row/column ownership and domains must already be present in the
custodian-signed capsule `vertical_cross_specs` contract.

## Usage

``` r
ds.vertChisqCross(
  data,
  var1 = NULL,
  var2 = NULL,
  correct = TRUE,
  fisher = FALSE,
  datasources = NULL,
  verbose = TRUE,
  simulations = 9999L,
  mc_confidence = 0.95
)
```

## Arguments

- data:

  One of the two signed dataset names, or an existing
  `ds.vertDPContingency` release.

- var1, var2:

  Row and column variables. When `data` is an existing release these are
  optional orientation assertions.

- correct:

  Apply the DP-aware Yates-style correction for a 2-by-2 table.

- fisher:

  Also compute the DP-aware conditional 2-by-2 calibration from the same
  release. No second DSI request is made.

- datasources:

  DataSHIELD connections. Omit for an existing release.

- verbose:

  Print one non-sensitive progress message.

- simulations:

  Monte Carlo replicates for each requested calibration.

- mc_confidence:

  Confidence level for its Monte Carlo interval.

## Value

A `ds.vertChisq` result. If `fisher=TRUE`, the same object also contains
`fisher`, `fisher_p`, and `source_dp_release`.
