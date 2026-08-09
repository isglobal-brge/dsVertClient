# DP-aware conditional test for a 2-by-2 contingency release

Test association by conditioning a plug-in latent table on deterministic
integer margins, then reproducing the signed capsule noise and clamp in
a Monte Carlo reference law. The input is either one already released
`ds.vertDPContingency` object (zero DSI calls) or data and variable
identifiers, in which case `ds.vertDPContingency` is called exactly
once.

## Usage

``` r
ds.vertFisher(
  data_name,
  var1 = NULL,
  var2 = NULL,
  server = NULL,
  alternative = c("two.sided", "greater", "less"),
  conf.int = TRUE,
  conf.level = 0.95,
  datasources = NULL,
  simulations = 9999L,
  mc_confidence = 0.95
)
```

## Arguments

- data_name:

  A released `ds.vertDPContingency` object, or the protected data-frame
  name.

- var1, var2:

  Row and column variables for a character `data_name`. They may be
  omitted for an existing release.

- server:

  Optional owner-peer assertion forwarded to `ds.vertDPContingency`. No
  separate column-discovery request is performed. It must be omitted for
  an existing release.

- alternative:

  One of `"two.sided"`, `"greater"`, or `"less"`. The directional
  alternatives refer to an odds ratio above or below one in the
  displayed row/column orientation.

- conf.int:

  Logical compatibility argument. A classical conditional odds-ratio
  interval is deliberately not returned from a DP release.

- conf.level:

  Compatibility confidence level recorded in the result.

- datasources:

  DataSHIELD connections. Omit for an existing release.

- simulations:

  Number of Monte Carlo replicates.

- mc_confidence:

  Confidence level of the exact binomial Monte Carlo interval.

## Value

An object of class `ds.vertFisher`. It retains the historical `p_value`,
`odds_ratio`, and `conf_int` fields, but the p-value is from the
explicitly labelled DP-aware plug-in conditional bootstrap and
`conf_int` is always `NULL`.

## Details

The projected DP table supplies a fitted contributing total and margins.
Those margins are rounded deterministically to a feasible positive
integer 2-by-2 configuration. Under the odds-ratio-one null, the
upper-left latent cell is drawn from its hypergeometric law. Each
simulated latent table then receives the same ideal one-draw or two-peer
discrete-Laplace law and public clamp as the signed release, after which
nuisance margins and the signed root-Pearson statistic are refitted.

This is not Fisher's finite-sample exact test for the confidential
table: noisy projected margins are neither the confidential margins nor
exact nuisance-sufficient statistics. The plug-in calibration is
asymptotic. Monte Carlo error and the signed finite-sampler
total-variation allowance are reported separately. Degenerate fitted
margins return a structured non-tested result. Only the signed Ring128
discrete-Laplace capsule is currently certified; other mechanisms fail
with a typed condition rather than being silently approximated.
