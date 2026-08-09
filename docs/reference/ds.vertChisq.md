# DP-aware independence test on a 2-way contingency table

Calibrate a categorical independence statistic against the sampling
distribution and the privacy mechanism. The input is either an already
released `ds.vertDPContingency` object (zero DSI calls) or the original
data/variable identifiers, in which case the one immutable signed DP
capsule is obtained through `ds.vertDPContingency`. Ordinary chi-square
or Fisher reference laws are never applied to the noisy cells.

## Usage

``` r
ds.vertChisq(
  data_name,
  var1 = NULL,
  var2 = NULL,
  server = NULL,
  correct = TRUE,
  datasources = NULL,
  simulations = 9999L,
  mc_confidence = 0.95
)
```

## Arguments

- data_name:

  A released `ds.vertDPContingency` object, or the name of the protected
  data frame.

- var1, var2:

  Row and column variables when `data_name` is a character string. They
  may be omitted for an already released object.

- server:

  Optional owner-peer assertion forwarded to `ds.vertDPContingency`. No
  separate column-discovery request is performed. It must be omitted for
  an existing release.

- correct:

  Logical. If TRUE (default for 2x2 tables), apply Yates' continuity
  correction to the statistic. Its null law is still obtained by the
  same DP-aware bootstrap rather than a chi-square approximation.

- datasources:

  DataSHIELD connections.

- simulations:

  Number of parametric-bootstrap replicates. The default gives p-value
  resolution `1 / 10000`.

- mc_confidence:

  Confidence level for the reported exact binomial interval around the
  Monte Carlo p-value.

## Value

An object of class `ds.vertChisq` with elements `statistic`, `df`,
`p_value`, `observed`, `expected`, `residuals`, `n`, `correct`, and the
bootstrap/mechanism uncertainty contract.

## Details

The test follows the Monte Carlo independence-testing principle of
Gaboardi et al. (2016), adapted to dsVert's signed mechanism. It
estimates the latent contributing sample size and null margins from the
released table, simulates multinomial tables under the fitted
independence model, adds the signed release's one exact-GC draw or two
independent peer discrete-Laplace draws per cell, applies the same
common-lattice clamp, and repeats the nuisance fit for every replicate.
This is a parametric-bootstrap test: it is asymptotically calibrated
under positive cell probabilities, but is not a finite-sample exact
conditional test. Tables whose fitted expected count is below five
return a structured non-tested result rather than an unreliable p-value.

The production sampler is a finite binary-geometric approximation to its
signed ideal one-draw or two-peer-convolution law. Its signed
total-variation certificate is propagated into a separate calibration
interval. Bootstrap randomness is deterministically derived from public
release commitments solely for reproducible post-processing; it is not
DP mechanism randomness and adds no privacy cost.
