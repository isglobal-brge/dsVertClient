# DP-aware categorical inference

## Scope

`ds.vertChisq()` is a client-only post-processing adapter for one validated
`ds.vertDPContingency` block in the immutable biomedical capsule. Passing the
release object performs zero DSI calls and spends no additional epsilon or
delta. Passing data and variable names obtains/replays that same capsule and
then executes the identical client path. The retired exact contingency endpoint
is not used.

The method tests independence of two fixed-domain categorical variables held
by the same owner or by two different custodians. A cross-owner pair must be
declared in the custodian-signed `vertical_cross_specs` contract with
`family = "categorical_pair"` and `version = "v2"`; the client never discovers
variable ownership or domains from protected data.

For a cross-owner pair, each source owner contributes a fixed-capacity,
public-domain one-hot vector as encrypted Ring128 shares to exactly two pinned
computation peers. Those peers execute one concatenated exact-GC
multiplication for all `K x L` cells, reduce public fixed-length segments
locally, and inject only additive table shares into the global joint-DP vector.
No one-hot row, validity mask, product, exact cell count, marginal, alignment
hash, or additive table share is returned through DSI. Missing, out-of-domain,
or inconsistent local categories contribute the all-zero vector under the
signed rule.

`ds.vertChisqCross()` obtains this contingency release once and applies the
same DP-aware chi-square calibration locally. With `fisher = TRUE`, its
DP-aware 2-by-2 conditional calibration reuses the identical release and makes
no second DSI request. The retired column-discovery, one-hot, guarded exact
table, and `stats::fisher.test()` routes are not reachable from this API.

## Released mechanism reproduced by the reference law

For each latent cell count `X[i,j]`, the productive release is

```
W[i,j] = clamp(X[i,j] + (L1[i,j] + L2[i,j]) / S, 0, U)
```

where:

- `S` is the signed common-lattice scale;
- `U` is the signed public coordinate capacity;
- `L1` and `L2` are independent complete peer draws;
- each `L` is the difference of two iid geometric variables with continuation
  probability `exp(-epsilon / sensitivity_steps)` in the ideal reference; and
- the clamp is coordinate-wise and is applied once after signed Ring128 decode.

The implementation uses a finite binary-geometric sampler. The release carries
an authenticated total-variation certificate for its distance from the ideal
law. The test therefore simulates the ideal two-draw convolution and widens the
reported Monte Carlo interval by a conservative bound covering both the
observed and reference mechanism distributions. It never pretends that the
finite sampler is bit-for-bit reproduced by R's pseudorandom generator.

This is the same ideal-law-plus-signed-TV convention used by the capsule's
published accuracy radii. The clamp itself, including its fractional common
lattice, is reproduced exactly by the client post-processing.

## Statistic and nuisance fit

The true number of units contributing a valid pair is not released as a
separate coordinate. It is estimated as the nearest feasible integer to the
noisy table total, bounded by `U`. The released table is projected in Euclidean
distance onto the non-negative simplex with this total. Row and column
probabilities are estimated from the projected table, and the fitted null mean
is their outer product multiplied by the fitted total.

The reported statistic is the Pearson-style `Q_D` distance used by the Monte
Carlo independence procedure. The privacy-noise variance is represented by the
simulated reference distribution rather than inserted into an ordinary
chi-square tail approximation. Yates' subtraction remains available for
2-by-2 tables as a choice of statistic, but it is never paired with an ordinary
chi-square reference law.

For every bootstrap replicate the implementation:

1. draws a multinomial table under the fitted independence model;
2. draws two complete discrete-Laplace noise contributions per cell;
3. applies the public `[0, U]` clamp;
4. repeats the sample-size, simplex and margin fit; and
5. recomputes the same statistic.

The p-value is `(1 + exceedances) / (B + 1)`. Ties are included in the upper
tail. A failed bootstrap nuisance fit is assigned an infinite statistic, which
implements the fail-to-reject behavior conservatively. An exact binomial
interval reports Monte Carlo uncertainty; a second interval additionally
includes the signed finite-sampler TV allowance.

The simulation seed is a SHA-256 derivation from public capsule, vector,
orientation and statistic commitments. It is not user-selectable. It makes
repeated client post-processing byte-stable without serving as mechanism
randomness, and the caller's R random-number state is restored unchanged.

## Validity claim

The method is a DP-aware parametric bootstrap with plug-in nuisance parameters.
Under fixed positive row and column probabilities and increasing contributing
sample size, the noisy total and projected margins are consistent, while the
privacy noise remains order one. This supports asymptotic calibration. It does
not establish a finite-sample exact level for the composite independence null.

The implementation returns a structured non-tested result, not a p-value, if a
margin is degenerate or a fitted expected cell count is below five. This is the
same conservative expected-count policy used in the Monte Carlo independence
procedure of Gaboardi, Lim, Rogers and Vadhan.

Primary methodological references:

- Gaboardi M, Lim H, Rogers R, Vadhan S. *Differentially Private Chi-Squared
  Hypothesis Testing: Goodness of Fit and Independence Testing*. ICML 2016.
  <https://proceedings.mlr.press/v48/rogers16.html>
- Wang Y, Lee J, Kifer D. *Revisiting Differentially Private Hypothesis Tests
  for Categorical Data*. 2015. <https://arxiv.org/abs/1511.03376>

## Claims that are explicitly not made

- No ordinary Pearson, Yates or chi-square p-value is computed from the noisy
  table.
- The method is not Fisher's exact conditional test.
- The fitted contributing total is not claimed to be the confidential exact
  total.
- The Monte Carlo interval quantifies simulation error and the reported TV
  allowance; it does not eliminate plug-in nuisance uncertainty.
- Reusing the sticky release adds no privacy cost, but it does not undo the
  privacy/utility trade-off already fixed by the capsule.
- A Gaussian-selected categorical release is not silently interpreted with a
  Laplace reference law: chi-square and Fisher post-processing reject it with
  a typed mechanism-not-certified condition until a formal Gaussian
  calibration is available.
- Capsule creation completes the two-peer cross-signed allocation lifecycle
  before source tickets, source materialization, or cross-owner private-input
  binding. Source owners independently verify both pinned opening signatures;
  the relay is never an allocation authority. This allocator gate burns the
  lifetime unit but does not select among sibling release-instance candidates.
  The first valid vector START at each designated peer is the irrevocable
  per-capsule anti-sibling boundary before sampling. Neither boundary is a
  request quota: replay and all later client-side inference remain unlimited
  post-processing of the same sticky release.

## DP-aware conditional 2-by-2 route

`ds.vertFisher()` and `ds.vert.fisher()` now consume the same single signed
`ds.vertDPContingency` release. They never call the retired exact contingency
endpoint and never apply `stats::fisher.test()` to noisy cells. Passing an
existing release performs no DSI request; the character front door obtains the
capsule once and performs all inference client-side.

For a 2-by-2 release, the method projects the noisy table onto the fitted
simplex and deterministically rounds its row and column margins to a feasible
positive integer configuration. Under the odds-ratio-one null, the upper-left
latent cell is simulated from its hypergeometric law. The implementation then
adds the same ideal two-peer discrete-Laplace convolution, applies the public
clamp, refits the nuisance margins, and recomputes a tail-oriented signed
root-Pearson statistic. Greater and less alternatives use the displayed
row/column orientation; the two-sided alternative uses its absolute value.

The p-value is the plus-one Monte Carlo p-value. Its exact binomial simulation
interval and a second interval widened by the authenticated finite-sampler TV
allowance are both reported. The simulation seed is derived from public capsule
commitments and the alternative, is not analyst-selectable, and does not alter
the caller's R random state. Reusing the sticky release has zero additional
privacy cost.

This route is deliberately labelled a **conditional hypergeometric plug-in
DP-aware bootstrap**, not a Fisher-exact test. Projected noisy margins are not
the confidential margins and are not exact nuisance-sufficient statistics.
Consequently the claim is asymptotic under positive margins; finite-sample
exactness is explicitly false. The displayed odds ratio is the cross-product
ratio of the continuous projected table, not Fisher's conditional maximum-
likelihood estimate, and no classical conditional confidence interval is
returned. Degenerate conditional support produces a structured non-tested
result. Only the signed Ring128 discrete-Laplace mechanism is certified for
this reference law; a Gaussian or other mechanism raises a typed
`dsvert_dp_fisher_mechanism_not_certified` condition until its own reference
law and TV certificate are available.
