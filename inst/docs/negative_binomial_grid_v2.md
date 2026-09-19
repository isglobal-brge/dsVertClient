# Same-owner negative-binomial grid v2

The accepted workload specification is `negative_binomial_grid_v2`; its
artifact is `bounded-negative-binomial-likelihood-grid-v2`. Both packages
must use these semantics. The v1 specification and artifact are **sealed and
defective**: they omitted `y * log(theta + mu)` from the NB2 negative log
likelihood. Existing v1 releases remain historical records and must not be
relabelled, re-signed, reused as v2, or recalibrated in place. Validators reject
v1 and mixed-version descriptors. A corrected signed workload has a distinct
semantic identity and must follow the ordinary authorized publication path.

## Corrected loss

For an admitted integer count `0 <= y <= M`, positive signed `theta`, and
`mu = exp(eta)`, let `s(z) = max(z, 0) + log1p(exp(-abs(z)))`. The v2 loss is

```
C(y, theta) = lgamma(theta) + lgamma(y + 1) - lgamma(y + theta)
L(y, eta, theta) = C(y, theta)
                + theta * s(eta - log(theta))
                + y * s(log(theta) - eta).
```

This equals `-dnbinom(y, size = theta, mu = exp(eta), log = TRUE)`.
The last term is the stable equivalent of
`y * log(theta + exp(eta)) - y * eta`; it avoids exponentiating `eta` and
subtracting large almost equal terms. At `theta = 2`, `y = 1`, `eta = 0`, the
loss is `1.21639532432449`, while v1 returned `0.117783035656384`.

## Per-row cap and sensitivity

Predictors are normalized to `[0, 1]`, with an intercept of one. For candidate
`j`, the public coefficient bounds imply `|eta| <= A_j = sum(abs(beta_j))`.
The omitted term is bounded above by
`M * log(theta + exp(A_j))`: the logarithm is increasing in `eta` and is
positive at `A_j >= 0`. Thus adding this amount to an old upper bound is a
conservative correction. The implemented bound instead recomputes the exact
maximum on the enclosing signed interval, which can be tighter:

```
L''(eta) = (theta + y) * theta * exp(eta) / (theta + exp(eta))^2 >= 0
B_j = max(0, max over y = 0,...,M and eta in {-A_j, A_j} of L).
```

Convexity places each fixed-count maximum at an endpoint. Enumerating every
admitted integer count therefore bounds every admissible row. This is an
interval bound; it need not be the tightest bound for a particular design.
Server and client independently recompute the same `B_j`, ordered by theta
then beta. It is invalid to retain v1 caps or sensitivities for the new loss.

At lattice scale `S = 2^numeric_grid_bits`, let `q_j = ceiling(S * B_j)`.
Every admitted row contributes an integer in `[0, q_j]`; with observation
capacity `N`, the signed statistic maximum is `N * q_j`. Add/remove
adjacency uses raw L1 sensitivity `sum(q_j)` and raw L2 sensitivity
`sqrt(sum(q_j^2))`. The existing conservative replace-one contract doubles
both. Natural sensitivities divide the raw values by `S`. This accounts for
the loss correction and lattice rounding before DP calibration.

For `M = 8`, beta rows `(0, 0), (0, 1)`, theta values `1, 2`, and `S = 256`,
the corrected lattice caps are `(1598, 2770, 1896, 3338)`. With capacity 20,
the maxima are `(31960, 55400, 37920, 66760)`, and add/remove raw L1
sensitivity is `9602`.

The coordinate loss, contribution-domain, sensitivity-basis, and estimation
scope identifiers are v2. Candidate ordering, source scaling, row admission,
and missingness policies keep their existing versions. The intercept-only
Frequency method-of-moments route is unaffected.
