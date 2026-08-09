# Formal GLM joint-inference artifact (internal review contract)

Status: internal specification and signed-release consumer only. It is not
exported, is not reachable from the public GLM front door, and has no DSI
producer. `production_ready` is therefore `FALSE`.

## Signed estimand

The artifact evaluates all protected contributions at the coefficient vector
from one already signed DP point release. That beta, its exact lattice
coordinates, release hash, privacy epoch, logical snapshot and materialization
root are inputs to the second signed release. A different beta, root or epoch is
a different adaptive release; marginal releases are never mixed.

For capacity `C`, coefficient vector `beta`, PWL mean `m`, signed left-knot
derivative `s`, patient design `x`, weight `w` and response `y`, the released
blocks are:

```text
H_total = sum_i w_i s(eta_i) x_i x_i' + C diag(ridge)
psi_i   = w_i x_i (y_i - m(eta_i))
J_HC0   = sum_i psi_i psi_i'
score   = sum_i psi_i
```

Optional scalar blocks are the bounded canonical-family log-likelihood at the
DP beta, the integrated PWL penalized surrogate loss and the active-patient
count. Their labels explicitly say that the canonical log-likelihood is not an
optimized PWL likelihood and that the count is DP noisy.

`H_total^-1` is exposed only as `penalized_curvature_inverse`. Because
`H_total` includes ridge, it is not called a frequentist model-based covariance
and does not produce sampling standard errors. A model-based covariance needs a
separate likelihood bread and model-variance meat. When `J_HC0` is certified,
the conditional research sampling estimator is
`H_total^-1 J_HC0 H_total^-1`; its HC0 standard errors remain explicitly
separate from DP mechanism uncertainty. Wald, LR, AIC and residual deviance are
not reconstructed. LR additionally requires one jointly calibrated nested
contrast artifact on the same snapshot; two marginal noisy log-likelihoods are
never subtracted.

## Server-owned bounds and joint sensitivity

The signed server plan is recomputed from the signed schema and fixed public DP
beta. The client has no sensitivity, bound or contribution-cap input.

For every selected coordinate the plan derives an exact rational upper bound
`b_j` on one admitted patient's absolute contribution. Public ridge terms have
zero sensitivity. For add/remove adjacency,

```text
Delta_2 <= Delta_1 <= sum_j b_j.
```

For replace-one adjacency the signed schema factor `gamma = 2` multiplies this
bound. Aggregate coordinate bounds use `C b_j`, add the deterministic ridge
curvature/loss terms, and are rounded outward to the output lattice. The
consumer recomputes the complete plan and rejects changes to the coordinate
order, beta, adjacency factor, bounds, sensitivity or lattice rounding.

The optional log-likelihood bounds use certified rational upper intervals for
`log`/`exp`; the surrogate loss uses the signed eta, outcome, mean and weight
bounds. If a block lacks such a derivation it must carry
`not_promoted_sensitivity_not_certified` and has no coordinate in the vector.

## DP and accuracy certificate

Each of the two designated pinned compute peers contributes one complete noise
draw. The per-peer Gaussian-core variance obeys the conservative signed rule

```text
sigma^2 >= Delta_2^2 * 2 log(1.25 / delta_core) / epsilon^2,
0 < epsilon <= 1.
```

The release delta is split exactly between the Gaussian core and finite-support
transfer. Independent full peer draws make the published variance upper bound
at least `2 sigma^2`; epsilon is not divided by peer count. The simultaneous
coordinate radius is checked by a distribution-independent union/Chebyshev
bound, `r^2 >= 20 d variance_upper`, for `d` released coordinates. This is
conservative but cannot advertise an implausibly narrow 95% region.

The inference release composes sequentially with the beta release by exact
rational addition of epsilon and delta. The durable ledger is accounting and
sticky memoization for this internal phase: retry, restart and rechunk preserve
the same signed coordinates and cannot be denied by request history. Any future
public capsule containing this release would still consume one unit at the
outer allocator commit; postpublication root rotation may not reroll it.

## Numerical and threat-model gates

The internal contract accepts only Ring128, a signed no-wrap certificate,
nonnegative shifted coordinate projection, hidden validity consumed before
noise, zero intermediate openings, one final vector opening, all-K manifest
admission, the complete K-peer pinset and two signatures over the identical
canonical release. Tests cover K=2, K=3 and K=5.

Covariance-like postprocessing additionally requires the complete simultaneous
DP region for beta to be strictly inside the coefficient box and
`lambda_min(H_hat) - (p r_H + deterministic_error) > 0`. Failure returns a
typed unavailable capability; there is no pseudoinverse, hidden ridge or
estimand change. A DP-noisy indefinite `J` is symmetrized and projected to PSD,
and the projection/eigenvalues remain visible.

This Ring128 artifact is not evidence that the generic Ring63/Ring127 planner,
exact comparison/truncation path or dynamic multiprecision fallback is
production-verified. No generic numerical-attestation flag may be promoted from
these fixture tests. Public E2E promotion still requires a real server-owned
materializer, exact protected execution, durable cross-signed release, DSI
source and public binomial/Poisson front door with no legacy fallback.
