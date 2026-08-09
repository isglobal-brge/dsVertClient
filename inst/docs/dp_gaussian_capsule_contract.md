# Bounded Gaussian capsule contract

`ds.vertDPGaussian()` is the formal same-owner Gaussian-regression route. It
retrieves one signed sticky joint-DP capsule and performs all model fitting as
deterministic client-side post-processing. It never calls an exact-statistic,
score, Hessian, row, residual or legacy GLM endpoint.

## Signed estimand

The custodian signs the dataset, outcome, predictor bounds and order, explicit
intercept choice, numeric grid, repeated-record collapse and complete-case
rule. Each finite value is clipped to its public bound and normalized to
`[0,1]`; repeated records are collapsed once per admitted patient before the
complete-case decision. A model is therefore a bounded finite-snapshot
estimand, not an unbounded ordinary least-squares fit to the original rows.

For `q` signed design terms the released block contains:

1. a noisy complete-case count `n`;
2. the noisy upper triangle of `X'X`, in column-major order;
3. noisy `X'y`, in signed design order; and
4. noisy `y'y`.

The block has `C = 1 + q(q+1)/2 + q + 1` coordinates. The server descriptor
binds the add/remove or replace-one sensitivity, coordinate order, maxima and
quantization grid. Cross-owner Gaussian products are
`reserved_not_materialized`; no legacy fallback is attempted.

When the design contains an intercept, `n` and `X'X[1,1]` are separate noisy
measurements. `n` is reported and used as the public upper clamp for moments;
the Gram coordinate governs the solve. They are not averaged. The returned
diagnostic reports their discrepancy and the signed reconciliation policy.

## Post-processing and numerical contract

The augmented second-moment matrix is symmetrized and negative eigenvalues are
clipped to zero before the solve. The result explicitly reports that this is
not an exact nearest projection subject to every moment constraint. Rank,
eigenvalue tolerance and condition diagnostics apply to the released,
projected normalized design—not to the unavailable exact design.

The default `ridge = 0` returns the unpenalized bounded-statistic estimand and
fails with `non_identifiable` when the released projected design is singular.
A positive `ridge` is an explicit predictor-only penalty on the normalized
design; the intercept is never penalized, and the result labels the changed
estimand. No pseudoinverse, ridge or identifiability repair is selected
silently.

Coefficients are returned both on normalized and original public-bound scales.
For outcome bounds `[L_y,U_y]` and predictor bounds `[L_j,U_j]`, slopes obey

`b_j(original) = (U_y-L_y) b_j(normalized) / (U_j-L_j)`.

With an intercept, its inverse transform also subtracts the predictor-lower-
bound offsets. Without a fitted normalized intercept, the reported original-
scale intercept entry is only the deterministic normalization offset and is
labelled as such.

## Accuracy and non-claims

The result includes a simultaneous 95% capsule-mechanism radius, a conservative
per-coordinate quantization enclosure, complete-case count bounds, mechanism
metadata and zero additional privacy cost after capsule retrieval. These are
finite-snapshot mechanism/quantization statements. The current method does not
claim nonlinear coefficient regions, classical standard errors, p-values,
population confidence intervals, individual fitted values or residuals.

`ds.vertGLM(..., family = "gaussian", dp_analysis_id = "...")` is an explicit
adapter to this contract only when the formula exactly matches the signed
same-owner artifact. Omitting `dp_analysis_id` preserves the pre-existing GLM
route and its separate maturity status. Binomial, Poisson, offsets, weights,
iterative controls and cross-owner models are not silently mapped to this
Gaussian capsule.
