# Formal binomial/Poisson DP result boundary

Status: internal, reviewable client boundary; no export, DataSHIELD call,
`AggregateMethod`, runtime command or production claim.

## Implemented result

`R/formal_glm_dp_result.R` consumes a canonical public coefficient release
only after it verifies all of the following:

- the local Phase-0 schema commitment and the complete trusted `K >= 2`
  pinset agree;
- the two designated computation peers signed the byte-identical release;
- their signed attestation binds the all-K manifest admission, Phase-1.8 v2 to
  Phase-1.9 v2 protected materializer, hidden execution-validity consumption,
  neutralisation before noise, zero intermediate openings and one final DP
  vector opening;
- bounds, contribution caps and sensitivity were derived server-side;
- the productive backend is the explicitly selected two-independent-full-draw
  discrete-Gaussian backend, with complete `(epsilon, delta)` at each peer,
  the signed variance/radius certificate, the finite-support transfer term and
  no automatic backend substitution;
- exact replay/post-processing of the already admitted outer capsule has no
  request limit or history-dependent denial, while any future promoted new
  capsule remains subject to the outer finite lifetime gate; and
- every released coordinate uses the bridge's non-negative shifted lattice,
  its even shifted upper bound decodes to an outward-rounded public box by
  less than one lattice step, and the signed coefficient reconstructed after
  subtracting the public midpoint is inside that box; the backend certificate
  names Ring128 and must report exactly 128 bits; the no-wrap numeric
  certificate's total error contains all listed components; and
- each Ed25519 signature has its fixed canonical encoded length before it is
  decoded, and the public release digest covers the complete signed artifact.

The returned object is deliberately not a `glm` or `ds.glm`. It exposes the DP
coefficient vector, odds ratios or incidence-rate ratios, a simultaneous
**mechanism-only** region, an upper bound on the conditional DP-noise
covariance, numeric/privacy certificates and release identity. Sticky
retry/restart/rechunk yields the same signed bytes. Losing or rotating a noise
root creates a visible new composed release identity and never creates a
history-based block.

Predictions accept caller-public covariates only. Numeric covariates follow the
signed clipping policy, categorical values must be in the signed domain, and
registered public offsets follow their signed bounds (a model with no offset
rejects one instead of interpreting it ambiguously). Prediction reports whether
any covariate or offset was clipped. Both the point and mechanism bounds remain
inside the certified PWL domain; out-of-domain or non-finite output fails
explicitly. The coefficient mechanism region is propagated through either the
committed PWL link or the canonical logit/log link. Prediction performs no DSI
call and cannot request protected fitted values, residuals or influence
diagnostics.

## Why classical inferential fields are not fabricated

The current formal mechanism releases only a noisy version of the deterministic
`beta_T` from a capacity-normalised, box-constrained, fully ridge-regularised,
fixed-iteration PWL surrogate. A coefficient draw alone does not contain the
active count, score, bread/information, score covariance, ordinary likelihood,
saturated likelihood or protected residuals.

Consequently:

- `mechanism_std_error_upper` and `mechanism_covariance_upper` describe DP
  perturbation conditional on the protected snapshot; they are not sampling
  standard errors or a population covariance;
- `log_likelihood`, `deviance` and `AIC` are unavailable rather than inferred
  from the released coefficients;
- Wald and likelihood-ratio p-values are unavailable; full ridge, a public
  box, a finite optimizer, a surrogate loss and DP noise invalidate the
  ordinary unpenalised chi-square shortcut; and
- separation, collinearity and the all-zero Poisson boundary remain explicit
  failures of the *ordinary comparator*, while the stated regularised formal
  target stays unique. No Firth, pseudoinverse or estimator-changing fallback
  is selected silently.

`.dsvert_formal_glm_statistic()` returns a typed
`requires_joint_dp_inference_artifact` state for these requests with zero
additional server calls. This is a scientific validity boundary, not a
privacy-budget or request-count denial.

## Contract for the separately reviewed inference artifact

`.dsvert_formal_glm_joint_inference_requirements()` is the machine-readable
handoff. It accepts only the signed compilation and derives the fixed vector
shape; there are no analyst sensitivity, bound or contribution-cap arguments.
For `p` coefficients, the proposed single stacked release contains:

1. the existing `p` DP coefficient coordinates;
2. `p(p+1)/2` capacity-normalised penalised-surrogate bread coordinates;
3. `p(p+1)/2` capacity-normalised bounded patient-score outer-product
   coordinates;
4. one integrated-PWL fitted surrogate loss; and
5. one pre-registered reference-model integrated-PWL surrogate loss.

The labels are normative. The loss is not ordinary `logLik`; the corresponding
penalised fixed-iteration contrast is not automatically a classical deviance,
LR statistic or AIC. A future method may report a penalised-surrogate sandwich
research estimate only after mechanism-aware coverage validation. Classical
Wald/LR names remain unavailable unless a separately stated estimand and null
calibration are proven.

The server implementation must derive one stacked L2 sensitivity from the
cross-signed schema and patient-level materializer. It must include the
data-dependent movement of `beta_T`, all coordinates, output quantisation and
numeric error. For a capacity-normalised Lipschitz statistic `S`, the generic
starting bound is

`Delta S <= gamma * row_range / C + L_beta * Delta beta`.

This shortcut cannot be applied blindly to the current PWL segment slope,
whose derivative is discontinuous at knots. The information block therefore
needs either a separately certified continuous curvature surrogate or a
server-derived global range bound. The implementation may not advertise a
tighter Lipschitz sensitivity that the committed arithmetic graph does not
satisfy.

The new artifact must reuse Phase-1.8 v2 and Phase-1.9 v2, consume hidden
validity before noise, mask invalid execution to the neutral tuple inside the
protected finalizer, perform no intermediate opening, and append the durable
sticky ledger before its one public vector opening. It supports K=2 and K>2
with exactly two pinned computation peers and at least one honest,
non-colluding peer. Malicious custodians and collusion of both computation
peers remain out of scope.

## Executable fixtures

`tests/testthat/test-formal-glm-dp-result.R` covers both families and K=2, K=3
and K=5. It checks central exact-rational oracle decoding, canonical-link public
prediction, separated binomial and all-zero Poisson comparators, signature and
pin tampering, cross-schema replay, hidden-validity/path substitution, quota
field substitution, sticky retry/restart/rechunk, visible root rotation, and
the typed non-promotion of covariance, log-likelihood, deviance, AIC, Wald, LR
and protected residuals.
