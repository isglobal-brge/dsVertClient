# dsVert numeric-result surface inventory

Audit date: 2026-07-31

This inventory separates three claims that must not be conflated:

1. a disclosure/release contract;
2. statistical identifiability and estimand validity; and
3. numerical execution certification.

No public estimator currently has an operation-complete numerical execution
certificate, and the method registry currently authorises no route to report
one. The historical Ring GLM appears below only as an internal candidate
inventory for a future audited adapter. The public `ds.vertGLM` can instead
delegate an explicitly authorised compatible Gaussian request, select a signed
finite binomial/Poisson likelihood-grid release, or read a completed formal
certificate; none claims operation-complete Ring execution certification.
Fresh formal computation fails before DSI. The `ds.vertDP*` methods have a separate
integer-mechanism/accountant contract and must not be described as
Ring-certified.

## GLM operation blockers

| Stage | Current operation | Covered now | Promotion blocker |
|---|---|---:|---|
| Raw and standardized input | `glmStandardizeDS`, offset and weight validation | Partial | Finite/domain, raw bounds, encoded bounds and fixed-point representability are attested, without releasing observed extrema. |
| Fixed-point input sharing | `k2ShareInputDS` and related K>=3 share assembly | Partial | Canonical records are checked, but there is no operation-complete workload attestation. |
| Linear predictor | `k2ComputeEtaShareDS` -> `handleK2ComputeEtaFP[127]` | No | Public beta times secret X uses `ScalarVectorProductPartyZero/One[127]`, including asymmetric local truncate logic. It must be migrated to an exact mul-truncate adapter or formally verified with its complete preconditions. |
| Nonlinear inverse link | Ring127 sigmoid/exp Chebyshev chains, Beaver vecmul, affine/local public scaling | No | Exact truncation of selected vecmul calls does not certify public-scalar/local-scale steps, every polynomial intermediate, raw-product capacity, or domain enforcement. |
| Protected weights | Beaver multiplication of weights and residuals | No | Raw product and truncation need exact execution binding and bounds. |
| Gradient | `k2-full-iter-r3`, `k2GradientR1DS/R2DS` | No | The accumulated raw Beaver matvec is followed by local asymmetric truncation; Ring63 and Ring127 both remain unverified. |
| Hessian/SE | Repeated perturbed score evaluations and client finite differences | No | Perturbation, accumulated arithmetic error, rank decision and estimator/covariance error are not covered by the one-path planner bound. |
| Deviance | identity, Ring127 softplus/exp and aggregate sums | No | This is a separate nonlinear and accumulation workload not yet bound by the certificate. |
| Optimizer and back-transform | client L-BFGS, coefficient rescaling, covariance algebra | Partial | Finite result postconditions exist, but convergence and estimator error are not numerically certified. |
| K>=3 adapter | fusion/reordering plus the same two-party compute core | No | The complete fan-out/reordering/workload path needs its own E2E fixture and attestation. |

Consequently, the Ring4096 primitive is not yet a GLM fallback adapter. The
remaining work crosses the record format, input sharing, eta, nonlinear link,
gradient, deviance and inference stages; it cannot be implemented correctly by
changing `effective_backend` or replacing one multiplication call. The current
typed `client_backend_adapter_unavailable` refusal is the required behavior for
`exact_gc`/`multiprecision` GLM selections until all of those stages are bound
to the selected ring and pass a whole-workload plaintext oracle.

The data-free preflight planner nevertheless has a complete representability
contract for a future adapter. It derives the minimum fractional precision
from the custodian error budget, retains the Ring63/f20 and Ring127/f50 fast
floors when they suffice, and then chooses the smallest jointly advertised
exact ring through Ring4096. Magnitude products are planned in log2 space, so
an overflowing IEEE-double diagnostic cannot hide the actual public ring
requirement. If a genuinely available and E2E-attested exact backend exhausts
its ring or precision domain, the client raises
`numeric_backend_unrepresentable` with `required_ring_bits`,
`required_frac_bits`, the attempted precision, ring ceiling, and error budget.
An absent or incomplete exact/multiprecision implementation remains
`numeric_backend_unavailable` and is never labelled as a selected backend.
This planner guarantee does not promote the current GLM execution path.

Other public-scalar uses with the same blocker include
`k2Ring127LocalScaleDS`, Ring63/Ring127 spline coefficient multiplication, and
the local scale stages used by reciprocal, exp, log, sigmoid and softplus
handlers. The mere presence of an exact-GC primitive is not evidence that any
of these complete workloads use it.

## Public method inventory

| Canonical route | Numeric route | Current disposition | Required work before a numeric claim |
|---|---|---|---|
| `ds.vertNumericPreflight` | Data-free policy/bound planner | Preflight only | Never describe it as execution certification. Bind a real workload adapter and operation-complete runtime attestations. |
| `ds.vertCor` | Client post-processing of one signed same-owner or exact-GC cross-owner complete-case Gaussian sufficient-statistic block | Formal sticky joint-DP capsule; signed numeric/quantization certificate, outward mechanism regions, explicit PSD projection, and no pairwise fallback | Population-sampling inference is deliberately unavailable; zero-variance or missing-intercept cases return typed non-identifiability. |
| `ds.vertPCA` | Client eigendecomposition of the complete-case DP correlation release | Zero-cost post-processing with outward spectral/eigenvalue and eigengap regions; pairwise artifacts and individual scores are rejected | Directions without a certified eigengap remain explicitly unresolved. |
| `ds.vertChisq`, `ds.vertFisher` | Validated signed common-lattice DP table plus client-only mechanism-aware bootstrap | No client Ring arithmetic; exact signed decode is certified by the release | Validate the signed capsule, clamp/lattice metadata, fitted nuisance contract, Monte Carlo error and sampler-TV allowance. |
| `ds.vertChisqCross` | Capacity-padded Ring128 one-hot shares, one fixed exact-GC multiplication, private segmented reduction, sticky joint-DP table | Operation-complete no-wrap certificate; exact intermediates unopened | Validate the signed cross-owner descriptor, allocation openings, private alignment contract, common lattice, clamp metadata, mechanism-aware bootstrap and sampler-TV allowance. Gaussian categorical calibration remains unavailable rather than falling back to an uncertified law. |
| `ds.vertGLM` | Compatible Gaussian signed Synopsis; signed finite binomial/Poisson likelihood grid; completed two-authority formal certificate; retained Ring63/127 candidate | Read-only/signed-artifact adapter; no public Ring execution claim | Fresh formal computation is sealed. Close every GLM blocker above before proposing an audited public iterative Ring adapter. |
| `ds.vertCoxProfileNonDisclosive`, `ds.vertCox` | Completed two-authority formal Cox public certificate; retained Ring127 prototype remains sealed | Read-only sticky certificate adapter | Fresh Cox computation, covariance, standard errors, p-values, baseline hazards and sampling inference remain unavailable until a production-attested protected runtime exists. |
| `ds.vertCoxDiscreteNonDisclosive` | Completed binomial formal-GLM lattice certificate over one custodian fixed pooled-logistic grid | Read-only two-authority public-certificate adapter; no fresh Ring execution | Keep the discrete-time estimand explicit; no Cox-PH claim, source expansion, covariance or sampling inference. |
| `ds.vertNBFullRegTheta` | Validated sticky bounded non-negative integer Frequency release, or same-owner signed finite `(beta, theta)` likelihood grid | Intercept-only log-mean/NB2 method-of-moments post-processing, or additive bare-predictor finite-grid coefficient/theta selection | The grid has bounded integer outcome and public predictor bounds. No free likelihood optimization, covariance or sampling inference. |
| `ds.vertMultinom`, `ds.vertMultinomJoint`, `ds.vertMultinomJointNewton` | Validated sticky categorical Frequency release, or same-owner signed finite softmax likelihood grid | Intercept-only Jeffreys-smoothed log-odds, or additive bare-predictor finite-grid coefficient selection | The grid binds the categorical domain and public predictor bounds. Interactions/transforms, free optimization, covariance and inference remain unavailable. |
| `ds.vertOrdinal`, `ds.vertOrdinalJointNewton` | Validated sticky categorical Frequency release, or same-owner signed finite cumulative-logit likelihood grid | Intercept-only cumulative-logit thresholds, or additive bare-predictor finite-grid coefficient/ordered-threshold selection | The grid binds the complete clinical order and public predictor bounds. Interactions/transforms, free optimization, covariance and inference remain unavailable. |
| `ds.vertGEE` | Signed finite binomial/Poisson grid or completed formal GLM certificate, matching signed Gaussian random-intercept Synopsis, or a signed Gaussian AR(1) beta/rho working-GLS grid with strict within-cluster order | Read-only independence point adapter, Gaussian exchangeable model-based working GLS, or Gaussian AR(1) working GLS | The independent score equals the corresponding GLM score; the Gaussian random-intercept covariance has the same exchangeable working form. AR(1) has its own signed finite-grid kernel. Fresh formal computation, cluster statistics, robust covariance, standard errors and inference remain unavailable. |
| `ds.vertGLMM` | Signed binary same-owner `outcome ~ 1` random-intercept moments, same-owner finite binary marginal-likelihood grid, or same-owner finite bounded-count Poisson random-intercept grid | Population-average log-odds/observed-scale ICC moment approximation, additive binary bare-covariate finite random-intercept/one-to-three-random-slope selection, or additive Poisson random-intercept selection | PQL/Laplace/ML/REML, interactions, Poisson random slopes, more than three binary random slopes, covariance and sampling inference remain unavailable. |
| `ds.vertLMM`, `ds.vertLMM.k3` | Signed same-owner random-intercept Synopsis moments, fixed-effect ML/REML finite variance-ratio grid, or finite Gaussian random-slope candidate grid | Promoted signed-artifact post-processing, including K>=3 compatibility | No unrestricted optimisation/covariance, standard errors or sampling inference. |
| `ds.vertIPW` | Validated sticky binary treatment-by-outcome table, or complete signed categorical stratum-treatment-by-outcome table | Exact intercept-only ATE, or one-categorical-stratum saturated IPW/g-formula ATE with public weights and ATT/ATC with release-bound target-arm weights; no individual weights or propensity fit | ATT/ATC mechanism regions jointly optimize the signed simultaneous count box but intentionally omit sampling inference. Continuous or multiple covariates, trimming, fitted propensity/outcome models, standard errors and p-values remain unavailable. |
| `ds.vertMI` | Signed sticky categorical Synopsis with strict-missing marginals or pairs | Deterministic categorical MCAR marginal completion, two-response joint completion, or opt-in conditional categorical star model, with no analyst seed | It returns neither source rows nor Rubin/sampling inference. Only one categorical conditioning column is accepted by the star model; chained equations and MAR/MNAR models remain unavailable. |
| `ds.vertLASSO`, `ds.vertLASSO1Step`, `ds.vertLASSOCV` | Client transformations of a supplied fit | Compatibility/inherits input | They cannot inherit numeric certification as new estimators; the legacy names/estimands must remain explicit. |
| `ds.vertLASSOProximal` | Client Gaussian normal-equation reconstruction and proximal solver on one validated same-owner Synopsis artifact | Promoted post-processing, inherits input; no sampling inference or cross-owner design | Condition-number/PSD and convergence certificate for the reconstructed objective; no silent curvature substitution. |
| `ds.vertLASSOIter` | Repeated protected score/Hessian GLMs plus client proximal updates | Quarantined; exact slope-binomial route fails closed before DSI | Materialize a signed `binomial_lasso_design_grams` majorant over exactly the raw standardized score design, then add a whole-path workload contract and KKT/reference validation. The clipped/normalised Gaussian capsule cannot certify that Gram. |
| Inference helpers (`Confint`, `Wald`, `Contrast`, `LR`) | Client algebra and, for LR, additional fits | Inherits supplied fits | Validate finite covariance/nesting and report inherited—not upgraded—numeric status. LR needs a composite certificate over both fits. |
| Epidemiology/standardization helpers | Client algebra over authorized aggregates | Inherits input | Denominator/zero/domain checks and explicit uncertainty scope; no new Ring claim. |
| `ds.vertDesc` | Compatibility data frame over `ds.vertDPDescribe` | Separate integer DP contract, provisional | Uses one explicit custodian-owned sticky capsule artifact. Fixed-grid mechanism regions and public grid resolution are reported; sampling uncertainty is excluded. |
| `ds.vertDP*` | Integer DP mechanisms and client post-processing | Separate DP contract | Keep Ring numeric and DP guarantees distinct; calibration covers mechanism noise, not sampling error. |
| PSI/security/status endpoints | No statistical estimator | Not applicable | Do not attach numeric certification fields. |

Aliases have exactly the status of their canonical route. The executable source
of truth is `ds.vertMethodStatus()`: every registered method currently reports
both `may_report_numerically_certified = FALSE` and
`currently_numerically_certified = FALSE`. Promoting a future workload would
require an explicit audited registry and implementation change.

## Fail-closed promotion rule

A method can move from `unattested_result` only when all of the following are
true for that exact route and deployment: authenticated custodian policy body;
workload-specific adapter probe; canonical encoding; checked raw products and
accumulators; exact or rigorously bounded truncation/comparison/public-scalar
operations; runtime input and intermediate attestations bound to
policy/session/data/variables; approximation-domain proof; finite output and
identifiability postconditions; a conservative aggregate and estimator error
bound; and E2E differential tests against an agreed plaintext reference. No
single primitive test, ring width, or input bound is sufficient by itself.
