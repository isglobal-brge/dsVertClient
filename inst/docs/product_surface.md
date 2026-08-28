# dsVertClient audited method surface

Status: development audit, 2026-08-11.

The release contract is available at runtime through
`ds.vertMethodStatus()`. It classifies every public analysis entry point as
`promoted`, `provisional`, `compatibility`, or `quarantine`. A status describes
both statistical maturity and disclosure behaviour; it is not a cryptographic
certification.

## Available, audited conditional routes

This table includes promoted and provisional routes; the runtime registry is
the authoritative maturity classification.

| Route | Defensible scope | Important condition |
|---|---|---|
| `ds.getIdentityPks()` | Public Ed25519 identity discovery | Deployment must validate persistence, rotation and pin mismatch handling |
| `ds.vertConfint(type = "mechanism")` | Certificate-revalidated client post-processing | A current signed same-owner Gaussian Synopsis has a finite inverse-norm margin; the result is a simultaneous DP-mechanism region, never a sampling interval |
| `ds.vertWald()`, `ds.vertContrast()` | Client-only algebra | The supplied fit must be converged, unpenalised and expose a valid covariance matrix |
| `ds.vertEpi2x2()` | Client-only measures from an authorised 2-by-2 table | Caller-supplied matrices have no DataSHIELD provenance |
| `ds.vertMantelHaenszel()` | Client-only common odds ratio and classical MH inference from authorised stratified 2-by-2 aggregates | Bare inputs are caller-attested; orientation, strata and common-odds-ratio assumptions must be valid |
| `ds.vertDirectStandardization()` | Client-only directly standardised rates | Inputs must already be disclosure-authorised aggregates with compatible strata |
| `ds.vertIndirectStandardization()` | Client-only SMR/SIR | Expected counts and target population must be scientifically valid |
| `ds.vertMethodStatus()` | Read-only maturity registry | Metadata is conservative guidance, not formal verification |
| `ds.vertDPStatus()` | Read-only joint-DP capsule-status v5/accounting handshake | All peers must agree on domain, fixed per-capsule privacy parameters, `N`, adjacency, complete pinset, designated peers and one stable accountant namespace; the allocator-commit history can deny a new capsule at `N`, while `request_limit` is false and exact replay remains free; namespace uniqueness across reconfiguration is a custodial assumption, not automatic enforcement or migration |
| `ds.vertDPCapsulePlan()` | Data-free signed dry-run of the immutable capsule workload | It accesses no protected snapshot, creates no release and has zero privacy cost; it cannot change the custodian workload |
| `ds.vertDPCalibrate()` | Data-free fixed-capsule noise calibration | Quantifies mechanism noise, not sampling error or a universally optimal epsilon; operation history never degrades the calibration |
| `ds.vertDPCount()` | Canonical signed current-snapshot privacy-unit Count: sticky exact-MPC add/remove DP or public fixed-cohort K-consensus | Privacy is per signed artifact and distinct artifacts compose; no finite global composition claim is made. Add/remove assumes signed bounds, secret persistent seeds and one non-colluding pinned authority; the fixed value has zero sensitivity |
| `ds.vertDPContingency()` | Fixed-domain one-contribution DP histogram | Ordinary chi-square/Fisher p-values are not calibrated for noisy counts |
| `ds.vertDPFrequency()` | One fixed-domain categorical marginal from a canonical signed sticky two-authority Ring128 artifact | The source owner is explicit; PSI-v4 keeps the first eligible source record per unit, missing/out-of-domain factor values contribute zero, and count/proportion regions cover mechanism noise only |
| `ds.vertDPFrequencyInference()` | Zero-call conservative population-proportion regions from a validated DP frequency release | Requires iid privacy units and scientifically ignorable exclusions; Bonferroni/Clopper--Pearson regions can be wide and provide no p-value |
| `ds.vertChisq()` | Client-only mechanism-aware independence bootstrap from one signed DP contingency release | Plug-in asymptotic calibration; not an ordinary chi-square reference law |
| `ds.vertFisher()` | Client-only 2-by-2 conditional hypergeometric plug-in bootstrap from one signed DP contingency release | Not Fisher-exact for confidential data; no conditional CMLE/CI; Gaussian calibration is not yet certified |
| `ds.vertDPMeanVar()` | Bounded DP moments | Clipping and privacy noise must be included in interpretation |
| `ds.vertDPQuantile()`, `ds.vertDPMedian()` | Client-only fixed-grid bins and intervals from one validated DP describe release | No exact sample quantile or within-bin interpolation is claimed; mechanism regions exclude sampling uncertainty |
| `ds.vertDPGaussian()` | Bounded same- or cross-owner complete-case Gaussian coefficients from one signed sticky capsule; cross-owner products use the fixed-shape exact-GC artifact | Population-sampling inference is unavailable; ridge is explicit and changes the estimand |
| `ds.vertDPDiagnostic2x2()` | Client-only diagnostic accuracy from an authorised disease-by-test DP 2-by-2 table | Positive disease row and test column must be explicit; regions cover mechanism noise, not sampling uncertainty |
| `ds.vertDPDiagnostic2x2Inference()` | Conservative joint DP-mechanism plus exact-binomial sampling confidence regions for diagnostic accuracy | Assumes iid privacy units from one joint disease/test population; regions can be wide or unbounded and provide no p-value |
| `ds.vertDPROC()` | Client-only threshold ROC curve and tie-adjusted finite-snapshot AUC from one ordered DP table | Bin order/direction must be public; simultaneous regions cover mechanism noise, not population sampling uncertainty |
| `ds.vertDPEpi2x2()` | Client-only risks, association measures, attributable fractions and typed NNB/NNH from one DP 2-by-2 table | Regions can be unbounded at zero denominators and exclude sampling uncertainty |
| `ds.vertDPEpi2x2Inference()` | Conservative joint DP-mechanism plus exact-binomial sampling confidence regions from the same authorised 2-by-2 release | Assumes iid privacy units under the declared joint exposure/outcome model; Bonferroni/Clopper-Pearson regions can be wide or unbounded and provide no p-value |
| `ds.vertDPPrevalenceRatio()` | Zero-call cross-sectional prevalence, prevalence-difference/ratio/odds-ratio, attributable-fraction and reciprocal-difference naming view of `ds.vertDPEpi2x2()` | Exposed and prevalent orientations are mandatory; all numbers and mechanism regions are unchanged, and study design is caller-declared rather than inferred from the table |
| `ds.vertDPPrevalenceRatioInference()` | The same cross-sectional naming view of the conservative joint DP-mechanism plus exact-binomial regions from `ds.vertDPEpi2x2Inference()` | It preserves the source coverage contract exactly; the iid cross-sectional sampling model must be scientifically justified and adds no causal claim |
| `ds.vertDPMantelHaenszel()` | Client-only common MH odds ratio from one fixed DP strata-by-four-cells table | Requires the one-global-cell add/remove contract (block L1 sensitivity 1) and public mapping; the mechanism region excludes sampling uncertainty and no CMH p-value is reported |
| `ds.vertDPDirectStandardization()` | Client-only directly standardized risk from one DP strata table and public weights | Public weights must match the fixed stratum domain; regions cover mechanism noise only |
| `ds.vertDPDirectStandardizationInference()` | Conservative joint DP-mechanism plus exact-binomial stratum sampling region for directly standardized risk | Conditions on confidential stratum sample sizes and treats public standard weights as fixed; no causal or transportability claim is added |
| `ds.vertDPCausalStandardization()` | Saturated public-stratum g-formula with treated/control risk contrasts and simultaneous DP-mechanism regions | Causal interpretation requires consistency, conditional exchangeability, positivity, no interference, correct public row mapping and valid fixed target weights; no propensity score is estimated |
| `ds.vertDPCausalStandardizationInference()` | The same g-formula with conservative joint DP-mechanism plus exact-binomial sampling regions | All identification assumptions remain external; Bonferroni/Clopper-Pearson regions can be wide or unbounded and no p-value is provided |
| `ds.vertDPIndirectStandardization()` | Client-only O/E ratio from an authorised DP strata table and public expected rates | Denominator-zero regions may be non-estimable or unbounded; no classical Garwood interval is used |
| `ds.vertDPIndirectStandardizationInference()` | Conservative joint DP-mechanism plus exact Poisson Garwood region for the O/E ratio | Requires an externally valid fixed expected total and Poisson total-count model; zero expected denominators return `[0, Inf]`, with no p-value, causal or transportability claim |
| `ds.vertDPRMST()` | Client-only fixed-grid RMST from one authorised DP survival release | Its simultaneous limits cover mechanism noise, not sampling uncertainty or continuous-time grid error |
| `ds.vertDPRMTL()` | Exact client-only restriction-width complement of RMST from the same authorised release | It makes no new DSI call or DP release; limits inherit the RMST mechanism coverage and exclude sampling uncertainty and continuous-time grid error |
| `ds.vertDPSurvivalContrast()`, `ds.vertDPRMSTContrast()` | Client-only survival/RMST differences and ratios from two compatible authorised releases | Distinct releases use the Bonferroni joint-confidence lower bound without assuming independence; zero denominators are typed, and regions exclude sampling uncertainty and continuous-time grid error |
| `ds.vertDPSurvivalQuantile()`, `ds.vertDPMedianSurvival()` | Client-only inversion of one authorised fixed-grid Kaplan--Meier curve and its simultaneous mechanism band | Beyond-horizon targets are typed without extrapolation; limits exclude sampling uncertainty and continuous-time grid error |

These functions make no new server disclosure beyond their stated input. The
epidemiological, frequency-inference, survival-quantile, survival-contrast,
RMST and RMTL helpers issue zero DSI
calls after their validated input release.

The fixed-grid `ds.vertDPDescribe()` and `ds.vertDPSurvival()` routes remain
provisional pending broader deployment/scientific validation. `ds.vertDesc()`
is a compatibility data-frame adapter over the former and shares that maturity
status; it requires an explicit custodian-owned `analysis_id` and cannot
request exact extrema or adaptive bins. Their precise
finite-snapshot estimands, mechanism-region semantics, replace-one
sensitivities, and non-claims are recorded in
[`dp_statistical_validity_20260801.md`](dp_statistical_validity_20260801.md).
The same-owner statistics remain local protected-data materializations. The
declared cross-owner categorical and Gaussian artifacts use purpose-bound,
fixed-transcript exact GC and enter the same one-opening joint-DP vector;
families without such a signed cross-owner artifact remain unavailable.

## Provisional routes

Formal capsule descriptives/correlation/PCA/contingency methods, the explicit
Gaussian GLM capsule adapter and selected client post-processing remain
provisional. Pinned fixed-capacity PSI
alignment and its count-free persistent attestation are promoted non-statistical
protocols. Unbounded/legacy categorical regression, legacy mutating MI and
the iterative LASSO route remain quarantined. Their
precise supported scope and limitation are returned by
`ds.vertMethodStatus()`.

Provisional does not mean unusable. It means that at least one production claim
is still open: live multi-host validation, fixed-point/truncation proof,
comparison-protocol hardening, output-composition control, complete inferential
coverage or a broader statistical validation matrix.

## Compatibility and quarantine

Compatibility names for client-only post-processing retain their documented
historical behaviour. Quarantined distributed names remain exported for API
discovery, but their legacy server endpoints are unregistered and unexported;
they fail locally with a stable migration error before connection discovery or
DSI submission. Current quarantined families include legacy Cox computation,
unbounded legacy NB2, multinomial and ordinal regression, AR(1) and robust
clustered GEE, GLMM-PQL, legacy mutating MI, iterative LASSO and the legacy
weight-column IPW contract. Use
`ds.vertMethodStatus(status = "quarantine")` for the authoritative list.

In particular:

- `ds.vertIPW()` is limited to the exact intercept-only identity or one
  categorical saturated-stratum IPW/g-formula identity over a signed
  treatment-by-outcome or stratum-treatment-by-outcome table. The categorical
  route supports ATE with fixed public target weights and ATT/ATC with target
  weights derived from that sticky signed table. ATT/ATC expose a
  DP-mechanism region only, not sampling inference. It never writes or
  releases weights, fits a propensity model, or supports continuous/multiple
  covariates or outcome regression.
- the signed LMM route is limited to its signed random-intercept, fixed-effect
  or finite random-slope artifacts; the signed binary GLMM route is limited to
  an `outcome ~ 1` moment projection or its signed additive finite
  random-intercept/one-or-two-random-slope grid.
  Gaussian exchangeable GEE is a model-based working-GLS post-processing of a
  matching signed random-intercept artifact; AR(1), robust/sandwich GEE and
  legacy PQL remain cluster-granular and quarantined;
- `ds.vertMI()` only post-processes strict-missing signed categorical Synopsis
  artifacts into MCAR marginals, a two-response joint pair, or an opt-in
  conditional categorical star model. It never writes source data, performs
  chained equations or supplies Rubin sampling inference;
- `ds.vertNBFullRegTheta(y ~ 1, frequency = ...)` post-processes one signed
  bounded count Frequency release into a no-inference NB2 method-of-moments
  result. With an `analysis_id`, it instead selects an additive bare-predictor
  coefficient/theta candidate from one same-owner signed finite likelihood
  grid; no free likelihood optimisation or inference is offered;
- multinomial and ordinal `analysis_id` routes likewise select one additive
  bare-predictor candidate from a same-owner signed finite likelihood grid.
  They never report covariance, standard errors, tests or sampling inference;
- interactions, transforms, cross-owner categorical designs and arbitrary
  candidate grids remain unavailable before DSI; a Gaussian correlation
  capsule is not substituted for their signed likelihood artifact;
- `ds.vertLASSOCV()` is a legacy name for information-criterion selection, not
  cross-validation and not a one-standard-error rule.

## Remote server surface

`dsVert/DESCRIPTION` is guarded by an exact allowlist test. The variable-shape
PSI closures, patient-index helper and exact `localCorDS` closure have been
deleted; the only registered PSI surface is the sixteen-method padded protocol.
Other unused primitives that could mutate analysis data, expose generic
session reductions, clear protected state or release cluster-granular legacy
summaries are neither registered nor exported. Functions retained internally
for source compatibility are not remotely invocable.

Registered low-level endpoints are purpose-bound phases of authenticated
protocol state machines: malformed, replay-conflicting or out-of-phase direct
calls fail instead of opening a statistic. Registration is not promotion of a
high-level estimator, but no legacy exact/adaptive analysis primitive remains
on the registered surface.

## Non-claims

The package does not claim any of the following:

- impossibility of reconstruction under unlimited adaptive exact queries;
- differential privacy for historical exact implementations (the installed
  gate keeps them unavailable; the guarantee is limited to the explicit
  purpose-bound DP release surface);
- malicious-peer security, robustness to host compromise, or collusion of all
  designated MPC peers;
- correctness against a custodian that lies about its own input: one malicious
  custodian can bias its contribution or abort, and the two designated compute
  peers can reconstruct their joint additive secrets if they collude;
- bit-exact real arithmetic for historical local fixed-point truncation, which
  remains only in unregistered/quarantined implementations;
- production validity of quarantined estimators.

Pinned peers remove analyst key substitution and prevent the relay from opening
peer ciphertext. They do not prevent reconstruction from a sufficiently rich
collection of authorised exact outputs. See `SECURITY.md` for the threat model
and claim boundaries.
