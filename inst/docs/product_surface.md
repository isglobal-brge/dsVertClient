# dsVertClient audited method surface

Status: implemented 1.4.2 surface audit, 2026-10-15.

The complete ownership inventory and remaining same-owner limitations are in
[the README](../../README.md#public-cross-owner-inventory-142).

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
| `ds.vertConfint(type = "mechanism")` | Certificate-revalidated client post-processing | A current authenticated same- or cross-owner Gaussian Synopsis has a finite inverse-norm margin; the result is a simultaneous DP-mechanism region, never a sampling interval |
| `ds.vertWald()`, `ds.vertContrast()`, `ds.vertLR()` | Current Gaussian `type="mechanism"` regions are certificate-revalidated post-processing for either ownership layout | No sampling inference is added; LR bounds nested RSS reduction. Preserved sampling routes still require their documented valid covariance/deviance inputs |
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
| `ds.vertCor()`, `ds.vertDPCor(artifact_family="gaussian_models")`, `ds.vertPCA()` and aliases | All variable pairs and PCA eigenstructure from one same- or cross-owner Gaussian Synopsis | An intercept, fixed registered outcome/anchor and whole-artifact complete cases are required; authenticated offline certificates make zero DSI calls; no scores are released |
| `dp_nb_grid()`, `dp_multinomial_grid()`, `dp_ordinal_grid()` | Promoted signed finite grids across K custodians with exactly two compute/noise authorities | Matching qualified formulas dispatch through `analysis_id` in safe `ds.vertNB`/`ds.vertNBFullRegTheta`, `ds.vertMultinom` and `ds.vertOrdinal`; no unconstrained optimizer or sampling inference |
| `dp_lasso_grid()`, non-Gaussian `ds.vertLASSOIter()` | Finite signed L1 candidate selection over one authenticated binomial/Poisson base loss release | The front door needs separately signed `cross_contract` and `base_contract`; lambda can only select exact signed levels; no second private release |
| `dp_cox_grid()`, `dp_lmm_grid()`, `dp_glmm_grid()`, `dp_gee_grid()` | Promoted cross-owner finite Breslow, random-intercept ML, GH5, and fixed-working-correlation grids respectively | Explicit signed contracts remain the cross entry points; historical front-door risk-set/REML/random-slope/covariance controls are not silently mapped onto different estimands |
| Gaussian `ds.vertLASSO*()` | KKT-certified paths or pseudo-information-criterion selection from validated current same- or cross-owner Gaussian certificates | Existing same-owner inputs remain valid; historic cross capsules and caller-supplied unverified moments are refused; the CV name does not mean cross-validation |
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
| `ds.vertDPLogRank()` | Client-only fixed-grid all-cause or cause-specific log-rank score from two compatible authorised releases | It returns a conservative DP-mechanism score region; its plug-in variance is descriptive only, the authenticated releases do not prove disjoint cohorts, and no test, p-value, sampling interval or PH effect is claimed |
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
fixed-transcript exact GC and enter the same one-opening joint-DP vector.
The promoted likelihood/grouped grids use their existing authenticated protected
release routes. The default pairwise correlation artifact, MI strict-missing
artifacts and time/event/entry survival histogram still require same-owner
inputs; no new producer is implied by 1.4.2.

## Provisional routes

Formal capsule descriptives/correlation/PCA/contingency methods, the explicit
Gaussian GLM capsule adapter and selected client post-processing remain
provisional. Pinned fixed-capacity PSI
alignment and its count-free persistent attestation are promoted non-statistical
protocols. Unbounded/legacy categorical regression, legacy mutating MI and
the retired iterative score-MPC LASSO route remain quarantined; the signed
Gaussian and finite-grid LASSO paths are admitted. Their
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
unbounded legacy NB2, multinomial and ordinal regression, arbitrary legacy
AR(1) and clustered-GEE solvers, GLMM-PQL, legacy mutating MI, retired score-MPC LASSO and the legacy
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
  or finite random-slope artifacts; the signed GLMM route is limited to a
  binary `outcome ~ 1` moment projection, a signed additive binary finite
  random-intercept/one-to-three-random-slope grid, or a signed additive
  bounded-count Poisson random-intercept/one-to-three-random-slope grid.
  Gaussian exchangeable GEE is a model-based working-GLS post-processing of a
  matching signed random-intercept artifact. A distinct signed finite
  beta/rho artifact supports Gaussian AR(1) working-GLS with a strict
  within-cluster order. A distinct bounded AR(1) artifact supports
  componentwise-clipped cluster-score sandwich covariance for at most three
  additive predictors when its DP bread remains positive definite. It has no
  standard errors, p-values or sampling inference. Matching same-owner
  binomial/Poisson robust-independence grids likewise support one categorical
  cluster id, no order column, at most three additive predictors and a
  clipped-score sandwich covariance only when DP bread is positive definite;
  legacy PQL and arbitrary clustered-GEE solvers remain quarantined. Cross-owner
  `dp_lmm_grid`, `dp_glmm_grid`, `dp_gee_grid` and `dp_cox_grid` expose their
  own signed finite-grid estimands; their historical front doors do not add
  automatic cross dispatch;
- `ds.vertMI()` only post-processes strict-missing signed categorical Synopsis
  artifacts into MCAR marginals, a two-response joint pair, or an opt-in
  conditional categorical star model. It never writes source data, performs
  chained equations or supplies Rubin sampling inference. Its current supported
  artifacts are same-owner: cross categorical exclusion rules do not satisfy the
  required strict missingness contract;
- `ds.vertNBFullRegTheta(y ~ 1, frequency = ...)` and safe alias `ds.vertNB`
  retain the signed bounded-count Frequency NB2 method-of-moments route.
  With `analysis_id`, bare formulas retain the existing same-owner finite
  likelihood grid; fully qualified additive formulas select the embedded
  authenticated cross-owner NB2 beta/theta contract through `dp_nb_grid`.
- Multinomial and ordinal primary `analysis_id` routes likewise admit the
  corresponding signed cross-owner grid with qualified additive formulas.
  Class/reference/order controls must exactly match the signed contract.
  `ds.vertOrdinalJointNewton` remains Frequency-only; it has no grid selector.
  No finite-grid route provides covariance, standard errors or sampling tests.
- Interactions, transforms, unsigned candidate grids and unrestricted
  optimizers remain unavailable. A Gaussian correlation artifact cannot replace
  a signed categorical likelihood contract.
- `ds.vertLASSOIter(..., cross_contract=..., base_contract=...)` exposes the
  separately signed binomial/Poisson finite L1 extension. Its base release owns
  protected computation and privacy spending; the public capacity-normalized
  penalty adds zero sensitivity. The v1.3.0 archive covers binomial; Poisson
  shares implementation/unit validation but has no separate archived campaign.
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
