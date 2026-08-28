# dsVert capsule method migration matrix

Audit date: 2026-08-11. The machine-readable source of truth is
`.dsvert_capsule_method_inventory()` in
`R/capsule_method_inventory.R`. This inventory is a migration and evidence
contract; it does not alter the public maturity registry in
`ds.vertMethodStatus()` and it does not claim that a planned capsule backend is
already available.

## Current route versus target contract

The inventory deliberately records two different facts:

- `same_capsule_replay_history_can_deny` is `FALSE` for every route. Replaying
  the same authenticated capsule release is post-processing and consumes no
  additional privacy unit.
- `new_capsule_reservation_history_can_deny` is `TRUE` exactly for the public
  formal joint-DP capsule routes. A new capsule burns one persistent lifetime
  unit at allocator commit, before protected-source access or sampling, and is
  rejected when the signed lifetime allowance is exhausted. Client-only
  post-processing and routes that do not yet use the formal capsule backend
  report `FALSE` on this axis.

The latter is not a statement that every planned statistical family is already
available. The workload has a deterministic server-local protected-data
materializer for admitted count, bounded numeric moments and histograms,
within-owner categorical marginals/pairs, fixed-grid survival coordinates, and
signed Gaussian sufficient statistics. The implemented cross-owner categorical
and Gaussian artifacts add confidential source-share transport, fixed-shape
exact-GC computation and one joint noisy opening; other cross-owner families
remain explicit reserved states. A mathematical non-identifiability, invalid
specification, failed attestation or infrastructure failure can still produce
a typed failure. The lifetime gate is a privacy-loss admission boundary, not a
request quota: unlimited authenticated replay of the same release remains
available.

## Machine-readable axes

- `current_route_status`: what the exported method does today and whether its
  output contract is exact/granular, allocator-backed DP, local-only, or a
  postprocessor of a legacy input. `formal_sticky_count_artifact` identifies
  the independent canonical Count release and is not a capsule status.
- `migration_feasibility`: whether migration needs a capsule artifact, a new
  secure protocol, an attested capsule input, or a capsule adapter around
  already implemented post-processing. `count_operation_implemented` records
  that Count has completed this migration independently of the capsule path.
- `artifact_implementation_state`: whether the artifact is only planned,
  reserved, absent, or waiting for an input adapter.
- `inference_implementation_state`: separates existing numerical/client
  algebra from inference that still needs a capsule backend, secure redesign,
  DP-aware null law or validated capsule adapter.
- `alias_kind`: distinguishes an exact compatibility alias, a semantic wrapper,
  a canonical entry point and a deprecated subroute.

In particular, `postprocess_implemented_requires_capsule_adapter` means that
the algebra exists but the function does not yet validate the new capsule
schema. It must not be read as “direct capsule support”.

## Migration matrix

| Family and preserved public names | Required target artifact and inferential contract | Honest migration state |
|---|---|---|
| Descriptive: `ds.vertDesc`, `ds.vert.desc` | admitted count, bounded numeric moments and fixed numeric histograms for quantiles from one explicit custodian-owned `analysis_id` | validated adapter over the implemented sticky joint DP capsule; variable selection and quantile probabilities are client-only post-processing; exact extrema and adaptive ranges/bins are rejected |
| Correlation/PCA: `ds.vertCor`, `ds.vert.cor`, `ds.vertPCA`, `ds.vert.pca` | one signed same-owner or cross-owner `gaussian_models` artifact with a joint secret complete-case mask; PSD projection and spectral-stability rules | validated capsule adapter implemented; one capsule workflow, no pairwise fallback, no PCA scores, and zero extra DSI or privacy cost after the release |
| Pearson/Yates compatibility surface: `ds.vertChisq`, `ds.vert.chisq` | one signed fixed-domain same-owner categorical capsule; mechanism-aware multinomial plug-in null and expected-cell policy | validated capsule adapter and DP-aware parametric bootstrap implemented; explicitly asymptotic rather than finite-sample conditional |
| Fisher compatibility surface: `ds.vertFisher`, `ds.vert.fisher` | one signed fixed-domain same-owner 2x2 capsule; hypergeometric odds-ratio-one plug-in null with the signed noise/clamp law and sampler-TV calibration | validated capsule adapter and DP-aware conditional bootstrap implemented; explicitly not finite-sample Fisher-exact, no conditional CMLE/CI, and Gaussian calibration remains uncertified |
| Cross-vertical categorical: `ds.vertChisqCross`, `ds.vert.chisq_cross` | signed fixed-domain cross-owner pair, fixed-capacity private one-hot inputs, one exact concatenated Ring128 multiplication, no-wrap certificate, and the released-mechanism-aware null | formal joint-DP capsule route implemented: the cross-signed allocator/registry gate precedes source access, exact shares remain private, and one-release chi-square/Fisher adapters include replay and tamper checks; inference remains asymptotic and the Gaussian categorical reference law is not yet certified |
| Bounded Gaussian: `ds.vertDPGaussian` | one signed same- or cross-owner complete-case model, public outcome/predictor bounds, explicit intercept and DP `n`, `X'X`, `X'y`, `y'y`; released-design diagnostics and mechanism certificate | formal sticky-capsule route implemented for both ownership layouts; cross-owner products use the signed fixed-shape exact-GC artifact, while nonlinear coefficient regions and population-sampling inference remain unavailable |
| GLM: `ds.vertGLM`, `ds.vert.glm` | explicit compatible Gaussian signed Synopsis; bounded additive binomial/Poisson finite likelihood grid; or a completed two-authority public certificate | an explicit Gaussian `dp_analysis_id` delegates to the same- or cross-owner `ds.vertDPGaussian` capsule. Binomial/Poisson `analysis_id` selects one signed finite grid and `formal_analysis_id` reads one completed certificate; default, legacy iterative and fresh-formal computation fail before DSI |
| Cox profile: `ds.vertCox`, `ds.vertCoxProfileNonDisclosive` | partial likelihood, risk sets, score/Hessian, baseline hazard and covariance under explicit ties, strata and delayed-entry contracts | new capsule backend required; current route is limited to its implemented time/status semantics |
| Cox dispatcher/discrete hazard: `ds.vert.cox`, `ds.vert.coxph`, `ds.vertCoxDiscreteNonDisclosive` | Cox artifacts above or a separately labelled fixed-grid pooled-logistic hazard artifact | a custodian-configured, completed two-authority binomial certificate bound to the `Surv()` formula and fixed grid is readable; it cannot start expansion or a release, has no sampling inference, and is never relabelled as Cox PH |
| Negative binomial: `ds.vertNBFullRegTheta`, `ds.vert.nb` | one validated sticky bounded non-negative integer Frequency release for an intercept-only log-mean and NB2 method-of-moments dispersion, or one same-owner signed finite `(beta, theta)` likelihood grid with bounded integer outcome and public predictor bounds | `y ~ 1` post-processing is available without a new DSI call; one signed grid selects additive bare-predictor coefficients and theta. Unconstrained likelihood optimization, covariance and sampling inference remain unavailable |
| Multinomial: `ds.vertMultinom`, `ds.vertMultinomJoint`, `ds.vertMultinomJointNewton`, `ds.vert.multinom` | one validated sticky categorical Frequency release for `y ~ 1`, or one same-owner signed finite softmax likelihood grid bound to a categorical domain, public predictor bounds and additive bare design terms | all listed names post-process the Frequency release for intercept-only log-odds or select one signed-grid additive fit. Interactions/transforms, unconstrained optimization, covariance and sampling inference remain unavailable |
| Ordinal: `ds.vertOrdinal`, `ds.vertOrdinalJointNewton`, `ds.vert.ordinal` | one validated sticky categorical Frequency release and a complete clinical order for `y ~ 1`, or one same-owner signed finite cumulative-logit grid with ordered category domain, public predictor bounds and additive bare design terms | listed names post-process the Frequency release for thresholds or select one signed-grid additive fit. Interactions/transforms, unconstrained optimization, covariance and sampling inference remain unavailable |
| LMM: `ds.vertLMM`, `ds.vert.lmm` | one signed same-owner random-intercept `outcome ~ 1` moment artifact, fixed-effect ML/REML finite variance-ratio grid, or finite Gaussian random-slope candidate grid | validated sticky Synopsis post-processing; no unrestricted optimisation/covariance, standard errors or sampling inference |
| LMM K>=3 compatibility: `ds.vertLMM.k3` | same signed LMM artifact and at least three peers | compatibility wrapper, tested through K=5; accepts the artifact's declared ML/REML profile |
| GEE: `ds.vertGEE`, `ds.vert.gee` | signed finite binomial/Poisson grid or completed formal GLM certificate for independence, or matching signed Gaussian random-intercept artifact for exchangeable working GLS; a future robust route needs contribution-bounded cluster score/meat and DP-aware covariance | independence and Gaussian model-based exchangeable point estimates are validated sticky-artifact post-processing. Fresh-formal computation, AR(1), sandwich covariance and inference still require a protected kernel |
| GLMM: `ds.vertGLMM`, `ds.vert.glmm` | one signed same-owner binary `outcome ~ 1` random-intercept moment artifact, or a same-owner signed finite marginal-likelihood grid bound to an additive bare design and zero, one to three named random slopes | validated sticky Synopsis post-processing for population-average log-odds/observed-scale ICC or one finite-grid fit. PQL/Laplace/ML/REML, interactions, more than three random slopes, covariance and sampling inference remain unavailable |
| IPW: `ds.vertIPW`, `ds.vert.ipw` | a signed binary treatment-by-outcome table for `outcome ~ treatment`, `treatment ~ 1`, or a signed categorical stratum-treatment-by-outcome table with complete public row mapping for `treatment ~ stratum` | the exact intercept-only ATE or one-categorical-stratum saturated IPW/g-formula ATE identity uses fixed public weights. ATT/ATC derive target-arm weights from the same sticky DP table and report a simultaneous-DP-mechanism outer region only; no weights or propensity fit are released. Continuous/multiple covariates, stabilization/clipping and outcome regression still require a purpose-bound secure protocol |
| Multiple imputation: `ds.vertMI`, `ds.vert.mi` | strict-missing categorical marginals or pairs in one signed Synopsis release, plus deterministic release-bound completion draws | validated categorical MCAR marginal, two-response joint-pair and opt-in conditional-star completion; the historical mutating route, chained equations, MAR/MNAR models and Rubin inference remain unavailable |
| Post-hoc LASSO: `ds.vertLASSO`, `ds.vert.lasso` | attested authorized coefficients | client algebra exists and inherits input limitations |
| Quadratic-surrogate LASSO: `ds.vertLASSO1Step`, `ds.vert.lasso_1step` | attested coefficients plus covariance/Fisher information | client solver exists; capsule input attestation absent |
| Gaussian proximal LASSO: `ds.vertLASSOProximal`, `ds.vert.lasso_proximal` | attested unpenalized Gaussian fit, admitted count, numeric moments and centered normal equations | client solver exists; capsule input attestation absent |
| Iterative LASSO: `ds.vertLASSOIter`, `ds.vert.lasso_iter` | one signed same-owner Gaussian Synopsis artifact with bounded moments, provenance validation and a deterministic KKT certificate for every L1 penalty | Gaussian is promoted through explicit `dp_analysis_id` and creates no second release; binomial and Poisson still fail closed pending their purpose-bound score-design artifacts and whole-path contracts |
| Information-criterion LASSO: `ds.vertLASSOCV`, `ds.vert.lasso_cv` | attested coefficients, covariance/Fisher information, n and model deviance | client selector exists; it is information-criterion selection, not cross-validation |
| Likelihood ratio: `ds.vertLR`, `ds.vert.lr` | two attested, nested, converged, unpenalized binomial/Poisson fits with canonical unweighted deviance, identical cohort/missingness/offset | client algebra exists and inherits input limitations |
| DP-mechanism interval: `ds.vertConfint`, `ds.vert.confint` | current signed same-owner Gaussian Synopsis, its simultaneous coordinate-error certificate and a finite released inverse-norm margin | certificate-revalidated simultaneous outer region implemented; it is not a sampling confidence interval |
| Wald inference: `ds.vertWald`, `ds.vert.wald`, `ds.vertContrast`, `ds.vert.contrast` | attested coefficients and valid compatible covariance | client algebra exists and inherits input limitations |
| Exact-input epidemiology: `ds.vertEpi2x2` | authorized 2x2 table; group/population risks, RR, OR, RD, exposed/population attributable fractions, and NNB/NNH | client algebra exists and inherits table provenance |
| Exact stratified epidemiology: `ds.vertMantelHaenszel` | authorised 2x2xK or strata-by-four-cells aggregate with explicit orientation; common MH odds ratio and conditional classical inference | client algebra exists and inherits input provenance; bare tables are caller-attested and inference remains model-dependent |
| Exact direct standardization: `ds.vertDirectStandardization` | authorized stratum-specific cases/person-time plus public reference weights; standardized rate | client algebra exists and inherits input provenance |
| Exact indirect standardization: `ds.vertIndirectStandardization` | authorized compatible strata aggregates; observed/expected standardized ratio | client algebra exists and inherits input provenance |
| DP Count: `ds.vertDPCount` | one canonical signed current-snapshot Count artifact; identity-seeded sticky scalar exact GC under add/remove adjacency, or an exact public K-consensus declaration under fixed-cohort replace-one adjacency | `count_operation_implemented`: independent of capsule allocation, lifetime accounting and SQLite; identical artifacts recompute one identical bounded release, distinct artifacts compose, and no finite global composition is claimed |
| Sticky Frequency producer: `ds.vertDPFrequency` | one canonical signed fixed-domain categorical vector from the explicit PSI-aligned source owner | `frequency_operation_implemented`: K peers compile and sign; only the source and one pinned secondary authority execute one sticky two-peer Ring128 release; K-2 witnesses retain no execution state |
| Capsule DP producers: `ds.vertDPContingency`, `ds.vertDPMeanVar`, `ds.vertDPCor`, `ds.vertDPDescribe`, `ds.vertDPGaussian`, `ds.vertDPSurvival` | one joint capsule vector covering categorical pairs, numeric moments/pairs, fixed histograms, signed Gaussian sufficient statistics, and fixed-grid survival as applicable | sticky-capsule adapters are implemented for the listed artifacts; declared cross-owner categorical tables use the pre-source cross-signed allocation gate and private exact-GC injection; unsupported cross-owner families remain explicit reserved states |
| DP frequency view: `ds.vertDPFrequencyInference` | one validated stateless sticky Frequency artifact from `ds.vertDPFrequency` | zero-cost client post-processing implemented; combines a simultaneous mechanism count box with conservative Bonferroni/Clopper--Pearson population-proportion regions and makes no p-value claim |
| DP describe views: `ds.vertDPQuantile`, `ds.vertDPMedian` | one validated capsule describe artifact with a fixed public histogram grid | zero-cost client post-processing implemented; returns only a bin and interval, with no exact-quantile or interpolation claim |
| DP survival views: `ds.vertDPKaplanMeier`, `ds.vertDPNelsonAalen`, `ds.vertDPCumulativeIncidence`, `ds.vertDPSurvivalQuantile`, `ds.vertDPMedianSurvival`, `ds.vertDPRMST`, `ds.vertDPRMTL`, `ds.vertDPSurvivalContrast`, `ds.vertDPRMSTContrast`, `ds.vertDPLogRank` | one validated capsule survival artifact, or two compatible artifacts for contrasts/scores | validated capsule adapters and zero-cost post-processing implemented; survival quantiles do not extrapolate, RMTL is the exact public-interval complement of RMST, and two-release contrasts/scores use a Bonferroni joint-confidence lower bound; log-rank exposes a mechanism-only score region and descriptive plug-in variance, never a p-value |
| DP epidemiology/diagnostic views: `ds.vertDPEpi2x2`, `ds.vertDPEpi2x2Inference`, `ds.vertDPMantelHaenszel`, `ds.vertDPDiagnostic2x2`, `ds.vertDPDiagnostic2x2Inference`, `ds.vertDPROC`, `ds.vertDPDirectStandardization`, `ds.vertDPDirectStandardizationInference`, `ds.vertDPCausalStandardization`, `ds.vertDPCausalStandardizationInference`, `ds.vertDPIndirectStandardization`, `ds.vertDPIndirectStandardizationInference` | one validated capsule 2x2, strata-by-four-cells, ordered diagnostic-bin or strata-treatment-by-outcome artifact; combined inference additionally requires its declared binomial or Poisson sampling model; causal standardisation requires fixed public target weights and explicit identification assumptions; DP MH requires the one-global-cell add/remove contract with block L1 sensitivity 1 | validated capsule adapters and zero-cost post-processing implemented; combined inference uses exact-binomial or Poisson-Garwood union-bound regions, while the DP MH route has no classical CMH p-value |

`ds.vertDPEpi2x2` returns group/population risks, risk difference, risk and odds
ratios, exposed and population attributable fractions, and a direction-typed
number needed to benefit/harm. `ds.vertDPEpi2x2Inference` adds conservative
joint mechanism-plus-sampling confidence regions without another server call
or privacy release. `ds.vertDPDiagnostic2x2Inference` does the same for the
explicitly oriented diagnostic metrics using six exact-binomial base regions.
`ds.vertDPDirectStandardizationInference` combines positive-weight stratum
intervals under fixed public standard-population weights.
`ds.vertDPIndirectStandardizationInference` envelopes an exact Poisson Garwood
interval over the signed integer count box and fixed public expected rates;
possible zero expected totals return the full non-negative parameter space.
`ds.vertDPCausalStandardization` and its inference companion apply the same
bounded-table logic to a saturated public stratum-treatment design, while
making every causal identification assumption and the absence of a
propensity/doubly-robust model explicit.
`ds.vertDPMantelHaenszel` returns a common MH
odds ratio plus a conservative simultaneous mechanism region from one fixed
strata-by-cells table. `ds.vertDPROC` returns a threshold ROC curve and
tie-adjusted finite-snapshot AUC for a fixed public score-bin order. Their
base finite-snapshot, MH and ROC regions quantify mechanism uncertainty
only; sampling intervals, tests or population confidence bands outside the
explicit combined 2x2, direct-standardisation and causal-standardisation
contracts require separate validation.

## Classified remote-call evidence

`legacy_remote_call_evidence` is a list-column of data frames with `call`,
`component`, `release_class` and a source-level evidence statement. The classes
prevent an encrypted peer relay from being conflated with an analyst-visible
exact aggregate:

- `plaintext_exact_aggregate`: exact moments, counts, tables or per-cluster
  summaries returned directly to the analyst.
- `client_reconstructed_exact_statistic`: a guarded exact final statistic
  reconstructed by the client.
- `share_reconstructed_by_client`: a response is a share, but the current route
  retrieves complementary outputs and combines them.
- `opaque_peer_ciphertext`: the client relays a blob sealed to a pinned peer and
  cannot inspect the protected payload.
- `public_metadata`: dimensions, receipts, domains or attestations.

Verified evidence includes the previously omitted exact outputs from
`glmStandardizeDS`; the NB bootstrap/profile and full share route (including
`dsvertLocalMomentsDS`, `dsvertNBProfileSumsDS` and
`dsvertNBPsiAggregateDS`); per-cluster `dsvertClusterZtZDS`; and the exact
`n_imputed`/`n_observed` values returned by `dsvertImputeColumnDS`. It also
records that `k2ShareWeightsDS` is recipient-pinned ciphertext and that the GEE
AR1 ordering payload is opaque to the relay, while separately classifying their
visible metadata. These are migration/audit facts, not an allowlist.

Any route classified as `legacy_granular_release_quarantine` is tested to have
a secure-redesign requirement, an unimplemented secure artifact and concrete
historical endpoint evidence. Those endpoints are now unregistered and
unexported, and their public frontdoors fail locally before DSI. The tests
intentionally do not freeze a particular number of
quarantined names: the classification must follow evidence and may change as a
route is replaced.

## Non-inferential exports

The remaining twelve public names are tracked explicitly by
`.dsvert_capsule_non_inference_exports()`: identity discovery, PSI alignment
and alignment checks, DP calibration/status/plan, security/numeric status and
the public method registry. Their union with the 95 inferential names is tested to
equal `NAMESPACE` exactly. PSI alignment remains outside the inferential matrix
because it returns no statistic. Its sole implementation
is now the promoted fixed-capacity, pinned-peer route: K and a public capacity
bucket determine the transcript, while exact counts and patient-derived
commitments remain server-local.
