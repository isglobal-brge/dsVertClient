# dsVert capsule method migration matrix

Current surface audit: 2026-10-15 (implemented 1.4.2). The machine-readable source of truth is
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
exact-GC computation and one joint noisy opening. Promoted finite likelihood,
LASSO, Cox and grouped grids use their existing signed protected release paths.
Strict-missing MI and survival time/event/entry artifacts still require one
owner; neither has a compatible cross-owner producer. A mathematical non-identifiability, invalid
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
| Correlation/PCA: `ds.vertCor`, `ds.vertDPCor`, `ds.vert.cor`, `ds.vertPCA`, `ds.vert.pca` | One authenticated current same- or cross-owner `gaussian_models` Synopsis with an intercept and a fixed outcome/anchor; complete cases cover the whole registered artifact | All variable pairs, PSD projection and eigenstructure are post-processing. DPCor explicitly opts in with `artifact_family="gaussian_models"`; its default same-owner pairwise artifact is unchanged. Offline certificates need external/session trust and full validation; PCA never returns scores |
| Pearson/Yates compatibility surface: `ds.vertChisq`, `ds.vert.chisq` | One signed fixed-domain same- or cross-owner categorical Synopsis; mechanism-aware multinomial plug-in null and expected-cell policy | Validated adapter and DP-aware parametric bootstrap; asymptotic rather than finite-sample conditional inference |
| Fisher compatibility surface: `ds.vertFisher`, `ds.vert.fisher` | One signed same- or cross-owner 2x2 Synopsis; hypergeometric odds-ratio-one plug-in null with signed noise/clamp law and sampler-TV calibration | Validated DP-aware conditional bootstrap; not finite-sample Fisher-exact, no conditional CMLE/CI, and Gaussian categorical calibration remains uncertified |
| Cross-vertical categorical: `ds.vertChisqCross`, `ds.vert.chisq_cross` | signed fixed-domain cross-owner pair, fixed-capacity private one-hot inputs, one exact concatenated Ring128 multiplication, no-wrap certificate, and the released-mechanism-aware null | formal joint-DP capsule route implemented: the cross-signed allocator/registry gate precedes source access, exact shares remain private, and one-release chi-square/Fisher adapters include replay and tamper checks; inference remains asymptotic and the Gaussian categorical reference law is not yet certified |
| Bounded Gaussian: `ds.vertDPGaussian` | Signed same- or cross-owner complete-case model with public bounds, explicit intercept and DP `n`, Gram, cross-products and outcome square | Both layouts use current authenticated sticky Synopsis certificates. Cross-owner products use exact GC; Gaussian mechanism regions are post-processing when identifiable, with no population-sampling inference |
| GLM: `ds.vertGLM`, `ds.vert.glm` | explicit compatible Gaussian signed Synopsis; bounded additive binomial/Poisson finite likelihood grid; or a completed two-authority public certificate | an explicit Gaussian `dp_analysis_id` delegates to the same- or cross-owner `ds.vertDPGaussian` capsule. Binomial/Poisson `analysis_id` selects one signed finite grid and `formal_analysis_id` reads one completed certificate; default, legacy iterative and fresh-formal computation fail before DSI |
| Cox profile: `ds.vertCox`, `ds.vertCoxProfileNonDisclosive` | same-owner signed bounded Breslow partial-likelihood finite grid over explicit time/event/predictor bounds, or a completed formal Cox certificate | validated sticky finite-grid selection or read-only certificate post-processing; no arbitrary optimiser, covariance, standard errors, p-values, baseline hazard, prediction, residuals, strata or delayed-entry support |
| Cox dispatcher/discrete hazard: `ds.vert.cox`, `ds.vert.coxph`, `ds.vertCoxDiscreteNonDisclosive` | the same signed Cox grid/certificate or a separately labelled fixed-grid pooled-logistic hazard artifact | `analysis_id` selects the signed Cox grid and `formal_analysis_id` reads its certificate; discrete time reads a completed two-authority binomial certificate. Neither starts expansion or a release, has sampling inference, nor relabels pooled logistic as Cox PH |
| Negative binomial: `ds.vertNB`, `ds.vertNBFullRegTheta`, `ds.vert.nb`, `dp_nb_grid` | Same-owner signed Frequency intercept-only moments or signed finite NB2 beta/theta grids; cross grids use the complete certified negative log-PMF | Safe public NB is the certified front-door alias, never the historical body. Qualified formulas plus `analysis_id` discover and validate embedded cross contracts; bare formulas preserve existing same-owner behavior; no unrestricted optimizer or sampling inference |
| Multinomial: `ds.vertMultinom`, `ds.vertMultinomJoint`, `ds.vertMultinomJointNewton`, `ds.vert.multinom`, `dp_multinomial_grid` | Signed Frequency intercept-only fit or finite same-/cross-owner softmax candidate grid with fixed class/reference domain | Qualified formulas plus `analysis_id` validate the embedded cross contract; bare formulas preserve same-owner grids. Interactions/transforms, unconstrained optimization and sampling inference remain unavailable |
| Ordinal: `ds.vertOrdinal`, `ds.vertOrdinalJointNewton`, `ds.vert.ordinal`, `dp_ordinal_grid` | Signed Frequency thresholds or same-/cross-owner finite cumulative-logit candidates with fixed clinical category order | Primary `ds.vertOrdinal`/alias dispatch qualified formulas with `analysis_id` through the embedded cross contract. JointNewton remains Frequency-only. No unconstrained optimization or sampling inference |
| LMM: `ds.vertLMM`, `ds.vert.lmm` | one signed same-owner random-intercept `outcome ~ 1` moment artifact, fixed-effect ML/REML finite variance-ratio grid, or finite Gaussian random-slope candidate grid | validated sticky Synopsis post-processing; no unrestricted optimisation/covariance, standard errors or sampling inference |
| Cross-owner Cox: `dp_cox_grid` | Signed finite bounded Breslow risk-set and coefficient contract | Explicit authenticated entry; `ds.vertCox` and its profile aliases retain same-owner/completed-certificate scopes because their controls do not transparently match this contract |
| Cross-owner LMM/GLMM/GEE: `dp_lmm_grid`, `dp_glmm_grid`, `dp_gee_grid` | Finite random-intercept ML, signed GH5 binomial/Poisson variance grids, and fixed-working-correlation GEE contracts | Existing promoted authenticated entries; no implicit dispatch from same-owner front doors with different REML/random-slope/covariance semantics; no unrestricted optimizer or sampling inference |
| LMM K>=3 compatibility: `ds.vertLMM.k3` | same signed LMM artifact and at least three peers | compatibility wrapper, tested through K=5; accepts the artifact's declared ML/REML profile |
| GEE: `ds.vertGEE`, `ds.vert.gee` | signed finite binomial/Poisson grid or completed formal GLM certificate for independence, matching signed Gaussian random-intercept artifact for exchangeable working GLS, signed Gaussian AR(1) beta/rho grid with a strict within-cluster order, or same-owner binomial/Poisson robust-independence bread/meat grid | all point routes are validated sticky-artifact post-processing. Gaussian AR(1) and binomial/Poisson robust-independence grids use componentwise-clipped cluster scores and expose sandwich covariance only when DP bread is positive definite; they provide no standard errors, p-values or sampling inference |
| GLMM: `ds.vertGLMM`, `ds.vert.glmm` | one signed same-owner binary `outcome ~ 1` random-intercept moment artifact, a same-owner signed finite binary marginal-likelihood grid bound to an additive bare design and zero to three named random slopes, or a same-owner signed finite bounded-count Poisson random-intercept/one-to-three-random-slope grid | validated sticky Synopsis post-processing for population-average log-odds/observed-scale ICC or one finite-grid fit. PQL/Laplace/ML/REML, interactions, more than three random slopes, covariance and sampling inference remain unavailable |
| IPW: `ds.vertIPW`, `ds.vert.ipw` | a signed same- or cross-owner binary treatment-by-outcome table for `outcome ~ treatment`, `treatment ~ 1`, or a signed categorical stratum-treatment-by-outcome table with complete public row mapping for `treatment ~ stratum` | the exact intercept-only ATE or one-categorical-stratum saturated IPW/g-formula ATE identity uses fixed public weights. ATT/ATC derive target-arm weights from the same sticky DP table and report a simultaneous-DP-mechanism outer region only; no weights or propensity fit are released. Continuous/multiple covariates, stabilization/clipping and outcome regression still require a purpose-bound secure protocol |
| Multiple imputation: `ds.vertMI`, `ds.vert.mi` | strict-missing same-owner categorical marginals or pairs in one signed Synopsis release, plus deterministic release-bound completion draws; cross categorical exclusion semantics do not satisfy this contract | validated categorical MCAR marginal, two-response joint-pair and opt-in conditional-star completion; the historical mutating route, chained equations, MAR/MNAR models and Rubin inference remain unavailable |
| Post-hoc LASSO: `ds.vertLASSO`, `ds.vert.lasso` | Validated signed bounded Gaussian moments and explicit L1 levels | Current same-/cross-owner Gaussian certificates are revalidated before KKT-certified post-processing; historic same-owner inputs remain supported and historic cross capsules are rejected |
| Gaussian one-step compatibility: `ds.vertLASSO1Step`, `ds.vert.lasso_1step` | Same validated Gaussian artifact as the main LASSO path | Delegates to the authenticated deterministic Gaussian L1 path, including current cross-owner certificates; no private optimization or second release |
| Gaussian proximal LASSO: `ds.vertLASSOProximal`, `ds.vert.lasso_proximal` | Authenticated unpenalized same-/cross-owner Gaussian moments and normal equations | Certificate-revalidated KKT solver; input limitations remain and no new source access occurs |
| Iterative LASSO: `ds.vertLASSOIter`, `ds.vert.lasso_iter`, `dp_lasso_grid` | Signed Gaussian moments, existing same-owner finite binomial/Poisson path, or a separately signed cross-owner LASSO extension plus base contract | Gaussian uses `dp_analysis_id`; same-owner finite grids keep `analysis_id`. Cross binomial/Poisson paths additionally require `cross_contract` and `base_contract`, qualified formulas and exact signed lambda levels. The extension adds zero sensitivity and no second release; retired score-MPC stays sealed |
| Information-criterion LASSO: `ds.vertLASSOCV`, `ds.vert.lasso_cv` | Authenticated same-/cross-owner Gaussian moments and deterministic L1 path | Validated pseudo-information-criterion selection, not cross-validation or a one-standard-error rule; no new source call |
| Likelihood ratio: `ds.vertLR`, `ds.vert.lr` | For `type="mechanism"`, current authenticated same-/cross-owner Gaussian certificates for compatible nested models; preserved sampling routes require their historical valid input contracts | Gaussian route bounds nested RSS reduction by post-processing and reports no likelihood-ratio test or p-value |
| DP-mechanism interval: `ds.vertConfint`, `ds.vert.confint` | Current authenticated same-/cross-owner Gaussian Synopsis with simultaneous coordinate errors and a finite inverse-norm margin | Revalidated simultaneous mechanism outer region; no sampling confidence claim |
| Wald/contrast inference: `ds.vertWald`, `ds.vert.wald`, `ds.vertContrast`, `ds.vert.contrast` | Current same-/cross-owner Gaussian certificates for `type="mechanism"`; historical sampling routes retain compatible covariance requirements | Certificate-revalidated mechanism regions do not introduce sampling p-values or sampling covariance |
| Exact-input epidemiology: `ds.vertEpi2x2` | authorized 2x2 table; group/population risks, RR, OR, RD, exposed/population attributable fractions, and NNB/NNH | client algebra exists and inherits table provenance |
| Exact stratified epidemiology: `ds.vertMantelHaenszel` | authorised 2x2xK or strata-by-four-cells aggregate with explicit orientation; common MH odds ratio and conditional classical inference | client algebra exists and inherits input provenance; bare tables are caller-attested and inference remains model-dependent |
| Exact direct standardization: `ds.vertDirectStandardization` | authorized stratum-specific cases/person-time plus public reference weights; standardized rate | client algebra exists and inherits input provenance |
| Exact indirect standardization: `ds.vertIndirectStandardization` | authorized compatible strata aggregates; observed/expected standardized ratio | client algebra exists and inherits input provenance |
| DP Count: `ds.vertDPCount` | one canonical signed current-snapshot Count artifact; identity-seeded sticky scalar exact GC under add/remove adjacency, or an exact public K-consensus declaration under fixed-cohort replace-one adjacency | `count_operation_implemented`: independent of capsule allocation, lifetime accounting and SQLite; identical artifacts recompute one identical bounded release, distinct artifacts compose, and no finite global composition is claimed |
| Sticky Frequency producer: `ds.vertDPFrequency` | one canonical signed fixed-domain categorical vector from the explicit PSI-aligned source owner | `frequency_operation_implemented`: K peers compile and sign; only the source and one pinned secondary authority execute one sticky two-peer Ring128 release; K-2 witnesses retain no execution state |
| Capsule DP producers: `ds.vertDPContingency`, `ds.vertDPMeanVar`, `ds.vertDPCor`, `ds.vertDPDescribe`, `ds.vertDPGaussian`, `ds.vertDPSurvival` | one joint capsule vector covering categorical pairs, numeric moments/pairs, fixed histograms, signed Gaussian sufficient statistics, and fixed-grid survival as applicable | sticky-capsule adapters are implemented for the listed artifacts; declared cross-owner categorical tables use the pre-source cross-signed allocation gate and private exact-GC injection; unsupported cross-owner families remain explicit reserved states |
| DP frequency view: `ds.vertDPFrequencyInference` | one validated stateless sticky Frequency artifact from `ds.vertDPFrequency` | zero-cost client post-processing implemented; combines a simultaneous mechanism count box with conservative Bonferroni/Clopper--Pearson population-proportion regions and makes no p-value claim |
| DP describe views: `ds.vertDPQuantile`, `ds.vertDPMedian` | one validated capsule describe artifact with a fixed public histogram grid | zero-cost client post-processing implemented; returns only a bin and interval, with no exact-quantile or interpolation claim |
| DP survival views: `ds.vertDPKaplanMeier`, `ds.vertDPNelsonAalen`, `ds.vertDPCumulativeIncidence`, `ds.vertDPSurvivalQuantile`, `ds.vertDPMedianSurvival`, `ds.vertDPRMST`, `ds.vertDPRMTL`, `ds.vertDPSurvivalContrast`, `ds.vertDPRMSTContrast`, `ds.vertDPLogRank` | one validated same-owner time/event/entry survival artifact, or two compatible artifacts for contrasts/scores; no split-owner survival producer | validated capsule adapters and zero-cost post-processing implemented; survival quantiles do not extrapolate, RMTL is the exact public-interval complement of RMST, and two-release contrasts/scores use a Bonferroni joint-confidence lower bound; log-rank exposes a mechanism-only score region and descriptive plug-in variance, never a p-value |
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

The non-inferential public names are tracked explicitly by
`.dsvert_capsule_non_inference_exports()`: identity discovery, PSI alignment
and alignment checks, DP calibration/status/plan, security/numeric status and
the public method registry. Their union with the inferential inventory is tested to
equal `NAMESPACE` exactly. PSI alignment remains outside the inferential matrix
because it returns no statistic. Its sole implementation
is now the promoted fixed-capacity, pinned-peer route: K and a public capacity
bucket determine the transcript, while exact counts and patient-derived
commitments remain server-local.

The [1.4.2 ownership inventory](../../README.md#public-cross-owner-inventory-142)
is the complete public family map. Four new grid exports expose previously
validated routes; they do not change wire capabilities, signed arithmetic, the
1.4.1 provenance epoch contract or Go runtime bytes. The v1.3.0 LASSO campaign
covered binomial; Poisson shares code/unit validation without a separate
archived campaign. Actual new connector evidence is reported independently.
