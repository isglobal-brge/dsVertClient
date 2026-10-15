# dsVertClient - DataSHIELD Client for Vertically Partitioned Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-1.4.1-blue.svg)](https://isglobal-brge.github.io/dsVertClient/news/index.html)

## Overview

**dsVertClient** provides R workflows for privacy-aware analysis of vertically
partitioned data across DataSHIELD servers. It combines ECDH-PSI, Ring63 /
Ring127 additive sharing, dealer-free IKNP OT-extension preprocessing,
peer-to-peer authenticated encryption and Ed25519 server identity pinning.

Method maturity is intentionally explicit. Run `ds.vertMethodStatus()` before
choosing a route: its method-specific release contract and limitation are
authoritative. A `promoted` route is available only within its declared
numerical, disclosure and biomedical scope; historical routes outside that
scope remain unavailable rather than silently falling back.

All connected servers must expose dsVert's single immutable
`disclosure_safe` profile. There is no client or server switch that enables an
exact legacy release plane; methods still awaiting migration return an
explicit unavailable/migration error instead of falling back silently.

Pair with the server-side companion package [**dsVert**](https://github.com/isglobal-brge/dsVert).

## 1.4.1 provenance epochs and migration

Upgrade dsVert and dsVertClient together. PSI, Count, Frequency and Synopsis
reject mixed provenance protocols before private materialization; analyst
function arguments are unchanged. These releases expose a `stable_receipt`
that binds public analysis semantics and all participating authorization
epochs, including Synopsis dependencies. A new custodian authorization changes
that identity even when the authorized bytes are identical.

**Legacy generic aligned imports require one custodian migration before use.**
In the trusted server R process, after approving a descriptor for the current
frame, run this once per imported dataset (replace `id` with its privacy-unit
column):

```r
dsVert::dsvertAuthorizeSource(
  data = aligned, descriptor = getOption("dsvert.dp.datasets")$aligned,
  event_id = "migration-1.4.1", kind = "dataset", patient_column = "id")
```

The command is administrative and is not a DataSHIELD analyst endpoint.
Registered material and committed numeric answers stay frozen across source
edits and restarts until explicit custodian reauthorization. Retain the server
authorization store with the existing identity/noise/replay state. See the
[provenance and custodian workflow](inst/docs/provenance_epochs_1.4.1.md) for
first-load rules, reauthorization, stable versus operational receipts and
privacy limits.

## Quick Start

```r
library(DSI); library(DSOpal); library(dsVertClient)

# Connect to servers
conns <- DSI::datashield.login(logins = builder$build())
DSI::datashield.assign.table(conns, "D", "project.data")

# 1. Align records without mutating the source table
ds.psiAlign("D", "patient_id", "DA", na.action = "na.omit",
            datasources = conns)

# 2. Inspect the immutable server policy and signed capsule plan
ds.vertDPStatus(conns)
ds.vertDPCapsulePlan(conns)

# 3. Obtain one sticky bounded Gaussian capsule and reuse it for analysis
gaussian <- ds.vertDPGaussian("DA", "gaussian_primary",
                              datasources = conns)

# 4. Correlation + PCA from the same authorised Gaussian artifact
# The admin-declared Gaussian artifact may be same-owner or cross-owner. All
# variables use its one signed, secret complete-case mask.
cor <- ds.vertCor(
  "DA", variables = list(site_a = c("age", "bmi", "glu")),
  analysis_id = "gaussian_primary", datasources = conns)
pca <- ds.vertPCA(cor_result = cor)

DSI::datashield.logout(conns)
```

## Differentially private releases

The package can inspect—never change—the immutable DP policy and run a public
utility preview. The capability check is data-free; a 1.4.1 Synopsis plan may
read and freeze authorized source material during first registration, without
releasing statistics. Count independently compiles one canonical signed
current-snapshot artifact and uses sticky exact-MPC noise or signed fixed-cohort
K-consensus. Frequency independently compiles one signed fixed-domain artifact
and executes one sticky two-authority Ring128 release; the source owner must be
named explicitly and the remaining K-2 peers are compile-only witnesses. The
signed biomedical capsule covers table, moment,
Gaussian-model and survival artifacts. Its deterministic selector uses
the exact-GC one-draw Laplace route for a scalar vector, the scalable
two-complete-draw Laplace convolution for wider vectors, or the formal
fixed-work dyadic discrete-Gaussian backend when its exact server plan strictly
improves the certified simultaneous radius. The preview does not impersonate
that decision: the signed server manifest is authoritative. The 1.4.0
[exact v4 Laplace fallback](inst/docs/exact_laplace_fallback_v4.md) uses
unbounded draws and modular Ring128 addition, accepts zero mechanism delta
under ideal independent bits, and carries a separate positive utility-wrap
certificate. Keyed replay assumes pseudorandomness; v3 artifacts retain their
original finite-sampler guarantees.

```r
ds.vertDPCalibrate(capsule_epsilon = c(1, 3, 5),
                   peer_count = length(conns),
                   sensitivity = 1)
ds.vertDPStatus(conns)
ds.vertDPCapsulePlan(conns) # signed preview; may capture authorization, no release
# Canonical stateless sticky artifacts:
ds.vertDPCount("DA", datasources = conns)
frequency <- ds.vertDPFrequency(
  "DA", "exposure", server = "site_a", datasources = conns)
ds.vertDPFrequencyInference(frequency)
# Sticky joint-DP capsule post-processing routes:
ds.vertDPContingency("DA", "exposure", "outcome",
                     datasources = conns)
ds.vertDPMeanVar("DA", "age", datasources = conns)
ds.vertDPDescribe("DA", "primary", datasources = conns)
gaussian <- ds.vertDPGaussian("DA", "gaussian_primary",
                              datasources = conns)
surv <- ds.vertDPSurvival("DA", "primary", datasources = conns)
ds.vertDPKaplanMeier(surv)
ds.vertDPSurvivalQuantile(surv, probabilities = c(0.25, 0.5, 0.75))
ds.vertDPMedianSurvival(surv)
ds.vertDPRMST(surv, tau = 365)
ds.vertDPRMTL(surv, tau = 365)
# Given a second compatible authorised release:
# ds.vertDPSurvivalContrast(surv_treated, surv_control)
# ds.vertDPRMSTContrast(surv_treated, surv_control, tau = 365)
# ds.vertDPLogRank(surv_treated, surv_control, cause = "event")
```

The table wrapper itself intentionally reports no ordinary chi-square or
Fisher p-value: their usual reference distributions ignore DP noise.
`ds.vertChisq()` provides a mechanism-aware independence bootstrap, while
`ds.vertFisher()` provides a signed-mechanism 2-by-2 conditional
hypergeometric plug-in bootstrap. The latter is explicitly not Fisher-exact
for the confidential table and returns no classical conditional odds-ratio
interval. Both can reuse the same sticky release with zero additional DSI
calls and zero additional privacy cost.

Its public repeated-record policy is
`consistent_cell_else_exclude_v1`: identical valid cells count once,
conflicting valid cells make the unit contribute zero, and missing or
out-of-domain rows do not create a conflict. No exact exclusion count is
released; studies should report the resulting concordant-unit estimand.

The public DP release plane is the signed Synopsis protocol, not an
epsilon ledger per peer or per user operation. Pinned peers cross-sign each
immutable canonical snapshot/workload artifact and its sticky publication
state. Retry, restart and lost acknowledgement must replay the same durable
artifact without resampling; corrupted or inconsistent durable state fails
closed. Server secrets and authenticated ciphertext spools remain in
owner-only Rock storage, and no analyst receives a root or derived seed.

There is no request, rate, catalog, or finite lifetime admission limit. A
supported canonical artifact has its own declared DP scope; a distinct
artifact is a distinct analysis and may compose with it. We do not claim a
finite global privacy bound across an unlimited set of different analyses.
`ds.vertDPStatus()` and `ds.vertDPCapsulePlan()` verify signed Synopsis
bootstraps: they bind the pinset, two designated authorities, authorization
dependencies and immutable workload. They may read and durably capture
authorized material on first registration, and create no statistical release. The retired
vector-lifecycle ABI remains internal compatibility code and is not a public
release route.

`ds.vertSecurityStatus()` consumes security-profile schema v5 and reports
readiness per route. Its top-level `ready` field, compatibility field
`biomedical_joint_dp_capsule_ready`, and server alias
`formal_dp_claim_eligible` refer only to the custodian-attested surface plus
the live policy/pinset/authenticated control-plane handshake. They do not prove
that a particular dataset passes admission, that its manifest can be built,
that a route's numeric runtime is available, or that a live release will
complete. The route map says
`execution_readiness = "not_evaluated_requires_route_specific_preflight"` and
separately reports formal GLM and formal Cox as not ready for new computation;
it also reports whether a completed formal GLM or Cox public certificate can
be read.

The surface token is connector-neutral and is provisioned only by server
administrators after an exact effective-inventory check. Opal stores it in the
dedicated server profile; Armadillo/Rock may set it in the server/container
environment. No client option or DSI argument carries an attestation, and an
installed package never auto-attests. See
`inst/docs/remote_surface_attestation_20260808.md` for the admin procedure.

## Public cross-owner inventory (1.4.2)

Cross-owner support means the listed route can use variables held by different
custodians under its signed contract. It does not mean every legacy scope is
available across owners. In particular, current MI and survival Synopsis
producers remain same-owner. The unrestricted statement “every family works
with columns held by several custodians” is therefore not true for this release.
K is the number of custodians, not a count of variables or a utility guarantee.
The runtime registry (`ds.vertMethodStatus()`) lists each exported entry point.

| Public method family (aliases share the stated route) | Public cross-owner entry point in 1.4.2 | Estimand and limits | Remaining same-owner scope and reason |
|---|---|---|---|
| Identity and alignment | `ds.getIdentityPks`, `ds.psiAlign`, `ds.isPsiAligned` | Pinned private record alignment across K peers; no public intersection size | None; alignment is not a statistical release |
| Status, planning, calibration and certificates | `ds.vertSecurityStatus`, `ds.vertMethodStatus`, `ds.vertNumericPreflight`, `ds.vertDPStatus`, `ds.vertDPCapsulePlan`, `ds.vertDPCalibrate`, `ds.validateDPGaussianCertificate` | Public contracts, utility planning and authenticated artifact checks; readiness is route-specific | `ds.validateDPLMMCertificate` validates the existing same-owner LMM certificate; cross LMM uses its own grid validation |
| Count | `ds.vertDPCount` | One signed aligned-cohort count artifact; distinct artifacts compose | No multivariate estimand |
| Frequency and frequency inference | `ds.vertDPFrequency`, `ds.vertDPFrequencyInference` | One explicit source owner, two noise authorities, other peers sign; fixed categorical domain | One variable naturally has one owner; not a joint-table producer |
| Description, moments, quantiles and medians | `ds.vertDesc`, `ds.vertDPMeanVar`, `ds.vertDPDescribe`, `ds.vertDPQuantile`, `ds.vertDPMedian` | Bounded per-column summaries across an aligned federation; quantiles are release-only fixed-bin post-processing | Each described column is owned locally; no extra cross moments are inferred |
| Correlation | `ds.vertCor`, `ds.vertDPCor(artifact_family="gaussian_models")`, `ds.vert.cor` | All pairs from one authenticated intercept-containing Gaussian Synopsis; full registered-set complete cases, fixed anchor, bounded quantized Pearson correlation | Default `ds.vertDPCor` pairwise-complete artifact remains same-owner; no cross pairwise producer |
| PCA | `ds.vertPCA`, `ds.vert.pca` | PSD-projected Gaussian correlation eigenvalues, loadings and variance; offline current certificates accepted | Pairwise artifacts rejected; individual scores are never released |
| Gaussian regression | `ds.vertDPGaussian`, `ds.vertGLM(family="gaussian", dp_analysis_id=...)` | Bounded normalized complete-case Gaussian moments; explicit ridge if required | Existing same-owner artifacts retained; no sampling inference |
| Binomial / Poisson GLM | `ds.vertGLM(..., analysis_id=...)` | Authenticated finite signed coefficient grid; fixed bounded design | Completed formal-certificate reader retained; fresh formal/legacy iterative routes remain sealed |
| Gaussian mechanism inference | `ds.vertConfint`, `ds.vertContrast`, `ds.vertWald`, `ds.vertLR`, `type="mechanism"` | Gaussian fits with current same- or cross-owner certificates; simultaneous mechanism regions, including nested RSS reduction | No new sampling covariance, p-values or population-inference claim; LR is not a likelihood-ratio test |
| Gaussian LASSO | `ds.vertLASSO`, `ds.vertLASSO1Step`, `ds.vertLASSOProximal`, `ds.vertLASSOCV`; Gaussian `ds.vertLASSOIter` | Deterministic Gaussian L1 path/KKT and pseudo-information-criterion selection from validated current cross Gaussian provenance | Historical same-owner inputs retained; historic cross capsules rejected; CV name is not cross-validation |
| Binomial / Poisson LASSO | `dp_lasso_grid`; `ds.vertLASSOIter(..., cross_contract=..., base_contract=...)` | Signed finite candidate and lambda grid over an authenticated base contract | Existing same-owner `analysis_id` path retained; no unrestricted optimizer; Poisson live-connector evidence is not claimed by unit tests |
| Negative binomial | `dp_nb_grid`; safe `ds.vertNB`, `ds.vertNBFullRegTheta`, `ds.vert.nb` with `analysis_id` | Bounded NB2 beta/theta grid with complete negative log-PMF, signed by source owners | Intercept-only Frequency method-of-moments retained; historical disclosive NB body is internal and never dispatched |
| Multinomial | `dp_multinomial_grid`; `ds.vertMultinom`, `ds.vertMultinomJoint`, `ds.vertMultinomJointNewton` with `analysis_id` | Signed finite additive softmax grid with fixed outcome levels | Intercept-only Frequency path retained; no transformed predictors or unconstrained optimization |
| Ordinal | `dp_ordinal_grid`; `ds.vertOrdinal`, `ds.vert.ordinal` with `analysis_id` | Signed finite proportional-odds coefficient/threshold grid, exact declared clinical order | Intercept-only Frequency path retained; no unconstrained optimization |
| Cox PH | `dp_cox_grid` | Signed finite Breslow coefficient grid within admitted capacity; no standard errors or baseline hazard | `ds.vertCox`, profile and aliases retain same-owner grid/completed-certificate routes; their risk-set controls do not map transparently to the cross contract |
| Discrete hazard | `ds.vertCoxDiscreteNonDisclosive`, `ds.vert.cox(method="discrete")` read a compatible completed certificate | Fixed pooled-logistic time-grid certificate; no fresh person-period computation | No new cross-owner producer; this is not Cox PH |
| Linear mixed models | `dp_lmm_grid` | Signed finite cross-owner random-intercept ML grid | `ds.vertDPLMM`, `ds.vertLMM`, `ds.vertLMM.k3` and aliases retain same-owner moments, ML/REML and random-slope scopes; different contracts prevent implicit dispatch |
| Generalized linear mixed models | `dp_glmm_grid` | Signed binomial/Poisson random-intercept finite grid with fixed quadrature | `ds.vertGLMM` and aliases retain same-owner moments and random-slope grids; different estimands/controls prevent implicit dispatch |
| GEE | `dp_gee_grid` | Signed binomial/Poisson finite grid with fixed working correlation; evaluated independence lane; no sandwich covariance | `ds.vertGEE` retains Gaussian exchangeable/AR(1), robust-independence and completed-certificate scopes; those require different covariance contracts |
| Contingency, chi-square and Fisher | `ds.vertDPContingency`, `ds.vertChisq`, `ds.vertChisqCross`, `ds.vertFisher` | Signed cross-owner categorical pair; concordant-unit estimand, mechanism-aware/plugin inference | No classical exact-test claim for a noisy confidential table |
| DP epidemiology, diagnosis, ROC and standardization | `ds.vertDPEpi2x2*`, `ds.vertDPPrevalenceRatio*`, `ds.vertDPMantelHaenszel`, `ds.vertDPDiagnostic2x2*`, `ds.vertDPROC`, `ds.vertDPDirectStandardization*`, `ds.vertDPIndirectStandardization*`, `ds.vertDPCausalStandardization*` | Post-processing of a compatible authenticated cross table and fixed public mappings/weights | A stratum or combined arm must be one registered categorical column; no new cross multiway producer; scientific sampling/causal assumptions remain external |
| IPW | `ds.vertIPW`, `ds.vert.ipw` using a signed cross contingency table | Binary intercept-only or one categorical saturated-stratum IPW/g-formula identities; no released individual weights | Stratum+treatment must already form one registered arm column; arbitrary propensity models unavailable |
| Multiple imputation | No current cross-owner producer for the supported MI contract | `ds.vertMI`, `ds.vert.mi`: categorical MCAR marginal, joint-pair or conditional-star completion summaries | Same-owner only: cross categorical exclusion semantics do not meet MI's strict missingness contract; chained/continuous MI and Rubin inference unavailable |
| Survival curves and contrasts | No cross-owner survival Synopsis producer | `ds.vertDPSurvival`, Kaplan–Meier, Nelson–Aalen, cumulative incidence, RMST/RMTL, survival/RMST contrasts, log-rank and survival quantile/median post-processors | Time/event/entry inputs must share an owner; the existing producer rejects split ownership |
| Epidemiology from supplied aggregates | `ds.vertEpi2x2`, `ds.vertMantelHaenszel`, `ds.vertDirectStandardization`, `ds.vertIndirectStandardization` | Client-only calculations on already-authorized tables, regardless of upstream ownership | Arbitrary matrices carry no DataSHIELD provenance; these helpers do not acquire cross data or establish a DP claim |

Correlation and PCA reuse the existing Gaussian DP artifact without another
release. See the [worked two- and three-owner custodian configurations](inst/docs/cross_owner_correlation_1.4.2.md)
and the [implemented design](DESIGN_CROSS_OWNER_CORRELATION.md). A server may be
1.4.1 if all required protocol capabilities validate; the coordinated packages
are 1.4.2. The provenance capability remains version 1.4.1 and runtime bytes
remain unchanged. Missing/incompatible capabilities fail closed.

## Validation evidence

For 1.4.2, [targeted validation results](inst/validation/v1.4.2/RESULTS.md)
record 8,671 client and 6,382 server assertions passing, one optional Cox harness
skip, unchanged same-owner golden/artifact/certificate/release bytes, and the
separate-process DSLite smoke. K2 Cor/PCA, exact replay, offline conversion with
zero DSI calls, binomial and Poisson LASSO, and K3 Cor/PCA over three owners
(n=64, four variables including a same-owner pair) completed. The K3 Cor release
took about 56 minutes, dominated by serial exact-GC multiplication stages; K5 was
not attempted. Two client lifecycle and certificate-validation defects present in
1.4.1 were repaired without changing server guards or release identities.
Archived promotion evidence below is historical and does not replace these checks.


The repository contains unit, adversarial, centralized-reference, real-DSLite,
connector-lifecycle and live-Opal harnesses. PSI additionally has a two-host
test with independent identity seeds. `K=2`, `K=3` and `K=5` statements for
the promoted Count, PSI and Synopsis routes,
identify topology coverage in the named test;
they do not imply that every topology completed on every connector. The local
Armadillo S4/httpuv smoke is a connector-lifecycle regression, not a live Rock
deployment, and cached Opal results are historical evidence rather than a
blanket release guarantee. Current deployment claims require a fresh,
artifact-exact live run with the method-specific tolerance and assumptions
recorded by the harness. Analytical calls use DSI; Opal reconciliation and
Armadillo TLS/session inspection remain deployment safeguards rather than
statistical-protocol dependencies.

Ordinary in-process DSLite cannot faithfully model separate persistent server
identities because logical servers share one R process. Identity-sensitive
security tests therefore use isolated host environments; live Opal validation
remains required before deployment.

## Security

The pinned-peer model prevents the analyst/relay from substituting a key and
opening peer ciphertext. PSI membership messages are signed, encrypted and
session-bound; completed alignment is verified through a secret-token-backed,
ordered-ID manifest rather than row counts. Session IDs and persistent identity
seeds use operating-system entropy, and unused high-risk remote primitives have
been removed from registration and export.

These controls do **not** imply that source data can never be reconstructed.
Unlimited adaptive exact outputs can be composed without breaking MPC.
`ds.vertDPCount()` is independent of the capsule allocator and accountant. It
compiles one canonical signed current-snapshot artifact; the same artifact
recomputes the same identity-seeded noise, while distinct artifacts compose
and have no finite global composition claim. Fixed-cohort Count instead uses
an exact public K-consensus declaration with zero sensitivity and no MPC noise.
`ds.vertDPFrequency()` likewise has no capsule allocator: one canonical signed
source/factor artifact derives one identity-seeded two-authority release;
replays are identical and distinct signed artifacts compose without a claimed
finite global budget.
The still-provisional table, moment, Gaussian-model and survival routes retain
their separate cross-signed biomedical-capsule contract while they are being
migrated. DP bounds indistinguishability; it is not a literal guarantee that no
original value can ever be guessed or reconstructed. The claim does not cover
collusion of every designated peer. Ordinary exact methods
remain outside any DP guarantee. Historical local fixed-point truncation and
legacy DCF comparison survive only inside unregistered/quarantined
implementations and retain their stated proof limitations. See
[SECURITY.md](SECURITY.md) for the full threat model and non-claims.

There is no request-count or lifetime admission limit on active public Synopsis
routes. Shape, byte, finite-range and storage caps remain separate resource and
arithmetic-safety controls, never privacy budgets or expiring replay rights.

## DSI communication

Persistent DSI connections need no generic application/request batching layer
and no per-peer wrapper serialization. The capsule-source route can carry a
byte-bounded consecutive window of unchanged signed inner artifacts in one
fetch/accept request after public capability negotiation. Its outer framing is
not independently signed and it never creates multiple in-flight DSI jobs.
Large opaque values still require frame
chunking because DSI connectors and parsers impose finite request sizes. One
named fan-out cycle carries the next dependency-ready immutable frame for every
peer, and the source-to-recipient pipeline completes `F` frames in `F + 2`
healthy phases, without pretending that dependent MPC rounds are parallel. The
frame path uses DSI's public job
primitives directly and avoids repeated deparsing of the same large expression
during polling. Stateless request-size probing is session-bound; a reconnect
revalidates the last connector-profile size before using it. The registered
typed transport uses absolute offsets, bounded private node spools, complete
SHA-256 verification and pinned-peer signed receipts. The analyst keeps only
one frame and no payload-sized integrity copy. The more general relay remains
unregistered until independent multi-host live-connector validation and
protocol migration prove that exposing it does not widen the remote surface.

Repeated analyses first refresh the complete pinned control-plane status and
rebuild the server-authoritative capsule manifest. An in-process LRU then
reuses an already validated public DP vector only when the current
`manifest_sha256`, `capsule_id`, policy, pinset, noise root and release domain
are identical. A hit skips the expensive source-sharing, sampling and vector
replay phases; it never skips snapshot/workload validation. The cache retains
no connection handles, is bounded to four entries and 64 MiB, and eviction
cannot deny an operation or generate different noise—the durable sticky replay
remains authoritative.

## Installation

```r
# From GitHub
devtools::install_github("isglobal-brge/dsVertClient")

# Or from source
R CMD build --no-build-vignettes .
R CMD INSTALL "$(ls -t dsVertClient_*.tar.gz | head -1)"
```

Requires the server-side package [dsVert](https://github.com/isglobal-brge/dsVert) installed on all DataSHIELD servers.

## Documentation

- [Reference index](https://isglobal-brge.github.io/dsVertClient/) (pkgdown site)
- [DP statistical validation evidence](https://github.com/isglobal-brge/dsVertClient/blob/main/inst/docs/benchmarks/dp_statistical_validation_20260802.md)
- [Historical capsule v6 large-scale audit; current workload is v7](https://github.com/isglobal-brge/dsVertClient/blob/main/inst/docs/capsule_v6_large_scale_audit_20260808.md)
- [Product surface audit](https://github.com/isglobal-brge/dsVertClient/blob/main/inst/docs/product_surface.md)
- [NEWS](https://isglobal-brge.github.io/dsVertClient/news/index.html)

## License

MIT - see [LICENSE](LICENSE.md).
