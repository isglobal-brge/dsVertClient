# dsVertClient - DataSHIELD Client for Vertically Partitioned Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)](https://isglobal-brge.github.io/dsVertClient/news/index.html)

## Overview

**dsVertClient** provides R workflows for privacy-aware analysis of vertically
partitioned data across DataSHIELD servers. It combines ECDH-PSI, Ring63 /
Ring127 additive sharing, dealer-free IKNP OT-extension preprocessing,
peer-to-peer authenticated encryption and Ed25519 server identity pinning.

Method maturity is intentionally explicit. Run `ds.vertMethodStatus()` before
choosing a route: only a small conditional surface is currently `promoted`;
the main distributed estimators remain `provisional`, `compatibility`, or
`quarantine` while their numerical, disclosure and biomedical contracts are
closed.

All connected servers must expose dsVert's single immutable
`disclosure_safe` profile. There is no client or server switch that enables an
exact legacy release plane; methods still awaiting migration return an
explicit unavailable/migration error instead of falling back silently.

Pair with the server-side companion package [**dsVert**](https://github.com/isglobal-brge/dsVert).

## Quick Start

```r
library(DSI); library(DSOpal); library(dsVertClient)

# Connect to servers
conns <- DSI::datashield.login(logins = builder$build())
DSI::datashield.assign.table(conns, "D", "project.data")

# 1. Align records without mutating the source table
ds.psiAlign("D", "patient_id", "DA", na.action = "na.omit",
            datasources = conns)

# 2. Inspect the immutable server policy and signed, data-free capsule plan
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

The package can inspect—never change—the immutable DP policy and run a public,
data-free utility preview. Count independently compiles one canonical signed
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
that decision: the signed server manifest is authoritative.

```r
ds.vertDPCalibrate(capsule_epsilon = c(1, 3, 5),
                   peer_count = length(conns),
                   sensitivity = 1)
ds.vertDPStatus(conns)
ds.vertDPCapsulePlan(conns) # signed dry-run; no data access or release
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
`ds.vertDPStatus()` and `ds.vertDPCapsulePlan()` are data-free signed Synopsis
bootstraps: they bind the pinset, two designated authorities and immutable
workload, but create neither a release nor a protected-source read. The retired
vector-lifecycle ABI remains internal compatibility code and is not a public
release route.

`ds.vertSecurityStatus()` consumes security-profile schema v4 and reports
readiness per route. Its top-level `ready` field, compatibility field
`biomedical_joint_dp_capsule_ready`, and server alias
`formal_dp_claim_eligible` refer only to the custodian-attested surface plus
the live policy/pinset/authenticated control-plane handshake. They do not prove
that a particular dataset passes admission, that its manifest can be built,
that a route's numeric runtime is available, or that a live release will
complete. The route map says
`execution_readiness = "not_evaluated_requires_route_specific_preflight"` and
separately reports formal GLM and formal Cox as not ready.

The surface token is connector-neutral and is provisioned only by server
administrators after an exact effective-inventory check. Opal stores it in the
dedicated server profile; Armadillo/Rock may set it in the server/container
environment. No client option or DSI argument carries an attestation, and an
installed package never auto-attests. See
`inst/docs/remote_surface_attestation_20260808.md` for the admin procedure.

## Function surface (v1.1.0 development audit)

This is a capability index, not a claim that every route has the same maturity.
Use the runtime registry for the authoritative status and limitation of each
entry point.

| Family | Functions |
|---|---|
| **Record alignment** | `ds.psiAlign()`, `ds.isPsiAligned()`, `ds.getIdentityPks()` |
| **Descriptive / 2nd-order** | `ds.vertDesc()`, `ds.vertCor()`, `ds.vertPCA()`, `ds.vertChisq()`, `ds.vertChisqCross()`, `ds.vertFisher()` |
| **Differential privacy** | `ds.vertDPStatus()`, `ds.vertDPCapsulePlan()`, `ds.vertDPCalibrate()`, `ds.vertDPCount()`, `ds.vertDPContingency()`, `ds.vertDPFrequency()`, `ds.vertDPFrequencyInference()`, `ds.vertDPMeanVar()`, `ds.vertDPDescribe()`, `ds.vertDPQuantile()`, `ds.vertDPMedian()`, `ds.vertDPGaussian()`, `ds.vertDPSurvival()`, `ds.vertDPSurvivalQuantile()`, `ds.vertDPMedianSurvival()`, `ds.vertDPSurvivalContrast()`, `ds.vertDPRMSTContrast()`, `ds.vertDPEpi2x2()`, `ds.vertDPEpi2x2Inference()`, `ds.vertDPMantelHaenszel()`, `ds.vertDPDiagnostic2x2()`, `ds.vertDPDiagnostic2x2Inference()`, `ds.vertDPROC()`, `ds.vertDPDirectStandardizationInference()`, `ds.vertDPIndirectStandardizationInference()`, `ds.vertDPCausalStandardization()`, `ds.vertDPCausalStandardizationInference()` and the remaining direct/indirect standardisation, competing-risk, RMST, and RMTL post-processors |
| **GLM** (gaussian / binomial / poisson) | `ds.vertGLM()` / `ds.vert.glm()` are promoted only for a compatible same- or cross-owner Gaussian request with a custodian-configured `dp_analysis_id`; it delegates to the signed sticky `ds.vertDPGaussian()` artifact. Binomial/Poisson and other iterative routes stay gated until their formal capsule backends are production-ready. |
| **Inference helpers** | `ds.vertConfint()`, `ds.vertContrast()`, `ds.vertWald()`, `ds.vertLR()` |
| **Epidemiology from authorised aggregates** | `ds.vertEpi2x2()`, `ds.vertMantelHaenszel()`, `ds.vertDirectStandardization()`, `ds.vertIndirectStandardization()` |
| **Survival** | Sticky fixed-grid DP survival artifacts and their Kaplan--Meier, Nelson--Aalen, competing-risk, survival-quantile, RMST and RMTL post-processors are available. Cox frontdoors are retained but quarantined and fail before DSI until the formal Cox capsule is production-ready. |
| **Negative binomial** | Historical names are retained but quarantined; no NB2 estimator is advertised as disclosure-safe until a complete bounded score/information capsule and validated joint inference are available. |
| **Multinomial** | `ds.vertMultinom()` / `ds.vertMultinomJointNewton()` are quarantined; slope models fail closed before DSI until a signed raw-design Gram workload exists |
| **Ordinal (proportional odds)** | Historical joint-Newton names are retained but quarantined and fail before DSI pending a formal bounded proportional-odds capsule. |
| **Mixed / clustered models** | Historical LMM, GEE and GLMM routes; every public frontdoor is currently quarantined and fails before DSI as reported by `ds.vertMethodStatus()` |
| **Causal / robustness** | `ds.vertIPW()` and `ds.vertMI()` are quarantined compatibility names; their retained implementations are respectively a known-weight workflow and server-local MI, not promoted end-to-end methods |
| **Penalised regression** | `ds.vertLASSOIter()` is quarantined; `ds.vertLASSO()` and `ds.vertLASSO1Step()` are compatibility helpers; the same-owner Gaussian Synopsis post-processors `ds.vertLASSOProximal()` and the information-criterion (not cross-validation) selector `ds.vertLASSOCV()` are promoted without sampling inference |

## K=2 vs K≥3 support

| Family | K=2 product route | K≥3 product route | Legacy / not offered |
|---|---|---|---|
| GLM (gauss / binom / poisson) | Compatible same- or cross-owner Gaussian signed-Synopsis adapter only; binomial/Poisson gated | same, tested through K=5 | legacy iterative routes are not production claims |
| Cox PH | Formal capsule under internal validation; public frontdoors gated | same | `legacy_rank`, `ds.vertCox.k3()` removed |
| Negative binomial | Quarantined pending bounded NB2 capsule and joint inference | same | disclosive `variant = "full_reg"` removed |
| Multinomial | Slope route unavailable pending signed `multinomial_design_grams` | same | unsafe local-moment/correlation Gram reconstruction and warm / OVR final-estimator route removed |
| Ordinal | Quarantined pending a formal proportional-odds capsule | same | warm final-estimator and patient-level joint reconstruction routes removed from the exported API |
| LMM | Historical K=2 route, quarantined | Historical K>=3 profile, quarantined; `ds.vertLMM.k3()` deprecated | confirmatory inference not supported |
| GEE | All working-correlation modes are quarantined and fail before DSI | same zero-DSI boundary | robust cluster meat must remain entirely in MPC before promotion |
| GLMM | Experimental PQL/Laplace routes, quarantined | same maturity boundary | valid final covariance/SE and safer cluster aggregation are incomplete |
| IPW / MI / LASSO / Cor / PCA / Chisq / Desc | mixed maturity; query `ds.vertMethodStatus()` | same | exact legacy outputs remain subject to composition risk; signed capsule routes retain their explicit formal-DP contract |

See `inst/docs/product_surface.md` for the audited maturity surface.

## Validation evidence

The repository contains unit, adversarial, centralized-reference, real-DSLite,
connector-lifecycle and live-Opal harnesses. PSI additionally has a two-host
test with independent identity seeds. `K=2`, `K=3` and `K=5` statements for
the promoted Count and PSI routes, and for provisional capsule routes,
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
