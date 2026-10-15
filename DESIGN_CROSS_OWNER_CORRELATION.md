# Cross-owner correlation and PCA — implemented 1.4.2 design

**Implemented after Claude consensus, 2026-10-15. Validation limits are recorded below; no push or publication.**

Baseline: dsVert `2dba3996764454268b03b82bb2eeb9c220e982bc`, dsVertClient `04a690341915ac59e2753a331dfb1826f8a01970`, both 1.4.1, branch `feature/cross-owner-correlation`. The phase-A baseline is retained below. Phase B changes only these two repositories; no pushes, issue operations, full R suite or full Go sweep.

## Decision

Reuse the existing signed cross-owner Gaussian sufficient-statistic Synopsis, with an intercept, for all-pair complete-case correlation and PCA. Preserve the existing same-owner pairwise artifact and Gaussian artifact unchanged. Complete and verify the client paths already present, repair cross-only verification inconsistencies, add explicit authenticated zero-call conversion, and publish the four promoted grid adapters. Do not introduce a dedicated correlation producer, new noise sampler, new DataSHIELD endpoint or new privacy budget.

This covers numeric columns split across K custodians, including several columns at one custodian and all predictor–predictor as well as predictor–outcome pairs. K=2,3,5 is topology coverage to establish in the targeted tests, not a claim of every connector/deployment having been tested. An existing scientific outcome may retain its role. Without one, custodians select the smallest owner-qualified reference in radix order as the fixed public storage outcome/anchor; no dummy variable or extra private outcome is required.

## Findings that correct the initial premise

At the 1.4.1 baseline, `ds.vertDPCor` was same-owner pairwise-only, but `ds.vertCor` has no runtime same-owner gate. Its implementation calls the Gaussian Synopsis reader, validates requested owners against the signed artifact and handles cross participants. Its own pairwise rejection message already directs callers to `ds.vertCor` for cross-owner use. Same-owner-only roxygen/registry statements are stale. See [R/ds.vertCor.R](R/ds.vertCor.R), especially the pairwise gate around 430, all-coordinate unpacking around 630, Gaussian constructor around 816, and stale adapter docs around 1008.

The test called “cross-owner complete-case Cor is quarantined for K2 and K3” forces the rejection in its retrieval mock (`test-dp-cor-pca.R`, helper around 494). It does not show that the production reader rejects cross artifacts. Existing public Gaussian tests also mock retrieval; they prove consumer behavior, not the complete current Synopsis lifecycle.

The baseline had a concrete defect (now fixed): correlation construction supplies the signed cross quantization allowance to `.dsvert_dp_cor_interval`, but `.dsvert_dp_cor_gaussian_certificate_match` reconstructs the interval using the same-owner default. PCA compares these values exactly and can reject a valid cross correlation. A bounded 80-unit, capacity-100, scale-256 diagnostic produces cross interval `[0.6818035,1]` versus verifier interval `[0.9052327,1]`. Broad `[-1,1]` test intervals can conceal the discrepancy. The implementation uses the same signed cross allowance in both paths, preserving the v1 arithmetic/order.

A second baseline defect (now fixed) was demonstrated for owner-qualified duplicate physical names. A current server-built artifact with outcome `peer_a$x` and predictor `peer_b$x` passes server construction but the client rejects it with “The signed cross-owner Gaussian predictor order changed”; the unique-name control succeeds. Parts of the Gaussian client compare canonical references directly to physical column names. The read-only diagnostic uses existing fixture signature/selector seams and is not federation evidence. Fix these cross-path inconsistencies; never drop owner/dataset binding or merge two owners' identically named columns. Distinguish canonical signed variable identifiers from physical column names throughout the new converter. Preserve unqualified same-owner outputs. Use the existing `owner$column` grammar and signed dataset binding; ambiguity when one owner has the same physical column in multiple datasets must fail closed unless separately designed.

The vertcor03 separate-process diagnostic found another production client defect
present in 1.4.1: the dedicated Gaussian computation session created a fresh UUID
without installing its Synopsis authorization before exact-MPC transport binding.
The actual endpoint error was `dsvert.dp.ledger_path must be one non-empty string`,
because the unauthorised session selected the legacy DP policy. The previously
recorded `dsvert_phase_not_ready` came from an earlier publication lookup, not
this binding failure. The client now sends the existing authenticated
`dsvertDPSynopsisPrepareDS` request to both designated authorities before creating
transport keys, using the same child UUID and manifest/claim/compilation evidence.
This matches the established cross-grid lifecycle. PREPARE validates the existing
signed authorization and installs artifact-derived execution state; it neither
starts a noise draw nor creates a new release identity. Server pin, signature,
attestation, provenance and fail-closed checks remain unchanged. Focused tests
cover ordering, exact context/session binding, rejected or malformed authorization
before transport, and the unchanged legacy route. Current federation outcomes
are recorded in [RESULTS.md](inst/validation/v1.4.2/RESULTS.md); this finding does
not establish that every historical cross-owner Gaussian deployment failed.

## Why the 26 August admission was reverted

There is **no documented author rationale** in `b83260a`, Git notes, or the searched tracked notes/NEWS. The revert says only that it reverts `766a663`; it followed 88 seconds later. Do not present a guessed motive as a historical fact.

There is direct evidence that the proposed promotion was premature:

- `766a663` changed documentation, method status and mocked tests, not executable admission/transport code.
- At that revision `dp_synopsis_vector_runner.R:343–347` admitted only categorical cross catalogs; Gaussian cross catalogs failed before execution.
- Its Gaussian certificate builder (`dp_gaussian_certificate.R:969–972`) still required closed same-owner Synopsis provenance. The changed tests bypassed this by mocking the release function.
- `b0a7c6e` later supplied the missing cross Gaussian Synopsis routing, purpose-bound Bind/Finalize, durable evidence/replay and certificate plumbing.
- `8e7a929` fixed the evidence shape (direct designated-peer map, not an analysis-id-nested map) and made correlation certificate matching expect the correct cross state. It did not add a same-owner gate.

Thus missing production transport/provenance is a code-supported explanation of why admission was not ready, not proof of the author's private motive. There is no evidence here of a failed DP proof or missing all-pairs coordinates. The new acceptance tests must exercise current signed Synopsis transport/provenance, not merely resurrect the old mocked promotion test.

## Existing artifact: exact release and estimand

Let m be the number of actual numeric variables, one outcome and m−1 predictors. With an intercept, the design has q=m terms. The descriptor is `bounded-normalized-gaussian-cross-sufficient-statistics-v1`, `spec_version="v2"`.

Each custodian signs finite bounds `[l_j,u_j]`. For each admitted patient, finite records are clipped **before** averaging once per variable; no finite value means invalid for that variable. The mean is normalized to `[0,1]`, then `a_j=round(S*z_j)` with ties-to-even and `S=2^b`, b from 8 through 18. Values and validity travel privately as padded encrypted Ring128 inputs to the two designated computation peers. Exact protected products form one joint validity mask across the full signed set, masked values and moments; products use `floor(a_i*a_j/S)`. Capacity/record-count overflow rejects the snapshot rather than truncating it.

The Gaussian block releases, in order: n, upper-triangular `X'X` column-major, `X'y` in design order, and `y'y`, with the existing joint noise and public clamps. With an intercept it contains every marginal sum and every second moment. For predictors x,z and outcome y, the natural-scale pre-noise order is:

```
n, n, sum(x), sum(x²), sum(z), sum(xz), sum(z²),
sum(y), sum(xy), sum(zy), sum(y²)
```

The separate n and intercept Gram entry are duplicate noisy masses with distinct existing coordinates; never average them. Gaussian block dimension is `d=q(q+1)/2+q+2=(m²+3m+4)/2`, namely 7,11,22 for m=2,3,5 (m is variable count, not custodian count K).

The signed cross quantization certificate bounds each summed monomial's absolute error by `capacity*(2+1/(4*S))` lattice steps. Cross v2 is not numerically identical to the same-owner v1 quantizer and does not claim to be. Ring safety requires `ceil(log2(capacity+1))+2*b+3 < 127` with the existing bounded-input/product checks.

The catalog additionally publishes the mandatory admitted-count coordinate and, for these Gaussian variables, three numeric marginal coordinates per variable (finite count, sum, sumsq). Those marginals can concern different valid-patient cohorts; do not use them to fill missing complete-case sums. The complete signed catalog, including any other configured families, determines sensitivity and noise.

Without an intercept, the Gaussian block omits the marginal sums: general centered Pearson correlation is not identifiable. Preserve `pearson_requires_intercept_marginal_sums`. Never substitute cosine similarity or unrelated marginals. A subset of a signed larger set retains the larger set's joint complete-case cohort; variable subsetting cannot silently change missingness.

Evidence: server [dpGaussianCrossDS.R](../dsVert/R/dpGaussianCrossDS.R), [dpCapsuleWorkload.R](../dsVert/R/dpCapsuleWorkload.R) around 2854–2915 and 5168–5368, and [dpPolicyDS.R](../dsVert/R/dpPolicyDS.R) around 2493 and 2626; client [ds.vertDPGaussian.R](R/ds.vertDPGaussian.R) around 297–620 and [ds.vertCor.R](R/ds.vertCor.R) around 630–679.

## Correlation, utility and PCA

After noisy release, the existing client clamps and PSD-projects the augmented moment matrix. The projection occurs on the client, not as an extra exact server release. Extract projected intercept mass h, marginal sums s and second moments Q. Compute `C=Q-s*s'/h` and `R_ij=C_ij/sqrt(C_ii*C_jj)`, preserving current symmetrization, eigenvalue clipping, diagonal rescaling and diagnostics. The projected intercept mass governs centering. Nonpositive/nonfinite mass or variances retain typed non-identifiability. A singular covariance matrix with positive individual variances remains valid for correlation/PCA; do not require successful unregularized Gaussian regression to obtain it.

Positive affine normalization preserves Pearson correlations of the clipped, patient-aggregated variables. All pairs use the same joint cohort and one existing vector. PCA returns eigenvalues, loadings and explained/cumulative variance, with existing mechanism-region/eigengap diagnostics. It never releases scores, source rows, fitted values or residuals. `scores=NULL` and `scores_available=FALSE` remain mandatory. Loading signs and repeated-eigenvalue subspaces are not identifiable as unique vectors; tests compare appropriate invariant quantities where needed.

Gaussian-block natural L1/L2 sensitivity is d/sqrt(d), respectively, under the existing add/remove-patient or replacement adjacency: each patient's vector lies in the nonnegative unit cube and zero/all-one complete-case endpoints are possible. This is not the deployed full-catalog sensitivity. For a minimal existing Gaussian catalog, `D=1+3*m+d`; add/remove natural L1 is D and L2 is sqrt(D). For fixed-cohort replacement the mandatory count is constant, giving D−1 and sqrt(D−1). Additional catalog entries retain their signed contributions.

For the discrete-Laplace convolution, each authority supplies an independent full-epsilon draw calibrated to the entire catalog. Ignoring modular wrap for this utility calculation, geometric decay parameter `r=exp(-epsilon/Delta1_steps)` gives normalized coordinate variance `4*r/((1-r)^2*S^2)`, approximately `4*(Delta1_natural/epsilon)^2`. Use the actual signed mechanism/accuracy plan; discrete Gaussian uses full-catalog L2 calibration and two full calibrated contributions. K does not divide epsilon or supply K noisy releases. All pair numerators inherit this common moment scale, but there is no constant additive noise scale for correlation itself: denominators, DP mass, clipping and PSD projection make error nonlinear. Near-zero variance can yield non-identifiability or vacuous `[-1,1]` regions. Retain simultaneous noise-plus-quantization regions; they are not sampling confidence intervals.

A dedicated symmetric artifact could release `1+m+m(m+1)/2=d−1` moments and potentially avoid the 3m auxiliary marginal coordinates. It could improve a newly registered workload's utility, but needs new producer, sensitivity/calibration, schema, signing, identity and lifecycle validation. Merely removing one coordinate from an otherwise identical catalog improves L1 scale only by `(D−1)/D`. Reusing an already released artifact incurs no further privacy cost and preserves its utility and bytes. This is the lower-risk 1.4.2 choice; no dedicated-artifact utility improvement is claimed.

## Analyst APIs and zero-call behavior

The following additive arguments are implemented and appended to preserve positional compatibility.

| API | Existing behavior retained | New supported behavior |
| --- | --- | --- |
| `ds.vertCor` | Current signed Gaussian complete-case reader and result structure | Verified cross owner map; optional `gaussian_result=NULL` accepting current Synopsis Gaussian result **or certificate**, plus `trusted_pinset=NULL`, for zero-call reconstruction |
| `ds.vertDPCor` | Default same-owner pairwise artifact, pairwise missingness, counts, outputs and refusal rules | Explicit `artifact_family="gaussian_models"` delegates to the complete-case reader; default stays `"correlation_artifacts"`. Optional authenticated input as above; inconsistent selector/input rejected |
| `ds.vertPCA` | Dataset route invokes Cor once; `cor_result` route performs no DSI calls and rejects pairwise/arbitrary matrices | Correct cross certificate reconstruction; append `trusted_pinset=NULL` for authenticated offline reload; optional direct Gaussian-input forwarding if kept unambiguous |
| `ds.vert.cor`, `ds.vert.pca` | Existing aliases/positional calls | Forward the additive arguments without changing result fields; bypass the current ds.vert.cor datasource resolution on the offline-input branch |

`ds.vertCor` may make `data_name=NULL` only for the new input-object branch, deriving dataset and analysis id from the certificate; supplied assertions must match exactly. All existing calls continue to require the same identifiers. Do not fall back between artifact families after a validation failure, concatenate different released pairwise matrices, infer a new private cohort or call retired exact kernels. An explicit complete-case selector makes the change in missingness visible; owner count alone must not silently change the estimand.

Factor only the existing verified-coordinate algebra into a shared internal conversion path. Online acquisition keeps its current certificate builder and verifier. Offline acquisition must validate the **current Synopsis** certificate, all roots/coordinates and authenticated pins, then use reconstructed moments. It must not trust a supplied correlation matrix, fit coefficients, metadata booleans, cached verification flag or self-consistent hash. Reject legacy capsule provenance through this new convenience input unless separately justified; preserve existing legacy APIs elsewhere.

An offline certificate can be integrity-valid but `unanchored`. Require `caller_anchored` via an external trusted pinset, or `session_transport_anchored` via the existing online trust cache. Forward trusted pins into validation; do not mutate the cache or label self-supplied certificate pins as trusted. A certificate is a valid input without first fitting a regression, avoiding singular-fit loss of valid PCA functionality.

“Zero calls” means conversion of an already available authenticated artifact. A dataset/id invocation may acquire/replay the existing Synopsis workflow. Preserve this distinction in fields/docs: zero additional calls after acquisition and zero additional epsilon/delta; do not promise an online cache read has no transport. New input branches need an explicit zero-discovery/zero-aggregate test.

## Custodian registration and signing

1. Every custodian configures its local patient column, dataset identity/version, numeric bounds, schema mapping, admission/capacity/adjacency and persistent peer identity. All agree the full pinset, designated two authorities, public epsilon/delta and closed catalog. K owners can feed the same two computation/noise authorities.
2. Use authorized frozen 1.4.1 sources with compatible private PSI alignment. A public equality of row counts is insufficient. Imported alignments lacking complete history follow the existing explicit authorization migration.
3. Exactly the outcome/anchor custodian registers one global `dsvert.dp.gaussian_specs` v2 fragment, using the canonical signed schema labels (bare when globally unique, qualified for collisions). For example, with globally unique x and z:

   ```r
   cor_primary = list(version="v2", dataset="anchor_data", outcome="anchor",
                      predictors=c("x", "z"), intercept=TRUE)
   ```

   This is a configuration fragment, not an analyst call. Other custodians register their owned columns/bounds and approve the combined schema; they do not duplicate the analysis-id fragment. Owner-qualified references must resolve to exactly the signed dataset/owner/column.
4. With a scientific outcome, keep that role; it may be included in the correlation set or omitted while retaining its complete-case restriction. Without one, choose a public canonical anchor from the requested set (e.g. smallest owner-qualified reference) during registration and freeze it. Relabelling an existing artifact's outcome/predictors is not free replay.
5. All peers sign the global deterministic schema/workload and the existing claims/compile receipts. Both designated authorities sign the protected computation/release evidence. Analyst input only selects an approved analysis and subset/order; it never supplies bounds, privacy budget or unsigned contracts that replace custodian policy.

See [dpCapsuleManifestDS.R](../dsVert/R/dpCapsuleManifestDS.R), particularly the outcome-custodian ownership rule around 496, and [dpSynopsisLifecycleDS.R](../dsVert/R/dpSynopsisLifecycleDS.R) around 1300–1577. Add a short two-owner correlation-only worked configuration and a three-owner mixed-pair example to current docs, labelling the wire role honestly.

## Privacy and replay argument

**Release boundary.** Existing protected products feed only the two-authority joint DP finalizer. The analyst receives the fixed bounded noised Synopsis vector and its public authenticated provenance, never private validity, masks, row IDs, exact moments, intermediate shares or PCA scores. New outputs are deterministic functions of that vector and public contracts, so post-processing adds `(0,0)` privacy cost.

**Mechanism and threat model.** Retain the pinned, protocol-compliant/semi-honest threat model allowing an analyst to collude with at most one of the two designated authorities, leaving at least one hidden authority contribution, authenticated source publication and intact private state. Each authority contributes a complete calibrated draw; conditioning on the adversary's permitted view leaves one complete hidden DP mechanism, and independent additional noise is post-processing. The current exact discrete-Laplace route has pure epsilon-DP under ideal bits (production uses keyed computational randomness); Gaussian/finite-sampler routes keep their declared `(epsilon,delta)` including implementation terms. Ring reduction, fixed clamping, PSD projection, normalization and eigendecomposition are post-processing. Do not claim malicious-peer security, whole-transcript DP or finite lifetime composition.

**Replay and validation.** Reuse exactly the signed catalog/source/layout/mechanism identity and current provenance epoch closure; no new private-content-derived identity, rerandomization or budget allocation. Identical authorized artifacts replay their existing release; distinct publications/queries compose. Certificate validation remains mandatory before all algebra, including offline use.

In particular preserve the 1.4.1 `dsvert-provenance-epochs-v1` check on every participant and E_star, the complete effective dataset dependency closure including membership-only PSI participants. Frozen authorized material is served when a live frame changes in place; a new explicit publication allocates a new random epoch even for identical bytes, and event retries reuse their epoch. Reauthorization changes affected artifact/authority streams; stable equality must not become a private-content-change oracle. Keep private seals, persistent UUID/state integrity, and fail-closed recovery. No fallback to 1.4.0 content-derived identities. These are conditional mechanism guarantees, not protection against private publication schedules or whole-state rollback. See [server PROVENANCE_EPOCHS.md](../dsVert/PROVENANCE_EPOCHS.md).

Retain the `8e7a929` direct-peer evidence shape. Validate exactly both designated pinned authority attestations and agreement on analysis/artifact/source/layout/block span/order, transcript and numeric-certificate hashes, ring/fraction bits and terminal no-exposure state. Revalidate manifest/schema signatures, compiler/release/replay proofs, coordinate roots/bytes and externally anchored identity. Recomputing an outer checksum must not repair a tampered inner field. See [dp_gaussian_cross_transport.R](R/dp_gaussian_cross_transport.R) around 576–644 and [dp_gaussian_certificate.R](R/dp_gaussian_certificate.R) around 1425–1509.

## No functionality loss

Do not rewrite or version existing same-owner artifacts, layouts, sensitivity formulas, masks, quantizers, noise inputs, replay identities, certificates or result structures. Default `ds.vertDPCor` remains pairwise-complete, with PCA ineligible. Same-owner Gaussian Cor/PCA remains joint complete-case. Preserve all existing aliases, fitting routes and documented refusal cases. Fix the quantization discrepancy in the cross branch only; a shared helper must return the exact v1 default and operation ordering.

Capture deterministic baseline serialization for same-owner Cor/PCA/aliases and compare byte-for-byte after implementation; additionally compare raw signed artifact/certificate/release bytes under retained fixed fixtures. JSON or numeric approximate equality is insufficient for this requirement. Existing tests that deliberately assert stale public inventory/quarantine labels need an explicit metadata-only expectation update; preserve their substantive security and same-owner assertions. “All existing tests green” cannot literally mean retaining assertions that four requested public exports must stay absent.

Do not rebuild/re-stamp an unchanged Go runtime merely to change the R package patch version. Runtime/source/key/protocol domain changes can change authenticated artifacts; preserve shipped runtime bytes/checksums and stable identity inputs where unchanged. If a required qualified-name fix changes a previously invalid cross descriptor, scope it to that case and do not redefine any already valid artifact.

## Four promoted grid families and LASSO front doors

Publish the existing authenticated `dp_nb_grid`, `dp_lasso_grid`, `dp_multinomial_grid`, `dp_ordinal_grid`: consistent roxygen, NAMESPACE, Rd, registry rows, examples and release flags. No unsafe historical optimizer becomes public. The promotion manifest records nine cells per simple family (K=2/3/5, epsilon=1/4/8, n=2000,p=6,two candidates) with the stated oracle/lifecycle gates. Retain its limits: archived LASSO cells are binomial; Poisson has shared code and unit coverage, not equivalent archived federation evidence. The old per-cell oracle-output digest gap remains disclosed. No mathematical blocker was found for publishing the authenticated adapters; a new focused Poisson smoke is required before claiming fresh Poisson federation evidence.

| Family | Public dispatch and preserved limitations |
| --- | --- |
| NB | Add cross dispatch to supported `ds.vertNBFullRegTheta` and `ds.vert.nb`. `ds.vertNB` is currently **not exported** and its historical body is unsafe to expose; implement an exported safe wrapper/alias to the certified route if retaining the requested name. Use the authenticated embedded NB contract; no unrestricted refinement or classical inference. |
| Multinomial | Dispatch supported qualified additive formulas in `ds.vertMultinom`/alias to the authenticated embedded multinomial grid. Preserve signed levels/reference category, parameter order and fixed candidate set; retain both existing same-owner Frequency/intercept-only and finite-grid behavior. |
| Ordinal | Dispatch supported qualified additive formulas in `ds.vertOrdinal`/alias to the authenticated embedded ordinal grid. Preserve signed level order/link/cutpoints/candidate set and current same-owner route. |
| LASSO binomial/Poisson | Export `dp_lasso_grid`. Extend the appropriate `ds.vertLASSOIter` branch with appended `cross_contract=NULL, base_contract=NULL` arguments for separately signed extension/base contracts, not guessed discovery from `analysis_id`. Preserve the signed lambda, family, beta grid and capacity-normalized public penalty/selection contract. `analysis_id` currently does not carry the LASSO extension. |

For NB/multinomial/ordinal, `analysis_id` may select the already embedded cross contract in the authenticated Synopsis; derive schema/policy from that verified bootstrap and match the owner-qualified formula before private execution. No implicit rewrite of bare formulas or unsigned grid creation. Low-level adapters keep their current explicit signed-contract API. Unsupported interactions/offsets/grids/families refuse locally when detectable; no fallback to exact kernels.

Gaussian LASSO is a separate existing route. `ds.vertLASSOCapsule.R`'s shared source verifier is ownership-agnostic and validates anchored Gaussian certificates; `ds.vertLASSOProximal` and `ds.vertLASSOCV` appear already to accept cross Gaussian fits. `ds.vertLASSO` has an explicit same-owner gate; `ds.vertLASSO1Step` and Gaussian Iter reach it. Admit only validated current cross Gaussian provenance through that gate and retain unchanged solver behavior, with direct tests for all affected front doors. Do not claim binomial/Poisson CV or arbitrary iterative optimization merely because a finite grid is public. Where a front door's mathematical operation does not match the signed artifact, document its narrower scope and retain refusal.

The `release_enabled`, `production_release_enabled`, and NB `production_release_registered` FALSE fields occur in both client and server registration metadata. Static call/hash inspection finds no production reader and no inclusion in signed computational descriptors. Replace these ambiguous unused flags with explicit descriptive promotion/production-route metadata and update stale comments/tests. In particular the server multinomial/ordinal descriptor `release` callbacks deliberately fail closed: retain that behavior and identify the real authenticated Synopsis lifecycle, rather than imply that those callbacks now release. Actual admission remains signatures, schema, mechanism and protocol capability; no signed artifact or release changes.

## Registry, documentation and versioning

Correct `ds.vertMethodStatus` for cross Gaussian, Cor/PCA and the four new exports, retaining one row per public name and precise supported scopes/limitations. Update safe NB alias and LASSO rows as applicable. Replace stale same-owner-only text in roxygen/Rd, README tables, product surface/migration docs and examples. Keep archived campaign documents immutable; label their evidence dates. Update both NEWS files with the actual partial-preexisting support, cross PCA fix, explicit complete-case/offline API, grid exports, safe frontdoor dispatch and provenance invariants. Do not claim all historical methods are unrestricted cross-owner estimators.

Both R packages are **1.4.2** because the four server register descriptors/documentation change as well as client APIs. No server arithmetic, sampler, endpoint allowlist or Go wire change is proposed. Preserve the exact provenance capability object, including its compatibility field `version="1.4.1"` (the current preflight compares the whole object), and existing `dsvert-provenance-epochs-v1`, cross Gaussian artifact, exact-GC, Synopsis and grid capabilities/versions. Do not invent a new capability string for client-only post-processing. A compatible 1.4.1 server may serve these unchanged contracts if it advertises every required capability; verify capability content, not a patch-number assumption. Deploying the paired 1.4.2 packages gives consistent registry metadata. Missing/unknown capabilities, incompatible descriptors and absent evidence fail closed before materialization; no legacy downgrade.

If implementation review finds a genuinely necessary server semantic change, revisit the wire-version/identity migration explicitly instead of hiding it in 1.4.2. Preserve old artifact readers and byte replay for retained old domains. No commits in phase A; phase B commits are small English changes with the requested Claude-Session trailer and are pushed only by Claude.

## Targeted validation plan

Every command has `timeout 900`; never run full suites or `go test ./...`. Record actual assertions, warnings, skips and deadlines, not only process exit codes. Phase-A commands/results are recorded below and in the isolated ENV.md.

1. Current signed Synopsis fixtures at K=2,3,5 for Cor/PCA and aliases; cover both K-peer federations with fewer data owners and variables genuinely distributed across all K owners (m>=K, plus an extra column for the mixed same-owner pair). Use at least three variables, including same-owner and cross-owner predictor–predictor/outcome pairs. Include K=5 with membership-only participants and computation peers distinct from some owners. Compare deterministic normalized bounded-moment references, matrix symmetry/PSD/diagonal, eigenvalue sum and reconstruction, rank-deficient valid cases, subset ordering and repeated-eigenvalue subspaces.
2. Missingness/quantization: repeated records, clipping before mean, nonfinite values, all-variable complete cases, subsets retaining the signed cohort, no-intercept refusal, zero variance/nonpositive DP mass, tight/nonvacuous cross intervals and broad intervals. Test duplicate physical names qualified by owner/dataset without dropping bindings.
3. Tamper individual coordinates, owners, bounds, canonical names, signatures/pins, source/layout/order/transcript/numeric hashes, provenance epochs and authority evidence (missing, nested incorrectly, inconsistent, wrong peer). Recompute outer hashes in attacks. Exercise real verifier paths; mocks only for isolated algebra tests.
4. Freeze all DSI discovery/aggregate hooks to error for offline Cor and supplied-result PCA. Test valid external pins after client restart, absent pins/unanchored failure, mismatched supplied identifiers, certificates without successful regression, and no rerandomization.
5. Replay: exact artifact, final-vector, release/certificate and same-owner result bytes with fixed fixtures. Retained authorized source mutation/reconnect/restart replays; new authorization (including membership-only dependency) changes affected identity. Test subset/order/rank choices do not create new private releases. Keep the 1.4.1 provenance tests.
6. Grid exports: check NAMESPACE and `getExportedValue` for four names; every registry export represented once; flags/comments consistent. Test authenticated frontdoor parity and zero unexpected calls for NB/multinomial/ordinal; LASSO binomial/Poisson contract mismatch and extension tampering; all current same-owner and Gaussian LASSO routes unchanged.
7. Server files: Gaussian cross, capsule workload/materializer (correlation descriptors included), relevant Synopsis source/attestation/replay/provenance files; four-family registration tests when flags change. Client files: current Cor/PCA/Gaussian/transport/certificate/method-status and matching cross-grid files. Add only focused tests for the changed behavior, not broad suite expansion.
8. **DSLite smoke:** use separate custodian R processes with owner-only persistent state outside TMPDIR/library, bounded public synthetic data and tiny signed grids. The original phase-B budget was at most ten fit/analysis invocations and 15 minutes total; its attempts and limits remain part of the historical validation record. For the vertcor03 debugging follow-up, the maintainer explicitly authorized **45 minutes cumulative smoke wall time across attempts**, with **every command still capped at 900 seconds**. Prioritize K2 Cor+PCA, then K3 Cor+PCA, K5 if time permits, one binomial and one Poisson LASSO grid release, exact Cor replay, and the offline Cor-to-PCA chain with zero DSI calls. Keep signing, the allowlist, source authorization, PSI and the production release path real; use no fixture seams. Record exact outcomes, timings and cumulative wall time in [RESULTS.md](inst/validation/v1.4.2/RESULTS.md). A deadline or partial run remains incomplete evidence; never label mocked algebra as federation validation. Other grid families retain targeted unit/adversarial tests and their named historical promotion evidence.

Identity-sensitive tests must not treat ordinary one-process DSLite as separate persistent custodians. A short DSLite smoke is a local connector check, not a live Opal/Rock release guarantee. Do not run the historical multi-hour campaign in this phase.

## Phase-A baseline record

Environment, exact commands, dependency/runtime versions, logs and machine-readable per-file results: `/workspace/dsflower/.tmp/vertcor/ENV.md`; source `/workspace/dsflower/.tmp/vertcor/env.sh`. Unmodified CRAN DSI 1.8.0 requires R >=4.3, so the isolated runtime uses user-space R 4.3.3 while system R 4.2.2 is untouched. Go 1.25.7 builds the required runtime stamp 1.4.0/API 1.2.0, byte-identical to the committed Linux binary. Both working-tree R packages installed as 1.4.1. No implementation changes preceded baseline.

Final corrected-environment baseline: **20 client files, 164 test blocks, 6,425 passing assertions, zero failures/errors/warnings, one optional Cox plaintext-harness skip** (survival and the tagged test binary were not provisioned). All four requested grid-family files execute and pass. **Eight server files, 126 test blocks, 1,523 passing assertions, zero failed assertions/warnings/skips, one error** in the existing capsule-materializer cross Gaussian alignment fixture. Commands and original environment failures/retries are retained; final per-file counts are in `.tmp/vertcor/BASELINE.md` and `baseline-summary.csv`. Every test batch completed below 900 seconds; the longest summed file time was 388.2 seconds.

The server error occurs constructing the second logical peer, before Gaussian materialization: the helper reuses one actual persistent identity and deterministic authorization/seal context for two different local frames. The 1.4.1 immutable-value seal correctly refuses this. Phase B isolates the two simulated identities during both creation and validation, with matching real-derived test pins and common PSI epoch vector. Do not vary the common token to evade the check, clear retained pins, disable seals or weaken production code. See `.tmp/vertcor/baseline-materializer-diagnosis.md`. The repaired test also checks the current v5 manifest, both value/snapshot digests and continued rejection of identity reuse with changed values.

Read-only Cor/PCA reproduction confirms the cross interval mismatch through installed public functions (retrieval/verifier fixture mocks). A separate current server/client descriptor diagnostic confirms the duplicate-name failure. Deterministic same-owner Cor/PCA/alias serialization was captured and repeated identically: `.tmp/vertcor/baseline-same-owner-cor-pca.rds`, serialized SHA256 `95ec566da9db7db10d1bf6f69f5203d847ed375e2df4021e62494d392102269a`. This is a client-output golden, not a signed server replay claim.

## Consensus and final implementation decisions

Claude approved Gaussian reuse, fixed public anchor, explicit Gaussian DPCor opt-in,
certificate-only offline Cor/PCA, safe NB alias, explicit LASSO extension contracts,
fixture identity repair, and coordinated 1.4.2 with unchanged native/protocol identities.
There is no additional external account of the August revert; the historical
explanation above is a code-supported readiness finding, not an asserted motive.

Implemented Cor/PCA reads only authenticated coordinates and recomputes matrices,
intervals and PSD projection. Cross-only quantization and canonical duplicate-name
repairs preserve every valid same-owner artifact. Live DSLite investigation also
exposed client v2 workload references being hashed before the server
canonicalization step. The client now resolves only those v2 references against
the signed global schema before hashing, matching existing server keys/bytes;
no producer identity or already-valid artifact is changed. The independent server
local-fragment equality check remains unchanged: raw registration uses bare
labels for globally unique columns and qualified labels for colliding names.
Qualified analyst references and the frozen owner-qualified anchor identity
remain supported. Raw qualified aliases for unique registration columns still
fail closed; this is documented rather than weakening the server comparator. Offline certificates require
caller/session anchoring, and no supplied matrix or validity flag establishes trust.
PCA returns loadings, eigenvalues and variance; no scores are produced.

NB, multinomial and ordinal front doors discover the authenticated embedded grid
when an owner-qualified formula and analysis_id select it. Existing bare-name and
Frequency routes retain their behavior. `ds.vertNB` is an exact alias of the safe
wrapper; the historical body is internal. LASSO extensions require both signed
cross/base contracts; Gaussian LASSO accepts current cross certificates only.
Confint/Contrast/Wald/LR already supported the authenticated cross moments; their
docs and direct tests now establish this scope without changing their arithmetic.

Cox/LMM/GLMM/GEE retain explicit `dp_*_grid` cross entry points: their richer
front-door risk-set, ML/REML/random-slope and covariance controls do not map exactly
to the promoted cross contracts. `ds.vertOrdinalJointNewton` remains Frequency-only.
Current MI and survival Synopsis producers remain same-owner for the concrete
missingness/producer reasons in the [complete public inventory](README.md#public-cross-owner-inventory-142).
Thus the unrestricted thesis statement about every family is not claimed.

The server changes only promotion metadata, documentation, version and the test
fixture. Go source/binary bytes, capabilities (including provenance version 1.4.1),
release coordinates, noise, identity domains and endpoint allowlists are unchanged.
The metadata no longer presents unused FALSE release flags as runtime gates.

See [worked custodian configurations](inst/docs/cross_owner_correlation_1.4.2.md)
and [validation results](inst/validation/v1.4.2/RESULTS.md). The current certificate
unit fixture signs real release/replay/authority evidence but uses explicit
bundle/compiler/calibration seams; it is not live-connector evidence. Any failed
DSLite run is recorded as incomplete, without weakening production validation.
