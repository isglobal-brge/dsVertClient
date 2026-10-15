# 1.4.2 targeted validation — 15 October 2026

## 1. Environment and command limits

Source `/workspace/dsflower/.tmp/vertcor/env.sh`; the complete environment record
is `.tmp/vertcor/ENV.md`. R 4.3.3 was built in user space because current CRAN DSI
requires R >= 4.3; no system installation was changed. Dependencies and Go 1.25.7
are isolated below `.tmp/vertcor`. All test/install commands are capped at 900
seconds. No full R suite, full Go sweep, GitHub issue operation or push was run.
Both working trees install as 1.4.2. The Linux native binary remains SHA256
`52cde05d6874d274b60e0772f37596017b344d4a1b90da86ed089995f1ab9784`.

## 2. Correlation, PCA and Gaussian tests

K=2,3,5 cover mixed same-/cross-owner variable pairs, all requested products,
qualified duplicate physical names, covariance/correlation reconstruction,
nonvacuous cross quantization intervals, singular-regression certificate-only
PCA, and rejection of no-intercept/zero-variance inputs. PCA never returns scores.
Client v2 canonical workload/schema tests compare unchanged server reconstruction
bytes, including duplicate names and owner/dataset ambiguity.

## 3. Authentication, privacy and replay

Tests reject release/replay/attestation/coordinate tampering and unanchored offline
inputs, check zero DSI/discovery/MPC calls for authenticated offline conversion,
and retain sticky output identity. Existing server provenance epoch, execution,
source and remote-surface tests pass unchanged. The repaired materializer fixture
uses distinct identity seeds and corresponding pins; reusing one identity with
changed authorized values still fails. Current certificate unit fixtures use real
Ed25519 release/replay and cross-attestation signatures, with explicitly documented
bundle/compiler/calibration seams. They are not live connector evidence.

## 4. Same-owner byte identity

The original Cor/PCA/alias golden remains exactly
`95ec566da9db7db10d1bf6f69f5203d847ed375e2df4021e62494d392102269a`.
The pre-edit four-LASSO golden remains
`326ade197332d6d7b40928f0daac5471f75f0e3396b94daeccd41e14fb9a9a8e`.

`timeout 900 Rscript tools/validate_same_owner_142.R <output_dir> 04a6903`
compares baseline definitions with current code on the same fixed signed fixture.
It confirms serialized identity for artifact (2,074 bytes), raw release (15,885),
certificate (11,680), Gaussian fit (32,109) and all eight consumer outputs (189,634).
It uses the documented current-fixture seams above; it does not claim a fresh
production noise draw or connector replay. Already-valid bare and qualified
homonym schema/workload bytes also match 04a6903 exactly after normalization.

## 5. Public grids and front-door parity

NB/LASSO/multinomial/ordinal exports and K=2,3,5 signed-contract dispatch tests pass.
The safe NB alias cannot dispatch the historical body. Categorical domain/order,
formula controls, signed base/extension equality and exact lambda subsets are
checked before release. Null extension arguments preserve same-owner dispatch.
Cox/LMM/GLMM/GEE keep their existing explicit grid entry points. The archived
LASSO promotion campaign is binomial; Poisson contract/unit coverage is new here,
but no fresh Poisson connector fit or separate promotion campaign is claimed.

## 6. Results and commands

Final unique-file totals: **8,595 client assertions / 31 files / 236 blocks** and
**6,382 server assertions / 11 files / 163 blocks**, zero failures/errors/warnings.
One client skip: the optional tagged Cox plaintext harness needs `survival`, which
is not installed. All 21 changed Rd files parse. Both package installations and
`git diff --check` are checked separately.

Each file below was run with the targeted helper:

```sh
source /workspace/dsflower/.tmp/vertcor/env.sh
timeout 900 Rscript /workspace/dsflower/.tmp/vertcor/run-targeted.R dsVertClient <test-file>...
timeout 900 Rscript /workspace/dsflower/.tmp/vertcor/run-targeted.R dsVert <test-file>...
```

The helper loads the package and its test helpers, then calls `testthat::test_file`
only for the named files. Per-file CSV/RDS and command logs are in
`.tmp/vertcor/logs`; `targeted-summary.csv` preserves the final unique-file counts.
Early red runs exposed stale registry counts and fixture assertions, plus test-only
alias metadata comparisons; their final reruns are counted once. Two grid batches
ended on mistakenly named nonexistent tests after real files passed; those missing
files are not counted (see `.tmp/vertcor/GRID_PHASE_B.md`).

| Package | Targeted file | Blocks | Pass | Fail | Error | Skip | Warn |
|---|---|---:|---:|---:|---:|---:|---:|
| dsVert | test-crossgrid-family-b.R | 13 | 1045 | 0 | 0 | 0 | 0 |
| dsVert | test-crossgrid-lasso.R | 11 | 118 | 0 | 0 | 0 | 0 |
| dsVert | test-crossgrid-nb.R | 13 | 3689 | 0 | 0 | 0 | 0 |
| dsVert | test-dp-capsule-materializer.R | 57 | 674 | 0 | 0 | 0 | 0 |
| dsVert | test-dp-capsule-workload.R | 28 | 384 | 0 | 0 | 0 | 0 |
| dsVert | test-dp-gaussian-cross.R | 6 | 39 | 0 | 0 | 0 | 0 |
| dsVert | test-dp-synopsis-exact-downstream.R | 8 | 172 | 0 | 0 | 0 | 0 |
| dsVert | test-dp-synopsis-execution-safety.R | 4 | 18 | 0 | 0 | 0 | 0 |
| dsVert | test-dp-synopsis-provenance-epochs.R | 6 | 56 | 0 | 0 | 0 | 0 |
| dsVert | test-dp-synopsis-remote-surface-red.R | 12 | 111 | 0 | 0 | 0 | 0 |
| dsVert | test-dp-synopsis-source-gate.R | 5 | 76 | 0 | 0 | 0 | 0 |
| dsVertClient | test-capsule-method-inventory.R | 13 | 966 | 0 | 0 | 0 | 0 |
| dsVertClient | test-cross-gaussian-schema.R | 5 | 54 | 0 | 0 | 0 | 0 |
| dsVertClient | test-cross-owner-cor-pca.R | 7 | 128 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-capsule-manifest-client.R | 4 | 96 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-cor-pca.R | 17 | 222 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-cox-grid-cross-e2e.R | 1 | 0 | 0 | 0 | 1 | 0 |
| dsVertClient | test-dp-cox-grid-cross.R | 10 | 190 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-gaussian-cross-transport.R | 9 | 97 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-gaussian-synopsis-release.R | 1 | 6 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-gaussian.R | 15 | 329 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-glm-grid-cross-contract.R | 20 | 360 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-glm-grid-cross-integer.R | 6 | 178 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-glm-grid-cross-sensitivity.R | 1 | 84 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-grid-cross-frontdoors.R | 7 | 205 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-grouped-grid-cross-registry.R | 2 | 24 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-grouped-grid-cross.R | 8 | 176 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-lasso-grid-cross.R | 13 | 122 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-multinomial-grid-cross-authorities.R | 4 | 48 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-multinomial-grid-cross.R | 7 | 318 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-nb-grid-cross.R | 15 | 3701 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-ordinal-grid-cross.R | 2 | 23 | 0 | 0 | 0 | 0 |
| dsVertClient | test-formal-multinom-frontdoor.R | 5 | 46 | 0 | 0 | 0 | 0 |
| dsVertClient | test-formal-nb2-frontdoor.R | 4 | 30 | 0 | 0 | 0 | 0 |
| dsVertClient | test-formal-ordinal-frontdoor.R | 5 | 41 | 0 | 0 | 0 | 0 |
| dsVertClient | test-gaussian-cross-consumers.R | 5 | 557 | 0 | 0 | 0 | 0 |
| dsVertClient | test-lasso-contract.R | 20 | 324 | 0 | 0 | 0 | 0 |
| dsVertClient | test-method-status.R | 3 | 167 | 0 | 0 | 0 | 0 |
| dsVertClient | test-pca-contract.R | 3 | 11 | 0 | 0 | 0 | 0 |
| dsVertClient | test-signed-correlation-consumers.R | 7 | 47 | 0 | 0 | 0 | 0 |
| dsVertClient | test-vertContrast.R | 6 | 14 | 0 | 0 | 0 | 0 |
| dsVertClient | test-vertInference.R | 11 | 31 | 0 | 0 | 0 | 0 |

## 7. Separate-process DSLite smoke

**Incomplete: zero completed Cor/PCA results; six analysis requests across seven
bounded setup/diagnostic attempts; 574.64430 cumulative runner seconds
(about 9m35s).** K=2 used independent persistent custodian R processes and real
PSI. K=3/K=5 were not reached. No fresh live Poisson LASSO fit was attempted.

The early fixture used an inadmissible epsilon=128; the unchanged safe policy
correctly rejected it (maximum 8). The corrected epsilon=8 run exposed differing
client/server v2 schema/workload hashes; the client canonicalization fix now matches
both captured server hashes and the new direct server-parity tests. A raw fragment
using qualified aliases for globally unique columns then hit the unchanged server
exact local-fragment check. Canonical registration labels resolve that setup issue.

The final attempt passed schema and authoritative manifest binding and failed at
`exact MPC pinned-peer binding`, with no completed correlation/PCA result. The
last recorded local protected condition was `dsvert_phase_not_ready`; that record
is not proof of the binding failure's cause. This transport/runtime stage was not
changed or bypassed. Live connector completion remains a release-verification
risk for Claude, and unit passes must not be presented as federation success.

The final bounded invocation was:

```sh
DSVERT_SMOKE_PRIOR_SECONDS=407.5 DSVERT_SMOKE_PRIOR_ANALYSES=5 \
  timeout 470 Rscript dsVertClient/tools/validate_cross_owner_142.R /workspace/dsflower
```

Earlier attempts used the same runner under timeout 900, 810, 750, 680, 610 and
530 seconds, exiting early. The prior-time/analysis counters bounded the remaining
attempts jointly; the runner was stopped after the failure, within both caps.
Logs: `.tmp/vertcor/logs/phase-b-dslite-smoke*.log`.

The synthetic-only diagnostic records public binding metadata and local fixture
errors, never real participant data. The committed runner is
`tools/validate_cross_owner_142.R`; it uses the existing process-isolated custodian
scaffold and the actual allowlist, source authorization, pinned PSI, signatures
and production release path. Unit fixture seams are not used in this runner.

## 8. Compatibility and release boundaries

Both R package versions are 1.4.2. Server production changes are descriptive grid
metadata only; the source/provenance/attestation comparators stay unchanged.
Native Go source/binary, capability objects (including provenance version 1.4.1),
endpoint allowlists, sampler parameters and existing identity domains are unchanged.
The client fixes v2 schema reconstruction to match the existing server normal form;
already-valid contracts retain exact bytes. Custodians register canonical labels
(bare globally unique names, qualified collisions); the server still rejects raw
noncanonical qualified aliases for unique registered columns. Qualified analyst
references and the fixed owner-qualified anchor identity remain supported.

The README and NEWS contain the full public family inventory. MI and survival
Synopsis producers retain their explicit same-owner restrictions; the release
does not claim that every possible family scope is cross-owner.
