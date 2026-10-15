# 1.4.2 targeted validation — 15 October 2026

## 1. Environment and command limits

Source `/workspace/dsflower/.tmp/vertcor/env.sh`; the complete environment record
is `.tmp/vertcor/ENV.md`. R 4.3.3 was built in user space because current CRAN DSI
requires R >= 4.3; no system installation was changed. Dependencies and Go 1.25.7
are isolated below `.tmp/vertcor`. All test/install commands are capped at 900
seconds. No full R suite, full Go sweep, GitHub issue operation or push was run.
Both working trees install as 1.4.2. The Linux native binary remains SHA256
`52cde05d6874d274b60e0772f37596017b344d4a1b90da86ed089995f1ab9784`.
After the follow-up fixes, both runner files parsed and both packages installed
again successfully (server 65.253471 s, client 50.836494 s). Installed versions
remain 1.4.2; active, source-tree and installed native binaries have the same
SHA256. Commands and results are appended to `ENV.md`.

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
LASSO promotion campaign is binomial; Poisson contract/unit coverage is new here.
The bounded live follow-up outcomes are reported separately below; a tiny smoke
does not establish a separate Poisson promotion campaign or statistical accuracy.

## 6. Results and commands

Final unique-file totals: **8,654 client assertions / 32 files / 239 blocks** and
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
| dsVertClient | test-dp-gaussian-cross-transport.R | 10 | 136 | 0 | 0 | 0 | 0 |
| dsVertClient | test-dp-gaussian-synopsis-provenance.R | 2 | 20 | 0 | 0 | 0 | 0 |
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

## 7. Separate-process DSLite smoke — vertcor03 follow-up

**K2 Cor/PCA, exact replay, offline zero-DSI conversion and binomial LASSO
passed. The full requested release gate is still incomplete: K3 Cor/PCA was
non-identifiable, Poisson did not finish before the cumulative deadline, and
optional K5 was not attempted.** No further smoke was run after the deadline.

### Production diagnosis and fixes

The original `exact MPC pinned-peer binding` error was a client production defect,
not a requirement to bypass the harness or server policy. Gaussian created a fresh
exact-computation UUID without first installing signed Synopsis authorization.
`exactGCBindPeersDS` therefore entered the legacy policy branch and failed with
`dsvert.dp.ledger_path must be one non-empty string` at both custodians. The
cross-grid transport already called PREPARE on its child session. Commit
`a957113` applies that same authenticated ordering to Gaussian, before transport
initialization, without changing pins, signatures, compilation or source checks.

Once the real exact computation completed, certificate validation exposed a
second client omission: reconstruction did not populate
`context$authorization_dependencies`, so the signed closure was compared with
`NULL`. Commit `f6f123c` derives the closure from signed federation evidence before
the unchanged trusted-bundle validator authenticates it. This also repairs
current same-owner Gaussian Synopsis certificates and the shared grid reader.
Neither change emits a new artifact, changes sticky identities or weakens epochs.

Both omissions are present in client 1.4.1 (`04a6903`). The standard path tested
here could not complete with that code. This is not evidence that every
historical deployment or legacy Gaussian route failed. The earlier
`dsvert_phase_not_ready` record was expected publication polling: it also occurs
in the successful K2 run. The schema normalization fix `1846e3e` from phase B
remains necessary; server canonical comparators are unchanged.

The worked v1.2.0/v1.3.0 custodian scaffolds are identical. Unlike the historical
promotion harness, this smoke uses no oracle pre-provisioning, native-output
capture wrapper, substituted compiler/calibration or mocked endpoint. Custodian
administration provisions identities, trusted peers, source authorization and
signed public contracts in each independent process. All analyst requests use
real PSI, the actual allowlist and production release endpoints. Diagnostic
observers record endpoint errors and already-public DP certificates without
replacing arguments, validation or results.

### Live outcomes

| Check | Outcome | Measured analysis seconds |
|---|---|---:|
| K2 Cor, 8 synthetic patients / 2 owners / 2 variables | PASS | 361.0614 |
| K2 PCA | PASS; loadings/variance, no scores | 2.0341 |
| Exact K2 Cor replay | PASS; serialized correlation/artifact/certificate and embedded release/replay JSON identical | 43.0975 |
| Offline Cor → PCA | PASS; externally supplied pins; zero aggregate, assign or connection-discovery calls | 4.1262 |
| K3 Cor, 4 synthetic patients / 3 owners / 3 variables | Signed DP release reached reconstruction; rejected non-positive DP complete-case count | 444.9315 |
| K3 PCA | Not run because Cor was non-identifiable | — |
| K2 binomial LASSO, 1 synthetic patient / 1 predictor | PASS via `ds.vertLASSOIter`; authenticated certificate and two penalty paths | 256.2188 |
| K2 Poisson LASSO | Signed registration passed (11.5835 s); public `dp_lasso_grid` call started, then deadline | approximately 39.84 in flight |
| Optional K5 Cor/PCA | Not attempted | — |

K2 startup plus real PSI took 51.2540 s; its complete command took **464.0009 s**.
K3 startup plus PSI took 65.4261 s; its command took **512.6773 s**. Do not interpret
K3's signed noisy count as a successful Cor/PCA result or replace it with a guessed
count. Its certificate is retained locally for inspection. Larger public
synthetic cohorts improve the chance of identifiability but require another
bounded run; the default K3 runner now uses eight rows, and `DSVERT_SMOKE_ROWS=4`
recreates the tested cohort shape (not a fresh noise draw).

The LASSO cases shared two custodian processes and one real PSI run (49.6051 s).
Binomial registration took 26.6559 s; its registration plus verified fit took
282.8852 s. Each family signs a distinct workload and requires its own current
certificate. Both planned releases use epsilon 4 (sum 8), delta 2^-100 each,
16-bit quantization and two beta candidates at penalties 0.5 and 0. No ledger,
identity, aligned data or authorization epoch is reset between workloads. This is
a tiny lifecycle check, not statistical-accuracy or Poisson promotion evidence.
The combined command stopped at **386.1519 s**, exit 124. No Poisson result or
certificate-validation success is claimed.

### Time accounting and reproducibility

| Follow-up attempt | Command wall seconds | Outcome |
|---|---:|---|
| `binding-before` | 163.6104 | Actual legacy-ledger error captured |
| `cor2-after-prepare` | 1.8673 | Runner's non-logical `stopifnot(length(cases))` rejected before federation |
| `cor2-prepare-2` | 476.3332 | Larger n64/three-variable diagnostic stopped during exact computation |
| `cor2-small` | 36.9762 | Runner configured PSI below its minimum capacity; rejected before release |
| `cor2-small-psi64` | 523.4151 | Real n16 computation completed; certificate closure omission exposed |
| `cor2-closure` | 464.0009 | K2 Cor/PCA/replay/offline PASS |
| `cor3` | 512.6773 | Released DP count non-positive; Cor/PCA incomplete |
| `lasso-combined` | 386.1519 | Binomial PASS; Poisson public call reached deadline |

Measured command wall time is **2,565.0322 s (42m45.03s)**. The stopped n64
attempt left native workers waiting after its R process exited; those workers
were explicitly terminated. Their last heartbeat and fixture creation times
bound that interval at about 124 s. A conservative **130 s cleanup allowance**
is charged without deducting its overlap with the next attempt. Total charged
smoke time is therefore **2,695.0322 s (44m55.03s)**, below the authorized 45 minutes.
Every command was below 900 s; no native worker remained after the final stop.

The measured records are committed as `vertcor03-smoke-ledger.jsonl` and
`vertcor03-analysis-events.jsonl`. The latter records completed analysis events;
termination can occur before an event is written. Logs are
`.tmp/vertcor/logs/vertcor03-<attempt>.log`; successful public outputs and the K3
certificate are `.tmp/vertcor/smoke-*.rds`. The bounded local launcher recorded
commands of this form (sequential invocations, cumulative remaining-time cap):

```sh
source /workspace/dsflower/.tmp/vertcor/env.sh
DSVERT_SMOKE_CASES=cor2 timeout 900 python3 ../.tmp/vertcor/run-smoke03.py \
  cor2-closure Rscript tools/validate_cross_owner_142.R /workspace/dsflower
DSVERT_SMOKE_CASES=cor3 timeout 900 python3 ../.tmp/vertcor/run-smoke03.py \
  cor3 Rscript tools/validate_cross_owner_142.R /workspace/dsflower
DSVERT_SMOKE_CASES=lasso timeout 900 python3 ../.tmp/vertcor/run-smoke03.py \
  lasso-combined Rscript tools/validate_cross_owner_142.R /workspace/dsflower
```

Run from the client root. The runner evolved during diagnosis; the effective
cohort sizes are recorded above. For a fresh isolated case use `binomial` or
`poisson`; `lasso` shares startup but still performs two independent signed
releases. No additional attempt is authorized by the exhausted follow-up budget.

The preceding phase-B smoke remains historical evidence: seven attempts, six
analysis requests, 574.64430 s, no completed Cor/PCA result. Its last unresolved
binding failure motivated this follow-up; it is not counted in the newly
authorized 45-minute budget.

### Focused regression checks after the production fixes

`timeout 900 Rscript ../.tmp/vertcor/run-targeted.R dsVertClient` was run with
`test-dp-gaussian-cross-transport.R` (136 passes / 10 blocks / 7.671469 s),
`test-dp-gaussian-synopsis-provenance.R` (20 / 2 / 7.207181 s), and
`test-dp-gaussian-synopsis-release.R` (6 / 1 / 0.144301 s). All **162 assertions**
passed, zero failures/errors/warnings/skips. New closure tests use real signed
bootstrap/manifest evidence, actual closure/trusted-bundle validation, K2/K3/K5,
changed epochs and missing epochs. The unique-file totals in section 6 include
these reruns once. No server production function changed.

The same-owner command reconfirmed all five byte comparisons in section 4:

`timeout 900 Rscript tools/validate_same_owner_142.R ../.tmp/vertcor/vertcor03-same-owner 04a6903`

Logs retain the SHA256 values after both fixes; the fixture seams
remain explicitly documented. These byte checks are separate from live evidence.

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
