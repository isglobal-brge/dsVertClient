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

Final unique-file totals: **8,671 client assertions / 33 files / 242 blocks** and
**6,382 server assertions / 11 files / 163 blocks**, zero failures/errors/warnings.
One client skip: the optional tagged Cox plaintext harness needs `survival`, which
is not installed. All 21 changed Rd files parse. Both package installations and
`git diff --check` are checked separately.

The production test files below were run with the targeted helper; the added
smoke-fixture file uses the direct targeted command in section 9:

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
| dsVertClient | test-cross-owner-smoke-fixture.R | 3 | 17 | 0 | 0 | 0 | 0 |
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

**At the end of vertcor03, K2 Cor/PCA, exact replay, offline zero-DSI conversion
and binomial LASSO passed. The requested release gate remained incomplete: K3 Cor/PCA was
non-identifiable, Poisson did not finish before the cumulative deadline, and
optional K5 was not attempted.** No further smoke was run under that exhausted
budget. Section 9 records the separately authorized vertcor04 follow-up.

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
bounded run. The runner at commit `8d6b83d` used eight rows by default, with
`DSVERT_SMOKE_ROWS=4` recreating the tested cohort shape. Section 9 describes
the subsequent four-variable n64 fixture.

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
releases. No additional attempt was made under that exhausted follow-up budget.

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

## 9. Sized live gates (`vertcor04`)

This is a separately authorized 2,700-second campaign. Every command is limited
to 900 seconds, including cleanup. The new tracked launcher
`tools/run_cross_owner_smoke.py` keeps an exclusive campaign lock, durable
start/completion records and conservative reservations for interrupted attempts.
Linux subreaper adoption tracks detached descendants; generation-checked pidfd
signals clean up only processes belonging to this attempt. Five focused Python
tests pass (3.612 s), including detached-child cleanup and preservation of an
unrelated process. The fixture test adds 17 passing R assertions in three blocks,
zero failures/errors/warnings/skips. These tests do not count as live federation
evidence or against the live-smoke budget.

### K3 cohort selection

The one bounded K3 attempt uses **n = unit capacity = 64**, epsilon 8, delta 0,
8-bit quantization, add/remove-patient adjacency and bounds [0,1]. The requested
variables are `site_a$a`, `site_a$b`, `site_b$x`, `site_c$x`: one same-owner pair,
all three data owners and cross-owner pairs. Balanced binary columns have true
centered sums of squares 16 and zero pairwise centered cross-products. The fixed
public anchor is `site_a$a`.

The existing compiler admits 16 Gaussian coordinates, 12 automatic numeric
moment coordinates and one global count. Total L1 sensitivity is 29 normalized
units (7,424 lattice steps); the two independent discrete-Laplace authorities
give pre-clamp additive-noise SD 7.24999965 per normalized coordinate at epsilon
8. Clamped/projected released moments need not retain that SD. No sensitivity,
epsilon, sampler, clipping, attestation or identity rule is changed.

Public planning through the unchanged Gaussian unpacker and Cor identifiability
gates gave 49,889/50,000 finite results at n32 and 50,000/50,000 at n64, n128 and
n256. At n64 the exact ideal nonpositive-count probability is 1.05705e-7; the
95% binomial interval for simulated finite output is [99.9926%, 100%]. Choosing
n64 provides more margin than n32 while limiting exact-computation cost.
This supports a functional gate, not a precision guarantee: the simultaneous
95% noise radius is 26.484375 over the 16 Gaussian sums, and balanced population
PCA directions are not uniquely identified. The n64 size was fixed before this
campaign from the earlier certificate's public calibration and public noise
simulation. It was not tuned to newly released values or retried for a favorable
draw. The reproducible planning script is
`tools/plan_cross_owner_142_noise.R`; exact counts, calibration and limitations
are in `vertcor04-noise-planning.json`. Public planning took 67.428 s, outside
the live-smoke budget.

### Earlier Poisson timeout diagnosis

The previous Poisson interruption came from the external cumulative smoke
budget, **not a 40-second deployment timeout**. The old launcher allocated
386.1197 seconds to the entire combined binomial/Poisson command; measured wall
time was 386.1519 seconds, exit 124. Only about 39.84 seconds remained after the
Poisson public call began. The client exact-GC inactivity default is 900 seconds,
capped by the signed server lease (default 180 seconds); authenticated progress
renews it. The server default maximum runtime is 21,600 seconds. The DSLite
scaffold's 1,000-ms polling interval is not an overall deadline. There is no
reason to raise a production timeout for this case.

The independent Poisson fixture provisions fresh K2 custodian state, one patient,
one predictor, epsilon 8, delta 2^-100 and one permitted artifact. It signs the
two-candidate base grid and two-penalty LASSO extension and uses the public
`dp_lasso_grid` entry point. Objective normalization uses public capacity, so it
does not require a positive noisy Gaussian count. The completed run establishes
protocol operation on this tiny fixture, not statistical utility or promotion
campaign coverage.

### K3 bounded outcome and runtime diagnosis

The single n64/m4 K3 command ran from 13:21:16.946751 to
13:36:14.083548 UTC, **897.1367984600365 seconds**, and stopped with exit 124
(`command_limit`). The launcher reserved three seconds of its 900-second
allowance for cleanup; all owned processes exited, including detached native
workers (`cleanup_complete=true`, no remaining PIDs). Setup and real PSI passed
in 159.1586 seconds. Cor entered its production release path, but neither a
completed Cor/PCA result nor a public certificate was returned before the limit.
The analysis event ledger contains completed events only; the timeout itself is
recorded in the command ledger.

Read-only process and public operation metadata showed forward progress, not an
inactivity timeout or a pin/signature/authorization rejection. Alignment-mask
execution completed at about 13:30:18, and subsequent exact-GC multiplications
ran. The first native operation was active by 13:28:29; the third had started by
13:34:24. Diagnostics inspected only process counters, public transport offsets,
file metadata and SQLite compilation completion bits, never private payloads.

The existing protocol has a concrete cost at this cohort size. In
`R/dp_gaussian_cross_transport.R:338–355`, validity has n entries, masked values
have n*m and moments have n*m*(m+1)/2. The fresh stage sequence at lines 760–790
requires m-1 validity stages. Server `R/dpGaussianCrossDS.R:1307–1318` fixes
Ring128; client `R/exact_gc_transport.R:457–464` caps direct multiplication chunks
at 64 (the separate ring127-OT alternative does not apply). Chunking at lines
973–1021 therefore requires
`3*ceil(64/64) + ceil(256/64) + ceil(640/64) = 17` multiplication operations.
This excludes alignment, R preparation/finalization and the final DP sampler.
The observed per-operation transport/coordination cost makes the 900-second
functional gate impractical here. This is a bounded-runtime result, not evidence
of a new security or correctness defect. No chunk limit or protocol guard was
changed, and neither reviewed production fix (`a957113`, `f6f123c`) was modified.

No smaller-cohort retry was made. The independently requested Poisson gate is
separate; optional K5 is not attempted after this K3 limit. The prepared K5
fixture has three data owners and two additional signed federation members,
each with an unrequested bounded local marker because empty local schemas are
rejected. An offline check of the unchanged server scope function confirms those
markers do not enter the requested block or its 29-coordinate calibration. This
planning is not K5 live evidence.

### Final gate outcomes and accounting

| Gate | Outcome | Measured seconds |
|---|---|---:|
| K3 Cor, n64 / four variables / mixed pairs | Incomplete at bounded command deadline; exit 124 | 897.1367984600365, entire command including cleanup |
| K3 PCA | Not reached | — |
| K2 Poisson LASSO via public `dp_lasso_grid` | PASS, including current certificate validation against independently obtained peer pins | 273.5162, analysis |
| K2 Poisson full command | PASS, exit 0, complete worker cleanup | 357.93558445200324 |
| Optional K5 Cor/PCA | Not attempted after the K3 stop condition | — |

Poisson setup/PSI took 54.9176 seconds and authenticated registration 26.9793
seconds. Registration plus the verified fit took 300.5216 seconds. Its certificate
passed integrity and anchored-authenticity checks; implementation state is
`cross_owner_exact_gc_materialized`, with two penalty paths and no standard
errors. The diagnostic observer retained `dsvert_phase_not_ready:v1` from normal
publication polling at both peers; it did not prevent completion. No fixture
seams, substituted samplers or relaxed checks were used. The actual public fit is
retained as `.tmp/vertcor/vertcor04-k2-poisson-lasso.rds` (37,766 bytes), SHA256
`dbfb8688654195adfa0a5771c8fbcee34eafa0eeed767feebd6471b55db2b84d`.

The two commands used exactly **1,255.0723829120398 seconds (20m55.072383s)** of
new live-command time, including cleanup, below 2,700 seconds. Their elapsed
campaign span was 1,270.874165572226 seconds (21m10.874166s), including the
15.80178266018629-second gap between commands. Each command was below 900 seconds.
Both completion records report `cleanup_complete=true` and zero remaining owned
PIDs; no R or native workers remained after the second command. These durations
are separate from the earlier vertcor03 budget. No new production defect was
found or changed; server source, runtime bytes and both fixes under independent
review remain unchanged. K3 remains a required live release gate; optional K5
has no new live evidence.

The exact durable records and complete short command logs are committed here:
`vertcor04-smoke-ledger.jsonl`, `vertcor04-analysis-events.jsonl`,
`vertcor04-k3-cor-pca-001.log`, `vertcor04-k2-poisson-lasso-002.log`.
Commands, from the client root after sourcing the environment, were:

```sh
source /workspace/dsflower/.tmp/vertcor/env.sh
DSVERT_SMOKE_CASES=cor3 DSVERT_SMOKE_ROWS=64 \
  python3 -B tools/run_cross_owner_smoke.py --root /workspace/dsflower \
  --label k3-cor-pca -- Rscript tools/validate_cross_owner_142.R /workspace/dsflower
DSVERT_SMOKE_CASES=poisson \
  python3 -B tools/run_cross_owner_smoke.py --root /workspace/dsflower \
  --label k2-poisson-lasso -- Rscript tools/validate_cross_owner_142.R /workspace/dsflower
```

The launcher defaults to campaign `vertcor04`, a 900-second command cap and a
2,700-second cumulative cap. Do not erase or reuse a campaign ledger to reset
its accounting. Targeted checks for the new harness files were:

```sh
timeout 120 Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-cross-owner-smoke-fixture.R", reporter="summary")'
python3 -B tools/test_run_cross_owner_smoke.py
```

Results: 17 R assertions in three blocks and five Python test cases passed. An
initial test-only comparison of a `table` with an unclassed array was corrected;
the final test run is the one included in `targeted-summary.csv`. No production
tests were repeated because this follow-up changes validation tools and records
only. No full suite or push was run.

## 10. Extended K3 run (maintainer-approved, one-off)

The maintainer approved one live K3 Cor/PCA run without the 900-second command cap. The fixture is unchanged from
section 9:
- n = unit capacity = 64;
- four variables `site_a$a`, `site_a$b`, `site_b$x`, `site_c$x`, giving one same-owner pair plus cross-owner pairs
  across all three data owners;
- epsilon 8, delta 0, 8-bit quantization and add/remove-patient adjacency;
- real PSI, signatures, allowlist and production release endpoints.

Only the time limits changed. The tracked launcher and runner were copied outside the repository, with exactly these
edits:
- `tools/run_cross_owner_smoke.py`: `MAX_COMMAND_SECONDS` 900 → 5400 and `MAX_CAMPAIGN_SECONDS` 2700 → 5400, with the
  two matching range messages.
- `tools/validate_cross_owner_142.R`: the per-analysis guard `< 2680` → `< 5300`, and `prior_seconds < 2700` →
  `< 5400`.

No package code, sensitivity, epsilon, sampler, chunk limit, attestation or identity rule changed.

```sh
source /workspace/dsflower/.tmp/vertcor/env.sh
DSVERT_SMOKE_CASES=cor3 DSVERT_SMOKE_ROWS=64 python3 -B <copy>/run_long.py --root /workspace/dsflower \
  --campaign k3long --label k3-cor-pca-long --timeout 5400 --budget 5400 -- \
  Rscript <copy>/validate_long.R /workspace/dsflower
```

| Gate | Outcome | Measured seconds |
|---|---|---:|
| K3 setup and real PSI | PASS | 71.5925 |
| K3 Cor, n64 / four variables / mixed pairs | **PASS** | 3,378.669 |
| K3 PCA (loadings/eigenvalues/variance, no scores) | **PASS** | 2.6024 |
| Whole command, including cleanup | exit 0, `cleanup_complete=true`, no remaining owned PIDs | 3,455.488 |

**The PCA result records:**
- `cross_owner_state = "exact_gc_to_joint_dp_vector_v1"`;
- `provenance_integrity = TRUE` and `provenance_authenticity = "session_transport_anchored"`;
- no scores, and no source or intermediate values exposed;
- zero additional privacy cost and zero server calls after the Synopsis;
- no legacy exact route called.

The eigenvalues are 1.746, 1.185, 0.819 and 0.250 (sum 4).

**This is a functional gate, not a precision result.** The fixture's population correlations are exactly zero, while
the released correlations range over |r| ≤ 0.655. That is consistent with the planned pre-clamp noise SD of 7.25 per
normalized coordinate against centered sums of squares of 16 (section 9).

**Runtime** is dominated by the 17 serial exact-GC multiplication operations that the cross Gaussian path needs at this
shape under the 64-element chunk cap: about 3 min each in this environment.

Outputs (local):
- `k3long-k3-cor.rds`, SHA256 `dd99a0578e3016c93b330d32d3f016dce7c941c48e73fcb74633621fc789b3b8`;
- `k3long-k3-pca.rds`, SHA256 `f23f281c394777d06dcacc8359deddbe3bfa8600405a3ffc8360085edf86dac7`.

The committed records are `k3long-smoke-ledger.jsonl`, `k3long-analysis-events.jsonl` and
`k3long-k3-cor-pca-long-001.log`. The synthetic last-endpoint `dsvert_phase_not_ready:v1` lines are the expected
publication polling seen in every passing run.

**The K3 live gate is closed.** Optional K5 was not attempted.
