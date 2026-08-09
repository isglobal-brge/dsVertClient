# Standard per-analysis DP and MPC plan

Status: proposed replacement architecture, 2026-08-09.

This document defines the intended product architecture for `dsVert` and
`dsVertClient`. It replaces the bounded-lifetime capsule/accountant model. It
does not describe the current shipped ABI until the implementation gates below have
passed.

## Decision

The product will use:

1. standard differential privacy for each distinct analytical artifact;
2. MPC or exact GC for confidential intermediate computation;
3. one globally sticky release for each canonical analysis;
4. unlimited calls, aliases and methods, with no per-user quota, request rate
   limit or privacy catalogue; and
5. explicit statistical-utility and performance gates before a method is
   promoted.

All DataSHIELD analysts may share one server principal. The design therefore
must not attribute a privacy budget to a login, connection or session.

The guarantee is deliberately scoped. Every new artifact has its declared
per-analysis `(epsilon, delta)` guarantee. Equivalent calls replay the same
artifact at zero additional privacy cost. An unlimited sequence of distinct,
overlapping analyses has no finite global composition claim and no absolute
non-reconstruction theorem. The status API and documentation must say this
directly.

One semantic analysis is the only privacy-accounting boundary. If its mechanism
contains several DP coordinates or internal rounds, their epsilon/delta shares
are calibrated and composed in memory inside that operation and are committed
in its semantic contract. They never update persistent state. The preferred
design keeps all rounds inside MPC and makes one final DP opening, avoiding
internal composition where possible. There is no accumulated budget per user,
dataset, resource, method, server or deployment.

## Product goals

- Preserve row values and all iteration-level statistics inside the sites or
  MPC workers.
- Return results that are numerically close to the corresponding centralized
  bounded estimand.
- Make the common path fast: one source preparation, one internal model job
  and one final opening. Equivalent calls use identical randomness even when
  they are recomputed.
- Make the client plug-and-play: ordinary formulas and method arguments, no
  ring selection, no manual release identifiers and no protocol phases.
- Support K=2 and the generic K>=3 source-owner topology. Exactly two pinned
  compute/noise peers remain the confidential-computation boundary.
- Reuse shared statistical engines so adding a method does not add a new
  disclosure protocol.

## Non-goals and threat boundary

- No lifetime DP claim across arbitrary distinct analyses.
- No query-count, rate-limit or catalogue enforcement.
- No malicious-host, all-peer-collusion or simultaneous total-state-rollback
  claim.
- No physical constant-time theorem for the host, scheduler, HTTP stack or
  WAN.
- No claim that arbitrary R objects are valid analytical datasets. The initial
  contract is a tabular privacy-unit dataset with supported atomic types,
  explicit missingness rules and custodian-owned bounds.

The supported threat model retains pinned identities, authenticated TLS,
private PSI and two semi-honest designated compute peers, with at least one of
those peers not colluding with the analyst.

## User-facing contract

The normal call remains a statistical call, for example:

```r
ds.vertGLM(y ~ age + treatment, family = binomial(), data = "D")
ds.vertCox(Surv(time, status) ~ age + treatment, data = "D")
ds.vertGEE(y ~ treatment + time, id = "patient", data = "D")
```

The user does not supply:

- an MPC ring;
- a noise seed;
- a sensitivity;
- a peer role;
- a release index; or
- a security-critical `analysis_id`.

The result should look like the corresponding familiar R fit and add a compact
`privacy` component containing the canonical artifact ID, per-analysis
epsilon/delta, bounds, mechanism, clipping diagnostics, K/pinset binding and
sticky-replay status.

Server policy owns privacy parameters and contribution bounds. A friendly
client option may request a stricter policy, but cannot weaken the server
policy or select an unsafe backend.

There is no privacy allowlist of operations. Every method implemented by the
installed package is callable with its supported dynamic formula and arguments.
Kernel names version executable semantics; they do not reserve, authorize or
charge a method. An unavailable method means that its implementation has not
yet passed its correctness and privacy gates, not that a catalogue or budget
denied the call.

## Canonical analysis artifact

Every one of the K custodians reconstructs and validates the same semantic
request. The sites sign opaque, site-keyed HMAC commitments to their immutable
snapshots; raw dataset or cohort digests are never published. The release
requires consensus over those commitments and normalized statistical semantics
before any result is opened.

The `semantic_release_contract` and its artifact key contain at least:

- stable owner-site identities and opaque snapshot/aligned-cohort commitments;
- privacy-unit and adjacency semantics;
- canonical formula AST and intercept;
- effective terms, contrasts, factor references and column order;
- offsets, weights and missingness policy;
- method family, link, tuning and convergence policy;
- custodian bounds, clipping and maximum contribution;
- the two stable Ed25519 identities that contribute private noise;
- scale, grid, rounding, overflow, sampler and output encoding;
- the complete server-compiled plan of stochastic lanes used by the operation;
- public result shape;
- statistical mechanism and semantic version; and
- server-owned epsilon and delta.

The initial contract derives the two noise authorities as the first two
canonical stable identity keys in the full owner pinset. A request cannot pick
another pair to obtain a fresh draw. Its randomness schema requires one audited
final-noise lane and may contain compiler-declared confidential internal lanes;
all are committed by the artifact key. This versions DP mechanisms, not a
catalogue of statistical operations.

A separate signed `execution_contract` contains build hashes, connection-name
to stable-identity mappings, chunk geometry, transport parameters and current
TLS details. Build, chunk and certificate details can change without rerolling
noise only when conformance tests prove that the public bytes remain bit-exact.
The physical ring is always execution metadata. If it cannot satisfy the
committed value domain and bit-exact numeric semantics, execution fails closed;
it does not mint another artifact for the same analysis. The two logical noise
authorities and identity seeds are fixed for one deployment. Replacing or
reassigning one starts a new deployment without a continuity claim. If an
implementation change alters the estimand, numeric result or DP mechanism, it
requires an explicit semantic-version change instead.

The key excludes the analyst login, session, connection, API alias, argument
ordering, formatting and pure post-processing choices. Semantically equivalent
aliases must resolve to the same primitive artifact. For example, one DP 2x2
table can serve epidemiology, diagnostic, chi-square and Fisher
post-processing without another release.

Unlike the current proposal object, method and effective arguments are part of
the artifact identity. No server accepts a caller-provided artifact key as
authoritative.

A public preflight can estimate shape, sensitivity and cost from policy and
declared schema. Exact factor levels, missingness and output shape are
canonicalized inside protected processing before the semantic contract is
signed and before sampling. Assignment/alignment creates an immutable logical
snapshot commitment that excludes the R symbol name; a changed binding is a
new snapshot. A service may cache the validated commitment in RAM, but it must
not rescan the full dataset on every call.

## Sticky randomness without a release database

Each designated peer already keeps one fixed owner-only `identity.seed`. It is
the only persistent secret. Domain-separated HMAC derivations create an
internal sticky-noise subkey and a separate snapshot-commitment subkey; the raw
master is never reused directly as noise. The sticky subkey and canonical
artifact key derive labelled pseudorandom streams for every stochastic choice.
The two private peer contributions are combined inside MPC to produce one valid
draw for the target standard DP mechanism.

The same semantic key reproduces the same random stream after a restart; a
different key obtains computationally independent material. The new route has
no SQLite release database, privacy ledger, allocator, reservation log, result
memoization or external storage dependency. Repeating a call may recompute the
model, but it cannot obtain independent noise.

Deterministic recomputation requires fixed-point computation, ordering,
serialization and signing; timestamps, attempt IDs and transport nonces are
excluded from the replayable artifact. Losing `identity.seed` requires
restoring the same bytes from authenticated backup or failing closed. The seed
never rotates inside one deployment and the service never silently generates a
replacement. Changing execution builds or TLS certificates must not alter the
derivation; changing a logical noise authority starts a new deployment without
continuity.

This is computational DP: the standard DP mechanism is driven by secret PRF
outputs that are computationally indistinguishable from independent random
coins. Seed custody, separation between noise/signing domains, owner-only
permissions and authenticated backup are therefore part of the threat model.

No state field records or uses the number of earlier analyses to deny, alter or
degrade a new distinct analysis.

Package `configure`, installation and `.onLoad()` never create the seed. The
first real service bootstrap runs only after the node's private persistent
volume is mounted, creates 32 bytes from the operating-system CSPRNG if the
file is absent, and commits them atomically with owner-only permissions. An
image must not contain a configured literal seed; different nodes use
different state volumes and therefore different identities.

## Confidential computation engine

The productive model backend will be one native MPC worker for the complete
analysis job, shared by method families. It need not be a cross-analysis daemon
unless benchmarks later prove that extra lifecycle complexity worthwhile. It
must:

- ingest bounded, privacy-unit-collapsed source shares once;
- use vectorized Ring127/Ring128 arithmetic and batched OT/Beaver operations;
- reserve exact GC for comparisons, clipping, selected products and other
  operations that actually require it;
- execute all model iterations inside the worker rather than through one DSI
  round trip per iteration;
- expose no per-site gradient, score, Hessian, risk set, cluster contribution,
  random effect, residual or individual prediction;
- return only the final privatized result and signed diagnostics; and
- retain no plaintext protected source outside the owner site.

Ring and scaling selection are implementation details. Numeric preflight must
choose the safe backend automatically and reject overflow or unsupported shape
before protected computation.

The current row-wise formal-GLM exact-GC pipeline is not the productive
backend: its work and wire volume scale poorly. Its mathematical oracles,
bounds and adversarial tests may be extracted into the new worker before the
phase-specific lifecycle is removed.

## Statistical mechanisms

There is no requirement that every family use the same DP construction.
Mechanisms should be standard, method-appropriate and selected for utility.

| Family | Intended confidential statistic | Intended release |
|---|---|---|
| Counts, tables, histograms | bounded per-unit cell contributions | discrete Laplace or Gaussian vector, once |
| Means, moments, correlation | clipped sums/products in MPC | noisy sufficient statistics, then post-processing |
| Gaussian regression/PCA/LASSO | bounded Gram and cross-product | one privatized moment artifact |
| Binomial/Poisson/NB | clipped score/information or strongly-convex objective | private final fit and covariance |
| Multinomial/ordinal | bounded class/threshold scores and information | private final fit and covariance |
| Cox | bounded risk-set score and information | private coefficients, covariance and safe diagnostics |
| GEE | bounded per-cluster bread and meat contributions | private coefficients and sandwich covariance |
| LMM/GLMM | bounded per-cluster likelihood/score contributions | private fixed effects, variance parameters and covariance |
| IPW/AIPW | propensity and outcome stages inside one MPC job, bounded weights | private target estimand and uncertainty |
| MI | fixed, non-rerollable imputations inside MPC | private pooled result only |

Where a direct sensitivity proof would be fragile, a standard
sample-and-aggregate construction may provide a generic first backend. High-use
families should then receive specialized bounded-score or objective-perturbation
backends when those materially improve speed or precision.

This table is a roadmap, not a privacy proof. Before promotion, every family
must specify adjacency, clipping, the sensitivity and mechanism for every
released coordinate and diagnostic, and any composition internal to that one
artifact. Fixed MI draws, confidential Cox fitting or MPC bread/meat alone do
not make the final output DP.

Exact public output is appropriate only for public constants or deterministic
post-processing of an already private artifact. Confidential model estimates
are not described as exact DP outputs even when their internal MPC calculation
matches the centralized fit.

## Utility policy

Before protected source access, a public preflight estimates:

- artifact family and a declared or bounded output dimension;
- L1/L2 sensitivity and clipping/contribution policy;
- candidate mechanisms and selected backend;
- certified mechanism radius in natural units;
- projected chunks, bytes and protocol phases; and
- the benchmark class used for the latency estimate.

Protected canonicalization then fixes the exact factor levels, missingness,
shape and opaque snapshot commitments before sampling. No data-dependent
detail discovered in that step is returned unless covered by the final DP
artifact.

Conditioning, event rarity and cluster adequacy may be unknown before the one
release. The result therefore includes private or post-processed diagnostics
without opening exact confidential counts.

Initial product acceptance targets, to be calibrated with committed evidence,
are:

- DP-excess RMSE no more than 0.25 times the non-private sampling SE in common
  biomedical scenarios;
- 0.25--0.5 times the sampling SE produces a visible utility warning;
- above 0.5 is labelled low utility rather than silently presented as precise;
- advertised interval coverage is validated with a predeclared Monte Carlo
  precision; and
- non-estimable, rare-event and ill-conditioned cases fail or return an honest
  typed result.

These are promotion targets, not claims about the current package.

## Performance policy

Performance is a release gate, not deferred cleanup.

- Use the minimum artifact needed by the analysis; remove universal
  `all_schema` expansion from the normal path.
- Batch vector operations and comparisons.
- Reuse site-local canonical preprocessing and safe in-memory shares while one
  analysis is active.
- Keep iteration inside one persistent worker.
- Derive one logical noise draw for each semantic key, even when a restart
  causes deterministic recomputation.
- Make concurrent equivalent calls derive identical randomness and results.
- Select one-draw versus independent-full/convolution backends using measured
  certified radius and p95 latency, not a hard-coded dimension guess.

The benchmark matrix covers at least K={2,3,5}, representative WAN/LAN and
local deployments, n={1k,10k,100k,1M}, method-appropriate p, rare events,
imbalance, missingness, collinearity and cluster-size distributions. It records
fresh/recomputed p50 and p95, CPU, peak RSS, temporary disk, wire bytes and DSI
phases.

Initial performance targets are:

- a restart may repeat protected computation but reproduces the same canonical
  result bytes and does not obtain fresh randomness;
- DP release overhead is no greater than 1.5x the same secure-MPC estimator
  without noise for common profiles; and
- no normal query performs O(history) validation because there is no release
  history or global accountant.

MI and large mixed models may publish explicit, evidence-backed exceptions.

## K topology

Every statistical family must pass K=2 and K=3 against the same centralized
estimand. Each shared backend must additionally pass K=5 to exercise generic
fanout and consensus. K=4 needs no separate implementation branch.

K counts source owners. It does not change the two-designated-peer MPC/noise
boundary. More owners must not receive a compute share or broaden the public
result.

## Cleanup map

### Keep

- identity, TLS pinning and signed peer sets;
- private padded PSI and authenticated aligned snapshots;
- typed streaming, adaptive DSI windows and bounded spools;
- Ring127/Ring128, IKNP/OT, Beaver and fixed-point kernels;
- exact-GC transport and guards where needed;
- custodian bounds, clipping and contribution policy;
- standard Laplace/Gaussian samplers and their distribution tests;
- canonical JSON, signed manifests and tamper checks;
- deterministic final-artifact/recomputation semantics; and
- current deterministic post-processing methods.

### Replace or simplify

- lifetime allocator/reservation -> remove; deterministic calls need no
  privacy charge or durable query claim;
- capsule registry and lifetime release ledger -> remove without replacement;
- privacy-accountant namespace and noise-root subsystem -> remove; derive
  domain-separated subkeys from the existing fixed `identity.seed`;
- universal capsule catalogue -> minimal automatic artifact contract;
- manual `analysis_id` -> optional display label, never authority;
- O(history) audit -> remove with the historical accountant;
- per-phase GLM/Cox lifecycle -> generic persistent model job.

### Remove after replacement is green

- `lifetime_max_distinct_capsules`, cumulative lifetime epsilon/delta and
  exhaustion errors;
- budget admission, reservation burns and remaining-unit telemetry;
- `capsuleRegistry` and allocator-registry integration;
- all SQLite/DBI code and both package dependencies. Any surviving transport
  spool uses bounded ephemeral files or memory, never a database;
- lifetime-only status fields, docs and tests;
- universal `catalog_v1/all_schema` normal operation;
- superseded formal GLM/Cox phase lifecycle and unreachable commands after
  useful kernels/tests are extracted;
- stale generated documentation and historical claims that imply current
  method availability; and
- superseded remote endpoints that expose iterative or site-level aggregates.

The large bounded-capsule commits must not be reverted wholesale because they
also contain valuable PSI, transport, exact-GC, numeric-policy, sampler and
biomedical-method work. Cleanup is dependency-aware and staged.

## Method implementation order

1. Existing DP summaries, categorical/epidemiology, Gaussian and survival
   artifacts move to the canonical per-analysis key and deterministic PRF.
2. Binomial and Poisson GLM, then Cox, use the common persistent model worker.
3. NB2, multinomial, ordinal and iterative LASSO reuse the score/information
   engine.
4. IPW/AIPW reuse the private binomial and outcome stages.
5. GEE introduces the private cluster bread/meat engine.
6. LMM and GLMM reuse the cluster engine with bounded random-effect models.
7. MI uses fixed internal draws and pooled-only release.
8. New methods such as log-rank, Fine--Gray and TMLE/AIPW are added only when
   they reuse a safe artifact or backend rather than adding granular remote
   endpoints.

Frontdoors remain fail-closed until their complete method gate passes. They are
not made available by deleting a quarantine check.

## Acceptance gates per method family

Security and replay:

- 10,000 equivalent calls, aliases, argument reorderings and concurrent
  sessions produce one semantic artifact ID, one randomness commitment and
  identical result bytes;
- restart and lost ACK may rerun source, sampler and MPC but must reproduce the
  same randomness commitment and canonical result bytes;
- loss, substitution or path drift of the persistent identity seed fails
  closed; immutable snapshot commitments prevent a silent in-place rewrite;
- no repeat can obtain an independent noise draw for the same semantic key;
- the new route creates no SQLite database and has no history-dependent gate;
  and
- result schemas contain no site-level or iteration-level intermediates.

Mathematics and statistics:

- zero-noise MPC matches the centralized bounded/clipped estimand within its
  numeric certificate;
- adjacent-dataset and distribution tests cover mechanism extremes;
- simulation reports bias, RMSE, DP-excess RMSE, coverage, type-I error, power,
  convergence and non-estimable rate; and
- factors, missingness, dates, offsets, weights, imbalance, rare events,
  collinearity and cluster edge cases are exercised where applicable.

Runtime and product:

- source and installed-package suites pass;
- K=2, K=3 and backend K=5 pass artifact-exact;
- DSLite and Opal benchmarks meet or explicitly document the targets;
- server/client ABI and registered surface are exact;
- no ring or protocol option is required from the user; and
- examples return familiar R objects with clear privacy and utility metadata.

## Implementation milestones

Each milestone is a separate reviewed commit pushed to `main` in both package
repositories where applicable.

1. Approve this architecture and freeze an evidence baseline.
2. Add the new semantic/execution contracts and deterministic sticky derivation
   without SQLite.
3. Migrate existing DP methods and prove canonical replay.
4. Remove lifetime admission, registry, accountant ABI and obsolete docs.
5. Optimize one-draw selection and publish utility/performance evidence.
6. Promote GLM and Cox through the new worker.
7. Promote the remaining model families in the order above.
8. Run complete source/build/check, K=2/K=3/K=5 and artifact-exact Opal gates.
9. Remove superseded formal code only after its replacement is proven.

This is a pre-release clean break. There is no persisted DP state to migrate or
replay. The release removes the old lifetime ABI and SQLite code rather than
shipping a transitional mode. A final source, dependency and artifact scan must
prove that the released packages neither create nor consult a DP SQLite
database.
