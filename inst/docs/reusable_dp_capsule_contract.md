# Reusable disclosure-safe analysis capsule

Status: normative migration contract, 2026-08-01.

## User-facing availability contract

For every method and argument combination advertised as supported for a
registered immutable snapshot:

- the number, order, or history of operations over the same authenticated
  capsule/release instance MUST NOT consume another privacy unit, advance the
  allocator, reduce accuracy, or deny replay;
- there MUST NOT be a request counter, rate limit or `queries_remaining` value;
  a previously unseen capsule is instead subject to the finite lifetime
  privacy-loss gate;
- retrying an operation, changing harmless syntax, or using a different method
  over the same capsule MUST reuse the same immutable private release; and
- mathematical non-identifiability, an invalid scientific specification, or a
  genuine infrastructure failure MUST be distinguished from privacy accounting
  and MUST NOT be reported as a query-limit failure.

The protected event is capsule creation, not an analysis operation. A capsule
is a single differentially private representation of one immutable logical
snapshot and one complete, versioned analysis-workload contract. All subsequent
answers are post-processing of that representation and therefore add no
privacy loss, regardless of how many are requested.

## Honest privacy accounting

Each capsule has fixed `capsule_epsilon` and `capsule_delta` parameters. They do
not decay with an allocation index. The custodian-owned
`lifetime_max_distinct_capsules = N` defaults to one. Exact decimal arithmetic
MUST prove `N * capsule_epsilon <= 8` and `N * capsule_delta < 1` before the
policy is accepted. This bound has status-v5 scope
`at_most_N_immutable_snapshot_workload_capsules_per_stable_privacy_accountant_namespace`.
The allocator commit consumes one persistent unit before
protected-source access or sampling; an abandoned commit is not refunded.
Authenticated durable state prevents rerolls, forked publications and duplicate
charging of the same capsule.

Publishing different capsules for overlapping cohorts or successive snapshots
does compose. The reservation registry MUST count every committed burn, whether
or not publication follows, and MUST deny a new capsule when `N` is reached.
Publication is recorded separately and MUST never exceed the reservation count.
This history gate does not reduce the fixed accuracy of an admitted capsule.
Reusing one exact capsule/release instance has zero additional privacy cost.

The fixed `[dsvert_dp_lifetime_budget_exhausted:v1]` condition is deliberately
an opaque terminal union: either global `N` is exhausted, or the requested
capsule cannot safely advance because its irrevocable release-instance claim or
sole publication slot is already bound and the exact instance cannot continue
or replay. It does not imply `remaining_distinct_capsules == 0`; separate public
causes would leak state and break the compatible failure boundary.

Joint-DP capsule status v5 therefore reports `capsule_epsilon`, `capsule_delta`,
`lifetime_max_distinct_capsules`, `capsules_created`,
`releases_published`, `remaining_distinct_capsules` and exact bounded
composition. It binds
`operation_accounting = "one_per_distinct_capsule_allocator_commit"`,
`operation_limit = true`, `request_limit = false`, and
`history_can_deny_operation = true`. “Remaining” denotes possible new capsules,
not end-user requests. It exposes no geometric `decay` or per-query epsilon
schedule.

Status v5 and result security claim v3 also bind the assumption that at least
one non-colluding designated peer retains and uses complete authenticated
monotonic history. The exact public value is
`authenticated_history_retention_assumption=at_least_one_noncolluding_designated_noise_peer_retains_and_uses_complete_authenticated_monotonic_history`.
The result claim additionally binds
`privacy_accountant_namespace_assumption=one_stable_unique_namespace_across_domain_cohort_policy_pinset_and_ledger_reconfiguration_per_protected_privacy_universe`.
Today this is a custodial deployment assumption: neither package discovers or
enforces a unique namespace globally, nor automatically migrates reservation
history across reconfiguration. They make no claim against simultaneous
rollback of both designated histories without an external linearizable CAS.
This is not protection against a compromised host or collusion of both
designated peers.

## Frozen vector ABI and rollout

The pre-release baseline is:

| Artifact | Exact version |
|---|---|
| Joint-DP control protocol | `dsvert-joint-dp-control-v3` |
| Canonical capsule identity | `dsvert-joint-dp-capsule-identity-v3` |
| PREPARE signed receipt | `dsvert-joint-dp-vector-prepare-v6` |
| START signed receipt | `dsvert-joint-dp-vector-chunk-start-v5` |
| RESULT signed receipt | `dsvert-joint-dp-vector-local-result-v5` |
| RELEASE signed receipt | `dsvert-joint-dp-vector-release-root-v5` |
| ACK signed receipt | `dsvert-joint-dp-vector-finalization-ack-v5` |
| Authenticated durable STORE schema | `dsvert-joint-dp-vector-store-v6` |
| Per-capsule instance claim | `dsvert-joint-dp-vector-instance-claim-v1` |
| Authenticated instance-claim set | `dsvert-joint-dp-vector-instance-claim-state-v1` |
| Replay response | `dsvert-joint-dp-vector-replay-v4` |

Every signed receipt and replay response MUST attest exactly
`history_gate=TRUE`, `request_limit=FALSE`, and `operation_limit=TRUE`.
STORE rows are separately HMAC-authenticated durable state.

Control v3 and capsule identity v3 bind the lifetime fields and biomedical
workload v7. STORE v6 makes both claim-v1 artifacts required durable
invariants. Their v2 artifacts and legacy v4/v5 stores receive no automatic
migration or re-signing. Rollout MUST use empty DP capsule state or a future
audited offline migration; any legacy artifact fails closed. Previous Opal
deployments had `POLICY_READY=FALSE`,
published no DP capsules and used ephemeral local K-site state, so this baseline
does not abandon a production publication history.

Exact COMMIT/RELEASE endpoint replay and sticky replay within the live session
are O(1). This does not imply O(1) cold replay after a process restart: the
end-to-end reconstruction re-enters through AllocationProof, which performs a
complete O(N) audit of the allocator journal before returning the proof. The
audit consumes no lifetime unit, invokes no sampler and reads no protected
source.

## Capsule identity and bindings

`capsule_id` is derived from a canonical, versioned contract containing:

- consortium, cohort, logical snapshot identifier and version;
- privacy epoch and noise-root identifier;
- complete pinned-peer map and designated computation peers;
- privacy unit, adjacency, contribution bounds, clipping, missingness and
  repeated-record rules;
- public categorical domains, numeric bounds and time grids;
- the full capsule workload/schema and every mechanism/circuit version; and
- data-independent capacity, padding, numeric-grid and result-contract rules.

It MUST NOT contain a public digest derived from patient identifiers or raw
values. Each custodian privately binds the public logical snapshot descriptor
to its local immutable data and ordered PSI/alignment attestation. Cross-signed
peer receipts bind those private checks to the same public capsule contract.

During the compatibility migration, an internal wire field named `query_id`
may carry exactly the same value as `capsule_id`. It MUST be validated as an
alias, MUST NOT include method names or arguments, and MUST NOT be interpreted
as an end-user query allocation.

### Server-authoritative manifest bootstrap

The relay has no API parameter for a dataset map, owner, bound, category,
version, snapshot or workload specification. Bootstrap is fixed to three
idempotent phases:

1. Each peer signs a data-free draft containing only its pin, public dataset
   identity/version, patient-key column name, numeric bounds, categorical
   domains, alignment protocol version and custodian-configured workload
   fragments. Snapshot/alignment hashes, row counts and patient-derived
   digests are forbidden.
2. The client verifies every pinned signature, rejects duplicate ownership or
   analysis IDs, and deterministically merges the drafts. Each workload entry
   retains the peer that signed it. The complete canonical workload contract
   is included in the logical-snapshot fingerprint.
3. Every peer compares its owned contract fragment byte-for-byte with fresh
   local policy, validates the complete schema/workload, and signs the same
   schema. Only after all pinned signatures are present may every peer build
   and durably memoize the manifest. The client accepts it only if all peers
   return identical manifest bytes, hash and capsule ID.

The source and sampler boundaries accept a manifest only when that exact byte
string exists in the server's authenticated local authority cache. Repeating
Draft, Sign or Build is unlimited and idempotent; it does not create a query
allocation or inspect request history. The cache also binds each public capsule
version to one HMAC-authenticated private snapshot/alignment policy. Changing a
private binding while reusing the same public dataset version is a typed
versioning conflict, not a privacy-budget denial; the custodian must publish a
new dataset version so changed data can never reuse the old capsule/noise.

## One-time creation protocol

Capsule creation is single-flight, crash-recoverable, and sticky:

1. Every peer validates the data-free immutable policy, bounds, admission
   contract, pinset, epoch and workload manifest.
2. A deterministic leader proposes one capsule sequence entry. At allocator
   commit, each designated peer durably burns the same lifetime unit and writes
   its dense reservation-journal row. Competing proposals cannot reserve
   divergent heads; no later abort refunds the unit.
3. Only after that reservation do peers resolve the protected snapshot and
   alignment and compute vertical measurements with purpose-bound MPC. Raw rows,
   exact local statistics, individual noisy shares, guard bits, pre-clamp values
   and data-dependent hashes are never opened to the relay.
4. Vector PREPARE persists an authenticated candidate without claiming the
   capsule, so sibling candidates may coexist. At each designated peer the first
   valid START atomically and irrevocably persists the per-capsule claim before
   local staged-source access or sampling. Source transport may already have
   staged encrypted protected material, but no noised share or public output
   exists at this boundary. Split sibling claims cannot form the matching
   bilateral receipts required for publication under the non-collusion model.
5. The two designated pinned peers contribute independent, full-parameter
   cryptographic noise inside the joint finalizer. Exact signed decoding,
   projection and the sole final clamp occur before the public opening.
6. The canonical capsule payload, mechanism certificate and cross-signed
   receipt are atomically committed before they can be served.
7. Lost acknowledgements, process crashes and restarts return the byte-identical
   committed capsule while its authenticated durable release remains available.
   They never resample merely because a request is retried.

Each designated peer also maintains an automatically generated, public,
data-independent release domain outside the package image. It is part of the
signed release-instance identity but never of the immutable `capsule_id`.
Rotation may select another candidate only before the first valid START claim.
Once claimed, that capsule has its sole allowed instance for life: missing state
must be authenticated and restored so the same instance can continue, or the
server fails closed. Once publication exists, only its exact bytes may be
replayed. It MUST NOT generate fresh noise or a second instance for
availability. Interrupted materialization may continue/retry its claimed
instance without creating a public differencing opportunity.

If every copy of release and reservation history is destroyed or rolled back
together, local software cannot prove continuity. That simultaneous-history
case is outside the package model unless an external linearizable anchor covers
it; missing history must never be presented as fresh allowance inside the
stated model.

If the capsule has not yet been materialised, the first supported operation may
trigger this protocol automatically and wait for the single shared build. It
must not require a server administrator to provision a seed manually.

## Analysis isolation

Once committed, a public analysis endpoint may read only:

- the immutable capsule payload and its public schema;
- public method arguments; and
- deterministic or explicitly public randomness used only for post-processing.

It MUST NOT retain a handle to the protected frame, invoke a legacy exact
statistic, or request fresh data-dependent randomness. This separation is the
enforcement boundary that makes unlimited analysis compatible with DP. Public
methods should preferably run client-side; large capsule slices may be fetched
on demand, but the server-side projection must be a pure function of the same
capsule.

## Statistical representation

A single naive synthetic table is not a sufficient biomedical inference
contract. The capsule is a hybrid, workload-aware representation:

- fixed-domain one-way and selected joint marginal measurements for counts,
  contingency, epidemiology, diagnostics and standardisation;
- bounded counts, sums, squared sums and a jointly noised cross-product matrix
  for descriptives, correlation, PCA and Gaussian regression;
- fixed public time-grid entry, risk, event and competing-risk measurements for
  survival summaries;
- mechanism covariance, quantisation/clipping bounds and consistency
  projections needed to propagate privacy uncertainty; and
- multiple noise-aware synthetic/posterior representations for supported
  nonlinear models, accompanied by method-specific validation rather than
  ordinary inference that pretends synthetic records are original observations.

The workload is complete for the published method registry and declared
columns at capsule creation. Data-dependent workload selection, if used, is
part of the one-time DP mechanism. A user operation never extends the workload
by querying the source data. New columns, bounds, domains, mechanisms or
snapshot contents create a new capsule version; they do not mutate an existing
one.

Method adapters MUST preserve historical high-level names where the scientific
estimand remains defensible. A method that cannot yet obtain valid inference
from the capsule may remain explicitly provisional while its adapter is built;
it MUST NOT silently call a legacy exact/disclosive server primitive. Removal
of a pre-existing method requires an explicit compatibility decision.

## Numerical and disclosure contracts

- Ring63 and Ring127 are fast paths selected only after public bound proofs.
- An exact multiprecision ring is the automatic fallback when either fast ring
  cannot certify the full operation.
- Exact GC/OT performs every chained signed truncation, comparison and final
  clamp. Modular wrap is never presented as a numerical result.
- Every result carries the capsule receipt, estimand, backend, scale, public
  magnitude bound, approximation domain and error/uncertainty certificate.
- Non-identifiable models return a typed mathematical state. Ridge, Firth,
  pseudoinverse or another estimand-changing remedy is used only when selected
  explicitly by the method contract.
- Small cells are handled by the same DP representation. No exact threshold or
  pass/fail guard bit is exposed.

Differential privacy bounds the change in an adversary's inference caused by
one privacy unit; it does not make a literal, distribution-free promise that an
adversary can never guess any original value. The package MUST state its exact
adjacency, parameters and pinned-peer/non-collusion assumptions and MUST NOT
market DP as deterministic impossibility of reconstruction.

## Performance and DSI contract

The expensive scan and vertical MPC occur once per capsule. Ordinary methods
then operate locally or on cached immutable slices, so repeated analysis avoids
both source-data scans and interactive MPC rounds.

Creation MUST use bounded-memory streaming and blockwise aggregation:

- O(frame) rather than O(dataset) transport buffers;
- direct scalar control messages and typed binary payloads with one unavoidable
  DSI encoding layer;
- asynchronous named fan-out, per-node forwarding state and bounded worker
  concurrency;
- exact short-write accounting, absolute offsets, acknowledgements,
  idempotent replay, bounded spool and backpressure; and
- fixed-shape or public-bucket semantic message geometry wherever non-output
  payload length could reveal protected state. DSI polling, retransmission and
  completion timing are not claimed constant-time by the R/DSI package path
  and require a separately scheduled gateway when they belong to the
  deployment threat model.

High-dimensional workload construction must exploit sparse domains,
block-vectorised marginal accumulation, packed bitsets where appropriate, and
parallel independent blocks. Pairwise work is intrinsically O(n p^2) if every
pair is requested; the capsule manifest therefore records the exact workload
and complexity certificate rather than hiding a combinatorial cost. Large
capsules are chunked, content-addressed and cached, not rebuilt for each method.

## Route-promotion acceptance tests

The package now exposes one disclosure-safe public mode; families that do not
meet its release contract remain provisional or fail before DSI. Promoting a
specific route within that mode requires the applicable checks below. The full
list remains the acceptance blueprint for complete package-family coverage:

1. Ten thousand repeated and syntactically varied operations over one capsule
   leave the ledger head, cumulative epsilon/delta and capsule bytes unchanged.
2. Different supported method families over one capsule never call the
   allocator, sampler or protected-frame resolver.
3. Configured `N` boundaries are exact: `N` committed capsules succeed, capsule
   `N+1` receives the fixed typed exhaustion condition, and an aborted commit
   remains charged. Boundary arithmetic covers `N * epsilon <= 8` and
   `N * delta < 1` without binary-floating ambiguity.
4. Concurrent different capsule proposals, reordered delivery, lost ACKs,
   crashes and restarts cannot fork or poison the durable head.
   Post-publication loss never produces a second release instance.
5. Neighboring-dataset mechanism tests, sampler tests and exact arithmetic
   tests cover domain extremes, NaN/Inf rejection and Ring63/127/multiprecision
   fallback without silent wrap.
6. Transcript tests prove that the relay sees no raw/pre-clamp value, exact
   small-cell bit, individual noise share, patient-derived public hash, or
   variable-size secret-dependent frame.
7. Statistical simulations cover bias, RMSE, coverage, Type-I error, power,
   imbalance, rare outcomes, missingness, repeated records, clipping, rank
   deficiency and non-identifiability for every promoted estimator.
8. Reconstruction, membership and attribute-inference attacks are included as
   adversarial empirical audits and interpreted against the formal DP bound;
   passing them is supporting evidence, not a replacement for the proof.
9. DSI benchmarks cover large row counts, wide schemas, reconnect/replay,
   backpressure, bounded memory/disk and byte-identical payload hashes.
10. Both R packages pass full tests and `R CMD check`; the MPC worker passes
    unit, fuzz, race, cross-platform build and two-peer multi-process E2E tests.

Obsolete or disclosive client-to-server primitives must remain absent from the
registered surface. A high-level method without a complete capsule adapter and
the applicable evidence above remains provisional or quarantined; it does not
create a second, less safe operating mode.

## Scientific basis for the hybrid design

- [NIST SP 800-226](https://doi.org/10.6028/NIST.SP.800-226) treats the query
  model, trust model, implementation correctness, side channels, utility and
  bias as necessary parts of evaluating a DP guarantee. It also records the
  uncertainty and subgroup-bias risks of synthetic data.
- [AIM](https://doi.org/10.14778/3551793.3551817) demonstrates the
  select-measure-generate pattern, workload-adaptive marginal selection and
  post-release per-query error bounds; it also explains why no synthetic
  representation is uniformly accurate for unrestricted workloads.
- [Noise-Aware Statistical Inference with Differentially Private Synthetic
  Data](https://proceedings.mlr.press/v206/raisa23a.html) shows that treating DP
  synthetic records as ordinary observations yields invalid uncertainty and
  develops noise-aware multiple-synthesis inference for supported discrete
  workloads.
- [Analysis of Differentially Private Synthetic Data: A Measurement Error
  Approach](https://doi.org/10.1609/aaai.v38i19.30114) gives another
  mechanism-aware route to bias correction and regression uncertainty rather
  than naive synthetic-data confidence intervals.
- [VertiMRF](https://doi.org/10.1145/3637528.3671771) and the later
  [MPC select-measure-generate framework](https://tpdp.journalprivacyconfidentiality.org/2025/pdf/fu.pdf)
  provide vertical-federated constructions and document the domain-size,
  communication and secure-marginal costs that the capsule implementation must
  benchmark rather than assume away.
- The biomedical evaluation [Does Differentially Private Synthetic Data Lead
  to Synthetic Discoveries?](https://doi.org/10.1055/a-2385-1355)
  reports inflated false-positive rates for several naive downstream
  combinations. This is why promotion is method-specific and requires Type-I
  error, power and coverage validation.
