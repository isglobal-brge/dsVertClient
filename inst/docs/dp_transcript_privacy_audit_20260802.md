# Joint-DP capsule transcript privacy audit

Status: implemented boundary and adversarial unit coverage, 2026-08-02.

## Claim boundary

The package makes a bounded lifetime `(epsilon, delta)` claim for at most `N`
released immutable capsule vectors per stable privacy accountant namespace.
For a successful publication, every authenticated semantic
protocol message observable by the analyst-relay is either:

1. determined by the signed public contract; or
2. deterministic post-processing of the released DP vector.

This is not a constant-time R, operating-system, HTTP or DSI claim. The
occurrence and timing of pre-release setup/admission failure, scheduling,
availability, DSI polling, physical retransmission and denial of service remain
outside the formal package claim. A snapshot rejected before publication is
not reported as a DP release. Covering those traffic observables requires a
separately validated scheduled/asynchronous deployment gateway.

The pre-release ABI fixes control at `dsvert-joint-dp-control-v3`, capsule
identity at `dsvert-joint-dp-capsule-identity-v3`, PREPARE at v6,
START/RESULT/RELEASE/ACK at v5, STORE at v6, the claim row/set at
`dsvert-joint-dp-vector-instance-claim-v1` and
`dsvert-joint-dp-vector-instance-claim-state-v1`, and replay at v4. Control v3
and identity v3 bind the lifetime fields and biomedical workload v7. Every signed public
phase receipt and replay response attests exactly `history_gate=TRUE`,
`request_limit=FALSE`, and `operation_limit=TRUE`; the durable STORE schema is
separately HMAC-authenticated. Legacy v2 control/identity artifacts and v4/v5
stores are neither migrated nor re-signed automatically and fail closed.
Rollout therefore requires empty DP capsule state or a future audited offline
migration. Earlier Opal deployments
had `POLICY_READY=FALSE`, no published DP capsule and only ephemeral local
K-site state.

Exact COMMIT/RELEASE endpoint replay and sticky replay in the live session are
O(1). End-to-end cold reconstruction after process restart instead reaches
AllocationProof, which performs a complete O(N) allocator-journal audit before
returning the proof. That difference is computational only: the cold path still
burns no lifetime unit, invokes no noise sampler and reads no protected source.

Here, *semantic message* means one authenticated protocol operation at a named
site. It does not mean one HTTP request, one DSI poll, one TCP write or one
retransmitted frame. Those physical counts can change with connector and
network state even when the protected dataset is unchanged.

## Relay-visible channel map

| Observable | Possible dependency | Enforced treatment |
|---|---|---|
| Logical site names, K and selected compute peers | Public pinset | Signed public contract; exactly two designated peers for every K >= 2 |
| Source-owner set | Signed schema ownership | Public and canonical before protected access |
| Coordinate and chunk counts | Signed workload and fixed chunk capacity | Public deterministic geometry |
| Source ciphertext length | Public coordinates in the chunk plus fixed-format encrypted header | No protected value controls it; final short chunk is publicly derivable |
| Cross-owner rounds | Signed artifact dimensions and padded unit capacity | Fixed public schedule; no match/event count is returned |
| Sampler rounds | Public mechanism plan and coordinate count | Exact-GC circuit shape or fixed-work public loop bounds |
| Final-vector JSON length | Noisy clamped vector | DP post-processing; the vector itself is released |
| Success receipts and progress flags | Public protocol state | Signed, canonical and replay-idempotent |
| Cache hit / replay fast path | Public capsule and release-instance identity | May be faster; never re-enters protected data or resamples |
| Retry/reconnect/poll count | Network and connector availability | Outside the claim; absolute offsets prevent semantic duplication |
| DSI expressions, method names, offsets and site mapping | Public protocol schedule | Parser-bounded and validated before submission |
| DSI polling and completion time | Host work, scheduling and transport | Explicitly excluded; no constant-time claim |
| Protected admission, snapshot, range or integrity error text | Could depend on protected state | Replaced by one public token: `[dsvert_dp_public_failure:v1] Protected capsule operation failed.` |
| Failure occurrence and phase reached | Can depend on pre-release setup | Excluded; valid DP use is conditional on the published bounded-admission domain |
| Resource oversize/backpressure | Public message geometry and process capacity | Remains typed for safe retry/terminal handling; never a request-count or privacy-budget gate |
| Lifetime/instance terminal refusal | Global `N` exhaustion, or an irrevocable instance-claim/publication binding that prevents the requested instance from safely continuing or replaying | Fixed detail-free `[dsvert_dp_lifetime_budget_exhausted:v1]`; opaque union that does not imply `remaining_distinct_capsules == 0`, because cause-specific tokens would reveal state |
| Unknown peer / pre-claim root refresh | Public authenticated identity state | Remains typed so administrators and the bounded automatic refresh path can recover before the first valid START; it never authorizes a sibling after claim |

The full ordered-PSI alignment hash in the current R cross-source path remains
inside each recipient-encrypted header. This audit claims confidentiality from
the analyst-relay, not recipient-specific secret-sharing of that hash or
malicious-peer security.

## Successful base cold-path schedule

Let `K` be the public number of custodians, `S` the public number of source
owners, `C = ceil(D / B)` the public number of vector chunks, `D` the signed
coordinate count, and `B` the fixed chunk capacity. The base cold successful
path has the following remote-site invocation multiplicities:

| Phase | Remote-site invocations |
|---|---:|
| Status handshake | K |
| Manifest draft, sign and build | 3 K |
| Allocation proof lookup | 2 |
| Cold allocation prepare, commit, authorize and open | 8 |
| Recipient ticket | 2 |
| Source prepare | S |
| Source ciphertext fetch | S C |
| Recipient acceptance | 2 S C |
| Vector prepare | 2 |
| Vector-result lookup | 2 |
| Vector start | 2 C |
| Vector-result commit | 2 |
| Vector-release lookup | 2 |
| Final-share fetch | 2 C |
| Vector-release commit | 2 |
| Final DP replay | 2 C |
| Durable finalization acknowledgement | K |

For `D = 16,385`, `B = 8,192`, `C = 3` and `S = K`, the tested totals are
respectively 70, 85 and 115 remote-site invocations for K = 2, 3 and 5. The
corresponding minimum sequential semantic dependency-round counts are 39, 45
and 57. Different protected values do not enter these formulae. A
byte-identical public release replay may skip the cold
materialization/sampling path.

Those totals deliberately describe the base control plane. Publicly declared
cross-artifact stages, exact-GC inner exchanges and typed-blob frame pumps add
work derived from their signed artifact/circuit/payload dimensions. DSI polls,
backpressure cycles and physical retransmissions are not included, because
they are transport-state observables rather than semantic protocol messages.

## Implementation changes

- `dsVert/R/dpTranscriptPrivacy.R` is internal and adds no remote endpoint.
- The source, vector, allocation, Gaussian-cross and categorical-cross public
  boundaries use one protected-failure representation.
- Peer-not-recognized, resource backpressure/oversize, independent-root and
  pre-claim current-release-root refresh conditions remain typed because their
  contents are public operational state and callers need the recovery action.
  Once START has claimed the capsule, root change never authorizes a sibling
  instance.
- The signed policy no longer says or implies that all transcript properties
  are equal. It states the successful semantic-message boundary and explicitly
  excludes physical timing, polling, retransmission and availability.
- No artificial sleeps, request-count limits or package-wide maximum padding
  were added. The authenticated lifetime-history gate is deliberately separate
  from transport chunking and backpressure, which remain bounded by public byte
  geometry.

## Verification

The focused server tests cover:

- byte-identical error text across unrelated protected failures;
- absence of secret markers and phase labels from public errors;
- preservation of the five typed public recovery/resource classes;
- exact sanitization and terminal propagation of the lifetime-exhaustion token;
- public base-control-plane invocation and dependency-round geometry for K =
  2, 3 and 5, including compact maximum-size geometry; and
- explicit negative assertions for R constant-time and traffic-flow DP.

Recorded focused results:

- `dsVert`: 57 transcript-privacy assertions passed;
- `dsVertClient`: 89 core-vector assertions passed;
- `dsVertClient`: 202 K-aware capsule-status assertions passed.

The broader server suite could not be interpreted during this audit because
agents were concurrently rebuilding `dsvert-mpc`: packaged binaries and
`SHA256SUMS` temporarily disagreed. In that run 117 assertions passed before
the first checksum failure. The complete suite must be rerun after the final
single binary rebuild/checksum step; a checksum bypass is not acceptable.

## What is and is not defensible

Defensible:

> Conditional on an admitted immutable snapshot, the declared adjacency and
> bounds, authenticated protocol-compliant peers, and at least one
> non-colluding designated noise peer retaining and using complete authenticated
> monotonic history,
> plus one stable, unique privacy accountant namespace across domain, cohort,
> policy, pinset and ledger reconfiguration for the protected privacy universe,
> at most the configured `N` successful released vectors satisfy the exact
> bounded lifetime composition. Successful authenticated semantic-message geometry is
> public, non-output-dependent message sizes are public, and output-dependent
> serialization is DP post-processing. Retry cannot reroll the release.

Simultaneous rollback of both designated histories is not covered without an
external linearizable CAS.

The namespace condition is a custodial deployment assumption, not automatic
cross-namespace enforcement or counter migration by either package.

Not defensible without another deployment component:

> The complete network transcript, response time, failure occurrence and
> availability of an R/DSI request are constant-time or differentially
> private.

Padding every payload to a package-wide maximum or sleeping for a guessed
duration would not establish that claim and would damage throughput. A real
traffic-flow claim needs a public service calendar, asynchronous enqueue/serve
boundary, fixed externally audited envelopes and an independently validated
gateway implementation.
