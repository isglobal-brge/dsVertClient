# dsVert security and disclosure contract

Status: development hardening, updated 2026-08-11.

## Threat model

The intended deployment has separately administered vertical data
holders. Peers are authenticated by a server-name-to-Ed25519-key pinset.
The analyst is an untrusted relay: it may reorder, replay, omit or
modify messages, but it must not be able to substitute its own peer key
or decrypt peer-to-peer payloads. The current MPC protocols are designed
for honest-but-curious, non-colluding computing peers; at least one
designated share holder must remain honest.

The model does not cover a compromised server host, malicious protocol
deviations by a peer, collusion of all required share holders,
compromised long-term identity seeds, traffic-analysis guarantees, or
denial of service.

## Properties implemented and tested

- Stable server identities use 32 bytes from the operating-system
  CSPRNG, private atomic storage, symlink rejection and concurrent
  initialisation locking.
- Trusted peers are pinned by logical server name and exact Ed25519
  identity; duplicate identities, missing peers and key/name
  substitutions fail closed.
- If a peer identity is regenerated, server and client raise the typed
  `dsvert_peer_not_recognized` error with only the logical name and
  public expected/observed SHA-256 fingerprints. A server administrator
  must obtain the new public key through a trusted administrative
  connection, verify its fingerprint independently at that site’s
  console, update the name-bound pin on every participating server,
  restart the affected services and retry. There is no rotation counter,
  but a replacement is never accepted from the analyst/relay or learned
  automatically from the failing handshake.
- PSI messages carrying matched/common membership are session-, sender-
  and recipient-bound, signed and peer-encrypted. Plaintext patient
  indices are not remotely exposed.
- A completed PSI run attaches a secret-token-backed manifest bound
  locally to the current ordered identifier vector. Equal row counts are
  insufficient; reorder, subset, ID change or missing manifest fails
  validation.
- Session identifiers use CSPRNG UUIDv4 values. Weak
  [`sample()`](https://rdrr.io/r/base/sample.html),
  [`runif()`](https://rdrr.io/r/stats/Uniform.html) and time-derived
  session IDs have been removed from client protocol paths.
- Session files are contained below validated private directories,
  written atomically with restrictive permissions, and cleaned with path
  checks.
- The registered DataSHIELD aggregate surface is checked against an
  exact allowlist. Retired patient-index, generic reduction and unsafe
  legacy primitives are not registered or namespace-exported.
- No request-count quota is imposed. Input dimensions, byte sizes,
  finite numeric ranges and storage paths are still bounded as
  denial-of-service and arithmetic-safety controls.
- The still-provisional capsule-backed DP releases use an immutable
  server-owned policy, fixed patient or row add-remove adjacency,
  bounded one-unit contribution, two independently rooted designated
  peers, sticky query memoisation and persistent SQLite composition
  ledgers. Their signed server plan
  selects either the exact discrete-Laplace construction or the
  fixed-work dyadic discrete-Gaussian mechanism; HMAC-SHA256 domain
  separation and ChaCha20 streams make retries byte-identical. No root
  or derived seed is returned to the client. Concurrent identical
  requests commit one noisy answer.
- The capsule contract binds a cohort identifier, every approved
  snapshot, the current order-authenticated PSI manifest and the
  complete logical peer-name-to-Ed25519 pinset. Exact snapshot and
  ordered-ID commitments never appear in disclosure-safe status
  responses. Each server validates them locally, then exposes only the
  same ex-ante cohort/dataset attestation; the client requires that
  stable attestation and the exact pinset to agree across peers.
- Capsule ledgers use a canonical path in an owner-only `0700`
  directory. The database, lock, WAL and SHM are created under
  `umask 077` and kept at `0600`; symbolic links and non-regular ledger
  files fail closed.
- Production DP uses one capsule-publication sequence cross-signed by
  the pinned peers, never one independent epsilon ledger per peer or per
  analysis operation. The canonical immutable snapshot/workload capsule,
  creation entry, opening and exact payload commitments are durable
  before a delivery token can exist; retries and lost acknowledgements
  replay that mapping without resampling. Methods and arguments are
  post-processing and never advance this sequence. Each designated peer
  signs an independently generated public release-domain identifier in
  addition to its secret noise-root identity. PREPARE persists an
  authenticated candidate but does not claim the capsule, so sibling
  candidates may coexist. At each designated peer, the first valid START
  atomically and irrevocably claims one instance before local
  staged-source access or sampling. Source transport may already have
  staged encrypted protected material, but no noised share or public
  output exists at that boundary. Domain rotation may select another
  candidate only before the claim. Afterwards only the claimed instance
  may progress, restore or replay. A split relay cannot turn sibling
  local claims into the matching bilateral receipts required for
  publication while at least one designated peer remains non-colluding
  and retains its history. Once authenticated history proves publication,
  loss of the corresponding vector must restore and replay that exact
  release instance or fail closed; it can never rotate, resample or
  publish a second release for the capsule. External linearizable CAS
  optionally hardens the publication head and cannot select the retired
  local allocation path; its current receipt does not cover the separate
  result spool.
- Joint-DP result receipts use a v2 domain-separated HMAC-SHA256 payload
  commitment under a server-private persistent ledger key. The relay
  cannot enumerate a low-entropy payload from that commitment; the
  former public-hash v1 receipts are rejected. The final public contract
  exposes neither payload nor commitment. This assumes the server key
  and host remain uncompromised and does not by itself authorise or make
  an exact statistic safe to release.
- The joint-DP DSI bridge validates the complete K-peer pinset and
  contacts only the exact two canonical `designated_noise_peers` chosen
  by the custodians. Its prepare input is a
  server-local-HMAC-authenticated proposal; later inputs are canonical
  Ed25519 receipts verified again client-side. It exposes no generic
  proposal minter, statistic/noise/seed input, result stager or payload
  getter. Count does not use this capsule bridge: its purpose-bound
  endpoints require a full-K signed analysis contract, authorize exactly
  two identity-bound peers and return one bounded release signed by the
  finalizer.
- Every client DP status, discovery and release call forces DSI
  `errors.print = FALSE` and rejects partial/NULL site results
  uniformly. Purpose-bound Count phases are session-idempotent and
  recompute the same artifact-bound noise from persistent identity
  seeds; capsule-source and joint control phases retain their separate
  durable ledger idempotency. The relay may replay only a byte-identical
  request after an ambiguous missing ACK. Retry has an inactivity
  deadline but no request-attempt quota and cannot reroll a Count
  artifact, capsule allocation, noise draw or payload. Detailed node
  errors remain in privileged server logs rather than becoming an
  analyst-visible disclosure channel.
- Each immutable capsule has fixed epsilon/delta parameters. There is no
  request-count quota, geometric decay or accuracy suppression caused by
  prior requests. Reusing the same capsule/release instance is zero-cost
  post-processing and is never denied by request count. A new capsule
  burns one non-refundable distinct-capsule reservation unit at
  allocator commit and is denied after the server-owned maximum `N` is
  exhausted. Policy is accepted only when exact arithmetic proves
  `N * epsilon <= 8` and `N * delta < 1`. Thus
  `history_can_deny_operation=TRUE`, `operation_limit=TRUE`,
  `request_limit=FALSE`, and the registry counts reservations separately
  from publication. The claim is scoped to
  `at_most_N_immutable_snapshot_workload_capsules_per_stable_privacy_accountant_namespace`
  and assumes
  `privacy_accountant_namespace_assumption=one_stable_unique_namespace_across_domain_cohort_policy_pinset_and_ledger_reconfiguration_per_protected_privacy_universe`.
  Namespace continuity is a custodial deployment obligation: neither
  package enforces global uniqueness or automatically migrates
  reservations across reconfiguration. The fixed, non-retryable
  `[dsvert_dp_lifetime_budget_exhausted:v1]` token is an opaque terminal
  union of global `N` exhaustion and a requested capsule whose
  irrevocable instance claim/publication binding prevents that instance
  from safely continuing or replaying; it does not imply
  `remaining_distinct_capsules == 0`.
- IKNP KOS consistency validation is mandatory; analyst and custodian
  options can no longer select unchecked OT. Unused participating-party
  dealer triple generators have been removed from registration and
  namespace export.

These properties have unit, adversarial and two-host PSI tests. They are
not a formal cryptographic proof or an external security audit.

## Output privacy boundary

MPC protects intermediate values while they are being computed. It
cannot make an exact released statistic reveal less than that statistic
mathematically contains. With unlimited adaptive queries, exact
histograms, recoded tables, scores at analyst-selected parameters,
cluster sums or related outputs can be combined to infer sensitive
values without breaking encryption.

Therefore dsVert cannot honestly guarantee all three of the following at
once:

1.  exact answers;
2.  arbitrary, unlimited adaptive analyses;
3.  impossibility of reconstructing the source data.

Pinned peers address the most dangerous relay attack—key substitution
followed by share opening—but do not change this output-composition
limitation. Strong non-reconstruction needs either a narrow
server-enforced catalogue of final analyses on immutable cohorts, or
differential privacy/composition controls.

The promoted Count route is separate from the capsule vector. Every peer
signs one canonical current-snapshot contract; exactly two
identity-bound authorities run the fixed scalar exact-GC mechanism and
only the finalizer opens one bounded signed value. The same canonical
artifact recomputes the same noise after a new session or process, while
distinct artifacts compose and have no finite global composition claim.
Under replace-one fixed-cohort adjacency, all K peers sign the same exact
public cohort size and no DP noise or MPC session is created.

The still-provisional capsule-backed table, moment, Gaussian-model and
survival routes select coordinates from their separate signed sticky
biomedical vector. Promotion of a client postprocessor does not by
itself promote the upstream release that it consumes. The vector’s
public deterministic selector uses the scalable
two-independent-complete-draw Laplace convolution or a formal fixed-work
dyadic discrete-Gaussian route when its exact server plan strictly
improves the certified simultaneous radius. The analyst cannot choose
epsilon, adjacency, categorical domains, clipping bounds, workload or
ledger state. Every supported method over the same immutable capsule
reuses its authenticated release, preventing averaging through retry or
a syntactic alias.
[`ds.vertDPCalibrate()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPCalibrate.md)
reports data-free floating-point utility previews, but never presents
them as the executed choice: only the signed server manifest contains
the exact selected mechanism and certificate. Operation count never
changes those values and calibration does not authorize a release.

These computational DP guarantees are scoped to each method’s published
adjacency and release contract after a custodian-approved snapshot has
passed setup validation. Count is conditional on secrecy and persistence
of the two authority identity seeds and correct pinned-peer operation;
capsule routes are additionally conditional on their documented
noise-root and authenticated ledger assumptions. dsVertClient itself
remains Windows-compatible.
Compromise of a noise root can reveal the keyed randomness for known
release contexts and is outside the guarantee.

For every successful capsule publication, authenticated
semantic-operation counts, chunk counts and non-output payload geometry
are functions of the signed public manifest: K, the fixed two designated
peers, the public source-owner set, coordinate count and chunk capacity.
Here a semantic operation is one authenticated protocol call at a named
site, not one HTTP request, DSI poll, TCP write or retransmitted frame.
Every protected pre-release failure has one relay-visible error token;
detailed admission, snapshot, range and integrity messages are not
returned. Final-vector text and its transport size may depend on the
already-DP released vector and are therefore post-processing. A memoized
replay may be faster only under the public capsule and release-instance
identity; it never resamples or consults source data.

This is not a host constant-time or traffic-flow claim. The occurrence
and timing of pre-release R admission/setup failures, process
scheduling, cache state, DSI/HTTP framing, polling, retransmission,
availability and denial of service are excluded. In particular, valid DP
use is conditional on the custodian snapshot already belonging to the
published bounded-admission domain; a rejected setup is not presented as
a DP release. Deployments requiring network-observer traffic-flow
privacy still need a separately validated scheduled/asynchronous
gateway. dsVert deliberately does not add long sleeps or pad every
transfer to the package-wide maximum: neither would make the R/DSI host
path constant-time, and both would severely reduce availability.

The guarantee also does not cover historical exact implementations
(which the installed disclosure-safe gate keeps unavailable),
compromised hosts, collusion of every designated peer, or simultaneous
rollback/corruption of all peer result ledgers. The optional v2 external
anchor currently protects only the allocation head. After irrecoverable
post-publication release-state loss, the capsule must restore and replay
its exact release instance or fail closed; automatic regeneration,
rotation, resampling or a second public instance is forbidden.
Destroying or rolling back every ledger can also destroy that local
accounting history and is a compromised-deployment case, not a
package-level rollback guarantee. The owner-only directory prevents
unprivileged filesystem substitution; an administrator able to act as
the service account or replace every retained ledger is a
compromised-host/deployment case. The claim requires at least one
honest, non-colluding designated peer to retain and use its complete
authenticated monotonic history; a relay cannot make it accept a forged
peer identity, receipt or payload. DP bounds inference risk; it is not
an absolute logical impossibility of reconstruction and there is no
universal epsilon that optimises every study.

## Numerical boundary

Ring63 and Ring127 encode fixed-point numbers. Their historical
local-share multiplication truncation is probabilistic even with proven
headroom and fresh uniform shares; the historical DCF/wide-spline
comparison also lacks a full-domain exact-comparison contract. Those
implementations remain only as unregistered, quarantined internal
migration code and are not used by a promoted analysis path.

Every promoted path that requires secret truncation or comparison uses
its purpose-bound exact-GC adapter, with canonical share encoding,
custodian-owned public bounds and fail-closed shape/range validation.
Unsupported bounds or an unavailable certified adapter produce a typed
failure instead of selecting a local-truncation or DCF fallback. This is
an exact integer/ring operation claim, not a proof of exact real
arithmetic, a whole-estimator error bound, or a constant-time claim for
the R/DSI host. Ring127 by itself only lowers ordinary quantisation
error; it does not make a legacy truncation chain exact.

## Method maturity

Run:

``` r

ds.vertMethodStatus()
ds.vertMethodStatus(status = "quarantine")
```

Only methods marked `promoted` should be described as release-ready, and
even then only inside the `safe_scope` and assumptions shown by the
registry. `provisional` means additional cryptographic, numerical or
statistical evidence is required. `compatibility` preserves a historical
estimand or name. `quarantine` means the current route should not
support confirmatory biomedical claims.

## Communication performance

DataSHIELD/DSI remains a request-response channel. Efficient portable
streaming uses one named fan-out across sites per dependency-ordered
pump cycle and one parser-bounded frame per target in that call. Before
any legacy/typed blob payload, a stateless public-ASCII probe selects
the largest common 8/4/2/1/0.625 MiB text geometry accepted by every
target, with a fixed expression-overhead reserve. The measured DSLite
profile starts at 640 KiB text/approximately 480 KiB raw; failure
descends through smaller public probes and never applies an unverified
opaque-payload fallback. Absolute offsets and exact replay make
ambiguous retries idempotent. Once a payload starts, geometry cannot be
re-probed or reduced. Persistent sessions and peer fan-out remove
avoidable setup and serial-site latency. Immutable frames use DSI’s
public job primitives directly, avoiding the high-level helper’s
repeated large-expression formatting during polling while preserving
connector serialization, named result mapping, keepalive and error
callbacks. The current DSI API does not expose a portable, ordered
multi-request window on one connection, so dsVert adds no generic
application/request batching layer. The capsule-source route has one
purpose-bound exception: a single fetch or accept request may carry a
byte-bounded consecutive window of unchanged inner artifacts after every
participant attests that capability. Its width is the minimum of the
source node’s verified response capacity, both designated recipients’
per-site request/expression capacities, and the signed 8 MiB application
caps. Response capacity is never inferred from request capacity. The
response probe reserves 12.5% plus 128 KiB for connector serialization,
and window construction reserves another 64 KiB for expression/framing
overhead; the hard chunk count remains eight. The outer window framing
is not signed; every inner ticket, ciphertext envelope and
acknowledgement remains independently signed and bound. This does not
create multiple outstanding DSI jobs. Independent peers use one named
fan-out; parser-bounded frame chunking remains dependency ordered and
the source pump is pipelined to `F + 2` healthy phases.

The signed adaptive negotiation is version 2. Within one synchronized
DSV1 client/server package generation, the pre-existing `byte-window-v1`
capability and omitted framing arguments retain their exact 768-KiB
source-response and 1-MiB recipient-request caps. Only a complete v2
attestation across both recipients and every source makes the client
attach the public v2 framing contract; a same-generation peer that lacks
v2 uses the byte-identical scalar route. DSV1 itself is a lockstep wire
change, not an old-client/new-server compatibility layer. The safe
rollout is server-first with release traffic paused: deploy, reconcile
and attest every server, verify that the framed data-free manifest
boundary is active, and then deploy the paired client. Any intermediate
mixed package generation fails before protected preparation.

DSI also exposes no portable cancellation for an asynchronous
`DSResult`. If a job remains unresolved at the monotonic deadline,
dsVert does not pretend that it stopped remotely: the authenticated
connection is poisoned for this client process and a fresh DSI
login/session is required before replay. This prevents Armadillo 3.0.x
`lastcommand`/`lastresult` state from being fetched as the result of a
later submission. dsVert never calls `dsDisconnect()` silently.

An expected cold miss in the allocation, result or release replay
lifecycle is different from an unresolved transport. The server emits
only the fixed `[dsvert_phase_not_ready:v1]` application tag, and the
client advances a prerequisite phase only after every relevant completed
result was fetched and carried that tag. A fetched terminal application
result is settled and leaves the DSI session reusable. Mixed, untyped,
partial, reset, timeout and ambiguous 500 paths remain terminal and
poison any handle whose result lifecycle cannot be proven settled; the
tag is never an automatic retry signal.

The probe endpoint reads no protected data, MPC state or DP ledger,
retains no state and caps request and response padding independently at
8 MiB. Its three-argument v1 acknowledgement is byte-identical; the
optional response extension returns only fixed public `R` padding, its
length and digest. A confirmed HTTP 413 may select a smaller public
candidate. An ambiguous singleton-result lifecycle poisons that exact
handle, and a present malformed acknowledgement is fatal. A server
without the extension uses scalar capsule source transfer. Probe results
are availability optimisations, not TLS or peer-identity attestations.
Exact session results are cached by authenticated handle. A
connector-profile cache only reorders later probes: every reconnect
revalidates its hint before opaque payload, and an unknown connector
version executes the complete ladder.

dsVert contains a signed, identity-bound, idempotent fan-out relay
substrate with replay/gap/conflict rejection and no request-count limit.
The first legacy blob replacement substrate is a registered
producer-ticket endpoint that has no caller-selected key, file or
purpose; its 20 producer/consumer classes are inventoried but are not
claimed as migrated until each method state machine is bound and
integration-tested. The more general relay remains internal until it is
validated across independent hosts with supported live DSI connectors.

For producer-spooled streaming, every frame is hash-checked at the
relay, the recipient verifies the complete ticket-bound SHA-256 before
atomic commit, and the producer verifies the pinned recipient’s
signature over hash and length. The analyst therefore stores no
redundant O(payload) integrity file; that file could not be a security
trust anchor because the analyst controls it.

Client-side TLS inspection is a fail-safe against benign deployment
mistakes, not a defence against a malicious analyst that controls the
client process. A recognized Armadillo connection with an inspectable
HTTPS endpoint is accepted automatically, matching dsFlowerClient; no
per-session attestation option is required. The Armadillo frontend or
reverse proxy remains authoritative for certificate and hostname
validation. Remote HTTP, a detected process-wide TLS downgrade, a
present but uninspectable global httr configuration, a missing endpoint,
and unknown generic connectors are rejected. An analyst-controlled
option cannot promote an unknown connector. The security boundary
against an untrusted relay is the server-side release contract plus
name-bound peer pinning and authenticated peer messages; outer TLS still
protects DataSHIELD credentials and ordinary responses in transit.

Surface readiness is separate from connector TLS recognition. The same
connector-neutral custodian token is used for Opal and Armadillo/Rock,
but it must be installed server-side only after connector-specific
administrative tooling verifies the exact callable inventory. The client
submits only `dsvertSecurityProfileDS()` with no arguments and cannot
forward an option, environment value or attestation. Missing, malformed,
conflicting or stale server assertions fail closed; dsVert never infers
readiness merely from an installed package. The Rock procedure is
documented in `inst/docs/remote_surface_attestation_20260808.md`.

Topology and connector claims are also separate. `K=2`, `K=3` and `K=5`
coverage for the promoted Count and PSI routes, and for provisional
capsule routes, means the named unit, adversarial, isolated-process or
real-DSLite harness unless a current live connector artifact is cited.
The local Armadillo S4/httpuv lifecycle regression is not a live Rock
service, and cached Opal runs are historical evidence only. The
analytical protocol constructs DSI calls; Opal reconciliation and
Armadillo TLS/session inspection are deployment safeguards. Fresh
artifact-exact Opal and Armadillo/Rock smoke remains part of the
deployment release gate.

## Reporting vulnerabilities

Do not include patient data, keys, tokens, raw shares or live service
credentials in an issue. Report the affected version, method, topology
and a synthetic reproducer through the repository’s private
security-reporting channel where available.
