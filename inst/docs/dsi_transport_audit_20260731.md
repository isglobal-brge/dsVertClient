# dsVert direct-DSI transport and scalability audit

Status: implementation audit and local DSLite benchmark, updated 2026-08-08.

## Scope and claim boundary

dsVert transports DataSHIELD calls directly through DSI.  It does not tunnel a
second RPC protocol such as Flower or gRPC through DSI.  Consequently there is
no inner message envelope to remove.  Opaque MPC values use one canonical,
unpadded Base64url encoding because an R/DSI expression cannot carry arbitrary
raw bytes portably.

The Flower tunnel's non-blocking `writeBin()`, short-write shim, TCP byte ACK
and local forwarder-port findings do not transfer literally: dsVert submits one
complete R call to each DSI connector and has no analyst-side socket stream.
Their equivalents here are pre-parse expression bounds, immutable absolute
application offsets, typed terminal receipts, bounded server spools and DSI
availability replay. Adding a second client spool or socket ACK layer would
duplicate work without making the connector's atomic HTTP/DSI submission more
correct.

This audit supports claims about request geometry, replay correctness, bounded
local resources and in-process DSLite behaviour.  It is not a production
multi-host latency or TLS benchmark.  Opal and Armadillo deployments still need
connector-specific load testing over their real HTTP/TLS path.

## Verification snapshot (2026-08-02)

- The exercised package versions were DSI 1.8.0, DSLite 1.4.1 and
  DSMolgenisArmadillo 3.0.2. A real `ArmadilloConnection` built with
  `httr::handle()` exposes its HTTPS URL and is accepted without
  `options(dsvert.dsi_tls_attested=...)` or any other dsVert-specific login
  input. The installed S4 methods for `dsAggregate`, `dsIsCompleted`, `dsFetch`
  and `dsKeepAlive`, and the connector's asynchronous aggregate capability,
  were exercised structurally. HTTPS downgrade, process-wide curl verification
  downgrade and malformed global httr configuration were rejected.
- The committed `test-dsi-armadillo-connector-smoke.R` localhost regression
  exercises those real
  `DSMolgenisArmadillo` S4 methods end to end through HTTP: both the direct
  fan-out cycle and the central strict aggregate submitted `/execute`, polled
  `/lastcommand`, fetched the serialized `/lastresult`, and recovered the
  expected typed value without any
  dsVert-specific connection option. This isolates connector integration and
  lifecycle logic; it is not evidence about a real Armadillo server, its TLS
  termination, authentication deployment or WAN performance.
- The client suites `dsi-security`, `dsi-fanout`, `dsi-transport-probe`,
  `dsi-sliding-window`, `dsi-benchmark-artifact`, `typed-blob-transport` and
  `dp-capsule-source-client-transport` all pass. The fan-out suite creates two
  actual DSLite servers/connections, maps distinct expressions by site, covers
  partial failure, and verifies launch-before-poll, fetch and keepalive order.
- The server suites `dsi-relay-substrate`, `transport-probe`,
  `typed-blob-transport`, `dp-capsule-source-transport`,
  `joint-dp-dsi-bridge` and `exact-gc-transport` all pass.
- A seven-repeat cold request-geometry benchmark on local DSLite 1.4.1 used
  the version-specific 640 KiB public probe and recorded medians of 466 ms for
  K=2 (range 461--846 ms) and 1,148 ms for K=5 (range 1,145--1,798 ms).
  Exactly one server probe ran per peer. Reapplying the authoritative geometry
  for the same authenticated session set over 500 iterations made zero server
  calls and averaged 0.242 ms for K=2 and 0.294 ms for K=5. Thus synchronous
  cold negotiation has O(K) connector cost, while the exact-session cache makes
  repeat use effectively client-local. These are in-process timing observations,
  not Opal, Armadillo, HTTP/TLS or WAN latency claims.
- A warmed, order-balanced two-peer, 1 MiB-per-peer DSLite A/B produced
  identical frame geometry, submitted bytes, replay/reconnect counts and
  terminal SHA-256 in every pair. Across four alternating pairs, median
  standard dispatch was 12.264 s (0.164 MiB/s) and median direct dispatch was
  8.825 s (0.234 MiB/s): 28.0% lower median wall time and 43.3% higher median
  useful throughput in this local run. All four recorded pairs favoured the
  direct path, but individual times remained variable; this is local DSLite
  evidence, not an Opal/Armadillo or WAN speed guarantee.
- A fresh 4 MiB typed source-stream smoke completed eight 512 KiB frames at
  7.87 MiB/s with one generation, two producer-frame replays, one recipient
  replay, one reconnect, both backpressure preflight rejections, zero
  controller payload-spool bytes and terminal `sha256_identical`. This
  validates the harness and bounded protocol, not remote-network throughput.
- An exploratory direct-dispatch scaling sweep on 2026-08-02 used K=2/3/5 and
  0.0625/0.5/2 MiB per peer. At 2 MiB it completed in 5.552/8.269/14.035 s,
  corresponding to 0.720/0.726/0.713 aggregate useful MiB/s. A paired 0.5 MiB
  standard-versus-direct run reduced wall time by 19.4%, 34.2% and 29.5% for
  K=2/3/5 respectively; submitted bytes, cycles, replay, reconnect, spool size
  and terminal SHA-256 were identical between dispatchers. Those earlier
  pairs did not alternate engine order and are retained as exploratory scaling
  observations rather than the primary dispatcher-effect estimate.
- Terminal retention stress now covers 50 sequential relay transfers, 100
  expired sessions and 128 durable SQLite store providers. Relay control state
  was exactly 12,864 bytes after transfer 1, 5 and 50; only one replay window
  remained, sender retained bytes were zero and recipient retained bytes were
  the fixed 4,160-byte payload-plus-metadata reservation. A 26-frame typed
  receipt occupied 10,336 bytes before consumption and 8,112 bytes after
  compaction, with one terminal replay frame and zero retained payload bytes.

No remote Armadillo service endpoint was available. Therefore this snapshot
does not claim a real Armadillo server, proxy, certificate-chain, or WAN
latency E2E. The installed connector contract is tested; a production remote
E2E remains an explicit deployment validation item.

## Implemented transport contract

- Independent site calls in one protocol phase are sent as one named,
  asynchronous DSI fan-out.  Results must contain every expected logical site,
  exactly once, and are reordered by site name before use.
- Independent assignment expressions use DSI's expression-assignment API as
  one named fan-out. Completion requires exactly one success callback for each
  logical site; missing, duplicate, unknown or mixed success/error callbacks
  reject the phase rather than accepting a partial assignment.
- The two reciprocal exact-GC validity receipts for each checked multiplication
  chunk are mapped by target and submitted as one named fan-out. A missing
  target can replay only that same idempotent phase; no chunk commits until both
  peers report the checked state.
- Each validated exact-GC envelope is serialized to its immutable delivery JSON
  once. Missing/partial DSI responses reuse those exact bytes until the target
  acknowledges the absolute offset, rather than repeatedly copying and
  serializing the payload during polling.
- Immutable frame/exchange fan-out uses DSI's public `dsAggregate`,
  `dsIsCompleted`, `dsFetch` and `dsKeepAlive` primitives directly. DSI's
  standard high-level helper deparses the complete expression again during
  polling and progress formatting; the direct path preserves lazy sessions,
  asynchronous launch-before-fetch, connector polling intervals, keepalive,
  callbacks and named `NULL` failures while leaving serialization to the
  connector once. Every production call that crosses the central strict
  aggregate boundary uses this dispatcher; injected aggregate functions remain
  available only to isolated tests. A remote statistical job now runs to
  completion by default instead of inheriting the five-minute idempotent-frame
  retry deadline. An operator may set the optional monotonic
  `dsvert.dsi.job_timeout_seconds` ceiling; reaching it returns a named
  unavailable result without claiming to cancel a still-pending job, because
  DSI exposes no portable cancellation primitive. Any authenticated
  session with an unresolved job is marked poisoned for the client process;
  replay requires fresh DSI login connections, preventing Armadillo's
  singleton `lastcommand`/`lastresult` lifecycle from being misassociated.
  With the default infinite job deadline, explicit user interruption remains
  the escape hatch for a connector that never settles; that interruption also
  poisons the ambiguous handle rather than silently reusing it.
- This mapping contract requires DSI >= 1.8.0 in both package metadata and a
  fail-closed runtime guard; an older DSI cannot silently fall back to positional
  or sequential routing.
- A transport exception, missing result, unexpected `NULL`, malformed
  acknowledgement or wrong result association cannot advance the phase. The
  server's error text is not reflected to the analyst.
- Protocol phases are not retried implicitly. The replayable paths are an
  audited allowlist of producer calls that cache their complete response under
  the identical request hash, plus absolute-offset blob/exchange operations.
  An absent/ambiguous response replays the identical expression or frame until
  a configurable monotonic inactivity/availability deadline; there is no
  attempt-count or request-count quota. The deadline is renewed only by a
  verified ACK, new immutable bytes/envelope or a terminal state transition;
  retry traffic and empty polling do not renew it. Producer/legacy callback
  errors, present rejections, malformed ACKs,
  signature failures and contract conflicts are fatal immediately. A
  purpose-bound typed store/receipt callback is terminal unless it carries the
  fixed public capacity-backpressure tag; semantic rejections also return a
  present typed result. Only an absent/ambiguous outer result or explicitly
  tagged backpressure can replay an identical typed request. A generic
  typed-transport deadline returns the
  structured `retry_deadline_exceeded`/`transport_unavailable` condition and
  preserves its resumable server state. If the deadline belongs to an
  unresolved asynchronous DSI job, replay first requires a fresh authenticated
  connection; the old handle is never silently reused. Exact-GC instead aborts and cleans the
  current private attempt on inactivity; any later retry must enter its
  purpose-bound retry contract and no partial result is returned.
- Producer retry eligibility is checked from the complete call form, not only
  its function name. IKNP input/output-to-legacy-blob arguments must be absent,
  `NULL` or empty. Transitional mutating branches therefore never inherit
  automatic retry from an active-first producer with the same name. Retired
  DCF comparison calls are neither registered nor present in this allowlist.
- Final `on.exit` cleanup is the sole best-effort exception.  It never advances
  protocol state, never exposes connector error details and will not cross an
  insecure outer transport.
- Real DSI calls accept DSLite in-process connections or verified HTTPS.
  Remote HTTP is rejected.  Opal must expose enabled certificate verification;
  a recognized Armadillo connection with an inspectable HTTPS endpoint is
  accepted automatically, matching dsFlowerClient, without a per-session
  option. Its frontend or reverse proxy remains authoritative for certificate
  and hostname validation. An unknown generic connector is rejected because
  an analyst-controlled option is not evidence of a connector TLS contract. A
  plaintext loopback exception exists only under explicit test mode.
- The portable default payload frame is 480 KiB raw.  It becomes at most
  640 KiB Base64url and remains below the measured DSI 1.8/DSLite 1.4.1
  parser boundary. Generic expressions are completely deparsed and measured
  before DSI sees them. Typed Base64url frames and normal legacy Base64 frames
  use conservative fixed-schema byte bounds instead, avoiding an O(frame)
  validation deparse while retaining a 1 KiB syntax reserve; regression tests
  prove both bounds exceed the actual deparsed size. An unusual legacy payload
  outside the normal Base64/storage alphabet falls back to exact generic
  deparsing rather than acquiring a compatibility rejection.
- Before the first legacy/typed blob payload, a stateless data-free probe tries
  8, 4, 2, 1, 0.625, 0.3125, 0.15625, 0.078125, 0.03125 and 0.015625 MiB of
  public ASCII text in descending order. The exact measured DSLite 1.4.1
  class/version starts at its recorded 0.625 MiB ceiling; unknown/new connector
  versions retain the full ladder. A size is selected only when every named
  target returns the exact typed acknowledgement. Missing or transport-failed
  probes may descend; any present malformed response is fatal. Probe calls are
  deliberately synchronous: they contain no data or protocol state, and this
  prevents an oversized candidate rejected after submission from leaving an
  unresolved async job that poisons the login before the next smaller public
  candidate. This one-time/cached negotiation trade-off does not change the
  named asynchronous fan-out used by statistical and MPC phases. If every probe
  fails, the client refuses to send the opaque payload rather than applying an
  unverified fallback. The resulting common geometry reserves 64 KiB for call
  metadata and is immutable after the first payload attempt.
- Relay inboxes and outboxes are disk-backed, private, size-bounded and subject
  to resource backpressure.  State retains descriptors and at most the selected
  frame rather than the complete payload.  Completed envelopes expose an
  idempotent absolute-offset reader capped by the exchange bound, so large
  internal consumers need not materialise the ciphertext in R.  There is no
  request-count or rate-based quota.
- The registered typed recipient independently reserves each admitted
  transfer's signed final length against an aggregate per-session capacity
  (`dsvert.typed_blob.spool_max_bytes`, 1 GiB by default). Completed blobs and
  legacy files are charged until consumed. This is byte backpressure, not an
  attempt or request quota; exact duplicate frames do not refresh retention.
- One process-wide byte manager (8 GiB default, operator-configurable) sums
  relay reservations, the authenticated typed `blobs`/`typed`/`typed_source`
  accounting head, the conservative two-direction exact-GC worker reservation,
  and every opened biomedical capsule-source SQLite store's MAC-verified
  `reserved_bytes`. Durable store owners are opaque path-derived hashes and are
  replaced, not accumulated, when the same store reopens. A transaction
  publishes its conservative in-memory head before SQLite commit: interruption
  can temporarily over-reserve, never under-count; reopening verifies the
  persisted row and reconciles rollback/restart without double charging.
  Biomedical capsule-source reservations are reclaimed after, and only after,
  the two pinned final release receipts have been verified. The source SQLite
  transaction atomically writes a MAC-authenticated compaction tombstone,
  deletes source shares/ciphertexts/aggregate and cross-method intermediates,
  and decrements `reserved_bytes`; a restart sees either the complete old state
  or the complete compacted state. The tombstone retains only final-release and
  final-chunk commitment hashes plus a 16 KiB accounting reservation. It blocks
  delayed ciphertexts or retries from recreating private state, while exact
  final-vector replay remains in the separate durable vector store. New stores
  use bounded incremental vacuum (at most 256 pages per compaction); request
  paths never run a full `VACUUM`.
  A single object or operation larger than its durable/local/process policy is
  tagged `[dsvert_resource_oversize:v1]`, has
  `code=resource_oversize`, `retryable=FALSE`, and is terminal end to end. The
  client recognises only that fixed tag and never retries it. Only competition
  with already retained bytes
  emits the fixed `[dsvert_resource_backpressure:v1]` retryable condition.
  This scheduler has no request counter, rate limit, identity penalty or query
  history gate.
- The retained legacy `mpcStoreBlobDS` compatibility path now enters the same
  process-wide admission decision before mutating memory or disk. Completed
  in-memory values, completed disk files and incomplete multipart bytes are
  counted; the first multipart frame conservatively reserves the complete
  public frame-metadata geometry. Distinct keys and distinct sessions therefore
  cannot evade capacity by opening many sparse transfers. An oversized frame,
  frame geometry or cumulative object is terminal and leaves every previously
  accepted chunk intact. Consumption or session cleanup reclaims the bytes;
  direct legacy producers also reserve their full replacement value before a
  temporary-file write or in-memory replacement. This preserves compatibility
  without introducing a request-count limit.
- Lifetime is split deliberately into admission and active leases. A signed
  typed ticket is valid for 24 hours only until its first accepted frame/read;
  exact-GC purpose manifests use the exact-GC TTL only until their first claim.
  These fixed admission deadlines prevent stockpiling old capabilities. Once
  admitted, neither deadline caps total transfer duration.
- Active typed source/recipient transfers enforce a 24-hour inactivity lease;
  the global session lifecycle reaper uses the same inactivity anchor rather
  than total age. Exact-GC uses a server-published, operator-set
  inactivity lease of 180 seconds by default; the client uses the shorter of
  its local monotonic availability setting and the two peers' published
  leases. The relay lease defaults to six hours and remains configurable up to
  24 hours. All are refreshed only by authenticated/durable forward progress.
  Exact duplicate frames, cached response replay, conflicting/invalid frames,
  repeated reads, retry attempts and empty polls cannot prolong retention.
  Thus a healthy transfer may run for longer than any one lease, while an
  abandoned spool is reclaimed after one full idle interval. These controls
  are byte/state retention bounds, not query or request quotas.
- Exact-GC v4 does not let empty relay polling impersonate computation
  progress. The worker emits a private HMAC-SHA256 heartbeat bound to its
  protocol-session nonce, actual PID and strictly increasing counter under a
  one-use 32-byte key removed from the worker config. Only a valid advancing
  counter, committed bytes/ACKs or a terminal transition renews inactivity;
  replayed/rolled-back/misbound heartbeats fail closed. A separate operator
  maximum runtime is a non-renewable total deadline, so a live but stuck or
  malicious worker cannot retain a spool indefinitely.
- Incoming capacity reserves payload plus frame metadata before accepting the
  first byte.  Outgoing streaming reserves capacity before publication, fixes
  frame geometry, validates exact duplicate bytes and seals only after full
  length and SHA-256 verification.
- Control-plane metadata has an independent bound from the disk spool.  Tiny
  frame geometry is rejected before allocating offset vectors, creating a
  spool file, decoding Base64 or verifying a signature; enlarging disk capacity
  therefore cannot induce unbounded R descriptor allocation.
- Absolute ACK offsets are monotonic.  An ACK advances only after a typed,
  operation-bound response; duplicated old frames cannot regress the cursor.
  Terminal receipts are signed and bound to peer, capability, operation, size,
  hash and final offset.
- A source-stream relay verifies each frame hash but does not copy the complete
  ciphertext to a second client file. The recipient computes the complete
  ticket-bound SHA-256 before atomic commit and signs the exact hash and length;
  the pinned producer verifies that receipt before releasing its source spool.
  The removed relay hash file was not a trust anchor—the analyst controls it—
  and cost one full payload write plus one full payload read.
- After the initial source read, that relay commits frame `i` at the recipient
  and reads immutable frame `i+1` from the producer in one named DSI fan-out.
  Missing either result causes an exact replay of both expressions with no
  cursor advance; tests cover a committed store with a lost read response and
  a successful read with a lost store response. The final store and signed
  receipt confirmation remain ordered, so `F` frames use `F + 2` healthy DSI
  phases instead of `2F + 1`.
- Pinned peer identities and capability strings are checked before exchange.
  The analyst/relay can forward peer ciphertext but cannot substitute a peer
  transport key or open the payload under the stated no-host-compromise threat
  model.

Ciphertext authentication does not hide transport metadata.  The relay still
observes public routing, frame count, total ciphertext length, round count,
timing and availability.  Current claims therefore treat the corresponding
geometry (notably aligned row count and method shape) as public protocol
metadata.  A deployment that treats any of it as sensitive must pad to fixed
custodian-declared buckets and use a fixed-round schedule; sticky output noise
alone cannot close a length/timing side channel.

Resource bounds are availability controls, not privacy-budget controls.  They
bound memory/disk amplification and apply storage/DSI backpressure; they do not
limit how many authorised requests an analyst may issue.  Output composition
and reconstruction resistance belong to the single release/privacy policy,
not to the transport scheduler.

The global byte manager is process-wide, not a host-wide distributed disk
quota. Separate Rserve/Opal worker processes must each receive the same policy
and the host/container still needs an operational filesystem quota. The
SQLite store retains its own authenticated durable cap across processes; an
opened store is reconciled into the local process manager before any new
reservation.

## Does dsVert need an application batching layer?

No generic layer is needed. There is no inner Flower envelope, no per-peer
wrapper serialization and no generic application/request batching layer to add
or tune. Independent calls in one protocol phase already travel as one named
DSI fan-out. The capsule-source route described below has a purpose-bound byte
window inside one fetch/accept request; it carries unchanged authenticated
inner artifacts under unsigned outer framing and does not create multiple
outstanding DSI jobs. DSI still pays one
aggregate/fetch lifecycle, HTTP/TLS handling and expression parse per submitted
phase, so dsVert optimises that lifecycle directly. Consequently:

1. independent scalar calls to different sites use one named fan-out rather
   than serial per-peer requests;
2. small scalar arguments travel directly rather than through JSON wrappers;
3. large opaque values still require parser-bounded immutable frame chunking;
   this is transport framing, not application batching. The source-to-recipient
   pump pipelines read `i+1` with commit `i`, giving `F + 2` healthy phases for
   `F` frames;
4. dependent MPC rounds remain separate, because combining them would change
   protocol semantics or make recovery ambiguous.

Increasing a frame beyond the minimum capability of every connector can make a
nominally faster path fail at the parser or proxy. The 480 KiB-raw value is the
starting geometry measured for the exact DSLite profile here; data-free
negotiation may raise or lower it and freezes a common bound only after every
target proves it.

## Data-free frame negotiation

`dsvertTransportProbeDS(nonce, padding, padding_sha256)` is a safe public
control-plane endpoint. It accepts only a fresh syntactically bounded nonce and
printable ASCII chosen by the client, verifies its SHA-256, enforces a server
hard cap of 8 MiB and returns only nonce, observed length, hash, protocol
version and server cap. It does not resolve a protected symbol, create/read an
MPC session, touch a DP ledger or persist state. Its cap is a memory/availability
bound, not a request-count or privacy-budget limit.

The client evaluates one named all-site request map per candidate across the
complete target set, normally descending through 8, 4, 2, 1, 0.625, 0.3125,
0.078125, 0.03125 and 0.015625 MiB of text. A transport exception, callback
failure, absent site or `NULL` may fall through to the next candidate because
no protocol state exists. A response that is present but malformed, misnamed,
wrongly hashed or associated with another nonce fails closed instead of being
interpreted as a smaller capacity. If no candidate verifies, the client fails
before sending any opaque payload; it never applies an untested 640 KiB
fallback. Unlike data-bearing phases, these public requests are synchronous,
so their cold connector cost is sequential in K; the exact-session cache
removes this cost from later transfers on the same login set.

The negotiated payload-text size is at most 8 MiB; the corresponding complete
expression receives a fixed 64 KiB metadata reserve and is still measured
before submission. Server-side legacy and typed blob stores independently cap
one frame at 8 MiB. Cache entries are bounded and keyed to the sorted logical
site set, connector class, sanitised endpoint and DSI session identity; Opal's
key is recomputed after lazy session creation, and Armadillo binds the
authoritative entry to the exact in-process libcurl handle rather than its
refreshable/reusable bearer token. The key combines the pointer address with a
random lifecycle id attached to that exact handle, so allocator address reuse
after finalization cannot inherit cached or poisoned state. Reconnects or a
changed login/profile cannot inherit another path's larger result. A second
bounded cache retains only the last successful connector-profile size as probe
ordering guidance. Every new session must re-prove that size with public
padding; failure descends safely. Thus a warm authenticated endpoint normally
needs one probe without trusting stale capacity. For the exactly benchmarked
DSLite 1.4.1 release, the class/version hint reduces cold-start padding from
15.625 MiB over five calls to 0.625 MiB in one call; another DSLite version is
not hard-coded and uses the ordinary extended ladder.

The first payload attempt locks geometry. From that point, an ambiguous result
may cause any number of byte-identical replays at the same absolute offset
until the monotonic availability deadline; it can never trigger a smaller
chunk, new chunk count or re-probe. Confirmed offsets are skipped on resume,
and a completed recipient can replay its signed terminal receipt from any
previously authenticated frame. Thus negotiation improves
availability/performance without weakening ACK, replay, backpressure or commit
semantics. It is not TLS attestation and no server security decision trusts the
reported capacity.

This replay statement assumes the same live DataSHIELD/R worker session. MPC
state is process-local and its temporary spools live below that worker's
`tempdir()`: an HTTP ambiguity or connector reconnection can resume while the
worker survives, but worker death, service restart or a new login cannot. Such
loss fails closed and the complete purpose-bound operation must restart with
fresh session identifiers; the client does not claim transparent transcript
recovery across process death.

DSOpal 1.5.0 deparses the call and posts it as
[`application/x-rscript`](https://github.com/datashield/DSOpal/blob/1.5.0/R/datashield.aggregate.r#L13-L26).
Neither DSI nor that connector publishes a request-body capability. Opal's
documented [`org.obiba.opal.maxFormContentSize`](https://opaldoc.obiba.org/en/latest/admin/configuration.html#miscelaneous-configuration)
applies to HTTP form posts, so it is not evidence for the DataSHIELD
`application/x-rscript` route. Reverse proxies can impose a still smaller
limit. Consequently 2--8 MiB may be used on Opal only when this exact deployed
path accepts the data-free probe; it is not enabled from connector class or an
analyst option alone.

## Sliding-window decision

A generic 2--4 request sliding window is not enabled. The capsule-source byte
window described below is one bounded request containing consecutive immutable
artifacts, not multiple outstanding DSI jobs. In DSI 1.8,
`datashield.aggregate()` owns the complete lifecycle: it submits at most one
aggregate job per connection, polls it, fetches it and returns only fetched
values.  No job handle escapes to the caller.  Its last-error registry is also
namespace-global, so running several high-level aggregate calls concurrently
would race error attribution.  DSLite explicitly advertises
`dsIsAsync(conn)$aggregate == FALSE`; its `dsAggregate(..., async=TRUE)` executes
synchronously before returning.  A regression test records both facts.

The lower-level DSI generic can return a remote job handle for connectors such
as Opal, but DSI exposes no contract for multiple outstanding jobs on the same
session, their execution order, cancellation, or a negotiated maximum
in-flight count.  Forking live connection objects is not a valid substitute,
and out-of-order execution would break endpoints that commit sequential state.
Therefore a window greater than one would be connector-specific speculation,
not a portable optimisation.

A future window is safe only after DSI/each connector advertises and tests an
explicit multi-in-flight capability.  Each operation must then carry immutable
absolute offsets and hashes, the server must accept out-of-order frames without
changing semantics, ACKs must be independently idempotent, and memory/disk must
reserve at most `window * frame_size`.  That bound is backpressure, not a query
quota, and it does not increase the per-expression byte ceiling.  Until those
conditions exist, the implemented portable optimum is concurrency across peers
inside one named fan-out, with one dependency-ordered frame per peer.

## Scaling model

For a payload of `B` raw bytes per peer and a frame of `F` bytes, the normal
transport requires `ceiling(B/F)` fan-out cycles.  Base64 expands the opaque
payload by approximately 4/3, plus a small expression overhead.  Fan-out removes
the avoidable multiplication of round trips by the number of peers, although
total transmitted bytes remain linear in peer count.  Some multi-party MPC
topologies can additionally require pairwise traffic.

For one producer-spooled source transfer, let `C = ceiling(B/F)`. The pipelined
producer-to-recipient path uses one prime read, `C - 1` paired read/store
fan-outs, one final store and one receipt confirmation: `C + 2` healthy phases.
This changes latency constants, not the unavoidable `O(B)` bytes or work.

Dataset row count does not automatically imply a large DSI transfer: methods
that compute sufficient statistics at the node communicate summaries whose
size depends mainly on features and model iterations.  Protocols that operate
on per-row secret shares still have unavoidable `O(n)` cryptographic work and
traffic.  Disk-backed framing prevents RAM blow-up and makes reconnect/replay
safe, but cannot remove that bandwidth or computation.

Repeated public capsule consumers use a separate in-process release cache.
Every call still refreshes and validates the complete pinned capsule status and
rebuilds the current server-authoritative manifest. Only then may a key over
the sorted site set, stable control-plane state, `manifest_sha256` and
`capsule_id` reuse the already public DP vector and coordinate layout. Thus a
snapshot, workload, pin, policy, noise-root or release-domain change cannot
serve stale numbers. The authenticated reservation count can deny a new
capsule, but does not cause a false miss for exact replay of the existing
capsule. A hit skips source sharing, allocation/sampling and vector replay, not
manifest authority. The LRU holds at most four entries and 64 MiB, strips the
live connection context before storage, never caches failures and never gates
an operation. Eviction merely causes the same durable sticky release to be
validated and replayed again.

Exact COMMIT/RELEASE endpoint replay and a sticky hit in the live session are
O(1). An end-to-end cold reconstruction after process restart is different: it
re-enters through AllocationProof and performs an O(N) audit of the complete
allocator journal before returning the proof. It remains replay rather than a
new release and consumes no lifetime unit, sampler invocation or protected-data
read.

The global constant-memory claim still does not apply to legacy statistical
producers: they return an already-materialised Base64 string in one DSI result
before the client validates its whole SHA-256 and slices parser-safe frames.
Their producer, client/coordinator and current consumer peaks remain
`O(payload)`. Recipient transport writes are `O(frame)` plus a private disk
spool with aggregate per-session byte reservation and backpressure before the
first frame; current statistical consumers may still read the committed object
in full. The global resource flags therefore remain
`producer_source_streaming=FALSE` and `client_source_streaming=FALSE`.
The retained generic/PSI compatibility accumulator is also O(payload), but its
per-frame bookkeeping is now O(1): an environment stores each indexed chunk
and digest once, and incremental byte/count state replaces repeated scans of
all prior frames. The first frame fixes `generic` versus `psi` purpose, so a
cross-purpose continuation cannot relabel an in-progress payload.

One deliberately data-free producer now implements the missing layers without
promoting a statistical capability. `mpcTypedSourceProbeDS` invokes the Go
runtime with a server-minted private path; the runtime writes Base64url output
directly to a bounded spool while computing length and SHA-256 incrementally.
The client receives only a purpose-bound ticket and calls
`mpcTypedBlobReadDS(ticket, offset, max_chars, session_id)`. That reader accepts
no path, key, purpose or recipient override, fixes geometry on the first read,
and permits byte-identical absolute-offset replay. The client verifies every
frame and forwards it to the existing absolute-ACK recipient. The recipient
verifies the complete ticket hash before commit; its signed terminal receipt
is verified by the producer and only then removes the producer spool. The
per-producer resource entry consequently marks producer/client source streaming
true only for `mpcTypedSourceProbeDS`, with
`statistical_producer=FALSE` and `recipient_consumer_streaming=FALSE`.

This pilot is the migration template, not evidence that `.callMpcTool` or any
statistical producer has become streaming. Each real producer still needs a
purpose-specific file-output runtime contract and each large consumer needs a
bounded file/offset input before its own capability flag can change.

The reproducible isolated-process typed-pump matrix is recorded in
`inst/docs/benchmarks/typed_blob_memory_20260801.csv`. With 640 KiB text
frames, 1/16/64/128 MiB inputs used 2/26/103/205 frames and every reconstructed
spool had an identical SHA-256. Maximum-RSS deltas over a separately measured
package-load process were 7,749,632; 52,002,816; 227,999,744; and 293,371,904
bytes respectively. The 64 MiB case used 104 DSI calls and 0.693 seconds in the
network-free mock; these figures confirm the `O(payload)` source claim and are
not remote-DSI throughput results. The fixed-schema byte bound and PCRE byte
validation keep the local pump fast enough that a real deployment should be
network/parser-bound rather than regex-bound.

Reproduce that measurement on macOS or Linux with:

```sh
Rscript inst/scripts/benchmark_typed_blob_memory.R \
  --sizes-mib=1,16,64,128 \
  --output=inst/docs/benchmarks/typed_blob_memory_20260801.csv
```

The committed typed-memory CSV SHA-256 is
`d32e214dcfca12e528062d08edd3cafdcb6b23e23d98b1bc8c6557c739de06c2`.

The new pilot's persistent three-process harness is recorded in
`inst/docs/benchmarks/typed_source_streaming_multiprocess_20260801.csv`. It
runs the real Go file-output producer and client frame pump with 512 KiB text
frames, discards one generation response, one source-frame response and one
recipient ACK, and closes/reopens both worker sockets halfway through. The
64/128/256 MiB cases completed 128/256/512 frames with identical terminal
SHA-256 values at 16.66/16.02/19.43 MiB/s under the recorded concurrent local
load. In the 128 MiB case, wall time was 7.990 seconds; exactly one source
generation occurred, both lost frame
responses caused one exact replay, one reconnect completed, over-capacity
preflights created no file, and the producer source was removed only after the
receipt. The Go producer peak RSS was 6.5 MiB. The controller created zero
payload-sized spool bytes in all three cases, avoiding an exact 2B of redundant
write-plus-final-hash-read I/O: 128/256/512 MiB respectively. R high-water
deltas remained frame-bounded rather than allocating the payload at the relay:
controller 81.2/82.1/96.9 MiB, producer 57.3/69.3/77.8 MiB and recipient
53.7/72.0/75.0 MiB. These include the R runtime/allocator high-water and show
bounded frame processing, not a 512 KiB total-process footprint. Wall time is
environment-sensitive; the security/resource claim rests on hashes, offsets,
process separation and spool bytes, not on these loopback throughput values.

Reproduce the multiprocess pilot after building the current runtime:

```sh
Rscript inst/scripts/benchmark_typed_source_streaming.R \
  --binary=/absolute/path/to/dsvert-mpc --payload-mib=128 --frame-kib=512
```

The committed multiprocess CSV SHA-256 is
`6452905e343c16f0f5db85bb619fd8242567457fdce9c0e7149098b10f5bcbf3`.
The harness uses persistent loopback sockets to prove process isolation,
reconnection and memory geometry; it is not an Opal/Armadillo HTTP/TLS latency
measurement.

A post-lease-change 256 MiB smoke run on 2026-08-01 again completed all 512
frames with one producer replay, one recipient replay, one reconnect, bounded
preflight rejection, source removal after the receipt and terminal
`sha256_identical`. It measured 25.884 seconds (9.89 MiB/s) under concurrent
test load; controller/producer/recipient RSS deltas were 49.3/70.5/84.1 MiB.
This deliberately noisy rerun is validation rather than a replacement for the
committed three-size matrix.

Dense second-order models eventually scale with `O(p^2)` score/Hessian traffic
and storage and approximately `O(p^3)` coordinator solves.  Very high-dimensional
workloads therefore need a separate matrix-free/sparse estimator (for example,
Hessian-vector products with iterative solves) rather than ever-larger DSI
messages.  That is a statistical/MPC method change and is not claimed as
implemented by this transport work.

Practical production scaling should use compiled/vectorised node kernels,
bounded row blocks for per-row MPC, immutable checkpointed offsets, concurrent
site fan-out and connector-specific negotiated frame sizes.  Backpressure must
slow producers instead of accumulating unbounded R objects or temporary files.

## Reproducible benchmark

The focused dispatch A/B is recorded in
`inst/docs/benchmarks/dsi_dispatch_ab_dslite_20260801.csv`. After one excluded
warm-up per engine, four timed pairs alternated which engine ran first. Every
case used two peers, 1 MiB per peer, 480 KiB raw frames, four fan-out cycles,
one deliberately replayed response and one reconnect. Submitted raw, Base64
and expression bytes, frame counts, spool bytes and terminal SHA-256 were
identical between engines. Standard `DSI::datashield.aggregate()` had a median
wall time of 12.264 s and useful throughput of 0.164 MiB/s; the direct immutable
transport path measured 8.825 s and 0.234 MiB/s. That is 28.0% lower median wall
time and 43.3% higher median useful throughput in this DSLite run, without
changing the wire payload or ACK semantics. Individual timings varied widely,
so this is evidence of removed local overhead rather than a portable speed
guarantee. The committed A/B CSV SHA-256 is
`8da0fee25a575ef7f98dbea82bb8b09a516f1dacf92560fad85cf7331a02acb7`.

Reproduce the paired case with:

```sh
Rscript inst/scripts/benchmark_dsi_transport.R --quick \
  --peers=2 --sizes-mib=1 --chunk-raw-kib=480 \
  --engines=direct,dsi-standard --balanced-ab --warmup --repeats=4 \
  --output=inst/docs/benchmarks/dsi_dispatch_ab_dslite_20260801.csv
```

This optimization removes coordinator-side formatting/copies. It cannot remove
the connector's necessary one-time serialization, Base64 expansion, HTTP/TLS
round trip, server parsing or statistical/MPC computation.

Run the current direct engine on the full grid without overwriting the archival
baseline:

```sh
Rscript inst/scripts/benchmark_dsi_transport.R \
  --engines=direct \
  --output=inst/docs/benchmarks/dsi_transport_direct_dslite.csv
```

The archival scaling matrix below predates the direct-dispatch optimization
and is retained as the standard-DSI baseline. It covers 2, 4 and 8 peers with
1, 16 and 64 MiB per peer. Each case
uses a 480 KiB raw frame, one Base64url encoding and one named DSI fan-out per
cycle.  It deliberately drops one complete response and reconnects all DSLite
connections halfway through.  Absolute offsets are replayed, and every peer's
disk spool must finish with a SHA-256 identical to the streaming source.

The CSV records wall time and useful throughput together with cycles, replay,
reconnection, raw/Base64/expression bytes, maximum expression size, process RSS,
source/spool disk use and terminal integrity.  DSLite may execute site work
serially in-process, so its peer-scaling numbers are a conservative local
implementation measurement rather than evidence of network concurrency.

Recorded results (useful bytes exclude the deliberately replayed frame):

| Peers | MiB/peer | Useful MiB | Fan-out cycles | Wall seconds | Useful MiB/s | RSS delta MiB | Spool peak MiB |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 2 | 4 | 10.224 | 0.1956 | 88.31 | 2 |
| 4 | 1 | 4 | 4 | 14.812 | 0.2701 | 29.56 | 4 |
| 8 | 1 | 8 | 4 | 32.264 | 0.2480 | 18.25 | 8 |
| 2 | 16 | 32 | 36 | 80.751 | 0.3963 | 0.00 | 32 |
| 4 | 16 | 64 | 36 | 243.485 | 0.2628 | 13.56 | 64 |
| 8 | 16 | 128 | 36 | 837.440 | 0.1528 | 0.00 | 128 |
| 2 | 64 | 128 | 138 | 1218.312 | 0.1051 | 14.30 | 128 |
| 4 | 64 | 256 | 138 | 1060.728 | 0.2413 | 40.25 | 256 |
| 8 | 64 | 512 | 138 | 1148.725 | 0.4457 | 16.23 | 512 |

All nine cases recorded exactly one ambiguous-response replay, one connection
rebuild, a maximum deparsed expression below 641 KiB and terminal
`sha256_identical`.  Base64 expansion was 1.333334 or lower before the small R
expression overhead.  The committed CSV SHA-256 is
`65cd02b1729e61e4ebafdcd171db542a5c71dd226e967bcc50b0747216a7d659`.
RSS deltas are sampled process high-water observations and should not be
compared as monotonic per-case allocations; cases run in one R process and R's
allocator retains memory.  The non-monotonic DSLite throughput likewise must
not be extrapolated to remote connector performance.

### Frame-size sweep

Run:

```sh
Rscript inst/scripts/benchmark_dsi_transport.R \
  --peers=1 --sizes-mib=16 \
  --chunk-raw-kib=512,1024,2048,4096 --continue-on-error \
  --output=inst/docs/benchmarks/dsi_transport_frame_sweep_dslite_20260731.csv
```

| Raw frame | Base64/expression scale | Result | Cycles | Seconds | Useful MiB/s |
|---:|---:|---|---:|---:|---:|
| 512 KiB | 682.95 KiB measured | SHA-256 identical | 33 | 63.619 | 0.2515 |
| 1 MiB | about 1.33 MiB | parser rejected | -- | 3.071 | -- |
| 2 MiB | about 2.67 MiB | parser rejected | -- | 7.256 | -- |
| 4 MiB | about 5.33 MiB | parser rejected | -- | 13.887 | -- |

Rejected candidates did not advance an ACK or commit payload bytes. The sweep
CSV SHA-256 is
`0458550ba3d948997d53d2ffa7603fd3790bb4d4b1a938974f95978a67ba4528`.
It shows why a global 1--4 MiB raw default would break DSLite, while the
data-free endpoint integration test descends to and verifies 640 KiB of text.
On a path that verifies the 8 MiB text candidate, a 64 MiB raw payload needs
about 11 normal payload cycles instead of 137; 4, 2 and 1 MiB text require
about 22, 43 and 86 respectively. Real Opal throughput still requires the
multi-host validation below.

## Legacy blob migration inventory

`dsVert/inst/docs/typed_blob_transport_inventory.json` is the machine-readable
map that must reach zero legacy consumers before `mpcStoreBlobDS` can be
unregistered.  It records 20 semantic transfer classes, each with family,
caller, destination slot, producer, recipient, payload shape, confidentiality
and migration stage.  A class can be shared by several method families, so the
family counts below intentionally overlap:

| Family | Transfer classes |
|---|---:|
| chi-square cross-site | 2 |
| correlation | 7 |
| Cox | 10 |
| GEE | 8 |
| GLM | 11 |
| GLMM | 6 |
| LMM | 3 |
| multinomial | 8 |
| negative binomial | 10 |
| ordinal | 8 |

The common replacement substrate is registered as
`mpcTypedBlobStoreDS(ticket, chunk, offset, session_id)`.  It deliberately has
no caller-selected key, path, recipient or purpose.  A server producer mints a
signed one-shot ticket binding the session UUID, random transfer ID,
allowlisted capability, pinned sender and recipient, complete identity and
transport-key manifest, validated shape context, payload length and SHA-256.
The recipient derives its storage slot from the server registry.  Absolute
ASCII offsets, exact duplicate checks, private disk spooling and atomic final
commit provide bounded, idempotent delivery.

The substrate is ready and adversarially unit-tested. The first-wave
input/gradient/weight/Beaver/IKNP/comparison payload classes now use
producer-minted typed tickets and recipient receipts. The unused mutable
correlation controls have no product caller and are no longer registered or
exported. The machine-readable inventory and source regressions enforce those
migrations. The generic endpoint remains only
for explicitly inventoried second-wave/quarantined callers; it can be removed
after those consumers are redesigned or retired and a source/surface test shows
no reachable use.

## Capsule-source byte window

The canonical Ring128 capsule source keeps its 8,192-coordinate encrypted
chunks, but may carry a consecutive prefix of those unchanged chunks through
one DSI fetch and one two-recipient acceptance phase. This is not a larger
cryptographic chunk and does not change a share, ciphertext, signature, hash,
ACK or DP draw. The wrapper is capped at 768 KiB and eight chunks; the actual
prefix stops before the byte cap. Those are per-message resource bounds, never
request or history quotas.

Rolling negotiation is fail-safe. The client first obtains the persistent
scalar v1 ticket. For capsules with at least two chunks it asks both designated
recipients synchronously for a signed wrapper around that exact ticket, then
asks every declared source owner synchronously for a signed capability bound to
the capsule and source contract. A missing method/argument, oversized optional
request or response, partial result, or old source selects the scalar path
before `PrepareDS` receives a wrapper. Identity/pin failures and interrupts
still propagate. This fallback requires evidence of a terminal application
response. A generic transport failure after a synchronous `execute` can leave
Armadillo's `/lastcommand` lifecycle ambiguous; the affected authenticated
handles are therefore marked unusable and the client requires fresh login
connections without issuing `PrepareDS`, `ChunkDS`, or a scalar retry on those
handles. The old-server path is exercised on real local DSLite, while the
ambiguous two-step Armadillo path is fault-injected against the installed
connector contract; it is not a live Armadillo-server E2E claim.

`experiments/capsule-source-window/benchmark.R` uses data-free production-like
132,000-byte ciphertexts per recipient. One canonical paired bundle measured
355,057 bytes, so two fit below the independent 786,432-byte response cap. At
that effective width, data phases change from `2*S*C` to
`2*S*ceiling(C/2)` plus two one-time public negotiation phases. Server
`dsAggregate` invocations change from `3*S*C` to
`3*S*ceiling(C/2) + 2 + S`: for `S=2,C=123`, 738 becomes 376; for
`S=3,C=123`, 1,107 becomes 563; for `S=5,C=123`, 1,845 becomes 937. Canonical
encoding of sixteen production-like bundles was about 0.08 seconds scalar and
0.15 seconds with wrappers on the recorded local run; this small CPU increase
buys the much larger reduction in remote DSI round trips and is not presented
as a network throughput benchmark.

Request and response geometry are tested separately. A 480 KiB raw Base64url
request remains below the portable expression parser boundary. A distinct real
DSLite regression returns a 768 KiB scalar response byte-for-byte; passing that
test does not raise the request-expression limit or claim equivalent Opal proxy
limits.

## Remaining production validation

- Repeat the same matrix on at least two real processes and hosts using every
  supported DSI connector, verified TLS and realistic proxy limits.
- Inject connector timeout, partial response, process restart, stale ACK,
  conflicting replay, disk-full and spool-compaction failures under load.
- Record p50/p95/p99 phase latency and node CPU in addition to coordinator RSS.
- Benchmark representative statistical protocols separately: row count,
  feature count, peer count and iteration count affect compute and transport in
  different ways.
- Migrate legacy whole-string statistical producers and consumers one
  capability at a time using the bounded file/offset pilot, then verify
  coordinator and node RSS under a large-`n` secret-share workload before
  changing that capability's resource flag.
- Treat the quarantined LMM/GEE/GLMM client orchestrators as legacy until their
  cluster-granular release contracts are replaced; the LMM raw DSI transport
  calls have been migrated, but that does not alter their release semantics.

## Raw-call inventory after migration

All active provisional routes and the retained GEE/GLMM/LMM orchestrators now
enter DSI through the strict aggregate/assignment, named fan-out, bounded blob
(typed where its producer exposes a typed transfer) or
best-effort-final-cleanup helpers. The central strict aggregate helper itself
dispatches real `DSI::datashield.aggregate` calls through the direct public-DSI
job primitives; dependency-injected test transports retain their original
call contract. The LMM transport-only migration removed its
60 direct runtime `DSI::datashield.aggregate()` references without changing the
public API or estimator:

| Legacy file | Direct runtime references after migration | Product status |
|---|---:|---|
| `R/ds.vertLMM.R` | 0 | quarantine |
| `R/ds.vertLMM.k3.R` | 0 | deprecated quarantine |
| `R/ds.vertLMM.closed_form.R` | 0 | quarantine/internal |

Independent two-peer phases now use named fan-out, repeated column discovery
inside each K=2/K>=3 implementation is one fan-out reused throughout that
implementation, and final cleanup uses the
only best-effort transport contract. Identity mismatch and unresolved DSI
session conditions are propagated unchanged, including through Beaver-policy
discovery. No LMM operation gained a retry unless its producer was already on
the global typed-idempotent allowlist.

This is deliberately not a statistical promotion. The LMM family still retains
cluster-granular outputs (including cluster sizes and reconstructible paired
shares of per-cluster residual moments), so routing its calls through a strict
transport boundary does not make the release contract disclosure-safe. It
remains quarantined pending a bounded, capsule-level redesign and validation.

Other textual references to `DSI::datashield.aggregate` are dependency-injected
default function arguments, wrapper implementations, tests/comments, DP
front-door plumbing or exact-GC transport injection points; they are not
unwrapped statistical protocol calls.
