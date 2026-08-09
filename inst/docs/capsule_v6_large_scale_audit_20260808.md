# Historical capsule v6 large-scale audit (2026-08-08)

## Scope and notation

This frozen audit covers the then-current capsule-v6 path from protected
snapshots to the final public vector. The active workload is v7 and adds the
authenticated bounded-lifetime contract; this file is performance/geometry
evidence, not an artifact-exact claim for the current release. Unless stated
otherwise, its byte and complexity observations remain descriptive rather than
current lifecycle/accounting assertions.

- `n`: protected rows at one source site.
- `p`: protected columns used by the declared workload at that site.
- `U`: public patient/unit capacity.
- `d`: final public release coordinates (`d <= 1,000,000`).
- `dt`: source-transport coordinates. Without a cross-owner artifact,
  `dt = d`; padded cross-owner layouts permit `dt <= 64 * 1024^2`.
- `S`: source owners participating in the capsule (`S <= K`).
- `K`: pinned peers.
- `B = 8,192`: productive source and final-vector chunk size.
- `Cs = ceiling(dt / B)` and `Cf = ceiling(d / B)`.
- `q`: columns in one within-owner Gaussian design, including an intercept.
- `J`: certified binary-geometric bits (`1..63`) and `W`: sampler word size
  (`128` or `256` bits).

The productive discrete-Laplace selector currently promotes exact GC only for
`d <= 1`. For `d > 1`, it selects the two-independent-complete-draw Ring128
convolution backend with `B = 8,192`; consequently `Cf <= 123`. The exact-GC
chunk ceiling of 128 is not a large-vector production geometry today.

## Complexity inventory

| Stage | Time | Peak/live storage | Scale observation |
|---|---:|---:|---|
| Snapshot validation and digest | `O(n p)` | snapshot plus one column block | Hashing is block based; it does not serialize the complete frame at once. |
| Admission, once per dataset/materialization | `O(n log n)` | `O(n + U)` | Canonical patient IDs are sorted, matched, and counted. |
| One bounded numeric variable | `O(n log n)` | `O(n + U)` transient, `O(U)` compact cache | Sorting by unit/value dominates. Across `p` numeric variables: `O(p n log n)` and up to `O(p U)` retained compact values. |
| Numeric moments | `O(U)` per variable | `O(U)` | Fixed three-coordinate output. |
| Numeric pair moments | `O(U)` per requested pair | `O(U)` | Six coordinates per pair. |
| Within-owner Gaussian statistics | `O(U q^2)` arithmetic after patient collapse | `O(U q + q^2)` | Output has `1 + q(q+1)/2 + q + 1` coordinates. The canonical lower-triangular output is preallocated, so coordinate assembly is linear in the `O(q^2)` emitted values; the Gram computation remains `O(U q^2)`. |
| Categorical marginal/pair | `O(n + U + cells)` per artifact | `O(n + U + cells)` | Pair domains and output cells are public and bounded by the coordinate contract. |
| Survival artifact | `O(n log n + U + T C)` | `O(n + U + T C)` | Patient-row selection sorts; `T` is public time-grid length and `C` outcome count. |
| Full local release material | `O(d)` | at least `8d` bytes for the R numeric vector | At `d = 10^6`, the bare vector is about 8 MB. Preparation/validation can temporarily hold another range vector. |
| Cross private source producer | `O(dt)` plus bounded-variable work | `O(d + U + B)`, not `O(dt)` | The private padded tail is read in ranges; the complete `dt` vector is not materialized. The public release prefix remains a full `d` vector. |
| Local convolution sampler | `O(d J W/64)` fixed-limb comparisons per designated peer | `O(B + J W/8)` | Each peer samples two geometrics per coordinate; both peers together account for four. Each peer expands `2 J W/8` private stream bytes per coordinate. |
| Exact GC | `O(m (512 + 4 J W))` circuit input bits for a chunk of `m` coordinates, plus gates/OT | circuit/transcript bounded to `m <= 128` | Currently used only at `d = 1`. Large-vector promotion cannot be inferred from the small scalar route. |
| Finalizer | `O(d)` | `O(B)` working data plus durable `O(d)` final chunks | Ring128 signed decode and fixed clamp are per coordinate. |
| Client replay and numeric conversion | `O(d + Cf)` after this change | public `O(d)` character vector and `O(d)` numeric result coexist during conversion | A caller requesting `d` coordinates necessarily receives an `O(d)` result. |
| Method post-processing | generally `O(block length)`; Gaussian solve `O(q^3)` and `O(q^2)` memory | copies the selected block/matrix | Repeated block extraction also revalidates the full `d` vector today, an `O(d)` scan per extraction (residual P2). |

The source producer computes the release commitment during preparation and
reconstructs the producer on the first uncached chunk. For a normal release
segment, that first chunk generates and durably stores every release chunk;
subsequent source fetches are SQLite replay rather than repeated protected-data
scans. This is approximately two materialization passes per new source
release, not `Cs` passes.

## DSI calls, bytes, and client memory

Recipient ticket creation and source preparation are asynchronous fan-outs.
Source delivery preserves canonical order but now negotiates a public,
byte-bounded consecutive window through the existing ticket, chunk and accept
methods. Let `Ws` be the effective window for one source (`1 <= Ws <= 8`). It
is derived before protected source access from the minimum of that source's
data-free response probe, both recipients' per-site request/expression probes,
the signed application capability and fixed hard caps. The response probe
keeps 12.5% plus 128 KiB of headroom, and window framing keeps another 64 KiB;
a request ceiling is never treated as evidence about a response. For each
source window the client performs one fetch followed by one two-recipient
acceptance fan-out. Excluding the fixed public capability negotiation, source
delivery costs:

The adaptive capability is a distinct v2 negotiation. Within one synchronized
DSV1 client/server package generation, the adaptive v1 capability and omitted
adaptive method arguments remain byte-exact and hard-limited to 768 KiB for
source responses and 1 MiB for recipient requests. The client adds the public
v2 adaptive contract only to a multichunk call after all signed v2 attestations
agree. DSV1 itself is a deliberate lockstep wire change: client and server
packages must be upgraded together, and a mixed package generation fails at
the first framed, data-free manifest phase before protected data access. The
safe staged rollout is server-first with release traffic paused: deploy,
reconcile and attest every server, verify this data-free boundary, and then
deploy the paired client. The intermediate mixed state is fail-closed downtime,
not a compatibility mode.

- `2 sum_s ceiling(Cs / Ws)` sequential DSI latency phases;
- `3 sum_s ceiling(Cs / Ws)` AggregateMethod invocations (one fetch and two
  accepts);
- `Theta(S dt)` bytes.

Ignoring small AEAD and JSON headers, one source has two 16-byte Ring128 shares
per coordinate. Fetching their Base64 bodies and then sending those same two
bodies to the recipients is approximately
`2 * (4/3) * 32 * dt = 85.33 * dt` application bytes through DSI. The common
8,192-coordinate shape has a 355,057-byte production-like canonical pair
bundle in the committed data-free benchmark. An 8 MiB response probe yields a
conservative 7,208,960-byte usable bound and therefore the hard maximum
`Ws=8`. A smaller or heterogeneous deployment chooses its public per-source
width independently; absence of the response extension uses `Ws=1` without
changing artifacts or correctness.

| Shape | `Cs` | Source relay phases (`Ws=8`) | Approx. source payload |
|---|---:|---:|---:|
| `d=10^6`, `S=2` | 123 | 64 | 171 MB |
| `d=10^6`, `S=10` | 123 | 320 | 853 MB |
| `dt=64*1024^2`, `S=10` | 8,192 | 20,480 | 57.3 GB |

These figures exclude HTTP/DSI framing, signatures, manifests, final-share
exchange, and public replay. They are lower-order for large cross layouts but
are material for ordinary `K=2` releases too. Direct scalar DSI calls remove an
extra application serialization layer; they do not make fixed chunks obsolete.
Chunks remain necessary for parser bounds, authenticated retries, durable
offsets, and backpressure.

The client retains only the current source bundle, so private source relay
memory is `O(B)` rather than `O(dt)`. It later retains the final public
character vector and numeric vector, both `O(d)`. Response fan-out also means
both designated peers return the same public replay chunk. The 64 MiB public
release cache simply declines to cache an oversized result; it does not block
the release.

## SQLite and spool behavior

Source state uses authenticated canonical JSON rows, WAL, and
`synchronous=FULL`. Its conservative public reservations are:

- recipient key/ticket per designated peer: three times its canonical row
  bytes plus 64 KiB of page, index and transaction overhead;
- outbound per source: `64 dt + 128 KiB Cs + 64 KiB`;
- inbound per designated peer: `24 dt + 16 KiB S Cs + 32 KiB`.

For `dt=10^6`, outbound is about 80.2 MB and inbound is about 28.1 MB at
`S=2` or 44.2 MB at `S=10`. For the maximum cross layout, the corresponding
figures are about 5.37 GB outbound and 1.88 GB (`S=2`) or 2.95 GB (`S=10`)
inbound. Source intermediates are compacted only after bilateral durable
publication; the compacting tombstone prevents recreation. The authenticated
inbound reservation now scales with every durable source/chunk receipt
`O(S Cs)`, so large consortia encounter byte backpressure before persistent
mutation rather than silently exceeding the source-store accounting. This is
not a request, session or release-history quota. Recipient tickets are also
reserved atomically before insertion, replay does not add a second reservation,
and authenticated legacy rows migrate under the same store lock.

Each designated vector store retains one Base64 noised share (`~21.33d`
characters) until finalization and enforces a 4 GiB intermediate cap. Final
public chunks are retained for sticky replay and grow `O(d)` per distinct
published release. Source and sampler intermediates are reclaimed. Under the
current v7 contract, at most `N` capsules and one public instance per capsule
bound this active-release growth by `O(N d)` for a fixed maximum shape; abandoned
allocator commits consume lifetime units without adding a public vector. An
authenticated archive/retention tier may still be operationally useful, but it
must preserve replay and rollback evidence.

The finalizer repeatedly opens/locks/initializes SQLite inside its chunk loop.
This is bounded by `Cf <= 123` today but is an avoidable P2 constant. Replay
also rebuilds the complete Merkle tree for each requested chunk: a complete
replay is `O(Cf^2)` server hashing. At the current maximum this is only about
`2 * 123^2 ~= 30,000` node/leaf hashes; with hypothetical 128-coordinate
large-vector exact-GC chunks it would become roughly 122 million hashes.

## Implemented byte-identical improvement

The client previously assembled replay with `scaled <- c(scaled, chunk)`.
For chunk sizes `b_i`, that copies approximately
`sum_i b_i * (Cf - i + 1)` vector references: `O(d Cf)`, quadratic in `d`
for a fixed chunk size. The new internal collector preallocates exactly `d`
character slots, validates every public chunk shape, and assigns it to the
signed offset. It changes neither bytes, order, network calls, API, DP, nor
release/cache identity.

Reproducible benchmark (`Rscript --vanilla
experiments/vector-replay-assembly/benchmark.R`, median of three runs):

| `d` | chunk | chunks | append | preallocated | estimated append refs | preallocated refs | identical |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 50,003 | 128 | 391 | 0.047 s | 0.001 s | 78,474,904 B | 400,024 B | yes |
| 100,003 | 128 | 782 | 0.185 s | 0.001 s | 313,499,928 B | 800,024 B | yes |
| 200,003 | 128 | 1,563 | 0.783 s | 0.003 s | 1,251,599,896 B | 1,600,024 B | yes |
| 1,000,000 | 8,192 | 123 | 0.526 s | 0.012 s | 499,716,608 B | 8,000,000 B | yes |

The reference-byte columns estimate vector-pointer traffic on a 64-bit R
build; they do not claim whole-process peak RSS. Network and cryptography are
excluded deliberately so the assembly regression is reproducible.

Relevant regression checks after the change:

- dsVertClient source/vector/exact-GC suites: 307 assertions;
- new replay-assembly suite: 4 assertions;
- dsVert materializer/source/cross suites: 675 assertions;
- dsVert vector/allocation/release-ledger/exact-GC suites: 531 assertions.

Total: 1,517 relevant assertions, with zero failures, errors, warnings, or
skips. This is not a real Opal/Armadillo throughput benchmark.

The adaptive window implementation was additionally checked with focused
client and server source-transport suites for `K=2,3,5`, effective widths one
through eight, heterogeneous node ceilings, adaptive-v1/v2 peers within one
synchronized DSV1 generation, exact
scalar/window hash equivalence, replay, reordering, tamper, oversize, HTTP 413,
reconnect/cache invalidation and ambiguous Armadillo failure injection. A real
local DSLite connector crossed the separate request/response limits. The
committed one-million-coordinate call-count artifact is
`inst/docs/benchmarks/capsule_source_adaptive_window_20260808.csv`. Its
`data_plus_incremental_v2_negotiation_*` columns add only the two incremental
v2 negotiation phases to the data phases; they deliberately exclude the
legacy ticket, protected Prepare, and outer status phases shared by the
comparison. The 492-to-64, 738-to-96, and 1230-to-160 reductions are explicitly
data-phase counts. These checks do not constitute a live multi-host Armadillo
throughput benchmark. Here `K=2,3,5` is topology regression coverage, not a
claim that those topologies ran on current live Opal or Armadillo artifacts.

## Residual priorities

### P1 residual: source-relay bytes

The bounded window removes repeated DSI round trips without changing any inner
ticket, ciphertext envelope, acknowledgement, hash or DP artifact. Both
designated recipients and every source owner must attest the same capability;
same-generation peers lacking adaptive v2 use the scalar route before
protected preparation; mixed DSV1 package generations fail at the earlier
framed manifest boundary. A transport
failure with ambiguous Armadillo singleton state poisons that handle and
requires a fresh login rather than falling back on the same connection. The
outer window is untrusted framing; every inner artifact remains independently
signed and context-bound.

Latency is now `O(sum_s ceiling(Cs/Ws))` with `Ws` derived from independently
verified request and response paths, but the analyst relay still carries
`Theta(S dt)` bytes. Very large padded cross-owner layouts therefore still need
a producer-spooled source-to-recipient stream (or a direct authenticated
server-to-server channel) to improve more than the constant window factor.
Such a redesign would require:

1. canonical source/chunk order unchanged inside each stream;
2. per-envelope signatures and hashes unchanged; outer framing must never be
   treated as an authenticated artifact;
3. atomic recipient commit with absolute offsets/ACK bitmap and byte-identical
   retry behavior after every crash point;
4. decode limits before allocation, bounded producer/recipient spool, and TCP
   backpressure; never aggregate all `dt` values in R;
5. corruption, duplicate, reordering, reconnect, and partial-write tests for
   `K=2,3,5+`;
6. real Opal and Armadillo E2E benchmarks over `d={8K,100K,1M}`, representative
   cross `dt`, and `K={2,5,10}`;
7. unchanged capsule/source/final hashes and unchanged DP output distribution.

### P2: public replay and durable-store constants

A versioned bounded `ReplayRangeDS` response with a Merkle multiproof (or an
authenticated streamed public-vector file) would reduce `Cf` replay rounds to
`ceiling(Cf/window)` and construct the Merkle tree once in `O(Cf)`. A persistent
authenticated node table is another option. Either requires retry/root/storage
tests and is not a surgical change to the current protocol.

One store handle/transaction scope per server-side finalization would remove
repeated SQLite open/schema/WAL work. It must retain the current per-chunk
durability checkpoints and crash recovery before promotion.

### Gate for large-vector exact GC

Do not raise the current `d <= 1` exact-GC promotion threshold by changing a
constant. At worst `J=63`, `W=256`, the circuit input budget is about 65,024
bits per coordinate before gates and OT. A one-million-coordinate vector at
128 coordinates per circuit would require 7,813 circuits. Promotion needs
measured compile time, garbling/OT bytes, peer RSS and spool, resumable durable
chunk semantics, deterministic retry tests, and an end-to-end utility/latency
comparison against the certified convolution route. A failed exact route must
not silently switch backend after private material is accessed.

## Defensible scaling claim

For signed workloads within the public coordinate and resource contracts, the
current route is bounded-memory at the client and for private cross-source
generation, has linear cryptographic/data volume, and never silently changes
backend or chunk geometry. It is not valid to claim constant latency in `K` or
arbitrary-dataset scalability: the windowed source relay still carries bytes
linear in `S * dt`, final output is necessarily `O(d)`, and durable release
history needs finite operational storage or an authenticated archival policy.
