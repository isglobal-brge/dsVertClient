# Single disclosure-safe mode: implementation and acceptance plan

Status: historical implementation and acceptance plan, updated 2026-08-08.
The authoritative current guarantees are the runtime method registry,
`SECURITY.md` and the installed product-surface document.  Unchecked targets
below are not claims about the current development tree.

## Target outcome

The released packages will expose one production behaviour, named
`disclosure_safe`. There will be no analyst- or custodian-selectable legacy
exact mode. Public client method names are retained, but low-level server
primitives are not part of the compatibility contract.

Every statistical result must be one of:

1. a fixed-shape, ledger-accounted DP release;
2. deterministic post-processing of an already authorised DP release; or
3. public protocol metadata that is independent of protected records.

MPC protects intermediate values and pinned peers prevent analyst key
substitution. Neither property makes an unrestricted exact output catalogue
safe. Exact scores, Hessians, cluster summaries, tables, counts, matched-index
sets and unnoised model outputs must therefore never cross the release plane.

The final claim is deliberately narrower than “reconstruction is logically
impossible”. For the declared adjacency, public invariants, clipping bounds,
privacy parameters and non-collusion threshold, released values and successful
authenticated semantic-message contents are either fixed by the public
contract or are post-processing of the DP vector. The R/DSI host path is not constant-time:
pre-release failure occurrence, timing, availability and traffic flow are not
inside that package-level DP claim. A deployment that needs those observables
covered must add an independently validated scheduled/asynchronous gateway.
Host compromise, loss of every noise key, collusion of all designated noise
peers, and data released through another system also remain outside the claim.

## Threat and trust contract

- The analyst/client/relay is untrusted. It may issue adaptive requests,
  reorder or replay messages, and inspect DSI sizes, timing and failures.
- Peer identities and destinations are pinned by the custodian. A request
  cannot substitute a peer, recipient, operation, dataset or result slot.
- At least one of the designated two noise/MPC peers must not collude with the
  analyst or the other peer. All-peer collusion is out of scope and is stated
  in every certificate.
- Peers are initially honest-but-curious. Malicious-peer robustness is a
  separate backend capability and must not be inferred from pinning.
- A protected dataset is an immutable logical snapshot. A mutation requires a
  new snapshot identifier and a newly composed release.
- External TLS verification, peer pinning and origin attestation are
  deployment requirements. Analyst-supplied attestation is not evidence of
  transport security.

## Two explicit adjacency contracts within the one safe mode

These are privacy-policy choices, not safe/unsafe execution modes.

### `replace_one_fixed_cohort`

The cohort cardinality and row capacity are declared public invariants. One
neighbour replaces all protected values for one patient. This avoids padding
when cohort membership/cardinality is already public, but it does not claim to
hide participation in that fixed cohort.

### `add_remove_patient`

Participation is protected. The policy declares a public row-capacity bucket;
the real row count remains inside MPC, and dummy rows plus a secret validity
mask make message count, payload shape and computation shape independent of
the real count within that bucket.

For repeated measures, one patient/cluster is the adjacency unit. The policy
also fixes `max_records_per_unit`, contribution norms and a public cluster/row
capacity. Records beyond those limits are rejected or clipped according to the
published estimand; they are never handled by an implicit fallback.

## Release state machine

1. **Canonical plan.** The client constructs a schema-validated descriptor of
   method, formula, public tuning parameters, public domains, requested output
   and logical snapshot. Equivalent syntax maps to one canonical query ID.
2. **Consortium validation.** Every peer checks the descriptor against its
   immutable policy, dataset binding, alignment attestation, pinned peer set,
   adjacency and public shape. Exact row hashes and counts are not returned.
3. **Durable allocation.** The two custodian-designated, identity-pinned
   allocator peers cross-sign and durably bind
   `(domain, cohort, snapshot, epoch, query_id)` to one global allocation index
   before any statistic can be released. At least one non-colluding designated
   peer must retain its ledger. A linearizable external CAS may additionally
   anchor the allocation and delivery heads, but it is hardening rather than a
   second allocator or an administrator prerequisite.
4. **Opaque capability.** Each peer mints an authenticated capability bound to
   the query, allowed state transitions, destination peers, public shape and
   result slots. The analyst can relay it but cannot broaden it.
5. **Fixed-shape computation.** Purpose-bound MPC endpoints accept only that
   capability. Intermediate replies are authenticated opaque bytes or
   fixed-shape public receipts; they are never exact statistical aggregates.
6. **Joint sticky noise.** The two pinned computation peers derive independent
   private seeds from their dedicated roots and the same committed query. An
   exact GC combines the seeds, samples one globally calibrated mechanism, adds
   it to the aggregate and opens only the noised result. Neither seed nor the
   realised noise is an output of the circuit. This gives the target mechanism,
   without the `sqrt(K)` utility penalty of summing complete independent draws,
   provided at least one of those two computation peers is honest and does not
   collude. Until that circuit is attested, the only safe fallback is for every
   designated noise peer to contribute a complete global draw; its additional
   convolution noise must be reported explicitly in the accuracy certificate.
7. **Commit before release.** The result digest, mechanism, allocation and
   certificate are durably committed. Crash recovery either reproduces the
   same payload or fails closed; it cannot draw again.
8. **Replay.** The same canonical query, logical snapshot and privacy epoch
   returns the byte-identical release without further privacy cost. A distinct
   query receives an independently derived stream and composes normally.

Every transition is idempotent by operation ID and transcript digest. Missing
DSI replies are retryable only for such transitions. A present malformed reply,
authentication failure, state conflict or typed terminal failure is fatal.

### Canonical-query rules

- Client aliases map to one canonical method and all defaults are materialised.
- Formulae use a validated AST, fixed contrasts/reference levels, explicit
  intercept/offset/weights and a canonical term/interaction order. Raw R code,
  environments and arbitrary predicates are not hash inputs or executable
  query descriptors.
- Dataset names map to custodian dataset/snapshot identifiers. Cohort selection
  uses approved named cohorts or a separately accounted fixed-domain query,
  never an analyst-provided filter expression.
- Numbers use one finite canonical representation; unordered named objects and
  fixed domains are sorted, while semantically ordered objects retain order.
- Output subsets, confidence levels and presentation choices are
  post-processing of the same maximal authorised release where possible. They
  cannot create a new noise draw for the same underlying statistic.
- The server recomputes and signs the canonical ID. It never trusts a query ID,
  allocation, epsilon, sensitivity or mechanism supplied by the analyst.

## Noise-root, identity and ledger contract

- One independent 256-bit DP root is generated automatically per peer and is
  cryptographically domain-separated for every privacy domain. It is never
  derived from the independently generated Ed25519 identity seed and neither
  secret is packaged in, or generated while building, the R library.
- `configure` creates only private state-directory scaffolding. The first real
  cryptographic service operation creates any missing roots outside the
  installed library, before it can return a response, under an inter-process
  lock, using a cryptographic operating-system random source, a same-directory
  temporary file, atomic rename and owner-only permissions. Package
  installation, `.onLoad()`, staged loading and `R CMD check` never mint
  deployment secrets, which prevents a container-image build from cloning one
  identity into every deployment.
- A file provider must use an absolute path outside the package, a regular
  non-symlink single-link file, correct owner, mode `0600`, parent mode `0700`,
  stable pre/post-read metadata checks, and location-policy revalidation after
  ancestor symlink resolution. A configured HSM/KMS provider may
  replace the file provider, but no administrator action is required for the
  secure default.
- Service restart reopens the same root. Once identity and noise root have
  both been validated, reciprocal encrypted-and-MACed recovery envelopes are
  maintained automatically: the identity is wrapped using domain-separated
  keys obtained from the noise-root HMAC interface (including HSM/KMS), and a
  file-backed noise root is wrapped using the independent identity. Loss of
  either primary is restored to the exact original value without changing
  pins, epoch or sticky noise. Malformed, permission-weakened or conflicting
  state is never silently replaced.
- Private, durably synced deployment receipts remain authoritative continuity
  markers. Existing pre-receipt installations are upgraded only after the
  candidate secret authenticates every surviving local or joint v1/v2 ledger;
  orphan non-empty WAL/SHM sidecars are incomplete history. If both independent
  primary roots disappear simultaneously, the reciprocal envelopes cannot
  authenticate each other. The server preserves the old evidence and generates
  a fresh identity, but every other server rejects it as
  `peer_not_recognized` until an administrator verifies its fingerprint out of
  band and updates the persistent name-bound pin. There is no analyst/relay
  autoaccept path.
- Derivation binds domain, cohort, method, canonical query, logical snapshot,
  mechanism version, allocation, coordinate/iteration and privacy epoch.
- A cryptographic deterministic generator and audited discrete sampler are
  used. R's `set.seed()` and non-injectable process-global randomness are not
  used for sticky releases.
- Root rotation before a capsule's first valid START changes its candidate
  release instance and begins a new privacy epoch; any new capsule still
  requires an available lifetime unit. After START claims an instance, rotation
  cannot rebind that capsule: it must continue/restore the exact claim or fail
  closed. A changed key ID with an unchanged epoch is a fatal configuration
  error. A published capsule is never rerolled under a replacement root.
- Encrypted, tested backup of the root and the cross-signed
  query-to-allocation/result receipts is operationally required. A root
  without durable receipts does not provide reinstall-safe replay. If every
  designated peer and every optional external anchor are rolled back together,
  no local process can distinguish that state from a genuinely new deployment.

### Bounds and data-dependent failure

Contribution, coefficient, time, weight and numeric bounds are provisioned by
the custodian or obtained through a separately accounted DP mechanism. They are
never inferred exactly from the protected snapshot during a query. DP adapters
clip contributions before ring encoding, so an outlying protected value cannot
select a visible success/failure branch. Finite/schema/alignment validation is
performed when the immutable snapshot is admitted; the query plane sees only a
stable attestation or one uniform unavailable response. A clipped-count or
data-quality diagnostic is a separate DP release, not an error-message field.

Snapshot admission extracts columns without invoking data-frame subclass
methods and freezes the accepted object as a base `data.frame`. Columns may use
base atomic storage, factor/ordered, Date, POSIXct or difftime. Unknown S3
classes and extra attributes are rejected with an explicit conversion error;
they are never silently interpreted through user-defined `[`/`as.numeric`/
`as.character` methods. The authenticated PSI manifest is the only retained
frame-level attribute and is revalidated against the frozen identifier order.
This restriction is a semantic-safety contract, not a data-size limitation.
Patient and PSI identifiers use one option-independent UTF-8 canonical form.
Character and factor IDs retain their stored labels; numeric IDs are accepted
only when finite, integer-valued and exactly representable within 53 bits, then
formatted without `digits`, `scipen`, locale or `OutDec` dependence. Fractional
or wider numeric identifiers are rejected and must be supplied as character
strings, preventing a display option from changing privacy-unit grouping.

Every informative publication is one complete reusable analysis capsule for an
immutable logical snapshot and versioned workload. Its creation uses one
consortium publication sequence replicated and anchored by the pinned
computation peers. The sequence orders capsules and prevents replay, fork and
rollback; it is not an operation counter and does not determine epsilon.
Methods and arguments never enter the durable capsule identity. Local and joint
measurements are materialised together under the capsule's fixed epsilon/delta
contract, avoiding both per-peer `epsilon/K` partitioning and per-query decay.
The retired local-ledger split remains only in internal compatibility tests.
Until the cross-signed DSI adapter passes E2E promotion, production fails before
creating a capsule, with or without an optional external CAS.

A `joint_mpc_single_opening` is a one-time coordinate finalizer during capsule
creation, not an allocation performed when an analyst calls a method.
The preferred two-computation-peer GC samples exactly one mechanism calibrated
to the full global `(epsilon, delta)`, not `(epsilon/K, delta/K)`. Its signed
commitments bind both sticky seeds, the capsule, snapshot, privacy epoch and
circuit digest before release. If the conservative independent-draw fallback
is used instead, every draw is calibrated to the same full global release and
the certificate covers their convolution. Opening any individual noised share
as a statistic, or opening more than one independently useful combination,
changes this accounting and is forbidden.

## Requests, composition and availability

There is no transport rate limit, arbitrary request-count quota or accuracy
decay caused by request history. Bounded spools, worker concurrency and
backpressure protect availability and do not alter privacy parameters. Ten or
ten thousand methods over one capsule are post-processing of the same bytes and
add zero privacy loss. The typed lifetime/publication token is an opaque
terminal union: either a previously unseen capsule reaches exhausted `N`, or
the requested capsule's irrevocable instance-claim/publication binding prevents
the exact instance from safely continuing or replaying. It does not reveal
which cause occurred or imply that the global remaining count is zero.

Distinct capsules for overlapping cohorts or successive snapshots compose.
Each uses the same fixed per-capsule parameters; authenticated exact composition
is bounded by `N * epsilon <= 8` and `N * delta < 1`, and capsule `N+1` is
denied. This is a finite bounded-lifetime contract, while exact replay remains
unlimited. A new epoch does not erase reservation or composition history.

## Method-family migration matrix

| Family | Target release | Required contribution contract | Promotion blocker |
|---|---|---|---|
| Identity, security/numeric status, method registry | Public fixed schema | No protected statistic | Unanimous schema/pin/policy validation |
| Alignment | Stable attestation only | Pre-aligned immutable cohort; no patient-derived digest/count | Fixed-shape padded PSI or custodian pre-alignment |
| Count | Sticky discrete scalar mechanism | One patient contribution | Dedicated root, recoverable receipt mapping and sampler audit |
| Contingency, epidemiology, standardisation | Sticky fixed-domain histogram; all effects are post-processing | One bounded cell contribution per patient; public domains | DP-aware inferential calibration; no ordinary Fisher/chi-square claim on noisy cells |
| Bounded mean/variance/descriptives | Sticky noisy count, sum and sum of squares | Public numeric bounds and one contribution per patient | Integer-grid overflow proof and sampling-uncertainty labelling |
| Correlation and PCA | Noisy clipped second-moment matrix followed by PSD projection/eigendecomposition | Public vector norm, dimension and centring contract | Joint cross-site covariance adapter; no individual scores; eigengap/utility certificate |
| Gaussian regression | Private stable/sufficient-statistic regression with explicit coefficient/design bounds | Clipped row vector and outcome | Positive-definite private solve, estimator label and DP-valid uncertainty |
| Binomial/logistic GLM | Objective/gradient perturbation with fixed public optimisation schedule | Clipped row gradient, bounded coefficient domain | Joint sticky vector noise; no exact score/Hessian oracle |
| Poisson GLM | Fixed-iteration clipped-gradient mechanism | Bounded outcome, feature norm, coefficient/eta domain | Tail/clipping estimand and approximation certificate |
| Negative binomial | Fixed-iteration clipped joint score for coefficients and dispersion | Poisson bounds plus public dispersion domain | DP-stable theta update and valid uncertainty |
| Multinomial and ordinal | Fixed-iteration clipped vector-gradient mechanism | Public class domain, feature norm and coefficient/threshold domain | Joint-estimator covariance and ordered-threshold postconditions |
| LASSO | DP proximal/coordinate optimisation; lambda path is public or privately selected | GLM contribution bound plus public lambda grid | KKT/reference validation; legacy post-hoc sparsifiers are not promoted |
| Cox PH | Clipped private partial-likelihood optimisation; baseline hazard is a separately budgeted release | Patient contribution, bounded features/time grid, public ties/strata policy | Exact risk-set adapter, robust covariance and DP-valid inference |
| Kaplan-Meier, Nelson-Aalen, competing risks | Fixed public time-grid event/at-risk histograms and post-processing | One bounded patient trajectory | New high-level methods and simultaneous-band validation |
| LMM, GEE and GLMM | Cluster-level clipped gradients/moments entirely inside MPC | Patient/cluster adjacency, max records and cluster norm | Remove client-visible cluster scores/meat; complete estimand and inference validation |
| IPW/causal | Composed private propensity and outcome stages with clipped weights | Public positivity/weight bounds and estimand | Current known-weight wrapper is not end-to-end IPW; add AIPW/g-computation only with separate validation |
| Multiple imputation | Private imputation model and releases, Rubin pooling as post-processing | Declared local/cross-site MAR scope and bounded variables | Cross-site predictors, keyed imputation randomness and non-zero between-imputation variance |
| LR/Wald/contrast/confidence helpers | Post-processing only of a DP-valid fit and DP-valid covariance/test release | Inherit parent fit contract | Naive non-private reference distributions must not be reused silently |

Where a family has no tight practical mechanism, a conservative bounded-output
mechanism may be used as an explicitly labelled fallback, but it is not called
“promoted” until utility tests show it is scientifically usable. Regularisation,
Firth correction, pseudoinverse, clipping and discretisation may be offered as
explicit estimators; none may silently replace the requested estimand.

## Numerical backend contract

- Inputs are finite and checked against custodian-owned bounds before encoding.
- The planner selects the smallest attested Ring63, Ring127 or dynamic `2^k`
  backend that satisfies raw-product capacity and the public error budget.
- Secret multiplication produces an untruncated Beaver share only after a
  server-owned no-wrap proof. Exact signed floor truncation then returns fresh
  random shares. Local share truncation is not a promoted fallback.
- Signed comparisons use the exact purpose-bound GC/OT adapter. Historical DCF
  comparison is removed from promoted paths.
- If no fixed ring suffices, the dynamic multi-limb backend is selected. It
  must support both the required magnitude and fractional precision; CRT is an
  optimisation, not a replacement for secure comparison/truncation.
- Approximate `exp`, sigmoid, reciprocal, softplus and log operations publish a
  validated domain and worst-case approximation error. Domain violation is a
  typed failure, not extrapolation.
- Every fit returns a numeric certificate containing backend, ring/precision,
  public bounds, operation counts, approximation error, aggregate error,
  estimator error and postcondition/identifiability status. Preflight alone is
  never labelled execution certification.

## DSI and large-data contract

- Use one named asynchronous DSI fan-out per protocol phase, with a distinct
  expression/payload for each node. Do not serialise a generic command batch.
- Use direct scalar arguments and one necessary Base64 layer. Outer JSON/Base64
  envelopes are removed after the purpose-bound transport is stable.
- Short/missing replies never advance offsets. Absolute offsets, exact ACKs,
  immutable chunk geometry and transcript digests make retries idempotent.
- Input sharing is streamed by row block from a server-minted operation;
  matrices, shares and ciphertexts are not assembled as multi-gigabyte R/JSON
  strings. The relay can pump only the fixed destination and operation.
- Memory is `O(block_rows * p + p^2)` for dense GLM sufficient statistics.
  High-dimensional fits use Hessian-vector/iterative algorithms rather than a
  dense `O(p^2)` Hessian when explicitly selected.
- Public circuits/manifests and base OT epochs may be cached. Triples, GC
  labels, OT-extension material, nonces, random shares and DP noise are never
  reused across operation IDs.
- Per-operation/session/node spool reservations, TTL cleanup and backpressure
  are resource controls, not request quotas. Terminal receipts allow cleanup
  while preserving idempotent retry.

## Mandatory acceptance gates

1. **Surface:** every client DS method literal is registered by the companion
   server; every registered DS method has a client consumer; retired primitives
   have neither registration nor caller.
2. **Single mode:** production rejects `legacy_exact_mpc` and `strict_dp`; no
   environment variable or option can enable an exact release plane.
3. **Coverage:** every exported high-level statistical method maps to a complete
   release contract. A method without one is unavailable in the release build,
   not silently routed through legacy code.
4. **Capabilities:** forged, cross-method, cross-snapshot, wrong-peer, expired,
   replay-conflicting and out-of-order capabilities fail before data access.
5. **Sticky DP:** semantic aliases, restart, reconnect and intact/restored local
   state reproduce bytes; distinct queries/epochs use independent streams.
   Missing roots are recovered where possible. A release domain may select
   another candidate only before the first valid START claim; afterwards state
   loss must restore and continue the exact instance or fail closed. After
   publication only exact replay/restore is allowed, never an old mapping under
   a replacement root or a second instance. A
   peer-identity replacement remains
   untrusted until administrators update its name-bound pin. Rollback and lost
   mappings are detected within the retained-ledger threat model.
6. **Composition:** every operation over one committed capsule is free
   post-processing and can never reduce accuracy or availability. Creating a
   genuinely different capsule burns one lifetime unit at allocator commit and
   adds its fixed declared `(epsilon, delta)` to authenticated bounded
   composition. That history may deny capsule `N+1`, but it is not a request
   quota and never affects exact replay.
   All coordinates, iterations, peers, model selection and inference products
   must already be included in the immutable capsule workload.
7. **Transcript:** successful authenticated semantic-operation counts and
   non-output payload geometry are derived only from the signed public
   contract; output-dependent serialization is post-processing of the DP
   vector. Protected failure text is uniform. Host timing, failure occurrence,
   DSI polling, retransmission, scheduling and availability are explicitly
   outside the package claim; a traffic-flow claim requires a separately
   validated scheduled/asynchronous gateway.
8. **Ring:** exhaustive small-ring and property tests cover negative values,
   boundaries, f=0, maximal fractions, no-wrap guards, dynamic widths and fresh
   output shares. No wrap or truncation error can be returned as a valid result.
9. **Statistics:** central-oracle simulations cover baseline, boundary,
   weighted, missingness, separation, non-identifiability and longitudinal
   cases. DP tests report bias, RMSE, coverage/type-I error and power over a
   declared `(n, p, epsilon, delta)` grid.
10. **Transport:** multiprocess large transfers cover partial writes, partial
    NULL replies, response loss, replay, reconnect, bounded spool/backpressure,
    tamper, terminal cleanup and byte-identical hashes.
11. **Scale:** benchmarks cover row streaming, large `n`, high `p`, many sites
    and constrained spool/memory. Performance claims include DSI phases, bytes,
    peak RSS, disk, CPU and wall time.
12. **Build:** all platform binaries are rebuilt from the final source with a
    checked runtime manifest and hashes; server/client tests, `R CMD check`, Go
    tests and real multi-process DSI E2E all pass from a clean installation.

## Transition rule (resolved)

The installed runtime now has one immutable `disclosure_safe` profile and no
selector that can enable the historical exact/adaptive release plane.  This did
not make incomplete estimators safe: retained public names without a reviewed
capsule route fail before DSI and remain explicitly `quarantine` in
`ds.vertMethodStatus()`.  A method is promoted independently only after its
own disclosure, numerical, statistical and end-to-end gates pass; compatibility
of a name never grants remote capability.

## Primary technical references

- [Differentially Private Empirical Risk Minimization](https://jmlr.org/papers/volume12/chaudhuri11a/chaudhuri11a.pdf)
- [Our Data, Ourselves: Privacy via Distributed Noise Generation](https://www.iacr.org/archive/eurocrypt2006/40040493/40040493.pdf)
- [The Distributed Discrete Gaussian Mechanism](https://proceedings.mlr.press/v139/kairouz21a.html)
- [Differentially Private Chi-Squared Hypothesis Testing](https://proceedings.mlr.press/v48/rogers16.html)
- [A Near-Optimal Algorithm for Differentially-Private PCA](https://jmlr.csail.mit.edu/beta/papers/v14/chaudhuri13a.html)
- [Differentially Private Survival Function Estimation](https://proceedings.mlr.press/v126/gondara20a.html)
- [SoK: Truncation Untangled](https://petsymposium.org/popets/2025/popets-2025-0135.php)
