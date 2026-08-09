# dsVertClient (development version)

### Estimators

* `ds.vertLMM` now dispatches on K internally (K=2 exact closed-form, K>=3
  variance-ratio profile), matching every other method; it previously used only
  the first peer at K>=3, silently producing an incomplete fit. `ds.vertLMM.k3`
  is retained as a deprecated passthrough, and the `ds.vert.lmm` alias forwards
  to `ds.vertLMM` for every K.
* Gaussian K=2 GLM is solved in one exact Gram step instead of the iterative
  L-BFGS loop (~11x faster, exact); `options(dsvert.gaussian_oneshot = FALSE)`
  restores the loop.

### Security

* Expected cold misses for the joint-DP allocation, result and release replay
  phases now cross DSI under one fixed `phase_not_ready` tag. The client
  reconstructs a non-retryable sanitized condition and advances a prerequisite
  phase only when every relevant failure carries that tag; mixed, untyped,
  poisoned, capacity and deadline failures remain terminal. A fetched tagged
  async result is terminal and consumed, so Opal or Armadillo sessions are not
  poisoned merely because a replay phase has not yet been committed.

* `ds.vertDPStatus()` now accepts only joint-DP capsule status v5 and prints the
  fixed per-capsule parameters, lifetime maximum, burned reservations,
  publications and remaining distinct-capsule units. Exact validation requires
  `N * epsilon <= 8` and `N * delta < 1`, so a vacuous delta state is rejected
  rather than warned about. `ds.vertSecurityStatus()` does not reinterpret
  exhaustion as consortium unreadiness. New capsule reservation is an
  operation/history gate, while `request_limit` is false and exact replay of
  the same capsule/release instance remains unlimited post-processing.

* Public DP results now carry `dsvert-capsule-security-claim-v3`. It binds the
  exact field
  `authenticated_history_retention_assumption=at_least_one_noncolluding_designated_noise_peer_retains_and_uses_complete_authenticated_monotonic_history`
  and
  `privacy_accountant_namespace_assumption=one_stable_unique_namespace_across_domain_cohort_policy_pinset_and_ledger_reconfiguration_per_protected_privacy_universe`.
  The namespace condition is a custodial deployment assumption, not client or
  server enforcement and not automatic accounting migration. The claim states
  explicitly that simultaneous rollback of both designated histories is not
  protected without an external linearizable CAS. The client recognizes the
  fixed, non-retryable,
  detail-free `[dsvert_dp_lifetime_budget_exhausted:v1]` condition and never
  treats it as a reason to issue another release attempt. The token is an
  opaque terminal union of global `N` exhaustion and a requested capsule whose
  irrevocable instance claim/publication binding prevents the exact instance
  from safely continuing or replaying; it does not imply
  `remaining_distinct_capsules == 0`, and separate tokens would reveal state.

* Froze the biomedical vector ABI as a pre-release baseline: signed PREPARE
  receipts are v6; signed START, RESULT, RELEASE and ACK receipts are v5; the
  authenticated STORE schema is v6; the claim row and authenticated claim-set
  state are `dsvert-joint-dp-vector-instance-claim-v1` and
  `dsvert-joint-dp-vector-instance-claim-state-v1`; replay is v4. Every signed
  receipt and replay response must attest exactly `history_gate=TRUE`,
  `request_limit=FALSE`, and `operation_limit=TRUE`. Neither package migrates or
  re-signs legacy v4/v5 stores automatically: rollout requires empty DP capsule
  state or a future audited offline migration, and legacy state fails closed.
  Earlier Opal deployments had `POLICY_READY=FALSE`, published no DP capsules,
  and used ephemeral local K-site state.
  Exact COMMIT/RELEASE and sticky replay in a live session are O(1). A cold
  end-to-end reconstruction after restart instead returns through
  AllocationProof and audits the complete O(N) allocator journal before the
  proof; it still consumes no lifetime unit, noise draw or protected data.
  The baseline also fixes `dsvert-joint-dp-control-v3` and
  `dsvert-joint-dp-capsule-identity-v3`, both binding the lifetime fields and
  biomedical workload v7. Their v2 artifacts are legacy, receive no automatic
  migration or re-signing, and fail closed; rollout requires empty DP capsule
  state or a future audited offline migration.

* The client now accepts only the v7 biomedical workload lifecycle contract:
  the signed manifest must bind every included coordinate to the registered
  secret-sharing, encrypted transport, joint sampler, confidential finalizer
  and durable replay path. It rejects missing, downgraded or internally
  contradictory lifecycle attestations before source transport. This does not
  assert that a live Opal or Armadillo deployment has executed successfully;
  connector availability remains a runtime preflight. Raw source tickets and
  ciphertext envelopes remain explicitly non-releasable.

* Repeated capsule consumers now keep a four-entry, 64 MiB in-process LRU of
  validated public DP vectors. Every call still refreshes pinned status and
  rebuilds the authoritative manifest; the cache key also binds its
  `manifest_sha256` and `capsule_id`, so a changed snapshot or workload cannot
  serve stale numbers. Hits omit only source sharing, sampling and vector
  replay. Connection handles and failures are never retained, and eviction
  cannot gate an operation or reroll sticky noise.

* Large capsule sources now relay consecutive authenticated encrypted chunks in
  a byte-bounded window over the existing DSI methods. The outer
  chunk and acknowledgement framing is not independently signed; every inner
  ticket/summary wrapper, ciphertext envelope and acknowledgement retains its
  existing signature and complete binding. A legacy scalar ticket is
  obtained first; capability negotiation is public and synchronous, and the
  client uses a window only after both designated recipients and every source
  owner attest the adaptive v2 contract. Within one synchronized DSV1 package
  generation, the adaptive v1 negotiation remains byte-exact and enforces its
  original 768 KiB/1 MiB caps; peers without adaptive v2 use the scalar route
  after a terminal application rejection, without poisoning an Armadillo
  session. DSV1 is a deliberate lockstep wire change, so client and server
  packages must be upgraded together; a mixed package generation fails at the
  first framed, data-free manifest phase before protected data access. A
  safe staged rollout is server-first with release traffic paused: deploy,
  reconcile and attest every server, verify the data-free framing boundary,
  and only then deploy the paired client. The intermediate mixed state is
  deliberate fail-closed downtime, not backward compatibility. A
  transport failure whose post-`execute`
  state cannot be proven terminal marks the affected authenticated handles
  unusable and requires fresh login connections rather than issuing a scalar
  request against an ambiguous singleton result lifecycle. This
  reduces DSI round trips while leaving ciphertexts, DP distribution and the
  no-request-quota contract unchanged.

* Capsule-source windows now derive their effective width from three separate
  public limits: the source node's data-free response probe, each designated
  recipient's independently observed request/expression geometry, and fixed
  8 MiB application caps with explicit response and framing headroom. The v1
  request probe remains byte-identical; an optional response-padding extension
  contains no protected data or protocol state. HTTP 413 descends only the
  public candidate ladder, reconnects re-probe the exact session handle, and
  an ambiguous singleton-result lifecycle poisons that handle. A server that
  lacks the extension within that synchronized DSV1 generation, or advertises
  the signed adaptive v1 768 KiB capability, uses the byte-identical scalar
  route. At one million coordinates, the
  committed data-free benchmark records `W=8` reducing source data phases from
  492 to 64 for two sources, 738 to 96 for three, and 1,230 to 160 for five;
  inner signed artifacts and their hashes are unchanged.

* Server-side capsule-source admission now accounts for the authenticated
  receipt retained for every source-by-chunk pair. The public capacity bound
  therefore scales with consortium size and applies byte backpressure before
  mutation. Recipient key/ticket rows receive a separate atomic reservation,
  with idempotent replay, authenticated migration and compaction release. Both
  bounds remain independent of the number of questions, sessions or releases.

* Remote-surface readiness now uses one connector-neutral custodian token over
  the exact alias-free dsVert method map. Opal stores it in the reconciled
  server profile; Armadillo/Rock may receive the same token only through the
  server/container environment after a native administrative inventory check.
  The client sends only the zero-argument profile query and cannot inject an
  option, environment value or attestation. Missing, malformed, conflicting or
  stale server assertions fail closed; no installed package auto-attests.

* Opal provisioning now uses a disabled, restricted, non-default `dsvert`
  profile; removes the complete prior aggregate/assign inventory; registers
  and re-verifies exactly the staged server tarball's alias-free allowlist; and
  stores its deterministic surface attestation before enablement. TLS
  verification and administrator/runner separation are mandatory. Provisioning
  resolves and validates the profile's actual Rock cluster, installs the server
  tarball's declared imports through the custodian-configured repository, and
  verifies package removal/installation against Rock inventory rather than an
  empty HTTP 200 response; a valid zero-row ACL inventory is handled explicitly.
  The runner
  receives only profile use plus Opal's DataSHIELD-capable dictionary/summary
  table permission, never individual-value, edit or administrative access.
  `ds.vertSecurityStatus()` cannot report ready when any server's custodian
  attestation is missing or stale, even if the joint-DP handshake succeeds.
  Its typed route map now separates biomedical joint-DP surface eligibility
  from runtime policy/consortium readiness and reports formal GLM and Cox as
  not ready, preventing the legacy eligibility flag from being read as a
  blanket model-readiness claim.

* Release evidence now distinguishes topology coverage from connector
  execution. `K=2`, `K=3` and `K=5` for the promoted capsule/PSI route denote
  unit, adversarial and isolated-process/DSLite coverage unless a current live
  harness artifact is cited explicitly. Analytical calls are constructed
  through DSI; Opal method reconciliation and Armadillo TLS/session inspection
  are deployment and lifecycle safeguards, not analytical protocol
  dependencies. Local S4/httpuv Armadillo smoke and cached Opal results are not
  current-artifact live E2E claims; artifact-exact Opal and Armadillo/Rock
  smoke remains a release gate.

* Removed the `allow_unattested_numeric` compatibility route from
  `ds.vertGLM()` and `ds.vertNumericPreflight()`. GLM now authenticates and
  validates every custodian numeric policy before schema discovery, row-count
  access, alignment inspection, or MPC setup; options and environment
  variables cannot restore an unattested certificate.

* `ds.vertSecurityStatus()` now accepts only the server's versioned single
  `disclosure_safe` profile and its joint-capsule readiness handshake. The
  former mode selector and `require_strict_ready` terminology were removed;
  `require_ready` controls only whether an unavailable consortium dependency
  is returned as status or raised as an error, never whether exact methods are
  enabled.

* Made `ds.vertDPCalibrate()` an honest client-only utility preview: both
  Laplace and fixed-work dyadic discrete-Gaussian backends are recognised as
  deployed, while the exact signed server plan remains the sole mechanism
  selection authority. Its retired `decay`, `release_indices`,
  `composition_partitions`, `total_epsilon` and `total_delta` arguments were
  removed so the public API cannot imply a request quota or let the analyst
  select the server-owned lifetime bound.

* The two remaining internal consumers of legacy correlation statistics are
  closed. Exact slope-binomial `ds.vertLASSOIter()` now requires an explicit
  signed id and fails with a typed `binomial_lasso_design_grams` availability
  condition before DSI; multinomial slope models do the same for
  `multinomial_design_grams`. Neither route substitutes `ds.vertCor`, pairwise
  statistics, or same-owner local moments whose clipping, mask, snapshot,
  scaling and design order are not bound to the score MPC.

* Joint-DP status and release-instance validation now include each designated
  peer's independent public release domain. Normal retries remain
  byte-identical; surviving authenticated release receipts and chunks restore
  missing domain metadata and replay the original vector without another
  sampler execution. PREPARE candidates may coexist, but the first valid START
  at each designated peer atomically and irrevocably claims one instance before
  local staged-source access or sampling. A release domain may select another
  candidate only before that claim. Afterwards only the claimed instance may
  progress, restore or replay. Matching bilateral receipts prevent split
  sibling claims from becoming a release. Once a capsule has one public release
  instance, irrecoverable loss must restore and replay that exact instance or
  fail closed; a fresh domain can never resample or publish a second instance
  for the capsule.

* `DSMolgenisArmadillo` 3.0.2 HTTPS connections are now validated directly
  from their real S4 handle and need no dsVert attestation option or extra
  client input. Remote HTTP, global curl verification downgrades, and present
  but uninspectable `httr_config` values fail closed. The real connector's
  aggregate/poll/fetch/keepalive S4 lifecycle is covered by regression tests.
* An asynchronous DSI timeout no longer implies nonexistent job cancellation.
  Because DSI has no portable cancellation primitive, a connection with an
  unresolved job is marked unusable and the retry must use a fresh login
  connection; this prevents Armadillo's singleton last-result lifecycle from
  being associated with a later request. No implicit `dsDisconnect()` occurs.
* Immutable DSI frame/exchange calls now launch and poll through DSI's public
  job primitives, avoiding repeated deparsing of large expressions while
  preserving named fan-out, keepalive and fail-closed callbacks. A warmed,
  order-balanced four-pair DSLite replay/reconnect benchmark recorded 28.0%
  lower median wall time with byte-identical submissions and terminal hashes;
  this is local evidence rather than an Opal/Armadillo speed guarantee.
* Request-size negotiation reuses only a connector-profile probe-order hint and
  revalidates it after reconnect. The exactly measured DSLite 1.4.1 path starts
  at its portable ceiling, reducing cold public padding from five candidates
  (15.625 MiB) to one (0.625 MiB). Unknown versions continue through an
  extended ladder down to 16 KiB, and failure of every candidate now stops
  before any opaque payload instead of applying an unverified fallback. These
  public, stateless probes run synchronously so an oversized candidate rejected
  by a remote parser cannot leave an ambiguous async job that blocks descent;
  statistical and MPC phases retain named asynchronous fan-out.
* Per-site assignments now use DSI's expression-assignment API and require one
  success callback for every named site. Direct asynchronous aggregate polling
  has a monotonic availability deadline instead of an unbounded connector
  poll; an unresolved job poisons that authenticated session because DSI does
  not expose portable cancellation.
* Reciprocal exact-GC validity receipts now share one named DSI fan-out per
  multiplication chunk, removing one sequential round trip while retaining
  target-specific peer blobs and fail-closed partial-result handling.
* Exact-GC serializes each authenticated outbound envelope once and reuses the
  frozen JSON for every ambiguous DSI retry until acknowledgement, avoiding a
  payload-sized allocation and parse preparation on every polling cycle.
* The exact-GC client now accepts the server's certified dynamic Ring63--4096
  contract, including 1024/2048/4096-bit wire containers and the matching
  64/16/4/1 direct-multiplication chunk limits. Oversized circuit shapes still
  fail before transport; the former hard Ring512 client ceiling no longer
  disables the multiprecision fallback.
* Producer-spooled streaming no longer writes a second payload-sized file at
  the analyst. The recipient's complete hash and pinned signed receipt remain
  mandatory and producer-verified, reducing client disk from O(payload) to
  O(frame) and avoiding two payload lengths of redundant integrity I/O.
* Producer-spooled transfers now pipeline recipient commit of frame `i` with
  producer read of frame `i+1` in one named two-site DSI fan-out. A partial
  reply replays the identical pair and advances neither cursor, covering both
  lost-read and lost-store ACKs while reducing a healthy `F`-frame transfer
  from `2F + 1` transport phases to `F + 2`.

* Added the internal client relay for the provisional global-DP control plane.
  It verifies the complete K-peer policy, derives the exact two
  custodian-designated noise peers, verifies their Ed25519 receipts and safely
  retries only byte-identical durable control phases after a missing ACK. No
  statistical method is promoted and payload delivery remains unavailable.
* The relay now accepts only v2 joint-DP result receipts carrying an opaque
  server-keyed HMAC commitment. Enumerable v1 public payload hashes fail
  closed, and lost result acknowledgements replay the exact durable receipt.
* NB full-regression and LMM K>=3 flows register the identity-verified peer set
  so their recipient-pinned share seals hold across all K.

### Differentially private statistical validity

* Added `ds.vertDPPrevalenceRatio()` and
  `ds.vertDPPrevalenceRatioInference()` as strict zero-call cross-sectional
  naming views of the existing DP 2-by-2 epidemiology results. Exposed and
  prevalent orientations are mandatory; points, regions, coverage, DP
  provenance, number-needed boundaries and zero-cost metadata are inherited
  without recalculation or narrowing. The result records that study design is
  caller-declared and cannot be inferred from the released table.
* Added `ds.vertDPFrequency()` for a signed fixed-domain categorical marginal
  from the reusable sticky capsule and `ds.vertDPFrequencyInference()` for
  zero-call conservative population-proportion regions. The latter combines
  the simultaneous DP count-box event with Bonferroni exact
  Clopper--Pearson intervals, validates capsule/owner/domain provenance, and
  adds neither a server call nor another privacy release.
* Added `ds.vertDPEpi2x2Inference()`,
  `ds.vertDPDiagnostic2x2Inference()` and
  `ds.vertDPDirectStandardizationInference()` as zero-call post-processing
  that combines the signed simultaneous DP count-box event with conservative
  exact-binomial sampling regions. Empty integer mechanism boxes return the
  full parameter space instead of a misleading narrow interval.
* Added `ds.vertDPCausalStandardization()` and its inference companion for a
  saturated public-stratum g-formula over one already-authorised
  stratum-by-treatment-by-outcome DP table and fixed public target weights.
  The output states consistency, conditional exchangeability, positivity, no
  interference and mapping/target-population assumptions explicitly; it does
  not claim propensity estimation or double robustness. Exhaustive integer-box
  tests and the reproducible 1,000-replicate statistical harness cover its
  identities, mechanism regions and boundary states.
* Fixed-work discrete-Gaussian releases now support arbitrary requested
  confidence levels from the signed rational variance, vector total-variation
  allowance and finite support. The published signed 95% radius remains the
  exact authority at 0.95; higher-confidence requests can never return a
  smaller radius, and an exhausted tail-transfer allowance falls back to the
  unconditional signed support bound.
* Added `ds.vertMantelHaenszel()` for a common odds ratio and explicitly
  conditional classical MH inference from already-authorised stratified 2x2
  aggregates. Added `ds.vertDPMantelHaenszel()` as zero-call, zero-additional-
  privacy-cost post-processing of one fixed strata-by-four-cells capsule. Its
  simultaneous mechanism region is exhaustively box-validated, types finite,
  zero, infinite and non-estimable boundaries, and never reports a classical
  CMH p-value for DP-noised counts.

* `ds.vertDPContingency()` and `ds.vertChisqCross()` now use the signed
  cross-owner categorical capsule artifact. The client relays only pinned
  allocation proofs and opaque fixed-shape transport, then performs the same
  DP-aware chi-square or conditional 2-by-2 plug-in bootstrap from one sticky
  release. The legacy exact-table, one-hot discovery, and classical Fisher
  routes are no longer reachable from these APIs.

* Added `ds.vertDPQuantile()` and `ds.vertDPMedian()` as zero-cost,
  release-only post-processing of one validated `ds.vertDPDescribe` capsule.
  They return a fixed public bin and interval plus a simultaneous mechanism
  region; they do not claim an exact sample quantile or interpolate within a
  bin, and issue no DSI calls.
* Added `ds.vertDPIndirectStandardization()` for zero-cost SMR/SIR-style
  observed-to-expected post-processing of one authorised DP strata table and
  public expected rates. It returns simultaneous mechanism-noise regions and
  never substitutes a classical Garwood interval on noisy counts.
* Added `ds.vertDPIndirectStandardizationInference()` as a separate zero-call
  Poisson/Garwood inference companion. It envelopes the exact interval over
  the signed integer DP count box, combines sampling and mechanism failure by
  union bound, preserves source provenance, and returns `[0, Inf]` when the
  integer box is empty or the expected denominator may be zero. It draws no
  randomness, creates no release, and makes no p-value or causal claim.
* Added `ds.vertDPRMST()` as zero-cost client-side post-processing of one DP
  survival release. It integrates the left-continuous product-limit curve on
  the public grid and propagates the simultaneous mechanism band without
  presenting it as sampling inference or continuous-time grid accuracy.
* Added `ds.vertDPRMTL()` as the exact zero-call complement of that RMST over
  the public release interval. It uses
  `(tau - time_lower_bound) - RMST`, reverses the existing mechanism limits
  without recalibration, preserves release provenance, and makes no sampling
  or continuous-time grid claim.
* Added `ds.vertDPSurvivalQuantile()` and `ds.vertDPMedianSurvival()` as
  zero-call inversion of the validated fixed-grid Kaplan--Meier curve and its
  simultaneous DP-mechanism band. Targets beyond the signed public horizon
  are represented by explicit estimability states rather than extrapolation;
  sampling uncertainty and continuous-time grid error remain excluded.
* Added pure post-processing `ds.vertDPDiagnostic2x2()` for sensitivity,
  specificity, PPV/VPP, NPV/VPN, prevalence, accuracy, balanced accuracy,
  F1, LR+, LR−, and diagnostic odds ratio. Disease-row and test-column
  orientation is explicit;
  simultaneous count-box regions type attainable zero, infinity, and
  non-estimability without continuity corrections or classical p-values.
* Added `ds.vertDPROC()` as zero-cost post-processing of one ordered
  disease-status by score-bin DP table. It reports the complete threshold ROC
  curve, tie-adjusted finite-snapshot AUC, and conservative simultaneous
  mechanism-noise regions without presenting them as population confidence
  bands.
* Added `ds.vertDPCapsulePlan()` as a signed, data-free dry-run of the
  server-authoritative capsule workload, coordinate/sensitivity cost and
  selected mechanism. It never materialises protected data or creates a DP
  release.
* DP contingency releases now validate and expose the sole custodian-owned
  `consistent_cell_else_exclude_v1` repeated-record rule and fail closed on
  protocol/payload drift; conflict or exclusion counts remain unreleased.
* `ds.vertDPDescribe()` now reports conservative simultaneous mean/variance
  regions obtained from the certified count/sum/sumsq noise box, bounded-moment
  geometry, and deterministic per-unit quantisation error. They are explicitly
  labelled as mechanism/grid regions rather than sampling confidence intervals.
* `ds.vertDPEpi2x2()` retains each individually estimable noisy-table effect
  when another group is empty, and high-confidence Laplace regions avoid
  cancellation by working directly with tail probabilities.
* Client survival validation now matches the server's distinct invalid-bin
  semantics with and without delayed entry.
* Added `inst/docs/dp_statistical_validity_20260801.md`, including the supported
  finite-snapshot estimands, non-claims, missing joint mechanisms, and legacy
  exact routes that cannot meet the DP reconstruction objective.

### Beaver preprocessing

* Reuse IKNP base-OT state within a DataSHIELD session for each sender,
  receiver and ring tuple. Each multiplication batch still uses a unique
  extension transcript key, which the server runtime uses for domain-separated
  PRG seeds.
* Added an OT-aware K>=3 LMM profile mode. Under IKNP preprocessing,
  `ds.vert.lmm()` avoids the repeated weighted-GLS
  golden-section profile and uses the protected moment + cluster-mean GLS
  route instead; set `options(dsvert.lmm_k3.profile_mode = "profile")` to
  force the exhaustive profile path.
* IKNP OT-extension is the sole product backend: dealer preprocessing has been
  removed because a participating-party dealer can reconstruct peer operands.

### Cleanup

* Removed generated package tarballs, vignette HTML, vignette caches,
  `demo.html`, and TeX logs from version control. Rmd sources remain as the
  editable validation documentation.
* Updated README and DESCRIPTION to describe the current product routes and
  the removed-route policy.
* Added `inst/docs/product_surface.md`, the K=2 / K>=3 route matrix used to
  distinguish product estimators from discarded historical paths.
* Removed disclosive or materially suboptimal user routes: Cox rank/K3
  prototypes, negative-binomial per-patient eta transport, multinomial OVR
  final estimator, ordinal patient-level reconstruction, and GLMM EM.
  Joint multinomial/ordinal routes still use non-exported warm starts
  internally.
* Made source installation instructions version-agnostic so they do not point
  at stale local tarball names.

## dsVertClient 1.1.0

Companion release to dsVert 1.1.0. Closes the v2 follow-up list (5/5
items shipped) and brings the client API to the J-BHI submission
state.

### New federated estimators

* **`ds.vertCox`** — federated Cox proportional-hazards regression
  (reverse-cumsum reformulation, Allison 1982 / Andreux 2020). Supports
  left-truncation (`tstart_col`), stratification (`strata_col`),
  Path B damped fixed-Fisher refinement, and `ring = 63 | 127`. Defaults
  to Ring127 (5/5 STRICT on Pima synthetic; ~2× faster Path B).
* Historical K=3 Cox Poisson-trick prototype with offset + baseline dummies
  was evaluated, then removed from the product API in favour of the
  non-disclosive profile route.
* **`ds.vertCoxDiscreteNonDisclosive`** — K=2 discrete-time
  pooled-logistic Cox with Aliasgari-Blanton 2013 share-mask gating.
* **`ds.vertMultinom`** — K-class one-vs-rest multinomial logistic.
* **`ds.vertMultinomJoint`** — compatibility wrapper for the joint Newton
  softmax route.
* **`ds.vertMultinomJointNewton`** — full softmax Newton via Ring127
  Chebyshev exp + reciprocal + Beaver vecmul on per-patient eta_k /
  mu_k shares (paper §V.A row, Bohning-bounded).
* **`ds.vertOrdinal`** — proportional-odds via cumulative-binomial
  BLUE pool + threshold correction (McCullagh-Agresti).
* **`ds.vertOrdinalJointNewton`** — Tutz 1990 §3.2 block-diagonal
  joint Newton with McCullagh §2.5 closed-form H_θθ + post-Newton
  refinement; sub-noise close-well at L2 K=2.
* **`ds.vertNB`** — negative binomial GLM via two-stage Poisson β + θ
  Newton-Raphson; `joint = TRUE` re-Poisson-fits with theta-adjusted
  weights.
* **`ds.vertNBMoMTheta`** — closed-form Method-of-Moments theta
  (Anscombe 1950 / Saha-Paul 2005), psi-free.
* **`ds.vertNBFullRegTheta`** — first-order corrected full-regression
  theta with Ring127 NR-LOG share-domain primitive (Goldschmidt 1964 +
  Pugh 2004); SUB-NOISE on K=2 (105×–109× σ-probe ratio).
* **`ds.vertLMM`** — Laird-Ware GLS closed-form with REML; supports
  random intercept + slopes via `random_slopes`; `ring = ring127`
  pipeline for STRICT closure on dense Gram matrices.
* **`ds.vertLMM.k3`** — REML 1-D profile LMM for K=3.
* **`ds.vertGLMM`** — aggregate GLMM-PQL binomial mixed model (single
  random intercept) with guarded share-domain cluster sufficient
  statistics; **`ds.vertGLMMLaplace`** / **`ds.vert.glmer`** provide the
  supported Laplace route against the `glmer(nAGQ = 1)` target.
* **`ds.vertIPW`** — propensity fit plus protected weighted outcome GLM
  using a server-side IPW column that is consumed but not returned.
* **`ds.vertLASSOProximal`** — proper LASSO via client-side
  proximal-gradient (FISTA-accelerated by default) on the normal
  equations; 4–10× tighter coefficient agreement vs post-hoc
  soft-threshold.
* **`ds.vertLASSOIter`** — standardized L1 path for Gaussian, binomial,
  and Poisson GLMs; binomial uses aggregate-score FISTA with a
  LASSO-specific 200-interval secure sigmoid, and Poisson uses
  aggregate score/Hessian prox-Newton.
* **`ds.vertLASSO1Step`** + **`ds.vertLASSOCV`** — one-step
  quadratic-surrogate LASSO with AIC / BIC / EBIC information-criterion
  λ selector.
* **`ds.vertGEE`** — generalised estimating equations (working
  exchangeable / AR1 correlation, sandwich SE) with per-call
  `binomial_sigmoid_intervals` precision control for protected binomial
  sigmoid evaluations.
* **`ds.vertMI`** — multiple-imputation wrapper with Rubin pooling.
* **`ds.vertChisq`** — two-way contingency χ² via Beaver dot product
  on one-hot shares.
* **`ds.vertDesc`** — descriptive aggregates (mean / SD / min / max /
  histogram-based quantiles).
* **`ds.vertDPCor`** — explicit same-owner, pairwise-complete bounded Pearson
  correlation from one signed sticky DP capsule vector. It is retained for
  research use but is marked ineligible for PCA.
* **`ds.vertCor` / `ds.vertPCA`** — joint complete-case correlation and
  client-only PCA from one signed `gaussian_models` artifact, including the
  cross-owner exact-GC artifact. There is no silent pairwise fallback, no
  individual score release, and no legacy Ring63 correlation route.
* **`ds.vertContrast`**, **`ds.vertWald`**, **`ds.vertConfint`**,
  **`ds.vertLR`** — inferential scaffolding (multivariate Wald,
  likelihood-ratio, contrast tests; CI helper).

### API changes to `ds.vertGLM`

* New args: `offset`, `weights` (column-name based), `ring = 63 | 127`,
  `binomial_sigmoid_intervals` (per-call binomial secure-sigmoid precision),
  `keep_session = TRUE` (session reuse for downstream LMM / GEE / Cox),
  `no_intercept`, `std_mode = "full" | "x_only"`, `data_name`, `y_var`
  (back-compat aliases).
* Auto-detects predictor → server mapping when `x_vars = NULL` (calls
  `dsvertColNamesDS` once).

### Documentation / quality

* All Rd warnings cleared: bracket-link traps, loose LaTeX macros,
  Rd parser END_OF_INPUT (`2\%`, `\%s_leq`), missing `\link{}` targets.
* All R sources ASCII-only (Greek / math symbols replaced).
* `call(name = "fn", ...)` adopted across 490 sites in 29 files
  (eliminates partial-arg-match NOTE).
* `R/zzz.R` adds `@importFrom stats / utils` for previously
  hidden-binding functions.
* LICENSE switched to DCF stub; full MIT in `LICENSE.md`
  (Rbuildignored).
* `.Rbuildignore` excludes vignette caches, scripts, docs.
* DESCRIPTION drops unused `digest` import.

### Testing

* Extended the reproducible DP statistical harness to Gaussian regression,
  complete-case correlation and PCA, including 1,000-replicate fidelity,
  mechanism-region and spectral-identity gates without inventing coefficient
  or loading confidence intervals.
* 184 client-side `testthat` checks PASS, 15 SKIP (require live DSI
  mocks or DSLite extensions covered by L3 opal-demo probes).
* `quick_impl_check.sh`: 8/8 L1 probes + Go test PASS in <3 min.
* `continuous_validation.sh fast`: L1 7/7 OK in 24 s.
* `continuous_validation.sh medium`: L1 7/7 + L2 3/3 OK (~33 min).
* Local DSLite K=2 smoke (Pima, n=132): max|Δβ| = 3.18e-04 vs lm
  (binomial GLM); max|Δβ| = 8.09e-05 vs `survival::coxph` (Cox PH);
  max|Δβ| = 8.32e-06 (Cox in LMM+Cox combined harness).

### R CMD check

* `Status: OK` (0 ERRORs, 0 WARNINGs, 0 NOTEs).

## dsVertClient 1.0.0

Initial public release.
