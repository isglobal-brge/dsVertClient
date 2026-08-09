# Formal DP design for vertically partitioned binomial and Poisson GLMs

Status: research design plus an isolated Phase-0 executable reference,
2026-08-02. Nothing in this document promotes the current legacy GLM path or
asserts that a protected GLM execution/release route is implemented.

## Phase-0 implementation boundary

The client package now contains an internal, side-effect-free Phase-0
reference in `R/formal_glm_rational.R`, `R/formal_glm_compiler.R` and
`R/formal_glm_oracle.R`. It is intentionally absent from the namespace and
does not call DSI, inspect a datasource, register a server endpoint or dispatch
to the legacy GLM implementation.

The implemented reference freezes and tests:

- a canonical signature-ready schema bound to the consortium, immutable
  snapshot, capsule, full pinset commitment and exactly two designated peers;
- verification of byte-identical Ed25519 signatures from both identities in a
  caller-supplied, already trusted pinset;
- binomial/logit and Poisson/log families, add/remove and replace-one
  adjacency, fixed padded capacity, a one-record-per-patient rule that maps
  duplicates/conflicts to the zero-weight tuple, complete-tuple zero-weight
  missingness, clipping, bounded weights, bounded offsets and bounded positive
  exposures transformed to a quantised log offset;
- a deterministic additive formula compiler, treatment coding over signed
  categorical domains, explicit intercept semantics, canonical coefficient
  order, full positive ridge, coefficient boxes, fixed zero start, fixed step
  and iteration count;
- arbitrary-precision signed rational arithmetic over OpenSSL BIGNUM, rigorous
  Taylor/geometric-tail intervals for `exp` and `log`, and generated monotone
  dyadic piecewise-linear link tables with rational continuity, slope and
  uniform-error certificates;
- an exact-rational projected full-gradient oracle and contraction/sensitivity
  certificate for the reference algorithm; and
- a deliberately separate ordinary central-GLM comparator that returns typed
  `non_identifiable` conditions for rank deficiency, separation/no finite fit
  and the all-zero Poisson boundary. Full ridge and the public box make the
  formal surrogate itself unique in these cases, so they are not silently
  replaced by Firth, pseudoinverse or another estimand.

Only unique additive registered columns are accepted. The compiler rejects
`.`, transformations, interactions, splines, data-dependent factor discovery,
unregistered/reference-invalid categorical domains, non-treatment contrasts,
unbounded fields, repeated-record aggregation (pending a separately proved
patient-level rule), non-positive ridge, unconstrained coefficients,
non-contractive optimisers, unsupported missingness/clipping/weight/offset
semantics and public numeric plans outside the Phase-0 certified domain.
Equivalent additive formula order and registry order canonicalise to identical
bytes and hashes.

The generated numeric plan is a conservative public **candidate** ring plan.
Its `rho = 0` statement applies only to the arbitrary-precision rational
reference. Every compiled object sets `production_release_ready = false`,
`production_fixed_point_certificate_pending_phase1 = true` and DP calibration
to `phase2_pending`. Phase 1 must still compile and verify the complete
fixed-point MPC DAG and its non-zero outward `rho`; Phase 2 must implement the
sticky one-global-draw protected coefficient release; Phase 3 must add the
typed client surface. No claim below should be read as completing those
phases.

## Decision

The first formally defensible migration should be **one-shot output
perturbation of a deterministic, fixed-iteration, projected full-gradient
algorithm**. The algorithm must fit a public, box-constrained, fully
ridge-regularised binomial or Poisson surrogate on one patient-level,
capacity-padded vertical cohort. Its complete coefficient vector is then
perturbed once, inside the two designated computation peers, and opened only
after noise and final public projection.

This is deliberately not the classical recipe “compute an exact unbounded MLE
and add noise”. The privacy proof covers the exact finite algorithm that dsVert
will run, including its public iteration count and a conservative fixed-point
error enclosure. It therefore avoids making privacy depend on an unproved
convergence tolerance or on a data-dependent stopping rule.

The choice is the smallest credible step from the current architecture:

- it uses one DP vector release rather than composing noise over iterations;
- it needs only first-order rowwise products, a bounded link approximation,
  exact comparisons, exact checked multiply/truncate and box projection;
- it is linear in cohort capacity and coefficient count per iteration, rather
  than requiring a private Hessian or Newton solve; and
- it can share the capsule identity, sticky ledger, pinned-peer control plane,
  typed source transport, dynamic-ring planner and exact-GC core already in
  the repository.

Objective perturbation is a worthwhile later utility project, not the first
promotion target. Classical objective-perturbation privacy depends on a
Jacobian argument and an exact minimiser; modern approximate-minimum variants
add another private release. DP-SGD is more general but requires per-example
clipping, noise and privacy accounting across many iterations. Polynomial or
moment releases are attractive for some specialised Bayesian approximations,
but logistic and Poisson models do not have a fixed finite sufficient statistic
for arbitrary continuous designs.

## Claim boundary after implementation

If every mandatory item and acceptance test in this document passes, dsVert
could accurately claim:

> For a registered immutable snapshot, signed patient-level bounds and a
> registered binomial or Poisson model, dsVert returns a sticky
> `(epsilon, delta)`-DP estimate of a stated box-constrained,
> ridge-regularised, quantised finite-iteration GLM surrogate. The protected
> vertical rows, intermediate scores, iterates and unnoised coefficients are
> not opened to the analyst relay. The numerical certificate proves the public
> no-wrap domain and bounds deterministic approximation error. This claim
> assumes two pinned, semi-honest, non-colluding computation peers.

It could not claim any of the following without additional work:

- that no adversary can ever guess or reconstruct any source value;
- malicious-secure MPC or protection when the two computation peers collude;
- an ordinary unbounded, unpenalised maximum-likelihood estimate;
- valid classical standard errors, Wald/LR tests, p-values or population
  confidence intervals;
- negligible clipping, ridge, box, quantisation, link-approximation or DP bias
  for every possible dataset; or
- production Armadillo/Opal WAN performance merely from local DSLite tests.

Differential privacy limits the change in the output distribution caused by
one protected unit. It is a formal inference-stability guarantee, not a
literal distribution-free impossibility of guessing a record.

## Threat and execution model

### Parties

- `K >= 2` custodians hold vertically partitioned columns for the aligned
  cohort. `K` is not the number of MPC parties.
- Exactly two distinct, custodian-designated, identity-pinned peers execute the
  arithmetic and joint release.
- The analyst/client is an untrusted opaque relay. It may reorder, replay,
  delay or drop messages but must not receive a plaintext input share,
  intermediate value, noise draw or unnoised coefficient.
- Each computation peer may know its own source columns. For every other
  source it sees at most one uniformly masked additive share and a
  purpose-bound semi-honest protocol view.

The confidentiality claim requires at least one of the two computation peers
not to collude with the other. Pinning removes the analyst relay from the
trusted computing base; it does not turn the existing Yao/OT implementation
into malicious-secure MPC. Authentication, signed purpose binding and KOS
consistency checks may cause an abort, but do not prove correct behaviour by a
malicious peer.

### Views that must be separated

There are two complementary guarantees:

1. MPC transcript privacy prevents the relay or one non-colluding computation
   peer from learning protected intermediate values.
2. Differential privacy constrains the one public coefficient release.

Neither guarantee substitutes for the other. Secret sharing an exact model
without DP does not make adaptive releases safe; adding DP noise after opening
the exact model is too late.

### Capsule semantics

Every supported GLM artifact is registered in the server-authoritative capsule
manifest before protected materialisation. Its model identifier binds the
family, formula grammar, term order, levels/interactions, bounds, contribution
rule, regulariser, box, numeric grid, link approximation, optimizer and DP
mechanism. An analyst cannot extend an immutable capsule with a new protected
formula.

If a capsule contains `M` registered GLMs, all coefficient blocks form one
joint vector query and are sampled/opened once. Repeating a fit, changing
equivalent formula syntax or requesting deterministic post-processing returns
the same sticky block and has zero additional privacy loss. A genuinely new
model workload, snapshot or privacy epoch is a different capsule release and
composes with earlier releases. Its outer allocator commit consumes one
lifetime unit; composition history can deny capsule `N+1` but is never a
request/operation quota for replay.

## Privacy unit, adjacency and contribution contract

### One aligned patient is one unit

All rows belonging to one patient across all `K` sources are collapsed before
the model contribution is formed. The signed repeated-record rule must map the
entire patient record set deterministically to at most one bounded tuple

`z_i = (a_i, x_i, y_i, o_i),  i = 1, ..., C`,

where `C` is the public capsule capacity:

- `a_i` is a non-negative analysis weight in `[0, W]`; `a_i = 0` represents a
  padded, missing, conflicting or otherwise excluded unit;
- `x_i in R^p` is the complete design row in signed coefficient order;
- `y_i in {0,1}` for binomial or `y_i in {0,...,Y_max}` for Poisson; and
- `o_i in [o_min,o_max]` is an optional bounded offset, with zero as the
  default.

`W`, `Y_max`, feature bounds and offset bounds are custodian-owned public
policy, not analyst-supplied sensitivity claims. The initial implementation
should default to `W = 1`. Supporting case/frequency weights is safe only when
the patient-level collapse and `W` are signed and included in every bound
below.

For repeated Poisson observations, a defensible rule is to sum the patient's
counts and exposures and then clip once to signed patient-level caps. For
binomial outcomes, a defensible rule is a signed concordance/tie-break policy.
The exact rule may differ by registered artifact, but it must be deterministic,
patient-level and invariant to row order. Missingness or disagreement maps to
`a_i = 0`; it must not select a visible error branch.

### Adjacency

The manifest chooses exactly one of:

- **add/remove in a fixed padded capacity:** adjacent inputs differ by toggling
  one complete patient tuple between its bounded value and the all-zero dummy;
  or
- **replace-one:** adjacent inputs have the same capacity and replace all
  columns/records of one patient by another admissible patient tuple.

The change is joint across all vertical owners. It is not one column at one
site. Let `gamma = 1` for add/remove and `gamma = 2` for the conservative
replace-one triangle bound. No certificate may silently reuse an add/remove
sensitivity for replace-one.

### Bounded design

For each coefficient coordinate, the signed manifest gives

`|x_ij| <= A_j` and `R = sqrt(sum_j A_j^2)`.

An intercept is represented by `x_i0 = 1`, so its contribution is included in
`R`. Categorical levels, reference cells, interactions and transformations
are a fixed, versioned design grammar. Data-dependent level discovery,
standardisation, knot selection or variable selection is not part of this
mechanism.

Features are clipped and quantised after patient-level collapse. The clipped,
quantised design is the protected finite-snapshot estimand. The certificate
must separately report its maximum pre-quantisation error; it must not describe
the result as a fit to unavailable unclipped values.

## Statistical target

### Public coefficient box and full regularisation

The coefficient domain is the public convex box

`B = product_j [-B_j, B_j]`.

The diagonal regulariser `Lambda` must satisfy

`lambda_min = min_j Lambda_jj > 0`.

For the first formal route this includes the intercept. Leaving any coordinate
unpenalised removes the public global strong-convexity bound unless a separate,
data-independent curvature lower bound is proved. The first version therefore
estimates a fully penalised target and must label it as such. Ridge and the box
make the target unique even under collinearity, logistic separation or sparse
Poisson outcomes; they do not establish classical unpenalised identifiability.

Public hyperparameters may be chosen from scientific prior knowledge or a
pre-registered policy. Selecting `Lambda`, the box, iteration count, step size
or link grid after inspecting protected unnoised fit quality is not free. A
pre-registered grid of models may be included as one joint capsule vector,
with its resulting joint sensitivity and dimension.

### Bounded linear predictor

Let

`O = max(abs(o_min), abs(o_max))`

and

`M_eta = O + sum_j A_j B_j`.

Then every admitted and dummy unit satisfies `|o_i + x_i' beta| <= M_eta` for
all `beta in B`. This public bound is essential for Poisson: the exponential
score and curvature are otherwise unbounded. It also supplies a finite domain
for certifying the binomial link approximation and ring headroom.

### Convex link surrogate

The first exact-MPC implementation should not use the old DCF spline path.
Instead, each family has a public continuous monotone piecewise-linear mean
function `mu_tilde` on `[-M_eta,M_eta]`:

- binomial approximates `sigmoid(eta)` and is range-confined to `[0,1]`;
- Poisson approximates `exp(eta)` and is range-confined to a certified
  `[L_mu,U_mu]`, with `0 < L_mu <= U_mu`.

All knots and dyadic coefficients are public and committed by hash. Exact-GC
signed comparisons and muxes choose a segment; checked multiplication and
exact truncation evaluate it. Define the surrogate row loss by integrating
the score,

`d l_tilde(y,eta) / d eta = mu_tilde(eta) - y`.

Because `mu_tilde` is continuous and non-decreasing, this defines a convex,
continuously differentiable loss even though its second derivative changes at
public knots. The fitted objective is

`F_D(beta) = (1/C) sum_i a_i l_tilde(y_i, o_i + x_i' beta)`

`            + (1/2) beta' Lambda beta`.

The public capacity `C`, not the private number of complete cases, is the
normaliser. This is equivalent to a capacity-scaled likelihood with a stated
ridge penalty. It avoids a private reciprocal and makes add/remove sensitivity
well defined. It may penalise more strongly when many units are excluded, so
the normalisation must be reported and must not be called an ordinary
complete-case average likelihood.

Let `S_mu` be the maximum non-negative segment slope. The objective is
`m`-strongly convex and `L`-smooth on `B`, with public bounds

`m = lambda_min`

`L = lambda_max + W * R^2 * S_mu`.

The maximum patient score norm is bounded by

`G_binomial = W * R`,

`G_poisson  = W * R * max(U_mu, Y_max)`.

The latter uses non-negativity: `|mu_tilde-y| <= max(U_mu,Y_max)`. The
fixed-beta gradient difference between adjacent objectives is therefore at
most

`s = gamma * G_family / C`.

Every quantity in these bounds is public, signed and machine-verifiable.

## Deterministic optimizer and sensitivity proof

### Fixed algorithm

Choose a public dyadic step size `alpha`, public iteration count `T >= 1` and
public starting point `beta_0 in B`, normally zero. Require

`0 < alpha <= 2 / (m + L)`.

For exactly `T` iterations execute

`beta_(t+1) = Project_B(beta_t - alpha * grad F_D(beta_t)).`

There is no early stopping, line search, data-dependent damping, convergence
bit, private iteration counter, adaptive precision or fallback estimator. Row
blocks may be streamed and reduced in a fixed public order, but that is
deterministic full-gradient evaluation, not stochastic minibatch training.

### Contraction

The ideal projected-gradient map is contractive in Euclidean norm. A valid
public contraction factor is

`q = max(abs(1 - alpha*m), abs(1 - alpha*L)) < 1`.

Projection onto `B` is non-expansive. The fixed-point implementation must
publish a uniform per-update enclosure `rho` such that, for every admissible
dataset and every lattice point in `B`, the decoded MPC update differs from
the ideal real-valued update by at most `rho` in L2 norm. `rho` is derived by
outward interval arithmetic over the complete operation graph; it is never an
empirical tolerance or an analyst assertion.

For two adjacent executions starting at the same `beta_0`, let `d_t` be their
decoded coefficient distance. Then

`d_(t+1) <= q*d_t + alpha*s + 2*rho`, with `d_0 = 0`.

Consequently the implemented deterministic query has certified global L2
sensitivity

`Delta2_real <= (alpha*s + 2*rho) * (1 - q^T) / (1 - q)`.

This recurrence is the primary phase-one privacy proof. It covers the finite
algorithm actually executed rather than assuming an exact minimiser. If the
common coefficient lattice has `S_beta` integer steps per coefficient unit,
the joint-DP planner must also account for the componentwise release
quantizer. A bare `ceil(S_beta * Delta2_real)` is not generally safe: several
coordinates can cross adjacent floor boundaries simultaneously even when the
scaled Euclidean distance is below one. For a signed-floor/downsampling map
from a finer integer lattice by divisor `d`, a simple machine-checkable bound
for a `p`-coordinate block is

`Delta2_steps <= ceil_up(Delta2_source_steps / d + sqrt(p))`.

The implementation may use a tighter exact integer-lattice optimisation, but
it must prove that bound for the committed quantizer. It may not omit the
rounding term. When no such proof is available, the translated public-box
diameter remains the safe fallback.

All model blocks in one capsule use the same coefficient lattice. For `M`
registered models that one patient can affect simultaneously, a safe stacked
bound is

`Delta2_capsule <= sqrt(sum_m Delta2_steps_m^2)`.

A future implementation may replace this with a tighter proved joint bound,
but may not take the maximum over blocks or pretend they are disjoint.

For a pure-DP coordinatewise discrete-Laplace alternative, a conservative
conversion is

`Delta1_capsule <= sqrt(P) * Delta2_capsule`,

where `P` is the total released coefficient dimension. This dimensional cost
makes the exact discrete Gaussian the preferred point-estimate mechanism when
a non-zero `delta` is available.

### Accuracy of the deterministic target

Let `D_B` be the public Euclidean diameter of `B`. The distance to the exact
surrogate minimiser after `T` iterations is bounded by

`optimization_error <= q^T * D_B + rho * (1 - q^T) / (1 - q)`.

If the link approximation satisfies the certified uniform bound

`sup_eta |mu_tilde(eta) - mu(eta)| <= e_mu`,

then the gradient perturbation caused by the link alone is at most
`W*R*e_mu`, and strong convexity gives the conservative minimiser displacement

`link_coefficient_error <= W * R * e_mu / m`.

The same pattern applies to a separately derived feature/outcome quantisation
gradient enclosure. There is no distribution-free bound from the regularised,
box-constrained, clipped target to an ordinary unbounded MLE when the latter
may not exist or may be ill-conditioned. Ridge, box and clipping effects must
therefore be reported as estimand changes, not folded into a misleading
“numerical error” number.

## One joint DP opening

### Required release flow

1. Every source uploads its fixed-capacity, fixed-schema patient contribution
   as two independently encrypted additive-share streams, one for each pinned
   computation peer. Slot alignment is privately attested; no patient
   identifier, patient-derived public hash or occupancy vector enters the
   manifest. The relay sees only fixed public framing and ciphertext.
2. The two peers execute the fixed optimizer on shares. No score, objective,
   complete-case count, Hessian, residual, convergence state, intermediate
   coefficient or validity bit is opened.
3. The two private noise-root contributions are combined inside a
   producer-bound joint sampler/finalizer to generate exactly one sticky
   coefficient-noise vector. Neither peer nor the relay learns the full noise
   vector or the unnoised coefficient vector.
4. Noise is added before reconstruction. The finalizer signed-decodes once,
   applies the committed public coefficient-box projection, and opens only the
   DP coefficient vector and public certificate.
5. Both peers sign byte-identical release bytes before the client accepts the
   object. Retries, reconnects and restarts replay those bytes and never draw
   again.

Final box projection is DP-safe post-processing and cannot increase Euclidean
error relative to a true deterministic coefficient vector already in `B`.
Projection metadata must be derivable from the released value; no hidden
pre-projection bit is returned.

### Mechanism

The preferred mechanism is the existing rationally planned, sticky discrete
Gaussian, calibrated to `Delta2_capsule` and the capsule's fixed
`(epsilon,delta)` allocation. Its exact/integer implementation and any bounded
support or total-variation transfer slack must be included in `delta_total`.
The primary reason to retain an integer sampler is that naive floating-point
sampling can invalidate a formal mechanism; Canonne, Kamath and Steinke give
an exact discrete-Gaussian construction and analyse its DP behaviour.

Using `rho_priv` for the privacy-conversion parameter (distinct from the
numeric update error `rho`), a conservative zCDP calibration chooses

`rho_priv + 2*sqrt(rho_priv*log(1/delta_core)) <= epsilon`

and

`sigma^2 >= Delta2_capsule^2 / (2*rho_priv)`.

The rational planner must round `rho_priv` downward, `sigma^2` upward and
verify the conversion inequality without binary floating-point shortcuts.
Sampler tail/TV transfer consumes `delta_impl`, with
`delta_core + delta_impl <= delta`. The certificate records all rational
numerators/denominators and the exact theorem/version used.

The current independent-full-draw fallback is not the utility target for GLM:
two full independent Gaussian draws preserve the one-hidden-seed privacy
argument but double nominal variance. The promoted GLM route should combine
the two private seeds in exact GC and emit one globally calibrated hidden draw.
Existing rational planner, exact-sampler reference code, sticky derivation and
certificate tests can be reused, but a producer-bound one-draw two-peer
finalizer is new work and must pass E2E tests before promotion.

If a deployment uses discrete Laplace, it must use the jointly certified
global L1 sensitivity and account for all implementation-tail slack. It cannot
add independently calibrated scalar noise using epsilon once per coefficient.

### Sticky identity

The canonical release identity must bind at least:

- consortium, immutable snapshot, patient capacity, privacy epoch and complete
  pinset;
- all registered model identifiers and their canonical design order;
- adjacency, patient collapse, missingness, weight, bounds and offset policy;
- `C`, `W`, `R`, `Y_max`, `M_eta`, box and full `Lambda`;
- link knots/coefficients/hash, `e_mu`, `S_mu` and integration convention;
- `alpha`, `T`, coefficient grid, ring plan, arithmetic DAG and `rho`;
- sensitivity proof version and exact rational sensitivity bound; and
- mechanism, `(epsilon,delta)`, sampler/circuit versions and final projection.

Equivalent client syntax maps to the same registered artifact. A noise-root
rotation may select a new candidate only before the first valid START claims an
instance and remains subject to the outer lifetime gate. After the claim, loss
must restore and continue the exact instance or fail closed; after publication
only replay/restore is allowed. It cannot create a new composed release for the
old capsule or block unrelated post-processing.

## Transcript and numerical requirements

### Transcript invariance

The following are public and fixed before source materialisation:

- `C`, total `p`, per-source column counts, row-block size and number of
  blocks;
- `T`, link segment count and comparison topology;
- ring/container/fractional widths and every circuit chunk shape; and
- the number and order of source uploads, peer phases, receipts and release
  bytes.

All-missing, all-zero, separated, dense and boundary-valued protected inputs
must induce the same logical phase/shape transcript. There is no visible
small-cell guard, overflow bit, convergence bit, line-search count, reciprocal
failure or private-domain error. Public schema/authentication failures may
abort. A private value at a signed bound is valid input and cannot choose a
different protocol path.

The random sampler may consume a data-independent random-bit stream, but its
externally visible framing and bounded-support certificate must follow the
committed mechanism. No random draw, seed-derived diagnostic or number of
rejection iterations is exposed.

### Ring planning

The producer compiles public interval bounds for every intermediate:

- quantised design, coefficient and offset products;
- linear-predictor accumulation;
- segment selection and link evaluation;
- residuals, per-row scores and capacity-wide gradient accumulation;
- regularisation, step update and projection; and
- coefficient plus complete supported noise range before final projection.

Ring63/Ring127 are fast paths only when those bounds prove all raw-product,
truncated-output and accumulator headroom. Otherwise the producer creates
fresh shares under the smallest certified dynamic ring, through Ring4096.
There is no reinterpretation of existing Ring127 shares and no fallback to
local truncation, DCF comparison, plaintext arithmetic or modular wrap. A
bound that cannot fit Ring4096 returns a typed public
`numeric_backend_unrepresentable` specification state before protected
execution; it is never presented as a model result.

Exact truncation still introduces a known fixed-point approximation to the
real target. The certificate reports the rounding convention, fractional
width, per-stage outward error and accumulated `rho`. “Exact GC” means exact
execution of the stated integer operation, not exact evaluation of sigmoid,
exponential or an ordinary real-valued GLM.

## Comparison of candidate mechanisms

| Route | Privacy object | Main requirements | Vertical-MPC cost | Statistical/utility issue | Decision |
|---|---|---|---|---|---|
| Classical exact-ERM output perturbation | one noisy minimiser | global strong convexity, bounded score and an exact or tightly certified solve | iterative GLM plus one sampler | simple sensitivity and one noise vector, but an approximate solver needs an additional proof | useful reference, but phase one proves the fixed algorithm directly |
| **Fixed-iteration algorithmic output perturbation** | one noisy deterministic `beta_T` | public box/full ridge, fixed `T`, contractive step, uniform arithmetic-error bound | `T` full gradients plus one sampler | conservative sensitivity; explicit finite-iteration and surrogate target | **selected phase-one route** |
| Objective perturbation | minimiser of a randomly perturbed objective | differentiability/Jacobian conditions, exact minimiser or an approximate-minimum correction, bounded Poisson gradient/curvature | similar private optimisation plus secure objective noise | can outperform output perturbation for convex linear learners, but proof/solver coupling is materially harder | phase-two research candidate |
| DP-GD/DP-SGD | composition of noisy clipped gradients | per-patient clipping, shared sampling, noise each step, accountant and private hyperparameter handling | many cryptographic rounds and random accesses | versatile; iteration noise and tuning can dominate a small GLM | not phase one |
| Polynomial/functional mechanism | noisy polynomial objective coefficients | fixed expansion degree, coefficient sensitivity, convexity repair | one release then local solve; coordinate count grows combinatorially | logistic/Poisson likelihood is approximated; noisy objectives may be non-convex/unbounded | only a separately named low-dimensional method |
| Approximate sufficient-statistic GLM | noisy finite moments, then noise-aware inference | declared moment basis/order, joint sensitivity and a distinct inferential model | cheap repeated post-processing after one potentially large release | not an exact sufficient-statistic representation for arbitrary continuous logistic/Poisson design | future Bayesian artifact, not a `ds.vertGLM` compatibility shortcut |

Chaudhuri, Monteleoni and Sarwate prove the classical regularised-ERM output
sensitivity pattern and objective-perturbation route under explicit convexity
and differentiability conditions. Redberg, Koskela and Wang show that improved
objective perturbation can be competitive for convex GLMs, while also noting
that exact minimisation is central and approximate minima require a corrective
private release. Song et al. establish strong results for DP gradient methods
on GLMs, including the fact that clipped gradients optimise a specifically
modified objective. These results motivate later work; none can be copied
verbatim without mapping its adjacency, constraint, optimisation and numeric
assumptions to this implementation.

Recent data-dependent sufficient-statistic perturbation reconstructs a
competitive **approximate** logistic objective from privately released
marginals. That is meaningful phase-five evidence, but it reinforces rather
than removes the distinction here: it is a separate marginal-release workload
and approximation, not an exact finite sufficient-statistic identity for the
ordinary logistic likelihood.

For a total-degree-`r` polynomial objective over `p` coefficients, the naive
monomial basis has `choose(p+r,r)` terms before exploiting symmetry or sparsity.
Second order is already quadratic in `p`, while accurate exponential/logistic
approximations can require higher order over a wide `M_eta` domain. This can be
excellent for a small fixed design and poor for a broad vertical formula
registry; its dimension and approximation bias must be benchmarked rather than
described as a universally faster sufficient-statistic route.

## Scientific output contract

The phase-one result object contains:

- the DP coefficient point vector in canonical term order and coefficient
  scale;
- family, link surrogate, patient-level estimand, capacity normalisation,
  complete-case rule, bounds, offset/weight policy, full ridge and box;
- capsule receipt, pinned peers, adjacency and fixed `(epsilon,delta)`;
- ring/backend, fractional grid, public magnitude proof, `rho`, `q`, `T`,
  deterministic optimization bound and link/quantisation bounds;
- global sensitivity and sampler/implementation-delta certificate; and
- a simultaneous mechanism-only coefficient-noise radius when supported by
  the exact sampler certificate.

It does not contain exact or merely hidden versions of `n`, log likelihood,
deviance, AIC, score, Hessian, fitted patient values, residuals, influence,
rank, separation flags, convergence traces or server error text. Any such
quantity must be either deterministic post-processing of the released
coefficient plus public data or a separately designed coordinate in the one
joint capsule mechanism.

The first promoted object is **point-estimate only**. A mechanism-noise radius
is not a sampling confidence interval. Ordinary GLM covariance computed as if
the released coefficient were an exact MLE is not valid. A later inferential
phase may jointly release a bounded DP Hessian/score-covariance artifact and
validate a noise-aware bootstrap or posterior. That phase changes the vector
sensitivity and scientific target and must not be inferred from this design.

For biomedical/epidemiological use, this phase can support bounded descriptive
effect estimation and prediction research where a penalised point target is
acceptable. It is not yet a confirmatory epidemiology suite: exposure-offset
semantics, overdispersion, robust/clustered covariance, missing-data methods,
hypothesis tests and calibration all require separate validation. Poisson does
not imply that overdispersion is absent, and this route does not silently
substitute negative binomial, quasi-Poisson, Firth or sandwich inference.

## Complexity and scaling

Let `p` be total coefficients across all vertical owners, `C` public capacity,
`J` link segments, `T` iterations and `K` source custodians. With an oblivious
balanced segment selector, the arithmetic work is approximately

`O(T * C * (p + log J))` time,

plus `O(p)` projection per iteration and one `O(p)` joint noise release. The
dominant terms are the `x*beta` and `x*residual` products. Persistent input
shares require `O(Cp)` storage per computation peer. Row-block streaming can
hold `O(Bp + p)` working arithmetic state for public block size `B`.

The realised valid complete-case count `n` satisfies `n <= C` but remains
protected. Padding deliberately makes runtime and transcript scale with `C`,
not with private `n`. For multiple registered models, replace `p` by the sum
of their compiled coefficient/design widths unless the compiler proves and
reuses an identical shared term block.

Source transport is `O(Cp)` protected values in total, encrypted once to each
of the two computation peers, plus `O(K)` manifest/fan-out control. Since
`sum_k p_k = p`, increasing `K` does not multiply the core optimizer by `K`;
it increases source fan-out, signatures and alignment work. The source upload
occurs once per capsule. Iterations run over persistent peer-local shares and
the peer channel, not through one DSI call per row, coefficient or iteration.

This is where “batching” needs precise language. Statistical minibatching is
unnecessary and undesirable for the selected sensitivity proof. Fixed public
row **chunking** is still necessary for bounded memory, exact-GC vector limits,
backpressure and efficient DSI source transfer. Chunking preserves exactly the
same full gradient and fixed transcript.

The design scales linearly in `C` and `p` per iteration, but cryptographic
constants are substantial. It cannot promise that every enormous dataset is
fast. Practical optimisation should prioritise:

- one source upload and persistent share store per capsule;
- direct named asynchronous DSI fan-out and immutable offset replay;
- long-lived authenticated peer sessions and circuit/OT setup reuse within the
  one purpose-bound job;
- Ring127/f50 vectorised OT multiplication only when public headroom proves it,
  with direct-wide chunks as the certified fallback;
- an oblivious binary PWL selector rather than a linear scan over segments;
- fused eta/residual/gradient row blocks with deterministic reduction order;
- bounded worker concurrency, spool/backpressure and no relay payload copy;
- public tuning of `J`, fractional precision and `T` against certified error
  rather than maximal precision by default; and
- published benchmark surfaces over `C`, `p`, `J`, `T`, ring width and `K`.

Compared with private Newton/IRLS, full gradient avoids `O(Cp^2)` Hessian work
and `O(p^3)` solves, at the cost of more iterations. Acceleration, quasi-Newton
updates or objective perturbation should be evaluated only after the reference
path is correct; their state and sensitivity proofs are more complex.

## Existing backend: reuse and missing work

| Component | Current evidence | GLM use | Required work before promotion |
|---|---|---|---|
| Capsule manifest, pinset, sticky ledger and cross-signed receipts | implemented architecture and tests | bind the complete registered model and replay one release | add the normative GLM artifact schema and stacked sensitivity certificate |
| Typed source streams and two-recipient encrypted chunks | implemented for current capsule source transport | upload fixed-capacity rowwise feature/outcome/offset/weight shares once | add a producer-owned row-block schema at the compiler-selected common ring; prove private alignment and fixed shape for `K=2..5+` |
| Dynamic-ring planner and canonical multiprecision encoding | exact core supports Rings 63–4096 | choose headroom for the complete GLM DAG | extend planning from isolated operations to eta/link/gradient/update/noise end to end |
| `compare-signed` exact-GC core | implemented/property-tested, but not advertised E2E | PWL segment selection and box projection | producer-minted purpose adapter, mux support, signed receipts and E2E capability |
| checked multiply/truncate | direct-wide and Ring127/f50 hybrid cores exist | eta, PWL and gradient products | producer-bound GLM handles, accumulator bounds, streaming reduction and E2E oracle |
| exact floor/ties-to-even truncation | core primitive exists | deterministic fixed-point semantics | choose one versioned convention and include every rounding term in `rho` |
| Legacy DCF sigmoid/exp/spline and monolithic GLM iteration | explicitly residual/not promoted | none | do not route the formal GLM through these paths |
| Sticky rational discrete-Gaussian planner/sampler code | implemented local mathematical component with certificate tests | one final L2-calibrated coefficient vector | one producer-bound, one-global-draw, two-peer hidden sampler/finalizer, dynamic-ring headroom and productive DSI E2E |
| DSI direct fan-out, immutable offsets, spool and backpressure | implemented/tested locally, with connector-specific claim limits | source/control transport and opaque relay | GLM multiprocess/reconnect tests plus real deployment Opal/Armadillo validation |

The current runtime capability explicitly reports core exact comparison but no
comparison E2E and `workload_glm_e2e_verified = false`. The repository's exact
arithmetic inventory also records local truncation and DCF residues in the
legacy GLM route. Therefore the existing GLM must not be relabelled DP-safe by
adding a client dispatch branch; the producer-bound variable-width path and
one joint opening are mandatory.

Repository evidence inspected for this boundary:

- `dsVert/inst/docs/exact_gc_residual_inventory.md`;
- `dsVert/inst/dsvert-mpc/k2_exact_gc_capability.go` and
  `dsVert/inst/dsvert-mpc/k2_exact_gc_core.go`;
- `dsVert/inst/dsvert-mpc/k2_joint_dp_discrete_gaussian_v1.go`;
- `dsVert/inst/docs/biomedical_capsule_local_materializer.md`;
- [reusable capsule contract](reusable_dp_capsule_contract.md); and
- [DSI transport audit](dsi_transport_audit_20260731.md).

## Implementation phases

### Phase 0 — freeze the theorem and artifact schema

1. Specify add/remove and replace-one artifacts separately.
2. Freeze patient collapse, feature grammar, coefficient order, offsets,
   weights, capacity normalisation, full ridge and box semantics.
3. Generate monotone dyadic PWL tables with high-precision error and slope
   certificates for both families.
4. Implement a pure public compiler that derives `R`, `M_eta`, `G`, `m`, `L`,
   `q`, the full interval DAG, ring plan, `rho` and capsule sensitivity.
5. Obtain independent mathematical review of the recurrence and discrete
   mechanism calibration before protected execution is connected.

Exit criterion: the same canonical manifest bytes deterministically produce
the same complete proof certificate, with no private or analyst-controlled
sensitivity fields.

### Phase 1 — plaintext oracle and producer-bound exact arithmetic

1. Implement an arbitrary-precision rational reference optimizer.
2. Add producer-minted row-block handles bound to capsule, snapshot, peer set,
   source digests, numeric plan and destination CAS.
3. Promote exact comparison/mux, checked multiplication/truncation, fixed-order
   accumulation and box projection for that handle only.
4. Add Ring63/Ring127 fast paths and automatic fresh-share replanning through
   Ring4096.
5. Keep every result sealed; no public coefficient opening exists in this
   phase.

Exit criterion: two-peer reconstruction inside the test harness agrees with
the reference within the published `rho` for every certified backend, and no
legacy arithmetic command appears in a GLM trace.

### Phase 2 — one sticky coefficient release

1. Stack every registered GLM coefficient block in canonical order.
2. Connect the exact rational L2 sensitivity to the discrete-Gaussian planner.
3. Implement the one-global-draw hidden two-peer sampler/finalizer, final
   projection and byte-identical cross-signed persistence.
4. Bind replay, pre-START candidate rotation, the irrevocable instance claim and
   privacy epoch semantics to the capsule ledger without a request counter; a
   new capsule still consumes one outer lifetime unit.
5. Open only the DP coefficient block and public certificate.

Exit criterion: crash/restart/reconnect returns identical bytes. A new epoch may
create a different candidate only before the first valid START claim; after the
claim, the exact instance continues/restores or fails closed, and after
publication only exact replay/restore is allowed.

### Phase 3 — client surface and compatibility

1. Add an explicit registered-artifact adapter for binomial/Poisson
   `ds.vertGLM` without changing the Gaussian capsule contract.
2. Return a new typed object that cannot be mistaken for ordinary `glm`.
3. Implement deterministic coefficient-scale conversion and public
   prediction helpers.
4. Reject unsupported offsets, weights, formula terms or inference requests
   before a server call; never fall back to the legacy GLM.
5. Remove only client/server GLM calls proven unused or disclosive after a
   separate compatibility inventory and deprecation decision.

Exit criterion: every advertised combination reaches only the capsule path,
and every unsupported combination gives a typed scientific-specification
message rather than an exact legacy result.

### Phase 4 — validation and performance promotion

Run the complete acceptance matrix below on all supported platforms, then
publish benchmark and statistical-validation artifacts. Production claims for
a connector require a real deployment smoke, not merely an installed S4 method
or localhost mock.

### Phase 5 — separately reviewed inference and utility research

Evaluate objective perturbation against the reference route. Separately design
a joint DP information/score-covariance artifact and mechanism-aware inference.
Consider a noise-aware approximate-sufficient-statistic Bayesian GLM only under
a distinct method name and validation contract. None is a silent update to the
phase-one estimand.

## Mandatory acceptance tests

### Formal sensitivity and mechanism

1. Exhaustively enumerate small fixed-point domains for both families,
   `C=1..4`, `p=1..3`, boundary weights/offsets/outcomes, add/remove and
   replace-one adjacency. The observed coefficient distance after every
   `T=1..small_T` is no greater than the exact rational certificate.
2. Property-test larger random adjacent pairs, including one patient changing
   simultaneously at all `K` sources and every bound endpoint.
3. Independently recompute `m`, `L`, `q`, `rho`, per-model and stacked
   sensitivity from canonical manifests; reject any downward-rounded field.
4. Verify the exact discrete-Gaussian distribution/planner against primary
   test vectors and exhaustive small parameters; verify the published
   likelihood/`delta` inequality and all sampler-tail/TV transfer terms.
5. Prove sticky replay for equivalent syntax, retries, process restarts and
   concurrent first use. Prove distinct model/snapshot/epoch domains do not
   reuse a random stream.
6. Verify that `M` stacked models use a joint norm bound, not max sensitivity,
   per-coordinate epsilon, or a hidden sequential release.
7. Exhaust signed-floor/downsampling boundaries in multiple coordinates. In
   particular, reject a certificate based only on
   `ceil(S_beta * Delta2_real)` when simultaneous boundary crossings require
   the additional integer-lattice rounding term.

### Numerical correctness

1. Check every PWL lattice point and every segment boundary against a
   high-precision sigmoid/exponential oracle; prove monotonicity, continuity,
   range, `S_mu` and `e_mu`.
2. Compare every MPC iteration with the arbitrary-precision rational oracle on
   no-noise fixtures and randomized boundary cases. Check eta, residual,
   gradient and beta only inside the test harness; none becomes a remote
   endpoint.
3. Exercise Ring63, Ring127/f50, direct-wide Ring128+, Ring513 and at least one
   multiprecision case above Ring2048. Force every headroom boundary and prove
   either certified equality/error or typed preflight rejection—never wrap.
4. Verify floor/ties-to-even negative values, segment knots, coefficient-box
   endpoints, large `C` accumulators, `Y_max`, extreme bounded offsets and
   post-noise projection.
5. Confirm the measured deterministic error never exceeds `rho`, optimizer
   error never exceeds its certificate and final projection never increases
   distance to an in-box deterministic target.

### Transcript and security

1. Compare phase count, frame count, public sizes, circuit shapes and receipt
   order for all-missing, all-zero, random, separated, collinear and every
   boundary-valued dataset with identical public schema.
2. Run two real worker processes for `K=2,3,4,5`, including computation peers
   that are non-adjacent in source order. Reorder delivery, duplicate frames,
   drop acknowledgements, reconnect and restart at every durable phase.
3. Tamper with capsule/model IDs, pinsets, purposes, numeric plans, source
   digests, offsets and terminal signatures. No partial result, private error
   detail or alternative backend may be returned.
4. Demonstrate that the relay cannot decode any source/noise payload and that
   either single additive share is uniform over its declared ring. Document
   two-peer collusion and malicious deviation as out of scope; do not simulate
   them away in marketing claims.
5. Assert that no GLM trace calls legacy DCF comparison, local truncation,
   analyst-selected slot binding, exact score/Hessian/deviance or plaintext
   finalisation endpoints.
6. Fuzz all parsers and canonical encoders; run race, crash-consistency, spool
   quota/backpressure and secret-file permission tests.

### Statistical validity

1. Compare four explicitly labelled targets on deterministic fixtures:
   ordinary high-precision GLM when it exists, bounded/clipped regularised
   exact-link target, PWL finite-iteration target and final DP release.
2. Cover rare binomial events, complete/quasi separation, collinearity, sparse
   and high Poisson counts, clipping, heavy missingness, offsets and weights.
3. Verify each deterministic difference is either within a published bound or
   labelled an estimand change with no universal bound.
4. In Monte Carlo, verify mechanism-only marginal/simultaneous radii against
   the implemented sampler. Do not call those intervals population confidence
   intervals.
5. Demonstrate coefficient scale and factor/reference coding against trusted R
   and independent high-precision references. No-noise parity is to the signed
   bounded surrogate first; closeness to ordinary `glm` is a reported empirical
   comparison, not a theorem.
6. Have a biomedical statistician review offset, weight, repeated-record,
   missingness and interpretability documentation before promotion.

### DSI and scaling

1. Run direct named fan-out and central aggregate paths over two real DSLite
   processes, then real Opal and Armadillo deployments with their normal login
   configuration. Verify there is no dsVert-specific Armadillo TLS-attestation
   option and no positional routing assumption.
2. Transfer large typed source blocks under bounded spool smaller than the
   payload, with backpressure, replay and reconnect; final hashes and source
   commitments must match byte for byte.
3. Benchmark a factorial grid over `C`, `p`, `J`, `T`, ring width and
   `K=2..5+`. Publish wall time, CPU, peak RSS, wire bytes, DSI calls, peer
   bytes and circuit-cache/OT setup costs with hardware and software versions.
4. Verify empirical source work follows total `Cp` plus `K` control overhead,
   and optimizer work follows `TC(p+log J)` within measured cryptographic
   constants. Investigate any unexplained super-linear regime before release.
5. Prove peak working memory is bounded by the public row block and worker
   concurrency; no implementation materialises `C*p*T`, a dense Hessian or a
   second relay-side payload copy.
6. Verify one DSI source upload per capsule and no analyst/DSI round trip per
   row, coefficient or optimization iteration.

### API and non-regression

1. Test canonical formulas, factor contrasts, intercept inclusion/exclusion,
   registered bounded offsets/weights and every explicitly unsupported term.
2. Test unlimited repeated client operations over the same capsule: output
   bytes and accuracy remain unchanged and no query counter, epsilon remaining
   or history gate appears.
3. Test noise-root regeneration as a visible new epoch and peer-key rotation
   as an unrecognised-peer error with administrative re-pinning instructions;
   neither is represented as privacy-budget exhaustion.
4. Run complete server/client unit, integration, multiprocess and package-check
   suites. Promotion requires zero skipped tests for available core security
   dependencies and an explicit list of environment-only deployment smokes.

## Primary sources

- Dwork, McSherry, Nissim and Smith, [*Calibrating Noise to Sensitivity in
  Private Data Analysis*](https://www.iacr.org/archive/tcc2006/38760266/38760266.pdf),
  TCC 2006. Primary sensitivity-calibrated DP mechanism and transcript
  formulation.
- Chaudhuri, Monteleoni and Sarwate, [*Differentially Private Empirical Risk
  Minimization*](https://jmlr.org/papers/v12/chaudhuri11a.html), JMLR 2011.
  Primary output- and objective-perturbation analysis for regularised ERM,
  including explicit convexity, differentiability and sensitivity assumptions.
- Abadi et al., [*Deep Learning with Differential
  Privacy*](https://research.google/pubs/deep-learning-with-differential-privacy/),
  CCS 2016. Primary DP-SGD clipping/noise/accounting route.
- Song, Steinke, Thakkar and Thakurta, [*Evading the Curse of Dimensionality in
  Unconstrained Private GLMs*](https://proceedings.mlr.press/v130/song21a.html),
  AISTATS 2021. Primary GLM-specific DP gradient analysis, including the
  modified objective induced by gradient clipping.
- Arora, Bassily, Guzmán, Menart and Ullah, [*Differentially Private
  Generalized Linear Models
  Revisited*](https://proceedings.neurips.cc/paper_files/paper/2022/hash/8d321ebb82b58987509b8624cbb85d65-Abstract-Conference.html),
  NeurIPS 2022. Primary modern analysis of private convex linear predictors and
  constrained regularised output perturbation.
- Redberg, Koskela and Wang, [*Improving the Privacy and Practicality of
  Objective Perturbation for Differentially Private Linear
  Learners*](https://papers.nips.cc/paper_files/paper/2023/hash/2ceda49041816da6d5a34eb3b612607f-Abstract-Conference.html),
  NeurIPS 2023. Primary modern objective-perturbation analysis and approximate
  minima treatment.
- Canonne, Kamath and Steinke, [*The Discrete Gaussian for Differential
  Privacy*](https://proceedings.neurips.cc/paper/2020/hash/b53b3a3d6ab90ce0268229151c9bde11-Abstract.html),
  NeurIPS 2020. Primary discrete-Gaussian privacy and exact sampling analysis.
- Zhang et al., [*Functional Mechanism: Regression Analysis under Differential
  Privacy*](https://www.vldb.org/pvldb/vol5/p1364_junzhang_vldb2012.pdf), PVLDB
  2012. Primary polynomial objective-perturbation proposal and discussion of
  invalid/noisy objective repair.
- Kulkarni, Jälkö, Koskela, Kaski and Honkela, [*Differentially Private
  Bayesian Inference for Generalized Linear
  Models*](https://proceedings.mlr.press/v139/kulkarni21a.html), ICML 2021.
  Primary logistic/Poisson inference from noisy approximate summary
  statistics; it is evidence for a distinct noise-aware Bayesian route, not
  for treating finite moments as exact GLM sufficient statistics.
- Ferrando and Sheldon, [*Private Regression via Data-Dependent Sufficient
  Statistic Perturbation*](https://openreview.net/forum?id=gtCfDKm9ME), TMLR
  2025. Primary current evidence that private marginals can support a
  competitive approximate logistic objective; it remains a separately
  specified approximation and release workload.
- Patra, Schneider, Suresh and Yalame, [*ABY2.0: Improved Mixed-Protocol
  Secure Two-Party Computation*](https://www.usenix.org/conference/usenixsecurity21/presentation/patra),
  USENIX Security 2021. Primary evidence that keeping dot products and matrix
  operations in a ring protocol and converting only comparison/nonlinear
  subgraphs can materially reduce online communication relative to GC-only
  training; adopting such a kernel would require a new reviewed backend, not
  a relaxation of dsVert's exact truncation contract.
- Kelkar, Le, Raykova and Seth, [*Secure Poisson
  Regression*](https://www.usenix.org/conference/usenixsecurity22/presentation/kelkar),
  USENIX Security 2022. Primary batched fixed-point exponentiation and
  correlated-matrix-multiplication design for secure Poisson regression. It
  motivates the post-E2E vector kernel and benchmark gate; its assumptions do
  not by themselves certify dsVert's DSI relay, sticky DP release or ring
  bounds.
