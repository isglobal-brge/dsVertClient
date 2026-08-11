# dsVert DP statistical-validity audit

Status: implementation audit, 2026-08-01.

This document defines the statistical claims supported by the current
`ds.vertDP*` surface. It is deliberately narrower than a claim that an output
is “95% accurate” without qualification. There are four different sources of
uncertainty or approximation:

1. DP mechanism noise, conditional on the fixed protected snapshot;
2. deterministic clipping, patient-level collapse, quantisation, and public
   grid discretisation;
3. sampling uncertainty for a target population or superpopulation; and
4. numerical approximation in a fitted model.

The base finite-snapshot accuracy certificates cover item 1.
`ds.vertDPDescribe()` also propagates the declared per-unit quantisation bound
for its mean and variance regions. Explicit methods whose names end in
`Inference` additionally combine the signed DP region with the sampling model
stated in their contract, using conservative exact Clopper--Pearson or Garwood
families; the corresponding base methods make no item-3 claim. The bounded
Gaussian route exposes a descriptive point fit and numerical identifiability
diagnostics, but deliberately does not claim nonlinear coefficient regions or
population-model inference for item 4.

## Defensible current surface

| Method | Finite-snapshot estimand | Contribution contract | Defensible output | Important non-claim |
|---|---|---|---|---|
| `ds.vertDPStatus()` | None | None | Read-only policy, accountant, pinset, epoch, and rollback-anchor handshake | It is not a statistical release or remote proof of an honest host |
| `ds.vertDPCalibrate()` | None | Public inputs only | Candidate mechanism radii and planning metrics | It does not select a universally correct epsilon or quantify sampling error |
| `ds.vertDPCount()` | Number of aligned privacy units in one canonical signed current-snapshot artifact; for replace-one, the signed fixed cohort size | One admitted unit under add/remove adjacency; zero sensitivity for the fixed public declaration | One bounded exact-MPC DP count with a conservative mechanism-noise radius, or an exact public K-consensus fixed-cohort count | No population-size sampling interval and no finite global composition claim across distinct Count artifacts |
| `ds.vertDPContingency()` | Fixed-domain cell counts on one server | At most one selected cell per unit | Non-negative DP table with marginal and simultaneous mechanism-noise radii | No ordinary Pearson, Fisher, or Yates p-value is valid merely by applying it to noisy cells |
| `ds.vertDPFrequency()` | Fixed-domain marginal counts and proportions after the signed repeated-record collapse | At most one public category per admitted unit; missing, out-of-domain and conflicting units are excluded | Non-negative DP counts plus simultaneous count and proportion mechanism regions from the same sticky capsule | Base regions exclude population sampling uncertainty and the excluded-unit mechanism is part of the estimand |
| `ds.vertDPFrequencyInference()` | Population category proportions under an iid multinomial privacy-unit model | No new contribution; client post-processing only | Conservative joint regions combining the signed simultaneous DP count box with Bonferroni exact Clopper--Pearson intervals | Exclusions must be scientifically ignorable; regions may be wide and no multinomial test or p-value is claimed |
| `ds.vertChisq()` | Independence in one fixed-domain same-owner categorical pair | No new contribution; client post-processing of one signed table | DP-aware multinomial plug-in bootstrap with Monte Carlo and sampler-TV calibration intervals | Asymptotic rather than finite-sample conditional calibration |
| `ds.vertFisher()` | Odds-ratio-one association null for one same-owner 2-by-2 pair | No new contribution; client post-processing of one signed table | DP-aware conditional hypergeometric plug-in bootstrap with Monte Carlo and sampler-TV calibration intervals | Not Fisher-exact for the confidential table; projected odds ratio is descriptive and no conditional CMLE/CI is returned |
| `ds.vertDPMeanVar()` | Mean and population central second moment of bounded patient-level values | Finite rows are clipped and averaged within unit; missing units are excluded | Bounded, internally consistent DP point estimates and simultaneous mechanism-noise regions obtained from the released noisy sufficient statistics | The regions exclude sampling uncertainty and can include a non-estimable boundary when the effective-count box reaches zero |
| `ds.vertDPDescribe()` | The same bounded patient-level moments plus fixed-grid empirical histograms | One bounded value per unit and variable | DP points; fixed-grid quantiles; simultaneous histogram quantile bands; simultaneous mean/variance regions that include mechanism noise and deterministic quantisation | Quantile bands include public bin width, but not sampling uncertainty; observed extrema are never estimated |
| `ds.vertDPGaussian()` | Bounded, patient-collapsed complete-case least-squares estimand from one signed Gaussian sufficient-statistic artifact | One clipped normalised monomial vector per admitted unit | Projected DP sufficient statistics, a typed point fit, residual second moment and numerical rank/conditioning diagnostics | No nonlinear coefficient mechanism region, classical standard error, p-value or population-sampling interval is claimed |
| `ds.vertCor()` | Bounded joint complete-case correlation matrix from the same signed Gaussian artifact | No new contribution; client post-processing of one sticky release | Raw complete-case DP correlations, outward mechanism/quantisation intervals and an explicitly labelled PSD projection | No pairwise fallback, hypothesis test or sampling confidence interval; the PSD point is not an exact reconstruction |
| `ds.vertPCA()` | Eigenstructure of the explicitly PSD-projected complete-case DP correlation matrix | No new contribution; client post-processing only | Eigenvalues/loadings, Weyl/Frobenius mechanism regions for eigenvalues and eigengap-based loading-stability diagnostics | No individual scores; loading signs are arbitrary and unresolved eigenspaces do not identify individual directions |
| `ds.vertDPSurvival()` | Fixed-grid event/censor/entry histogram after the declared patient tie-break | One exit cell and, for delayed entry, one entry cell per unit | Kaplan–Meier, Nelson–Aalen, Aalen–Johansen competing-risk curves, conservative simultaneous mechanism bands, and zero-cost fixed-grid survival-quantile, median, RMST and RMTL post-processing | The bands are not sampling confidence bands; beyond-horizon quantiles are not extrapolated, no p-value or proportional-hazards effect is inferred, and time-grid error remains separate and unquantified |
| `ds.vertDPEpi2x2()` | Risks, risk difference, risk ratio, and odds ratio of one released DP 2-by-2 table | No new contribution; client post-processing only | Per-estimand typed points and simultaneous mechanism-noise regions, including zero-denominator boundaries | Regions can be unbounded and are not epidemiological sampling confidence intervals |
| `ds.vertDPEpi2x2Inference()` | The same risks and derived effects under an iid joint exposure/outcome superpopulation model | No new contribution; client post-processing only | Conservative joint regions formed from a signed simultaneous DP count box and three Bonferroni exact Clopper-Pearson intervals | The sampling model must be scientifically appropriate; regions can be wide or unbounded and no hypothesis test or p-value is claimed |
| `ds.vertDPPrevalenceRatio()` | The identical finite-snapshot effects relabelled as exposed/unexposed/population prevalence, prevalence difference, prevalence ratio and prevalence odds ratio under a caller-declared cross-sectional interpretation | No new contribution; exact naming post-processing of `ds.vertDPEpi2x2()` | Numerically identical typed points, mechanism regions, attributable prevalence fractions and reciprocal prevalence-difference summary | The table cannot establish a cross-sectional design, temporality, representativeness or causality; mechanism regions exclude sampling uncertainty |
| `ds.vertDPPrevalenceRatioInference()` | The identical iid effects and regions relabelled under a caller-declared cross-sectional population model | No new contribution; exact naming post-processing of `ds.vertDPEpi2x2Inference()` | Numerically identical conservative joint DP-mechanism and exact-binomial sampling regions | The cross-sectional iid sampling model is an external scientific assumption; no design or causal claim is inferred from the table |
| `ds.vertDPDiagnostic2x2()` | Sensitivity, specificity, PPV/VPP, NPV/VPN, prevalence, accuracy, LR+, LR−, and diagnostic odds ratio from one disease-by-test DP 2-by-2 table | No new contribution; client post-processing only | Explicitly oriented typed points and simultaneous mechanism-noise regions with zero, infinity, and non-estimability flags | Rows must encode disease status and columns test result; regions exclude sampling uncertainty and may be unbounded |
| `ds.vertDPDiagnostic2x2Inference()` | The same diagnostic measures under an iid joint disease/test superpopulation model | No new contribution; client post-processing only | Conservative joint regions formed from a signed simultaneous DP count box and six Bonferroni exact Clopper-Pearson intervals | The sampling model must be scientifically appropriate; zero-denominator states remain explicit, regions may be unbounded, and no hypothesis test or p-value is claimed |
| `ds.vertDPDirectStandardization()` | Directly standardised risk from one released DP strata-by-binary-outcome table and public weights | No new contribution; client post-processing only | Point estimate and a simultaneous mechanism-noise region | Public weights define the target standard population; no sampling interval or causal interpretation is added |
| `ds.vertDPDirectStandardizationInference()` | The same directly standardised risk under conditionally binomial outcomes within public strata | No new contribution; client post-processing only | Conservative joint region formed from a signed simultaneous DP count box and Bonferroni exact Clopper-Pearson intervals for positive-weight strata | Confidential stratum sample sizes are conditioned on and standard weights are fixed; no causal effect or population-transportability claim is added |
| `ds.vertDPCausalStandardization()` | Saturated treated/control risks standardised over fixed public strata and target weights | No new contribution; client post-processing only | Treated/control risks, risk difference, risk ratio, odds ratio and number-needed regions from one simultaneous DP count box | A causal interpretation requires consistency, conditional exchangeability, positivity, no interference, correct public mapping and valid target weights; no propensity or outcome nuisance model is estimated |
| `ds.vertDPCausalStandardizationInference()` | The same public-stratum g-formula under conditionally binomial outcomes within each stratum-treatment arm | No new contribution; client post-processing only | Conservative joint regions formed from a signed simultaneous DP count box and Bonferroni exact Clopper-Pearson arm intervals | Identification assumptions remain external scientific requirements; intervals can be wide or unbounded and no p-value or doubly robust claim is made |
| `ds.vertDPIndirectStandardization()` | Observed-to-expected event ratio from one released DP strata-by-binary-outcome table and public expected rates | No new contribution; client post-processing only | Typed O/E estimate and an outward-bracketed simultaneous mechanism-noise region | Denominator-zero boxes can be non-estimable or unbounded; no Garwood or sampling confidence interval is implied |
| `ds.vertDPIndirectStandardizationInference()` | The same O/E ratio under a Poisson total-count model with fixed externally valid expected rates | No new contribution; client post-processing only | Conservative joint region formed from the signed simultaneous DP count box and an exact Garwood family over compatible integer tables | Expected-denominator boxes containing zero are vacuous `[0, Inf]`; no p-value, causal effect or transportability claim is added |

Every accuracy region is simultaneous only over the coordinates stated in its
release certificate. Laplace regions use the exact granular integer mechanism
confidence calculation and a union bound. Gaussian regions subtract the
published sampler total-variation bound from the tail allocation before the
same union-bound logic. Cellwise projection to non-negative integers does not
widen error relative to a non-negative true count.

The deployed Gaussian implementation slack is not an empirical tolerance. Its
server plan uses exact rational zCDP calibration, a certified finite-support
tail bound and an outward-enclosed dyadic CDF approximation. Tail-truncation
and sampler total variation are both allocated explicitly, transferred to DP
through the certified `(1 + exp(epsilon))` factor, and bound into the signed
manifest. Two pinned noise peers each generate one complete independent draw;
the privacy accounting is the maximum per-peer guarantee, while the reported
simultaneous 95% radius reserves the two-draw total-variation transfer before
applying its vector union bound. If the rational plan, finite support or
accuracy reserve is not representable, Gaussian is unavailable and the public
pre-data selector uses the certified Laplace route. The fixed-work transcript
has public geometry, but dsVert makes no host constant-time claim.

## Current vertical-analysis boundary

The present DP release surface is statistically useful but is not yet a full
cross-peer vertical inference suite:

- fixed-cohort Count uses one public signed K-consensus declaration; add/remove
  Count uses its independent canonical current-snapshot exact-MPC artifact and
  is not a coordinate of the capsule vector;
- contingency supports a declared same-owner table or the signed cross-owner
  exact-GC categorical artifact;
- mean/variance releases one bounded variable from one peer; and
- each Describe or Survival specification executes on one peer;
- correlation/PCA consumes one signed same-owner or cross-owner Gaussian
  artifact whose complete-case mask is shared across all variables.

Consequently, `ds.vertDPEpi2x2()` and `ds.vertDPDiagnostic2x2()` can analyse
two variables only after one authorised same-peer DP contingency release. They
do not infer a cross-owner table unless the workload contains the explicit
signed cross-owner categorical artifact. No combination of separate
univariate DP releases can recover a joint distribution. The implemented
cross-owner categorical and Gaussian artifacts therefore use purpose-bound
fixed-transcript exact GC before their one joint-DP opening; other cross-owner
families still require a dedicated signed MPC/DP artifact.

## Estimand details that must be reported

- The privacy unit is controlled by the server policy. With multiple records
  per patient, numeric descriptives average clipped finite rows within patient.
- A contingency table deterministically contributes at most one category pair
  per patient under public policy
  `consistent_cell_else_exclude_v1`. Repeated copies of the same valid cell
  count once. Two or more distinct valid cells make that patient contribute
  zero; missing and out-of-domain pairs are ignored and do not themselves
  create a conflict. No exact conflict/exclusion count is released. Analyses
  must report that the table targets concordant units and consider bias when
  disagreement across repeated records may be informative. The policy is
  custodian-owned, protocol-versioned, and bound into both policy and canonical
  query hashes.
- Diagnostic post-processing requires rows to be disease status and columns to
  be test result, with both positive levels selected explicitly. Sensitivity,
  specificity, likelihood ratios, and diagnostic odds ratio describe the
  finite released cohort; PPV/VPP, NPV/VPN, prevalence, and accuracy also
  inherit its case mix. All inherit the concordant-unit selection above, so
  neither population transportability nor a sampling model is implied.
- Indirect standardisation treats expected stratum rates as public and fixed.
  Its O/E region optimises the linear-fractional estimand over the continuous
  simultaneous count box, which conservatively covers every integer table;
  it does not reuse a Poisson/Garwood interval on DP-noised events.
- Mean and variance use custodian bounds and the population denominator `n`,
  not the `n - 1` sample-variance convention. The released point is a
  constrained post-processing of noisy quantised moments and can be biased by
  clipping, quantisation, noise, and consistency projection.
- Describe quantiles estimate a fixed-grid empirical quantile of valid bounded
  patient values. Their region covers every histogram in the simultaneous DP
  coordinate box and the public interval represented by a selected bin.
- Survival times are clipped to public limits and mapped to the first upper
  grid endpoint. A patient's first event wins; without an event, its latest
  censoring record wins. The point curves also apply a public consistency
  projection to noisy risk-flow counts and jointly rescale cause hazards if
  their sum exceeds one. `ds.vertDPRMST()` integrates the left-continuous
  product-limit step curve from the public lower bound to public `tau`.
  `ds.vertDPRMTL()` is exactly
  `(tau - time_lower_bound) - RMST`; the usual `tau - RMST` shorthand is valid
  only when the public lower bound is zero. Both are fixed-grid estimands, not
  assertions about unobserved within-bin event times.
  `ds.vertDPSurvivalQuantile()` reports the first public endpoint where the
  product-limit curve crosses the requested event-distribution probability;
  it inverts the simultaneous mechanism band and types a target as beyond the
  public horizon instead of extrapolating. `ds.vertDPMedianSurvival()` is its
  exact probability-one-half view.
- DP mechanism regions are probability statements over the mechanism for a
  fixed snapshot. They are not confidence intervals over repeated samples from
  a population and must not be labelled simply “95% CI”.

## Methods intentionally not inferred from current releases

The following would be statistically misleading and are therefore not added as
client post-processing:

- ordinary chi-square, Fisher, or Yates p-values from DP-noised tables (the
  separately labelled `ds.vertChisq()` and `ds.vertFisher()` routes instead
  reproduce the signed mechanism in explicit plug-in bootstrap laws);
- Wald or log-ratio sampling intervals from DP-noised 2-by-2 cells;
- Greenwood, log-rank, or standard Cox inference from the fixed-grid DP
  survival histogram without a mechanism-aware calibration;
- a separate landmark-analysis API: selecting an already released exact grid
  row is ordinary Kaplan--Meier subsetting, whereas a genuine post-landmark
  risk-set/cohort analysis changes the estimand and needs its own predeclared
  artifact and scientific assumptions;
- correlations or PCA reconstructed from separate univariate DP moments; and
- GLM, Cox, negative-binomial, multinomial, ordinal, mixed, GEE, or penalised
  model inference assembled from unrelated local releases.

## Concrete prerequisites and implemented boundaries

1. **DP contingency tests:** one purpose-bound table release and a
   mechanism-aware null calibration that includes nuisance estimation and the
   exact post-processing rule. Simulation must be driven only by public inputs
   and the already released DP object. Naive reference distributions remain
   forbidden.
2. **Correlation and PCA (implemented):** one signed `gaussian_models`
   release contains bounded count, Gram, outcome-cross and outcome-square
   coordinates under one complete-case mask, explicit joint L1/L2 sensitivity,
   PSD projection and clipping/quantisation regions. Cross-owner coordinates
   are produced by fixed-transcript exact GC before the joint noisy opening.
   Separate pairwise releases remain ineligible for PCA.
3. **GLM and related regression:** a purpose-bound joint MPC/DP mechanism for
   the final estimator or a formally accounted sequence of clipped gradients
   and Hessians. It needs public design/outcome bounds, contribution limits,
   convergence and identifiability states, and mechanism-aware covariance or
   an explicit point-estimate-only contract.
4. **Survival sampling inference:** the current simultaneous bands propagate
   the DP histogram mechanism only. Population confidence bands, contrasts or
   fitted survival models still need a predeclared target and a separately
   justified sampling model; mechanism bands must not be presented as
   sampling confidence bands.
5. **Missing data and repeated/clustered data:** a declared privacy unit and
   joint contribution bound before imputation or cluster aggregation. No
   patient- or cluster-granular intermediate may be opened to the analyst.

## Existing legacy surface that does not meet the DP reconstruction objective

No historical method was removed by this audit. The following distinction is
essential:

- Remaining legacy exact contingency variants, GLM, Cox, negative-binomial,
  ordinal, and non-capsule LASSO variants are
  outside the DP accountant. Even where a single invocation has a disclosure
  guard or uses MPC internally, adaptive exact outputs can be differenced and
  composed. They cannot support a general “source rows cannot be reconstructed”
  claim.
- `ds.vertCor()` and `ds.vertPCA()` now consume one signed complete-case
  Gaussian capsule; they no longer expose an exact or pairwise compatibility
  fallback. The exact slope-binomial iterative-LASSO and slope-bearing
  multinomial routes cannot reuse that artifact because their score MPC uses a
  different raw design. Those routes fail with typed signed-workload conditions
  before DSI until their purpose-bound design-Gram artifacts exist.
- `ds.vertDesc` is now a compatibility view over one explicit
  custodian-owned `ds.vertDPDescribe` artifact. It performs no exact column,
  moment, extrema or adaptive-histogram release; its legacy variable and
  quantile arguments are only local post-processing of the immutable capsule.
- LMM, clustered GEE, GLMM, and the current IPW contract are already
  quarantined. In addition to missing confirmatory statistical validation,
  some historical paths use cluster-granular intermediates and are especially
  unsuitable for a production non-reconstruction claim.
- Legacy client-only epidemiology, standardisation, and inference helpers make
  no server call, but inherit the disclosure and statistical validity of their
  input. An arbitrary caller-supplied matrix has no DataSHIELD or DP
  provenance.

These methods may remain for compatibility while the package is migrated, but
they must not be described as DP-safe or promoted into the single safe release
surface without the joint mechanisms above.

## Validation performed in this audit

- Reference tests for no-noise bounded moments, fixed-grid quantiles, KM,
  Nelson–Aalen, competing-risk cumulative incidence, fixed-grid RMST, and the
  exact RMTL complement identity (including nonzero public lower bounds).
- Deterministic finite-snapshot simulations of the real Gaussian solve,
  complete-case Cor projection/interval arithmetic, PCA eigenvalue regions and
  eigendirection stability. Gaussian coefficients and loading angles are
  recorded as point-only diagnostics because their API makes no region claim.
- Exhaustive enumeration of integer tables inside selected simultaneous-noise
  regression boxes, plus deterministic randomized property cases, for risk
  difference, risk ratio, odds ratio, direct and indirect standardisation, and
  all diagnostic-accuracy estimands. Diagnostic tests
  additionally exhaust every lower/upper cell box with endpoints from zero to
  two and verify exact flags for attainable zero, infinity, and
  non-estimability.
- Exhaustive small-grid and 250 deterministic property cases showing that the
  new Describe mean/variance regions contain the true bounded values whenever
  all count/sum/sumsq noises lie inside their simultaneous radii, including
  coordinate projection to zero and per-unit quantisation error.
- Separate add/remove and replace-one sensitivity checks for contingency,
  Describe histograms, and survival with and without delayed entry.
- Cross-package contract correction for the survival invalid-quality bin in
  both delayed-entry modes.
- High-confidence, million-coordinate Laplace radius regression avoiding
  floating-point cancellation in `1 - (1 - level) / cells`.

This evidence supports the finite-snapshot, mechanism-aware claims above. It
does not constitute clinical validation, a universal utility guarantee, or a
proof against malicious hosts or collusion outside the declared threat model.
