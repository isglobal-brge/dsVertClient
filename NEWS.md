# dsVertClient (development version)

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
* **`ds.vertCox.k3`** — Cox K=3 via the Allison-1982 Poisson trick
  with offset + 100 baseline dummies; L-BFGS step-damping + backtracking
  guard (CGD).
* **`ds.vertCoxDiscreteNonDisclosive`** — K=2 discrete-time
  pooled-logistic Cox with Aliasgari-Blanton 2013 share-mask gating.
* **`ds.vertMultinom`** — K-class one-vs-rest multinomial logistic.
* **`ds.vertMultinomJoint`** — covariance-rescaled OVR (`full_irls`
  toggle).
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
  statistics; the older EM/Laplace-style route is diagnostics-only.
* **`ds.vertIPW`** — propensity-weighted two-stage GLM wrapper.
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
* **`ds.vertCor`** — Pearson correlation (Ring63 Beaver).
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
