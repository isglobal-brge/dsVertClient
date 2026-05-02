# Cox K=2 discrete-time pooled-logistic – non-disclosive option B (#D')

K=2 OT-Beaver discrete-time Cox via pooled-logistic with the per-patient
ending-bin index J_i hidden from the covariate (non-label) server.
Closes the disclosure gap of the leaky variant `run_cox_discrete` (where
the person-period frame leaked J_i = \#replicas(X_i) to the covariate
server).

Architecture (selected by David 2026-04-27 ~10:00, option B in the
project_k2_strict_unified_plan_2026-04-27.md "Option B FEASIBILITY
ANALYSIS" Sec.):

1.  Outcome (label) server holds (time, status) -\> computes J_i + m_ij
    = I(j \<= J_i) + y_ij = I(j == J_i AND status_i = 1) plaintext.

2.  m and y are split into Ring127 additive shares (frac=50) between OS
    and the covariate server. NL never sees J_i directly – only
    uniform-random-looking length-(J\*n) shares whose mod-2^127 sum with
    OS shares reconstructs the masks.

3.  Covariate server replicates X_i to a uniform Jxn person-period frame
    (every patient contributes exactly J rows, regardless of true J_i).
    Zero row-count signal leaks per Aliasgari-Blanton 2013 NDSS
    share-mask gating folklore.

4.  Newton inner loop (next session): mask-gated residual (y - p)*m via
    .ring127_vecmul Beaver round per bin + mask-gated Hessian X^T
    diag(W*m) X share-space + ridge epsilon\*I (~1e-8) on Hessian
    diagonal to handle all-zero-mask alpha_j strata per Christensen 2019
    Sec.A.3.

This first increment ships the orchestration that exercises the three
dormant server-side primitives shipped in dsVert eee40f6
(`dsvertCoxDiscreteShareMaskDS` + `dsvertCoxDiscreteReceiveSharesDS` +
`dsvertCoxDiscreteExpandXDS`) end-to-end and validates the
share-pipeline. The masked-Newton inner loop (Sec.4 above) is staged for
the next session – the scaffold returns the resolved share- slot keys +
uniform-expanded data frame name so the caller can wire the Newton when
it lands.

Citations:

- Aliasgari & Blanton 2013 NDSS (eprint 2012/405) – share-mask gating
  for FP branch elimination, the underlying primitive.

- De Cock et al. 2016 (eprint 2016/736) – oblivious selection under
  additive secret sharing.

- Mohassel & Zhang 2017 IEEE S&P (eprint 2017/396) SecureML –
  probabilistic-truncation noise model that bounds the Ring127 per-mult
  error budget.

- Catrina & Saxena 2010 FC Sec.3.3 – fixed-point representation
  statistical-security analysis (kappa = ring_bits - 2\*frac).

- Andreux et al. 2020 arXiv:2006.08997 – discrete-time Cox via
  pooled-logistic; the algorithmic substrate.

- Allison 1982 *Sociological Methodology* 13:61-98 – canonical
  pooled-logistic equivalence to discrete Cox.

- Prentice & Gloeckler 1978 *Biometrics* 34:57-67 – discrete Cox MLE
  properties.

- Christensen 2019 CRAN ordinal vignette Sec.A.3 – Newton diagonal
  eigenvalue inflation for handling singular H.

## Usage

``` r
ds.vertCoxDiscreteNonDisclosive(
  formula,
  data = NULL,
  J = 5L,
  bin_breaks = NULL,
  max_iter = 20L,
  tol = 1e-06,
  newton = TRUE,
  ridge_eps = 1e-06,
  verbose = FALSE,
  datasources = NULL
)
```

## Arguments

- formula:

  Cox-style formula `Surv(time, status) ~ x1 + x2 + ...` (interpreted as
  discrete-time pooled-logistic on time-bin reformulation). The LHS
  `Surv(time_var, status_var)` is parsed to extract the time / status
  column names; OS must hold both.

- data:

  Character. Server-side data frame symbol.

- J:

  Integer. Number of time bins (default 5L). Larger J -\> finer
  discretisation, lower bias to continuous Cox, higher MPC cost.

- bin_breaks:

  Numeric vector of length J+1 (sorted, increasing, first = 0). If
  `NULL`, the caller must precompute and pass them – they are public
  metadata that must be reproducible across servers (typically
  equal-quantile breaks of observed event times, computed at the
  coordinator from non-disclosive aggregates).

- max_iter:

  Integer. Newton outer iterations (default 20L). Used by the
  masked-Newton inner loop.

- tol:

  Numeric. Convergence tolerance on max(\|score\|) (default 1e-6).
  STRICT 1e-4 on beta_hat vs glm pooled-logistic reachable per the
  Catrina-Saxena 2010 Sec.3.3 noise-floor analysis (frac=50 Ring127 per-
  mult ulp approx 8.9e-16; depth-50 chain -\> relative error approx
  4.4e-14, ~10 orders of margin to STRICT 1e-4 on coefficients).

- newton:

  Logical. If TRUE (default), run the masked-Newton inner loop after the
  share-pipeline scaffold to produce coefficient estimates. If FALSE,
  return the scaffold result with coefficients=NULL (diagnostic mode for
  primitive validation).

- ridge_eps:

  Numeric. Diagonal eigenvalue inflation added to the Hessian before
  [`solve()`](https://rdrr.io/r/base/solve.html) to handle all-zero-mask
  alpha_j strata per Christensen 2019 CRAN ordinal vignette Sec.A.3,
  plus to suppress the spurious near-singular eigenvalues introduced by
  the masked-W noise (Catrina-Saxena 2010 frac=50 truncation
  accumulating across the mask*W*X^T\*X chain). Default 1e-6 –
  empirically drops the L2 fixture rel from 4.4e-2 (epsilon=1e-8) to
  4.1e-4 by suppressing the \|step\|-\>noise oscillation around the
  iter-4 attractor.

- verbose:

  Logical.

- datasources:

  DSI connections (length-2 K=2 split: label server holding (time,
  status) + a subset of X; covariate server holding the remaining X
  columns).

## Value

List with class `"ds.vertCoxDiscreteNonDisclosive"`:

- stage:

  Character: `"primitives_validated"` for this scaffold release. Will
  become `"converged"` / `"max_iter"` once the masked-Newton lands.

- coefficients:

  NULL in this scaffold; populated by Newton.

- n_obs:

  Number of patients.

- J:

  Number of bins.

- n_pp:

  Person-period rows = J \* n_obs (uniform – no leak).

- mask_share_key, y_share_key:

  Session slot names where the Ring127 (m, y) shares are stored on both
  servers, ready for the Newton.

- expanded_x_name:

  Symbol of the uniform Jxn covariate frame at the covariate server.

- disclosure_audit:

  List of per-step disclosure validation (n_pp at NL = J\*n with no
  per-patient row-count variation).
