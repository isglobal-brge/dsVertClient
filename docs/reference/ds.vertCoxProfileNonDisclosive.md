# Cox proportional hazards via non-disclosive Breslow profile Newton

Convenience wrapper for `ds.vertCoxDiscreteNonDisclosive` with
`target = "cox_profile"`. Event times, event indicators, and risk-set
membership remain local or in Ring127 share domain; the analyst receives
only the usual slope coefficients and aggregate Newton traces.

## Usage

``` r
ds.vertCoxProfileNonDisclosive(
  formula,
  data = NULL,
  max_event_times = NULL,
  max_iter = 20L,
  tol = 1e-06,
  newton = TRUE,
  ridge_eps = 1e-06,
  debug_trace = FALSE,
  verbose = FALSE,
  datasources = NULL
)
```

## Arguments

- formula:

  Cox-style `Surv(time, status) ~ x1 + ...` formula.

- data:

  Character server-side data symbol.

- max_event_times:

  Integer runtime guard for the number of distinct observed event times.

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
  share-domain setup to produce coefficient estimates. If FALSE, return
  primitive/setup audit metadata with coefficients=NULL (diagnostic mode
  for primitive validation).

- ridge_eps:

  Numeric. Diagonal eigenvalue inflation added to the Hessian before
  [`solve()`](https://rdrr.io/r/base/solve.html) to handle all-zero-mask
  alpha_j strata per Christensen 2019 CRAN ordinal vignette Sec.A.3,
  plus to suppress the spurious near-singular eigenvalues introduced by
  the masked-W noise (Catrina-Saxena 2010 frac=50 truncation
  accumulating across the `mask * W * X^T * X` chain). Default 1e-6 –
  empirically drops the L2 fixture rel from 4.4e-2 (epsilon=1e-8) to
  4.1e-4 by suppressing the \|step\|-\>noise oscillation around the
  iter-4 attractor.

- debug_trace:

  Logical. If `TRUE`, retain per-iteration beta and gradient traces plus
  diagnostic bin summaries for local debugging. Defaults to `FALSE`;
  setting it to `TRUE` requires
  `options(dsvert.allow_cox_debug_trace = TRUE)` or the
  `DSVERT_ALLOW_COX_DEBUG_TRACE=true` environment variable.

- verbose:

  Logical.

- datasources:

  DSI connections. The outcome server holds `time_var/status_var`. For
  K\>=3, the outcome server and one selected fusion server are the DCF
  parties; other servers contribute encrypted additive shares of their
  uniform person-period covariate frames.
