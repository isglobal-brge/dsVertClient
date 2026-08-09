# Client-side information-criterion selection for Gaussian LASSO

Select an L1 penalty entirely client-side. With a `ds.vertDPGaussian`
input, each candidate is solved from the same validated projected DP
moments and ranked by an explicitly labelled DP projected
pseudo-AIC/BIC/EBIC. With a historical `ds.glm` input, the existing
quadratic-surrogate selector is retained. Despite the historical
function/class suffix `CV`, neither route performs cross-validation or
resampling.

For a DP source, the score uses the projected noisy residual moment and
noisy effective count. It is therefore a deterministic model-selection
heuristic, not a classical likelihood information criterion; selection
is returned as unavailable when those released quantities cannot define
it. For a legacy source, the score is the quadratic-surrogate misfit
plus the requested degrees-of-freedom penalty. The preferred metadata
name for the more parsimonious alternative is `lambda.parsimonious`; the
historical `lambda.1se` slot is retained as an exact alias, but it is
not a one-standard-error rule. It uses a relative IC tolerance set by
`se_threshold` and involves no estimated sampling standard error.

Candidate paths are reusable post-processing and never trigger K-fold
refitting or repeated private releases.

## Usage

``` r
ds.vertLASSOCV(
  fit,
  lambda_grid = NULL,
  criterion = c("BIC", "AIC", "EBIC"),
  ebic_gamma = 0.5,
  keep_intercept = TRUE,
  se_threshold = 0.02
)
```

## Arguments

- fit:

  Preferably a `ds.vertDPGaussian` object. A historical `ds.glm` object
  with `fit$covariance` remains accepted.

- lambda_grid:

  Numeric vector of candidate lambda values (default: a 50-point
  log-spaced grid from `lambda_max` to `lambda_max / 1000`).

- criterion:

  One of `"BIC"` (default), `"AIC"`, or `"EBIC"`.

- ebic_gamma:

  Extended-BIC gamma parameter (default 0.5; effective only when
  `criterion = "EBIC"`).

- keep_intercept:

  Never penalise the intercept.

- se_threshold:

  Relative IC tolerance for the parsimonious selection: retain the
  sparsest lambda whose IC is no more than `abs(IC_min) * se_threshold`
  above `IC_min`. The historical argument name is retained for
  compatibility; it is not a standard error.

## Value

A `ds.vertLASSOCV` object: `lambda`, `ic`, `df`, `lambda.min`,
`lambda.parsimonious`, compatibility aliases `lambda.1se`/`beta.1se`,
and metadata explicitly identifying information-criterion selection
without cross-validation.
