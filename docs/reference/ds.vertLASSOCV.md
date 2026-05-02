# Information-criterion lambda selection for one-step LASSO

Select the L1 penalty on the one-step quadratic-surrogate LASSO path
(`ds.vertLASSO1Step`) by minimising AIC / BIC / extended-BIC of the
quadratic-surrogate criterion. Because the surrogate uses only the
already-computed coefficient vector and the full covariance matrix,
selection is entirely client-side: no new MPC rounds are spent beyond
the single `ds.vertGLM` fit that produced `fit`.

For each lambda we solve the one-step LASSO, then score IC(lambda) =
surrogate_misfit(beta_lambda) + penalty \* df where `df` is the number
of nonzero coefficients and `penalty` is 2 (AIC), log(n) (BIC, default),
or log(n) + 2 gamma log(p) (extended BIC). The selected `lambda.min` is
the global minimiser; `lambda.1se` is reported as a more parsimonious
alternative that preserves at least `(1 - se_threshold) * IC_min` of the
fit.

This is the standard selector for one-step / SCAD-style penalised
maximum likelihood and is defensible in large samples without the need
for full K-fold refitting (which would require rerunning the MPC GLM
loop K times).

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

  A `ds.glm` object (with `fit$covariance`).

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

  For `lambda.1se`, retain the sparsest lambda whose IC is within this
  fraction of `IC_min` (default 0.02, i.e. 2 percent).

## Value

A `ds.vertLASSOCV` object: `lambda`, `ic`, `df`, `lambda.min`,
`lambda.1se`, `beta.min`, `beta.1se`, the original fit.
