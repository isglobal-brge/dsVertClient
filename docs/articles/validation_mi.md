# Multiple imputation validation

## What is validated

Functions:
[`ds.vertMI()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMI.md).

Missing covariates are imputed server-side for each imputation round,
then GLM fits are pooled client-side with Rubin’s rules.

## Mathematical target

Rubin pooling uses beta_bar = mean_m beta_m and T = W + (1 + 1/M)B for
total variance.

## Fixture and reference

Fixture: Synthetic Gaussian regression fixture with deterministic
missingness in one covariate.

Centralized reference: Central mean-imputation regression used as a
compact reproducibility check.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Imputed columns stay server-side. The client sees only per-imputation
beta/covariance and pooled model summaries.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("mi", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertMI | synthetic missing-covariate fixture | central mean-imputation reference | pooled_coef_abs_delta | 0.009884 | 0.02 | strict-practical | PASS | 8.4 |
| K\>=3 | ds.vertMI | synthetic missing-covariate fixture | central mean-imputation reference | pooled_coef_abs_delta | 0.009884 | 0.02 | strict-practical | PASS | 9.6 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
