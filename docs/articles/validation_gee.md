# GEE validation

## What is validated

Functions:
[`ds.vertGEE()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGEE.md).

The GEE route computes estimating-equation aggregates and robust
sandwich summaries for clustered data.

## Mathematical target

For independence working correlation, the estimating equation reduces to
sum_i D_i^T V_i^{-1}(y_i - mu_i)=0. The validation compares beta to
geepack.

## Fixture and reference

Fixture: MASS::Pima.tr clustered fixture with vertical predictors.

Centralized reference: geepack::geeglm(corstr=‘independence’) on the
pooled fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Returned values are regression and sandwich-level aggregates; row scores
are not returned.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("gee", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertGEE(corstr=‘independence’) | MASS::Pima.tr clustered fixture | geepack::geeglm | coef_max_abs_delta | 6.32e-05 | 0.01 | strict-practical | PASS | 12.2 |
| K\>=3 | ds.vertGEE(corstr=‘independence’) | MASS::Pima.tr clustered fixture | geepack::geeglm | coef_max_abs_delta | 6.16e-05 | 0.01 | strict-practical | PASS | 13.3 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
