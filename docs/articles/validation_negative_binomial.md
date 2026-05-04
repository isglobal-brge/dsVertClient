# Negative binomial validation

## What is validated

Functions:
[`ds.vertNBMoMTheta()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNBMoMTheta.md).

The product route fits a negative-binomial mean model and estimates
overdispersion through guarded aggregate moments.

## Mathematical target

The model has Var(Y\|X) = mu + mu^2/theta and log(mu)=X beta. The
validation compares beta to MASS::glm.nb().

## Fixture and reference

Fixture: Synthetic overdispersed count regression fixture with vertical
predictors.

Centralized reference: MASS::glm.nb() on the pooled fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Returned values are beta and scalar theta/diagnostics. Per-patient
means, residuals, and score terms are not returned.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("negative_binomial", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertNBMoMTheta | synthetic NB fixture | MASS::glm.nb | coef_max_abs_delta | 0.005339 | 0.02 | strict-practical | PASS | 63.8 |
| K\>=3 | ds.vertNBMoMTheta | synthetic NB fixture | MASS::glm.nb | coef_max_abs_delta | 0.005339 | 0.02 | strict-practical | PASS | 64.1 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
