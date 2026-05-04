# LASSO validation

## What is validated

Functions:
[`ds.vertLASSOProximal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOProximal.md).

The proximal LASSO solver consumes the distributed GLM fit and optimizes
the penalized quadratic objective client-side.

## Mathematical target

The checked identity is the lambda = 0 limit: argmin 1/2 (beta -
beta_hat)^T H (beta - beta_hat) equals beta_hat.

## Fixture and reference

Fixture: MASS::Pima.tr Gaussian GLM fixture.

Centralized reference: The prime ds.vertGLM coefficients at lambda = 0.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

The proximal step uses only model-level beta/Hessian information and
returns coefficients; it does not request row-level data.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("lasso", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertLASSOProximal(lambda=0) | MASS::Pima.tr fixture | OLS limit of ds.vertGLM | lambda0_coef_abs_delta | 0 | 0 | strict-precise | PASS | 7.7 |
| K\>=3 | ds.vertLASSOProximal(lambda=0) | MASS::Pima.tr fixture | OLS limit of ds.vertGLM | lambda0_coef_abs_delta | 0 | 0 | strict-precise | PASS | 8.7 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
