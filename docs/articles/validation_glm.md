# GLM validation

## What is validated

Functions:
[`ds.vertGLM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md).

The distributed route evaluates score and Fisher/Hessian aggregates with
Ring MPC. This compact vignette uses a Gaussian GLM so the full
validation runs quickly while exercising the same vertical design
assembly.

## Mathematical target

For mean mu_i = g^({-1}(x_i)T beta), beta solves X^T(y - mu) = 0 with
family-specific weights. The Gaussian fixture is compared to lm().

## Fixture and reference

Fixture: MASS::Pima.tr fixture with diabetes coded 0/1 and predictors
split vertically.

Centralized reference: stats::lm() on the pooled fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

The analyst receives model-level coefficients, covariance/standard
errors, and scalar diagnostics. Eta, residuals, row scores, and weights
remain internal.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("glm", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertGLM | MASS::Pima.tr fixture | stats::lm | coef_max_abs_delta | 0.0000607 | 0.001 | strict-practical | PASS | 7.4 |
| K\>=3 | ds.vertGLM | MASS::Pima.tr fixture | stats::lm | coef_max_abs_delta | 0.0001753 | 0.001 | strict-practical | PASS | 8.4 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
