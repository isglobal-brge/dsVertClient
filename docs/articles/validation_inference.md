# GLM inference validation

## What is validated

Functions:
`ds.vertConfint(), ds.vertWald(), ds.vertContrast(), ds.vertLR()`.

Inference helpers post-process the beta, covariance, deviance, and
degrees of freedom already returned by ds.vertGLM.

## Mathematical target

Wald tests use z = (K beta - m) / sqrt(K V K^T). LR tests use
D_reduced - D_full against a chi-square reference.

## Fixture and reference

Fixture: The same MASS::Pima.tr Gaussian GLM fixture used for GLM
validation.

Centralized reference: Manual algebra on ds.vertGLM output.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

No new server aggregate is required beyond the already accepted
model-level GLM outputs.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("inference", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR | MASS::Pima.tr fixture | manual algebra on ds.vertGLM output | algebra_max_abs_delta | 0 | 0 | strict-precise | PASS | 15.3 |
| K\>=3 | ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR | MASS::Pima.tr fixture | manual algebra on ds.vertGLM output | algebra_max_abs_delta | 0 | 0 | strict-precise | PASS | 16.1 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
