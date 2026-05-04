# GLMM validation

## What is validated

Functions:
[`ds.vertGLMM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLMM.md).

The product route is aggregate PQL for a binomial random-intercept GLMM.
The compact executable vignette checks the no-outer-step smoke path
because full PQL is intentionally heavier.

## Mathematical target

With zero PQL outer updates, the route reduces to the protected binomial
GLM prime. The validation checks the coefficient envelope of that public
smoke path.

## Fixture and reference

Fixture: Balanced binomial cluster fixture with cluster sizes above
privacy thresholds.

Centralized reference: Central binomial glm no-random-effect limit.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

The route returns fixed effects and scalar variance diagnostics only;
per-cluster BLUPs and row probabilities are not returned.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("glmm", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertGLMM(smoke validation) | balanced binomial cluster fixture | central glm no-random-effect limit | coef_max_abs_delta | 0 | 0.05 | strict-smoke | PASS | 68.1 |
| K\>=3 | ds.vertGLMM(smoke validation) | balanced binomial cluster fixture | central glm no-random-effect limit | coef_max_abs_delta | 0 | 0.05 | strict-smoke | PASS | 69.0 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
