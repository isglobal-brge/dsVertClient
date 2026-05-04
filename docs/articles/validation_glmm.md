# GLMM validation

## What is validated

Functions:
[`ds.vertGLMM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLMM.md).

The product route is aggregate PQL for a binomial random-intercept GLMM.
The compact executable vignette runs one protected PQL outer update and
checks that the fixed-effect symmetry target and variance trace are
finite.

## Mathematical target

A PQL step alternates protected binomial working-response updates with
aggregate weighted mixed-model normal equations. The validation compares
fixed effects to the central symmetric glm target and asserts that the
PQL trace is populated.

## Fixture and reference

Fixture: Paired balanced binomial cluster fixture with cluster sizes
above privacy thresholds.

Centralized reference: Central binomial glm fixed-effect symmetry plus a
one-step PQL trace check.

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
| K=2 | ds.vertGLMM(max_outer=1) | paired balanced binomial cluster fixture | central glm fixed-effect symmetry plus PQL trace check | coef_max_abs_delta_and_pql_trace | 0 | 1e-04 | strict-pql | PASS | 89.4 |
| K\>=3 | ds.vertGLMM(max_outer=1) | paired balanced binomial cluster fixture | central glm fixed-effect symmetry plus PQL trace check | coef_max_abs_delta_and_pql_trace | 0 | 1e-04 | strict-pql | PASS | 90.6 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
