# IPW validation

## What is validated

Functions:
[`ds.vertIPW()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertIPW.md).

IPW runs a protected propensity GLM and then a protected weighted
outcome GLM using an already server-side weight column.

## Mathematical target

The checked outcome target solves X^T W(y - X beta)=0 with
inverse-probability weights W.

## Fixture and reference

Fixture: Synthetic confounded treatment/outcome fixture with known IPW
column.

Centralized reference: Central weighted lm using the same weights.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Weights remain a server column; the analyst receives only propensity and
outcome model-level fits.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("ipw", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertIPW | synthetic confounded IPW fixture | central weighted lm using same weights | weighted_outcome_coef_abs_delta | 0.0001155 | 0.001 | strict-practical | PASS | 10.7 |
| K\>=3 | ds.vertIPW | synthetic confounded IPW fixture | central weighted lm using same weights | weighted_outcome_coef_abs_delta | 0.0001152 | 0.001 | strict-practical | PASS | 11.6 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
