# Ordinal validation

## What is validated

Functions: `ds.vertOrdinal(), ds.vertOrdinalJointNewton()`.

The supported route is joint proportional-odds Newton over Ring127
shares. Historical cumulative-binomial output is internal warm start
only.

## Mathematical target

P(Y \<= k \| X)=sigmoid(theta_k - X beta). The validation compares
cumulative probabilities to MASS::polr().

## Fixture and reference

Fixture: Balanced synthetic three-level ordered fixture with vertical
predictors.

Centralized reference: MASS::polr cumulative probabilities.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Cumulative probabilities, residuals, and row scores stay in shares; the
client receives thresholds, slopes, and scalar optimizer diagnostics.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("ordinal", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertOrdinalJointNewton | balanced synthetic 3-level ordinal fixture | MASS::polr cumulative probabilities | cumulative_probability_max_abs_delta | 0.02872 | 0.15 | strict-practical | PASS | 177.0 |
| K\>=3 | ds.vertOrdinalJointNewton | balanced synthetic 3-level ordinal fixture | MASS::polr cumulative probabilities | cumulative_probability_max_abs_delta | 0.02870 | 0.15 | strict-practical | PASS | 176.8 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
