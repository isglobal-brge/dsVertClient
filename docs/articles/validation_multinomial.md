# Multinomial validation

## What is validated

Functions: `ds.vertMultinom(), ds.vertMultinomJointNewton()`.

The supported route is joint softmax Newton over Ring127 shares.
Historical one-vs-rest output is kept only as an internal warm start.

## Mathematical target

P(Y=c\|X)=exp(eta_c)/sum_l exp(eta_l). The validation compares predicted
class probabilities to nnet::multinom().

## Fixture and reference

Fixture: Balanced soft-signal synthetic three-class fixture with
vertical predictors.

Centralized reference: nnet::multinom predicted probabilities.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Softmax probabilities, residuals, and row scores remain Ring127 shares.
The client receives coefficients and scalar optimizer diagnostics.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("multinomial", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertMultinomJointNewton | balanced soft-signal synthetic 3-class fixture | nnet::multinom probabilities | class_probability_max_abs_delta | 0.008508 | 0.02 | strict-practical | PASS | 152.3 |
| K\>=3 | ds.vertMultinomJointNewton | balanced soft-signal synthetic 3-class fixture | nnet::multinom probabilities | class_probability_max_abs_delta | 0.008520 | 0.02 | strict-practical | PASS | 152.4 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
