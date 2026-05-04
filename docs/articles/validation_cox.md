# Cox PH validation

## What is validated

Functions:
[`ds.vertCox()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCox.md).

The current product route evaluates Cox partial-likelihood score terms
over discretised event-time risk sets without releasing risk-set rows.

## Mathematical target

The target solves sum_i delta_i (x_i - weighted_risk_mean(t_i; beta)) =
0. The fixture compares beta to coxph(ties=‘breslow’).

## Fixture and reference

Fixture: Synthetic survival fixture with event times discretised into
guarded bins.

Centralized reference: survival::coxph(…, ties = ‘breslow’) on the
pooled fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Risk-set contributions stay in the share domain; the client receives
coefficients and scalar convergence diagnostics.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("cox", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertCox | synthetic discretised survival fixture | survival::coxph(ties=‘breslow’) | coef_max_abs_delta | 0.0001888 | 0.001 | strict-practical | PASS | 88.6 |
| K\>=3 | ds.vertCox | synthetic discretised survival fixture | survival::coxph(ties=‘breslow’) | coef_max_abs_delta | 0.0001888 | 0.001 | strict-practical | PASS | 90.1 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
