# Contingency tests validation

## What is validated

Functions: `ds.vertChisq(), ds.vertFisher(), ds.vertChisqCross()`.

Same-server categorical tests use guarded local contingency counts.
Cross-server tests build one-hot shares and aggregate cell counts with
MPC before applying chi-square/Fisher tests.

## Mathematical target

The validation compares O_ab counts and X^2 = sum_ab (O_ab - E_ab)^2 /
E_ab with centralized R.

## Fixture and reference

Fixture: Categorised MASS::Pima.tr variables with all margins above
disclosure thresholds.

Centralized reference: chisq.test() and fisher.test() on the pooled
fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

The released object contains guarded table counts and scalar test
statistics. Small cells fail closed.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("contingency", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertChisq / ds.vertFisher / ds.vertChisqCross | MASS::Pima.tr categorical fixture | chisq.test / fisher.test | max(count_or_chisq_delta) | 0 | 0 | strict-precise | PASS | 2.2 |
| K\>=3 | ds.vertChisq / ds.vertFisher / ds.vertChisqCross | MASS::Pima.tr categorical fixture | chisq.test / fisher.test | max(count_or_chisq_delta) | 0 | 0 | strict-precise | PASS | 2.0 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
