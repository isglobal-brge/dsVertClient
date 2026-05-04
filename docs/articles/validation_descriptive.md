# Descriptive statistics validation

## What is validated

Functions:
[`ds.vertDesc()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDesc.md).

Each server computes local guarded moments for its own variables.
Quantiles are estimated from guarded histograms rather than releasing
sorted observations.

## Mathematical target

For a variable x, the exact checked moments are n, mean(x), and sd(x).
Histogram quantiles approximate F^{-1}(p) from binned counts under the
disclosure threshold.

## Fixture and reference

Fixture: MASS::Pima.tr numeric fixture split vertically.

Centralized reference: Central mean and standard deviation on the pooled
fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Only scalar summaries and guarded histogram information are returned; no
row-level values or order statistics are disclosed.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("descriptive", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertDesc | MASS::Pima.tr fixture | central mean/sd | max(mean_sd_abs_delta) | 0 | 0 | strict-precise | PASS | 0.6 |
| K\>=3 | ds.vertDesc | MASS::Pima.tr fixture | central mean/sd | max(mean_sd_abs_delta) | 0 | 0 | strict-precise | PASS | 0.7 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
