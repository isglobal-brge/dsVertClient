# Correlation validation

## What is validated

Functions:
[`ds.vertCor()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md).

The method releases a low-dimensional correlation matrix from
server-local moments and MPC cross-products.

## Mathematical target

For variables x_j and x_l, r_jl = cov(x_j, x_l) / (sd(x_j) sd(x_l)). The
validation compares the full matrix to stats::cor().

## Fixture and reference

Fixture: MASS::Pima.tr numeric fixture split vertically.

Centralized reference: stats::cor() on the pooled fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

The disclosure surface is the guarded p by p correlation matrix, the
same aggregate tier used by downstream PCA.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("correlation", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertCor | MASS::Pima.tr fixture | stats::cor | correlation_max_abs_delta | 1.24e-05 | 1e-04 | strict-practical | PASS | 2.0 |
| K\>=3 | ds.vertCor | MASS::Pima.tr fixture | stats::cor | correlation_max_abs_delta | 1.24e-05 | 1e-04 | strict-practical | PASS | 2.8 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
