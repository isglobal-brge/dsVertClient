# PCA validation

## What is validated

Functions:
[`ds.vertPCA()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertPCA.md).

PCA is computed from the validated correlation matrix, not from
patient-level scores.

## Mathematical target

The target is the eigendecomposition R = V Lambda V^T. The validation
compares eigenvalues and sign-aligned loadings.

## Fixture and reference

Fixture: MASS::Pima.tr numeric fixture split vertically.

Centralized reference: eigen(cor(X)) on the pooled fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Only eigenvalues and loadings derived from the aggregate correlation
matrix are returned; no individual PCA scores are released.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("pca", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertPCA | MASS::Pima.tr fixture | eigen(cor(X)) | max(eigen_loading_abs_delta) | 3.30e-05 | 1e-04 | strict-practical | PASS | 0 |
| K\>=3 | ds.vertPCA | MASS::Pima.tr fixture | eigen(cor(X)) | max(eigen_loading_abs_delta) | 3.09e-05 | 1e-04 | strict-practical | PASS | 0 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
