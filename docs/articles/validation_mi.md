# Multiple imputation validation

## What is validated

Functions:
[`ds.vertMI()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMI.md)

For $`m=1,\ldots,M`$, missing values are imputed server-side, the
analysis GLM returns $`\hat\beta_m`$ and $`\hat\Sigma_m`$, and Rubin’s
rules compute $`\bar\beta=M^{-1}\sum_m\hat\beta_m`$ and
$`T=W+(1+M^{-1})B`$, where $`W`$ is the mean within-imputation
covariance and $`B`$ is the between-imputation covariance.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fits <- lapply(seq_len(3), function(m) {
  imputed <- same_server_imputer(pooled, seed = 100 + m)
  glm(y ~ x1 + x2 + x3, data = imputed)
})
rubin_pool(fits)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "x1")],
  s2 = pooled[c("patient_id", "x2", "x3", "y")]
)
```

``` r


fit_mi <- dsVertClient::ds.vertMI(
  y ~ x1 + x2 + x3,
  data = "DA",
  impute_columns = "x1",
  m = 3L,
  family = "gaussian",
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_mi.R
```

## Disclosure review

Imputed values never leave the server that owns the column. The client
receives only per-imputation model-scale coefficient/covariance
summaries and Rubin pooled summaries. Exact extrema in imputation
diagnostics are suppressed.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("mi")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertMI | synthetic_gaussian_missing_x1 | same-imputer centralized Rubin pooling | pooled_coef_abs | 1.11e-05 | 1e-04 | strict-precise | PASS | mi_dslite_20260503-003321.rds |
| K\>=3 | ds.vertMI | synthetic_gaussian_missing_x1 | same-imputer centralized Rubin pooling | pooled_coef_abs | 1.08e-05 | 1e-04 | strict-precise | PASS | mi_dslite_20260503-003321.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
