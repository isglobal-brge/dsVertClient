# PCA validation

## What is validated

Functions:
[`ds.vertPCA()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertPCA.md)

PCA on standardized variables is the eigendecomposition of the
correlation matrix $`R = V \Lambda V^T`$.
[`ds.vertPCA()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertPCA.md)
reuses
[`ds.vertCor()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md)
and performs the eigendecomposition client-side. The validated target is
therefore `eigen(cor(X))`.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


X <- MASS::Pima.tr[seq_len(120), c("age", "bmi", "ped", "glu", "bp")]
ref <- eigen(stats::cor(X), symmetric = TRUE)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = data.frame(patient_id = sprintf("P%03d", seq_len(nrow(X))),
                  X[c("age", "bmi", "ped")]),
  s2 = data.frame(patient_id = sprintf("P%03d", seq_len(nrow(X))),
                  X[c("glu", "bp")])
)
```

``` r


cor_fit <- dsVertClient::ds.vertCor("DA", variables = variables,
                                    datasources = conns)
pca_fit <- dsVertClient::ds.vertPCA(cor_result = cor_fit,
                                    n_components = 3)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_descriptive.R
```

## Disclosure review

Loadings and eigenvalues are aggregate outputs derived from the released
correlation matrix. Individual component scores are not returned because
scores would be row-level projections.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("pca")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertPCA | MASS::Pima.tr | eigen(cor(X)) | loading_max_abs | 4.41e-05 | 1e-04 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |
| K\>=3 | ds.vertPCA | MASS::Pima.tr | eigen(cor(X)) | loading_max_abs | 5.92e-05 | 1e-04 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
