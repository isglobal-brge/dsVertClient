# Descriptive statistics validation

## What is validated

Functions:
[`ds.vertDesc()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDesc.md)

Means and standard deviations are computed from local scalar moments
$`\sum x`$, $`\sum x^2`$, and $`n`$. Quantiles are approximated by
guarded histograms: after safe bucket counts are released, the requested
quantile is linearly interpolated inside the bucket containing $`p n`$.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


pooled <- MASS::Pima.tr[seq_len(120), ]
summary(pooled[c("age", "bmi", "glu", "bp")])
stats::quantile(pooled$age, probs = c(.25, .5, .75))
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "age", "bmi", "ped")],
  s2 = pooled[c("patient_id", "npreg", "glu", "bp", "skin")]
)
```

``` r


dsVertClient::ds.psiAlign("D", "patient_id", "DA",
                          datasources = conns, verbose = FALSE)
desc <- dsVertClient::ds.vertDesc(
  "DA",
  variables = list(s1 = c("age", "bmi", "ped"),
                   s2 = c("npreg", "glu", "bp", "skin")),
  probs = c(.25, .5, .75),
  exact_extrema = FALSE,
  datasources = conns)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_descriptive.R
```

## Disclosure review

Exact extrema are suppressed by default. Histogram buckets with positive
counts below `datashield.privacyLevel` fail closed or are coarsened
before release. The quantile error is therefore a controlled
disclosure/approximation tradeoff, not an MPC reconstruction error.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("descriptive")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertDesc | MASS::Pima.tr | central summary/quantile | quantile_max_abs | 2.233 | 3 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |
| K\>=3 | ds.vertDesc | MASS::Pima.tr | central summary/quantile | quantile_max_abs | 2.233 | 3 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
