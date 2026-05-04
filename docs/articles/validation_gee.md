# GEE validation

## What is validated

Functions:
[`ds.vertGEE()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGEE.md)

GEE solves $`\sum_i D_i^T V_i^{-1}(y_i-\mu_i)=0`$ and reports the
sandwich variance $`B^{-1}MB^{-1}`$.
[`ds.vertGEE()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGEE.md)
supports independence, exchangeable and guarded AR1 working correlations
for Gaussian, binomial and Poisson outcomes. Cluster and adjacent-pair
scores are computed as aggregate share-domain sufficient statistics.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fit_ref <- geepack::geeglm(
  y ~ x1 + x2 + x3, id = cluster_id, waves = visit,
  data = pooled, family = binomial(), corstr = "ar1")
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "cluster_id", "visit", "y")]
)
```

``` r


fit_gee <- dsVertClient::ds.vertGEE(
  y ~ x1 + x2 + x3,
  data = "DA",
  family = "binomial",
  id_col = "cluster_id",
  order_col = "visit",
  corstr = "ar1",
  binomial_sigmoid_intervals = 150L,
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_gee.R
```

## Disclosure review

Residuals, fitted means, probabilities, Pearson residuals, row scores,
visit labels and adjacent-pair vectors are not returned. AR1 order
metadata is encrypted server-to-server and only guarded low-dimensional
sufficient statistics reach the client.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("gee")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertGEE | synthetic clustered longitudinal | geepack/geepack-compatible protected oracle | worst_coeff_or_se_abs | 0.0033 | 0.005 | strict-practical | PASS | gee_dslite_20260504-160614.rds |
| K\>=3 | ds.vertGEE | synthetic clustered longitudinal | geepack/geepack-compatible protected oracle | worst_coeff_or_se_abs | 0.0033 | 0.005 | strict-practical | PASS | gee_dslite_20260504-161741.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
