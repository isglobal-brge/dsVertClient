# Linear mixed model validation

## What is validated

Functions:
[`ds.vertLMM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLMM.md),
[`ds.vertLMM.k3()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLMM.k3.md)

The random-intercept LMM is $`y_{ij}=x_{ij}^T\beta+b_i+\epsilon_{ij}`$,
with $`b_i\sim N(0,\sigma_b^2)`$ and
$`\epsilon_{ij}\sim N(0,\sigma^2)`$. The estimator profiles variance
components and solves GLS normal equations. K=2 uses the closed-form
route; K\>=3 uses share-domain residual REML and exact cluster-mean GLS
transforms.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fit_ref <- nlme::lme(y ~ x1 + x2 + x3,
                     random = ~ 1 | cluster_id,
                     data = pooled, method = "REML")
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "cluster_id", "y")]
)
```

``` r


fit_lmm_k2 <- dsVertClient::ds.vertLMM(
  y ~ x1 + x2 + x3, data = "DA", cluster_col = "cluster_id",
  datasources = conns, verbose = FALSE)
fit_lmm_k3 <- dsVertClient::ds.vertLMM.k3(
  y ~ x1 + x2 + x3, data = "DA", cluster_col = "cluster_id",
  datasources = conns, verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_lmm.R
```

## Disclosure review

Cluster membership is method-level metadata and is sent only through
encrypted server-to-server channels. Original cluster labels are not
returned, small clusters fail closed, and patient-level residuals or
BLUPs are not returned.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("lmm")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertLMM | synthetic_balanced_random_intercept | nlme::lme / lme4::lmer | max(fixed_effect_abs, variance_abs) | 0.0007200 | 0.001 | strict-practical | PASS | lmm_dslite_20260503-030040.rds |
| K\>=3 | ds.vertLMM.k3 | synthetic_balanced_random_intercept | nlme::lme / lme4::lmer | max(fixed_effect_abs, variance_abs) | 0.0002815 | 0.001 | strict-practical | PASS | lmm_dslite_20260503-214223.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
