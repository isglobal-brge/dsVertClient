# GLMM validation

## What is validated

Functions:
[`ds.vertGLMM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLMM.md)

The supported GLMM route is binomial random-intercept PQL. The model is
$`\mathrm{logit}\{P(y_{ij}=1|b_i)\}=x_{ij}^T\beta+b_i`$ with
$`b_i\sim N(0,\sigma_b^2)`$. PQL iterates working responses and weights,
but the product implementation keeps these quantities in Ring127 shares
and returns only fixed effects, scalar variance components and
diagnostics.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fit_ref <- MASS::glmmPQL(
  y_bin ~ x1 + x2 + x3,
  random = ~ 1 | cluster_id,
  family = binomial(),
  data = pooled,
  verbose = FALSE)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "cluster_id", "y_bin")]
)
```

``` r


fit_glmm <- dsVertClient::ds.vertGLMM(
  y_bin ~ x1 + x2 + x3,
  data = "DA",
  cluster_col = "cluster_id",
  ring = 127L,
  compute_se = FALSE,
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_glmm.R
```

## Disclosure review

The removed method-selection path is not user-invokable. Working
responses, weights, probabilities, residuals and row scores remain in
shares. The public object is audited for no length-n and no
cluster-length vectors.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("glmm")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertGLMM | synthetic_clustered_binomial | MASS::glmmPQL | coef_max_abs_vs_glmmPQL | 0.000672 | 0.002 | strict-practical | PASS | glmm_dslite_20260504-082705.rds |
| K\>=3 | ds.vertGLMM | synthetic_clustered_binomial | MASS::glmmPQL | coef_max_abs_vs_glmmPQL | 0.000672 | 0.002 | strict-practical | PASS | glmm_dslite_20260504-082705.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
