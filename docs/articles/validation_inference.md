# GLM inference helper validation

## What is validated

Functions:
[`ds.vertConfint()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertConfint.md),
[`ds.vertWald()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertWald.md),
[`ds.vertContrast()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertContrast.md),
[`ds.vertLR()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLR.md)

The helpers are deterministic transformations of a `ds.glm` result: Wald
intervals use $`\hat\beta_j \pm z_{\alpha/2}\mathrm{SE}_j`$, univariate
Wald tests use $`(\hat\beta_j-\beta_{0j})/\mathrm{SE}_j`$, contrasts use
$`(K\hat\beta-m)^T(K\hat\Sigma K^T)^{-1}(K\hat\beta-m)`$, and LR tests
use the difference in returned deviances.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fit_full <- glm(glu ~ age + bmi + ped + bp + skin + npreg,
                data = MASS::Pima.tr[seq_len(80), ])
fit_reduced <- update(fit_full, . ~ . - skin)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "age", "bmi", "ped")],
  s2 = pooled[c("patient_id", "bp", "skin", "npreg", "glu")]
)
```

``` r


fit_full <- dsVertClient::ds.vertGLM(
  glu ~ age + bmi + ped + bp + skin + npreg,
  data = "DA", family = "gaussian", lambda = 0,
  datasources = conns, verbose = FALSE)
ci <- dsVertClient::ds.vertConfint(fit_full)
wald <- dsVertClient::ds.vertWald(fit_full, parm = "age")
contrast <- dsVertClient::ds.vertContrast(fit_full, K = "age = bmi")
lr <- dsVertClient::ds.vertLR(fit_reduced, fit_full)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_glm_wrappers.R
```

## Disclosure review

No server call is made by these helpers. They add no disclosure beyond
beta, covariance, and deviance already returned by the validated GLM
route.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("inference")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR | MASS::Pima.tr | manual algebra on ds.glm | algebra_delta | 0 | 0 | strict-precise | PASS | glm_wrappers_dslite_20260504-091506.rds |
| K\>=3 | ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR | MASS::Pima.tr | manual algebra on ds.glm | algebra_delta | 0 | 0 | strict-precise | PASS | glm_wrappers_dslite_20260504-091506.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
