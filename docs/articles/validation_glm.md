# GLM validation

## What is validated

Functions:
[`ds.vertGLM()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md)

For a GLM with mean $`\mu_i=g^{-1}(x_i^T\beta)`$, the maximum-likelihood
estimate solves $`X^T(y-\mu)=0`$ (with the usual variance weights by
family). The distributed route evaluates aggregate scores,
Hessian/Fisher summaries, and deviance with Ring63/Ring127 MPC. The
validation table uses the binomial precision route with 150 secure
sigmoid intervals.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


pooled <- MASS::Pima.tr[seq_len(50), ]
pooled$diabetes <- as.integer(pooled$type == "Yes")
fit_ref <- glm(diabetes ~ age + bmi + ped + glu + bp + skin + npreg,
               data = pooled, family = binomial())
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "age", "bmi", "ped")],
  s2 = pooled[c("patient_id", "glu", "bp", "skin")],
  s3 = pooled[c("patient_id", "npreg", "diabetes")]
)
```

``` r


fit_ds <- dsVertClient::ds.vertGLM(
  diabetes ~ age + bmi + ped + glu + bp + skin + npreg,
  data = "DA",
  family = "binomial",
  lambda = 0,
  binomial_sigmoid_intervals = 150L,
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_glm.R --family binomial --binomial-sigmoid-intervals 150
```

## Disclosure review

The analyst receives model-scale coefficients, covariance/standard
errors, and scalar diagnostics. Per-patient eta, probabilities,
residuals, weights, and row scores stay in the share domain. Increasing
sigmoid intervals changes only numerical approximation, not the
disclosure surface.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("glm")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertGLM(binomial_sigmoid_intervals=150) | MASS::Pima.tr | stats::glm | coef_max_abs | 0.004733 | 0.01 | strict-practical | PASS | glm_dslite_20260504-182112.rds |
| K\>=3 | ds.vertGLM(binomial_sigmoid_intervals=150) | MASS::Pima.tr | stats::glm | coef_max_abs | 0.004816 | 0.01 | strict-practical | PASS | glm_dslite_20260504-182112.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
