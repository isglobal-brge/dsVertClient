# Negative binomial validation

## What is validated

Functions:
[`ds.vertNB()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNB.md),
[`ds.vertNBMoMTheta()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNBMoMTheta.md),
[`ds.vertNBFullRegTheta()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNBFullRegTheta.md)

For NB2 log-link regression,
$`\mathrm{Var}(Y_i|X_i)=\mu_i+\mu_i^2/\theta`$. The full-regression
theta score contains terms such as $`\sum_i \log(\mu_i+\theta)`$ and
$`\sum_i (y_i+\theta)/(\mu_i+\theta)`$. The product route evaluates
those terms in Ring127 shares and alternates theta Newton updates with
aggregate Fisher beta updates.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fit_ref <- MASS::glm.nb(y ~ x1 + x2 + x3, data = pooled)
coef(fit_ref)
fit_ref$theta
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "y")]
)
```

``` r


fit_nb <- dsVertClient::ds.vertNBFullRegTheta(
  y ~ x1 + x2 + x3,
  data = "DA",
  variant = "full_reg_nd",
  theta_max_iter = 8L,
  beta_max_iter = 2L,
  compute_covariance = TRUE,
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_nb.R --variant full_reg_nd
```

## Disclosure review

The removed `variant='full_reg'` route disclosed non-label eta in
plaintext. The `full_reg_nd` route keeps eta, mu, logarithms,
reciprocals, score residuals, Fisher weights, and weighted covariates in
Ring127 additive shares. Only scalar theta sums, aggregate beta scores,
and aggregate Fisher matrices are reconstructed.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("negative_binomial")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertNBFullRegTheta(variant=‘full_reg_nd’) | Synthetic NB regression fixture | MASS::glm.nb | max(beta_abs, theta_abs) | 0.0006369 | 0.001 | strict-practical | PASS | nb_dslite_20260503-200529.rds |
| K\>=3 | ds.vertNBFullRegTheta(variant=‘full_reg_nd’) | Synthetic NB regression fixture | MASS::glm.nb | max(beta_abs, theta_abs) | 0.0004528 | 0.001 | strict-practical | PASS | nb_dslite_20260503-201831.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
