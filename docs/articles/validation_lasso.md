# LASSO validation

## What is validated

Functions:
[`ds.vertLASSO()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO.md),
[`ds.vertLASSO1Step()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO1Step.md),
[`ds.vertLASSOCV()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOCV.md),
[`ds.vertLASSOIter()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOIter.md),
[`ds.vertLASSOProximal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOProximal.md)

The LASSO target is
$`\arg\min_\beta n^{-1}\ell(\beta)+\lambda\|\beta_{-0}\|_1`$. Gaussian
routes solve the normal-equation LASSO from aggregate covariance/Hessian
objects. Binomial routes use secure aggregate-score proximal gradient
with FISTA restart and a guarded aggregate-Gram Lipschitz bound. Poisson
routes use aggregate score/Hessian probes in a damped proximal-Newton
update.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


X <- model.matrix(diabetes ~ age + bmi + ped + glu, data = pooled)
y <- pooled$diabetes
fit_ref <- glmnet::glmnet(X[, -1], y, family = "binomial",
                          lambda = 0.02, standardize = TRUE)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "age", "bmi")],
  s2 = pooled[c("patient_id", "ped", "glu", "diabetes")]
)
```

``` r


lasso <- dsVertClient::ds.vertLASSOIter(
  diabetes ~ age + bmi + ped + glu,
  data = "DA",
  family = "binomial",
  lambda = 0.02,
  exact_non_gaussian = TRUE,
  binomial_sigmoid_intervals = 200L,
  lipschitz = "auto",
  fista_restart = TRUE,
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_lasso_binomial.R --sigmoid-intervals 200
```

## Disclosure review

The LASSO routes open only quantities already accepted for GLM,
inference, and guarded correlation: beta/covariance/Hessian,
p-dimensional aggregate scores, or p-by-p aggregate Hessians/Gram
bounds. No eta, probability, residual, objective-by-row, or row score
vector is returned.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("lasso")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertLASSOIter | MASS::Pima.tr | glmnet / centralized standardized solver | coef_max_abs_vs_glmnet | 0.001401 | 0.005 | strict-practical | PASS | lasso_binomial_dslite_20260504-171300.rds |
| K\>=3 | ds.vertLASSOIter | MASS::Pima.tr | glmnet / centralized standardized solver | coef_max_abs_vs_glmnet | 0.001331 | 0.005 | strict-practical | PASS | lasso_binomial_dslite_20260504-173214.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
