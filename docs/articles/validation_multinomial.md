# Multinomial validation

## What is validated

Functions:
[`ds.vertMultinom()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinom.md),
[`ds.vertMultinomJointNewton()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinomJointNewton.md)

The softmax model uses $`P(Y=c|x)=\exp(\eta_c)/\sum_h\exp(\eta_h)`$,
with one reference class. The old one-vs-rest route is only an internal
warm start. The product route optimizes the joint softmax objective with
Ring127 share-domain probabilities, residuals, aggregate gradients and a
Bohning/aggregate-Gram Hessian surrogate.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fit_ref <- nnet::multinom(class ~ age + lwt + smoke + ht,
                          data = pooled, trace = FALSE)
predict(fit_ref, type = "probs")
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "age", "lwt")],
  s2 = pooled[c("patient_id", "smoke", "ht", "class_ind_low", "class_ind_mid")]
)
```

``` r


fit_mn <- dsVertClient::ds.vertMultinom(
  class ~ age + lwt + smoke + ht,
  data = "DA",
  classes = c("low", "mid"),
  reference = "high",
  indicator_template = "class_ind_%s",
  max_outer = 30L,
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_multinom_joint.R --K 2 && Rscript scripts/validate_method_multinom_joint.R --K 3
```

## Disclosure review

Softmax probabilities and residuals remain Ring127 shares. The client
sees only low-dimensional aggregate gradients/Gram information and
coefficients. Warm OVR final estimation is not exported and is excluded
from the paper path.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("multinomial")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertMultinom / ds.vertMultinomJointNewton | MASS::birthwt tertiles | nnet::multinom | class_probability_max_abs | 0.0007369 | 0.001 | strict-precise | PASS | multinom_joint_dslite_20260502-214939.rds |
| K\>=3 | ds.vertMultinom / ds.vertMultinomJointNewton | MASS::birthwt tertiles | nnet::multinom | class_probability_max_abs | 0.0007658 | 0.001 | strict-precise | PASS | multinom_joint_dslite_20260503-042440.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
