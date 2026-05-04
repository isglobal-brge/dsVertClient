# Ordinal validation

## What is validated

Functions:
[`ds.vertOrdinal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertOrdinal.md),
[`ds.vertOrdinalJointNewton()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertOrdinalJointNewton.md)

The proportional-odds model is $`P(Y\le k|x)=F(\theta_k-x^T\beta)`$,
with ordered thresholds $`\theta_1<\cdots<\theta_{K-1}`$. The product
route evaluates the joint score for beta and thresholds with class
masks, cumulative probabilities, reciprocals and finite-difference
Hessian probes held in Ring127 shares.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fit_ref <- MASS::polr(y_ord ~ x1 + x2 + x3,
                      data = pooled, method = "logistic", Hess = TRUE)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "low_leq", "mid_leq", "y_ord")]
)
```

``` r


fit_ord <- dsVertClient::ds.vertOrdinal(
  y_ord ~ x1 + x2 + x3,
  data = "DA",
  levels_ordered = c("low", "mid", "high"),
  cumulative_template = "%s_leq",
  max_outer = 2L,
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_ordinal_joint.R --K 2 && Rscript scripts/validate_method_ordinal_joint.R --K 3
```

## Disclosure review

Class masks, probabilities, weights and residual-like terms remain
encrypted shares. The client receives guarded class counts and aggregate
gradient/parameter summaries only. The warm cumulative-binomial route is
not a user-facing final estimator.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("ordinal")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertOrdinal / ds.vertOrdinalJointNewton | synthetic proportional-odds fixture | MASS::polr cumulative probabilities | cumulative_probability_max_abs | 0.0004281 | 0.001 | strict-precise | PASS | ordinal_joint_dslite_k2_20260503-160920.rds |
| K\>=3 | ds.vertOrdinal / ds.vertOrdinalJointNewton | synthetic proportional-odds fixture | MASS::polr cumulative probabilities | cumulative_probability_max_abs | 0.0005697 | 0.001 | strict-precise | PASS | ordinal_joint_dslite_k3_20260503-155749.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
