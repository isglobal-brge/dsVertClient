# Cox PH validation

## What is validated

Functions:
[`ds.vertCox()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCox.md),
[`ds.vertCoxProfileNonDisclosive()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCoxProfileNonDisclosive.md),
[`ds.vertCoxDiscreteNonDisclosive()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCoxDiscreteNonDisclosive.md)

The Cox PH model has hazard $`h(t|x)=h_0(t)\exp(x^T\beta)`$. With
Breslow ties, the profile partial-likelihood score is built from
event-time risk-set sums $`\sum_{i:t_i \ge t_j}\exp(x_i^T\beta)x_i`$.
The product route hides risk-set and event masks as Ring127 shares and
opens only slope score/Hessian aggregates.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


fit_ref <- survival::coxph(
  survival::Surv(time, event) ~ x1 + x2 + x3,
  data = pooled, ties = "breslow")
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "time", "event")]
)
```

``` r


fit_cox <- dsVertClient::ds.vertCox(
  survival::Surv(time, event) ~ x1 + x2 + x3,
  data = "DA",
  max_iter = 5L,
  max_event_times = 50L,
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_cox.R --target cox_profile
```

## Disclosure review

Legacy rank/permutation and person-time Cox routes are not offered.
Event times, event indicators, event ranks, risk-set membership,
baseline dummies, and per-person period rows are not returned. Debug
traces and bin summaries are gated diagnostics only.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("cox")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertCox / ds.vertCoxProfileNonDisclosive | synthetic_event_time_survival | survival::coxph(ties=‘breslow’) | slope_max_abs | 3.51e-05 | 1e-04 | strict-precise | PASS | cox_dslite_k2_20260504-144510.rds |
| K\>=3 | ds.vertCox / ds.vertCoxProfileNonDisclosive | synthetic_event_time_survival | survival::coxph(ties=‘breslow’) | slope_max_abs | 3.51e-05 | 1e-04 | strict-precise | PASS | cox_dslite_k3_20260504-144829.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
