# IPW validation

## What is validated

Functions:
[`ds.vertIPW()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertIPW.md)

IPW fits a propensity model $`e(x)=P(T=1|X=x)`$, then an outcome model
with weights $`w_i=T_i/e(x_i)+(1-T_i)/(1-e(x_i))`$. The product route
validates the two-stage propensity plus weighted GLM workflow while
keeping weights and sqrt-weights server-side or in additive shares.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


prop_ref <- glm(treat ~ x1 + x2, data = pooled, family = binomial())
pooled$ipw <- with(pooled, ifelse(treat == 1, 1 / fitted(prop_ref),
                                  1 / (1 - fitted(prop_ref))))
out_ref <- glm(y ~ treat + x1 + x2, data = pooled,
               family = binomial(), weights = ipw)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "treat", "ipw", "y")]
)
```

``` r


fit_ipw <- dsVertClient::ds.vertIPW(
  outcome_formula = y ~ treat + x1 + x2,
  propensity_formula = treat ~ x1 + x2,
  data = "DA",
  weights_column = "ipw",
  outcome_family = "binomial",
  datasources = conns,
  verbose = FALSE)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_ipw.R
```

## Disclosure review

Patient weights are not transported plaintext to DCF peers. Weighted
gradients/deviance use Beaver products of weight shares and residual
shares. The widened secure sigmoid fixed the previous stress accuracy
gap without changing disclosure.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("ipw")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertIPW | Synthetic confounded IPW fixture | central propensity GLM + weighted outcome GLM | weighted_outcome_coef_abs | 0.002159 | 0.003 | strict-practical | PASS | ipw_dslite_20260504-010636.rds |
| K\>=3 | ds.vertIPW | Synthetic confounded IPW fixture | central propensity GLM + weighted outcome GLM | weighted_outcome_coef_abs | 0.002177 | 0.003 | strict-practical | PASS | ipw_dslite_20260504-010636.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
