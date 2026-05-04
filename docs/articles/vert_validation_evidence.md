# dsVert validation evidence

## Scope

This page is the validation index for the current dsVertClient product
surface. The historical vignettes under `vignettes/methods/` were
discarded because they mixed product methods with legacy or diagnostic
routes. The pages linked here document the current supported route for
each method family, with K=2 and K\>=3 evidence from deterministic
DSLite harnesses.

Each method vignette contains:

- the statistical target and mathematical estimating equation;
- the centralized reference construction;
- the vertical DSLite split and product function call;
- a disclosure review describing what the analyst and peer servers can
  see;
- an executed evidence check against the cached DSLite result table.

## Summary table

``` r

all_rows <- load_validation_summary()
assert_validation(all_rows)
knitr::kable(all_rows[, c("method_name", "k_mode", "function_route",
                          "primary_metric", "observed", "tolerance",
                          "tier", "status", "cache")])
```

| method_name | k_mode | function_route | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| PSI alignment | K=2 | ds.psiAlign | intersection_count_delta | 0.0000000 | 0e+00 | strict-precise | PASS | psi_dslite_20260503-003316.rds |
| PSI alignment | K\>=3 | ds.psiAlign | intersection_count_delta | 0.0000000 | 0e+00 | strict-precise | PASS | psi_dslite_20260503-003316.rds |
| Descriptive statistics | K=2 | ds.vertDesc | quantile_max_abs | 2.2327101 | 3e+00 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |
| Descriptive statistics | K\>=3 | ds.vertDesc | quantile_max_abs | 2.2327101 | 3e+00 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |
| Contingency tests | K=2 | ds.vertChisq / ds.vertFisher / ds.vertChisqCross | count_or_pvalue_delta | 0.0000000 | 0e+00 | strict-precise | PASS | tables_dslite_20260504-113158.rds |
| Contingency tests | K\>=3 | ds.vertChisq / ds.vertFisher / ds.vertChisqCross | count_or_pvalue_delta | 0.0000000 | 0e+00 | strict-precise | PASS | tables_dslite_20260504-113158.rds |
| Correlation | K=2 | ds.vertCor | correlation_max_abs | 0.0000078 | 1e-04 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |
| Correlation | K\>=3 | ds.vertCor | correlation_max_abs | 0.0000093 | 1e-04 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |
| PCA | K=2 | ds.vertPCA | loading_max_abs | 0.0000441 | 1e-04 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |
| PCA | K\>=3 | ds.vertPCA | loading_max_abs | 0.0000592 | 1e-04 | strict-practical | PASS | descriptive_dslite_20260502-224817.rds |
| GLM | K=2 | ds.vertGLM(binomial_sigmoid_intervals=150) | coef_max_abs | 0.0047327 | 1e-02 | strict-practical | PASS | glm_dslite_20260504-182112.rds |
| GLM | K\>=3 | ds.vertGLM(binomial_sigmoid_intervals=150) | coef_max_abs | 0.0048161 | 1e-02 | strict-practical | PASS | glm_dslite_20260504-182112.rds |
| Inference helpers | K=2 | ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR | algebra_delta | 0.0000000 | 0e+00 | strict-precise | PASS | glm_wrappers_dslite_20260504-091506.rds |
| Inference helpers | K\>=3 | ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR | algebra_delta | 0.0000000 | 0e+00 | strict-precise | PASS | glm_wrappers_dslite_20260504-091506.rds |
| LASSO | K=2 | ds.vertLASSOIter | coef_max_abs_vs_glmnet | 0.0014013 | 5e-03 | strict-practical | PASS | lasso_binomial_dslite_20260504-171300.rds |
| LASSO | K\>=3 | ds.vertLASSOIter | coef_max_abs_vs_glmnet | 0.0013311 | 5e-03 | strict-practical | PASS | lasso_binomial_dslite_20260504-173214.rds |
| Negative binomial | K=2 | ds.vertNBFullRegTheta(variant=‘full_reg_nd’) | max(beta_abs, theta_abs) | 0.0006369 | 1e-03 | strict-practical | PASS | nb_dslite_20260503-200529.rds |
| Negative binomial | K\>=3 | ds.vertNBFullRegTheta(variant=‘full_reg_nd’) | max(beta_abs, theta_abs) | 0.0004528 | 1e-03 | strict-practical | PASS | nb_dslite_20260503-201831.rds |
| Cox PH | K=2 | ds.vertCox / ds.vertCoxProfileNonDisclosive | slope_max_abs | 0.0000351 | 1e-04 | strict-precise | PASS | cox_dslite_k2_20260504-144510.rds |
| Cox PH | K\>=3 | ds.vertCox / ds.vertCoxProfileNonDisclosive | slope_max_abs | 0.0000351 | 1e-04 | strict-precise | PASS | cox_dslite_k3_20260504-144829.rds |
| LMM | K=2 | ds.vertLMM | max(fixed_effect_abs, variance_abs) | 0.0007200 | 1e-03 | strict-practical | PASS | lmm_dslite_20260503-030040.rds |
| LMM | K\>=3 | ds.vertLMM.k3 | max(fixed_effect_abs, variance_abs) | 0.0002815 | 1e-03 | strict-practical | PASS | lmm_dslite_20260503-214223.rds |
| GEE | K=2 | ds.vertGEE | worst_coeff_or_se_abs | 0.0033000 | 5e-03 | strict-practical | PASS | gee_dslite_20260504-160614.rds |
| GEE | K\>=3 | ds.vertGEE | worst_coeff_or_se_abs | 0.0033000 | 5e-03 | strict-practical | PASS | gee_dslite_20260504-161741.rds |
| GLMM | K=2 | ds.vertGLMM | coef_max_abs_vs_glmmPQL | 0.0006720 | 2e-03 | strict-practical | PASS | glmm_dslite_20260504-082705.rds |
| GLMM | K\>=3 | ds.vertGLMM | coef_max_abs_vs_glmmPQL | 0.0006720 | 2e-03 | strict-practical | PASS | glmm_dslite_20260504-082705.rds |
| IPW | K=2 | ds.vertIPW | weighted_outcome_coef_abs | 0.0021589 | 3e-03 | strict-practical | PASS | ipw_dslite_20260504-010636.rds |
| IPW | K\>=3 | ds.vertIPW | weighted_outcome_coef_abs | 0.0021774 | 3e-03 | strict-practical | PASS | ipw_dslite_20260504-010636.rds |
| Multiple imputation | K=2 | ds.vertMI | pooled_coef_abs | 0.0000111 | 1e-04 | strict-precise | PASS | mi_dslite_20260503-003321.rds |
| Multiple imputation | K\>=3 | ds.vertMI | pooled_coef_abs | 0.0000108 | 1e-04 | strict-precise | PASS | mi_dslite_20260503-003321.rds |
| Multinomial | K=2 | ds.vertMultinom / ds.vertMultinomJointNewton | class_probability_max_abs | 0.0007369 | 1e-03 | strict-precise | PASS | multinom_joint_dslite_20260502-214939.rds |
| Multinomial | K\>=3 | ds.vertMultinom / ds.vertMultinomJointNewton | class_probability_max_abs | 0.0007658 | 1e-03 | strict-precise | PASS | multinom_joint_dslite_20260503-042440.rds |
| Ordinal | K=2 | ds.vertOrdinal / ds.vertOrdinalJointNewton | cumulative_probability_max_abs | 0.0004281 | 1e-03 | strict-precise | PASS | ordinal_joint_dslite_k2_20260503-160920.rds |
| Ordinal | K\>=3 | ds.vertOrdinal / ds.vertOrdinalJointNewton | cumulative_probability_max_abs | 0.0005697 | 1e-03 | strict-precise | PASS | ordinal_joint_dslite_k3_20260503-155749.rds |

## Method pages

- [PSI
  alignment](https://isglobal-brge.github.io/dsVertClient/articles/validation_psi.md)
- [Descriptive
  statistics](https://isglobal-brge.github.io/dsVertClient/articles/validation_descriptive.md)
- [Contingency
  tests](https://isglobal-brge.github.io/dsVertClient/articles/validation_contingency.md)
- [Correlation](https://isglobal-brge.github.io/dsVertClient/articles/validation_correlation.md)
- [PCA](https://isglobal-brge.github.io/dsVertClient/articles/validation_pca.md)
- [GLM](https://isglobal-brge.github.io/dsVertClient/articles/validation_glm.md)
- [GLM inference
  helpers](https://isglobal-brge.github.io/dsVertClient/articles/validation_inference.md)
- [LASSO](https://isglobal-brge.github.io/dsVertClient/articles/validation_lasso.md)
- [Negative
  binomial](https://isglobal-brge.github.io/dsVertClient/articles/validation_negative_binomial.md)
- [Cox
  PH](https://isglobal-brge.github.io/dsVertClient/articles/validation_cox.md)
- [LMM](https://isglobal-brge.github.io/dsVertClient/articles/validation_lmm.md)
- [GEE](https://isglobal-brge.github.io/dsVertClient/articles/validation_gee.md)
- [GLMM](https://isglobal-brge.github.io/dsVertClient/articles/validation_glmm.md)
- [IPW](https://isglobal-brge.github.io/dsVertClient/articles/validation_ipw.md)
- [Multiple
  imputation](https://isglobal-brge.github.io/dsVertClient/articles/validation_mi.md)
- [Multinomial](https://isglobal-brge.github.io/dsVertClient/articles/validation_multinomial.md)
- [Ordinal](https://isglobal-brge.github.io/dsVertClient/articles/validation_ordinal.md)

## Removed legacy routes

The package should not expose routes that are disclosive or materially
less accurate than an available product route. This check is executed
during rendering.

``` r

library(dsVertClient)
legacy <- load_legacy_removed()
checks <- assert_legacy_routes_removed()
knitr::kable(merge(legacy, checks, by.x = "route", by.y = "route", all.x = TRUE))
```

| route | product_replacement | reason | expected_state | pass |
|:---|:---|:---|:---|:---|
| ds.vertCox.k3 | ds.vertCoxProfileNonDisclosive | Person-time/legacy K\>=3 Cox path is not the product Cox PH claim. | absent export | TRUE |
| ds.vertCox(method=…) | ds.vertCoxProfileNonDisclosive | Rank/permutation metadata route was removed from the user wrapper. | unused argument | TRUE |
| ds.vertGLMM(method=…) | ds.vertGLMM aggregate PQL default | Only aggregate PQL is user-facing; method selection to legacy EM is removed. | unused argument | TRUE |
| ds.vertMultinom(method=‘warm’) | ds.vertMultinomJointNewton | Warm OVR is an initializer only and does not target the softmax MLE. | unused argument | TRUE |
| ds.vertNBFullRegTheta(variant=‘full_reg’) | ds.vertNBFullRegTheta(variant=‘full_reg_nd’) | Plain non-label eta transport is disclosive; share-domain theta/beta route replaces it. | variant rejected | TRUE |
| ds.vertOrdinal(method=‘warm’) | ds.vertOrdinalJointNewton | Warm cumulative-binomial is an initializer only and does not target the joint PO MLE. | unused argument | TRUE |

## Interpretation

All rows in the current table are `PASS`, non-disclosive under the
project standard, and cover both K=2 and K\>=3. Some methods are exact
to their centralized target; others are strict-practical because the
remaining gap is a documented secure approximation floor, such as
histogram quantile bucketization, secure sigmoid spline resolution, PQL
versus GLMM ML, or penalized optimizer depth.
