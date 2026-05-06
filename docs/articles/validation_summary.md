# Validation evidence summary

## Purpose

This article is an index and audit summary for the executable validation
vignettes. It does not run a DSLite analysis itself. The method-specific
articles linked below are the executable evidence: each one generates
its dataset in the vignette, partitions it vertically, starts an
in-memory DSLite deployment, aligns records, runs the dsVert route, runs
the centralized R reference on the same pooled data, computes the
numerical delta, and asserts the result.

The current evidence snapshot covers 17 method blocks, both `K=2` and
`K>=3` modalities, and 38 route-level validation rows. All rendered rows
are `PASS`.

## Validation contract

Each method vignette builds a `rows` data frame with one row per
validated product route. The common columns are:

| Column | Meaning |
|----|----|
| `k_mode` | Whether the route is validated with two peers (`K=2`) or at least three peers (`K>=3`). |
| `function_route` | User-facing `ds.vert.*` route or route option being exercised. |
| `dataset` | The generated fixture used for the comparison. |
| `reference_target` | The centralized R reference used as the comparator. |
| `primary_metric` | The numerical error criterion, usually a max absolute coefficient, probability, loading, or algebra delta. |
| `observed` | The observed distributed-vs-centralized delta. |
| `tolerance` | The accepted validation envelope for the method and route. |
| `tier` | Accuracy/disclosure class used in the validation matrix. |
| `non_disclosive` | Route-level disclosure assertion. Accepted evidence rows must remain `TRUE`. |
| `status` | `PASS` when `observed` is finite and within `tolerance`; otherwise `FAIL`. |

`assert_validation(rows)` is the executable assertion layer at the end
of each method vignette. It does not compute, tune, or change any
result. It stops the render if any row has `status != "PASS"` or if any
row is marked `non_disclosive == FALSE`. Therefore, a rendered method
article is evidence that the route stayed within the declared numerical
envelope and passed the route-level disclosure check for that fixture.

## Disclosure envelope

The product evidence accepts only route outputs that are consistent with
the project disclosure standard: model-scale outputs, scalar
diagnostics, and approved aggregate summaries may be returned;
patient-level identifiers, matched row indices, row scores, fitted
probabilities, BLUPs, cluster vectors, rank primitives, or other
observation-level reconstruction aids must not be returned to the
analyst.

Working quantities used by the algorithms can remain encrypted,
secret-shared, or server-side. The disclosure question in these
vignettes is what the product route returns to the analyst.

The DSLite examples disable trusted-peer pinning because DSLite is an
in-memory validation backend, not an Opal/Rock deployment with
persistent peer identity. That does not relax the method-level
disclosure checks.

## Method coverage

| Method | User-facing route | K=2 evidence | K\>=3 evidence | Current product path |
|----|----|----|----|----|
| PSI alignment | [`ds.vert.align()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_psi_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_psi_kge3.md) | Strict exact alignment counts/status only. |
| Descriptive statistics | [`ds.vert.desc()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_descriptive_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_descriptive_kge3.md) | Strict aggregate means and standard deviations. |
| Contingency tests | [`ds.vert.chisq()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md), [`ds.vert.fisher()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md), [`ds.vert.chisq_cross()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_contingency_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_contingency_kge3.md) | Strict count/test aggregate route. |
| Correlation | [`ds.vert.cor()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_correlation_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_correlation_kge3.md) | Practical strict correlation route. |
| PCA | [`ds.vert.pca()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_pca_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_pca_kge3.md) | Practical strict eigen/loading route. |
| GLM | [`ds.vert.glm()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_glm_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_glm_kge3.md) | Practical strict encrypted-label GLM route. |
| Inference helpers | [`ds.vert.confint()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md), [`ds.vert.wald()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md), [`ds.vert.contrast()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md), [`ds.vert.lr()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_inference_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_inference_kge3.md) | Exact algebra on accepted GLM output. |
| LASSO | [`ds.vert.lasso_proximal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_lasso_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_lasso_kge3.md) | Strict proximal route; validation checks the lambda=0 OLS limit. |
| Negative binomial | [`ds.vert.nb()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_negative_binomial_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_negative_binomial_kge3.md) | Accurate full-regression route and fast MoM approximation. |
| Cox PH | [`ds.vert.cox()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_cox_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_cox_kge3.md) | Practical non-disclosive profile/discrete route. |
| LMM | [`ds.vert.lmm()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_lmm_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_lmm_kge3.md) | Practical random-intercept route. |
| GEE | [`ds.vert.gee()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_gee_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_gee_kge3.md) | Practical aggregate estimating-equation route. |
| GLMM | [`ds.vert.glmm()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_glmm_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_glmm_kge3.md) | PQL fast route plus Laplace practical route. |
| IPW | [`ds.vert.ipw()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_ipw_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_ipw_kge3.md) | High-precision weighted outcome route. |
| MI | [`ds.vert.mi()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_mi_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_mi_kge3.md) | Practical pooled coefficient route. |
| Multinomial | [`ds.vert.multinom()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_multinomial_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_multinomial_kge3.md) | Joint softmax Newton route. |
| Ordinal | [`ds.vert.ordinal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md) | [K=2](https://isglobal-brge.github.io/dsVertClient/articles/validation_ordinal_k2.md) | [K\>=3](https://isglobal-brge.github.io/dsVertClient/articles/validation_ordinal_kge3.md) | Proportional-odds cumulative probability route. |

## Route-level result snapshot

The following table is copied from the current rendered validation
outputs. It is a static summary; the linked method articles contain the
executable code and the rendered DSLite logs.

| Method | K | Article | Route | Metric | Observed | Tolerance | Tier | Status |
|----|----|----|----|----|---:|---:|----|----|
| PSI | K=2 | [validation_psi_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_psi_k2.md) | `ds.vert.align` | `max(count_delta, correlation_delta)` | 0 | 0 | strict-precise | PASS |
| PSI | K\>=3 | [validation_psi_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_psi_kge3.md) | `ds.vert.align` | `max(count_delta, correlation_delta)` | 0 | 0 | strict-precise | PASS |
| Descriptive | K=2 | [validation_descriptive_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_descriptive_k2.md) | `ds.vert.desc` | `max(mean_sd_abs_delta)` | 0 | 0 | strict-precise | PASS |
| Descriptive | K\>=3 | [validation_descriptive_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_descriptive_kge3.md) | `ds.vert.desc` | `max(mean_sd_abs_delta)` | 0 | 0 | strict-precise | PASS |
| Contingency | K=2 | [validation_contingency_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_contingency_k2.md) | `ds.vert.chisq / ds.vert.fisher / ds.vert.chisq_cross` | `max(count_or_chisq_delta)` | 0 | 0 | strict-precise | PASS |
| Contingency | K\>=3 | [validation_contingency_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_contingency_kge3.md) | `ds.vert.chisq / ds.vert.fisher / ds.vert.chisq_cross` | `max(count_or_chisq_delta)` | 0 | 0 | strict-precise | PASS |
| Correlation | K=2 | [validation_correlation_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_correlation_k2.md) | `ds.vert.cor` | `correlation_max_abs_delta` | 1.24e-05 | 1e-04 | strict-practical | PASS |
| Correlation | K\>=3 | [validation_correlation_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_correlation_kge3.md) | `ds.vert.cor` | `correlation_max_abs_delta` | 1.24e-05 | 1e-04 | strict-practical | PASS |
| PCA | K=2 | [validation_pca_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_pca_k2.md) | `ds.vert.pca` | `max(eigen_loading_abs_delta)` | 3.29e-05 | 1e-04 | strict-practical | PASS |
| PCA | K\>=3 | [validation_pca_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_pca_kge3.md) | `ds.vert.pca` | `max(eigen_loading_abs_delta)` | 3.09e-05 | 1e-04 | strict-practical | PASS |
| GLM | K=2 | [validation_glm_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_glm_k2.md) | `ds.vert.glm` | `coef_max_abs_delta` | 6.24e-05 | 0.001 | strict-practical | PASS |
| GLM | K\>=3 | [validation_glm_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_glm_kge3.md) | `ds.vert.glm` | `coef_max_abs_delta` | 5.93e-05 | 0.001 | strict-practical | PASS |
| Inference | K=2 | [validation_inference_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_inference_k2.md) | `ds.vert.confint / ds.vert.wald / ds.vert.contrast / ds.vert.lr` | `algebra_max_abs_delta` | 0 | 0 | strict-precise | PASS |
| Inference | K\>=3 | [validation_inference_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_inference_kge3.md) | `ds.vert.confint / ds.vert.wald / ds.vert.contrast / ds.vert.lr` | `algebra_max_abs_delta` | 0 | 0 | strict-precise | PASS |
| LASSO | K=2 | [validation_lasso_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_lasso_k2.md) | `ds.vert.lasso_proximal(lambda=0)` | `lambda0_coef_abs_delta` | 0 | 0 | strict-precise | PASS |
| LASSO | K\>=3 | [validation_lasso_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_lasso_kge3.md) | `ds.vert.lasso_proximal(lambda=0)` | `lambda0_coef_abs_delta` | 0 | 0 | strict-precise | PASS |
| Negative binomial | K=2 | [validation_negative_binomial_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_negative_binomial_k2.md) | `ds.vert.nb(method='accurate')` | `coef_max_abs_delta` | 0.0000277 | 0.02 | strict-practical | PASS |
| Negative binomial | K=2 | [validation_negative_binomial_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_negative_binomial_k2.md) | `ds.vert.nb(method='fast')` | `max(beta_delta, theta_delta)` | 0.0010410 | 0.02 | fast-approximation | PASS |
| Negative binomial | K\>=3 | [validation_negative_binomial_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_negative_binomial_kge3.md) | `ds.vert.nb(method='accurate')` | `coef_max_abs_delta` | 0.0000317 | 0.02 | strict-practical | PASS |
| Negative binomial | K\>=3 | [validation_negative_binomial_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_negative_binomial_kge3.md) | `ds.vert.nb(method='fast')` | `max(beta_delta, theta_delta)` | 0.0010410 | 0.02 | fast-approximation | PASS |
| Cox | K=2 | [validation_cox_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_cox_k2.md) | `ds.vert.cox(method='profile')` | `coef_max_abs_delta` | 4.87e-05 | 0.001 | strict-practical | PASS |
| Cox | K\>=3 | [validation_cox_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_cox_kge3.md) | `ds.vert.cox(method='profile')` | `coef_max_abs_delta` | 4.87e-05 | 0.001 | strict-practical | PASS |
| LMM | K=2 | [validation_lmm_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_lmm_k2.md) | `ds.vert.lmm` | `fixed_effect_max_abs_delta` | 0.0006448 | 0.02 | strict-practical | PASS |
| LMM | K\>=3 | [validation_lmm_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_lmm_kge3.md) | `ds.vert.lmm` | `fixed_effect_max_abs_delta` | 0.0001614 | 0.02 | strict-practical | PASS |
| GEE | K=2 | [validation_gee_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_gee_k2.md) | `ds.vert.gee(corstr='independence')` | `coef_max_abs_delta` | 6.07e-05 | 0.01 | strict-practical | PASS |
| GEE | K\>=3 | [validation_gee_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_gee_kge3.md) | `ds.vert.gee(corstr='independence')` | `coef_max_abs_delta` | 6e-05 | 0.01 | strict-practical | PASS |
| GLMM | K=2 | [validation_glmm_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_glmm_k2.md) | `ds.vert.glmm(method='laplace')` | `max_fixed_effect_or_sigma_b2_abs_delta` | 0.0407200 | 0.060 | laplace-practical | PASS |
| GLMM | K=2 | [validation_glmm_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_glmm_k2.md) | `ds.vert.glmm(method='pql')` | `fixed_effect_max_abs_delta_and_pql_quality` | 0.0005547 | 0.005 | strict-pql | PASS |
| GLMM | K\>=3 | [validation_glmm_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_glmm_kge3.md) | `ds.vert.glmm(method='laplace')` | `max_fixed_effect_or_sigma_b2_abs_delta` | 0.0407200 | 0.060 | laplace-practical | PASS |
| GLMM | K\>=3 | [validation_glmm_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_glmm_kge3.md) | `ds.vert.glmm(method='pql')` | `fixed_effect_max_abs_delta_and_pql_quality` | 0.0005547 | 0.005 | strict-pql | PASS |
| IPW | K=2 | [validation_ipw_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_ipw_k2.md) | `ds.vert.ipw(precision='high')` | `weighted_outcome_coef_abs_delta` | 0.0001785 | 0.001 | strict-practical | PASS |
| IPW | K\>=3 | [validation_ipw_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_ipw_kge3.md) | `ds.vert.ipw(precision='high')` | `weighted_outcome_coef_abs_delta` | 0.0001785 | 0.001 | strict-practical | PASS |
| MI | K=2 | [validation_mi_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_mi_k2.md) | `ds.vert.mi` | `pooled_coef_abs_delta` | 0.01797 | 0.02 | strict-practical | PASS |
| MI | K\>=3 | [validation_mi_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_mi_kge3.md) | `ds.vert.mi` | `pooled_coef_abs_delta` | 5.51e-05 | 0.02 | strict-practical | PASS |
| Multinomial | K=2 | [validation_multinomial_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_multinomial_k2.md) | `ds.vert.multinom` | `class_probability_max_abs_delta` | 3.23e-05 | 0.005 | strict-practical | PASS |
| Multinomial | K\>=3 | [validation_multinomial_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_multinomial_kge3.md) | `ds.vert.multinom` | `class_probability_max_abs_delta` | 3.28e-05 | 0.005 | strict-practical | PASS |
| Ordinal | K=2 | [validation_ordinal_k2](https://isglobal-brge.github.io/dsVertClient/articles/validation_ordinal_k2.md) | `ds.vert.ordinal` | `cumulative_probability_max_abs_delta` | 1.61e-05 | 0.001 | strict-practical | PASS |
| Ordinal | K\>=3 | [validation_ordinal_kge3](https://isglobal-brge.github.io/dsVertClient/articles/validation_ordinal_kge3.md) | `ds.vert.ordinal` | `cumulative_probability_max_abs_delta` | 3.98e-05 | 0.001 | strict-practical | PASS |
