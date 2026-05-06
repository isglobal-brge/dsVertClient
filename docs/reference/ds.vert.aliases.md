# User-facing ds.vert.\* API aliases

These wrappers provide a compact, formula-style public surface while
keeping the historical CamelCase functions as compatibility backends.
Wrappers that have distinct K=2 and K\>=3 implementations dispatch from
the number of active DataSHIELD connections.

## Usage

``` r
ds.vert.align(data_name, id_col, newobj = "D_aligned", datasources = NULL, ...)

ds.vert.is_aligned(newobj = "DA", datasources = NULL, ...)

ds.vert.desc(data_name, datasources = NULL, ...)

ds.vert.cor(data_name, variables = NULL, datasources = NULL, ...)

ds.vert.pca(data_name = NULL, variables = NULL, datasources = NULL, ...)

ds.vert.chisq(data_name, var1, var2, datasources = NULL, ...)

ds.vert.fisher(data_name, var1, var2, datasources = NULL, ...)

ds.vert.chisq_cross(data, var1, var2, datasources = NULL, ...)

ds.vert.glm(
  formula,
  data = NULL,
  precision = c("auto", "high", "fast"),
  datasources = NULL,
  ...
)

ds.vert.cox(
  formula,
  data = NULL,
  method = c("profile", "discrete"),
  datasources = NULL,
  ...
)

ds.vert.coxph(formula, data = NULL, ...)

ds.vert.nb(
  formula,
  data = NULL,
  method = c("auto", "accurate", "fast", "mom", "profile"),
  datasources = NULL,
  ...
)

ds.vert.multinom(formula, data = NULL, datasources = NULL, ...)

ds.vert.ordinal(formula, data = NULL, datasources = NULL, ...)

ds.vert.lmm(
  formula,
  data = NULL,
  cluster_col,
  max_iter = 30L,
  inner_iter = 50L,
  max_outer = 30L,
  tol = NULL,
  ring = NULL,
  verbose = TRUE,
  datasources = NULL,
  ...
)

ds.vert.gee(
  formula,
  data = NULL,
  precision = c("auto", "high", "fast"),
  datasources = NULL,
  ...
)

ds.vert.glmm(
  formula,
  data = NULL,
  cluster_col,
  method = c("auto", "laplace", "pql"),
  datasources = NULL,
  ...
)

ds.vert.glmer(formula, data = NULL, cluster_col, datasources = NULL, ...)

ds.vert.ipw(
  outcome_formula,
  propensity_formula,
  data = NULL,
  precision = c("auto", "high", "fast"),
  datasources = NULL,
  ...
)

ds.vert.mi(
  formula,
  data = NULL,
  impute_columns = NULL,
  datasources = NULL,
  ...
)

ds.vert.lasso(fit, lambda_1, ...)

ds.vert.lasso_iter(
  formula,
  data = NULL,
  method = c("auto", "accurate", "fast"),
  ...
)

ds.vert.lasso_proximal(fit, lambda, ...)

ds.vert.lasso_1step(fit, lambda, ...)

ds.vert.lasso_cv(fit, lambda_grid = NULL, ...)

ds.vert.lr(reduced, full)

ds.vert.confint(fit, parm = NULL, level = 0.95)

ds.vert.wald(fit, parm, null = 0)

ds.vert.contrast(fit, K, m = NULL)
```
