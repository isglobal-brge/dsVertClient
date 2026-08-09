# Public ds.vert.\* API aliases

These wrappers preserve the compact formula-style API names.
Availability is defined by
[`ds.vertMethodStatus()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)
and the selected signed capsule, not by the presence of an alias.
Aliases for quarantined families raise a typed
`dsvert_route_unavailable` condition before any DSI call.
`ds.vert.glm()` can delegate an explicit Gaussian `dp_analysis_id` to
[`ds.vertDPGaussian()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPGaussian.md);
its no-id legacy route is unavailable, and the reserved binomial/Poisson
`formal_analysis_id` route also fails before DSI in this release. No
alias re-enables a retired remote endpoint or weakens the
signed-artifact and custodian-owned policy gates of an available
backend.

## Usage

``` r
ds.vert.align(data_name, id_col, newobj = "D_aligned", datasources = NULL, ...)

ds.vert.is_aligned(newobj = "DA", datasources = NULL, ...)

ds.vert.desc(data_name, datasources = NULL, ...)

ds.vert.cor(data_name, variables = NULL, datasources = NULL, ...)

ds.vert.pca(data_name = NULL, variables = NULL, datasources = NULL, ...)

ds.vert.chisq(data_name, var1 = NULL, var2 = NULL, datasources = NULL, ...)

ds.vert.fisher(data_name, var1 = NULL, var2 = NULL, datasources = NULL, ...)

ds.vert.chisq_cross(data, var1 = NULL, var2 = NULL, datasources = NULL, ...)

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

ds.vert.coxph(formula, data = NULL, method = "profile", ...)

ds.vert.nb(formula, data = NULL, method = "accurate", datasources = NULL, ...)

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
  method = "pql",
  datasources = NULL,
  ...
)

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

## Arguments

- data_name, data, formula, datasources:

  Aligned data-frame symbol (or model `formula`) and the DataSHIELD
  connections.

- id_col, newobj:

  Record-identifier column and output symbol for alignment.

- ...:

  Further arguments forwarded to the backend.

- variables, var1, var2:

  Column selections for descriptive / bivariate routes.

- precision, method, ring, verbose:

  Binomial-sigmoid precision preset, estimator/route selector,
  fixed-point ring, and progress flag.

- cluster_col:

  Grouping column for the mixed-model routes.

- max_iter, inner_iter, max_outer, tol:

  Iteration caps and convergence tolerance for the iterative fits.

- outcome_formula, propensity_formula:

  Outcome and propensity models (IPW).

- impute_columns, m:

  Columns to impute and number of imputations (MI).

- fit, reduced, full, parm, level, null, K:

  Inference-helper inputs (fitted object, nested models, parameter,
  confidence level, null value, class count).

- lambda, lambda_1, lambda_grid:

  Penalty or penalty grid for the LASSO routes.
