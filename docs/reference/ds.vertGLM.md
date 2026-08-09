# DP-capsule GLM compatibility frontdoor

This public frontdoor has one available analysis route: an explicit
`dp_analysis_id` with `family = "gaussian"` delegates to
[`ds.vertDPGaussian()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPGaussian.md)
and returns that signed, contribution-bounded sticky joint-DP capsule
estimand. A call without `dp_analysis_id` raises a typed
`dsvert_route_unavailable` condition before any DSI call. An explicit
`formal_analysis_id` for binomial or Poisson also fails before DSI until
its durable worker and single common joint-DP opening are promoted. The
retained iterative Ring/Beaver code below the local gate is unreachable
through this frontdoor and carries no public disclosure, accuracy, or
availability claim.

## Usage

``` r
ds.vertGLM(
  formula,
  data = NULL,
  x_vars = NULL,
  y_server = NULL,
  family = "gaussian",
  max_iter = 100,
  tol = 1e-04,
  lambda = 0,
  log_n = 12,
  offset = NULL,
  weights = NULL,
  ring = 63L,
  binomial_sigmoid_intervals = NULL,
  verbose = TRUE,
  datasources = NULL,
  eta_privacy = "auto",
  keep_session = FALSE,
  no_intercept = FALSE,
  std_mode = "full",
  start = NULL,
  compute_se = TRUE,
  compute_deviance = TRUE,
  gradient_only = FALSE,
  data_name = NULL,
  y_var = NULL,
  missing = "fail",
  numeric_backend = "auto",
  dp_analysis_id = NULL,
  formal_analysis_id = NULL
)
```

## Arguments

- formula, data, x_vars, y_server:

  Additive model specification, aligned data name and optional
  signed-artifact ownership checks for the Gaussian capsule route.

- family:

  Must be `"gaussian"` with `dp_analysis_id`. Binomial and Poisson are
  not available through the public frontdoor.

- max_iter, tol, log_n, offset, weights, ring,
  binomial_sigmoid_intervals, eta_privacy, keep_session, std_mode,
  start, compute_se, compute_deviance, gradient_only, numeric_backend:

  Retained legacy arguments. They are rejected when explicitly supplied
  to the Gaussian capsule adapter, and the no-id legacy route is
  unavailable.

- lambda, no_intercept, data_name, y_var, missing:

  Gaussian capsule estimand selectors. `lambda` is the explicit
  non-negative ridge penalty; `missing`, when supplied, must be
  `"complete_case_capsule"`.

- verbose, datasources:

  Progress flag and DataSHIELD connections used only after the Gaussian
  signed-artifact request has passed local checks.

- dp_analysis_id:

  Custodian-configured signed bounded Gaussian artifact id. This is
  required for the available route.

- formal_analysis_id:

  Reserved custodian-configured binomial/Poisson selector. Supplying it
  returns a typed `dsvert_formal_glm_frontdoor_unavailable` condition
  before DSI.

## Value

With a valid Gaussian `dp_analysis_id`, a `ds.vertDPGaussian` object
containing bounded noisy sufficient- statistic regression output and no
classical standard errors, p-values, individual fitted values, residuals
or scores. All other routes raise a typed condition before DSI and
return no fitted object.

## Details

**Available route.** Supply an additive formula, an aligned data name,
`family = "gaussian"` and a custodian-configured `dp_analysis_id`. The
signed artifact owns clipping bounds, the complete-case cohort,
contribution caps, privacy parameters and variable ownership. The
adapter accepts only arguments that describe that bounded Gaussian
estimand; it never falls back to the retired iterative GLM.

**Unavailable routes.** Default/no-id calls and all legacy iterative
routes stop locally with zero DSI calls. The `formal_analysis_id`
selector is reserved for binomial/logit and Poisson/log models but is
also sealed locally in this release. No binomial or Poisson fit is
therefore returned by this function.

## See also

[`ds.vertDPGaussian`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertDPGaussian.md),
[`ds.vertMethodStatus`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMethodStatus.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- ds.vertGLM(
  y ~ x1 + x2, data = "D", family = "gaussian",
  dp_analysis_id = "custodian-gaussian-analysis")
} # }
```
