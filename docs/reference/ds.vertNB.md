# Federated negative binomial regression with dispersion estimate

Fit a negative binomial GLM on vertically partitioned DataSHIELD data.
Two-stage:

1.  Fit \\\hat\beta\\ via the dsVert Poisson GLM (identical score under
    canonical log link so the point estimate is the same; only the
    covariance differs).

2.  Estimate \\\hat\theta\\ by client-side Newton-Raphson on the NB
    profile log-likelihood, evaluated at each candidate \\\theta\\
    through `dsvertNBProfileSumsDS` which returns
    \\\sum\psi(y_i+\theta)\\, \\\sum\psi_1(y_i+\theta)\\, \\n\\, and
    \\\bar y\\ as scalar aggregates on the outcome server. The
    Anscombe/Lawless score
    \\s(\theta)=\sum\psi(y_i+\theta)-n\psi(\theta)+n\log(\theta/(\bar
    y+\theta))\\ is zero at the MLE; derivative
    \\s'(\theta)=\sum\psi_1(y_i+\theta)-n\psi_1(\theta)+n\[1/\theta-1/(\bar
    y+\theta)\]\\ supplies the Newton step. Initial value comes from
    method-of-moments \\\hat\theta_0=\bar y^2/(s_y^2-\bar y)\\.

3.  Rescale the Poisson SE by \\\sqrt{1 + \bar y / \hat\theta}\\ so the
    reported z-stats reflect NB variance inflation.

All aggregates are scalar sums over y; no per-patient disclosure.

## Usage

``` r
ds.vertNB(
  formula,
  data = NULL,
  theta = NULL,
  joint = TRUE,
  theta_max_iter = 5L,
  theta_tol = 0.001,
  verbose = TRUE,
  datasources = NULL,
  ...
)
```

## Arguments

- formula:

  Model formula for the count outcome (LHS) and linear predictor (RHS).

- data:

  Aligned data-frame name on each server.

- theta:

  Optional fixed dispersion parameter; if supplied, the theta refinement
  step is skipped and only the Poisson beta path runs with this theta
  plugged into the SE rescaling.

- joint:

  Logical. If TRUE (default), iterate between the Poisson \\\hat\beta\\
  update and a \\\hat\theta\\ update until both converge. If FALSE,
  return the Poisson \\\hat\beta\\ with a single one-shot MoM
  \\\hat\theta\\.

- theta_max_iter:

  Outer iterations for the joint update (default 5). Each iteration
  refits the Poisson GLM with the current theta-adjusted mean estimate.

- theta_tol:

  Convergence tolerance on the relative change in theta during the
  Newton refinement.

- verbose:

  Logical. Print stage-by-stage progress.

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

- ...:

  Extra arguments forwarded to `ds.vertGLM`.
