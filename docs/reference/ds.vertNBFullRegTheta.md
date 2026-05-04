# Federated NB regression with full-regression theta refinement

Extends `ds.vertNB` (which uses the iid-mu profile MLE for theta:
assumes mu_i == ybar when evaluating the profile score) with a
variance-corrected refinement that accounts for mu_i variation across
patients without requiring per-patient MPC reveals.

## Usage

``` r
ds.vertNBFullRegTheta(
  formula,
  data = NULL,
  theta = NULL,
  joint = TRUE,
  theta_max_iter = 5L,
  theta_tol = 0.001,
  variant = "full_reg_nd",
  beta_max_iter = 2L,
  beta_tol = 1e-04,
  compute_covariance = TRUE,
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

- variant:

  Character. `"full_reg_nd"` (default) runs the non-disclosive
  share-domain full-regression theta refinement. `"iid_mu"` returns the
  unmodified `ds.vertNB` result. `"corrected"` applies the aggregate
  variance correction described in Details. The legacy disclosive
  `"full_reg"` eta-transport path has been removed.

- beta_max_iter:

  Integer. Maximum beta refinements for the non-disclosive
  full-regression theta variant.

- beta_tol:

  Numeric. Relative convergence tolerance for beta refinements in the
  non-disclosive full-regression theta variant.

- compute_covariance:

  Logical. If `TRUE`, request covariance and standard-error diagnostics
  where the selected beta path supports them.

- verbose:

  Logical. Print stage-by-stage progress.

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

- ...:

  Extra arguments forwarded to `ds.vertGLM`.

## Value

Object of class `c("ds.vertNBFullRegTheta", "ds.vertNB")`. Fields as
`ds.vertNB`, plus `$theta_iid` (original iid-mu estimate) and
`$variance_correction` (the \\\hat V\_\mu\\ used). For the
non-disclosive full-regression variant, the object also contains
`$theta_trace`, `$theta_iter`, and `$theta_converged`.

## Details

The NB(mu_i, theta) log-likelihood score for theta is \$\$s(\theta) =
\sum_i \psi(y_i + \theta) - n \psi(\theta) + n \log \theta - \sum_i
\log(\mu_i + \theta).\$\$ The iid-mu approximation used in `ds.vertNB`
replaces the last term by \\n \log(\bar y + \theta)\\. For homogeneous
cohorts (small \\\text{Var}(\mu)\\) this is tight; for regression-rich
settings the bias on theta can reach ~16% (quine, overdispersed counts)
relative to [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

A first-order correction uses the aggregate marginal variance of y
decomposed via the NB law of total variance: \\\text{Var}(y) =
E\[\mu\] + E\[\mu^2\]/\theta + \text{Var}(\mu)\\. With \\\bar y\\ and
\\s_y^2\\ (scalar aggregates from `dsvertLocalMomentsDS`) and the iid
theta_0 as seed, we refine via Brent root-finding on the corrected score
\\s\_{\text{corr}}(\theta) = s\_{\text{iid}}(\theta) - \frac{1}{2}
\frac{n \\ \hat{V}\_\mu}{(\bar y + \theta)^2}\\ where \\\hat V\_\mu\\ is
the aggregate estimate of \\\text{Var}(\mu)\\ and the second term is the
Taylor correction to \\\sum_i \log(\mu_i + \theta)\\ around \\\bar y\\.

The aggregate \\\hat V\_\mu\\ is computed as \\\hat V\_\mu = \max(0,\\
s_y^2 - \bar y - \bar y^2 / \hat\theta_0)\\ – the portion of total y
variance not explained by NB conditional variance \\\mu +
\mu^2/\theta\\. All quantities are scalar aggregates; no per-patient
disclosure.

Full-share-space Clenshaw evaluation of \\\sum_i \log(\mu_i + \theta)\\
using the shipped `Ring127LogShiftPlaintext` Chebyshev primitive + DCF
argument reduction is a stricter variant scheduled separately; the
first-order correction here closes the bulk of the iid-mu bias (on
quine: 16% -\> 4-5%) without any new MPC machinery.

## See also

[`ds.vertNB`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNB.md)
