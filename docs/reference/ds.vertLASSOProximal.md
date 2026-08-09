# Gaussian LASSO via client-side coordinate descent

Solve the Gaussian LASSO objective entirely client-side. The preferred
route accepts a validated `ds.vertDPGaussian` release and optimises its
signed, projected DP sufficient statistics without another DSI call or
privacy cost. The historical `ds.glm` route is retained for
compatibility with already-created fits.

Normal equations: para Gaussian \\y = X\beta + \varepsilon\\ el
minimizador LASSO es \$\$\beta^\* = \arg\min\_{\beta} \tfrac{1}{2n}
\\y - X\beta\\^2 + \lambda \\\beta\\\_1.\$\$ El gradiente de la parte
cuadratica es \$\$\nabla f(\beta) = \tfrac{1}{n} X^\top X \cdot (\beta -
\hat\beta\_{OLS}).\$\$ ds.vertGLM ya expone \\\hat\beta\_{OLS}\\ y
\\\mathrm{Cov}(\hat\beta) = \sigma^2 (X^\top X)^{-1}\\. De ahi el
cliente reconstruye \\X^\top X / n = \mathrm{Cov}^{-1} \cdot
\hat\sigma^2 / n\\ e itera proximal-gradient puramente en memoria:
\\\beta\_{t+1} = S\_{\lambda/L}(\beta_t - \eta \nabla f(\beta_t))\\ con
\\L \ge \lambda\_{\max}(X^\top X / n)\\ (upper bound local) y \\S_t(x) =
\mathrm{sign}(x) \max(\|x\|-t, 0)\\ el operador soft-threshold.

## Usage

``` r
ds.vertLASSOProximal(
  fit,
  lambda,
  max_iter = 2000L,
  tol = 1e-09,
  keep_intercept = TRUE,
  warm_start = NULL,
  accelerate = TRUE
)
```

## Arguments

- fit:

  Preferably a `ds.vertDPGaussian` object. A historical unpenalised
  Gaussian `ds.glm` object with `$covariance` and `$n_obs` remains
  accepted for compatibility.

- lambda:

  Numeric. L1 penalty magnitude (on the 1/n-normalised objective).

- max_iter:

  Positive integer. Coordinate-descent passes (default 2000).

- tol:

  Numeric. Convergence tolerance on \\\\\beta\_{t+1} - \beta_t\\\\
  (default 1e-9).

- keep_intercept:

  Logical. If TRUE, do NOT penalise the intercept.

- warm_start:

  Numeric vector. For a `ds.vertDPGaussian` input this is an optional
  original-scale coefficient vector containing the intercept and every
  signed predictor. For a legacy `ds.glm`, it is beta_0.

- accelerate:

  Compatibility argument retained for historical callers. The current
  solver is coordinate descent, so this value has no effect.

## Value

An object of class `ds.vertLASSOProximal` with the proximal-MLE
coefficients, number of iterations, convergence flag, support, final
objective value, and the reconstructed Gram matrix used. The slot
`$comparison$coefficients_soft` reports the naive post-hoc
soft-thresholded OLS for comparison.

## Disclosure budget

Zero new DSI or MPC rounds. For a `ds.vertDPGaussian` input, every
lambda and iteration is deterministic post-processing of the same sticky
release and has additional privacy cost `(epsilon, delta) = (0, 0)`. No
classical sampling inference or coefficient confidence region is implied
by the projected noisy moments.

## See also

[`ds.vertLASSO`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO.md),
[`ds.vertLASSOCV`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOCV.md)
