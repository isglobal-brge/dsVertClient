# Proper LASSO via client-side proximal gradient on normal equations

Cierra el gap OLS-soft-threshold -\> proper-LASSO usando solo cantidades
ya expuestas por ds.vertGLM (beta, covariance, n), sin anadir ninguna
ronda MPC adicional.

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

  A `ds.glm` object from `ds.vertGLM(family="gaussian")`. Must have
  `$covariance` and `$n_obs` populated.

- lambda:

  Numeric. L1 penalty magnitude (on the 1/n-normalised objective).

- max_iter:

  Integer. Outer proximal-gradient iterations (default 200).

- tol:

  Numeric. Convergence tolerance on \\\\\beta\_{t+1} - \beta_t\\\\
  (default 1e-7).

- keep_intercept:

  Logical. If TRUE, do NOT penalise the intercept.

- warm_start:

  Numeric vector. Optional beta_0 (default = beta_OLS).

- accelerate:

  Logical (default TRUE). Use Beck-Teboulle FISTA acceleration on the
  proximal-gradient inner loop. FALSE falls back to plain ISTA.

## Value

An object of class `ds.vertLASSOProximal` with the proximal-MLE
coefficients, number of iterations, convergence flag, support, final
objective value, and the reconstructed Gram matrix used. The slot
`$comparison$coefficients_soft` reports the naive post-hoc
soft-thresholded OLS for comparison.

## P3 disclosure budget

Zero new MPC rounds beyond the initial `ds.vertGLM` call. All iteration
is client-side on quantities already in the `fit` object. The
optimisation volume argument is therefore moot – no additional reveals
happen regardless of iteration count.

## See also

[`ds.vertLASSO`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSO.md),
[`ds.vertLASSOCV`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertLASSOCV.md)
