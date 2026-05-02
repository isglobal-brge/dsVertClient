# Federated multiple imputation with Rubin pooling

Fit a GLM on vertically partitioned DataSHIELD data with multiple
imputation of missing values. Imputations stay on the server holding the
missing variable (via `dsvertImputeColumnDS`); the client only ever sees
the \\M\\ pooled coefficient vectors and covariance matrices and applies
Rubin's rules client-side.

Protocol for each of \\m = 1..M\\: 1. On every server holding a column
with missingness, call `dsvertImputeColumnDS` with seed \\s_m\\. The
server draws a local imputation using a Bayesian-ridge model conditional
on the other complete-case columns available on that server. The imputed
column is written back into the aligned data frame under a per-round
name (e.g. `__dsvert_imp_<var>_<m>`). 2. Run
[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md)
on the imputed data, collect `beta_m` and `Cov(beta_m)`. 3. Client
accumulates \\(\beta_m, \mathrm{Cov}\_m)\\.

Rubin's pooling rules are applied client-side: \$\$\bar\beta =
\frac{1}{M}\sum_m \beta_m\$\$ \$\$W = \frac{1}{M}\sum_m
\mathrm{Cov}\_m\$\$ \$\$B = \frac{1}{M-1} \sum_m
(\beta_m-\bar\beta)(\beta_m-\bar\beta)^T\$\$ \$\$T = W + (1 + 1/M) B\$\$

## Usage

``` r
ds.vertMI(
  formula,
  data = NULL,
  impute_columns = NULL,
  m = 20L,
  family = "gaussian",
  max_iter = 50L,
  tol = 1e-04,
  lambda = 1e-04,
  verbose = TRUE,
  datasources = NULL,
  seed = 1L
)
```

## Arguments

- formula:

  Model formula.

- data:

  Aligned data-frame name.

- impute_columns:

  Character vector of column names with missingness that should be
  imputed (on whichever server holds them). Per-server column presence
  is auto-detected.

- m:

  Number of imputations (default 20).

- family:

  GLM family.

- max_iter:

  Inner `ds.vertGLM` `max_iter`.

- tol:

  Convergence tolerance for inner fits.

- lambda:

  L2 regularisation for inner fits.

- verbose:

  Print progress.

- datasources:

  DataSHIELD connection object.

- seed:

  RNG seed (default 1L). Per-round seed = `seed + m`.

## Value

A `ds.vertMI` object with fields `coefficients`, `covariance` (Rubin
total variance T), `std_errors`, `within`, `between`, `fmi` (fraction of
missing information), `m`, `family`, `fits` (list of the M inner
`ds.glm` fits).
