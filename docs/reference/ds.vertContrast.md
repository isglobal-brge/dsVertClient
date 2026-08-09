# Multi-coefficient Wald test via linear contrast K\*beta

Test H0: K \* beta = m using the multivariate Wald statistic W = (K \*
beta_hat - m)^T inv(K \* Cov \* K^T) (K \* beta_hat - m), using F = W /
rank(K) for Gaussian fits with residual degrees of freedom, and the
asymptotic chi-square reference otherwise. Requires the fit's full
covariance matrix.

## Usage

``` r
ds.vertContrast(fit, K, m = NULL)
```

## Arguments

- fit:

  A ds.glm object with a non-NULL `covariance` slot.

- K:

  Contrast matrix: numeric matrix with ncol equal to the number of
  coefficients. Rows define the contrasts under test. Alternatively a
  named-coef character vector (treated as indicator rows) or a character
  RHS parsed against the coefficient names.

- m:

  Null vector (length nrow(K)); default zero.

## Value

A list of class ds.vertContrast with estimates, variance, reference
`statistic`, raw `wald_statistic`, `distribution`, degrees of freedom
and `p_value`.
