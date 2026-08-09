# Univariate Wald test for a single ds.vertGLM coefficient

Test H0: beta_j = null against a two-sided alternative using the
diagonal statistic (estimate - null) / SE. Gaussian fits use Student t
with residual degrees of freedom; binomial and Poisson fits use the
asymptotic normal reference.

## Usage

``` r
ds.vertWald(fit, parm, null = 0)
```

## Arguments

- fit:

  A ds.glm object.

- parm:

  Single coefficient name.

- null:

  Null value (default 0).

## Value

List with estimate, SE, statistic, distribution, p_value and null;
Gaussian results include `t` and `df`, other families include `z`.
