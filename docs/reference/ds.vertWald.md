# Univariate Wald test for a single ds.vertGLM coefficient

Test H0: beta_j = null against a two-sided alternative using the
diagonal Wald statistic (estimate - null) / SE. A convenience wrapper
over the z_values / p_values already stored in a ds.glm object, extended
to non-zero nulls.

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

List with estimate, SE, z, p_value, null.
