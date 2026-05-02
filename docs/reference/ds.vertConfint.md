# Wald confidence intervals for ds.vertGLM coefficients

Return Wald-type confidence intervals using the finite-difference
standard errors already stored in a ds.glm object. Observation-level
quantities are never touched; this is a scalar client-side
transformation.

## Usage

``` r
ds.vertConfint(fit, parm = NULL, level = 0.95)
```

## Arguments

- fit:

  A ds.glm object.

- parm:

  Optional character vector of coefficient names to report; default all.

- level:

  Confidence level, default 0.95.

## Value

A data frame with columns `estimate`, `std_error`, `lower`, `upper`; row
names are coefficient names.
