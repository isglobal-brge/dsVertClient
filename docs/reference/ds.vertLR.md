# Likelihood-ratio test on two nested ds.vertGLM fits

Compute the standard LR chi-square statistic between a reduced and a
full ds.vertGLM fit. Both fits must be on the same cohort and the
reduced model must be nested within the full model.

## Usage

``` r
ds.vertLR(reduced, full)
```

## Arguments

- reduced:

  ds.glm object with fewer coefficients.

- full:

  ds.glm object with more coefficients.

## Value

A list with:

- `statistic`: 2 \* (reduced\$deviance - full\$deviance)

- `df`: full\$n_vars - reduced\$n_vars

- `p_value`: upper-tail chi-square p-value

- `deviance_reduced`, `deviance_full`

## Details

This is a pure client-side computation over the scalar deviances already
returned by ds.vertGLM; no additional MPC round is performed and no
observation-level information is exposed.
