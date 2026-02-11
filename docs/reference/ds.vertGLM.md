# Generalized Linear Model for Vertically Partitioned Data

Client-side function that fits a Generalized Linear Model across
vertically partitioned data using Block Coordinate Descent.

## Usage

``` r
ds.vertGLM(
  data_name,
  y_var,
  x_vars,
  family = "gaussian",
  max_iter = 100,
  tol = 1e-06,
  lambda = 1e-04,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- data_name:

  Character string. Name of the (aligned) data frame on each server.

- y_var:

  Character string. Name of the response variable (must exist on ALL
  servers).

- x_vars:

  A named list where each name corresponds to a server name and each
  element is a character vector of predictor variable names from that
  server.

- family:

  Character string. GLM family: "gaussian", "binomial", "poisson",
  "Gamma", or "inverse.gaussian". Default is "gaussian".

- max_iter:

  Integer. Maximum number of iterations. Default is 100.

- tol:

  Numeric. Convergence tolerance. Default is 1e-6.

- lambda:

  Numeric. L2 regularization parameter. Default is 1e-4.

- verbose:

  Logical. Print iteration progress. Default is TRUE.

- datasources:

  DataSHIELD connection object or list of connections. If NULL, uses all
  available connections.

## Value

A list with class "ds.glm" containing:

- `coefficients`: Named vector of coefficient estimates

- `iterations`: Number of iterations until convergence

- `converged`: Logical indicating convergence

- `family`: Family used

- `n_obs`: Number of observations

- `n_vars`: Number of predictor variables

- `deviance`: Residual deviance of the fitted model

- `null_deviance`: Null deviance (intercept-only model)

- `pseudo_r2`: McFadden's pseudo R-squared

- `aic`: Akaike Information Criterion

- `call`: The matched call

## Details

This function implements the Block Coordinate Descent (BCD) algorithm
for privacy-preserving GLM fitting on vertically partitioned data.

The algorithm iteratively:

1.  For each partition i:

    - Compute eta_other = sum of X_j \* beta_j for all j != i

    - Update beta_i using IRLS with eta_other

    - Share new eta_i = X_i \* beta_i (NOT raw data or coefficients)

2.  Check convergence (sum of \|beta_new - beta_old\| \< tol)

3.  Repeat until converged or max_iter reached

Privacy is preserved because:

- Only linear predictor contributions (eta) are shared

- Raw data never leaves servers

- Final coefficients are computed locally on each server

## References

van Kesteren, E.J. et al. (2019). Privacy-preserving generalized linear
models using distributed block coordinate descent. arXiv:1911.03183.

## See also

[`glm`](https://rdrr.io/r/stats/glm.html) for standard GLM fitting

## Examples

``` r
if (FALSE) { # \dontrun{
# Define predictor variables per server
x_vars <- list(
  server1 = c("age", "weight"),
  server2 = c("height", "bmi"),
  server3 = c("glucose", "cholesterol")
)

# Fit Gaussian GLM (linear regression)
model <- ds.vertGLM("D_aligned", "outcome", x_vars, family = "gaussian")
print(model)

# Fit logistic regression
model_logit <- ds.vertGLM("D_aligned", "binary_outcome", x_vars,
                          family = "binomial")
} # }
```
