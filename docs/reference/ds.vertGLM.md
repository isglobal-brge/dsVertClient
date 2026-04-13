# Generalized Linear Model for Vertically Partitioned Data

Fits a GLM across vertically partitioned data using Ring63 Beaver MPC
with DCF wide spline for the link function. The system auto-detects
which server holds each variable. Only p-dimensional aggregate gradients
are revealed to the client per iteration. No observation-level data is
ever disclosed.

## Usage

``` r
ds.vertGLM(
  formula,
  data = NULL,
  x_vars = NULL,
  y_server = NULL,
  family = "gaussian",
  max_iter = 100,
  tol = 1e-04,
  lambda = 1e-04,
  log_n = 12,
  verbose = TRUE,
  datasources = NULL,
  eta_privacy = "auto",
  data_name = NULL,
  y_var = NULL
)
```

## Arguments

- formula:

  A formula (e.g. `npreg ~ age + bmi + glu`) or character string
  (`"npreg ~ age + bmi + glu"`). Can also be a data_name string for
  backward compatibility.

- data:

  Character string. Name of the (aligned) data frame on each server.

- x_vars:

  Optional. Character vector of predictor names, or a named list mapping
  server names to variable vectors. If NULL and formula is used,
  extracted from the formula. If NULL and no formula, all available
  columns (minus y and IDs) are used.

- y_server:

  Character string. Name of the server holding the response.

- family:

  Character string. GLM family: "gaussian", "binomial", or "poisson".
  Default is "gaussian".

- max_iter:

  Integer. Maximum L-BFGS iterations. Default is 100.

- tol:

  Numeric. Convergence tolerance on coefficient change. Default is 1e-4.

- lambda:

  Numeric. L2 regularization parameter. Default is 1e-4.

- log_n:

  Integer. Legacy parameter (ignored). Kept for backward compatibility.

- verbose:

  Logical. Print progress messages. Default is TRUE.

- datasources:

  DataSHIELD connection object or list of connections.

- eta_privacy:

  Character. `"auto"` (default) selects `"k2_beaver"` for K=2 or
  `"secure_agg"` for K\>=3.

## Value

A list with class "ds.glm" containing:

- `coefficients`: Named coefficient vector (original scale)

- `std_errors`: Standard errors (finite-difference Hessian)

- `z_values`: z-statistics (coef / SE)

- `p_values`: Two-sided p-values

- `iterations`: Number of iterations

- `converged`: Logical

- `family`: Family used

- `n_obs`: Number of observations

- `deviance`: Residual sum of squares

- `pseudo_r2`: 1 - deviance/null_deviance

## Details

### Protocol

All computation uses Ring63 fixed-point arithmetic with Beaver MPC:

1.  **Transport keys**: X25519 keypairs on all servers (~0.5s)

2.  **Standardize**: Each server standardizes its features

3.  **Input sharing**: Features split into additive Ring63 shares
    between 2 DCF parties. Non-DCF servers contribute shares.

4.  **L-BFGS loop**: Per iteration:

    - Compute eta shares (Ring63 matrix-vector)

    - DCF wide spline for sigmoid/exp (binomial/Poisson) or identity
      link (Gaussian)

    - Beaver matvec for gradient (server-generated triples)

    - Client aggregates Ring63 shares -\> p gradient scalars

    - L-BFGS quasi-Newton update

5.  **SE**: p+1 gradient evaluations (finite-difference Hessian)

6.  **Deviance**: Beaver dot-product for residual sum of squares

7.  **Unstandardize**: Coefficients + SE via Jacobian transform

### Security

No observation-level data is disclosed. The client sees only
p-dimensional aggregate gradients per iteration. Beaver triples are
generated server-side (never seen by client). Dealer rotation for K\>=4
ensures the analyst must compromise (K-1)/K servers to extract data.

## References

van Kesteren, E.J. et al. (2019). Privacy-preserving generalized linear
models using distributed block coordinate descent. arXiv:1911.03183.

Mouchet, C. et al. (2021). "Multiparty Homomorphic Encryption from
Ring-Learning-With-Errors". *Proceedings on Privacy Enhancing
Technologies*.

## See also

[`ds.vertCor`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCor.md)
for correlation analysis,
[`ds.vertPCA`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertPCA.md)
for PCA analysis

## Examples

``` r
if (FALSE) { # \dontrun{
# Simplest: formula interface (auto-detects everything)
model <- ds.vertGLM(npreg ~ age + bmi + glu + bp + skin,
                     data = "DA", family = "gaussian")

# String formula
model <- ds.vertGLM("diabetes ~ age + bmi + glu",
                     data = "DA", family = "binomial")

# Auto-detect all features
model <- ds.vertGLM("DA", "npreg", family = "poisson")

# Manual server mapping (legacy)
model <- ds.vertGLM("DA", "npreg",
  list(s1 = c("age", "bmi"), s2 = c("glu", "bp")),
  y_server = "s2", family = "gaussian")
} # }
```
