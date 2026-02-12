# Generalized Linear Model for Vertically Partitioned Data

Client-side function that fits a Generalized Linear Model across
vertically partitioned data using Block Coordinate Descent with
encrypted labels. The response variable y only needs to exist on ONE
server (the "label server"). Non-label servers compute gradient updates
using y encrypted under the MHE collective public key, and only the
aggregated p_k-length gradient is revealed via threshold decryption.

## Usage

``` r
ds.vertGLM(
  data_name,
  y_var,
  x_vars,
  y_server = NULL,
  family = "gaussian",
  max_iter = 100,
  tol = 1e-04,
  lambda = 1e-04,
  log_n = 12,
  log_scale = 40,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- data_name:

  Character string. Name of the (aligned) data frame on each server.

- y_var:

  Character string. Name of the response variable (must exist on the
  label server specified by `y_server`).

- x_vars:

  A named list where each name corresponds to a server name and each
  element is a character vector of predictor variable names from that
  server.

- y_server:

  Character string. Name of the server holding the response variable.
  This server uses plaintext IRLS; all other servers use the encrypted
  gradient protocol.

- family:

  Character string. GLM family: "gaussian", "binomial", or "poisson".
  Default is "gaussian".

- max_iter:

  Integer. Maximum number of BCD iterations. Default is 100.

- tol:

  Numeric. Convergence tolerance on coefficient change. Default is 1e-4
  (accounts for CKKS approximation noise).

- lambda:

  Numeric. L2 regularization parameter. Default is 1e-4.

- log_n:

  Integer. CKKS ring dimension parameter (12, 13, or 14). Default is 12
  (2048 slots, supports up to 2048 observations).

- log_scale:

  Integer. CKKS scale parameter. Default is 40.

- verbose:

  Logical. Print progress messages. Default is TRUE.

- datasources:

  DataSHIELD connection object or list of connections. If NULL, uses all
  available connections.

## Value

A list with class "ds.glm" containing:

- `coefficients`: Named vector of coefficient estimates (on original
  scale, including intercept)

- `iterations`: Number of iterations until convergence

- `converged`: Logical indicating convergence

- `family`: Family used

- `n_obs`: Number of observations

- `n_vars`: Number of predictor variables (including intercept)

- `lambda`: Regularization parameter used

- `deviance`: Residual deviance of the fitted model

- `null_deviance`: Null deviance (intercept-only model)

- `pseudo_r2`: McFadden's pseudo R-squared

- `aic`: Akaike Information Criterion

- `y_server`: Name of the label server

- `call`: The matched call

## Details

### Feature Standardization

Features are automatically standardized (centered and scaled) on each
server before BCD to ensure fast convergence. For Gaussian family, the
response is also standardized. Coefficients are transformed back to the
original scale after convergence, and an intercept is computed.

### Encrypted-Label BCD-IRLS Protocol

The response variable y resides on a single "label server". Non-label
servers never see y in plaintext. The protocol proceeds as:

1.  **MHE Key Setup**: All servers generate key shares and combine them
    into a Collective Public Key (CPK) with Galois keys.

2.  **Standardize**: Each server standardizes its features.

3.  **Encrypt y**: The label server encrypts (standardized) y under the
    CPK and distributes the ciphertext to non-label servers.

4.  **BCD Loop**: For each iteration, each server updates its block of
    coefficients on the standardized scale.

5.  **Unstandardize**: Coefficients are transformed back to the original
    scale and an intercept is computed.

6.  **Deviance**: Computed on the label server using plaintext y and the
    final linear predictor (original scale).

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
x_vars <- list(
  server1 = c("age", "bmi"),
  server2 = c("glucose"),
  server3 = c("cholesterol", "heart_rate")
)

# Gaussian GLM (bp on server2)
model <- ds.vertGLM("D_aligned", "bp", x_vars,
                     y_server = "server2", family = "gaussian")
print(model)
} # }
```
