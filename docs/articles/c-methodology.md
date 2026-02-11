# Methodology

This vignette explains the mathematical methodology behind
dsVertClient’s privacy-preserving algorithms. Understanding these
methods helps users interpret results correctly and appreciate the
privacy guarantees.

## Overview of Methods

dsVertClient implements three main distributed algorithms:

1.  **Block SVD** - For correlation and PCA
2.  **Block Coordinate Descent (BCD)** - For GLM coefficient estimation
3.  **Iteratively Reweighted Least Squares (IRLS)** - For non-Gaussian
    GLMs

All methods share a key principle: **decompose computations so each
institution only sees its own data, while sharing only aggregate
statistics**.

------------------------------------------------------------------------

## Block Singular Value Decomposition (Block SVD)

### The Problem

Given data matrix $`X`$ split column-wise across $`K`$ institutions:

``` math
X = [X_1 | X_2 | ... | X_K]
```

where institution $`k`$ holds $`X_k`$ (its subset of columns), we want
to compute:

- Correlation matrix: $`R = X^T X / (n-1)`$ (for standardized data)
- Principal components: From SVD of $`X`$

But no institution can see the full $`X`$.

### The Solution

The key insight is that we can compute the combined SVD by:

1.  Each institution $`k`$ computes SVD of its own data:
    $`X_k = U_k D_k V_k^T`$
2.  Each institution returns $`U_k D_k`$ (the “scores” component)
3.  The client stacks these: $`[U_1 D_1 | U_2 D_2 | ... | U_K D_K]`$
4.  A final SVD on this combined matrix yields the global result

**Why this preserves privacy:**

- $`U_k D_k`$ is an aggregate transformation of the data
- It has the same number of rows as observations but typically far fewer
  columns
- Individual data values cannot be recovered from $`U_k D_k`$

### Mathematical Details

For correlation analysis on standardized data (mean=0, sd=1):

``` math
\text{Cor}(X) = \frac{X^T X}{n-1} = \frac{V D^2 V^T}{n-1}
```

The Block SVD approach computes this as:

1.  Stack the partial results:
    $`M = [U_1 D_1 | U_2 D_2 | ... | U_K D_K]`$
2.  Compute SVD of $`M`$: $`M = \tilde{U} \tilde{D} \tilde{V}^T`$
3.  The correlation matrix is:
    $`R = \tilde{V} \tilde{D}^2 \tilde{V}^T / (n-1)`$

For PCA, the components are:

- **PC scores**: $`\tilde{U} \tilde{D}`$
- **PC loadings**: $`\tilde{V}`$
- **Variance explained**: $`\tilde{D}^2 / (n-1)`$

### Reference

Iwen, M. & Ong, B.W. (2016). [A distributed and incremental SVD
algorithm for agglomerative data analysis on large
networks](https://doi.org/10.1137/16M1058467). *SIAM Journal on Matrix
Analysis and Applications*.

------------------------------------------------------------------------

## Block Coordinate Descent (BCD) for GLMs

### The Standard GLM Framework

A Generalized Linear Model relates response $`y`$ to predictors $`X`$
through:

``` math
g(\mu) = \eta = X\beta
```

where:

- $`g(\cdot)`$ is the link function
- $`\mu = E[y]`$ is the expected response
- $`\eta`$ is the linear predictor
- $`\beta`$ are the coefficients to estimate

### The Distributed Setting

With vertically partitioned data:

``` math
X = [X_1 | X_2 | ... | X_K]
```
``` math
\beta = [\beta_1^T, \beta_2^T, ..., \beta_K^T]^T
```

The linear predictor decomposes as:

``` math
\eta = X_1\beta_1 + X_2\beta_2 + ... + X_K\beta_K
```

### The BCD Algorithm

Block Coordinate Descent iteratively updates each block while holding
others fixed:

    Initialize: β₁ = β₂ = ... = βₖ = 0

    Repeat until convergence:
        For each institution k = 1, ..., K:

            1. Compute contribution from OTHER institutions:
               η₋ₖ = Σⱼ≠ₖ Xⱼβⱼ  (shared as aggregate)

            2. Institution k solves locally:
               βₖ = argmin L(y; η₋ₖ + Xₖβₖ) + λ||βₖ||²

            3. Compute and share new contribution:
               ηₖ = Xₖβₖ  (shared as aggregate)

        Check convergence: Σₖ||βₖⁿᵉʷ - βₖᵒˡᵈ|| < tolerance

**Privacy guarantee:** Only the linear predictor contributions
$`\eta_k`$ are shared, never the raw data $`X_k`$ or intermediate
coefficient values during iteration.

### IRLS for Non-Gaussian Families

For non-Gaussian families, each local update uses Iteratively Reweighted
Least Squares (IRLS). The update at institution $`k`$ is:

``` math
\beta_k^{new} = (X_k^T W X_k + \lambda I)^{-1} X_k^T W (z - \eta_{-k})
```

where:

- $`W`$ is the weight matrix (depends on family)
- $`z`$ is the working response (depends on family)
- $`\lambda`$ is the L2 regularization parameter
- $`\eta_{-k}`$ is the contribution from other institutions

### Family-Specific Formulas

For link function $`g(\mu) = \eta`$, variance function $`V(\mu)`$, and
response $`y`$:

| Family | Link | Mean $`\mu`$ | Weight $`w`$ | Working Response $`z`$ |
|----|----|----|----|----|
| Gaussian | Identity | $`\eta`$ | $`1`$ | $`y`$ |
| Binomial | Logit | $`\frac{1}{1+e^{-\eta}}`$ | $`\mu(1-\mu)`$ | $`\eta + \frac{y-\mu}{\mu(1-\mu)}`$ |
| Poisson | Log | $`e^\eta`$ | $`\mu`$ | $`\eta + \frac{y-\mu}{\mu}`$ |
| Gamma | Log | $`e^\eta`$ | $`1`$ | $`\eta + \frac{y-\mu}{\mu}`$ |
| Inv. Gaussian | Log | $`e^\eta`$ | $`\frac{1}{\mu}`$ | $`\eta + \frac{y-\mu}{\mu}`$ |

### Convergence Properties

The BCD algorithm for GLMs has the following properties:

- **Monotonic descent**: Each iteration decreases (or maintains) the
  objective
- **Linear convergence**: Converges at a geometric rate for
  well-conditioned problems
- **Regularization**: The $`\lambda`$ parameter ensures numerical
  stability

Typical convergence requires 20-100 iterations depending on data
characteristics.

### Reference

van Kesteren, E.J. et al. (2019). [Privacy-preserving generalized linear
models using distributed block coordinate
descent](https://arxiv.org/abs/1911.03183). *arXiv:1911.03183*.

------------------------------------------------------------------------

## Privacy-Preserving Record Alignment

### The Hashing Approach

To align records across institutions without revealing identifiers:

1.  Each institution computes SHA-256 hashes of its identifiers
2.  Only hashes are shared (irreversible)
3.  Institutions reorder data to match a reference hash order

**SHA-256 properties:**

- One-way function: Cannot recover input from output
- Collision resistant: Extremely unlikely for two inputs to produce same
  hash
- Deterministic: Same input always produces same hash

### Example

| Original ID     | →   | SHA-256 Hash                          |
|-----------------|-----|---------------------------------------|
| `PATIENT_00042` | →   | `a3f8c9d2e1b7...` (64 hex characters) |

Given only the hash, recovering “PATIENT_00042” would require:

- Trying all possible inputs (computationally infeasible)
- Even knowing it’s a patient ID doesn’t significantly reduce the search
  space

------------------------------------------------------------------------

## Model Diagnostics

### Deviance

Deviance measures how well the model fits, compared to a saturated
model:

``` math
D = 2 \times [\log L(\text{saturated}) - \log L(\text{fitted})]
```

For each family:

| Family | Deviance Formula |
|----|----|
| Gaussian | $`\sum_i (y_i - \mu_i)^2`$ |
| Binomial | $`2\sum_i [y_i \log\frac{y_i}{\mu_i} + (1-y_i)\log\frac{1-y_i}{1-\mu_i}]`$ |
| Poisson | $`2\sum_i [y_i \log\frac{y_i}{\mu_i} - (y_i - \mu_i)]`$ |
| Gamma | $`2\sum_i [-\log\frac{y_i}{\mu_i} + \frac{y_i - \mu_i}{\mu_i}]`$ |
| Inv. Gaussian | $`\sum_i \frac{(y_i - \mu_i)^2}{\mu_i^2 y_i}`$ |

### Pseudo R-squared

McFadden’s pseudo R-squared:

``` math
R^2_{pseudo} = 1 - \frac{D_{model}}{D_{null}}
```

where $`D_{null}`$ is the deviance of an intercept-only model.

Interpretation:

- 0: Model no better than intercept only
- 1: Perfect fit (deviance = 0)
- Typical values: 0.2-0.4 considered good for logistic regression

### AIC

Akaike Information Criterion:

``` math
AIC = D + 2k
```

where $`k`$ is the number of parameters. Lower AIC indicates better
model (balancing fit and complexity).

------------------------------------------------------------------------

## Limitations and Considerations

### Semi-honest Security

The current implementation assumes **semi-honest** (honest-but-curious)
parties:

- Institutions follow the protocol correctly
- But may try to infer information from received messages
- More sophisticated attacks (malicious adversaries) are not addressed

### Numerical Considerations

- **Convergence**: BCD may require many iterations for ill-conditioned
  data
- **Regularization**: Small $`\lambda`$ (default 1e-4) prevents singular
  matrices
- **Scaling**: Very large or very small values may cause numerical
  issues

### When to Use Each Method

| Scenario                    | Recommended Approach |
|-----------------------------|----------------------|
| Exploratory analysis        | Correlation → PCA    |
| Continuous outcome          | Gaussian GLM         |
| Binary outcome              | Binomial GLM         |
| Count data                  | Poisson GLM          |
| Positive continuous (costs) | Gamma GLM            |
| High-variance positive data | Inverse Gaussian GLM |

------------------------------------------------------------------------

## Summary

dsVertClient enables privacy-preserving analysis through:

1.  **Block SVD**: Combines partial SVD results to compute global
    correlation/PCA
2.  **Block Coordinate Descent**: Iteratively updates coefficients,
    sharing only linear predictors
3.  **Cryptographic Hashing**: Aligns records without revealing
    identifiers

These methods provide strong privacy guarantees while enabling
meaningful statistical analysis across vertically partitioned data.
