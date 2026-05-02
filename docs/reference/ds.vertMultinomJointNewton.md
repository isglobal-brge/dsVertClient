# Federated joint-softmax multinomial logistic regression via Ring127 MPC-orchestrated Newton iteration

Full per-patient softmax Newton: orchestrates K-1 parallel exp(eta_k)
shares, sums to denominator D = 1+Sumexp(eta_k), computes 1/D via
Ring127 Chebyshev + Newton-Raphson, multiplies per class to get p_k(x_i)
share per patient, builds residual y_k - p_k on outcome server, and
aggregates X^T(y_k - p_k) via existing Beaver matvec pipeline for each
class. Client-side Bohning-Hessian-bounded Newton step on stacked beta.

All per-patient quantities stay as Ring127 additive shares; only the
final p(K-1)-dim aggregate gradient per iter is revealed – same privacy
class as the single-class gradient of ds.vertGLM. **P3 delta: zero new
reveal types.**

## Usage

``` r
ds.vertMultinomJointNewton(
  formula,
  data = NULL,
  levels,
  indicator_template = "%s_ind",
  max_outer = 8L,
  tol = 1e-04,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  R formula with categorical outcome on LHS.

- data:

  Aligned data frame name.

- levels:

  Character vector of outcome levels (first = reference).

- indicator_template:

  sprintf template for class-indicator columns on the outcome server,
  e.g. `"%s_ind"`. Must exist server-side.

- max_outer:

  Outer Newton iterations (default 8).

- tol:

  Convergence tolerance on max \|Deltabeta\| (default 1e-4).

- verbose:

  Logical.

- datasources:

  DataSHIELD connections.

## Value

`ds.vertMultinomJointNewton` object.
