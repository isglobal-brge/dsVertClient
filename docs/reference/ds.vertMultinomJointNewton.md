# Federated joint-softmax multinomial logistic regression via Ring127 MPC-orchestrated Newton iteration

Full softmax Newton path for vertical splits: orchestrates
class-specific exp(eta) shares, shared softmax denominators, shared
residuals, and Beaver matvec score aggregation. K=2 uses both servers as
DCF parties; K\>=3 selects the outcome server plus one fusion DCF party
and has the other servers contribute encrypted additive shares. The
client performs a Bohning-Hessian-bounded Newton step on stacked
coefficients using only aggregate gradients and low-dimensional
Gram/Hessian objects.

All per-patient probabilities and residuals remain Ring127 additive
shares. The raw design Gram is built from scalar local moments and the
federated correlation matrix, which is the same aggregate-disclosure
tier as `ds.vertCor`.

## Usage

``` r
ds.vertMultinomJointNewton(
  formula,
  data = NULL,
  levels,
  indicator_template = "%s_ind",
  max_outer = 8L,
  tol = 1e-04,
  warm_max_iter = NULL,
  warm_tol = NULL,
  binomial_sigmoid_intervals = NULL,
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

- warm_max_iter:

  Optional maximum iterations for each internal binomial warm-start GLM.

- warm_tol:

  Optional tolerance for each internal binomial warm-start GLM.

- binomial_sigmoid_intervals:

  Optional DCF spline interval count for internal binomial warm-start
  GLMs.

- verbose:

  Logical.

- datasources:

  DataSHIELD connections.

## Value

`ds.vertMultinomJointNewton` object.
