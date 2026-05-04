# Federated joint-softmax multinomial logistic regression

Compatibility wrapper for
[`ds.vertMultinomJointNewton`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinomJointNewton.md).
The historical one-vs-rest/covariance-rescale implementation was removed
from the product package because it did not close the softmax MLE
accuracy gap. Calls now always dispatch to the non-disclosive joint
Newton route for K \>= 3 classes.

## Usage

``` r
ds.vertMultinomJoint(
  formula,
  data = NULL,
  levels = NULL,
  max_iter = 30L,
  tol = 1e-04,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  R formula with the categorical outcome on the LHS.

- data:

  Aligned data-frame name.

- levels:

  Optional character vector of outcome levels (first is reference). If
  NULL, inferred from the outcome server.

- max_iter:

  Outer Newton iterations.

- tol:

  Convergence tolerance on max \|delta beta\|.

- verbose:

  Print progress.

- datasources:

  DataSHIELD connections.

## Value

A `ds.vertMultinomJointNewton` object.
