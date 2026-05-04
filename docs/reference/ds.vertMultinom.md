# Federated multinomial logistic regression

User-facing multinomial wrapper. Dispatches to
[`ds.vertMultinomJointNewton`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertMultinomJointNewton.md),
the paper-safe joint softmax Newton route for K \>= 3. The historical
one-vs-rest approximation is no longer exposed as a user-facing
estimator; it remains only as an internal warm start for the joint
Newton route.

## Usage

``` r
ds.vertMultinom(
  formula,
  data = NULL,
  classes = NULL,
  reference = NULL,
  indicator_template = "%s_ind",
  max_iter = NULL,
  max_outer = 8L,
  tol = NULL,
  verbose = TRUE,
  datasources = NULL,
  ...
)
```

## Arguments

- formula:

  R formula with the class indicator on the LHS. The class column must
  be a factor with K levels OR a pre-existing set of binary indicator
  columns named `paste0(class_col, "_is_", level_name)` (one per
  non-reference level) on a single server.

- data:

  Name of the aligned data frame on all servers.

- classes:

  Optional character vector specifying which levels to fit (default: all
  non-reference). The reference level is excluded and its probability is
  computed as \\1 - \sum p_k\\ client-side for any subsequent
  prediction.

- reference:

  Optional name of the reference level.

- indicator_template:

  String format with "%s" replaced by each class name to construct
  indicator column names on the server. Default "\\ indicator columns
  must already exist server-side.

- max_iter:

  Optional alias for `max_outer`.

- max_outer:

  Maximum outer Newton iterations for the joint route.

- tol:

  Convergence tolerance for the joint route.

- verbose:

  Logical (default TRUE). Print per-class fit progress.

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

- ...:

  Reserved for future extensions.

## Value

ds.vertMultinom object: a list with per-class `ds.glm` fits, the level
vector, the reference, and a consolidated coefficient matrix (rows =
coefficients, columns = non-reference classes).
