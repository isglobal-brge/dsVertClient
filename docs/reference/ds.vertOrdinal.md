# Federated ordinal logistic regression

User-facing ordinal wrapper. Dispatches to
[`ds.vertOrdinalJointNewton`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertOrdinalJointNewton.md),
the paper-safe joint proportional-odds route for K \>= 3. The historical
cumulative-binomial approximation is no longer exposed as a user-facing
estimator; it remains only as an internal warm start for the joint
route.

## Usage

``` r
ds.vertOrdinal(
  formula,
  data = NULL,
  levels_ordered,
  cumulative_template = "%s_leq",
  max_iter = NULL,
  max_outer = 8L,
  tol = NULL,
  warm_max_iter = NULL,
  warm_tol = NULL,
  binomial_sigmoid_intervals = NULL,
  verbose = TRUE,
  datasources = NULL,
  ...
)
```

## Arguments

- formula:

  R formula with the ORDERED outcome on the LHS (passed through as a
  factor level name in the per-threshold formulas).

- data:

  Name of the aligned data frame on each server.

- levels_ordered:

  Character vector of the ordered levels, smallest-to-largest.

- cumulative_template:

  String format (`sprintf`-style) used to build cumulative indicator
  column names; for instance the default `"%s_leq"` produces
  `<level>_leq`, a 0/1 column that is 1 when the patient's outcome is
  at-or-below that level. Columns must already exist server-side.

- max_iter:

  Optional alias for `max_outer`.

- max_outer:

  Maximum outer Newton iterations for the joint route.

- tol:

  Convergence tolerance for the joint route.

- warm_max_iter:

  Optional maximum iterations for each internal binomial warm-start GLM.

- warm_tol:

  Optional tolerance for each internal binomial warm-start GLM.

- binomial_sigmoid_intervals:

  Optional DCF spline interval count for internal binomial warm-start
  GLMs.

- verbose:

  Logical (default TRUE). Print per-threshold fit progress.

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

- ...:

  Reserved for future extensions.

## Value

`ds.vertOrdinal` object with (among other fields): `thresholds`
\\\alpha_k\\ (intercepts of the K-1 cumulative binomial fits) and
`beta_po` \\\gamma\\ (BLUE-pooled slope coefficients from the K-1 fits).
Both are in the CUMULATIVE-BINOMIAL GLM convention, i.e.\\ the fit form
is

      \eqn{P(Y \leq k | X) = \mathrm{sigmoid}(\alpha_k + X^\top \gamma)},

    NOT the \code{MASS::polr} convention
    \eqn{P(Y \leq k | X) = \mathrm{sigmoid}(\theta_k - X^\top \beta)}.
    The two agree under \eqn{\theta_k = \alpha_k} and \eqn{\beta = -\gamma}.
    Therefore a caller comparing against \code{coef(polr)} must flip the
    sign of \code{beta_po} (or equivalently evaluate predictions with
    \eqn{\mathrm{sigmoid}(\theta_k + X^\top \gamma)} on the \code{ds.vertOrdinal}
    outputs). Empirically the cumulative probabilities agree with polr
    to max \eqn{|\Delta P| \approx 5 \times 10^{-2}} on the housing
    subset once the convention is honoured (probe_ordinal_harness.R,
    2026-04-21).
