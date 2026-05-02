# Federated ordinal logistic regression (naive cumulative binomials)

Fit a K-category ordinal logistic regression by training K-1 cumulative
binomial logistics of the form \\P(Y \leq k \| X) = \sigma(\theta_k -
X^\top \beta)\\, one for each threshold level. This is the NAIVE
approach: each cumulative binomial is fit independently, giving
per-level beta estimates that may differ across thresholds. Proper
proportional-odds fitting with a single shared beta and K-1 threshold
parameters jointly optimised requires a dedicated joint L-BFGS objective
(Month 2 follow-on).

The naive version is useful for exploratory analyses and for testing the
proportional-odds assumption by comparing the beta estimates across
threshold levels.

## Usage

``` r
ds.vertOrdinal(
  formula,
  data = NULL,
  levels_ordered,
  cumulative_template = "%s_leq",
  ring = 63L,
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

- ring:

  Integer (63 or 127). Ring selector for the underlying binomial
  sub-fits. Currently defaults to 63L while the Ring127 binomial
  wide-spline path is wired up.

- verbose:

  Logical (default TRUE). Print per-threshold fit progress.

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

- ...:

  Passed to each underlying `ds.vertGLM` call.

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
