# Post-hoc soft-thresholded GLM coefficients (naive LASSO)

Apply the LASSO soft-threshold operator \\\mathrm{sign}(\beta_j) \cdot
\max(\|\beta_j\| - \alpha \lambda_1, 0)\\ to the coefficient vector of a
fitted `ds.glm` object, without re-running the MPC iteration loop. This
is the simplest possible approximation to LASSO: it zeros coefficients
whose magnitude falls below the threshold but does not iteratively
re-optimise under the sparsity constraint.

Useful for quick variable-selection sketches and for checking whether a
subsequent proper proximal-gradient implementation would pay off. Proper
LASSO / elastic net via client-side proximal gradient with per-iteration
MPC gradient calls is a planned Month 3 deliverable (see
V2_PROGRESS.md).

## Usage

``` r
ds.vertLASSO(
  fit,
  lambda_1,
  alpha_grid = c(1, 0.5, 0.25, 0.125, 0.0625),
  keep_intercept = TRUE
)
```

## Arguments

- fit:

  A `ds.glm` object from `ds.vertGLM`.

- lambda_1:

  L1 penalty magnitude.

- alpha_grid:

  Step-size multipliers to sweep (default 1, 0.5, 0.25, 0.125, 0.0625).

- keep_intercept:

  If TRUE (default) never threshold the intercept.

## Value

List of class `ds.vertLASSO` containing the thresholded coefficient
vectors indexed by effective lambda = alpha \* lambda_1.
