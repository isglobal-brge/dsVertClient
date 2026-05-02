# Federated Cox proportional-hazards via Allison-1982 Poisson trick (K=3)

Discrete-time Allison-1982 / Whitehead-1980 / Prentice-Gloeckler 1978
Sec.2 equivalence: the partial-likelihood Cox estimator is the
maximum-likelihood Poisson regression on person-time-expanded data with
one row per (subject, at-risk-interval) pair, log-link, and a
piecewise-constant baseline hazard \\\alpha_t\\ as a factor covariate.
The slopes \\\beta\\ of the Poisson fit are the Cox PH slopes
(asymptotically).

Implementation reduces to a single
[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md)
call on the K=3 secure-aggregation path with `family="poisson"`,
`offset=log(at_risk * width)`, and the baseline factor `factor(t_bin)`
included as design columns. No new MPC primitives – full inheritance of
the K=3 GLM threat model (D-INV-1..5 preserved; honest-majority
secure-agg envelope; OT-Beaver dishonest-majority at the DCF-pair
level).

The expanded design is expected to live on the cluster under a
user-supplied symbol (default "DA_pt"), constructed by the caller prior
to invocation. The validation probe in
`scripts/run_opal_demo_probes.R::probe_cox_k3` pre-expands the lung K=3
split client-side and uploads `lung_pt_s1/s2/s3` to the test cluster;
production deployments would do the expansion server-side via a non-
disclosive helper (see `coxDiscreteShareDS.R` for the K=2 share-mask
precedent – the K=3 generalisation is straightforward 3-party additive
sharing).

## Usage

``` r
ds.vertCox.k3(
  formula,
  data,
  event_col,
  offset_col,
  baseline_col,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  One-sided formula listing the slope covariates only – *not* a
  Surv(...) formula. Example: `~ age + sex_num + bmi`. The outcome
  column is taken to be the user-supplied `event_col`.

- data:

  Aligned data-frame name on each server (e.g. "DA_pt").

- event_col:

  Outcome column on the outcome server holding the per-row 0/1 event
  indicator (1 = event occurred in this person-time row).

- offset_col:

  Offset column on the outcome server holding \\\log(\Delta t)\\ where
  \\\Delta t\\ is the at-risk width of the row. For not-at-risk rows the
  row is excluded prior to upload (or a weights column is set to 0).

- baseline_col:

  Column name on the outcome server holding the discrete-time bin index
  (e.g. `"t_bin"`). Will be included via `factor(t_bin)` in the design.

- verbose:

  Print progress.

- datasources:

  DataSHIELD K=3 connections.

## Value

list of class `ds.vertCox.k3` with `coefficients` (slopes only –
baseline factor levels dropped), `n_pp` (number of person-time rows),
and the underlying `ds.glm` object as `$fit`.

## References

Allison, P. D. (1982). Discrete-time methods for the analysis of event
histories. *Sociological Methodology*, 13, 61-98. Prentice, R. L. &
Gloeckler, L. A. (1978). Regression analysis of grouped survival data
with application to breast cancer data. *Biometrics*, 34(1), 57-67.
Whitehead, J. (1980). Fitting Cox's regression model to survival data
using GLIM. *Applied Statistics*, 29(3), 268-275. Andreux, M. *et al.*
(2020). Federated survival analysis with discrete-time Cox models.
arXiv:2006.08997.

## See also

[`ds.vertGLM`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertGLM.md),
[`ds.vertCox`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertCox.md)
(K=2 masked-Newton path).
