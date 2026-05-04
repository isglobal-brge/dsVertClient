# Federated joint proportional-odds ordinal regression via Ring127 MPC-orchestrated Newton iteration

Strict K=2 share-domain proportional-odds Newton route by default.
Patient-level linear predictors, class probabilities, reciprocals, and
score terms remain Ring127 shares; the client receives only guarded
class counts, aggregate score/Hessian probes, and final model
parameters. The older per-patient reconstruction route has been removed
from the product package.

For a K-level ordered outcome the PO log-likelihood is \$\$\ell(\beta,
\theta) = \sum_i \log\[F(\theta\_{y_i} - \eta_i) - F(\theta\_{y_i-1} -
\eta_i)\],\$\$ with \\F = \mathrm{sigmoid}\\, \\\eta_i = X_i \beta\\,
\\\theta_0 = -\infty\\, \\\theta_K = +\infty\\. Score: \$\$\partial \ell
/ \partial \beta_j = -\sum_i x\_{ij} \cdot
\frac{f(\theta\_{y_i}-\eta_i) - f(\theta\_{y_i-1}-\eta_i)}
{F(\theta\_{y_i}-\eta_i) - F(\theta\_{y_i-1}-\eta_i)},\$\$ with \\f(u) =
F(u)(1-F(u))\\.

MPC pipeline per outer Newton iter:

1.  Compute \\\eta\\ share via `k2ComputeEtaShareDS`.

2.  For each threshold \\k\\: compute \\\theta_k - \eta\\ share via
    affine-combine (plaintext \\\theta_k\\ public).

3.  Apply exp + recip on \\\theta_k - \eta\\ -\> \\F_k\\ share (\\F(u) =
    1/(1+\exp(-u)) = e^u/(1+e^u)\\ – evaluated as
    `exp(u) * (1/(1+exp(u)))` via existing primitives).

4.  \\f_k = F_k (1 - F_k)\\ via Beaver vecmul.

5.  Per-patient residual numerator/denominator built from
    indicator-weighted differences (done on outcome server since it
    holds \\y_i\\ plaintext).

6.  `.ring127_recip_round_keyed` on the F-difference share.

7.  Beaver vecmul (f-diff) \* (1/F-diff) -\> \\T_i\\ share.

8.  Beaver matvec \\X^\top T\\ -\> aggregate score for \\\beta\\.

9.  Client Newton on stacked \\(\beta, \theta)\\ using an aggregate
    finite-difference Hessian over share-domain scores.

## Usage

``` r
ds.vertOrdinalJointNewton(
  formula,
  data = NULL,
  levels_ordered,
  cumulative_template = "%s_leq",
  max_outer = 8L,
  tol = 1e-04,
  verbose = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  Ordered outcome on LHS.

- data:

  Aligned data name.

- levels_ordered:

  Character vector of ordered levels (low -\> high).

- cumulative_template:

  e.g. `"%s_leq"` for Y \<= k indicator.

- max_outer:

  Outer Newton iterations.

- tol:

  Convergence tolerance on \\\\\Delta (\beta, \theta)\\\_\infty\\.

- verbose:

  Logical.

- datasources:

  DataSHIELD connections.

## Value

`ds.vertOrdinalJointNewton` object.
