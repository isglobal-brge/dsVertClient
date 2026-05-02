# NB regression with Method-of-Moments theta-estimator (K=2-safe)

Negative-binomial GLM with Anscombe 1950 / Saha-Paul 2005
Method-of-Moments dispersion estimator. Replaces the digamma- based MLE
theta-Newton
([`ds.vertNB`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNB.md))
with the closed-form sample-moment estimator
\$\$\hat\theta\_{\mathrm{MoM}} = \bar y^2 / (s^2 - \bar y)\$\$ where
\\s^2 = (\sum y_i^2 - n\bar y^2)/(n-1)\\ is the bias-corrected sample
variance. Under the iid-mu approximation (mu == ybar) this is the
reduction of the regression Saha-Paul 2005 Sec.3 moment equation, and is
psi-free – so it has ZERO new MPC primitive cost beyond the four
y-aggregates already in the existing iid-mu disclosure budget.

## Usage

``` r
ds.vertNBMoMTheta(
  formula,
  data = NULL,
  verbose = TRUE,
  datasources = NULL,
  ...
)
```

## Arguments

- formula:

  Model formula for the count outcome (LHS) and linear predictor (RHS).

- data:

  Aligned data-frame name on each server.

- verbose:

  Logical. Print stage-by-stage progress.

- datasources:

  DataSHIELD connections; if NULL, uses
  [`DSI::datashield.connections_find()`](https://datashield.github.io/DSI/reference/datashield.connections_find.html).

- ...:

  Extra arguments forwarded to `ds.vertGLM`.

## Value

Object of class `c("ds.vertNBMoMTheta", "ds.vertNB", "ds.glm")`
compatible with
[`ds.vertNB`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNB.md)
consumers. `$theta` carries the MoM estimate; `$theta_method = "mom"`;
`$theta_mom_underdispersed` flags pathological s^2 \<= ybar cases.

## Details

Disclosure-equivalent to
[`ds.vertNB`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNB.md):
only y aggregates (Sumy, Sumy^2, n, ybar) are revealed at the
coordinator, which is identical to the `dsvertNBProfileSumsDS` cost. No
per-patient mu disclosure (which would require eta-reveal at OS,
currently outside the K=2-safe budget). Trade-off vs MLE: ~5-10%
asymptotic efficiency loss for moderate overdispersion (Lloyd-Smith 2007
PLoS ONE 2(2):e180); but consistent under correct model specification
(Saha-Paul 2005 Sec.3 Theorem 1) and computable in closed form (no
Newton iteration on theta).

Honesty note: under the iid-mu approximation, this estimator propagates
the same heteroscedasticity bias as iid-mu MLE-theta; the structural
full-regression MoM (Saha-Paul Method 2 with per-patient mu) requires
eta at OS plaintext, which is the same disclosure pattern that was
disabled in ord_joint K=2-safe. Future work: share-space mu_i Beaver
vecmul to recover the full-regression form without eta-reveal.

References:

- Anscombe 1950 *Biometrika* 37(3-4):358-382 (NB moment estimators).

- Saha & Paul 2005 *Biometrics* 61(1):179-185 Sec.3 (bias-corrected
  regression MoM).

- Lloyd-Smith 2007 *PLoS ONE* 2(2):e180 (MLE-vs-MoM efficiency
  comparison for overdispersed counts).

## See also

[`ds.vertNB`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNB.md),
[`ds.vertNBFullRegTheta`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertNBFullRegTheta.md)
