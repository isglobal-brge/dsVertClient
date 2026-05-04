# LMM validation

## What is validated

Functions: `ds.vertLMM(), ds.vertLMM.k3()`.

Random-intercept LMM uses guarded cluster metadata and MPC
cross-products to fit fixed effects and variance components.

## Mathematical target

The model is y = X beta + b_cluster + epsilon with b ~ N(0, sigma_b^2).
The validation compares fixed effects to lme4::lmer().

## Fixture and reference

Fixture: Synthetic balanced random-intercept fixture, with cluster sizes
above privacy thresholds.

Centralized reference: lme4::lmer(… + (1 \| cluster)) on the pooled
fixture.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

Cluster membership is used internally at the accepted LMM tier.
Per-cluster residuals and BLUP vectors are not returned.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("lmm", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vertLMM | synthetic random-intercept fixture | lme4::lmer | fixed_effect_max_abs_delta | 0.0006365 | 0.02 | strict-practical | PASS | 26.7 |
| K\>=3 | ds.vertLMM.k3 | synthetic random-intercept fixture | lme4::lmer | fixed_effect_max_abs_delta | 0.0001614 | 0.02 | strict-practical | PASS | 169.5 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
