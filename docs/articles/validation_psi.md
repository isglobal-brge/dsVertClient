# PSI alignment validation

## What is validated

Functions: `ds.psiAlign(), ds.isPsiAligned()`.

Private set intersection aligns vertically split tables on a shared
identifier. Each server blinds its identifiers, the protocol computes
the intersection, and the aligned data object keeps only matched rows in
a common order.

## Mathematical target

The statistical target is the set I = intersection_k I_k. The validation
checks that every server has \|I\| rows after alignment and that a
correlation of deterministic id-derived columns is 1 after alignment.

## Fixture and reference

Fixture: Synthetic overlapping patient_id ranges with shuffled row
order.

Centralized reference: Deterministic set intersection in centralized R.

The executable chunk below calls `run_validation()` from
`vignettes/validation_helpers.R`. That helper constructs the fixture,
opens a DSLite server, performs PSI alignment, runs the dsVertClient
product route for K=2 and K=3, computes the centralized reference, and
compares both results. No RDS or result table outside this package is
required; if a local `vignettes/validation-cache/` file exists it is a
cache produced by this same execution path.

## Disclosure review

The analyst receives counts and status fields. Patient identifiers,
matched row indices, and length-n index vectors are not returned; the
legacy index reveal aggregate is checked as blocked.

The fixture keeps `datashield.privacyLevel = 5` and is sized so the
standard disclosure guards remain active. Only
`dsvert.require_trusted_peers` is disabled for DSLite because there is
no real Opal/Rock deployment in this local validation context.

## Executed evidence

``` r

rows <- run_validation("psi", force = force_run)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.psiAlign | synthetic ID intersection | deterministic set intersection | max(count_delta, correlation_delta) | 0 | 0 | strict-precise | PASS | 1.0 |
| K\>=3 | ds.psiAlign | synthetic ID intersection | deterministic set intersection | max(count_delta, correlation_delta) | 0 | 0 | strict-precise | PASS | 0.8 |

## Verdict

The vignette fails during rendering if either K=2 or K\>=3 leaves the
accepted numerical envelope or is marked as disclosive.
