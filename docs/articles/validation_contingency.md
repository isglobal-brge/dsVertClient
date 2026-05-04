# Contingency table validation

## What is validated

Functions:
[`ds.vertChisq()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertChisq.md),
[`ds.vertFisher()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertFisher.md),
[`ds.vertChisqCross()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertChisqCross.md)

For an observed table $`O_{ab}`$, Pearson’s statistic is
$`\sum_{ab}(O_{ab}-E_{ab})^2/E_{ab}`$, where $`E_{ab}=O_{a+}O_{+b}/n`$.
Same-server tables release guarded aggregate counts. Cross-server tables
construct one-hot indicators and compute $`\sum_i 1\{A_i=a\}1\{B_i=b\}`$
with Beaver dot products before the exact table is opened.

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact
seed, split and reference object in the cache listed below.

``` r


pooled <- MASS::Pima.tr[seq_len(120), ]
pooled$age_group <- cut(pooled$age, breaks = c(0, 30, Inf))
pooled$glu_group <- cut(pooled$glu, breaks = c(0, 120, Inf))
tab <- table(pooled$age_group, pooled$glu_group)
stats::chisq.test(tab, correct = TRUE)
stats::fisher.test(tab)
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with
[`ds.psiAlign()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.psiAlign.md).

``` r


tables <- list(
  s1 = pooled[c("patient_id", "age_group")],
  s2 = pooled[c("patient_id", "glu_group")]
)
```

``` r


same <- dsVertClient::ds.vertChisq(
  "DA", "age_group", "glu_group", server = "s2", datasources = conns)
exact <- dsVertClient::ds.vertFisher(
  "DA", "age_group", "glu_group", server = "s2", datasources = conns)
cross <- dsVertClient::ds.vertChisqCross(
  "DA", "age_group", "glu_group", fisher = TRUE, datasources = conns)
```

To reproduce the cache from the repository root:

``` r
Rscript scripts/validate_method_tables.R
```

## Disclosure review

Positive cells, row margins, and column margins below the configured
privacy threshold block the same-server release. In the cross-server
route, DCF threshold checks reveal only pass/fail before exact cell
counts are released. No one-hot patient vector is returned to the
analyst.

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K
mode is not marked non-disclosive, is not `PASS`, or exceeds its
accepted tolerance.

``` r

rows <- validation_rows("contingency")
assert_validation(rows)
display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | cache |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|:---|
| K=2 | ds.vertChisq / ds.vertFisher / ds.vertChisqCross | MASS::Pima.tr categorical fixture | chisq.test/fisher.test | count_or_pvalue_delta | 0 | 0 | strict-precise | PASS | tables_dslite_20260504-113158.rds |
| K\>=3 | ds.vertChisq / ds.vertFisher / ds.vertChisqCross | MASS::Pima.tr categorical fixture | chisq.test/fisher.test | count_or_pvalue_delta | 0 | 0 | strict-precise | PASS | tables_dslite_20260504-113158.rds |

## Verdict

Both K=2 and K\>=3 validation rows are inside their accepted numerical
envelope and use the current non-disclosive product route. Any legacy
route mentioned in the package is excluded from this evidence path.
