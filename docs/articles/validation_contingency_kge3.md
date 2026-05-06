# Contingency tests validation (K\>=3)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_pima), collapse = "\n"))
#> function (n = 60L) 
#> {
#>     validation_require("MASS")
#>     pdf <- as.data.frame(MASS::Pima.tr)
#>     pdf$patient_id <- sprintf("P%03d", seq_len(nrow(pdf)))
#>     pdf$diabetes <- as.integer(pdf$type == "Yes")
#>     pdf$type <- NULL
#>     pdf <- stats::na.omit(pdf)
#>     pdf <- pdf[seq_len(min(nrow(pdf), n)), , drop = FALSE]
#>     pdf[, c("patient_id", "npreg", "glu", "bp", "skin", "bmi", 
#>         "ped", "age", "diabetes")]
#> }
```

## Scope

Functions:
[`ds.vert.chisq_cross()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K>=3` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

Cross-server contingency testing builds guarded cell counts across
vertically separated categorical variables before applying the usual
chi-square/Fisher calculations.

## Mathematical target

The validation compares the released observed table and X^2 statistic
with table() and chisq.test() on the pooled fixture.

## Disclosure review

Exact counts are released only after small-cell checks pass. Sparse
positive cells fail closed before release.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 3L

pooled <- build_pima(80L)
pooled$diabetes <- factor(ifelse(pooled$diabetes == 1, "yes", "no"),
                          levels = c("no", "yes"))
pooled$age_grp <- factor(ifelse(pooled$age >= median(pooled$age),
                                "older", "younger"),
                         levels = c("younger", "older"))
pooled$bp_grp <- factor(ifelse(pooled$bp >= median(pooled$bp),
                               "high_bp", "low_bp"),
                        levels = c("low_bp", "high_bp"))

tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "age_grp"), drop = FALSE],
       s2 = pooled[, c("patient_id", "bp_grp", "diabetes"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "age_grp"), drop = FALSE],
       s2 = pooled[, c("patient_id", "bp_grp"), drop = FALSE],
       s3 = pooled[, c("patient_id", "diabetes"), drop = FALSE])
}

knitr::kable(utils::head(pooled[, c("patient_id", "age_grp", "bp_grp",
                                    "diabetes")]))
```

| patient_id | age_grp | bp_grp  | diabetes |
|:-----------|:--------|:--------|:---------|
| P001       | younger | low_bp  | no       |
| P002       | older   | high_bp | yes      |
| P003       | older   | high_bp | no       |
| P004       | younger | high_bp | no       |
| P005       | younger | low_bp  | no       |
| P006       | older   | high_bp | yes      |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns              |
|:----|:-------|----:|:---------------------|
| s1  | s1     |  80 | patient_id, age_grp  |
| s2  | s2     |  80 | patient_id, bp_grp   |
| s3  | s3     |  80 | patient_id, diabetes |

``` r


tab <- table(pooled$age_grp, pooled$diabetes)
chi_ref <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.chisq_cross(
      "DA", "age_grp", "diabetes", correct = FALSE, fisher = TRUE,
      verbose = validation_demo_verbose(), datasources = conns)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})
#> 
#> Logging into the collaborating servers
#> [NA handling] Per-server na.omit (before PSI)...
#> === ECDH-PSI Record Alignment (Blind Relay) ===
#> Reference: s1, Targets: s2, s3
#> [Phase 0] Transport key exchange...
#>   Transport keys exchanged.
#> [Phase 1] Reference server masking IDs...
#>   s1: 80 IDs masked (stored server-side)
#> [Phase 2] s1: exporting encrypted points for s2...
#> [Phase 2] s1: exporting encrypted points for s3...
#> [Phase 3] s2, s3: processing (parallel)...
#>   s2: 80 IDs masked
#>   s3: 80 IDs masked
#> [Phase 4-5] s2: double-masking via s1 (blind relay)...
#> [Phase 6-7] s2: matching and aligning (blind relay)...
#> [Phase 4-5] s3: double-masking via s1 (blind relay)...
#> [Phase 6-7] s3: matching and aligning (blind relay)...
#> [Phase 7] s1: self-aligning...
#> [Phase 8] Computing encrypted multi-server intersection...
#>   Common records: 80
#> Server 's1': 80 of 80 records matched (100.0%)
#> Server 's2': 80 of 80 records matched (100.0%)
#> Server 's3': 80 of 80 records matched (100.0%)
#> PSI alignment complete.
#> [ds.vertChisqCross] age_grp@s1 x diabetes@s3, session=85a6f458
#>   n[older,no] = 21
#>   n[older,yes] = 20
#>   n[younger,no] = 32
#>   n[younger,yes] = 7

cross <- result$value
ref_mat <- unclass(tab)[rownames(cross$observed),
                        colnames(cross$observed), drop = FALSE]
observed <- max(max(abs(cross$observed - ref_mat)),
                abs(cross$chisq - unname(chi_ref$statistic)))

rows <- row_result(
  "contingency", "Contingency tests", K,
  "ds.vert.chisq / ds.vert.fisher / ds.vert.chisq_cross",
  "MASS::Pima.tr categorical fixture", "chisq.test / fisher.test",
  "max(count_or_chisq_delta)", observed, 1e-8, "strict-precise",
  "Releases guarded table counts and test statistics; small cells fail closed.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K\>=3 | ds.vert.chisq / ds.vert.fisher / ds.vert.chisq_cross | MASS::Pima.tr categorical fixture | chisq.test / fisher.test | max(count_or_chisq_delta) | 0 | 0 | strict-precise | PASS | 3.1 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
