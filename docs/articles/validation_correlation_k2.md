# Correlation validation (K=2)

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
[`ds.vert.cor()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

The method combines server-local moments and MPC cross-products to
release a guarded low-dimensional correlation matrix.

## Mathematical target

For variables x_j and x_l, r_jl = cov(x_j, x_l) / (sd(x_j) sd(x_l)).
Validation compares the full matrix to stats::cor().

## Disclosure review

The disclosure surface is the aggregate p by p correlation matrix, with
no row-level cross-products returned.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 2L

pooled <- build_pima(60L)
tables <- if (K == 2L) {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi", "ped"), drop = FALSE],
    s2 = pooled[, c("patient_id", "npreg", "glu", "bp", "skin"),
                drop = FALSE])
} else {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi"), drop = FALSE],
    s2 = pooled[, c("patient_id", "ped", "skin"), drop = FALSE],
    s3 = pooled[, c("patient_id", "npreg", "glu", "bp"), drop = FALSE])
}
vars <- lapply(tables, function(x) setdiff(names(x), "patient_id"))

knitr::kable(utils::head(pooled))
```

| patient_id | npreg | glu |  bp | skin |  bmi |   ped | age | diabetes |
|:-----------|------:|----:|----:|-----:|-----:|------:|----:|---------:|
| P001       |     5 |  86 |  68 |   28 | 30.2 | 0.364 |  24 |        0 |
| P002       |     7 | 195 |  70 |   33 | 25.1 | 0.163 |  55 |        1 |
| P003       |     5 |  77 |  82 |   41 | 35.8 | 0.156 |  35 |        0 |
| P004       |     0 | 165 |  76 |   43 | 47.9 | 0.259 |  26 |        0 |
| P005       |     0 | 107 |  60 |   25 | 26.4 | 0.133 |  23 |        0 |
| P006       |     5 |  97 |  76 |   27 | 35.6 | 0.378 |  52 |        1 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns                          |
|:----|:-------|----:|:---------------------------------|
| s1  | s1     |  60 | patient_id, age, bmi, ped        |
| s2  | s2     |  60 | patient_id, npreg, glu, bp, skin |

``` r


result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.cor(
      "DA", variables = vars,
      verbose = validation_demo_verbose(), datasources = conns)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})
#> 
#> Logging into the collaborating servers
#> [NA handling] Per-server na.omit (before PSI)...
#> === ECDH-PSI Record Alignment (Blind Relay) ===
#> Reference: s1, Targets: s2
#> [Phase 0] Transport key exchange...
#>   Transport keys exchanged.
#> [Phase 1] Reference server masking IDs...
#>   s1: 60 IDs masked (stored server-side)
#> [Phase 2] s1: exporting encrypted points for s2...
#> [Phase 3] s2: processing (parallel)...
#>   s2: 60 IDs masked
#> [Phase 4-5] s2: double-masking via s1 (blind relay)...
#> [Phase 6-7] s2: matching and aligning (blind relay)...
#> [Phase 7] s1: self-aligning...
#> [Phase 8] Computing encrypted multi-server intersection...
#>   Common records: 60
#> Server 's1': 60 of 60 records matched (100.0%)
#> Server 's2': 60 of 60 records matched (100.0%)
#> PSI alignment complete.
#> === Ring63 Correlation (p=7, 2 servers) ===
#> [Phase 0] Transport keys...
#> [Phase 1] Standardizing...
#> [Phase 2] Local correlations...
#> [Phase 3] Input sharing...
#> [Phase 4] Beaver cross-products (7 columns)...
#>   Column 1/7 done
#>   Column 2/7 done
#>   Column 3/7 done
#>   Column 4/7 done
#>   Column 5/7 done
#>   Column 6/7 done
#>   Column 7/7 done
#> Correlation complete (2s)

cor_ds <- result$value
cor_ref <- stats::cor(pooled[cor_ds$var_names])
observed <- max(abs(cor_ds$correlation - cor_ref))

rows <- row_result(
  "correlation", "Correlation", K, "ds.vert.cor",
  "MASS::Pima.tr fixture", "stats::cor",
  "correlation_max_abs_delta", observed, 1e-4, "strict-practical",
  "Releases the low-dimensional correlation matrix only.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.cor | MASS::Pima.tr fixture | stats::cor | correlation_max_abs_delta | 1.24e-05 | 1e-04 | strict-practical | PASS | 2.9 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
