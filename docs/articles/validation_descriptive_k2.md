# Descriptive statistics validation (K=2)

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
[`ds.vert.desc()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

Each server computes guarded scalar summaries for its own variables.
Histogram quantiles are guarded by disclosure thresholds, but this
vignette validates the exact mean/sd surface used in the matrix.

## Mathematical target

For each variable x, the checked target is max(\|mean_DS - mean_local\|,
\|sd_DS - sd_local\|).

## Disclosure review

Only scalar summaries and guarded histogram information are returned; no
patient-level values or order statistics are disclosed.

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
all_vars <- unlist(vars, use.names = FALSE)

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


ref_mean <- vapply(pooled[all_vars], mean, numeric(1), na.rm = TRUE)
ref_sd <- vapply(pooled[all_vars], stats::sd, numeric(1), na.rm = TRUE)

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.desc(
      "DA", variables = vars, n_buckets = 40L,
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
#> [ds.vertDesc] s1:age
#> [ds.vertDesc] s1:bmi
#> [ds.vertDesc] s1:ped
#> [ds.vertDesc] s2:npreg
#> [ds.vertDesc] s2:glu
#> [ds.vertDesc] s2:bp
#> [ds.vertDesc] s2:skin

desc <- result$value
ds_mean <- stats::setNames(desc$mean, desc$variable)
ds_sd <- stats::setNames(desc$sd, desc$variable)
observed <- max(max_named_delta(ds_mean, ref_mean),
                max_named_delta(ds_sd, ref_sd))

rows <- row_result(
  "descriptive", "Descriptive statistics", K, "ds.vert.desc",
  "MASS::Pima.tr fixture", "central mean/sd",
  "max(mean_sd_abs_delta)", observed, 1e-8, "strict-precise",
  "Returns guarded scalar summaries and histogram-based quantiles, not rows.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.desc | MASS::Pima.tr fixture | central mean/sd | max(mean_sd_abs_delta) | 0 | 0 | strict-precise | PASS | 1.8 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
