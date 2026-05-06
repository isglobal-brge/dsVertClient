# PCA validation (K=2)

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
[`ds.vert.pca()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

PCA is computed from the federated correlation matrix, not from
individual-level scores.

## Mathematical target

The target is R = V Lambda V^T. Validation compares eigenvalues and
sign-aligned loadings with eigen(cor(X)).

## Disclosure review

Only eigenvalues and loadings derived from the aggregate correlation
matrix are returned; no individual PCA scores are released.

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
    cor_ds <- dsVertClient::ds.vert.cor(
      "DA", variables = vars,
      verbose = validation_demo_verbose(), datasources = conns)
    pca_ds <- dsVertClient::ds.vert.pca(
      cor_result = cor_ds,
      verbose = validation_demo_verbose(), datasources = conns)
    list(cor_ds = cor_ds, pca_ds = pca_ds)
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
#> Using provided correlation matrix (skipping Ring63 protocol)...
#> Performing PCA via eigen decomposition...
#> PCA complete: 7 components extracted.

cor_ds <- result$value$cor_ds
pca_ds <- result$value$pca_ds
cor_ref <- stats::cor(pooled[cor_ds$var_names])
eig_ref <- eigen(cor_ref, symmetric = TRUE)
load_ref <- eig_ref$vectors[, seq_len(ncol(pca_ds$loadings)), drop = FALSE]
rownames(load_ref) <- rownames(pca_ds$loadings)
load_ds <- pca_ds$loadings
for (j in seq_len(ncol(load_ds))) {
  if (sum(load_ds[, j] * load_ref[, j]) < 0) load_ds[, j] <- -load_ds[, j]
}
observed <- max(max(abs(pca_ds$eigenvalues - eig_ref$values)),
                max(abs(load_ds - load_ref)))

rows <- row_result(
  "pca", "PCA", K, "ds.vert.pca",
  "MASS::Pima.tr fixture", "eigen(cor(X))",
  "max(eigen_loading_abs_delta)", observed, 1e-4,
  "strict-practical",
  "Reuses the correlation matrix; no score or loading per patient is returned.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.pca | MASS::Pima.tr fixture | eigen(cor(X)) | max(eigen_loading_abs_delta) | 3.29e-05 | 1e-04 | strict-practical | PASS | 3.4 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
