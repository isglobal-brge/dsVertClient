# Ordinal validation (K\>=3)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_ordinal), collapse = "\n"))
#> function (n = 60L, seed = 55L) 
#> {
#>     set.seed(seed)
#>     d <- data.frame(patient_id = sprintf("O%03d", seq_len(n)), 
#>         x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
#>     y <- factor(rep(c("low", "med", "high"), each = n/3L), levels = c("low", 
#>         "med", "high"), ordered = TRUE)
#>     d$y_ord <- y
#>     for (k in levels(y)) {
#>         d[[paste0(k, "_leq")]] <- as.integer(as.integer(y) <= 
#>             match(k, levels(y)))
#>     }
#>     d
#> }
```

## Scope

Functions:
[`ds.vert.ordinal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K>=3` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

The supported route is joint proportional-odds Newton. Historical
cumulative-binomial approximations are only warm starts.

## Mathematical target

P(Y\<=k\|X)=sigmoid(theta_k-X beta). Validation compares cumulative
probabilities with MASS::polr().

## Disclosure review

Cumulative probabilities, residuals, and row scores stay in shares; the
client receives thresholds, slopes, and scalar optimizer diagnostics.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 3L

validation_require("MASS")
old_fd <- getOption("dsvert.ord_strict_fd_max_dim", NULL)
options(dsvert.ord_strict_fd_max_dim = 0L)

pooled <- build_ordinal()
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x3", "y_ord",
                       "low_leq", "med_leq", "high_leq"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "x3", "y_ord",
                       "low_leq", "med_leq", "high_leq"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
```

| patient_id |         x1 |         x2 |         x3 | y_ord | low_leq | med_leq | high_leq |
|:-----------|-----------:|-----------:|-----------:|:------|--------:|--------:|---------:|
| O001       |  0.1201391 |  1.5654063 | -0.1466108 | low   |       1 |       1 |        1 |
| O002       | -1.8123769 |  0.0529084 | -0.4223145 | low   |       1 |       1 |        1 |
| O003       |  0.1515830 |  1.4382989 |  0.5396566 | low   |       1 |       1 |        1 |
| O004       | -1.1192210 |  0.5064401 |  0.3629795 | low   |       1 |       1 |        1 |
| O005       |  0.0019082 | -0.7341779 | -1.8250143 | low   |       1 |       1 |        1 |
| O006       |  1.1885185 | -0.0462339 | -1.0021388 | low   |       1 |       1 |        1 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns                                           |
|:----|:-------|----:|:--------------------------------------------------|
| s1  | s1     |  60 | patient_id, x1                                    |
| s2  | s2     |  60 | patient_id, x2                                    |
| s3  | s3     |  60 | patient_id, x3, y_ord, low_leq, med_leq, high_leq |

``` r


ref <- MASS::polr(y_ord ~ x1 + x2 + x3, data = pooled, Hess = TRUE)
ref_cum <- sapply(ref$zeta, function(th) {
  stats::plogis(th - as.numeric(as.matrix(pooled[, names(coef(ref))]) %*%
                                  coef(ref)))
})

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.ordinal(
      y_ord ~ x1 + x2 + x3, data = "DA",
      levels_ordered = c("low", "med", "high"),
      cumulative_template = "%s_leq",
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
#>   s1: 60 IDs masked (stored server-side)
#> [Phase 2] s1: exporting encrypted points for s2...
#> [Phase 2] s1: exporting encrypted points for s3...
#> [Phase 3] s2, s3: processing (parallel)...
#>   s2: 60 IDs masked
#>   s3: 60 IDs masked
#> [Phase 4-5] s2: double-masking via s1 (blind relay)...
#> [Phase 6-7] s2: matching and aligning (blind relay)...
#> [Phase 4-5] s3: double-masking via s1 (blind relay)...
#> [Phase 6-7] s3: matching and aligning (blind relay)...
#> [Phase 7] s1: self-aligning...
#> [Phase 8] Computing encrypted multi-server intersection...
#>   Common records: 60
#> Server 's1': 60 of 60 records matched (100.0%)
#> Server 's2': 60 of 60 records matched (100.0%)
#> Server 's3': 60 of 60 records matched (100.0%)
#> PSI alignment complete.
#> [ds.vertOrdinal] dispatching to ds.vertOrdinalJointNewton
#> [OrdinalJointStrict] session ordinalStrict_1778061337_1206657760 n=60 p=3 class_counts=[20,20,20]
#> [OrdinalJointStrict eval 1] |g|_inf=6.964e-03
#> [OrdinalJointStrict eval 2] |g|_inf=3.983e-03
#> [OrdinalJointStrict eval 3] |g|_inf=1.422e-03
#> [OrdinalJointStrict eval 4] |g|_inf=3.550e-04
#> [OrdinalJointStrict eval 5] |g|_inf=8.358e-04
#> [OrdinalJointStrict eval 6] |g|_inf=2.427e-04
#> [OrdinalJointStrict eval 7] |g|_inf=2.000e-05
#> [OrdinalJointStrict eval 8] |g|_inf=1.633e-05

fit <- result$value
ds_cum <- ordinal_cumprob(fit, pooled)
colnames(ds_cum) <- colnames(ref_cum)
observed <- max(abs(ds_cum - ref_cum))

rows <- row_result(
  "ordinal", "Ordinal", K, "ds.vert.ordinal",
  "balanced synthetic 3-level ordinal fixture",
  "MASS::polr cumulative probabilities",
  "cumulative_probability_max_abs_delta", observed, 0.001,
  "strict-practical",
  "Class probabilities/residuals stay Ring127 shares; no row probabilities are returned.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K\>=3 | ds.vert.ordinal | balanced synthetic 3-level ordinal fixture | MASS::polr cumulative probabilities | cumulative_probability_max_abs_delta | 3.98e-05 | 0.001 | strict-practical | PASS | 674 |

``` r

assert_validation(rows)
options(dsvert.ord_strict_fd_max_dim = old_fd)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
