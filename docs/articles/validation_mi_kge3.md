# Multiple imputation validation (K\>=3)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_mi), collapse = "\n"))
#> function (n = 60L, seed = 99L, missing_every = 7L, beta = c(1, 
#>     0.4, -0.2, 0.1), noise_sd = 0.2) 
#> {
#>     set.seed(seed)
#>     x1 <- rnorm(n)
#>     x2 <- rnorm(n)
#>     x3 <- rnorm(n)
#>     y <- beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3 + 
#>         rnorm(n, 0, noise_sd)
#>     d <- data.frame(patient_id = sprintf("MI%03d", seq_len(n)), 
#>         x1, x2, x3, y)
#>     d$x2[seq(5L, n, by = missing_every)] <- NA
#>     d
#> }
```

## Scope

Functions:
[`ds.vert.mi()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K>=3` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

Missing covariates are imputed server-side for each round, then GLM fits
are pooled client-side.

## Mathematical target

The compact reference uses central mean imputation, and the distributed
route is checked by pooled coefficient delta.

## Disclosure review

Imputed columns stay server-side. The client receives pooled
beta/covariance summaries, not imputed patient-level values.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 3L

pooled <- build_mi()
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2", "x3", "y"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "x3", "y"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
```

| patient_id |         x1 |         x2 |         x3 |         y |
|:-----------|-----------:|-----------:|-----------:|----------:|
| MI001      |  0.2139625 | -0.1629179 |  0.5943204 | 1.0253829 |
| MI002      |  0.4796581 | -0.1142152 |  0.3793513 | 0.8646986 |
| MI003      |  0.0878287 | -0.4446594 | -0.9439510 | 0.9415261 |
| MI004      |  0.4438585 |  0.2569592 | -1.6630335 | 0.7977042 |
| MI005      | -0.3628379 |         NA |  0.0289771 | 0.8345062 |
| MI006      |  0.1226740 | -1.3365759 |  0.7939753 | 1.0702498 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns           |
|:----|:-------|----:|:------------------|
| s1  | s1     |  60 | patient_id, x1    |
| s2  | s2     |  60 | patient_id, x2    |
| s3  | s3     |  60 | patient_id, x3, y |

``` r


ref_data <- pooled
ref_data$x2[is.na(ref_data$x2)] <- mean(ref_data$x2, na.rm = TRUE)
ref <- coef(stats::lm(y ~ x1 + x2 + x3, data = ref_data))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns, na.action = "none")
    dsVertClient::ds.vert.mi(
      y ~ x1 + x2 + x3, data = "DA", impute_columns = "x2",
      family = "gaussian",
      verbose = validation_demo_verbose(), datasources = conns, seed = 12L)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})
#> 
#> Logging into the collaborating servers
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
#> [ds.vertMI] M=20 imputations, variables: x2
#> [ds.vertMI] Imputation round 1/20
#> [ds.vertMI] Imputation round 2/20
#> [ds.vertMI] Imputation round 3/20
#> [ds.vertMI] Imputation round 4/20
#> [ds.vertMI] Imputation round 5/20
#> [ds.vertMI] Imputation round 6/20
#> [ds.vertMI] Imputation round 7/20
#> [ds.vertMI] Imputation round 8/20
#> [ds.vertMI] Imputation round 9/20
#> [ds.vertMI] Imputation round 10/20
#> [ds.vertMI] Imputation round 11/20
#> [ds.vertMI] Imputation round 12/20
#> [ds.vertMI] Imputation round 13/20
#> [ds.vertMI] Imputation round 14/20
#> [ds.vertMI] Imputation round 15/20
#> [ds.vertMI] Imputation round 16/20
#> [ds.vertMI] Imputation round 17/20
#> [ds.vertMI] Imputation round 18/20
#> [ds.vertMI] Imputation round 19/20
#> [ds.vertMI] Imputation round 20/20

fit <- result$value
observed <- max_named_delta(fit$coefficients, ref)

rows <- row_result(
  "mi", "Multiple imputation", K, "ds.vert.mi",
  "synthetic missing-covariate fixture",
  "central mean-imputation reference",
  "pooled_coef_abs_delta", observed, 0.02, "strict-practical",
  "Imputed columns stay server-side; client pools beta/covariance only.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K\>=3 | ds.vert.mi | synthetic missing-covariate fixture | central mean-imputation reference | pooled_coef_abs_delta | 5.51e-05 | 0.02 | strict-practical | PASS | 82 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
