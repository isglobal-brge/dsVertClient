# LMM validation (K=2)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_lmm), collapse = "\n"))
#> function (ncl = 12L, m = 5L, seed = 101L, b_sd = 0.8, eps_sd = 0.4, 
#>     beta = c(1, 0.5, -0.3, 0.2)) 
#> {
#>     set.seed(seed)
#>     n <- ncl * m
#>     cluster <- rep(seq_len(ncl), each = m)
#>     x1 <- rnorm(n)
#>     x2 <- rnorm(n)
#>     x3 <- rnorm(n)
#>     b <- rnorm(ncl, sd = b_sd)[cluster]
#>     y <- beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3 + 
#>         b + rnorm(n, sd = eps_sd)
#>     data.frame(patient_id = sprintf("L%03d", seq_len(n)), cluster, 
#>         x1, x2, x3, y)
#> }
```

## Scope

Functions:
[`ds.vert.lmm()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

Random-intercept LMM uses guarded cluster metadata and protected
cross-products to fit fixed effects and variance components.

## Mathematical target

The model is y = X beta + b_cluster + epsilon with b ~ N(0,sigma_b^2).
Validation compares fixed effects to lme4::lmer().

## Disclosure review

Cluster membership is used internally at the accepted method tier.
Per-cluster residuals and BLUP vectors are not returned.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 2L

validation_require("lme4")
pooled <- build_lmm()
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x3", "y", "cluster"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "x3", "y", "cluster"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
```

| patient_id | cluster |         x1 |         x2 |         x3 |         y |
|:-----------|--------:|-----------:|-----------:|-----------:|----------:|
| L001       |       1 | -0.3260365 | -0.2598021 | -0.6616053 | 2.0959015 |
| L002       |       1 |  0.5524619 | -1.4111730 | -0.7724177 | 2.4229857 |
| L003       |       1 | -0.6749438 | -0.6413576 | -2.0184735 | 0.8332413 |
| L004       |       1 |  0.2143595 |  0.1124575 | -0.5335854 | 1.6932413 |
| L005       |       1 |  0.3107692 |  0.4226043 |  0.4347283 | 1.5580919 |
| L006       |       2 |  1.1739663 |  0.3868353 | -0.7711673 | 1.3529767 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns                    |
|:----|:-------|----:|:---------------------------|
| s1  | s1     |  60 | patient_id, x1, x2         |
| s2  | s2     |  60 | patient_id, x3, y, cluster |

``` r


ref <- lme4::fixef(lme4::lmer(
  y ~ x1 + x2 + x3 + (1 | cluster), data = pooled, REML = TRUE))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.lmm(
      y ~ x1 + x2 + x3, data = "DA", cluster_col = "cluster",
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
#> [ds.vertLMM] 12 clusters, n_total=60 (sizes: median=5.0, max=5)
#> [ds.vertLMM] non-outcome-server predictors (x1,x2) are absorbed into the intercept for residual SS; ICC estimate is on the outcome-server projection only
#> [LMM] closed-form OK (lambda range 0.000; max |coef| 1.19)
#> [LMM] sigma^2 client-side refit: sigma^2=0.125609 (was 0.136651 MoM); sigma_b^2 kept at 0.269692
#> [LMM] iter 1  sigma^2=0.1256  sigma_b^2=0.2697  rho=0.6822
#> [LMM] closed-form OK (lambda range 0.000; max |coef| 1.19)
#> [LMM] sigma^2 client-side refit: sigma^2=0.1165 (was 0.110601 MoM); sigma_b^2 kept at 0.318845
#> [LMM] iter 2  sigma^2=0.1165  sigma_b^2=0.3188  rho=0.7324
#> [LMM] closed-form OK (lambda range 0.000; max |coef| 1.19)
#> [LMM] sigma^2 client-side refit: sigma^2=0.116364 (was 0.11045 MoM); sigma_b^2 kept at 0.320473
#> [LMM] iter 3  sigma^2=0.1164  sigma_b^2=0.3205  rho=0.7336
#> [LMM] closed-form OK (lambda range 0.000; max |coef| 1.19)
#> [LMM] sigma^2 client-side refit: sigma^2=0.11637 (was 0.110453 MoM); sigma_b^2 kept at 0.320499
#> [LMM] iter 4  sigma^2=0.1164  sigma_b^2=0.3205  rho=0.7336

fit <- result$value
observed <- max_named_delta(fit$coefficients, ref)

rows <- row_result(
  "lmm", "LMM", K, "ds.vert.lmm",
  "synthetic random-intercept fixture", "lme4::lmer",
  "fixed_effect_max_abs_delta", observed, 0.02, "strict-practical",
  "Cluster membership is used internally; per-cluster residuals/BLUPs are not returned.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.lmm | synthetic random-intercept fixture | lme4::lmer | fixed_effect_max_abs_delta | 0.0006448 | 0.02 | strict-practical | PASS | 36.9 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
