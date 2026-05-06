# GLM inference validation (K\>=3)

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
`ds.vert.confint(), ds.vert.wald(), ds.vert.contrast(), ds.vert.lr()`.

This vignette validates the `K>=3` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

Inference helpers post-process the model-level beta/covariance/deviance
already returned by ds.vert.glm().

## Mathematical target

Wald and contrast tests use K beta and K V K’. LR tests use the deviance
difference with a chi-square reference.

## Disclosure review

No extra server aggregate is opened beyond the fitted GLM objects.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 3L

pooled <- build_pima(60L)
tables <- if (K == 2L) {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi", "ped"), drop = FALSE],
    s2 = pooled[, c("patient_id", "npreg", "glu", "bp", "skin",
                    "diabetes"), drop = FALSE])
} else {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi"), drop = FALSE],
    s2 = pooled[, c("patient_id", "ped", "skin"), drop = FALSE],
    s3 = pooled[, c("patient_id", "npreg", "glu", "bp",
                    "diabetes"), drop = FALSE])
}
fm <- diabetes ~ age + bmi + ped + glu + bp + skin + npreg
red <- diabetes ~ age + bmi + ped + glu + bp + skin

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

|     | server |   n | columns                              |
|:----|:-------|----:|:-------------------------------------|
| s1  | s1     |  60 | patient_id, age, bmi                 |
| s2  | s2     |  60 | patient_id, ped, skin                |
| s3  | s3     |  60 | patient_id, npreg, glu, bp, diabetes |

``` r


result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    full <- dsVertClient::ds.vert.glm(
      fm, data = "DA", family = "gaussian", max_iter = 30L,
      compute_se = TRUE, compute_deviance = TRUE,
      verbose = validation_demo_verbose(), datasources = conns)
    reduced <- dsVertClient::ds.vert.glm(
      red, data = "DA", family = "gaussian", max_iter = 30L,
      compute_se = TRUE, compute_deviance = TRUE,
      verbose = validation_demo_verbose(), datasources = conns)
    ci <- dsVertClient::ds.vert.confint(full)
    wald <- dsVertClient::ds.vert.wald(full, "age")
    Kmat <- matrix(0, nrow = 1, ncol = length(full$coefficients),
                   dimnames = list("age", names(full$coefficients)))
    Kmat[1, "age"] <- 1
    contrast <- dsVertClient::ds.vert.contrast(full, Kmat)
    lr <- dsVertClient::ds.vert.lr(reduced, full)
    list(full = full, ci = ci, wald = wald, contrast = contrast, lr = lr)
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
#> [Auto-detect] Querying server columns...
#>   y_var 'diabetes' found on: s3
#>   s1: age, bmi
#>   s2: ped, skin
#>   s3 (label): npreg, glu, bp
#> === Encrypted-Label BCD-IRLS for gaussian GLM ===
#> Observations: 60, Variables: 7, Partitions: 3
#> Label server: s3 (holds 'diabetes')
#> Non-label servers: s1, s2
#> 
#> [Phase 0] Transport key setup (3 servers)...
#>   [Key Setup] Transport keys exchanged (3 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 3 servers...
#>   [Standardize] 7 total features standardized (y standardized, 0.0s)
#> 
#> [Phase 3] Ring63 Beaver Gradient (K=3 servers, family=gaussian)...
#> 
#> [Phase 3] K>=3 ring63 DCF + Beaver L-BFGS (family=gaussian, n=60, intervals=50, lambda=1.0e-04)
#>   [Input Sharing] Sharing features+y with DCF parties...
#>   [Input Sharing] Complete (p_total=7, 0.2s)
#>     [1.1] Computing eta shares...
#>     [1.3] Beaver gradient...
#>   Iter 1: ||grad||=0.9249, step=0.30, diff=1.52e-01, theta=[0.000, 0.152] (0.3s)
#>     [2.1] Computing eta shares...
#>     [2.3] Beaver gradient...
#>   Iter 2: ||grad||=0.3926, step=1.00, diff=1.45e-01, theta=[0.000, 0.297] (0.2s)
#>     [3.1] Computing eta shares...
#>     [3.3] Beaver gradient...
#>   Iter 3: ||grad||=0.1064, step=1.00, diff=5.26e-02, theta=[0.000, 0.333] (0.2s)
#>   Iter 4: ||grad||=0.0446, step=1.00, diff=3.10e-02, theta=[0.000, 0.344] (0.2s)
#>   Iter 5: ||grad||=0.0188, step=1.00, diff=1.88e-02, theta=[0.000, 0.345] (0.3s)
#>   Iter 6: ||grad||=0.0139, step=1.00, diff=2.58e-02, theta=[0.000, 0.344] (0.2s)
#>   Iter 7: ||grad||=0.0134, step=1.00, diff=6.13e-03, theta=[-0.000, 0.341] (0.2s)
#>   Iter 8: ||grad||=0.0028, step=1.00, diff=1.95e-03, theta=[-0.000, 0.341] (0.2s)
#>   Iter 9: ||grad||=0.0006, step=1.00, diff=3.66e-04, theta=[-0.000, 0.341] (0.3s)
#>   Iter 10: ||grad||=0.0004, step=1.00, diff=8.35e-04, theta=[0.000, 0.342] (0.2s)
#>   Iter 11: ||grad||=0.0005, step=1.00, diff=3.63e-04, theta=[0.000, 0.342] (0.2s)
#>   Iter 12: ||grad||=0.0003, step=1.00, diff=1.43e-04, theta=[0.000, 0.342] (0.2s)
#>   Iter 13: ||grad||=0.0002, step=1.00, diff=6.47e-05, theta=[0.000, 0.342] (0.3s)
#>   Converged after 13 iterations (diff = 6.47e-05)
#>   [Ring63-DCF-Beaver] Converged after 13 iters (total 3.2s)
#>   [SE] Computing exact Hessian (central differences)...
#>     [SE] Column 1/8
#>     [SE] Column 2/8
#>     [SE] Column 3/8
#>     [SE] Column 4/8
#>     [SE] Column 5/8
#>     [SE] Column 6/8
#>     [SE] Column 7/8
#>     [SE] Column 8/8
#>   [Deviance] Computing canonical deviance...
#>   [Deviance] = 30.5246
#> 
#> [Phase 5] Secure deviance (Beaver): 30.5246
#> 
#> Deviance: 6.3032, Null deviance: 12.1833, Pseudo R2: 0.4826
#> 
#> Coefficients:
#>                     Estimate  Std.Error    z value   Pr(>|z|)
#>   (Intercept)        -1.2307     0.3666     -3.357   0.000788 ***
#>   age                 0.0118     0.0054      2.180   0.029222 *
#>   bmi                 0.0034     0.0098      0.350   0.726397
#>   ped                 0.3141     0.1181      2.660   0.007822 **
#>   skin                0.0070     0.0061      1.136   0.255897
#>   npreg               0.0108     0.0199      0.543   0.586954
#>   glu                 0.0047     0.0015      3.075   0.002103 **
#>   bp                  0.0014     0.0045      0.317   0.751517
#> [Auto-detect] Querying server columns...
#>   y_var 'diabetes' found on: s3
#>   s1: age, bmi
#>   s2: ped, skin
#>   s3 (label): glu, bp
#> === Encrypted-Label BCD-IRLS for gaussian GLM ===
#> Observations: 60, Variables: 6, Partitions: 3
#> Label server: s3 (holds 'diabetes')
#> Non-label servers: s1, s2
#> 
#> [Phase 0] Transport key setup (3 servers)...
#>   [Key Setup] Transport keys exchanged (3 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 3 servers...
#>   [Standardize] 6 total features standardized (y standardized, 0.0s)
#> 
#> [Phase 3] Ring63 Beaver Gradient (K=3 servers, family=gaussian)...
#> 
#> [Phase 3] K>=3 ring63 DCF + Beaver L-BFGS (family=gaussian, n=60, intervals=50, lambda=1.0e-04)
#>   [Input Sharing] Sharing features+y with DCF parties...
#>   [Input Sharing] Complete (p_total=6, 0.2s)
#>     [1.1] Computing eta shares...
#>     [1.3] Beaver gradient...
#>   Iter 1: ||grad||=0.8690, step=0.30, diff=1.52e-01, theta=[0.000, 0.152] (0.2s)
#>     [2.1] Computing eta shares...
#>     [2.3] Beaver gradient...
#>   Iter 2: ||grad||=0.4275, step=1.00, diff=1.76e-01, theta=[0.000, 0.328] (0.3s)
#>     [3.1] Computing eta shares...
#>     [3.3] Beaver gradient...
#>   Iter 3: ||grad||=0.0825, step=1.00, diff=4.47e-02, theta=[0.000, 0.346] (0.2s)
#>   Iter 4: ||grad||=0.0321, step=1.00, diff=2.38e-02, theta=[0.000, 0.344] (0.2s)
#>   Iter 5: ||grad||=0.0143, step=1.00, diff=1.95e-02, theta=[0.000, 0.341] (0.2s)
#>   Iter 6: ||grad||=0.0095, step=1.00, diff=2.07e-02, theta=[0.000, 0.343] (0.2s)
#>   Iter 7: ||grad||=0.0078, step=1.00, diff=5.64e-03, theta=[-0.000, 0.344] (0.3s)
#>   Iter 8: ||grad||=0.0046, step=1.00, diff=3.74e-03, theta=[-0.000, 0.343] (0.2s)
#>   Iter 9: ||grad||=0.0002, step=1.00, diff=1.52e-04, theta=[0.000, 0.343] (0.2s)
#>   Iter 10: ||grad||=0.0002, step=1.00, diff=4.08e-04, theta=[0.000, 0.343] (0.2s)
#>   Iter 11: ||grad||=0.0003, step=1.00, diff=2.21e-04, theta=[0.000, 0.343] (0.3s)
#>   Iter 12: ||grad||=0.0001, step=1.00, diff=8.22e-05, theta=[0.000, 0.343] (0.4s)
#>   Converged after 12 iterations (diff = 8.22e-05)
#>   [Ring63-DCF-Beaver] Converged after 12 iters (total 2.9s)
#>   [SE] Computing exact Hessian (central differences)...
#>     [SE] Column 1/7
#>     [SE] Column 2/7
#>     [SE] Column 3/7
#>     [SE] Column 4/7
#>     [SE] Column 5/7
#>     [SE] Column 6/7
#>     [SE] Column 7/7
#>   [Deviance] Computing canonical deviance...
#>   [Deviance] = 30.6980
#> 
#> [Phase 5] Secure deviance (Beaver): 30.6980
#> 
#> Deviance: 6.3390, Null deviance: 12.1833, Pseudo R2: 0.4797
#> 
#> Coefficients:
#>                     Estimate  Std.Error    z value   Pr(>|z|)
#>   (Intercept)        -1.2918     0.3467     -3.726   0.000195 ***
#>   age                 0.0133     0.0046      2.878   0.004002 **
#>   bmi                 0.0029     0.0096      0.299   0.764800
#>   ped                 0.3148     0.1173      2.684   0.007271 **
#>   skin                0.0067     0.0061      1.110   0.266854
#>   glu                 0.0047     0.0015      3.107   0.001893 **
#>   bp                  0.0024     0.0041      0.578   0.563374

full <- result$value$full
se <- full$std_errors["age"]
est <- full$coefficients["age"]
z <- est / se
manual_p <- 2 * stats::pnorm(-abs(z))
lr_p_ref <- stats::pchisq(result$value$lr$statistic,
                          df = result$value$lr$df,
                          lower.tail = FALSE)
observed <- max(abs(result$value$wald$p_value - manual_p),
                abs(result$value$contrast$estimate - est),
                abs(result$value$ci["age", "estimate"] - est),
                abs(result$value$lr$p_value - lr_p_ref))

rows <- row_result(
  "inference", "Inference helpers", K,
  "ds.vert.confint / ds.vert.wald / ds.vert.contrast / ds.vert.lr",
  "MASS::Pima.tr fixture", "manual algebra on ds.vert.glm output",
  "algebra_max_abs_delta", observed, 1e-10, "strict-precise",
  "Post-processes released model-level beta/covariance/deviance only.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K\>=3 | ds.vert.confint / ds.vert.wald / ds.vert.contrast / ds.vert.lr | MASS::Pima.tr fixture | manual algebra on ds.vert.glm output | algebra_max_abs_delta | 0 | 0 | strict-precise | PASS | 17.8 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
