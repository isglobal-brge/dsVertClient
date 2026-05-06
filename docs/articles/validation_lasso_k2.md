# LASSO validation (K=2)

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
[`ds.vert.lasso_proximal()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

The proximal LASSO path consumes the distributed GLM fit. The validation
uses lambda = 0, where the penalized solution must equal the prime GLM
estimate.

## Mathematical target

At lambda = 0, the proximal quadratic argmin is beta_hat, so the checked
target is max \|beta_lasso0 - beta_glm\|.

## Disclosure review

The proximal step uses only model-level beta/Hessian information and
returns coefficients; it does not request row-level data.

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

|     | server |   n | columns                                    |
|:----|:-------|----:|:-------------------------------------------|
| s1  | s1     |  60 | patient_id, age, bmi, ped                  |
| s2  | s2     |  60 | patient_id, npreg, glu, bp, skin, diabetes |

``` r


result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    fit <- dsVertClient::ds.vert.glm(
      fm, data = "DA", family = "gaussian", max_iter = 30L,
      verbose = validation_demo_verbose(), datasources = conns)
    prox0 <- dsVertClient::ds.vert.lasso_proximal(
      fit, lambda = 0, max_iter = 1000L, tol = 1e-8)
    list(fit = fit, prox0 = prox0)
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
#> [Auto-detect] Querying server columns...
#>   y_var 'diabetes' found on: s2
#>   s1: age, bmi, ped
#>   s2 (label): npreg, glu, bp, skin
#> === Encrypted-Label BCD-IRLS for gaussian GLM ===
#> Observations: 60, Variables: 7, Partitions: 2
#> Label server: s2 (holds 'diabetes')
#> Non-label servers: s1
#> 
#> [Phase 0] Transport key setup (2 servers)...
#>   [Key Setup] Transport keys exchanged (2 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 2 servers...
#>   [Standardize] 7 total features standardized (y standardized, 0.0s)
#> 
#> [Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...
#>   Gaussian one-shot: X^T X + X^T y via Beaver, direct solve
#>   [Input Sharing] Generating additive secret shares (2 servers)...
#>   [Input Sharing] Complete: p_coord=4, p_nl=3 (0.1s)
#> 
#> [Phase 3] K=2 L-BFGS iterations (p=7, n=60, lambda=1.0e-04)
#>   [L-BFGS] sum_res=-0.0002, ||grad||=0.924871, grad_range=[-0.5051, -0.0000]
#>   Iter 1: ||grad||=0.9249, step=0.30, diff=1.52e-01, theta=[0.000, 0.152] (0.2s)
#>   [L-BFGS] sum_res=-0.0005, ||grad||=0.392576, grad_range=[-0.2405, -0.0000]
#>   Iter 2: ||grad||=0.3926, step=1.00, diff=1.45e-01, theta=[0.000, 0.297] (0.2s)
#>   [L-BFGS] sum_res=-0.0002, ||grad||=0.106406, grad_range=[-0.0682, 0.0502]
#>   Iter 3: ||grad||=0.1064, step=1.00, diff=5.26e-02, theta=[0.000, 0.333] (0.3s)
#>   Iter 4: ||grad||=0.0446, step=1.00, diff=3.10e-02, theta=[0.000, 0.344] (0.2s)
#>   Iter 5: ||grad||=0.0188, step=1.00, diff=1.88e-02, theta=[0.000, 0.345] (0.2s)
#>   Iter 6: ||grad||=0.0139, step=1.00, diff=2.58e-02, theta=[0.000, 0.344] (0.2s)
#>   Iter 7: ||grad||=0.0134, step=1.00, diff=6.13e-03, theta=[-0.000, 0.341] (0.2s)
#>   Iter 8: ||grad||=0.0028, step=1.00, diff=1.95e-03, theta=[-0.000, 0.341] (0.3s)
#>   Iter 9: ||grad||=0.0006, step=1.00, diff=3.66e-04, theta=[-0.000, 0.341] (0.2s)
#>   Iter 10: ||grad||=0.0004, step=1.00, diff=8.34e-04, theta=[0.000, 0.342] (0.2s)
#>   Iter 11: ||grad||=0.0005, step=1.00, diff=3.63e-04, theta=[0.000, 0.342] (0.2s)
#>   Iter 12: ||grad||=0.0003, step=1.00, diff=1.43e-04, theta=[0.000, 0.342] (0.3s)
#>   Iter 13: ||grad||=0.0002, step=1.00, diff=6.50e-05, theta=[0.000, 0.342] (0.2s)
#>   Converged after 13 iterations (diff = 6.50e-05)
#>   [SE] Computing Hessian (central differences)...
#>     [SE] Column 1/8
#>     [SE] Column 2/8
#>     [SE] Column 3/8
#>     [SE] Column 4/8
#>     [SE] Column 5/8
#>     [SE] Column 6/8
#>     [SE] Column 7/8
#>     [SE] Column 8/8
#>   [Deviance] Computing canonical deviance...
#>   [Deviance] = 30.5305
#> 
#> [Phase 5] Secure deviance (K=2 Beaver): 30.5305
#> 
#> Deviance: 6.3045, Null deviance: 12.1833, Pseudo R2: 0.4825
#> 
#> Coefficients:
#>                     Estimate  Std.Error    z value   Pr(>|z|)
#>   (Intercept)        -1.2307     0.3667     -3.357   0.000789 ***
#>   age                 0.0118     0.0054      2.180   0.029238 *
#>   bmi                 0.0034     0.0098      0.350   0.726423
#>   ped                 0.3141     0.1181      2.659   0.007828 **
#>   npreg               0.0108     0.0199      0.543   0.586996
#>   glu                 0.0047     0.0015      3.075   0.002105 **
#>   bp                  0.0014     0.0045      0.317   0.751536
#>   skin                0.0070     0.0061      1.136   0.255947

observed <- max_named_delta(result$value$prox0$coefficients,
                            result$value$fit$coefficients)

rows <- row_result(
  "lasso", "LASSO", K, "ds.vert.lasso_proximal(lambda=0)",
  "MASS::Pima.tr fixture", "OLS limit of ds.vert.glm",
  "lambda0_coef_abs_delta", observed, 1e-8, "strict-precise",
  "Consumes only model-level GLM aggregates; no extra patient-level values are returned.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.lasso_proximal(lambda=0) | MASS::Pima.tr fixture | OLS limit of ds.vert.glm | lambda0_coef_abs_delta | 0 | 0 | strict-precise | PASS | 8.5 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
