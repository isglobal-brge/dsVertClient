# IPW validation (K=2)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_ipw), collapse = "\n"))
#> function (n = 60L, seed = 88L, beta_p = c(-0.2, 0.4, -0.3), beta_y = c(1, 
#>     2, 0.5, -0.2), noise_sd = 0.2) 
#> {
#>     set.seed(seed)
#>     w1 <- rnorm(n)
#>     w2 <- rnorm(n)
#>     p <- stats::plogis(beta_p[1] + beta_p[2] * w1 + beta_p[3] * 
#>         w2)
#>     tr <- stats::rbinom(n, 1, p)
#>     ipw <- ifelse(tr == 1, 1/p, 1/(1 - p))
#>     y <- beta_y[1] + beta_y[2] * tr + beta_y[3] * w1 + beta_y[4] * 
#>         w2 + rnorm(n, 0, noise_sd)
#>     data.frame(patient_id = sprintf("I%03d", seq_len(n)), w1, 
#>         w2, tr, ipw, y)
#> }
```

## Scope

Functions:
[`ds.vert.ipw()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

IPW runs a protected propensity GLM and then a protected weighted
outcome GLM using a server-side weight column.

## Mathematical target

The checked outcome target solves X’W(y-X beta)=0 with the same IPW
weights used by the centralized lm().

## Disclosure review

Weights remain server-side; the analyst receives only propensity and
outcome model-level fits.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 2L

pooled <- build_ipw()
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "w1", "w2"), drop = FALSE],
       s2 = pooled[, c("patient_id", "tr", "ipw", "y"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "w1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "w2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "tr", "ipw", "y"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
```

| patient_id |         w1 |         w2 |  tr |      ipw |          y |
|:-----------|-----------:|-----------:|----:|---------:|-----------:|
| I001       | -0.2262557 |  0.4949529 |   1 | 2.551137 |  3.0325474 |
| I002       |  0.6465698 | -0.5307771 |   1 | 1.804237 |  3.2010731 |
| I003       |  2.3460583 |  0.0803406 |   0 | 3.042799 |  2.1916986 |
| I004       | -1.8456187 | -0.2881362 |   0 | 1.426643 | -0.0567660 |
| I005       |  0.4596911 | -0.3813339 |   1 | 1.906398 |  3.3391627 |
| I006       |  0.1237811 |  0.5091922 |   0 | 1.738418 |  0.9774241 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns                |
|:----|:-------|----:|:-----------------------|
| s1  | s1     |  60 | patient_id, w1, w2     |
| s2  | s2     |  60 | patient_id, tr, ipw, y |

``` r


ref <- coef(stats::lm(y ~ tr + w1 + w2, data = pooled, weights = ipw))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.ipw(
      y ~ tr + w1 + w2, tr ~ w1 + w2, data = "DA",
      precision = "high",
      outcome_family = "gaussian",
      compute_se = FALSE, compute_deviance = FALSE,
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
#> [ds.vertIPW] Stage 1: propensity model
#> [Auto-detect] Querying server columns...
#>   y_var 'tr' found on: s2
#>   s1: w1, w2
#>   s2 (label): (response only)
#> === Encrypted-Label BCD-IRLS for binomial GLM ===
#> Observations: 60, Variables: 2, Partitions: 2
#> Label server: s2 (holds 'tr')
#> Non-label servers: s1
#> 
#> [Phase 0] Transport key setup (2 servers)...
#>   [Key Setup] Transport keys exchanged (2 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 2 servers...
#>   [Standardize] 2 total features standardized (y raw, 0.0s)
#> 
#> [Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...
#>   L-BFGS: binomial, 150-interval spline (6 rounds/iter, no Fisher)
#>   [Input Sharing] Generating additive secret shares (2 servers)...
#>   [Input Sharing] Complete: p_coord=0, p_nl=2 (0.1s)
#>   [DCF] Server s1 generating keys (n=60, 150 intervals)...
#>   [DCF] Keys distributed (47.0s)
#> 
#> [Phase 3] K=2 L-BFGS iterations (p=2, n=60, lambda=1.0e-04)
#>   [L-BFGS] sum_res=6.0000, ||grad||=0.168099, grad_range=[-0.1159, 0.1000]
#>   Iter 1: ||grad||=0.1681, step=0.30, diff=3.48e-02, theta=[-0.030, 0.035] (3.2s)
#>   [L-BFGS] sum_res=5.5505, ||grad||=0.154513, grad_range=[-0.1066, 0.0925]
#>   Iter 2: ||grad||=0.1545, step=1.00, diff=3.98e-01, theta=[-0.400, 0.433] (3.2s)
#>   [L-BFGS] sum_res=0.4062, ||grad||=0.012532, grad_range=[-0.0097, 0.0067]
#>   Iter 3: ||grad||=0.0125, step=1.00, diff=4.01e-02, theta=[-0.430, 0.473] (3.5s)
#>   Iter 4: ||grad||=0.0019, step=1.00, diff=7.63e-03, theta=[-0.435, 0.481] (3.3s)
#>   Iter 5: ||grad||=0.0000, step=1.00, diff=1.94e-04, theta=[-0.435, 0.481] (3.1s)
#>   Iter 6: ||grad||=0.0000, step=1.00, diff=4.44e-05, theta=[-0.435, 0.481] (3.3s)
#>   Converged after 6 iterations (diff = 4.44e-05)
#> 
#> [Phase 5] Secure deviance (K=2 Beaver): NA
#> 
#> Deviance: NA, Null deviance: NA, Pseudo R2: NA
#> [ds.vertIPW] Stage 2: outcome model (weighted)
#> [Auto-detect] Querying server columns...
#>   y_var 'y' found on: s2
#>   s1: w1, w2
#>   s2 (label): tr
#> === Encrypted-Label BCD-IRLS for gaussian GLM ===
#> Observations: 60, Variables: 3, Partitions: 2
#> Label server: s2 (holds 'y')
#> Non-label servers: s1
#> 
#> [Phase 0] Transport key setup (2 servers)...
#>   [Key Setup] Transport keys exchanged (2 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 2 servers...
#>   [Standardize] 3 total features standardized (y standardized, 0.0s)
#> Registering weights 'ipw' on server s2
#> 
#> [Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...
#>   Gaussian one-shot: X^T X + X^T y via Beaver, direct solve
#>   [Input Sharing] Generating additive secret shares (2 servers)...
#>   [Input Sharing] Complete: p_coord=1, p_nl=2 (0.2s)
#> 
#> [Phase 3] K=2 L-BFGS iterations (p=3, n=60, lambda=1.0e-04)
#>   [L-BFGS] sum_res=-12.8466, ||grad||=1.884134, grad_range=[-1.5982, 0.3591]
#>   Iter 1: ||grad||=1.8841, step=0.30, diff=4.79e-01, theta=[-0.108, 0.479] (0.4s)
#>   [L-BFGS] sum_res=3.1533, ||grad||=0.667044, grad_range=[-0.5939, 0.0632]
#>   Iter 2: ||grad||=0.6670, step=1.00, diff=2.92e-01, theta=[-0.125, 0.772] (0.4s)
#>   [L-BFGS] sum_res=1.4765, ||grad||=0.027053, grad_range=[-0.0088, 0.0246]
#>   Iter 3: ||grad||=0.0271, step=1.00, diff=1.20e-02, theta=[-0.121, 0.774] (0.4s)
#>   Iter 4: ||grad||=0.0033, step=1.00, diff=1.24e-03, theta=[-0.120, 0.773] (0.4s)
#>   Iter 5: ||grad||=0.0004, step=1.00, diff=1.93e-04, theta=[-0.120, 0.773] (0.4s)
#>   Iter 6: ||grad||=0.0001, step=1.00, diff=7.91e-05, theta=[-0.120, 0.773] (0.4s)
#>   Converged after 6 iterations (diff = 7.91e-05)
#> 
#> [Phase 5] Secure deviance (K=2 Beaver): NA
#> 
#> Deviance: NA, Null deviance: 98.9560, Pseudo R2: NA

fit <- result$value
observed <- max_named_delta(fit$outcome$coefficients, ref)

rows <- row_result(
  "ipw", "IPW", K, "ds.vert.ipw(precision='high')",
  "synthetic confounded IPW fixture",
  "central weighted lm using same weights",
  "weighted_outcome_coef_abs_delta", observed, 1e-3,
  "strict-practical",
  "Uses protected propensity/outcome GLM aggregates; only model-level fits are returned.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.ipw(precision=‘high’) | synthetic confounded IPW fixture | central weighted lm using same weights | weighted_outcome_coef_abs_delta | 0.0001785 | 0.001 | strict-practical | PASS | 71.6 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
