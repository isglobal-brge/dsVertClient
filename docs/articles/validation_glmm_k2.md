# GLMM validation (K=2)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_balanced_glmm), collapse = "\n"))
#> function (seed = 77L, ncl = 12L, m = 8L, b_sd = 0.5, beta = c(-0.2, 
#>     0.45, -0.3)) 
#> {
#>     set.seed(seed)
#>     n <- ncl * m
#>     cluster <- rep(seq_len(ncl), each = m)
#>     x1 <- rnorm(n)
#>     x2 <- rnorm(n)
#>     b <- rnorm(ncl, sd = b_sd)[cluster]
#>     eta <- beta[1] + beta[2] * x1 + beta[3] * x2 + b
#>     data.frame(patient_id = sprintf("G%03d", seq_len(n)), cluster = cluster, 
#>         x1 = x1, x2 = x2, y = stats::rbinom(n, 1L, stats::plogis(eta)))
#> }
```

## Scope

Functions:
[`ds.vert.glmm()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

This vignette shows both available product routes: PQL as the lower-cost
approximation and Laplace as the more accurate route against
lme4::glmer.

## Mathematical target

PQL alternates binomial working responses with aggregate mixed-model
equations. Laplace approximates the random-intercept marginal likelihood
and checks both fixed effects and scalar random-effect variance.

## Disclosure review

Both routes return fixed effects and scalar variance/quality diagnostics
only. Per-patient probabilities, row scores, BLUPs, cluster labels, and
cluster vectors are not returned.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 2L

validation_require(c("MASS", "nlme", "lme4"))
pooled <- build_balanced_glmm(seed = 1L, ncl = 10L, m = 5L, b_sd = 0.6)
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2", "cluster", "y"),
                   drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "cluster", "y"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
```

| patient_id | cluster |         x1 |         x2 |   y |
|:-----------|--------:|-----------:|-----------:|----:|
| G001       |       1 | -0.6264538 |  0.3981059 |   0 |
| G002       |       1 |  0.1836433 | -0.6120264 |   0 |
| G003       |       1 | -0.8356286 |  0.3411197 |   0 |
| G004       |       1 |  1.5952808 | -1.1293631 |   1 |
| G005       |       1 |  0.3295078 |  1.4330237 |   1 |
| G006       |       2 | -0.8204684 |  1.9803999 |   0 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns                    |
|:----|:-------|----:|:---------------------------|
| s1  | s1     |  50 | patient_id, x1             |
| s2  | s2     |  50 | patient_id, x2, cluster, y |

``` r


pooled_ref <- pooled
pooled_ref$cluster <- factor(pooled_ref$cluster)
ref_pql <- suppressWarnings(MASS::glmmPQL(
  y ~ x1 + x2, random = ~1 | cluster, family = binomial(),
  data = pooled_ref, verbose = validation_demo_verbose()))
#> iteration 1
#> iteration 2
#> iteration 3
ref_pql_beta <- nlme::fixef(ref_pql)
ref_laplace <- suppressWarnings(lme4::glmer(
  y ~ x1 + x2 + (1 | cluster), data = pooled_ref, family = binomial(),
  control = lme4::glmerControl(
    optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    run_pql <- elapsed(fit_pql <- dsVertClient::ds.vert.glmm(
      y ~ x1 + x2, data = "DA", cluster_col = "cluster",
      method = "pql", compute_se = FALSE,
      verbose = validation_demo_verbose(), datasources = conns))
    run_laplace <- elapsed(fit_laplace <- dsVertClient::ds.vert.glmm(
      y ~ x1 + x2, data = "DA", cluster_col = "cluster",
      method = "laplace", max_outer = 1L, mode_max_iter = 1L,
      prime_iter = 5L, compute_se = FALSE,
      verbose = validation_demo_verbose(), datasources = conns))
    list(fit_pql = fit_pql, fit_laplace = fit_laplace,
         runtime_pql = run_pql$runtime_s,
         runtime_laplace = run_laplace$runtime_s)
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
#>   s1: 50 IDs masked (stored server-side)
#> [Phase 2] s1: exporting encrypted points for s2...
#> [Phase 3] s2: processing (parallel)...
#>   s2: 50 IDs masked
#> [Phase 4-5] s2: double-masking via s1 (blind relay)...
#> [Phase 6-7] s2: matching and aligning (blind relay)...
#> [Phase 7] s1: self-aligning...
#> [Phase 8] Computing encrypted multi-server intersection...
#>   Common records: 50
#> Server 's1': 50 of 50 records matched (100.0%)
#> Server 's2': 50 of 50 records matched (100.0%)
#> PSI alignment complete.
#> [ds.vertGLMM] prime: binomial ds.vertGLM
#> [Auto-detect] Querying server columns...
#>   y_var 'y' found on: s2
#>   s1: x1
#>   s2 (label): x2
#> === Encrypted-Label BCD-IRLS for binomial GLM ===
#> Observations: 50, Variables: 2, Partitions: 2
#> Label server: s2 (holds 'y')
#> Non-label servers: s1
#> 
#> [Phase 0] Transport key setup (2 servers)...
#>   [Key Setup] Transport keys exchanged (2 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 2 servers...
#>   [Standardize] 2 total features standardized (y raw, 0.0s)
#> 
#> [Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...
#>   L-BFGS: binomial, 100-interval spline (6 rounds/iter, no Fisher)
#>   [Input Sharing] Generating additive secret shares (2 servers)...
#>   [Input Sharing] Complete: p_coord=1, p_nl=1 (0.1s)
#>   [DCF] Server s1 generating keys (n=50, 100 intervals)...
#>   [DCF] Keys distributed (76.0s)
#> 
#> [Phase 3] K=2 L-BFGS iterations (p=2, n=50, lambda=0.0e+00)
#>   [L-BFGS] sum_res=-2.0000, ||grad||=0.141719, grad_range=[-0.0594, 0.1223]
#>   Iter 1: ||grad||=0.1417, step=0.30, diff=3.67e-02, theta=[-0.037, 0.018] (4.3s)
#>   [L-BFGS] sum_res=-1.8503, ||grad||=0.131014, grad_range=[-0.0547, 0.1131]
#>   Iter 2: ||grad||=0.1310, step=1.00, diff=4.54e-01, theta=[-0.491, 0.225] (4.3s)
#>   [L-BFGS] sum_res=-0.1350, ||grad||=0.007778, grad_range=[-0.0029, 0.0067]
#>   Iter 3: ||grad||=0.0078, step=1.00, diff=2.85e-02, theta=[-0.519, 0.236] (4.4s)
#>   Iter 4: ||grad||=0.0009, step=1.00, diff=3.53e-03, theta=[-0.523, 0.238] (4.1s)
#>   Iter 5: ||grad||=0.0000, step=1.00, diff=8.63e-05, theta=[-0.523, 0.238] (4.1s)
#>   Converged after 5 iterations (diff = 8.63e-05)
#>   [Deviance] Computing canonical deviance...
#>   [Deviance] = 65.2801
#> 
#> [Phase 5] Secure deviance (K=2 Beaver): 65.2801
#> 
#> Deviance: 65.2801, Null deviance: NA, Pseudo R2: NA
#> [GLMM] share-domain cluster moments: K=10
#> [GLMM PQL] outer 1  sigma_b^2=0.1071  phi=0.9765  tau=0.1097  beta_delta=2.550e-02
#> [GLMM] share-domain cluster moments: K=10
#> [GLMM PQL] outer 2  sigma_b^2=0.1158  phi=0.9706  tau=0.1193  beta_delta=1.313e-03
#> [GLMM] share-domain cluster moments: K=10
#> [GLMM PQL] outer 3  sigma_b^2=0.1165  phi=0.971  tau=0.1199  beta_delta=7.863e-05
#> [GLMM] share-domain cluster moments: K=10
#> [GLMM PQL] outer 4  sigma_b^2=0.1165  phi=0.971  tau=0.12  beta_delta=6.770e-06
#> [ds.vertGLMMLaplace] prime: binomial ds.vertGLM
#> [Auto-detect] Querying server columns...
#>   y_var 'y' found on: s2
#>   s1: x1
#>   s2 (label): x2
#> === Encrypted-Label BCD-IRLS for binomial GLM ===
#> Observations: 50, Variables: 2, Partitions: 2
#> Label server: s2 (holds 'y')
#> Non-label servers: s1
#> 
#> [Phase 0] Transport key setup (2 servers)...
#>   [Key Setup] Transport keys exchanged (2 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 2 servers...
#>   [Standardize] 2 total features standardized (y raw, 0.0s)
#> 
#> [Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...
#>   L-BFGS: binomial, 100-interval spline (6 rounds/iter, no Fisher)
#>   [Input Sharing] Generating additive secret shares (2 servers)...
#>   [Input Sharing] Complete: p_coord=1, p_nl=1 (0.1s)
#>   [DCF] Server s1 generating keys (n=50, 100 intervals)...
#>   [DCF] Keys distributed (75.7s)
#> 
#> [Phase 3] K=2 L-BFGS iterations (p=2, n=50, lambda=0.0e+00)
#>   [L-BFGS] sum_res=-2.0000, ||grad||=0.141719, grad_range=[-0.0594, 0.1223]
#>   Iter 1: ||grad||=0.1417, step=0.30, diff=3.67e-02, theta=[-0.037, 0.018] (4.2s)
#>   [L-BFGS] sum_res=-1.8503, ||grad||=0.131014, grad_range=[-0.0547, 0.1131]
#>   Iter 2: ||grad||=0.1310, step=1.00, diff=4.54e-01, theta=[-0.491, 0.225] (4.1s)
#>   [L-BFGS] sum_res=-0.1350, ||grad||=0.007778, grad_range=[-0.0029, 0.0067]
#>   Iter 3: ||grad||=0.0078, step=1.00, diff=2.85e-02, theta=[-0.519, 0.236] (4.0s)
#>   Iter 4: ||grad||=0.0009, step=1.00, diff=3.53e-03, theta=[-0.523, 0.238] (4.1s)
#>   Iter 5: ||grad||=0.0000, step=1.00, diff=8.63e-05, theta=[-0.523, 0.238] (4.2s)
#>   Converged after 5 iterations (diff = 8.63e-05)
#> 
#> [Phase 5] Secure deviance (K=2 Beaver): NA
#> 
#> Deviance: NA, Null deviance: NA, Pseudo R2: NA
#> [ds.vertGLMMLaplace] adaptive sd starts: 0.2292
#> [ds.vertGLMMLaplace] Laplace start 1/1: sd=0.2292
#> [ds.vertGLMMLaplace] eval 1: nll=32.6219 sigma_b=0.2292 mode_step=0.121

fit_pql <- result$value$fit_pql
fit_laplace <- result$value$fit_laplace
fixed_delta <- max_named_delta(fit_pql$coefficients, ref_pql_beta)
pql_ran <- isTRUE(fit_pql$iterations >= 1L) &&
  is.data.frame(fit_pql$trace) && nrow(fit_pql$trace) >= 1L &&
  is.finite(fit_pql$sigma_b2) &&
  identical(fit_pql$quality$status, "ok")
observed_pql <- max(fixed_delta, if (pql_ran) 0 else Inf)

fixed_laplace <- max_named_delta(fit_laplace$coefficients,
                                 lme4::fixef(ref_laplace))
sigma_laplace <- abs(fit_laplace$sigma_b2 -
                       as.numeric(lme4::VarCorr(ref_laplace)$cluster)[1L])
safe_return <- identical(fit_laplace$quality$status, "ok") &&
  isFALSE(fit_laplace$disclosure$patient_level_returned) &&
  isFALSE(fit_laplace$disclosure$random_effects_returned) &&
  isFALSE(fit_laplace$disclosure$cluster_vectors_returned)
observed_laplace <- max(fixed_laplace, sigma_laplace,
                        if (safe_return) 0 else Inf)

rows <- rbind(
  row_result(
    "glmm", "GLMM PQL", K, "ds.vert.glmm(method='pql')",
    "synthetic mixed binomial random-intercept fixture",
    "MASS::glmmPQL",
    "fixed_effect_max_abs_delta_and_pql_quality", observed_pql, 0.005,
    "strict-pql",
    "PQL route returns fixed effects and scalar variance diagnostics only.",
    result$value$runtime_pql),
  row_result(
    "glmm", "GLMM Laplace", K, "ds.vert.glmm(method='laplace')",
    "synthetic mixed binomial random-intercept fixture",
    "lme4::glmer",
    "max_fixed_effect_or_sigma_b2_abs_delta", observed_laplace, 0.06,
    "laplace-practical",
    paste("Returns fixed effects, scalar variance, and scalar quality only;",
          "no BLUPs, cluster labels, cluster vectors, or row scores."),
    result$value$runtime_laplace)
)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.glmm(method=‘pql’) | synthetic mixed binomial random-intercept fixture | MASS::glmmPQL | fixed_effect_max_abs_delta_and_pql_quality | 0.0005547 | 0.005 | strict-pql | PASS | 541.7 |
| K=2 | ds.vert.glmm(method=‘laplace’) | synthetic mixed binomial random-intercept fixture | lme4::glmer | max_fixed_effect_or_sigma_b2_abs_delta | 0.0407200 | 0.060 | laplace-practical | PASS | 318.2 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
