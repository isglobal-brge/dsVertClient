# Negative binomial validation (K=2)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_nb), collapse = "\n"))
#> function (n = 60L, seed = 303L, theta = 3, beta = c(1, 0.2, -0.15, 
#>     0.1)) 
#> {
#>     set.seed(seed)
#>     x1 <- rnorm(n)
#>     x2 <- rnorm(n)
#>     x3 <- rnorm(n)
#>     mu <- exp(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * 
#>         x3)
#>     y <- stats::rnbinom(n, size = theta, mu = mu)
#>     data.frame(patient_id = sprintf("N%03d", seq_len(n)), x1, 
#>         x2, x3, y)
#> }
```

## Scope

Functions:
[`ds.vert.nb()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

This vignette shows the accurate non-disclosive full-regression route
and the fast MoM route. The accurate route targets MASS::glm.nb(); the
fast route targets a central Poisson beta plus closed-form MoM theta
approximation.

## Mathematical target

NB2 uses Var(Y\|X)=mu+mu^2/theta and log(mu)=X beta. The accurate route
evaluates beta/theta score pieces in shares; the fast route uses
theta_MoM = ybar^(2/(s_y)2-ybar).

## Disclosure review

The accurate route keeps eta, mu, reciprocal, score, and Fisher pieces
in shares and returns beta plus scalar theta. The fast route opens only
outcome-level scalar moments and Poisson-model aggregates.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 2L

validation_require("MASS")
pooled <- build_nb(60L)
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x3", "y"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "x3", "y"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
```

| patient_id |         x1 |         x2 |         x3 |   y |
|:-----------|-----------:|-----------:|-----------:|----:|
| N001       | -0.4969391 | -0.3490398 |  1.4723099 |   3 |
| N002       |  0.2788383 | -1.6389588 |  0.3217739 |   1 |
| N003       |  0.3844678 | -0.0714228 | -0.4890526 |   6 |
| N004       |  0.0685822 |  0.5115725 |  0.6384081 |   3 |
| N005       | -0.5782101 | -1.0443054 | -0.6071904 |   4 |
| N006       |  1.3138932 | -0.4724609 |  1.0910059 |   3 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns            |
|:----|:-------|----:|:-------------------|
| s1  | s1     |  60 | patient_id, x1, x2 |
| s2  | s2     |  60 | patient_id, x3, y  |

``` r


ref_accurate <- coef(MASS::glm.nb(y ~ x1 + x2 + x3, data = pooled))
ref_fast_beta <- coef(stats::glm(y ~ x1 + x2 + x3, data = pooled,
                                 family = poisson()))
ybar <- mean(pooled$y)
yvar <- stats::var(pooled$y)
ref_fast_theta <- if (is.finite(yvar) && yvar > ybar) {
  ybar * ybar / (yvar - ybar)
} else {
  Inf
}

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    run_accurate <- elapsed(accurate <- dsVertClient::ds.vert.nb(
      y ~ x1 + x2 + x3, data = "DA", method = "accurate",
      max_iter = 25L,
      verbose = validation_demo_verbose(), datasources = conns))
    run_fast <- elapsed(fast <- dsVertClient::ds.vert.nb(
      y ~ x1 + x2 + x3, data = "DA", method = "fast",
      max_iter = 25L,
      verbose = validation_demo_verbose(), datasources = conns))
    list(accurate = accurate, fast = fast,
         runtime_accurate = run_accurate$runtime_s,
         runtime_fast = run_fast$runtime_s)
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
#> [ds.vertNB] Stage 1: Poisson GLM point estimator
#>   [NB theta-MLE] iter 1  theta=3.86515  s=-5.621e-02  sp=-2.247e-01  step=2.502e-01
#>   [NB theta-MLE] iter 2  theta=3.89426  s=8.657e-03  sp=-2.974e-01  step=-2.911e-02
#>   [NB theta-MLE] iter 3  theta=3.89475  s=1.400e-04  sp=-2.878e-01  step=-4.866e-04
#>   [NB theta-MLE] iter 4  theta=3.89475  s=3.824e-08  sp=-2.876e-01  step=-1.330e-07
#> [ds.vertNB] joint IRLS (initial theta=3.895)
#>   joint iter 1  theta=3.895  delta=0
#> [NBFullRegND] K>=3 Ring127 setup OK n=60 DCF=(s2,s1)
#> [NBFullRegND] theta=7.7895  score=-1.3686e-01  deriv=-3.5109e-03  n=60  Sumlog(mu+theta)=141.492  Sum1/(theta+mu)=5.683  Sum(y+theta)/(theta+mu)=60.013  Sum(y+theta)/(theta+mu)^2=5.685  Sumpsi=137.434  Sumpsi_1=6.194
#> [NBFullRegND] theta=1.9474  score=+2.7435e+00  deriv=-3.6449e+00  n=60  Sumlog(mu+theta)=93.000  Sum1/(theta+mu)=12.823  Sum(y+theta)/(theta+mu)=60.031  Sum(y+theta)/(theta+mu)^2=12.833  Sumpsi=79.081  Sumpsi_1=18.372
#> [NBFullRegND] theta=2.8662  score=+8.6900e-01  deriv=-1.0650e+00  n=60  Sumlog(mu+theta)=103.746  Sum1/(theta+mu)=10.698  Sum(y+theta)/(theta+mu)=60.024  Sum(y+theta)/(theta+mu)^2=10.704  Sumpsi=93.571  Sumpsi_1=13.695
#> [NBFullRegND] theta=3.8102  score=+2.4803e-01  deriv=-3.8255e-01  n=60  Sumlog(mu+theta)=113.075  Sum1/(theta+mu)=9.146  Sum(y+theta)/(theta+mu)=60.020  Sum(y+theta)/(theta+mu)^2=9.150  Sumpsi=105.127  Sumpsi_1=11.005
#> [NBFullRegND] theta=4.5170  score=+5.3271e-02  deriv=-1.9217e-01  n=60  Sumlog(mu+theta)=119.212  Sum1/(theta+mu)=8.252  Sum(y+theta)/(theta+mu)=60.018  Sum(y+theta)/(theta+mu)^2=8.254  Sumpsi=112.398  Sumpsi_1=9.635
#> [NBFullRegND] theta=4.8029  score=+5.0290e-03  deriv=-1.4681e-01  n=60  Sumlog(mu+theta)=121.526  Sum1/(theta+mu)=7.938  Sum(y+theta)/(theta+mu)=60.017  Sum(y+theta)/(theta+mu)^2=7.940  Sumpsi=115.086  Sumpsi_1=9.179
#> [NBFullRegND] Fisher column 1/4
#> [NBFullRegND] Fisher column 2/4
#> [NBFullRegND] Fisher column 3/4
#> [NBFullRegND] Fisher column 4/4
#> [NBFullRegND] theta=4.8029  score=+4.5170e-03  deriv=-1.4678e-01  n=60  Sumlog(mu+theta)=121.543  Sum1/(theta+mu)=7.935  Sum(y+theta)/(theta+mu)=60.001  Sum(y+theta)/(theta+mu)^2=7.935  Sumpsi=115.086  Sumpsi_1=9.179
#> [NBFullRegND] theta=4.8338  score=-1.6957e-04  deriv=-1.4257e-01  n=60  Sumlog(mu+theta)=121.788  Sum1/(theta+mu)=7.903  Sum(y+theta)/(theta+mu)=60.001  Sum(y+theta)/(theta+mu)^2=7.902  Sumpsi=115.369  Sumpsi_1=9.132
#> [NBFullRegND] theta=4.8326  score=+8.5011e-05  deriv=-1.4274e-01  n=60  Sumlog(mu+theta)=121.778  Sum1/(theta+mu)=7.904  Sum(y+theta)/(theta+mu)=60.001  Sum(y+theta)/(theta+mu)^2=7.904  Sumpsi=115.358  Sumpsi_1=9.134
#> [NBFullRegND] Fisher column 1/4
#> [NBFullRegND] Fisher column 2/4
#> [NBFullRegND] Fisher column 3/4
#> [NBFullRegND] Fisher column 4/4
#> [NBFullRegND] theta=4.8326  score=+7.4391e-05  deriv=-1.4274e-01  n=60  Sumlog(mu+theta)=121.779  Sum1/(theta+mu)=7.904  Sum(y+theta)/(theta+mu)=60.000  Sum(y+theta)/(theta+mu)^2=7.904  Sumpsi=115.358  Sumpsi_1=9.134
#> [NBFullRegND] theta=4.8331  score=+1.6755e-04  deriv=-1.4269e-01  n=60  Sumlog(mu+theta)=121.783  Sum1/(theta+mu)=7.903  Sum(y+theta)/(theta+mu)=60.000  Sum(y+theta)/(theta+mu)^2=7.903  Sumpsi=115.363  Sumpsi_1=9.133
#> [NBFullRegND] Fisher column 1/4
#> [NBFullRegND] Fisher column 2/4
#> [NBFullRegND] Fisher column 3/4
#> [NBFullRegND] Fisher column 4/4
#> [ds.vertNBMoMTheta] Stage 1: Poisson GLM point estimator
#> [ds.vertNBMoMTheta] theta_MoM = ybar^2/(s^2-ybar) = 2.8^2/(4.705-2.8) = 4.115 (n=60)

accurate_delta <- max_named_delta(result$value$accurate$coefficients,
                                  ref_accurate)
fast_beta_delta <- max_named_delta(result$value$fast$coefficients,
                                   ref_fast_beta)
fast_theta_delta <- abs(result$value$fast$theta - ref_fast_theta)
fast_delta <- max(fast_beta_delta, fast_theta_delta)

rows <- rbind(
  row_result(
    "negative_binomial", "Negative binomial accurate", K,
    "ds.vert.nb(method='accurate')",
    "synthetic NB fixture", "MASS::glm.nb",
    "coef_max_abs_delta", accurate_delta, 0.02,
    "strict-practical",
    "Full-regression beta/theta score components remain aggregate.",
    result$value$runtime_accurate),
  row_result(
    "negative_binomial", "Negative binomial fast", K,
    "ds.vert.nb(method='fast')",
    "synthetic NB fixture", "central Poisson beta + MoM theta",
    "max(beta_delta, theta_delta)", fast_delta, 0.02,
    "fast-approximation",
    "Fast path uses outcome-only scalar moments for theta and Poisson beta.",
    result$value$runtime_fast)
)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.nb(method=‘accurate’) | synthetic NB fixture | MASS::glm.nb | coef_max_abs_delta | 0.0000277 | 0.02 | strict-practical | PASS | 1072.8 |
| K=2 | ds.vert.nb(method=‘fast’) | synthetic NB fixture | central Poisson beta + MoM theta | max(beta_delta, theta_delta) | 0.0010410 | 0.02 | fast-approximation | PASS | 79.3 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
