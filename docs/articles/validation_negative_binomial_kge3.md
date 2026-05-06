# Negative binomial validation (K\>=3)

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

This vignette validates the `K>=3` modality. It creates the local pooled
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

K <- 3L

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

|     | server |   n | columns           |
|:----|:-------|----:|:------------------|
| s1  | s1     |  60 | patient_id, x1    |
| s2  | s2     |  60 | patient_id, x2    |
| s3  | s3     |  60 | patient_id, x3, y |

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
#> [ds.vertNB] Stage 1: Poisson GLM point estimator
#>   [NB theta-MLE] iter 1  theta=3.86515  s=-5.621e-02  sp=-2.247e-01  step=2.502e-01
#>   [NB theta-MLE] iter 2  theta=3.89426  s=8.657e-03  sp=-2.974e-01  step=-2.911e-02
#>   [NB theta-MLE] iter 3  theta=3.89475  s=1.400e-04  sp=-2.878e-01  step=-4.866e-04
#>   [NB theta-MLE] iter 4  theta=3.89475  s=3.824e-08  sp=-2.876e-01  step=-1.330e-07
#> [ds.vertNB] joint IRLS (initial theta=3.895)
#>   joint iter 1  theta=3.895  delta=0
#> [NBFullRegND] K>=3 Ring127 setup OK n=60 DCF=(s3,s1)
#> [NBFullRegND] theta=7.7895  score=-1.3686e-01  deriv=-3.5108e-03  n=60  Sumlog(mu+theta)=141.490  Sum1/(theta+mu)=5.684  Sum(y+theta)/(theta+mu)=60.014  Sum(y+theta)/(theta+mu)^2=5.685  Sumpsi=137.434  Sumpsi_1=6.194
#> [NBFullRegND] theta=1.9474  score=+2.7435e+00  deriv=-3.6449e+00  n=60  Sumlog(mu+theta)=92.997  Sum1/(theta+mu)=12.824  Sum(y+theta)/(theta+mu)=60.034  Sum(y+theta)/(theta+mu)^2=12.835  Sumpsi=79.081  Sumpsi_1=18.372
#> [NBFullRegND] theta=2.8662  score=+8.6900e-01  deriv=-1.0650e+00  n=60  Sumlog(mu+theta)=103.743  Sum1/(theta+mu)=10.699  Sum(y+theta)/(theta+mu)=60.027  Sum(y+theta)/(theta+mu)^2=10.705  Sumpsi=93.571  Sumpsi_1=13.695
#> [NBFullRegND] theta=3.8102  score=+2.4799e-01  deriv=-3.8255e-01  n=60  Sumlog(mu+theta)=113.072  Sum1/(theta+mu)=9.147  Sum(y+theta)/(theta+mu)=60.022  Sum(y+theta)/(theta+mu)^2=9.150  Sumpsi=105.127  Sumpsi_1=11.005
#> [NBFullRegND] theta=4.5169  score=+5.2938e-02  deriv=-1.9215e-01  n=60  Sumlog(mu+theta)=119.209  Sum1/(theta+mu)=8.252  Sum(y+theta)/(theta+mu)=60.020  Sum(y+theta)/(theta+mu)^2=8.255  Sumpsi=112.396  Sumpsi_1=9.635
#> [NBFullRegND] theta=4.8010  score=+4.9478e-03  deriv=-1.4703e-01  n=60  Sumlog(mu+theta)=121.509  Sum1/(theta+mu)=7.940  Sum(y+theta)/(theta+mu)=60.019  Sum(y+theta)/(theta+mu)^2=7.943  Sumpsi=115.068  Sumpsi_1=9.182
#> [NBFullRegND] Fisher column 1/4
#> [NBFullRegND] Fisher column 2/4
#> [NBFullRegND] Fisher column 3/4
#> [NBFullRegND] Fisher column 4/4
#> [NBFullRegND] theta=4.8010  score=+4.4364e-03  deriv=-1.4700e-01  n=60  Sumlog(mu+theta)=121.528  Sum1/(theta+mu)=7.937  Sum(y+theta)/(theta+mu)=60.001  Sum(y+theta)/(theta+mu)^2=7.937  Sumpsi=115.068  Sumpsi_1=9.182
#> [NBFullRegND] theta=4.8312  score=+6.3880e-04  deriv=-1.4297e-01  n=60  Sumlog(mu+theta)=121.767  Sum1/(theta+mu)=7.905  Sum(y+theta)/(theta+mu)=60.001  Sum(y+theta)/(theta+mu)^2=7.905  Sumpsi=115.345  Sumpsi_1=9.136
#> [NBFullRegND] theta=4.8357  score=-2.3423e-04  deriv=-1.4234e-01  n=60  Sumlog(mu+theta)=121.803  Sum1/(theta+mu)=7.901  Sum(y+theta)/(theta+mu)=60.001  Sum(y+theta)/(theta+mu)^2=7.900  Sumpsi=115.386  Sumpsi_1=9.129
#> [NBFullRegND] Fisher column 1/4
#> [NBFullRegND] Fisher column 2/4
#> [NBFullRegND] Fisher column 3/4
#> [NBFullRegND] Fisher column 4/4
#> [NBFullRegND] theta=4.8357  score=-2.4485e-04  deriv=-1.4234e-01  n=60  Sumlog(mu+theta)=121.803  Sum1/(theta+mu)=7.901  Sum(y+theta)/(theta+mu)=60.000  Sum(y+theta)/(theta+mu)^2=7.900  Sumpsi=115.386  Sumpsi_1=9.129
#> [NBFullRegND] theta=4.8340  score=-1.5709e-04  deriv=-1.4255e-01  n=60  Sumlog(mu+theta)=121.790  Sum1/(theta+mu)=7.903  Sum(y+theta)/(theta+mu)=60.000  Sum(y+theta)/(theta+mu)^2=7.902  Sumpsi=115.371  Sumpsi_1=9.132
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
| K\>=3 | ds.vert.nb(method=‘accurate’) | synthetic NB fixture | MASS::glm.nb | coef_max_abs_delta | 0.0000317 | 0.02 | strict-practical | PASS | 1121.7 |
| K\>=3 | ds.vert.nb(method=‘fast’) | synthetic NB fixture | central Poisson beta + MoM theta | max(beta_delta, theta_delta) | 0.0010410 | 0.02 | fast-approximation | PASS | 81.6 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
