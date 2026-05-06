# Cox PH validation (K=2)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_cox), collapse = "\n"))
#> function (n = 60L, seed = 202L, beta = c(0.4, -0.3, 0.2), n_bins = 8L) 
#> {
#>     set.seed(seed)
#>     x1 <- rnorm(n)
#>     x2 <- rnorm(n)
#>     x3 <- rnorm(n)
#>     eta <- beta[1] * x1 + beta[2] * x2 + beta[3] * x3
#>     t0 <- stats::rexp(n, rate = exp(eta)/20)
#>     cens <- stats::rexp(n, rate = 0.02)
#>     event <- as.integer(t0 <= cens)
#>     time <- pmin(t0, cens)
#>     br <- unique(stats::quantile(time, probs = seq(0, 1, length.out = n_bins + 
#>         1L), names = FALSE))
#>     time <- as.integer(cut(time, br, include.lowest = TRUE))
#>     data.frame(patient_id = sprintf("C%03d", seq_len(n)), x1, 
#>         x2, x3, time, event)
#> }
```

## Scope

Functions:
[`ds.vert.cox()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

The product Cox route is the non-disclosive Breslow profile route. The
discrete-time route is a different estimand, not the fast approximation
to this profile target, so this vignette validates the route used in the
product matrix.

## Mathematical target

The Cox PH score is sum_i delta_i \[x_i - E_beta{x \| t_i in risk
set}\]. Validation compares beta to survival::coxph(ties=‘breslow’).

## Disclosure review

Risk-set/event-time score terms remain shared. The client receives
coefficients and scalar convergence diagnostics, not event-time vectors
or patient risk-set membership.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 2L

validation_require("survival")
pooled <- build_cox(60L)
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x3", "time", "event"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "x3", "time", "event"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
```

| patient_id |         x1 |         x2 |         x3 | time | event |
|:-----------|-----------:|-----------:|-----------:|-----:|------:|
| C001       | -1.1317843 |  1.0204538 |  0.2248802 |    7 |     1 |
| C002       | -0.4404935 |  0.9188791 |  0.7585108 |    1 |     1 |
| C003       | -0.3364002 |  0.0838143 | -0.8405271 |    5 |     1 |
| C004       | -0.8470490 | -1.8174928 | -1.3290907 |    8 |     0 |
| C005       | -0.1460274 |  0.0751841 | -0.5502363 |    6 |     0 |
| C006       | -1.4305493 | -0.2279426 |  0.8314930 |    7 |     1 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns                     |
|:----|:-------|----:|:----------------------------|
| s1  | s1     |  60 | patient_id, x1, x2          |
| s2  | s2     |  60 | patient_id, x3, time, event |

``` r


ref <- coef(survival::coxph(
  survival::Surv(time, event) ~ x1 + x2 + x3,
  data = pooled, ties = "breslow"))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.cox(
      survival::Surv(time, event) ~ x1 + x2 + x3, data = "DA",
      method = "profile",
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
#> [#D' non-disclosive] K=2 OS=s2 fusion/NL=s1 session=cox_nd_20260506095519_1547559884 target=cox_profile
#> [#D' non-disclosive] Phase 1: OS share-mask
#> [#D' non-disclosive] Phase 2: NL receive shares
#> [#D' non-disclosive] Phase 3: fusion uniform expansion
#> [#D' Cox-profile iter 1] |g|_max=1.813e+01 |g|_L2=2.174e+01
#> [#D' Cox-profile iter 1] |step|_max=5.460e-01 |beta|_max=5.460e-01 cap=1.000
#> [#D' Cox-profile iter 2] |g|_max=1.670e+00 |g|_L2=2.031e+00
#> [#D' Cox-profile iter 2] |step|_max=4.221e-02 |beta|_max=5.038e-01 cap=1.000
#> [#D' Cox-profile iter 3] |g|_max=5.887e-03 |g|_L2=6.928e-03
#> [#D' Cox-profile iter 3] |step|_max=2.208e-04 |beta|_max=5.038e-01 cap=1.000
#> [#D' Cox-profile iter 4] |g|_max=2.561e-03 |g|_L2=3.035e-03
#> [#D' Cox-profile iter 4] |step|_max=7.914e-05 < 1.000e-04; stopping

fit <- result$value
observed <- max_named_delta(fit$coefficients, ref)

rows <- row_result(
  "cox", "Cox PH", K, "ds.vert.cox(method='profile')",
  "synthetic discretised survival fixture",
  "survival::coxph(ties='breslow')",
  "coef_max_abs_delta", observed, 1e-3, "strict-practical",
  "Event-time score terms remain shared; client receives beta and scalar diagnostics.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.cox(method=‘profile’) | synthetic discretised survival fixture | survival::coxph(ties=‘breslow’) | coef_max_abs_delta | 4.87e-05 | 0.001 | strict-practical | PASS | 123.6 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
