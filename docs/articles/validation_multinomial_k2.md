# Multinomial validation (K=2)

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

``` r

cat(paste(deparse(build_multinomial), collapse = "\n"))
#> function (n = 60L, seed = 44L, shift = 0.15) 
#> {
#>     set.seed(seed)
#>     stopifnot(n%%3L == 0L)
#>     base_n <- n/3L
#>     base <- data.frame(x1 = rnorm(base_n), x2 = rnorm(base_n), 
#>         x3 = rnorm(base_n))
#>     high <- transform(base, x1 = x1 - shift, x2 = x2 + 0.5 * 
#>         shift, y_cls = "high")
#>     low <- transform(base, x1 = x1 + shift, y_cls = "low")
#>     med <- transform(base, x2 = x2 - shift, y_cls = "med")
#>     d <- rbind(high, low, med)
#>     d <- cbind(patient_id = sprintf("M%03d", seq_len(n)), d)
#>     d$y_cls <- factor(d$y_cls, levels = c("high", "low", "med"))
#>     for (k in levels(d$y_cls)) d[[paste0(k, "_ind")]] <- as.integer(d$y_cls == 
#>         k)
#>     d
#> }
```

## Scope

Functions:
[`ds.vert.multinom()`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vert.aliases.md).

This vignette validates the `K=2` modality. It creates the local pooled
fixture, partitions it vertically, starts an in-memory DSLite server,
aligns records with PSI, runs the distributed dsVert route, runs the
centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering
the vignette executes the workflow from dataset generation to assertion.

## Method

The supported route is joint softmax Newton. Historical one-vs-rest
approximations are not exposed as the product evidence path.

## Mathematical target

P(Y=c\|X)=exp(eta_c)/sum_l exp(eta_l). Validation compares class
probabilities with nnet::multinom().

## Disclosure review

Softmax probabilities, residuals, and row scores remain Ring127 shares;
no row probabilities are returned.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend,
not an Opal/Rock deployment. That does not weaken the method-level check
of what the product route returns to the analyst.

## Executed validation

``` r

K <- 2L

validation_require("nnet")
pooled <- build_multinomial()
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x3", "y_cls",
                       "high_ind", "low_ind", "med_ind"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "x3", "y_cls",
                       "high_ind", "low_ind", "med_ind"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
```

| patient_id |         x1 |         x2 |         x3 | y_cls | high_ind | low_ind | med_ind |
|:-----------|-----------:|-----------:|-----------:|:------|---------:|--------:|--------:|
| M001       |  0.5039183 |  0.1308794 | -0.1749227 | high  |        1 |       0 |       0 |
| M002       | -0.1309477 |  0.4894249 | -0.3390776 | high  |        1 |       0 |       0 |
| M003       | -1.9995040 | -0.8544092 | -0.1263074 | high  |        1 |       0 |       0 |
| M004       | -0.2827633 | -0.5823493 |  0.9804672 | high  |        1 |       0 |       0 |
| M005       | -1.3488182 |  0.0409361 |  2.4742333 | high  |        1 |       0 |       0 |
| M006       | -1.4797415 | -2.1832039 | -0.5472331 | high  |        1 |       0 |       0 |

``` r

knitr::kable(partition_summary(tables))
```

|     | server |   n | columns                                           |
|:----|:-------|----:|:--------------------------------------------------|
| s1  | s1     |  60 | patient_id, x1, x2                                |
| s2  | s2     |  60 | patient_id, x3, y_cls, high_ind, low_ind, med_ind |

``` r


ref <- nnet::multinom(y_cls ~ x1 + x2 + x3, data = pooled,
                      trace = FALSE, maxit = 200L)
ref_prob <- stats::predict(ref, pooled, type = "probs")

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.multinom(
      y_cls ~ x1 + x2 + x3, data = "DA",
      classes = c("high", "low", "med"),
      indicator_template = "%s_ind",
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
#> [ds.vertMultinom] dispatching to ds.vertMultinomJointNewton
#> [ds.vertMultinom] fitting class 'low' vs rest (indicator 'low_ind')
#> [Auto-detect] Querying server columns...
#>   y_var 'low_ind' found on: s2
#>   s1: x1, x2
#>   s2 (label): x3
#> === Encrypted-Label BCD-IRLS for binomial GLM ===
#> Observations: 60, Variables: 3, Partitions: 2
#> Label server: s2 (holds 'low_ind')
#> Non-label servers: s1
#> 
#> [Phase 0] Transport key setup (2 servers)...
#>   [Key Setup] Transport keys exchanged (2 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 2 servers...
#>   [Standardize] 3 total features standardized (y raw, 0.0s)
#> 
#> [Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...
#>   L-BFGS: binomial, 100-interval spline (6 rounds/iter, no Fisher)
#>   [Input Sharing] Generating additive secret shares (2 servers)...
#>   [Input Sharing] Complete: p_coord=1, p_nl=2 (0.1s)
#>   [DCF] Server s1 generating keys (n=60, 100 intervals)...
#>   [DCF] Keys distributed (28.8s)
#> 
#> [Phase 3] K=2 L-BFGS iterations (p=3, n=60, lambda=1.0e-04)
#>   [L-BFGS] sum_res=10.0000, ||grad||=0.172104, grad_range=[-0.0420, 0.1667]
#>   Iter 1: ||grad||=0.1721, step=0.30, diff=5.00e-02, theta=[-0.050, 0.013] (2.4s)
#>   [L-BFGS] sum_res=9.2516, ||grad||=0.159226, grad_range=[-0.0389, 0.1542]
#>   Iter 2: ||grad||=0.1592, step=1.00, diff=6.18e-01, theta=[-0.668, 0.171] (2.2s)
#>   [L-BFGS] sum_res=0.4043, ||grad||=0.007906, grad_range=[-0.0041, 0.0067]
#>   Iter 3: ||grad||=0.0079, step=1.00, diff=2.82e-02, theta=[-0.696, 0.189] (2.3s)
#>   Iter 4: ||grad||=0.0009, step=1.00, diff=3.16e-03, theta=[-0.699, 0.191] (2.2s)
#>   Iter 5: ||grad||=0.0000, step=1.00, diff=1.10e-04, theta=[-0.699, 0.191] (2.2s)
#>   Iter 6: ||grad||=0.0000, step=1.00, diff=2.78e-05, theta=[-0.699, 0.191] (2.3s)
#>   Converged after 6 iterations (diff = 2.78e-05)
#>   [SE] Computing Hessian (central differences)...
#>     [SE] Column 1/4
#>     [SE] Column 2/4
#>     [SE] Column 3/4
#>     [SE] Column 4/4
#>   [Deviance] Computing canonical deviance...
#>   [Deviance] = 75.9647
#> 
#> [Phase 5] Secure deviance (K=2 Beaver): 75.9647
#> 
#> Deviance: 75.9647, Null deviance: NA, Pseudo R2: NA
#> 
#> Coefficients:
#>                     Estimate  Std.Error    z value   Pr(>|z|)
#>   (Intercept)        -0.6763     0.2815     -2.403   0.016276 *
#>   x1                  0.1607     0.2318      0.693   0.488095
#>   x2                  0.0380     0.3001      0.126   0.899368
#>   x3                  0.0127     0.2051      0.062   0.950569
#> [ds.vertMultinom] fitting class 'med' vs rest (indicator 'med_ind')
#> [Auto-detect] Querying server columns...
#>   y_var 'med_ind' found on: s2
#>   s1: x1, x2
#>   s2 (label): x3
#> === Encrypted-Label BCD-IRLS for binomial GLM ===
#> Observations: 60, Variables: 3, Partitions: 2
#> Label server: s2 (holds 'med_ind')
#> Non-label servers: s1
#> 
#> [Phase 0] Transport key setup (2 servers)...
#>   [Key Setup] Transport keys exchanged (2 servers, 0.1s)
#> 
#> [Phase 1] Standardizing features across 2 servers...
#>   [Standardize] 3 total features standardized (y raw, 0.0s)
#> 
#> [Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...
#>   L-BFGS: binomial, 100-interval spline (6 rounds/iter, no Fisher)
#>   [Input Sharing] Generating additive secret shares (2 servers)...
#>   [Input Sharing] Complete: p_coord=1, p_nl=2 (0.1s)
#>   [DCF] Server s1 generating keys (n=60, 100 intervals)...
#>   [DCF] Keys distributed (28.8s)
#> 
#> [Phase 3] K=2 L-BFGS iterations (p=3, n=60, lambda=1.0e-04)
#>   [L-BFGS] sum_res=10.0000, ||grad||=0.172242, grad_range=[-0.0000, 0.1667]
#>   Iter 1: ||grad||=0.1722, step=0.30, diff=5.00e-02, theta=[-0.050, 0.000] (2.5s)
#>   [L-BFGS] sum_res=9.2516, ||grad||=0.159364, grad_range=[-0.0005, 0.1542]
#>   Iter 2: ||grad||=0.1594, step=1.00, diff=6.19e-01, theta=[-0.669, 0.025] (2.2s)
#>   [L-BFGS] sum_res=0.4000, ||grad||=0.008387, grad_range=[-0.0013, 0.0066]
#>   Iter 3: ||grad||=0.0084, step=1.00, diff=2.81e-02, theta=[-0.697, 0.028] (2.3s)
#>   Iter 4: ||grad||=0.0012, step=1.00, diff=3.83e-03, theta=[-0.700, 0.032] (2.2s)
#>   Iter 5: ||grad||=0.0001, step=1.00, diff=5.38e-04, theta=[-0.700, 0.032] (2.3s)
#>   Iter 6: ||grad||=0.0001, step=1.00, diff=1.30e-04, theta=[-0.700, 0.032] (2.3s)
#>   Iter 7: ||grad||=0.0000, step=1.00, diff=3.29e-05, theta=[-0.700, 0.032] (2.2s)
#>   Converged after 7 iterations (diff = 3.29e-05)
#>   [SE] Computing Hessian (central differences)...
#>     [SE] Column 1/4
#>     [SE] Column 2/4
#>     [SE] Column 3/4
#>     [SE] Column 4/4
#>   [Deviance] Computing canonical deviance...
#>   [Deviance] = 75.9351
#> 
#> [Phase 5] Secure deviance (K=2 Beaver): 75.9351
#> 
#> Deviance: 75.9351, Null deviance: NA, Pseudo R2: NA
#> 
#> Coefficients:
#>                     Estimate  Std.Error    z value   Pr(>|z|)
#>   (Intercept)        -0.6859     0.2818     -2.434   0.014936 *
#>   x1                  0.0100     0.2332      0.043   0.965843
#>   x2                 -0.2128     0.2937     -0.725   0.468717
#>   x3                  0.0230     0.2054      0.112   0.910833
#> [MultinomJointNewton] session multinomJoint_1778057179_1504136423  n=60  classes=3  DCF=(s2,s1)  p=4
#> [MnlJoint] iter 1 |g|_L2=3.280e-02 |g|_max=2.566e-02 |step|_pre=8.053e-02 |step|_post=8.053e-02 beta_max=2.393e-01 dt=24.5s
#> [MnlJoint] iter 2 |g|_L2=1.164e-02 |g|_max=9.235e-03 |step|_pre=2.847e-02 |step|_post=2.847e-02 beta_max=2.496e-01 dt=24.8s
#> [MnlJoint] iter 3 |g|_L2=4.241e-03 |g|_max=3.317e-03 |step|_pre=1.019e-02 |step|_post=1.019e-02 beta_max=2.535e-01 dt=25.4s
#> [MnlJoint] iter 4 |g|_L2=1.653e-03 |g|_max=1.260e-03 |step|_pre=3.945e-03 |step|_post=3.945e-03 beta_max=2.549e-01 dt=24.9s
#> [MnlJoint] iter 5 |g|_L2=6.319e-04 |g|_max=4.602e-04 |step|_pre=1.537e-03 |step|_post=1.537e-03 beta_max=2.555e-01 dt=25.9s
#> [MnlJoint] iter 6 |g|_L2=2.582e-04 |g|_max=1.989e-04 |step|_pre=6.571e-04 |step|_post=6.571e-04 beta_max=2.557e-01 dt=24.4s
#> [MnlJoint] iter 7 |g|_L2=1.211e-04 |g|_max=9.479e-05 |step|_pre=3.073e-04 |step|_post=3.073e-04 beta_max=2.558e-01 dt=25.2s
#> [MnlJoint] iter 8 |g|_L2=3.870e-05 |g|_max=2.609e-05 |step|_pre=1.556e-04 |step|_post=1.556e-04 beta_max=2.559e-01 dt=24.9s
#> [MnlJoint] iter 9 |g|_L2=4.307e-05 |g|_max=3.629e-05 |step|_pre=8.226e-05 |step|_post=8.226e-05 beta_max=2.559e-01 dt=25.0s
#> [MnlJoint] iter 10 |g|_L2=4.697e-05 |g|_max=3.778e-05 |step|_pre=4.981e-05 |step|_post=4.981e-05 beta_max=2.558e-01 dt=26.0s

fit <- result$value
ds_prob <- softmax_prob(fit$coefficients, pooled)
ds_prob <- ds_prob[, colnames(ref_prob), drop = FALSE]
observed <- max(abs(ds_prob - ref_prob))

rows <- row_result(
  "multinomial", "Multinomial", K, "ds.vert.multinom",
  "balanced soft-signal synthetic 3-class fixture",
  "nnet::multinom probabilities",
  "class_probability_max_abs_delta", observed, 0.001,
  "strict-practical",
  "Softmax probabilities and residuals remain Ring127 shares; no row probabilities are returned.",
  result$runtime_s)

display_validation(rows)
```

| k_mode | function_route | dataset | reference_target | primary_metric | observed | tolerance | tier | status | runtime_s |
|:---|:---|:---|:---|:---|---:|---:|:---|:---|---:|
| K=2 | ds.vert.multinom | balanced soft-signal synthetic 3-class fixture | nnet::multinom probabilities | class_probability_max_abs_delta | 3.23e-05 | 0.001 | strict-practical | PASS | 427.4 |

``` r

assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
