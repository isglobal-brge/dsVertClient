
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, warning = FALSE,
  message = TRUE, collapse = TRUE, comment = "#>")
helper_candidates <- c("validation_helpers.R",
  file.path("..", "validation_helpers.R"),
  file.path("vignettes", "validation_helpers.R"))
source(helper_candidates[file.exists(helper_candidates)][1])
validation_load_packages()
options(dsvert.beaver_preprocessing = "auto")
partition_summary <- function(tables) {
  data.frame(
    server = names(tables),
    n = vapply(tables, nrow, integer(1)),
    columns = vapply(tables, function(x) paste(names(x), collapse = ", "), character(1)),
    stringsAsFactors = FALSE)
}
```

## Dataset generator

The fixture is generated inside this vignette run. The helper body used
to create the pooled local data is printed here before the analysis.

```{r dataset-generator}
cat(paste(deparse(build_lmm), collapse = "\n"))
```

## Scope

Functions: `ds.vert.lmm()`.

This vignette validates the `K=2` modality. It creates the local
pooled fixture, partitions it vertically, starts an in-memory DSLite
server, aligns records with PSI, runs the distributed dsVert route, runs
the centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering the
vignette executes the workflow from dataset generation to assertion.

## Method

Random-intercept LMM uses guarded cluster metadata and protected cross-products to fit fixed effects and variance components.

## Mathematical target

The model is y = X beta + b_cluster + epsilon with b ~ N(0,sigma_b^2). Validation compares fixed effects to lme4::lmer().

## Disclosure review

Cluster membership is used internally at the accepted method tier. Per-cluster residuals and BLUP vectors are not returned.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend, not
an Opal/Rock deployment. That does not weaken the method-level check of
what the product route returns to the analyst.

## Executed validation

```{r execution}
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
knitr::kable(partition_summary(tables))

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

fit <- result$value
observed <- max_named_delta(fit$coefficients, ref)

rows <- row_result(
  "lmm", "LMM", K, "ds.vert.lmm",
  "synthetic random-intercept fixture", "lme4::lmer",
  "fixed_effect_max_abs_delta", observed, 0.001, "strict-practical",
  "Cluster membership is used internally; per-cluster residuals/BLUPs are not returned.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
