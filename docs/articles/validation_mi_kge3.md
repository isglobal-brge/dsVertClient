
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
cat(paste(deparse(build_mi), collapse = "\n"))
```

## Scope

Functions: `ds.vert.mi()`.

This vignette validates the `K>=3` modality. It creates the local
pooled fixture, partitions it vertically, starts an in-memory DSLite
server, aligns records with PSI, runs the distributed dsVert route, runs
the centralized R reference on the same pooled fixture, and computes the
numerical delta.

No external RDS, cache, or pre-loaded server table is used. Rendering the
vignette executes the workflow from dataset generation to assertion.

## Method

Missing covariates are imputed server-side for each round with the
same server-local Bayesian-ridge helper used by the product route.
The completed-data GLM fits are then pooled client-side using Rubin's
rules.

## Mathematical target

The reference mirrors the route centrally on the same vertical split:
it applies `dsvertImputeColumnDS()` to the split table that owns the
missing column, fits one centralized GLM per imputation round and applies
Rubin pooling. This checks the distributed orchestration, fitting and
pooling against the same imputation model rather than comparing it with
a different imputation method.

## Disclosure review

Imputed columns stay server-side. The client receives pooled beta/covariance summaries, not imputed patient-level values.

The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer
pinning is disabled only because DSLite is an in-memory test backend, not
an Opal/Rock deployment. That does not weaken the method-level check of
what the product route returns to the analyst.

## Executed validation

```{r execution}
K <- 3L

pooled <- build_mi()
tables <- split_mi(pooled, K)

knitr::kable(utils::head(pooled))
knitr::kable(partition_summary(tables))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns, na.action = "none")
    dsVertClient::ds.vert.mi(
      y ~ x1 + x2 + x3, data = "DA", impute_columns = "x2",
      family = "gaussian", lambda = 0,
      verbose = validation_demo_verbose(), datasources = conns, seed = 12L)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

fit <- result$value
ref <- central_mi_same_imputer(
  pooled, K, y ~ x1 + x2 + x3, impute_columns = "x2",
  m = fit$m, seed = 12L)
observed <- max_named_delta(fit$coefficients, ref$coefficients)
se_observed <- max_named_delta(fit$std_errors, ref$std_errors)

knitr::kable(ref$imputation_log)
cat(sprintf("Maximum pooled coefficient delta: %.6g\n", observed))
cat(sprintf("Maximum Rubin standard-error delta: %.6g\n", se_observed))

rows <- row_result(
  "mi", "Multiple imputation", K, "ds.vert.mi",
  "synthetic missing-covariate fixture",
  "same-imputer centralized Rubin pooling",
  "pooled_coef_abs_delta", observed, 0.001, "strict-practical",
  "Imputed columns stay server-side; client pools beta/covariance only.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
```

## Verdict

Rendering fails if the distributed result leaves the accepted numerical
envelope or if the method row is marked disclosive.
