#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)[1]
script_path <- if (length(file_arg) && !is.na(file_arg)) {
  sub("^--file=", "", file_arg)
} else {
  "tools/build_validation_vignettes.R"
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(root)

dir.create("vignettes", showWarnings = FALSE, recursive = TRUE)
unlink(Sys.glob("vignettes/validation_*.Rmd"), force = TRUE)
unlink(Sys.glob("vignettes/validation_*.html"), force = TRUE)
unlink("vignettes/vert_validation_evidence.Rmd", force = TRUE)
unlink("vignettes/vert_validation_evidence.html", force = TRUE)

literal <- function(x) {
  x <- sub("^\\n", "", x)
  x <- sub("\\n$", "", x)
  strsplit(x, "\n", fixed = TRUE)[[1]]
}

chunk <- function(label, code, include = TRUE) {
  c(
    sprintf("```{r %s%s}", label, if (include) "" else ", include=FALSE"),
    literal(code),
    "```"
  )
}

setup_chunk <- function(builder = NULL) {
  src_lines <- if (is.null(builder)) {
    character()
  } else {
    c(
      "",
      "## Dataset generator",
      "",
      "The fixture is generated inside this vignette run. The helper body used",
      "to create the pooled local data is printed here before the analysis.",
      "",
      "```{r dataset-generator}",
      sprintf("cat(paste(deparse(%s), collapse = \"\\n\"))", builder),
      "```"
    )
  }
  c(
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = TRUE, eval = TRUE, warning = FALSE,",
    "  message = TRUE, collapse = TRUE, comment = \"#>\")",
    "helper_candidates <- c(\"validation_helpers.R\",",
    "  file.path(\"..\", \"validation_helpers.R\"),",
    "  file.path(\"vignettes\", \"validation_helpers.R\"))",
    "source(helper_candidates[file.exists(helper_candidates)][1])",
    "validation_load_packages()",
    "partition_summary <- function(tables) {",
    "  data.frame(",
    "    server = names(tables),",
    "    n = vapply(tables, nrow, integer(1)),",
    "    columns = vapply(tables, function(x) paste(names(x), collapse = \", \"), character(1)),",
    "    stringsAsFactors = FALSE)",
    "}",
    "```",
    src_lines
  )
}

write_case <- function(id, K, title, functions, method, math, disclosure,
                       builder, code) {
  file <- file.path("vignettes", sprintf("validation_%s_%s.Rmd", id,
                                         if (K == 2L) "k2" else "kge3"))
  k_label <- if (K == 2L) "K=2" else "K>=3"
  lines <- c(
    "---",
    sprintf("title: \"%s (%s)\"", title, k_label),
    "author: \"David Sarrat Gonzalez\"",
    "date: \"`r format(Sys.Date(), '%Y-%m-%d')`\"",
    "output:",
    "  rmarkdown::html_vignette:",
    "    toc: true",
    "    toc_depth: 2",
    "vignette: >",
    sprintf("  %%\\VignetteIndexEntry{%s (%s)}", title, k_label),
    "  %\\VignetteEngine{knitr::rmarkdown}",
    "  %\\VignetteEncoding{UTF-8}",
    "---",
    "",
    setup_chunk(builder),
    "",
    "## Scope",
    "",
    sprintf("Functions: `%s`.", functions),
    "",
    sprintf("This vignette validates the `%s` modality. It creates the local",
            k_label),
    "pooled fixture, partitions it vertically, starts an in-memory DSLite",
    "server, aligns records with PSI, runs the distributed dsVert route, runs",
    "the centralized R reference on the same pooled fixture, and computes the",
    "numerical delta.",
    "",
    "No external RDS, cache, or pre-loaded server table is used. Rendering the",
    "vignette executes the workflow from dataset generation to assertion.",
    "",
    "## Method",
    "",
    method,
    "",
    "## Mathematical target",
    "",
    math,
    "",
    "## Disclosure review",
    "",
    disclosure,
    "",
    "The DSLite validation uses `datashield.privacyLevel = 5`. Trusted-peer",
    "pinning is disabled only because DSLite is an in-memory test backend, not",
    "an Opal/Rock deployment. That does not weaken the method-level check of",
    "what the product route returns to the analyst.",
    "",
    "## Executed validation",
    "",
    chunk("execution", code),
    "",
    "## Verdict",
    "",
    "Rendering fails if the distributed result leaves the accepted numerical",
    "envelope or if the method row is marked disclosive."
  )
  writeLines(lines, file)
  file
}

codes <- list()

codes$psi <- r"--(
K <- VALIDATION_K

contains_ids <- function(x) {
  if (is.character(x)) return(any(grepl("^P[0-9]{4}$", x)))
  if (is.list(x)) return(any(vapply(x, contains_ids, logical(1))))
  FALSE
}

make_table <- function(ids, server_index) {
  id_num <- as.integer(sub("^P", "", ids))
  ord <- sample(seq_along(ids))
  data.frame(patient_id = ids[ord], id_num = id_num[ord],
             value = id_num[ord] + server_index / 10)
}

set.seed(9100 + K)
id_sets <- if (K == 2L) {
  list(s1 = sprintf("P%04d", 1:80),
       s2 = sprintf("P%04d", 11:90))
} else {
  list(s1 = sprintf("P%04d", 1:85),
       s2 = sprintf("P%04d", 11:95),
       s3 = sprintf("P%04d", 6:80))
}
common_ids <- Reduce(intersect, id_sets)
tables <- Map(make_table, id_sets, seq_along(id_sets))

knitr::kable(partition_summary(tables))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi <- psi_align(conns)
    counts <- DSI::datashield.aggregate(
      conns, call(name = "getObsCountDS", data_name = "DA"))
    n_by_server <- vapply(counts, function(x) as.integer(x$n_obs), integer(1))
    variables <- stats::setNames(lapply(names(conns), function(.x) "id_num"),
                                 names(conns))
    cor_fit <- dsVertClient::ds.vert.cor(
      "DA", variables = variables,
      verbose = validation_demo_verbose(), datasources = conns)
    legacy_blocked <- tryCatch({
      DSI::datashield.aggregate(conns[1],
        call(name = "psiGetMatchedIndicesDS"))
      FALSE
    }, error = function(e) TRUE)
    list(psi = psi, n_by_server = n_by_server, cor_fit = cor_fit,
         legacy_blocked = legacy_blocked)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

offdiag <- result$value$cor_fit$correlation[
  upper.tri(result$value$cor_fit$correlation)]
observed <- max(abs(result$value$n_by_server - length(common_ids)),
                abs((result$value$psi$n_common %||%
                       result$value$psi[[1]]$n_matched) - length(common_ids)),
                max(abs(offdiag - 1)))
audit_ok <- !contains_ids(result$value$psi) &&
  isTRUE(result$value$legacy_blocked)

rows <- row_result(
  "psi", "PSI alignment", K, "ds.vert.align",
  "synthetic ID intersection", "deterministic set intersection",
  "max(count_delta, correlation_delta)", observed, 1e-8,
  "strict-precise",
  "Returns match counts/status only; matched IDs and row indices are not returned.",
  result$runtime_s,
  if (audit_ok) "index reveal blocked" else "audit warning")

display_validation(rows)
assert_validation(rows)
)--"

codes$descriptive <- r"--(
K <- VALIDATION_K

pooled <- build_pima(60L)
tables <- if (K == 2L) {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi", "ped"), drop = FALSE],
    s2 = pooled[, c("patient_id", "npreg", "glu", "bp", "skin"),
                drop = FALSE])
} else {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi"), drop = FALSE],
    s2 = pooled[, c("patient_id", "ped", "skin"), drop = FALSE],
    s3 = pooled[, c("patient_id", "npreg", "glu", "bp"), drop = FALSE])
}
vars <- lapply(tables, function(x) setdiff(names(x), "patient_id"))
all_vars <- unlist(vars, use.names = FALSE)

knitr::kable(utils::head(pooled))
knitr::kable(partition_summary(tables))

ref_mean <- vapply(pooled[all_vars], mean, numeric(1), na.rm = TRUE)
ref_sd <- vapply(pooled[all_vars], stats::sd, numeric(1), na.rm = TRUE)

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.desc(
      "DA", variables = vars, n_buckets = 40L,
      verbose = validation_demo_verbose(), datasources = conns)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

desc <- result$value
ds_mean <- stats::setNames(desc$mean, desc$variable)
ds_sd <- stats::setNames(desc$sd, desc$variable)
observed <- max(max_named_delta(ds_mean, ref_mean),
                max_named_delta(ds_sd, ref_sd))

rows <- row_result(
  "descriptive", "Descriptive statistics", K, "ds.vert.desc",
  "MASS::Pima.tr fixture", "central mean/sd",
  "max(mean_sd_abs_delta)", observed, 1e-8, "strict-precise",
  "Returns guarded scalar summaries and histogram-based quantiles, not rows.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$contingency <- r"--(
K <- VALIDATION_K

pooled <- build_pima(80L)
pooled$diabetes <- factor(ifelse(pooled$diabetes == 1, "yes", "no"),
                          levels = c("no", "yes"))
pooled$age_grp <- factor(ifelse(pooled$age >= median(pooled$age),
                                "older", "younger"),
                         levels = c("younger", "older"))
pooled$bp_grp <- factor(ifelse(pooled$bp >= median(pooled$bp),
                               "high_bp", "low_bp"),
                        levels = c("low_bp", "high_bp"))

tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "age_grp"), drop = FALSE],
       s2 = pooled[, c("patient_id", "bp_grp", "diabetes"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "age_grp"), drop = FALSE],
       s2 = pooled[, c("patient_id", "bp_grp"), drop = FALSE],
       s3 = pooled[, c("patient_id", "diabetes"), drop = FALSE])
}

knitr::kable(utils::head(pooled[, c("patient_id", "age_grp", "bp_grp",
                                    "diabetes")]))
knitr::kable(partition_summary(tables))

tab <- table(pooled$age_grp, pooled$diabetes)
chi_ref <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.chisq_cross(
      "DA", "age_grp", "diabetes", correct = FALSE, fisher = TRUE,
      verbose = validation_demo_verbose(), datasources = conns)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

cross <- result$value
ref_mat <- unclass(tab)[rownames(cross$observed),
                        colnames(cross$observed), drop = FALSE]
observed <- max(max(abs(cross$observed - ref_mat)),
                abs(cross$chisq - unname(chi_ref$statistic)))

rows <- row_result(
  "contingency", "Contingency tests", K,
  "ds.vert.chisq / ds.vert.fisher / ds.vert.chisq_cross",
  "MASS::Pima.tr categorical fixture", "chisq.test / fisher.test",
  "max(count_or_chisq_delta)", observed, 1e-8, "strict-precise",
  "Releases guarded table counts and test statistics; small cells fail closed.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$correlation <- r"--(
K <- VALIDATION_K

pooled <- build_pima(60L)
tables <- if (K == 2L) {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi", "ped"), drop = FALSE],
    s2 = pooled[, c("patient_id", "npreg", "glu", "bp", "skin"),
                drop = FALSE])
} else {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi"), drop = FALSE],
    s2 = pooled[, c("patient_id", "ped", "skin"), drop = FALSE],
    s3 = pooled[, c("patient_id", "npreg", "glu", "bp"), drop = FALSE])
}
vars <- lapply(tables, function(x) setdiff(names(x), "patient_id"))

knitr::kable(utils::head(pooled))
knitr::kable(partition_summary(tables))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.cor(
      "DA", variables = vars,
      verbose = validation_demo_verbose(), datasources = conns)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

cor_ds <- result$value
cor_ref <- stats::cor(pooled[cor_ds$var_names])
observed <- max(abs(cor_ds$correlation - cor_ref))

rows <- row_result(
  "correlation", "Correlation", K, "ds.vert.cor",
  "MASS::Pima.tr fixture", "stats::cor",
  "correlation_max_abs_delta", observed, 1e-4, "strict-practical",
  "Releases the low-dimensional correlation matrix only.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$pca <- r"--(
K <- VALIDATION_K

pooled <- build_pima(60L)
tables <- if (K == 2L) {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi", "ped"), drop = FALSE],
    s2 = pooled[, c("patient_id", "npreg", "glu", "bp", "skin"),
                drop = FALSE])
} else {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi"), drop = FALSE],
    s2 = pooled[, c("patient_id", "ped", "skin"), drop = FALSE],
    s3 = pooled[, c("patient_id", "npreg", "glu", "bp"), drop = FALSE])
}
vars <- lapply(tables, function(x) setdiff(names(x), "patient_id"))

knitr::kable(utils::head(pooled))
knitr::kable(partition_summary(tables))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    cor_ds <- dsVertClient::ds.vert.cor(
      "DA", variables = vars,
      verbose = validation_demo_verbose(), datasources = conns)
    pca_ds <- dsVertClient::ds.vert.pca(
      cor_result = cor_ds,
      verbose = validation_demo_verbose(), datasources = conns)
    list(cor_ds = cor_ds, pca_ds = pca_ds)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

cor_ds <- result$value$cor_ds
pca_ds <- result$value$pca_ds
cor_ref <- stats::cor(pooled[cor_ds$var_names])
eig_ref <- eigen(cor_ref, symmetric = TRUE)
load_ref <- eig_ref$vectors[, seq_len(ncol(pca_ds$loadings)), drop = FALSE]
rownames(load_ref) <- rownames(pca_ds$loadings)
load_ds <- pca_ds$loadings
for (j in seq_len(ncol(load_ds))) {
  if (sum(load_ds[, j] * load_ref[, j]) < 0) load_ds[, j] <- -load_ds[, j]
}
observed <- max(max(abs(pca_ds$eigenvalues - eig_ref$values)),
                max(abs(load_ds - load_ref)))

rows <- row_result(
  "pca", "PCA", K, "ds.vert.pca",
  "MASS::Pima.tr fixture", "eigen(cor(X))",
  "max(eigen_loading_abs_delta)", observed, 1e-4,
  "strict-practical",
  "Reuses the correlation matrix; no score or loading per patient is returned.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$glm <- r"--(
K <- VALIDATION_K

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
knitr::kable(partition_summary(tables))

ref <- coef(stats::lm(fm, pooled))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.glm(
      fm, data = "DA", family = "gaussian", max_iter = 30L,
      verbose = validation_demo_verbose(), datasources = conns)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

fit <- result$value
observed <- max_named_delta(fit$coefficients, ref)

rows <- row_result(
  "glm", "GLM", K, "ds.vert.glm",
  "MASS::Pima.tr fixture", "stats::lm",
  "coef_max_abs_delta", observed, 1e-3, "strict-practical",
  "Uses secure score aggregates and returns model-level coefficients/covariance.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$inference <- r"--(
K <- VALIDATION_K

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
knitr::kable(partition_summary(tables))

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
assert_validation(rows)
)--"

codes$lasso <- r"--(
K <- VALIDATION_K

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
knitr::kable(partition_summary(tables))

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

observed <- max_named_delta(result$value$prox0$coefficients,
                            result$value$fit$coefficients)

rows <- row_result(
  "lasso", "LASSO", K, "ds.vert.lasso_proximal(lambda=0)",
  "MASS::Pima.tr fixture", "OLS limit of ds.vert.glm",
  "lambda0_coef_abs_delta", observed, 1e-8, "strict-precise",
  "Consumes only model-level GLM aggregates; no extra patient-level values are returned.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$negative_binomial <- r"--(
K <- VALIDATION_K

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
knitr::kable(partition_summary(tables))

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
assert_validation(rows)
)--"

codes$cox <- r"--(
K <- VALIDATION_K

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
knitr::kable(partition_summary(tables))

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
assert_validation(rows)
)--"

codes$lmm <- r"--(
K <- VALIDATION_K

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
  "fixed_effect_max_abs_delta", observed, 0.02, "strict-practical",
  "Cluster membership is used internally; per-cluster residuals/BLUPs are not returned.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$gee <- r"--(
K <- VALIDATION_K

validation_require("geepack")
pooled <- build_pima(60L)
pooled$cluster <- as.integer((seq_len(nrow(pooled)) - 1L) %/% 10L) + 1L
tables <- if (K == 2L) {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi", "ped"), drop = FALSE],
    s2 = pooled[, c("patient_id", "npreg", "glu", "bp", "skin",
                    "diabetes", "cluster"), drop = FALSE])
} else {
  list(
    s1 = pooled[, c("patient_id", "age", "bmi"), drop = FALSE],
    s2 = pooled[, c("patient_id", "ped", "skin"), drop = FALSE],
    s3 = pooled[, c("patient_id", "npreg", "glu", "bp",
                    "diabetes", "cluster"), drop = FALSE])
}
fm <- diabetes ~ age + bmi + ped + glu + bp + skin + npreg

knitr::kable(utils::head(pooled))
knitr::kable(partition_summary(tables))

ref <- coef(geepack::geeglm(
  fm, data = pooled, id = cluster, family = gaussian(),
  corstr = "independence"))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.gee(
      fm, data = "DA", family = "gaussian", id_col = "cluster",
      corstr = "independence", max_iter = 30L,
      verbose = validation_demo_verbose(), datasources = conns)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

fit <- result$value
observed <- max_named_delta(fit$coefficients, ref)

rows <- row_result(
  "gee", "GEE", K, "ds.vert.gee(corstr='independence')",
  "MASS::Pima.tr clustered fixture", "geepack::geeglm",
  "coef_max_abs_delta", observed, 0.01, "strict-practical",
  "Returns regression and sandwich-level aggregates, not row scores.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$glmm <- r"--(
K <- VALIDATION_K

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
knitr::kable(partition_summary(tables))

pooled_ref <- pooled
pooled_ref$cluster <- factor(pooled_ref$cluster)
ref_pql <- suppressWarnings(MASS::glmmPQL(
  y ~ x1 + x2, random = ~1 | cluster, family = binomial(),
  data = pooled_ref, verbose = validation_demo_verbose()))
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
assert_validation(rows)
)--"

codes$ipw <- r"--(
K <- VALIDATION_K

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
knitr::kable(partition_summary(tables))

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
assert_validation(rows)
)--"

codes$mi <- r"--(
K <- VALIDATION_K

pooled <- build_mi()
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2", "x3", "y"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "x3", "y"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
knitr::kable(partition_summary(tables))

ref_data <- pooled
ref_data$x2[is.na(ref_data$x2)] <- mean(ref_data$x2, na.rm = TRUE)
ref <- coef(stats::lm(y ~ x1 + x2 + x3, data = ref_data))

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns, na.action = "none")
    dsVertClient::ds.vert.mi(
      y ~ x1 + x2 + x3, data = "DA", impute_columns = "x2",
      family = "gaussian",
      verbose = validation_demo_verbose(), datasources = conns, seed = 12L)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

fit <- result$value
observed <- max_named_delta(fit$coefficients, ref)

rows <- row_result(
  "mi", "Multiple imputation", K, "ds.vert.mi",
  "synthetic missing-covariate fixture",
  "central mean-imputation reference",
  "pooled_coef_abs_delta", observed, 0.02, "strict-practical",
  "Imputed columns stay server-side; client pools beta/covariance only.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$multinomial <- r"--(
K <- VALIDATION_K

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
knitr::kable(partition_summary(tables))

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

fit <- result$value
ds_prob <- softmax_prob(fit$coefficients, pooled)
ds_prob <- ds_prob[, colnames(ref_prob), drop = FALSE]
observed <- max(abs(ds_prob - ref_prob))

rows <- row_result(
  "multinomial", "Multinomial", K, "ds.vert.multinom",
  "balanced soft-signal synthetic 3-class fixture",
  "nnet::multinom probabilities",
  "class_probability_max_abs_delta", observed, 0.005,
  "strict-practical",
  "Softmax probabilities and residuals remain Ring127 shares; no row probabilities are returned.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
)--"

codes$ordinal <- r"--(
K <- VALIDATION_K

validation_require("MASS")
old_fd <- getOption("dsvert.ord_strict_fd_max_dim", NULL)
options(dsvert.ord_strict_fd_max_dim = 0L)

pooled <- build_ordinal()
tables <- if (K == 2L) {
  list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x3", "y_ord",
                       "low_leq", "med_leq", "high_leq"), drop = FALSE])
} else {
  list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
       s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
       s3 = pooled[, c("patient_id", "x3", "y_ord",
                       "low_leq", "med_leq", "high_leq"), drop = FALSE])
}

knitr::kable(utils::head(pooled))
knitr::kable(partition_summary(tables))

ref <- MASS::polr(y_ord ~ x1 + x2 + x3, data = pooled, Hess = TRUE)
ref_cum <- sapply(ref$zeta, function(th) {
  stats::plogis(th - as.numeric(as.matrix(pooled[, names(coef(ref))]) %*%
                                  coef(ref)))
})

result <- elapsed({
  conns <- connect_dslite(tables)
  tryCatch({
    psi_align(conns)
    dsVertClient::ds.vert.ordinal(
      y_ord ~ x1 + x2 + x3, data = "DA",
      levels_ordered = c("low", "med", "high"),
      cumulative_template = "%s_leq",
      verbose = validation_demo_verbose(), datasources = conns)
  }, finally = try(DSI::datashield.logout(conns), silent = TRUE))
})

fit <- result$value
ds_cum <- ordinal_cumprob(fit, pooled)
colnames(ds_cum) <- colnames(ref_cum)
observed <- max(abs(ds_cum - ref_cum))

rows <- row_result(
  "ordinal", "Ordinal", K, "ds.vert.ordinal",
  "balanced synthetic 3-level ordinal fixture",
  "MASS::polr cumulative probabilities",
  "cumulative_probability_max_abs_delta", observed, 0.001,
  "strict-practical",
  "Class probabilities/residuals stay Ring127 shares; no row probabilities are returned.",
  result$runtime_s)

display_validation(rows)
assert_validation(rows)
options(dsvert.ord_strict_fd_max_dim = old_fd)
)--"

meta <- list(
  psi = list(
    title = "PSI alignment validation",
    functions = "ds.vert.align()",
    builder = NULL,
    method = "ECDH-PSI aligns vertically split tables on a shared identifier. Each server masks identifiers, the protocol computes the common set, and the aligned object keeps only matched rows in a common order.",
    math = "The target is the set intersection I = intersection_k I_k. Validation checks matched counts and that deterministic id-derived columns align perfectly after PSI.",
    disclosure = "The analyst receives counts and status fields only. Patient IDs and matched row-index vectors are not returned; the legacy matched-index aggregate is checked as blocked."
  ),
  descriptive = list(
    title = "Descriptive statistics validation",
    functions = "ds.vert.desc()",
    builder = "build_pima",
    method = "Each server computes guarded scalar summaries for its own variables. Histogram quantiles are guarded by disclosure thresholds, but this vignette validates the exact mean/sd surface used in the matrix.",
    math = "For each variable x, the checked target is max(|mean_DS - mean_local|, |sd_DS - sd_local|).",
    disclosure = "Only scalar summaries and guarded histogram information are returned; no patient-level values or order statistics are disclosed."
  ),
  contingency = list(
    title = "Contingency tests validation",
    functions = "ds.vert.chisq_cross()",
    builder = "build_pima",
    method = "Cross-server contingency testing builds guarded cell counts across vertically separated categorical variables before applying the usual chi-square/Fisher calculations.",
    math = "The validation compares the released observed table and X^2 statistic with table() and chisq.test() on the pooled fixture.",
    disclosure = "Exact counts are released only after small-cell checks pass. Sparse positive cells fail closed before release."
  ),
  correlation = list(
    title = "Correlation validation",
    functions = "ds.vert.cor()",
    builder = "build_pima",
    method = "The method combines server-local moments and MPC cross-products to release a guarded low-dimensional correlation matrix.",
    math = "For variables x_j and x_l, r_jl = cov(x_j, x_l) / (sd(x_j) sd(x_l)). Validation compares the full matrix to stats::cor().",
    disclosure = "The disclosure surface is the aggregate p by p correlation matrix, with no row-level cross-products returned."
  ),
  pca = list(
    title = "PCA validation",
    functions = "ds.vert.pca()",
    builder = "build_pima",
    method = "PCA is computed from the federated correlation matrix, not from individual-level scores.",
    math = "The target is R = V Lambda V^T. Validation compares eigenvalues and sign-aligned loadings with eigen(cor(X)).",
    disclosure = "Only eigenvalues and loadings derived from the aggregate correlation matrix are returned; no individual PCA scores are released."
  ),
  glm = list(
    title = "GLM validation",
    functions = "ds.vert.glm()",
    builder = "build_pima",
    method = "The public frontdoor routes to the protected GLM engine. This vignette uses the Gaussian fixture from the matrix, where the centralized target is exact least squares.",
    math = "For the Gaussian identity-link case, beta solves X'X beta = X'y. The distributed route reconstructs only aggregate score/Hessian information needed for beta.",
    disclosure = "The analyst receives coefficients and model-level covariance/diagnostics. Eta, residuals, row scores, and fitted values are not returned."
  ),
  inference = list(
    title = "GLM inference validation",
    functions = "ds.vert.confint(), ds.vert.wald(), ds.vert.contrast(), ds.vert.lr()",
    builder = "build_pima",
    method = "Inference helpers post-process the model-level beta/covariance/deviance already returned by ds.vert.glm().",
    math = "Wald and contrast tests use K beta and K V K'. LR tests use the deviance difference with a chi-square reference.",
    disclosure = "No extra server aggregate is opened beyond the fitted GLM objects."
  ),
  lasso = list(
    title = "LASSO validation",
    functions = "ds.vert.lasso_proximal()",
    builder = "build_pima",
    method = "The proximal LASSO path consumes the distributed GLM fit. The validation uses lambda = 0, where the penalized solution must equal the prime GLM estimate.",
    math = "At lambda = 0, the proximal quadratic argmin is beta_hat, so the checked target is max |beta_lasso0 - beta_glm|.",
    disclosure = "The proximal step uses only model-level beta/Hessian information and returns coefficients; it does not request row-level data."
  ),
  negative_binomial = list(
    title = "Negative binomial validation",
    functions = "ds.vert.nb()",
    builder = "build_nb",
    method = "This vignette shows the accurate non-disclosive full-regression route and the fast MoM route. The accurate route targets MASS::glm.nb(); the fast route targets a central Poisson beta plus closed-form MoM theta approximation.",
    math = "NB2 uses Var(Y|X)=mu+mu^2/theta and log(mu)=X beta. The accurate route evaluates beta/theta score pieces in shares; the fast route uses theta_MoM = ybar^2/(s_y^2-ybar).",
    disclosure = "The accurate route keeps eta, mu, reciprocal, score, and Fisher pieces in shares and returns beta plus scalar theta. The fast route opens only outcome-level scalar moments and Poisson-model aggregates."
  ),
  cox = list(
    title = "Cox PH validation",
    functions = "ds.vert.cox()",
    builder = "build_cox",
    method = "The product Cox route is the non-disclosive Breslow profile route. The discrete-time route is a different estimand, not the fast approximation to this profile target, so this vignette validates the route used in the product matrix.",
    math = "The Cox PH score is sum_i delta_i [x_i - E_beta{x | t_i in risk set}]. Validation compares beta to survival::coxph(ties='breslow').",
    disclosure = "Risk-set/event-time score terms remain shared. The client receives coefficients and scalar convergence diagnostics, not event-time vectors or patient risk-set membership."
  ),
  lmm = list(
    title = "LMM validation",
    functions = "ds.vert.lmm()",
    builder = "build_lmm",
    method = "Random-intercept LMM uses guarded cluster metadata and protected cross-products to fit fixed effects and variance components.",
    math = "The model is y = X beta + b_cluster + epsilon with b ~ N(0,sigma_b^2). Validation compares fixed effects to lme4::lmer().",
    disclosure = "Cluster membership is used internally at the accepted method tier. Per-cluster residuals and BLUP vectors are not returned."
  ),
  gee = list(
    title = "GEE validation",
    functions = "ds.vert.gee()",
    builder = "build_pima",
    method = "The GEE route computes protected estimating-equation aggregates and robust sandwich-level summaries for clustered data.",
    math = "For independence working correlation, sum_i D_i' V_i^{-1}(y_i-mu_i)=0. Validation compares beta to geepack::geeglm().",
    disclosure = "The client receives regression and sandwich-level aggregates, not row scores, fitted values, or residual vectors."
  ),
  glmm = list(
    title = "GLMM validation",
    functions = "ds.vert.glmm()",
    builder = "build_balanced_glmm",
    method = "This vignette shows both available product routes: PQL as the lower-cost approximation and Laplace as the more accurate route against lme4::glmer.",
    math = "PQL alternates binomial working responses with aggregate mixed-model equations. Laplace approximates the random-intercept marginal likelihood and checks both fixed effects and scalar random-effect variance.",
    disclosure = "Both routes return fixed effects and scalar variance/quality diagnostics only. Per-patient probabilities, row scores, BLUPs, cluster labels, and cluster vectors are not returned."
  ),
  ipw = list(
    title = "IPW validation",
    functions = "ds.vert.ipw()",
    builder = "build_ipw",
    method = "IPW runs a protected propensity GLM and then a protected weighted outcome GLM using a server-side weight column.",
    math = "The checked outcome target solves X'W(y-X beta)=0 with the same IPW weights used by the centralized lm().",
    disclosure = "Weights remain server-side; the analyst receives only propensity and outcome model-level fits."
  ),
  mi = list(
    title = "Multiple imputation validation",
    functions = "ds.vert.mi()",
    builder = "build_mi",
    method = "Missing covariates are imputed server-side for each round, then GLM fits are pooled client-side.",
    math = "The compact reference uses central mean imputation, and the distributed route is checked by pooled coefficient delta.",
    disclosure = "Imputed columns stay server-side. The client receives pooled beta/covariance summaries, not imputed patient-level values."
  ),
  multinomial = list(
    title = "Multinomial validation",
    functions = "ds.vert.multinom()",
    builder = "build_multinomial",
    method = "The supported route is joint softmax Newton. Historical one-vs-rest approximations are not exposed as the product evidence path.",
    math = "P(Y=c|X)=exp(eta_c)/sum_l exp(eta_l). Validation compares class probabilities with nnet::multinom().",
    disclosure = "Softmax probabilities, residuals, and row scores remain Ring127 shares; no row probabilities are returned."
  ),
  ordinal = list(
    title = "Ordinal validation",
    functions = "ds.vert.ordinal()",
    builder = "build_ordinal",
    method = "The supported route is joint proportional-odds Newton. Historical cumulative-binomial approximations are only warm starts.",
    math = "P(Y<=k|X)=sigmoid(theta_k-X beta). Validation compares cumulative probabilities with MASS::polr().",
    disclosure = "Cumulative probabilities, residuals, and row scores stay in shares; the client receives thresholds, slopes, and scalar optimizer diagnostics."
  )
)

order_ids <- c("psi", "descriptive", "contingency", "correlation", "pca",
               "glm", "inference", "lasso", "negative_binomial", "cox",
               "lmm", "gee", "glmm", "ipw", "mi", "multinomial", "ordinal")

files <- character()
for (id in order_ids) {
  for (K in c(2L, 3L)) {
    code <- gsub("VALIDATION_K", paste0(K, "L"), codes[[id]], fixed = TRUE)
    files <- c(files, write_case(
      id = id,
      K = K,
      title = meta[[id]]$title,
      functions = meta[[id]]$functions,
      method = meta[[id]]$method,
      math = meta[[id]]$math,
      disclosure = meta[[id]]$disclosure,
      builder = meta[[id]]$builder,
      code = code))
  }
}

cat("Wrote", length(files), "validation vignettes\n")
cat(paste(files, collapse = "\n"), "\n")
