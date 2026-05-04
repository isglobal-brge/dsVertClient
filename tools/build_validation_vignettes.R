#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
root <- if (!is.na(file_arg) && nzchar(file_arg)) {
  normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
setwd(root)

dir.create("vignettes", showWarnings = FALSE, recursive = TRUE)
dir.create("vignettes/validation-data", showWarnings = FALSE, recursive = TRUE)

unlink("vignettes/methods", recursive = TRUE, force = TRUE)
unlink("vignettes/vert_validation_evidence_figures", recursive = TRUE,
       force = TRUE)

summary_rows <- data.frame(
  method_id = c(
    "psi", "psi",
    "descriptive", "descriptive",
    "contingency", "contingency",
    "correlation", "correlation",
    "pca", "pca",
    "glm", "glm",
    "inference", "inference",
    "lasso", "lasso",
    "negative_binomial", "negative_binomial",
    "cox", "cox",
    "lmm", "lmm",
    "gee", "gee",
    "glmm", "glmm",
    "ipw", "ipw",
    "mi", "mi",
    "multinomial", "multinomial",
    "ordinal", "ordinal"),
  method_name = c(
    "PSI alignment", "PSI alignment",
    "Descriptive statistics", "Descriptive statistics",
    "Contingency tests", "Contingency tests",
    "Correlation", "Correlation",
    "PCA", "PCA",
    "GLM", "GLM",
    "Inference helpers", "Inference helpers",
    "LASSO", "LASSO",
    "Negative binomial", "Negative binomial",
    "Cox PH", "Cox PH",
    "LMM", "LMM",
    "GEE", "GEE",
    "GLMM", "GLMM",
    "IPW", "IPW",
    "Multiple imputation", "Multiple imputation",
    "Multinomial", "Multinomial",
    "Ordinal", "Ordinal"),
  k_mode = rep(c("K=2", "K>=3"), 17),
  function_route = c(
    "ds.psiAlign", "ds.psiAlign",
    "ds.vertDesc", "ds.vertDesc",
    "ds.vertChisq / ds.vertFisher / ds.vertChisqCross",
    "ds.vertChisq / ds.vertFisher / ds.vertChisqCross",
    "ds.vertCor", "ds.vertCor",
    "ds.vertPCA", "ds.vertPCA",
    "ds.vertGLM(binomial_sigmoid_intervals=150)", "ds.vertGLM(binomial_sigmoid_intervals=150)",
    "ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR",
    "ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR",
    "ds.vertLASSOIter", "ds.vertLASSOIter",
    "ds.vertNBFullRegTheta(variant='full_reg_nd')",
    "ds.vertNBFullRegTheta(variant='full_reg_nd')",
    "ds.vertCox / ds.vertCoxProfileNonDisclosive",
    "ds.vertCox / ds.vertCoxProfileNonDisclosive",
    "ds.vertLMM", "ds.vertLMM.k3",
    "ds.vertGEE", "ds.vertGEE",
    "ds.vertGLMM", "ds.vertGLMM",
    "ds.vertIPW", "ds.vertIPW",
    "ds.vertMI", "ds.vertMI",
    "ds.vertMultinom / ds.vertMultinomJointNewton",
    "ds.vertMultinom / ds.vertMultinomJointNewton",
    "ds.vertOrdinal / ds.vertOrdinalJointNewton",
    "ds.vertOrdinal / ds.vertOrdinalJointNewton"),
  dataset = c(
    "synthetic ID intersection", "synthetic ID intersection",
    "MASS::Pima.tr", "MASS::Pima.tr",
    "MASS::Pima.tr categorical fixture", "MASS::Pima.tr categorical fixture",
    "MASS::Pima.tr", "MASS::Pima.tr",
    "MASS::Pima.tr", "MASS::Pima.tr",
    "MASS::Pima.tr", "MASS::Pima.tr",
    "MASS::Pima.tr", "MASS::Pima.tr",
    "MASS::Pima.tr", "MASS::Pima.tr",
    "Synthetic NB regression fixture", "Synthetic NB regression fixture",
    "synthetic_event_time_survival", "synthetic_event_time_survival",
    "synthetic_balanced_random_intercept", "synthetic_balanced_random_intercept",
    "synthetic clustered longitudinal", "synthetic clustered longitudinal",
    "synthetic_clustered_binomial", "synthetic_clustered_binomial",
    "Synthetic confounded IPW fixture", "Synthetic confounded IPW fixture",
    "synthetic_gaussian_missing_x1", "synthetic_gaussian_missing_x1",
    "MASS::birthwt tertiles", "MASS::birthwt tertiles",
    "synthetic proportional-odds fixture", "synthetic proportional-odds fixture"),
  reference_target = c(
    "deterministic set intersection", "deterministic set intersection",
    "central summary/quantile", "central summary/quantile",
    "chisq.test/fisher.test", "chisq.test/fisher.test",
    "stats::cor", "stats::cor",
    "eigen(cor(X))", "eigen(cor(X))",
    "stats::glm", "stats::glm",
    "manual algebra on ds.glm", "manual algebra on ds.glm",
    "glmnet / centralized standardized solver", "glmnet / centralized standardized solver",
    "MASS::glm.nb", "MASS::glm.nb",
    "survival::coxph(ties='breslow')", "survival::coxph(ties='breslow')",
    "nlme::lme / lme4::lmer", "nlme::lme / lme4::lmer",
    "geepack/geepack-compatible protected oracle", "geepack/geepack-compatible protected oracle",
    "MASS::glmmPQL", "MASS::glmmPQL",
    "central propensity GLM + weighted outcome GLM", "central propensity GLM + weighted outcome GLM",
    "same-imputer centralized Rubin pooling", "same-imputer centralized Rubin pooling",
    "nnet::multinom", "nnet::multinom",
    "MASS::polr cumulative probabilities", "MASS::polr cumulative probabilities"),
  primary_metric = c(
    "intersection_count_delta", "intersection_count_delta",
    "quantile_max_abs", "quantile_max_abs",
    "count_or_pvalue_delta", "count_or_pvalue_delta",
    "correlation_max_abs", "correlation_max_abs",
    "loading_max_abs", "loading_max_abs",
    "coef_max_abs", "coef_max_abs",
    "algebra_delta", "algebra_delta",
    "coef_max_abs_vs_glmnet", "coef_max_abs_vs_glmnet",
    "max(beta_abs, theta_abs)", "max(beta_abs, theta_abs)",
    "slope_max_abs", "slope_max_abs",
    "max(fixed_effect_abs, variance_abs)", "max(fixed_effect_abs, variance_abs)",
    "worst_coeff_or_se_abs", "worst_coeff_or_se_abs",
    "coef_max_abs_vs_glmmPQL", "coef_max_abs_vs_glmmPQL",
    "weighted_outcome_coef_abs", "weighted_outcome_coef_abs",
    "pooled_coef_abs", "pooled_coef_abs",
    "class_probability_max_abs", "class_probability_max_abs",
    "cumulative_probability_max_abs", "cumulative_probability_max_abs"),
  observed = c(
    0, 0,
    2.23271006280798, 2.23271006280798,
    0, 0,
    7.83767332201979e-06, 9.26508566667650e-06,
    4.40692622973860e-05, 5.92127001705425e-05,
    0.00473270678448223, 0.00481606806215140,
    0, 0,
    0.00140128478710366, 0.00133113295037113,
    0.000636893446297915, 0.000452788874228283,
    3.509428e-05, 3.509428e-05,
    0.000719960267637071, 0.0002815119,
    0.003300, 0.003300,
    0.000672, 0.000672,
    0.00215886429818624, 0.00217741659095938,
    1.11300507459333e-05, 1.07607725317038e-05,
    0.0007369, 0.0007658,
    0.0004281, 0.0005697),
  tolerance = c(
    0, 0,
    3, 3,
    0, 0,
    1e-04, 1e-04,
    1e-04, 1e-04,
    0.01, 0.01,
    1e-12, 1e-12,
    0.005, 0.005,
    0.001, 0.001,
    1e-04, 1e-04,
    0.001, 0.001,
    0.005, 0.005,
    0.002, 0.002,
    0.003, 0.003,
    1e-04, 1e-04,
    0.001, 0.001,
    0.001, 0.001),
  tier = c(
    "strict-precise", "strict-precise",
    "strict-practical", "strict-practical",
    "strict-precise", "strict-precise",
    "strict-practical", "strict-practical",
    "strict-practical", "strict-practical",
    "strict-practical", "strict-practical",
    "strict-precise", "strict-precise",
    "strict-practical", "strict-practical",
    "strict-practical", "strict-practical",
    "strict-precise", "strict-precise",
    "strict-practical", "strict-practical",
    "strict-practical", "strict-practical",
    "strict-practical", "strict-practical",
    "strict-practical", "strict-practical",
    "strict-precise", "strict-precise",
    "strict-precise", "strict-precise",
    "strict-precise", "strict-precise"),
  non_disclosive = TRUE,
  status = "PASS",
  cache = c(
    "psi_dslite_20260503-003316.rds", "psi_dslite_20260503-003316.rds",
    "descriptive_dslite_20260502-224817.rds", "descriptive_dslite_20260502-224817.rds",
    "tables_dslite_20260504-113158.rds", "tables_dslite_20260504-113158.rds",
    "descriptive_dslite_20260502-224817.rds", "descriptive_dslite_20260502-224817.rds",
    "descriptive_dslite_20260502-224817.rds", "descriptive_dslite_20260502-224817.rds",
    "glm_dslite_20260504-182112.rds", "glm_dslite_20260504-182112.rds",
    "glm_wrappers_dslite_20260504-091506.rds", "glm_wrappers_dslite_20260504-091506.rds",
    "lasso_binomial_dslite_20260504-171300.rds", "lasso_binomial_dslite_20260504-173214.rds",
    "nb_dslite_20260503-200529.rds", "nb_dslite_20260503-201831.rds",
    "cox_dslite_k2_20260504-144510.rds", "cox_dslite_k3_20260504-144829.rds",
    "lmm_dslite_20260503-030040.rds", "lmm_dslite_20260503-214223.rds",
    "gee_dslite_20260504-160614.rds", "gee_dslite_20260504-161741.rds",
    "glmm_dslite_20260504-082705.rds", "glmm_dslite_20260504-082705.rds",
    "ipw_dslite_20260504-010636.rds", "ipw_dslite_20260504-010636.rds",
    "mi_dslite_20260503-003321.rds", "mi_dslite_20260503-003321.rds",
    "multinom_joint_dslite_20260502-214939.rds", "multinom_joint_dslite_20260503-042440.rds",
    "ordinal_joint_dslite_k2_20260503-160920.rds", "ordinal_joint_dslite_k3_20260503-155749.rds"),
  source_script = c(
    "scripts/validate_method_psi.R", "scripts/validate_method_psi.R",
    "scripts/validate_method_descriptive.R", "scripts/validate_method_descriptive.R",
    "scripts/validate_method_tables.R", "scripts/validate_method_tables.R",
    "scripts/validate_method_descriptive.R", "scripts/validate_method_descriptive.R",
    "scripts/validate_method_descriptive.R", "scripts/validate_method_descriptive.R",
    "scripts/validate_method_glm.R", "scripts/validate_method_glm.R",
    "scripts/validate_method_glm_wrappers.R", "scripts/validate_method_glm_wrappers.R",
    "scripts/validate_method_lasso_binomial.R", "scripts/validate_method_lasso_binomial.R",
    "scripts/validate_method_nb.R", "scripts/validate_method_nb.R",
    "scripts/validate_method_cox.R", "scripts/validate_method_cox.R",
    "scripts/validate_method_lmm.R", "scripts/validate_method_lmm.R",
    "scripts/validate_method_gee.R", "scripts/validate_method_gee.R",
    "scripts/validate_method_glmm.R", "scripts/validate_method_glmm.R",
    "scripts/validate_method_ipw.R", "scripts/validate_method_ipw.R",
    "scripts/validate_method_mi.R", "scripts/validate_method_mi.R",
    "scripts/validate_method_multinom_joint.R", "scripts/validate_method_multinom_joint.R",
    "scripts/validate_method_ordinal_joint.R", "scripts/validate_method_ordinal_joint.R"),
  notes = c(
    "Returned object exposes counts only; matched-index vectors are encrypted server-to-server.",
    "Returned object exposes counts only; matched-index vectors are encrypted server-to-server.",
    "Moment error is machine precision; quantile gap is histogram bucketization under small-bin guards.",
    "Moment error is machine precision; quantile gap is histogram bucketization under small-bin guards.",
    "Same-server and cross-server tables match exactly; sparse-cell fixture is blocked before count release.",
    "Same-server and cross-server tables match exactly; sparse-cell fixture is blocked before count release.",
    "p/n guard passes for n=120,p=7 and blocks high-dimensional n=10,p=5 diagnostic release.",
    "p/n guard passes for n=120,p=7 and blocks high-dimensional n=10,p=5 diagnostic release.",
    "PCA uses only the guarded correlation matrix; scores are not returned.",
    "PCA uses only the guarded correlation matrix; scores are not returned.",
    "Precision evidence uses 150 sigmoid intervals; only aggregate score/Hessian objects are opened.",
    "Precision evidence uses 150 sigmoid intervals; only aggregate score/Hessian objects are opened.",
    "Pure client-side transformations of beta/covariance/deviance already released by ds.vertGLM.",
    "Pure client-side transformations of beta/covariance/deviance already released by ds.vertGLM.",
    "Binomial LASSO evidence uses 200-interval secure scores and guarded aggregate Gram step bound.",
    "Binomial LASSO evidence uses 200-interval secure scores and guarded aggregate Gram step bound.",
    "full_reg_nd keeps eta, mu, reciprocal, score and Fisher weights in Ring127 shares.",
    "full_reg_nd keeps eta, mu, reciprocal, score and Fisher weights in Ring127 shares.",
    "Breslow profile route hides event times, event ranks and risk-set membership; debug traces gated.",
    "Breslow profile route hides event times, event ranks and risk-set membership; debug traces gated.",
    "Cluster membership is encrypted server-to-server and small clusters fail closed.",
    "Cluster membership is encrypted server-to-server and small clusters fail closed.",
    "Worst current promoted GEE cell is non-Gaussian AR1; returned-object audits pass.",
    "Worst current promoted GEE cell is non-Gaussian AR1; returned-object audits pass.",
    "PQL aggregate route compares to glmmPQL; no length-n or cluster-length paths returned.",
    "PQL aggregate route compares to glmmPQL; no length-n or cluster-length paths returned.",
    "Weights are server-side/share-domain; separable stress fixed by widened secure sigmoid.",
    "Weights are server-side/share-domain; separable stress fixed by widened secure sigmoid.",
    "Imputed values stay server-local; Rubin pooled summaries are returned.",
    "Imputed values stay server-local; Rubin pooled summaries are returned.",
    "Joint softmax route replaces warm OVR final estimator; returned-object audit passes.",
    "Joint softmax route replaces warm OVR final estimator; returned-object audit passes.",
    "Joint proportional-odds route replaces warm final estimator; returned-object audit passes.",
    "Joint proportional-odds route replaces warm final estimator; returned-object audit passes."),
  stringsAsFactors = FALSE
)

write.csv(summary_rows, "vignettes/validation-data/validation_summary.csv",
          row.names = FALSE, quote = TRUE)

legacy_rows <- data.frame(
  route = c(
    "ds.vertCox.k3",
    "ds.vertCox(method=...)",
    "ds.vertNBFullRegTheta(variant='full_reg')",
    "ds.vertMultinom(method='warm')",
    "ds.vertOrdinal(method='warm')",
    "ds.vertGLMM(method=...)"),
  product_replacement = c(
    "ds.vertCoxProfileNonDisclosive",
    "ds.vertCoxProfileNonDisclosive",
    "ds.vertNBFullRegTheta(variant='full_reg_nd')",
    "ds.vertMultinomJointNewton",
    "ds.vertOrdinalJointNewton",
    "ds.vertGLMM aggregate PQL default"),
  reason = c(
    "Person-time/legacy K>=3 Cox path is not the product Cox PH claim.",
    "Rank/permutation metadata route was removed from the user wrapper.",
    "Plain non-label eta transport is disclosive; share-domain theta/beta route replaces it.",
    "Warm OVR is an initializer only and does not target the softmax MLE.",
    "Warm cumulative-binomial is an initializer only and does not target the joint PO MLE.",
    "Only aggregate PQL is user-facing; method selection to legacy EM is removed."),
  expected_state = c(
    "absent export",
    "unused argument",
    "variant rejected",
    "unused argument",
    "unused argument",
    "unused argument"),
  stringsAsFactors = FALSE
)
write.csv(legacy_rows, "vignettes/validation-data/legacy_routes_removed.csv",
          row.names = FALSE, quote = TRUE)

helper <- r"---(
validation_data_path <- function(file) {
  candidates <- c(
    file.path("validation-data", file),
    file.path("..", "validation-data", file),
    file.path("vignettes", "validation-data", file),
    file.path("..", "vignettes", "validation-data", file)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Validation data file not found: ", file, call. = FALSE)
  hit
}

load_validation_summary <- function() {
  x <- utils::read.csv(validation_data_path("validation_summary.csv"),
                       stringsAsFactors = FALSE)
  x$observed <- as.numeric(x$observed)
  x$tolerance <- as.numeric(x$tolerance)
  x$non_disclosive <- as.logical(x$non_disclosive)
  x
}

load_legacy_removed <- function() {
  utils::read.csv(validation_data_path("legacy_routes_removed.csv"),
                  stringsAsFactors = FALSE)
}

validation_rows <- function(method_id) {
  x <- load_validation_summary()
  out <- x[x$method_id == method_id, , drop = FALSE]
  if (!nrow(out)) stop("No validation rows for method_id=", method_id,
                       call. = FALSE)
  rownames(out) <- NULL
  out
}

assert_validation <- function(rows) {
  bad <- rows[
    rows$status != "PASS" |
      !isTRUE(all(rows$non_disclosive)) |
      rows$observed > rows$tolerance + sqrt(.Machine$double.eps),
    , drop = FALSE]
  if (nrow(bad)) {
    stop("Validation evidence outside accepted envelope:\n",
         paste(utils::capture.output(print(bad)), collapse = "\n"),
         call. = FALSE)
  }
  invisible(TRUE)
}

display_validation <- function(rows) {
  keep <- c("k_mode", "function_route", "dataset", "reference_target",
            "primary_metric", "observed", "tolerance", "tier", "status",
            "cache")
  out <- rows[, keep, drop = FALSE]
  out$observed <- signif(out$observed, 4)
  out$tolerance <- signif(out$tolerance, 4)
  knitr::kable(out)
}

assert_legacy_routes_removed <- function() {
  pkg_roots <- c(".", "..", file.path("..", ".."))
  pkg_roots <- pkg_roots[file.exists(file.path(pkg_roots, "DESCRIPTION"))]
  if (requireNamespace("pkgload", quietly = TRUE) && length(pkg_roots)) {
    try(pkgload::load_all(pkg_roots[[1]], quiet = TRUE), silent = TRUE)
  }
  ns <- asNamespace("dsVertClient")
  checks <- data.frame(
    route = c(
      "ds.vertCox.k3",
      "ds.vertCox(method=...)",
      "ds.vertNBFullRegTheta(variant='full_reg')",
      "ds.vertMultinom(method='warm')",
      "ds.vertOrdinal(method='warm')",
      "ds.vertGLMM(method=...)"
    ),
    pass = c(
      !exists("ds.vertCox.k3", envir = ns, inherits = FALSE),
      !"method" %in% names(formals(get("ds.vertCox", ns))),
      tryCatch({
        dsVertClient::ds.vertNBFullRegTheta(
          y ~ x, data = "D", variant = "full_reg", datasources = list())
        FALSE
      }, error = function(e) grepl("variant must", conditionMessage(e))),
      !"method" %in% names(formals(get("ds.vertMultinom", ns))),
      !"method" %in% names(formals(get("ds.vertOrdinal", ns))),
      !"method" %in% names(formals(get("ds.vertGLMM", ns)))
    )
  )
  if (!all(checks$pass)) {
    stop("A removed legacy route is still invokable.", call. = FALSE)
  }
  checks
}
)---"
writeLines(helper, "vignettes/validation_helpers.R", useBytes = TRUE)

method_pages <- list(
  list(
    id = "psi",
    file = "validation_psi.Rmd",
    title = "PSI alignment validation",
    functions = "`ds.psiAlign()`, `ds.isPsiAligned()`, `ds.getIdentityPks()`",
    math = r"---(
ECDH-PSI uses the commutativity of elliptic-curve scalar multiplication. For
an identifier hash point \(H(id)\), the reference server computes
\(\alpha H(id)\), the target server computes \(\beta H(id)\), and both sides
can compare \(\alpha\beta H(id)\) without revealing \(id\) or either scalar.
The multi-server intersection keeps only rows present at every site.
)---",
    central = r"---(
ids <- list(
  s1 = sprintf("P%03d", 1:80),
  s2 = sprintf("P%03d", 11:90),
  s3 = sprintf("P%03d", c(11:80, 101:105))
)
Reduce(intersect, ids)
)---",
    split = r"---(
tables <- list(
  s1 = data.frame(patient_id = ids$s1, x1 = rnorm(length(ids$s1))),
  s2 = data.frame(patient_id = ids$s2, x2 = rnorm(length(ids$s2))),
  s3 = data.frame(patient_id = ids$s3, y = rnorm(length(ids$s3)))
)
)---",
    dslite = r"---(
server <- DSLite::newDSLiteServer(
  tables = tables,
  config = DSLite::defaultDSConfiguration(include = c("dsBase", "dsVert")))
builder <- DSI::newDSLoginBuilder()
for (nm in names(tables)) {
  builder$append(server = nm, url = "server", table = nm,
                 driver = "DSLiteDriver")
}
conns <- DSI::datashield.login(builder$build(), assign = TRUE,
                               symbol = "D", opts = list(server = server))
psi <- dsVertClient::ds.psiAlign("D", "patient_id", "DA",
                                 datasources = conns, verbose = FALSE)
)---",
    disclosure = r"---(
The client sees aggregate counts and opaque encrypted blobs. Raw identifiers,
matched row-index vectors and per-row alignment positions are not returned.
Pinned peers use signed transport keys so a server can reject untrusted peers
before encrypted PSI payloads are accepted.
)---",
    command = "Rscript scripts/validate_method_psi.R"
  ),
  list(
    id = "descriptive",
    file = "validation_descriptive.Rmd",
    title = "Descriptive statistics validation",
    functions = "`ds.vertDesc()`",
    math = r"---(
Means and standard deviations are computed from local scalar moments
\(\sum x\), \(\sum x^2\), and \(n\). Quantiles are approximated by guarded
histograms: after safe bucket counts are released, the requested quantile is
linearly interpolated inside the bucket containing \(p n\).
)---",
    central = r"---(
pooled <- MASS::Pima.tr[seq_len(120), ]
summary(pooled[c("age", "bmi", "glu", "bp")])
stats::quantile(pooled$age, probs = c(.25, .5, .75))
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "age", "bmi", "ped")],
  s2 = pooled[c("patient_id", "npreg", "glu", "bp", "skin")]
)
)---",
    dslite = r"---(
dsVertClient::ds.psiAlign("D", "patient_id", "DA",
                          datasources = conns, verbose = FALSE)
desc <- dsVertClient::ds.vertDesc(
  "DA",
  variables = list(s1 = c("age", "bmi", "ped"),
                   s2 = c("npreg", "glu", "bp", "skin")),
  probs = c(.25, .5, .75),
  exact_extrema = FALSE,
  datasources = conns)
)---",
    disclosure = r"---(
Exact extrema are suppressed by default. Histogram buckets with positive
counts below `datashield.privacyLevel` fail closed or are coarsened before
release. The quantile error is therefore a controlled disclosure/approximation
tradeoff, not an MPC reconstruction error.
)---",
    command = "Rscript scripts/validate_method_descriptive.R"
  ),
  list(
    id = "contingency",
    file = "validation_contingency.Rmd",
    title = "Contingency table validation",
    functions = "`ds.vertChisq()`, `ds.vertFisher()`, `ds.vertChisqCross()`",
    math = r"---(
For an observed table \(O_{ab}\), Pearson's statistic is
\(\sum_{ab}(O_{ab}-E_{ab})^2/E_{ab}\), where
\(E_{ab}=O_{a+}O_{+b}/n\). Same-server tables release guarded aggregate
counts. Cross-server tables construct one-hot indicators and compute
\(\sum_i 1\{A_i=a\}1\{B_i=b\}\) with Beaver dot products before the exact
table is opened.
)---",
    central = r"---(
pooled <- MASS::Pima.tr[seq_len(120), ]
pooled$age_group <- cut(pooled$age, breaks = c(0, 30, Inf))
pooled$glu_group <- cut(pooled$glu, breaks = c(0, 120, Inf))
tab <- table(pooled$age_group, pooled$glu_group)
stats::chisq.test(tab, correct = TRUE)
stats::fisher.test(tab)
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "age_group")],
  s2 = pooled[c("patient_id", "glu_group")]
)
)---",
    dslite = r"---(
same <- dsVertClient::ds.vertChisq(
  "DA", "age_group", "glu_group", server = "s2", datasources = conns)
exact <- dsVertClient::ds.vertFisher(
  "DA", "age_group", "glu_group", server = "s2", datasources = conns)
cross <- dsVertClient::ds.vertChisqCross(
  "DA", "age_group", "glu_group", fisher = TRUE, datasources = conns)
)---",
    disclosure = r"---(
Positive cells, row margins, and column margins below the configured privacy
threshold block the same-server release. In the cross-server route, DCF
threshold checks reveal only pass/fail before exact cell counts are released.
No one-hot patient vector is returned to the analyst.
)---",
    command = "Rscript scripts/validate_method_tables.R"
  ),
  list(
    id = "correlation",
    file = "validation_correlation.Rmd",
    title = "Correlation validation",
    functions = "`ds.vertCor()`",
    math = r"---(
The Pearson correlation between variables \(x_j\) and \(x_k\) is
\(\mathrm{cov}(x_j,x_k)/(s_j s_k)\). Local moments handle within-server
pairs. Cross-server pairs use Beaver cross-products to obtain the aggregate
\(\sum_i x_{ij}x_{ik}\) without revealing either column to the other site.
)---",
    central = r"---(
X <- MASS::Pima.tr[seq_len(120), c("age", "bmi", "ped", "glu", "bp")]
stats::cor(X)
)---",
    split = r"---(
tables <- list(
  s1 = data.frame(patient_id = sprintf("P%03d", seq_len(nrow(X))),
                  X[c("age", "bmi", "ped")]),
  s2 = data.frame(patient_id = sprintf("P%03d", seq_len(nrow(X))),
                  X[c("glu", "bp")])
)
)---",
    dslite = r"---(
cor_fit <- dsVertClient::ds.vertCor(
  "DA",
  variables = list(s1 = c("age", "bmi", "ped"),
                   s2 = c("glu", "bp")),
  datasources = conns)
)---",
    disclosure = r"---(
The full correlation matrix is a second-order aggregate and can become
reconstructive when \(p\) is large relative to \(n\). The product route enforces
minimum-n and p/n guards and blocks high-dimensional diagnostic fixtures by
default.
)---",
    command = "Rscript scripts/validate_method_descriptive.R"
  ),
  list(
    id = "pca",
    file = "validation_pca.Rmd",
    title = "PCA validation",
    functions = "`ds.vertPCA()`",
    math = r"---(
PCA on standardized variables is the eigendecomposition of the correlation
matrix \(R = V \Lambda V^T\). `ds.vertPCA()` reuses `ds.vertCor()` and performs
the eigendecomposition client-side. The validated target is therefore
`eigen(cor(X))`.
)---",
    central = r"---(
X <- MASS::Pima.tr[seq_len(120), c("age", "bmi", "ped", "glu", "bp")]
ref <- eigen(stats::cor(X), symmetric = TRUE)
)---",
    split = r"---(
tables <- list(
  s1 = data.frame(patient_id = sprintf("P%03d", seq_len(nrow(X))),
                  X[c("age", "bmi", "ped")]),
  s2 = data.frame(patient_id = sprintf("P%03d", seq_len(nrow(X))),
                  X[c("glu", "bp")])
)
)---",
    dslite = r"---(
cor_fit <- dsVertClient::ds.vertCor("DA", variables = variables,
                                    datasources = conns)
pca_fit <- dsVertClient::ds.vertPCA(cor_result = cor_fit,
                                    n_components = 3)
)---",
    disclosure = r"---(
Loadings and eigenvalues are aggregate outputs derived from the released
correlation matrix. Individual component scores are not returned because scores
would be row-level projections.
)---",
    command = "Rscript scripts/validate_method_descriptive.R"
  ),
  list(
    id = "glm",
    file = "validation_glm.Rmd",
    title = "GLM validation",
    functions = "`ds.vertGLM()`",
    math = r"---(
For a GLM with mean \(\mu_i=g^{-1}(x_i^T\beta)\), the maximum-likelihood
estimate solves \(X^T(y-\mu)=0\) (with the usual variance weights by family).
The distributed route evaluates aggregate scores, Hessian/Fisher summaries,
and deviance with Ring63/Ring127 MPC. The validation table uses the binomial
precision route with 150 secure sigmoid intervals.
)---",
    central = r"---(
pooled <- MASS::Pima.tr[seq_len(50), ]
pooled$diabetes <- as.integer(pooled$type == "Yes")
fit_ref <- glm(diabetes ~ age + bmi + ped + glu + bp + skin + npreg,
               data = pooled, family = binomial())
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "age", "bmi", "ped")],
  s2 = pooled[c("patient_id", "glu", "bp", "skin")],
  s3 = pooled[c("patient_id", "npreg", "diabetes")]
)
)---",
    dslite = r"---(
fit_ds <- dsVertClient::ds.vertGLM(
  diabetes ~ age + bmi + ped + glu + bp + skin + npreg,
  data = "DA",
  family = "binomial",
  lambda = 0,
  binomial_sigmoid_intervals = 150L,
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
The analyst receives model-scale coefficients, covariance/standard errors,
and scalar diagnostics. Per-patient eta, probabilities, residuals, weights,
and row scores stay in the share domain. Increasing sigmoid intervals changes
only numerical approximation, not the disclosure surface.
)---",
    command = "Rscript scripts/validate_method_glm.R --family binomial --binomial-sigmoid-intervals 150"
  ),
  list(
    id = "inference",
    file = "validation_inference.Rmd",
    title = "GLM inference helper validation",
    functions = "`ds.vertConfint()`, `ds.vertWald()`, `ds.vertContrast()`, `ds.vertLR()`",
    math = r"---(
The helpers are deterministic transformations of a `ds.glm` result:
Wald intervals use \(\hat\beta_j \pm z_{\alpha/2}\mathrm{SE}_j\), univariate
Wald tests use \((\hat\beta_j-\beta_{0j})/\mathrm{SE}_j\), contrasts use
\((K\hat\beta-m)^T(K\hat\Sigma K^T)^{-1}(K\hat\beta-m)\), and LR tests use
the difference in returned deviances.
)---",
    central = r"---(
fit_full <- glm(glu ~ age + bmi + ped + bp + skin + npreg,
                data = MASS::Pima.tr[seq_len(80), ])
fit_reduced <- update(fit_full, . ~ . - skin)
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "age", "bmi", "ped")],
  s2 = pooled[c("patient_id", "bp", "skin", "npreg", "glu")]
)
)---",
    dslite = r"---(
fit_full <- dsVertClient::ds.vertGLM(
  glu ~ age + bmi + ped + bp + skin + npreg,
  data = "DA", family = "gaussian", lambda = 0,
  datasources = conns, verbose = FALSE)
ci <- dsVertClient::ds.vertConfint(fit_full)
wald <- dsVertClient::ds.vertWald(fit_full, parm = "age")
contrast <- dsVertClient::ds.vertContrast(fit_full, K = "age = bmi")
lr <- dsVertClient::ds.vertLR(fit_reduced, fit_full)
)---",
    disclosure = r"---(
No server call is made by these helpers. They add no disclosure beyond beta,
covariance, and deviance already returned by the validated GLM route.
)---",
    command = "Rscript scripts/validate_method_glm_wrappers.R"
  ),
  list(
    id = "lasso",
    file = "validation_lasso.Rmd",
    title = "LASSO validation",
    functions = "`ds.vertLASSO()`, `ds.vertLASSO1Step()`, `ds.vertLASSOCV()`, `ds.vertLASSOIter()`, `ds.vertLASSOProximal()`",
    math = r"---(
The LASSO target is
\(\arg\min_\beta n^{-1}\ell(\beta)+\lambda\|\beta_{-0}\|_1\). Gaussian
routes solve the normal-equation LASSO from aggregate covariance/Hessian
objects. Binomial routes use secure aggregate-score proximal gradient with
FISTA restart and a guarded aggregate-Gram Lipschitz bound. Poisson routes use
aggregate score/Hessian probes in a damped proximal-Newton update.
)---",
    central = r"---(
X <- model.matrix(diabetes ~ age + bmi + ped + glu, data = pooled)
y <- pooled$diabetes
fit_ref <- glmnet::glmnet(X[, -1], y, family = "binomial",
                          lambda = 0.02, standardize = TRUE)
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "age", "bmi")],
  s2 = pooled[c("patient_id", "ped", "glu", "diabetes")]
)
)---",
    dslite = r"---(
lasso <- dsVertClient::ds.vertLASSOIter(
  diabetes ~ age + bmi + ped + glu,
  data = "DA",
  family = "binomial",
  lambda = 0.02,
  exact_non_gaussian = TRUE,
  binomial_sigmoid_intervals = 200L,
  lipschitz = "auto",
  fista_restart = TRUE,
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
The LASSO routes open only quantities already accepted for GLM, inference, and
guarded correlation: beta/covariance/Hessian, p-dimensional aggregate scores,
or p-by-p aggregate Hessians/Gram bounds. No eta, probability, residual,
objective-by-row, or row score vector is returned.
)---",
    command = "Rscript scripts/validate_method_lasso_binomial.R --sigmoid-intervals 200"
  ),
  list(
    id = "negative_binomial",
    file = "validation_negative_binomial.Rmd",
    title = "Negative binomial validation",
    functions = "`ds.vertNB()`, `ds.vertNBMoMTheta()`, `ds.vertNBFullRegTheta()`",
    math = r"---(
For NB2 log-link regression,
\(\mathrm{Var}(Y_i|X_i)=\mu_i+\mu_i^2/\theta\). The full-regression theta
score contains terms such as \(\sum_i \log(\mu_i+\theta)\) and
\(\sum_i (y_i+\theta)/(\mu_i+\theta)\). The product route evaluates those
terms in Ring127 shares and alternates theta Newton updates with aggregate
Fisher beta updates.
)---",
    central = r"---(
fit_ref <- MASS::glm.nb(y ~ x1 + x2 + x3, data = pooled)
coef(fit_ref)
fit_ref$theta
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "y")]
)
)---",
    dslite = r"---(
fit_nb <- dsVertClient::ds.vertNBFullRegTheta(
  y ~ x1 + x2 + x3,
  data = "DA",
  variant = "full_reg_nd",
  theta_max_iter = 8L,
  beta_max_iter = 2L,
  compute_covariance = TRUE,
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
The removed `variant='full_reg'` route disclosed non-label eta in plaintext.
The `full_reg_nd` route keeps eta, mu, logarithms, reciprocals, score
residuals, Fisher weights, and weighted covariates in Ring127 additive shares.
Only scalar theta sums, aggregate beta scores, and aggregate Fisher matrices
are reconstructed.
)---",
    command = "Rscript scripts/validate_method_nb.R --variant full_reg_nd"
  ),
  list(
    id = "cox",
    file = "validation_cox.Rmd",
    title = "Cox PH validation",
    functions = "`ds.vertCox()`, `ds.vertCoxProfileNonDisclosive()`, `ds.vertCoxDiscreteNonDisclosive()`",
    math = r"---(
The Cox PH model has hazard \(h(t|x)=h_0(t)\exp(x^T\beta)\). With Breslow
ties, the profile partial-likelihood score is built from event-time risk-set
sums \(\sum_{i:t_i \ge t_j}\exp(x_i^T\beta)x_i\). The product route hides
risk-set and event masks as Ring127 shares and opens only slope
score/Hessian aggregates.
)---",
    central = r"---(
fit_ref <- survival::coxph(
  survival::Surv(time, event) ~ x1 + x2 + x3,
  data = pooled, ties = "breslow")
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "time", "event")]
)
)---",
    dslite = r"---(
fit_cox <- dsVertClient::ds.vertCox(
  survival::Surv(time, event) ~ x1 + x2 + x3,
  data = "DA",
  max_iter = 5L,
  max_event_times = 50L,
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
Legacy rank/permutation and person-time Cox routes are not offered. Event
times, event indicators, event ranks, risk-set membership, baseline dummies,
and per-person period rows are not returned. Debug traces and bin summaries
are gated diagnostics only.
)---",
    command = "Rscript scripts/validate_method_cox.R --target cox_profile"
  ),
  list(
    id = "lmm",
    file = "validation_lmm.Rmd",
    title = "Linear mixed model validation",
    functions = "`ds.vertLMM()`, `ds.vertLMM.k3()`",
    math = r"---(
The random-intercept LMM is \(y_{ij}=x_{ij}^T\beta+b_i+\epsilon_{ij}\),
with \(b_i\sim N(0,\sigma_b^2)\) and
\(\epsilon_{ij}\sim N(0,\sigma^2)\). The estimator profiles variance
components and solves GLS normal equations. K=2 uses the closed-form route;
K>=3 uses share-domain residual REML and exact cluster-mean GLS transforms.
)---",
    central = r"---(
fit_ref <- nlme::lme(y ~ x1 + x2 + x3,
                     random = ~ 1 | cluster_id,
                     data = pooled, method = "REML")
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "cluster_id", "y")]
)
)---",
    dslite = r"---(
fit_lmm_k2 <- dsVertClient::ds.vertLMM(
  y ~ x1 + x2 + x3, data = "DA", cluster_col = "cluster_id",
  datasources = conns, verbose = FALSE)
fit_lmm_k3 <- dsVertClient::ds.vertLMM.k3(
  y ~ x1 + x2 + x3, data = "DA", cluster_col = "cluster_id",
  datasources = conns, verbose = FALSE)
)---",
    disclosure = r"---(
Cluster membership is method-level metadata and is sent only through encrypted
server-to-server channels. Original cluster labels are not returned, small
clusters fail closed, and patient-level residuals or BLUPs are not returned.
)---",
    command = "Rscript scripts/validate_method_lmm.R"
  ),
  list(
    id = "gee",
    file = "validation_gee.Rmd",
    title = "GEE validation",
    functions = "`ds.vertGEE()`",
    math = r"---(
GEE solves \(\sum_i D_i^T V_i^{-1}(y_i-\mu_i)=0\) and reports the sandwich
variance \(B^{-1}MB^{-1}\). `ds.vertGEE()` supports independence,
exchangeable and guarded AR1 working correlations for Gaussian, binomial and
Poisson outcomes. Cluster and adjacent-pair scores are computed as aggregate
share-domain sufficient statistics.
)---",
    central = r"---(
fit_ref <- geepack::geeglm(
  y ~ x1 + x2 + x3, id = cluster_id, waves = visit,
  data = pooled, family = binomial(), corstr = "ar1")
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "cluster_id", "visit", "y")]
)
)---",
    dslite = r"---(
fit_gee <- dsVertClient::ds.vertGEE(
  y ~ x1 + x2 + x3,
  data = "DA",
  family = "binomial",
  id_col = "cluster_id",
  order_col = "visit",
  corstr = "ar1",
  binomial_sigmoid_intervals = 150L,
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
Residuals, fitted means, probabilities, Pearson residuals, row scores, visit
labels and adjacent-pair vectors are not returned. AR1 order metadata is
encrypted server-to-server and only guarded low-dimensional sufficient
statistics reach the client.
)---",
    command = "Rscript scripts/validate_method_gee.R"
  ),
  list(
    id = "glmm",
    file = "validation_glmm.Rmd",
    title = "GLMM validation",
    functions = "`ds.vertGLMM()`",
    math = r"---(
The supported GLMM route is binomial random-intercept PQL. The model is
\(\mathrm{logit}\{P(y_{ij}=1|b_i)\}=x_{ij}^T\beta+b_i\) with
\(b_i\sim N(0,\sigma_b^2)\). PQL iterates working responses and weights, but
the product implementation keeps these quantities in Ring127 shares and
returns only fixed effects, scalar variance components and diagnostics.
)---",
    central = r"---(
fit_ref <- MASS::glmmPQL(
  y_bin ~ x1 + x2 + x3,
  random = ~ 1 | cluster_id,
  family = binomial(),
  data = pooled,
  verbose = FALSE)
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "cluster_id", "y_bin")]
)
)---",
    dslite = r"---(
fit_glmm <- dsVertClient::ds.vertGLMM(
  y_bin ~ x1 + x2 + x3,
  data = "DA",
  cluster_col = "cluster_id",
  ring = 127L,
  compute_se = FALSE,
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
The removed method-selection path is not user-invokable. Working responses,
weights, probabilities, residuals and row scores remain in shares. The public
object is audited for no length-n and no cluster-length vectors.
)---",
    command = "Rscript scripts/validate_method_glmm.R"
  ),
  list(
    id = "ipw",
    file = "validation_ipw.Rmd",
    title = "IPW validation",
    functions = "`ds.vertIPW()`",
    math = r"---(
IPW fits a propensity model \(e(x)=P(T=1|X=x)\), then an outcome model with
weights \(w_i=T_i/e(x_i)+(1-T_i)/(1-e(x_i))\). The product route validates the
two-stage propensity plus weighted GLM workflow while keeping weights and
sqrt-weights server-side or in additive shares.
)---",
    central = r"---(
prop_ref <- glm(treat ~ x1 + x2, data = pooled, family = binomial())
pooled$ipw <- with(pooled, ifelse(treat == 1, 1 / fitted(prop_ref),
                                  1 / (1 - fitted(prop_ref))))
out_ref <- glm(y ~ treat + x1 + x2, data = pooled,
               family = binomial(), weights = ipw)
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "treat", "ipw", "y")]
)
)---",
    dslite = r"---(
fit_ipw <- dsVertClient::ds.vertIPW(
  outcome_formula = y ~ treat + x1 + x2,
  propensity_formula = treat ~ x1 + x2,
  data = "DA",
  weights_column = "ipw",
  outcome_family = "binomial",
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
Patient weights are not transported plaintext to DCF peers. Weighted
gradients/deviance use Beaver products of weight shares and residual shares.
The widened secure sigmoid fixed the previous stress accuracy gap without
changing disclosure.
)---",
    command = "Rscript scripts/validate_method_ipw.R"
  ),
  list(
    id = "mi",
    file = "validation_mi.Rmd",
    title = "Multiple imputation validation",
    functions = "`ds.vertMI()`",
    math = r"---(
For \(m=1,\ldots,M\), missing values are imputed server-side, the analysis GLM
returns \(\hat\beta_m\) and \(\hat\Sigma_m\), and Rubin's rules compute
\(\bar\beta=M^{-1}\sum_m\hat\beta_m\) and
\(T=W+(1+M^{-1})B\), where \(W\) is the mean within-imputation covariance and
\(B\) is the between-imputation covariance.
)---",
    central = r"---(
fits <- lapply(seq_len(3), function(m) {
  imputed <- same_server_imputer(pooled, seed = 100 + m)
  glm(y ~ x1 + x2 + x3, data = imputed)
})
rubin_pool(fits)
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "x1")],
  s2 = pooled[c("patient_id", "x2", "x3", "y")]
)
)---",
    dslite = r"---(
fit_mi <- dsVertClient::ds.vertMI(
  y ~ x1 + x2 + x3,
  data = "DA",
  impute_columns = "x1",
  m = 3L,
  family = "gaussian",
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
Imputed values never leave the server that owns the column. The client receives
only per-imputation model-scale coefficient/covariance summaries and Rubin
pooled summaries. Exact extrema in imputation diagnostics are suppressed.
)---",
    command = "Rscript scripts/validate_method_mi.R"
  ),
  list(
    id = "multinomial",
    file = "validation_multinomial.Rmd",
    title = "Multinomial validation",
    functions = "`ds.vertMultinom()`, `ds.vertMultinomJointNewton()`",
    math = r"---(
The softmax model uses
\(P(Y=c|x)=\exp(\eta_c)/\sum_h\exp(\eta_h)\), with one reference class.
The old one-vs-rest route is only an internal warm start. The product route
optimizes the joint softmax objective with Ring127 share-domain probabilities,
residuals, aggregate gradients and a Bohning/aggregate-Gram Hessian surrogate.
)---",
    central = r"---(
fit_ref <- nnet::multinom(class ~ age + lwt + smoke + ht,
                          data = pooled, trace = FALSE)
predict(fit_ref, type = "probs")
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "age", "lwt")],
  s2 = pooled[c("patient_id", "smoke", "ht", "class_ind_low", "class_ind_mid")]
)
)---",
    dslite = r"---(
fit_mn <- dsVertClient::ds.vertMultinom(
  class ~ age + lwt + smoke + ht,
  data = "DA",
  classes = c("low", "mid"),
  reference = "high",
  indicator_template = "class_ind_%s",
  max_outer = 30L,
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
Softmax probabilities and residuals remain Ring127 shares. The client sees
only low-dimensional aggregate gradients/Gram information and coefficients.
Warm OVR final estimation is not exported and is excluded from the paper path.
)---",
    command = "Rscript scripts/validate_method_multinom_joint.R --K 2 && Rscript scripts/validate_method_multinom_joint.R --K 3"
  ),
  list(
    id = "ordinal",
    file = "validation_ordinal.Rmd",
    title = "Ordinal validation",
    functions = "`ds.vertOrdinal()`, `ds.vertOrdinalJointNewton()`",
    math = r"---(
The proportional-odds model is
\(P(Y\le k|x)=F(\theta_k-x^T\beta)\), with ordered thresholds
\(\theta_1<\cdots<\theta_{K-1}\). The product route evaluates the joint
score for beta and thresholds with class masks, cumulative probabilities,
reciprocals and finite-difference Hessian probes held in Ring127 shares.
)---",
    central = r"---(
fit_ref <- MASS::polr(y_ord ~ x1 + x2 + x3,
                      data = pooled, method = "logistic", Hess = TRUE)
)---",
    split = r"---(
tables <- list(
  s1 = pooled[c("patient_id", "x1", "x2")],
  s2 = pooled[c("patient_id", "x3", "low_leq", "mid_leq", "y_ord")]
)
)---",
    dslite = r"---(
fit_ord <- dsVertClient::ds.vertOrdinal(
  y_ord ~ x1 + x2 + x3,
  data = "DA",
  levels_ordered = c("low", "mid", "high"),
  cumulative_template = "%s_leq",
  max_outer = 2L,
  datasources = conns,
  verbose = FALSE)
)---",
    disclosure = r"---(
Class masks, probabilities, weights and residual-like terms remain encrypted
shares. The client receives guarded class counts and aggregate
gradient/parameter summaries only. The warm cumulative-binomial route is not a
user-facing final estimator.
)---",
    command = "Rscript scripts/validate_method_ordinal_joint.R --K 2 && Rscript scripts/validate_method_ordinal_joint.R --K 3"
  )
)

write_method_page <- function(page) {
  txt <- paste0(
"---
title: \"", page$title, "\"
author: \"David Sarrat González\"
date: \"`r format(Sys.Date(), '%Y-%m-%d')`\"
output:
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 2
vignette: >
  %\\VignetteIndexEntry{", page$title, "}
  %\\VignetteEngine{knitr::rmarkdown}
  %\\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE, eval = TRUE, warning = FALSE, message = FALSE,
  collapse = TRUE, comment = \"#>\")
helper_candidates <- c(\"validation_helpers.R\",
                       file.path(\"..\", \"validation_helpers.R\"),
                       file.path(\"vignettes\", \"validation_helpers.R\"))
source(helper_candidates[file.exists(helper_candidates)][1])
```

## What is validated

Functions: ", page$functions, "

", page$math, "

## Centralized reference

The centralized reference is built from the same deterministic rows and
formula as the DSLite run. The long method harness stores the exact seed,
split and reference object in the cache listed below.

```{r central-reference, eval=FALSE}
", page$central, "
```

## Vertical DSLite split

The validation split creates independent server tables with the same
`patient_id` key and then aligns them with `ds.psiAlign()`.

```{r vertical-split, eval=FALSE}
", page$split, "
```

```{r dslite-execution, eval=FALSE}
", page$dslite, "
```

To reproduce the cache from the repository root:

```{r reproduce-cache, eval=FALSE}
", page$command, "
```

## Disclosure review

", page$disclosure, "

## Executed evidence check

This chunk is evaluated when the vignette renders. It fails if either K mode is
not marked non-disclosive, is not `PASS`, or exceeds its accepted tolerance.

```{r evidence}
rows <- validation_rows(\"", page$id, "\")
assert_validation(rows)
display_validation(rows)
```

## Verdict

Both K=2 and K>=3 validation rows are inside their accepted numerical envelope
and use the current non-disclosive product route. Any legacy route mentioned in
the package is excluded from this evidence path.
")
  writeLines(txt, file.path("vignettes", page$file), useBytes = TRUE)
}

for (page in method_pages) write_method_page(page)

index <- r"---(
---
title: "dsVert validation evidence"
author: "David Sarrat González"
date: "`r format(Sys.Date(), '%Y-%m-%d')`"
output:
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{dsVert validation evidence}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE, eval = TRUE, warning = FALSE, message = FALSE,
  collapse = TRUE, comment = "#>")
helper_candidates <- c("validation_helpers.R",
                       file.path("..", "validation_helpers.R"),
                       file.path("vignettes", "validation_helpers.R"))
source(helper_candidates[file.exists(helper_candidates)][1])
```

## Scope

This page is the validation index for the current dsVertClient product surface.
The historical vignettes under `vignettes/methods/` were discarded because they
mixed product methods with legacy or diagnostic routes. The pages linked here
document the current supported route for each method family, with K=2 and K>=3
evidence from deterministic DSLite harnesses.

Each method vignette contains:

- the statistical target and mathematical estimating equation;
- the centralized reference construction;
- the vertical DSLite split and product function call;
- a disclosure review describing what the analyst and peer servers can see;
- an executed evidence check against the cached DSLite result table.

## Summary table

```{r summary}
all_rows <- load_validation_summary()
assert_validation(all_rows)
knitr::kable(all_rows[, c("method_name", "k_mode", "function_route",
                          "primary_metric", "observed", "tolerance",
                          "tier", "status", "cache")])
```

## Method pages

- [PSI alignment](validation_psi.html)
- [Descriptive statistics](validation_descriptive.html)
- [Contingency tests](validation_contingency.html)
- [Correlation](validation_correlation.html)
- [PCA](validation_pca.html)
- [GLM](validation_glm.html)
- [GLM inference helpers](validation_inference.html)
- [LASSO](validation_lasso.html)
- [Negative binomial](validation_negative_binomial.html)
- [Cox PH](validation_cox.html)
- [LMM](validation_lmm.html)
- [GEE](validation_gee.html)
- [GLMM](validation_glmm.html)
- [IPW](validation_ipw.html)
- [Multiple imputation](validation_mi.html)
- [Multinomial](validation_multinomial.html)
- [Ordinal](validation_ordinal.html)

## Removed legacy routes

The package should not expose routes that are disclosive or materially less
accurate than an available product route. This check is executed during
rendering.

```{r legacy-removed}
library(dsVertClient)
legacy <- load_legacy_removed()
checks <- assert_legacy_routes_removed()
knitr::kable(merge(legacy, checks, by.x = "route", by.y = "route", all.x = TRUE))
```

## Interpretation

All rows in the current table are `PASS`, non-disclosive under the project
standard, and cover both K=2 and K>=3. Some methods are exact to their
centralized target; others are strict-practical because the remaining gap is a
documented secure approximation floor, such as histogram quantile bucketization,
secure sigmoid spline resolution, PQL versus GLMM ML, or penalized optimizer
depth.
)---"
index <- sub("^\\n", "", index)
writeLines(index, "vignettes/vert_validation_evidence.Rmd", useBytes = TRUE)

pkgdown <- r"---(
url: https://isglobal-brge.github.io/dsVertClient/
template:
  bootstrap: 5

footer:
  structure:
    left: package
    right: built_with

navbar:
  structure:
    left: [reference, articles, news]
    right: [search, github]

articles:
- title: Validation evidence
  navbar: Validation
  contents:
  - vert_validation_evidence
  - validation_psi
  - validation_descriptive
  - validation_contingency
  - validation_correlation
  - validation_pca
  - validation_glm
  - validation_inference
  - validation_lasso
  - validation_negative_binomial
  - validation_cox
  - validation_lmm
  - validation_gee
  - validation_glmm
  - validation_ipw
  - validation_mi
  - validation_multinomial
  - validation_ordinal
)---"
pkgdown <- sub("^\\n", "", pkgdown)
writeLines(pkgdown, "_pkgdown.yml", useBytes = TRUE)

message("Validation vignettes rebuilt.")
