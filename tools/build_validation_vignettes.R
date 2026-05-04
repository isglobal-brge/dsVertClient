#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

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
unlink("vignettes/methods", recursive = TRUE, force = TRUE)
unlink("vignettes/validation-data", recursive = TRUE, force = TRUE)

method_pages <- list(
  psi = list(
    title = "PSI alignment validation",
    functions = "ds.psiAlign(), ds.isPsiAligned()",
    method = "Private set intersection aligns vertically split tables on a shared identifier. Each server blinds its identifiers, the protocol computes the intersection, and the aligned data object keeps only matched rows in a common order.",
    math = "The statistical target is the set I = intersection_k I_k. The validation checks that every server has |I| rows after alignment and that a correlation of deterministic id-derived columns is 1 after alignment.",
    fixture = "Synthetic overlapping patient_id ranges with shuffled row order.",
    reference = "Deterministic set intersection in centralized R.",
    disclosure = "The analyst receives counts and status fields. Patient identifiers, matched row indices, and length-n index vectors are not returned; the legacy index reveal aggregate is checked as blocked."
  ),
  descriptive = list(
    title = "Descriptive statistics validation",
    functions = "ds.vertDesc()",
    method = "Each server computes local guarded moments for its own variables. Quantiles are estimated from guarded histograms rather than releasing sorted observations.",
    math = "For a variable x, the exact checked moments are n, mean(x), and sd(x). Histogram quantiles approximate F^{-1}(p) from binned counts under the disclosure threshold.",
    fixture = "MASS::Pima.tr numeric fixture split vertically.",
    reference = "Central mean and standard deviation on the pooled fixture.",
    disclosure = "Only scalar summaries and guarded histogram information are returned; no row-level values or order statistics are disclosed."
  ),
  contingency = list(
    title = "Contingency tests validation",
    functions = "ds.vertChisq(), ds.vertFisher(), ds.vertChisqCross()",
    method = "Same-server categorical tests use guarded local contingency counts. Cross-server tests build one-hot shares and aggregate cell counts with MPC before applying chi-square/Fisher tests.",
    math = "The validation compares O_ab counts and X^2 = sum_ab (O_ab - E_ab)^2 / E_ab with centralized R.",
    fixture = "Categorised MASS::Pima.tr variables with all margins above disclosure thresholds.",
    reference = "chisq.test() and fisher.test() on the pooled fixture.",
    disclosure = "The released object contains guarded table counts and scalar test statistics. Small cells fail closed."
  ),
  correlation = list(
    title = "Correlation validation",
    functions = "ds.vertCor()",
    method = "The method releases a low-dimensional correlation matrix from server-local moments and MPC cross-products.",
    math = "For variables x_j and x_l, r_jl = cov(x_j, x_l) / (sd(x_j) sd(x_l)). The validation compares the full matrix to stats::cor().",
    fixture = "MASS::Pima.tr numeric fixture split vertically.",
    reference = "stats::cor() on the pooled fixture.",
    disclosure = "The disclosure surface is the guarded p by p correlation matrix, the same aggregate tier used by downstream PCA."
  ),
  pca = list(
    title = "PCA validation",
    functions = "ds.vertPCA()",
    method = "PCA is computed from the validated correlation matrix, not from patient-level scores.",
    math = "The target is the eigendecomposition R = V Lambda V^T. The validation compares eigenvalues and sign-aligned loadings.",
    fixture = "MASS::Pima.tr numeric fixture split vertically.",
    reference = "eigen(cor(X)) on the pooled fixture.",
    disclosure = "Only eigenvalues and loadings derived from the aggregate correlation matrix are returned; no individual PCA scores are released."
  ),
  glm = list(
    title = "GLM validation",
    functions = "ds.vertGLM()",
    method = "The distributed route evaluates score and Fisher/Hessian aggregates with Ring MPC. This compact vignette uses a Gaussian GLM so the full validation runs quickly while exercising the same vertical design assembly.",
    math = "For mean mu_i = g^{-1}(x_i^T beta), beta solves X^T(y - mu) = 0 with family-specific weights. The Gaussian fixture is compared to lm().",
    fixture = "MASS::Pima.tr fixture with diabetes coded 0/1 and predictors split vertically.",
    reference = "stats::lm() on the pooled fixture.",
    disclosure = "The analyst receives model-level coefficients, covariance/standard errors, and scalar diagnostics. Eta, residuals, row scores, and weights remain internal."
  ),
  inference = list(
    title = "GLM inference validation",
    functions = "ds.vertConfint(), ds.vertWald(), ds.vertContrast(), ds.vertLR()",
    method = "Inference helpers post-process the beta, covariance, deviance, and degrees of freedom already returned by ds.vertGLM.",
    math = "Wald tests use z = (K beta - m) / sqrt(K V K^T). LR tests use D_reduced - D_full against a chi-square reference.",
    fixture = "The same MASS::Pima.tr Gaussian GLM fixture used for GLM validation.",
    reference = "Manual algebra on ds.vertGLM output.",
    disclosure = "No new server aggregate is required beyond the already accepted model-level GLM outputs."
  ),
  lasso = list(
    title = "LASSO validation",
    functions = "ds.vertLASSOProximal()",
    method = "The proximal LASSO solver consumes the distributed GLM fit and optimizes the penalized quadratic objective client-side.",
    math = "The checked identity is the lambda = 0 limit: argmin 1/2 (beta - beta_hat)^T H (beta - beta_hat) equals beta_hat.",
    fixture = "MASS::Pima.tr Gaussian GLM fixture.",
    reference = "The prime ds.vertGLM coefficients at lambda = 0.",
    disclosure = "The proximal step uses only model-level beta/Hessian information and returns coefficients; it does not request row-level data."
  ),
  negative_binomial = list(
    title = "Negative binomial validation",
    functions = "ds.vertNBMoMTheta()",
    method = "The product route fits a negative-binomial mean model and estimates overdispersion through guarded aggregate moments.",
    math = "The model has Var(Y|X) = mu + mu^2/theta and log(mu)=X beta. The validation compares beta to MASS::glm.nb().",
    fixture = "Synthetic overdispersed count regression fixture with vertical predictors.",
    reference = "MASS::glm.nb() on the pooled fixture.",
    disclosure = "Returned values are beta and scalar theta/diagnostics. Per-patient means, residuals, and score terms are not returned."
  ),
  cox = list(
    title = "Cox PH validation",
    functions = "ds.vertCox()",
    method = "The current product route evaluates Cox partial-likelihood score terms over discretised event-time risk sets without releasing risk-set rows.",
    math = "The target solves sum_i delta_i (x_i - weighted_risk_mean(t_i; beta)) = 0. The fixture compares beta to coxph(ties='breslow').",
    fixture = "Synthetic survival fixture with event times discretised into guarded bins.",
    reference = "survival::coxph(..., ties = 'breslow') on the pooled fixture.",
    disclosure = "Risk-set contributions stay in the share domain; the client receives coefficients and scalar convergence diagnostics."
  ),
  lmm = list(
    title = "LMM validation",
    functions = "ds.vertLMM(), ds.vertLMM.k3()",
    method = "Random-intercept LMM uses guarded cluster metadata and MPC cross-products to fit fixed effects and variance components.",
    math = "The model is y = X beta + b_cluster + epsilon with b ~ N(0, sigma_b^2). The validation compares fixed effects to lme4::lmer().",
    fixture = "Synthetic balanced random-intercept fixture, with cluster sizes above privacy thresholds.",
    reference = "lme4::lmer(... + (1 | cluster)) on the pooled fixture.",
    disclosure = "Cluster membership is used internally at the accepted LMM tier. Per-cluster residuals and BLUP vectors are not returned."
  ),
  gee = list(
    title = "GEE validation",
    functions = "ds.vertGEE()",
    method = "The GEE route computes estimating-equation aggregates and robust sandwich summaries for clustered data.",
    math = "For independence working correlation, the estimating equation reduces to sum_i D_i^T V_i^{-1}(y_i - mu_i)=0. The validation compares beta to geepack.",
    fixture = "MASS::Pima.tr clustered fixture with vertical predictors.",
    reference = "geepack::geeglm(corstr='independence') on the pooled fixture.",
    disclosure = "Returned values are regression and sandwich-level aggregates; row scores are not returned."
  ),
  glmm = list(
    title = "GLMM validation",
    functions = "ds.vertGLMM()",
    method = "The product route is aggregate PQL for a binomial random-intercept GLMM. The compact executable vignette runs one protected PQL outer update and checks that the fixed-effect symmetry target and variance trace are finite.",
    math = "A PQL step alternates protected binomial working-response updates with aggregate weighted mixed-model normal equations. The validation compares fixed effects to the central symmetric glm target and asserts that the PQL trace is populated.",
    fixture = "Paired balanced binomial cluster fixture with cluster sizes above privacy thresholds.",
    reference = "Central binomial glm fixed-effect symmetry plus a one-step PQL trace check.",
    disclosure = "The route returns fixed effects and scalar variance diagnostics only; per-cluster BLUPs and row probabilities are not returned."
  ),
  ipw = list(
    title = "IPW validation",
    functions = "ds.vertIPW()",
    method = "IPW runs a protected propensity GLM and then a protected weighted outcome GLM using an already server-side weight column.",
    math = "The checked outcome target solves X^T W(y - X beta)=0 with inverse-probability weights W.",
    fixture = "Synthetic confounded treatment/outcome fixture with known IPW column.",
    reference = "Central weighted lm using the same weights.",
    disclosure = "Weights remain a server column; the analyst receives only propensity and outcome model-level fits."
  ),
  mi = list(
    title = "Multiple imputation validation",
    functions = "ds.vertMI()",
    method = "Missing covariates are imputed server-side for each imputation round, then GLM fits are pooled client-side with Rubin's rules.",
    math = "Rubin pooling uses beta_bar = mean_m beta_m and T = W + (1 + 1/M)B for total variance.",
    fixture = "Synthetic Gaussian regression fixture with deterministic missingness in one covariate.",
    reference = "Central mean-imputation regression used as a compact reproducibility check.",
    disclosure = "Imputed columns stay server-side. The client sees only per-imputation beta/covariance and pooled model summaries."
  ),
  multinomial = list(
    title = "Multinomial validation",
    functions = "ds.vertMultinom(), ds.vertMultinomJointNewton()",
    method = "The supported route is joint softmax Newton over Ring127 shares. Historical one-vs-rest output is kept only as an internal warm start.",
    math = "P(Y=c|X)=exp(eta_c)/sum_l exp(eta_l). The validation compares predicted class probabilities to nnet::multinom().",
    fixture = "Balanced soft-signal synthetic three-class fixture with vertical predictors.",
    reference = "nnet::multinom predicted probabilities.",
    disclosure = "Softmax probabilities, residuals, and row scores remain Ring127 shares. The client receives coefficients and scalar optimizer diagnostics."
  ),
  ordinal = list(
    title = "Ordinal validation",
    functions = "ds.vertOrdinal(), ds.vertOrdinalJointNewton()",
    method = "The supported route is joint proportional-odds Newton over Ring127 shares. Historical cumulative-binomial output is internal warm start only.",
    math = "P(Y <= k | X)=sigmoid(theta_k - X beta). The validation compares cumulative probabilities to MASS::polr().",
    fixture = "Balanced synthetic three-level ordered fixture with vertical predictors.",
    reference = "MASS::polr cumulative probabilities.",
    disclosure = "Cumulative probabilities, residuals, and row scores stay in shares; the client receives thresholds, slopes, and scalar optimizer diagnostics."
  )
)

write_page <- function(id, meta) {
  file <- file.path("vignettes", paste0("validation_", id, ".Rmd"))
  title <- meta$title
  lines <- c(
    "---",
    paste0("title: \"", title, "\""),
    "author: \"David Sarrat Gonzalez\"",
    "date: \"`r format(Sys.Date(), '%Y-%m-%d')`\"",
    "output:",
    "  rmarkdown::html_vignette:",
    "    toc: true",
    "    toc_depth: 2",
    "vignette: >",
    paste0("  %\\VignetteIndexEntry{", title, "}"),
    "  %\\VignetteEngine{knitr::rmarkdown}",
    "  %\\VignetteEncoding{UTF-8}",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = TRUE, eval = TRUE, warning = FALSE,",
    "  message = FALSE, collapse = TRUE, comment = \"#>\")",
    "helper_candidates <- c(\"validation_helpers.R\",",
    "  file.path(\"..\", \"validation_helpers.R\"),",
    "  file.path(\"vignettes\", \"validation_helpers.R\"))",
    "source(helper_candidates[file.exists(helper_candidates)][1])",
    "force_run <- identical(tolower(Sys.getenv(\"DSVERT_VALIDATION_FORCE\", \"false\")), \"true\")",
    "```",
    "",
    "## What is validated",
    "",
    paste0("Functions: `", meta$functions, "`."),
    "",
    meta$method,
    "",
    "## Mathematical target",
    "",
    meta$math,
    "",
    "## Fixture and reference",
    "",
    paste0("Fixture: ", meta$fixture),
    "",
    paste0("Centralized reference: ", meta$reference),
    "",
    "The executable chunk below calls `run_validation()` from",
    "`vignettes/validation_helpers.R`. That helper constructs the fixture,",
    "opens a DSLite server, performs PSI alignment, runs the dsVertClient",
    "product route for K=2 and K=3, computes the centralized reference, and",
    "compares both results. No RDS or result table outside this package is",
    "required; if a local `vignettes/validation-cache/` file exists it is a",
    "cache produced by this same execution path.",
    "",
    "## Disclosure review",
    "",
    meta$disclosure,
    "",
    "The fixture keeps `datashield.privacyLevel = 5` and is sized so the",
    "standard disclosure guards remain active. Only `dsvert.require_trusted_peers`",
    "is disabled for DSLite because there is no real Opal/Rock deployment in",
    "this local validation context.",
    "",
    "## Executed evidence",
    "",
    "```{r evidence}",
    paste0("rows <- run_validation(\"", id, "\", force = force_run)"),
    "display_validation(rows)",
    "```",
    "",
    "```{r assertion, include=FALSE}",
    "assert_validation(rows)",
    "```",
    "",
    "## Verdict",
    "",
    "The vignette fails during rendering if either K=2 or K>=3 leaves the",
    "accepted numerical envelope or is marked as disclosive.")
  writeLines(lines, file)
}

for (id in names(method_pages)) write_page(id, method_pages[[id]])

index_lines <- c(
  "---",
  "title: \"dsVert validation evidence\"",
  "author: \"David Sarrat Gonzalez\"",
  "date: \"`r format(Sys.Date(), '%Y-%m-%d')`\"",
  "output:",
  "  rmarkdown::html_vignette:",
  "    toc: true",
  "    toc_depth: 2",
  "vignette: >",
  "  %\\VignetteIndexEntry{dsVert validation evidence}",
  "  %\\VignetteEngine{knitr::rmarkdown}",
  "  %\\VignetteEncoding{UTF-8}",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = TRUE, eval = TRUE, warning = FALSE,",
  "  message = FALSE, collapse = TRUE, comment = \"#>\")",
  "helper_candidates <- c(\"validation_helpers.R\",",
  "  file.path(\"..\", \"validation_helpers.R\"),",
  "  file.path(\"vignettes\", \"validation_helpers.R\"))",
  "source(helper_candidates[file.exists(helper_candidates)][1])",
  "force_run <- identical(tolower(Sys.getenv(\"DSVERT_VALIDATION_FORCE\", \"false\")), \"true\")",
  "```",
  "",
  "## Scope",
  "",
  "This is the executable validation index for the current dsVertClient",
  "product surface. It covers every supported method in K=2 and K>=3 mode.",
  "Each method page documents the statistical target, fixture, centralized",
  "reference, DSLite execution path, and disclosure surface.",
  "",
  "No external RDS files or private repository caches are required. On a clean",
  "checkout, rendering this page builds all fixtures and runs the DSLite",
  "validations. Local cache files under `vignettes/validation-cache/` are",
  "generated by this same code path and can be deleted safely.",
  "",
  "## Summary table",
  "",
  "```{r summary}",
  "all_rows <- run_all_validations(force = force_run)",
  "assert_validation(all_rows)",
  "knitr::kable(all_rows[, c(\"method_name\", \"k_mode\", \"function_route\",",
  "  \"primary_metric\", \"observed\", \"tolerance\", \"tier\", \"status\",",
  "  \"runtime_s\")])",
  "```",
  "",
  "## Method pages",
  "",
  paste0("- [", vapply(method_pages, `[[`, character(1), "title"), "](",
         paste0("validation_", names(method_pages), ".html"), ")"),
  "",
  "## Removed legacy routes",
  "",
  "The package should not expose routes that were disclosive or materially",
  "less accurate than an available product route. This check is executed",
  "during rendering.",
  "",
  "```{r legacy-removed}",
  "legacy <- assert_legacy_routes_removed()",
  "knitr::kable(legacy)",
  "```")
writeLines(index_lines, "vignettes/vert_validation_evidence.Rmd")
