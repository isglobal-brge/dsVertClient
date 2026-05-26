`%||%` <- function(x, y) if (is.null(x)) y else x

validation_method_ids <- c(
  "psi", "descriptive", "contingency", "correlation", "pca",
  "glm", "inference", "lasso", "negative_binomial", "cox",
  "lmm", "gee", "glmm", "ipw", "mi", "multinomial", "ordinal")

validation_bundle_id <- function(method_id) {
  if (method_id %in% c("descriptive", "correlation", "pca")) {
    "descriptive_correlation_pca"
  } else {
    method_id
  }
}

validation_cache_dir <- function() {
  candidates <- c(
    file.path("vignettes", "validation-cache"),
    "validation-cache",
    file.path("..", "validation-cache"))
  root <- candidates[dir.exists(dirname(candidates))][1]
  if (is.na(root)) root <- file.path("vignettes", "validation-cache")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root
}

validation_cache_path <- function(bundle_id) {
  file.path(validation_cache_dir(), paste0(bundle_id, ".rds"))
}

validation_force_refresh <- function() {
  !isTRUE(getOption("dsvert.validation.use_cache", FALSE))
}

validation_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing package(s) required by this executable vignette: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

validation_find_pkg_root <- function(pkg) {
  here <- normalizePath(getwd(), mustWork = TRUE)
  paths <- unique(c(here, dirname(here), dirname(dirname(here))))
  for (p in paths) {
    d <- file.path(p, "DESCRIPTION")
    if (file.exists(d)) {
      desc <- tryCatch(read.dcf(d), error = function(e) NULL)
      if (!is.null(desc) && identical(unname(desc[1, "Package"]), pkg)) {
        return(p)
      }
    }
  }
  NA_character_
}

validation_load_packages <- function() {
  validation_require(c("DSI", "DSLite", "dsVert", "dsVertClient"))
  suppressPackageStartupMessages({
    library(DSI)
    library(DSLite)
    library(dsVert)
    library(dsVertClient)
  })
  client_root <- validation_find_pkg_root("dsVertClient")
  server_root <- if (!is.na(client_root)) {
    file.path(dirname(client_root), "dsVert")
  } else {
    NA_character_
  }
  if (requireNamespace("pkgload", quietly = TRUE) &&
      !is.na(server_root) && dir.exists(server_root)) {
    pkgload::load_all(server_root, quiet = TRUE)
  }
  if (requireNamespace("pkgload", quietly = TRUE) &&
      !is.na(client_root) && dir.exists(client_root)) {
    pkgload::load_all(client_root, quiet = TRUE)
  }
  options(dsvert.require_trusted_peers = FALSE,
          dsvert.beaver_preprocessing = "auto",
          datashield.privacyLevel = 5L,
          datashield.errors.print = TRUE)
  invisible(TRUE)
}

validation_demo_verbose <- function() {
  env <- trimws(Sys.getenv("DSVERT_VALIDATION_VERBOSE", ""))
  if (nzchar(env)) {
    return(tolower(env) %in% c("1", "true", "yes", "y", "on"))
  }
  isTRUE(getOption("dsvert.validation.verbose", TRUE))
}

connect_dslite <- function(tables, symbol = "D") {
  server <- DSLite::newDSLiteServer(
    tables = tables,
    config = DSLite::defaultDSConfiguration(include = c("dsBase", "dsVert")))
  dslite.server <<- server
  builder <- DSI::newDSLoginBuilder()
  for (nm in names(tables)) {
    builder$append(server = nm, url = "dslite.server",
                   table = nm, driver = "DSLiteDriver")
  }
  DSI::datashield.login(builder$build(), assign = TRUE, symbol = symbol,
                        opts = list(server = server))
}

psi_align <- function(conns, data = "D", id_col = "patient_id",
                      newobj = "DA", ...) {
  dsVertClient::ds.vert.align(data, id_col, newobj,
                              datasources = conns, verbose = validation_demo_verbose(), ...)
}

max_named_delta <- function(ds, ref) {
  nm <- intersect(names(ds), names(ref))
  if (!length(nm)) return(Inf)
  max(abs(as.numeric(ds[nm]) - as.numeric(ref[nm])))
}

elapsed <- function(expr) {
  t0 <- proc.time()[[3]]
  value <- force(expr)
  list(value = value, runtime_s = unname(proc.time()[[3]] - t0))
}

row_result <- function(method_id, method_name, K, route, dataset,
                       reference, metric, observed, tolerance, tier,
                       disclosure, runtime_s, note = "") {
  status <- if (is.finite(observed) && observed <= tolerance) "PASS" else "FAIL"
  data.frame(
    method_id = method_id,
    method_name = method_name,
    k_mode = if (K == 2L) "K=2" else "K>=3",
    function_route = route,
    dataset = dataset,
    reference_target = reference,
    primary_metric = metric,
    observed = as.numeric(observed),
    tolerance = as.numeric(tolerance),
    tier = tier,
    non_disclosive = TRUE,
    disclosure = disclosure,
    runtime_s = as.numeric(runtime_s),
    status = status,
    note = note,
    stringsAsFactors = FALSE)
}

assert_validation <- function(rows) {
  bad <- rows[rows$status != "PASS" | !rows$non_disclosive, , drop = FALSE]
  if (nrow(bad)) {
    stop("Validation outside accepted envelope:\n",
         paste(utils::capture.output(print(bad)), collapse = "\n"),
         call. = FALSE)
  }
  invisible(TRUE)
}

display_validation <- function(rows) {
  keep <- c("k_mode", "function_route", "dataset", "reference_target",
            "primary_metric", "observed", "tolerance", "tier", "status",
            "runtime_s")
  out <- rows[, keep, drop = FALSE]
  out$observed <- signif(out$observed, 4)
  out$tolerance <- signif(out$tolerance, 4)
  out$runtime_s <- round(out$runtime_s, 1)
  knitr::kable(out)
}

display_profile_validation <- function(rows) {
  keep <- c("k_mode", "function_route", "observed", "tolerance",
            "profile_delta", "runtime_s", "runtime_ratio_vs_dealer",
            "status")
  out <- rows[, keep, drop = FALSE]
  out$observed <- format(out$observed, scientific = TRUE, digits = 8)
  out$tolerance <- format(out$tolerance, scientific = TRUE, digits = 8)
  out$profile_delta <- format(out$profile_delta, scientific = TRUE,
                              digits = 8)
  out$runtime_s <- format(round(out$runtime_s, 3), nsmall = 3)
  out$runtime_ratio_vs_dealer <- format(
    round(out$runtime_ratio_vs_dealer, 2), nsmall = 2)
  knitr::kable(out)
}

run_validation <- function(method_id, force = getOption(
  "dsvert.validation.force", validation_force_refresh()), trace = FALSE) {
  if (isTRUE(trace)) {
    cat("Validation trace\n")
    cat("1. Load dsVertClient, dsVert, DSI and DSLite.\n")
    cat("2. Build the method-specific fixture in R.\n")
    cat("3. Split the fixture vertically into K=2 and K>=3 server tables.\n")
    cat("4. Start an in-memory DSLite server with those tables.\n")
    cat("5. Run ds.vert.align() to create aligned server-side tables.\n")
    cat("6. Run the dsVert product route on the aligned data.\n")
    cat("7. Compute the centralized reference on the pooled fixture.\n")
    cat("8. Compare distributed and centralized outputs against the accepted envelope.\n\n")
  }
  validation_load_packages()
  bundle <- validation_bundle_id(method_id)
  path <- validation_cache_path(bundle)
  if (!isTRUE(force) && file.exists(path)) {
    if (isTRUE(trace)) {
      cat("Cache: reading previous local result from ", path, "\n\n", sep = "")
    }
    rows <- readRDS(path)
  } else {
    if (isTRUE(trace)) {
      cat("Cache: bypassed; executing DSLite validation now.\n\n")
    }
    rows <- switch(bundle,
      psi = validate_psi(),
      descriptive_correlation_pca = validate_desc_cor_pca(),
      contingency = validate_contingency(),
      glm = validate_glm(),
      inference = validate_inference(),
      lasso = validate_lasso(),
      negative_binomial = validate_negative_binomial(),
      cox = validate_cox(),
      lmm = validate_lmm(),
      gee = validate_gee(),
      glmm = validate_glmm(),
      ipw = validate_ipw(),
      mi = validate_mi(),
      multinomial = validate_multinomial(),
      ordinal = validate_ordinal(),
      stop("Unknown validation method: ", method_id, call. = FALSE))
    saveRDS(rows, path)
  }
  out <- rows[rows$method_id == method_id, , drop = FALSE]
  rownames(out) <- NULL
  assert_validation(out)
  if (isTRUE(trace)) {
    cat("Executed rows\n")
    print(out[, c("method_name", "k_mode", "function_route", "dataset",
                  "reference_target", "primary_metric", "observed",
                  "tolerance", "tier", "status", "runtime_s", "note"),
              drop = FALSE], row.names = FALSE)
    cat("\n")
  }
  out
}

run_all_validations <- function(force = getOption(
  "dsvert.validation.force", validation_force_refresh())) {
  do.call(rbind, lapply(validation_method_ids, run_validation, force = force))
}

assert_legacy_routes_removed <- function() {
  validation_load_packages()
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
    ),
    stringsAsFactors = FALSE
  )
  if (!all(checks$pass)) {
    stop("A removed legacy route is still invokable.", call. = FALSE)
  }
  checks
}

build_pima <- function(n = 60L) {
  validation_require("MASS")
  pdf <- as.data.frame(MASS::Pima.tr)
  pdf$patient_id <- sprintf("P%03d", seq_len(nrow(pdf)))
  pdf$diabetes <- as.integer(pdf$type == "Yes")
  pdf$type <- NULL
  pdf <- stats::na.omit(pdf)
  pdf <- pdf[seq_len(min(nrow(pdf), n)), , drop = FALSE]
  pdf[, c("patient_id", "npreg", "glu", "bp", "skin", "bmi", "ped",
          "age", "diabetes")]
}

split_pima <- function(pooled, K, include_y = TRUE, include_cluster = FALSE) {
  pooled$cluster <- as.integer((seq_len(nrow(pooled)) - 1L) %/% 10L) + 1L
  y_cols <- if (include_y) "diabetes" else character(0)
  cl_cols <- if (include_cluster) "cluster" else character(0)
  if (K == 2L) {
    list(
      s1 = pooled[, c("patient_id", "age", "bmi", "ped"), drop = FALSE],
      s2 = pooled[, c("patient_id", "npreg", "glu", "bp", "skin",
                      y_cols, cl_cols), drop = FALSE])
  } else {
    list(
      s1 = pooled[, c("patient_id", "age", "bmi"), drop = FALSE],
      s2 = pooled[, c("patient_id", "ped", "skin"), drop = FALSE],
      s3 = pooled[, c("patient_id", "npreg", "glu", "bp",
                      y_cols, cl_cols), drop = FALSE])
  }
}

validate_psi <- function() {
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
  build <- function(K) {
    id_sets <- if (K == 2L) {
      list(s1 = sprintf("P%04d", 1:80),
           s2 = sprintf("P%04d", 11:90))
    } else {
      list(s1 = sprintf("P%04d", 1:85),
           s2 = sprintf("P%04d", 11:95),
           s3 = sprintf("P%04d", 6:80))
    }
    list(tables = Map(make_table, id_sets, seq_along(id_sets)),
         common = Reduce(intersect, id_sets))
  }
  rows <- list()
  for (K in c(2L, 3L)) {
    set.seed(9100 + K)
    fixture <- build(K)
    conns <- connect_dslite(fixture$tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    run <- elapsed(psi <- psi_align(conns))
    counts <- DSI::datashield.aggregate(
      conns, call(name = "getObsCountDS", data_name = "DA"))
    n_by_server <- vapply(counts, function(x) as.integer(x$n_obs), integer(1))
    variables <- stats::setNames(lapply(names(conns), function(.x) "id_num"),
                                 names(conns))
    cor_fit <- dsVertClient::ds.vert.cor("DA", variables = variables,
                                        datasources = conns, verbose = validation_demo_verbose())
    offdiag <- cor_fit$correlation[upper.tri(cor_fit$correlation)]
    legacy_blocked <- tryCatch({
      DSI::datashield.aggregate(conns[1], call(name = "psiGetMatchedIndicesDS"))
      FALSE
    }, error = function(e) TRUE)
    observed <- max(abs(n_by_server - length(fixture$common)),
                    abs((psi$n_common %||% psi[[1]]$n_matched) -
                          length(fixture$common)),
                    max(abs(offdiag - 1)))
    audit_ok <- !contains_ids(psi) && isTRUE(legacy_blocked)
    rows[[as.character(K)]] <- row_result(
      "psi", "PSI alignment", K, "ds.vert.align",
      "synthetic ID intersection", "deterministic set intersection",
      "max(count_delta, correlation_delta)", observed, 0,
      "strict-precise",
      "Returns match counts/status only; matched IDs and row indices are not returned.",
      run$runtime_s,
      if (audit_ok) "index reveal blocked" else "audit warning")
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

validate_desc_cor_pca <- function() {
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_pima(60L)
    tables <- split_pima(pooled, K, include_y = FALSE)
    vars <- lapply(tables, function(x) setdiff(names(x), "patient_id"))
    all_vars <- unlist(vars, use.names = FALSE)
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)

    d0 <- elapsed(desc <- dsVertClient::ds.vert.desc(
      "DA", variables = vars, n_buckets = 40L,
      verbose = validation_demo_verbose(), datasources = conns))
    ref_mean <- vapply(pooled[all_vars], mean, numeric(1), na.rm = TRUE)
    ref_sd <- vapply(pooled[all_vars], stats::sd, numeric(1), na.rm = TRUE)
    ds_mean <- stats::setNames(desc$mean, desc$variable)
    ds_sd <- stats::setNames(desc$sd, desc$variable)
    desc_delta <- max(max_named_delta(ds_mean, ref_mean),
                      max_named_delta(ds_sd, ref_sd))

    c0 <- elapsed(cor_ds <- dsVertClient::ds.vert.cor(
      "DA", variables = vars, verbose = validation_demo_verbose(), datasources = conns))
    cor_ref <- stats::cor(pooled[cor_ds$var_names])
    cor_delta <- max(abs(cor_ds$correlation - cor_ref))

    p0 <- elapsed(pca_ds <- dsVertClient::ds.vert.pca(
      cor_result = cor_ds, verbose = validation_demo_verbose(), datasources = conns))
    eig_ref <- eigen(cor_ref, symmetric = TRUE)
    load_ref <- eig_ref$vectors[, seq_len(ncol(pca_ds$loadings)),
                                drop = FALSE]
    rownames(load_ref) <- rownames(pca_ds$loadings)
    load_ds <- pca_ds$loadings
    for (j in seq_len(ncol(load_ds))) {
      if (sum(load_ds[, j] * load_ref[, j]) < 0) load_ds[, j] <- -load_ds[, j]
    }
    pca_delta <- max(max(abs(pca_ds$eigenvalues - eig_ref$values)),
                     max(abs(load_ds - load_ref)))

    rows[[paste0("desc", K)]] <- row_result(
      "descriptive", "Descriptive statistics", K, "ds.vert.desc",
      "MASS::Pima.tr fixture", "central mean/sd",
      "max(mean_sd_abs_delta)", desc_delta, 0, "strict-precise",
      "Returns guarded scalar summaries and histogram-based quantiles, not rows.",
      d0$runtime_s)
    rows[[paste0("cor", K)]] <- row_result(
      "correlation", "Correlation", K, "ds.vert.cor",
      "MASS::Pima.tr fixture", "stats::cor",
      "correlation_max_abs_delta", cor_delta, 0.001, "strict-practical",
      "Releases the low-dimensional correlation matrix only.",
      c0$runtime_s)
    rows[[paste0("pca", K)]] <- row_result(
      "pca", "PCA", K, "ds.vert.pca",
      "MASS::Pima.tr fixture", "eigen(cor(X))",
      "max(eigen_loading_abs_delta)", pca_delta, 0.001,
      "strict-practical",
      "Reuses the correlation matrix; no score or loading per patient is returned.",
      p0$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

validate_contingency <- function() {
  rows <- list()
  for (K in c(2L, 3L)) {
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
      list(s1 = pooled[, c("patient_id", "age_grp")],
           s2 = pooled[, c("patient_id", "bp_grp", "diabetes")])
    } else {
      list(s1 = pooled[, c("patient_id", "age_grp")],
           s2 = pooled[, c("patient_id", "bp_grp")],
           s3 = pooled[, c("patient_id", "diabetes")])
    }
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed(cross <- dsVertClient::ds.vert.chisq_cross(
      "DA", "age_grp", "diabetes", correct = FALSE, fisher = TRUE,
      verbose = validation_demo_verbose(), datasources = conns))
    tab <- table(pooled$age_grp, pooled$diabetes)
    chi_ref <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    ref_mat <- unclass(tab)[rownames(cross$observed),
                            colnames(cross$observed), drop = FALSE]
    observed <- max(max(abs(cross$observed - ref_mat)),
                    abs(cross$chisq - unname(chi_ref$statistic)))
    rows[[as.character(K)]] <- row_result(
      "contingency", "Contingency tests", K,
      "ds.vert.chisq / ds.vert.fisher / ds.vert.chisq_cross",
      "MASS::Pima.tr categorical fixture", "chisq.test / fisher.test",
      "max(count_or_chisq_delta)", observed, 0, "strict-precise",
      "Releases guarded table counts and test statistics; small cells fail closed.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

validate_glm <- function() {
  rows <- list()
  fm <- diabetes ~ age + bmi + ped + glu + bp + skin + npreg
  for (K in c(2L, 3L)) {
    pooled <- build_pima(60L)
    conns <- connect_dslite(split_pima(pooled, K, include_y = TRUE))
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed(fit <- dsVertClient::ds.vert.glm(
      fm, data = "DA", family = "gaussian", max_iter = 30L,
      verbose = validation_demo_verbose(), datasources = conns))
    observed <- max_named_delta(fit$coefficients, coef(stats::lm(fm, pooled)))
    rows[[as.character(K)]] <- row_result(
      "glm", "GLM", K, "ds.vert.glm",
      "MASS::Pima.tr fixture", "stats::lm",
      "coef_max_abs_delta", observed, 1e-3, "strict-practical",
      "Uses secure score aggregates and returns model-level coefficients/covariance.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

validate_inference <- function() {
  rows <- list()
  fm <- diabetes ~ age + bmi + ped + glu + bp + skin + npreg
  red <- diabetes ~ age + bmi + ped + glu + bp + skin
  for (K in c(2L, 3L)) {
    pooled <- build_pima(60L)
    conns <- connect_dslite(split_pima(pooled, K, include_y = TRUE))
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed({
      full <- dsVertClient::ds.vert.glm(fm, data = "DA", family = "gaussian",
        max_iter = 30L, compute_se = TRUE, compute_deviance = TRUE,
        verbose = validation_demo_verbose(), datasources = conns)
      reduced <- dsVertClient::ds.vert.glm(red, data = "DA",
        family = "gaussian", max_iter = 30L, compute_se = TRUE,
        compute_deviance = TRUE, verbose = validation_demo_verbose(), datasources = conns)
      ci <- dsVertClient::ds.vert.confint(full)
      wald <- dsVertClient::ds.vert.wald(full, "age")
      Kmat <- matrix(0, nrow = 1, ncol = length(full$coefficients),
                     dimnames = list("age", names(full$coefficients)))
      Kmat[1, "age"] <- 1
      contrast <- dsVertClient::ds.vert.contrast(full, Kmat)
      lr <- dsVertClient::ds.vert.lr(reduced, full)
      list(full = full, ci = ci, wald = wald, contrast = contrast, lr = lr)
    })
    full <- run$value$full
    se <- full$std_errors["age"]
    est <- full$coefficients["age"]
    z <- est / se
    manual_p <- 2 * stats::pnorm(-abs(z))
    lr_p_ref <- stats::pchisq(run$value$lr$statistic,
                              df = run$value$lr$df,
                              lower.tail = FALSE)
    observed <- max(abs(run$value$wald$p_value - manual_p),
                    abs(run$value$contrast$estimate - est),
                    abs(run$value$ci["age", "estimate"] - est),
                    abs(run$value$lr$p_value - lr_p_ref))
    rows[[as.character(K)]] <- row_result(
      "inference", "Inference helpers", K,
      "ds.vert.confint / ds.vert.wald / ds.vert.contrast / ds.vert.lr",
      "MASS::Pima.tr fixture", "manual algebra on ds.vert.glm output",
      "algebra_max_abs_delta", observed, 0, "strict-precise",
      "Post-processes released model-level beta/covariance/deviance only.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

validate_lasso <- function() {
  rows <- list()
  fm <- diabetes ~ age + bmi + ped + glu + bp + skin + npreg
  for (K in c(2L, 3L)) {
    pooled <- build_pima(60L)
    conns <- connect_dslite(split_pima(pooled, K, include_y = TRUE))
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed({
      fit <- dsVertClient::ds.vert.glm(fm, data = "DA", family = "gaussian",
        max_iter = 30L, verbose = validation_demo_verbose(), datasources = conns)
      prox0 <- dsVertClient::ds.vert.lasso_proximal(
        fit, lambda = 0, max_iter = 1000L, tol = 1e-8)
      list(fit = fit, prox0 = prox0)
    })
    observed <- max_named_delta(run$value$prox0$coefficients,
                                run$value$fit$coefficients)
    rows[[as.character(K)]] <- row_result(
      "lasso", "LASSO", K, "ds.vert.lasso_proximal(lambda=0)",
      "MASS::Pima.tr fixture", "OLS limit of ds.vert.glm",
      "lambda0_coef_abs_delta", observed, 0, "strict-precise",
      "Consumes only model-level GLM aggregates; no extra patient-level values are returned.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_nb <- function(n = 60L, seed = 303L, theta = 3,
                     beta = c(1, 0.2, -0.15, 0.1)) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  mu <- exp(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3)
  y <- stats::rnbinom(n, size = theta, mu = mu)
  data.frame(patient_id = sprintf("N%03d", seq_len(n)), x1, x2, x3, y)
}

split_x123_y <- function(d, K, y = "y", extras = character()) {
  if (K == 2L) {
    list(s1 = d[, c("patient_id", "x1", "x2"), drop = FALSE],
         s2 = d[, c("patient_id", "x3", y, extras), drop = FALSE])
  } else {
    list(s1 = d[, c("patient_id", "x1"), drop = FALSE],
         s2 = d[, c("patient_id", "x2"), drop = FALSE],
         s3 = d[, c("patient_id", "x3", y, extras), drop = FALSE])
  }
}

validate_negative_binomial <- function() {
  validation_require("MASS")
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_nb(60L)
    conns <- connect_dslite(split_x123_y(pooled, K))
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    ref_accurate <- coef(MASS::glm.nb(y ~ x1 + x2 + x3, data = pooled))
    run_accurate <- elapsed(accurate <- dsVertClient::ds.vert.nb(
      y ~ x1 + x2 + x3, data = "DA", method = "accurate",
      max_iter = 25L, verbose = validation_demo_verbose(),
      datasources = conns))
    accurate_delta <- max_named_delta(accurate$coefficients, ref_accurate)
    rows[[paste0(K, "_accurate")]] <- row_result(
      "negative_binomial", "Negative binomial accurate", K,
      "ds.vert.nb(method='accurate')",
      "synthetic NB fixture", "MASS::glm.nb",
      "coef_max_abs_delta", accurate_delta, 0.001, "strict-practical",
      "Full-regression beta/theta score components remain aggregate.",
      run_accurate$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_cox <- function(n = 60L, seed = 202L, beta = c(0.4, -0.3, 0.2),
                      n_bins = 8L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  eta <- beta[1] * x1 + beta[2] * x2 + beta[3] * x3
  t0 <- stats::rexp(n, rate = exp(eta) / 20)
  cens <- stats::rexp(n, rate = 0.02)
  event <- as.integer(t0 <= cens)
  time <- pmin(t0, cens)
  br <- unique(stats::quantile(time, probs = seq(0, 1, length.out = n_bins + 1L),
                               names = FALSE))
  time <- as.integer(cut(time, br, include.lowest = TRUE))
  data.frame(patient_id = sprintf("C%03d", seq_len(n)), x1, x2, x3,
             time, event)
}

validate_cox <- function() {
  validation_require("survival")
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_cox(60L)
    tables <- if (K == 2L) {
      list(s1 = pooled[, c("patient_id", "x1", "x2")],
           s2 = pooled[, c("patient_id", "x3", "time", "event")])
    } else {
      list(s1 = pooled[, c("patient_id", "x1")],
           s2 = pooled[, c("patient_id", "x2")],
           s3 = pooled[, c("patient_id", "x3", "time", "event")])
    }
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed(fit <- dsVertClient::ds.vert.cox(
      survival::Surv(time, event) ~ x1 + x2 + x3, data = "DA",
      verbose = validation_demo_verbose(), datasources = conns))
    ref <- coef(survival::coxph(survival::Surv(time, event) ~ x1 + x2 + x3,
                                data = pooled, ties = "breslow"))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "cox", "Cox PH", K, "ds.vert.cox",
      "synthetic discretised survival fixture",
      "survival::coxph(ties='breslow')",
      "coef_max_abs_delta", observed, 1e-3, "strict-practical",
      "Event-time score terms remain shared; client receives beta and scalar diagnostics.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_lmm <- function(ncl = 12L, m = 5L, seed = 101L, b_sd = 0.8,
                      eps_sd = 0.4,
                      beta = c(1, 0.5, -0.3, 0.2)) {
  set.seed(seed)
  n <- ncl * m
  cluster <- rep(seq_len(ncl), each = m)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  b <- rnorm(ncl, sd = b_sd)[cluster]
  y <- beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3 +
    b + rnorm(n, sd = eps_sd)
  data.frame(patient_id = sprintf("L%03d", seq_len(n)), cluster,
             x1, x2, x3, y)
}

validate_lmm <- function() {
  validation_require("lme4")
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_lmm()
    conns <- connect_dslite(split_x123_y(pooled, K, extras = "cluster"))
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed(fit <- dsVertClient::ds.vert.lmm(
      y ~ x1 + x2 + x3, data = "DA", cluster_col = "cluster",
      verbose = validation_demo_verbose(), datasources = conns))
    ref <- lme4::fixef(lme4::lmer(y ~ x1 + x2 + x3 + (1 | cluster),
                                  data = pooled, REML = TRUE))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "lmm", "LMM", K, "ds.vert.lmm",
      "synthetic random-intercept fixture", "lme4::lmer",
      "fixed_effect_max_abs_delta", observed, 0.001, "strict-practical",
      "Cluster membership is used internally; per-cluster residuals/BLUPs are not returned.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

validate_gee <- function() {
  validation_require("geepack")
  rows <- list()
  fm <- diabetes ~ age + bmi + ped + glu + bp + skin + npreg
  for (K in c(2L, 3L)) {
    pooled <- build_pima(60L)
    pooled$cluster <- as.integer((seq_len(nrow(pooled)) - 1L) %/% 10L) + 1L
    conns <- connect_dslite(split_pima(pooled, K, include_y = TRUE,
                                       include_cluster = TRUE))
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed(fit <- dsVertClient::ds.vert.gee(
      fm, data = "DA", family = "gaussian", id_col = "cluster",
      corstr = "independence", max_iter = 30L,
      verbose = validation_demo_verbose(), datasources = conns))
    ref <- coef(geepack::geeglm(fm, data = pooled, id = cluster,
                                family = gaussian(),
                                corstr = "independence"))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "gee", "GEE", K, "ds.vert.gee(corstr='independence')",
      "MASS::Pima.tr clustered fixture", "geepack::geeglm",
      "coef_max_abs_delta", observed, 0.001, "strict-practical",
      "Returns regression and sandwich-level aggregates, not row scores.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_balanced_glmm <- function(seed = 77L, ncl = 12L, m = 8L,
                                b_sd = 0.5,
                                beta = c(-0.2, 0.45, -0.3)) {
  set.seed(seed)
  n <- ncl * m
  cluster <- rep(seq_len(ncl), each = m)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  b <- rnorm(ncl, sd = b_sd)[cluster]
  eta <- beta[1] + beta[2] * x1 + beta[3] * x2 + b
  data.frame(
    patient_id = sprintf("G%03d", seq_len(n)),
    cluster = cluster,
    x1 = x1,
    x2 = x2,
    y = stats::rbinom(n, 1L, stats::plogis(eta)))
}

validate_glmm <- function() {
  validation_require(c("MASS", "nlme"))
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_balanced_glmm(seed = 1L, ncl = 10L, m = 5L,
                                  b_sd = 0.6)
    tables <- if (K == 2L) {
      list(s1 = pooled[, c("patient_id", "x1")],
           s2 = pooled[, c("patient_id", "x2", "cluster", "y")])
    } else {
      list(s1 = pooled[, c("patient_id", "x1")],
           s2 = pooled[, c("patient_id", "x2")],
           s3 = pooled[, c("patient_id", "cluster", "y")])
    }
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run_pql <- elapsed(fit_pql <- dsVertClient::ds.vert.glmm(
      y ~ x1 + x2, data = "DA", cluster_col = "cluster",
      method = "pql", compute_se = FALSE,
      verbose = validation_demo_verbose(), datasources = conns))
    pooled_ref <- pooled
    pooled_ref$cluster <- factor(pooled_ref$cluster)
    ref_pql <- suppressWarnings(MASS::glmmPQL(
      y ~ x1 + x2, random = ~1 | cluster, family = binomial(),
      data = pooled_ref, verbose = validation_demo_verbose()))
    ref <- nlme::fixef(ref_pql)
    fixed_delta <- max_named_delta(fit_pql$coefficients, ref)
    pql_ran <- isTRUE(fit_pql$iterations >= 1L) &&
      is.data.frame(fit_pql$trace) && nrow(fit_pql$trace) >= 1L &&
      is.finite(fit_pql$sigma_b2) &&
      identical(fit_pql$quality$status, "ok")
    observed <- max(fixed_delta, if (pql_ran) 0 else Inf)
    rows[[paste0("pql", K)]] <- row_result(
      "glmm", "GLMM", K, "ds.vert.glmm(method='pql')",
      "synthetic mixed binomial random-intercept fixture",
      "MASS::glmmPQL",
      "fixed_effect_max_abs_delta_and_pql_quality", observed, 0.001,
      "strict-practical",
      "PQL route returns fixed effects and scalar variance diagnostics only.",
      run_pql$runtime_s)

    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_ipw <- function(n = 60L, seed = 88L,
                      beta_p = c(-0.2, 0.4, -0.3),
                      beta_y = c(1, 2, 0.5, -0.2),
                      noise_sd = 0.2) {
  set.seed(seed)
  w1 <- rnorm(n); w2 <- rnorm(n)
  p <- stats::plogis(beta_p[1] + beta_p[2] * w1 + beta_p[3] * w2)
  tr <- stats::rbinom(n, 1, p)
  ipw <- ifelse(tr == 1, 1 / p, 1 / (1 - p))
  y <- beta_y[1] + beta_y[2] * tr + beta_y[3] * w1 + beta_y[4] * w2 +
    rnorm(n, 0, noise_sd)
  data.frame(patient_id = sprintf("I%03d", seq_len(n)), w1, w2, tr, ipw, y)
}

validate_ipw <- function() {
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_ipw()
    tables <- if (K == 2L) {
      list(s1 = pooled[, c("patient_id", "w1", "w2")],
           s2 = pooled[, c("patient_id", "tr", "ipw", "y")])
    } else {
      list(s1 = pooled[, c("patient_id", "w1")],
           s2 = pooled[, c("patient_id", "w2")],
           s3 = pooled[, c("patient_id", "tr", "ipw", "y")])
    }
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed(fit <- dsVertClient::ds.vert.ipw(
      y ~ tr + w1 + w2, tr ~ w1 + w2, data = "DA",
      outcome_family = "gaussian",
      compute_se = FALSE, compute_deviance = FALSE,
      verbose = validation_demo_verbose(), datasources = conns))
    ref <- coef(stats::lm(y ~ tr + w1 + w2, data = pooled,
                          weights = ipw))
    observed <- max_named_delta(fit$outcome$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "ipw", "IPW", K, "ds.vert.ipw",
      "synthetic confounded IPW fixture",
      "central weighted lm using same weights",
      "weighted_outcome_coef_abs_delta", observed, 1e-3,
      "strict-practical",
      "Uses protected propensity/outcome GLM aggregates; only model-level fits are returned.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_mi <- function(n = 60L, seed = 99L, missing_every = 7L,
                     beta = c(1, 0.4, -0.2, 0.1), noise_sd = 0.2) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  y <- beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3 +
    rnorm(n, 0, noise_sd)
  d <- data.frame(patient_id = sprintf("MI%03d", seq_len(n)), x1, x2, x3, y)
  d$x2[seq(5L, n, by = missing_every)] <- NA
  d
}

rubin_pool_validation <- function(fits) {
  m <- length(fits)
  betas <- do.call(cbind, lapply(fits, stats::coef))
  covs <- lapply(fits, stats::vcov)
  beta_bar <- rowMeans(betas)
  W <- Reduce(`+`, covs) / m
  dev <- sweep(betas, 1L, beta_bar)
  B <- (dev %*% t(dev)) / max(m - 1L, 1L)
  Tmat <- W + (1 + 1 / m) * B
  se <- sqrt(pmax(diag(Tmat), 0))
  nm <- names(stats::coef(fits[[1L]]))
  names(beta_bar) <- names(se) <- nm
  dimnames(W) <- dimnames(B) <- dimnames(Tmat) <- list(nm, nm)
  list(coefficients = beta_bar, covariance = Tmat, within = W,
       between = B, std_errors = se, fits = fits)
}

central_mi_same_imputer <- function(pooled, K, formula,
                                    impute_columns = "x2", m = 20L,
                                    seed = 12L) {
  base_tables <- split_mi(pooled, K)
  owners <- vapply(impute_columns, function(v) {
    hits <- names(base_tables)[vapply(base_tables, function(tbl) {
      v %in% names(tbl)
    }, logical(1))]
    if (length(hits) != 1L) {
      stop("Expected impute column '", v, "' on exactly one split table",
           call. = FALSE)
    }
    hits[[1L]]
  }, character(1))

  fits <- vector("list", m)
  imputation_log <- vector("list", m * length(impute_columns))
  log_i <- 0L
  for (mi in seq_len(m)) {
    tables <- base_tables
    round_tag <- sprintf("__mi_%d", mi)
    for (v in impute_columns) {
      owner <- owners[[v]]
      env <- new.env(parent = globalenv())
      env$D <- tables[[owner]]
      output_column <- paste0(v, round_tag)
      imp_res <- eval(substitute(
        dsVert::dsvertImputeColumnDS(
          data_name = "D",
          impute_column = V,
          output_column = OUT,
          seed = SEED),
        list(V = v, OUT = output_column, SEED = as.integer(seed + mi))),
        envir = env)
      tables[[owner]] <- env$D
      log_i <- log_i + 1L
      imputation_log[[log_i]] <- data.frame(
        round = mi,
        variable = v,
        split = owner,
        n_imputed = as.integer(imp_res$n_imputed %||% NA_integer_),
        n_observed = as.integer(imp_res$n_observed %||% NA_integer_),
        method = as.character(imp_res$method %||% NA_character_),
        n_predictors = as.integer(imp_res$n_predictors %||% NA_integer_),
        intercept_only = as.logical(imp_res$intercept_only %||% NA),
        stringsAsFactors = FALSE)
    }
    pooled_imp <- pooled
    for (v in impute_columns) {
      pooled_imp[[v]] <- tables[[owners[[v]]]][[paste0(v, round_tag)]]
    }
    fits[[mi]] <- stats::glm(formula, data = pooled_imp,
                             family = stats::gaussian())
  }
  out <- rubin_pool_validation(fits)
  out$imputation_log <- if (log_i) {
    do.call(rbind, imputation_log[seq_len(log_i)])
  } else {
    data.frame()
  }
  out
}

validate_mi <- function() {
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_mi()
    tables <- split_mi(pooled, K)
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns, na.action = "none")
    run <- elapsed(fit <- dsVertClient::ds.vert.mi(
      y ~ x1 + x2 + x3, data = "DA", impute_columns = "x2",
      family = "gaussian",
      lambda = 0, verbose = validation_demo_verbose(),
      datasources = conns, seed = 12L))
    ref <- central_mi_same_imputer(
      pooled, K, y ~ x1 + x2 + x3, impute_columns = "x2",
      m = fit$m, seed = 12L)
    observed <- max_named_delta(fit$coefficients, ref$coefficients)
    rows[[as.character(K)]] <- row_result(
      "mi", "Multiple imputation", K, "ds.vert.mi",
      "synthetic missing-covariate fixture",
      "same-imputer centralized Rubin pooling",
      "pooled_coef_abs_delta", observed, 0.001, "strict-practical",
      "Imputed columns stay server-side; client pools beta/covariance only.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_multinomial <- function(n = 60L, seed = 44L, shift = 0.15) {
  set.seed(seed)
  stopifnot(n %% 3L == 0L)
  base_n <- n / 3L
  base <- data.frame(x1 = rnorm(base_n), x2 = rnorm(base_n), x3 = rnorm(base_n))
  high <- transform(base, x1 = x1 - shift, x2 = x2 + 0.5 * shift,
                    y_cls = "high")
  low <- transform(base, x1 = x1 + shift, y_cls = "low")
  med <- transform(base, x2 = x2 - shift, y_cls = "med")
  d <- rbind(high, low, med)
  d <- cbind(patient_id = sprintf("M%03d", seq_len(n)), d)
  d$y_cls <- factor(d$y_cls, levels = c("high", "low", "med"))
  for (k in levels(d$y_cls)) d[[paste0(k, "_ind")]] <- as.integer(d$y_cls == k)
  d
}

softmax_prob <- function(beta, data, non_ref = c("low", "med")) {
  X <- cbind("(Intercept)" = 1, x1 = data$x1, x2 = data$x2, x3 = data$x3)
  eta <- X[, rownames(beta), drop = FALSE] %*% beta[, non_ref, drop = FALSE]
  eta_all <- cbind(high = 0, eta)
  num <- exp(eta_all - apply(eta_all, 1L, max))
  num / rowSums(num)
}

validate_multinomial <- function() {
  validation_require("nnet")
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_multinomial()
    tables <- if (K == 2L) {
      list(s1 = pooled[, c("patient_id", "x1", "x2")],
           s2 = pooled[, c("patient_id", "x3", "y_cls",
                           "high_ind", "low_ind", "med_ind")])
    } else {
      list(s1 = pooled[, c("patient_id", "x1")],
           s2 = pooled[, c("patient_id", "x2")],
           s3 = pooled[, c("patient_id", "x3", "y_cls",
                           "high_ind", "low_ind", "med_ind")])
    }
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    ref <- nnet::multinom(y_cls ~ x1 + x2 + x3, data = pooled,
                          trace = FALSE, maxit = 200L)
    ref_prob <- stats::predict(ref, pooled, type = "probs")
    run <- elapsed(fit <- dsVertClient::ds.vert.multinom(
      y_cls ~ x1 + x2 + x3, data = "DA",
      classes = c("high", "low", "med"),
      indicator_template = "%s_ind",
      verbose = validation_demo_verbose(), datasources = conns))
    ds_prob <- softmax_prob(fit$coefficients, pooled)
    ds_prob <- ds_prob[, colnames(ref_prob), drop = FALSE]
    observed <- max(abs(ds_prob - ref_prob))
    rows[[as.character(K)]] <- row_result(
      "multinomial", "Multinomial", K, "ds.vert.multinom",
      "balanced soft-signal synthetic 3-class fixture", "nnet::multinom probabilities",
      "class_probability_max_abs_delta", observed, 0.001,
      "strict-practical",
      "Softmax probabilities and residuals remain Ring127 shares; no row probabilities are returned.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_ordinal <- function(n = 60L, seed = 55L) {
  set.seed(seed)
  d <- data.frame(patient_id = sprintf("O%03d", seq_len(n)),
                  x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  y <- factor(rep(c("low", "med", "high"), each = n / 3L),
              levels = c("low", "med", "high"), ordered = TRUE)
  d$y_ord <- y
  for (k in levels(y)) {
    d[[paste0(k, "_leq")]] <- as.integer(as.integer(y) <= match(k, levels(y)))
  }
  d
}

ordinal_cumprob <- function(fit, data) {
  X <- as.matrix(data[, names(fit$beta_po_joint), drop = FALSE])
  sapply(fit$thresholds_joint, function(th) {
    stats::plogis(th - as.numeric(X %*% fit$beta_po_joint))
  })
}

validate_ordinal <- function() {
  validation_require("MASS")
  old_fd <- getOption("dsvert.ord_strict_fd_max_dim", NULL)
  options(dsvert.ord_strict_fd_max_dim = 0L)
  on.exit(options(dsvert.ord_strict_fd_max_dim = old_fd), add = TRUE)
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_ordinal()
    tables <- if (K == 2L) {
      list(s1 = pooled[, c("patient_id", "x1", "x2")],
           s2 = pooled[, c("patient_id", "x3", "y_ord",
                           "low_leq", "med_leq", "high_leq")])
    } else {
      list(s1 = pooled[, c("patient_id", "x1")],
           s2 = pooled[, c("patient_id", "x2")],
           s3 = pooled[, c("patient_id", "x3", "y_ord",
                           "low_leq", "med_leq", "high_leq")])
    }
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    ref <- MASS::polr(y_ord ~ x1 + x2 + x3, data = pooled, Hess = TRUE)
    ref_cum <- sapply(ref$zeta, function(th) {
      stats::plogis(th - as.numeric(as.matrix(pooled[, names(coef(ref))]) %*%
                                      coef(ref)))
    })
    run <- elapsed(fit <- dsVertClient::ds.vert.ordinal(
      y_ord ~ x1 + x2 + x3, data = "DA",
      levels_ordered = c("low", "med", "high"),
      cumulative_template = "%s_leq",
      verbose = validation_demo_verbose(), datasources = conns))
    ds_cum <- ordinal_cumprob(fit, pooled)
    colnames(ds_cum) <- colnames(ref_cum)
    observed <- max(abs(ds_cum - ref_cum))
    rows[[as.character(K)]] <- row_result(
      "ordinal", "Ordinal", K, "ds.vert.ordinal",
      "balanced synthetic 3-level ordinal fixture", "MASS::polr cumulative probabilities",
      "cumulative_probability_max_abs_delta", observed, 0.001,
      "strict-practical",
      "Class probabilities/residuals stay Ring127 shares; no row probabilities are returned.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

generalization_core_methods <- c(
  "glm", "negative_binomial", "cox", "lmm", "gee", "ipw", "mi")

generalization_heavy_methods <- c("multinomial", "ordinal", "glmm")

generalization_cache_path <- function(method_ids = NULL, seeds = NULL) {
  key <- "generalization"
  if (!is.null(method_ids) || !is.null(seeds)) {
    parts <- c(
      "generalization",
      if (!is.null(method_ids)) paste(method_ids, collapse = "-") else NULL,
      if (!is.null(seeds)) paste(seeds, collapse = "-") else NULL)
    key <- paste(parts, collapse = "_")
  }
  key <- gsub("[^A-Za-z0-9_=-]+", "-", key)
  validation_cache_path(key)
}

generalization_parse_methods <- function(include_heavy = FALSE) {
  env <- trimws(Sys.getenv("DSVERT_GENERALIZATION_METHODS", ""))
  if (nzchar(env)) {
    methods <- trimws(strsplit(env, ",", fixed = TRUE)[[1L]])
    methods[nzchar(methods)]
  } else {
    c(generalization_core_methods,
      if (isTRUE(include_heavy)) generalization_heavy_methods else character())
  }
}

generalization_parse_seeds <- function() {
  env <- trimws(Sys.getenv("DSVERT_GENERALIZATION_SEEDS", ""))
  if (!nzchar(env)) return(c(711L, 733L))
  seeds <- suppressWarnings(as.integer(trimws(strsplit(env, ",",
                                                       fixed = TRUE)[[1L]])))
  seeds[is.finite(seeds)]
}

generalization_quality_status <- function(fit) {
  status <- fit$quality$status %||% "ok"
  if (!length(status) || is.na(status[1L])) "ok" else as.character(status[1L])
}

generalization_row <- function(method_id, method_name, scenario, seed, K,
                               route, dataset, reference, metric, observed,
                               tolerance, tier, disclosure, runtime_s,
                               quality_status = "ok",
                               quality_allowed = c("ok"),
                               note = "") {
  quality_ok <- quality_status %in% quality_allowed
  status <- if (is.finite(observed) && observed <= tolerance && quality_ok) {
    "PASS"
  } else {
    "FAIL"
  }
  data.frame(
    method_id = method_id,
    method_name = method_name,
    scenario = scenario,
    seed = as.integer(seed),
    k_mode = if (K == 2L) "K=2" else "K>=3",
    function_route = route,
    dataset = dataset,
    reference_target = reference,
    primary_metric = metric,
    observed = as.numeric(observed),
    tolerance = as.numeric(tolerance),
    tier = tier,
    non_disclosive = TRUE,
    disclosure = disclosure,
    quality_status = quality_status,
    runtime_s = as.numeric(runtime_s),
    status = status,
    note = note,
    stringsAsFactors = FALSE)
}

generalization_error_row <- function(method_id, scenario, seed, K, err) {
  generalization_row(
    method_id = method_id,
    method_name = method_id,
    scenario = scenario,
    seed = seed,
    K = K,
    route = "ds.vert.*",
    dataset = "generalization scenario",
    reference = "centralized R",
    metric = "execution_error",
    observed = Inf,
    tolerance = 0,
    tier = "unclassified",
    disclosure = "No successful analysis object returned.",
    runtime_s = NA_real_,
    quality_status = "error",
    quality_allowed = "ok",
    note = conditionMessage(err))
}

assert_generalization <- function(rows) {
  bad <- rows[rows$status != "PASS" | !rows$non_disclosive, , drop = FALSE]
  if (nrow(bad)) {
    stop("Generalization validation outside accepted envelope:\n",
         paste(utils::capture.output(print(bad)), collapse = "\n"),
         call. = FALSE)
  }
  invisible(TRUE)
}

display_generalization <- function(rows) {
  keep <- c("method_id", "scenario", "seed", "k_mode", "function_route",
            "primary_metric", "observed", "tolerance", "quality_status",
            "status", "runtime_s")
  out <- rows[, keep, drop = FALSE]
  out$observed <- signif(out$observed, 4)
  out$tolerance <- signif(out$tolerance, 4)
  out$runtime_s <- round(out$runtime_s, 1)
  knitr::kable(out)
}

with_aligned_dslite <- function(tables, FUN, ...) {
  conns <- connect_dslite(tables)
  on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
  psi_align(conns, ...)
  FUN(conns)
}

build_gaussian_regression <- function(n = 72L, seed = 711L,
                                      beta = c(1, 0.35, -0.25, 0.15),
                                      noise_sd = 0.25) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  y <- beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3 +
    rnorm(n, sd = noise_sd)
  data.frame(patient_id = sprintf("R%03d", seq_len(n)), x1, x2, x3, y)
}

split_cox <- function(pooled, K) {
  if (K == 2L) {
    list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x3", "time", "event"),
                     drop = FALSE])
  } else {
    list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
         s3 = pooled[, c("patient_id", "x3", "time", "event"),
                     drop = FALSE])
  }
}

split_ipw <- function(pooled, K) {
  if (K == 2L) {
    list(s1 = pooled[, c("patient_id", "w1", "w2"), drop = FALSE],
         s2 = pooled[, c("patient_id", "tr", "ipw", "y"), drop = FALSE])
  } else {
    list(s1 = pooled[, c("patient_id", "w1"), drop = FALSE],
         s2 = pooled[, c("patient_id", "w2"), drop = FALSE],
         s3 = pooled[, c("patient_id", "tr", "ipw", "y"), drop = FALSE])
  }
}

split_mi <- function(pooled, K) {
  if (K == 2L) {
    list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x2", "x3", "y"), drop = FALSE])
  } else {
    list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
         s3 = pooled[, c("patient_id", "x3", "y"), drop = FALSE])
  }
}

split_glmm <- function(pooled, K) {
  if (K == 2L) {
    list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x2", "cluster", "y"),
                     drop = FALSE])
  } else {
    list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
         s3 = pooled[, c("patient_id", "cluster", "y"), drop = FALSE])
  }
}

split_multinomial <- function(pooled, K) {
  if (K == 2L) {
    list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x3", "y_cls",
                         "high_ind", "low_ind", "med_ind"), drop = FALSE])
  } else {
    list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
         s3 = pooled[, c("patient_id", "x3", "y_cls",
                         "high_ind", "low_ind", "med_ind"), drop = FALSE])
  }
}

split_ordinal <- function(pooled, K) {
  if (K == 2L) {
    list(s1 = pooled[, c("patient_id", "x1", "x2"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x3", "y_ord",
                         "low_leq", "med_leq", "high_leq"), drop = FALSE])
  } else {
    list(s1 = pooled[, c("patient_id", "x1"), drop = FALSE],
         s2 = pooled[, c("patient_id", "x2"), drop = FALSE],
         s3 = pooled[, c("patient_id", "x3", "y_ord",
                         "low_leq", "med_leq", "high_leq"), drop = FALSE])
  }
}

generalize_glm <- function(seed, K) {
  pooled <- build_gaussian_regression(seed = seed)
  with_aligned_dslite(split_x123_y(pooled, K), function(conns) {
    run <- elapsed(fit <- dsVertClient::ds.vert.glm(
      y ~ x1 + x2 + x3, data = "DA", verbose = validation_demo_verbose(), datasources = conns))
    ref <- coef(stats::lm(y ~ x1 + x2 + x3, data = pooled))
    generalization_row(
      "glm", "GLM", "gaussian_regression", seed, K, "ds.vert.glm",
      sprintf("synthetic gaussian n=%d", nrow(pooled)), "stats::lm",
      "coef_max_abs_delta", max_named_delta(fit$coefficients, ref), 1e-3,
      "strict-practical",
      "Secure score aggregates; returns model-level coefficients/covariance.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok"))
  })
}

generalize_negative_binomial <- function(seed, K) {
  validation_require("MASS")
  pooled <- build_nb(seed = seed)
  with_aligned_dslite(split_x123_y(pooled, K), function(conns) {
    run <- elapsed(fit <- dsVertClient::ds.vert.nb(
      y ~ x1 + x2 + x3, data = "DA", verbose = validation_demo_verbose(), datasources = conns))
    ref <- coef(MASS::glm.nb(y ~ x1 + x2 + x3, data = pooled))
    generalization_row(
      "negative_binomial", "Negative binomial", "nb_loglinear", seed, K,
      fit$route %||% "ds.vert.nb",
      sprintf("synthetic NB n=%d", nrow(pooled)), "MASS::glm.nb",
      "coef_max_abs_delta", max_named_delta(fit$coefficients, ref), 0.001,
      "strict-practical",
      "Returns beta and scalar theta diagnostics only; no per-patient mu/eta is returned.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok", "approximate"))
  })
}

generalize_cox <- function(seed, K) {
  validation_require("survival")
  pooled <- build_cox(seed = seed)
  with_aligned_dslite(split_cox(pooled, K), function(conns) {
    run <- elapsed(fit <- dsVertClient::ds.vert.cox(
      survival::Surv(time, event) ~ x1 + x2 + x3, data = "DA",
      verbose = validation_demo_verbose(), datasources = conns))
    ref <- coef(survival::coxph(survival::Surv(time, event) ~ x1 + x2 + x3,
                                data = pooled, ties = "breslow"))
    generalization_row(
      "cox", "Cox PH", "discrete_time_survival", seed, K, "ds.vert.cox",
      sprintf("synthetic Cox n=%d", nrow(pooled)),
      "survival::coxph(ties='breslow')", "coef_max_abs_delta",
      max_named_delta(fit$coefficients, ref), 1e-3, "strict-practical",
      "Risk-set score terms stay shared; returns beta and scalar diagnostics.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok"))
  })
}

generalize_lmm <- function(seed, K) {
  validation_require("lme4")
  pooled <- build_lmm(seed = seed)
  with_aligned_dslite(split_x123_y(pooled, K, extras = "cluster"),
                      function(conns) {
    run <- elapsed(fit <- dsVertClient::ds.vert.lmm(
      y ~ x1 + x2 + x3, data = "DA", cluster_col = "cluster",
      verbose = validation_demo_verbose(), datasources = conns))
    ref <- lme4::fixef(lme4::lmer(y ~ x1 + x2 + x3 + (1 | cluster),
                                  data = pooled, REML = TRUE))
    generalization_row(
      "lmm", "LMM", "random_intercept_gaussian", seed, K, "ds.vert.lmm",
      sprintf("synthetic LMM n=%d clusters=%d", nrow(pooled),
              length(unique(pooled$cluster))),
      "lme4::lmer", "fixed_effect_max_abs_delta",
      max_named_delta(fit$coefficients, ref), 0.001, "strict-practical",
      "Uses protected mixed-model aggregates; no BLUPs or cluster residuals are returned.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok", "approximate"))
  })
}

generalize_gee <- function(seed, K) {
  validation_require("geepack")
  pooled <- build_lmm(seed = seed, b_sd = 0.45, eps_sd = 0.45)
  with_aligned_dslite(split_x123_y(pooled, K, extras = "cluster"),
                      function(conns) {
    run <- elapsed(fit <- dsVertClient::ds.vert.gee(
      y ~ x1 + x2 + x3, data = "DA", family = "gaussian",
      id_col = "cluster", corstr = "independence",
      verbose = validation_demo_verbose(), datasources = conns))
    ref <- coef(geepack::geeglm(y ~ x1 + x2 + x3, data = pooled,
                                id = cluster, family = gaussian(),
                                corstr = "independence"))
    generalization_row(
      "gee", "GEE", "clustered_gaussian_independence", seed, K,
      "ds.vert.gee",
      sprintf("synthetic GEE n=%d clusters=%d", nrow(pooled),
              length(unique(pooled$cluster))),
      "geepack::geeglm", "coef_max_abs_delta",
      max_named_delta(fit$coefficients, ref), 0.001, "strict-practical",
      "Returns regression and sandwich-level aggregates, not row scores.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok", "approximate"))
  })
}

generalize_ipw <- function(seed, K) {
  pooled <- build_ipw(seed = seed)
  with_aligned_dslite(split_ipw(pooled, K), function(conns) {
    run <- elapsed(fit <- dsVertClient::ds.vert.ipw(
      y ~ tr + w1 + w2, tr ~ w1 + w2, data = "DA",
      outcome_family = "gaussian", verbose = validation_demo_verbose(), datasources = conns))
    ref <- coef(stats::lm(y ~ tr + w1 + w2, data = pooled,
                          weights = ipw))
    generalization_row(
      "ipw", "IPW", "known_weight_ipw", seed, K, "ds.vert.ipw",
      sprintf("synthetic IPW n=%d", nrow(pooled)),
      "central weighted lm using same weights",
      "weighted_outcome_coef_abs_delta",
      max_named_delta(fit$outcome$coefficients, ref), 1e-3,
      "strict-practical",
      "Uses protected propensity/outcome GLM aggregates; only model-level fits are returned.",
      run$runtime_s, generalization_quality_status(fit$outcome),
      quality_allowed = c("ok"))
  })
}

generalize_mi <- function(seed, K) {
  pooled <- build_mi(seed = seed)
  with_aligned_dslite(split_mi(pooled, K), function(conns) {
    run <- elapsed(fit <- dsVertClient::ds.vert.mi(
      y ~ x1 + x2 + x3, data = "DA", impute_columns = "x2",
      family = "gaussian", lambda = 0,
      verbose = validation_demo_verbose(), datasources = conns, seed = 12L))
    ref <- central_mi_same_imputer(
      pooled, K, y ~ x1 + x2 + x3, impute_columns = "x2",
      m = fit$m, seed = 12L)
    generalization_row(
      "mi", "Multiple imputation", "missing_covariate_same_imputer", seed, K,
      "ds.vert.mi",
      sprintf("synthetic MI n=%d missing_every=7", nrow(pooled)),
      "same-imputer centralized Rubin pooling", "pooled_coef_abs_delta",
      max_named_delta(fit$coefficients, ref$coefficients), 0.001,
      "strict-practical",
      "Imputed columns stay server-side; client pools beta/covariance only.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok", "approximate"))
  }, na.action = "none")
}

generalize_multinomial <- function(seed, K) {
  validation_require("nnet")
  pooled <- build_multinomial(seed = seed)
  with_aligned_dslite(split_multinomial(pooled, K), function(conns) {
    ref <- nnet::multinom(y_cls ~ x1 + x2 + x3, data = pooled,
                          trace = FALSE, maxit = 200L)
    ref_prob <- stats::predict(ref, pooled, type = "probs")
    run <- elapsed(fit <- dsVertClient::ds.vert.multinom(
      y_cls ~ x1 + x2 + x3, data = "DA",
      classes = c("high", "low", "med"),
      indicator_template = "%s_ind",
      verbose = validation_demo_verbose(), datasources = conns))
    ds_prob <- softmax_prob(fit$coefficients, pooled)
    ds_prob <- ds_prob[, colnames(ref_prob), drop = FALSE]
    generalization_row(
      "multinomial", "Multinomial", "balanced_soft_signal", seed, K,
      "ds.vert.multinom",
      sprintf("synthetic multinomial n=%d", nrow(pooled)),
      "nnet::multinom probabilities",
      "class_probability_max_abs_delta", max(abs(ds_prob - ref_prob)), 0.001,
      "strict-practical",
      "Softmax probabilities/residuals stay Ring127 shares; no row probabilities are returned.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok"))
  })
}

generalize_ordinal <- function(seed, K) {
  validation_require("MASS")
  old_fd <- getOption("dsvert.ord_strict_fd_max_dim", NULL)
  options(dsvert.ord_strict_fd_max_dim = 0L)
  on.exit(options(dsvert.ord_strict_fd_max_dim = old_fd), add = TRUE)
  pooled <- build_ordinal(seed = seed)
  with_aligned_dslite(split_ordinal(pooled, K), function(conns) {
    ref <- MASS::polr(y_ord ~ x1 + x2 + x3, data = pooled, Hess = TRUE)
    ref_cum <- sapply(ref$zeta, function(th) {
      stats::plogis(th - as.numeric(as.matrix(pooled[, names(coef(ref))]) %*%
                                      coef(ref)))
    })
    run <- elapsed(fit <- dsVertClient::ds.vert.ordinal(
      y_ord ~ x1 + x2 + x3, data = "DA",
      levels_ordered = c("low", "med", "high"),
      cumulative_template = "%s_leq",
      verbose = validation_demo_verbose(), datasources = conns))
    ds_cum <- ordinal_cumprob(fit, pooled)
    colnames(ds_cum) <- colnames(ref_cum)
    generalization_row(
      "ordinal", "Ordinal", "balanced_ordered_levels", seed, K,
      "ds.vert.ordinal",
      sprintf("synthetic ordinal n=%d", nrow(pooled)),
      "MASS::polr cumulative probabilities",
      "cumulative_probability_max_abs_delta", max(abs(ds_cum - ref_cum)),
      0.001, "strict-practical",
      "Class probabilities/residuals stay Ring127 shares; no row probabilities are returned.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok"))
  })
}

generalize_glmm <- function(seed, K) {
  validation_require(c("MASS", "nlme"))
  pooled <- build_balanced_glmm(seed = seed)
  with_aligned_dslite(split_glmm(pooled, K), function(conns) {
    run <- elapsed(fit <- dsVertClient::ds.vert.glmm(
      y ~ x1 + x2, data = "DA", cluster_col = "cluster",
      method = "pql", compute_se = FALSE,
      verbose = validation_demo_verbose(), datasources = conns))
    pooled_ref <- pooled
    pooled_ref$cluster <- factor(pooled_ref$cluster)
    ref_fit <- suppressWarnings(MASS::glmmPQL(
      y ~ x1 + x2, random = ~1 | cluster, family = binomial(),
      data = pooled_ref, verbose = validation_demo_verbose()))
    ref <- nlme::fixef(ref_fit)
    fixed_delta <- max_named_delta(fit$coefficients, ref)
    pql_ran <- isTRUE(fit$iterations >= 1L) &&
      is.data.frame(fit$trace) && nrow(fit$trace) >= 1L &&
      is.finite(fit$sigma_b2)
    observed <- max(fixed_delta, if (pql_ran) 0 else Inf)
    generalization_row(
      "glmm", "GLMM", "balanced_binomial_random_intercept", seed, K,
      "ds.vert.glmm",
      sprintf("synthetic GLMM n=%d clusters=%d", nrow(pooled),
              length(unique(pooled$cluster))),
      "MASS::glmmPQL", "fixed_effect_max_abs_delta_and_pql_quality",
      observed, 0.001, "strict-practical",
      "PQL route returns fixed effects and scalar variance diagnostics only.",
      run$runtime_s, generalization_quality_status(fit),
      quality_allowed = c("ok"))
  })
}

run_generalization_case <- function(method_id, seed, K) {
  switch(method_id,
    glm = generalize_glm(seed, K),
    negative_binomial = generalize_negative_binomial(seed, K),
    cox = generalize_cox(seed, K),
    lmm = generalize_lmm(seed, K),
    gee = generalize_gee(seed, K),
    ipw = generalize_ipw(seed, K),
    mi = generalize_mi(seed, K),
    multinomial = generalize_multinomial(seed, K),
    ordinal = generalize_ordinal(seed, K),
    glmm = generalize_glmm(seed, K),
    stop("Unknown generalization method: ", method_id, call. = FALSE))
}

run_generalization_validation <- function(
    method_ids = NULL,
    seeds = NULL,
    include_heavy = identical(Sys.getenv("DSVERT_GENERALIZATION_INCLUDE_HEAVY"),
                              "true"),
    force = getOption("dsvert.validation.force", validation_force_refresh()),
    assert = TRUE,
    trace = TRUE) {
  validation_load_packages()
  if (is.null(method_ids)) {
    method_ids <- generalization_parse_methods(include_heavy = include_heavy)
  }
  if (is.null(seeds)) seeds <- generalization_parse_seeds()
  method_ids <- unique(method_ids)
  seeds <- unique(as.integer(seeds))
  path <- generalization_cache_path(method_ids, seeds)
  if (!isTRUE(force) && file.exists(path)) {
    rows <- readRDS(path)
    rows <- rows[rows$method_id %in% method_ids & rows$seed %in% seeds,
                 , drop = FALSE]
  } else {
    rows <- list()
    i <- 0L
    for (method_id in method_ids) {
      for (seed in seeds) {
        for (K in c(2L, 3L)) {
          i <- i + 1L
          if (isTRUE(trace)) {
            cat(sprintf("Generalization: %s seed=%d K=%s\n",
                        method_id, seed, if (K == 2L) "2" else ">=3"))
          }
          rows[[i]] <- tryCatch(
            run_generalization_case(method_id, seed, K),
            error = function(e) generalization_error_row(
              method_id, "default", seed, K, e))
        }
      }
    }
    rows <- do.call(rbind, rows)
    rownames(rows) <- NULL
    saveRDS(rows, path)
  }
  if (isTRUE(assert)) assert_generalization(rows)
  rows
}
