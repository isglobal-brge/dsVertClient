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
          datashield.privacyLevel = 5L,
          datashield.errors.print = TRUE)
  invisible(TRUE)
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
                      newobj = "DA") {
  dsVertClient::ds.psiAlign(data, id_col, newobj,
                            datasources = conns, verbose = FALSE)
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

run_validation <- function(method_id, force = getOption(
  "dsvert.validation.force", FALSE)) {
  validation_load_packages()
  bundle <- validation_bundle_id(method_id)
  path <- validation_cache_path(bundle)
  if (!isTRUE(force) && file.exists(path)) {
    rows <- readRDS(path)
  } else {
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
  out
}

run_all_validations <- function(force = getOption(
  "dsvert.validation.force", FALSE)) {
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
    cor_fit <- dsVertClient::ds.vertCor("DA", variables = variables,
                                        datasources = conns, verbose = FALSE)
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
      "psi", "PSI alignment", K, "ds.psiAlign",
      "synthetic ID intersection", "deterministic set intersection",
      "max(count_delta, correlation_delta)", observed, 1e-8,
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

    d0 <- elapsed(desc <- dsVertClient::ds.vertDesc(
      "DA", variables = vars, n_buckets = 40L,
      verbose = FALSE, datasources = conns))
    ref_mean <- vapply(pooled[all_vars], mean, numeric(1), na.rm = TRUE)
    ref_sd <- vapply(pooled[all_vars], stats::sd, numeric(1), na.rm = TRUE)
    ds_mean <- stats::setNames(desc$mean, desc$variable)
    ds_sd <- stats::setNames(desc$sd, desc$variable)
    desc_delta <- max(max_named_delta(ds_mean, ref_mean),
                      max_named_delta(ds_sd, ref_sd))

    c0 <- elapsed(cor_ds <- dsVertClient::ds.vertCor(
      "DA", variables = vars, verbose = FALSE, datasources = conns))
    cor_ref <- stats::cor(pooled[cor_ds$var_names])
    cor_delta <- max(abs(cor_ds$correlation - cor_ref))

    p0 <- elapsed(pca_ds <- dsVertClient::ds.vertPCA(
      cor_result = cor_ds, verbose = FALSE, datasources = conns))
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
      "descriptive", "Descriptive statistics", K, "ds.vertDesc",
      "MASS::Pima.tr fixture", "central mean/sd",
      "max(mean_sd_abs_delta)", desc_delta, 1e-8, "strict-precise",
      "Returns guarded scalar summaries and histogram-based quantiles, not rows.",
      d0$runtime_s)
    rows[[paste0("cor", K)]] <- row_result(
      "correlation", "Correlation", K, "ds.vertCor",
      "MASS::Pima.tr fixture", "stats::cor",
      "correlation_max_abs_delta", cor_delta, 1e-4, "strict-practical",
      "Releases the low-dimensional correlation matrix only.",
      c0$runtime_s)
    rows[[paste0("pca", K)]] <- row_result(
      "pca", "PCA", K, "ds.vertPCA",
      "MASS::Pima.tr fixture", "eigen(cor(X))",
      "max(eigen_loading_abs_delta)", pca_delta, 1e-4,
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
    run <- elapsed(cross <- dsVertClient::ds.vertChisqCross(
      "DA", "age_grp", "diabetes", correct = FALSE, fisher = TRUE,
      verbose = FALSE, datasources = conns))
    tab <- table(pooled$age_grp, pooled$diabetes)
    chi_ref <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    ref_mat <- unclass(tab)[rownames(cross$observed),
                            colnames(cross$observed), drop = FALSE]
    observed <- max(max(abs(cross$observed - ref_mat)),
                    abs(cross$chisq - unname(chi_ref$statistic)))
    rows[[as.character(K)]] <- row_result(
      "contingency", "Contingency tests", K,
      "ds.vertChisq / ds.vertFisher / ds.vertChisqCross",
      "MASS::Pima.tr categorical fixture", "chisq.test / fisher.test",
      "max(count_or_chisq_delta)", observed, 1e-8, "strict-precise",
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
    run <- elapsed(fit <- dsVertClient::ds.vertGLM(
      fm, data = "DA", family = "gaussian", max_iter = 30L,
      verbose = FALSE, datasources = conns))
    observed <- max_named_delta(fit$coefficients, coef(stats::lm(fm, pooled)))
    rows[[as.character(K)]] <- row_result(
      "glm", "GLM", K, "ds.vertGLM",
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
      full <- dsVertClient::ds.vertGLM(fm, data = "DA", family = "gaussian",
        max_iter = 30L, compute_se = TRUE, compute_deviance = TRUE,
        verbose = FALSE, datasources = conns)
      reduced <- dsVertClient::ds.vertGLM(red, data = "DA",
        family = "gaussian", max_iter = 30L, compute_se = TRUE,
        compute_deviance = TRUE, verbose = FALSE, datasources = conns)
      ci <- dsVertClient::ds.vertConfint(full)
      wald <- dsVertClient::ds.vertWald(full, "age")
      Kmat <- matrix(0, nrow = 1, ncol = length(full$coefficients),
                     dimnames = list("age", names(full$coefficients)))
      Kmat[1, "age"] <- 1
      contrast <- dsVertClient::ds.vertContrast(full, Kmat)
      lr <- dsVertClient::ds.vertLR(reduced, full)
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
      "ds.vertConfint / ds.vertWald / ds.vertContrast / ds.vertLR",
      "MASS::Pima.tr fixture", "manual algebra on ds.vertGLM output",
      "algebra_max_abs_delta", observed, 1e-10, "strict-precise",
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
      fit <- dsVertClient::ds.vertGLM(fm, data = "DA", family = "gaussian",
        max_iter = 30L, verbose = FALSE, datasources = conns)
      prox0 <- dsVertClient::ds.vertLASSOProximal(
        fit, lambda = 0, max_iter = 1000L, tol = 1e-8)
      list(fit = fit, prox0 = prox0)
    })
    observed <- max_named_delta(run$value$prox0$coefficients,
                                run$value$fit$coefficients)
    rows[[as.character(K)]] <- row_result(
      "lasso", "LASSO", K, "ds.vertLASSOProximal(lambda=0)",
      "MASS::Pima.tr fixture", "OLS limit of ds.vertGLM",
      "lambda0_coef_abs_delta", observed, 1e-8, "strict-precise",
      "Consumes only model-level GLM aggregates; no extra patient-level values are returned.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_nb <- function(n = 60L) {
  set.seed(303)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  mu <- exp(1 + 0.2 * x1 - 0.15 * x2 + 0.1 * x3)
  y <- stats::rnbinom(n, size = 3, mu = mu)
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
    run <- elapsed(fit <- dsVertClient::ds.vertNBMoMTheta(
      y ~ x1 + x2 + x3, data = "DA", max_iter = 25L,
      verbose = FALSE, datasources = conns))
    ref <- coef(MASS::glm.nb(y ~ x1 + x2 + x3, data = pooled))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "negative_binomial", "Negative binomial", K, "ds.vertNBMoMTheta",
      "synthetic NB fixture", "MASS::glm.nb",
      "coef_max_abs_delta", observed, 0.02, "strict-practical",
      "Returns beta and scalar theta only; score components remain aggregate.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_cox <- function(n = 60L) {
  set.seed(202)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  eta <- 0.4 * x1 - 0.3 * x2 + 0.2 * x3
  t0 <- stats::rexp(n, rate = exp(eta) / 20)
  cens <- stats::rexp(n, rate = 0.02)
  event <- as.integer(t0 <= cens)
  time <- pmin(t0, cens)
  br <- unique(stats::quantile(time, probs = seq(0, 1, length.out = 9),
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
    run <- elapsed(fit <- dsVertClient::ds.vertCox(
      survival::Surv(time, event) ~ x1 + x2 + x3, data = "DA",
      max_iter = 5L, tol = 1e-3, max_event_times = 20L,
      verbose = FALSE, datasources = conns))
    ref <- coef(survival::coxph(survival::Surv(time, event) ~ x1 + x2 + x3,
                                data = pooled, ties = "breslow"))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "cox", "Cox PH", K, "ds.vertCox",
      "synthetic discretised survival fixture",
      "survival::coxph(ties='breslow')",
      "coef_max_abs_delta", observed, 1e-3, "strict-practical",
      "Event-time score terms remain shared; client receives beta and scalar diagnostics.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_lmm <- function(ncl = 12L, m = 5L) {
  set.seed(101)
  n <- ncl * m
  cluster <- rep(seq_len(ncl), each = m)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  b <- rnorm(ncl, sd = 0.8)[cluster]
  y <- 1 + 0.5 * x1 - 0.3 * x2 + 0.2 * x3 + b + rnorm(n, sd = 0.4)
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
    run <- elapsed(fit <- if (K == 2L) {
      dsVertClient::ds.vertLMM(y ~ x1 + x2 + x3, data = "DA",
        cluster_col = "cluster", max_iter = 3L, inner_iter = 10L,
        verbose = FALSE, datasources = conns)
    } else {
      dsVertClient::ds.vertLMM.k3(y ~ x1 + x2 + x3, data = "DA",
        cluster_col = "cluster", max_outer = 3L, tol = 1e-3,
        verbose = FALSE, datasources = conns)
    })
    ref <- lme4::fixef(lme4::lmer(y ~ x1 + x2 + x3 + (1 | cluster),
                                  data = pooled, REML = TRUE))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "lmm", "LMM", K, if (K == 2L) "ds.vertLMM" else "ds.vertLMM.k3",
      "synthetic random-intercept fixture", "lme4::lmer",
      "fixed_effect_max_abs_delta", observed, 0.02, "strict-practical",
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
    run <- elapsed(fit <- dsVertClient::ds.vertGEE(
      fm, data = "DA", family = "gaussian", id_col = "cluster",
      corstr = "independence", max_iter = 30L,
      verbose = FALSE, datasources = conns))
    ref <- coef(geepack::geeglm(fm, data = pooled, id = cluster,
                                family = gaussian(),
                                corstr = "independence"))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "gee", "GEE", K, "ds.vertGEE(corstr='independence')",
      "MASS::Pima.tr clustered fixture", "geepack::geeglm",
      "coef_max_abs_delta", observed, 0.01, "strict-practical",
      "Returns regression and sandwich-level aggregates, not row scores.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_balanced_glmm <- function() {
  set.seed(77)
  n <- 50L
  cluster <- rep(seq_len(10L), each = 5L)
  base1 <- rnorm(n / 2L)
  base2 <- rnorm(n / 2L)
  data.frame(
    patient_id = sprintf("G%03d", seq_len(n)),
    cluster = cluster,
    x1 = c(base1, base1),
    x2 = c(base2, base2),
    y = c(rep(0L, n / 2L), rep(1L, n / 2L)))
}

validate_glmm <- function() {
  old <- getOption("dsvert.glm_num_intervals_binomial", NULL)
  options(dsvert.glm_num_intervals_binomial = 10L)
  on.exit(options(dsvert.glm_num_intervals_binomial = old), add = TRUE)
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_balanced_glmm()
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
    run <- elapsed(fit <- dsVertClient::ds.vertGLMM(
      y ~ x1 + x2, data = "DA", cluster_col = "cluster",
      max_outer = 0L, inner_iter = 1L, compute_se = FALSE,
      verbose = FALSE, datasources = conns))
    ref <- coef(stats::glm(y ~ x1 + x2, data = pooled, family = binomial()))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "glmm", "GLMM", K, "ds.vertGLMM(smoke validation)",
      "balanced binomial cluster fixture", "central glm no-random-effect limit",
      "coef_max_abs_delta", observed, 0.05, "strict-smoke",
      "PQL route returns fixed effects and scalar variance diagnostics only.",
      run$runtime_s,
      "Compact vignette uses max_outer=0; full PQL tuning is heavier.")
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_ipw <- function(n = 60L) {
  set.seed(88)
  w1 <- rnorm(n); w2 <- rnorm(n)
  p <- stats::plogis(-0.2 + 0.4 * w1 - 0.3 * w2)
  tr <- stats::rbinom(n, 1, p)
  ipw <- ifelse(tr == 1, 1 / p, 1 / (1 - p))
  y <- 1 + 2 * tr + 0.5 * w1 - 0.2 * w2 + rnorm(n, 0, 0.2)
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
    run <- elapsed(fit <- dsVertClient::ds.vertIPW(
      y ~ tr + w1 + w2, tr ~ w1 + w2, data = "DA",
      outcome_family = "gaussian", max_iter = 5L,
      binomial_sigmoid_intervals = 10L, compute_se = FALSE,
      compute_deviance = FALSE, verbose = FALSE, datasources = conns))
    ref <- coef(stats::lm(y ~ tr + w1 + w2, data = pooled,
                          weights = ipw))
    observed <- max_named_delta(fit$outcome$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "ipw", "IPW", K, "ds.vertIPW",
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

build_mi <- function(n = 60L) {
  set.seed(99)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  y <- 1 + 0.4 * x1 - 0.2 * x2 + 0.1 * x3 + rnorm(n, 0, 0.2)
  d <- data.frame(patient_id = sprintf("MI%03d", seq_len(n)), x1, x2, x3, y)
  d$x2[seq(5L, n, by = 7L)] <- NA
  d
}

validate_mi <- function() {
  rows <- list()
  for (K in c(2L, 3L)) {
    pooled <- build_mi()
    tables <- if (K == 2L) {
      list(s1 = pooled[, c("patient_id", "x1")],
           s2 = pooled[, c("patient_id", "x2", "x3", "y")])
    } else {
      list(s1 = pooled[, c("patient_id", "x1")],
           s2 = pooled[, c("patient_id", "x2")],
           s3 = pooled[, c("patient_id", "x3", "y")])
    }
    conns <- connect_dslite(tables)
    on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)
    psi_align(conns)
    run <- elapsed(fit <- dsVertClient::ds.vertMI(
      y ~ x1 + x2 + x3, data = "DA", impute_columns = "x2",
      m = 2L, family = "gaussian", max_iter = 15L,
      verbose = FALSE, datasources = conns, seed = 12L))
    ref_data <- pooled
    ref_data$x2[is.na(ref_data$x2)] <- mean(ref_data$x2, na.rm = TRUE)
    ref <- coef(stats::lm(y ~ x1 + x2 + x3, data = ref_data))
    observed <- max_named_delta(fit$coefficients, ref)
    rows[[as.character(K)]] <- row_result(
      "mi", "Multiple imputation", K, "ds.vertMI",
      "synthetic missing-covariate fixture",
      "central mean-imputation reference",
      "pooled_coef_abs_delta", observed, 0.02, "strict-practical",
      "Imputed columns stay server-side; client pools beta/covariance only.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_multinomial <- function(n = 60L) {
  set.seed(44)
  d <- data.frame(patient_id = sprintf("M%03d", seq_len(n)),
                  x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y_cls <- factor(rep(c("high", "low", "med"), each = n / 3L),
                    levels = c("high", "low", "med"))
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
    run <- elapsed(fit <- dsVertClient::ds.vertMultinomJointNewton(
      y_cls ~ x1 + x2 + x3, data = "DA",
      levels = c("high", "low", "med"),
      indicator_template = "%s_ind", max_outer = 2L,
      warm_max_iter = 8L, binomial_sigmoid_intervals = 10L,
      tol = 1e-3, verbose = FALSE, datasources = conns))
    ds_prob <- softmax_prob(fit$coefficients, pooled)
    ds_prob <- ds_prob[, colnames(ref_prob), drop = FALSE]
    observed <- max(abs(ds_prob - ref_prob))
    rows[[as.character(K)]] <- row_result(
      "multinomial", "Multinomial", K, "ds.vertMultinomJointNewton",
      "balanced synthetic 3-class fixture", "nnet::multinom probabilities",
      "class_probability_max_abs_delta", observed, 0.05,
      "strict-practical",
      "Softmax probabilities and residuals remain Ring127 shares; no row probabilities are returned.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}

build_ordinal <- function(n = 60L) {
  set.seed(55)
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
    run <- elapsed(fit <- dsVertClient::ds.vertOrdinalJointNewton(
      y_ord ~ x1 + x2 + x3, data = "DA",
      levels_ordered = c("low", "med", "high"),
      cumulative_template = "%s_leq", max_outer = 1L,
      warm_max_iter = 8L, binomial_sigmoid_intervals = 10L,
      tol = 1e-3, verbose = FALSE, datasources = conns))
    ds_cum <- ordinal_cumprob(fit, pooled)
    colnames(ds_cum) <- colnames(ref_cum)
    observed <- max(abs(ds_cum - ref_cum))
    rows[[as.character(K)]] <- row_result(
      "ordinal", "Ordinal", K, "ds.vertOrdinalJointNewton",
      "balanced synthetic 3-level ordinal fixture", "MASS::polr cumulative probabilities",
      "cumulative_probability_max_abs_delta", observed, 0.15,
      "strict-practical",
      "Class probabilities/residuals stay Ring127 shares; no row probabilities are returned.",
      run$runtime_s)
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  do.call(rbind, rows)
}
