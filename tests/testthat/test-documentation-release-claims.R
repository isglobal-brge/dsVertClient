quarantined_doc_sources <- c(
  "ds.vertMultinomJoint.R",
  "ds.vertMultinomJointNewton.R",
  "ds.vertOrdinalJointNewton.R",
  "ds.vertLMM.R",
  "ds.vertLMM.k3.R",
  "ds.vertGLMM.R",
  "ds.vertIPW.R",
  "ds.vertMI.R")

quarantined_doc_topics <- c(
  "ds.vertMultinomJoint" = "ds.vertMultinomJoint.R",
  "ds.vertMultinomJointNewton" = "ds.vertMultinomJointNewton.R",
  "ds.vertOrdinalJointNewton" = "ds.vertOrdinalJointNewton.R",
  "ds.vertLMM" = "ds.vertLMM.R",
  "ds.vertLMM.k3" = "ds.vertLMM.k3.R",
  "ds.vertGLMM" = "ds.vertGLMM.R",
  "ds.vertIPW" = "ds.vertIPW.R",
  "ds.vertMI" = "ds.vertMI.R")

.dsvert_public_roxygen_text <- function(filename, function_name) {
  lines <- readLines(.dsvert_client_source_file(filename), warn = FALSE)
  index <- which(startsWith(lines, paste0(function_name, " <- function")))
  if (!length(index)) {
    testthat::fail(paste("public function not found:", function_name))
  }
  end <- index[[1L]] - 1L
  start <- end
  while (start >= 1L && startsWith(lines[[start]], "#'")) start <- start - 1L
  block <- lines[seq.int(start + 1L, end)]
  paste(sub("^#'[[:space:]]?", "", block), collapse = "\n")
}

test_that("quarantined frontdoors describe their release status candidly", {
  for (filename in quarantined_doc_sources) {
    text <- paste(readLines(.dsvert_client_source_file(filename), warn = FALSE),
                  collapse = "\n")
    expect_match(
      text, "quarantin", ignore.case = TRUE,
      info = paste(filename, "must identify the route as quarantined"))
  }
})

test_that("generated help preserves quarantine warnings", {
  package_root <- dirname(.dsvert_client_source_root())
  topics <- names(quarantined_doc_topics)
  paths <- file.path(package_root, "man", paste0(topics, ".Rd"))
  if (!all(file.exists(paths))) {
    skip("generated Rd source is unavailable for the release-claim audit")
  }
  for (i in seq_along(paths)) {
    text <- paste(readLines(paths[[i]], warn = FALSE), collapse = "\n")
    expect_match(
      text, "quarantin", ignore.case = TRUE,
      info = paste(topics[[i]], "generated help lost its quarantine warning"))
  }
})

test_that("quarantined Roxygen and Rd describe only the zero-DSI frontdoor", {
  package_root <- dirname(.dsvert_client_source_root())
  forbidden <- paste(c(
    "No observation-level data is (ever )?disclosed",
    "[Tt]he client (receives|sees) only",
    "[Pp]er-patient quantities never",
    "[Pp]er-patient data never"), collapse = "|")

  for (topic in names(quarantined_doc_topics)) {
    source <- .dsvert_public_roxygen_text(
      quarantined_doc_topics[[topic]], topic)
    rd_path <- file.path(package_root, "man", paste0(topic, ".Rd"))
    if (!file.exists(rd_path)) {
      testthat::fail(paste("generated help is unavailable:", topic))
    }
    rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

    for (text in list(source, rd)) {
      expect_match(text, "dsvert_route_unavailable", fixed = TRUE,
                   info = paste(topic, "must name its typed condition"))
      expect_match(text, "before any DSI\\s+call", perl = TRUE,
                   info = paste(topic, "must promise the local zero-DSI gate"))
      expect_match(text, "No fitted object", fixed = TRUE,
                   info = paste(topic, "must not document a legacy return"))
      expect_false(grepl(forbidden, text, ignore.case = TRUE, perl = TRUE),
                   info = paste(topic, "contains an active legacy claim"))
    }
  }
})

test_that("GLM help documents the Synopsis route matrix instead of legacy MPC", {
  package_root <- dirname(.dsvert_client_source_root())
  source <- .dsvert_public_roxygen_text("ds.vertGLM.R", "ds.vertGLM")
  rd_path <- file.path(package_root, "man", "ds.vertGLM.Rd")
  if (!file.exists(rd_path)) {
    testthat::fail("generated ds.vertGLM help is unavailable")
  }
  rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

  forbidden <- paste(c(
    "Fits a GLM across vertically partitioned data",
    "No observation-level data is ever disclosed",
    "Only p-dimensional aggregate gradients are revealed"), collapse = "|")
  for (text in list(source, rd)) {
    expect_match(text, "two available routes",
                 ignore.case = TRUE, perl = TRUE)
    expect_match(text, "dp_analysis_id", fixed = TRUE)
    expect_match(text, "ds.vertDPGaussian", fixed = TRUE)
    expect_match(text, "formal_analysis_id", fixed = TRUE)
    expect_match(text, "completed public", ignore.case = TRUE)
    expect_match(text, "coefficient-only", fixed = TRUE)
    expect_match(text,
                 "without[\\s\\S]*signed-analysis id[\\s\\S]*before any DSI call",
                 ignore.case = TRUE, perl = TRUE)
    expect_false(grepl(forbidden, text, ignore.case = TRUE, perl = TRUE))
  }
  glm_status <- ds.vertMethodStatus(c("ds.vertGLM", "ds.vert.glm"))
  expect_true(all(glm_status$status == "promoted"))
  expect_match(glm_status$principal_limitation[[1L]],
               "Binomial and Poisson", fixed = TRUE)
})

test_that("multinomial help limits Frequency post-processing to y ~ 1", {
  package_root <- dirname(.dsvert_client_source_root())
  source <- .dsvert_public_roxygen_text("ds.vertMultinom.R", "ds.vertMultinom")
  rd_path <- file.path(package_root, "man", "ds.vertMultinom.Rd")
  if (!file.exists(rd_path)) {
    testthat::fail("generated ds.vertMultinom help is unavailable")
  }
  rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

  for (text in list(source, rd)) {
    expect_match(text, "ds.vertDPFrequency", fixed = TRUE)
    expect_match(text, "y ~ 1", fixed = TRUE)
    expect_match(text, "Jeffreys", fixed = TRUE)
    expect_match(text, "never starts a new analysis", fixed = TRUE)
    expect_match(text, "joint softmax", ignore.case = TRUE)
    expect_match(text, "standard errors", ignore.case = TRUE)
  }
})

test_that("NB2 help limits Frequency post-processing to bounded y ~ 1 counts", {
  package_root <- dirname(.dsvert_client_source_root())
  source <- .dsvert_public_roxygen_text(
    "ds.vertNBFullRegTheta.R", "ds.vertNBFullRegTheta")
  rd_path <- file.path(package_root, "man", "ds.vertNBFullRegTheta.Rd")
  if (!file.exists(rd_path)) {
    testthat::fail("generated ds.vertNBFullRegTheta help is unavailable")
  }
  rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

  for (text in list(source, rd)) {
    expect_match(text, "ds.vertDPFrequency", fixed = TRUE)
    expect_match(text, "y ~ 1", fixed = TRUE)
    expect_match(text, "non-negative integer", fixed = TRUE)
    expect_match(text, "method-of-moments", fixed = TRUE)
    expect_match(text, "never starts a new analysis", fixed = TRUE)
    expect_match(text, "covariance", ignore.case = TRUE)
    expect_match(text, "standard errors", ignore.case = TRUE)
  }
})

test_that("ordinal help requires an explicit order and limits output to thresholds", {
  package_root <- dirname(.dsvert_client_source_root())
  source <- .dsvert_public_roxygen_text("ds.vertOrdinal.R", "ds.vertOrdinal")
  rd_path <- file.path(package_root, "man", "ds.vertOrdinal.Rd")
  if (!file.exists(rd_path)) {
    testthat::fail("generated ds.vertOrdinal help is unavailable")
  }
  rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

  for (text in list(source, rd)) {
    expect_match(text, "ds.vertDPFrequency", fixed = TRUE)
    expect_match(text, "y ~ 1", fixed = TRUE)
    expect_match(text, "complete clinical", ignore.case = TRUE)
    expect_match(text, "cumulative-logit", fixed = TRUE)
    expect_match(text, "never starts a new analysis", fixed = TRUE)
    expect_match(text, "covariates", ignore.case = TRUE)
    expect_match(text, "standard errors", ignore.case = TRUE)
  }
})

test_that("discrete-hazard help limits the completed formal release", {
  package_root <- dirname(.dsvert_client_source_root())
  source <- .dsvert_public_roxygen_text(
    "ds.vertCoxDiscreteNonDisclosive.R",
    "ds.vertCoxDiscreteNonDisclosive")
  rd_path <- file.path(package_root, "man",
                     "ds.vertCoxDiscreteNonDisclosive.Rd")
  if (!file.exists(rd_path)) {
    testthat::fail("generated discrete-hazard help is unavailable")
  }
  rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")
  for (text in list(source, rd)) {
    expect_match(text, "two-authority-signed", fixed = TRUE)
    expect_match(text, "fixed time grid", fixed = TRUE)
    expect_match(text, "distinct from Cox proportional hazards",
                 ignore.case = TRUE, perl = TRUE)
    expect_match(text, "never starts expansion", fixed = TRUE)
    expect_match(text, "Covariance, standard errors, p-values", fixed = TRUE)
  }
})

test_that("Describe help matches the cold Synopsis execution and replay route", {
  package_root <- dirname(.dsvert_client_source_root())
  source <- .dsvert_public_roxygen_text("ds.vertDPDescribe.R",
                                        "ds.vertDPDescribe")
  rd_path <- file.path(package_root, "man", "ds.vertDPDescribe.Rd")
  if (!file.exists(rd_path)) {
    testthat::fail("generated ds.vertDPDescribe help is unavailable")
  }
  rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

  for (text in list(source, rd)) {
    normalized <- gsub("[[:space:]]+", " ", text)
    expect_match(normalized, "cold exact-GC Synopsis execution performs",
                 fixed = TRUE)
    expect_match(normalized, "durable publication fast path", fixed = TRUE)
    expect_false(grepl("fails before Claim or Compile", normalized,
                       fixed = TRUE))
  }
})

test_that("security-status docs cannot promote sealed routes via top-level readiness", {
  source <- .dsvert_public_roxygen_text(
    "ds.vertSecurityStatus.R", "ds.vertSecurityStatus")
  package_root <- dirname(.dsvert_client_source_root())
  readme <- paste(readLines(file.path(package_root, "README.md"), warn = FALSE),
                  collapse = "\n")

  for (text in list(source, readme)) {
    expect_match(text, "security-profile schema v4", fixed = TRUE)
    expect_match(text, "top-level `ready`", fixed = TRUE)
    expect_match(text, "`formal_dp_claim_eligible`", fixed = TRUE)
    expect_match(
      text,
      paste0("(never\\s+promote[sd]?\\s+formal GLM or formal Cox|",
             "separately\\s+reports formal GLM and formal Cox as not ready)"),
      perl = TRUE)
  }
})

test_that("alias help cannot imply that a name makes a route available", {
  package_root <- dirname(.dsvert_client_source_root())
  source_lines <- readLines(
    .dsvert_client_source_file("ds.vert.aliases.R"), warn = FALSE)
  end <- which(!startsWith(source_lines, "#'"))[[1L]] - 1L
  source <- paste(sub("^#'[[:space:]]?", "", source_lines[seq_len(end)]),
                  collapse = "\n")
  rd_path <- file.path(package_root, "man", "ds.vert.aliases.Rd")
  if (!file.exists(rd_path)) {
    testthat::fail("generated alias help is unavailable")
  }
  rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

  for (text in list(source, rd)) {
    expect_match(text, "ds.vertMethodStatus", fixed = TRUE)
    expect_match(text, "dsvert_route_unavailable", fixed = TRUE)
    expect_match(text, "before any\\s+DSI call", perl = TRUE)
    expect_false(grepl(
      "dispatch from the number of active|compatibility backends",
      text, ignore.case = TRUE, perl = TRUE))
  }
})

test_that("rendered reference pages preserve the audited route contracts", {
  package_root <- dirname(.dsvert_client_source_root())
  reference_root <- file.path(package_root, "docs", "reference")
  if (!dir.exists(reference_root)) {
    skip("pkgdown reference output is absent from the package tarball")
  }

  for (topic in names(quarantined_doc_topics)) {
    paths <- file.path(reference_root, paste0(topic, c(".md", ".html")))
    expect_true(all(file.exists(paths)),
                info = paste(topic, "rendered help is incomplete"))
    for (path in paths) {
      text <- paste(readLines(path, warn = FALSE), collapse = "\n")
      expect_match(text, "dsvert_route_unavailable", fixed = TRUE)
      expect_match(text, "before any DSI\\s+call", perl = TRUE)
      expect_match(text, "no\\s+fitted object", ignore.case = TRUE,
                   perl = TRUE)
    }
  }

  glm_paths <- file.path(reference_root,
                         paste0("ds.vertGLM", c(".md", ".html")))
  alias_paths <- file.path(reference_root,
                           paste0("ds.vert.aliases", c(".md", ".html")))
  for (path in glm_paths) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(text, "(one|only) available analysis route",
                 ignore.case = TRUE, perl = TRUE)
    expect_match(text, "dp_analysis_id", fixed = TRUE)
    expect_match(text, "No binomial or Poisson fit", fixed = TRUE)
  }
  for (path in alias_paths) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(text, "ds.vertMethodStatus", fixed = TRUE)
    expect_match(text, "dsvert_route_unavailable", fixed = TRUE)
    expect_match(text, "before any\\s+DSI call", perl = TRUE)
  }
})

test_that("retired validation vignettes cannot become website evidence", {
  package_root <- dirname(.dsvert_client_source_root())
  pkgdown <- file.path(package_root, "_pkgdown.yml")
  archive <- file.path(package_root, "_archive_vignettes")
  if (!file.exists(pkgdown) || !dir.exists(archive)) {
    skip("repository-only website sources are absent from the package tarball")
  }

  config <- paste(readLines(pkgdown, warn = FALSE), collapse = "\n")
  expect_false(grepl("^[[:space:]]*articles:", config, perl = TRUE))
  expect_false(grepl("validation_[a-z0-9_]", config, perl = TRUE))

  rmd <- list.files(archive, pattern = "[.]Rmd$", full.names = TRUE)
  expect_gt(length(rmd), 0L)
  for (path in rmd) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(
      text, "Archived historical", fixed = TRUE,
      info = paste(basename(path), "must carry the historical-evidence banner"))
  }
})

test_that("Gaussian ownership claims match the implemented capsule registry", {
  package_root <- dirname(.dsvert_client_source_root())
  paths <- c(
    file.path(package_root, "README.md"),
    file.path(package_root, "inst", "docs", "product_surface.md"),
    file.path(package_root, "inst", "docs",
              "capsule_method_migration_matrix.md"))
  if (!all(file.exists(paths))) {
    skip("repository documentation is absent from the package tarball")
  }
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)),
                collapse = "\n")

  expect_match(text, "same- or cross-owner", fixed = TRUE)
  expect_false(grepl(
    "Cross-owner products and sampling inference are unavailable",
    text, fixed = TRUE))
  expect_false(grepl(
    "binomial/Poisson/cross-owner routes still require migration",
    text, fixed = TRUE))
})

test_that("README maturity and numeric claims match the runtime registry", {
  package_root <- dirname(.dsvert_client_source_root())
  paths <- c(
    file.path(package_root, "README.md"),
    file.path(package_root, "inst", "docs", "numeric_surface_inventory.md"))
  if (!all(file.exists(paths))) {
    skip("repository maturity documents are absent from the installed package")
  }
  readme <- paste(readLines(paths[[1L]], warn = FALSE), collapse = "\n")
  numeric <- paste(readLines(paths[[2L]], warn = FALSE), collapse = "\n")

  gee <- ds.vertMethodStatus(c("ds.vertGEE", "ds.vert.gee"))
  expect_true(all(gee$status == "promoted"))
  expect_match(
    readme,
    "independence-working GEE point estimate",
    fixed = TRUE)
  expect_match(
    numeric,
    "Read-only independence-working point adapter",
    fixed = TRUE)

  penalised <- ds.vertMethodStatus(c(
    "ds.vertLASSOIter", "ds.vertLASSO", "ds.vertLASSO1Step",
    "ds.vertLASSOProximal", "ds.vertLASSOCV"))
  expect_identical(
    stats::setNames(penalised$status, penalised$method)[c(
      "ds.vertLASSOIter", "ds.vertLASSO", "ds.vertLASSO1Step",
      "ds.vertLASSOProximal", "ds.vertLASSOCV")],
    c(
      ds.vertLASSOIter = "promoted",
      ds.vertLASSO = "promoted",
      ds.vertLASSO1Step = "promoted",
      ds.vertLASSOProximal = "promoted",
      ds.vertLASSOCV = "promoted"))
  expect_match(
    readme,
    paste0("`ds.vertLASSOIter()`, `ds.vertLASSO()` and ",
           "`ds.vertLASSO1Step()` are promoted Gaussian Synopsis paths ",
           "when an explicit signed `dp_analysis_id` selects the ",
           "same-owner artifact"),
    fixed = TRUE)

  registry <- ds.vertMethodStatus()
  expect_false(any(registry$may_report_numerically_certified))
  expect_false(any(registry$currently_numerically_certified))
  expect_match(
    numeric, "may_report_numerically_certified = FALSE", fixed = TRUE)
  expect_match(
    numeric, "currently_numerically_certified = FALSE", fixed = TRUE)
  expect_false(grepl(
    "permitted to report `numerically_certified`", numeric, fixed = TRUE))
})

test_that("formal GEE help states its narrow certified-reader contract", {
  package_root <- dirname(.dsvert_client_source_root())
  source <- .dsvert_public_roxygen_text("ds.vertGEE.R", "ds.vertGEE")
  rd_path <- file.path(package_root, "man", "ds.vertGEE.Rd")
  expect_true(file.exists(rd_path))
  rd <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

  for (text in list(source, rd)) {
    expect_match(text, "formal_analysis_id", fixed = TRUE)
    expect_match(text, "independence", fixed = TRUE)
    expect_match(text, "two-authority", fixed = TRUE)
    expect_match(text, "sandwich covariance", fixed = TRUE)
    expect_match(text, "failing locally before DSI", fixed = TRUE)
  }
})

test_that("the installed single-mode plan cannot describe a live dual profile", {
  package_root <- dirname(.dsvert_client_source_root())
  path <- file.path(
    package_root, "inst", "docs", "disclosure_safe_single_mode_plan.md")
  if (!file.exists(path)) {
    skip("single-mode implementation record is absent from the package tarball")
  }
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "historical implementation and acceptance plan",
               ignore.case = TRUE)
  expect_match(text, "one immutable `disclosure_safe` profile", fixed = TRUE)
  expect_match(text, "fail before DSI", fixed = TRUE)
  expect_false(grepl("current dual-profile development gate", tolower(text),
                     fixed = TRUE))
})

test_that("installed Cox and LMM notes retain their exact release boundaries", {
  package_root <- dirname(.dsvert_client_source_root())
  paths <- file.path(
    package_root, "inst", "docs", "disclosure_budget", c("cox.md", "lmm.md"))
  paths <- paths[file.exists(paths)]
  expect_gte(length(paths), 2L)

  forbidden <- "\\b(shipping|shipped|currently|current|proof|PASS(_PRACTICAL)?)\\b"
  for (path in paths) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_false(grepl(forbidden, text, ignore.case = TRUE, perl = TRUE),
                 info = paste(path, "contains a release-status overclaim"))
    if (identical(basename(path), "cox.md")) {
      expect_match(text, "Security-profile schema v4", fixed = TRUE)
      expect_match(text, "`route_claims$formal_cox_ready = FALSE`",
                   fixed = TRUE)
      expect_match(text, "read-only", ignore.case = TRUE)
      expect_match(text, "before any DSI\\s+call", perl = TRUE)
    } else {
      expect_match(text, "quarantined", ignore.case = TRUE,
                   info = paste(basename(path), "must state quarantine"))
      expect_match(text, "before any DSI call", fixed = TRUE,
                   info = paste(basename(path), "must state the zero-DSI gate"))
    }
  }
})

test_that("installed legacy Cox and LMM experiments are marked archival", {
  package_root <- dirname(.dsvert_client_source_root())
  roots <- c(package_root, file.path(dirname(package_root), "dsVert"))
  roots <- roots[dir.exists(roots)]
  paths <- unlist(lapply(roots, function(root) file.path(
    root, "inst", "docs", c(
      "acceptance/path_b_targets.md",
      "error_bounds/cox_newton_onestep.md"))), use.names = FALSE)
  installed <- vapply(c(
    "acceptance/path_b_targets.md",
    "error_bounds/cox_newton_onestep.md"), function(relative_path) {
      system.file("docs", relative_path, package = "dsVertClient")
    }, character(1L))
  paths <- unique(c(paths[file.exists(paths)],
                    installed[nzchar(installed) & file.exists(installed)]))
  expect_gte(length(paths), 2L)
  for (path in paths) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(text, "Archived historical design record", fixed = TRUE)
    expect_match(text, "fails? before any DSI call")
  }
})

test_that("installed and prominent historical scripts have no unsafe defaults", {
  package_root <- dirname(.dsvert_client_source_root())
  installed <- list.files(file.path(package_root, "inst", "scripts"),
                          pattern = "[.]R$", full.names = TRUE)

  repo_root <- normalizePath(file.path(package_root, ".."), mustWork = FALSE)
  prominent <- file.path(repo_root, c(
    "demo_dsvert.R",
    "register_methods.R",
    "demo/run_opal_demo.R",
    "demo/run_pima_demo.R",
    "scripts/opal_demo_provision_install.R",
    "scripts/opal_demo_provision_datasets.R",
    "scripts/opal_demo_provision_register.R",
    "scripts/opal_demo_validate_e2e.R",
    "scripts/run_opal_demo_probes.R",
    "scripts/publish_methods.R",
    "vignettes/methods/_opal_demo_setup.R",
    "vignettes/methods/_opal_demo_helpers.R"))
  prominent <- prominent[file.exists(prominent)]
  repository_launchers <- character(0)
  if (file.exists(file.path(repo_root, "demo_dsvert.R"))) {
    directories <- c(repo_root,
                     file.path(repo_root, c("demo", "scripts",
                                             "vignettes/methods")))
    directories <- directories[dir.exists(directories)]
    repository_launchers <- unlist(lapply(
      directories,
      function(directory) list.files(
        directory, pattern = "[.](R|sh)$", full.names = TRUE,
        recursive = FALSE)), use.names = FALSE)
  }
  paths <- unique(c(installed, prominent, repository_launchers))

  forbidden <- c(
    "admin123",
    "opal-demo[.]obiba[.]org",
    "user[[:space:]]*=[[:space:]]*['\"]administrator['\"]",
    "password[[:space:]]*(=|<-)[[:space:]]*['\"][^'\"]+['\"]",
    "ssl_verify(host|peer)[[:space:]]*=[[:space:]]*(0|FALSE)",
    "opal[.]verifyssl[[:space:]]*=[[:space:]]*FALSE",
    "setSSLVerifyPeer[[:space:]]*[(][[:space:]]*FALSE",
    "curl[^\\n]*(--insecure|[[:space:]]-k([[:space:]]|$))",
    "wget[^\\n]*--no-check-certificate",
    "datashield[.]privacyLevel[[:space:]]*=[[:space:]]*0",
    "base::(c|list|numeric|character)",
    "(localCorDS|psi(GetMatchedIndices|ExportCommonIndices|ExportMasked|ExportMatchedIndices)DS)" )
  pattern <- paste(forbidden, collapse = "|")

  for (path in paths) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_false(grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
                 info = paste(path, "contains an unsafe executable default"))
  }
})

test_that("repository validation vignettes are explicitly archival", {
  package_root <- dirname(.dsvert_client_source_root())
  repo_root <- normalizePath(file.path(package_root, ".."), mustWork = FALSE)
  if (!dir.exists(file.path(repo_root, "vignettes"))) {
    skip("repository-only historical vignettes are absent")
  }
  paths <- c(
    file.path(repo_root, "vignettes", c(
      "method_validation.Rmd", "method_validation_supplementary.Rmd")),
    list.files(file.path(repo_root, "vignettes", "methods"),
               pattern = "[.]Rmd$", full.names = TRUE))
  paths <- paths[file.exists(paths)]
  expect_gt(length(paths), 2L)
  for (path in paths) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(text, "Archived historical validation record", fixed = TRUE,
                 info = paste(path, "lacks the archive warning"))
  }

  rendered <- c(
    list.files(file.path(repo_root, "vignettes"),
               pattern = "[.]html$", full.names = TRUE,
               recursive = FALSE),
    list.files(file.path(repo_root, "vignettes", "methods"),
               pattern = "[.]html$", full.names = TRUE,
               recursive = FALSE))
  expect_gt(length(rendered), 2L)
  for (path in rendered) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(text, "Archived historical validation record", fixed = TRUE,
                 info = paste(path, "is a stale rendered validation report"))
    expect_false(grepl("live verdict|opal-demo[.]obiba[.]org", text,
                       ignore.case = TRUE, perl = TRUE),
                 info = paste(path, "still advertises retired live evidence"))
  }
})
