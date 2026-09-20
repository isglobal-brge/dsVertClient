test_that("family reference harness remains source-only and explicitly tagged", {
  package_path <- getNamespaceInfo(asNamespace("dsVertClient"), "path")
  if (file.exists(file.path(package_path, "Meta", "package.rds"))) {
    expect_identical(system.file("tools", "validate_fama_reference_dslite.R",
                                package = "dsVertClient"), "")
  } else {
    # pkgload's system.file shim deliberately resolves source-only tools too.
    expect_false(file.exists(file.path(package_path, "inst", "tools",
                                       "validate_fama_reference_dslite.R")))
  }
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
  script <- file.path(root, "tools", "validate_fama_reference_dslite.R")
  skip_if_not(file.exists(script), "source-tree harness unavailable after install")
  source <- paste(readLines(script, warn = FALSE), collapse = "\n")
  expect_match(source, '"test", "-c", "-tags", "dsvert_family_reference"', fixed = TRUE)
  expect_match(source, "production_release_validated = FALSE", fixed = TRUE)
  expect_match(source, "for (epsilon in c(1,4,8))", fixed = TRUE)
  expect_match(source, "finally = unlink(work, recursive = TRUE)", fixed = TRUE)
  help <- processx::run(file.path(R.home("bin"), "Rscript"), c(script, "--help"))
  expect_identical(help$status, 0L)
  expect_match(help$stdout, "no production fallback", fixed = TRUE)
})

test_that("two-custodian synthetic reference compares signed grid selections", {
  skip_if_not(identical(Sys.getenv("DSVERT_RUN_FAMA_REFERENCE"), "1"),
    "run source harness explicitly; this is not production release evidence")
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  script <- file.path(root, "tools", "validate_fama_reference_dslite.R")
  skip_if_not(file.exists(script), "source-tree harness unavailable after install")
  output <- Sys.getenv("DSVERT_FAMA_REFERENCE_OUTPUT", "")
  if (!nzchar(output)) {
    output <- tempfile(fileext = ".json")
    on.exit(unlink(output), add = TRUE)
  }
  result <- processx::run(file.path(R.home("bin"), "Rscript"),
    c(script, paste0("--output=", output)), error_on_status = FALSE, timeout = 3600)
  expect_identical(result$status, 0L, info = paste(result$stdout, result$stderr))
  if (result$status != 0L) return(invisible(NULL))
  evidence <- jsonlite::read_json(output, simplifyVector = FALSE)
  expect_true(evidence$reference_only)
  expect_false(evidence$production_release_validated)
  expect_identical(evidence$validation_context$f50_json_encoding,
                   "canonical_decimal_strings")
  expect_length(evidence$cases, 12L)
  expect_setequal(vapply(evidence$cases, function(case)
    paste(case$family, case$epsilon, sep = ":"), character(1L)),
    as.vector(outer(c("nb", "binomial", "poisson", "gaussian"), c(1, 4, 8),
                    paste, sep = ":")))
  for (case in evidence$cases) {
    expect_identical(case$agrees_with_exact, case$selected == case$exact_best)
    expect_true(case$production_failed_closed)
    expect_true(case$missing_signature_failed_closed)
    expect_true(case$replay_identical)
    expect_true(all(is.finite(unlist(case$central_coefficients))))
    expect_true(is.finite(case$central_objective))
    expect_true(is.finite(case$selected_exact_objective))
  }
})
