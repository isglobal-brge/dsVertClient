.dp_statistical_validation_environment <- new.env(parent = globalenv())
.dp_statistical_validation_script <- system.file(
  "scripts", "validate_dp_statistical_methods.R",
  package = "dsVertClient")
if (!nzchar(.dp_statistical_validation_script)) {
  .dp_statistical_validation_script <- testthat::test_path(
    "..", "..", "inst", "scripts",
    "validate_dp_statistical_methods.R")
}
sys.source(
  .dp_statistical_validation_script,
  envir = .dp_statistical_validation_environment)

test_that("DP statistical validation harness is deterministic and gated", {
  set.seed(8675309)
  caller_seed <- .Random.seed

  first <- .dp_statistical_validation_environment$
    dsvert_run_dp_statistical_validation(8L)
  expect_identical(.Random.seed, caller_seed)
  second <- .dp_statistical_validation_environment$
    dsvert_run_dp_statistical_validation(8L)
  expect_identical(.Random.seed, caller_seed)

  stable_columns <- setdiff(
    names(first$summary), "elapsed_seconds_family")
  expect_identical(
    first$summary[stable_columns], second$summary[stable_columns])
  expect_identical(first$contracts, second$contracts)
  expect_identical(first$gates, second$gates)
  expect_identical(first$edge_cases, second$edge_cases)
  expect_true(all(first$gates$pass))
  expect_true(all(first$edge_cases$pass))
  expect_true(all(first$contracts$certified_coverage_lower >= 0.95))
  expect_setequal(
    unique(first$summary$family),
    c("correlation", "describe", "epidemiology", "gaussian", "pca",
      "mantel_haenszel", "survival"))
  expect_setequal(
    unique(first$summary$method),
    c(
      "ds.vertDPDescribe", "ds.vertDPQuantile", "ds.vertDPMedian",
      "ds.vertDPEpi2x2", "ds.vertDPDiagnostic2x2", "ds.vertDPROC",
      "ds.vertDPDirectStandardization",
      "ds.vertDPCausalStandardization",
      "ds.vertDPIndirectStandardization", "ds.vertDPMantelHaenszel",
      "ds.vertDPSurvival",
      "ds.vertDPRMST", "ds.vertDPSurvivalQuantile",
      "ds.vertDPGaussian", "ds.vertCor",
      "ds.vertPCA"))
  point_only <- first$summary$family == "gaussian" |
    first$summary$estimand == "PC1_sine_angle"
  expect_true(all(!first$summary$mechanism_region_available[point_only]))
  expect_true(all(is.na(first$summary$mechanism_coverage[point_only])))
  expect_true(all(first$summary$estimate_finite_rate[point_only] == 1))
  expect_true(all(first$summary$mechanism_region_available[!point_only]))
  mh_gate <- first$gates[
    first$gates$gate ==
      "mantel_haenszel_zero_cost_no_classical_dp_inference", , drop = FALSE]
  expect_equal(nrow(mh_gate), 1L)
  expect_true(mh_gate$pass)
})

test_that("DP statistical validation writer records its scope", {
  result <- .dp_statistical_validation_environment$
    dsvert_run_dp_statistical_validation(2L)
  output_dir <- withr::local_tempdir()
  paths <- .dp_statistical_validation_environment$
    .dv_write_validation(result, output_dir)

  expect_length(paths, 4L)
  expect_true(all(file.exists(paths)))
  report <- readLines(paths[[4L]], warn = FALSE)
  expect_true(any(grepl("not DSI E2E", report, fixed = TRUE)))
  expect_true(any(grepl(
    "No population or sampling confidence interval", report,
    fixed = TRUE)))
  expect_true(any(grepl(
    "HMAC/HKDF/ChaCha20, Ring128, peer pinning, DSI transport", report,
    fixed = TRUE)))
})
