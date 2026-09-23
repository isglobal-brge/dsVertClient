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
    dsvert_run_dp_statistical_validation(8L, sampler = "ideal")
  expect_identical(.Random.seed, caller_seed)
  second <- .dp_statistical_validation_environment$
    dsvert_run_dp_statistical_validation(8L, sampler = "ideal")
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
    dsvert_run_dp_statistical_validation(2L, sampler = "ideal")
  output_dir <- withr::local_tempdir()
  paths <- .dp_statistical_validation_environment$
    .dv_write_validation(result, output_dir)

  expect_length(paths, 5L)
  expect_true(all(file.exists(paths)))
  report <- readLines(paths[[4L]], warn = FALSE)
  expect_true(any(grepl("not DSI E2E", report, fixed = TRUE)))
  expect_true(any(grepl(
    "No population or sampling confidence interval", report,
    fixed = TRUE)))
  expect_true(any(grepl(
    "peer pinning, DSI transport", report, fixed = TRUE)))
  replay <- jsonlite::read_json(paths[[5L]], simplifyVector = FALSE)
  expect_identical(replay$metadata$sampler_mode, "ideal")
  expect_length(replay$sampler_records$draws, 12L)
  expect_true(all(vapply(replay$sampler_records$draws, function(draw) {
    exact <- as.numeric(unlist(draw$exact_coordinate_integer))
    noise <- as.numeric(unlist(draw$noise_integer))
    released <- as.numeric(unlist(draw$released_coordinate_integer))
    identical(released, pmin(as.numeric(draw$clamp_upper_integer),
                             pmax(0, exact + noise)))
  }, logical(1L))))
})


test_that("the battery defaults to production and records replayable requests", {
  env <- .dp_statistical_validation_environment
  binary <- Sys.getenv("DSVERT_MPC_BINARY", unset = "")
  skip_if(!nzchar(binary), "set DSVERT_MPC_BINARY to run the production oracle")
  result <- env$dsvert_run_dp_statistical_validation(1L, sampler_binary = binary)
  repeated <- env$dsvert_run_dp_statistical_validation(1L, sampler_binary = binary)
  expect_identical(result$sampler_records, repeated$sampler_records)
  expect_identical(result$metadata$sampler_mode, "production")
  expect_identical(result$metadata$randomness, "keyed-stream-computational")
  expect_identical(result$metadata$guarantee, "pure-dp-under-ideal-bits")
  expect_true(all(vapply(result$sampler_records$oracle_plans, function(record) {
    identical(record$request$allocated_delta, "0") &&
      identical(record$response$plan$implementation_delta_numerator, "0") &&
      identical(record$response$plan$wrap_bound_certified, TRUE) &&
      !identical(record$response$plan$representability_bound, "0")
  }, logical(1L))))
  expect_length(result$sampler_records$oracle_plans, 5L)
  expect_length(result$sampler_records$oracle_batches, 6L)
  expect_length(result$sampler_records$draws, 6L)
  expect_true(all(result$gates$pass))
  batches <- result$sampler_records$oracle_batches
  expect_true(all(vapply(batches, function(batch) {
    request <- batch$request
    all(grepl("^[0-9a-f]{64}$", c(
      request$garbler_seed, request$evaluator_seed,
      request$release_contract_hash, request$transcript_hash))) &&
      request$garbler_seed != request$evaluator_seed &&
      identical(batch$response$plan_version, env$.dv_laplace_plan)
  }, logical(1L))))
  expect_true(all(vapply(result$sampler_records$draws, function(draw) {
    exact <- as.numeric(draw$exact_coordinate_integer)
    noise <- as.numeric(draw$noise_integer)
    released <- as.numeric(draw$released_coordinate_integer)
    identical(released, pmin(as.numeric(draw$clamp_upper_integer),
                             pmax(0, exact + noise)))
  }, logical(1L))))
  changed <- batches[[1L]]$response
  changed$guarantee <- "approximate-dp-under-ideal-bits"
  expect_error(env$.dv_check_oracle(changed, batches[[1L]]$request),
               "inconsistent privacy metadata")
  expect_null(env$.dv_sampler_state)
})

test_that("the ideal comparison bypasses production and unknown modes fail", {
  env <- .dp_statistical_validation_environment
  original <- env$.dv_call_sampler
  withr::defer(assign(".dv_call_sampler", original, envir = env))
  env$.dv_call_sampler <- function(...) stop("production called")
  result <- env$dsvert_run_dp_statistical_validation(1L, sampler = "ideal")
  expect_identical(result$metadata$sampler_mode, "ideal")
  expect_length(result$sampler_records$oracle_batches, 0L)
  expect_error(env$dsvert_run_dp_statistical_validation(
    1L, sampler = "approximate"), "arg")
})

test_that("the battery rejects stale production sampler output", {
  env <- .dp_statistical_validation_environment
  original <- env$.dv_call_sampler
  withr::defer(assign(".dv_call_sampler", original, envir = env))
  env$.dv_call_sampler <- function(...) list(
    version = "dsvert-dp-statistical-noise-oracle-output-v1",
    sampler = "hkdf-sha256-chacha20-independent-full-draw-binary-geometric-tv-v3")
  expect_error(env$dsvert_run_dp_statistical_validation(
    1L, sampler_binary = "synthetic-test-binary"), "incompatible contract")
  expect_null(env$.dv_sampler_state)
})
