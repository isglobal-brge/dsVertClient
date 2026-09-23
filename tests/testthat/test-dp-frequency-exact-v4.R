test_that("Frequency v4 production vectors preserve zero delta and modular certificates", {
  config <- jsonlite::read_json(test_path("fixtures", "frequency-v4-exact.json"),
                               simplifyVector = FALSE)
  plan <- .dsvert_dp_frequency_client_plan_v1(config)
  expect_identical(plan$version, "dsvert-frequency-plan-summary-v2")
  expect_identical(plan$implementation_delta, list(numerator = "0", denominator = "1"))
  expect_identical(plan$maximum_noise_per_peer, "unbounded")
  expect_null(plan$no_wrap_sha256)
  expect_true(plan$wrap_certificate$wrap_bound_certified)
  expect_identical(plan$wrap_certificate$guarantee, "pure-dp-under-ideal-bits")
  expect_identical(plan$wrap_certificate$randomness, "keyed-stream-computational")
  expect_identical(plan$wrap_certificate$representability_bound,
    "1e-56713727820156410577229101238628035241")
  expect_identical(plan$wrap_certificate$sum_wrap_bound,
    "1e-28356863910078205288614550619314017603")
  expect_equal(.dsvert_dp_frequency_accuracy_v1(config, .95)$radius, 18)
  expect_equal(.dsvert_dp_frequency_accuracy_v1(config, .99)$radius, 24)
  expect_identical(.dsvert_dp_frequency_client_plan_v1(config), plan)
  for (field in c("representability_bound", "wrap_bound", "guarantee", "randomness")) {
    tampered <- config
    tampered$backend_selection$selected_plan[[field]] <- "0"
    expect_error(.dsvert_dp_frequency_client_plan_v1(tampered), "certificate")
  }
  tampered <- config
  tampered$backend_selection$selected_plan$no_wrap_headroom_certified <- TRUE
  expect_error(.dsvert_dp_frequency_client_plan_v1(tampered), "certificate")
  tampered <- config
  tampered$backend_selection$selected_accuracy_certificate$simultaneous_95_abs <- "17"
  expect_error(.dsvert_dp_frequency_client_plan_v1(tampered), "certificate")
})

test_that("Frequency v4 selection retains a positive-delta Gaussian route", {
  config <- jsonlite::read_json(test_path("fixtures", "frequency-v4-positive-1000.json"),
                               simplifyVector = FALSE)
  plan <- .dsvert_dp_frequency_client_plan_v1(config)
  expect_match(plan$physical_plan_version, "gaussian")
  expect_false(identical(plan$core_delta$numerator, "0"))
  expect_false(identical(plan$implementation_delta$numerator, "0"))
  expect_match(plan$no_wrap_sha256, "^[0-9a-f]{64}$")
  expect_null(plan$wrap_certificate)
  expect_equal(.dsvert_dp_frequency_accuracy_v1(config, .95)$radius, 28)
  config$privacy$delta <- 0
  config$calibration$implementation_delta <- 0
  expect_error(.dsvert_dp_frequency_client_plan_v1(config), "backend selection")
})
