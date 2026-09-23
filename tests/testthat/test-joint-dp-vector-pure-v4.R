.vector_pure_fixture <- function() jsonlite::read_json(testthat::test_path(
  "fixtures", "joint-dp-vector-convolution-v4.json"), simplifyVector = FALSE)

test_that("exact production v4 plan binds pure privacy and a nonzero utility bound", {
  fixture <- .vector_pure_fixture()
  plan <- fixture$output$plan
  profile <- .dsvert_vector_profile("discrete-laplace",
    backend = .DSVERT_CLIENT_VECTOR_PURE_BACKEND)
  expect_true(profile$pure)
  expect_false(profile$exact_gc)
  expect_silent(.dsvert_vector_plan_validate(plan, .dsvert_vector_hash(plan),
    profile, 3, "256"))
  expect_identical(.dsvert_vector_implementation_delta(list(
    profile = profile, mechanism_plan = plan)), "0/1")
  expect_identical(plan$guarantee, "pure-dp-under-ideal-bits")
  expect_identical(plan$randomness, "keyed-stream-computational")
  expect_identical(plan$wrap_bound, plan$representability_bound)
  expect_match(plan$representability_bound, "^1e-[1-9][0-9]*$")
  expect_gt(.dsvert_vector_exact_probability_upper(plan$wrap_bound), 0)
  expect_identical(plan$maximum_noise_magnitude, "unbounded")
  expect_false(any(grepl("no_wrap", names(plan))))
  expect_identical(.dsvert_dp_vector_sampler_tv_upper(plan, FALSE), 0)
  legacy <- .dsvert_vector_profile("discrete-laplace",
    backend = .DSVERT_CLIENT_VECTOR_BACKEND)
  expect_false(isTRUE(legacy$pure))
  expect_match(legacy$plan_version, "v3$")
})

test_that("exact v4 readers reject misbound ranges and privacy certificates", {
  plan <- .vector_pure_fixture()$output$plan
  mutations <- list(
    guarantee = "approximate-dp-under-ideal-bits", randomness = "ideal",
    noise_support = "bounded_integer", noise_commitment = "clamp",
    representability_bound = "0", wrap_bound = "0", wrap_bound_certified = FALSE,
    wrap_bound_event = "no_wrap", maximum_noise_magnitude = "1000",
    per_peer_implementation_delta_numerator = "1",
    per_peer_implementation_delta_bound = "1e-9",
    epsilon_effective_upper_numerator = "10001",
    epsilon_effective_upper_denominator = "1267650600228229401496703205376")
  for (field in names(mutations)) {
    altered <- plan
    altered[[field]] <- mutations[[field]]
    expect_false(.dsvert_vector_exact_plan_valid(altered, "256", 3), info = field)
  }
  altered <- plan
  altered$no_wrap_headroom_certified <- TRUE
  expect_false(.dsvert_vector_exact_plan_valid(altered, "256", 3))
  for (dimension in c(0, 1000001, 1.5)) {
    expect_false(.dsvert_vector_exact_plan_valid(plan, "256", dimension))
  }
  for (sensitivity in c("0", "-1", "1.5", "512")) {
    expect_false(.dsvert_vector_exact_plan_valid(plan, sensitivity, 3))
  }
})

test_that("exact v4 rational tail envelopes the analytic convolution", {
  plan <- .vector_pure_fixture()$output$plan
  tail <- .dsvert_dp_vector_dyadic_tail_context(plan)
  for (radius in c(0, 1, 128, 1024, 10000)) {
    expect_gte(.dsvert_dp_vector_plan_log_tail_upper(radius, tail, TRUE),
      .dsvert_dp_vector_convolution_log_tail(radius, 1, 256))
  }
  certificate <- .dsvert_vector_exact_sum_certificate(plan,
    openssl::bignum("9007199254740991"), 3)
  expect_identical(certificate$sum_wrap_threshold,
    "85070591730234615865839148258314682368")
  expect_match(certificate$sum_wrap_bound, "^1e-[1-9][0-9]*$")
})

test_that("production v4 known-answer draws replay exactly", {
  binary <- Sys.getenv("DSVERT_MPC_BINARY", unset = "")
  skip_if(!nzchar(binary), "set DSVERT_MPC_BINARY for production replay")
  fixture <- .vector_pure_fixture()
  input <- tempfile()
  output <- tempfile()
  withr::defer(unlink(c(input, output)))
  jsonlite::write_json(fixture$input, input, auto_unbox = TRUE)
  run <- function() {
    status <- system2(binary, "joint-dp-vector-convolution-oracle-v1",
      stdin = input, stdout = output)
    expect_identical(status, 0L)
    jsonlite::read_json(output, simplifyVector = FALSE)
  }
  first <- run()
  expect_identical(first, fixture$output)
  expect_identical(run(), first)
})

test_that("public v4 selection permits zero delta without a finite-plan prerequisite", {
  hex <- strrep("a", 64)
  for (policy in c("dsvert-joint-dp-vector-exact-gc-cost-policy-v2",
                   "dsvert-joint-dp-vector-pure-laplace-policy-v1")) {
    pure_policy <- grepl("pure-laplace", policy, fixed = TRUE)
    assessment <- list(version = "dsvert-joint-dp-vector-exact-gc-assessment-v2",
      manifest_sha256 = hex, representable = TRUE,
      exact_gc_capability_id = "joint_dp_biomedical_vector_exact_gc_v1",
      plan_sha256 = hex, maximum_chunk_coordinates = 8192L,
      cost_policy_version = policy, total_coordinate_count = 9000L,
      maximum_promoted_coordinates = if (pure_policy) 0L else 1L,
      promoted = FALSE, selection_reason = if (pure_policy) {
        "zero_delta_requires_exact_unbounded_laplace"
      } else "above_public_exact_gc_cost_ceiling",
      private_material_accessed = FALSE, runtime_failure_consulted = FALSE)
    assessment$assessment_sha256 <-
      .dsvert_joint_dp_vector_exact_gc_client_hash(assessment)
    expect_silent(.dsvert_joint_dp_vector_exact_gc_client_assessment(assessment, hex))
    selection <- list(version = "dsvert-joint-dp-vector-backend-selection-v2",
      manifest_sha256 = hex, backend = .DSVERT_CLIENT_VECTOR_PURE_BACKEND,
      one_draw = FALSE, cost_policy_version = policy,
      total_coordinate_count = assessment$total_coordinate_count,
      maximum_promoted_coordinates = assessment$maximum_promoted_coordinates,
      selection_reason = assessment$selection_reason,
      assessment_sha256 = assessment$assessment_sha256,
      exact_gc_plan_sha256 = hex, exact_gc_maximum_chunk_coordinates = 8192L,
      selected_before_private_material = TRUE, retry_may_change_backend = FALSE)
    selection$selection_sha256 <- .dsvert_joint_dp_vector_exact_gc_client_hash(selection)
    expect_silent(.dsvert_joint_dp_vector_exact_gc_client_selection(selection, hex, FALSE))
    selection$backend <- .DSVERT_CLIENT_VECTOR_BACKEND
    selection$selection_sha256 <- .dsvert_joint_dp_vector_exact_gc_client_hash(
      selection[setdiff(names(selection), "selection_sha256")])
    expect_error(.dsvert_joint_dp_vector_exact_gc_client_selection(selection, hex), "invalid")
  }
})


test_that("zero-delta staged sources fail before selecting an unsupported gate", {
  manifest <- list(workload = list(mechanism_selection = list(allocated_delta = 0),
    families = list(gaussian_models = list(artifacts = list(cox = list(family = "cox"))))))
  expect_error(.dsvert_dp_glm_grid_cross_noise_policy(manifest),
    "Staged source requires positive delta")
  manifest$workload$families$gaussian_models$artifacts <- list()
  expect_identical(.dsvert_dp_glm_grid_cross_noise_policy(manifest),
    "dsvert-joint-dp-vector-pure-laplace-policy-v1")
})
