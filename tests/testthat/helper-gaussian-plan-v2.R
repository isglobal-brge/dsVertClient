.client_complete_gaussian_plan_v2 <- function(plan) {
  stopifnot(is.list(plan), !is.null(names(plan)), !anyDuplicated(names(plan)))
  defaults <- list(
    reference = "CKS2020-theorem14-equation16-corollary17",
    epsilon_numerator = "1", epsilon_denominator = "1",
    allocated_delta_numerator = "1",
    allocated_delta_denominator = "1000000000000000000000000000000",
    core_delta_numerator = "1",
    core_delta_denominator = "2000000000000000000000000000000",
    tail_delta_numerator = "1",
    tail_delta_denominator = "4000000000000000000000000000000",
    l2_sensitivity_numerator = "1", l2_sensitivity_denominator = "1",
    rho_numerator = "1", rho_denominator = "2",
    zcdp_log_upper_integer = "100",
    zcdp_conversion_exponent_numerator = "1",
    zcdp_conversion_exponent_denominator = "2",
    sigma_squared_numerator = "25", sigma_squared_denominator = "1",
    proposal_scale = "5", maximum_noise_magnitude_per_peer = "1048576",
    maximum_noise_magnitude_two_peers = "2097152",
    tail_proof_exponent_numerator = "100",
    tail_proof_exponent_denominator = "1",
    tail_proof_target_numerator = "1",
    tail_proof_target_denominator = "1000000000000000000000000000000",
    sampler_full_scan_steps = as.numeric(plan$sampler_magnitude_count),
    sampler_cdf_table_bytes = as.numeric(plan$sampler_magnitude_count) *
      as.numeric(plan$sampler_random_bytes_per_coordinate),
    simultaneous_95_abs = "25",
    accuracy_accounting = "fixed public Gaussian fixture accounting",
    accounting = "fixed public Gaussian fixture DP accounting",
    privacy_theorem = "CKS2020 Theorem 14 + Equation 16 + Corollary 17",
    nominal_variance_multiplier = 2,
    nominal_standard_deviation_factor =
      "sqrt(2)_relative_to_one_full_draw",
    at_least_one_honest_noise_peer = TRUE,
    maximum_colluding_noise_peers = 1,
    adversary_view = paste0(
      "analyst_plus_at_most_one_designated_noise_peer_including_",
      "its_seed_draw_source_share_and_protocol_transcript"),
    adversary_view_privacy_argument = paste0(
      "conditioned_on_a_simulatable_own_share_and_fixed_corrupt_",
      "peer_view_the_other_independent_complete_epsilon_full_",
      "sensitivity_draw_is_an_epsilon_delta_DP_mechanism; own_",
      "draw_translation_second_draw_signed_decode_and_public_",
      "clamp_are_post_processing; release_delta_is_max_of_the_",
      "two_symmetric_conditional_guarantees"),
    source_share_hiding_precondition = paste0(
      "each_single_pre_noise_aggregate_share_is_computationally_",
      "simulatable_without_the_protected_query_under_",
      "authenticated_semi_honest_fanin"),
    sampler_branches_on_private_randomness = FALSE,
    transcript_dp_claim = TRUE,
    logical_transcript_fixed_shape = TRUE,
    physical_timing_dp_claim = FALSE,
    unavailable_reason = "")
  for (field in names(defaults)) {
    if (is.null(plan[[field]])) plan[[field]] <- defaults[[field]]
  }
  plan$sampler_branches_on_private_randomness <- FALSE
  plan$transcript_dp_claim <- TRUE
  plan$logical_transcript_fixed_shape <- TRUE
  plan$physical_timing_dp_claim <- FALSE
  plan[.dsvert_vector_gaussian_plan_fields()]
}

.client_fresh_go_vector_plan <- local({
  binary <- NULL
  function(request, command) {
    go <- Sys.which("go")
    if (is.null(binary) || !file.exists(binary)) {
      source <- Sys.getenv("DSVERT_SERVER_SOURCE", unset = "")
      if (nzchar(source)) source <- file.path(source, "inst", "dsvert-mpc")
      if (!dir.exists(source)) {
        source <- testthat::test_path(
          "..", "..", "..", "dsVert", "inst", "dsvert-mpc")
      }
      if (dir.exists(source) && nzchar(go)) {
        source <- normalizePath(source, mustWork = TRUE)
        binary <<- tempfile("dsvert-client-fresh-go-")
        output <- withr::with_dir(source, system2(
          go, c("build", "-o", binary, "."), stdout = TRUE, stderr = TRUE))
        status <- attr(output, "status")
        testthat::skip_if(!is.null(status) && status != 0L,
                          paste(output, collapse = "\n"))
        Sys.chmod(binary, mode = "0700")
      } else {
        installed <- tryCatch(
          dsVert:::.findMpcBinary(), error = identity)
        if (inherits(installed, "error")) {
          testthat::skip(paste(
            "the installed dsVert MPC runtime is unavailable:",
            conditionMessage(installed)))
        }
        binary <<- installed
      }
    }
    input <- tempfile("vector-plan-input-")
    output <- tempfile("vector-plan-output-")
    error <- tempfile("vector-plan-error-")
    on.exit(unlink(c(input, output, error), force = TRUE), add = TRUE)
    writeLines(.dsvert_joint_dp_client_json(request), input, useBytes = TRUE)
    status <- system2(
      binary, command, stdin = input,
      stdout = output, stderr = error)
    if (!identical(status, 0L)) {
      stop(paste(readLines(c(output, error), warn = FALSE), collapse = "\n"),
           call. = FALSE)
    }
    jsonlite::read_json(output, simplifyVector = TRUE)
  }
})

.client_fresh_go_gaussian_plan_v2 <- function(request) {
  .client_fresh_go_vector_plan(
    request, "joint-dp-vector-gaussian-plan-v2")
}

.client_fresh_go_laplace_plan_v3 <- function(request) {
  .client_fresh_go_vector_plan(
    request, "joint-dp-vector-laplace-plan-v3")
}
