.security_profile_fixture <- function(
    attested = TRUE,
    state = if (isTRUE(attested)) "verified_custodian_attestation" else
      "missing_custodian_attestation") {
  list(
    schema_version = 4L,
    release_mode = "disclosure_safe",
    exact_adaptive_releases_enabled = FALSE,
    disclosure_safe_gate_active = TRUE,
    deployment_surface_attested = attested,
    remote_surface_profile = "dsvert",
    remote_surface_attestation_authority =
      "custodian_owned_deployment_assertion",
    remote_surface_attestation_state = state,
    formal_dp_claim_eligible = attested,
    route_claims = list(
      schema_version = 1L,
      biomedical_joint_dp_capsule_profile_surface_eligible = attested,
      biomedical_joint_dp_capsule_runtime_readiness =
        "not_evaluated_requires_client_joint_dp_status_handshake",
      formal_glm_ready = FALSE,
      formal_glm_state =
        "sealed_no_registered_r_dsi_joint_dp_release_lifecycle",
      formal_cox_ready = FALSE,
      formal_cox_state = paste0(
        "sealed_no_recipient_encrypted_r_dsi_lifecycle_or_",
        "end_to_end_numeric_certificate")),
    unconditional_non_reconstruction_guarantee = FALSE,
    mpc_transport_is_opaque_to_analyst = TRUE,
    caveat = "test caveat")
}

test_that("security status requires the single profile and DP readiness", {
  conns <- list(a = structure(list(), class = "mock"),
                b = structure(list(), class = "mock"))
  aggregate <- function(conns, expr, error = NULL, errors.print = TRUE, ...) {
    expect_false(errors.print)
    stats::setNames(rep(list(.security_profile_fixture()), length(conns)),
                    names(conns))
  }
  dp_status <- function(datasources, .aggregate) {
    structure(list(a = list(enabled = TRUE), b = list(enabled = TRUE)),
              class = c("ds.vertDPStatus", "list"))
  }
  result <- .dsvert_security_status_impl(
    conns, TRUE, aggregate, dp_status)
  expect_s3_class(result, "ds.vertSecurityStatus")
  expect_true(result$ready)
  expect_true(result$surface_attested)
  expect_true(
    result$route_readiness[[
      "biomedical_joint_dp_capsule_profile_surface_eligible"]])
  expect_true(
    result$route_readiness[[
      "biomedical_joint_dp_capsule_runtime_policy_consortium_ready"]])
  expect_true(
    result$route_readiness$biomedical_joint_dp_capsule_ready)
  expect_identical(result$route_readiness$schema_version, 2L)
  expect_identical(
    result$route_readiness$execution_readiness,
    "not_evaluated_requires_route_specific_preflight")
  expect_false(result$route_readiness$formal_glm_ready)
  expect_false(result$route_readiness$formal_cox_ready)
  expect_identical(result$release_mode, "disclosure_safe")
  expect_null(result$warning)
  printed <- capture.output(print(result))
  expect_true(any(grepl(
    "biomedical policy/control-plane: yes", printed, fixed = TRUE)))
  expect_true(any(grepl(
    "Execution readiness: not_evaluated_requires_route_specific_preflight",
    printed, fixed = TRUE)))
  expect_false(any(grepl(
    "biomedical joint-DP capsules: yes", printed, fixed = TRUE)))
  expect_true(any(grepl("formal GLM: no", printed, fixed = TRUE)))
  expect_true(any(grepl("formal Cox: no", printed, fixed = TRUE)))
})

test_that("Armadillo readiness carries no client attestation input", {
  conns <- list(
    rock_a = structure(list(), class = "ArmadilloConnection"),
    rock_b = structure(list(), class = "ArmadilloConnection"))
  submitted <- NULL
  aggregate <- function(conns, expr, error = NULL, errors.print = TRUE, ...) {
    submitted <<- expr
    stats::setNames(rep(list(.security_profile_fixture()), length(conns)),
                    names(conns))
  }
  dp_status <- function(datasources, .aggregate) {
    structure(list(rock_a = list(enabled = TRUE),
                   rock_b = list(enabled = TRUE)),
              class = c("ds.vertDPStatus", "list"))
  }
  withr::local_envvar(c(
    DSVERT_REMOTE_SURFACE_ATTESTATION = "analyst-must-not-be-forwarded"))
  withr::local_options(list(
    dsvert.remote_surface_attestation = "analyst-must-not-be-forwarded"))

  result <- .dsvert_security_status_impl(
    conns, TRUE, aggregate, dp_status)
  expect_true(result$ready)
  expect_identical(submitted, call(name = "dsvertSecurityProfileDS"))
  expect_length(as.list(submitted), 1L)
})

test_that("security status rejects obsolete or contradictory profiles", {
  conns <- list(a = structure(list(), class = "mock"),
                b = structure(list(), class = "mock"))
  obsolete <- function(conns, expr, error = NULL, errors.print = TRUE, ...) {
    old <- .security_profile_fixture()
    old$schema_version <- 1L
    old$release_mode <- "legacy_exact_mpc"
    list(a = old, b = old)
  }
  expect_error(
    .dsvert_security_status_impl(conns, FALSE, obsolete, identity),
    "invalid security profile")

  old_v3 <- function(
      conns, expr, error = NULL, errors.print = TRUE, ...) {
    value <- .security_profile_fixture()
    value$schema_version <- 3L
    list(a = value, b = value)
  }
  expect_error(
    .dsvert_security_status_impl(conns, FALSE, old_v3, identity),
    "invalid security profile")

  bad <- function(conns, expr, error = NULL, errors.print = TRUE, ...) {
    value <- .security_profile_fixture()
    value$unconditional_non_reconstruction_guarantee <- TRUE
    list(a = value, b = value)
  }
  expect_error(
    .dsvert_security_status_impl(conns, FALSE, bad, identity),
    "contradictory")

  inconsistent <- function(
      conns, expr, error = NULL, errors.print = TRUE, ...) {
    value <- .security_profile_fixture(FALSE)
    value$formal_dp_claim_eligible <- TRUE
    list(a = value, b = value)
  }
  expect_error(
    .dsvert_security_status_impl(conns, FALSE, inconsistent, identity),
    "contradictory")

  overstated_route <- function(
      conns, expr, error = NULL, errors.print = TRUE, ...) {
    value <- .security_profile_fixture()
    value$route_claims$formal_glm_ready <- TRUE
    list(a = value, b = value)
  }
  expect_error(
    .dsvert_security_status_impl(conns, FALSE, overstated_route, identity),
    "contradictory")

  stale_claim <- function(
      conns, expr, error = NULL, errors.print = TRUE, ...) {
    value <- .security_profile_fixture()
    value$remote_surface_attestation_state <-
      "mismatched_custodian_attestation"
    list(a = value, b = value)
  }
  expect_error(
    .dsvert_security_status_impl(conns, FALSE, stale_claim, identity),
    "contradictory")

  tampered_schema <- function(
      conns, expr, error = NULL, errors.print = TRUE, ...) {
    value <- .security_profile_fixture()
    value$schema_version <- 3.5
    value$remote_surface_profile <- "default"
    list(a = value, b = value)
  }
  expect_error(
    .dsvert_security_status_impl(conns, FALSE, tampered_schema, identity),
    "invalid security profile")

  unknown_state <- function(
      conns, expr, error = NULL, errors.print = TRUE, ...) {
    value <- .security_profile_fixture(FALSE)
    value$remote_surface_attestation_state <- "analyst_asserted"
    list(a = value, b = value)
  }
  expect_error(
    .dsvert_security_status_impl(conns, FALSE, unknown_state, identity),
    "invalid security profile")

  for (state in c(
      "invalid_custodian_attestation",
      "conflicting_custodian_attestation")) {
    expect_identical(
      .dsvert_security_profile_validate(
        .security_profile_fixture(FALSE, state), "rock"),
      .security_profile_fixture(FALSE, state))
  }
})

test_that("missing or stale custodian surface attestation prevents readiness", {
  conns <- list(a = structure(list(), class = "mock"),
                b = structure(list(), class = "mock"))
  aggregate <- function(conns, expr, error = NULL, errors.print = TRUE, ...) {
    list(
      a = .security_profile_fixture(FALSE),
      b = .security_profile_fixture(
        FALSE, "mismatched_custodian_attestation"))
  }
  dp_status <- function(datasources, .aggregate) {
    structure(list(a = list(enabled = TRUE), b = list(enabled = TRUE)),
              class = c("ds.vertDPStatus", "list"))
  }

  result <- .dsvert_security_status_impl(
    conns, FALSE, aggregate, dp_status)
  expect_false(result$ready)
  expect_false(result$surface_attested)
  expect_false(
    result$route_readiness[[
      "biomedical_joint_dp_capsule_profile_surface_eligible"]])
  expect_true(
    result$route_readiness[[
      "biomedical_joint_dp_capsule_runtime_policy_consortium_ready"]])
  expect_false(
    result$route_readiness$biomedical_joint_dp_capsule_ready)
  expect_false(result$route_readiness$formal_glm_ready)
  expect_false(result$route_readiness$formal_cox_ready)
  expect_s3_class(result$dp_status, "ds.vertDPStatus")
  expect_match(
    result$warning, "custodian-owned exclusive DataSHIELD surface")
  expect_match(result$warning, "a \\(missing_custodian_attestation\\)")
  expect_match(result$warning, "b \\(mismatched_custodian_attestation\\)")
  expect_match(result$warning, "Reconcile and re-attest")

  expect_error(
    .dsvert_security_status_impl(conns, TRUE, aggregate, dp_status),
    "not ready.*missing or stale")
})

test_that("readiness can be inspected without weakening the gate", {
  conns <- list(a = structure(list(), class = "mock"))
  aggregate <- function(conns, expr, error = NULL, errors.print = TRUE, ...) {
    list(a = .security_profile_fixture())
  }
  unavailable <- function(...) stop("DP ledger unavailable", call. = FALSE)
  expect_error(
    .dsvert_security_status_impl(conns, TRUE, aggregate, unavailable),
    "not ready")
  result <- .dsvert_security_status_impl(
    conns, FALSE, aggregate, unavailable)
  expect_false(result$ready)
  expect_true(
    result$route_readiness[[
      "biomedical_joint_dp_capsule_profile_surface_eligible"]])
  expect_false(
    result$route_readiness[[
      "biomedical_joint_dp_capsule_runtime_policy_consortium_ready"]])
  expect_false(
    result$route_readiness$biomedical_joint_dp_capsule_ready)
  expect_false(result$route_readiness$formal_glm_ready)
  expect_false(result$route_readiness$formal_cox_ready)
  expect_match(result$warning, "DP ledger unavailable")
})

test_that("lifetime exhaustion does not redefine consortium readiness", {
  conns <- list(a = structure(list(), class = "mock"))
  aggregate <- function(conns, expr, error = NULL, errors.print = TRUE, ...) {
    list(a = .security_profile_fixture())
  }
  dp_status <- function(...) {
    structure(list(a = list(
      composition_telemetry = list(remaining_distinct_capsules = 0),
      release_instance_telemetry = list(
        remaining_distinct_capsules = 0))),
    class = c("ds.vertDPStatus", "list"))
  }

  result <- .dsvert_security_status_impl(
    conns, TRUE, aggregate, dp_status)
  output <- capture.output(returned <- print(result))

  expect_true(result$ready)
  expect_true(result$route_readiness$biomedical_joint_dp_capsule_ready)
  expect_null(result$warning)
  expect_identical(
    result$dp_status$a$composition_telemetry$remaining_distinct_capsules, 0)
  expect_identical(returned, result)
  expect_false(any(grepl("Warning:", output, fixed = TRUE)))
})

test_that("the public security status uses only the joint capsule handshake", {
  namespace <- asNamespace("dsVertClient")
  walk <- function(name, seen = character()) {
    if (name %in% seen || !exists(name, namespace, inherits = FALSE)) {
      return(seen)
    }
    value <- get(name, namespace, inherits = FALSE)
    if (!is.function(value)) return(seen)
    seen <- c(seen, name)
    globals <- tryCatch(
      unique(unlist(codetools::findGlobals(value, merge = FALSE),
                    use.names = FALSE)),
      error = function(error) character())
    for (global in intersect(globals, ls(namespace, all.names = TRUE))) {
      seen <- walk(global, seen)
    }
    seen
  }

  reachable <- walk("ds.vertSecurityStatus")
  expect_true(".dsvert_joint_dp_capsule_status_impl" %in% reachable)
  expect_false(".dsvert_dp_status_impl" %in% reachable)
})
