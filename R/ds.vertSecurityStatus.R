.dsvert_security_route_claims_valid <- function(value) {
  required <- c(
    "schema_version",
    "biomedical_joint_dp_capsule_profile_surface_eligible",
    "biomedical_joint_dp_capsule_runtime_readiness",
    "formal_glm_ready", "formal_glm_state",
    "formal_glm_public_result_ready", "formal_glm_public_result_state",
    "formal_cox_ready", "formal_cox_state",
    "formal_cox_public_result_ready", "formal_cox_public_result_state")
  is.list(value) && identical(sort(names(value)), sort(required)) &&
    is.numeric(value$schema_version) &&
    length(value$schema_version) == 1L &&
    is.finite(value$schema_version) &&
    identical(as.numeric(value$schema_version), 1) &&
    is.logical(
      value$biomedical_joint_dp_capsule_profile_surface_eligible) &&
    length(
      value$biomedical_joint_dp_capsule_profile_surface_eligible) == 1L &&
    !is.na(value$biomedical_joint_dp_capsule_profile_surface_eligible) &&
    identical(
      value$biomedical_joint_dp_capsule_runtime_readiness,
      "not_evaluated_requires_client_joint_dp_status_handshake") &&
    is.logical(value$formal_glm_ready) &&
    length(value$formal_glm_ready) == 1L &&
    !is.na(value$formal_glm_ready) &&
    identical(
      value$formal_glm_state,
      "sealed_no_registered_r_dsi_joint_dp_release_lifecycle") &&
    is.logical(value$formal_glm_public_result_ready) &&
    length(value$formal_glm_public_result_ready) == 1L &&
    !is.na(value$formal_glm_public_result_ready) &&
    identical(value$formal_glm_public_result_state,
              paste0("read_only_completed_two_authority_signed_public_",
                     "certificate")) &&
    is.logical(value$formal_cox_ready) &&
    length(value$formal_cox_ready) == 1L &&
    !is.na(value$formal_cox_ready) &&
    identical(
      value$formal_cox_state,
      paste0("sealed_no_recipient_encrypted_r_dsi_lifecycle_or_",
             "end_to_end_numeric_certificate")) &&
    is.logical(value$formal_cox_public_result_ready) &&
    length(value$formal_cox_public_result_ready) == 1L &&
    !is.na(value$formal_cox_public_result_ready) &&
    identical(value$formal_cox_public_result_state,
              paste0("read_only_completed_two_authority_signed_sticky_",
                     "opening_certificate"))
}

.dsvert_security_profile_validate <- function(value, server) {
  required <- c(
    "schema_version", "release_mode", "exact_adaptive_releases_enabled",
    "disclosure_safe_gate_active", "deployment_surface_attested",
    "remote_surface_profile", "remote_surface_attestation_authority",
    "remote_surface_attestation_state", "formal_dp_claim_eligible",
    "route_claims",
    "unconditional_non_reconstruction_guarantee",
    "mpc_transport_is_opaque_to_analyst", "caveat")
  logical_fields <- c(
    "exact_adaptive_releases_enabled", "disclosure_safe_gate_active",
    "deployment_surface_attested", "formal_dp_claim_eligible",
    "unconditional_non_reconstruction_guarantee",
    "mpc_transport_is_opaque_to_analyst")
  valid_states <- c(
    "verified_custodian_attestation",
    "missing_custodian_attestation",
    "invalid_custodian_attestation",
    "conflicting_custodian_attestation",
    "mismatched_custodian_attestation",
    "server_contract_unavailable")
  valid <- is.list(value) && identical(sort(names(value)), sort(required)) &&
    is.numeric(value$schema_version) && length(value$schema_version) == 1L &&
    is.finite(value$schema_version) &&
    identical(as.numeric(value$schema_version), 5) &&
    is.character(value$release_mode) && length(value$release_mode) == 1L &&
    identical(value$release_mode, "disclosure_safe") &&
    all(vapply(logical_fields, function(field) {
      is.logical(value[[field]]) && length(value[[field]]) == 1L &&
        !is.na(value[[field]])
    }, logical(1L))) &&
    is.character(value$remote_surface_profile) &&
    length(value$remote_surface_profile) == 1L &&
    identical(value$remote_surface_profile, "dsvert") &&
    is.character(value$remote_surface_attestation_authority) &&
    length(value$remote_surface_attestation_authority) == 1L &&
    identical(
      value$remote_surface_attestation_authority,
      "custodian_owned_deployment_assertion") &&
    is.character(value$remote_surface_attestation_state) &&
    length(value$remote_surface_attestation_state) == 1L &&
    !is.na(value$remote_surface_attestation_state) &&
    value$remote_surface_attestation_state %in% valid_states &&
    .dsvert_security_route_claims_valid(value$route_claims) &&
    is.character(value$caveat) && length(value$caveat) == 1L &&
    !is.na(value$caveat) && nzchar(value$caveat)
  if (!valid) {
    stop("Server '", server, "' returned an invalid security profile",
         call. = FALSE)
  }
  if (!identical(value$disclosure_safe_gate_active, TRUE) ||
      !identical(value$exact_adaptive_releases_enabled, FALSE) ||
      !identical(value$formal_dp_claim_eligible,
                 value$deployment_surface_attested) ||
      !identical(
        value$route_claims[[
          "biomedical_joint_dp_capsule_profile_surface_eligible"]],
        value$deployment_surface_attested) ||
      !identical(
        value$formal_dp_claim_eligible,
        value$route_claims[[
          "biomedical_joint_dp_capsule_profile_surface_eligible"]]) ||
      isTRUE(value$route_claims$formal_glm_ready) ||
      isTRUE(value$route_claims$formal_cox_ready) ||
      !isTRUE(value$route_claims$formal_glm_public_result_ready) ||
      !isTRUE(value$route_claims$formal_cox_public_result_ready) ||
      !identical(
        value$deployment_surface_attested,
        identical(value$remote_surface_attestation_state,
                  "verified_custodian_attestation")) ||
      isTRUE(value$unconditional_non_reconstruction_guarantee) ||
      !identical(value$mpc_transport_is_opaque_to_analyst, TRUE)) {
    stop("Server '", server, "' returned a contradictory security profile",
         call. = FALSE)
  }
  value
}

.dsvert_security_status_impl <- function(
    datasources = NULL, require_ready = TRUE,
    .aggregate = DSI::datashield.aggregate,
    .dp_status = .dsvert_dp_status_impl) {
  datasources <- .dsvert_dp_datasources(datasources)
  if (!is.logical(require_ready) || length(require_ready) != 1L ||
      is.na(require_ready)) {
    stop("require_ready must be TRUE or FALSE", call. = FALSE)
  }
  profiles <- .dsvert_aggregate_strict(
    conns = datasources, expr = call(name = "dsvertSecurityProfileDS"),
    operation = "security profile handshake", .aggregate = .aggregate)
  if (!is.list(profiles) ||
      !identical(sort(names(profiles)), sort(names(datasources)))) {
    stop("The security-profile handshake returned an invalid server set",
         call. = FALSE)
  }
  profiles <- stats::setNames(lapply(names(profiles), function(server) {
    .dsvert_security_profile_validate(profiles[[server]], server)
  }), names(datasources))
  modes <- unique(vapply(profiles, `[[`, character(1L), "release_mode"))
  if (length(modes) != 1L) {
    stop("Connected servers disagree on the disclosure release mode",
         call. = FALSE)
  }
  if (!identical(modes, "disclosure_safe")) {
    stop("Connected servers do not expose dsVert's disclosure-safe profile",
         call. = FALSE)
  }
  dp <- tryCatch(
    .dp_status(datasources, .aggregate),
    error = function(e) e)
  surface_ready <- all(vapply(
    profiles, `[[`, logical(1L), "deployment_surface_attested"))
  dp_ready <- !inherits(dp, "error")
  biomedical_profile_eligible <- all(vapply(profiles, function(profile) {
    profile$route_claims[[
      "biomedical_joint_dp_capsule_profile_surface_eligible"]]
  }, logical(1L)))
  route_readiness <- list(
    schema_version = 3L,
    biomedical_joint_dp_capsule_profile_surface_eligible =
      biomedical_profile_eligible,
    biomedical_joint_dp_capsule_runtime_policy_consortium_ready = dp_ready,
    biomedical_joint_dp_capsule_ready =
      biomedical_profile_eligible && dp_ready,
    execution_readiness =
      "not_evaluated_requires_route_specific_preflight",
    formal_glm_ready = FALSE,
    formal_glm_state = profiles[[1L]]$route_claims$formal_glm_state,
    formal_glm_public_result_ready = TRUE,
    formal_glm_public_result_state =
      profiles[[1L]]$route_claims$formal_glm_public_result_state,
    formal_cox_ready = FALSE,
    formal_cox_state = profiles[[1L]]$route_claims$formal_cox_state,
    formal_cox_public_result_ready = TRUE,
    formal_cox_public_result_state =
      profiles[[1L]]$route_claims$formal_cox_public_result_state)
  ready <- route_readiness$biomedical_joint_dp_capsule_ready
  warnings <- character()
  if (!surface_ready) {
    unavailable <- names(profiles)[!vapply(
      profiles, `[[`, logical(1L), "deployment_surface_attested")]
    states <- vapply(
      profiles[unavailable], `[[`, character(1L),
      "remote_surface_attestation_state")
    warnings <- c(warnings, paste0(
      "The custodian-owned exclusive DataSHIELD surface attestation is missing ",
      "or stale on: ",
      paste0(unavailable, " (", states, ")", collapse = ", "),
      ". Reconcile and re-attest the dedicated dsvert surface with ",
      "connector-specific admin tooling after verifying the effective ",
      "callable inventory."))
  }
  if (!dp_ready) warnings <- c(warnings, conditionMessage(dp))
  warning <- if (length(warnings)) paste(warnings, collapse = "; ") else NULL
  if (isTRUE(require_ready) && !ready) {
    stop("The disclosure-safe consortium DP contract is not ready: ",
         warning, call. = FALSE)
  }
  result <- list(
    ready = ready,
    release_mode = modes,
    surface_attested = surface_ready,
    route_readiness = route_readiness,
    profiles = profiles,
    dp_status = if (inherits(dp, "error")) NULL else dp,
    warning = warning)
  class(result) <- c("ds.vertSecurityStatus", "list")
  result
}

#' Verify the consortium disclosure profile
#'
#' Checks that every connected server enforces dsVert's single disclosure-safe
#' profile, a matching custodian-owned attestation of the dedicated logical
#' dsVert surface and the coherent signed Synopsis policy and pinset. The
#' surface
#' attestation is an
#' administrative assertion, not live introspection: it must be renewed after
#' an administrator changes the effective callable surface. Opal and
#' Armadillo/Rock use the same logical contract and token; connector-specific
#' admin tooling provisions it in the server profile or server/container
#' environment only. This client sends no attestation option, token or claim in
#' the DSI expression. Under security-profile schema v5, the
#' returned route map reports biomedical joint-DP profile eligibility and
#' authenticated control-plane readiness separately. It explicitly does not
#' evaluate route-specific dataset admission, manifest construction, numeric
#' runtime capabilities, or a live release. The top-level `ready` value and
#' server compatibility alias `formal_dp_claim_eligible` apply only to that
#' biomedical route; they never promote formal GLM or formal Cox compute,
#' whose route-specific `ready` fields remain false while sealed. It does
#' report the separate read-only route for a completed two-authority GLM or
#' Cox certificate.
#'
#' @param datasources DataSHIELD connections; active connections by default.
#' @param require_ready Fail if the custodian-attested remote surface or
#'   current sticky-artifact consortium policy is not ready. Historical
#'   capsule lifetime telemetry is not an admission gate. This is not a
#'   route-specific execution preflight and never
#'   promotes formal GLM or Cox compute.
#' @return A `ds.vertSecurityStatus` object with explicit route readiness.
#' @export
ds.vertSecurityStatus <- function(datasources = NULL,
                                  require_ready = TRUE) {
  .dsvert_security_status_impl(
    datasources, require_ready,
    DSI::datashield.aggregate, .dsvert_dp_status_impl)
}

#' @export
print.ds.vertSecurityStatus <- function(x, ...) {
  cat("dsVert release mode:", x$release_mode,
      " | policy/control-plane ready:",
      if (isTRUE(x$ready)) "yes" else "no", "\n")
  routes <- x$route_readiness
  cat(
    "Routes: biomedical policy/control-plane:",
    if (isTRUE(routes$biomedical_joint_dp_capsule_ready)) "yes" else "no",
    "| formal GLM:", if (isTRUE(routes$formal_glm_ready)) "yes" else "no",
    "| formal GLM completed public result:",
    if (isTRUE(routes$formal_glm_public_result_ready)) "yes" else "no",
    "| formal Cox:", if (isTRUE(routes$formal_cox_ready)) "yes" else "no",
    "| formal Cox completed public result:",
    if (isTRUE(routes$formal_cox_public_result_ready)) "yes" else "no",
    "\n")
  cat("Execution readiness:", routes$execution_readiness, "\n")
  if (!is.null(x$warning)) cat("Warning:", x$warning, "\n")
  invisible(x)
}
