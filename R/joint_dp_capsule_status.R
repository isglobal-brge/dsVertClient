# Client validation for the reusable joint-DP capsule status.  This is kept
# separate from the retired per-release accountant schema so joint methods can
# never confuse remaining distinct-capsule reservation units with request
# slots or decay.

.DSVERT_CLIENT_JOINT_DP_CAPSULE_STATUS_VERSION <-
  "dsvert-joint-dp-capsule-status-v5"

.dsvert_joint_dp_lifetime_decimal_parts <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^(0|[1-9][0-9]*)(\\.[0-9]+)?(e[+-][0-9]+)?$",
             value, perl = TRUE)) return(NULL)
  match <- regexec(
    "^([0-9]+)(?:\\.([0-9]+))?(?:e([+-])([0-9]+))?$",
    value, perl = TRUE)
  parts <- regmatches(value, match)[[1L]]
  if (!length(parts)) return(NULL)
  fractional <- if (length(parts) >= 3L) parts[[3L]] else ""
  exponent <- 0
  if (length(parts) >= 5L && nzchar(parts[[4L]])) {
    exponent <- suppressWarnings(as.numeric(paste0(
      parts[[4L]], parts[[5L]])))
  }
  if (length(exponent) != 1L || is.na(exponent) || !is.finite(exponent) ||
      exponent != floor(exponent) || abs(exponent) > 4096) return(NULL)
  coefficient <- sub("^0+", "", paste0(parts[[2L]], fractional))
  if (!nzchar(coefficient)) {
    return(list(coefficient = "0", exponent = 0))
  }
  trailing <- nchar(coefficient, type = "bytes") -
    nchar(sub("0+$", "", coefficient), type = "bytes")
  if (trailing > 0L) {
    coefficient <- substr(
      coefficient, 1L, nchar(coefficient, type = "bytes") - trailing)
  }
  list(
    coefficient = coefficient,
    exponent = exponent - nchar(fractional, type = "bytes") + trailing)
}

.dsvert_joint_dp_lifetime_decimal_compose <- function(
    coefficient, exponent) {
  coefficient <- if (inherits(coefficient, "bignum")) {
    as.character(coefficient)
  } else coefficient
  if (!is.character(coefficient) || length(coefficient) != 1L ||
      !grepl("^[0-9]+$", coefficient) || !is.numeric(exponent) ||
      length(exponent) != 1L || is.na(exponent) || !is.finite(exponent) ||
      exponent != floor(exponent)) return(NULL)
  coefficient <- sub("^0+", "", coefficient)
  if (!nzchar(coefficient)) return("0")
  trailing <- nchar(coefficient, type = "bytes") -
    nchar(sub("0+$", "", coefficient), type = "bytes")
  if (trailing > 0L) {
    coefficient <- substr(
      coefficient, 1L, nchar(coefficient, type = "bytes") - trailing)
  }
  exponent <- exponent + trailing
  scientific_exponent <- exponent + nchar(coefficient, type = "bytes") - 1
  mantissa <- if (nchar(coefficient, type = "bytes") == 1L) {
    coefficient
  } else paste0(
    substr(coefficient, 1L, 1L), ".",
    substr(coefficient, 2L, nchar(coefficient, type = "bytes")))
  if (scientific_exponent == 0) return(mantissa)
  paste0(mantissa, "e", if (scientific_exponent > 0) "+" else "-",
         sprintf("%.0f", abs(scientific_exponent)))
}

.dsvert_joint_dp_lifetime_exact_total <- function(value, count) {
  if (!.dsvert_dp_is_number(value, 0) ||
      !.dsvert_dp_is_integer(count, 1, 2^53 - 1)) return(NULL)
  text <- sub(
    getOption("OutDec", "."), ".",
    format(as.numeric(value), digits = 17L,
           scientific = TRUE, trim = TRUE),
    fixed = TRUE)
  parts <- .dsvert_joint_dp_lifetime_decimal_parts(text)
  if (is.null(parts)) return(NULL)
  coefficient <- openssl::bignum(parts$coefficient) *
    openssl::bignum(sprintf("%.0f", count))
  .dsvert_joint_dp_lifetime_decimal_compose(
    coefficient, parts$exponent)
}

.dsvert_joint_dp_lifetime_exact_compare <- function(left, right) {
  parts <- lapply(list(left, right),
                  .dsvert_joint_dp_lifetime_decimal_parts)
  if (any(vapply(parts, is.null, logical(1L)))) return(NA_integer_)
  exponent <- min(vapply(parts, `[[`, numeric(1L), "exponent"))
  values <- lapply(parts, function(part) {
    value <- openssl::bignum(part$coefficient)
    shift <- part$exponent - exponent
    if (shift > 0) value <- value * openssl::bignum(10)^as.integer(shift)
    value
  })
  if (values[[1L]] < values[[2L]]) return(-1L)
  if (values[[1L]] > values[[2L]]) return(1L)
  0L
}

.dsvert_joint_dp_validate_release_domain <- function(value) {
  expected <- c(
    "version", "generation", "domain_id", "rotation_count",
    "automatic_generation", "automatic_rotation", "snapshot_derived",
    "key_material_exposed")
  .dsvert_dp_has_exact_names(value, expected) &&
    identical(value$version, "dsvert-joint-dp-release-domain-v1") &&
    .dsvert_dp_is_integer(value$generation, 1, 2^53 - 1) &&
    .dsvert_dp_is_string(value$domain_id) &&
    grepl("^rd_[0-9a-f]{64}$", value$domain_id) &&
    .dsvert_dp_is_integer(value$rotation_count, 0, 2^53 - 1) &&
    identical(as.numeric(value$rotation_count),
              as.numeric(value$generation) - 1) &&
    identical(value$automatic_generation, TRUE) &&
    identical(value$automatic_rotation, TRUE) &&
    identical(value$snapshot_derived, FALSE) &&
    identical(value$key_material_exposed, FALSE)
}

.dsvert_joint_dp_validate_capsule_policy <- function(policy, server) {
  invalid <- function() {
    stop("Server '", server,
         "' returned an invalid reusable joint-DP capsule policy",
         call. = FALSE)
  }
  expected <- c(
    "contract", "domain", "cohort_id", "peer_name", "own_identity_pk",
    "peer_pinset", "peer_pinset_sha256", "peer_count",
    "designated_noise_peers", "capsule_epsilon", "capsule_delta",
    "lifetime_max_distinct_capsules", "lifetime_epsilon_upper_bound",
    "lifetime_delta_upper_bound",
    "adjacency", "patient_column", "unit_capacity",
    "max_records_per_unit", "overflow_policy", "sampler")
  expected_lifetime_epsilon <- .dsvert_joint_dp_lifetime_exact_total(
    policy$capsule_epsilon,
    policy$lifetime_max_distinct_capsules)
  expected_lifetime_delta <- .dsvert_joint_dp_lifetime_exact_total(
    policy$capsule_delta,
    policy$lifetime_max_distinct_capsules)
  epsilon_comparison <- .dsvert_joint_dp_lifetime_exact_compare(
    policy$lifetime_epsilon_upper_bound, "8")
  delta_comparison <- .dsvert_joint_dp_lifetime_exact_compare(
    policy$lifetime_delta_upper_bound, "1")
  if (!.dsvert_dp_has_exact_names(policy, expected) ||
      !identical(policy$contract, "immutable_reusable_capsule_v1") ||
      !.dsvert_dp_is_string(policy$domain) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", policy$domain) ||
      !.dsvert_dp_is_string(policy$cohort_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", policy$cohort_id) ||
      !.dsvert_dp_is_string(policy$peer_name) ||
      !.dsvert_dp_validate_pinset(policy$peer_pinset) ||
      !policy$peer_name %in% names(policy$peer_pinset) ||
      !.dsvert_dp_is_string(policy$own_identity_pk) ||
      !identical(policy$own_identity_pk,
                 unname(policy$peer_pinset[[policy$peer_name]])) ||
      !.dsvert_dp_is_string(policy$peer_pinset_sha256) ||
      !grepl("^[0-9a-f]{64}$", policy$peer_pinset_sha256) ||
      !identical(.dsvert_dp_pinset_hash(policy$peer_pinset),
                 policy$peer_pinset_sha256) ||
      !.dsvert_dp_is_integer(policy$peer_count, 2,
                            .Machine$integer.max) ||
      policy$peer_count != length(policy$peer_pinset) ||
      !is.character(policy$designated_noise_peers) ||
      length(policy$designated_noise_peers) != 2L ||
      anyNA(policy$designated_noise_peers) ||
      anyDuplicated(policy$designated_noise_peers) ||
      !identical(policy$designated_noise_peers,
                 sort(policy$designated_noise_peers, method = "radix")) ||
      !all(policy$designated_noise_peers %in% names(policy$peer_pinset)) ||
      !.dsvert_dp_is_number(
        policy$capsule_epsilon, 0, .DSVERT_DP_MAXIMUM_EPSILON, TRUE) ||
      !.dsvert_dp_is_number(policy$capsule_delta, 0, 1) ||
      policy$capsule_delta >= 1 ||
      !.dsvert_dp_is_integer(
        policy$lifetime_max_distinct_capsules, 1, 2^53 - 1) ||
      !identical(
        policy$lifetime_epsilon_upper_bound,
        expected_lifetime_epsilon) ||
      !identical(
        policy$lifetime_delta_upper_bound,
        expected_lifetime_delta) ||
      !isTRUE(epsilon_comparison <= 0L) ||
      !isTRUE(delta_comparison < 0L) ||
      !policy$adjacency %in% c(
        "add_remove_patient", "replace_one_fixed_cohort") ||
      !.dsvert_dp_is_string(policy$patient_column) ||
      !.dsvert_dp_is_integer(
        policy$unit_capacity, 1, .DSVERT_DP_MAX_COORDINATES) ||
      !.dsvert_dp_is_integer(
        policy$max_records_per_unit, 1, .DSVERT_DP_MAX_COORDINATES) ||
      !identical(policy$overflow_policy, "reject_snapshot") ||
      !identical(policy$sampler,
                 "two_private_hmac_seeds_one_gc_sample_v1")) {
    invalid()
  }
  invisible(TRUE)
}

.dsvert_joint_dp_validate_capsule_privacy_contract <- function(value) {
  expected <- c(
    "definition", "scope", "adversary_model", "assumptions",
    "simultaneous_designated_history_rollback_protection",
    "transcript_security", "malicious_security", "operation_accounting",
    "privacy_budget_gate", "operation_limit", "request_limit",
    "history_can_deny_operation",
    "release_instance_accounting", "accuracy_depends_on_request_history",
    "reuse", "new_capsules", "lifetime_max_distinct_capsules",
    "lifetime_epsilon_upper_bound", "lifetime_delta_upper_bound")
  .dsvert_dp_has_exact_names(value, expected) &&
    identical(value$definition,
              "bounded_lifetime_epsilon_delta_dp") &&
    identical(value$scope,
              paste0(
                "at_most_N_immutable_snapshot_workload_capsules_per_",
                "stable_privacy_accountant_namespace")) &&
    identical(value$adversary_model,
              "authenticated_semi_honest_noncollusion") &&
    identical(
      value$assumptions,
      paste0("declared_adjacency_bounds_immutable_snapshot_protocol_",
             "compliant_peers_at_least_one_noncolluding_designated_",
             "noise_peer_retains_and_uses_complete_authenticated_monotonic_",
             "history_stable_unique_privacy_accountant_namespace_per_",
             "protected_privacy_",
             "universe")) &&
    identical(
      value$simultaneous_designated_history_rollback_protection,
      "not_claimed_without_external_linearizable_cas") &&
    identical(value$transcript_security,
              "computational_mpc_and_csprng") &&
    identical(value$malicious_security, FALSE) &&
    identical(value$operation_accounting,
              "one_per_distinct_capsule_allocator_commit") &&
    identical(value$privacy_budget_gate, TRUE) &&
    identical(value$operation_limit, TRUE) &&
    identical(value$request_limit, FALSE) &&
    identical(value$history_can_deny_operation, TRUE) &&
    identical(value$release_instance_accounting,
              "one_public_release_instance_per_capsule_id") &&
    identical(value$accuracy_depends_on_request_history, FALSE) &&
    identical(value$reuse, "unlimited_sticky_postprocessing") &&
    identical(value$new_capsules,
              "allowed_until_authenticated_lifetime_bound") &&
    .dsvert_dp_is_integer(
      value$lifetime_max_distinct_capsules, 1, 2^53 - 1) &&
    .dsvert_dp_is_string(value$lifetime_epsilon_upper_bound) &&
    .dsvert_dp_is_string(value$lifetime_delta_upper_bound)
}

.dsvert_joint_dp_validate_release_instance_telemetry <- function(
    value, policy, server) {
  invalid <- function() {
    stop("Server '", server,
         "' returned invalid release-instance composition telemetry",
         call. = FALSE)
  }
  expected <- c(
    "releases_published", "release_epsilon", "release_delta",
    "lifetime_max_distinct_capsules", "remaining_distinct_capsules",
    "cumulative_epsilon_upper_bound", "cumulative_delta_upper_bound",
    "cumulative_delta_vacuous", "composition_role",
    "release_accounting", "replay_accounting", "rotation_accounting",
    "privacy_budget_gate", "operation_limit", "request_limit",
    "history_can_deny_operation", "admission_role")
  if (!.dsvert_dp_has_exact_names(value, expected) ||
      !.dsvert_dp_is_integer(value$releases_published) ||
      !identical(as.numeric(value$lifetime_max_distinct_capsules),
                 as.numeric(policy$lifetime_max_distinct_capsules)) ||
      !.dsvert_dp_is_integer(value$remaining_distinct_capsules, 0,
                            policy$lifetime_max_distinct_capsules) ||
      value$releases_published > policy$lifetime_max_distinct_capsules ||
      value$remaining_distinct_capsules >
        policy$lifetime_max_distinct_capsules - value$releases_published ||
      !.dsvert_dp_num_equal(value$release_epsilon,
                            policy$capsule_epsilon) ||
      !.dsvert_dp_num_equal(value$release_delta,
                            policy$capsule_delta) ||
      !.dsvert_dp_is_number(value$cumulative_epsilon_upper_bound, 0) ||
      !.dsvert_dp_is_number(value$cumulative_delta_upper_bound, 0) ||
      !.dsvert_dp_num_equal(
        value$cumulative_epsilon_upper_bound,
        value$releases_published * value$release_epsilon, 2048) ||
      !.dsvert_dp_num_equal(
        value$cumulative_delta_upper_bound,
        value$releases_published * value$release_delta, 2048) ||
      !is.logical(value$cumulative_delta_vacuous) ||
      length(value$cumulative_delta_vacuous) != 1L ||
      is.na(value$cumulative_delta_vacuous) ||
      !identical(value$cumulative_delta_vacuous, FALSE) ||
      value$cumulative_delta_upper_bound >= 1 ||
      !identical(value$composition_role,
                 "basic_composition_authenticated_lifetime_bound") ||
      !identical(value$release_accounting,
                 "one_public_release_instance_per_capsule_id") ||
      !identical(value$replay_accounting,
                 "same_instance_replay_not_recounted") ||
      !identical(value$rotation_accounting,
                 paste0("prepublication_rotation_only_postpublication_",
                        "replay_or_fail_closed")) ||
      !identical(value$privacy_budget_gate, TRUE) ||
      !identical(value$operation_limit, TRUE) ||
      !identical(value$request_limit, FALSE) ||
      !identical(value$history_can_deny_operation, TRUE) ||
      !identical(value$admission_role,
                 "authenticated_lifetime_gate_before_sampler")) {
    invalid()
  }
  invisible(TRUE)
}

.dsvert_joint_dp_validate_capsule_telemetry <- function(
    value, policy, server) {
  invalid <- function() {
    stop("Server '", server,
         "' returned invalid capsule composition telemetry", call. = FALSE)
  }
  expected <- c(
    "capsules_created", "capsule_epsilon", "capsule_delta",
    "lifetime_max_distinct_capsules", "remaining_distinct_capsules",
    "cumulative_epsilon_upper_bound", "cumulative_delta_upper_bound",
    "cumulative_delta_vacuous", "composition_role",
    "registration_policy", "privacy_budget_gate", "operation_limit",
    "request_limit", "history_can_deny_new_release", "admission_role")
  if (!.dsvert_dp_has_exact_names(value, expected) ||
      !.dsvert_dp_is_integer(value$capsules_created) ||
      !identical(as.numeric(value$lifetime_max_distinct_capsules),
                 as.numeric(policy$lifetime_max_distinct_capsules)) ||
      !.dsvert_dp_is_integer(value$remaining_distinct_capsules, 0,
                            policy$lifetime_max_distinct_capsules) ||
      value$capsules_created > policy$lifetime_max_distinct_capsules ||
      !identical(
        as.numeric(value$remaining_distinct_capsules),
        as.numeric(policy$lifetime_max_distinct_capsules -
          value$capsules_created)) ||
      !.dsvert_dp_num_equal(value$capsule_epsilon,
                            policy$capsule_epsilon) ||
      !.dsvert_dp_num_equal(value$capsule_delta, policy$capsule_delta) ||
      !.dsvert_dp_is_number(value$cumulative_epsilon_upper_bound, 0) ||
      !.dsvert_dp_is_number(value$cumulative_delta_upper_bound, 0) ||
      !.dsvert_dp_num_equal(
        value$cumulative_epsilon_upper_bound,
        value$capsules_created * policy$capsule_epsilon, 2048) ||
      !.dsvert_dp_num_equal(
        value$cumulative_delta_upper_bound,
        value$capsules_created * policy$capsule_delta, 2048) ||
      !is.logical(value$cumulative_delta_vacuous) ||
      length(value$cumulative_delta_vacuous) != 1L ||
      is.na(value$cumulative_delta_vacuous) ||
      !identical(value$cumulative_delta_vacuous, FALSE) ||
      value$cumulative_delta_upper_bound >= 1 ||
      !identical(
        value$composition_role,
        "basic_composition_authenticated_lifetime_bound") ||
      !identical(
        value$registration_policy,
        "allocator_admitted_distinct_capsules_up_to_lifetime_limit") ||
      !identical(value$privacy_budget_gate, TRUE) ||
      !identical(value$operation_limit, TRUE) ||
      !identical(value$request_limit, FALSE) ||
      !identical(value$history_can_deny_new_release, TRUE) ||
      !identical(value$admission_role,
                 "allocator_reservation_before_protected_access")) {
    invalid()
  }
  invisible(TRUE)
}

.dsvert_joint_dp_validate_capsule_status <- function(status, server) {
  invalid <- function() {
    stop("Server '", server,
         "' does not expose the required reusable joint-DP capsule contract",
         call. = FALSE)
  }
  expected <- c(
    "version", "enabled", "privacy_contract", "policy", "noise_root",
    "release_domain", "role", "composition_telemetry",
    "release_instance_telemetry")
  if (!.dsvert_dp_has_exact_names(status, expected) ||
      !identical(status$version,
                 .DSVERT_CLIENT_JOINT_DP_CAPSULE_STATUS_VERSION) ||
      !identical(status$enabled, TRUE) ||
      !.dsvert_joint_dp_validate_capsule_privacy_contract(
        status$privacy_contract) ||
      !.dsvert_dp_validate_noise_root(status$noise_root) ||
      !.dsvert_dp_has_exact_names(
        status$role, c("designated_noise_peer", "allocator")) ||
      !is.logical(status$role$designated_noise_peer) ||
      length(status$role$designated_noise_peer) != 1L ||
      is.na(status$role$designated_noise_peer)) {
    invalid()
  }
  .dsvert_joint_dp_validate_capsule_policy(status$policy, server)
  if (!identical(
        as.numeric(status$privacy_contract$lifetime_max_distinct_capsules),
        as.numeric(status$policy$lifetime_max_distinct_capsules)) ||
      !identical(status$privacy_contract$lifetime_epsilon_upper_bound,
                 status$policy$lifetime_epsilon_upper_bound) ||
      !identical(status$privacy_contract$lifetime_delta_upper_bound,
                 status$policy$lifetime_delta_upper_bound)) {
    invalid()
  }
  designated <- server %in% status$policy$designated_noise_peers
  if (!identical(status$role$designated_noise_peer, designated) ||
      (designated && !identical(
        status$role$allocator, "authenticated_ready")) ||
      (!designated && !identical(
        status$role$allocator, "not_applicable_policy_attestor"))) {
    invalid()
  }
  if (designated) {
    if (!.dsvert_joint_dp_validate_release_domain(
          status$release_domain)) invalid()
    .dsvert_joint_dp_validate_capsule_telemetry(
      status$composition_telemetry, status$policy, server)
    .dsvert_joint_dp_validate_release_instance_telemetry(
      status$release_instance_telemetry, status$policy, server)
    if (status$release_instance_telemetry$releases_published >
          status$composition_telemetry$capsules_created ||
        !identical(
          as.numeric(
            status$release_instance_telemetry$remaining_distinct_capsules),
          as.numeric(
            status$composition_telemetry$remaining_distinct_capsules))) {
      invalid()
    }
  } else if (!is.null(status$release_domain) ||
             !is.null(status$composition_telemetry) ||
             !is.null(status$release_instance_telemetry)) {
    invalid()
  }
  invisible(TRUE)
}

.dsvert_joint_dp_capsule_status_consensus <- function(result, datasources) {
  servers <- sort(names(datasources), method = "radix")
  if (!is.list(result) || is.null(names(result)) || anyNA(names(result)) ||
      anyDuplicated(names(result)) || !setequal(names(result), servers)) {
    stop("The capsule status handshake returned an invalid server set",
         call. = FALSE)
  }
  result <- result[servers]
  for (server in servers) {
    .dsvert_joint_dp_validate_capsule_status(result[[server]], server)
    policy <- result[[server]]$policy
    if (!identical(policy$peer_name, server) ||
        !identical(names(policy$peer_pinset), servers)) {
      stop("Server '", server,
           "' does not match its pinned capsule identity", call. = FALSE)
    }
  }

  # Diagnose a regenerated identity before collapsing the disagreement into a
  # generic policy mismatch. Every other server's name-bound map supplies the
  # expected public pin; the affected server supplies only its observed public
  # identity. The relay is never allowed to update either value.
  for (server in servers) {
    observed_pk <- result[[server]]$policy$own_identity_pk
    observers <- setdiff(servers, server)
    configured <- vapply(observers, function(observer) {
      server %in% names(result[[observer]]$policy$peer_pinset)
    }, logical(1L))
    expected <- unique(vapply(observers[configured], function(observer) {
      unname(result[[observer]]$policy$peer_pinset[[server]])
    }, character(1L)))
    expected_fingerprint <- if (!all(configured) || !length(expected)) {
      "unconfigured"
    } else if (length(expected) > 1L) {
      "inconsistent"
    } else {
      .dsvert_client_identity_fingerprint(expected[[1L]])
    }
    if (!identical(expected_fingerprint,
                   .dsvert_client_identity_fingerprint(observed_pk))) {
      .dsvert_client_stop_peer_not_recognized(
        server,
        .dsvert_client_identity_fingerprint(observed_pk),
        expected_fingerprint)
    }
  }

  reference <- result[[1L]]
  common_policy <- setdiff(
    names(reference$policy), c("peer_name", "own_identity_pk"))
  for (server in servers) {
    current <- result[[server]]
    if (!identical(current$policy[common_policy],
                   reference$policy[common_policy])) {
      stop("The connected servers disagree on the reusable capsule policy",
           call. = FALSE)
    }
  }
  designated <- reference$policy$designated_noise_peers
  if (!identical(
        unname(vapply(result, function(value) {
          isTRUE(value$role$designated_noise_peer)
        }, logical(1L))), servers %in% designated)) {
    stop("The connected servers disagree on the designated noise peers",
         call. = FALSE)
  }
  if (!identical(
        result[[designated[[1L]]]]$composition_telemetry,
        result[[designated[[2L]]]]$composition_telemetry)) {
    stop("The designated peers disagree on capsule composition telemetry",
         call. = FALSE)
  }
  if (!identical(
        result[[designated[[1L]]]]$release_instance_telemetry,
        result[[designated[[2L]]]]$release_instance_telemetry)) {
    stop("The designated peers disagree on release-instance telemetry",
         call. = FALSE)
  }
  class(result) <- c(
    "ds.vertDPStatus", "ds.vertJointDPCapsuleStatus", "list")
  result
}

.dsvert_joint_dp_capsule_status_impl <- function(
    datasources = NULL, .aggregate = DSI::datashield.aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  result <- .dsvert_aggregate_strict(
    conns = datasources,
    expr = call(name = "dsvertJointDPCapsuleStatusDS"),
    operation = "reusable joint-DP capsule status",
    .aggregate = .aggregate)
  .dsvert_joint_dp_capsule_status_consensus(result, datasources)
}
