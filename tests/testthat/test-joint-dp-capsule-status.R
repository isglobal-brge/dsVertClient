.dsvert_joint_dp_capsule_status_impl <- function(...) {
  impl <- get(
    ".dsvert_joint_dp_capsule_status_impl",
    envir = asNamespace("dsVertClient"), inherits = FALSE)
  testthat::with_mocked_bindings(
    impl(...),
    .dsvert_block_retired_remote_route = function(...) invisible(FALSE),
    .package = "dsVertClient")
}

.capsule_status_pk <- function(offset) {
  chartr("+/", "-_", sub(
    "=+$", "",
    jsonlite::base64_enc(as.raw((seq_len(32L) + offset) %% 256L)),
    perl = TRUE))
}

.capsule_status_fixture <- function(
    servers = c("site_a", "site_b", "site_c"),
    designated = c("site_a", "site_c"), capsules = 7,
    lifetime_max = 8, capsule_epsilon = 1, capsule_delta = 2^-100,
    releases = 2) {
  pins <- stats::setNames(vapply(seq_along(servers), function(index) {
    .capsule_status_pk((index - 1L) * 40L)
  }, character(1L)), servers)
  pin_hash <- .dsvert_dp_pinset_hash(pins)
  lifetime_epsilon <- .dsvert_joint_dp_lifetime_exact_total(
    capsule_epsilon, lifetime_max)
  lifetime_delta <- .dsvert_joint_dp_lifetime_exact_total(
    capsule_delta, lifetime_max)
  values <- lapply(servers, function(server) {
    is_designated <- server %in% designated
    telemetry <- if (is_designated) list(
      capsules_created = as.numeric(capsules),
      lifetime_max_distinct_capsules = as.numeric(lifetime_max),
      remaining_distinct_capsules = as.numeric(lifetime_max - capsules),
      capsule_epsilon = capsule_epsilon,
      capsule_delta = capsule_delta,
      cumulative_epsilon_upper_bound = capsules * capsule_epsilon,
      cumulative_delta_upper_bound = capsules * capsule_delta,
      cumulative_delta_vacuous = FALSE,
      composition_role = "basic_composition_authenticated_lifetime_bound",
      registration_policy =
        "allocator_admitted_distinct_capsules_up_to_lifetime_limit",
      privacy_budget_gate = TRUE, operation_limit = TRUE,
      request_limit = FALSE, history_can_deny_new_release = TRUE,
      admission_role = "allocator_reservation_before_protected_access") else NULL
    release_telemetry <- if (is_designated) list(
      releases_published = as.numeric(releases),
      lifetime_max_distinct_capsules = as.numeric(lifetime_max),
      remaining_distinct_capsules = as.numeric(lifetime_max - capsules),
      release_epsilon = capsule_epsilon,
      release_delta = capsule_delta,
      cumulative_epsilon_upper_bound = releases * capsule_epsilon,
      cumulative_delta_upper_bound = releases * capsule_delta,
      cumulative_delta_vacuous = FALSE,
      composition_role = "basic_composition_authenticated_lifetime_bound",
      release_accounting = "one_public_release_instance_per_capsule_id",
      replay_accounting = "same_instance_replay_not_recounted",
      rotation_accounting = paste0(
        "prepublication_rotation_only_postpublication_replay_or_",
        "fail_closed"),
      privacy_budget_gate = TRUE, operation_limit = TRUE,
      request_limit = FALSE, history_can_deny_operation = TRUE,
      admission_role = "authenticated_lifetime_gate_before_sampler") else NULL
    list(
      version = "dsvert-joint-dp-capsule-status-v5",
      enabled = TRUE,
      privacy_contract = list(
        definition = "bounded_lifetime_epsilon_delta_dp",
        scope = paste0(
          "at_most_N_immutable_snapshot_workload_capsules_per_stable_",
          "privacy_accountant_namespace"),
        adversary_model = "authenticated_semi_honest_noncollusion",
        assumptions = paste0(
          "declared_adjacency_bounds_immutable_snapshot_protocol_",
          "compliant_peers_at_least_one_noncolluding_designated_noise_peer_",
          "retains_and_uses_complete_authenticated_monotonic_history_",
          "stable_unique_privacy_accountant_namespace_per_protected_",
          "privacy_universe"),
        simultaneous_designated_history_rollback_protection =
          "not_claimed_without_external_linearizable_cas",
        transcript_security = "computational_mpc_and_csprng",
        malicious_security = FALSE,
        operation_accounting =
          "one_per_distinct_capsule_allocator_commit",
        privacy_budget_gate = TRUE,
        operation_limit = TRUE,
        request_limit = FALSE,
        history_can_deny_operation = TRUE,
        release_instance_accounting =
          "one_public_release_instance_per_capsule_id",
        accuracy_depends_on_request_history = FALSE,
        reuse = "unlimited_sticky_postprocessing",
        new_capsules = "allowed_until_authenticated_lifetime_bound",
        lifetime_max_distinct_capsules = as.numeric(lifetime_max),
        lifetime_epsilon_upper_bound = lifetime_epsilon,
        lifetime_delta_upper_bound = lifetime_delta),
      policy = list(
        contract = "immutable_reusable_capsule_v1",
        domain = "capsule-study", cohort_id = "capsule-cohort",
        peer_name = server, own_identity_pk = unname(pins[[server]]),
        peer_pinset = pins, peer_pinset_sha256 = pin_hash,
        peer_count = as.integer(length(servers)),
        designated_noise_peers = designated,
        capsule_epsilon = capsule_epsilon, capsule_delta = capsule_delta,
        lifetime_max_distinct_capsules = as.numeric(lifetime_max),
        lifetime_epsilon_upper_bound = lifetime_epsilon,
        lifetime_delta_upper_bound = lifetime_delta,
        adjacency = "add_remove_patient", patient_column = "patient_id",
        unit_capacity = 1000, max_records_per_unit = 2,
        overflow_policy = "reject_snapshot",
        sampler = "two_private_hmac_seeds_one_gc_sample_v1"),
      noise_root = list(
        protocol = "dsvert-dp-noise-root-v1",
        provider_id = paste0("provider-", server),
        key_id = paste0("key-", server), privacy_epoch = 1,
        external = FALSE, storage = "owner_only_file",
        automatic_generation = TRUE, automatic_recovery = TRUE,
        automatic_rotation = FALSE, rotation_count = 0,
        key_material_exposed = FALSE),
      release_domain = if (is_designated) list(
        version = "dsvert-joint-dp-release-domain-v1",
        generation = 1,
        domain_id = paste0("rd_", digest::digest(
          paste0("capsule-status/", server),
          "sha256", serialize = FALSE)),
        rotation_count = 0,
        automatic_generation = TRUE, automatic_rotation = TRUE,
        snapshot_derived = FALSE, key_material_exposed = FALSE) else NULL,
      role = list(
        designated_noise_peer = is_designated,
        allocator = if (is_designated) "authenticated_ready" else
          "not_applicable_policy_attestor"),
      composition_telemetry = telemetry,
      release_instance_telemetry = release_telemetry)
  })
  stats::setNames(values, servers)
}

.capsule_status_aggregate <- function(values, calls = NULL) {
  force(values)
  function(conns, expr, error = NULL, errors.print = TRUE, ...) {
    expect_false(errors.print)
    expect_identical(as.character(expr[[1L]]),
                     "dsvertJointDPCapsuleStatusDS")
    if (!is.null(calls)) calls[[length(calls) + 1L]] <- expr
    values[names(conns)]
  }
}

test_that("K=2 through K=5 handshakes expose a lifetime gate, not a request quota", {
  for (k in 2:5) {
    servers <- paste0("site_", letters[seq_len(k)])
    designated <- if (k == 2L) servers else c("site_a", "site_c")
    values <- .capsule_status_fixture(
      servers = servers, designated = designated)
    conns <- lapply(
      names(values), function(value) structure(1, class = "fake"))
    names(conns) <- names(values)
    calls <- list()
    status <- .dsvert_joint_dp_capsule_status_impl(
      conns, .capsule_status_aggregate(values, calls))

    expect_s3_class(status, "ds.vertJointDPCapsuleStatus")
    expect_identical(names(status), servers)
    expect_true(all(vapply(status[designated], function(value) {
      isTRUE(value$role$designated_noise_peer)
    }, logical(1L))))
    for (server in setdiff(servers, designated)) {
      expect_false(status[[server]]$role$designated_noise_peer)
      expect_null(status[[server]]$composition_telemetry)
      expect_null(status[[server]]$release_instance_telemetry)
    }
    expect_identical(status[[designated[[1L]]]]$composition_telemetry,
                     status[[designated[[2L]]]]$composition_telemetry)
    for (value in status) {
      expect_identical(
        value$privacy_contract$operation_accounting,
        "one_per_distinct_capsule_allocator_commit")
      expect_true(value$privacy_contract$privacy_budget_gate)
      expect_true(value$privacy_contract$operation_limit)
      expect_false(value$privacy_contract$request_limit)
      expect_true(value$privacy_contract$history_can_deny_operation)
      expect_identical(
        value$privacy_contract$
          simultaneous_designated_history_rollback_protection,
        "not_claimed_without_external_linearizable_cas")
      expect_identical(
        value$privacy_contract$release_instance_accounting,
        "one_public_release_instance_per_capsule_id")
      expect_false(
        value$privacy_contract$accuracy_depends_on_request_history)
      expect_identical(value$privacy_contract$reuse,
                       "unlimited_sticky_postprocessing")
      expect_identical(
        as.numeric(value$privacy_contract$lifetime_max_distinct_capsules), 8)
    }
    expect_identical(
      status[[designated[[1L]]]]$composition_telemetry$
        remaining_distinct_capsules, 1)
    expect_identical(
      status[[designated[[1L]]]]$release_instance_telemetry$
        remaining_distinct_capsules, 1)
    encoded <- jsonlite::toJSON(status, auto_unbox = TRUE, null = "null")
    expect_false(any(vapply(
      c("decay", "allocation_slots", "queries_remaining", "request_quota"),
      grepl, logical(1L), x = encoded, fixed = TRUE)))
  }
})

test_that("relay connection order cannot choose the two designated noise peers", {
  servers <- paste0("site_", letters[1:5])
  designated <- c("site_a", "site_c")
  values <- .capsule_status_fixture(
    servers = servers, designated = designated)
  relay_order <- c("site_e", "site_b", "site_d", "site_c", "site_a")
  conns <- stats::setNames(lapply(relay_order, function(value) {
    structure(1, class = "fake")
  }), relay_order)

  status <- .dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(values))

  expect_identical(names(status), servers)
  expect_identical(
    status[[relay_order[[1L]]]]$policy$designated_noise_peers,
    designated)
  selected <- names(Filter(
    function(value) isTRUE(value$role$designated_noise_peer), status))
  expect_length(selected, 2L)
  expect_setequal(selected, designated)
  expect_true(all(vapply(status[setdiff(servers, designated)],
                         function(value) is.null(value$composition_telemetry),
                         logical(1L))))
})

test_that("a regenerated peer identity gives actionable repinning guidance", {
  values <- .capsule_status_fixture()
  conns <- lapply(names(values), function(value) structure(1, class = "fake"))
  names(conns) <- names(values)
  replacement <- .capsule_status_pk(200L)

  values$site_b$policy$own_identity_pk <- replacement
  values$site_b$policy$peer_pinset[["site_b"]] <- replacement
  values$site_b$policy$peer_pinset_sha256 <-
    .dsvert_dp_pinset_hash(values$site_b$policy$peer_pinset)

  unrecognized <- tryCatch(
    .dsvert_joint_dp_capsule_status_impl(
      conns, .capsule_status_aggregate(values)),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(unrecognized, "dsvert_peer_not_recognized")
  expect_identical(unrecognized$code, "peer_not_recognized")
  expect_identical(unrecognized$peer_name, "site_b")
  expect_match(unrecognized$expected_fingerprint_sha256, "^[0-9a-f]{64}$")
  expect_match(unrecognized$observed_fingerprint_sha256, "^[0-9a-f]{64}$")
  expect_match(conditionMessage(unrecognized), "ds.getIdentityPks", fixed = TRUE)
  expect_match(conditionMessage(unrecognized), "verify the observed fingerprint")
  expect_match(conditionMessage(unrecognized), "dsvert.trusted_peers", fixed = TRUE)
  expect_match(conditionMessage(unrecognized), "each other participating server")
  expect_match(conditionMessage(unrecognized), "must not pin its own identity")
  expect_false(grepl(
    "identity_sk|private key|identity.seed|noise_root",
    conditionMessage(unrecognized), ignore.case = TRUE))

  # Repinning changes no lifetime allowance. Once all custodians update the
  # verified name-bound pin, consensus resumes.
  for (server in names(values)) {
    values[[server]]$policy$peer_pinset[["site_b"]] <- replacement
    values[[server]]$policy$peer_pinset_sha256 <-
      .dsvert_dp_pinset_hash(values[[server]]$policy$peer_pinset)
  }
  expect_s3_class(
    .dsvert_joint_dp_capsule_status_impl(
      conns, .capsule_status_aggregate(values)),
    "ds.vertJointDPCapsuleStatus")
})

test_that("partial repinning remains a typed unrecognized-peer error", {
  values <- .capsule_status_fixture()
  conns <- lapply(names(values), function(value) structure(1, class = "fake"))
  names(conns) <- names(values)
  replacement <- .capsule_status_pk(201L)

  values$site_b$policy$own_identity_pk <- replacement
  values$site_b$policy$peer_pinset[["site_b"]] <- replacement
  values$site_b$policy$peer_pinset_sha256 <-
    .dsvert_dp_pinset_hash(values$site_b$policy$peer_pinset)
  values$site_c$policy$peer_pinset[["site_b"]] <- replacement
  values$site_c$policy$peer_pinset_sha256 <-
    .dsvert_dp_pinset_hash(values$site_c$policy$peer_pinset)

  condition <- tryCatch(
    .dsvert_joint_dp_capsule_status_impl(
      conns, .capsule_status_aggregate(values)),
    dsvert_peer_not_recognized = identity)
  expect_s3_class(condition, "dsvert_peer_not_recognized")
  expect_identical(condition$peer_name, "site_b")
  expect_identical(
    condition$expected_fingerprint_sha256, "inconsistent")
  expect_match(conditionMessage(condition), "inconsistent name-bound")
})

test_that("vacuous lifetime policies fail before becoming client admission", {
  values <- .capsule_status_fixture(lifetime_max = 8, capsule_delta = 0.5)
  for (server in names(values)) {
    expect_identical(values[[server]]$policy$lifetime_delta_upper_bound, "4")
  }
  conns <- lapply(names(values), function(value) structure(1, class = "fake"))
  names(conns) <- names(values)

  expect_error(
    .dsvert_joint_dp_capsule_status_impl(
      conns, .capsule_status_aggregate(values)),
    "invalid reusable")
})

test_that("DP status print summarizes the non-vacuous lifetime gate", {
  values <- .capsule_status_fixture(capsules = 7)
  conns <- lapply(names(values), function(value) structure(1, class = "fake"))
  names(conns) <- names(values)
  status <- .dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(values))

  output <- capture.output(returned <- print(status))

  expect_identical(returned, status)
  expect_true(any(grepl("peers: 3", output, fixed = TRUE)))
  expect_true(any(grepl(
    "designated noise peers: site_a, site_c", output, fixed = TRUE)))
  expect_true(any(grepl(
    "allocator-committed reservation units: 7e+00", output, fixed = TRUE)))
  expect_true(any(grepl(
    "published release instances: 2e+00", output, fixed = TRUE)))
  expect_false(any(grepl("WARNING", output, fixed = TRUE)))
  expect_true(any(grepl(
    "allocator-committed reservation units remaining: 1e+00",
    output, fixed = TRUE)))
  expect_true(any(grepl(
    "request quota: none", output, fixed = TRUE)))
})

test_that("lifetime status uses exact decimal composition boundaries", {
  conns <- stats::setNames(lapply(1:3, function(value) {
    structure(value, class = "fake")
  }), c("site_a", "site_b", "site_c"))

  valid <- .capsule_status_fixture(
    capsules = 3, releases = 3, lifetime_max = 3,
    capsule_epsilon = 8 / 3)
  expect_identical(
    valid$site_a$policy$lifetime_epsilon_upper_bound,
    "7.9999999999999995")
  expect_s3_class(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(valid)),
    "ds.vertJointDPCapsuleStatus")

  rounded_over <- .capsule_status_fixture(
    capsules = 1, releases = 1, lifetime_max = 10,
    capsule_epsilon = 0.8)
  expect_identical(
    rounded_over$site_a$policy$lifetime_epsilon_upper_bound,
    "8.0000000000000004")
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(rounded_over)), "invalid reusable")

  delta_one <- .capsule_status_fixture(
    capsules = 1, releases = 1, lifetime_max = 2,
    capsule_delta = 0.5)
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(delta_one)), "invalid reusable")

  withr::local_options(list(OutDec = ","))
  expect_identical(
    .dsvert_joint_dp_lifetime_exact_total(0.8, 10),
    "8.0000000000000004")
})

test_that("the public DP status uses only the signed Synopsis bootstrap", {
  conns <- list(site_a = "a", site_b = "b")
  sentinel <- structure(list(ok = TRUE), class = c(
    "ds.vertDPStatus", "dsvert_synopsis_bootstrap_v1", "list"))
  observed <- NULL
  value <- testthat::with_mocked_bindings(
    ds.vertDPStatus(conns),
    .dsvert_dp_status_impl = function(
        datasources, .aggregate) {
      observed <<- datasources
      sentinel
    },
    .package = "dsVertClient")
  expect_identical(value, sentinel)
  expect_identical(observed, conns)
})

test_that("capsule handshake requires gates and consistent telemetry", {
  conns <- lapply(c("site_a", "site_b", "site_c"), function(value) {
    structure(1, class = "fake")
  })
  names(conns) <- c("site_a", "site_b", "site_c")

  gated <- .capsule_status_fixture()
  gated$site_a$privacy_contract$operation_limit <- FALSE
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(gated)), "required reusable")

  old_schema <- .capsule_status_fixture()
  old_schema$site_a$version <- "dsvert-joint-dp-capsule-status-v4"
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(old_schema)), "required reusable")

  overstated_rollback <- .capsule_status_fixture()
  overstated_rollback$site_a$privacy_contract$
    simultaneous_designated_history_rollback_protection <-
      "guaranteed_without_external_linearizable_cas"
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(overstated_rollback)),
    "required reusable")

  inaccurate <- .capsule_status_fixture()
  inaccurate$site_a$privacy_contract$accuracy_depends_on_request_history <-
    TRUE
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(inaccurate)), "required reusable")

  overpublished <- .capsule_status_fixture(capsules = 1, releases = 2)
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(overpublished)),
    "invalid release-instance composition telemetry")

  malformed_rotation <- .capsule_status_fixture()
  malformed_rotation$site_a$noise_root$automatic_rotation <- 1
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(malformed_rotation)),
    "required reusable")

  malformed_rotation_count <- .capsule_status_fixture()
  malformed_rotation_count$site_a$noise_root$rotation_count <- 0.5
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(malformed_rotation_count)),
    "required reusable")

  impossible_external_rotation <- .capsule_status_fixture()
  root <- impossible_external_rotation$site_a$noise_root
  root$external <- TRUE
  root$storage <- "hsm_kms_provider"
  root$automatic_generation <- FALSE
  root$automatic_recovery <- FALSE
  root$automatic_rotation <- TRUE
  root$rotation_count <- 1
  impossible_external_rotation$site_a$noise_root <- root
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(impossible_external_rotation)),
    "required reusable")

  divergent <- .capsule_status_fixture()
  divergent$site_c$composition_telemetry$capsules_created <- 8
  divergent$site_c$composition_telemetry$remaining_distinct_capsules <- 0
  divergent$site_c$composition_telemetry$cumulative_epsilon_upper_bound <- 8
  divergent$site_c$composition_telemetry$cumulative_delta_upper_bound <-
    8 * 2^-100
  divergent$site_c$release_instance_telemetry$
    remaining_distinct_capsules <- 0
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(divergent)),
    "disagree on capsule composition telemetry")

  wrong_role <- .capsule_status_fixture()
  wrong_role$site_b$role$designated_noise_peer <- TRUE
  wrong_role$site_b$role$allocator <- "authenticated_ready"
  wrong_role$site_b$composition_telemetry <-
    wrong_role$site_a$composition_telemetry
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(wrong_role)),
    "required reusable")

  substituted_designation <- .capsule_status_fixture()
  substituted_designation$site_b$policy$designated_noise_peers <-
    c("site_a", "site_b")
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(substituted_designation)),
    "disagree|required reusable|invalid reusable", ignore.case = TRUE)

  changed_epoch <- .capsule_status_fixture()
  changed_epoch$site_c$noise_root$privacy_epoch <- 2
  changed_epoch$site_c$release_instance_telemetry$releases_published <- 3
  changed_epoch$site_c$release_instance_telemetry$
    remaining_distinct_capsules <- 1
  changed_epoch$site_c$release_instance_telemetry$
    cumulative_epsilon_upper_bound <- 3
  changed_epoch$site_c$release_instance_telemetry$
    cumulative_delta_upper_bound <- 3 * 2^-100
  expect_error(.dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(changed_epoch)),
    "disagree on release-instance telemetry")

  changed_epoch$site_c$release_instance_telemetry <-
    changed_epoch$site_a$release_instance_telemetry
  rotated <- .dsvert_joint_dp_capsule_status_impl(
    conns, .capsule_status_aggregate(changed_epoch))
  expect_identical(rotated$site_a$noise_root$privacy_epoch, 1)
  expect_identical(rotated$site_c$noise_root$privacy_epoch, 2)
})
