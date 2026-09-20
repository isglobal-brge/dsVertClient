test_that("fixed-rho GEE reaches the exact transport only with its Ring128 family binding", {
  servers <- c("site_a", "site_b")
  peers <- setNames(paste0("dsv1_", c(strrep("a", 64), strrep("b", 64))), servers)
  operation <- "grouped-gee-fixed-rho-staged-v1"
  purpose <- paste0(operation, "/", strrep("c", 64))
  states <- setNames(lapply(servers, function(server) list(
    capability_id = "exact_gc_v1", peer_id = peers[[server]],
    peer_peer_id = peers[[setdiff(servers, server)]],
    role = if (server == "site_a") "garbler" else "evaluator",
    context_hash = strrep("d", 64), operation = operation,
    output_kind = "grouped-gee-fixed-rho-staged-ring128-share-v1",
    purpose = purpose, source_producer = "dp.gee-fixed-rho-grid-cross.staged-v1",
    ring_bits = 128L, frac_bits = 0L, vector_len = 42L,
    threshold = "", chunk_bytes = 65536L, ttl_seconds = 10L,
    max_runtime_seconds = 120L, worker_heartbeat = 1,
    state = "running", stored = FALSE)), servers)
  conns <- setNames(lapply(servers, function(...) structure(list(), class = "mock")), servers)
  calls <- 0L
  aggregate <- function(conns, expr, ...) {
    calls <<- calls + 1L
    expect_true(is.list(expr) && !is.call(expr))
    expect_true(all(vapply(expr, function(x)
      identical(as.character(x[[1L]]), "exactGCExchangeDS"), logical(1L))))
    setNames(lapply(names(conns), function(server) list(
      capability_id = "exact_gc_v1", peer_id = peers[[server]],
      state = "complete", stored = TRUE, inbound_size = 0,
      outbound = NULL, worker_heartbeat = 2)), names(conns))
  }
  args <- list(datasources = conns, server_names = servers, servers = 1:2,
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = paste0("op_", strrep("3", 32)),
    source_key = paste0("exact_gc_in_", strrep("3", 32)),
    output_key = paste0("exact_gc_out_", strrep("3", 32)),
    operation = operation, ring = 128L, frac_bits = 0L, vector_len = 42L,
    purpose = purpose, transport_ready = TRUE, initialized = states,
    timeout_seconds = 1, .aggregate = aggregate)
  result <- do.call(.dsvert_exact_gc_run, args)
  expect_identical(result$operation, operation)
  expect_identical(result$ring_bits, 128L)
  expect_gt(calls, 0L)
  expect_false(any(grepl("share|payload|release|result", names(result))))
  successful_calls <- calls
  for (mutation in list(
      function(x) { x$ring <- 64L; x },
      function(x) { x$frac_bits <- 16L; x },
      function(x) { x$purpose <- paste0("grouped-lmm-staged-v1/", strrep("c", 64)); x },
      function(x) { x$purpose <- paste0(operation, "/", strrep("c", 63)); x },
      function(x) { x$operation <- "grouped-gee-estimated-alpha-staged-v1"; x })) {
    expect_error(do.call(.dsvert_exact_gc_run, mutation(args)), "operation (contract|shape)")
  }
  expect_identical(calls, successful_calls)
  validate <- function(value) .dsvert_exact_gc_validate_init(value, servers,
    operation, 128L, 0L, 42L, purpose)
  for (mutation in list(
      function(x) { x[[1L]]$output_kind <- "grouped-lmm-staged-ring128-share-v1"; x },
      function(x) { x[[2L]]$operation <- "grouped-glmm-staged-v1"; x },
      function(x) { x[[1L]]$vector_len <- 21L; x },
      function(x) { x[[2L]]$purpose <- paste0(operation, "/", strrep("e", 64)); x })) {
    expect_error(validate(mutation(states)), "operation contract")
  }
})

.gee_sampler_client_fixture <- function(family) {
  f <- .grouped_cross_client_fixture(family)
  f$raw$parameters$composition <- "staged_fixed_rho_v1"
  spec <- .dsvert_dp_grouped_grid_cross_spec(f$raw, f$policy, f$schema)
  artifact <- .dsvert_dp_grouped_grid_cross_artifact(spec)
  contract <- f$sign(list(version = f$contract$version, spec = spec,
    artifact = artifact,
    source_contract = .dsvert_dp_grouped_grid_cross_source_contract(spec, artifact)))
  admitted <- .dsvert_dp_glm_grid_profile_admit(contract, f$policy,
    f$schema_manifest)
  workload <- .dsvert_dp_grouped_cross_workload_artifact(admitted)
  list(artifact = list(artifact_key = strrep("a", 64), semantic = list(
    catalog_projection = list(catalog = list(families = list(
      gaussian_models = list(artifacts = list(gee = workload))))))),
    layout = list(coordinate_count = 43L), physical = list(
      backend_selection = list(policy_version =
        "dsvert-lmm-grid-exact-gc-cost-policy-v1"),
      full_plan = list(total_coordinate_count = 43L,
        maximum_chunk_coordinates = 36L, epsilon = 8, delta = 2^-100,
        sensitivity_steps = "123456789", complete_epsilon_per_peer = TRUE)))
}

test_that("staged fixed-rho GEE schedules 16 coordinates without changing its full plan", {
  execution <- list(geometry = list(coordinate_count = 43L,
    public_chunk_coordinates = 8192L))
  for (family in c("binomial_gee", "poisson_gee")) {
    compiled <- .gee_sampler_client_fixture(family)
    before <- serialize(compiled, NULL)
    size <- .dsvert_dp_synopsis_runner_exact_chunk_size(compiled, execution)
    expect_equal(size, 16L, info = family)
    expect_equal(.dsvert_dp_synopsis_runner_exact_chunk_size(compiled, execution),
      size, info = family)
    expect_equal(pmin(size, 43L - (0:2) * size), c(16L, 16L, 11L), info = family)
    expect_identical(serialize(compiled, NULL), before, info = family)
    smaller <- compiled
    smaller$physical$full_plan$maximum_chunk_coordinates <- 8L
    expect_equal(.dsvert_dp_synopsis_runner_exact_chunk_size(smaller, execution),
      8L, info = family)
    for (control in c("lmm", "binomial_glmm", "poisson_glmm", "legacy_gee",
                      "wrong_version")) {
      changed <- compiled
      artifact <- changed$artifact$semantic$catalog_projection$catalog$families$
        gaussian_models$artifacts$gee
      if (control == "legacy_gee") {
        artifact$composition <- NULL
      } else if (control == "wrong_version") {
        artifact$version <- "bounded-lmm-cross-grid-v1"
      } else {
        artifact$family <- control
        artifact$version <- paste0("bounded-", gsub("_", "-", control),
          "-cross-grid-v1")
        artifact$spec_version <- paste0(control, "_grid_cross_v1")
      }
      changed$artifact$semantic$catalog_projection$catalog$families$
        gaussian_models$artifacts$gee <- artifact
      expect_equal(.dsvert_dp_synopsis_runner_exact_chunk_size(changed, execution),
        36L, info = paste(family, control))
    }
  }
})

test_that("signed staged GEE START receipts bind the 16-coordinate sampler windows", {
  peers <- c("site_a", "site_b")
  keys <- setNames(lapply(peers, function(peer) openssl::ed25519_keygen()), peers)
  b64 <- function(value) sub("=+$", "", chartr("+/", "-_",
    gsub("[\r\n]", "", jsonlite::base64_enc(value))))
  pins <- vapply(keys, function(key) b64(tail(as.raw(as.list(key)$pubkey), 32L)),
    character(1L))
  trusted <- list(context = list(pinset = pins))
  execution <- list(execution_id = strrep("b", 64), geometry = list(
    coordinate_count = 43L, public_chunk_coordinates = 8192L))
  for (family in c("binomial_gee", "poisson_gee")) {
    compiled <- .gee_sampler_client_fixture(family)
    responses_for <- function(index, offset, count) setNames(lapply(peers,
      function(peer) {
        unsigned <- list(version = .DSVERT_CLIENT_SYNOPSIS_EXACT_START_VERSION,
          phase = "synopsis_exact_gc_initialized",
          execution_id = execution$execution_id,
          artifact_key = compiled$artifact$artifact_key,
          contract_sha256 = strrep("c", 64), attempt_sha256 = strrep("d", 64),
          source_contract_sha256 = strrep("e", 64), local_authority = list(
            peer_name = peer, identity_pk = unname(pins[[peer]]),
            role = c("primary_noise_authority", "secondary_noise_authority")[[
              match(peer, peers)]]),
          chunk_index = index, coordinate_offset = offset, coordinate_count = count,
          backend_selection_sha256 = strrep("f", 64),
          worker_contract_sha256 = strrep("1", 64), binding_sha256 = strrep("2", 64),
          operation_id = paste0("op_", strrep("3", 32)),
          purpose = "joint-dp-vector-laplace-v3/fixture",
          local_chunk_durable = FALSE, intermediate_payload_exposed = FALSE,
          source_share_exposed = FALSE, private_seed_exposed = FALSE,
          preclamp_values_exposed = FALSE)
        signature <- b64(openssl::ed25519_sign(charToRaw(paste0(
          .DSVERT_CLIENT_SYNOPSIS_EXACT_START_DOMAIN,
          .dsvert_joint_dp_client_json(unsigned))), keys[[peer]]))
        .dsvert_joint_dp_client_json(list(
          version = .DSVERT_CLIENT_SYNOPSIS_EXACT_START_RESPONSE_VERSION,
          receipt = c(unsigned, list(signature = signature)),
          initialization = list(state = "running", stored = FALSE)))
      }), peers)
    for (index in 0:2) {
      count <- c(16L, 16L, 11L)[[index + 1L]]
      accepted <- .dsvert_dp_synopsis_runner_exact_start_set(
        responses_for(index, index * 16L, count), peers, trusted, compiled,
        execution, index)
      expect_false(accepted$complete, info = family)
      expect_equal(accepted$receipts[[1L]]$coordinate_count, count, info = family)
    }
    for (old in list(c(0L, 0L, 36L), c(1L, 36L, 7L))) {
      expect_error(.dsvert_dp_synopsis_runner_exact_start_set(
        responses_for(old[[1L]], old[[2L]], old[[3L]]), peers, trusted,
        compiled, execution, old[[1L]]), "misbound exact-GC START", info = family)
    }
  }
})
