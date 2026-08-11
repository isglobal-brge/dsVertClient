.frequency_execution_b64url <- function(value) {
  sub("=+$", "", chartr("+/", "-_", gsub(
    "[\r\n]", "", jsonlite::base64_enc(value))), perl = TRUE)
}

.frequency_execution_client_fixture <- function() {
  source <- "site_source"
  secondary <- "site_secondary"
  witness <- "site_witness"
  private <- openssl::read_ed25519_key(as.raw(seq.int(0L, 31L)))
  secondary_pk <- .frequency_execution_b64url(as.list(private)$pubkey$data)
  source_pk <- .frequency_execution_b64url(as.raw(seq.int(32L, 63L)))
  hex <- function(value) digest::digest(value, "sha256", serialize = FALSE)
  roles <- list(
    source_owner = source_pk,
    secondary_noise_authority = secondary_pk)
  profile <- list(
    mechanism = "discrete_laplace_convolution",
    sampler = "binary_geometric_convolution_laplace_v3",
    output_transform = "clamp_each_coordinate_to_public_bounds_v1")
  worker <- list(
    selected_primitive = "independent_full_global_draw_convolution_ring128_v3",
    selected_profile = profile, ring_bits = 128L, frac_bits = 0L,
    d = 3L, chunk_coordinates = 3L, chunk_count = 1L,
    raw_bound = list(lower = "0", upper = "1000", scale = 0L),
    authority_roles = roles, release_contract_hash = hex("release"))
  contract <- list(
    artifact_key = hex("artifact"),
    semantic = list(analysis = list(effective_arguments = list(
      sampler_plan = list(coordinate_order_sha256 = hex("coordinates"))))))
  compiled <- list(
    claim = list(version = "claim"),
    config = list(factor_domain = list(
      variable_name = "category", levels = as.list(c("a", "b", "c")))),
    receipts = list(), contract = contract,
    contract_sha256 = hex("contract"), binding = list(sha256 = hex("binding")),
    worker_static = worker, worker_static_sha256 = hex("worker"))
  authorities <- c(source, secondary)
  datasources <- stats::setNames(lapply(
    c(authorities, witness), function(peer) structure(list(peer), class = "mock")),
    c(authorities, witness))
  context <- list(
    compiled = compiled, session_id = "00000000-0000-4000-8000-000000000001",
    authorities = authorities, source = source, secondary = secondary,
    conns = datasources[authorities], authorizations = stats::setNames(
      list(list(role = "source"), list(role = "secondary")), authorities),
    transport = stats::setNames(list("source-key", "secondary-key"), authorities),
    geometry = .dsvert_dp_frequency_execution_geometry_v1(worker),
    authorization_set_sha256 = hex("authorizations"))
  operation <- "op_00000000000000000000000000000001"
  window <- .dsvert_dp_frequency_execution_window_v1(context$geometry, 0L)
  expected_context <- .dsvert_dp_frequency_execution_expected_context_v1(
    context, operation, window)
  ciphertext <- strrep("A", .DSVERT_CLIENT_DP_FREQUENCY_CIPHERTEXT_CHARS)
  transfer <- list(
    ticket = "ticket", transfer_id = "tb_00000000000000000000000000000001",
    capability_id = .DSVERT_CLIENT_DP_FREQUENCY_SOURCE_CAPABILITY,
    sender_name = source, recipient_name = secondary,
    payload_chars = as.numeric(nchar(ciphertext, type = "bytes")),
    payload_sha256 = digest::digest(ciphertext, "sha256", serialize = FALSE))
  issued <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_SOURCE_VERSION, state = "issued",
    artifact_key = contract$artifact_key, window_index = 0L, window_count = 1L,
    context = expected_context, ciphertext_chars = ciphertext,
    transfer = transfer, intermediate_values_exposed = FALSE))
  delivered <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_SOURCE_VERSION, state = "delivered",
    artifact_key = contract$artifact_key, window_index = 0L, window_count = 1L,
    intermediate_values_exposed = FALSE))
  values <- c(4, 4, 3)
  chunk_hash <- .dsvert_dp_frequency_execution_chunk_hash_v1(values, 0L, 0L)
  window_core <- list(
    version = "dsvert-dp-frequency-final-window-v1", window_index = 0L,
    coordinate_offset = 0L, coordinate_count = 3L,
    chunk_hashes = list(chunk_hash))
  window_hash <- .dsvert_dp_frequency_client_hash_v1(
    .DSVERT_CLIENT_DP_FREQUENCY_WINDOW_DOMAIN, window_core)
  release_core <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_RELEASE_VERSION,
    artifact_key = contract$artifact_key,
    contract_sha256 = compiled$contract_sha256,
    analysis_binding_sha256 = compiled$binding$sha256,
    worker_static_sha256 = compiled$worker_static_sha256,
    authorization_set_sha256 = context$authorization_set_sha256,
    release_contract_hash = worker$release_contract_hash,
    primitive = worker$selected_primitive, mechanism = profile$mechanism,
    sampler = profile$sampler, d = 3L, chunk_coordinates = 8192L,
    chunk_count = 1L, window_count = 1L,
    coordinate_order_sha256 = contract$semantic$analysis$effective_arguments$
      sampler_plan$coordinate_order_sha256,
    bounds = worker$raw_bound, authority_roles = roles,
    final_chunk_hashes = list(chunk_hash),
    final_window_hashes = list(window_hash),
    final_vector_root = .dsvert_dp_frequency_execution_merkle_v1(chunk_hash),
    postprocessing = profile$output_transform,
    intermediate_values_exposed = FALSE, public_openings = 1L))
  release_sha256 <- .dsvert_dp_frequency_client_hash_v1(
    .DSVERT_CLIENT_DP_FREQUENCY_RELEASE_DOMAIN, release_core)
  signed <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    release_core, list(release_sha256 = release_sha256)))
  signature <- .frequency_execution_b64url(openssl::ed25519_sign(
    charToRaw(paste0(.DSVERT_CLIENT_DP_FREQUENCY_RELEASE_SIGNATURE_DOMAIN,
      .dsvert_joint_dp_client_json(signed))), private))
  release <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    signed, list(signature = signature)))
  final <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_FINAL_VERSION,
    state = "release_committed", artifact_key = contract$artifact_key,
    release = release, intermediate_values_exposed = FALSE))
  public_window <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = "dsvert-dp-frequency-final-window-v1", window_index = 0L,
    coordinate_offset = 0L, coordinate_count = 3L,
    chunks = list(list(
      version = "dsvert-dp-frequency-final-chunk-v1", chunk_index = 0L,
      coordinate_offset = 0L, coordinate_count = 3L,
      values = as.list(as.character(values)), chunk_sha256 = chunk_hash)),
    window_sha256 = window_hash))
  replay <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_FINAL_VERSION,
    state = "release_replay", release = release, window = public_window,
    intermediate_values_exposed = FALSE))
  prepared <- list(
    session_id = context$session_id, contract = contract,
    transport = context$transport)
  list(
    context = context, datasources = datasources, prepared = prepared,
    operation = operation, window = window, expected_context = expected_context,
    issued = issued, delivered = delivered, final = final, replay = replay,
    release = release, values = values)
}

test_that("Frequency execution endpoints are retry-safe typed producers", {
  endpoints <- c(
    "dsvertDPFrequencySourceWindowDS",
    "dsvertDPFrequencyFinalizeWindowDS",
    "dsvertDPFrequencyReplayDS")
  expect_true(exists(
    ".dsvert_dp_frequency_execute_v1", mode = "function", inherits = TRUE))
  expect_true(all(endpoints %in% .DSVERT_IDEMPOTENT_TYPED_PRODUCERS))
})

test_that("Frequency release and replay validators match the signed wire", {
  fixture <- .frequency_execution_client_fixture()
  expect_identical(.dsvert_dp_frequency_execution_source_v1(
    fixture$issued, fixture$context, fixture$expected_context,
    fixture$window), fixture$issued)
  expect_identical(.dsvert_dp_frequency_execution_source_v1(
    fixture$delivered, fixture$context, fixture$expected_context,
    fixture$window, issued = FALSE), fixture$delivered)
  expect_identical(.dsvert_dp_frequency_execution_release_v1(
    fixture$release, fixture$context), fixture$release)
  replay <- .dsvert_dp_frequency_execution_replay_v1(
    fixture$replay, fixture$context, fixture$window, fixture$release)
  expect_identical(replay$values, fixture$values)
  expect_identical(
    .dsvert_dp_frequency_execution_chunk_hash_v1(c(4, 4, 3), 0L, 0L),
    "98f83d6d48d205ce7f60bea8ad8c79b99ff0384c1113a3c95daf9a68877c23e1")
  expect_identical(.dsvert_dp_frequency_execution_merkle_v1(c(
    strrep("0", 64L), strrep("1", 64L), strrep("a", 64L))),
    "4e51a8fe50d056988e7f8a3564a02cb1c9a125680a3340655fa3797b4908d8eb")
})

test_that("Frequency execution uses only its two roles and cleans in order", {
  fixture <- .frequency_execution_client_fixture()
  state <- new.env(parent = emptyenv())
  state$source_calls <- 0L
  state$calls <- list()
  state$cleanups <- character()
  aggregate <- function(conns, expr, ...) {
    command <- as.character(expr[[1L]])
    state$calls[[length(state$calls) + 1L]] <- list(
      command = command, peers = names(conns))
    value <- switch(command,
      dsvertDPFrequencySourceWindowDS = {
        state$source_calls <- state$source_calls + 1L
        if (state$source_calls < 3L) fixture$issued else fixture$delivered
      },
      dsvertDPFrequencyFinalizeWindowDS = fixture$final,
      dsvertDPFrequencyReplayDS = fixture$replay,
      stop("unexpected endpoint: ", command))
    stats::setNames(list(value), names(conns))
  }
  retry <- function(attempt, classify, ...) {
    result <- classify(attempt())
    if (identical(result$state, "missing")) result <- classify(attempt())
    result
  }
  stored <- 0L
  result <- .dsvert_dp_frequency_execute_v1(
    fixture$prepared, "D", fixture$datasources, .aggregate = aggregate,
    .new_context = function() list(operation_id = fixture$operation),
    .store_typed = function(blob, transfer, conn, producer_conn, ...) {
      stored <<- stored + 1L
      expect_identical(blob, fixture$issued$ciphertext_chars)
      expect_identical(transfer, fixture$issued$transfer)
      invisible(TRUE)
    }, .retry = retry,
    .execution_context = function(...) fixture$context,
    .frequency_cleanup = function(conns, ...) {
      state$cleanups <- c(state$cleanups, "frequency")
      expect_identical(names(conns), fixture$context$authorities)
      invisible(TRUE)
    }, .exact_cleanup = function(conns, ...) {
      state$cleanups <- c(state$cleanups, "exact")
      expect_identical(names(conns), fixture$context$authorities)
      invisible(TRUE)
    })
  expect_identical(result$values, fixture$values)
  expect_identical(stored, 1L)
  expect_identical(state$cleanups, c("frequency", "exact"))
  expect_true(all(vapply(state$calls, function(call) {
    identical(call$peers, fixture$context$source) ||
      identical(call$peers, fixture$context$secondary)
  }, logical(1L))))
  expect_false(any(c("transport", "ciphertext_chars") %in%
    names(result$proof)))
  expect_identical(result$proof$release, fixture$release)
})

test_that("Frequency execution rejects tampering and oversized shapes early", {
  fixture <- .frequency_execution_client_fixture()
  bad <- fixture$release
  bad$release_sha256 <- strrep("f", 64L)
  expect_error(.dsvert_dp_frequency_execution_release_v1(
    bad, fixture$context), "corrupt|misbound")
  bad <- fixture$release
  bad$signature <- paste0("A", substring(bad$signature, 2L))
  expect_error(.dsvert_dp_frequency_execution_release_v1(
    bad, fixture$context), "signature|verification")
  bad <- fixture$replay
  bad$window$chunks[[1L]]$values[[1L]] <- "5"
  expect_error(.dsvert_dp_frequency_execution_replay_v1(
    bad, fixture$context, fixture$window, fixture$release), "corrupt")
  bad <- fixture$replay
  bad$window$chunks[[1L]]$values <- rep(list("0"), 8193L)
  expect_error(.dsvert_dp_frequency_execution_replay_v1(
    bad, fixture$context, fixture$window, fixture$release), "corrupt")
  expect_error(.dsvert_dp_frequency_execution_integer_v1(
    strrep("1", 10000L), "index"), "Invalid")
  expect_error(.dsvert_dp_frequency_execution_cleanup_capability_v1(
    strrep("A", 16385L), "site", fixture$context$session_id,
    strrep("a", 64L), c(site = strrep("A", 43L))), "Invalid")
  geometry <- .dsvert_dp_frequency_execution_geometry_v1(within(
    fixture$context$compiled$worker_static, {
      d <- 65537L; chunk_coordinates <- 8192L; chunk_count <- 9L
    }))
  expect_identical(geometry$window_count, 2L)
  expect_identical(geometry$chunk_count, 9L)
})
