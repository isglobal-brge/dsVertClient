.frequency_public_fixture_code <- parse(testthat::test_path(
  "test-dp-frequency-analysis.R"))
for (.frequency_public_fixture_expression in .frequency_public_fixture_code) {
  if (is.call(.frequency_public_fixture_expression) && identical(
      .frequency_public_fixture_expression[[1L]], quote(`<-`))) {
    .frequency_public_fixture_name <- as.character(
      .frequency_public_fixture_expression[[2L]])
    if (startsWith(.frequency_public_fixture_name, ".frequency_client_"))
      eval(.frequency_public_fixture_expression)
  }
}
rm(.frequency_public_fixture_code, .frequency_public_fixture_expression,
   .frequency_public_fixture_name)

.frequency_public_execution <- function(fixture, values = c(4, 4, 3)) {
  compiled <- fixture$compiled
  authorities <- compiled$authorities
  worker <- compiled$worker_static
  geometry <- .dsvert_dp_frequency_execution_geometry_v1(worker)
  chunks <- vapply(seq_len(geometry$chunk_count) - 1L, function(index) {
    offset <- index * .DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES
    count <- min(.DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES,
                 geometry$d - offset)
    .dsvert_dp_frequency_execution_chunk_hash_v1(
      values[seq.int(offset + 1L, offset + count)], index, offset)
  }, character(1L))
  windows <- vapply(seq_len(geometry$window_count) - 1L, function(index) {
    window <- .dsvert_dp_frequency_execution_window_v1(geometry, index)
    range <- seq.int(window$first_chunk + 1L,
                     window$first_chunk + window$chunks)
    .dsvert_dp_frequency_client_hash_v1(
      .DSVERT_CLIENT_DP_FREQUENCY_WINDOW_DOMAIN, list(
        version = "dsvert-dp-frequency-final-window-v1",
        window_index = index, coordinate_offset = window$offset,
        coordinate_count = window$count,
        chunk_hashes = as.list(chunks[range])))
  }, character(1L))
  authorization_hash <- .dsvert_dp_frequency_client_hash_v1(
    .DSVERT_CLIENT_DP_FREQUENCY_AUTH_SET_DOMAIN,
    unname(fixture$authorizations[authorities]))
  profile <- worker$selected_profile
  core <- .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_RELEASE_VERSION,
    artifact_key = compiled$contract$artifact_key,
    contract_sha256 = compiled$contract_sha256,
    analysis_binding_sha256 = compiled$binding$sha256,
    worker_static_sha256 = compiled$worker_static_sha256,
    authorization_set_sha256 = authorization_hash,
    release_contract_hash = worker$release_contract_hash,
    primitive = worker$selected_primitive, mechanism = profile$mechanism,
    sampler = profile$sampler, d = worker$d, chunk_coordinates = 8192L,
    chunk_count = worker$chunk_count, window_count = geometry$window_count,
    coordinate_order_sha256 = compiled$contract$semantic$analysis$
      effective_arguments$sampler_plan$coordinate_order_sha256,
    bounds = worker$raw_bound, authority_roles = worker$authority_roles,
    final_chunk_hashes = as.list(chunks),
    final_window_hashes = as.list(windows),
    final_vector_root = .dsvert_dp_frequency_execution_merkle_v1(chunks),
    postprocessing = profile$output_transform,
    intermediate_values_exposed = FALSE, public_openings = 1L))
  release_hash <- .dsvert_dp_frequency_client_hash_v1(
    .DSVERT_CLIENT_DP_FREQUENCY_RELEASE_DOMAIN, core)
  signed <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    core, list(release_sha256 = release_hash)))
  finalizer <- authorities[[2L]]
  release <- .dsvert_dp_analysis_client_canonical_value_v1(c(
    signed, list(signature = .frequency_client_sign(charToRaw(paste0(
      .DSVERT_CLIENT_DP_FREQUENCY_RELEASE_SIGNATURE_DOMAIN,
      .dsvert_joint_dp_client_json(signed))),
      fixture$keys$private[[finalizer]]))))
  list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_RESULT_VERSION,
    operation_id = "op_00000000000000000000000000000001",
    variable_name = compiled$config$factor_domain$variable_name,
    levels = unlist(compiled$config$factor_domain$levels, use.names = FALSE),
    values = unname(values), source_owner = authorities[[1L]],
    finalizer_peer = finalizer,
    proof = list(
      session_id = fixture$session_id, claim = fixture$claim,
      config = fixture$config, receipts = fixture$receipts,
      contract = compiled$contract, worker_static = worker,
      authorities = authorities, authorizations = fixture$authorizations,
      release = release))
}

.frequency_public_run <- function(fixture, execution = NULL) {
  if (is.null(execution)) execution <- .frequency_public_execution(fixture)
  calls <- character()
  datasources <- stats::setNames(lapply(names(fixture$pins), function(peer) {
    structure(list(peer), class = "mock")
  }), names(fixture$pins))
  result <- .dsvert_dp_frequency_impl(
    "D", "status", fixture$source_peer, datasources,
    function(...) stop("unexpected raw DSI call", call. = FALSE),
    .prepare = function(data_name, variable_name, source_owner, ...) {
      calls <<- c(calls, "prepare")
      expect_identical(c(data_name, variable_name, source_owner),
                       c("D", "status", fixture$source_peer))
      list(prepared = TRUE)
    }, .execute = function(prepared, data_name, datasources, ...) {
      calls <<- c(calls, "execute")
      expect_true(isTRUE(prepared$prepared))
      execution
    })
  list(result = result, calls = calls)
}

test_that("Frequency requires an explicit source before any DSI call", {
  calls <- 0L
  expect_error(.dsvert_dp_frequency_impl(
    "D", "status", NULL, NULL, function(...) {
      calls <<- calls + 1L
    }), "explicit source")
  expect_identical(calls, 0L)
})

test_that("Frequency routes K2/K3/K5 and both physical backends", {
  cases <- data.frame(
    k = c(2L, 3L, 5L),
    kind = c("convolution", "gaussian", "convolution"))
  for (index in seq_len(nrow(cases))) {
    fixture <- .frequency_client_fixture(cases$k[[index]], cases$kind[[index]])
    run <- .frequency_public_run(fixture)
    result <- run$result
    expect_identical(run$calls, c("prepare", "execute"))
    expect_s3_class(result, "ds.vertDPFrequency")
    expect_identical(result$levels, c("case", "control", "other"))
    expect_identical(result$counts, c(case = 4, control = 4, other = 3))
    expect_equal(sum(result$proportions), 1)
    expect_identical(result$source_owner, fixture$source_peer)
    expect_identical(length(result$proof$receipts), cases$k[[index]])
    expect_match(result$primitive, cases$kind[[index]])
    certificate <- fixture$config$backend_selection$selected_accuracy_certificate
    expect_identical(result$accuracy_implementation_tv_upper_bound, min(
      1, .dsvert_dp_vector_fraction_upper(
        certificate$release_tv_upper_numerator,
        certificate$release_tv_upper_denominator)))
    expect_false(any(c("transport", "ciphertext_chars") %in%
                       names(result$proof)))
    expect_false(any(c("capsule_id", "history_gate", "operation_limit") %in%
                       names(result)))
  }
})

test_that("Frequency recompiles and authenticates every public proof layer", {
  fixture <- .frequency_client_fixture(3L, "convolution")
  execution <- .frequency_public_execution(fixture)
  mutations <- list(
    claim = function(x) {
      x$proof$claim$signature <- .frequency_client_flip(
        x$proof$claim$signature); x
    }, config = function(x) {
      x$proof$config$coordinate_upper_bound <- 999; x
    }, receipt = function(x) {
      x$proof$receipts[[1L]]$signature <- .frequency_client_flip(
        x$proof$receipts[[1L]]$signature); x
    }, authorization = function(x) {
      x$proof$authorizations[[1L]]$signature <- .frequency_client_flip(
        x$proof$authorizations[[1L]]$signature); x
    }, release = function(x) {
      x$proof$release$signature <- .frequency_client_flip(
        x$proof$release$signature); x
    }, values = function(x) {
      x$values[[1L]] <- x$values[[1L]] + 1; x
    })
  for (mutate in mutations) {
    expect_error(.dsvert_dp_frequency_public_execution_v1(
      mutate(execution)), "Frequency|frequency")
  }
})

test_that("Frequency caps proof peers before compiling envelopes", {
  proof_fields <- c(
    "session_id", "claim", "config", "receipts", "contract",
    "worker_static", "authorities", "authorizations", "release")
  proof <- stats::setNames(vector("list", length(proof_fields)), proof_fields)
  peers <- sprintf("site_%04d", seq_len(4097L))
  proof$receipts <- stats::setNames(rep(list(NULL), length(peers)), peers)
  execution <- list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_RESULT_VERSION,
    operation_id = "op_00000000000000000000000000000001",
    variable_name = "status", levels = "a", values = 0,
    source_owner = "site_0001", finalizer_peer = "site_0002", proof = proof)
  expect_error(.dsvert_dp_frequency_public_execution_v1(execution),
               "Invalid public Frequency proof")
})

test_that("frequency regions are O(d), sharp, and safe for infeasible boxes", {
  counts <- c(a = 3, b = 5, c = 2)
  regions <- .dsvert_dp_frequency_regions(counts, radius = 2, capacity = 10)
  grids <- expand.grid(a = 1:5, b = 3:7, c = 0:4)
  grids <- grids[rowSums(grids) <= 10 & rowSums(grids) > 0, , drop = FALSE]
  for (index in seq_len(nrow(grids))) {
    value <- unlist(grids[index, ], use.names = FALSE)
    proportion <- value / sum(value)
    expect_true(all(proportion >= regions$proportion[, "lower"] - 1e-15))
    expect_true(all(proportion <= regions$proportion[, "upper"] + 1e-15))
  }
  infeasible <- .dsvert_dp_frequency_regions(c(a = 9, b = 9), 0, 10)
  expect_identical(unname(infeasible$proportion),
                   matrix(c(0, 0, 1, 1), ncol = 2L))
  body_text <- paste(deparse(body(.dsvert_dp_frequency_regions)),
                     collapse = "\n")
  quadratic <- grepl("setdiff|vapply\\(seq_along", body_text)
  expect_false(quadratic)
  if (!quadratic) {
    large <- .dsvert_dp_frequency_regions(rep(1, 1e6), 1, 1e6)
    expect_identical(nrow(large$proportion), 1000000L)
  }
})

test_that("Frequency inference is zero-call, vectorized, and confidence-aware", {
  for (kind in c("convolution", "gaussian")) {
    fixture <- .frequency_client_fixture(3L, kind)
    frequency <- .frequency_public_run(fixture)$result
    signed <- as.numeric(fixture$config$backend_selection$
      selected_accuracy_certificate$simultaneous_95_abs)
    expect_identical(.dsvert_dp_frequency_radius(frequency, 0.95),
                     min(signed, frequency$coordinate_maximum))
    expect_lte(.dsvert_dp_frequency_radius(frequency, 0.9),
               .dsvert_dp_frequency_radius(frequency, 0.99))
    inference <- ds.vertDPFrequencyInference(
      frequency, level = 0.95, dp_fraction = 0.5)
    expect_s3_class(inference, "ds.vertDPFrequencyInference")
    expect_identical(inference$additional_server_calls, 0L)
    expect_gte(inference$mechanism_radius,
               frequency$accuracy_simultaneous_95_abs)
    expect_true(all(inference$intervals >= 0 & inference$intervals <= 1))
  }
  inference_body <- paste(deparse(body(ds.vertDPFrequencyInference)),
                          collapse = "\n")
  cp_body <- paste(deparse(body(.dsvert_dp_frequency_cp_regions_v1)),
                   collapse = "\n")
  expect_false(grepl("DSI|aggregate|setdiff|vapply\\(seq_along",
                     paste(inference_body, cp_body)))
  wide <- .dsvert_dp_frequency_cp_regions_v1(
    rep(0, 1e6), rep(2, 1e6), 2e6, 0.01)
  expect_identical(dim(wide), c(1000000L, 2L))
})

test_that("Frequency accuracy accepts an exhausted signed release-TV bound", {
  testthat::local_mocked_bindings(
    .dsvert_dp_frequency_client_plan_v1 = function(...) list(
      full_plan_sha256 = strrep("a", 64L)),
    .dsvert_dp_analysis_frequency_profile_v1 = function(primitive) list(
      gaussian = identical(primitive, "gaussian")),
    .dsvert_dp_vector_gaussian_accuracy_steps = function(plan, ...) {
      expect_identical(plan$maximum_noise_magnitude_two_peers, "100")
      list(steps = 20, finite_support = FALSE)
    },
    .package = "dsVertClient")
  config <- list(
    coordinate_upper_bound = 100,
    factor_domain = list(dimension = 3L),
    backend_selection = list(
      summary = list(selected_primitive = "convolution"),
      selected_plan = list(),
      selected_accuracy_certificate = list(
        method = "finite_support_v1",
        simultaneous_95_abs = "18446744073709551615",
        absolute_support = "18446744073709551615",
        release_tv_upper_numerator = "1",
        release_tv_upper_denominator = "1")))
  signed <- .dsvert_dp_frequency_accuracy_v1(config, 0.95)
  exhausted <- .dsvert_dp_frequency_accuracy_v1(config, 0.999)
  expect_identical(signed$radius, 100)
  expect_identical(signed$implementation_tv_upper_bound, 1)
  expect_identical(exhausted$radius, 100)
  expect_true(exhausted$finite_support_fallback)
  gaussian <- config
  gaussian$backend_selection$summary$selected_primitive <- "gaussian"
  gaussian$backend_selection$selected_plan <- list(
    maximum_noise_magnitude_two_peers = "18446744073709551615")
  gaussian$backend_selection$selected_accuracy_certificate$method <-
    "gaussian_plan_v2_subgaussian_mgf_tv_transfer"
  gaussian$backend_selection$selected_accuracy_certificate$
    simultaneous_95_abs <- "7"
  expect_identical(.dsvert_dp_frequency_accuracy_v1(
    gaussian, 0.99)$radius, 20)
})

test_that("Frequency inference rejects tampered returned values", {
  fixture <- .frequency_client_fixture(2L, "convolution")
  frequency <- .frequency_public_run(fixture)$result
  frequency$counts[[1L]] <- frequency$counts[[1L]] + 1
  expect_error(ds.vertDPFrequencyInference(frequency),
               "released, validated ds.vertDPFrequency")
  frequency <- .frequency_public_run(fixture)$result
  frequency$repeated_record_policy <- "legacy_policy"
  expect_error(ds.vertDPFrequencyInference(frequency),
               "released, validated ds.vertDPFrequency")
  frequency <- .frequency_public_run(fixture)$result
  frequency$uncertainty_scope <- "mechanism and sampling uncertainty"
  expect_error(ds.vertDPFrequencyInference(frequency),
               "released, validated ds.vertDPFrequency")
})
