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
