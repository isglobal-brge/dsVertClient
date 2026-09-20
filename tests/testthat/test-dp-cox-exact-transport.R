test_that("Cox staged transport accepts only its typed server-minted operation", {
  servers <- c("site_a", "site_b")
  peers <- setNames(paste0("dsv1_", c(strrep("a", 64L), strrep("b", 64L))), servers)
  purpose <- paste0("cox-loss-staged-v1/", strrep("c", 64L))
  initialized <- setNames(lapply(servers, function(server) list(
    capability_id = "exact_gc_v1", peer_id = peers[[server]],
    peer_peer_id = peers[[setdiff(servers, server)]],
    role = if (server == "site_a") "garbler" else "evaluator",
    context_hash = strrep("d", 64L), operation = "cox-loss-staged-v1",
    output_kind = "cox-loss-staged-ring128-share-v1", purpose = purpose,
    source_producer = "dp.cox-grid-cross.staged-v1",
    ring_bits = 128L, frac_bits = 0L, vector_len = 2L,
    threshold = "", chunk_bytes = 65536L,
    ttl_seconds = 10L, max_runtime_seconds = 120L,
    worker_heartbeat = 1, state = "running", stored = FALSE)), servers)
  conns <- setNames(rep(list(structure(list(), class = "mock")), 2L), servers)
  exchanges <- 0L
  aggregate <- function(conns, expr, ...) {
    if (is.call(expr)) return(setNames(rep(list(TRUE), length(conns)), names(conns)))
    exchanges <<- exchanges + 1L
    expect_true(all(vapply(expr, function(call)
      identical(as.character(call[[1L]]), "exactGCExchangeDS"), logical(1L))))
    setNames(lapply(names(conns), function(server) list(
      capability_id = "exact_gc_v1", peer_id = peers[[server]],
      state = "complete", stored = TRUE, inbound_size = 0,
      outbound = NULL, worker_heartbeat = 2)), names(conns))
  }
  args <- list(datasources = conns, server_names = servers, servers = 1:2,
    session_id = "12345678-1234-4234-9234-123456789abc",
    operation_id = "op_33333333333333333333333333333333",
    source_key = "exact_gc_in_33333333333333333333333333333333",
    output_key = "exact_gc_out_33333333333333333333333333333333",
    operation = "cox-loss-staged-v1", ring = 128L, frac_bits = 0L,
    vector_len = 2L, purpose = purpose, transport_ready = TRUE,
    initialized = initialized, timeout_seconds = 1, .aggregate = aggregate)
  run <- function(changes = list()) do.call(.dsvert_exact_gc_run, modifyList(args, changes))
  result <- run()
  expect_identical(result$operation, "cox-loss-staged-v1")
  expect_identical(result$ring_bits, 128L)
  expect_identical(exchanges, 1L)
  expect_false(any(grepl("share|payload|release|result", names(result), ignore.case = TRUE)))
  expect_error(run(list(ring = 127L)), "staged Cox operation shape")
  expect_error(run(list(frac_bits = 1L)), "Invalid exact MPC operation contract")
  for (bad in c("cox-loss-staged-v1", paste0("grouped-lmm-staged-v1/", strrep("c", 64L)))) {
    expect_error(run(list(purpose = bad)), "staged Cox operation shape")
  }
  for (size in c(0L, 4097L)) {
    expect_error(run(list(vector_len = size)), "Invalid exact MPC operation contract")
  }
  for (kind in c("grouped-lmm-staged-ring128-share-v1", "cross-grid-ring128-share-v2", "xor-bit-share")) {
    changed <- initialized
    changed[[1L]]$output_kind <- kind
    expect_error(run(list(initialized = changed)), "changed the operation contract")
  }
  expect_identical(exchanges, 1L)
  # Cox remains available only through its purpose-specific signed producer.
  advertised <- .dsvert_exact_gc_capability_contract()
  expect_false("cox-loss-staged-v1" %in% advertised$operations)
  expect_false("cox-loss-staged-v1" %in% advertised$core_operations)
})
