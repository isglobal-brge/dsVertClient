.count_synopsis_e2e_code <- parse(testthat::test_path(
  "test-dp-synopsis-describe-rock-e2e.R"))
for (.count_synopsis_e2e_expression in .count_synopsis_e2e_code) {
  if (is.call(.count_synopsis_e2e_expression) && identical(
      .count_synopsis_e2e_expression[[1L]], quote(`<-`))) {
    .count_synopsis_e2e_name <- as.character(
      .count_synopsis_e2e_expression[[2L]])
    if (.count_synopsis_e2e_name %in% c(
          ".synopsis_describe_real_e2e_server",
          ".synopsis_describe_real_e2e_fixture",
          ".synopsis_describe_real_e2e_session",
          ".synopsis_describe_real_e2e_dispatch")) {
      eval(.count_synopsis_e2e_expression)
    }
  }
}
rm(.count_synopsis_e2e_code, .count_synopsis_e2e_expression,
   .count_synopsis_e2e_name)

test_that("real Synopsis Count is plausible and Rock-replayable at K=2/3/5", {
  server_ns <- .synopsis_describe_real_e2e_server()
  count <- get(".dsvert_dp_count_impl", asNamespace("dsVertClient"),
               inherits = FALSE)
  for (k in c(2L, 3L, 5L)) {
    fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- count("data_peer_a", "peer_b", conns, dispatch)
    expect_s3_class(first, "ds.vertDPCount")
    expect_true(isTRUE(first$released))
    expect_true(isTRUE(first$one_joint_draw))
    expect_identical(first$implementation,
                     .DSVERT_CLIENT_VECTOR_EXACT_BACKEND)
    expect_identical(first$backend_selection$backend,
                     .DSVERT_CLIENT_VECTOR_EXACT_BACKEND)
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(first$source_owner, "peer_a")
    expect_identical(first$server, "peer_b")
    expect_identical(first$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(first$release_provenance$ordered_peer_pinset, k)
    expect_true(is.finite(first$value))
    expect_gte(first$value, 0)
    expect_lte(first$value, 200)
    expect_true(is.finite(first$accuracy_95_abs))
    expect_gte(first$accuracy_95_abs, 0)
    expect_lte(first$accuracy_95_abs, 200)
    true_count <- nrow(fixture$snapshots$peer_a$data_peer_a$data)
    expect_lte(abs(first$value - true_count),
               1.5 * max(first$accuracy_95_abs, 1))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- count("data_peer_a", "peer_b", conns, dispatch)
    expect_identical(replay$value, first$value)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})
