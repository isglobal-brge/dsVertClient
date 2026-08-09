.cross_production_guard <- function() {
  environment <- new.env(parent = baseenv())
  environment$.dsvert_identity_test_mode <- function() TRUE
  package_root <- normalizePath(
    file.path(.dsvert_client_source_root(), ".."), mustWork = FALSE)
  source <- file.path(
    package_root, "..", "dsVert", "R", "securityProfileDS.R")
  if (!file.exists(source)) {
    testthat::skip(
      "the companion dsVert source tree is unavailable for this cross-package audit")
  }
  sys.source(source, envir = environment)
  environment$.dsvert_enforce_release_mode
}

test_that("K=2/3/4/5 checked cross runner traverses the real production guard", {
  guard <- .cross_production_guard()
  expected <- c(
    "exactGCVecmulClaimInputsDS", "exactGCVecmulStartDS",
    "exactGCVecmulValidityDS", "exactGCVecmulValidityReceiveDS",
    "exactGCVecmulCommitDS")
  expect_error(
    guard("exactGCVecmulBindInputsDS"),
    "single disclosure-safe profile", fixed = TRUE)
  expect_error(
    guard("exactGCGLMSoftplusPrepareDS"),
    "single disclosure-safe profile", fixed = TRUE)

  for (k in 2:5) {
    server_names <- paste0("site_", letters[seq_len(k)])
    datasources <- stats::setNames(as.list(seq_len(k)), server_names)
    selected <- server_names[1:2]
    commands <- character()
    policy_id <- strrep("1", 64L)
    plan_id <- strrep("2", 64L)
    context_hash <- strrep("3", 64L)
    max_chunk <- .dsvert_exact_gc_direct_mul_max_chunk(128L)
    manifests <- stats::setNames(lapply(seq_along(selected), function(index) {
      list(
        manifest_handle = paste0(
          if (index == 1L) "A" else "B", strrep(as.character(index), 42L)),
        total_n = 2L)
    }), selected)
    aggregate <- function(conns, expr, ...) {
      expressions <- if (is.list(expr) && !is.call(expr)) {
        expr[names(conns)]
      } else {
        stats::setNames(rep(list(expr), length(conns)), names(conns))
      }
      stats::setNames(lapply(names(conns), function(server) {
        expression <- expressions[[server]]
        command <- as.character(expression[[1L]])
        guard(command)
        commands <<- c(commands, command)
        switch(command,
          exactGCVecmulClaimInputsDS = list(
            stored = TRUE, state = "claimed", capability_id = "exact_gc_v1",
            context_hash = context_hash, policy_id = policy_id,
            plan_id = plan_id, backend = "direct-wide",
            bound_x = "256", bound_y = "256", ring_bits = 128L,
            frac_bits = 8L, max_chunk = max_chunk),
          exactGCVecmulStartDS = list(
            capability_id = "exact_gc_v1", state = "running"),
          exactGCVecmulValidityDS = list(
            capability_id = "exact_gc_v1", state = "sealed",
            peer_blob = paste0("opaque_", server)),
          exactGCVecmulValidityReceiveDS = list(
            stored = TRUE, state = "checked"),
          exactGCVecmulCommitDS = list(
            stored = TRUE, state = "committed"),
          stop("unexpected guarded endpoint: ", command))
      }), names(conns))
    }

    result <- testthat::with_mocked_bindings(
      .dsvert_exact_gc_vecmul_run(
        datasources, server_names = server_names, servers = 1:2,
        session_id = paste0(
          "0000000", k, "-0000-4000-8000-00000000000", k),
        total_n = 2L, input_manifests = manifests,
        transport_ready = TRUE, .aggregate = aggregate),
      .dsvert_exact_gc_run = function(...) invisible(list()),
      .package = "dsVertClient")
    expect_identical(result$capability_id, "exact_gc_v1")
    expect_setequal(unique(commands), expected)
    expect_false(any(commands %in% c(
      "exactGCVecmulBindInputsDS", "exactGCGLMSoftplusPrepareDS",
      "mpcCleanupDS")))
  }
})
