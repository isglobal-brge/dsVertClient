 .chisq_cross_manifest_fixture <- function() {
  artifact <- list(
    version = .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION,
    analysis_id = "cross_table",
    left = list(dataset = "left_data", column = "left_cat",
                owner_peer = "site_b", levels = list("A", "B")),
    right = list(dataset = "right_data", column = "right_cat",
                 owner_peer = "site_c", levels = list("X", "Y", "Z")),
    participating_peers = list("site_b", "site_c"),
    computation_peers = list("site_a", "site_b"),
    coordinate_count = 6,
    repeated_record_policy = "consistent_level_else_zero_v1",
    missingness_policy = "missing_is_zero_v1",
    selected_l1_sensitivity = 1,
    selected_l2_sensitivity = 1,
    numeric_certificate = list(frac_bits = 8),
    transcript = list(padded_units = 2))
  manifest <- list(
    admission = list(unit_capacity = 2),
    workload = list(coordinate_count = 7, families = list(
      admitted_count = list(owner_peer = "site_a", dataset = "count_data"),
      numeric_moments = list(artifacts = list()),
      numeric_pair_moments = list(artifacts = list()),
      gaussian_models = list(artifacts = list()),
      fixed_numeric_histograms = list(artifacts = list()),
      categorical_marginals = list(artifacts = list()),
      categorical_pairs = list(
        sets = list(), cross_artifacts = list(cross_table = artifact)),
      correlation_artifacts = list(), describe_artifacts = list(),
      survival_artifacts = list())),
    capsule_identity = list(capsule_id = strrep("a", 64L)))
  conns <- stats::setNames(lapply(c("site_a", "site_b", "site_c"),
    function(peer) structure(list(peer = peer), class = "fake")),
    c("site_a", "site_b", "site_c"))
  list(manifest = manifest, conns = conns)
}

test_that("cross chi-square surface uses only the Synopsis contingency route", {
  surface <- paste(deparse(body(ds.vertChisqCross)), collapse = "\n")
  retired <- c(
    "dsvertColNamesDS", "dsvertOneHotDS", "exactGCChisq",
    "k2ChisqCross", "glmRing63TransportInitDS", "fisher.test")
  expect_false(any(vapply(
    retired, grepl, logical(1L), x = surface, fixed = TRUE)))
  expect_true(grepl("release <- ds.vertDPContingency", surface,
                    fixed = TRUE))
  expect_match(surface, ".dsvert_dp_chisq_from_release", fixed = TRUE)
  expect_match(surface, ".dsvert_dp_fisher_from_release", fixed = TRUE)
})

test_that("cross chi-square and Fisher reuse the identical DP release", {
  release <- structure(list(
    row_var = "row", col_var = "column", cross_owner = TRUE,
    servers = c("site_a", "site_b")),
    class = c("ds.vertDPContingency", "list"))
  chisq_seen <- fisher_seen <- NULL
  result <- testthat::with_mocked_bindings(
    ds.vertChisqCross(
      release, fisher = TRUE, verbose = FALSE,
      simulations = 10L, mc_confidence = 0.9),
    .dsvert_dp_chisq_from_release = function(x, ...) {
      chisq_seen <<- x
      structure(list(p_value = 0.5), class = c("ds.vertChisq", "list"))
    },
    .dsvert_dp_fisher_from_release = function(x, ...) {
      fisher_seen <<- x
      structure(list(p_value = 0.25), class = c("ds.vertFisher", "list"))
    },
    .package = "dsVertClient")
  expect_identical(chisq_seen, release)
  expect_identical(fisher_seen, release)
  expect_identical(result$fisher_p, 0.25)
  expect_identical(result$source_dp_release, release)
  expect_true(result$cross_owner)
})

test_that("new cross-owner categorical requests use one Synopsis release", {
  fixture <- .chisq_cross_manifest_fixture()
  calls <- 0L
  release <- structure(list(
    row_var = "left_cat", col_var = "right_cat", cross_owner = TRUE,
    servers = c("site_b", "site_c")),
    class = c("ds.vertDPContingency", "list"))
  result <- testthat::with_mocked_bindings(
    ds.vertChisqCross(
      "right_data", "left_cat", "right_cat",
      datasources = fixture$conns, verbose = FALSE, simulations = 10L),
    ds.vertDPContingency = function(
        data_name, row_var, col_var, server, datasources) {
      calls <<- calls + 1L
      expect_identical(data_name, "right_data")
      expect_identical(c(row_var, col_var), c("left_cat", "right_cat"))
      expect_null(server)
      expect_identical(datasources, fixture$conns)
      release
    }, .dsvert_dp_chisq_from_release = function(x, ...) {
      expect_identical(x, release)
      structure(list(p_value = 0.5), class = c("ds.vertChisq", "list"))
    }, .package = "dsVertClient")
  expect_identical(calls, 1L)
  expect_identical(result$source_dp_release, release)
  expect_true(result$cross_owner)
})
