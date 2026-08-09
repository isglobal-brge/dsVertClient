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

test_that("cross chi-square surface contains no retired exact/discovery route", {
  surface <- paste(deparse(body(ds.vertChisqCross)), collapse = "\n")
  retired <- c(
    "dsvertColNamesDS", "dsvertOneHotDS", "exactGCChisq",
    "k2ChisqCross", "glmRing63TransportInitDS", "fisher.test")
  expect_false(any(vapply(
    retired, grepl, logical(1L), x = surface, fixed = TRUE)))
  expect_match(surface, "ds.vertDPContingency", fixed = TRUE)
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

test_that("cross contingency maps the signed orientation without discovery", {
  fixture <- .chisq_cross_manifest_fixture()
  layout <- .dsvert_dp_capsule_vector_layout(fixture$manifest)
  context <- list(
    layout = layout, manifest = fixture$manifest,
    release = structure(list(), class = "test-release"),
    adjacency = "add_remove_patient")
  values <- c(1, 2, 0, 1, 2, 0)
  run <- function(data_name, row_var, col_var, server = NULL) {
    testthat::with_mocked_bindings(
      .dsvert_dp_contingency_impl(
        data_name, row_var, col_var, server = server,
        datasources = fixture$conns,
        .aggregate = function(...) stop("no DSI call expected")),
      .dsvert_dp_datasources = function(datasources) datasources,
      .dsvert_dp_capsule_vector_run = function(...) list(),
      .dsvert_dp_vector_context = function(...) context,
      .dsvert_dp_capsule_vector_values = function(release, block) values,
      .dsvert_dp_vector_public_metadata = function(...) list(
        mechanism = "discrete-laplace"),
      .dsvert_dp_vector_accuracy_radius = function(...) list(
        radius = 1, confidence = 0.95, method = "test",
        implementation_tv_upper_bound = 0,
        additional_privacy_cost = 0),
      .package = "dsVertClient")
  }
  canonical <- run("right_data", "left_cat", "right_cat", "site_b")
  expect_true(canonical$cross_owner)
  expect_identical(canonical$servers, c("site_b", "site_c"))
  expect_identical(canonical$datasets, c("left_data", "right_data"))
  expect_equal(unname(canonical$table), matrix(values, nrow = 2L))

  transposed <- run("left_data", "right_cat", "left_cat", "site_c")
  expect_equal(unname(transposed$table), t(matrix(values, nrow = 2L)))
  expect_error(run("unknown", "left_cat", "right_cat"),
               "does not contain exactly one")
})
