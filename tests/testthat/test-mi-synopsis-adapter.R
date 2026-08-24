test_that("categorical MI completion is deterministic and preserves the signed total", {
  first <- dsVertClient:::.dsvert_mi_complete_counts_v1(
    observed_counts = c(control = 30, case = 10), admitted_count = 60,
    m = 8L, release_sha256 = paste(rep("a", 64L), collapse = ""))
  second <- dsVertClient:::.dsvert_mi_complete_counts_v1(
    observed_counts = c(control = 30, case = 10), admitted_count = 60,
    m = 8L, release_sha256 = paste(rep("a", 64L), collapse = ""))

  expect_identical(first, second)
  expect_identical(first$missing_count_dp, 20)
  expect_identical(first$completed_count_dp, 60)
  expect_identical(dim(first$completed_counts), c(8L, 2L))
  expect_true(all(rowSums(first$completed_counts) == 60))
  expect_true(all(first$completed_counts >= 0))
  expect_equal(colMeans(first$completed_counts) / 60,
               first$pooled_probabilities)

  changed_root <- dsVertClient:::.dsvert_mi_complete_counts_v1(
    observed_counts = c(control = 30, case = 10), admitted_count = 60,
    m = 8L, release_sha256 = paste(rep("b", 64L), collapse = ""))
  expect_false(identical(first$completed_counts, changed_root$completed_counts))
})

test_that("categorical MI completion accepts lattice counts and rejects invalid inputs", {
  root <- paste(rep("a", 64L), collapse = "")
  lattice <- dsVertClient:::.dsvert_mi_complete_counts_v1(
    observed_counts = c(control = 30.5, case = 10), admitted_count = 60.25,
    m = 8L, release_sha256 = root)
  expect_equal(rowSums(lattice$completed_counts),
               rep(60.25, 8L), tolerance = 1e-12)
  expect_error(dsVertClient:::.dsvert_mi_complete_counts_v1(
    observed_counts = c(control = -1, case = 10), admitted_count = 60,
    m = 8L, release_sha256 = root), "non-negative")
  expect_error(dsVertClient:::.dsvert_mi_complete_counts_v1(
    observed_counts = c(control = 30, case = 10), admitted_count = 60,
    m = 1L, release_sha256 = root), "m")
  expect_error(dsVertClient:::.dsvert_mi_complete_counts_v1(
    observed_counts = c(control = 30, case = 10), admitted_count = 60,
    m = 8L, release_sha256 = "not-a-root"), "release")
})

.mi_synopsis_adapter_fixture <- function() {
  root <- paste(rep("c", 64L), collapse = "")
  binding <- list(
    artifact_key = paste0("artifact_", root),
    execution_id = paste0("execution_", root),
    manifest_sha256 = root, contract_sha256 = root,
    attempt_sha256 = root, source_contract_sha256 = root,
    result_set_sha256 = root, final_vector_root = root,
    coordinate_order_sha256 = root, release_sha256 = root)
  state <- new.env(parent = emptyenv())
  state$runs <- 0L
  run <- function(...) {
    state$runs <- state$runs + 1L
    list(root = root)
  }
  count <- function(data_name, server, datasources, .aggregate, .run) {
    .run(datasources, .aggregate = .aggregate)
    c(binding, list(
      value = 60, source_owner = "peer_a", dataset = data_name))
  }
  frequency <- function(data_name, variable, source, datasources,
                        .aggregate, .run) {
    .run(datasources, .aggregate = .aggregate)
    c(binding, list(
      source_owner = source, variable = variable,
      levels = c("control", "case"), counts = c(control = 30, case = 10),
      missingness_policy =
        dsVertClient:::.dsvert_mi_strict_missingness_policy_v1,
      coordinate_descriptor = list(
        missingness_policy =
          dsVertClient:::.dsvert_mi_strict_missingness_policy_v1)))
  }
  list(state = state, run = run, count = count, frequency = frequency)
}

test_that("MI frontdoor consumes one strict signed Synopsis release", {
  fixture <- .mi_synopsis_adapter_fixture()
  fit <- dsVertClient:::.dsvert_mi_synopsis_result_v1(
    outcome ~ 1, "protected", NULL, 6L, "auto", list(peer_a = NULL),
    identity, .run = fixture$run, .count = fixture$count,
    .frequency = fixture$frequency)

  expect_s3_class(fit, "ds.vertMI")
  expect_identical(fixture$state$runs, 1L)
  expect_identical(fit$status, "ok")
  expect_identical(fit$family, "binomial")
  expect_identical(fit$missing_count_dp, 20)
  expect_identical(fit$completed_count_dp, 60)
  expect_true(is.finite(fit$coefficients[["(Intercept)"]]))
  expect_null(fit$standard_errors)
  expect_null(fit$covariance)
  expect_false("completed_counts" %in% names(fit))
  expect_identical(fit$additional_privacy_cost, c(epsilon = 0, delta = 0))
})

test_that("MI frontdoor rejects an unbound or non-strict marginal", {
  fixture <- .mi_synopsis_adapter_fixture()
  bad_frequency <- function(...) {
    value <- fixture$frequency(...)
    value$missingness_policy <- "missing_or_out_of_domain_rows_are_ignored"
    value$coordinate_descriptor$missingness_policy <- value$missingness_policy
    value
  }
  expect_error(dsVertClient:::.dsvert_mi_synopsis_result_v1(
    outcome ~ 1, "protected", NULL, 6L, "auto", list(peer_a = NULL),
    identity, .run = fixture$run, .count = fixture$count,
    .frequency = bad_frequency), "strict missingness")
})
