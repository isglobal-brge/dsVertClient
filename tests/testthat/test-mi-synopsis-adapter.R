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
    artifact_key = root, execution_id = root,
    manifest_sha256 = root, contract_sha256 = root,
    attempt_sha256 = root, source_contract_sha256 = root,
    result_set_sha256 = root, final_vector_root = root,
    coordinate_order_sha256 = root, release_sha256 = root)
  state <- new.env(parent = emptyenv())
  state$runs <- 0L
  state$contingencies <- 0L
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
    marginal <- switch(variable,
      outcome = list(levels = c("control", "case"),
                     counts = c(control = 30, case = 10)),
      exposure = list(levels = c("unexposed", "exposed", "unknown"),
                      counts = c(unexposed = 24, exposed = 12, unknown = 4)),
      region = list(levels = c("north", "south", "west"),
                    counts = c(north = 20, south = 12, west = 8)),
      stop("unknown fixture marginal", call. = FALSE))
    c(binding, list(
      source_owner = source, variable = variable,
      levels = marginal$levels, counts = marginal$counts,
      missingness_policy =
        dsVertClient:::.dsvert_mi_strict_missingness_policy_v1,
      coordinate_descriptor = list(
        missingness_policy =
          dsVertClient:::.dsvert_mi_strict_missingness_policy_v1)))
  }
  contingency <- function(data_name, row_var, col_var, server, datasources,
                          .aggregate) {
    state$contingencies <- state$contingencies + 1L
    if (!identical(data_name, "protected") || !identical(row_var, "outcome") ||
        !col_var %in% c("exposure", "region") || !is.null(server)) {
      stop("unexpected fixture joint pair", call. = FALSE)
    }
    table <- switch(col_var,
      exposure = matrix(c(18, 6, 8, 2, 5, 1), nrow = 2L,
                        dimnames = list(c("control", "case"),
                                        c("unexposed", "exposed", "unknown"))),
      region = matrix(c(14, 6, 9, 3, 7, 1), nrow = 2L,
                      dimnames = list(c("control", "case"),
                                      c("north", "south", "west"))))
    c(binding, list(
      row_var = row_var, col_var = col_var,
      missingness_policy =
        dsVertClient:::.dsvert_mi_strict_joint_missingness_policy_v1,
      admitted_count_dp = 60,
      table = table))
  }
  list(state = state, run = run, count = count, frequency = frequency,
       contingency = contingency)
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

test_that("MI completes one strict categorical pair jointly", {
  fixture <- .mi_synopsis_adapter_fixture()
  first <- dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure) ~ 1, "protected", NULL, 6L, "auto",
    list(peer_a = NULL), identity, .run = fixture$run,
    .count = fixture$count, .frequency = fixture$frequency,
    .contingency = fixture$contingency)
  second <- dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure) ~ 1, "protected", NULL, 6L, "auto",
    list(peer_a = NULL), identity, .run = fixture$run,
    .count = fixture$count, .frequency = fixture$frequency,
    .contingency = fixture$contingency)

  expect_s3_class(first, "ds.vertMI")
  expect_identical(first$method, "signed_categorical_mcar_joint_pair_v3")
  expect_identical(first$joint_model,
                   "strict_missing_signed_joint_pair_completion_v1")
  expect_identical(names(first$variables), c("outcome", "exposure"))
  expect_true(all(is.finite(first$variables$outcome$probabilities)))
  expect_true(all(is.finite(first$variables$exposure$probabilities)))
  expect_equal(sum(first$variables$outcome$probabilities), 1, tolerance = 1e-12)
  expect_equal(sum(first$variables$exposure$probabilities), 1, tolerance = 1e-12)
  expect_equal(sum(first$joint_probabilities), 1, tolerance = 1e-12)
  expect_identical(first$completed_draws_sha256, second$completed_draws_sha256)
  expect_identical(first$observed_joint_table_dp, second$observed_joint_table_dp)
  expect_identical(first$additional_privacy_cost, c(epsilon = 0, delta = 0))
  expect_identical(fixture$state$runs, 0L)
  expect_identical(fixture$state$contingencies, 2L)

  bad_pair <- function(...) {
    value <- fixture$contingency(...)
    value$missingness_policy <- "missing_or_out_of_domain_rows_are_ignored"
    value
  }
  expect_error(dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure) ~ 1, "protected", NULL, 6L, "auto",
    list(peer_a = NULL), identity, .contingency = bad_pair),
    "strict missingness")
  expect_error(dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure) ~ 1, "protected", c("exposure", "outcome"),
    6L, "auto", list(peer_a = NULL), identity, .run = fixture$run,
    .count = fixture$count, .frequency = fixture$frequency,
    .contingency = fixture$contingency),
    "formula order")
  expect_error(dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure) ~ 1, "protected", NULL, 6L, "binomial",
    list(peer_a = NULL), identity, .run = fixture$run,
    .count = fixture$count, .frequency = fixture$frequency,
    .contingency = fixture$contingency),
    "family = 'auto'")
})

test_that("MI completes three categorical marginals without claiming a joint model", {
  fixture <- .mi_synopsis_adapter_fixture()
  fit <- dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure, region) ~ 1, "protected", NULL, 6L, "auto",
    list(peer_a = NULL), identity, .run = fixture$run,
    .count = fixture$count, .frequency = fixture$frequency)
  expect_identical(fit$method,
                   "signed_categorical_mcar_independent_marginals_v2")
  expect_identical(names(fit$variables), c("outcome", "exposure", "region"))
  expect_identical(fixture$state$runs, 1L)
  expect_error(dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure, region) ~ 1, "protected",
    c("exposure", "outcome", "region"), 6L, "auto",
    list(peer_a = NULL), identity, .run = fixture$run,
    .count = fixture$count, .frequency = fixture$frequency), "formula order")
})

test_that("MI builds a sticky multivariable categorical star model from signed pairs", {
  fixture <- .mi_synopsis_adapter_fixture()
  star_pair <- function(run, data_name, row_var, col_var) {
    fixture$contingency(
      data_name, row_var, col_var, NULL, list(peer_a = NULL), identity)
  }
  first <- dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure, region) ~ 1, "protected", NULL, 6L, "auto",
    list(peer_a = NULL), identity, .run = fixture$run,
    .count = fixture$count, .frequency = fixture$frequency,
    .star_pair = star_pair, dependence = "star")
  second <- dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure, region) ~ 1, "protected", NULL, 6L, "auto",
    list(peer_a = NULL), identity, .run = fixture$run,
    .count = fixture$count, .frequency = fixture$frequency,
    .star_pair = star_pair, dependence = "star")

  expect_s3_class(first, "ds.vertMI")
  expect_identical(first$method, "signed_categorical_mcar_star_joint_v1")
  expect_identical(first$joint_model,
                   "strict_missing_signed_pairwise_star_completion_v1")
  expect_identical(first$root_column, "outcome")
  expect_identical(names(first$conditional_probabilities), c("exposure", "region"))
  expect_equal(sum(first$root_probabilities), 1, tolerance = 1e-12)
  expect_true(all(vapply(first$conditional_probabilities, function(value) {
    all(is.finite(value)) && all(value >= 0) &&
      all(abs(rowSums(value) - 1) < 1e-12)
  }, logical(1L))))
  expect_identical(first$completed_draws_sha256, second$completed_draws_sha256)
  expect_identical(fixture$state$runs, 2L)
  expect_identical(fixture$state$contingencies, 4L)

  inconsistent <- function(...) {
    release <- fixture$contingency(...)
    if (identical(release$col_var, "region")) {
      release$final_vector_root <- strrep("d", 64L)
    }
    release
  }
  expect_error(dsVertClient:::.dsvert_mi_synopsis_result_v1(
    cbind(outcome, exposure, region) ~ 1, "protected", NULL, 6L, "auto",
    list(peer_a = NULL), identity, .run = fixture$run,
    .star_pair = function(run, data_name, row_var, col_var) {
      inconsistent(
        data_name, row_var, col_var, NULL, list(peer_a = NULL), identity)
    },
    dependence = "star"), "same signed vector release")
})

test_that("the MI aliases forward one multivariable categorical request", {
  seen <- NULL
  response <- structure(list(
    status = "ok",
    method = "signed_categorical_mcar_independent_marginals_v2",
    variables = list(), source_values_exposed = FALSE,
    intermediate_values_exposed = FALSE, production_ready = FALSE),
    class = c("ds.vertMI", "list"))
  testthat::local_mocked_bindings(
    .dsvert_datasources = function(value) value,
    .dsvert_federation_argument = function(value, datasources) {
      list(value = value, datasources = datasources)
    },
    .dsvert_mi_synopsis_result_v1 = function(
        formula, data_name, impute_columns, m, family, datasources, .aggregate,
        dependence = c("independent", "star")) {
      seen <<- list(
        formula = formula, data_name = data_name,
        impute_columns = impute_columns, m = m, family = family,
        dependence = dependence, datasources = datasources)
      response
    },
    .package = "dsVertClient")

  fit <- ds.vert.mi(
    cbind(outcome, exposure) ~ 1, data = "protected",
    impute_columns = c("outcome", "exposure"), m = 6L,
    family = "auto", datasources = list(site_a = structure(list(), class = "mock")))

  expect_s3_class(fit, "ds.vertMI")
  expect_identical(fit$frontdoor, "ds.vert.mi")
  expect_identical(seen$data_name, "protected")
  expect_identical(seen$impute_columns, c("outcome", "exposure"))
  expect_identical(seen$m, 6L)
  expect_identical(seen$family, "auto")
  expect_identical(seen$dependence, "independent")
  expect_identical(seen$formula, cbind(outcome, exposure) ~ 1)
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
